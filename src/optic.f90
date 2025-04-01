subroutine linear_optic
    use wmpi
    use para
    implicit none
    
   !------------------------------------------------------------------!
   !> This subroutine is to calculate the linear optical conductivity !
   !>                                                                 !
   !>                                                                 !
   !> Use TBA approximation                                           !
   !> Ref. : PHYSICAL REVIEW B 97, 245143 (2018)                      !
   !------------------------------------------------------------------!
    
    integer :: ik, ikx, iky, ikz, knv3, ifreq, i, j, m, n, index, ierr

    real(dp) :: k(3), time_start, time_end, x, fac_H, fac_AH
    complex(dp) :: cmplx_i, cmplx_1, cmplx_0, omega

    complex(dp), allocatable :: Freq_array(:)

    !> Hermitian and anti-Hermitian part of conductivity tensor
    !> dimention: 3*3*FreqNum
    complex(dp), allocatable :: sigma_kubo_H(:, :, :)
    complex(dp), allocatable :: sigma_kubo_H_mpi(:, :, :)
    complex(dp), allocatable :: sigma_kubo_AH(:, :, :)
    complex(dp), allocatable :: sigma_kubo_AH_mpi(:, :, :)

    !> eigen value of H
    real(dp), allocatable :: W(:)
    complex(dp), allocatable :: Hamk_bulk(:, :)
    complex(dp), allocatable :: UU(:, :)
    real(dp), allocatable :: occ(:)

    !> V^ham_mn = UU_mj*dHdk*UU_jn
    complex(dp), allocatable :: V_ham(:, :, :)
    complex(dp), allocatable :: dHdk(:, :, :)

    !> D^ham_mn=V^ham_mn/(En-Em) for m!=n
    !> D_nn=0 
    complex(dp), allocatable :: D_ham(:, :, :)

    !> AA_mn = i*D_mn
    complex(dp), allocatable :: AA(:, :, :)

    !> delta function
    real(dp), external :: delta

    !> Re(sigma_S_xx), Im(sigma_S_xx), Re(sigma_S_yy), Im(sigma_S_yy), ...
    !> Re(sigma_A_yz), Im(sigma_A_yz), Re(sigma_A_zx), Im(sigma_A_zx), ...
    real(dp) :: sigma_S(12), sigma_A(6)
    integer :: alpha_S(6), beta_S(6), alpha_A(3), beta_A(3)

    !> 1 <--> xx
    !> 2 <--> yy
    !> 3 <--> zz
    !> 4 <--> xy
    !> 5 <--> xz
    !> 6 <--> yz
    alpha_S = (/1, 2, 3, 1, 1, 2/)
    beta_S = (/1, 2, 3, 2, 3, 3/)
    !> 1 <--> (y,z)
    !> 2 <--> (z,x)
    !> 3 <--> (x,y)
    alpha_A = (/2, 3, 1/)
    beta_A = (/3, 1, 2/)

    cmplx_i = (0.0d0, 1.0d0)
    cmplx_1 = (1.0d0, 0.0d0)
    cmplx_0 = (0.0d0, 0.0d0)

    allocate( W (Num_wann))
    allocate( Hamk_bulk(Num_wann, Num_wann))
    allocate( UU(Num_wann, Num_wann))
    allocate( occ(Num_wann))

    allocate(D_ham(Num_wann, Num_wann, 3))
    allocate(V_ham(Num_wann, Num_wann, 3))
    allocate(dHdk(Num_wann, Num_wann, 3))
    allocate(AA(Num_wann, Num_wann, 3))

    allocate(sigma_kubo_H(3, 3, FreqNum))
    allocate(sigma_kubo_H_mpi(3, 3, FreqNum))
    allocate(sigma_kubo_AH(3, 3, FreqNum))
    allocate(sigma_kubo_AH_mpi(3, 3, FreqNum))

    allocate(Freq_array(FreqNum))

    W = 0.0_dp
    Hamk_bulk = cmplx_0
    UU = cmplx_0
    AA = cmplx_0
    sigma_kubo_H = cmplx_0
    sigma_kubo_H_mpi = cmplx_0
    sigma_kubo_AH = cmplx_0
    sigma_kubo_AH_mpi = cmplx_0
    Freq_array = 0.0_dp

    if (FreqNum==1) then
        Freq_array(1)= FreqMin + cmplx_i*eta_smr_fixed
    else
        do i = 1, FreqNum
            Freq_array(i)= FreqMin+ (FreqMax-FreqMin)* (i-1d0)/dble(FreqNum-1) + cmplx_i*eta_smr_fixed
        enddo ! i
    endif

    knv3= Nk1*Nk2*Nk3

    call now(time_start) 
    do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0.and. mod(ik/num_cpu, 100).eq.0) then
            call now(time_end) 
            write(stdout, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
            ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/100d0
            time_start= time_end
        endif

        ikx= (ik-1)/(nk2*nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
        ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
        k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
            + K3D_vec2_cube*(iky-1)/dble(nk2)  &
            + K3D_vec3_cube*(ikz-1)/dble(nk3)

        ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
        call ham_bulk_atomicgauge(k, Hamk_bulk)
        ! call ham_bulk_latticegauge(k, Hamk_bulk)

        !> diagonalization by call zheev in lapack
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)

        occ = 0.0_dp
        do i = 1, Num_wann
            if (W(i) < 0.0_dp) occ(i) = 1.0_dp
        end do

        !> get velocity operator in Hamiltonian basis
        call dHdk_atomicgauge_Ham(k, UU, V_Ham)
        ! dHdk = cmplx_0
        ! call dHdk_latticegauge_wann(k, dHdk)
        ! call dHdk_latticegauge_Ham2(UU, dHdk, V_Ham)

        call get_Dmn_Ham(W, V_Ham, D_Ham)

        !> here we use TBA approximation 
        AA = cmplx_i*D_Ham

        do ifreq = 1, FreqNum
            omega = Freq_array(ifreq)
            do m = 1, Num_wann
                do n = 1, Num_wann
                    if (n == m) cycle
                    x = W(m)-W(n)-real(omega)
                    fac_H = delta(eta_smr_fixed, x)*(occ(m)-occ(n))*(W(m)-W(n))
                    fac_AH = (occ(m)-occ(n))*(W(m)-W(n))/real((W(m)-W(n)-omega))
                    do j = 1, 3
                        do i = 1, 3
                            sigma_kubo_H_mpi(i, j, ifreq) = sigma_kubo_H_mpi(i, j, ifreq) &
                                                      + AA(n, m, i)*AA(m, n, j)* fac_H
                            sigma_kubo_AH_mpi(i, j, ifreq) = sigma_kubo_AH_mpi(i, j, ifreq) &
                                                       + AA(n, m, i)*AA(m, n, j)*fac_AH
                        enddo 
                    enddo 
                enddo 
            enddo 
        enddo 

    enddo ! ik

#if defined (MPI)
    call mpi_allreduce(sigma_kubo_H_mpi, sigma_kubo_H, size(sigma_kubo_H), mpi_dc,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(sigma_kubo_AH_mpi, sigma_kubo_AH, size(sigma_kubo_AH), mpi_dc,mpi_sum,mpi_cmw,ierr)
#else
    sigma_kubo_H = sigma_kubo_H_mpi
    sigma_kubo_AH = sigma_kubo_AH_mpi
#endif

   ! ------------------------------------------------------------------------
   ! At this point 
   !
   ! sigma_kubo_H =N*V_c*int dk/(2pi)^3 (f_m-f_n)*(E_m-E_n)*A_nm*A_mn*delta(Em-En-omega)
   ! sigma_kubo_AH=N*V_c*int dk/(2pi)^3 (f_m-f_n)*Re[(E_m-E_n)/(E_m-E_n-omega-i*eta)]*A_nm*A_mn
   !
   ! (N is the number of kpoints, V_c is the cell volume). We want
   !
   ! sigma_kubo_H =-pi*e^2/hbar* int dk/(2pi)^3 (f_m-f_n)*(E_m-E_n)*A_nm*A_mn*delta(Em-En-omega)
   ! sigma_kubo_AH=i*e^2/hbar* int dk/(2pi)^3 (f_m-f_n)*Re[(E_m-E_n)/(E_m-E_n-omega-i*eta)]*A_nm*A_mn
   ! 
   ! --------------------------------------------------------------------

      
    !> in the latest version, we use the atomic unit
    sigma_kubo_H = -sigma_kubo_H/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume&
                  *pi*Echarge**2/hbar/Bohr_radius/100
    sigma_kubo_AH = sigma_kubo_AH/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume&
                   *Echarge**2/hbar*cmplx_i/Bohr_radius/100

    if (cpuid.eq.0) then
        outfileindex= outfileindex+ 1
        open(unit=outfileindex, file='sigma_kubo_symm.dat')
        write(outfileindex, '("#",10a)')' the symmetric part of linear conductivity'
        write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 13)
        write(outfileindex, '("#",a13, 20a16)')'Frequency (eV)', 'Re[\sigma_xx]', 'Im[\sigma_xx]', 'Re[\sigma_yy]', &
                            'Im[\sigma_yy]', 'Re[\sigma_zz]', 'Im[\sigma_zz]', 'Re[\sigma_xy]', 'Im[\sigma_xy]', &
                            'Im[\sigma_xz]', 'Re[\sigma_xz]', 'Im[\sigma_yz]', 'Re[\sigma_yz]'
        do ifreq=1, FreqNum
            do index = 1, 6
                i = alpha_S(index)
                j = beta_S(index)
                sigma_S(index*2-1) = real(0.5_dp*(sigma_kubo_H(i, j, ifreq) + sigma_kubo_H(j, i, ifreq)), dp)
                sigma_S(index*2)   = aimag(0.5_dp*(sigma_kubo_AH(i, j, ifreq) + sigma_kubo_AH(j, i, ifreq)))
            enddo
            write(outfileindex, '(200E16.8)') real(Freq_array(ifreq), dp)/eV2Hartree, sigma_S
        enddo 
        close(outfileindex)

        outfileindex= outfileindex+ 1
        open(unit=outfileindex, file='sigma_kubo_asymm.dat')
        write(outfileindex, '("#",10a)')' the antisymmetric part of linear conductivity'
        write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 7)
        write(outfileindex, '("#",a13, 20a16)')'Frequency (eV)', 'Re[\sigma_yz]', 'Im[\sigma_yz]', 'Re[\sigma_zx]', &
        'Im[\sigma_zx]', 'Re[\sigma_xy]', 'Im[\sigma_xy]'
        do ifreq=1, FreqNum
            do index = 1, 3
                i = alpha_A(index)
                j = beta_A(index)
                sigma_A(index*2-1) = real(0.5_dp*(sigma_kubo_AH(i, j, ifreq) - sigma_kubo_AH(j, i, ifreq)), dp)
                sigma_A(index*2)   = aimag(0.5_dp*(sigma_kubo_H(i, j, ifreq) - sigma_kubo_H(j, i, ifreq)))
            enddo
            write(outfileindex, '(200E16.8)') real(Freq_array(ifreq), dp)/eV2Hartree, sigma_A
        enddo ! ie
        close(outfileindex)
    endif

end subroutine linear_optic

subroutine bulk_photovoltaic
    use wmpi
    use para
    implicit none
   !------------------------------------------------------------------!
   !> This subroutine is to calculate the bulk photovotaic effect     !
   !> including shift current and injection current                   !
   !>                                                                 !
   !>                                                                 !
   !> Use TBA approximation                                           !
   !>                                                                 !
   !------------------------------------------------------------------!

    integer :: ik, ikx, iky, ikz, knv3, ifreq, i, j, a, b, c, m, n, index, ierr

    real(dp) :: k(3), time_start, time_end, xplus, xminus
    complex(dp) :: cmplx_i, cmplx_1, cmplx_0, omega

    complex(dp), allocatable :: Freq_array(:)

    !> dimention: 3*6*FreqNum
    real(dp), allocatable :: lshiftcur(:, :, :)
    real(dp), allocatable :: cshiftcur(:, :, :)
    real(dp), allocatable :: linjectcur(:, :, :)
    real(dp), allocatable :: cinjectcur(:, :, :)
    real(dp), allocatable :: lshiftcur_mpi(:, :, :)
    real(dp), allocatable :: cshiftcur_mpi(:, :, :)
    real(dp), allocatable :: linjectcur_mpi(:, :, :)
    real(dp), allocatable :: cinjectcur_mpi(:, :, :)
    

    !> eigen value of H
    real(dp), allocatable :: W(:)
    complex(dp), allocatable :: Hamk_bulk(:, :)
    complex(dp), allocatable :: UU(:, :)
    real(dp), allocatable :: occ(:)

    !> V^ham_mn = UU_mj*dHdk*UU_jn
    complex(dp), allocatable :: V_ham(:, :, :)
    complex(dp), allocatable :: dHdk(:, :, :)
    complex(dp), allocatable :: dHdkdk(:, :, :, :)

    !> D^ham_mn=V^ham_mn/(En-Em) for m!=n
    !> D_nn=0 
    complex(dp), allocatable :: D_ham(:, :, :)

    !> AA_mn = i*D_mn
    complex(dp), allocatable :: AA(:, :, :)
    !> Wmn_ham = U^dag*dHdkdk*U
    complex(dp), allocatable :: Wmn_ham(:, :, :, :)
    !> general derivative
    complex(dp), allocatable :: gen_der_r(:, :, :, :)

    !> delta function
    real(dp), external :: delta

    !> Re(sigma_S_xx), Im(sigma_S_xx), Re(sigma_S_yy), Im(sigma_S_yy), ...
    !> Re(sigma_A_yz), Im(sigma_A_yz), Re(sigma_A_zx), Im(sigma_A_zx), ...
    real(dp) :: sigma_S(12), sigma_A(6)
    integer :: alpha_S(6), beta_S(6), alpha_A(3), beta_A(3)

    !> 1 <--> xx
    !> 2 <--> xy
    !> 3 <--> xz
    !> 4 <--> yy
    !> 5 <--> yz
    !> 6 <--> zz
    alpha_S = (/1, 1, 1, 2, 2, 3/)
    beta_S = (/1, 2, 3, 2, 3, 3/)

    cmplx_i = (0.0d0, 1.0d0)
    cmplx_1 = (1.0d0, 0.0d0)
    cmplx_0 = (0.0d0, 0.0d0)

    allocate( W (Num_wann))
    allocate( Hamk_bulk(Num_wann, Num_wann))
    allocate( UU(Num_wann, Num_wann))
    allocate( occ(Num_wann))

    allocate(D_ham(Num_wann, Num_wann, 3))
    allocate(V_ham(Num_wann, Num_wann, 3))
    allocate(dHdk(Num_wann, Num_wann, 3))
    allocate(AA(Num_wann, Num_wann, 3))
    allocate(dHdkdk(Num_wann, Num_wann, 3, 3))
    allocate(Wmn_ham(Num_wann, Num_wann, 3, 3))
    allocate(gen_der_r(Num_wann, Num_wann, 3, 3))

    allocate(lshiftcur(3, 6, FreqNum))
    allocate(linjectcur(3, 6, FreqNum))
    allocate(cshiftcur(3, 6, FreqNum))
    allocate(cinjectcur(3, 6, FreqNum))
    allocate(lshiftcur_mpi(3, 6, FreqNum))
    allocate(cshiftcur_mpi(3, 6, FreqNum))
    allocate(linjectcur_mpi(3, 6, FreqNum))
    allocate(cinjectcur_mpi(3, 6, FreqNum))

    allocate(Freq_array(FreqNum))

    W = 0.0_dp
    Hamk_bulk = cmplx_0
    UU = cmplx_0
    AA = cmplx_0
    lshiftcur = 0.0_dp
    cshiftcur = 0.0_dp
    lshiftcur_mpi = 0.0_dp
    cshiftcur_mpi = 0.0_dp
    linjectcur = 0.0_dp
    cinjectcur = 0.0_dp
    linjectcur_mpi = 0.0_dp
    cinjectcur_mpi = 0.0_dp
    Freq_array = 0.0_dp

    if (FreqNum==1) then
        Freq_array(1)= FreqMin + cmplx_i*eta_smr_fixed
    else
        do i = 1, FreqNum
            Freq_array(i)= FreqMin+ (FreqMax-FreqMin)* (i-1d0)/dble(FreqNum-1) + cmplx_i*eta_smr_fixed
        enddo ! i
    endif

    knv3= Nk1*Nk2*Nk3

    call now(time_start) 
    do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0.and. mod(ik/num_cpu, 100).eq.0) then
            call now(time_end) 
            write(stdout, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
            ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/100d0
            time_start= time_end
        endif

        ikx= (ik-1)/(nk2*nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
        ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
        k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
            + K3D_vec2_cube*(iky-1)/dble(nk2)  &
            + K3D_vec3_cube*(ikz-1)/dble(nk3)

        ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
        call ham_bulk_atomicgauge(k, Hamk_bulk)
        ! call ham_bulk_latticegauge(k, Hamk_bulk)

        !> diagonalization by call zheev in lapack
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)

        occ = 0.0_dp
        do m = 1, Num_wann
            if (W(m) < 0.0_dp) occ(m) = 1.0_dp
        end do
        
        !> get velocity operator in Hamiltonian basis
        ! call dHdk_atomicgauge_Ham(k, UU, V_Ham)
        dHdk = cmplx_0
        ! call dHdk_latticegauge_wann(k, dHdk)
        ! call dHdk_latticegauge_Ham2(UU, dHdk, V_Ham)
        call dHdk_atomicgauge_Ham(k, UU, V_Ham)
        call get_Dmn_Ham(W, V_Ham, D_Ham)
        !> here we use TBA approximation 
        AA = cmplx_i*D_Ham

        !> get d^2H/dk^2 
        dHdkdk = cmplx_0
        Wmn_ham = cmplx_0
        ! call dHdkdk_latticegauge_wann(k, dHdkdk)
        ! !> eq(29c)
        ! call dHdkdk_latticegauge_Ham(UU, dHdkdk, Wmn_Ham)

        !YANYI ATTENTION!!!!!!! keep wann, no wann is a test!!!
        !call d2Hdk2_atomicgauge_wann(k, dHdkdk) 
        call d2Hdk2_atomicgauge(k, dHdkdk)

        !YANYI ATTENTION!!!!!!! keep Ham, no Ham is a test!!!
        !call d2Hdk2_atomicgauge_Ham(UU, dHdkdk, Wmn_ham)
        call d2Hdk2_atomicgauge(UU, dHdkdk, Wmn_ham)

        call generalderivative(W, V_Ham, D_Ham, Wmn_Ham, gen_der_r)

        do ifreq = 1, FreqNum
            omega = Freq_array(ifreq)
            do m = 1, Num_wann
                do n = 1, Num_wann
                    if (n == m) cycle
                    if (abs(occ(m)-occ(n)) < eps9) cycle
                    xminus = W(m)-W(n)-real(omega)
                    xplus = W(m)-W(n)+real(omega)
                    do a = 1, 3
                        do i = 1, 6
                            b = alpha_S(i)
                            c = beta_S(i)
                            lshiftcur_mpi(a, i, ifreq) = lshiftcur_mpi(a, i, ifreq) &
                                                        + (occ(m)-occ(n))*aimag((gen_der_r(m, n, c, a)*AA(n, m, b) &
                                                        + AA(n, m, c)*gen_der_r(m, n, b, a))*(delta(eta_smr_fixed, xminus) &
                                                        + delta(eta_smr_fixed, xplus)))
                            ! lshiftcur_mpi(a, i, ifreq) = lshiftcur_mpi(a, i, ifreq) &
                            !                             + (occ(m)-occ(n))*aimag((gen_der_r(m, n, c, a) &
                            !                             + gen_der_r(m, n, b, a))*(delta(eta_smr_fixed, xminus) &
                            !                             + delta(eta_smr_fixed, xplus)))
                            cshiftcur_mpi(a, i, ifreq) = cshiftcur_mpi(a, i, ifreq) &
                                                        + (occ(m)-occ(n))*real((gen_der_r(m, n, c, a)*AA(n, m, b) &
                                                        - AA(n, m, c)*gen_der_r(m, n, b, a))*(delta(eta_smr_fixed, xminus) &
                                                        - delta(eta_smr_fixed, xplus)))
                            linjectcur_mpi(a, i, ifreq) = linjectcur_mpi(a, i, ifreq) &
                                                        + (occ(m)-occ(n))*real((V_ham(m, m, a)-V_ham(n, n, a)) &
                                                        * (AA(n, m, b)*AA(m, n, c)+AA(n, m, c)*AA(m, n, b)) &
                                                        * delta(eta_smr_fixed, xminus))
                            cinjectcur_mpi(a, i, ifreq) = cinjectcur_mpi(a, i, ifreq) &
                                                        + (occ(m)-occ(n))*aimag((V_ham(m, m, a)-V_ham(n, n, a)) &
                                                        * (AA(n, m, b)*AA(m, n, c)-AA(n, m, c)*AA(m, n, b)) &
                                                        * delta(eta_smr_fixed, xminus))
                        enddo
                    enddo 
                enddo
            enddo  ! n
        enddo ! ifreq
    enddo ! ik

#if defined (MPI)
    call mpi_allreduce(lshiftcur_mpi, lshiftcur, size(lshiftcur), mpi_dp,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(cshiftcur_mpi, cshiftcur, size(cshiftcur), mpi_dp,mpi_sum,mpi_cmw,ierr) 
    call mpi_allreduce(linjectcur_mpi, linjectcur, size(linjectcur), mpi_dp,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(cinjectcur_mpi, cinjectcur, size(cinjectcur), mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
    lshiftcur = lshiftcur_mpi
    cshiftcur = cshiftcur_mpi
    linjectcur = linjectcur_mpi
    cinjectcur = cinjectcur_mpi
#endif

   ! ------------------------------------------------------------------------
   ! At this point 
   !
   !
   ! --------------------------------------------------------------------
      
    !> in the latest version, we use the atomic unit
    

    if (cpuid.eq.0) then
        outfileindex= outfileindex+ 1
        open(unit=outfileindex, file='linear_shift.dat')
        write(outfileindex, '("#",10a)')' the linear shift conductivity'
        write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 19)
        write(outfileindex, '("#",a13, 20a16)')'Frequency (eV)', 'xxx', 'xxy', 'xxz', 'xyy', 'xyz', 'xzz', &
                                                'yxx', 'yxy', 'yxz', 'yyy', 'yyz', 'yzz', 'zxx', 'zxy', 'zxz', &
                                                'zyy', 'zyz', 'zzz'
        do ifreq=1, FreqNum
            write(outfileindex, '(200E16.8)') real(Freq_array(ifreq), dp)/eV2Hartree, lshiftcur(1, :, ifreq),&
                                              lshiftcur(2, :, ifreq), lshiftcur(3, :, ifreq)
        enddo 
        close(outfileindex)

        outfileindex= outfileindex+ 1
        open(unit=outfileindex, file='circular_shift.dat')
        write(outfileindex, '("#",10a)')' the circular shift conductivity'
        write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 19)
        write(outfileindex, '("#",a13, 20a16)')'Frequency (eV)', 'xxx', 'xxy', 'xxz', 'xyy', 'xyz', 'xzz', &
                                                'yxx', 'yxy', 'yxz', 'yyy', 'yyz', 'yzz', 'zxx', 'zxy', 'zxz', &
                                                'zyy', 'zyz', 'zzz'
        do ifreq=1, FreqNum
            write(outfileindex, '(200E16.8)') real(Freq_array(ifreq), dp)/eV2Hartree, cshiftcur(1, :, ifreq),&
                                              cshiftcur(2, :, ifreq), cshiftcur(3, :, ifreq)
        enddo 
        close(outfileindex)

        outfileindex= outfileindex+ 1
        open(unit=outfileindex, file='circular_inject.dat')
        write(outfileindex, '("#",10a)')' the circular inject conductivity'
        write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 19)
        write(outfileindex, '("#",a13, 20a16)')'Frequency (eV)', 'xxx', 'xxy', 'xxz', 'xyy', 'xyz', 'xzz', &
                                                'yxx', 'yxy', 'yxz', 'yyy', 'yyz', 'yzz', 'zxx', 'zxy', 'zxz', &
                                                'zyy', 'zyz', 'zzz'
        do ifreq=1, FreqNum
            write(outfileindex, '(200E16.8)') real(Freq_array(ifreq), dp)/eV2Hartree, cinjectcur(1, :, ifreq),&
                                              cinjectcur(2, :, ifreq), cinjectcur(3, :, ifreq)
        enddo 
        close(outfileindex)

        outfileindex= outfileindex+ 1
        open(unit=outfileindex, file='linear_inject.dat')
        write(outfileindex, '("#",10a)')' the circular shift conductivity'
        write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 19)
        write(outfileindex, '("#",a13, 20a16)')'Frequency (eV)', 'xxx', 'xxy', 'xxz', 'xyy', 'xyz', 'xzz', &
                                                'yxx', 'yxy', 'yxz', 'yyy', 'yyz', 'yzz', 'zxx', 'zxy', 'zxz', &
                                                'zyy', 'zyz', 'zzz'
        do ifreq=1, FreqNum
            write(outfileindex, '(200E16.8)') real(Freq_array(ifreq), dp)/eV2Hartree, linjectcur(1, :, ifreq),&
                                              linjectcur(2, :, ifreq), linjectcur(3, :, ifreq)
        enddo 
        close(outfileindex)
    endif

end subroutine bulk_photovoltaic

subroutine generalderivative(W, V_Ham, D_Ham,  Wmn_Ham, gen_der_r)
    use para, only :Num_wann, zi, pi2zi, dp, stdout
    use wmpi
    implicit none

    real(dp), intent(in) :: W(Num_wann)
    complex(dp), intent(in) :: V_Ham(Num_wann, Num_wann, 3)
    complex(dp), intent(in) :: D_Ham(Num_wann, Num_wann, 3)
    complex(dp), intent(in) :: Wmn_Ham(Num_wann, Num_wann, 3, 3)
    complex(dp), intent(out) :: gen_der_r(Num_wann, Num_wann, 3, 3)

    integer :: m, n, i, j, p
    complex(dp) :: cmplx_i, cmplx_0

    cmplx_i = (0.0d0, 1.0d0)
    cmplx_0 = (0.0d0, 0.0d0)
    gen_der_r = cmplx_0
    do m = 1, num_wann
        do n = 1, num_wann
            ! if (n == m) cycle
            do i = 1, 3
                do j = 1, 3
                    do p = 1, Num_wann
                        if (p == m .or. p == n) cycle
                        gen_der_r(n, m, i, j) = gen_der_r(n, m, i, j) - V_ham(n, p, i)*D_Ham(p, m, j) &
                                                + D_ham(n, p, j)*V_ham(p, m, i)
                    enddo
                    gen_der_r(n, m, i, j) = gen_der_r(n, m, i, j) + (V_ham(n, m, i)*(V_ham(n, n, j)-V_ham(m, m, j)) &
                                                + (V_ham(n, n, i)-V_ham(m, m, i))*V_ham(n, m, j))/(W(n)-W(m)) - Wmn_ham(n, m, i, j)
                    gen_der_r(n, m, i, j) = cmplx_i*gen_der_r(n, m, i, j)/(W(n)-W(m))
                enddo
            enddo
        enddo
    enddo
end subroutine generalderivative

subroutine energygap
    use wmpi
    use para
    implicit none
    
   !------------------------------------------------------------------!
   !> This subroutine is to calculate the direct gap of semiconductor !
   !>                                                                 !
   !>                                                                 !
   !------------------------------------------------------------------!
    
    integer :: ik, ikx, iky, ikz, knv3, ifreq, i, j, m, n, index, ierr

    real(dp) :: k(3), time_start, time_end, x, fac_H, fac_AH
    real(dp) :: gap_min, VBM_max, CBM_min, indirect_gap
    complex(dp) :: cmplx_i, cmplx_1, cmplx_0, omega


    !> eigen value of H
    real(dp), allocatable :: W(:), occ(:)
    real(dp), allocatable :: gap(:), VBM(:), CBM(:)
    real(dp), allocatable :: gap_mpi(:), VBM_mpi(:), CBM_mpi(:)
    complex(dp), allocatable :: Hamk_bulk(:, :)
    complex(dp), allocatable :: UU(:, :)

    allocate( W (Num_wann))
    allocate( Hamk_bulk(Num_wann, Num_wann))
    allocate( UU(Num_wann, Num_wann))


    W = 0.0_dp
    Hamk_bulk = cmplx_0
    UU = cmplx_0
    knv3= Nk1*Nk2*Nk3

    allocate(occ(Num_wann))
    allocate(gap(knv3))
    allocate(CBM(knv3))
    allocate(VBM(knv3))
    allocate(gap_mpi(knv3))
    allocate(CBM_mpi(knv3))
    allocate(VBM_mpi(knv3))

    gap_mpi = 0.0_dp

    call now(time_start) 
    do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0.and. mod(ik/num_cpu, 100).eq.0) then
            call now(time_end) 
            write(stdout, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
            ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/100d0
            time_start= time_end
        endif

        ikx= (ik-1)/(nk2*nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
        ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
        k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
            + K3D_vec2_cube*(iky-1)/dble(nk2)  &
            + K3D_vec3_cube*(ikz-1)/dble(nk3)

        ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
        ! call ham_bulk_atomicgauge(k, Hamk_bulk)
        call ham_bulk_latticegauge(k, Hamk_bulk)

        !> diagonalization by call zheev in lapack
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)

        occ = 0.0_dp
        do i = 1, Num_wann
            if (W(i) < 0.0_dp) occ(i) = 1.0_dp
        end do

        do i = 1, Num_wann-1
            if (occ(i) == 1.0_dp .and. occ(i+1) == 0.0_dp) then
                gap_mpi(ik) = W(i+1) - W(i)
                exit
            endif
        enddo

        VBM_mpi(ik) = W(Numoccupied)
        CBM_mpi(ik) = W(Numoccupied+1)

    enddo ! ik

#if defined (MPI)
    call mpi_allreduce(gap_mpi, gap, size(gap), mpi_dp,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(VBM_mpi, VBM, size(VBM), mpi_dp,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(CBM_mpi, CBM, size(CBM), mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
    gap = gap_mpi
    VBM = VBM_mpi
    CBM = CBM_mpi
#endif
    
    !> get the minimal of the gap
    gap_min = minval(gap)
    !> get the maximal of the VBM
    VBM_max = maxval(VBM)
    !> get the minimal of the CBM
    CBM_min = minval(CBM)

    if (CBM_min > VBM_max) then
        indirect_gap = CBM_min - VBM_max
        if (cpuid .eq. 0 .and. indirect_gap .lt. gap_min) &
        write(stdout, '(a, E10.3, a)') 'The indirect gap is ', indirect_gap/eV2Hartree, ' eV'
    else
        indirect_gap = 0.0_dp
    endif

    if (cpuid .eq. 0) &
        write(stdout, '(a, E10.3, a)') 'The optical direct gap is ', gap_min/eV2Hartree, ' eV'
    
    
    outfileindex= outfileindex+ 1
    if (cpuid.eq.0) then
        open(unit=outfileindex, file='energygap.dat')
        write(outfileindex, "('#kx   ', 'ky   ', 'kz   ','gap  ')")
        do ik = 1, knv3
            ikx= (ik-1)/(nk2*nk3)+1
            iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
            ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
            k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
                + K3D_vec2_cube*(iky-1)/dble(nk2)  &
                + K3D_vec3_cube*(ikz-1)/dble(nk3)
        
            write(outfileindex, '(3f5.3, E12.3)') k(1), k(2), k(3), gap(ik)
        enddo
        close(outfileindex)
    endif

end subroutine energygap

