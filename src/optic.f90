subroutine linear_optic
    use wmpi
    use para
    implicit none

   !------------------------------------------------------------------!
   !> This subroutine calculates the interband optical conductivity   !
   !> within a conservative diagonal-TBA implementation.               !
   !> It does not include off-diagonal position matrix elements, a     !
   !> complete Wannier-gauge optical matrix, or a degenerate-subspace  !
   !> treatment.                                                       !
   !------------------------------------------------------------------!

    integer :: ik, ikx, iky, ikz, knv3, ifreq, i, j, m, n, index, ierr
    real(dp) :: k(3), time_start, time_end, x, fac_H, fac_AH
    real(dp) :: delta_E, k_weight, degen_tol, sigma_scale_H, sigma_scale_AH
    complex(dp) :: cmplx_i, cmplx_0, omega, kernel_AH

    complex(dp), allocatable :: Freq_array(:)

    !> Hermitian and anti-Hermitian parts of the interband conductivity tensor.
    complex(dp), allocatable :: sigma_kubo_H(:, :, :)
    complex(dp), allocatable :: sigma_kubo_H_mpi(:, :, :)
    complex(dp), allocatable :: sigma_kubo_AH(:, :, :)
    complex(dp), allocatable :: sigma_kubo_AH_mpi(:, :, :)

    real(dp), allocatable :: W(:)
    complex(dp), allocatable :: Hamk_bulk(:, :)
    complex(dp), allocatable :: UU(:, :)
    real(dp), allocatable :: occ(:)

    complex(dp), allocatable :: V_ham(:, :, :)
    complex(dp), allocatable :: D_ham(:, :, :)
    complex(dp), allocatable :: AA(:, :, :)

    real(dp), external :: delta

    real(dp) :: sigma_S(12), sigma_A(6)
    integer :: alpha_S(6), beta_S(6), alpha_A(3), beta_A(3)

    alpha_S = (/1, 2, 3, 1, 1, 2/)
    beta_S = (/1, 2, 3, 2, 3, 3/)
    alpha_A = (/2, 3, 1/)
    beta_A = (/3, 1, 2/)

    cmplx_i = (0.0d0, 1.0d0)
    cmplx_0 = (0.0d0, 0.0d0)
    degen_tol = eps6

    allocate(W(Num_wann))
    allocate(Hamk_bulk(Num_wann, Num_wann))
    allocate(UU(Num_wann, Num_wann))
    allocate(occ(Num_wann))
    allocate(D_ham(Num_wann, Num_wann, 3))
    allocate(V_ham(Num_wann, Num_wann, 3))
    allocate(AA(Num_wann, Num_wann, 3))
    allocate(sigma_kubo_H(3, 3, FreqNum))
    allocate(sigma_kubo_H_mpi(3, 3, FreqNum))
    allocate(sigma_kubo_AH(3, 3, FreqNum))
    allocate(sigma_kubo_AH_mpi(3, 3, FreqNum))
    allocate(Freq_array(FreqNum))

    W = 0.0_dp
    Hamk_bulk = cmplx_0
    UU = cmplx_0
    occ = 0.0_dp
    D_ham = cmplx_0
    V_ham = cmplx_0
    AA = cmplx_0
    sigma_kubo_H = cmplx_0
    sigma_kubo_H_mpi = cmplx_0
    sigma_kubo_AH = cmplx_0
    sigma_kubo_AH_mpi = cmplx_0
    Freq_array = cmplx_0

    if (FreqNum == 1) then
        Freq_array(1) = FreqMin + cmplx_i*eta_smr_fixed
    else
        do i = 1, FreqNum
            Freq_array(i) = FreqMin + (FreqMax-FreqMin)*(i-1d0)/dble(FreqNum-1) + cmplx_i*eta_smr_fixed
        end do
    endif

    knv3 = Nk1*Nk2*Nk3

    call now(time_start)
    do ik = 1+cpuid, knv3, num_cpu
        if (cpuid.eq.0 .and. mod(ik/num_cpu, 100).eq.0) then
            call now(time_end)
            write(stdout, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
                ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/100d0
            time_start = time_end
        endif

        ikx = (ik-1)/(nk2*nk3)+1
        iky = ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
        ikz = (ik-(iky-1)*Nk3-(ikx-1)*Nk2*Nk3)
        k = K3D_start_cube + K3D_vec1_cube*(ikx-1)/dble(nk1) &
            + K3D_vec2_cube*(iky-1)/dble(nk2) &
            + K3D_vec3_cube*(ikz-1)/dble(nk3)

        call ham_bulk_atomicgauge(k, Hamk_bulk)

        UU = Hamk_bulk
        call eigensystem_c('V', 'U', Num_wann, UU, W)

        occ = 0.0_dp
        do i = 1, Num_wann
            if (W(i) < 0.0_dp) occ(i) = 1.0_dp
        end do

        call dHdk_atomicgauge_Ham(k, UU, V_Ham)
        call get_Dmn_Ham_safe(W, V_Ham, degen_tol, D_Ham)
        AA = cmplx_i*D_Ham

        do ifreq = 1, FreqNum
            omega = Freq_array(ifreq)
            do m = 1, Num_wann
                do n = 1, Num_wann
                    if (n == m) cycle
                    if (abs(occ(m)-occ(n)) < eps9) cycle

                    delta_E = W(m) - W(n)
                    x = delta_E - real(omega, dp)
                    fac_H = delta(eta_smr_fixed, x)*(occ(m)-occ(n))*delta_E

                    kernel_AH = cmplx(delta_E, 0.0_dp, dp) / (cmplx(delta_E, 0.0_dp, dp) - omega)
                    fac_AH = (occ(m)-occ(n))*real(kernel_AH, dp)

                    do j = 1, 3
                        do i = 1, 3
                            sigma_kubo_H_mpi(i, j, ifreq) = sigma_kubo_H_mpi(i, j, ifreq) &
                                + AA(n, m, i)*AA(m, n, j)*fac_H
                            sigma_kubo_AH_mpi(i, j, ifreq) = sigma_kubo_AH_mpi(i, j, ifreq) &
                                + AA(n, m, i)*AA(m, n, j)*fac_AH
                        end do
                    end do
                end do
            end do
        end do
    end do

#if defined (MPI)
    call mpi_allreduce(sigma_kubo_H_mpi, sigma_kubo_H, size(sigma_kubo_H), mpi_dc, mpi_sum, mpi_cmw, ierr)
    call mpi_allreduce(sigma_kubo_AH_mpi, sigma_kubo_AH, size(sigma_kubo_AH), mpi_dc, mpi_sum, mpi_cmw, ierr)
#else
    sigma_kubo_H = sigma_kubo_H_mpi
    sigma_kubo_AH = sigma_kubo_AH_mpi
#endif

    !> The accumulated sums are normalized to the sampled k volume and then
    !> converted from atomic units to the interband optical conductivity in S/cm.
    k_weight = kCubeVolume / Origin_cell%ReciprocalCellVolume / dble(knv3)
    sigma_kubo_H = sigma_kubo_H * k_weight
    sigma_kubo_AH = sigma_kubo_AH * k_weight

    sigma_scale_H = -pi * Echarge**2 / (hbar * Bohr_radius * Origin_cell%CellVolume) / 100d0
    sigma_scale_AH = Echarge**2 / (hbar * Bohr_radius * Origin_cell%CellVolume) / 100d0
    sigma_kubo_H = sigma_kubo_H * sigma_scale_H
    sigma_kubo_AH = cmplx_i * sigma_kubo_AH * sigma_scale_AH

    call validate_complex_tensor_3d('sigma_kubo_H', sigma_kubo_H, 3, 3, FreqNum)
    call validate_complex_tensor_3d('sigma_kubo_AH', sigma_kubo_AH, 3, 3, FreqNum)

    if (cpuid.eq.0) then
        outfileindex = outfileindex + 1
        open(unit=outfileindex, file='sigma_kubo_symm.dat')
        write(outfileindex, '("#",a)') 'Symmetric interband optical conductivity in S/cm'
        write(outfileindex, '("#",a)') 'Diagonal TBA; 0 K step occupation; no intraband/Drude contribution'
        write(outfileindex, '("#",a)') '3D bulk normalization; for slab/supercell models the result depends on the vacuum thickness'
        write(outfileindex, '("#",a)') 'No off-diagonal position matrix elements or degenerate-subspace treatment are included'
        write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 13)
        write(outfileindex, '("#",a16, 20a16)') 'Frequency (eV)', 'Re[\sigma_xx]', 'Im[\sigma_xx]', 'Re[\sigma_yy]', &
            'Im[\sigma_yy]', 'Re[\sigma_zz]', 'Im[\sigma_zz]', 'Re[\sigma_xy]', 'Im[\sigma_xy]', &
            'Re[\sigma_xz]', 'Im[\sigma_xz]', 'Re[\sigma_yz]', 'Im[\sigma_yz]'
        do ifreq = 1, FreqNum
            do index = 1, 6
                i = alpha_S(index)
                j = beta_S(index)
                sigma_S(index*2-1) = real(0.5_dp*(sigma_kubo_H(i, j, ifreq) + sigma_kubo_H(j, i, ifreq)), dp)
                sigma_S(index*2) = aimag(0.5_dp*(sigma_kubo_AH(i, j, ifreq) + sigma_kubo_AH(j, i, ifreq)))
            end do
            write(outfileindex, '(200E16.8)') real(Freq_array(ifreq), dp)/eV2Hartree, sigma_S
        end do
        close(outfileindex)

        outfileindex = outfileindex + 1
        open(unit=outfileindex, file='sigma_kubo_asymm.dat')
        write(outfileindex, '("#",a)') 'Antisymmetric interband optical conductivity in S/cm'
        write(outfileindex, '("#",a)') 'Diagonal TBA; 0 K step occupation; no intraband/Drude contribution'
        write(outfileindex, '("#",a)') '3D bulk normalization; for slab/supercell models the result depends on the vacuum thickness'
        write(outfileindex, '("#",a)') 'No off-diagonal position matrix elements or degenerate-subspace treatment are included'
        write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 7)
        write(outfileindex, '("#",a16, 20a16)') 'Frequency (eV)', 'Re[\sigma_yz]', 'Im[\sigma_yz]', 'Re[\sigma_zx]', &
            'Im[\sigma_zx]', 'Re[\sigma_xy]', 'Im[\sigma_xy]'
        do ifreq = 1, FreqNum
            do index = 1, 3
                i = alpha_A(index)
                j = beta_A(index)
                sigma_A(index*2-1) = real(0.5_dp*(sigma_kubo_AH(i, j, ifreq) - sigma_kubo_AH(j, i, ifreq)), dp)
                sigma_A(index*2) = aimag(0.5_dp*(sigma_kubo_H(i, j, ifreq) - sigma_kubo_H(j, i, ifreq)))
            end do
            write(outfileindex, '(200E16.8)') real(Freq_array(ifreq), dp)/eV2Hartree, sigma_A
        end do
        close(outfileindex)
    endif

    if (cpuid.eq.0) write(stdout, '(a)') ' <<< The linear optic conductivity calculation finished'
    if (cpuid.eq.0) write(stdout, '(a)') ' '

    return

end subroutine linear_optic

subroutine bulk_photovoltaic
    use wmpi
    use para
    implicit none
   !------------------------------------------------------------------!
   !> This subroutine is to calculate the bulk photovotaic effect     !
   !> including shift current and injection current                   !
   !> Use TBA approximation                                           !
   !>                                                                 !
   !> Written by Hanqi Pi (hqpi1999@gmail.com)                        !
   !> Modified by Yi Yan (yanyi200207@gmail.com)                      ! 
   !>                                                                 !
   !> References : [1] PHYSICAL REVIEW B 97, 245143 (2018)            !
   !>              [2] Quantum Front 2, 6 (2023)                      !
   !------------------------------------------------------------------!

    integer :: ik, ikx, iky, ikz, knv3, ifreq, i, a, b, c, m, n, ierr

    real(dp) :: k(3), time_start, time_end, xplus, xminus
    real(dp) :: k_weight, degen_tol, shift_scale, inject_scale
    complex(dp) :: cmplx_i, cmplx_0, omega

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

    integer :: alpha_S(6), beta_S(6)

    alpha_S = (/1, 1, 1, 2, 2, 3/)
    beta_S = (/1, 2, 3, 2, 3, 3/)

    cmplx_i = (0.0d0, 1.0d0)
    cmplx_0 = (0.0d0, 0.0d0)
    degen_tol = eps6

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
        
        call dHdk_atomicgauge_Ham(k, UU, V_Ham)
        call get_Dmn_Ham_safe(W, V_Ham, degen_tol, D_Ham)
        AA = cmplx_i*D_Ham

        dHdkdk = cmplx_0
        Wmn_ham = cmplx_0
        call d2Hdk2_atomicgauge_wann(k, dHdkdk)
        call d2Hdk2_atomicgauge_Ham(UU, dHdkdk, Wmn_ham)
        call generalderivative(W, V_Ham, D_Ham, Wmn_Ham, degen_tol, gen_der_r)

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

    !> The internal energy axis is Hartree. For shift current the energy delta
    !> contributes one factor of the atomic unit of time (Time_atomic = hbar/Eh),
    !> while injection current uses a current-rate kernel and keeps the extra 1/s.
    !> All responses below remain 3D bulk-normalized; slab/supercell values
    !> therefore depend on the vacuum thickness.
    k_weight = kCubeVolume / Origin_cell%ReciprocalCellVolume / dble(knv3)
    lshiftcur = lshiftcur * k_weight
    cshiftcur = cshiftcur * k_weight
    linjectcur = linjectcur * k_weight
    cinjectcur = cinjectcur * k_weight

    shift_scale = Time_atomic * pi * Echarge**3 / (4d0*hbar**2*Origin_cell%CellVolume)
    inject_scale = pi * Echarge**3 / (2d0*hbar**2*Origin_cell%CellVolume)

    lshiftcur = lshiftcur * shift_scale
    cshiftcur = cshiftcur * shift_scale
    linjectcur = linjectcur * inject_scale
    cinjectcur = cinjectcur * inject_scale

    call validate_real_tensor_3d('lshiftcur', lshiftcur, 3, 6, FreqNum)
    call validate_real_tensor_3d('cshiftcur', cshiftcur, 3, 6, FreqNum)
    call validate_real_tensor_3d('linjectcur', linjectcur, 3, 6, FreqNum)
    call validate_real_tensor_3d('cinjectcur', cinjectcur, 3, 6, FreqNum)


    if (cpuid.eq.0) then
        outfileindex= outfileindex+ 1
        open(unit=outfileindex, file='linear_shift.dat')
        write(outfileindex, '("#",a)') 'Linear shift current coefficient in A/V^2'
        write(outfileindex, '("#",a)') 'Interband only; diagonal TBA; 0 K step occupation; 3D bulk normalization'
        write(outfileindex, '("#",a)') 'For slab/supercell models this bulk response depends on the vacuum thickness'
        write(outfileindex, '("#",a)') 'No off-diagonal position matrix elements or degenerate-subspace treatment are included'
        write(outfileindex, '("#",a)') 'Near true degeneracies or strong near-degeneracies the BPVE can remain gauge/mesh sensitive'
        write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 19)
        write(outfileindex, '("#",a16, 20a16)')'Frequency (eV)', 'xxx', 'xxy', 'xxz', 'xyy', 'xyz', 'xzz', &
                                                'yxx', 'yxy', 'yxz', 'yyy', 'yyz', 'yzz', 'zxx', 'zxy', 'zxz', &
                                                'zyy', 'zyz', 'zzz'
        do ifreq=1, FreqNum
            write(outfileindex, '(200E16.8)') real(Freq_array(ifreq), dp)/eV2Hartree, lshiftcur(1, :, ifreq),&
                                              lshiftcur(2, :, ifreq), lshiftcur(3, :, ifreq)
        enddo 
        close(outfileindex)

        outfileindex= outfileindex+ 1
        open(unit=outfileindex, file='circular_shift.dat')
        write(outfileindex, '("#",a)') 'Circular shift current coefficient in A/V^2'
        write(outfileindex, '("#",a)') 'Interband only; diagonal TBA; 0 K step occupation; 3D bulk normalization'
        write(outfileindex, '("#",a)') 'For slab/supercell models this bulk response depends on the vacuum thickness'
        write(outfileindex, '("#",a)') 'No off-diagonal position matrix elements or degenerate-subspace treatment are included'
        write(outfileindex, '("#",a)') 'Near true degeneracies or strong near-degeneracies the BPVE can remain gauge/mesh sensitive'
        write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 19)
        write(outfileindex, '("#",a16, 20a16)')'Frequency (eV)', 'xxx', 'xxy', 'xxz', 'xyy', 'xyz', 'xzz', &
                                                'yxx', 'yxy', 'yxz', 'yyy', 'yyz', 'yzz', 'zxx', 'zxy', 'zxz', &
                                                'zyy', 'zyz', 'zzz'
        do ifreq=1, FreqNum
            write(outfileindex, '(200E16.8)') real(Freq_array(ifreq), dp)/eV2Hartree, cshiftcur(1, :, ifreq),&
                                              cshiftcur(2, :, ifreq), cshiftcur(3, :, ifreq)
        enddo 
        close(outfileindex)

        outfileindex= outfileindex+ 1
        open(unit=outfileindex, file='circular_inject_rate.dat')
        write(outfileindex, '("#",a)') 'Circular injection current rate in A/(V^2 s)'
        write(outfileindex, '("#",a)') 'Interband only; diagonal TBA; 0 K step occupation; 3D bulk normalization'
        write(outfileindex, '("#",a)') 'No relaxation-time factor is included; multiply by a model-specific tau for a steady-state estimate'
        write(outfileindex, '("#",a)') 'For slab/supercell models this bulk response depends on the vacuum thickness'
        write(outfileindex, '("#",a)') 'No off-diagonal position matrix elements or degenerate-subspace treatment are included'
        write(outfileindex, '("#",a)') 'Near true degeneracies or strong near-degeneracies the BPVE can remain gauge/mesh sensitive'
        write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 19)
        write(outfileindex, '("#",a16, 20a16)') 'Frequency (eV)', 'xxx', 'xxy', 'xxz', 'xyy', 'xyz', 'xzz', &
                                                'yxx', 'yxy', 'yxz', 'yyy', 'yyz', 'yzz', 'zxx', 'zxy', 'zxz', &
                                                'zyy', 'zyz', 'zzz'
        do ifreq=1, FreqNum
            write(outfileindex, '(200E16.8)') real(Freq_array(ifreq), dp)/eV2Hartree, cinjectcur(1, :, ifreq),&
                                              cinjectcur(2, :, ifreq), cinjectcur(3, :, ifreq)
        enddo 
        close(outfileindex)

        outfileindex= outfileindex+ 1
        open(unit=outfileindex, file='linear_inject_rate.dat')
        write(outfileindex, '("#",a)') 'Linear injection current rate in A/(V^2 s)'
        write(outfileindex, '("#",a)') 'Interband only; diagonal TBA; 0 K step occupation; 3D bulk normalization'
        write(outfileindex, '("#",a)') 'No relaxation-time factor is included; multiply by a model-specific tau for a steady-state estimate'
        write(outfileindex, '("#",a)') 'For slab/supercell models this bulk response depends on the vacuum thickness'
        write(outfileindex, '("#",a)') 'No off-diagonal position matrix elements or degenerate-subspace treatment are included'
        write(outfileindex, '("#",a)') 'Near true degeneracies or strong near-degeneracies the BPVE can remain gauge/mesh sensitive'
        write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 19)
        write(outfileindex, '("#",a16, 20a16)') 'Frequency (eV)', 'xxx', 'xxy', 'xxz', 'xyy', 'xyz', 'xzz', &
                                                'yxx', 'yxy', 'yxz', 'yyy', 'yyz', 'yzz', 'zxx', 'zxy', 'zxz', &
                                                'zyy', 'zyz', 'zzz'
        do ifreq=1, FreqNum
            write(outfileindex, '(200E16.8)') real(Freq_array(ifreq), dp)/eV2Hartree, linjectcur(1, :, ifreq),&
                                              linjectcur(2, :, ifreq), linjectcur(3, :, ifreq)
        enddo 
        close(outfileindex)
    endif

    outfileindex = outfileindex + 1
    !> write script for gnuplot
    if (cpuid == 0) then
       open(unit=outfileindex, file='bpve.gnu')
       write(outfileindex, '(a)') '# BPVE Spectrum: Type and Component'
       write(outfileindex, '(a)') 'type = "LS"    # LS, LI, CS, CI'
       write(outfileindex, '(a)') 'col = 11        # See below '
    
       write(outfileindex, '(a)') ''
       write(outfileindex, '(a)') 'component = ""'
       write(outfileindex, '(a)') 'if (col == 2) component = "xxx"'
       write(outfileindex, '(a)') 'if (col == 3) component = "xxy"'
       write(outfileindex, '(a)') 'if (col == 4) component = "xxz"'
       write(outfileindex, '(a)') 'if (col == 5) component = "xyy"'
       write(outfileindex, '(a)') 'if (col == 6) component = "xyz"'
       write(outfileindex, '(a)') 'if (col == 7) component = "xzz"'
       write(outfileindex, '(a)') 'if (col == 8) component = "yxx"'
       write(outfileindex, '(a)') 'if (col == 9) component = "yxy"'
       write(outfileindex, '(a)') 'if (col == 10) component = "yxz"'
       write(outfileindex, '(a)') 'if (col == 11) component = "yyy"'
       write(outfileindex, '(a)') 'if (col == 12) component = "yyz"'
       write(outfileindex, '(a)') 'if (col == 13) component = "yzz"'
       write(outfileindex, '(a)') 'if (col == 14) component = "zxx"'
       write(outfileindex, '(a)') 'if (col == 15) component = "zxy"'
       write(outfileindex, '(a)') 'if (col == 16) component = "zxz"'
       write(outfileindex, '(a)') 'if (col == 17) component = "zyy"'
       write(outfileindex, '(a)') 'if (col == 18) component = "zyz"'
       write(outfileindex, '(a)') 'if (col == 19) component = "zzz"'

       write(outfileindex, '(a)') ''
       write(outfileindex, '(a)') 'if (type eq "LS") filename = "linear_shift.dat"'
       write(outfileindex, '(a)') 'if (type eq "LI") filename = "linear_inject_rate.dat"'
       write(outfileindex, '(a)') 'if (type eq "CS") filename = "circular_shift.dat"'
        write(outfileindex, '(a)') 'if (type eq "CI") filename = "circular_inject_rate.dat"'

       write(outfileindex, '(a)') ''
       write(outfileindex, '(a)') 'outfile = sprintf("sigma_%s_%s.pdf", type, component)'
       write(outfileindex, '(a)') ''
       write(outfileindex, '(a)') 'set terminal pdfcairo enhanced font "Arial,30" size 6,6'
       write(outfileindex, '(a)') 'set output outfile'
       write(outfileindex, '(a)') ''
       write(outfileindex, '(a)') 'ylabel = "Response"'
       write(outfileindex, '(a)') 'if (type eq "LS" || type eq "CS") ylabel = "Shift coefficient (A/V^{2})"'
       write(outfileindex, '(a)') 'if (type eq "LI" || type eq "CI") ylabel = "Injection rate (A/V^{2}/s)"'
       write(outfileindex, '(a)') 'set title sprintf("BPVE Response (%s-%s)", type, component)'
       write(outfileindex, '(a)') 'set xlabel "Frequency (eV)"'
       write(outfileindex, '(a)') 'set ylabel ylabel'
       write(outfileindex, '(a)') 'set grid'
       write(outfileindex, '(a)') ''
       write(outfileindex, '(a)') 'plot filename using 1:col title component with lines lw 5'
    endif

    if (cpuid.eq.0) write(stdout, '(a)') ' <<< The Bulk Photovoltaic Effect calculation finished'
    if (cpuid.eq.0) write(stdout, '(a)') ' '

    return

end subroutine bulk_photovoltaic

subroutine get_Dmn_Ham_safe(W, velocity_Ham, degen_tol, Dmn_Ham)
    use para, only : dp, Num_wann, zzero
    implicit none

    real(dp), intent(in) :: W(Num_wann)
    complex(dp), intent(in) :: velocity_Ham(Num_wann, Num_wann, 3)
    real(dp), intent(in) :: degen_tol
    complex(dp), intent(out) :: Dmn_Ham(Num_wann, Num_wann, 3)

    integer :: m, n, i

    Dmn_Ham = zzero
    do i = 1, 3
        do n = 1, Num_wann
            do m = 1, Num_wann
                if (m == n) cycle
                if (abs(W(n)-W(m)) <= degen_tol) cycle
                Dmn_Ham(m, n, i) = velocity_Ham(m, n, i)/(W(n)-W(m))
            end do
        end do
    end do
end subroutine get_Dmn_Ham_safe

subroutine validate_complex_tensor_3d(label, tensor, n1, n2, n3)
    use para, only : dp, stdout
    use wmpi, only : cpuid
    implicit none

    character(*), intent(in) :: label
    integer, intent(in) :: n1, n2, n3
    complex(dp), intent(in) :: tensor(n1, n2, n3)

    integer :: i1, i2, i3
    real(dp) :: re_part, im_part

    do i3 = 1, n3
        do i2 = 1, n2
            do i1 = 1, n1
                re_part = real(tensor(i1, i2, i3), dp)
                im_part = aimag(tensor(i1, i2, i3))
                if ((re_part /= re_part) .or. (im_part /= im_part)) then
                    if (cpuid.eq.0) write(stdout, '(a, 3i6)') 'ERROR: NaN detected in '//trim(label)//' at', i1, i2, i3
                    stop 1
                endif
                if (abs(re_part) > huge(re_part) .or. abs(im_part) > huge(im_part)) then
                    if (cpuid.eq.0) write(stdout, '(a, 3i6)') 'ERROR: Inf detected in '//trim(label)//' at', i1, i2, i3
                    stop 1
                endif
            end do
        end do
    end do
end subroutine validate_complex_tensor_3d

subroutine validate_real_tensor_3d(label, tensor, n1, n2, n3)
    use para, only : dp, stdout
    use wmpi, only : cpuid
    implicit none

    character(*), intent(in) :: label
    integer, intent(in) :: n1, n2, n3
    real(dp), intent(in) :: tensor(n1, n2, n3)

    integer :: i1, i2, i3

    do i3 = 1, n3
        do i2 = 1, n2
            do i1 = 1, n1
                if (tensor(i1, i2, i3) /= tensor(i1, i2, i3)) then
                    if (cpuid.eq.0) write(stdout, '(a, 3i6)') 'ERROR: NaN detected in '//trim(label)//' at', i1, i2, i3
                    stop 1
                endif
                if (abs(tensor(i1, i2, i3)) > huge(tensor(i1, i2, i3))) then
                    if (cpuid.eq.0) write(stdout, '(a, 3i6)') 'ERROR: Inf detected in '//trim(label)//' at', i1, i2, i3
                    stop 1
                endif
            end do
        end do
    end do
end subroutine validate_real_tensor_3d

subroutine generalderivative(W, V_Ham, D_Ham, Wmn_Ham, degen_tol, gen_der_r)
    use para, only : Num_wann, dp, zzero
    implicit none

    real(dp), intent(in) :: W(Num_wann)
    complex(dp), intent(in) :: V_Ham(Num_wann, Num_wann, 3)
    complex(dp), intent(in) :: D_Ham(Num_wann, Num_wann, 3)
    complex(dp), intent(in) :: Wmn_Ham(Num_wann, Num_wann, 3, 3)
    real(dp), intent(in) :: degen_tol
    complex(dp), intent(out) :: gen_der_r(Num_wann, Num_wann, 3, 3)

    integer :: m, n, i, j, p
    real(dp) :: delta_nm
    complex(dp) :: cmplx_i

    cmplx_i = (0.0d0, 1.0d0)
    gen_der_r = zzero
    do m = 1, Num_wann
        do n = 1, Num_wann
            if (n == m) cycle
            delta_nm = W(n) - W(m)
            if (abs(delta_nm) <= degen_tol) cycle
            do i = 1, 3
                do j = 1, 3
                    do p = 1, Num_wann
                        if (p == m .or. p == n) cycle
                        gen_der_r(n, m, i, j) = gen_der_r(n, m, i, j) - V_ham(n, p, i)*D_Ham(p, m, j) &
                                                + D_ham(n, p, j)*V_ham(p, m, i)
                    end do
                    gen_der_r(n, m, i, j) = gen_der_r(n, m, i, j) + (V_ham(n, m, i)*(V_ham(n, n, j)-V_ham(m, m, j)) &
                                                + (V_ham(n, n, i)-V_ham(m, m, i))*V_ham(n, m, j))/delta_nm - Wmn_ham(n, m, i, j)
                    gen_der_r(n, m, i, j) = cmplx_i*gen_der_r(n, m, i, j)/delta_nm
                end do
            end do
        end do
    end do
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
    
    integer :: ik, ikx, iky, ikz, knv3, i, ierr
    real(dp) :: k(3), time_start, time_end
    real(dp) :: gap_min, VBM_max, CBM_min, indirect_gap

    !> eigen value of H
    real(dp), allocatable :: W(:), occ(:)
    real(dp), allocatable :: gap(:), VBM(:), CBM(:)
    real(dp), allocatable :: gap_mpi(:), VBM_mpi(:), CBM_mpi(:)
    complex(dp), allocatable :: Hamk_bulk(:, :)
    complex(dp), allocatable :: UU(:, :)

    knv3 = Nk1*Nk2*Nk3

    allocate( W (Num_wann))
    allocate( Hamk_bulk(Num_wann, Num_wann))
    allocate( UU(Num_wann, Num_wann))
    allocate(occ(Num_wann))
    allocate(gap(knv3))
    allocate(CBM(knv3))
    allocate(VBM(knv3))
    allocate(gap_mpi(knv3))
    allocate(CBM_mpi(knv3))
    allocate(VBM_mpi(knv3))

    W = 0.0_dp
    Hamk_bulk = (0.0_dp, 0.0_dp)
    UU = (0.0_dp, 0.0_dp)
    occ = 0.0_dp
    gap = 0.0_dp
    CBM = 0.0_dp
    VBM = 0.0_dp
    gap_mpi = 0.0_dp
    CBM_mpi = 0.0_dp
    VBM_mpi = 0.0_dp

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
        write(stdout, '(a, E10.3, a)') 'The minimum direct single-particle gap is ', gap_min/eV2Hartree, ' eV'
    
    
    outfileindex= outfileindex+ 1
    if (cpuid.eq.0) then
        open(unit=outfileindex, file='energygap.dat')
        write(outfileindex, '("#",a)') 'Minimum direct single-particle gap on the sampled k mesh'
        write(outfileindex, '("#",a)') 'This is a direct single-particle gap only; no optical matrix elements or selection rules are included'
        write(outfileindex, '("#",a)') 'k_x  k_y  k_z  gap(eV)'
        do ik = 1, knv3
            ikx= (ik-1)/(nk2*nk3)+1
            iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
            ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
            k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
                + K3D_vec2_cube*(iky-1)/dble(nk2)  &
                + K3D_vec3_cube*(ikz-1)/dble(nk3)

            write(outfileindex, '(3f10.6, 1x, E16.8)') k(1), k(2), k(3), gap(ik)/eV2Hartree
        enddo
        close(outfileindex)
    endif

end subroutine energygap

