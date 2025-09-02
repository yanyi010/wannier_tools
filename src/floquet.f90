subroutine floquet_hamiltonian_atomicgauge(k, H_floquet, dimF)
    !=========================================================
    ! Construct Floquet Hamiltonian using Peierls substitution
    !=========================================================
    use wmpi
    use para
    implicit none

    ! Input
    real(dp), intent(in) :: k(3)
    integer,  intent(in) :: dimF

  
    ! dimF = (2*N_Floquet+1)*Num_wann
    complex(Dp), intent(out) :: H_floquet(dimF, dimF)

    ! -------------------- parameters & locals --------------------
    real(dp) :: omega, T, dt, times
    real(dp) :: A0(3), kshift(3)
    real(dp) :: Ax_t, Ay_t, Az_t
    integer  :: i, j, m, n, p
    integer  :: idx_m, idx_n

    complex(dp), allocatable :: Ham_t(:,:,:)
    complex(dp) :: Hk(Num_wann, Num_wann)
    ! 注意：p = m-n 的范围是 [-2N, 2N]
    complex(dp), allocatable :: Hdiff(:,:,:)

    complex(dp) :: phase

    ! -------------------- field parameters -----------------------

    A0    = (/ A0_x, A0_y, A0_z /)  ! 线偏振-x
    omega = Floquet_omega
    T     = (twopi / omega) 
    dt    = T / real(Nt, dp)

    !write(stdout,*)'omega',omega, Floquet_omega
    !write(stdout,*)'T',T
    !write(stdout,*)'dt',dt

    ! -------------------- init -----------------------

    H_floquet = (0.0_dp, 0.0_dp)
    allocate(Ham_t(Num_wann, Num_wann, Nt))
    allocate(Hdiff(Num_wann, Num_wann, -2*N_Floquet:2*N_Floquet))
    Hdiff = (0.0_dp, 0.0_dp)

    ! -------------------- Step 1: build H(k+A(t)) over time (Peierls Substitution)------

    ! Two Problem: (1) Which gauge? (2) The substitution? e-charge and hbar?

    do i = 1, Nt
        times  = (i-1) * dt
        ! 圆偏振：A(t) = (A0_x cos(ωt), A0_y sin(ωt), A0_z cos(ωt))
        Ax_t   = A0(1) * cos(omega * times)
        Ay_t   = A0(2) * sin(omega * times)
        Az_t   = A0(3) * cos(omega * times)
        kshift(1) = k(1) + Ax_t
        kshift(2) = k(2) + Ay_t
        kshift(3) = k(3) + Az_t
        call ham_bulk_atomicgauge(kshift, Hk)
        !call ham_bulk_latticegauge(kshift, Hk)
        Ham_t(:,:,i) = Hk
    end do

    ! -------------------- Step 2: Fourier comps H_p (p=m-n) integration ------
    ! H_p = (1/T) ∫_0^T H(t) e^{i p ω t} dt

    do p = -2*N_Floquet, 2*N_Floquet  !> p= m-n, m, n ∈ [-N_F, N_F]. so p ∈ [-2*N_F, 2*N_F]
        Hdiff(:,:,p) = (0.0_dp, 0.0_dp) !> zero complex
        do i = 1, Nt !> Time-zone integration
            times = (i-1) * dt
            phase = exp((0.0_dp,1.0_dp) * real(p,dp) * omega * times)
            Hdiff(:,:,p) = Hdiff(:,:,p) + Ham_t(:,:,i) * phase * (dt / T)
        end do
    end do

    ! -------------------- Step 3: assemble big Floquet H ---------
    ! 块索引规则：m, n ∈ [-N_F, N_F]
    ! block(m,n) = H_{m-n}，且对角块再加 m*ω I
    do m = -N_Floquet, N_Floquet
        idx_m =  m + N_Floquet    ! 0-based块序号
        do n = -N_Floquet, N_Floquet
            idx_n =  n + N_Floquet
            ! 目标矩阵切片位置（1-based元素下标）
            H_floquet( idx_m*Num_wann + 1 : idx_m*Num_wann + Num_wann, &
                        idx_n*Num_wann + 1 : idx_n*Num_wann + Num_wann ) = &
                        Hdiff(:,:, m - n )
            if (m == n) then
                do j = 1, Num_wann
                    H_floquet(idx_m*Num_wann + j, idx_n*Num_wann + j) = &
                        H_floquet(idx_m*Num_wann + j, idx_n*Num_wann + j) + &
                        cmplx( m*omega, 0.0_dp, kind=dp )
                end do
            end if
        end do
    end do

    ! -------------------- cleanup --------------------------------
    deallocate(Ham_t, Hdiff)

end subroutine floquet_hamiltonian_atomicgauge


subroutine floquet_band
    !> This subroutine is used for calculating band structure
    !> from floquet hamiltonian

    use wmpi
    use para
    implicit none 

    integer :: ik, il, ig, io, i, j, knv3, ierr, dimF
    integer :: m0_start, m0_end
    real(dp) :: emin,  emax,  k(3)
    character*40 :: filename

    !> eigenvalues of H_floquet
    real(Dp), allocatable :: W(:)

    ! Floquet Hamiltonian
    complex(Dp), allocatable :: H_floquet(:, :)

    ! eigenectors of H
    real(dp), allocatable :: eigv(:,:), eigv_mpi(:,:)
    real(dp), allocatable :: weight(:,:,:), weight_mpi(:,:,:), weight_sum(:,:)

    dimF = (2*N_Floquet+1)*Num_wann

    knv3= nk3_band
    allocate(W(dimF))
    allocate(H_floquet(dimF, dimF))
    allocate( eigv    (dimF, knv3))
    allocate( eigv_mpi(dimF, knv3))
    allocate( weight    (NumberofSelectedOrbitals_groups,dimF, knv3))
    allocate( weight_mpi(NumberofSelectedOrbitals_groups,dimF, knv3))
    allocate( weight_sum(dimF, knv3))
    W       = 0d0; H_floquet = 0d0
    eigv    = 0d0; eigv_mpi= 0d0
    weight  = 0d0; weight_sum = 0d0; weight_mpi = 0d0

    do ik= 1+cpuid, knv3, num_cpu

       k = kpath_3d(:, ik)

       ! 计算 Floquet 哈密顿量（为每个 k 构造 dimF×dimF 的 Floquet 大矩阵）
       H_floquet= 0d0
       call floquet_hamiltonian_atomicgauge(k, H_floquet, dimF)
       write(stdout,*)'ik',ik

       ! 对角化
       W= 0d0
       call eigensystem_c('V', 'U', dimF ,H_floquet, W)
       eigv(:, ik)= W

       do j= 1, dimF  !> band
          do ig= 1, NumberofSelectedOrbitals_groups
             do i= 1, NumberofSelectedOrbitals(ig)
                io= Selected_WannierOrbitals(ig)%iarray(i)
                if (io<=dimF) then
                   weight(ig, j, ik)= weight(ig, j, ik)+ abs(H_floquet(io, j))**2 
                endif
             enddo !i
          enddo !ig
       enddo !j

       ! m=0 Floquet 子块投影权重（静态组分）
       m0_start = N_Floquet*Num_wann + 1
       m0_end   = (N_Floquet+1)*Num_wann
       do j = 1, dimF
          weight_sum(j, ik) = 0d0
          do i = m0_start, m0_end
             weight_sum(j, ik) = weight_sum(j, ik) + abs(H_floquet(i, j))**2
          enddo
       enddo
    enddo !ik

#if defined (MPI)
    call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
       mpi_dp,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(weight, weight_mpi,size(weight),&
       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
    eigv_mpi= eigv
    weight_mpi= weight
#endif

#if defined (MPI)
    call mpi_allreduce(MPI_IN_PLACE, weight_sum, size(weight_sum), &
       mpi_dp, mpi_sum, mpi_cmw, ierr)
#endif

    if (index(Particle,'phonon')/=0) then
       eigv_mpi = eigv_mpi - MINVAL(eigv_mpi)
    endif

    !> deal with phonon system
    if (index(Particle,'phonon')/=0) then
       do ik=1, knv3
          do j=1, dimF
             eigv_mpi(j, ik)= sqrt(abs(eigv_mpi(j, ik)))*sign(1d0, eigv_mpi(j, ik))
          enddo
       enddo
    endif
    eigv_mpi= eigv_mpi/eV2Hartree

    outfileindex= outfileindex+ 1
    if (cpuid==0)then
       do il= 1, nk3lines
          if (il<10) then
             write(filename, '(a,i1)')'floquetband.dat-segment', il
          elseif (il>=10.and.il<100) then
             write(filename, '(a,i2)')'floquetband.dat-segment', il
          elseif (il>=100.and.il<1000) then
             write(filename, '(a,i3)')'floquetband.dat-segment', il
          endif
          write(filename, '(5a)')'floquetband.dat-', trim(adjustl(k3line_name(il))), '-', trim(adjustl(k3line_name(il+1)))
          outfileindex= outfileindex+ 1
          open(unit=outfileindex, file=filename)
          write(outfileindex, '(a, a6, a, a6)')'# segment between', k3line_name(il), ' and ',  k3line_name(il+1)
          write(outfileindex, '(a, 3f10.6)')'# kstart (fractional units)', k3line_start(:, il)
          write(outfileindex, '(a, 3f10.6)')'# kend   (fractional units)', k3line_end(:, il)

          do i=1, dimF
             do ik=1+(il-1)*Nk1, il*Nk1
                write(outfileindex, '(2f19.9, 10000i5)')(k3len(ik)*Angstrom2atomic-k3len((il-1)*Nk1+1)*Angstrom2atomic),eigv_mpi(i, ik)
             enddo
             write(outfileindex, *)' '
          enddo ! i
          close(outfileindex)
       enddo ! il
    endif

    ! 输出 m=0 Floquet 组分权重
    outfileindex= outfileindex+ 1
    if (cpuid==0)then
       open(unit=outfileindex, file='floquet_m0_weight.dat')
       write(outfileindex, '(2a19)')'% klen', 'w_m0(n)'
       do ik=1, knv3
          write(outfileindex, '(2000f19.9)')k3len(ik)*Angstrom2atomic, weight_sum(:, ik)
       enddo
       close(outfileindex)
    endif

    outfileindex= outfileindex+ nk3lines+1
    if (cpuid==0)then
       open(unit=outfileindex, file='floquetband.dat')

       write(outfileindex, &
          "('#', a12, a14, 1X, '| projection', 100(3X,'|group', i2, ': A '))")&
          'klen', 'E', (i, i=1, NumberofSelectedOrbitals_groups)
       write(outfileindex, "('#column', i5, 200i16)")(i, i=1, 2+NumberofSelectedOrbitals_groups)
       do i=1, dimF
          do ik=1, knv3
             write(outfileindex, '(200f16.9)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik), &
                weight_mpi(:, i, ik)
          enddo
          write(outfileindex, *)' '
       enddo
       close(outfileindex)
    endif

    outfileindex= outfileindex+ nk3lines+1
    if (cpuid==0)then
       open(unit=outfileindex, file='floquetband.dat-matlab')
       write(outfileindex, '(2a19)')'% klen', 'E(n)'

       do ik=1, knv3
          write(outfileindex, '(2000f19.9)')k3len(ik)*Angstrom2atomic,eigv_mpi(:, ik)
       enddo
       close(outfileindex)
    endif

    !> minimum and maximum value of energy bands
    emin=  minval(eigv_mpi)-0.5d0
    emax=  maxval(eigv_mpi)+0.5d0

    call generate_ek_kpath_gnu('floquetband.dat', 'floquetband.gnu', 'floquetband.pdf', &
                                  emin, emax, knv3, Nk3lines, &
                                  k3line_name, k3line_stop, k3len)

    deallocate(W)
    deallocate(H_floquet)
    deallocate( eigv    )
    deallocate( eigv_mpi)
    deallocate( weight    )
    deallocate( weight_mpi)
    deallocate( weight_sum)

end subroutine floquet_band