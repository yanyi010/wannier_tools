subroutine floquet_hamiltonian_atomicgauge(k,H_floquet)
    !> This subroutine is used for calculating floquet hamiltonian
    !> 1. get hamiltonian from subroutine ham_bulk_atomicgauge
    !> 2. add vector potential to hamiltonian
    !> 3. get floquet hamiltonian

   
    use wmpi
    use para
    implicit none 

    real(dp), intent(in) :: k(3)
    complex(dp), intent(out) :: H_floquet(Num_wann, Num_wann)
    complex(dp) :: Hamk_bulk(Num_wann, Num_wann)

    !> 1. get hamiltonian from subroutine ham_bulk_atomicgauge
    call ham_bulk_atomicgauge(k, Hamk_bulk)

    !> 2. Peierls substitution(Herein we use the form of A_0 * cos(omega * t))
    



    !> 3. 当前先将 Floquet 哈密顿量设为静态哈密顿量
    H_floquet = Hamk_bulk



end subroutine floquet_hamiltonian_atomicgauge


subroutine floquet_band
    !> This subroutine is used for calculating band structure
    !> from floquet hamiltonian

    use wmpi
    use para
    implicit none 

    integer :: ik, il, ig, io, i, j, knv3, ierr
    real(dp) :: emin,  emax,  k(3)
    character*40 :: filename

    !> eigenvalues of H_floquet
    real(Dp), allocatable :: W(:)

    ! Floquet Hamiltonian
    complex(Dp), allocatable :: H_floquet(:, :)

    ! eigenectors of H
    real(dp), allocatable :: eigv(:,:), eigv_mpi(:,:)
    real(dp), allocatable :: weight(:,:,:), weight_mpi(:,:,:), weight_sum(:,:)

    knv3= nk3_band
    allocate(W(Num_wann))
    allocate(H_floquet(Num_wann, Num_wann))
    allocate( eigv    (Num_wann, knv3))
    allocate( eigv_mpi(Num_wann, knv3))
    allocate( weight    (NumberofSelectedOrbitals_groups,Num_wann, knv3))
    allocate( weight_mpi(NumberofSelectedOrbitals_groups,Num_wann, knv3))
    allocate( weight_sum(Num_wann, knv3))
    W       = 0d0; H_floquet = 0d0
    eigv    = 0d0; eigv_mpi= 0d0
    weight  = 0d0; weight_sum = 0d0; weight_mpi = 0d0

    do ik= 1+cpuid, knv3, num_cpu

       k = kpath_3d(:, ik)

       ! 计算 Floquet 哈密顿量
       H_floquet= 0d0
       call floquet_hamiltonian_atomicgauge(k, H_floquet)

       ! 对角化
       W= 0d0
       call eigensystem_c('V', 'U', Num_wann ,H_floquet, W)
       eigv(:, ik)= W

       do j= 1, Num_wann  !> band
          do ig= 1, NumberofSelectedOrbitals_groups
             do i= 1, NumberofSelectedOrbitals(ig)
                io= Selected_WannierOrbitals(ig)%iarray(i)
                weight(ig, j, ik)= weight(ig, j, ik)+ abs(H_floquet(io, j))**2 
             enddo !i
          enddo !ig
       enddo !j
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

    if (index(Particle,'phonon')/=0) then
       eigv_mpi = eigv_mpi - MINVAL(eigv_mpi)
    endif

    !> deal with phonon system
    if (index(Particle,'phonon')/=0) then
       do ik=1, knv3
          do j=1, Num_wann
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

          do i=1, Num_wann
             do ik=1+(il-1)*Nk1, il*Nk1
                write(outfileindex, '(2f19.9, 10000i5)')(k3len(ik)*Angstrom2atomic-k3len((il-1)*Nk1+1)*Angstrom2atomic),eigv_mpi(i, ik)
             enddo
             write(outfileindex, *)' '
          enddo ! i
          close(outfileindex)
       enddo ! il
    endif

    outfileindex= outfileindex+ nk3lines+1
    if (cpuid==0)then
       open(unit=outfileindex, file='floquetband.dat')

       write(outfileindex, &
          "('#', a12, a14, 1X, '| projection', 100(3X,'|group', i2, ': A '))")&
          'klen', 'E', (i, i=1, NumberofSelectedOrbitals_groups)
       write(outfileindex, "('#column', i5, 200i16)")(i, i=1, 2+NumberofSelectedOrbitals_groups)
       do i=1, Num_wann
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