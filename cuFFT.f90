!! compile: pgf90 -acc -Minfo=accel -O2 -ta=tesla,cc75 -Mcudalib=cufft cuFFT.f90
program cufft2dTest
    use cufft
    use openacc
    integer      , parameter :: nx = 128           ! Grid number in x direction
    integer      , parameter :: ny = 128            ! Grid number in y direction
    integer      , parameter :: mnx = nx/2+1            
    real, parameter :: pi=3.14159265
    complex(kind(0d0)), allocatable  :: b(:,:)
    double precision, allocatable     :: r(:,:),q(:,:)
    integer :: iplan1, iplan2, iplan3,i,j, ierr
    real :: x,y,dx,dy,ceff_norm
    character*80 :: fname

    allocate(r(nx,ny),q(nx,ny))
    allocate(b(mnx,ny))
    ierr=0
    ceff_norm = 1.0/(nx*ny)

    x = 2*pi
    y = 2*pi
    dx = x/nx
    dy = y/ny
    do i=1,nx
        do j=1,ny
            r(i,j) = sin(2*dx*i+3*dy*j)
        enddo
    enddo
!    a = 1; r = 1

!    ierr = cufftPlan2D(iplan1,m,n,CUFFT_C2C)
!!$acc data copyin(a)  create(b) copyout(c)
!    !$acc host_data use_device(a,b,c)
!    ierr = ierr + cufftExecZ2Z(iplan1,a,b,CUFFT_FORWARD)
!    ierr = ierr + cufftExecZ2Z(iplan1,b,c,CUFFT_INVERSE)
!    !$acc end host_data
!
!    ! scale c
!    !$acc kernels
!    c = c / (m*n)
!    !$acc end kernels
!!$acc end data
!
!! Check forward answer
!    write(*,*) 'Max error C2C FWD: ', cmplx(maxval(real(b)) - sum(real(b)), &
!                                        maxval(imag(b)))
!! Check inverse answer
!    write(*,*) 'Max error C2C INV: ', maxval(abs(a-c))

! Real transform
    ierr = ierr + cufftPlan2D(iplan2,ny,nx,CUFFT_D2Z)
    write(*,*) 'plan2 Real transform:', ierr
    ierr = ierr + cufftPlan2D(iplan3,ny,nx,CUFFT_Z2D)
    write(*,*) 'plan3 Real transform:', ierr

!$acc data copy(r,b)
    !$acc host_data use_device(r,b)
    ierr = ierr + cufftExecD2Z(iplan2,r,b)
    write(*,*) 'plan2 D2Z transform:', ierr
    !$acc end host_data
!$acc end data

    write(fname,'("test_cuFFT_kk.dat")')
    open(20,file=fname)

    do i=1,mnx
        do j=1,ny
            write(20,'(4e15.7)') i*dx,j*dy,abs(b(i,j))*ceff_norm
        enddo
        write(20,*)
    enddo
    close(20)

!$acc data copy(b,q)
    !$acc host_data use_device(b,q)
    ierr = ierr + cufftExecZ2D(iplan3,b,q)
    write(*,*) 'plan3 Z2D transform:', ierr
    !$acc end host_data
!$acc end data

    write(fname,'("test_cuFFT_xy.dat")')
    open(20,file=fname)

    do i=1,nx
        do j=1,ny
            write(20,'(4e15.7)') i*dx,j*dy,r(i,j),q(i,j)*ceff_norm
        enddo
        write(20,*)
    enddo
    close(20)


!    ierr = ierr + cufftDestroy(iplan1)
    ierr = ierr + cufftDestroy(iplan2)
    ierr = ierr + cufftDestroy(iplan3)




    if (ierr.eq.0) then
        print *,"test PASSED"
    else
        print *,"test FAILED"
    endif

end program cufft2dTest
