

subroutine jacobi2d_sub(n, tol, maxiter, iter, dumax, nthreads)
    ! Implements 2d jacobi iteration on n by n grid
    ! Inputs:  n, tol, maxiter, nthreads
    ! Outputs: iter, dumax
    
    use omp_lib ! #include <omp.h>
    implicit none
    integer, intent(in) :: n, maxiter, nthreads
    integer, intent(out) :: iter

    ! Timing variables
    integer :: clock_start, clock_finish, clock_rate
    real(kind=8) :: elapsed_time


    real(kind=8), intent(in) :: tol
    real(kind=8), intent(out) :: dumax

    ! local variables:
    integer, parameter :: nprint = 100000
    real(kind=8), dimension(:,:), allocatable :: u,uold,f
    ! two dimension array
    real(kind=8), dimension(:), allocatable :: x,y
    ! one dimension array
    real(kind=8) :: h
    integer :: i,j

    ! allocate storage for boundary points too:
    allocate(x(0:n+1), y(0:n+1), u(0:n+1, 0:n+1), uold(0:n+1, 0:n+1), &
             f(0:n+1, 0:n+1))

    open(unit=20, file="solution.txt", status="unknown")
    ! where 'status' may be 'OLD', 'NEW', 'REPLACE', 'SCRATCH' or 'UNKNOWN'. If 'OLD' is specified the file must exist; if 'NEW' the file must not exist; if 'REPLACE' and the file exists it will be deleted before a new file is created; and if 'SCRATCH' the file will be deleted when closed. In general use 'OLD' for input and 'NEW' for output. 'UNKNOWN', to connect an existing file, or to create and connect a new file. If the file exists, it is connected as OLD. If the file does not exist, it is connected as NEW.

    ! 1x1 size square
    ! grid spacing:
    h = 1.d0 / (n+1.d0)

    do i=0,n+1
        ! grid points in x:
        x(i) = i*h

        ! boundary conditions:
        ! bottom boundary y=0:
        u(i,0) = 0.d0      
        ! top boundary y=1:
        if ((x(i) <= 0.3d0) .or. (x(i) >= 0.7d0)) then
           u(i,n+1) = 0.d0
        else
           u(i,n+1) = 1.d0
        endif
    enddo

    do j=0,n+1
        ! grid points in y:
        y(j) = j*h

        ! boundary conditions:
        ! left boundary x=0
        u(0,j) = 0.d0      
        ! right boundary x=1
        if ((y(j) <= 0.3d0) .or. (y(j) >= 0.7d0)) then
           u(n+1,j) = 0.d0
        else
           u(n+1,j) = 1.d0
        endif
        enddo

    do j=1,n
        do i=1,n
            ! source term:
            f(i,j) = 0.d0
            ! initial guess:
            u(i,j) = 0.d0
            enddo
        enddo

    ! tolerance and max number of iterations:
    !print *, "Convergence tolerance: tol = ",tol
    !print *, "Maximum number of iterations: maxiter = ",maxiter

    !$ call omp_set_num_threads(nthreads)
    !$ if (nthreads == 1) then
    !$    print "('Using only 1 thread')"
    !$ else
    !$    print "('Using OpenMP with ',i3,' threads')", nthreads
    !$ endif

    call system_clock(clock_start,clock_rate)

    ! Jacobi iteratation:

    do iter=1,maxiter
        uold = u         ! old values
        dumax = 0.d0
        !$omp parallel do private(i) reduction(max : dumax)
        ! private - only be accessed by a single thread.
        do j=1,n
            do i=1,n
                u(i,j) = 0.25d0*(uold(i-1,j) + uold(i+1,j) + &
                                 uold(i,j-1) + uold(i,j+1) + h**2*f(i,j))
                dumax = max(dumax, abs(u(i,j)-uold(i,j)))
            enddo
        enddo

        if (mod(iter,nprint)==0) then
            print 203, iter, dumax
203         format("After ",i8," iterations, dumax = ",d16.6,/)
            endif
        ! check for convergence:
        if (dumax .lt. tol) exit

        enddo

    call system_clock(clock_finish,clock_rate)
    elapsed_time = real(clock_finish - clock_start,kind=8) / &
                    real(clock_rate,kind=8)

    print "('Converged in ',i6,' iterations')", iter
    ! i6 - 6 characters wide
    print "('Elapsed time: ', f10.3,' seconds')", elapsed_time
    ! f10.3 - floating number, 10 characters wide, 3 decimal places.

    ! print out solution to file:

    do j=0,n+1
        do i=0,n+1
            write(20,222) x(i), y(j), u(i,j)
            enddo
        enddo
222 format(3e20.10)
! three 

    !print *, "Solution is in solution.txt"
    print *, ' '

    close(20)

end subroutine jacobi2d_sub
