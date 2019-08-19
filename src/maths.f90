module maths

  ! Useful mathematics functions for double precision reals
  ! Some complex functionality is built in
  ! Removed functions which depend on LAPACK and FFT libraries

  use types

  implicit none

  public

contains

  function cross(a,b)
    ! Generic cross product function in double precision
    real(dp), dimension(3) :: cross
    real(dp), dimension(3), intent(in) :: a,b
    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  end function cross
  ! ----------------------------------------------------------------------------

  subroutine RotMatrix(angle,vector,matrix)
    ! Generate rotation matrix for <angle> about <vector>
    ! Uses Rodrigues rotation formula
    ! Respective to right hand rule (i.e clockwise looking down vector axis)
    real(dp), intent(in) :: angle,vector(3)                 !Angle and axis
    real(dp), intent(out) :: matrix(3,3)                    !Rotaiton Matrix
    real(dp) :: uvec(3)                                     !Unit vector||axis
    real(dp) :: identity(3,3)                               !Identity matrix
    real(dp) :: W(3,3)                                      !Internal matrix for calculation

    uvec = vector/norm2(vector)
    identity =reshape(dble([1,0,0,0,1,0,0,0,1]),shape(identity))
    W = reshape([0._dp,uvec(3),-1*uvec(2),&
         -1*uvec(3),0._dp,uvec(1),&
         uvec(2),-1*uvec(1),0._dp],shape(W))
    matrix = identity + W*sin(angle)+matmul(W,W)*(1-cos(angle))
  end subroutine RotMatrix
  ! -------------------------------------------------------------------------------

  recursive function factorial(n) result(r)
    ! Calculates N!
    integer,intent(in) :: n                                 !Number N
    integer :: r                                            !N factorial
    if (n == 0) then
       r = 1
    else if (n>0) then
       r = n*factorial(n-1)
    else
       r=0
       write(*,"(A)") 'Negative number passed to factorial function. Zero retuned.'
    end if
  end function factorial
  ! ------------------------------------------------------------------------------

  function linear_interp(x,y,x0) result(r)
    ! Linearly interpolates a funciton y(x) for a value of y0 at the point x0
    ! Assumes the x,y data has been sorted
    real(dp), intent(in) :: x(:),y(:)                       !Discrete function y(x)
    real(dp) :: x0                                          !Desired point x0
    real(dp) :: r                                           !Desired value y(x0)
    integer :: i                                            !Dummy integer

    do i=1,size(x)
       if (x(i)<=x0 .and. x(i+1) >= x0) then
          r = y(i) + (y(i+1)-y(i))*(x0-x(i))/(x(i+1)-x(i))
          exit
       end if
    end do
  end function linear_interp
  ! -------------------------------------------------------------------------------

  function bin_data(x,y,n_bin) result(r)
    ! Recasts discretized function on evenly spaced bins using linear interpolation
    ! Result is the same size as x,y data
    real(dp), intent(in) :: x(:),y(:)                       !Discrete function y(x)
    real(dp) :: r(size(x),2)                                !Return funciton
    integer, intent(in) :: n_bin                            !Number of desired bins
    real(dp) :: x_min                                       !Minimum x value
    real(dp) :: x_max                                       !Maximum x value
    real(dp) :: bin_size                                    !Bin width
    integer :: i                                            !Dummy integer

    r=0.
    x_max = maxval(x)
    x_min = minval(x)
    bin_size = (x_max-x_min)/n_bin

    do i=1,size(x)
       r(i,1) = x_min+bin_size/2.0_dp
       x_min = x_min + bin_size
       r(i,2) = linear_interp(x,y,r(i,1))
    end do
  end function bin_data
  ! -----------------------------------------------------------------------------
end module maths
