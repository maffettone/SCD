module types

  ! The following contains global constants and paramters for use in programs
  ! It also contains some simple functions, which would be universally useful


  implicit none

  public
  ! Double Precision
  integer, parameter, public :: dp = kind(1.d0)
  ! Quadrouple Precision
  integer, parameter, public :: qp = selected_real_kind (32)

  ! Fundemental Constants
  ! Pi, double precision
  real(dp), parameter, public :: pi = acos(-1.d0)
  ! N_av, Avagadros number
  real(dp), parameter, public :: N_av = 6.0221413d23
  ! Boltzmans constant [J/K]
  real(dp), parameter, public :: k_b = 1.3806488d-23
  ! Planks constant, h/2pi [Js]
  real(dp), parameter, public :: hbar = 1.05459d-34
  ! Electron charge [C]
  real(dp), parameter, public :: electron = 1.60219d-19
  ! Dielectric permittivity of free space [F/m]
  real(dp), parameter, public :: epsilon0 = 8.8542d-12
  ! Imaginary unit
  complex :: imag = (0,1)

contains
  function strlen(st) result (r)
    ! Calculates number of occupied elements in a string with trailing blanks
    integer :: r, i                                         !Return and dummy
    character(len=*), intent(in) :: st                      !String
    character(len=1) :: c
    i=len(st)
    c=st(i:i)
    do while(c .eq. ' ')
       i = i-1
       c=st(i:i)
    end do
    r = i
  end function strlen
  ! ----------------------------------------------------------------------------

  elemental subroutine lower_case(word)
    ! convert entire word to lower case
    character (len=*) , intent(in out) :: word
    integer                            :: i,ic,nlen
    nlen = len(word)
    do i=1,nlen
       ic = ichar(word(i:i))
       if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
    end do
  end subroutine lower_case
  ! -----------------------------------------------------------------------------

  elemental subroutine upper_case(word)
    ! convert a word to upper case
    character (len=*), intent (in out) :: word
    integer :: i, ic, nlen
    nlen = len(word)
    do i=1,nlen
       ic = ichar(word(i:i))
       if (ic >= iachar("a") .and. ic <=iachar("z")) word(i:i) = char(ic-32)
    end do
  end subroutine upper_case
  ! ---------------------------------------------------------------------------

  logical function alphabetical(a) result(r)
    ! Determines if a character is a memeber of the alphabet
    character(len=1), intent(in) :: a
    character(len=26) :: alpha, alpha2
    integer :: i
    alpha = 'abcdefghijklmnopqrstuvwxyz'
    alpha2 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    r= .false.
    do i=1,26
       if (a==alpha(i:i) .or. a==alpha2(i:i)) then
          r = .true.
       end if
    end do
  end function alphabetical
  ! -----------------------------------------------------------------------------

end module types
