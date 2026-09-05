!
! The error type the grid modules report through.
!
! Those modules came from metalquicha, where `error_t` carries a code, a
! message, and a call stack. This is the subset they use -- `set`,
! `has_error`, `get_message`, `clear` -- with the same names, so the code
! could move with nothing changed but its module prefix. The stack is not
! here because nothing in terco pushes onto it.
!
! terco's own API has no error path: it is handed a basis and a density
! and returns integrals, and a wrong argument is a programming error. The
! grid is the exception because a grid request is user input -- a level, a
! count, a scheme -- and refusing a bad one is part of the contract.
!
module trc_error
   implicit none
   private

   public :: error_t
   public :: SUCCESS, ERROR_GENERIC, ERROR_VALIDATION

   integer, parameter :: SUCCESS = 0
   integer, parameter :: ERROR_GENERIC = 1
   integer, parameter :: ERROR_VALIDATION = 4  ! metalquicha's value, kept

   type :: error_t
      integer :: code = SUCCESS
      character(len=:), allocatable :: message
   contains
      procedure :: has_error => error_has_error
      procedure :: set => error_set
      procedure :: clear => error_clear
      procedure :: get_code => error_get_code
      procedure :: get_message => error_get_message
   end type error_t

contains

   pure function error_has_error(this) result(has_err)
      class(error_t), intent(in) :: this
      logical :: has_err
      has_err = (this%code /= SUCCESS)
   end function error_has_error

   pure subroutine error_set(this, code, message)
      class(error_t), intent(inout) :: this
      integer, intent(in) :: code
      character(len=*), intent(in) :: message
      this%code = code
      this%message = trim(message)
   end subroutine error_set

   pure subroutine error_clear(this)
      class(error_t), intent(inout) :: this
      this%code = SUCCESS
      if (allocated(this%message)) deallocate (this%message)
   end subroutine error_clear

   pure function error_get_code(this) result(code)
      class(error_t), intent(in) :: this
      integer :: code
      code = this%code
   end function error_get_code

   pure function error_get_message(this) result(message)
      class(error_t), intent(in) :: this
      character(len=:), allocatable :: message
      if (allocated(this%message)) then
         message = this%message
      else
         message = ""
      end if
   end function error_get_message

end module trc_error
