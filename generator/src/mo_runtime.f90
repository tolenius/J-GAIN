module mo_runtime 
implicit none 
  
public :: abort 
private    
  interface abort
     module procedure abortImpl
  end interface abort

contains
!> @brief stops the program with error message.
!!
subroutine abortImpl(msg)

implicit none
  ! Arguments  
  character(len=*)  , intent(in) :: msg
#ifdef MPI
  include 'mpif.h' 
  integer :: ierror
#endif  
  ! Body 
  write(*,*) "Abort, Message:"
  write(*,*) trim(msg)
  write(*,*) "Program will stop..."

  
#ifdef MPI
  call MPI_ABORT(MPI_COMM_WORLD,ierror)
#endif 

  stop 1  
end subroutine abortImpl



end module mo_runtime 