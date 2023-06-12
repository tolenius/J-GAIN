module mo_runtime 
implicit none 
  
public :: abort 
private    
  interface abort
     module procedure abortImpl
  end interface abort

contains
!> @brief stop the application with message.
!!
subroutine abortImpl(msg)

implicit none
  ! Arguments  
  character(len=*)  , intent(in) :: msg
  ! Body 
  write(*,*) "Abort, Message:"
  write(*,*) trim(msg)
  write(*,*) "Program will stop..."

  stop 1  
end subroutine abortImpl



end module mo_runtime 