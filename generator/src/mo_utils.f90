module mo_utils
implicit none 

  public :: advanceLoopsBy
  public :: getNewUnit
private 


  interface advanceLoopsBy
     module procedure advanceLoopsByOneImpl
     module procedure advanceLoopsByImpl
  end interface advanceLoopsBy 
  integer :: unitCounter = 3412
  interface getNewUnit
     module procedure getNewUnitImpl
  end interface getNewUnit  
contains 



subroutine advanceLoopsByOneImpl(loopsMax,loopVars)
! It also do the initial step 
use mo_debug 
use mo_runtime, only : abort 
implicit none 
! Arguments
integer, dimension(:), intent(in) :: loopsMax
integer, dimension(:), intent(inout) :: loopVars
! Local Variables 
integer :: n,i 
logical :: hold 
  
  if ( all(loopVars >= loopsMax)) then 
   	call abort("Loop is Over")  
  end if

  
  if (any(loopVars<1)) then 
    loopVars  = 1 
	if (debug) write(*,*) "init loop"
	return 
  end if 
  
  
  n = size(loopsMax,1)
 
  do i=n,1,-1
   hold = .false.
   loopVars(i) = loopVars(i) + 1 
   if (loopVars(i) > loopsMax(i) ) then 
     loopVars(i) = 1 
	 hold = .true.
   end if  
   if ( .not. hold) return 
  end do 
end subroutine advanceLoopsByOneImpl

subroutine advanceLoopsByImpl(loopsMax,by,loopVars)
! It also do the initial step 
use mo_debug 
use mo_runtime, only : abort 
implicit none 
! Arguments
integer, dimension(:), intent(in) :: loopsMax
integer, intent(in) :: by
integer, dimension(:), intent(inout) :: loopVars
! Local Variables 
integer :: n,i,j,lby
logical :: hold 

  
  if ( all(loopVars >= loopsMax)) then 
   	call abort("Loop is Over")  
  end if
  lby = by 
  if (any(loopVars<1)) then 
    loopVars  = 1 
    lby  = by - 1 
  end if 

  n = size(loopsMax,1)
  do j=1,lby
     do i=n,1,-1
      hold = .false.
      loopVars(i) = loopVars(i) + 1 
      if (loopVars(i) > loopsMax(i) ) then 
        loopVars(i) = 1 
	    hold = .true.
      end if  
      if ( .not. hold) exit  
     end do
     if ( all(loopVars >= loopsMax)) return 	 
  end do 	 
end subroutine advanceLoopsByImpl

subroutine getNewUnitImpl(u)
implicit none 
integer, intent(out) :: u 
  u = unitCounter
  unitCounter = unitCounter + 1
end subroutine getNewUnitImpl 

end module mo_utils
