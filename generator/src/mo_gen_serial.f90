module mo_gen_serial 
  use mo_lookup_config, only: gen_config_ranges_type
  use mo_kind, only : rk4,rk8
  use get_acdc_J, only : acdc_plugin
  implicit none 
  
  public genLookupTable_
  
  private 
  interface genLookupTable_
    module procedure genLookupTableImp
  end interface genLookupTable_
  
contains
 

subroutine genLookupTableImp(gen_config_ranges,table)
implicit none 
  ! Arguments 
  type(gen_config_ranges_type), intent(in) :: gen_config_ranges
  real(kind=rk4), allocatable,  dimension(:,:), intent(out) :: table 

    call genLookupTableFunc1(gen_config_ranges,table)

end subroutine genLookupTableImp

!> @brief calculates the lookup table values in serial mode
!!
!! @param[in]   gen_config_ranges   defines  the independent variables  sapace 
!! @param[out]   table   the resulting table of values. here it is Nucleation Rate.
subroutine genLookupTableFunc1(gen_config_ranges,table)
  use mo_utils, only: advanceLoopsBy,getNewUnit 
  use mo_debug 
  use mo_runtime, only: abort

  implicit none 
  ! Arguments 
  type(gen_config_ranges_type), intent(in) :: gen_config_ranges
  real(kind=rk4), allocatable,  dimension(:,:), intent(out) :: table 
  
  !Local Variables 
  real(kind=rk8) :: t_sim = 60 
  real(kind=rk8) :: t_tot = 0    
  real(kind=rk8) :: diameter_acdc
  real(kind=rk8) :: nuc_by_charge(3)
  real(kind=rk8) :: j_acdc
  real(kind=rk8) :: CS,T,IPR,RH
  real(kind=rk8), allocatable, dimension(:) :: vapours_consentrations
   
  integer, allocatable, dimension(:) :: loopVals
  integer :: k,tu
  integer :: i
  !
  integer :: rc ! error code 
  character(len=500) :: file_name ! temp file name   
  ! Body 
  if (.not. allocated(table)) allocate(table(gen_config_ranges%totalCount,1))
 
  ! 
  if(diagOut .and. gen_config_ranges%totalCount > 1.E7)  call abort("Very big totalCount. set diagOut=.false. in the namelist.")

  allocate(loopVals(gen_config_ranges%dimCount))
  loopVals = -1 
  k = 0 
  
  if(diagOut) then 
    call getNewUnit(tu)
	file_name = trim(gen_config_ranges%outputDirectory)//"/"//trim(gen_config_ranges%outputFileBase)//".diag"
    open(unit = tu, file = trim(file_name),&
	& status="unknown", iostat=rc)	
	if (rc /= 0) stop "Error: failed to open file: "//trim(file_name)
  end if 
  

  allocate(vapours_consentrations(gen_config_ranges%vapours_count))
  do while(.true.)
    	
    call advanceLoopsBy(gen_config_ranges%dimensions,loopVals)
	k = k + 1
	CS   = gen_config_ranges%varSteps(gen_config_ranges%idxCS  )%steps(loopVals(gen_config_ranges%idxCS  ))
	T    = gen_config_ranges%varSteps(gen_config_ranges%idxT   )%steps(loopVals(gen_config_ranges%idxT   ))
	IPR  = gen_config_ranges%varSteps(gen_config_ranges%idxIPR )%steps(loopVals(gen_config_ranges%idxIPR ))
	if (gen_config_ranges%idxRH > 0) then 
	   RH   = gen_config_ranges%varSteps(gen_config_ranges%idxRH  )%steps(loopVals(gen_config_ranges%idxRH  ))
	else 
	   RH = -1
	end if  
    
    do i=1,gen_config_ranges%vapours_count
	   vapours_consentrations(i) = gen_config_ranges%varSteps(gen_config_ranges%vapours_indices(i))%steps(loopVals(gen_config_ranges%vapours_indices(i)))
	end do 
	
    call acdc_plugin(gen_config_ranges%vapours_names,vapours_consentrations,CS,T,RH,IPR,t_sim,t_tot,j_acdc,diameter_acdc)
	!    acdc_plugin(names_vapor,                     c_vapor,         cs_ref,temp,rh,ipr,t_sim,t_tot,j_acdc,diameter_acdc)
    table(k,1) = real(j_acdc)
	if(debug  .and. (mod(k,1000) == 0 .or. all(loopVals >= gen_config_ranges%dimensions) )) then 
	  write(*,*) k , loopVals
	  write(*,*) k ,vapours_consentrations,CS,T,IPR
    end if 
    if(diagOut) write(tu,'(I9,7E17.8)') k ,vapours_consentrations ,CS  ,T   ,IPR, j_acdc 	
    
	if (all(loopVals >= gen_config_ranges%dimensions)) exit 
	
  end do 
   if(diagOut) close(tu) 
   deallocate(vapours_consentrations)
end subroutine genLookupTableFunc1

end module mo_gen_serial 