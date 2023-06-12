module mo_gen_profile
  use mo_lookup_config, only: gen_config_profile_type
  use mo_kind, only : rk4,rk8
  use get_acdc_J, only : acdc_plugin
  implicit none 
  
  public genAndWriteProfile
  
  private 
  interface genAndWriteProfile
    module procedure genAndWriteProfileImp
  end interface genAndWriteProfile
  
contains
 
 
subroutine genAndWriteProfileImp(gen_config_profile)
implicit none 
  ! Arguments 
  type(gen_config_profile_type), intent(in) :: gen_config_profile
    call genAndWriteProfileFunc1(gen_config_profile)

end subroutine genAndWriteProfileImp

!> @brief the main function that calls every thing else. it generates the lookup table
subroutine genAndWriteProfileFunc1(gen_config_profile)
  use mo_utils, only: advanceLoopsBy,getNewUnit 
  use mo_debug 
  use mo_runtime, only: abort

  implicit none 
  ! Arguments 
  type(gen_config_profile_type), intent(in) :: gen_config_profile
  
  !Local Variables 
  real(kind=rk8) :: t_sim = 60 
  real(kind=rk8) :: t_tot = 0    
  real(kind=rk8) :: diameter_acdc
  real(kind=rk8) :: nuc_by_charge(3)
  real(kind=rk8) :: j_acdc
  real(kind=rk8) :: CS,T,IPR,RH
  real(kind=rk8), allocatable, dimension(:) :: vapours_consentrations

  
  real(kind=rk8) , allocatable, dimension(:) :: allvars
  integer :: allvarsCount
  integer :: k,inu,outu, readstat
  integer :: i, shift 
  character(len=200) :: header 
  !
  integer :: rc ! error code  
      
  ! Body 
 
 
  ! 
  k = 0 
  
  call getNewUnit(inu)
  call getNewUnit(outu)
  open(unit = inu, file = trim(gen_config_profile%inputProfilePath), status="old", iostat=rc)
  if (rc /= 0) stop ("Error: failed to open file: "//trim(gen_config_profile%inputProfilePath))

  open(unit = outu, file = trim(gen_config_profile%outputProfilePath), status="unknown", iostat=rc)
  if (rc /= 0) stop ("Error: failed to open file: "//trim(gen_config_profile%outputProfilePath))

  
  read(inu,'(A)',iostat = readstat)  header 
  write(outu,*)  trim(header) , " j_acdc"
  shift = 1 
  allvarsCount = gen_config_profile%dimCount + shift ! this shift is the first column time for example 
  
 
  if (gen_config_profile%idxRH > 0) then 
     allvarsCount = allvarsCount + 1 
  end if  
  
  allocate(allvars(allvarsCount))
  
  allocate(vapours_consentrations(gen_config_profile%vapours_count))
  do while( readstat == 0 )
    	
    read(inu,*, iostat = readstat) allvars
	
	k = k + 1
	CS   =  allvars(gen_config_profile%idxCS  + shift )
	T    =  allvars(gen_config_profile%idxT   + shift )
	IPR  =  allvars(gen_config_profile%idxIPR + shift )
	if (gen_config_profile%idxRH > 0) then 
	   RH   = allvars(gen_config_profile%idxRH + shift )
	else 
	   RH = -1
	end if  
	
    do i=1,gen_config_profile%vapours_count
	   vapours_consentrations(i) = allvars(gen_config_profile%vapours_indices(i) + shift )
	end do 
    
    call acdc_plugin(gen_config_profile%vapours_names,vapours_consentrations,CS,T,RH,IPR,t_sim,t_tot,j_acdc,diameter_acdc)
    write(outu,*) allvars,j_acdc
	
    

  end do 
  close(inu) 
  close(outu) 
  deallocate(vapours_consentrations)
  deallocate(allvars)
end subroutine genAndWriteProfileFunc1

end module mo_gen_profile 







