module mo_lookup_config
  use mo_kind, only: dp
  use mo_debug
implicit none 


public gen_config_ranges_type
public initConf

public FORMAT_BIN
public FORMAT_NC3


type steps_type
    real(dp), dimension(:), allocatable :: steps ! index is 1, dimensions(i) in gen_config_ranges_type
end type steps_type

integer, parameter :: fillValue = -999999


type gen_config_profile_type
    character(len=20), dimension(:), allocatable :: varNames
	logical, dimension(:), allocatable :: isVapour
    character(len=20), dimension(:), allocatable :: outVarNames
    integer :: idxCS  	
    integer :: idxT   	
    integer :: idxRH
    integer :: idxIPR
	
    integer :: depCount  = -1
    integer :: vapours_count
	integer, dimension(:), allocatable :: vapours_indices
	character(len=11), dimension(:), allocatable :: vapours_names
	
    ! from shared namelist 
    integer :: dimCount  = -1
	
	logical :: isActiveMode

    character(len=500) :: inputProfilePath
    character(len=500) :: outputProfilePath
	
	
end  type gen_config_profile_type


type gen_config_ranges_type

    ! from lookup_gen namelist independent 
    character(len=20), dimension(:), allocatable :: varNames
    character(len=20), dimension(:), allocatable :: units
    integer, dimension(:), allocatable :: dimensions
    logical, dimension(:), allocatable :: isLog10    
	logical, dimension(:), allocatable :: isVapour

    real(dp), dimension(:), allocatable :: minVals
    real(dp), dimension(:), allocatable :: maxVals
    type(steps_type), dimension(:), allocatable  :: varSteps ! index is 1,dimCount
	integer :: outputFormat  = 0 ! 0: bin, 1: NC3 
    character(len=200) :: outputFileBase
    character(len=500) :: outputDirectory

    ! from shared namelist 
    integer :: dimCount  = -1
	! 
	! Calculated Fields - Dependent 
    integer*8 :: totalCount  = -1
    integer :: depCount  = -1
    integer :: vapours_count
	integer, dimension(:), allocatable :: vapours_indices
	character(len=11), dimension(:), allocatable :: vapours_names
	
    character(len=20), dimension(:), allocatable :: outVarNames
    character(len=20), dimension(:), allocatable :: outUnits
	
    integer :: idxCS  	
    integer :: idxT   	
    integer :: idxRH
    integer :: idxIPR
	! Profile 
	logical :: isActiveMode

	
end  type gen_config_ranges_type

   integer, parameter :: FORMAT_BIN=0
   integer, parameter :: FORMAT_NC3=1
contains

!> @brief initlizaes the generator from namelist
!! 
!! @param[out]   gen_config_ranges    ranges and nature of the independent variables 
!! @param[out]   gen_config_profile   informatioj about the position of the independent variables in the namelist 
subroutine initConf(namelistFile, gen_config_ranges, gen_config_profile)
use mo_runtime, only: abort
use mo_utils, only: getNewUnit
use mo_config, only: canRead
implicit none 
    
	
	
    character(len=*) , intent(in)  :: namelistFile
    type(gen_config_ranges_type) , intent(inout) :: gen_config_ranges
    type(gen_config_profile_type), intent(inout) :: gen_config_profile 
    ! Local Variables 
    integer:: indepCount, u 
    character(len=20), dimension(:), allocatable, target :: varNames
    character(len=20), dimension(:), allocatable, target :: units
    integer, dimension(:), allocatable, target :: dimensions
    logical, dimension(:), allocatable, target :: isLog10
    logical, dimension(:), allocatable, target :: isVapour
    real(dp), dimension(:), allocatable, target :: minVals
    real(dp), dimension(:), allocatable, target :: maxVals
	
	logical :: isProfileMode
    character(len=200) :: inputProfilePath
    character(len=200) :: outputProfilePath
	
    ! from shared namelist 
    integer :: dimCount=1 
    integer :: depCount = 1
	! 
	! dep 
	integer :: outputFormat
    character(len=200) :: outputFileBase
    character(len=500) :: outputDirectory	
    character(len=20), dimension(:), allocatable :: outVarNames
    character(len=20), dimension(:), allocatable :: outUnits
	!
	integer :: i,j,k
	real(dp) :: step
	!
	integer :: rc ! error code  
    
    ! NameLists
    namelist /shared/  indepCount,depCount,debug,diagOut,isProfileMode
    
 
	 
    namelist /independent/    &
     & varNames              ,&   
     & units                 ,&   
     & dimensions            ,&   
     & isLog10               ,&   
     & isVapour              ,&   
     & minVals               ,&   
     & maxVals      
    
    namelist /dependent/     &
     & outVarNames          ,&   
     & outUnits             ,&    
     & outputFormat         ,&   
     & outputFileBase       ,&   
     & outputDirectory      

	namelist /profile/       &
	  & varNames            ,&
	  & isVapour            ,&
	  & inputProfilePath    ,&
	  & outputProfilePath   
	  
    ! Body 
	
    !if (.not. canRead) return  allow all to read 


    call finalize(gen_config_ranges)
	
	 
    call getNewUnit(u)
	write (*,*) "Opening file: "//trim(namelistFile)
    open(unit=u, file=trim(namelistFile), status="old", iostat=rc)
	if (rc /= 0) stop 'Error: failed to open file: ' // trim(namelistFile)

    isProfileMode = .false. 
    rewind(u)
    read(u, nml=shared)    
    gen_config_ranges%dimCount      = indepCount
    gen_config_ranges%depCount      = depCount
    gen_config_profile%dimCount      = indepCount
    gen_config_profile%depCount      = depCount
	gen_config_profile%isActiveMode = .false.
	gen_config_ranges%isActiveMode  = .false.	
    
    ! At least 4 independent variable CS,T,IPR and one vapour
	if (indepCount < 4 ) call abort("less than 4 variables. At least 4 independent variable CS,T,IPR,... and one vapour")	
	
	
	
   
	
	if (isProfileMode) then 
	!!!!!!!!!!!!!!!!!!!!!!! Profile Mode configs !!!!!!!!!!!!!!!!!!!!!!!!!
        gen_config_profile%isActiveMode      = .True.
        
		allocate(varNames(indepCount))
        allocate(isVapour(indepCount))        
		
		rewind(u)
        read(u, nml=profile) 

        gen_config_profile%varNames            = varNames   
        gen_config_profile%isVapour            = isVapour  
		gen_config_profile%inputProfilePath    = trim(inputProfilePath)
        gen_config_profile%outputProfilePath   = trim(outputProfilePath)

       ! 
        allocate(outVarNames(depCount))
       ! 

	    
  
	    

        rewind(u)
        read(u, nml=dependent)
        close(u)
        gen_config_profile%outVarNames      = outVarNames      
	    
	    ! Find Indexes in varNames 
        gen_config_profile%idxCS       = -1 
        gen_config_profile%idxT   	 = -1 
        gen_config_profile%idxRH       = -1 
        gen_config_profile%idxIPR 	 = -1 
	    do i=1,indepCount
	       if (trim(varNames(i)) == "CS"  ) gen_config_profile%idxCS   = i
	       if (trim(varNames(i)) == "T"   ) gen_config_profile%idxT    = i
	       if (trim(varNames(i)) == "RH " ) gen_config_profile%idxRH   = i
	       if (trim(varNames(i)) == "IPR" ) gen_config_profile%idxIPR  = i
	    end do 
		
	    ! Check mandatory variables 
	    if ( gen_config_profile%idxCS == -1 )  call abort("CS is not one of the variables" )
	    if ( gen_config_profile%idxT == -1 )   call abort("Temperature (T) is not one of the variables" )
	    if ( gen_config_profile%idxIPR == -1 ) call abort("idxIPR is not one of the variables" )
	    
        ! Count number of vapours
        gen_config_profile%vapours_count = 0
        do i=1,gen_config_profile%dimCount
            if (gen_config_profile%isVapour(i)) then 
  	            gen_config_profile%vapours_count = gen_config_profile%vapours_count + 1
  	        end if
        end do 
		
        ! find vapurs indices 
        allocate (gen_config_profile%vapours_indices(gen_config_profile%vapours_count))
        allocate (gen_config_profile%vapours_names(gen_config_profile%vapours_count))
        !allocate (vapours_consentrations(vapours_count))
        k = 0
        do i=1,gen_config_profile%dimCount
           if (gen_config_profile%isVapour(i)) then 
  	        k = k + 1
  	        gen_config_profile%vapours_indices(k)  = i
  	    	gen_config_profile%vapours_names(k)(:) = gen_config_profile%varNames(i)
  	       end if
        end do 	
	end if 
	
	
	
	
	if (.not. isProfileMode) then 
	!!!!!!!!!!!!!!!!!!!!!!! Not Profile Mode configs !!!!!!!!!!!!!!!!!!!!!!!!!

	    gen_config_ranges%isActiveMode = .True.
        ! 
        allocate(varNames(indepCount))
        allocate(units(indepCount))
        allocate(dimensions(indepCount))
        allocate(isLog10(indepCount))
        allocate(isVapour(indepCount))
        allocate(minVals(indepCount))
        allocate(maxVals(indepCount))
	    
        allocate(outVarNames(depCount))
        allocate(outUnits(depCount))
	    
	    
	    
        ! 
        rewind(u)
        read(u, nml=independent) 
	    
        gen_config_ranges%varNames        = varNames   
        gen_config_ranges%units           = units      
        gen_config_ranges%dimensions      = dimensions 
        gen_config_ranges%isLog10         = isLog10    
        gen_config_ranges%isVapour        = isVapour    
        gen_config_ranges%minVals         = minVals    
        gen_config_ranges%maxVals         = maxVals    
	    
	    
	    
        rewind(u)
        read(u, nml=dependent)
        close(u)
        gen_config_ranges%outVarNames      = outVarNames      
        gen_config_ranges%outUnits         = outUnits      
        gen_config_ranges%outputFormat     = outputFormat      
        gen_config_ranges%outputFileBase   = outputFileBase  
        gen_config_ranges%outputDirectory  = outputDirectory 
	    
        
	    
	    
	    ! Find Indexes in varNames 
        gen_config_ranges%idxCS       = -1 
        gen_config_ranges%idxT   	 = -1 
        gen_config_ranges%idxRH       = -1 
        gen_config_ranges%idxIPR 	 = -1 
	    do i=1,indepCount
	       if (trim(varNames(i)) == "CS"  ) gen_config_ranges%idxCS   = i
	       if (trim(varNames(i)) == "T"   ) gen_config_ranges%idxT    = i
	       if (trim(varNames(i)) == "RH" ) gen_config_ranges%idxRH   = i
	       if (trim(varNames(i)) == "IPR" ) gen_config_ranges%idxIPR  = i
	    end do 
	    ! Check mandatory variables 
	    if ( gen_config_ranges%idxCS == -1 ) call abort("CS is not one of the variables" )
	    if ( gen_config_ranges%idxT == -1 ) call abort("Temperature (T) is not one of the variables" )
	    if ( gen_config_ranges%idxIPR == -1 ) call abort("idxIPR is not one of the variables" )
	    
	    
	    ! check variables 
	    do i=1,indepCount
	       if (maxVals(i) < minVals(i)) call abort("Check min and max values of the variables. min vale is larger !!!" )
	       if ((maxVals(i) == minVals(i)) .and. dimensions(i) > 1 )  call abort("dimension have to be 1 when min and max values are equal" )
	       if ((maxVals(i) > minVals(i)) .and. dimensions(i) < 2 )   call abort("dimension have to be more then 2 when variables minVal < maxVal " )
	    end do 
	    
        
        
        gen_config_ranges%totalCount = 1 
        do i=1,indepCount
           gen_config_ranges%totalCount = gen_config_ranges%totalCount * dimensions(i)
        end do 
        
   !     Generating steps 
        allocate(gen_config_ranges%varSteps(indepCount))
	    
        do i=1,indepCount
	       allocate(gen_config_ranges%varSteps(i)%steps(dimensions(i)))
	       if (dimensions(i) > 1) then 
               step = (gen_config_ranges%maxVals(i) - gen_config_ranges%minVals(i))/(dimensions(i)-1)
	       else
	           step = 0 
	       end if 
	       gen_config_ranges%varSteps(i)%steps(1) = gen_config_ranges%minVals(i)
	       do j=2,dimensions(i)
	          gen_config_ranges%varSteps(i)%steps(j) = gen_config_ranges%varSteps(i)%steps(j-1) + step 
	       end do 
	       if (gen_config_ranges%isLog10(i)) then 
	       
	       do j=1,dimensions(i)
	          gen_config_ranges%varSteps(i)%steps(j) = (10.0_dp)**gen_config_ranges%varSteps(i)%steps(j) 
	       end do 
	       end if 
	       
        end do 
        ! Count number of vapours
        gen_config_ranges%vapours_count = 0
        do i=1,gen_config_ranges%dimCount
           if (gen_config_ranges%isVapour(i)) then 
  	        gen_config_ranges%vapours_count = gen_config_ranges%vapours_count + 1
  	       end if
        end do 
        ! find vapurs indices 
        allocate (gen_config_ranges%vapours_indices(gen_config_ranges%vapours_count))
        allocate (gen_config_ranges%vapours_names(gen_config_ranges%vapours_count))
        !allocate (vapours_consentrations(vapours_count))
        k = 0
        do i=1,gen_config_ranges%dimCount
           if (gen_config_ranges%isVapour(i)) then 
  	        k = k + 1
  	        gen_config_ranges%vapours_indices(k)  = i
  	    	gen_config_ranges%vapours_names(k)(:) = gen_config_ranges%varNames(i)
  	     end if
        end do 	
	end if 

  
end subroutine initConf 

subroutine finalize(gen_config_ranges)
implicit none 
    ! Arguments
    type(gen_config_ranges_type) , intent(inout) :: gen_config_ranges
    ! Local Variables 
    integer :: i 	
    ! body 
    if (allocated(gen_config_ranges%varNames))           deallocate(gen_config_ranges%varNames)
    if (allocated(gen_config_ranges%units))              deallocate(gen_config_ranges%units)
    if (allocated(gen_config_ranges%dimensions))         deallocate(gen_config_ranges%dimensions)
    if (allocated(gen_config_ranges%isLog10))            deallocate(gen_config_ranges%isLog10)
    if (allocated(gen_config_ranges%isVapour))           deallocate(gen_config_ranges%isVapour)
    if (allocated(gen_config_ranges%minVals))            deallocate(gen_config_ranges%minVals)
    if (allocated(gen_config_ranges%maxVals))            deallocate(gen_config_ranges%maxVals)
    if (allocated(gen_config_ranges%vapours_indices))     deallocate(gen_config_ranges%vapours_indices)
    if (allocated(gen_config_ranges%vapours_names))      deallocate(gen_config_ranges%vapours_names)
    
 
    do i=1,gen_config_ranges%dimCount
       if (allocated(gen_config_ranges%varSteps(i)%steps)) then 
          deallocate(gen_config_ranges%varSteps(i)%steps)
	   end if 
    end do 
    
    if (allocated(gen_config_ranges%varSteps))  deallocate(gen_config_ranges%varSteps)


    gen_config_ranges%dimCount = -1 
    

end subroutine finalize 



end  module mo_lookup_config


