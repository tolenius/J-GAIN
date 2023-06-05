module mo_gen_parallel
  use mo_lookup_config, only: gen_config_ranges_type
  use mo_kind, only : dp,rk4,ik8,ik4,rk8
  use get_acdc_J, only : acdc_plugin  
  use mo_mpi
  implicit none 

  public genLookupTable_
  
  private 
  interface genLookupTable_
    module procedure genLookupTableImp
  end interface genLookupTable_
  
contains
 
 
subroutine genLookupTableImp(gen_config_ranges,table)
use mo_runtime, only: abort
implicit none 
  ! Arguments 
  type(gen_config_ranges_type), intent(in) :: gen_config_ranges
  real(kind=rk4), allocatable,  dimension(:,:), intent(out) :: table 
   
   call genLookupTableFunc1(gen_config_ranges,table)

end subroutine genLookupTableImp


subroutine genLookupTableFunc1(gen_config_ranges,table)
  use mo_utils, only: advanceLoopsBy,getNewUnit 
  use mo_debug 
  use mo_config
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
  integer(kind=ik8) :: totalCount,counter,lcounter
  integer :: Dest
  integer :: ierror  
  integer :: source
  integer :: request,dbletype
  integer :: i
  !
  integer :: rc ! error code 
  character(len=500) :: file_name ! temp file name
  
  
  include 'mpif.h'
 

 call MPI_TYPE_CREATE_F90_REAL(rk8, MPI_UNDEFINED, dbletype, ierror) 
  
  totalCount = gen_config_ranges%totalCount
  
  !write(*,*) process_rank
  ! call bcast(totalCount)
  !write(*,*) process_rank, totalCount
  
  allocate(vapours_consentrations(gen_config_ranges%vapours_count))

  if(isMaster) then 
      
  
      if (.not. allocated(table)) allocate(table(totalCount,1))
      

      if(diagOut .and. gen_config_ranges%totalCount > 1.E8)  call abort("Very big totalCount. set diagOut=.false.  in the namelist.")
      
	  if(diagOut) then 
        call getNewUnit(tu)
		file_name = trim(gen_config_ranges%outputDirectory)//"/"//trim(gen_config_ranges%outputFileBase)//".diag"
        open(unit = tu, file = trim(file_name), status="unknown", iostat=rc)
		if (rc /= 0) stop "Error: failed to open file: "//trim(file_name)

      end if       
	  allocate(loopVals(gen_config_ranges%dimCount))
      loopVals = -1 
      k = 0 
      Dest = 0
	  counter = 0 
      do while(.true.)
    	 counter = counter + 1;
		 Dest = Dest + 1
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
	     !write(*,*) "send start " , counter  
         call MPI_SEND(vapours_consentrations, gen_config_ranges%vapours_count, dbletype, Dest, 100, MPI_COMM_WORLD, REQUEST, IERROR)
         call MPI_SEND(CS  , 1, dbletype, Dest, 102, MPI_COMM_WORLD, REQUEST, IERROR)
         call MPI_SEND(T   , 1, dbletype, Dest, 103, MPI_COMM_WORLD, REQUEST, IERROR)
         call MPI_SEND(RH  , 1, dbletype, Dest, 104, MPI_COMM_WORLD, REQUEST, IERROR)
         call MPI_SEND(IPR , 1, dbletype, Dest, 105, MPI_COMM_WORLD, REQUEST, IERROR)
         call MPI_SEND(counter, 1, MPI_INTEGER8, Dest, 106, MPI_COMM_WORLD, REQUEST, IERROR) 
		 !write(*,*) "S ", gen_config_ranges%vapours_count, vapours_consentrations, CS,T,RH,IPR
		 !write(*,*) "send end " , counter    
         if(Dest == (size_of_cluster-1) .or. counter == totalCount) then 
		 
           ! Recieve 
           do source = 1, Dest                                      
             call MPI_RECV(j_acdc, 1, dbletype, source, 110, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
             call MPI_RECV(lcounter, 1, MPI_INTEGER8, source, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
			 table(lcounter,1) =  real(j_acdc)                    
           end do 
           DEST = 0
		   if (debug) then 
	         if(mod(counter,1000) == 0 .or. all(loopVals >= gen_config_ranges%dimensions) ) then 
	           write(*,*) k , loopVals
	           write(*,*) k ,vapours_consentrations ,CS  ,T   ,IPR
             end if 
          end if 
          !if(diagOut) write(tu,'(I9,7E17.8)') k ,vapours_consentrations,CS  ,T   ,IPR , j_acdc
          if(diagOut) write(tu,*) k ,vapours_consentrations,CS  ,T   ,IPR , j_acdc
		     
		   
		 end if 

		 
 
 	
         if (all(loopVals >= gen_config_ranges%dimensions)) exit 
      end do 
	  
	  deallocate(loopVals)
      if(diagOut) close(tu)
	  
    else ! if slave 
		!write(*,*) "TTTT", process_rank
	    counter = -100000 
        do while ( counter + size_of_cluster - 1 <= totalCount)
           call MPI_RECV(vapours_consentrations, gen_config_ranges%vapours_count, dbletype, 0, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
		   call MPI_RECV(CS, 1, dbletype, 0, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
           call MPI_RECV(T,  1, dbletype, 0, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
           call MPI_RECV(RH, 1, dbletype, 0, 104, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
           call MPI_RECV(IPR,1, dbletype, 0, 105, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
           call MPI_RECV(counter, 1, MPI_INTEGER8, 0, 106, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
           !write(*,*) "R ", gen_config_ranges%vapours_count, vapours_consentrations, CS,T,RH,IPR
		   call acdc_plugin(gen_config_ranges%vapours_names,vapours_consentrations,CS,T,RH,IPR,t_sim,t_tot,j_acdc,diameter_acdc)


           call MPI_SEND(j_acdc, 1, dbletype,0, 110, MPI_COMM_WORLD, REQUEST, IERROR)
           call MPI_SEND(counter, 1, MPI_INTEGER8,0, 111, MPI_COMM_WORLD, REQUEST, IERROR)
        end do 
	
    end if 
	
	deallocate(vapours_consentrations)
end subroutine genLookupTableFunc1

end module mo_gen_parallel 