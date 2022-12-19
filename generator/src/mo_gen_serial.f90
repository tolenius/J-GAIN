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
  ! Body 
  if (.not. allocated(table)) allocate(table(gen_config_ranges%totalCount,1))
 
  ! 
  if(diagOut .and. gen_config_ranges%totalCount > 1.E7)  call abort("Very big totalCount. set diagOut=.false. in the namelist.")

  allocate(loopVals(gen_config_ranges%dimCount))
  loopVals = -1 
  k = 0 
  
  if(diagOut) then 
    call getNewUnit(tu)
    open(unit = tu, file = trim(gen_config_ranges%outputDirectory)//"/"//trim(gen_config_ranges%outputFileBase)//".diag",status="unknown")
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











!!!!!!!1234 subroutine create_lookup
!!!!!!!1234 USE second_Precision,  ONLY : dp    ! KPP Numerical type
!!!!!!!1234 use mo_lookup_config
!!!!!!!1234 
!!!!!!!1234 implicit none
!!!!!!!1234 
!!!!!!!1234 logical :: solve_ss = .true.
!!!!!!!1234 real(dp) :: dt = 0 
!!!!!!!1234 real(dp) :: c_acid,c_base,c_org,CS_H2SO4,T
!!!!!!!1234 real(dp),allocatable :: ipr(:)
!!!!!!!1234 real(dp),allocatable :: Jnucl(:)
!!!!!!!1234 real(dp) :: diameter_acdc
!!!!!!!1234 real(dp) :: Nuc_by_charge(3)
!!!!!!!1234 real(dp) :: c_acid_log10,c_base_log10,CS_H2SO4_log10
!!!!!!!1234 real(dp) :: c_acid_d,c_base_d,CS_H2SO4_d,T_d,ipr_d
!!!!!!!1234 integer :: n1,n2,n3,n4,n5
!!!!!!!1234 integer :: i1,i2,i3,i4,i5
!!!!!!!1234 real(dp) :: min1,min2,min3,min4,min5
!!!!!!!1234 real(dp) :: max1,max2,max3,max4,max5
!!!!!!!1234 
!!!!!!!1234 integer :: dimCount = 5 
!!!!!!!1234 integer :: process_Rank, size_Of_Cluster
!!!!!!!1234 integer*8 :: counter,rec_counter,localcounter,counter_l
!!!!!!!1234 integer :: REQUEST
!!!!!!!1234 integer :: IERROR
!!!!!!!1234 integer :: STATUS
!!!!!!!1234 integer :: DEST,SOURCE
!!!!!!!1234 real(dp) c_acid_l  
!!!!!!!1234 real(dp) c_base_l 
!!!!!!!1234 real(dp) CS_H2SO4_l
!!!!!!!1234 real(dp) T_l 
!!!!!!!1234 real(dp) ipr_l 
!!!!!!!1234 real(dp) :: Jnucl_l 
!!!!!!!1234 
!!!!!!!1234 real(dp), allocatable, dimension(:) :: steps 
!!!!!!!1234 integer :: doubletype  
!!!!!!!1234 integer total_count
!!!!!!!1234 
!!!!!!!1234 type(gen_config_ranges_type) :: gen_config_ranges  
!!!!!!!1234   
!!!!!!!1234   
!!!!!!!1234 
!!!!!!!1234 write(*,*) "Create lookup table - Serial"
!!!!!!!1234 
!!!!!!!1234 call init("namelist.gen" , gen_config_ranges)
!!!!!!!1234 
!!!!!!!1234  
!!!!!!!1234 
!!!!!!!1234  
!!!!!!!1234 
!!!!!!!1234  
!!!!!!!1234 c_acid_d    = (max1 - min1)/(n1-1)
!!!!!!!1234 c_base_d    = (max2 - min2)/(n2-1)
!!!!!!!1234 CS_H2SO4_d  = (max3 - min3)/(n3-1)
!!!!!!!1234 T_d         = (max4 - min4)/(n4-1)
!!!!!!!1234 ipr_d       = (max5 - min5)/(n5-1)
!!!!!!!1234  
!!!!!!!1234 
!!!!!!!1234 if (process_Rank .eq.0) then 
!!!!!!!1234    
!!!!!!!1234    write(*,*) "Started: Create descriptor file"
!!!!!!!1234    OPEN(1102,FILE='lookup_M7.descriptor',STATUS='unknown',ACTION='WRITE') 
!!!!!!!1234    
!!!!!!!1234    write(1102,*)    "VARS_COUNT     " 
!!!!!!!1234    write(1102,'(i8)')     1
!!!!!!!1234    write(1102,*)     "Jnucl_log10    " 
!!!!!!!1234    write(1102,*)    "DIMS_COUNT     "  
!!!!!!!1234    write(1102,'(i8)')     dimCount
!!!!!!!1234    write(1102,*)    "DIMS           "
!!!!!!!1234    write(1102,'(i8)')  n1
!!!!!!!1234    write(1102,'(i8)')  n2
!!!!!!!1234    write(1102,'(i8)')  n3
!!!!!!!1234    write(1102,'(i8)')  n4
!!!!!!!1234    write(1102,'(i8)')  n5
!!!!!!!1234    write(1102,*)          "MINVALS,MAXVALS,STEPS" 
!!!!!!!1234    write(1102,'(3f18.7)') min1,max1,c_acid_d
!!!!!!!1234    write(1102,'(3f18.7)') min2,max2,c_base_d
!!!!!!!1234    write(1102,'(3f18.7)') min3,max3,CS_H2SO4_d
!!!!!!!1234    write(1102,'(3f18.7)') min4,max4,T_d
!!!!!!!1234    write(1102,'(3f18.7)') min5,max5,ipr_d
!!!!!!!1234    write(1102,*) 'IS_LOG10'
!!!!!!!1234    write(1102,*) .TRUE.
!!!!!!!1234    write(1102,*) .TRUE.
!!!!!!!1234    write(1102,*) .TRUE.
!!!!!!!1234    write(1102,*) .FALSE.
!!!!!!!1234    write(1102,*) .FALSE.
!!!!!!!1234    write(1102,*) 'BIN FILE NAME'
!!!!!!!1234    write(1102,*) 'lookup_M7.bin'
!!!!!!!1234    write(1102,*)  "VALUES_COUNT   "
!!!!!!!1234    write(1102,'(i9)')   n1*n2*n3*n4*n5
!!!!!!!1234    flush(1102)
!!!!!!!1234    close(1102)
!!!!!!!1234    write(*,*) "Finished: Create descriptor file"
!!!!!!!1234 
!!!!!!!1234 end if 
!!!!!!!1234 !call MPI_FINALIZE(ierror) 
!!!!!!!1234 !stop  
!!!!!!!1234 !!!!!!!!!!!!!!!!!! header section
!!!!!!!1234       
!!!!!!!1234 
!!!!!!!1234 if (process_Rank .eq.0) OPEN(1104,FILE="lookup_M7.txt",STATUS='unknown') 
!!!!!!!1234 if (process_Rank .eq.0) open(1001,file="lookup_M7.bin",status="unknown",form="unformatted",access="direct", recl = 4)
!!!!!!!1234  
!!!!!!!1234 
!!!!!!!1234  
!!!!!!!1234 if (process_Rank .eq. 0 ) then
!!!!!!!1234 ! Master
!!!!!!!1234 counter = 0 
!!!!!!!1234 rec_counter = 0 
!!!!!!!1234 DEST = 0
!!!!!!!1234 
!!!!!!!1234 
!!!!!!!1234 ipr(1) = min5 
!!!!!!!1234 do i5 =2,n5 
!!!!!!!1234     ipr(i5) = ipr(i5-1)  + ipr_d
!!!!!!!1234 end do 
!!!!!!!1234 
!!!!!!!1234 
!!!!!!!1234 c_acid_log10   = min1 - c_acid_d
!!!!!!!1234 do i1= 1,n1
!!!!!!!1234    c_acid_log10 = c_acid_log10 + c_acid_d
!!!!!!!1234    c_acid = 10.d0**c_acid_log10
!!!!!!!1234    !write(*,*)  'c_acid 1' , c_acid
!!!!!!!1234 
!!!!!!!1234    
!!!!!!!1234    c_base_log10   = min2 - c_base_d   
!!!!!!!1234    do i2= 1,n2
!!!!!!!1234       c_base_log10 = c_base_log10 + c_base_d
!!!!!!!1234       c_base = 10.d0**c_base_log10  
!!!!!!!1234       CS_H2SO4_log10 = min3 - CS_H2SO4_d       
!!!!!!!1234       do i3= 1,n3
!!!!!!!1234          CS_H2SO4_log10 = CS_H2SO4_log10 + CS_H2SO4_d
!!!!!!!1234          CS_H2SO4 = 10.d0**CS_H2SO4_log10 
!!!!!!!1234          T = min4 - T_d
!!!!!!!1234          do i4 = 1,n4 
!!!!!!!1234             T = T + T_d
!!!!!!!1234             do i5 =1,n5 
!!!!!!!1234                 counter = counter + 1;
!!!!!!!1234                 !write(*,*) "RRRRRR", counter, i1,i2,i3,i4,i5
!!!!!!!1234 
!!!!!!!1234                 if (counter > 0 ) then 
!!!!!!!1234                   DEST = DEST + 1 
!!!!!!!1234                   !if (mod(counter,1000) .eq. 0 ) write(*,*) "QQQQQ", i1,i2,i3,i4,i5 
!!!!!!!1234                   ! use Isend is better but it makes error on bi nowadays 
!!!!!!!1234                   call MPI_SEND(c_acid, 1, doubletype, DEST, 100, MPI_COMM_WORLD, REQUEST, IERROR)
!!!!!!!1234                   call MPI_SEND(c_base, 1, doubletype, DEST, 101, MPI_COMM_WORLD, REQUEST, IERROR)
!!!!!!!1234                   call MPI_SEND(CS_H2SO4, 1, doubletype, DEST, 102, MPI_COMM_WORLD, REQUEST, IERROR)
!!!!!!!1234                   call MPI_SEND(T, 1, doubletype, DEST, 103, MPI_COMM_WORLD, REQUEST, IERROR)
!!!!!!!1234                   call MPI_SEND(ipr(i5), 1, doubletype, DEST, 104, MPI_COMM_WORLD, REQUEST, IERROR)
!!!!!!!1234                   call MPI_SEND(counter, 1, MPI_INTEGER8, DEST, 105, MPI_COMM_WORLD, REQUEST, IERROR)
!!!!!!!1234                   
!!!!!!!1234                   if(DEST == (size_Of_Cluster-1) .or. counter == total_count) then 
!!!!!!!1234                    ! Recieve 
!!!!!!!1234                     write(100,*) "QQQQQ", DEST, counter, i1,i2,i3,i4,i5
!!!!!!!1234                     do source = 1, DEST                                      
!!!!!!!1234                       call MPI_RECV(Jnucl_l, 1, doubletype, source, 110, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
!!!!!!!1234                       call MPI_RECV(localcounter, 1, MPI_INTEGER8, source, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)                      
!!!!!!!1234                       write(1104,'(i15,E15.7)') localcounter,Jnucl_l ! counter,real(Jnucl)
!!!!!!!1234                       rec_counter = rec_counter + 1
!!!!!!!1234                       write(1001,rec=rec_counter) real(Jnucl_l)
!!!!!!!1234                       flush(1104)
!!!!!!!1234                       flush(1001)                      
!!!!!!!1234                     end do 
!!!!!!!1234                     DEST = 0
!!!!!!!1234                   end if
!!!!!!!1234                 end if    
!!!!!!!1234             end do               
!!!!!!!1234          end do            
!!!!!!!1234       end do 
!!!!!!!1234    end do 
!!!!!!!1234 end do 
!!!!!!!1234   flush(1104)
!!!!!!!1234   flush(1001)
!!!!!!!1234   close(1104)
!!!!!!!1234   close(1001)  
!!!!!!!1234 end if 
!!!!!!!1234 
!!!!!!!1234  
!!!!!!!1234 counter_l=0
!!!!!!!1234 if (process_Rank > 0) then 
!!!!!!!1234    do while ( counter_l + size_Of_Cluster - 1 <= total_count)
!!!!!!!1234    ! SLAVE OR Master 
!!!!!!!1234    call MPI_RECV(c_acid_l, 1, doubletype, 0, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
!!!!!!!1234    !write(*,*) 'recieved ',c_acid_l
!!!!!!!1234    call MPI_RECV(c_base_l, 1, doubletype, 0, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
!!!!!!!1234    call MPI_RECV(CS_H2SO4_l, 1, doubletype, 0, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
!!!!!!!1234    call MPI_RECV(T_l, 1, doubletype, 0, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
!!!!!!!1234    call MPI_RECV(ipr_l, 1, doubletype, 0, 104, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
!!!!!!!1234    call MPI_RECV(counter_l, 1, MPI_INTEGER8, 0, 105, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
!!!!!!!1234    !if (process_Rank .eq. 1) write(*,'(5E15.7)') c_acid_l,c_base_l,CS_H2SO4_l,T_l,ipr_l
!!!!!!!1234    CALL  get_acdc_J(c_acid_l,c_base_l,c_org,CS_H2SO4_l,T_l,ipr_l,dt,solve_ss,Jnucl_l,diameter_acdc,Nuc_by_charge)
!!!!!!!1234    call MPI_SEND(Jnucl_l, 1, doubletype,0, 110, MPI_COMM_WORLD, REQUEST, IERROR)
!!!!!!!1234    call MPI_SEND(counter_l, 1, MPI_INTEGER8,0, 111, MPI_COMM_WORLD, REQUEST, IERROR)
!!!!!!!1234    end do 
!!!!!!!1234 end if  
!!!!!!!1234   deallocate(ipr)
!!!!!!!1234   deallocate(Jnucl) 
!!!!!!!1234   write(*,*) "END    counter" , counter_l , "in process " , process_Rank
!!!!!!!1234   call MPI_Barrier(MPI_COMM_WORLD,IERROR)
!!!!!!!1234   call MPI_FINALIZE(ierror)
!!!!!!!1234   stop 
!!!!!!!1234   
!!!!!!!1234   
!!!!!!!1234 end subroutine create_lookup