subroutine create_lookup_mpi
USE second_Precision,  ONLY : dp    ! KPP Numerical type
use config
implicit none

logical :: solve_ss
real(dp) :: dt
real(dp) :: c_acid,c_base,c_org,CS_H2SO4,T
real(dp),allocatable :: ipr(:)
real(dp),allocatable :: Jnucl(:)
real(dp) :: diameter_acdc
real(dp) :: Nuc_by_charge(3)
real(dp) :: c_acid_log10,c_base_log10,CS_H2SO4_log10
real(dp) :: c_acid_d,c_base_d,CS_H2SO4_d,T_d,ipr_d
integer :: n1,n2,n3,n4,n5
integer :: i1,i2,i3,i4,i5
real(dp) :: min1,min2,min3,min4,min5
real(dp) :: max1,max2,max3,max4,max5

integer :: dimCount = 5 
integer :: process_Rank, size_Of_Cluster
integer*8 :: counter,rec_counter,localcounter,counter_l
integer :: REQUEST
integer :: IERROR
integer :: STATUS
integer :: DEST,SOURCE
real(dp) c_acid_l  
real(dp) c_base_l 
real(dp) CS_H2SO4_l
real(dp) T_l 
real(dp) ipr_l 
real(dp) :: Jnucl_l 
integer :: doubletype  
integer total_count
include 'mpif.h'
 

write(*,*) "Create lookup table"

call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size_Of_Cluster, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, process_Rank, ierror)


CALL MPI_Type_create_f90_real(14,300, doubletype , ierror )

!!!!!!!!!!!!!!!!! 



solve_ss = .true. 
dt = 0
c_org = 0.0

!      n1,      n2,         n3,              
!Log10 c_acid	c_base	    CS_H2SO4	    
!min   0009.919	0014.867	-0012.536		 
!max   0014.018	0016.795	-0001.470		 


n1 = 50
n2 = 50
n3 = 50
n4 = 50  
n5 = 50
min1 = 10.D0; min2 =  7.D0; min3 =  -6.D0
max1 = 14.D0; max2 = 18.D0; max3 =  0.0D0
 
total_count=n1*n2*n3*n4*n5
 
! T ipr 
min4 = 245_dp  ; min5 = 1500000_dp 
max4 = 300_dp  ; max5 = 3500000_dp 

allocate(ipr(n5))
allocate(Jnucl(n5))
 
c_acid_d    = (max1 - min1)/(n1-1)
c_base_d    = (max2 - min2)/(n2-1)
CS_H2SO4_d  = (max3 - min3)/(n3-1)
T_d         = (max4 - min4)/(n4-1)
ipr_d       = (max5 - min5)/(n5-1)
 

if (process_Rank .eq.0) then 
   
   write(*,*) "Started: Create descriptor file"
   OPEN(1102,FILE='lookup_M7.descriptor',STATUS='unknown',ACTION='WRITE') 
   
   write(1102,*)    "VARS_COUNT     " 
   write(1102,'(i8)')     1
   write(1102,*)     "Jnucl_log10    " 
   write(1102,*)    "DIMS_COUNT     "  
   write(1102,'(i8)')     dimCount
   write(1102,*)    "DIMS           "
   write(1102,'(i8)')  n1
   write(1102,'(i8)')  n2
   write(1102,'(i8)')  n3
   write(1102,'(i8)')  n4
   write(1102,'(i8)')  n5
   write(1102,*)          "MINVALS,MAXVALS,STEPS" 
   write(1102,'(3f18.7)') min1,max1,c_acid_d
   write(1102,'(3f18.7)') min2,max2,c_base_d
   write(1102,'(3f18.7)') min3,max3,CS_H2SO4_d
   write(1102,'(3f18.7)') min4,max4,T_d
   write(1102,'(3f18.7)') min5,max5,ipr_d
   write(1102,*) 'IS_LOG10'
   write(1102,*) .TRUE.
   write(1102,*) .TRUE.
   write(1102,*) .TRUE.
   write(1102,*) .FALSE.
   write(1102,*) .FALSE.
   write(1102,*) 'BIN FILE NAME'
   write(1102,*) 'lookup_M7.bin'
   write(1102,*)  "VALUES_COUNT   "
   write(1102,'(i9)')   n1*n2*n3*n4*n5
   flush(1102)
   close(1102)
   write(*,*) "Finished: Create descriptor file"

end if 
!call MPI_FINALIZE(ierror) 
!stop  
!!!!!!!!!!!!!!!!!! header section
      

if (process_Rank .eq.0) OPEN(1104,FILE="lookup_M7.txt",STATUS='unknown') 
if (process_Rank .eq.0) open(1001,file="lookup_M7.bin",status="unknown",form="unformatted",access="direct", recl = 4)
 

 
if (process_Rank .eq. 0 ) then
! Master
counter = 0 
rec_counter = 0 
DEST = 0


ipr(1) = min5 
do i5 =2,n5 
    ipr(i5) = ipr(i5-1)  + ipr_d
end do 


c_acid_log10   = min1 - c_acid_d
do i1= 1,n1
   c_acid_log10 = c_acid_log10 + c_acid_d
   c_acid = 10.d0**c_acid_log10
   !write(*,*)  'c_acid 1' , c_acid

   
   c_base_log10   = min2 - c_base_d   
   do i2= 1,n2
      c_base_log10 = c_base_log10 + c_base_d
      c_base = 10.d0**c_base_log10  
      CS_H2SO4_log10 = min3 - CS_H2SO4_d       
      do i3= 1,n3
         CS_H2SO4_log10 = CS_H2SO4_log10 + CS_H2SO4_d
         CS_H2SO4 = 10.d0**CS_H2SO4_log10 
         T = min4 - T_d
         do i4 = 1,n4 
            T = T + T_d
            do i5 =1,n5 
                counter = counter + 1;
                !write(*,*) "RRRRRR", counter, i1,i2,i3,i4,i5

                if (counter > 0 ) then 
                  DEST = DEST + 1 
                  !if (mod(counter,1000) .eq. 0 ) write(*,*) "QQQQQ", i1,i2,i3,i4,i5 
                  ! use Isend is better but it makes error on bi nowadays 
                  call MPI_SEND(c_acid, 1, doubletype, DEST, 100, MPI_COMM_WORLD, REQUEST, IERROR)
                  call MPI_SEND(c_base, 1, doubletype, DEST, 101, MPI_COMM_WORLD, REQUEST, IERROR)
                  call MPI_SEND(CS_H2SO4, 1, doubletype, DEST, 102, MPI_COMM_WORLD, REQUEST, IERROR)
                  call MPI_SEND(T, 1, doubletype, DEST, 103, MPI_COMM_WORLD, REQUEST, IERROR)
                  call MPI_SEND(ipr(i5), 1, doubletype, DEST, 104, MPI_COMM_WORLD, REQUEST, IERROR)
                  call MPI_SEND(counter, 1, MPI_INTEGER8, DEST, 105, MPI_COMM_WORLD, REQUEST, IERROR)
                  
                  if(DEST == (size_Of_Cluster-1) .or. counter == total_count) then 
                   ! Recieve 
                    write(100,*) "QQQQQ", DEST, counter, i1,i2,i3,i4,i5
                    do source = 1, DEST                                      
                      call MPI_RECV(nucl, 1, doubletype, source, 110, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
                      call MPI_RECV(localcounter, 1, MPI_INTEGER8, source, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)                      
                      write(1104,'(i15,E15.7)') localcounter,Jnucl_l ! counter,real(Jnucl)
                      rec_counter = rec_counter + 1
                      write(1001,rec=rec_counter) real(Jnucl_l)
                      flush(1104)
                      flush(1001)                      
                    end do 
                    DEST = 0
                  end if
                end if    
            end do               
         end do            
      end do 
   end do 
end do 
  flush(1104)
  flush(1001)
  close(1104)
  close(1001)  
end if 

 
counter_l=0
if (process_Rank > 0) then 
   do while ( counter_l + size_Of_Cluster - 1 <= total_count)
   ! SLAVE OR Master 
   call MPI_RECV(c_acid_l, 1, doubletype, 0, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
   !write(*,*) 'recieved ',c_acid_l
   call MPI_RECV(c_base_l, 1, doubletype, 0, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
   call MPI_RECV(CS_H2SO4_l, 1, doubletype, 0, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
   call MPI_RECV(T_l, 1, doubletype, 0, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
   call MPI_RECV(ipr_l, 1, doubletype, 0, 104, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
   call MPI_RECV(counter_l, 1, MPI_INTEGER8, 0, 105, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)
   !if (process_Rank .eq. 1) write(*,'(5E15.7)') c_acid_l,c_base_l,CS_H2SO4_l,T_l,ipr_l
   CALL  get_acdc_J(c_acid_l,c_base_l,c_org,CS_H2SO4_l,T_l,ipr_l,dt,solve_ss,Jnucl_l,diameter_acdc,Nuc_by_charge)
   call MPI_SEND(Jnucl_l, 1, doubletype,0, 110, MPI_COMM_WORLD, REQUEST, IERROR)
   call MPI_SEND(counter_l, 1, MPI_INTEGER8,0, 111, MPI_COMM_WORLD, REQUEST, IERROR)
   end do 
end if  
  deallocate(ipr)
  deallocate(Jnucl) 
  write(*,*) "END    counter" , counter_l , "in process " , process_Rank
  call MPI_Barrier(MPI_COMM_WORLD,IERROR)
  call MPI_FINALIZE(ierror)
  stop 
  
  
end subroutine create_lookup_mpi