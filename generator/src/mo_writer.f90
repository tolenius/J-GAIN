module mo_writer
  use mo_lookup_config, only: gen_config_ranges_type,FORMAT_BIN,FORMAT_NC3
  use mo_kind, only : rk4 
  use mo_runtime, only: abort 
  use mo_debug, only: debug 
  use mo_utils, only: getNewUnit

implicit none 
  
  public writeTable
 
  
  private 
  interface writeTable
    module procedure writeTableImpl
  end interface writeTable
  
  
contains 

subroutine writeTableImpl(gen_config_ranges,table)
  use mo_config, only : canWrite 
implicit none 
  ! Arguments 
  type(gen_config_ranges_type), intent(in) :: gen_config_ranges
  real(kind=rk4), allocatable,  dimension(:,:), intent(in) :: table 
  ! Body 
  
  if( .not. canWrite ) return 
  
  select case (gen_config_ranges%outputFormat )
     case (FORMAT_BIN)
        call writeBin(gen_config_ranges,table)
     case (FORMAT_NC3)
        call abort("NetCDF is not implemented yet")
     case default
        call abort("Not suported File Format")
  
  end select


  
end subroutine writeTableImpl
  
subroutine writeBin(gen_config_ranges,table)
  implicit none 
  ! Arguments 
  type(gen_config_ranges_type), intent(in) :: gen_config_ranges
  real(kind=rk4), allocatable,  dimension(:,:), intent(in) :: table 
  
  integer :: i
  integer :: du,bu  
  logical, allocatable , dimension(:)  :: dims_gt_1
  !
  integer :: rc ! error code 
  character(len=500) :: file_name ! temp file name
  
  
  dims_gt_1 = gen_config_ranges%dimensions > 1 
  ! first of all write the desciptor  file 
  
  if (debug) write(*,*) "Started: Create descriptor file"
  call getNewUnit(du)
  file_name = trim(gen_config_ranges%outputDirectory)//"/"//trim(gen_config_ranges%outputFileBase)//".desc"
  OPEN(unit = du ,file = trim(file_name), &
     & STATUS='unknown',ACTION='WRITE', iostat=rc) 
  if (rc /= 0) stop "Error: failed to open file: "//trim(file_name)

  write(du,*) "Dep Vars Count (depCount)" 
  write(du,*) size(table,2)
  write(du,*) "Var Names" 
  do i=1,gen_config_ranges%depCount
    write(du,*) trim(gen_config_ranges%outVarNames(i))
  end do 
  write(du,*) "Var Units" 
  do i=1,gen_config_ranges%depCount
    write(du,*) trim(gen_config_ranges%outUnits(i))
  end do 
  write(du,*) "Indep Vars Count (dimCount)"  
  write(du,*)  count(dims_gt_1) ! gen_config_ranges%dimCount 
 
  write(du,*) "Dims"
  do i=1,gen_config_ranges%dimCount
    if (dims_gt_1(i)) write(du,*) gen_config_ranges%dimensions(i)
  end do 
  
  write(du,*) "Indep Vars Names"  
  
  do i=1,gen_config_ranges%dimCount 
    if (dims_gt_1(i))  write(du,*) trim(gen_config_ranges%varNames(i))
  end do
  
  write(du,*) "Indep Vars Units"  
  do i=1,gen_config_ranges%dimCount
    if (dims_gt_1(i)) write(du,*) trim(gen_config_ranges%units(i))
  end do 
  
  write(du,*)"minVals,maxVals" 
  do i=1,gen_config_ranges%dimCount
    if (dims_gt_1(i)) write(du,'(2e20.8)') gen_config_ranges%minVals(i),gen_config_ranges%maxVals(i)
  end do 
   
  write(du,*) 'isLog10'
  
  do i=1,gen_config_ranges%dimCount
    if (dims_gt_1(i))  write(du,*) gen_config_ranges%isLog10(i) 
  end do 
   
  write(du,*) 'isVapour'
  
  do i=1,gen_config_ranges%dimCount
    if (dims_gt_1(i))  write(du,*) gen_config_ranges%isVapour(i) 
  end do 

  
  write(du,*) 'BinFile'
  write(du,*) trim(gen_config_ranges%outputDirectory)//"/"//trim(gen_config_ranges%outputFileBase)//".bin"
  
  write(du,*)  "totalCount"
  write(du,*) gen_config_ranges%totalCount
  flush(du)
  close(du)
  
  if (debug) write(*,*) "Finished: Create descriptor file"
  ! 
  
  ! Writing the binary file 
  if (debug) write(*,*) "Started writing the binary file"
  file_name = trim(gen_config_ranges%outputDirectory)//"/"//trim(gen_config_ranges%outputFileBase)//".bin"
  open(unit = bu, file = trim(file_name),&
    &status="unknown",form="unformatted", access = "direct", recl = 4*size(table), iostat=rc)
  if (rc /= 0) stop "Error: failed to open file: "//trim(file_name)

	
  write(bu,rec=1) table
  close(bu)
  
  if (debug) write(*,*) "Started writing the text file"
  file_name = trim(gen_config_ranges%outputDirectory)//"/"//trim(gen_config_ranges%outputFileBase)//".txt"
  open(unit = bu, file = trim(file_name),&
    &status="unknown", iostat=rc)
  if (rc /= 0) stop "Error: failed to open file: "//trim(file_name) 
  do i=1,gen_config_ranges%totalCount
    write(bu,*) i, table(i,1)
   ! write(bu,*) i,real(table(i,1))
  end do 
 
  close(bu)
  
  
  if (debug) write(*,*) "Finshed writing to the file"
end subroutine writeBin
  
end module mo_writer