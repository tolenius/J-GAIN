program lookup_gen
  use mo_lookup_config
  use mo_gen
  use mo_writer
  use mo_kind, only: rk4
  use mo_config, only: canWrite, isMaster
  use mo_gen_profile
#ifdef MPI
  use mo_mpi, only: init_mpi,finalize_mpi
#endif
  implicit none
  
#ifdef MPI
  include 'mpif.h'
#endif

  type(gen_config_ranges_type):: gen_config_ranges
  type(gen_config_profile_type):: gen_config_profile
  real(kind=rk4), allocatable,  dimension(:,:) :: table 
  
#ifdef MPI
  call init_mpi() 
#endif

   call initConf("namelist.gen", gen_config_ranges,gen_config_profile) 
#ifdef MPI
   call mpi_barrier(MPI_COMM_WORLD) 
#endif
   if (gen_config_profile%isActiveMode) then 
#ifdef MPI 
        write(*,*) 'Profile mode is not implemented to run in parallel. Only one process (Master) will be used.'  
        if(isMaster) call genAndWriteProfile(gen_config_profile)
#else 
        call genAndWriteProfile(gen_config_profile)  
#endif	   
   else if (gen_config_ranges%isActiveMode) then 
       call genLookupTable(gen_config_ranges,table) 

#ifdef MPI
       call mpi_barrier(MPI_COMM_WORLD) 
#endif

       if(canWrite) then 
         call writeTable(gen_config_ranges,table)
       end if 
   
       write(*,*) '===== Normal End ====='
	   
   end if 
   
   


#ifdef MPI
  call finalize_mpi() 
#endif

end program lookup_gen
 