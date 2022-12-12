module mo_gen  
  use mo_lookup_config, only: gen_config_ranges_type
  use mo_kind, only : rk4 
#ifdef MPI
  use mo_gen_parallel
#else
  use mo_gen_serial
#endif
  implicit none 
  
  
  interface genLookupTable
    module procedure genLookupTableImp
  end interface genLookupTable 
  
contains
 
 
subroutine genLookupTableImp(gen_config_ranges,table)
implicit none 
  ! Arguments 
  type(gen_config_ranges_type), intent(in) :: gen_config_ranges
  real(kind=rk4), allocatable,  dimension(:,:), intent(out) :: table 
  
   call genLookupTable_(gen_config_ranges,table)
end subroutine genLookupTableImp

end module mo_gen 
