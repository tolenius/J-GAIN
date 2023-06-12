module mo_serial_lookup
    use mo_kind
    use mo_table_lookup
    use mo_debug
    
    implicit none
    
    public init
    public lookup
    public type_serial_lookup
 

    private
    type type_serial_lookup
        real(kind=rk),    allocatable,dimension(:)   :: lowWghts
        real(kind=rk),    allocatable,dimension(:)   :: interpWghts
        integer(kind=ik), allocatable,dimension(:)   :: lowIdxs
        integer(kind=ik), allocatable,dimension(:)   :: upIdxs
        integer(kind=ik), allocatable,dimension(:)   :: tblIdxs
        integer(kind=ik), allocatable,dimension(:,:) :: comIdxs
        integer, allocatable,dimension(:)   :: inrange
        
        type(type_table_lookup) :: table_lookup
    end type type_serial_lookup
    
    interface init
      module procedure type_serial_lookup_init_from_bin_file
    end interface init
   
    interface finalize 
      module procedure type_serial_lookup_finalize
    end interface finalize
    
    
    interface lookup
      module procedure type_serial_lookup_lookupImpl
    end interface lookup
    
    
contains    

!> @brief initilize a serial lookup structure from a desc and a binary file
!!
!! @param[in]   descriptor_file_path     path to the descrriptor file 
!! @param[in]   table_file_path          path to the binary file 
!! @param[out]  serial_lookup            the lookup structure 
subroutine type_serial_lookup_init_from_bin_file(descriptor_file_path,table_file_path,serial_lookup)
    use mo_table_lookup, only : load
    implicit none
    type(type_serial_lookup),intent(inout) :: serial_lookup
    !    
    character(len=*), intent(in) :: descriptor_file_path
    character(len=*), intent(in) :: table_file_path
    !
    call type_serial_lookup_finalize(serial_lookup)
    call load(descriptor_file_path,table_file_path,serial_lookup%table_lookup)
	 
end subroutine type_serial_lookup_init_from_bin_file

!> @brief do the lookup job by calling the necessary subroutines.
!!
!! @param[in]   serial_lookup    the lookup structure. 
!! @param[in]   lookfor          values of the independet variables.
!! @param[out]  interpVals       the interpolated output values.            
subroutine type_serial_lookup_lookupImpl(serial_lookup,lookfor,interpVals)
    use mo_table_lookup 
    implicit none
    !
    type(type_serial_lookup),intent(inout) :: serial_lookup    
    !
    real(kind=rk),allocatable,dimension(:), intent(in) :: lookfor
    !
    real(kind=rk),allocatable,dimension(:), intent(out) :: interpVals
    !
    call type_table_lookup_find_surrounding_indices(                     &
	    serial_lookup%table_lookup,lookfor,serial_lookup%lowIdxs,         &
        serial_lookup%upIdxs,serial_lookup%lowWghts,serial_lookup%inrange)  

    !if (.not. all(serial_lookup%inrange))  return 
    
    call type_table_lookup_find_interpolation_weights(                  &
	    serial_lookup%lowIdxs,serial_lookup%upIdxs,serial_lookup%lowWghts,&
        serial_lookup%comIdxs,serial_lookup%interpWghts)
       
    call type_table_lookup_find_table_indices(                           &
	    serial_lookup%table_lookup,serial_lookup%comIdxs,serial_lookup%tblIdxs)
    
    call interpolate(serial_lookup%table_lookup,serial_lookup%tblIdxs,serial_lookup%   &
       interpWghts,interpVals)
      
end subroutine type_serial_lookup_lookupImpl

!> @brief releases the allocated arrays.
!!
subroutine type_serial_lookup_finalize(serial_lookup)
    use mo_table_lookup, only : tble_finalize => finalize 
    implicit none 
    type(type_serial_lookup),intent(inout) :: serial_lookup
    !
    if (allocated(serial_lookup%lowWghts   )) deallocate(serial_lookup%lowWghts   )
    if (allocated(serial_lookup%interpWghts)) deallocate(serial_lookup%interpWghts)
    if (allocated(serial_lookup%lowIdxs    )) deallocate(serial_lookup%lowIdxs    )
    if (allocated(serial_lookup%upIdxs     )) deallocate(serial_lookup%upIdxs     )
    if (allocated(serial_lookup%tblIdxs    )) deallocate(serial_lookup%tblIdxs    )
    if (allocated(serial_lookup%comIdxs    )) deallocate(serial_lookup%comIdxs    )
    if (allocated(serial_lookup%inrange    )) deallocate(serial_lookup%inrange    )
    call tble_finalize(serial_lookup%table_lookup)
end subroutine type_serial_lookup_finalize
end module mo_serial_lookup