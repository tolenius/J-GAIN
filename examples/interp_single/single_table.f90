program single_table
    use mo_kind, only : rk 
    use mo_serial_lookup, only : type_serial_lookup, serial_lookup_lookup => lookup , serial_lookup_init => init
    !
    implicit none 
    !
    real(kind=rk), allocatable,dimension(:) :: lookfor,interpVals
    type(type_serial_lookup) :: serial_lookup
	integer :: i
	
	!
	call serial_lookup_init("../tables/lookup.desc","../tables/lookup.bin",serial_lookup)
	!
    allocate(lookfor(serial_lookup%table_lookup%dimsCount))
    allocate(interpVals(serial_lookup%table_lookup%varsCount))

	
    ! fill the actual values 
    lookfor = (/1.D11,1.D14,298.15D0/)  
  
    call serial_lookup_lookup(serial_lookup,lookfor,interpVals)
                
    !
    write(*,*) 
    write(*,*) "interpolated values are:" 
	do i=1,serial_lookup%table_lookup%varsCount
	    write(*,*) trim(serial_lookup%table_lookup%depVarNames(i))//": ", interpVals(i)
	end do 
    
	!check weather all values are in range serial_lookup
    write(*,*) "Are Independent variables in-range?" 
    write(*,*) "0: In-range" 
    write(*,*) "-1: Less than the minimum value." 
    write(*,*) "+1: Greater than the maximum value." 
	
	do i=1,serial_lookup%table_lookup%dimsCount
	    write(*,*) trim(serial_lookup%table_lookup%indepVarNames(i))//": ", serial_lookup%inrange(i)
	end do 
    
    deallocate(interpVals)
    deallocate(lookfor)

end program single_table
    
    