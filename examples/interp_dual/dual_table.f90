program single_table
    use mo_kind, only : rk 
    use mo_serial_lookup, only : type_serial_lookup, serial_lookup_lookup => lookup , serial_lookup_init => init
    !
    implicit none 
    !
    real(kind=rk), allocatable,dimension(:) :: lookfor,interpVals,lookfor_2nd,interpVals_2nd
    type(type_serial_lookup) :: serial_lookup,serial_lookup_2nd
	integer :: i
	logical                 ::  interp_table_log
	logical                 ::  interp_indep_log
	
	
	!True: interplate log10 of the table values (adviced for ACDC). 
	!False: interpolate the actual values 
    interp_table_log  =.TRUE.
    !
	! TRUE: consider the values 
	interp_indep_log = .TRUE. 
	
	!
	call serial_lookup_init("../tables/lookup.desc","../tables/lookup.bin",serial_lookup)
	call serial_lookup_init("../tables/lookup.2.desc","../tables/lookup.2.bin",serial_lookup_2nd)
	!
    allocate(lookfor(serial_lookup%table_lookup%dimsCount))
    allocate(interpVals(serial_lookup%table_lookup%varsCount))

    allocate(lookfor_2nd(serial_lookup_2nd%table_lookup%dimsCount))
    allocate(interpVals_2nd(serial_lookup_2nd%table_lookup%varsCount))


	!
	if(interp_indep_log) then 
		serial_lookup%table_lookup%isLog10 = .FALSE. 
		serial_lookup_2nd%table_lookup%isLog10 = .FALSE. 
    end if 
	
	
    !
	if (interp_table_log) then 
	    serial_lookup%table_lookup%tbl = dlog10(serial_lookup%table_lookup%tbl)
		serial_lookup%table_lookup%isLog10 = .FALSE.
	    serial_lookup_2nd%table_lookup%tbl = dlog10(serial_lookup_2nd%table_lookup%tbl)
		serial_lookup_2nd%table_lookup%isLog10 = .FALSE.		
	end if 
	

	
    !
    lookfor = (/13.D0,15.D0,298.15D0/)  
	
    lookfor_2nd = (/13.D0,14.D0,298.15D0/)  
  
    call serial_lookup_lookup(serial_lookup,lookfor,interpVals)
    call serial_lookup_lookup(serial_lookup_2nd,lookfor_2nd,interpVals_2nd)
                
 
    
	!check weather all values are in range serial_lookup
    write(*,*) "check whether independent variables are in-range:" 
    write(*,*) "0: In-range" 
    write(*,*) "-1: Less than the minimum value." 
    write(*,*) "+1: Greater than the maximum value." 
	write(*,*) "First table:" 
	do i=1,serial_lookup%table_lookup%dimsCount
	    write(*,*) trim(serial_lookup%table_lookup%indepVarNames(i))//": ", serial_lookup%inrange(i)
		if(serial_lookup%inrange(i) == -1 .and. serial_lookup%table_lookup%isVapour(i)) then
		   write(*,*) trim(serial_lookup%table_lookup%indepVarNames(i))//" is a vapour and lower than the lowest limit. Setting the generation rate to zero fro the firest table."
		   interpVals = 0;
		end if 
	end do 
	!
    write(*,*) "Second table:" 
	do i=1,serial_lookup_2nd%table_lookup%dimsCount
	    write(*,*) trim(serial_lookup_2nd%table_lookup%indepVarNames(i))//": ", serial_lookup_2nd%inrange(i)
		if(serial_lookup_2nd%inrange(i) == -1 .and. serial_lookup_2nd%table_lookup%isVapour(i)) then
		   write(*,*) trim(serial_lookup_2nd%table_lookup%indepVarNames(i))//" is a vapour and lower than the lowest limit. Setting the generation rate to zero for the second table."
		   interpVals_2nd = 0;
		end if 
	end do 
	
	
    !
    write(*,*) 
    write(*,*) "interpolated values are (Asuming the output is of the same dimensions):" 
	do i=1,serial_lookup%table_lookup%varsCount
	    if (interp_table_log) then 
	        write(*,*) trim(serial_lookup%table_lookup%depVarNames(i))//": ", 10.**interpVals(i) + 1.**interpVals_2nd(i)
		else 
		    write(*,*) trim(serial_lookup%table_lookup%depVarNames(i))//": ", interpVals(i) + interpVals_2nd(i)
	    end if 		
	end do

	
    deallocate(interpVals)
    deallocate(lookfor)

    deallocate(interpVals_2nd)
    deallocate(lookfor_2nd)
	
end program single_table
    
    