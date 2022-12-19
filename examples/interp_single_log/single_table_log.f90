program single_table_log
    use mo_kind, only : rk 
    use mo_serial_lookup, only : type_serial_lookup, serial_lookup_lookup => lookup , serial_lookup_init => init
    !
    implicit none 
    !
    real(kind=rk), allocatable,dimension(:) :: lookfor,interpVals
    type(type_serial_lookup) :: serial_lookup
    logical                 ::  interp_table_log
    logical                 ::  interp_indep_log
    integer :: i 
    !True: interplate log10 of the table values (adviced for ACDC). 
    !False: interpolate the actual values 
    interp_table_log  =.TRUE.
    !
    ! TRUE: consider the values 
    interp_indep_log = .TRUE. 
    
    !
    call serial_lookup_init("../tables/lookup.desc","../tables/lookup.bin",serial_lookup)
    !
    allocate(lookfor(serial_lookup%table_lookup%dimsCount))
    allocate(interpVals(serial_lookup%table_lookup%varsCount))

    !
    if(interp_indep_log) then 
        serial_lookup%table_lookup%isLog10 = .FALSE. 
    end if 
    
    !
    if (interp_table_log) then 
        where (serial_lookup%table_lookup%tbl.eq.0.0)
            serial_lookup%table_lookup%tbl = 1.d-20
        end where        
        serial_lookup%table_lookup%tbl = dlog10(serial_lookup%table_lookup%tbl)
    end if 
    
    ! 
    if(interp_indep_log)
        lookfor = (/11.D0,14D0,298.15D0/) ! Compare to  (/1.D11,1.D14,298.15D0/) in single_table.f90
    else 
        lookfor = (/1.D11,1.D14,298.15D0/)
    end if 
    
    call serial_lookup_lookup(serial_lookup,lookfor,interpVals)
     

    write(*,*) 
    write(*,*) "interpolated values are:" 
    do i=1,serial_lookup%table_lookup%varsCount
        if (interp_table_log) then 
            write(*,*) trim(serial_lookup%table_lookup%depVarNames(i))//": ", 10.**interpVals(i)
        else 
            write(*,*) trim(serial_lookup%table_lookup%depVarNames(i))//": ", interpVals(i)
        end if 
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

end program single_table_log
    
    