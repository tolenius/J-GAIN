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
    
    
    ! In this example, 2 separate tables are used to determine the formation rate corresponding to each table, and
    ! the rates are summed to yield the total rate. In an actual implementation, these would correspond to different
    ! chemical pathways (but here only dummy tables are used for demonstration purposes).
    
    ! Also, the interpolation accuracy is optimised by interpolating along logarithmic axes for both the formation rate, and
    ! the independent parameters that are defined as 'logarithmic' upon table generation (as formation rates typically
    ! behave in this way on such log scales).
    
    
    ! Set log interpolation axis for the formation rate
    ! TRUE: interpolate log10 of the table values (recommended)
    ! FALSE: interpolate the actual values on linear axis
    interp_table_log  =.TRUE.
    
    ! Set log interpolation axis for the 'logarithmic' independent parameters
    ! TRUE: interpolate log10 of the parameter values (recommended); note that here the input values must be given as log10 of the actual values
    ! FALSE: interpolate the actual values on linear axis for all independent parameters
    interp_indep_log = .TRUE.
    
    
    call serial_lookup_init("../tables/lookup.desc","../tables/lookup.bin",serial_lookup)
    call serial_lookup_init("../tables/lookup.2.desc","../tables/lookup.2.bin",serial_lookup_2nd)
    
    allocate(lookfor(serial_lookup%table_lookup%dimsCount))
    allocate(interpVals(serial_lookup%table_lookup%varsCount))

    allocate(lookfor_2nd(serial_lookup_2nd%table_lookup%dimsCount))
    allocate(interpVals_2nd(serial_lookup_2nd%table_lookup%varsCount))


    ! Settings needed when using log axes
    if(interp_indep_log) then 
        serial_lookup%table_lookup%isLog10 = .FALSE. 
        serial_lookup_2nd%table_lookup%isLog10 = .FALSE. 
    end if 
    
    if(interp_table_log) then 
        where (serial_lookup%table_lookup%tbl.eq.0.0)
            serial_lookup%table_lookup%tbl = 1.d-20
        end where
        serial_lookup%table_lookup%tbl = dlog10(serial_lookup%table_lookup%tbl)
        !
        where (serial_lookup_2nd%table_lookup%tbl.eq.0.0)
            serial_lookup_2nd%table_lookup%tbl = 1.d-20
        end where        
        serial_lookup_2nd%table_lookup%tbl = dlog10(serial_lookup_2nd%table_lookup%tbl)
    end if

    
    ! Set input values for which the formation rate is determined; note that for 'log' parameters,
    ! log10 of the actual value needs to be used when using logarithmic axes
    if(interp_indep_log) then
        lookfor = (/13.D0,15.D0,298.15D0/)
        lookfor_2nd = (/13.D0,14.D0,298.15D0/)
    else
        lookfor = (/1.D13,1.D15,298.15D0/)
        lookfor_2nd = (/1.D13,1.D14,298.15D0/)
    end if
  
    call serial_lookup_lookup(serial_lookup,lookfor,interpVals)
    call serial_lookup_lookup(serial_lookup_2nd,lookfor_2nd,interpVals_2nd)
    
    if(interp_table_log) then
        interpVals=10.**interpVals
        interpVals_2nd=10.**interpVals_2nd
    end if
 
    
    ! Check whether all values are in range of serial_lookup
    write(*,*) "check whether independent variables are in-range:" 
    write(*,*) "0: In-range" 
    write(*,*) "-1: Less than the minimum value." 
    write(*,*) "+1: Greater than the maximum value." 
    write(*,*) "First table:" 
    do i=1,serial_lookup%table_lookup%dimsCount
        write(*,*) trim(serial_lookup%table_lookup%indepVarNames(i))//": ", serial_lookup%inrange(i)
        if(serial_lookup%inrange(i) == -1 .and. serial_lookup%table_lookup%isVapour(i)) then
           write(*,*) trim(serial_lookup%table_lookup%indepVarNames(i))//" is a vapour and lower than the lowest limit. Setting the formation rate to zero for the first table."
           interpVals = 0;
        end if 
    end do 
    !
    write(*,*) "Second table:" 
    do i=1,serial_lookup_2nd%table_lookup%dimsCount
        write(*,*) trim(serial_lookup_2nd%table_lookup%indepVarNames(i))//": ", serial_lookup_2nd%inrange(i)
        if(serial_lookup_2nd%inrange(i) == -1 .and. serial_lookup_2nd%table_lookup%isVapour(i)) then
           write(*,*) trim(serial_lookup_2nd%table_lookup%indepVarNames(i))//" is a vapour and lower than the lowest limit. Setting the formation rate to zero for the second table."
           interpVals_2nd = 0;
        end if 
    end do 

    
    ! Sum the rates (assuming that the number of output parameters is the same, i.e. here 1 as the formation rate is the only output)
    write(*,*) 
    write(*,*) "interpolated values are:" 
    do i=1,serial_lookup%table_lookup%varsCount
        write(*,*) trim(serial_lookup%table_lookup%depVarNames(i))//": ", interpVals(i) + interpVals_2nd(i)
    end do

    
    deallocate(interpVals)
    deallocate(lookfor)

    deallocate(interpVals_2nd)
    deallocate(lookfor_2nd)
    
end program single_table
    
    