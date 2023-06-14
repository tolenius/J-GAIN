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
    
    
    ! This example can be used as a subroutine template for implementing formation rate look-up tables in a host model
    
    ! In this example, 2 separate tables are used to determine the formation rate corresponding to each table, and
    ! the rates are summed to yield the total rate. In an actual implementation, these would correspond to different
    ! chemical pathways (but here only dummy tables are used for demonstration purposes).
    
    ! Also, the interpolation accuracy is optimised by interpolating along logarithmic axes for both the formation rate, and
    ! the independent parameters that are defined as 'logarithmic' upon table generation (as formation rates typically
    ! behave smoothly on such log scales).
    
    !------------------------------------------------------------------------------------------
    ! User-defined settings for interpolation along linear or logarithmic axes
    !------------------------------------------------------------------------------------------
    
    ! Set logarithmic interpolation axis for the formation rate
    ! TRUE: interpolate log10 of the table values (recommended)
    ! FALSE: interpolate the actual values on linear axis
    interp_table_log  =.TRUE.
    
    ! Set logarithmic interpolation axis for the 'logarithmic' independent parameters
    ! TRUE: interpolate log10 of the parameter values (recommended); note that here the input values must be given as log10 of the actual values
    ! FALSE: interpolate the actual values on linear axis for all independent parameters
    interp_indep_log = .TRUE.
    
    !------------------------------------------------------------------------------------------
    ! Initialization at the start of a model run
    ! (This is not modified by the user)
    !------------------------------------------------------------------------------------------

    ! Initialize the interpolation procedures (for each table): load the tables and parameter information
    call serial_lookup_init("../tables/lookup.desc","../tables/lookup.bin",serial_lookup)
    call serial_lookup_init("../tables/lookup.2.desc","../tables/lookup.2.bin",serial_lookup_2nd)
    
    ! Allocate the sizes of the table input/output, i.e. the numbers of independent parameters (and dependent parameters, which is here only one)
    allocate(lookfor(serial_lookup%table_lookup%dimsCount))
    allocate(interpVals(serial_lookup%table_lookup%varsCount))

    allocate(lookfor_2nd(serial_lookup_2nd%table_lookup%dimsCount))
    allocate(interpVals_2nd(serial_lookup_2nd%table_lookup%varsCount))

    ! Settings needed when using logarithmic axes
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

    !------------------------------------------------------------------------------------------
    ! Start the model run: call the interpolator with given input for each table; the user needs
    ! to set the input values of the independent parameters according to their model application
    !
    ! Here we perform one call for each table for demonstration, but in an actual application
    ! the interpolator is called repeatedly with input ambient conditions from the host model
    !------------------------------------------------------------------------------------------
    
    ! Set input values for the ambient conditions at which the formation rate is determined; note that for
    ! 'logarithmic' parameters, log10 of the actual value needs to be used when using logarithmic axes
    if(interp_indep_log) then
        lookfor = (/13.D0,15.D0,298.15D0/)
        lookfor_2nd = (/13.D0,14.D0,298.15D0/)
    else
        lookfor = (/1.D13,1.D15,298.15D0/)
        lookfor_2nd = (/1.D13,1.D14,298.15D0/)
    end if
    
    ! Call the interpolator to obtain the formation rate
    call serial_lookup_lookup(serial_lookup,lookfor,interpVals)
    call serial_lookup_lookup(serial_lookup_2nd,lookfor_2nd,interpVals_2nd)
    
    ! If the formation rate is interpolated on a logarithmic axis, convert the log10 of the rate to the actual rate
    if(interp_table_log) then
        interpVals=10.**interpVals
        interpVals_2nd=10.**interpVals_2nd
    end if
    
    ! Check whether all input values are within the ranges of the tables
    ! If a vapor concentration is below the lowest value of the table (in practice too low for particle formation),
    ! it is recommended to set the formation rate to zero to avoid artefact rates
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
    
    ! Sum the rates obtained from the different tables to obtain the total formation rate from all chemical pathways
    ! (assuming that the number of output parameters is the same, i.e. here one as the formation rate is the only output)
    write(*,*) 
    write(*,*) "interpolated values are:" 
    do i=1,serial_lookup%table_lookup%varsCount
        write(*,*) trim(serial_lookup%table_lookup%depVarNames(i))//": ", interpVals(i) + interpVals_2nd(i)
    end do
    
    !------------------------------------------------------------------------------------------
    ! End of the model run; no more calls to the interpolator
    !------------------------------------------------------------------------------------------
    
    ! Deallocate, if needed - included here only for demonstration, in an actual model implementation
    ! the input/output sizes should not be deallocated between calls to the table interpolation routine
    deallocate(interpVals)
    deallocate(lookfor)

    deallocate(interpVals_2nd)
    deallocate(lookfor_2nd)
    
end program single_table
