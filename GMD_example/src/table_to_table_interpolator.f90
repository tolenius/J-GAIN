program compare_main
    use mo_kind, only : rk 
    use mo_serial_lookup, only : type_serial_lookup, serial_lookup_lookup => lookup , serial_lookup_init => init
    use mo_utils, only : getNewUnit,advanceLoopsBy	
    !
    implicit none 
    !
    real(kind=rk), allocatable,dimension(:) :: ref_lookfor,ref_interpVals,exp_lookfor,exp_interpVals
    type(type_serial_lookup) :: ref_serial_lookup,exp_serial_lookup
    integer :: u,i,k
	logical :: loop_over
    integer, allocatable, dimension(:) :: loopVals
    real(kind=rk), allocatable, dimension(:) :: minVals  ! independent 
    real(kind=rk), allocatable, dimension(:) :: steps    ! independent 
        
    character (len = 500 )  :: ref_dir
    character (len = 500 )  :: ref_desc
    character (len = 500 )  :: ref_bin
    character (len = 500 )  :: exp_dir
    character (len = 500 )  :: exp_desc
    character (len = 500 )  :: exp_bin
    character (len = 500 )  :: result_file_path
	logical                 :: write_ref
	logical                 :: write_exp
	logical                 :: benchmark, linear_log_intep
	logical, allocatable    :: islog10Temp(:)
	integer                 :: jump_over
    real :: start, finish

    integer :: jump_counter
	integer :: c1,c2,cr,cm 
	real ::  rate
	
    namelist /interp/  ref_dir,ref_desc,ref_bin,exp_dir,exp_desc,exp_bin,result_file_path, write_ref, write_exp, jump_over, benchmark, linear_log_intep
	                    
    write_ref = .true.
	write_exp  =.true.
	benchmark  =.false.
	linear_log_intep  =.false.
	jump_over = 10
	
	call getNewUnit(u)
    
	open(unit=u, file = "namelist.t2t", status='old')

    rewind(u)
    read(u, nml=interp) 
    close(u)
    
 
    if (benchmark) then 
	
        write(*,*) 'comparing: '  
        write(*,*)  trim(ref_dir)//"/"//trim(ref_desc)
        write(*,*) 'with: ' 
        write(*,*) trim(exp_dir)//"/"//trim(exp_desc)
	    
	    
	    
        call serial_lookup_init(trim(ref_dir)//"/"//trim(ref_desc), &
        & trim(ref_dir)//"/"//trim(ref_bin),ref_serial_lookup)
		
 
	    
        allocate(ref_lookfor(ref_serial_lookup%table_lookup%dimsCount))
        allocate(ref_interpVals(ref_serial_lookup%table_lookup%varsCount))
        
        
        call serial_lookup_init(trim(exp_dir)//"/"//trim(exp_desc), &
        & trim(exp_dir)//"/"//trim(exp_bin),exp_serial_lookup)
        
		
		if (linear_log_intep) then 
		    ref_serial_lookup%table_lookup%tbl = dlog10(ref_serial_lookup%table_lookup%tbl)
		    ref_serial_lookup%table_lookup%isLog10 = .False.
		    exp_serial_lookup%table_lookup%tbl = dlog10(exp_serial_lookup%table_lookup%tbl)
			exp_serial_lookup%table_lookup%isLog10 = .False.
		end if 
		
        allocate(exp_lookfor(exp_serial_lookup%table_lookup%dimsCount))
        allocate(exp_interpVals(exp_serial_lookup%table_lookup%varsCount))
	    
        
        
        allocate(loopVals(ref_serial_lookup%table_lookup%dimsCount))
        loopVals = -1 
        
        allocate(minVals(ref_serial_lookup%table_lookup%dimsCount))
        allocate(steps(ref_serial_lookup%table_lookup%dimsCount))
        minVals = ref_serial_lookup%table_lookup%minVals
        steps   = ref_serial_lookup%table_lookup%steps
	    
        if (write_exp .or. write_ref) then 
	        call getNewUnit(u)
            open(unit=u, file=trim(result_file_path), status='unknown')
	    end if 
	
	    loop_over = .true.
	    k = 0 
		write(*,*) "ignoring jump_over for benchmarking"
		call system_clock(count_rate=cr)
        call system_clock(count_max=cm)
		rate = REAL(cr)
		CALL SYSTEM_CLOCK(c1)
		call cpu_time(start)
        do while(loop_over)
	        k = k + 1 
            call advanceLoopsBy(ref_serial_lookup%table_lookup%idims,loopVals)
             
            do i=1,ref_serial_lookup%table_lookup%dimsCount
               ref_lookfor(i) = minVals(i) + (loopVals(i)-1)*steps(i)
               if (ref_serial_lookup%table_lookup%isLog10(i)) ref_lookfor(i) = 10**ref_lookfor(i)
            end do 
                           
            call serial_lookup_lookup(exp_serial_lookup,ref_lookfor,exp_interpVals)
	    	    
            if (linear_log_intep) then
			    exp_interpVals = 10.D0**exp_interpVals
            end if

			
	    	jump_counter = jump_counter + 1
            if (all(loopVals >= ref_serial_lookup%table_lookup%idims)) loop_over = .false. 
          
        end do 
              
        call cpu_time(finish)
		CALL SYSTEM_CLOCK(c2)
        print '("XXXXXXRRR  ", f8.5,f8.5,i20,i20)', finish-start, (c2 - c1)/rate , product(exp_serial_lookup%table_lookup%idims), product(ref_serial_lookup%table_lookup%idims)
		
	else 
	
        write(*,*) 'comparing: '  
        write(*,*)  trim(ref_dir)//"/"//trim(ref_desc)
        write(*,*) 'with: ' 
        write(*,*) trim(exp_dir)//"/"//trim(exp_desc)
	    
	    
	    
        call serial_lookup_init(trim(ref_dir)//"/"//trim(ref_desc), &
        & trim(ref_dir)//"/"//trim(ref_bin),ref_serial_lookup)
	    
        allocate(ref_lookfor(ref_serial_lookup%table_lookup%dimsCount))
        allocate(islog10Temp(ref_serial_lookup%table_lookup%dimsCount))
		
       islog10Temp(:) = ref_serial_lookup%table_lookup%isLog10(:)



        allocate(ref_interpVals(ref_serial_lookup%table_lookup%varsCount))


		
        call serial_lookup_init(trim(exp_dir)//"/"//trim(exp_desc), &
        & trim(exp_dir)//"/"//trim(exp_bin),exp_serial_lookup)
        
		
		if (linear_log_intep) then 
		    ref_serial_lookup%table_lookup%tbl = dlog10(ref_serial_lookup%table_lookup%tbl)
			ref_serial_lookup%table_lookup%isLog10 = .False.
		    exp_serial_lookup%table_lookup%tbl = dlog10(exp_serial_lookup%table_lookup%tbl)
			exp_serial_lookup%table_lookup%isLog10 = .False.
		end if 

		
        allocate(exp_lookfor(exp_serial_lookup%table_lookup%dimsCount))
        allocate(exp_interpVals(exp_serial_lookup%table_lookup%varsCount))
	    
        
        
        allocate(loopVals(ref_serial_lookup%table_lookup%dimsCount))
        loopVals = -1 
        
        allocate(minVals(ref_serial_lookup%table_lookup%dimsCount))
        allocate(steps(ref_serial_lookup%table_lookup%dimsCount))
        minVals = ref_serial_lookup%table_lookup%minVals
        steps   = ref_serial_lookup%table_lookup%steps
	    
        if (write_exp .or. write_ref) then 
	        call getNewUnit(u)
            open(unit=u, file=trim(result_file_path), status='unknown')
	    end if 
	    
        call cpu_time(start)
	    	
	    loop_over = .true.
	    k = 0 
	    jump_counter = jump_over
        do while(loop_over)
	        k = k + 1 
            call advanceLoopsBy(ref_serial_lookup%table_lookup%idims,loopVals)
            if (jump_counter == jump_over) then 
                do i=1,ref_serial_lookup%table_lookup%dimsCount
                   ref_lookfor(i) = minVals(i) + (loopVals(i)-1)*steps(i)
                   if (ref_serial_lookup%table_lookup%isLog10(i)) ref_lookfor(i) = 10**ref_lookfor(i)
                end do 
                
                exp_lookfor = ref_lookfor
                
                call serial_lookup_lookup(ref_serial_lookup,ref_lookfor,ref_interpVals)
                call serial_lookup_lookup(exp_serial_lookup,exp_lookfor,exp_interpVals)
				! Recalc lookfor values 
	            do i=1,ref_serial_lookup%table_lookup%dimsCount
                   ref_lookfor(i) = minVals(i) + (loopVals(i)-1)*steps(i)
                   if (islog10Temp(i)) ref_lookfor(i) = 10**ref_lookfor(i)
                end do 
				exp_lookfor = ref_lookfor
				
				!
                if (linear_log_intep) then
			       exp_interpVals = 10.D0**exp_interpVals
			       ref_interpVals = 10.D0**ref_interpVals
                end if
			
	    	    if (write_ref .and. write_exp ) then 
                   write(u,*) k, (ref_lookfor(i), i=1,size(ref_lookfor,1)), &
	    	         (ref_interpVals(i), i=1,size(ref_interpVals,1)), &
	    	         (exp_interpVals(i), i=1,size(exp_interpVals,1))
	    	    else if (write_ref .and. .not. write_exp ) then 
                    write(u,*) k,(ref_lookfor(i), i=1,size(ref_lookfor,1)), &
	    	    	(ref_interpVals(i), i=1,size(ref_interpVals,1)), &
	    	        (ref_interpVals(i), i=1,size(ref_interpVals,1)) 
	    	    else if (.not. write_ref .and. write_exp ) then 
                    write(u,*) k, (ref_lookfor(i), i=1,size(ref_lookfor,1)), &
	    	    	(exp_interpVals(i), i=1,size(exp_interpVals,1)), &
	    	        (exp_interpVals(i), i=1,size(exp_interpVals,1))
                end if 		
	    	    jump_counter = 0 
	    	end if 
	    	jump_counter = jump_counter + 1
            if (all(loopVals >= ref_serial_lookup%table_lookup%idims)) loop_over = .false. 
          
        end do 
              ! put code to test here
        call cpu_time(finish)
        print '("Time = ",f6.3," seconds.")', finish-start
	    if (write_exp .or. write_ref) then
	        flush(u)
            close(u)
        end if 
	end if 
end program compare_main
    
    