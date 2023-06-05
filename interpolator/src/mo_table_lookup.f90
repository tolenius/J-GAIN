module mo_table_lookup
    use mo_kind
    use mo_debug    
    implicit none
    
    private
    
    public :: load
    public :: type_table_lookup
    public :: interpolate
    public :: finalize

    public :: type_table_lookup_find_table_index
    public :: type_table_lookup_find_nearest_index
    public :: type_table_lookup_find_surrounding_indices
    public :: type_table_lookup_find_table_indices

    public :: type_table_lookup_find_interpolation_weights
	
	
! dv: dependent varaibels  
! module types  
    type type_table_lookup
        real(kind=rk), allocatable, dimension(:,:) :: tbl  ! first index is row number and the second index is the vraible index 
        real(kind=rk), allocatable, dimension(:) :: minVals, maxVals ! independent 
        real(kind=rk), allocatable, dimension(:) :: steps    ! independent 
        integer(kind=ik), allocatable, dimension(:) :: idims ! independent 
        logical, allocatable, dimension(:) :: isLog10
        logical, allocatable, dimension(:) :: isVapour
        integer(kind=ik) :: dimsCount ! independent 
        character(len=20), dimension(:), allocatable :: indepVarNames ! independent 
        character(len=20), dimension(:), allocatable :: depVarNames ! dependent 
        character(len=20), dimension(:), allocatable :: units	! independent   
        integer(kind=ik) :: varsCount ! dependent 
        integer(kind=ik) :: valuesCount ! mult idims      
	    real(kind=rk), allocatable, dimension(:) :: maxval_plus_epsilon ! independent 
	    real(kind=rk), allocatable, dimension(:) :: minval_minus_epsilon ! independent 
    end type type_table_lookup    
! module variables 

! module interfaces
    interface load
        module procedure load_from_bin_file
    end interface load
    



    
    interface interpolate
       module procedure type_table_lookup_interpolate    
    end interface interpolate

    interface finalize
       module procedure type_table_lookup_finialize    
    end interface finalize
    
contains

! Type table_lookup Operations
subroutine type_table_lookup_finialize(table_lookup)
    implicit none
    ! Arguments
    type(type_table_lookup), intent(inout) :: table_lookup
    !
    if (debug) write(*,*) '==== Entering type_table_lookup_finialize ===='
    !
    if(allocated(table_lookup%tbl))       deallocate(table_lookup%tbl)
    if(allocated(table_lookup%minVals))   deallocate(table_lookup%minVals)
    if(allocated(table_lookup%maxVals))   deallocate(table_lookup%maxVals)
    if(allocated(table_lookup%steps))     deallocate(table_lookup%steps)
    if(allocated(table_lookup%idims))     deallocate(table_lookup%idims)
    if(allocated(table_lookup%isLog10))   deallocate(table_lookup%isLog10) 
    if(allocated(table_lookup%isVapour))   deallocate(table_lookup%isVapour) 
    if(allocated(table_lookup%depVarNames))   deallocate(table_lookup%depVarNames) 
    if(allocated(table_lookup%indepVarNames))   deallocate(table_lookup%indepVarNames) 
	
	if(allocated(table_lookup%maxval_plus_epsilon ))  deallocate(table_lookup%maxval_plus_epsilon )
	if(allocated(table_lookup%minval_minus_epsilon))  deallocate(table_lookup%minval_minus_epsilon)
end subroutine type_table_lookup_finialize

      
subroutine type_table_lookup_find_table_index(table_lookup,indices,tblIdx)
    ! Find the index in the table type_table_lookup%tbl
    implicit none
    type(type_table_lookup), intent(in) :: table_lookup
    integer(kind=ik), dimension(:) , intent(in) :: indices ! independent 
    integer(kind=ik), intent(out) :: tblIdx ! dependent
    ! 
    integer :: i,ml
    !
    if (debug) write(*,*) '==== Entering type_table_lookup_Find_table_excact_value ===='
    !
    tblIdx = indices(table_lookup%dimsCount)
    ml  = table_lookup%idims(table_lookup%dimsCount) 

    do i = table_lookup%dimsCount-1,1,-1
        tblIdx = tblIdx + (indices(i)-1)*ml
        ml = ml*table_lookup%idims(i)    
    enddo
    
    !if (debug) write(*,*) "Caclated table Index ", tblIdx
    
    !depVals = table_lookup%tbl(tblIdx,:)
end subroutine type_table_lookup_find_table_index


subroutine type_table_lookup_find_table_indices(table_lookup,comIdxs,tblIds)
    ! Find multiple indices in the table type_table_lookup%tbl
    implicit none
    type(type_table_lookup), intent(in) :: table_lookup
    integer(kind=ik), dimension(:,:) , intent(in) :: comIdxs ! independent first index is table_lookup%dimsCount the second is 2**table_lookup%dimsCount
    integer(kind=ik), allocatable, dimension(:), intent(out) :: tblIds ! dependent
    ! 
    integer :: i,ml,k,n
    !
    if (debug) write(*,*) '==== Entering type_table_lookup_find_table_indices ===='   
    !
    n = size(comIdxs,2)
    !
    if (.not. allocated(tblIds)) allocate(tblIds(n))
    do k = 1,n
        tblIds(k) = comIdxs(table_lookup%dimsCount,k)
        ml  = table_lookup%idims(table_lookup%dimsCount) 
        do i = table_lookup%dimsCount-1,1,-1
            tblIds(k) = tblIds(k) + (comIdxs(i,k)-1)*ml
            ml = ml*table_lookup%idims(i)    
        enddo
    enddo
    
    if (debug) then 
        write(*,*) 'Found table indices   ',tblIds   
        write(*,*) '===='
    endif 
   
end subroutine type_table_lookup_find_table_indices



subroutine type_table_lookup_find_nearest_index(table_lookup,lookfor,indices,inrange)
    implicit none 
    type(type_table_lookup), intent(in) :: table_lookup   
    real(kind=rk), dimension(:) , intent(in) :: lookfor ! independent
    integer, dimension(:), allocatable, intent(out) :: indices ! independent
    integer, dimension(:), allocatable, intent(out) :: inrange ! independent
    !
    integer :: i 
    real(kind=rk) :: zlookfor   
    !
    if (debug) write(*,*) '==== Entering type_table_lookup_find_nearest_index ===='
    ! 
!#ifdef CHECKALLOC TODO
    if (.not. allocated(indices)) allocate(indices(table_lookup%dimsCount))
    if (.not. allocated(inrange)) allocate(inrange(table_lookup%dimsCount))
!#endif   
 
 
    do i=1,table_lookup%dimsCount
	
	    if (table_lookup%isLog10(i)) then 
          zlookfor = dlog10(lookfor(i))
        else
          zlookfor = lookfor(i)
        end if 
	    
	    if (table_lookup%idims(i) > 1 ) then 
	    
            if (table_lookup%minval_minus_epsilon(i) >  zlookfor ) then 
                inrange(i) = -1 
                indices(i) = 1
   	        elseif ((table_lookup%maxval_plus_epsilon(i) <  zlookfor )) then 
	            inrange(i) = +1
		  	    indices(i) = table_lookup%idims(i)
	        else 
	            inrange(i) = 0
		  	    indices(i) = NINT((zlookfor-table_lookup%minVals(i))/table_lookup%steps(i)) + 1
	        end if 
        else 
	        inrange(i) = 0
	        indices(i) = 1
	    end if 
	  
	  
    enddo
   
end subroutine type_table_lookup_find_nearest_index 

subroutine type_table_lookup_find_surrounding_indices(table_lookup,lookfor,lowIdxs,upIdxs,lowWghts,inrange)
    implicit none 
    type(type_table_lookup), intent(in) :: table_lookup      
    real(kind=rk), dimension(:), intent(in) :: lookfor ! independent is the actual value  looking for without log10 even if minVals, maxVals and steps  is log10
    integer(kind=ik), dimension(:), allocatable, intent(out) :: lowIdxs,upIdxs ! independent
    real(kind=rk), dimension(:), allocatable, intent(out) :: lowWghts ! independent ! in [0,1] normalized distance from  table_lookup%minVals(lowIdxs)
    integer, dimension(:), allocatable, intent(out) :: inrange ! independent
    !
    integer :: i 
    real(kind=rk) :: zlookfor   
    real(kind=rk) :: zlowloc,znextlloc  
    !
    if (debug) write(*,*) '==== Entering type_table_lookup_find_surrounding_indices ===='
    !   
!#ifdef CHECKALLOC TODO  
    if (.not. allocated(lowIdxs )) allocate(lowIdxs(table_lookup%dimsCount))
    if (.not. allocated(upIdxs  )) allocate(upIdxs(table_lookup%dimsCount))
    if (.not. allocated(lowWghts)) allocate(lowWghts(table_lookup%dimsCount))
    if (.not. allocated(inrange )) allocate(inrange(table_lookup%dimsCount))
!#endif  
    if (debug) then 
        write(*,*) '==== Find ===='   
        write(*,*) 'minVals     ',table_lookup%minVals
        write(*,*) 'maxVals     ',table_lookup%maxVals
        write(*,*) 'steps       ',table_lookup%steps    
        write(*,*) 'lookfor     ',lookfor   
        write(*,*) '===='
    
    endif 
   
    do i=1,table_lookup%dimsCount
       
	    if (table_lookup%isLog10(i)) then 
            zlookfor = dlog10(lookfor(i))
        else
            zlookfor = lookfor(i)
        end if    


	    if (table_lookup%idims(i) > 1 ) then 
            if (table_lookup%minval_minus_epsilon(i) >  zlookfor) then 
                inrange(i)  = -1 
	  		    lowIdxs(i)  =  1
	  		    upIdxs(i)   =  1
	  		    lowWghts(i) =  1.
            elseif ((table_lookup%maxval_plus_epsilon(i) <  zlookfor )) then 
	            inrange(i)  = +1
	  		    lowIdxs(i)  = table_lookup%idims(i)
	  		    upIdxs(i)   = table_lookup%idims(i)
	  		    lowWghts(i) = 1.			
	        else 
	            inrange(i) = 0
                lowIdxs(i) = FLOOR((zlookfor-table_lookup%minVals(i))/table_lookup%steps(i)) + 1
                if (lowIdxs(i) < table_lookup%idims(i) .and. lowIdxs(i) > 0 ) then 
                    if (table_lookup%isLog10(i)) then 
                        zlowloc     = 10**(table_lookup%minVals(i) + ((lowIdxs(i) - 1)*table_lookup%steps(i)))
                        znextlloc   = 10**(table_lookup%minVals(i) + ((lowIdxs(i))*table_lookup%steps(i)))
                        lowWghts(i) = 1. - (lookfor(i) - zlowloc)/(znextlloc - zlowloc)
                    else
                        lowWghts(i) = 1.                                                                         &
                        - (lookfor(i) - (table_lookup%minVals(i) + ((lowIdxs(i) - 1)*table_lookup%steps(i))))    &
                        / (table_lookup%steps(i))
                    end if               
                    upIdxs(i) = lowIdxs(i) + 1
                elseif ( lowIdxs(i) < 1 ) then 
				    lowIdxs(i) = 1
	                upIdxs(i)  =  1   
                    lowWghts(i) = 1
                else 				
                    lowIdxs(i) =  table_lookup%idims(i)  
                    upIdxs(i)  =  table_lookup%idims(i)  
                    lowWghts(i) = 1
                endif
	        end if 
        else 
	        inrange(i)  = 0
	  	    lowIdxs(i)  = 1
	  	    upIdxs(i)   = 1
	  	    lowWghts(i) = 1.
	    end if 
	  
    end do
    if (debug) then 
        write(*,*) 'lowIdxs    ',lowIdxs 
        write(*,*) 'upIdxs     ',upIdxs
        write(*,*) 'lowWghts   ',lowWghts    
        write(*,*) 'inrange    ',inrange   
        write(*,*) '===='
    endif 
   
end subroutine type_table_lookup_find_surrounding_indices 


! With  descriptor file 
subroutine load_from_bin_file(descriptor_file_path,table_bin_file_path,table_lookup)
   implicit none 
   character(len=*), intent(in):: descriptor_file_path
   character(len=*), intent(in):: table_bin_file_path
   type(type_table_lookup), intent(out) :: table_lookup
   !
   integer :: i,j,max_dim,imax_dim,irec,ivar,chunk_size,itemp 
   real :: f
   character(len=200) :: bin_file_path_in_desc
   character(len=200) :: ctemp
   real*4 , allocatable :: temp(:)
   !
   integer :: rc ! error code  
   !
   if (debug) write(*,*) '==== Entering load_from_text_file ===='
   !
   open(100,file=trim(descriptor_file_path),status="old", iostat=rc)
   if (rc /= 0) stop "Error: failed to open file: "//trim(descriptor_file_path)

   
   call type_table_lookup_finialize(table_lookup)
 
   read(100,*)    !"Dep Vars Count (depCount)" 
   read(100,*)     table_lookup%varsCount
   read(100,*)    !"Var Names     " 
   
   allocate(table_lookup%depVarNames(table_lookup%varsCount))

   do i=1,table_lookup%varsCount
     read(100,*)  table_lookup%depVarNames(i) 
   end do 
  
   read(100,*)     !"Var Units" 
   do i=1,table_lookup%varsCount
 	 read(100,*)  
   end do 
   
  
   
   read(100,*)    !"Indep Vars Count (dimCount)"  
   read(100,*)     table_lookup%dimsCount
   
  allocate(table_lookup%idims(table_lookup%dimsCount))

  
   read(100,*)    !"Dims           "
   do i=1,table_lookup%dimsCount
      read(100,*)  table_lookup%idims(i)
   enddo
   

   allocate(table_lookup%minVals(table_lookup%dimsCount))
   allocate(table_lookup%maxVals(table_lookup%dimsCount))
   allocate(table_lookup%steps(table_lookup%dimsCount))
   allocate(table_lookup%isLog10(table_lookup%dimsCount))
   allocate(table_lookup%isVapour(table_lookup%dimsCount))
   allocate(table_lookup%indepVarNames(table_lookup%dimsCount))
   allocate(table_lookup%units(table_lookup%dimsCount))
   
   allocate(table_lookup%maxval_plus_epsilon(table_lookup%dimsCount))
   allocate(table_lookup%minval_minus_epsilon(table_lookup%dimsCount))
   
   
   read(100,*)    !"Indep Vars Names"
   do i=1,table_lookup%dimsCount
      read(100,*)  table_lookup%indepVarNames(i)
   enddo
   
   read(100,*)    !"Indep Vars Units"  
   do i=1,table_lookup%dimsCount
      read(100,*)  table_lookup%units(i)
   enddo
   
   
   read(100,*)         ! "minVals,maxVals" 
   do i=1,table_lookup%dimsCount
      read(100,*) table_lookup%minVals(i),table_lookup%maxVals(i)
   enddo
   
   
    do i=1,table_lookup%dimsCount
       if ( table_lookup%minVals(i) > table_lookup%maxVals(i)) then 
           write(*,*) "min: ", table_lookup%minVals(i),  "   max: table_lookup%maxVals(i)"
           write(*,*) "Min value is larger than maximum value (see above). Will abort."
           stop 1
       end if 

       if ( table_lookup%idims(i) > 1 ) then 
	      table_lookup%steps(i) = (table_lookup%maxVals(i) - table_lookup%minVals(i))/(table_lookup%idims(i)-1)
	   else 
	       write(*,*) "Diminsion should be larger than 1. Will abort."
           stop 1 
	   end if 
	   
	   table_lookup%maxval_plus_epsilon(i)  = table_lookup%maxVals(i) + abs(table_lookup%maxVals(i))/10000.
	   table_lookup%minval_minus_epsilon(i) = table_lookup%minVals(i) - abs(table_lookup%minVals(i))/10000.
    enddo

   
   read(100,*)  !"isLog10"
   
   do i=1,table_lookup%dimsCount
      read(100,*) table_lookup%isLog10(i)
   enddo
   
   read(100,*)  !"isVapour"
   
   do i=1,table_lookup%dimsCount
      read(100,*) table_lookup%isVapour(i)
   enddo
   
   
   read(100,*) !'BIN FILE NAME'
   read(100,*) bin_file_path_in_desc
   
   if ( trim(bin_file_path_in_desc) /= trim(table_bin_file_path)) then   
      write(*,*) "bin file name mismatch. Please check ", trim(bin_file_path_in_desc) , "   ",trim(table_bin_file_path)
   end if 
   
   read(100,*)  !"totalCount   "
   read(100,*)  table_lookup%valuesCount
       
   allocate(table_lookup%tbl(table_lookup%valuesCount,table_lookup%varsCount))
 
 
   max_dim = maxVal(table_lookup%idims)
   chunk_size = (table_lookup%valuesCount/max_dim)/table_lookup%varsCount
   allocate(temp(chunk_size))
   ! Daniel Todo check how it works 
   
   
   ! TODO: Remove this commented line
   !allocate(table_lookup%tbl(table_lookup%valuesCount,table_lookup%varsCount))
   !allocate(temp(table_lookup%valuesCount,table_lookup%varsCount))
 
   open(101,file=trim(table_bin_file_path),status="old",form="unformatted",access="direct",&
   recl = 4*chunk_size, iostat=rc)
   if (rc /= 0) stop "Error: failed to open file: "//trim(table_bin_file_path)

   irec = 0 
   itemp = 0 
   do ivar =1,table_lookup%varsCount
       do imax_dim=1,max_dim
	     irec = irec + 1 
         read(101,rec=irec) temp    
		 table_lookup%tbl((1+(irec-1)*chunk_size):(irec*chunk_size),ivar) = temp(:)
      	 itemp = itemp + chunk_size
       end do
   
   end do 
   ! TODO: Remove this commented line 
   !write(*,*) "Hello", max_dim , table_lookup%valuesCount , itemp,  mod(table_lookup%valuesCount, max_dim)
 
   
   deallocate(temp)
    
   !
   if (debug) then 
       write(*,*) '==== Loaded ===='
       write(*,*) 'bin file    ',trim(bin_file_path_in_desc)
       write(*,*) 'dimsCount   ',table_lookup%dimsCount
       write(*,*) 'depVars Names   ',table_lookup%depVarNames
       write(*,*) 'indepVars Names   ',table_lookup%indepVarNames
       write(*,*) 'indepVars Units   ',table_lookup%units    	   
       write(*,*) 'idims       ',table_lookup%idims
       write(*,*) 'varsCount   ',table_lookup%varsCount
       write(*,*) 'valuesCount ',table_lookup%valuesCount
       write(*,*) 'minVals     ',table_lookup%minVals
       write(*,*) 'maxVals     ',table_lookup%maxVals
       write(*,*) 'isLog10     ',table_lookup%isLog10    
       write(*,*) 'isVapour    ',table_lookup%isVapour    
       write(*,*) '===='
   endif 
   
   close(100)
   close(101)
   
end subroutine load_from_bin_file

subroutine type_table_lookup_find_interpolation_weights(lowIdxs,upIdxs,lowWghts,comIdxs,interpWghts)
    implicit none
    ! Note: This subroutine assumes that all indices and wights are correct   
    integer(kind=ik), dimension(:), intent(in) :: lowIdxs,upIdxs ! independent
    real(kind=rk), dimension(:), intent(in) :: lowWghts ! independent
    !
    integer(kind=ik), allocatable, dimension(:,:), intent(out) :: comIdxs ! independent Thse first index is for dependent variable   
    real(kind=rk), allocatable, dimension(:), intent(out) :: interpWghts    
    !
    integer :: i,k,n,nout,nh,p,q
    logical :: islow    
    !
    if (debug) write(*,*) '==== Entering type_table_lookup_find_interpolation_weights ===='
    !
    n = size(lowIdxs,1)

    nout = 2**n
    if (.not. allocated(interpWghts))  allocate(interpWghts(nout))
    if (.not. allocated(comIdxs)) allocate(comIdxs(n,nout))
    
    interpWghts  = 1._rk
    
    nh = nout/2  
    do i=1,n
        isLow = .false.
        k = 0
        do p=1,nout/nh
            isLow = .not. isLow
            do q=1,nh       
                k = k + 1
                if (isLow) then 
                    comIdxs(i,k) = lowIdxs(i)
                    interpWghts(k) = interpWghts(k)*lowWghts(i)
                else
                    comIdxs(i,k) = upIdxs(i)
                    interpWghts(k) = interpWghts(k)*(1._rk - lowWghts(i))
                endif
            enddo 
        enddo
        nh = nh / 2 
    enddo
    if (debug) then 
	    write(*,*) 'interploation weights:',interpWghts 
        !write(*,*)
        !write(*,*) 'comIdxs(i,k): '
		!
	    !do i=1,n
	    !    do k=1,nout
	    !        write(*,*) i, k , comIdxs(i,k)
		!    end do 
		!end do 
	end if 
    		
end subroutine type_table_lookup_find_interpolation_weights


subroutine type_table_lookup_interpolate(table_lookup,tblIdxs,interpWghts,interpVals)
    ! This function uses the output from calcInterpolationWeights
    implicit none
    type(type_table_lookup), intent(in) :: table_lookup      
    integer(kind=ik), dimension(:), intent(in) :: tblIdxs ! dependent calculated from type_table_lookup_find_table_indices  
    real(kind=rk), dimension(:), intent(in) :: interpWghts 
    real(kind=rk), allocatable, dimension(:), intent(inout) :: interpVals ! interpolated output values for all variables
    !
    integer :: k,nw,p
    !
    if (debug) write(*,*) '==== Entering type_table_lookup_interpolate ===='
    !    
    if (.not. allocated(interpVals)) allocate(interpVals(table_lookup%varsCount))
    nw = 2**table_lookup%dimsCount
    interpVals = 0.0 
    do p=1,table_lookup%varsCount
        do k=1,nw
            if (debug) then 
               write(*,*) "var", p, "point", k, "location", tblIdxs(k),"value",table_lookup%tbl(tblIdxs(k),p)
            end if 
            interpVals(p) = interpVals(p) + interpWghts(k)*table_lookup%tbl(tblIdxs(k),p)         
        end do 
    enddo
    
    
end subroutine type_table_lookup_interpolate


end module mo_table_lookup

 