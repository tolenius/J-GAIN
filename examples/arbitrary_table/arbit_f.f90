program arbitrary
implicit none

integer, parameter :: indepCount = 4 ! independent variables count
integer, parameter :: depCount = 3 ! dependent variables count
integer, parameter:: d1 = 5, d2 = 9, d3 = 17 , d4 = 33  ! independent dimensions 
integer :: dims
real(kind=4) :: table(d1*d2*d3*d4,depCount)
integer:: k ! index counter
integer::i1,i2,i3,i4 
integer::j ! dependent variables 
!
integer:: rc 
do j = 1, depCount
  k = 0
  do i1=1,d1  ! slowest indep  
    do i2=1,d2
      do i3=1,d3
        do i4=1,d4  ! fastest indep
		 k = k + 1 
  	     table(k,j) = j + 10*i1 + 100*i2 + 1000*i3 + 10000*i4  ! a facked formula 
        end do 
      end do 
    end do 
  end do 
end do 

open(unit = 1000, file = 'kookup_3dep_4indep_fortran.bin',&
    &status="unknown",form="unformatted", access = "direct", recl = 4*size(table), iostat=rc)
if (rc /= 0) stop "Error: failed to open file to write"

write(1000,rec=1) table

close(1000)

end program arbitrary