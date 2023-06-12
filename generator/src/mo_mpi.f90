!> @brief functions that implement MPI interfaces for different data types 
!!
module mo_mpi
  use mo_kind 
implicit none 

   public :: init_mpi
   public :: finalize_mpi
   public :: size_of_cluster
   public :: process_rank
   public :: bcast
private 

  interface init_mpi
     module procedure init
  end interface init_mpi
  
  
  interface finalize_mpi
     module procedure finalize 
  end interface finalize_mpi
  
 
  interface bcast
     module procedure bcastInt4
     module procedure bcastInt8
  end interface bcast
  
  integer :: size_of_cluster
  integer :: process_rank

  include 'mpif.h'	
contains 
subroutine bcastInt4(n)
  use mo_config, only: canWrite,isMaster,canRead
  implicit none 
   integer(kind=ik4), intent(inout) :: n 
   integer  ::  mpistat(mpi_status_size), ierror
      call MPI_BCAST(n, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD,ierror);
end subroutine bcastInt4

subroutine bcastInt8(n)
  use mo_config, only: canWrite,isMaster,canRead
  implicit none 
   integer(kind=ik8), intent(inout) :: n 
   integer  ::  mpistat(mpi_status_size), ierror
      call MPI_BCAST(n, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD,ierror)
end subroutine bcastInt8

subroutine finalize
  implicit none 
  integer:: ierror
  ! Body
  call MPI_FINALIZE(ierror)
  call process_mpi_code(ierror)

end subroutine finalize
subroutine init
  use mo_config, only: canWrite,isMaster,canRead
implicit none 
  
  ! Local Vatriables 
  integer :: ierror
  call MPI_INIT(ierror)
  call process_mpi_code(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size_of_cluster, ierror)
  call process_mpi_code(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, process_rank, ierror)
  call process_mpi_code(ierror)
  
  canRead  = process_rank == 0 
  canWrite = process_rank == 0 
  isMaster = process_rank == 0 

end subroutine init

subroutine process_mpi_code(code)

implicit none 
 
  integer, intent(in) :: code 
  integer  :: ierror 
  integer errorcode, resultlen
  character(len=300) string
  if (code /= MPI_SUCCESS) then 
    call MPI_ERROR_STRING(errorcode, string, resultlen, ierror)	
	if (ierror == MPI_SUCCESS ) write(*,*) trim(string(1:resultlen))
    call MPI_ABORT(MPI_COMM_WORLD,ierror)
    stop 
  end if 
end subroutine process_mpi_code
end module mo_mpi