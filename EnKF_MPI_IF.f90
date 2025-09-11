module EnKF_MPI_IF
    use mpi
    USE,INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER       !用于邻居通信的指针
    implicit none

    private

    integer                :: prcnumb, myprcid
    integer                :: local_communicator
    logical                :: OnMonitor

    ! public interfaces
    public :: EnKF_MPI_init, EnKF_MPI_Finalize

    ! public variables
    public :: local_communicator, prcnumb, myprcid, OnMonitor

    contains

    subroutine EnKF_MPI_init()
       integer :: ierr

       call MPI_Init(ierr)
       call MPI_Comm_rank(MPI_COMM_WORLD, myprcid, ierr)
       call MPI_Comm_size(MPI_COMM_WORLD, prcnumb, ierr)

       local_communicator = MPI_COMM_WORLD
       
       if (myprcid .eq. 0) then
          OnMonitor = .true.
       else
          OnMonitor = .false.
       end if

    end subroutine EnKF_MPI_init

    subroutine EnKF_MPI_Finalize(ierr)
       integer, intent(out) :: ierr

       call MPI_Finalize(ierr)

    end subroutine EnKF_MPI_Finalize
    
end module EnKF_MPI_IF