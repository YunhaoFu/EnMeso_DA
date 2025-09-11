module enkf_mpi
   use mpi
   use enkf_obs
   implicit none

   private

   public :: enkf_bcast, enkf_bcast_init
   public :: enkf_bcast_check

contains

   subroutine enkf_bcast_init(pe_rank,pe_root,obs)
      integer             , intent(in)                :: pe_rank
      integer             , intent(in)                :: pe_root
      type(obs_instrument), intent(inout)             :: obs

      integer                                         :: ierr

      call mpi_bcast(obs%insn             ,len(obs%insn) ,mpi_char   ,pe_root,mpi_comm_world,ierr)
      call mpi_bcast(obs%nbatch           ,1             ,mpi_integer,pe_root,mpi_comm_world,ierr)
      call mpi_bcast(obs%nobs             ,1             ,mpi_integer,pe_root,mpi_comm_world,ierr)
      call mpi_bcast(obs%nens             ,1             ,mpi_integer,pe_root,mpi_comm_world,ierr)
     !call mpi_bcast(obs%nallocated       ,1             ,mpi_integer,pe_root,mpi_comm_world,ierr)
      call mpi_bcast(obs%nobs_without_qc  ,1             ,mpi_integer,pe_root,mpi_comm_world,ierr)
      call mpi_bcast(obs%avail            ,1             ,mpi_logical,pe_root,mpi_comm_world,ierr)
      if (pe_rank /= pe_root) obs%nallocated = 0
   endsubroutine enkf_bcast_init

   subroutine enkf_bcast(pe_rank,pe_root,obs,size_read,step,nstep)
      integer             , intent(in)                :: pe_rank
      integer             , intent(in)                :: pe_root
      type(obs_instrument), intent(inout)             :: obs
      integer             , intent(in)                :: size_read
      integer             , intent(in)                :: step
      integer             , intent(in)                :: nstep

      integer                                         :: ierr

      integer                                         :: nallocated

      nallocated = obs%nallocated
      call mpi_bcast(obs%nallocated,1,mpi_integer,pe_root,mpi_comm_world,ierr)

      ! expand for other cpus in the first call
      if (nallocated == 0) then
         allocate(obs%lon  (obs%nallocated))
         allocate(obs%lat  (obs%nallocated))
         allocate(obs%alt  (obs%nallocated))
         allocate(obs%oid  (obs%nallocated))
         allocate(obs%err  (obs%nallocated))
         allocate(obs%qc   (obs%nallocated,obs%nens))
         allocate(obs%hx   (obs%nallocated,obs%nens))
         allocate(obs%inno (obs%nallocated))
      endif

      nallocated = ubound(obs%lon, dim=1)

      if (obs%nobs+size_read > nallocated) then
         ! expand for other cpus

         if (step < nstep) then
            obs%nallocated = ceiling(real(obs%nallocated) * real(nstep) / real(step))
         else
            obs%nallocated = obs%nobs+size_read
         endif

         call expand(obs%lon   ,obs%nallocated)
         call expand(obs%lat   ,obs%nallocated)
         call expand(obs%alt   ,obs%nallocated)
         call expand(obs%oid   ,obs%nallocated)
         call expand(obs%err   ,obs%nallocated)
         call expand(obs%qc    ,obs%nallocated)
         call expand(obs%hx    ,obs%nallocated)
         call expand(obs%inno  ,obs%nallocated)

      endif

      call mpi_bcast(obs%lon  (obs%nobs+1:obs%nobs+size_read)  ,size_read,          mpi_real   ,pe_root,mpi_comm_world,ierr)
      call mpi_bcast(obs%lat  (obs%nobs+1:obs%nobs+size_read)  ,size_read,          mpi_real   ,pe_root,mpi_comm_world,ierr)
      call mpi_bcast(obs%alt  (obs%nobs+1:obs%nobs+size_read)  ,size_read,          mpi_real   ,pe_root,mpi_comm_world,ierr)
      call mpi_bcast(obs%oid  (obs%nobs+1:obs%nobs+size_read)  ,size_read,          mpi_integer,pe_root,mpi_comm_world,ierr)
      call mpi_bcast(obs%err  (obs%nobs+1:obs%nobs+size_read)  ,size_read,          mpi_real   ,pe_root,mpi_comm_world,ierr)
      call mpi_bcast(obs%qc   (obs%nobs+1:obs%nobs+size_read,:),size_read*obs%nens, mpi_integer,pe_root,mpi_comm_world,ierr)
      call mpi_bcast(obs%hx   (obs%nobs+1:obs%nobs+size_read,:),size_read*obs%nens, mpi_real   ,pe_root,mpi_comm_world,ierr)
     !call mpi_bcast(obs%inno (obs%nobs+1:obs%nobs+size_read)  ,size_read,          mpi_real   ,pe_root,mpi_comm_world,ierr)
   endsubroutine enkf_bcast

   subroutine enkf_bcast_check(pe_rank,nins,nens,radius,infl,infl_rtpp,verbose,obs)
      integer             , intent(in) :: pe_rank
      integer             , intent(in) :: nins
      integer             , intent(in) :: nens
      real                , intent(in) :: radius
      real                , intent(in) :: infl,infl_rtpp
      logical             , intent(in) :: verbose
      type(obs_instrument), intent(in) :: obs(nins)

      character(len=64)                :: filename
      integer                          :: fileunit

      integer                          :: i, j, ierr, n1

      call mpi_barrier(mpi_comm_world, ierr)
      write(filename,'(A,I0)') 'output_XXX_file_', pe_rank
      fileunit = 100 + pe_rank
      open(unit=fileunit,file=trim(filename),status='replace',action='write',iostat=ierr)
      write(fileunit,*) nins, nens, radius, infl, infl_rtpp, verbose
      do i=1,nins
         write(fileunit,*) obs(i)%insn, obs(i)%nens, obs(i)%nobs, obs(i)%avail
      enddo
      close(fileunit)

      call mpi_barrier(mpi_comm_world, ierr)
      write(filename,'(A,I0)') 'output_lon_file_', pe_rank
      fileunit = 100 + pe_rank
      open(unit=fileunit,file=trim(filename),status='replace',action='write',iostat=ierr)
      write(fileunit,*)
      do i=1,nins
         write(fileunit,*) obs(i)%lon(1:obs(i)%nobs)
      enddo
      close(fileunit)

      call mpi_barrier(mpi_comm_world, ierr)
      write(filename,'(A,I0)') 'output_lat_file_', pe_rank
      fileunit = 100 + pe_rank
      open(unit=fileunit,file=trim(filename),status='replace',action='write',iostat=ierr)
      do i=1,nins
         write(fileunit,*) obs(i)%lat(1:obs(i)%nobs)
      enddo
      close(fileunit)

      call mpi_barrier(mpi_comm_world, ierr)
      write(filename,'(A,I0)') 'output_alt_file_', pe_rank
      fileunit = 100 + pe_rank
      open(unit=fileunit,file=trim(filename),status='replace',action='write',iostat=ierr)
      do i=1,nins
         write(fileunit,*) obs(i)%alt(1:obs(i)%nobs)
      enddo
      close(fileunit)

      call mpi_barrier(mpi_comm_world, ierr)
      write(filename,'(A,I0)') 'output_hx_file_', pe_rank
      fileunit = 100 + pe_rank
      open(unit=fileunit,file=trim(filename),status='replace',action='write',iostat=ierr)
      do i=1,nins
         write(fileunit,*) obs(i)%hx(1:obs(i)%nobs,1:nens)
      enddo
      close(fileunit)

      call mpi_barrier(mpi_comm_world, ierr)
      write(filename,'(A,I0)') 'output_err_file_', pe_rank
      fileunit = 100 + pe_rank
      open(unit=fileunit,file=trim(filename),status='replace',action='write',iostat=ierr)
      do i=1,nins
         write(fileunit,*) obs(i)%err(1:obs(i)%nobs)
      enddo
      close(fileunit)

      call mpi_barrier(mpi_comm_world, ierr)
      write(filename,'(A,I0)') 'output_inno_file_', pe_rank
      fileunit = 100 + pe_rank
      open(unit=fileunit,file=trim(filename),status='replace',action='write',iostat=ierr)
      do i=1,nins
         write(fileunit,*) obs(i)%inno(1:obs(i)%nobs)
      enddo
      close(fileunit)
   endsubroutine enkf_bcast_check

end module enkf_mpi
