program EnKF_Meso
    use mpi
    use kdtree
    use enkf_obs
    use enkf_mpi
    use EnKF_MPI_IF
    use EnKF_DM
    use EnKF_IO
    use letkf_driver
    implicit none

    integer                            :: nins        ! instrument count
    integer                            :: batch_size
    integer                            :: nens        ! ensemble size
    real                               :: radius      ! localization r [km]
    real                               :: infl        ! inflation >= 1.
    real                               :: infl_rtpp
    logical                            :: verbose     !
    logical                            :: write_y     !
    type(obs_instrument), allocatable  :: obs(:)      ! insc
    type(kd_root)                      :: bkg_kdtree

    integer                            :: begin_read, size_read
    integer                            :: nobs_min, nobs_max
    integer                            :: i, batch, root = 0
    integer                            :: ierr
    double precision                   :: t0, t1
    
    call EnKF_MPI_init()

    call enkf_cfg("EnKF_Meso.nml", obs, nins, batch_size, nens, radius, infl, infl_rtpp, verbose, write_y, myprcid)

    call grid%model_init("EnKF_Meso.nml")

    call DM_Init()

    call DM_Alloc(xs=grid%xs, xe=grid%xe, ys=grid%ys, ye=grid%ye, zs=grid%zs, ze=grid%ze)

    call EnKF_IO_init()

    t0 = MPI_Wtime()

    call kd_init(bkg_kdtree, reshape(grid%lon, (/grid%nx * grid%ny/)), reshape(grid%lat, (/grid%nx * grid%ny/)))

    do i=1,nins
       if (OnMonitor) then
          write(*, "('read instrument(', I0, '): ', A)") i, trim(obs(i)%insn)
          call obs(i)%obs_init(batch_size)
       endif

       call enkf_bcast_init(myprcid,root,obs(i))

       if (obs(i)%avail) then
          do batch = 1, obs(i)%nbatch
             begin_read = (batch-1)*batch_size + 1
             size_read  = min(begin_read + batch_size - 1, obs(i)%nobs_without_qc) - begin_read + 1
             if (OnMonitor) then
                call obs(i)%obs_read(begin_read,size_read,batch_size,batch)
             endif

             call enkf_bcast(myprcid,root,obs(i),size_read,batch,obs(i)%nbatch)
             call obs_partition(obs(i),bkg_kdtree,radius*1000,myprcid,size_read) ! will update size_read
             call obs(i)%obs_qc(myprcid,size_read)                               ! will update obs(i)%nobs
          enddo
       endif

       obs(i)%avail = obs(i)%nobs > 0

       ! shrink
       call expand(obs(i)%lon  ,obs(i)%nobs)
       call expand(obs(i)%lat  ,obs(i)%nobs)
       call expand(obs(i)%alt  ,obs(i)%nobs)
       call expand(obs(i)%oid  ,obs(i)%nobs)
       call expand(obs(i)%err  ,obs(i)%nobs)
       call expand(obs(i)%qc   ,obs(i)%nobs)
       call expand(obs(i)%hx   ,obs(i)%nobs)
       call expand(obs(i)%inno ,obs(i)%nobs)

       if (write_y) call obs(i)%obs_write(myprcid)

       call MPI_Reduce(obs(i)%nobs, nobs_min, 1, MPI_INTEGER, MPI_MIN, root, local_communicator, ierr)
       call MPI_Reduce(obs(i)%nobs, nobs_max, 1, MPI_INTEGER, MPI_MAX, root, local_communicator, ierr)
       if (OnMonitor) write(*, "('obs(', I0 ,')%nobs between ', I0, ' and ', I0, ' nobs_without_qc=', I0)") i, nobs_min, nobs_max, obs(i)%nobs_without_qc
    enddo

    if (verbose) call enkf_bcast_check(myprcid,nins,nens,radius,infl,infl_rtpp,verbose,obs)

    ! build obs tree
    do i=1,nins
       call obs(i)%kdtree_init()
    enddo

    t1 = MPI_Wtime() - t0
    call MPI_Reduce(t1, t0, 1, MPI_DOUBLE_PRECISION, MPI_MAX, root, local_communicator, ierr)
    if (OnMonitor) write(*, "('processing observations took ', F10.3, ' seconds')") t0

    t0 = MPI_Wtime()
    call EnKF_IO_Read()
    t1 = MPI_Wtime() - t0
    call MPI_Reduce(t1, t0, 1, MPI_DOUBLE_PRECISION, MPI_MAX, root, local_communicator, ierr)
    if (OnMonitor) write(*, "('reading background took ', F10.3, ' seconds')") t0

    t0 = MPI_Wtime()
    call letkf(nlon=grid%nx, nlat=grid%ny, lon=grid%lon, lat=grid%lat, nlev=grid%nz, nens=nens, nins=nins, y=obs, radius=radius*1000., infl=infl, var=buf4d_assimilated)
    t1 = MPI_Wtime() - t0
    call MPI_Reduce(t1, t0, 1, MPI_DOUBLE_PRECISION, MPI_MAX, root, local_communicator, ierr)
    if (OnMonitor) write(*, "('letkf took ', F10.3, ' seconds')") t0

    t0 = MPI_Wtime()
    call EnKF_IO_Write()
    t1 = MPI_Wtime() - t0
    call MPI_Reduce(t1, t0, 1, MPI_DOUBLE_PRECISION, MPI_MAX, root, local_communicator, ierr)
    if (OnMonitor) write(*, "('writing analysis took ', F10.3, ' seconds')") t0

    call kd_free(bkg_kdtree)

    do i=1,nins
        call obs(i)%kdtree_final()
        call obs(i)%obs_final()
    enddo
    call enkf_final(obs)

    call DM_Finalize()

    call EnKF_MPI_Finalize(ierr)

    if (ierr == MPI_SUCCESS .and. OnMonitor) then
        write(*, *) "EnKF_Meso successfully finished."
    endif

end program EnKF_Meso
