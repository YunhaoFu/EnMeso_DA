module EnKF_IO

    use mpi
    use EnKF_MPI_IF
    use EnKF_DM

    implicit none
    
    private

    ! public interfaces
    public :: EnKF_IO_init, EnKF_IO_Read, EnKF_IO_Write

    ! public variables
    public :: buf4d_assimilated

    ! parameters
    integer, parameter :: i_params = 14

    ! local variables

    integer :: ids, ide, jds, jde, kds, kde
    integer :: ims, ime, jms, jme, kms, kme
    integer :: its, ite, jts, jte, kts, kte

    integer :: e_vert, num_soil_layers
    integer :: e_we, e_sn
    real(kind=4) :: xs_we, ys_sn, xd, yd

    integer :: total_lev, levs, leve
    integer,allocatable :: is(:),ie(:),js(:),je(:),level_start(:),level_end(:), levels(:)

    integer(kind=mpi_offset_kind), protected :: disp
    real(kind=4),allocatable :: buf3d(:,:,:),buf_send(:),buf_recv(:),buf1d(:)
    real(kind=4),allocatable :: buf_hrlon(:,:,:), buf_hrztl(:,:,:)
    integer,allocatable        :: send_displ_1(:), send_count_1(:), recv_displ_1(:), recv_count_1(:) 
    integer,allocatable        :: send_displ_2(:), send_count_2(:), recv_displ_2(:), recv_count_2(:) 
    integer,dimension(3)   ::substar, globsize, subsize
    integer          :: filetype
    integer          :: isize
    integer :: kp

logical, dimension(:), allocatable :: assimilated

    ! MPI variables
    integer :: ierr, fh, stats
    integer :: root = 0

    real(kind=4),allocatable :: buf4d(:,:,:,:)

    ! pubclic variables
    real(kind=4),allocatable :: buf4d_assimilated(:,:,:,:)

contains

    subroutine EnKF_IO_init()

        character(len=256) :: inputfn
        integer :: ihead(i_params)
        integer :: i, j, ens
        real(kind=4), allocatable, dimension(:) :: rhead

        call get_ijk_from_grid (    grid ,                           &
                                    ids, ide, jds, jde, kds, kde,    &
                                    ims, ime, jms, jme, kms, kme,    &
                                    its, ite, jts, jte, kts, kte     )

        ! print *, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, myprcid

        ! print *, its, ite, jts, jte

        ! allocate(grid%lon_globe(ids:ide,jds:jde))
        ! allocate(grid%lat_globe(ids:ide,jds:jde))
        allocate(grid%lon(its:ite,jts:jte))
        allocate(grid%lat(its:ite,jts:jte))

        grid%nx = ite - its + 1
        grid%ny = jte - jts + 1

        ens = 1
        write(inputfn,'(A,I3.3)') "model/xa.", ens

        disp=0
        if (myprcid_b == 0)then
          call MPI_FILE_OPEN(nscomm,trim(inputfn),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
          call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER,MPI_INTEGER,'native',MPI_INFO_NULL,ierr)
        endif

        !ihead
        if (OnMonitor) then
            call MPI_FILE_READ(fh,ihead,i_params,MPI_INTEGER,stats,ierr)
        endif

        call MPI_BCAST ( ihead,i_params,mpi_integer,0,local_communicator,ierr)

        e_we   = ihead(10)
        e_sn   = ihead(11)
        e_vert = ihead(12)
        num_soil_layers = ihead(14)

        !rhead
        disp = i_params*4
        kp = 5+(e_vert+1)+2*(kme-kms+1)+2*num_soil_layers
        allocate(rhead(kp))

        if (myprcid_b == 0) call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,MPI_REAL,'native',MPI_INFO_NULL,ierr)
        if (myprcid == 0) then
            call MPI_FILE_READ(fh,rhead,kp,MPI_REAL,stats,ierr )
        endif
        call MPI_BCAST ( rhead,kp,MPI_REAL,root,local_communicator,ierr)

        xs_we = rhead(1)
        ys_sn = rhead(2)
        xd = rhead(3)
        yd = rhead(4)

        deallocate(rhead)

        do i = its, ite
            grid%lon(i,jts:jte) = xs_we + (i-1)*xd
        end do

        do j = jts, jte
            grid%lat(its:ite,j) = ys_sn + (j-1)*yd
        end do

        ! print *, myprcid , its ,minval(grid%lon(:,jts)), ite, maxval(grid%lon(:,jts)), jts, minval(grid%lat(its,:)), jte, maxval(grid%lat(its,:))

        if (myprcid_b == 0)then
            call MPI_FILE_CLOSE(fh,ierr)
        end if

    end subroutine EnKF_IO_init

    subroutine EnKF_IO_Read()
        character(len=256) :: inputfn
        integer :: ihead(i_params)
        integer :: l, m, n, i, j, k, ens
        real(kind=4), allocatable, dimension(:) :: rhead
        integer :: nens

        nens = grid%nens

        disp=disp + kp*4_MPI_OFFSET_KIND
        total_lev = 14*(kme-kms+1)+2*num_soil_layers   !14 3D  of kms~kme and 2 3D of soil
        total_lev = total_lev + 10
        total_lev = total_lev+3

        allocate(assimilated(total_lev))
        assimilated = .false.

        assimilated(289:289+72-1)= .true.  ! u
        assimilated(361:361+72-1)= .true.  ! v
        assimilated(505:505+72-1) = .true. ! th
        assimilated(577:577+72-1) = .true. ! qv
        assimilated(1020:1021) = .true.    ! u10, v10
        assimilated(1022:1023) = .true.    ! t2, q2

        grid%nz = count(assimilated)

        if (OnMonitor) print *, "total_lev = ", total_lev, " grid%nz = ", grid%nz

        allocate(level_start(prcnumb_a))
        allocate(level_end(prcnumb_a))
        allocate(levels(prcnumb_a))
        allocate(is(prcnumb_b))
        allocate(ie(prcnumb_b))
        allocate(js(prcnumb_a))
        allocate(je(prcnumb_a))

!evenly divide total vertical layers by NS direction
        DO n=1, prcnumb_a
          IF(n <= mod(total_lev,prcnumb_a)) THEN
             levels(n) = total_lev/prcnumb_a + 1
          ELSE
             levels(n) = total_lev/prcnumb_a
          ENDIF
        ENDDO
        level_start(1) = 1
        level_end(1) = levels(1)
        DO n=2, prcnumb_a
          level_start(n) = level_end(n-1) + 1
          level_end(n)   = level_start(n) + levels(n) - 1
        ENDDO
   ! start/end index of this processs level
        n = myprcid_a + 1
        levs = level_start(n)
        leve = level_end(n)

!***************************************************

        call MPI_ALLGATHER(its,1,MPI_INTEGER,is,1,MPI_INTEGER,ewcomm,ierr)
        call MPI_ALLGATHER(ite,1,MPI_INTEGER,ie,1,MPI_INTEGER,ewcomm,ierr)
        call MPI_ALLGATHER(jts,1,MPI_INTEGER,js,1,MPI_INTEGER,nscomm,ierr)
        call MPI_ALLGATHER(jte,1,MPI_INTEGER,je,1,MPI_INTEGER,nscomm,ierr)

        allocate(buf1d((ite-its+1)*(jde-jds+1)*(leve-levs+1)))
        allocate(buf_recv((ite-its+1)*(jte-jts+1)*total_lev))
        allocate(buf4d(its:ite,jts:jte,1:total_lev,nens))
        allocate(buf4d_assimilated(its:ite,jts:jte,grid%nz,nens))

        allocate(send_count_1(prcnumb_b))
        allocate(send_displ_1(prcnumb_b))
        allocate(recv_count_1(1))

        allocate(send_count_2(prcnumb_a))
        allocate(send_displ_2(prcnumb_a))
        allocate(recv_count_2(prcnumb_a))
        allocate(recv_displ_2(prcnumb_a))

        if (myprcid_b == 0) then
            allocate(buf3d(ids:ide,jds:jde,levs:leve))
            allocate(buf_send((ide-ids+1)*(jde-jds+1)*(leve-levs+1)))
        else
            allocate(buf_send(1))
        end if

        if (myprcid_b == 0)then
            globsize(1)=ide-ids+1
            globsize(2)=jde-jds+1
            globsize(3)=total_lev

            subsize(1)=ide-ids+1
            subsize(2)=jde-jds+1
            subsize(3)=leve-levs+1

            substar(1)=0
            substar(2)=0
            substar(3)=levs-1

            call mpi_type_create_subarray(3,globsize,subsize,substar,MPI_ORDER_FORTRAN,MPI_REAL,filetype,ierr)
            call mpi_type_commit(filetype,ierr)

            isize=subsize(1)*subsize(2)*subsize(3)
            call MPI_FILE_CLOSE(fh,ierr)
        end if

        recv_count_1(1) = (ite-its+1)*(jde-jds+1)*(leve-levs+1)
        do n=1,prcnumb_b
           send_count_1(n) = (ie(n)-is(n)+1)*(jde-jds+1)*(leve-levs+1)
        enddo
        send_displ_1(1) = 0
        if (prcnumb_b >=2) then
            do n=2,prcnumb_b
               send_displ_1(n) = send_displ_1(n-1) + send_count_1(n-1)
            enddo
        endif

        do n=1, prcnumb_a
          send_count_2(n) = (ite-its+1)*(je(n)-js(n)+1)*(leve-levs+1)
          recv_count_2(n) = (ite-its+1)*(jte-jts+1)*levels(n)
        enddo
        send_displ_2(1) = 0
        recv_displ_2(1) = 0
        if (prcnumb_a >=2)then
            do n=2, prcnumb_a
               send_displ_2(n) = send_displ_2(n-1) + send_count_2(n-1)
               recv_displ_2(n) = recv_displ_2(n-1) + recv_count_2(n-1)
            enddo
        endif

        call MPI_Barrier(local_communicator,ierr)

        do ens = 1, nens
            if (myprcid_b == 0)then
                write(inputfn,'(A,I3.3)') "model/xa.", ens
                print *, "Going to read ", trim(inputfn)
                call MPI_FILE_OPEN(nscomm,trim(inputfn),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
                call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,filetype,'native',MPI_INFO_NULL,ierr)
                call MPI_FILE_READ_ALL(fh,buf3d,isize,MPI_REAL,stats,ierr)
                call MPI_FILE_CLOSE(fh,ierr)

                l = 0
                do n=1, prcnumb_b
                  do m=1, prcnumb_a
                    do k=levs, leve
                      do j=js(m), je(m)
                        do i=is(n), ie(n)
                           l=l+1
                           buf_send(l) = buf3d(i,j,k)
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
            end if

            call MPI_SCATTERV(buf_send,send_count_1,send_displ_1,MPI_REAL,buf1d,recv_count_1(1),MPI_REAL,root,ewcomm,ierr)
            call MPI_ALLTOALLV(buf1d,send_count_2,send_displ_2,MPI_REAL,buf_recv,recv_count_2,recv_displ_2,MPI_REAL,nscomm,ierr)

            l = 0
            do n=1,prcnumb_a
              do k=level_start(n),level_end(n)
                 do j=jts,jte
                    do i=its,ite
                       l = l + 1
                       buf4d(i,j,k,ens) = buf_recv(l)
                    enddo
                 enddo
              enddo
            enddo

            l = 0
            do k = 1, total_lev
                if ( .not. assimilated(k) ) cycle
                l = l + 1
                buf4d_assimilated(:,:,l,ens) = buf4d(:,:,k,ens)
            end do

        end do

        call MPI_Barrier(local_communicator,ierr)

        deallocate(buf1d)
        deallocate(buf_recv)

        deallocate(send_count_1)
        deallocate(send_displ_1)
        deallocate(recv_count_1)

        deallocate(send_count_2)
        deallocate(send_displ_2)
        deallocate(recv_count_2)
        deallocate(recv_displ_2)

        if (myprcid_b == 0) deallocate(buf3d)
        deallocate(buf_send)

    end subroutine EnKF_IO_Read

    subroutine EnKF_IO_Write()

        character(len=256) :: inputfn
        integer :: l, m, n, i, j, k, ens
        integer :: nens

        nens = grid%nens

        allocate ( send_count_2(prcnumb_a))
        allocate ( recv_count_2(prcnumb_a))
        allocate ( send_displ_2(prcnumb_a))
        allocate ( recv_displ_2(prcnumb_a))

        allocate ( send_count_1(1        ) )
        allocate ( recv_count_1(prcnumb_b) )
        allocate ( recv_displ_1(prcnumb_b) )

        allocate( buf_recv((jde-jds+1)*(ite-its+1)*(leve-levs+1) ) )
        allocate( buf_hrlon(its:ite,jds:jde,levs:leve) )
        if (myprcid_b == 0 ) allocate ( buf_hrztl(ids:ide,jds:jde,levs:leve) )

        if (myprcid_b == 0 ) then
            allocate ( buf1d( (ide-ids+1)*(jde-jds+1)*(leve-levs+1)) )
        else
            allocate ( buf1d(1) )
        endif

        do n=1, prcnumb_a
            send_count_2(n) = (ite-its+1)*(jte-jts+1)*levels(n)
            recv_count_2(n) = (je(n)-js(n)+1)*(ite-its+1)*(leve-levs+1)
        end do

        send_displ_2(1)=0
        recv_displ_2(1)=0
        if (prcnumb_a >= 2) then
          do n=2, prcnumb_a
            send_displ_2(n)=send_displ_2(n-1)+send_count_2(n-1)
            recv_displ_2(n)=recv_displ_2(n-1)+recv_count_2(n-1)
          end do
        end if

        send_count_1 = (jde-jds+1)*(ite-its+1)*(leve-levs+1)
        do n=1,prcnumb_b
            recv_count_1(n)= (jde-jds+1)*(ie(n)-is(n)+1)*(leve-levs+1)
        end do
        recv_displ_1(1)=0
        if ( prcnumb_b >= 2 ) then
            do n=2,prcnumb_b
                recv_displ_1(n)= recv_displ_1(n-1) + recv_count_1(n-1 )
            end do
        end if

        call MPI_Barrier(local_communicator,ierr)

        do ens = 1, nens

            l = 0
            do k = 1, total_lev
                if ( .not. assimilated(k) ) cycle
                l = l + 1
                buf4d(:,:,k,ens) = buf4d_assimilated(:,:,l,ens)
            end do

            call MPI_ALLTOALLV(buf4d(:,:,:,ens),send_count_2,send_displ_2,MPI_REAL,buf_recv,recv_count_2,recv_displ_2,MPI_REAL,nscomm,ierr )

            l=0
            do n=1, prcnumb_a
                do k=levs,leve
                    do j=js(n),je(n)
                        do i=its,ite
                            l=l+1
                            buf_hrlon(i,j,k) = buf_recv( l )
                        end do
                    end do
                end do
            end do

            call MPI_GATHERV(buf_hrlon,send_count_1,MPI_REAL,buf1d,recv_count_1,recv_displ_1,MPI_REAL,root,ewcomm,ierr)

            if (myprcid_b == 0) then
                l=0
                do n=1, prcnumb_b
                    do k=levs,leve
                        do j=jds,jde
                            do i=is(n),ie(n)
                                l=l+1
                                buf_hrztl(i,j,k) = buf1d( l )
                            end do
                        end do
                    end do
                end do

                write(inputfn,'(A,I3.3)') "model/xa.", ens
                print *, "Going to write ", trim(inputfn)
                call MPI_FILE_OPEN(nscomm,trim(inputfn),MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
                call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,filetype,'native',MPI_INFO_NULL,ierr)
                CALL MPI_FILE_WRITE_ALL(fh,buf_hrztl,isize,MPI_REAL,stats,ierr)
                call MPI_FILE_CLOSE(fh,ierr)
            end if

        end do

        deallocate(send_count_2)
        deallocate(recv_count_2)
        deallocate(send_displ_2)
        deallocate(recv_displ_2)

        deallocate(send_count_1)
        deallocate(recv_count_1)
        deallocate(recv_displ_1)

        deallocate(buf_recv)
        deallocate(buf_hrlon)
        if (myprcid_b == 0) deallocate(buf_hrztl)

        deallocate(buf1d)

        deallocate(buf4d)
        deallocate(buf4d_assimilated)

        deallocate(level_start)
        deallocate(level_end)
        deallocate(levels)
        deallocate(is)
        deallocate(ie)
        deallocate(js)
        deallocate(je)

        ! deallocate(grid%lon_globe)
        ! deallocate(grid%lat_globe)
        deallocate(grid%lon)
        deallocate(grid%lat)

        deallocate(assimilated)

    end subroutine EnKF_IO_Write
    
end module EnKF_IO
