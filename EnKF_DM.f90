module EnKF_DM
    use MPI
    use EnKF_MPI_IF
    use, intrinsic :: iso_fortran_env, only : stderr => error_unit
    implicit none

    private

    ! public interfaces

    public :: DM_Init, DM_Finalize
    public :: DM_Alloc
    public :: get_ijk_from_grid

    ! public variables

    public :: grid, myprcid_a, myprcid_b
    public :: nscomm, ewcomm
    public :: prcnumb_a, prcnumb_b

    integer :: prcnumb_a ! Y
    integer :: prcnumb_b ! X
    integer :: myprcid_a ! Y
    integer :: myprcid_b ! X
    integer, dimension(:,:), allocatable :: procmap, cordmap   ! cordmap和procmap相反，以进程号为索引，存储x、y坐标（具体参考cordmap的内存分配）
    integer, dimension(:), allocatable :: ewmember
    integer, dimension(:), allocatable :: nsmember
    integer                            :: ewgroup, ewcomm
    integer                            :: nsgroup, nscomm
    integer                            :: ewcommsize, nscommsize
    integer                            :: world_group
    integer, dimension(:), allocatable :: ewcomm_lon, ewgroup_lon
    integer, dimension(:), allocatable :: nscomm_lat, nsgroup_lat
    integer, dimension(:,:),allocatable:: istedp, jstedp      ! 1 for i_start(j_start), 2 for i_end(j_end)

    integer, parameter ::  shw_control = 4
    integer, parameter ::  bdyzone_x = 4
    integer, parameter ::  bdyzone_y = 4

    type :: domain

        ! real(kind=8), allocatable, dimension(:,:) :: lon_globe, lat_globe
        real(kind=4), allocatable, dimension(:,:) :: lon      , lat

        integer     ::  nx, ny, nz

        integer     ::  sd31,   ed31,   sd32,   ed32,   sd33,   ed33,         &
                        sd21,   ed21,   sd22,   ed22,                         &
                        sd11,   ed11

        integer     ::  sp31,   ep31,   sp32,   ep32,   sp33,   ep33,         &
                        sp21,   ep21,   sp22,   ep22,                         &
                        sp11,   ep11,                                         &
                        sm31,   em31,   sm32,   em32,   sm33,   em33,         &
                        sm21,   em21,   sm22,   em22,                         &
                        sm11,   em11                                         

        integer     ::  nens
        integer     ::  nvar
        character(len=30), dimension(:), allocatable :: varn_list
        integer     ::  xs, xe, ys, ye, zs, ze

        contains

            procedure :: model_init => read_nml

    end type domain

    type(domain) :: grid

contains

    subroutine DM_Init(nx, ny)

        integer, optional, intent(in) :: nx, ny
        integer :: nproc_x, nproc_y
        integer :: ierr
        integer :: n, i, j
        integer :: errorcode
        integer, dimension(1:2) :: coords
        integer,dimension(:,:),allocatable      :: procmap_for_mpi, cordmap_for_mpi
        integer :: size_of_procmap, size_of_cordmap

        ! calculate the number of processes in x and y direction
        nproc_x = 0
        nproc_y = 0
        if( present(nx) )then
            nproc_x = nx
        end if
        if( present(ny) )then
            nproc_y = ny
        end if

        call verifynxny( nproc_x, nproc_y, prcnumb, ierr )

        if( ierr == 0 )then
            call sqrtnxny( nproc_x, nproc_y, prcnumb )
            n = 0
            n = nproc_x
            nproc_x = nproc_y
            nproc_y = n
        end if

        if( nproc_y <= 1 )then
            !南北向进程数应大于1，除非启动串行模式（串行模式下不应调用本函数）
            if(OnMonitor) write(*,*)"the ny is: ",nproc_y,". it's too small."
            call MPI_Abort( local_communicator, errorcode, ierr )
        end if

        ! print *, "DM_Init: nproc_x, nproc_y = ", nproc_x, nproc_y, " prcnumb = ", prcnumb, " OnMonitor = ", OnMonitor
        if (OnMonitor) write(*,'(A, 2I, A, I, A)') "DM_Init: nproc_x, nproc_y = ", nproc_x, nproc_y, " prcnumb = ", prcnumb
        prcnumb_a = nproc_y
        prcnumb_b = nproc_x

        ! calculate the number of processes in x and y direction, DONE

        ! calculate procmap and cordmap
        allocate( procmap_for_mpi( 0:nproc_y-1, 0:nproc_x-1 ) )
        allocate( cordmap_for_mpi( 1:2, 0:prcnumb-1 ) )

        allocate( procmap( 0:nproc_y-1, 0:nproc_x-1 ) )
        allocate( cordmap( 1:2, 0:prcnumb-1 ) )
        procmap = 0
        procmap_for_mpi = 0
        cordmap = 0
        cordmap_for_mpi = 0

        coords(1) = myprcid/nproc_x
        coords(2) = mod(myprcid,nproc_x)

        !全通信域内，procmap与cordmap的交换，
        !采用MPI_ALLREDUCE函数的原因，是调用简单。
        myprcid_a = coords(1)
        myprcid_b = coords(2)
        procmap_for_mpi( coords(1), coords(2) ) = myprcid
        cordmap_for_mpi( 1, myprcid ) = coords(1)
        cordmap_for_mpi( 2, myprcid ) = coords(2)
        size_of_procmap = nproc_y*nproc_x
        size_of_cordmap = 2*prcnumb
        CALL MPI_ALLREDUCE( procmap_for_mpi, procmap, size_of_procmap, MPI_INTEGER, MPI_SUM, local_communicator, IERR )
        CALL MPI_ALLREDUCE( cordmap_for_mpi, cordmap, size_of_cordmap, MPI_INTEGER, MPI_SUM, local_communicator, IERR )
        
    ! block
    !     integer :: i, j
    ! if (OnMonitor) then
    !     print *, 'procmap:'
    !     do i = 0, nproc_y-1
    !         print *, procmap(i, :)
    !     end do

    !     print *, 'cordmap:'
    !     do i = 0, prcnumb-1
    !         print *, cordmap(:, i)
    !     end do
    ! endif
    ! end block

        deallocate( procmap_for_mpi )
        deallocate( cordmap_for_mpi )

        ! calculate procmap and cordmap, DONE

        ! build communicator for NSWE

        !***********************************************************************************
        !生成用于东西向通信或南北向通信的通信域
        allocate( nsmember      ( 0:nproc_y-1 ) )
        allocate( ewmember      ( 0:nproc_x-1 ) )
        allocate( nscomm_lat    ( 0:nproc_x-1 ) )
        allocate( ewcomm_lon    ( 0:nproc_y-1 ) )
        allocate( nsgroup_lat   ( 0:nproc_x-1 ) )
        allocate( ewgroup_lon   ( 0:nproc_y-1 ) )

        call mpi_comm_group( local_communicator, world_group, ierr )

        do j=0,nproc_y-1
            do i=0,nproc_x-1
                ewmember(i) = procmap(j,i)
            end do
            call mpi_group_incl( world_group, nproc_x, ewmember, ewgroup_lon(j), ierr )
            call mpi_comm_create( local_communicator, ewgroup_lon(j), ewcomm_lon(j), ierr )
        end do
        ewgroup = ewgroup_lon(coords(1))
        ewcomm = ewcomm_lon(coords(1))

        do j=0,nproc_x-1
            do i=0,nproc_y-1
                nsmember(i) = procmap(i,j)
            end do
            call mpi_group_incl( world_group, nproc_y, nsmember, nsgroup_lat(j),ierr )
            call mpi_comm_create( local_communicator, nsgroup_lat(j), nscomm_lat(j), ierr )
        end do
        nsgroup = nsgroup_lat(coords(2))
        nscomm = nscomm_lat(coords(2))

        call mpi_comm_size( ewcomm, ewcommsize, ierr )
        call mpi_comm_size( nscomm, nscommsize, ierr )
        !生成用于东西向通信或南北向通信的通信域
        !***********************************************************************************

        ! build communicator for NSWE, DONE

    end subroutine DM_Init

    subroutine read_nml(self, nmlfn)

        class(domain), intent(inout) :: self
        character(*) , intent(in)    :: nmlfn
        integer                      :: nens
        integer                      :: nvar
        character(len=30), dimension(:), allocatable :: varn_list
        integer                      :: xs
        integer                      :: xe
        integer                      :: ys
        integer                      :: ye
        integer                      :: zs
        integer                      :: ze

        integer                      :: i

        logical                      :: exist
        integer                      :: status, file_unit
        integer                      :: ierr


        if (OnMonitor) then

            namelist /MESO_FIXED/ nens, nvar, xs, xe, ys, ye, zs, ze
            namelist /MESO_ALLOC/ varn_list

            inquire (file=nmlfn, exist=exist)

            if (.not. exist) then
                print *, __LINE__ , __FILE__
                write (stderr, '("Error: input file ", a, " does not exist")') nmlfn
                stop 30
            end if

            open (action='read', file=nmlfn, iostat=status, newunit=file_unit)
            if (status /= 0) write (stderr, '("Error: open Namelist file failed")') nmlfn
            read (nml=MESO_FIXED, iostat=status, unit=file_unit)
            if (status /= 0) write (stderr, '("Error: invalid Namelist format in MESO_FIXED")')

            self%nens = nens
            self%nvar = nvar
            self%xs = xs
            self%xe = xe
            self%ys = ys
            self%ye = ye
            self%zs = zs
            self%ze = ze

            allocate(self%varn_list(nvar))
            allocate(varn_list(nvar))

            read (nml=MESO_ALLOC, iostat=status, unit=file_unit)
            if (status /= 0) write (stderr, '("Error: invalid Namelist format in MESO_ALLOC")')

            do i = 1, nvar
                self%varn_list(i) = trim(varn_list(i))
                ! print *, i, trim(self%varn_list(i))
            enddo
            
            deallocate(varn_list)

        endif

        call MPI_BCAST ( self%nens, 1, MPI_INTEGER, 0, local_communicator, ierr)
        ! call MPI_BCAST ( self%nvar, 1, MPI_INTEGER, 0, local_communicator, ierr)
        call MPI_BCAST ( self%xs, 1, MPI_INTEGER, 0, local_communicator, ierr)
        call MPI_BCAST ( self%xe, 1, MPI_INTEGER, 0, local_communicator, ierr)
        call MPI_BCAST ( self%ys, 1, MPI_INTEGER, 0, local_communicator, ierr)
        call MPI_BCAST ( self%ye, 1, MPI_INTEGER, 0, local_communicator, ierr)
        call MPI_BCAST ( self%zs, 1, MPI_INTEGER, 0, local_communicator, ierr)
        call MPI_BCAST ( self%ze, 1, MPI_INTEGER, 0, local_communicator, ierr)

        ! print *, "GR: ", myprcid, " read_nml: xs, xe, ys, ye, zs, ze, nens = ", self%xs, self%xe, self%ys, self%ye, self%zs, self%ze, self%nens

    endsubroutine read_nml

    subroutine DM_Alloc(xs, xe, ys, ye, zs, ze)

        integer , intent(in) :: xs, xe, ys, ye, zs, ze

        integer                     :: sd1 , ed1 , sp1 , ep1 , sm1 , em1
        integer                     :: sd2 , ed2 , sp2 , ep2 , sm2 , em2
        integer                     :: sd3 , ed3 , sp3 , ep3 , sm3 , em3

        ! initialization of domain variables
        sd1 = xs
        ed1 = xe
        sd2 = zs
        ed2 = ze
        sd3 = ys
        ed3 = ye
        ! initialization of domain variables, DONE

        ! print *, "GR: ", myprcid, " DM_Alloc: xs, xe, ys, ye, zs, ze = ", xs, xe, ys, ye, zs, ze

        sd2 = sd2-1
        ed2 = ed2+1

        call grapes_patch_domain( sd1 , ed1 , sp1 , ep1 , sm1 , em1 , &     ! z-xpose dims
                                  sd2 , ed2 , sp2 , ep2 , sm2 , em2 , &     ! (standard)
                                  sd3 , ed3 , sp3 , ep3 , sm3 , em3 , &
                                  bdyzone_x  , bdyzone_y , shw_control)
                            
        ! if (OnMonitor) then
        !     write(*,*) "DM_Alloc: sd1, ed1, sd2, ed2, sd3, ed3 = ", sd1, ed1, sd2, ed2, sd3, ed3
        !     write(*,*) "DM_Alloc: sp1, ep1, sp2, ep2, sp3, ep3 = ", sp1, ep1, sp2, ep2, sp3, ep3
        !     write(*,*) "DM_Alloc: sm1, em1, sm2, em2, sm3, em3 = ", sm1, em1, sm2, em2, sm3, em3
            ! end if

        ! call alloc_space_field (    new_grid, domain_id ,                   &
        !                             sd1, ed1, sd2, ed2, sd3, ed3,           &
        !                             sm1,  em1,  sm2,  em2,  sm3,  em3       )
        
        grid%sd31                            = sd1 
        grid%ed31                            = ed1
        grid%sp31                            = sp1 
        grid%ep31                            = ep1 
        grid%sm31                            = sm1 
        grid%em31                            = em1
        grid%sd32                            = sd2 
        grid%ed32                            = ed2
        grid%sp32                            = sp2 
        grid%ep32                            = ep2 
        grid%sm32                            = sm2 
        grid%em32                            = em2
        grid%sd33                            = sd3 
        grid%ed33                            = ed3
        grid%sp33                            = sp3 
        grid%ep33                            = ep3 
        grid%sm33                            = sm3 
        grid%em33                            = em3

    end subroutine DM_Alloc

    subroutine grapes_patch_domain( sd1 , ed1 , sp1 , ep1 , sm1 , em1 , &
                                    sd2 , ed2 , sp2 , ep2 , sm2 , em2 , &
                                    sd3 , ed3 , sp3 , ep3 , sm3 , em3 , &
                                    bdx , bdy , shw                     )
   
        integer, intent(in)   :: sd1 , ed1 , sd2 , ed2 , sd3 , ed3 , bdx , bdy ,shw
        integer, intent(out)  :: sp1  , ep1  , sp2  , ep2  , sp3  , ep3  , &  ! z-xpose (std)
                                 sm1  , em1  , sm2  , em2  , sm3  , em3

        call dm_decomposition_reg( sd1 , ed1+1 , sp1 , ep1 , sm1 , em1 , &
                                   sd2 , ed2 , sp2 , ep2 , sm2 , em2 , &
                                   sd3 , ed3+1 , sp3 , ep3 , sm3 , em3 , &
                                   bdx , bdy , shw                       )

        if(ep3>ed3) then
        ep3=ep3-1
        em3=em3-1
        endif

        if(ep1>ed1) then
        ep1=ep1-1
        em1=em1-1
        endif

    end subroutine grapes_patch_domain

    subroutine dm_decomposition_reg( sd1 , ed1 , sp1 , ep1 , sm1 , em1 , &
                                     sd2 , ed2 , sp2 , ep2 , sm2 , em2 , &
                                     sd3 , ed3 , sp3 , ep3 , sm3 , em3 , &
                                     bdx , bdy , shw                     )

        integer,intent(in)     :: sd1 , ed1 ,sd2 , ed2 ,sd3 , ed3, bdx , bdy, shw
        integer,intent(out)    :: sp1 , ep1 ,sp2 , ep2 ,sp3 , ep3
        integer,intent(out)    :: sm1 , em1 ,sm2 , em2 ,sm3 , em3

        !local vars:
        integer   :: bdwith

        !bdwith = max(shw,bdx,bdy)

        bdwith =shw

        !write(0,*) "for x decomposition*******************d=",sd1 , ed1, bdwith," prcnumb_b=",prcnumb_b
        allocate(istedp(2,0:(prcnumb_b-1)))
        call promap( prcnumb_b , myprcid_b, bdwith,                &
                      sd1 , ed1 , sp1 , ep1 , sm1 , em1 , istedp )

        !write(0,*) "for y decomposition*******************d=",sd3 , ed3," prcnumb_a=",prcnumb_a
        allocate(jstedp(2,0:(prcnumb_a-1)))
        call promap( prcnumb_a , myprcid_a, bdwith,                &
                      sd3 , ed3 , sp3 , ep3 , sm3 , em3 , jstedp )

        !write(0,*) "for z decomposition*************************d=",sd2,ed2
        sp2 = sd2+1
        ep2 = ed2-1
        !sp2 = sd2
        !ep2 = ed2
        sm2 = sd2
        em2 = ed2
        !write(0,*) "for z decomposition*************************p/m=",sp2,ep2,sm2,em2

    end subroutine dm_decomposition_reg

    subroutine promap(nprcnum , mypid, bdw,                  &
                    sd , ed , sp , ep , sm , em ,   &
                    ijstedp )

        integer ,intent(in)    :: nprcnum, mypid, sd , ed , bdw
        integer,intent(out)    :: sp , ep , sm , em
        integer,dimension(2,0:nprcnum-1),intent(out)    :: ijstedp

        integer :: ix
        integer :: n, np, nm, nr

!*************************************************************

        ix = ed-sd+1
        if(nprcnum > 2) then               !-----------------------------------

            np = ix / nprcnum              !org
            if( np < 0) then
               write(0,*) "patch on altitude is too small,please set a new cpu number !"
            endif

!**************
            ijstedp(1,0) = 1
            ijstedp(2,0) = np   !org ijstedp(2,0) = np

            nr = ix - np*nprcnum !org

            nm = nprcnum-nr                  ! "-" number processor
            do n = 1, nm/2
               ijstedp(1,n) = ijstedp(2,n-1) + 1
               ijstedp(2,n) = ijstedp(2,n-1) + np
            enddo
            do n = 1+nm/2, nm/2+nr               ! nm/2+1 is the begin position for "+" processor
               ijstedp(1,n) = ijstedp(2,n-1) + 1
               ijstedp(2,n) = ijstedp(2,n-1) + np + 1
            enddo
            do n = 1+nm/2+nr, nprcnum-1
               ijstedp(1,n) = ijstedp(2,n-1) + 1
               ijstedp(2,n) = ijstedp(2,n-1) + np
            enddo
            if(ijstedp(2,nprcnum-1) /= ix) then
                 write(0,*) 'in promap, ijstedp(2,nprcnum-1)=',ijstedp(2,nprcnum-1),ix
                 write(0,*) "something error in promap!!!"
            endif

        else                                 !-----------------------------------
!**************
            np = (ix+1) / nprcnum        
            if( nprcnum ==1) then
               ijstedp(1, 0) = 1
               ijstedp(2, 0) = ix
            else
               ijstedp(1, 0) = 1
               ijstedp(2, 0) = np
               ijstedp(1, 1) = np+1
               ijstedp(2, 1) = ix
            endif
        endif                                !-----------------------------------

        sp = sd+ijstedp(1, mypid)-1
        ep = sd+ijstedp(2, mypid)-1
        sm = sp - bdw - 1
        em = ep + bdw + 1

    end subroutine promap

    subroutine DM_Finalize()

        if ( allocated(grid%varn_list) ) then
            deallocate(grid%varn_list)
        end if

        if( allocated(procmap) )then
            deallocate( procmap )
        end if
        if( allocated(cordmap) )then
            deallocate( cordmap )
        end if
        if( allocated(nsmember) )then
            deallocate( nsmember )
        end if
        if( allocated(ewmember) )then
            deallocate( ewmember )
        end if
        if( allocated(nscomm_lat) )then
            deallocate( nscomm_lat )
        end if
        if( allocated(ewcomm_lon) )then
            deallocate( ewcomm_lon )
        end if
        if( allocated(nsgroup_lat) )then
            deallocate( nsgroup_lat )
        end if
        if( allocated(ewgroup_lon) )then
            deallocate( ewgroup_lon )
        end if
        if( allocated(istedp) )then
            deallocate( istedp )
        end if
        if( allocated(jstedp) )then
            deallocate( jstedp )
        end if

    end subroutine DM_Finalize

    subroutine get_ijk_from_grid (  input_grid ,                     &
                                    ids, ide, jds, jde, kds, kde,    &
                                    ims, ime, jms, jme, kms, kme,    &
                                    ips, ipe, jps, jpe, kps, kpe     )

        type( domain ), intent (in)  :: input_grid
        integer, intent(out) ::                                 &
                               ids, ide, jds, jde, kds, kde,    &
                               ims, ime, jms, jme, kms, kme,    &
                               ips, ipe, jps, jpe, kps, kpe

            ids             = input_grid%sd31 
            ide             = input_grid%ed31 
            jds             = input_grid%sd33 
            jde             = input_grid%ed33 
!             kds             = input_grid%sd32
            kds             = input_grid%sd32+1 
            kde             = input_grid%ed32 
            ims             = input_grid%sm31 
            ime             = input_grid%em31 
            jms             = input_grid%sm33 
            jme             = input_grid%em33 
            kms             = input_grid%sm32
            kme             = input_grid%em32
            ips             = input_grid%sp31 
            ipe             = input_grid%ep31 
            jps             = input_grid%sp33 
            jpe             = input_grid%ep33 
            kps             = input_grid%sp32 
            kpe             = input_grid%ep32 


    end subroutine get_ijk_from_grid

    subroutine verifynxny( nproc_x, nproc_y, ntasks, ierr )

        integer,intent(inout)       :: nproc_x, nproc_y
        integer,intent(in)          :: ntasks
        integer,intent(out)         :: ierr

        ! print *, nproc_x, nproc_y, ierr, ntasks, OnMonitor

        ierr = 0
        if( nproc_x > 0 .and. nproc_y > 0 .and. nproc_x*nproc_y == ntasks)then
            ierr = 1
        else if( nproc_x == 0 .and. nproc_y == 0 )then
            if(OnMonitor) write(*,*)"the nx ny will be calculated by dm_initialize。"
            ierr = 0
        else
            if(OnMonitor) write(*,*)"the nproc_x,nproc_y or prcnumb is wrong, ", nproc_x, nproc_y, ntasks
            if(OnMonitor) write(*,*)"the condition is: nproc_x > 0 .and. nproc_y > 0 .and. nproc_x*nproc_y == prcnumb"
            if(OnMonitor) write(*,*)"the nx ny will be calculated by dm_initialize。"
            ierr = 0
        end if

    end subroutine verifynxny

    subroutine sqrtnxny( nproc_x, nproc_y, ntasks )

        integer,intent(in)          :: ntasks
        integer,intent(inout)       :: nproc_x, nproc_y
        integer                     :: m, n, errorcode, ierr

        !将给定的整数ntasks拆分为nproc_x*nproc_y的形式，
        !nproc_x和nproc_y尽量接近，并优先保证nproc_x为偶数。
        m = sqrt( real(ntasks) ) + 1.5
        if( mod(m,2) /= 0 )then
            m = m - 1
        end if

        do n=m,1,-2
            if( mod(ntasks,n) == 0 )then
                nproc_x = n
                nproc_y = ntasks/n
                exit
            end if
            nproc_x = 1
            nproc_y = ntasks
        end do

    end subroutine sqrtnxny
    
end module EnKF_DM