module enkf_obs
   USE kdtree
   USE precision
   USE netcdf
   use, intrinsic :: iso_fortran_env, only : stderr => error_unit

   implicit none
   private

   public :: enkf_cfg, enkf_final, obs_partition, expand

   interface expand
      module procedure expand_1D_int, expand_1D_real, expand_2D_int, expand_2D_real
   end interface expand

!----------------------------------------------
!-----------  obs instrument unit -------------
   type, public :: obs_instrument
      character(len=64)            :: insn               ! instrument name
      integer                      :: nbatch             ! time to read obs
      integer                      :: nobs               ! number of obs passed qc
      integer                      :: nens               ! ensemble size
      integer                      :: nallocated         !
      integer                      :: nobs_without_qc    !
      logical                      :: avail              ! has obs or not
      real           , allocatable :: lon(:)             ! longitude of obs
      real           , allocatable :: lat(:)             ! latitude of obs
      real           , allocatable :: alt(:)             ! altitude (can be geometric height, geopotential height, or pressure) of obs (not used now)
      integer        , allocatable :: oid(:)             ! obs ID (not used for DA)
      real           , allocatable :: err(:)             ! error of obs
      integer        , allocatable :: qc(:,:)            ! qc flag
      real           , allocatable :: hx(:,:)            ! HX deviation
      real           , allocatable :: inno(:)            ! innovation mean
      type(kd_root)                :: obs_tree           ! kdtree built in this instrument

      contains
         procedure :: obs_init        => enkf_obs_init
         procedure :: obs_read        => enkf_obs_read
         procedure :: obs_final       => enkf_obs_final
         procedure :: obs_qc          => enkf_obs_qc
         procedure :: obs_write       => enkf_obs_write
         procedure :: kdtree_init     => enkf_obs_kdinit
         procedure :: kdtree_final    => enkf_obs_kdfinal
   endtype obs_instrument
!----------------------------------------------


contains

   subroutine enkf_cfg(nmlfn, obs, nins, batch, nens, radius, infl, infl_rtpp, verbose, write_y, mype)
      character(*)                     , intent(in)  :: nmlfn
      type(obs_instrument), allocatable, intent(out) :: obs(:)
      integer                          , intent(out) :: nins
      integer                          , intent(out) :: batch
      integer                          , intent(out) :: nens
      real                             , intent(out) :: radius
      real                             , intent(out) :: infl
      real                             , intent(out) :: infl_rtpp
      logical                          , intent(out) :: verbose
      logical                          , intent(out) :: write_y
      integer                          , intent(in)  :: mype


      character(len=64), allocatable                 :: insn_list(:)
      integer                                        :: status, file_unit
      integer                                        :: i
      logical                                        :: exist

      namelist /STATE_INFO/ nins, batch, nens, radius, infl, infl_rtpp, verbose, write_y
      namelist /ALLOC_INFO/ insn_list

      inquire (file=nmlfn, exist=exist)

      if (.not. exist) then
          print *, __LINE__ , __FILE__
          write (stderr, '("Error: input file ", a, " does not exist")') nmlfn
          stop 30
      end if

      open (action='read', file=nmlfn, iostat=status, newunit=file_unit)
      if (status /= 0) write (stderr, '("Error: open Namelist file failed")') nmlfn
      read (nml=STATE_INFO, iostat=status, unit=file_unit)
      if (status /= 0) write (stderr, '("Error: invalid Namelist format in STATE_INFO")')

      allocate(insn_list(nins))
      allocate(obs      (nins))

      open (action='read', file=nmlfn, iostat=status, newunit=file_unit)
      if (status /= 0) write (stderr, '("Error: open Namelist file failed")') nmlfn
      read (nml=ALLOC_INFO, iostat=status, unit=file_unit)
      if (status /= 0) write (stderr, '("Error: invalid Namelist format in ALLOC_INFO")')

      do i = 1, nins
         obs(i)%insn = insn_list(i)
         obs(i)%nens = nens
      enddo

      deallocate(insn_list)
   endsubroutine enkf_cfg

   subroutine enkf_obs_init(self,batch_size)
      class(obs_instrument), intent(inout) :: self
      integer              , intent(in)    :: batch_size

      integer                              :: i, ens
      character(len=256)                   :: ncfn
      integer                              :: status
      integer                              :: ncid, dimid, varid
      logical                              :: exist

      ens = 1
      write(ncfn,"('obs/',A,I0.3)") trim(self%insn)//".", ens
      inquire (file=ncfn, exist=exist)

      if (.not. exist) then
        print *, __LINE__ , __FILE__
        write (stderr, '("Error: input file ", a, " does not exist")') ncfn
        stop 30
      end if

      call check( nf90_open(trim(ncfn), nf90_nowrite, ncid), 'open '//trim(ncfn))
      call check( nf90_inq_dimid(ncid, "nobs", dimid), 'reading dimension "obs"')
      call check( nf90_inquire_dimension(ncid, dimid, len=self%nobs_without_qc) )
      call check(nf90_close(ncid))

      if (self%nobs_without_qc <= 0) then
         print *, "nobs_without_qc <= 0 for instrument "//trim(self%insn)//" for ncfn = "//trim(ncfn)
         self%avail = .false.
      else
         self%avail  = .true.
         self%nobs   = 0
         self%nbatch = ceiling(real(self%nobs_without_qc) / real(batch_size))
         if (self%nobs_without_qc > batch_size) then
            self%nallocated = batch_size
         else
            self%nallocated = self%nobs_without_qc
         endif
         allocate(self%lon  (self%nallocated))
         allocate(self%lat  (self%nallocated))
         allocate(self%alt  (self%nallocated))
         allocate(self%oid  (self%nallocated))
         allocate(self%err  (self%nallocated))
         allocate(self%qc   (self%nallocated,self%nens))
         allocate(self%hx   (self%nallocated,self%nens))
         allocate(self%inno (self%nallocated))
      endif
   end subroutine enkf_obs_init

   subroutine enkf_obs_read(self,begin_read,size_read,batch_size,step)
      class(obs_instrument), intent(inout) :: self
      integer              , intent(in)    :: begin_read, size_read
      integer              , intent(in)    :: batch_size
      integer              , intent(in)    :: step
      integer                              :: nstep

      integer                              :: i, ens
      character(len=256)                   :: ncfn
      integer                              :: status
      integer                              :: ncid, dimid, varid

      real, allocatable                    :: tmpr(:)
      integer, allocatable                 :: tmpi(:)

      logical                              :: exist

      nstep = self%nbatch

      if (self%nobs+size_read > self%nallocated) then

        if (step < nstep) then
            self%nallocated = ceiling(real(self%nallocated) * real(nstep) / real(step))
        else
            self%nallocated = self%nallocated + batch_size
        endif

         call expand(self%lon   ,self%nallocated)
         call expand(self%lat   ,self%nallocated)
         call expand(self%alt   ,self%nallocated)
         call expand(self%oid   ,self%nallocated)
         call expand(self%err   ,self%nallocated)
         call expand(self%qc    ,self%nallocated)
         call expand(self%hx    ,self%nallocated)
         call expand(self%inno  ,self%nallocated)
      endif

      allocate(tmpr(size_read), tmpi(size_read))

      do ens = 1, self%nens
         write(ncfn,"('obs/',A,I0.3)") trim(self%insn)//".", ens
         inquire (file=ncfn, exist=exist)

         if (.not. exist) then
            print *, __LINE__ , __FILE__
            write (stderr, '("Error: input file ", a, " does not exist")') ncfn
            stop 30
         end if

         call check( nf90_open(trim(ncfn), nf90_nowrite, ncid), 'open '//trim(ncfn))

         if (ens == 1) then
            ! lon
            call check( nf90_inq_varid(ncid, "lon", varid),'read lon')
            call check( nf90_get_var  (ncid, varid, tmpr, start=(/begin_read/),count=(/size_read/)),'get lon')
            do i = 1, size_read
               self%lon(self%nobs+i) = tmpr(i)
            enddo

            ! lat
            call check( nf90_inq_varid(ncid, "lat", varid),'read lat')
            call check( nf90_get_var  (ncid, varid, tmpr, start=(/begin_read/),count=(/size_read/)),'get lat')
            do i = 1, size_read
               self%lat(self%nobs+i) = tmpr(i)
            enddo

            ! alt
            call check( nf90_inq_varid(ncid, "alt", varid),'read alt')
            call check( nf90_get_var  (ncid, varid, tmpr, start=(/begin_read/),count=(/size_read/)),'get alt')
            do i = 1, size_read
               self%alt(self%nobs+i) = tmpr(i)
            enddo

            ! oid
            call check( nf90_inq_varid(ncid, "oid", varid),'read oid')
            call check( nf90_get_var  (ncid, varid, tmpi, start=(/begin_read/),count=(/size_read/)),'get oid')
            do i = 1, size_read
               self%oid(self%nobs+i) = tmpi(i)
            enddo

            ! err
            call check( nf90_inq_varid(ncid, "err", varid),'read err')
            call check( nf90_get_var  (ncid, varid, tmpr, start=(/begin_read/),count=(/size_read/)),'get err')
            do i = 1, size_read
               self%err(self%nobs+i) = tmpr(i)
            enddo
         endif

         ! omb
         call check( nf90_inq_varid(ncid, "omb", varid),'read omb')
         call check( nf90_get_var  (ncid, varid, tmpr, start=(/begin_read/),count=(/size_read/)),'get omb')
         do i = 1, size_read
            self%hx(self%nobs+i,ens) = tmpr(i)
         enddo

         ! qc
         call check( nf90_inq_varid(ncid, "qc", varid),'read qc')
         call check( nf90_get_var  (ncid, varid, tmpi, start=(/begin_read/),count=(/size_read/)),'get qc')
         do i = 1, size_read
            self%qc(self%nobs+i,ens) = tmpi(i)
         enddo

         call check(nf90_close(ncid))
      enddo

      deallocate(tmpr, tmpi)
   endsubroutine enkf_obs_read

   subroutine enkf_obs_qc(self,pe_rank,size_read)
      class(obs_instrument), intent(inout) :: self
      integer              , intent(in)    :: pe_rank
      integer              , intent(in)    :: size_read

      integer                              :: i, j, k
      integer                              :: nens

      if (self%avail) then
         i = self%nobs
         do j = self%nobs+1, self%nobs+size_read
            nens = 0
            do k = 1, self%nens
               if (self%qc(j,k) /= 0) exit
               nens = nens + 1
            enddo
            if (nens == self%nens) then
               i = i + 1
               self%lon(i)   = self%lon(j)
               self%lat(i)   = self%lat(j)
               self%alt(i)   = self%alt(j)
               self%oid(i)   = self%oid(j)
               self%err(i)   = self%err(j) * self%err(j)            ! R
               self%qc(i,:)  = self%qc(j,:)
               self%inno(i)  = sum(self%hx(j,:)) / real(nens)             ! mean(Hx)
               self%hx(i,:)  = self%inno(i) - self%hx(j,:)          ! Hx - mean(Hx)

            endif
         enddo
         self%nobs = i
      endif
   endsubroutine enkf_obs_qc

   subroutine enkf_obs_write(self,pe_rank)
      class(obs_instrument), intent(in) :: self
      integer              , intent(in) :: pe_rank
      integer                           :: i
      character(len=256)                :: ncfn
      real  , allocatable               :: tmp_r(:,:)
      integer , allocatable             :: tmp_i(:)
      integer   ::   ncid, dimid_ns, dimid_ens
      integer   ::   varid_lon, varid_lat, varid_alt, varid_oid, varid_R, varid_omb

      if (.not. self%avail) return

      allocate(tmp_r(self%nobs,self%nens))
      allocate(tmp_i(self%nobs))

      write(ncfn,"(A,'_',I0.4)") trim(self%insn)//"_qc", pe_rank
      ! print *, "obsfn = "//trim(ncfn)

      call check( nf90_create(trim(ncfn), nf90_64bit_offset, ncid), 'write '//trim(ncfn))
      call check( nf90_def_dim  (ncid, "nobs", self%nobs,dimid_ns))
      call check( nf90_def_dim  (ncid, "nens", self%nens,dimid_ens))

      call check(nf90_def_var(ncid,'lon',nf90_real,(/dimid_ns/),varid_lon))
      call check(nf90_def_var(ncid,'lat',nf90_real,(/dimid_ns/),varid_lat))
      call check(nf90_def_var(ncid,'alt',nf90_real,(/dimid_ns/),varid_alt))
      call check(nf90_def_var(ncid,'oid',nf90_int,(/dimid_ns/),varid_oid))
      call check(nf90_def_var(ncid,'R'  ,nf90_real,(/dimid_ns/),varid_R))
      call check(nf90_def_var(ncid,'omb',nf90_real,(/dimid_ns,dimid_ens/),varid_omb))

      call check(nf90_enddef(ncid))

      do i=1,self%nobs
         tmp_r(i,1) = self%lon(i)
      enddo
      call check(nf90_put_var(ncid,varid_lon,tmp_r(:,1)))

      do i=1,self%nobs
         tmp_r(i,1) = self%lat(i)
      enddo
      call check(nf90_put_var(ncid,varid_lat,tmp_r(:,1)))

      do i=1,self%nobs
         tmp_r(i,1) = self%alt(i)
      enddo
      call check(nf90_put_var(ncid,varid_alt,tmp_r(:,1)))

      do i=1,self%nobs
         tmp_i(i) = self%oid(i)
      enddo
      call check(nf90_put_var(ncid,varid_oid,tmp_i(:)))

      do i=1,self%nobs
         tmp_r(i,1) = self%err(i)
      enddo
      call check(nf90_put_var(ncid,varid_R,tmp_r(:,1)))

      do i=1,self%nobs
         tmp_r(i,:) = self%inno(i) - self%hx(i,:)
      enddo
      call check(nf90_put_var(ncid,varid_omb,tmp_r))

      call check(nf90_close(ncid))

      deallocate(tmp_i)
      deallocate(tmp_r)

   endsubroutine enkf_obs_write

   subroutine enkf_obs_kdinit(self)
      class(obs_instrument), intent(inout) :: self
      integer                              :: nobs_tmp

      nobs_tmp = self%nobs

      if (self%avail) then
         call kd_init(self%obs_tree,self%lon(1:nobs_tmp),self%lat(1:nobs_tmp))
      endif

   endsubroutine enkf_obs_kdinit

   subroutine enkf_obs_kdfinal(self)
      class(obs_instrument), intent(inout) :: self

      if (self%avail) then
         call kd_free(self%obs_tree)
      endif

   endsubroutine enkf_obs_kdfinal

   subroutine enkf_obs_final(self)
      class(obs_instrument), intent(inout) :: self

      if (allocated(self%lon)  ) deallocate(self%lon  )
      if (allocated(self%lat)  ) deallocate(self%lat  )
      if (allocated(self%alt)  ) deallocate(self%alt  )
      if (allocated(self%oid)  ) deallocate(self%oid  )
      if (allocated(self%err)  ) deallocate(self%err  )
      if (allocated(self%qc)   ) deallocate(self%qc   )
      if (allocated(self%hx)   ) deallocate(self%hx   )
      if (allocated(self%inno) ) deallocate(self%inno )

   endsubroutine enkf_obs_final

   subroutine enkf_final(obs)
      type(obs_instrument), allocatable, intent(inout):: obs(:)

      if(allocated(obs)) deallocate(obs)
   endsubroutine enkf_final

   subroutine obs_partition(obs,bkg_tree,radius,pe_rank,size_read)
      type(obs_instrument), intent(inout)   :: obs
      type(kd_root)       , intent(in)      :: bkg_tree
      real                , intent(in)      :: radius
      integer             , intent(in)      :: pe_rank
      integer             , intent(inout)   :: size_read

      integer     ::  i, j, k, m
      logical     ::  exist

      m = 0
      i = obs%nobs
      j = obs%nobs
      do k = 1, size_read
         j = j + 1
         call kd_search_radius_TF(bkg_tree, obs%lon(j), obs%lat(j), radius, exist, .FALSE.)
         if (.not. exist) cycle

         m = m + 1
         i = i + 1
         obs%lon  (i)   = obs%lon  (j)
         obs%lat  (i)   = obs%lat  (j)
         obs%alt  (i)   = obs%alt  (j)
         obs%oid  (i)   = obs%oid  (j)
         obs%err  (i)   = obs%err  (j)
         obs%qc   (i,:) = obs%qc   (j,:)
         obs%hx   (i,:) = obs%hx   (j,:)
        !obs%inno (i)   = obs%inno (j)
      enddo
      size_read = m
   endsubroutine obs_partition

   subroutine check(status, str)
      integer, intent(in) :: status
      character(*), optional, intent(in) :: str

      if(status /= 0) then
         print *, status
         if (present(str)) then
            write (*,*) trim(nf90_strerror(status)), ": ",str
         else
            write (*,*) trim(nf90_strerror(status))
         end if
         stop 30 ! todo
      end if
   end subroutine check

   subroutine expand_1D_int(v,n)
      integer, dimension(:), allocatable, intent(inout) :: v
      integer                        , intent(in)       :: n

      integer                                           :: bn
      integer, dimension(:), allocatable                :: vb

      bn = size(v)
      allocate(vb(bn))
      vb = v
      if (allocated(v)) deallocate(v)
      allocate(v(n))
      if (n<bn) then
         v(1:n) = vb(1:n)
      else
         v(1:bn) = vb(1:bn)
      endif
      deallocate(vb)
   endsubroutine expand_1D_int

   subroutine expand_2D_int(v,n)
      integer, dimension(:,:), allocatable, intent(inout) :: v
      integer                        , intent(in)         :: n
      integer                                             :: bn1, bn2
      integer, dimension(:,:), allocatable                :: vb

      bn1 = size(v,dim=1)
      bn2 = size(v,dim=2)
      allocate(vb(bn1,bn2))
      vb = v
      if (allocated(v)) deallocate(v)
      allocate(v(n,bn2))
      if (n<bn1) then
         v(1:n,:) = vb(1:n,:)
      else
         v(1:bn1,:) = vb(1:bn1,:)
      endif
      deallocate(vb)
   endsubroutine expand_2D_int

   subroutine expand_1D_real(v,n)
      real, dimension(:), allocatable, intent(inout) :: v
      integer                        , intent(in)    :: n
      integer                                        :: bn
      real, dimension(:), allocatable                :: vb

      bn = size(v)
      allocate(vb(bn))
      vb = v
      if (allocated(v)) deallocate(v)
      allocate(v(n))
      if (n<bn) then
         v(1:n) = vb(1:n)
      else
         v(1:bn) = vb(1:bn)
      endif
      deallocate(vb)
   endsubroutine expand_1D_real

   subroutine expand_2D_real(v,n)
      real, dimension(:,:), allocatable, intent(inout) :: v
      integer                        , intent(in)      :: n
      integer                                          :: bn1, bn2
      real, dimension(:,:), allocatable                :: vb

      bn1 = size(v,dim=1)
      bn2 = size(v,dim=2)
      allocate(vb(bn1,bn2))
      vb = v
      if (allocated(v)) deallocate(v)
      allocate(v(n,bn2))
      if (n<bn1) then
         v(1:n,:) = vb(1:n,:)
      else
         v(1:bn1,:) = vb(1:bn1,:)
      endif
      deallocate(vb)
   endsubroutine expand_2D_real

end module enkf_obs
