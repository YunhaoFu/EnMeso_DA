module letkf_driver
  use mpi
  use enkf_obs
  use letkf_core
  use kdtree, only: kd_search_radius
  use precision
  implicit none

  private

  public :: letkf, letkf_loc_gc

contains
   subroutine debug_check(hdxb, dep, rdiag, rloc, bkg_raw, bkg_mean, bkg, tmp, trans)
      real, intent(in) :: hdxb(:,:)
      real, intent(in) :: dep(:)
      real, intent(in) :: rdiag(:)
      real, intent(in) :: rloc(:)
      real, intent(in) :: bkg_raw(:,:)
      real, intent(in) :: bkg_mean(:)
      real, intent(in) :: bkg(:,:)
      real, intent(in) :: tmp(:,:)
      real, intent(in) :: trans(:,:)

      integer :: nobs, nlev, nens
      integer :: i, j, k
      integer :: rank, ierr, fu
      character(len=128)  :: fn

      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

      nens = size(hdxb,1)
      nobs = size(hdxb,2)
      nlev = size(bkg,1)

      fu = 100 + rank
      write(fn, "('check.',I3.3)") rank
      open(fu, file=fn, status='replace')

      write(fu, "(A)") 'hdxb='
      do j=1,nobs
         write(fu, "(I3,': ')", advance='no') j
         do i=1,nens
            write(fu, "(F30.15)", advance='no') hdxb(i,j)
         enddo
         write(fu, *)
      enddo

      write(fu, "(A)") 'dep='
      do j=1,nobs
         write(fu, "(I3,': ',F30.15)") j, dep(j)
      enddo

      write(fu, "(A)") 'rdiag='
      do j=1,nobs
         write(fu, "(I3,': ',F30.15)") j, rdiag(j)
      enddo

      write(fu, "(A)") 'rloc='
      do j=1,nobs
         write(fu, "(I3,': ',F30.15)") j, rloc(j)
      enddo

      write(fu, "(A)") 'bkg_raw='
      do k=1,nlev
         write(fu, "(I3,': ')", advance='no') k
         do i=1,nens
            write(fu, "(F30.15)", advance='no') bkg_raw(k,i)
         enddo
         write(fu, *)
      enddo

      write(fu, "(A)") 'bkg_mean='
      do k=1,nlev
         write(fu, "(I3,': 'F30.15)") k, bkg_mean(k)
      enddo

      write(fu, "(A)") 'bkg='
      do k=1,nlev
         write(fu, "(I3,': ')", advance='no') k
         do i=1,nens
            write(fu, "(F30.15)", advance='no') bkg(k,i)
         enddo
         write(fu, *)
      enddo

      write(fu, "(A)") 'tmp='
      do k=1,nlev
         write(fu, "(I3,': ')", advance='no') k
         do i=1,nens
            write(fu, "(F30.15)", advance='no') tmp(k,i)
         enddo
         write(fu, *)
      enddo

      write(fu, "(A)") 'trans='
      do j=1,nens
         do i=1,nens
            write(fu, "(F30.15)", advance='no') trans(j,i)
         enddo
         write(fu, *)
      enddo

      close(fu)

      call MPI_Finalize(ierr)
      stop
   end subroutine debug_check

   subroutine output(lon, lat, tq, us, vs, k, prefix)
      real, intent(in)             :: lon(:,:,:)
      real, intent(in)             :: lat(:,:,:)
      real, intent(in)             :: tq(:,:,:,:,:,:)
      real, intent(in)             :: us(:,:,:,:,:)
      real, intent(in)             :: vs(:,:,:,:,:)
      integer, intent(in)          :: k
      character(len=*), intent(in) :: prefix

      integer, parameter :: word=4
      integer            :: ntile, nlon, nlat, nlev, nens, nvar
      integer            :: rank, ierr
      integer            :: tile, ens
      character(len=128) :: fname

      ntile = size(tq,1)
      nlon = size(tq,2)
      nlat = size(tq,3)
      nlev = size(tq,4)
      nens = size(tq,5)
      nvar = size(tq,6)

      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

      ! lon
      write(fname, "('lon.P',I3.3)") rank
      open(100, file=fname, access='direct', form='UNFORMATTED', recl=word*nlon*nlat, status='replace')
      do tile=1,ntile
         write(100, rec=tile) lon(tile, :, :)
      enddo
      close(100)
      ! lat
      write(fname, "('lat.P',I3.3)") rank
      open(100, file=fname, access='direct', form='UNFORMATTED', recl=word*nlon*nlat, status='replace')
      do tile=1,ntile
         write(100, rec=tile) lat(tile, :, :)
      enddo
      close(100)

      ! T
      do ens=1,nens
         write(fname, "(A, 'T',I2.2, '.E',I2.2, '.P',I3.3)") prefix, k, ens, rank
         open(100, file=fname, access='direct', form='UNFORMATTED', recl=word*nlon*nlat, status='replace')
         do tile=1,ntile
            write(100, rec=tile) tq(tile, :, :, k, ens, 1)
         enddo
         close(100)
      enddo

      ! Q
      do ens=1,nens
         write(fname, "(A, 'Q',I2.2, '.E',I2.2, '.P',I3.3)") prefix, k, ens, rank
         open(100, file=fname, access='direct', form='UNFORMATTED', recl=word*nlon*nlat, status='replace')
         do tile=1,ntile
            write(100, rec=tile) tq(tile, :, :, k, ens, 2)
         enddo
         close(100)
      enddo

      ! U
      do ens=1,nens
         write(fname, "(A, 'U',I2.2, '.E',I2.2, '.P',I3.3)") prefix, k, ens, rank
         open(100, file=fname, access='direct', form='UNFORMATTED', recl=word*nlon*(nlat+1), status='replace')
         do tile=1,ntile
            write(100, rec=tile) us(tile, :, :, k, ens)
         enddo
         close(100)
      enddo

      ! V
      do ens=1,nens
         write(fname, "(A, 'V',I2.2, '.E',I2.2, '.P',I3.3)") prefix, k, ens, rank
         open(100, file=fname, access='direct', form='UNFORMATTED', recl=word*(nlon+1)*nlat, status='replace')
         do tile=1,ntile
            write(100, rec=tile) vs(tile, :, :, k, ens)
         enddo
         close(100)
      enddo
   end subroutine output

   subroutine letkf(nlon, nlat, lon, lat, nlev, nens, nins, y, radius, infl, var)
      integer , intent(in)             ::    nlon
      integer , intent(in)             ::    nlat
      real    , intent(in)             ::    lon(nlon,nlat)
      real    , intent(in)             ::    lat(nlon,nlat)
      integer , intent(in)             ::    nlev
      integer , intent(in)             ::    nens
      integer , intent(in)             ::    nins
      type(obs_instrument),intent(in)  ::    y(nins)
      real    , intent(in)             ::    radius
      real    , intent(in)             ::    infl
      real    , intent(inout)          ::    var(nlon,nlat,nlev,nens)

      integer                          ::    ins
      integer, allocatable             ::    idx(:)
      real, allocatable                ::    dis(:)

      real, allocatable                ::    hdxb(:,:)
      real, allocatable                ::    dep(:)
      real, allocatable                ::    rdiag(:)
      real, allocatable                ::    dist(:)
      real, allocatable                ::    rloc(:)

      integer, parameter               ::    nobs_limit = 1000000
      integer                          ::    nobs_loc, nobs_max, nobs, ny
      real                             ::    loc_tmp
      real                             ::    cutoff

      real, allocatable                ::    trans(:, :), bkg_mean(:), bkg(:,:), tmp(:,:)

      integer :: i, j, k, e, l

      call letkf_core_init(nens)

      cutoff = radius / (sqrt(40.0/3.0))

      nobs_max = 0
      nobs = 0
      do ins=1,nins
         nobs_max = max(nobs_max, y(ins)%nobs)
         nobs = nobs + y(ins)%nobs
      enddo
      nobs_max = min(nobs_max, nobs_limit)
      nobs = min(nobs, nobs_limit)

      allocate(idx(nobs_max))
      allocate(dis(nobs_max))

      allocate(hdxb(nens,nobs))
      allocate(dep(nobs))
      allocate(rdiag(nobs))
      allocate(dist(nobs))
      allocate(rloc(nobs))

      allocate(trans(nens,nens))
      allocate(bkg_mean(nlev))
      allocate(bkg(nlev,nens))
      allocate(tmp(nlev,nens))

      do j=1,nlat
      do i=1,nlon
         nobs = 0
         do ins=1,nins
            if (.not. y(ins)%avail) cycle
            ny = 0
            call kd_search_radius(y(ins)%obs_tree, lon(i,j), lat(i,j), radius, idx, dis, ny, .FALSE.)

            do l = 1, ny
               if (nobs < nobs_limit) then
                  nobs = nobs + 1
                  hdxb(:,nobs) = y(ins)%hx(idx(l),:)
                  dep(nobs)    = y(ins)%inno(idx(l))
                  rdiag(nobs)  = y(ins)%err(idx(l))
                  dist(nobs)   = dis(l)
               else
                  write(*, "('observations within the radius of ', F20.10, ' km: nobs=', I0, ' nobs_limit=', I0)") radius, nobs, nobs_limit
                  exit
               endif
            enddo
         enddo
         if (nobs <= 0) cycle

         nobs_loc = 0
         do l=1,nobs
            loc_tmp = letkf_loc_gc(dist(l), cutoff)
            if (loc_tmp <= 1.0e-6) cycle
            nobs_loc = nobs_loc + 1
            hdxb(:,nobs_loc)  = hdxb(:,l)
            dep(nobs_loc)     = dep(l)
            rdiag(nobs_loc)   = rdiag(l)
            rloc(nobs_loc)    = loc_tmp
         enddo
         if (nobs_loc <= 0) cycle

         call letkf_core_solve(nobs_loc, hdxb(:,:nobs_loc),       &
                               rdiag(:nobs_loc), rloc(:nobs_loc), &
                               dep(:nobs_loc), infl, trans)

         !*****************************************
         ! variables at grid center such as T and Q
         !*****************************************
         do k=1,nlev
            bkg_mean(k) = sum(var(i,j,k,:)) / nens
            do e=1,nens
               bkg(k,e) = var(i,j,k,e) - bkg_mean(k)
            enddo
         enddo
         call sgemm('n', 'n', nlev, nens, nens, 1.0, bkg, nlev, trans, nens, 0.0, tmp, nlev)

         do e=1,nens
         do k=1,nlev
            var(i,j,k,e) = bkg_mean(k) + tmp(k,e)
         enddo
         enddo

      enddo
      enddo

      deallocate(idx)
      deallocate(dis)

      deallocate(hdxb)
      deallocate(dep)
      deallocate(rdiag)
      deallocate(dist)
      deallocate(rloc)

      deallocate(trans)
      deallocate(bkg_mean)
      deallocate(bkg)
      deallocate(tmp)

   endsubroutine letkf

   !================================================================================
   !> Gaspari-Cohn localization function.
   !! Possibly faster than the Gaussian function, depending on computer architecture.
   !! Similar shape to Gaussian, except it is compact, goes to 0 at 2L*sqrt(0.3)
   !--------------------------------------------------------------------------------
   PURE FUNCTION letkf_loc_gc(z, L) RESULT(res)
      REAL, INTENT(in) :: z  !< value to localize
      REAL, INTENT(in) :: L  !< the equivalent to the Gaussian standard deviation
      REAL :: res
      REAL(8) :: c
      REAL(8) :: abs_z, z_c

      c = L / SQRT(0.3d0)
      abs_z = ABS(z)
      z_c = abs_z / c

      IF (abs_z >= 2*c) THEN
         res = 0.0
      ELSEIF (abs_z > c) THEN
         res = &
               0.08333d0 * z_c**5 &
             - 0.50000d0 * z_c**4 &
             + 0.62500d0 * z_c**3 &
             + 1.66667d0 * z_c**2 &
             - 5.00000d0 * z_c &
             + 4d0 &
             - 0.66667d0 * c/abs_z
      ELSE
         res = &
              -0.25000d0 * z_c**5 &
             + 0.50000d0 * z_c**4 &
             + 0.62500d0 * z_c**3 &
             - 1.66667d0 * z_c**2 &
             + 1d0
      END IF
   END FUNCTION letkf_loc_gc
   !================================================================================

end module letkf_driver
