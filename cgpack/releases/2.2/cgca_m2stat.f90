!$Id: cgca_m2stat.f90 14 2014-12-01 10:14:16Z mexas $

!*robodoc*m* CGPACK/cgca_m2stat
!  NAME
!    cgca_m2stat
!  SYNOPSIS
!    module cgca_m2stat
!  DESCRIPTION
!    Module dealing with statistics
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See CGPACK_Copyright
!  CONTAINS
!    cgca_gv, cgca_gvl
!  USES
!    cgca_m1co
!  USED BY
!    cgca
!  SOURCE

module cgca_m2stat
use cgca_m1co
implicit none

private
public :: cgca_gv, cgca_gvl

contains

!*roboend*


!*robodoc*s* cgca_m2stat/cgca_gv
!  NAME
!    cgca_gv
!  SYNOPSIS

subroutine cgca_gv(coarray,gv)

!  INPUTS

integer(kind=iarr),allocatable,intent(in) :: coarray(:,:,:,:)[:,:,:]
integer(kind=ilrg),allocatable,intent(inout) :: gv(:)[:,:,:]

!  SIDE EFFECTS
!    The state of gv array changes
!  DESCRIPTION
!    This routine does grain volume calculation.
!    For each cell (i,j,k) in coarray, the routine increments
!    gv(coarray(i,j,k)).
!  NOTES
!    All images must call this routine!
!
!    There are several SYNC ALL barriers, because all images must
!    get the updated gv array.
!    It is possible (probable?) that there's too much syncronisation,
!    leading to poor performance.
!    This should be investigated at depth.
!  SOURCE

integer(kind=ilrg),allocatable :: gvimg1(:)
integer(kind=ilrg) :: imagevol
integer :: errstat, i1, i2, i3, &
  lbr(4)      ,& ! lower bounds of the "real" coarray, lower virtual+1
  ubr(4)      ,& ! uppper bounds of the "real" coarray, lower virtual+1
  lcob_coar(3),& ! lower cobounds of the coarray
  ucob_coar(3),& ! upper cobounds of the coarray
  lcob_gv(3)  ,& ! lower cobounds of gv
  ucob_gv(3)  ,& ! upper cobounds of gv
  nimages        ! to store num_images() output
logical(kind=ldef) :: image1

!**********************************************************************73
! checks
!**********************************************************************73

if (.not. allocated(coarray)) then
  write (*,'(a,i0)') "ERROR in cgca_gv: image: ", this_image()
  write (*,'(a)') "ERROR in cgca_gv: coarray is not allocated"
  error stop
end if

if (.not. allocated(gv)) then
  write (*,'(a,i0)') "ERROR in cgca_gv: image: ", this_image()
  write (*,'(a)') "ERROR in cgca_gv: gv is not allocated"
  error stop
end if

! make sure coarray and gv have the same cobounds

lcob_coar=lcobound(coarray)
ucob_coar=ucobound(coarray)
lcob_gv=lcobound(gv)
ucob_gv=ucobound(gv)

if ( any (lcob_coar .ne. lcob_gv .or. ucob_coar .ne. ucob_gv)) then
  write (*,'(a,i0)') "ERROR in cgca_gv: image: ", this_image()
  write (*,'(a)') &
    "ERROR in cgca_gv: codimensions of coarray and gv do not match"
  error stop
end if

!**********************************************************************73
! end of checks
!**********************************************************************73

! initialise few variables
errstat = 0
nimages = num_images()

! set image1
image1 = .false.
if (this_image() .eq. 1) image1 = .true.

! Assume the coarray has halos. Ignore those
lbr=lbound(coarray)+1
ubr=ubound(coarray)-1

! zero gv on every image
gv = 0_ilrg

! each image calculates its gv
do i3=lbr(3),ubr(3)
do i2=lbr(2),ubr(2)
do i1=lbr(1),ubr(1)
 gv(coarray(i1,i2,i3,cgca_state_type_grain)) = &
  gv(coarray(i1,i2,i3,cgca_state_type_grain)) + 1_ilrg
end do
end do
end do

! image volume
imagevol = int( size( coarray(lbr(1):ubr(1), lbr(2):ubr(2), &
 lbr(3):ubr(3), cgca_state_type_grain)), kind=ilrg)

! local check on each image: sum(gv) must equal the coarray volume
if (sum(gv) .ne. imagevol) then
  write (*,'(a,i0)') "ERROR in cgca_gv: image: ", this_image()
  write (*,'(a)') "ERROR in cgca_gv: sum(gv) does not match coarray volume"
  error stop
end if

! cannot proceed further until all images
! finish calculating their volumes
sync all

! image1 adds to its own volume volumes from all other images

if (image1) then

 ! preserve gv from image 1
 allocate( gvimg1(size(gv)), stat=errstat)
 if (errstat .ne. 0) then
  write (*,'(a)') "ERROR in cgca_gv: cannot allocate gvimg1"
  write (*,'(a,i0)') "ERROR in cgca_gv: error code: ", errstat
  error stop
 end if

 gvimg1 = gv

 do i3=lcob_gv(3),ucob_gv(3)
 do i2=lcob_gv(2),ucob_gv(2)
 do i1=lcob_gv(1),ucob_gv(1)
  ! image 1 will be counted twice! So need to subtract its
  ! preserved value, gvimg1, from the total
  gv(:) = gv(:) + gv(:)[i1,i2,i3]
 end do
 end do
 end do

 gv = gv - gvimg1

end if
sync all

! get the global volume from image 1
gv(:) = gv(:)[lcob_gv(1),lcob_gv(2),lcob_gv(3)]

! global check: sum(gv) must equal the model volume
if (sum(gv) .ne. imagevol*nimages) then
  write (*,'(a,i0)') "ERROR in cgca_gv: image: ", this_image()
  write (*,'(2(a,i0))') "ERROR in cgca_gv: sum(gv): ", sum(gv), &
    " does not match model volume: ", imagevol*nimages
  error stop
end if

! sync before leaving 
sync all

end subroutine cgca_gv

!*roboend*


!*robodoc*s* cgca_m2stat/cgca_gvl
!  NAME
!    cgca_gvl
!  SYNOPSIS

subroutine cgca_gvl(coarray,gv)

!  INPUTS

integer(kind=iarr),allocatable,intent(in) :: coarray(:,:,:,:)[:,:,:]
integer(kind=ilrg),allocatable,intent(inout) :: gv(:)[:,:,:]

!  SIDE EFFECTS
!    The state of gv array changes
!  DESCRIPTION
!    This routine does grain volume calculation on every image.
!    For each cell (i,j,k) in coarray, the routine increments
!    gv(coarray(i,j,k)). The intention is that after a call
!    to this routine a collective routine is called, e.g. CO_SUM,
!    to calculate grain volumes across all images.
!  NOTES
!    All images must call this routine!
!  SOURCE

integer(kind=ilrg) :: imagevol
integer :: i1, i2, i3, &
  lbr(4)      ,& ! lower bounds of the "real" coarray, lower virtual+1
  ubr(4)      ,& ! uppper bounds of the "real" coarray, lower virtual+1
  lcob_coar(3),& ! lower cobounds of the coarray
  ucob_coar(3),& ! upper cobounds of the coarray
  lcob_gv(3)  ,& ! lower cobounds of gv
  ucob_gv(3)     ! upper cobounds of gv

!**********************************************************************73
! checks
!**********************************************************************73

if (.not. allocated(coarray)) then
  write (*,'(a,i0)') "ERROR in cgca_gvl: image: ", this_image()
  write (*,'(a)') "ERROR in cgca_gvl: coarray is not allocated"
  error stop
end if

if (.not. allocated(gv)) then
  write (*,'(a,i0)') "ERROR in cgca_gvl: image: ", this_image()
  write (*,'(a)') "ERROR in cgca_gvl: gv is not allocated"
  error stop
end if

! make sure coarray and gv have the same cobounds

lcob_coar=lcobound(coarray)
ucob_coar=ucobound(coarray)
lcob_gv=lcobound(gv)
ucob_gv=ucobound(gv)

if ( any (lcob_coar .ne. lcob_gv .or. ucob_coar .ne. ucob_gv)) then
  write (*,'(a,i0)') "ERROR in cgca_gvl: image: ", this_image()
  write (*,'(a)') &
    "ERROR in cgca_gvl: codimensions of coarray and gv do not match"
  error stop
end if

!**********************************************************************73
! end of checks
!**********************************************************************73

! Assume the coarray has halos. Ignore those
lbr=lbound(coarray)+1
ubr=ubound(coarray)-1

! zero gv
gv = 0_ilrg

! each image calculates its gv
do i3=lbr(3),ubr(3)
do i2=lbr(2),ubr(2)
do i1=lbr(1),ubr(1)
 gv(coarray(i1,i2,i3,cgca_state_type_grain)) = &
  gv(coarray(i1,i2,i3,cgca_state_type_grain)) + 1_ilrg
end do
end do
end do

! image volume
imagevol = int( size( coarray(lbr(1):ubr(1), lbr(2):ubr(2), &
 lbr(3):ubr(3), cgca_state_type_grain)), kind=ilrg)

! local check on each image: sum(gv) must equal the coarray volume
if (sum(gv) .ne. imagevol) then
  write (*,'(a,i0)') "ERROR in cgca_gvl: image: ", this_image()
  write (*,'(a)') "ERROR in cgca_gvl: sum(gv) does not match coarray volume"
  error stop
end if

end subroutine cgca_gvl

!*roboend*

end module cgca_m2stat
