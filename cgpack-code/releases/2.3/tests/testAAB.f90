!$Id: testAAB.f90 8 2014-12-01 09:13:54Z mexas $

!*robodoc*u* tests/testAAB
!  NAME
!    testAAB
!  SYNOPSIS

program testAAB

!  PURPOSE
!    Checking: cgca_as, cgca_ds, cgca_swci
!  DESCRIPTION
!    Writing all coarrays to a file in order, as if it were a single
!    large 3D array (super array).
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See CGPACK_Copyright
!  USES
!    cgca testaux
!  USED BY
!    Part of CGPACK test suite
!  SOURCE

use testaux

implicit none

 integer(kind=idef),parameter :: size1=10, size2=10, size3=10
 integer(kind=idef) :: nimages
 integer(kind=idef) :: codim(3)[*]
 integer(kind=iarr),allocatable :: space1(:,:,:,:)[:,:,:]
logical(kind=ldef) :: image1

!**********************************************************************73
! first executable statement

nimages=num_images()
image1 = .false.
if (this_image() .eq. 1) image1 = .true.

! do a check on image 1
if (image1) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("AAB")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
 write (*,*) "codim:", codim
end if

sync all

if (image1) then
  write (*,*) "cgca kinds:"
  write (*,*) "default integer:", idef
  write (*,*) "coarray integer:", iarr
  write (*,*) "default logical:", ldef
end if

codim(:) = codim(:)[1]

if (this_image() .eq. 2) call system("sleep 1")

call cgca_as(1,size1,1,size2,1,size3,1,codim(1),1,codim(2),1,2,space1)

if (allocated(space1)) then
  write (*,'(a,i0,a)')"Image:",this_image(), " space1 allocated"
  write (*,'(a,i0,a,3(i0,tr1),a)') &
    "Image: ",this_image()," is ",this_image(space1)," in the grid"
end if

space1(:,:,:,cgca_state_type_grain) = this_image()
space1(:,:,:,cgca_state_type_frac) = 0
sync all

if (this_image() .eq. 1) write (*,*) "coarrays defined, calling cgca_swci"

call cgca_swci(space1, cgca_state_type_grain, 10, 'z.raw')
sync all

if (this_image() .eq. 3) call system("sleep 2")

call cgca_ds(space1)
 
if (.not. allocated(space1)) &
  write (*,'(a,i0,a)')"Image:",this_image(), " space1 not allocated"

end program testAAB

!*roboend*
