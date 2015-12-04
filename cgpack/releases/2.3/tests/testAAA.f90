!$Id: testAAA.f90 8 2014-12-01 09:13:54Z mexas $

!*robodoc*u* tests/testAAA
!  NAME
!    testAAA
!  SYNOPSIS

program testAAA

!  PURPOSE
!    Checking: getcodim, cgca_as, cgca_ds 
!  DESCRIPTION
!    Testing allocating and deallocating a coarray
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

!**********************************************************************73
! first executable statement

 nimages=num_images()

! do a check on image 1

! do a check on image 1
if (this_image() .eq. 1) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("AAA")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
 write (*,*) "codim:", codim
end if

sync all

codim(:) = codim(:)[1]

if (this_image() .eq. 2) call system("sleep 1")

call cgca_as(1,size1,1,size2,1,size3,1,codim(1),1,codim(2),1,2,space1)

if (allocated(space1)) then
  write (*,'(a,i0,a)')"Image:",this_image(), " space1 allocated"
  write (*,*) "My array is:"
  write (*,'(a,4(i0,tr1))') "Shape:" , shape(space1)
  write (*,'(a,4(i0,tr1))') "lbound:", lbound(space1)
  write (*,'(a,4(i0,tr1))') "ubound:", ubound(space1)
  write (*,'(a,i0)')        "size:"  , size(space1)
  write (*,'(a,3(i0,tr1))') "lcobound:", lcobound(space1)
  write (*,'(a,3(i0,tr1))') "ucobound:", ucobound(space1)
end if

sync all

write (*,'(a,i0,a,3(i0,tr1),a)') &
  "Image: ",this_image()," is ",this_image(space1)," in the grid"

if (this_image() .eq. 3) call system("sleep 2")

call cgca_ds(space1)
 
if (.not. allocated(space1)) &
  write (*,'(a,i0,a)')"Image:",this_image(), " space1 not allocated"

end program testAAA 

!*roboend*
