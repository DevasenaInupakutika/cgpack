!$Id: testABT.f90 8 2014-12-01 09:13:54Z mexas $

!*robodoc*u* tests/testABT
!  NAME
!    testABT
!  SYNOPSIS

program testABT

!  PURPOSE
!    Checking: cgca_gdim, cgca_cadim
!  DESCRIPTION
!    cgca_gdim finds the optimum coarray grid layout for a
!    given total number of images. It also reports the
!    quality of this optimum, from 0 - worst, to 1 - best.
!    cgca_cadim then calculates the coarray dimensions
!    the new updated box size.
!  NOTE
!    Both cgca_gdim and cgca_cadim are serial routines.
!    It makes no sence to run this test at high numbers of images.
!    A single image is enough to test the routines.
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

integer( kind=idef ) :: n, ir(3), nimgs, lres, &
 ng                       ! number of grains in the whole model
integer( kind=iarr ) :: c(3)  ! coarray dimensions
logical( kind=ldef ) :: image1
real( kind=rdef ) ::    &
 qual,                  & ! quality
 bsz0(3),               & ! the given "box" size
 bsz(3),                & ! updated "box" size
 dm,                    & ! mean grain size, linear dim, phys units
 res                      ! resolutions, cells per grain

!*********************************************************************72
! first executable statement

! physical dimensions of the box, assume mm
bsz0 = (/ 10, 20, 30 /)

! mean grain size, also mm
dm = 1.0e-1

! resolution
res = 1.0e5 

 nimgs = num_images()
image1 = .false.
if (this_image() .eq. 1) image1 = .true.

! do a check on image 1
if (image1) then

 ! print a banner
 call banner("ABT")

 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimgs, " images in a 3D grid"

 ! calculate the coarray grid dimensions
 do n = 1, 2**14
   call cgca_gdim( n, ir, qual )

   ! subroutine cgca_cadim( bsz, res, dm, ir, c, lres, ng )
   bsz = bsz0
   call cgca_cadim( bsz, res, dm, ir, c, lres, ng )

   write ( *, "(8(i0,a),g0,tr1,i0,3(a,g0),a)" )                &
    n, " (", c(1), ",", c(2), ",", c(3), ")[" ,                &
    ir(1), ",", ir(2), ",", ir(3), "] ", ng, " ",              &
    qual, lres,                                                &
    " (", bsz(1), ",", bsz(2), ",", bsz(3), ")"
 end do

end if

sync all

end program testABT


!*roboend*
