!$Id: testABC.f90 8 2014-12-01 09:13:54Z mexas $

!*robodoc*u* tests/testABC
!  NAME
!    testABC
!  SYNOPSIS

program testABC

!  PURPOSE
!    Checking: cgca_pdmp
!  DESCRIPTION
!    Dump the global CGPACK parameters to stdout.
!    Any or all images can call this routine.
!    Here only image 1 does it.
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

integer(kind=idef) :: nimages,codim(3)[*]
logical(kind=ldef) :: image1

!**********************************************************************73
! first executable statement

nimages=num_images()
image1=.false.
if (this_image() .eq. 1) image1 = .true.

! do a check on image 1
if (image1) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("ABC")
 ! print parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "Running on ", nimages, " images in a 3D grid"
end if

end program testABC

!*roboend*
