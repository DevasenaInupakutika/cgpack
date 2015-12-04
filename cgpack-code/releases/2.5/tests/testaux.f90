!$Id: testaux.f90 8 2014-12-01 09:13:54Z mexas $

!*robodoc*m* tests/testaux
!  NAME
!    testaux
!  SYNOPSIS

module testaux

!  DESCRIPTION
!    A helper module to use with the tests.
!    This is not a part of CGCA.
!    Use of this module is not required.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See CGPACK_Copyright
!  CONTAINS
!    getcodim, banner
!  USES
!    cgca
!  USED BY
!    Part of CGPACK test suite
!  SOURCE

use cgca
implicit none

contains

!*roboend*


!*robodoc*s* testaux/getcodim
!  NAME
!    getcodim
!  SYNOPSIS

subroutine getcodim(nimages,codim)

!  INPUT

integer( kind=idef ), intent( in ) :: nimages

!  OUTPUTS

integer( kind=idef ), intent( out ) :: codim(3)

!  DESCRIPTION
!    Subroutine to get coarray codimensions from
!    command line. Call this routine from all tests.
!  USES
!    none
!  USED BY
!    all tests
!  SOURCE

integer, parameter :: maxarglen=20
character( len=maxarglen ) :: value,fmt
integer :: arglen,errstat,i 

errstat = 0

do i = 1, 2
   call get_command_argument (i,value,arglen,errstat)

   if (errstat .gt. 0) then
     write (*,'(a,i0,a)') "ERROR: argument ", i, " cannot be retrieved"
     error stop
   elseif (errstat .eq. -1) then
     write (*,'(a,i0,a)') "ERROR: argument ", i, &
       " length is longer than the string to store it"
     error stop
   elseif (errstat .lt. -1) then
     write (*,'(a,i0)') &
       "ERROR: get_command_argument(", i, " returned ",errstat
     write (*,'(a)') "ERROR: unknown error, should never end up here"
     error stop 
   end if

   write(fmt,'(a,i0,a)') "(i", arglen, ")"
   read (value,fmt) codim(i)
end do

if (codim(1) .le. 0 .or. codim(2) .le. 0) &
  error stop "ERROR: one or both arguments are zero or negative"

codim(3)=nimages/(codim(1)*codim(2))

if (codim(1)*codim(2)*codim(3) .ne. nimages) &
  error stop "ERROR: ended up with codimension 3 not a positive integer"

end subroutine getcodim

!*roboend*


!*robodoc*s* testaux/banner
!  NAME
!    banner
!  SYNOPSIS

subroutine banner( test )

!  INPUT

character( len=3 ), intent( in ) :: test

!  SIDE EFFECTS
!    A banner with the test name is printed on sdtout
!  USES
!    none
!  USED BY
!    all tests
!  SOURCE

write (*,'(a)') "******************************************************"
write (*,'(a)') "******                                          ******"
write (*,'(a)') "******                 TEST "//test//"                 ******"
write (*,'(a)') "******                                          ******"
write (*,'(a)') "******************************************************"

end subroutine banner

!*roboend*

end module testaux
