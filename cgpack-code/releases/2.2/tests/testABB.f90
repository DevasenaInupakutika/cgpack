!$Id: testABB.f90 8 2014-12-01 09:13:54Z mexas $

!*robodoc*u* tests/testABB
!  NAME
!    testABB
!  SYNOPSIS

program testABB

!  PURPOSE
!    Checking: cgca_clvgp, cgca_clvgsd
!  DESCRIPTION
!    Checking cleavage propagation across grain boundary
!    with many grains. Still a single cleavage nucleus.
!    Put it at the centre of one of the faces of the model.
!  NOTE
!    With no grain boundary smoothing, crack
!    finds it very hard to propagate across a grain boundary.
!    This is because the GB is locally very irregular, and
!    it is likely the first cell in the new grain will find
!    itself in some sort of corner or a tunnel, from where
!    it cannot see enough neighbours of the same grains
!    to propagate into. Hence use cgca_gbs, at least once,
!    possibly multiple iterations.
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

logical( kind=ldef ), parameter :: yesdebug=.true., nodebug=.false., &
 periodicbc=.true.,  noperiodicbc=.false.

real, parameter :: gigabyte=real(2**30), resolution=1.0e-5,          &
! cleavage stress on 100, 110, 111 planes for BCC,
! see the manual for derivation.
 scrit(3) = (/ 1.05e4, 1.25e4, 4.90e4 /)

real( kind=rdef ), allocatable :: grt(:,:,:)[:,:,:]
real( kind=rdef ) :: t(3,3)   ! stress tensor

integer( kind=idef ) :: l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,  &
 cou3, nuc, nimages, codim(3)[*], iter
integer( kind=iarr ), allocatable :: space(:,:,:,:)[:,:,:]
integer( kind=ilrg ) :: icells,mcells

real :: image_storage

logical( kind=ldef ) :: solid, image1
character(6) :: image

!*********************************************************************72
! first executable statement

nimages = num_images()
 image1 = .false.
if (this_image() .eq. 1) image1 = .true.

! do a check on image 1
if (image1) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("ABB")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
end if

sync all

codim(:) = codim(:)[1]

l1 = 1
l2 = l1
l3 = l1

u1 = 2**5 ! 32
u2 = u1
u3 = u1

col1 = 1
cou1 = codim(1) - col1 + 1
col2 = 1
cou2 = codim(2) - col2 + 1
col3 = 1
cou3 = codim(3) - col3 + 1

! total number of cells in a coarray
icells = int( u1-l1+1, kind=ilrg ) *    &
         int( u2-l2+1, kind=ilrg ) *    &
         int( u3-l3+1, kind=ilrg )

! total number of cells in the model
mcells = icells * int( nimages, kind=ilrg )

! total number of nuclei
! nuc should not exceed resolution*mcells
nuc = 20

100 format( a,3(i0,":",i0,a) )

if (image1) then
  write (*,100) "bounds: (", l1,u1,",", l2,u2,",", l3,u3, ")"
  write (*,100) "cobounds: [", &
    col1,cou1,",", col2,cou2,",", col3,cou3, "]"

  ! An absolute minimum of storage, in GB, per image.
  ! A factor of 2 is used because will call _sld, which
  ! allocates another array of the same size and kind as
  ! coarray.
  image_storage = real(2 * icells*storage_size(space)/8)/gigabyte

  write (*,'(a,i0,a)') "Each image has ",icells, " cells"
  write (*,'(a,i0,a)') "The model has ", mcells, " cells"
  write (*,'(a,i0,a)') "The model has ", nuc, " nuclei"
  write (*,'(a,es9.2,a)') "Each image will use at least ",       &
                           image_storage, " GB memory"
end if

! initialise random number seed
call cgca_irs( nodebug )

! allocate space with two layers
call cgca_as( l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,2,space )

! allocate rotation tensors
call cgca_art( 1, nuc, col1, cou1, col2, cou2, col3, grt )

! initialise space
space( :, :, :, cgca_state_type_grain ) = cgca_liquid_state
space( :, :, :, cgca_state_type_frac  ) = cgca_intact_state

! nuclei, sync all inside
call cgca_nr( space, nuc, yesdebug )

! assign rotation tensors, sync all inside
call cgca_rt( grt )

! set a single crack nucleus in the centre of the front face
space( u1/2, u2/2, u3, cgca_state_type_frac )[ cou1/2, cou2/2, cou3] = &
 cgca_clvg_state_100_edge

! solidify, implicit sync all inside
call cgca_sld( space, noperiodicbc, 0, 10, solid )

! smoothen the GB, several iterations, sync needed
do iter=1,2
  call cgca_gbs( space )
  sync all

  ! halo exchange, following smoothing
  call cgca_hxi( space )
  sync all
end do

! update grain connectivity, local routine, no sync needed
call cgca_gcu( space )

! dump grain connectivity to files, local routine, no sync needed
write ( image, "(i0)" ) this_image()
call cgca_gcp( ounit=10, fname="z_gc_1_"//image )

if ( image1 ) write (*,*) "dumping model to files"

! dump space arrays to files, only image 1 does it, all others
! wait at the barrier, hence sync needed
call cgca_swci( space, cgca_state_type_grain, 10, "zg1.raw" )
call cgca_swci( space, cgca_state_type_frac,  10, "zf1.raw" )

if ( image1 ) write (*,*) "finished dumping model to files"

sync all

! set the stress tensor
t(1,1) = 1.0e6
t(2,2) = -1.0e6

! propagate cleavage, sync inside
! subroutine cgca_clvgp( coarray, rt, t, scrit, sub, periodicbc,    &
!                        iter, heartbeat, debug )
call cgca_clvgp( space, grt, t, scrit, cgca_clvgsd, noperiodicbc,   &
                 30, 10, yesdebug )

! dump grain connectivity to files, local routine, no sync needed.
write ( image, "(i0)" ) this_image()
call cgca_gcp( ounit=10, fname="z_gc_2_"//image )

! dump the fracture space array to files, only image 1 does it,
! all others wait at the barrier, hence sync needed
call cgca_swci( space, cgca_state_type_frac, 10, "zf2.raw" )

! However, since there's nothing more to do, no sync is needed.

! deallocate all arrays, implicit sync all.
call cgca_ds( space )
call cgca_dgc
call cgca_drt( grt )

end program testABB

!*roboend*
