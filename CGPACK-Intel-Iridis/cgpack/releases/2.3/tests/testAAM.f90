!$Id: testAAM.f90 8 2014-12-01 09:13:54Z mexas $

!*robodoc*u* tests/testAAM
!  NAME
!    testAAM
!  SYNOPSIS

program testAAM

!  PURPOSE
!    Checking: cgca_av, cgca_dv
!  DESCRIPTION
!    Checking grain volume array alloc/deallocation.
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

 logical(kind=ldef),parameter :: yesdebug=.true., nodebug=.false., &
   periodicbc=.true.
 real,parameter :: gigabyte=real(2**30), resolution=1.0e-5

 integer(kind=idef) :: l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3, &
   cou3, &
   nuc,    & ! number of nuclei in the model
   nimages,codim(3)[*]
 integer(kind=iarr),allocatable :: space1(:,:,:,:)[:,:,:]
 integer(kind=ilrg),allocatable :: grainvol(:)[:,:,:]
 integer(kind=ilrg) :: icells,mcells

 logical(kind=ldef) :: image1

 real :: image_storage

!**********************************************************************73
! first executable statement

nimages=num_images()
image1=.false.
if (this_image().eq.1) image1=.true.

! do a check on image 1
if (image1) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("AAM")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
 write (*,*) "codim:", codim
end if

sync all

codim(:) = codim(:)[1]

if (image1) call system("sleep 1")

l1=1
l2=l1
l3=l1

! The array size is only controlled by this value
! in this program.
u1=50
u2=u1
u3=u1

col1=1
cou1=codim(1)-col1+1
col2=1
cou2=codim(2)-col2+1
col3=1
cou3=codim(3)-col3+1

! total number of cells in a coarray
icells = int(u1-l1+1,kind=ilrg) * int(u2-l2+1,kind=ilrg) * &
  int(u3-l3+1,kind=ilrg)

! total number of cells in the model
mcells = icells * int(codim(1),kind=ilrg) * int(codim(2),kind=ilrg) * &
  int(codim(3),kind=ilrg)

! total number of nuclei
nuc = int(resolution*mcells)

if (image1) then
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "bounds: (",l1,u1,l2,u2,l3,u3
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "cobounds: (",col1,cou1,col2,cou2,col3,cou3

  ! An absolute minimum of storage, in GB, per image.
  ! A factor of 2 is used because will call _sld, which
  ! allocates another array of the same size and kind as
  ! coarray.
  image_storage = real(2 * icells*storage_size(space1)/8)/gigabyte 

  write (*,'(a,i0,a)') "Each image has ",icells, " cells"
  write (*,'(a,i0,a)') "The model has ", mcells, " cells"
  write (*,'(a,i0,a)') "The model has ", nuc, " nuclei"
  write (*,'(a,es9.2,a)') "Each image will use at least ", &
    image_storage, " GB memory"
end if

! initialise random number seed
call cgca_irs(nodebug)

! allocate coarray
call cgca_as(l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,1,space1)

! initialise coarray to liquid
space1(:,:,:,cgca_state_type_grain) = cgca_liquid_state

! allocate grain volume
call cgca_av(0,nuc,col1,cou1,col2,cou2,col3,grainvol)

if (image1) write (*,'(a)') "grain volume coarray allocated"

! deallocate all arrays
call cgca_ds(space1)
call cgca_dv(grainvol)
if (image1) write (*,'(a)') "grain volume coarray deallocated"

end program testAAM

!*roboend*
