!$Id: cgca_m2pfem.f90 82 2015-03-30 09:09:37Z mexas $

!*********************************************************************72

!*robodoc*m* CGPACK/cgca_m2pfem
!  NAME
!    cgca_m2pfem
!  SYNOPSIS

module cgca_m2pfem

!  DESCRIPTION
!    Module dealing with interfacing CGPACK with ParaFEM.
!  AUTHOR
!    Anton Shterenlikht, Luis Cebamanos
!  COPYRIGHT
!    See CGPACK_Copyright
!  CONTAINS
!    Public coarray variables of derived types:
!    cgca_pfem_centroid_tmp,
!    cgca_pfem_integrity,
!    cgca_pfem_stress.
!
!    Public *local*, non-coarray, variable:
!    cgca_pfem_enew.
!
!    Private *local*, non-coarray, variables:
!    lcentr.
!
!    Routines:
!    cgca_pfem_cenc,
!    cgca_pfem_cendmp,
!    cgca_pfem_ctalloc, cgca_pfem_ctdalloc,
!    cgca_pfem_ealloc, cgca_pfem_edalloc,
!    cgca_pfem_intcalc1,
!    cgca_pfem_integalloc, cgca_pfem_integdalloc,
!    cgca_pfem_salloc, cgca_pfem_sdalloc,
!    cgca_pfem_sdmp,
!    cgca_pfem_simg,
!    cgca_pfem_uym
!  USES
!    cgca_m1co
!  USED BY
!    end user?
!  SOURCE

use :: cgca_m1co, only : iarr, idef, ilrg, rdef
implicit none

private
public :: &
! routines
  cgca_pfem_cenc, cgca_pfem_cendmp, &
  cgca_pfem_ctalloc, cgca_pfem_ctdalloc,        &
  cgca_pfem_ealloc, cgca_pfem_edalloc, &
  cgca_pfem_integalloc, cgca_pfem_integdalloc, &
  cgca_pfem_intcalc1, &
  cgca_pfem_salloc,  cgca_pfem_sdalloc, &
  cgca_pfem_sdmp, &
  cgca_pfem_simg, cgca_pfem_uym, &
! variables 
  cgca_pfem_centroid_tmp, &
  cgca_pfem_enew, &
  cgca_pfem_integrity, &
  cgca_pfem_stress

! corresponds to typical double precision real.
integer, parameter :: cgca_pfem_iwp = selected_real_kind(15,300)

!*roboend*


!*robodoc*d* cgca_m2pfem/lcentr
!  NAME
!    lcentr
!  SYNOPSIS

type mcen
  integer( kind=idef ) :: image
  integer( kind=idef ) :: elnum
  real( kind=cgca_pfem_iwp ) :: centr(3)
end type mcen
type( mcen ), allocatable :: lcentr(:)

!  DESCRIPTION
!    A *private* *local* allocatable array of derived type with 3
!    components: (1) image number (2) the local element number on that
!    image and (3) centroid coordinates. Each entry in this array
!    corresponds to an FE with centroid coordinates within the coarray
!    "box" on this image.
!
!    Assumption!! This is a 3D problems, so the centroid is
!    defined by 3 coordinates, hence centr(3).
!
!    MCEN stands for Mixed CENtroid data type.
!    LCENTR stands for *Local* array of CENTRoids.
!  NOTE
!    This is *private* array, hence the name does not start
!    with "cgca_pfem".
!  USED BY
!    Many routines of this module.
!*roboend*


!*robodoc*d* cgca_m2pfem/cgca_pfem_centroid_tmp
!  NAME
!    cgca_pfem_centroid_tmp
!  SYNOPSIS
  
type rca
  real( kind=cgca_pfem_iwp ), allocatable :: r(:,:)
end type rca
type( rca ) :: cgca_pfem_centroid_tmp[*]

!  DESCRIPTION
!    RCA stands for Rugged CoArray.
!    cgca_pfem_centroid_tmp is a temporary scalar *coarray* of derived
!    type with allocatable array component, storing centroids of ParaFEM
!    finite elements on this image.
!    The array might be of different length on different images,
!    so have to use an allocatable component of a coarray variable
!    of derived type.
!  USED BY
!    routines of this module + end user
!*roboend*


!*robodoc*d* cgca_m2pfem/cgca_pfem_stress
!  NAME
!    cgca_pfem_stress
!  SYNOPSIS

type type_stress
  real( kind=cgca_pfem_iwp ), allocatable :: stress(:,:,:)
end type type_stress
type( type_stress ) :: cgca_pfem_stress[*]

!  DESCRIPTION
!    This type has a single allocatable array to store
!    all stress components for all integration points
!    for all elements on an image.
!    This data will be read by all images.
!*roboend*

  
!*robodoc*d* cgca_m2pfem/cgca_pfem_integrity
!  NAME
!    cgca_pfem_integrity
!  SYNOPSIS

type cgca_pfem_integ_type
 real, allocatable :: i(:)
end type cgca_pfem_integ_type
type( cgca_pfem_integ_type) :: cgca_pfem_integrity[*]

!  DESCRIPTION
!    A derived type is needed because the length of the integrity
!    array will differ from image to image. So this is a scalar coarray
!    of derived type with a single component: allocatable array of
!    integrity, i. i=1 means to damage, i=0 means no remaining load
!    bearing capacity.
!    This data will be used to update the Young's modulus
!  USED BY
!    cgca_pfem_uym + end user?
!*roboend*
  
  
!*robodoc*d* cgca_m2pfem/cgca_pfem_enew
!  NAME
!    cgca_pfem_enew
!  SYNOPSIS
  
real( kind=cgca_pfem_iwp ), allocatable :: cgca_pfem_enew(:,:)

!  DESCRIPTION
!    Naming: E New as in new Young's modulus. This *local* array
!    stores Young's moduli for each integration point of each
!    FE on this image.
!  USED BY
!    cgca_pfem_uym + end user
!*roboend*


contains


!*robodoc*s* cgca_m2pfem/cgca_pfem_integalloc
!  NAME
!    cgca_pfem_integalloc
!  SYNOPSIS

subroutine cgca_pfem_integalloc( nels_pp )

!  INPUT
!    nels_pp - elements per MPI process (per image).

integer, intent(in) :: nels_pp

!  SIDE EFFECTS
!    Allocatable array component cgca_pfem_integrity%i becomes allocated
!  DESCRIPTION
!    This routine allocates cgca_pfem_integrity%i on this image.
!    This is a *local*, non-coarray, array. Hence this routine can be
!    called by any or all images. It should be called by all images,
!    of course.
!  USES
!    cgca_pfem_integrity via host association.
!  USED BY
!    end user?
!  SOURCE

integer :: errstat=0

allocate( cgca_pfem_integrity%i( nels_pp ),  stat=errstat )
if ( errstat .ne. 0 ) error stop                                       &
  "ERROR: cgca_pfem_integalloc: allocate( cgca_pfem_integrity%i )"

end subroutine cgca_pfem_integalloc

!*roboend*
 

!*robodoc*s* cgca_m2pfem/cgca_pfem_integdalloc
!  NAME
!    cgca_pfem_integdalloc
!  SYNOPSIS

subroutine cgca_pfem_integdalloc

!  SIDE EFFECTS
!    Allocatable array component of cgca_pfem_integrity coarray becomes
!    deallocated.
!  DESCRIPTION
!    This routine deallocates allocatable array component of integrity:
!    cgca_pfem_integrity%i.
!  USES
!    cgca_pfem_integrity via host association
!  USED BY
!    end user?
!  SOURCE

integer :: errstat=0

deallocate( cgca_pfem_integrity%i, stat=errstat )
if ( errstat .ne. 0 ) error stop                                       &
  "ERROR: cgca_pfem_integalloc: deallocate( cgca_pfem_integrity%i )"

end subroutine cgca_pfem_integdalloc

!*roboend*


!*robodoc*s* cgca_m2pfem/cgca_pfem_ealloc
!  NAME
!    cgca_pfem_ealloc
!  SYNOPSIS

subroutine cgca_pfem_ealloc( nip, nels_pp )
    
!  INPUTS
!    nip - integer, number of integration points
!    nels_pp - elements per MPI process (per image).

integer, intent( in ) :: nip, nels_pp
    
!  SIDE EFFECTS
!    Allocatable *local* array enew becomes allocated
!  DESCRIPTION
!    This routine allocates an allocatable *local* array.
!    The allocatable array stores the Young's modulus
!    per FE element and integration point.
!  USED BY
!    end user?
!  SOURCE

integer :: errstat=0

allocate( cgca_pfem_enew( nip, nels_pp ), stat=errstat )
if ( errstat .ne. 0 ) error stop                                       &
  "ERROR: cgca_pfem_ealloc: allocate( cgca_pfem_enew )"

end subroutine cgca_pfem_ealloc

!*roboend*

  
!*robodoc*s* cgca_m2pfem/cgca_pfem_edalloc
!  NAME
!    cgca_pfem_edalloc
!  SYNOPSIS

subroutine cgca_pfem_edalloc
    
!  SIDE EFFECTS
!    Allocatable *local* array cgca_pfem_enew becomes deallocated.
!  DESCRIPTION
!    This routine deallocates an allocatable *local* array used to
!    store the Young's modulus per FE element and integration point
!  USES
!    cgca_pfem_enew via host association.
!  USED BY
!    end user?
!  SOURCE

integer :: errstat=0

deallocate( cgca_pfem_enew, stat=errstat )
if ( errstat .ne. 0 ) error stop                                       &
  "ERROR: cgca_pfem_ealloc: deallocate( cgca_pfem_enew )"
    
end subroutine cgca_pfem_edalloc

!*roboend*


!*robodoc*s* cgca_m2pfem/cgca_pfem_uym
!  NAME
!    cgca_pfem_uym
!  SYNOPSIS

subroutine cgca_pfem_uym( e_orig, nels_pp )

!  INPUTS
!     e_orig - *real* is the original Young's modulus.
!    For now assume a single value, i.e. all int points
!    have identical original value.
!     nels_pp - number of FEs for this image.

real( kind=cgca_pfem_iwp ), intent(in) :: e_orig
integer, intent(in) :: nels_pp

!  SIDE EFFECTS
!    The Young's modulus gets updated with integrity
!  DESCRIPTION
!    UYM stands for Update Young's Modulus
!    This routine updates the value of the Young's modulus, e,
!    e = e_original * integrity.
!    Integrity - integer, cell integrity (from 0.0 to 1.0)
!  NOTES
!    Purely local routine, not coarray operations.
!  USES
!    cgca_pfem_enew, cgca_pfem_integrity, all via host association.
!  USED BY
!    end user?
!  SOURCE

integer :: fe

do fe = 1, nels_pp
 cgca_pfem_enew( : , fe ) =  e_orig * cgca_pfem_integrity%i( fe )

 ! It seems the Young's modulus of 0 causes instability.
 ! So don't let it get to 0, use a small factor, e.g. 1.0e-3 instead.
 cgca_pfem_enew(:,fe) = max( cgca_pfem_enew( : , fe ) , 1.0e-3*e_orig )

end do

end subroutine cgca_pfem_uym

!*roboend*
  

!*robodoc*s* cgca_m2pfem/cgca_pfem_ctalloc
!  NAME
!    cgca_pfem_ctalloc
!  SYNOPSIS

subroutine cgca_pfem_ctalloc( ndim, nels_pp )

!  INPUTS

integer, intent( in ) :: ndim, nels_pp

!    ndim - integer, number of DOF per node.
!    nels_pp - elements per MPI process (per image).
!  SIDE EFFECTS
!    Allocatable array component cgca_pfem_centroid_tmp%r becomes
!    allocated.
!  DESCRIPTION
!    CTA stands for Centroids Temporary Allocate.
!    This routine allocates an allocatable array component of scalar
!    coarray cgca_pfem_centroid_tmp.
!    The allocatable array stores FE centroid coordinates
!    together with their numbers and MPI ranks where these are stored.
!  NOTES
!    The array component can of different length on different images.
!  USES
!  USED BY
!  SOURCE

integer :: errstat=0

allocate( cgca_pfem_centroid_tmp%r( ndim, nels_pp ), stat=errstat )
if ( errstat .ne. 0 ) error stop                                       &
  "ERROR: cgca_pfem_ctalloc: allocate( cgca_pfem_centroid_tmp%r )"

end subroutine cgca_pfem_ctalloc

!*roboend*


!*robodoc*s* cgca_m2pfem/cgca_pfem_ctdalloc
!  NAME
!    cgca_pfem_ctdalloc
!  SYNOPSIS

subroutine cgca_pfem_ctdalloc

!  SIDE EFFECTS
!    Allocatable array component of cgca_pfem_centroid_tmp becomes
!    deallocate
!  DESCRIPTION
!    CTD stands for Centroids Temporary Allocate.
!    This routine deallocates an allocatable array component of coar.
!    This must be done only after all images copied the contents of
!    type( rca ) :: cgca_pfem_centroid_tmp[*] into their local,
!    *not* coarray centroid arrays.
!  USES
!  USED BY
!  SOURCE

integer :: errstat=0

deallocate( cgca_pfem_centroid_tmp%r, stat=errstat )
if ( errstat .ne. 0 ) error stop                                       &
  "ERROR: cgca_pfem_ctd: deallocate( cgca_pfem_centroid_tmp%r )"

end subroutine cgca_pfem_ctdalloc

!*roboend*


!*robodoc*s* cgca_m2pfem/cgca_pfem_cenc
!  NAME
!    cgca_pfem_cenc
!  SYNOPSIS

subroutine cgca_pfem_cenc( origin, rot, bcol, bcou )

!  INPUTS

real( kind=rdef ), intent( in ) ::                                     &
 origin(3),        & ! origin of the "box" cs, in FE cs
 rot(3,3),         & ! rotation tensor *from* FE cs *to* CA cs
 bcol(3),          & ! lower phys. coords of the coarray on image
 bcou(3)             ! upper phys. coords of the coarray on image

!  SIDE EFFECTS
!    Array lcentr is changed.
!  DESCRIPTION
!    CENC stands for CENtroids Collection.
!    This routine reads centroids of all elements from all MPI
!    processes and adds those with centroids within its
!    CA "box" to its cgca_pfem_centroid array.
!  NOTES
!    This routine must be called only after coarray
!    cgca_pfem_centroid_tmp has been established on all images.
!    This routine accesses coarrays on other images, hence
!    sync must be used before calling this routine.
!  USES
!    lcentr via host association.
!  SOURCE

integer :: errstat, i, j

! centroid coords in CA cs
real( kind=cgca_pfem_iwp ), allocatable :: cen_ca(:)

! allocate the CA cs centroid array for a single point
! set to zero initially
allocate( cen_ca( size(cgca_pfem_centroid_tmp%r, dim=1) ),             &
         source = 0.0_cgca_pfem_iwp, stat=errstat )
if ( errstat .ne. 0 ) error stop                                       &
  "ERROR: cgca_pfem_cenc: cannot allocate( cen_ca )"

! Allocate to size 0. Any new value added to the array
! automatically increases its size.
allocate( lcentr( 0 ), stat=errstat )
if ( errstat .ne. 0 ) error stop                                       &
  "ERROR: cgca_pfem_cenc: cannot allocate( lcentr )"

! loop over all images
images: do i=1, num_images()

  ! loop over all elements on that image
  elements: do j = 1, size( cgca_pfem_centroid_tmp[i]%r, dim=2 )

    ! Convert centroid coordinates from FE cs to CA cs
    ! cgca_pfem_centroid_tmp[i]        - variable on image i
    ! cgca_pfem_centroid_tmp[i]%r      - component that is the centroids
    !                                    real array
    ! cgca_pfem_centroid_tmp[i]%r(:,j) - take finite element j, and all
    !                                    centroid coordinates for it.
    cen_ca = matmul( rot, (cgca_pfem_centroid_tmp[i]%r(:,j) - origin) )

    ! Check whether CA cs centroid is within the box.
    ! If all CA cs centroid coordinates are greater or equal to
    ! the lower bound of the box, and all of them are also
    ! less of equal to the upper bound of the box, then the centroid
    ! is inside. Then add the new entry. The allocatable array is
    ! automaticaly enlarged.
    if ( all( cen_ca .ge. bcol ) .and. all( cen_ca .le. bcou ) )       &
      lcentr = (/ lcentr, mcen( i, j, cen_ca ) /)

  end do elements
end do images

end subroutine cgca_pfem_cenc

!*roboend*


!*robodoc*s* cgca_m2pfem/cgca_pfem_cendmp
!  NAME
!    cgca_pfem_cendmp
!  SYNOPSIS

subroutine cgca_pfem_cendmp

!  SIDE EFFECTS
!    Dumps some data to stdout
!  DESCRIPTION
!    CENDMP stands for CENtroids array dump.
!    This routine dumps cgca_pfem_centroid to stdout.
!  NOTES
!    Must call from all images.
!  SOURCE

integer :: i, img

img = this_image()

do i = 1, size( lcentr )
  write (*,*) "CA on img"      , img             ,                     &
              "<-> FE"         , lcentr(i)%elnum ,                     &
              "on img"         , lcentr(i)%image ,                     &
              "centr. in CA cs", lcentr(i)%centr
end do

end subroutine cgca_pfem_cendmp

!*roboend*


!*robodoc*s* cgca_m2pfem/cgca_pfem_salloc
!  NAME
!    cgca_pfem_salloc
!  SYNOPSIS

subroutine cgca_pfem_salloc( nels_pp, intp, comp )

!  INPUTS
!    nels_pp - number of elements on this image
!    intp - number of integration points per element
!    comp - number of stress tensor components

integer, intent( in ) :: nels_pp, intp, comp

!  SIDE EFFECTS
!    Allocatable component array cgca_pfem_stress%stress becomes
!    allocated
!  DESCRIPTION
!    SALLOC stands for Allocate Stress tensor array.
!    This routine allocates an allocatable array component of coar.
!    The allocatable array stores all stress tensor components,
!    for all integration points on all elements on an image.
!  USES
!    cgca_pfem_iwp, host association
!  USED BY
!    end user
!  SOURCE

integer :: errstat=0

allocate( cgca_pfem_stress%stress( nels_pp, intp, comp ),              &
          source=0.0_cgca_pfem_iwp, stat=errstat )
if ( errstat .ne. 0 ) error stop                                       &
  "ERROR: cgca_pfem_salloc: allocate( cgca_pfem_stress%stress )"

end subroutine cgca_pfem_salloc

!*roboend*


!*robodoc*s* cgca_m2pfem/cgca_pfem_sdalloc
!  NAME
!    cgca_pfem_sdalloc
!  SYNOPSIS

subroutine cgca_pfem_sdalloc

!  SIDE EFFECTS
!    allocatable array cgca_pfem_stress%stress become deallocated
!  DESCRIPTION
!    SDALLOC stands for Deallocate Stress tensor array.
!    This routine deallocates allocatable array component of coar.
!    This routine should be called only when the analysis is complete.
!    Any and every image can call this routine.
!  USES
!    cgca_pfem_stress%stress, host association
!  USED BY
!  SOURCE

integer :: errstat=0

deallocate( cgca_pfem_stress%stress, stat=errstat )
if ( errstat .ne. 0 ) error stop                                       &
  "ERROR: cgca_pfem_ctd: deallocate( cgca_pfem_stress%stress )"

end subroutine cgca_pfem_sdalloc

!*roboend*


  !*robodoc*s* cgca_m2pfem/cgca_pfem_sdmp
  !  NAME
  !    cgca_pfem_sdmp
  !  SYNOPSIS

  subroutine cgca_pfem_sdmp

    !  SIDE EFFECTS
    !    Dumps some data to stdout
    !  DESCRIPTION
    !    SDMP stands for Stress tensor dump.
    !    This routine dumps stress tensors to stdout.
    !  NOTES
    !    Must call from all images.
    !  SOURCE

    integer :: img, nel, nintp, el, intp

    img = this_image()
    nel = size( cgca_pfem_stress%stress, dim=1 )
    nintp = size( cgca_pfem_stress%stress, dim=2 )

    do   el = 1, nel
       do intp = 1, nintp
          write (*,*) "img", img, "FE", el, "int p.", intp, "stress",          &
               cgca_pfem_stress%stress( el, intp, : )
       end do
    end do

  end subroutine cgca_pfem_sdmp
  !*roboend*


!*robodoc*s* cgca_m2pfem/cgca_pfem_simg
!  NAME
!    cgca_pfem_simg
!  SYNOPSIS

subroutine cgca_pfem_simg( simg )

!  OUTPUT
!    simg - mean stress tensor over all integration points on all
!    finite elements linked to CA on this image.
!    Note that I use CGPACK kind, because this var will be input to
!    a CGPACK routine.

real( kind=rdef ), intent(out) :: simg(3,3)

!  DESCRIPTION
!    SIMG stands for mean Stress on an Image.
!    The routine reads all stress tensors from all integration
!    points for all elements which are linked to CA on this image
!    and calculates the mean value.
!    This value is then used to pass to the cleavage routine.
!  SOURCE

! Running total stress array.
! The number of components is read from cgca_pfem_stress%stress
! Note I use the CGPACK real kind!
real( kind=rdef ) :: stot( size( cgca_pfem_stress%stress, dim=3 ) )

integer :: el, nel, rel, intp, nintp, rimg

! Total number of elements linked to CA model on this image
! Total number of int. points per element
  nel = size( lcentr )
nintp = size( cgca_pfem_stress%stress, dim=2 )

! Add all stress tensors together. Loop over all elements linked
! to CA on this image and over all int. points.
stot = 0.0_rdef
do el=1, nel
  do intp=1, nintp

    ! Calculate the image and the element numbers to read the stress
    ! data from.
    rimg = lcentr(el)%image
     rel = lcentr(el)%elnum

    ! Add the new tensor values to the running total.
    ! Convert from ParaFEM real kind to CGPACK.
    ! Hopefully, any potential loss of accuracy is of no importance.
     stot = stot +                                                     &
       real( cgca_pfem_stress[ rimg ]%stress( rel, intp, :), kind=rdef )

  end do
end do

! Construct a (3,3) matrix from (6) vector.
! Observe the component order of ParaFEM
!   sx=stress(1)
!   sy=stress(2)
!   sz=stress(3)
!   txy=stress(4)
!   tyz=stress(5)
!   tzx=stress(6)
!   sigm=(sx+sy+sz)/three
!https://code.google.com/p/parafem/source/browse/trunk/parafem/src/modules/shared/new_library.f90
simg(1,1) = stot(1)
simg(2,2) = stot(2)
simg(3,3) = stot(3)
simg(1,2) = stot(4)
simg(2,3) = stot(5)
simg(3,1) = stot(6)
simg(2,1) = simg(1,2)
simg(3,2) = simg(2,3)
simg(1,3) = simg(3,1)

! calculate the mean
simg = simg / real( nel*nintp, kind=rdef )

end subroutine cgca_pfem_simg

!*roboend*


!*robodoc*s* cgca_m2pfem/cgca_pfem_intcalc1
!  NAME
!    cgca_pfem_intcalc1
!  SYNOPSIS

subroutine cgca_pfem_intcalc1( arrsize, fracvol )

!  INPUTS
!     arrsize - contains the 3 sizes of the space coarray.
!    Using the coarray sizes, the characteristic coarray
!    area is calculated.
!     fracvol - *real*, the number of failed (fractured) cells for each
!    image. It is calculated by cgca_fv, which *must* be called prior to
!    calling this routine.

integer( kind=iarr ), intent( in ) :: arrsize(3)
real( kind=rdef), intent( in ) :: fracvol

!  SIDE EFFECTS
!    cgca_pfem_integrity array changes
!  DESCRIPTION
!    All FEs linked to this image get the same value of integrity.
!    These are all FEs in lcentr array. For entry i in this
!    array this is FE cgca_pfem_integrity( lcentr(i)%elnum )
!    on image lcentr(i)%image.
!
!    The integrity is 1 minus the ratio of number of
!    failed cells to the cracteristic coarray area.
!    If integrity < 0, set it to zero.
!  USES
!    lcentr via host association
!  USED BY
!    end user?
!  SOURCE

real :: carea ! characteristic area
integer, parameter :: kind_integ=kind( cgca_pfem_integrity%i )
real( kind=kind_integ ), parameter :: one=1_kind_integ
integer( kind=idef ) :: i

! Volume, in cells, is the product of 3 coarray sizes.
! Don't forget to remove the halo cells! 
! Characteristic area is volume ** 2/3
carea = product( real( arrsize-2 ) ) ** 0.66666666666666666667

do i = 1, size( lcentr )

 ! integrity is calculate as: i = 1 - min(1,f),
 ! Integrity has the range [1..0], where i=1 for f=0, i=0 for f=1.  
 ! f=fracvol / carea   - Fraction of failed cells, 0 if no fracture,
 !                       1 or above when I consider the CA to have no
 !                       load bearing capacity. 
 ! min( 1, fraction) - To make sure fraction is [0..1].
 cgca_pfem_integrity[ lcentr(i)%image ] % i( lcentr(i)%elnum ) =       &
          one - min( one, fracvol / carea )
end do

end subroutine cgca_pfem_intcalc1

  !*roboend*

!*********************************************************************72

end module cgca_m2pfem
