!$Id: cgca_m1co.f90 14 2014-12-01 10:14:16Z mexas $

!*robodoc*m* CGPACK/cgca_m1co
!  NAME
!    cgca_m1co
!  SYNOPSIS

module cgca_m1co

!  DESCRIPTION
!    Lowest level module, contains named global constants, e.g. kinds.
!    Contains routines which do not use modules (level 1 routines).
!  COPYRIGHT
!    See CGPACK_Copyright
!  CONTAINS
!    Various constants (parameters)
!  USED BY
!    Probably all higher level modules:
!    cgca_m2alloc, cgca_m2gb, cgca_m2out, cgca_m2rot
!  SOURCE

implicit none

! These 2008 functions have not been implemented yet in ifort 12.0
! But they might be useful in future.
!
!use iso_fortran_env, only: compiler_options, compiler_version
!private compiler_options, compiler_version
!character(*),parameter :: compiled_by = compiler_version()
!character(*),parameter :: compiled_with = compiler_options()

!*roboend*

!*robodoc*p* cgca_m1co/rdef
!  PARAMETER

integer, parameter :: rdef = selected_real_kind( 6, 30 )

!  DESCRIPTION
!    Default real kind
!*roboend*

!*robodoc*p* cgca_m1co/rlrg
!  PARAMETER

integer, parameter :: rlrg = selected_real_kind( 15, 300 )

!  DESCRIPTION
!    High precision real kind, most likely will be double
!    precision
!*roboend*

!*robodoc*p* cgca_m1co/idef
!  PARAMETER

integer, parameter :: idef = selected_int_kind( 8 )

!  DESCRIPTION
!    Default integer kind
!*roboend*

!*robodoc*p* cgca_m1co/iarr
!  PARAMETER

integer, parameter :: iarr = selected_int_kind( 5 )

!  DESCRIPTION
!    Integer kind for cellular arrays
!*roboend*

!*robodoc*p* cgca_m1co/ilrg
!  PARAMETER

integer, parameter :: ilrg = selected_int_kind( 10 )

! DESCRIPTION
!   Integer kind for large numbers, e.g. volumes, total
!   number of cells, etc.
!*roboend*

!*robodoc*p* cgca_m1co/ldef
!  PARAMETER

integer, parameter :: ldef = kind( .true. )

!  DESCRIPTION
!    Default logical kind
!*roboend*

!*robodoc*p* cgca_m1co/pi
!  PARAMETER

real( kind=rdef ), parameter :: pi = 3.14159265358979323846264338_rdef

!  DESCRIPTION
!    pi
!*roboend*

!*robodoc*p* cgca_m1co/cgca_state_type_grain
!  PARAMETER

integer( kind=idef ), parameter :: cgca_state_type_grain = 1_idef

! DESCRIPTION
!   Cell state type for grains
!*roboend*

!*robodoc*p* cgca_m1co/cgca_state_type_frac
!  PARAMETER

integer( kind=idef ), parameter :: cgca_state_type_frac = 2_idef

! DESCRIPTION
!   Cell state type for fractures
!*roboend*

!*robodoc*p* cgca_m1co/cgca_liquid_state
!  PARAMETER

integer( kind=iarr ), parameter :: cgca_liquid_state = 0_iarr

!  DESCRIPTION
!    Liquid phase, cell state of type cgca_state_type_grain.
!    All states of the same type must be unique.
!*roboend*

!*robodoc*p* cgca_m1co/cgca_intact_state
!  PARAMETER

integer( kind=iarr ), parameter :: cgca_intact_state = huge( 0_iarr )

!  DESCRIPTION
!    Intact state for fracture array, cell state of type cgca_state_type_frac.
!    All states of the same type must be unique.
!*roboend*

!*robodoc*p* cgca_m1co/cgca_gb_state_intact
!  PARAMETER

integer( kind=iarr ), parameter :: cgca_gb_state_intact = 2_iarr

!  DESCRIPTION
!    Intact grain boundary, cell state of type cgca_state_type_frac.
!    All states of the same type must be unique.
!*roboend*

!*robodoc*p* cgca_m1co/cgca_gb_state_fractured
!  PARAMETER

integer( kind=iarr ), parameter :: cgca_gb_state_fractured = 1_iarr

!  DESCRIPTION
!    Fractured grain boundary, cell state of type cgca_state_type_frac.
!    All states of the same type must be unique.
!*roboend*

!*robodoc*p* cgca_m1co/cgca_clvg_state_100_flank
!  PARAMETER

integer( kind=iarr ), parameter :: cgca_clvg_state_100_flank = 0_iarr

!  DESCRIPTION
!    Flanks of a cleavage crack of {100} family.
!    Cell state of type cgca_state_type_frac.
!    All states of the same type must be unique.
!*roboend*

!*robodoc*p* cgca_m1co/cgca_clvg_state_110_flank
!  PARAMETER

integer( kind=iarr ), parameter :: cgca_clvg_state_110_flank = -1_iarr

!  DESCRIPTION
!    Flanks of a cleavage crack of {110} family.
!    Cell state of type cgca_state_type_frac.
!    All states of the same type must be unique.
!*roboend*

!*robodoc*p* cgca_m1co/cgca_clvg_state_111_flank
!  PARAMETER

integer( kind=iarr ), parameter :: cgca_clvg_state_111_flank = -2_iarr

!  DESCRIPTION
!    Flanks of a cleavage crack of {111} family.
!    Cell state of type cgca_state_type_frac.
!    All states of the same type must be unique.
!*roboend*

!*robodoc*p* cgca_m1co/cgca_clvg_state_100_edge
!  PARAMETER

integer( kind=iarr ), parameter :: cgca_clvg_state_100_edge = -3_iarr

!  DESCRIPTION
!    Edges of a cleavage crack of {100} family.
!    Cell state of type cgca_state_type_frac.
!    All states of the same type must be unique.
!*roboend*

!*robodoc*p* cgca_m1co/cgca_clvg_state_110_edge
!  PARAMETER

integer( kind=iarr ), parameter :: cgca_clvg_state_110_edge = -4_iarr

!  DESCRIPTION
!    Edges of a cleavage crack of {110} family.
!    Cell state of type cgca_state_type_frac.
!    All states of the same type must be unique.
!*roboend*

!*robodoc*p* cgca_m1co/cgca_clvg_state_111_edge
!  PARAMETER

integer( kind=iarr ), parameter :: cgca_clvg_state_111_edge = -5_iarr

!  DESCRIPTION
!    Edges of a cleavage crack of {111} family.
!    Cell state of type cgca_state_type_frac.
!    All states of the same type must be unique.
!*roboend*

!*robodoc*p* cgca_m1co/cgca_clvg_lowest_state
!  PARAMETER

integer( kind=iarr ), parameter ::                             &
 cgca_clvg_lowest_state = cgca_clvg_state_111_edge

!  DESCRIPTION
!    The the lowest cleavage state, used for sizing the lower
!    bound of e.g. the grain volume array.
!    Cell state of type cgca_state_type_frac.
!    All states of the same type must be unique.
!*roboend*

!*robodoc*p* cgca_m1co/cgca_clvg_states_flank
!  PARAMETER

integer(kind=iarr), parameter ::        &
 cgca_clvg_states_flank(3) =            &
  (/ cgca_clvg_state_100_flank,         &
     cgca_clvg_state_110_flank,         &
     cgca_clvg_state_111_flank /)

!  DESCRIPTION
!    Array to store all flank cleavage states.
!    Cell state of type cgca_state_type_frac.
!*roboend*

!*robodoc*p* cgca_m1co/cgca_clvg_states_edge
!  PARAMETER

integer(kind=iarr), parameter ::        &
 cgca_clvg_states_edge(3)  =            &
  (/ cgca_clvg_state_100_edge,          &
     cgca_clvg_state_110_edge,          & 
     cgca_clvg_state_111_edge /)

!  DESCRIPTION
!    Array to store all edge cleavage states.
!    Cell state of type cgca_state_type_frac.
!*roboend*

!*robodoc*p* cgca_m1co/cgca_clvg_states
!  PARAMETER

integer(kind=iarr), parameter ::                      &
 cgca_clvg_states( size(cgca_clvg_states_flank) +     &
                   size(cgca_clvg_states_edge) ) =   &
  (/ cgca_clvg_states_flank, cgca_clvg_states_edge /)

!  DESCRIPTION
!    Array to store all cleavage states, flanks and edges.
!    Cell state of type cgca_state_type_frac.
!*roboend*

!*robodoc*p* cgca_m1co/cgca_lowest_state
!  PARAMETER

integer(kind=iarr), parameter :: cgca_lowest_state = -huge(0_iarr)

!  DESCRIPTION
!    Lowest possible state in the model
!*roboend*

end module cgca_m1co
