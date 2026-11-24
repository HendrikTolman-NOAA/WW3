!> @file
!> @brief Contains module W3STR1MD.
!>
!> @author A. J. van der Westhuysen @date 13-Jan-2013
!> @author A. Roland @date 22-Feb-2023
!>

#include "w3macros.h"
!/ ------------------------------------------------------------------- /
!>
!> @brief Module for inclusion of triad nonlinear interaction
!>  according to Eldeberky's (1996) Lumped Triad Interaction (LTA)
!>  source term.
!>
!> @author A. J. van der Westhuysen @date 13-Jan-2013
!>
MODULE W3STR1MD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III           NOAA/NCEP |
  !/                  |     A. J. van der Westhuysen      |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         13-Jan-2013 |
  !/                  +-----------------------------------+
  !/
  !/    13 Jan-2013 : Origination, based on SWAN v40.91 code ( version 4.08 )
  !/
  !/    Copyright 2009 National Weather Service (NWS),
  !/       National Oceanic and Atmospheric Administration.  All rights
  !/       reserved.  WAVEWATCH III is a trademark of the NWS.
  !/       No unauthorized use without permission.
  !/
  !  1. Purpose :
  !
  !     Module for inclusion of triad nonlinear interaction according to
  !     Eldeberky's (1996) Lumped Triad Interaction (LTA) source term.
  !
  !  2. Variables and types :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !     ----------------------------------------------------------------
  !
  !  3. Subroutines and functions :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !      W3STR1    Subr. Public   User supplied triad interactions.
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines and functions used :
  !
  !     See subroutine documentation.
  !
  !  5. Remarks :
  !
  !     WAVEWATCH III is designed as a highly plug-compatible code.
  !     Source term modules can be included as self-contained modules,
  !     with limited changes needed to the interface of routine calls
  !     in W3SRCE, and in the point postprocessing programs only.
  !     Codes submitted for inclusion in WAVEWATCH III should be
  !     self-contained in the way described below, and might be
  !     provided with distributions fully integrated in the data
  !     structure, or as an optional version of this module to be
  !     included by the user.
  !
  !     Rules for preparing a module to be included in or distributed
  !     with WAVEWATCH III :
  !
  !      - Fully document the code following the outline given in this
  !        file, and according to all other WAVEWATCH III routines.
  !      - Provide a file with necessary modifications to W3SRCE and
  !        all other routines that require modification.
  !      - Provide a test case with expected results.
  !      - It is strongly recommended that the programming style used
  !        in WAVEWATCH III is followed, in particular
  !          a) for readability, write as if in fixed FORTRAN format
  !             regarding column use, even though all files are F90
  !             free format.
  !          b) I prefer upper case programming for permanent code,
  !             as I use lower case in debugging and temporary code.
  !
  !     This module needs to be self-contained in the following way.
  !
  !      a) All saved variables connected with this source term need
  !         to be declared in the module header. Upon acceptance as
  !         permanent code, they will be converted to the WAVEWATCH III
  !         dynamic data structure.
  !      b) Provide a separate computation and initialization routine.
  !         In the submission, the initialization should be called
  !         from the computation routine upon the first call to the
  !         routine. Upon acceptance as permanent code, the
  !         initialization routine will be moved to a more appropriate
  !         location in the code (i.e., being absorbed in ww3_grid or
  !         being moved to W3IOGR).
  !
  !     See notes in the file below where to add these elements.
  !
  !  6. Switches :
  !
  !     None.
  !
  !  7. Source code :
  !/
  !/ ------------------------------------------------------------------- /
  !/
  !     *****************************************
  !     ***    Declare saved variables here   ***
  !     ***  public or private as appropriate ***
  !     *****************************************
  !
  PUBLIC
  !/
CONTAINS
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Triad interaction source term computed using the Lumped
  !>  Triad Appproximation (LTA) of Eldeberky (1996).
  !>
  !> @verbatim
  !>     The parametrized biphase is given by:
  !>
  !>                                  0.2
  !>     beta = - pi/2 + pi/2 tanh ( ----- )
  !>                                   Ur
  !>
  !>     where Ur is the Ursell number.
  !>
  !>     The source term as function of frequency p is:
  !>
  !>             +      -
  !>     S(p) = S(p) + S(p)
  !>
  !>     in which
  !>
  !>      +
  !>     S(p) = alpha Cp Cg,p (R(p/2,p/2))**2 sin (|beta|) ( E(p/2)**2 -2 E(p) E(p/2) )
  !>
  !>      -          +
  !>     S(p) = - 2 S(2p)
  !>
  !>     with alpha a tunable coefficient and R(p/2,p/2) is the interaction
  !>     coefficient of which the expression can be found in Eldeberky (1996).
  !>
  !>     Note that a slightly adapted formulation of the LTA is used in
  !>     in the SWAN model:
  !>
  !>     - Only positive contributions to higher harmonics are considered
  !>       here (no energy is transferred to lower harmonics).
  !>
  !>     - The mean frequency in the expression of the Ursell number
  !>       is calculated according to the first order moment over the
  !>       zeroth order moment (personal communication, Y.Eldeberky, 1997).
  !>
  !>     - The interactions are calculated up to 2.5 times the mean
  !>       frequency only.
  !>
  !>     - Since the spectral grid is logarithmically distributed in frequency
  !>       space, the interactions between central bin and interacting bin
  !>       are interpolated such that the distance between these bins is
  !>       factor 2 (nearly).
  !>
  !>     - The interactions are calculated in terms of energy density
  !>       instead of action density. So the action density spectrum
  !>       is firstly converted to the energy density grid, then the
  !>       interactions are calculated and then the spectrum is converted
  !>       to the action density spectrum back.
  !> @endverbatim
  !>
  !> @param[in] A Action density spectrum (1-D)
  !> @param[in] CG Group velocities.
  !> @param[in] WN Wavenumbers.
  !> @param[in] DEPTH Mean water depth.
  !> @param[in] IX
  !> @param[out] S Source term (1-D version).
  !> @param[out] D Diagonal term of derivative (1-D version).
  !>
  !> @author A. J. van der Westhuysen  @date 13-Jan-2013
  !>
  SUBROUTINE W3STR1 (A, CG, WN, DEPTH, IX, S, D)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |     A. J. van der Westhuysen      |
    !/                  |     A. Roland                     |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         13-Jan-2013 |
    !/                  +-----------------------------------+
    !/
    !/    13 Jan-2013 : Origination, based on SWAN v40.91 code ( version 4.08 )
    !/    05 Oct-2016 : Avoiding divide by zero for EMEAN      ( version 5.15 )
    !/    28 Feb-2023 : Improvement of efficiency and stability ( version 7.xx)
    !/
    !  1. Purpose :
    !
    !     Triad interaction source term computed using the Lumped Triad
    !     Appproximation (LTA) of Eldeberky (1996).
    !
    !  2. Method :
    !
    !     (Taken from SWAN v40.91, based on code by Marcel Zijlema, TU Delft)
    !
    !     The parametrized biphase is given by:
    !
    !                                  0.2
    !     beta = - pi/2 + pi/2 tanh ( ----- )
    !                                   Ur
    !
    !     where Ur is the Ursell number.
    !
    !     The source term as function of frequency p is:
    !
    !             +      -
    !     S(p) = S(p) + S(p)
    !
    !     in which
    !
    !      +
    !     S(p) = alpha Cp Cg,p (R(p/2,p/2))**2 sin (|beta|) ( E(p/2)**2 -2 E(p) E(p/2) )
    !
    !      -          +
    !     S(p) = - 2 S(2p)
    !
    !     with alpha a tunable coefficient and R(p/2,p/2) is the interaction
    !     coefficient of which the expression can be found in Eldeberky (1996).
    !
    !     Note that a slightly adapted formulation of the LTA is used in
    !     in the SWAN model:
    !
    !     - Only positive contributions to higher harmonics are considered
    !       here (no energy is transferred to lower harmonics).
    !
    !     - The mean frequency in the expression of the Ursell number
    !       is calculated according to the first order moment over the
    !       zeroth order moment (personal communication, Y.Eldeberky, 1997).
    !
    !     - The interactions are calculated up to 2.5 times the mean
    !       frequency only.
    !
    !     - Since the spectral grid is logarithmically distributed in frequency
    !       space, the interactions between central bin and interacting bin
    !       are interpolated such that the distance between these bins is
    !       factor 2 (nearly).
    !
    !     - The interactions are calculated in terms of energy density
    !       instead of action density. So the action density spectrum
    !       is firstly converted to the energy density grid, then the
    !       interactions are calculated and then the spectrum is converted
    !       to the action density spectrum back.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       A       R.A.  I   Action density spectrum (1-D)
    !       CG      R.A.  I   Group velocities.
    !       WN      R.A.  I   Wavenumbers.
    !       DEPTH   Real  I   Mean water depth.
    !       FMEAN   Real  I   Mean wave frequency.
    !       S       R.A.  O   Source term (1-D version).
    !       D       R.A.  O   Diagonal term of derivative (1-D version).
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !     None.
    !
    !  5. Called by :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      W3SRCE    Subr. W3SRCEMD Source term integration.
    !      W3EXPO    Subr.   N/A    Point output post-processor.
    !      GXEXPO    Subr.   N/A    GrADS point output post-processor.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !       None.
    !
    !  7. Remarks :
    !
    !     - The source term is computed by the C++ routine W3STR1_CPP
    !       (TR1_ww3.cpp, declared in TR1_ww3.h).
    !     - IX is no longer used, but is kept in the parameter list so
    !       that the calling routines are unchanged.
    !     - As before, S and D are not set if the spectrum is zero, and
    !       are not set for the highest frequency band if the
    !       interactions are computed. W3SRCE sets them to zero first.
    !
    !  8. Structure :
    !
    !     Determine resonance condition and the maximum discrete freq.
    !     for which the interactions are calculated.
    !
    !     If Ursell number larger than prescribed value compute interactions
    !        Calculate biphase
    !        Do for each direction
    !           Convert action density to energy density
    !           Do for all frequencies
    !             Calculate interaction coefficient and interaction factor
    !             Compute interactions and store results in matrix
    !
    !  9. Switches :
    !
    !     None.
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /
    USE CONSTANTS, ONLY: GRAV
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, DDEN, FTE, FTF
    USE, INTRINSIC :: ISO_C_BINDING
    !/
    IMPLICIT NONE
    !/
    INTERFACE
      SUBROUTINE W3STR1_CPP(A, CG, WN, DEPTH, S, D, NTH, NK, SIG, DDEN, FTE, FTF, GRAV) &
           BIND(C, name="w3str1_cpp")
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(C_INT), VALUE, INTENT(IN) :: NTH, NK
        REAL(C_FLOAT), VALUE, INTENT(IN) :: DEPTH, FTE, FTF, GRAV
        REAL(C_FLOAT), INTENT(IN) :: A(NTH*NK), CG(NK), WN(NK), SIG(NK), DDEN(NK)
        REAL(C_FLOAT), INTENT(INOUT) :: S(NTH*NK), D(NTH*NK)
      END SUBROUTINE W3STR1_CPP
    END INTERFACE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    REAL, INTENT(IN)        :: CG(NK), WN(NK), DEPTH, A(NSPEC)
    INTEGER, INTENT(IN)     :: IX
    REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
    !/
    !/ ------------------------------------------------------------------- /
    !/
    !
    ! 1.  Compute source term via C++ ------------------------------------- *
    !
    ! SIG(1:NK) is passed as SIG is dimensioned 0:NK+1.
    !
    CALL W3STR1_CPP(A, CG, WN, REAL(DEPTH, C_FLOAT), S, D,                &
                    INT(NTH, C_INT), INT(NK, C_INT), SIG(1:NK), DDEN(1:NK), &
                    REAL(FTE, C_FLOAT), REAL(FTF, C_FLOAT), REAL(GRAV, C_FLOAT))
    !
    RETURN
    !/
    !/ End of W3STR1 ----------------------------------------------------- /
    !/
  END SUBROUTINE W3STR1
  !/ ------------------------------------------------------------------- /
END MODULE W3STR1MD
