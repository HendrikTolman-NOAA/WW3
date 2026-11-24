!> @file
!> @brief SHOWEX bottom friction source term (Ardhuin et al 2003).
!>
!> @author F. Ardhuin
!> @author J. Lepesqueur
!> @date   14-Mar-2012
!>

#include "w3macros.h"
!/ ------------------------------------------------------------------- /
!>
!> @brief SHOWEX bottom friction source term (Ardhuin et al. 2003).
!>
!> @details Using a subgrid depth parameterization based on Tolman (CE 1995).
!>
!> @author F. Ardhuin
!> @author J. Lepesqueur
!> @date   14-Mar-2012
!>
!> @copyright Copyright 2009-2022 National Weather Service (NWS),
!>       National Oceanic and Atmospheric Administration.  All rights
!>       reserved.  WAVEWATCH III is a trademark of the NWS.
!>       No unauthorized use without permission.
!>
MODULE W3SBT4MD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III           NOAA/NCEP |
  !/                  |    F. Ardhuin and J. Lepesqueur   |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         14-Mar-2012 |
  !/                  +-----------------------------------+
  !/
  !/    20-Dec-2004 : Origination.                        ( version 3.06 )
  !/    23-Jun-2006 : Formatted for submitting code for   ( version 3.09 )
  !/                  inclusion in WAVEWATCH III.
  !/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
  !/    14-Mar-2012 : Preparing distribution version.     ( version 4.05 )
  !/
  !/    Copyright 2009 National Weather Service (NWS),
  !/       National Oceanic and Atmospheric Administration.  All rights
  !/       reserved.  WAVEWATCH III is a trademark of the NWS.
  !/       No unauthorized use without permission.
  !/
  !  1. Purpose :
  !
  !     SHOWEX bottom friction source term (Ardhuin et al. 2003),
  !     using a subgrid depth parameterization based on Tolman (CE 1995).
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
  !      W3SBT4    Subr. Public   SHOWEX bottom friction (movable bed)
  !      INSBT4    Subr. Public   Corresponding initialization routine.
  !      TABU_ERF  Subr. Public   Tabulation of ERF function
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines and functions used :
  !
  !     See subroutine documentation.
  !
  !  5. Remarks :
  !
  !     The computations are done by the C++ routines TABU_ERF_CPP and
  !     W3SBT4_CPP (BT4_ww3.cpp, declared in BT4_ww3.h). The ERF table
  !     is kept in this module.
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
  !

  PUBLIC
  !
  ! Parameters for ERF function
  !
  INTEGER, PARAMETER      :: SIZEERFTABLE=300
  REAL                    :: ERFTABLE(0:SIZEERFTABLE)
  REAL                    :: DELXERF
  REAL,    PARAMETER      :: XERFMAX =  4. ! number of stdev
  !/
CONTAINS


  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Initialization for bottom friction source term routine.
  !>
  !> @author F. Ardhuin
  !> @date   14-Mar-2012
  !>
  SUBROUTINE INSBT4
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |                         SHOM      |
    !/                  |            F. Ardhuin             |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         14-Mar-2012 |
    !/                  +-----------------------------------+
    !/
    !/    14-Mar-2012 : Origination.                        ( version 4.05 )
    !
    !  1. Purpose :
    !
    !     Initialization for bottom friction source term routine.
    !
    !  2. Method :
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      TABU_ERF  Subr. W3SBT4MD Tabulation of ERF function.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      W3SBT4    Subr. W3SRC3MD Corresponding source term.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !       None.
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    !     None.
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /
    !
    IMPLICIT NONE
    !/
    !/ ------------------------------------------------------------------- /
    !
    ! 1.  .... ----------------------------------------------------------- *
    !
    CALL TABU_ERF   !tabulates ERF function
    !/
    !/ End of INSBT4 ----------------------------------------------------- /
    !/
  END SUBROUTINE INSBT4
  ! ----------------------------------------------------------------------

  !>
  !> @brief Tabulation of ERF function, which is used in bottom friction subgrid modeling.
  !>
  !> @details Initialization for source term routine.
  !>
  !> @author J. Lepesqueur
  !> @date   14-Mar-2012
  !>
  SUBROUTINE TABU_ERF
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |        J. Lepesqueur              |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         14-Mar-2012 |
    !/                  +-----------------------------------+
    !/
    !/    14-Mar-2012 : Origination.                        ( version 3.13 )
    !/
    !  1. Purpose :
    !     Tabulation of ERF function, which is used in bottom friction subgrid modelling
    !
    !     Initialization for source term routine.
    !
    !  2. Method :
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
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
    !      W3SIN3    Subr. W3SRC3MD Corresponding source term.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !       None.
    !
    !  7. Remarks :
    !
    !     The table is computed by the C++ routine TABU_ERF_CPP
    !     (BT4_ww3.cpp, declared in BT4_ww3.h).
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    !     None.
    !
    ! 10. Source code :
    !
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    !
    INTERFACE
      FUNCTION TABU_ERF_CPP(ERFTABLE, SIZEERFTABLE, XERFMAX) BIND(C, name="tabu_erf_cpp") &
           RESULT(DELXERF)
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(C_INT), VALUE, INTENT(IN) :: SIZEERFTABLE
        REAL(C_FLOAT), VALUE, INTENT(IN) :: XERFMAX
        REAL(C_FLOAT), INTENT(OUT) :: ERFTABLE(0:SIZEERFTABLE)
        REAL(C_FLOAT) :: DELXERF
      END FUNCTION TABU_ERF_CPP
    END INTERFACE
    !
    DELXERF = TABU_ERF_CPP(ERFTABLE, INT(SIZEERFTABLE, C_INT), REAL(XERFMAX, C_FLOAT))
    RETURN
    !/ ------------------------------------------------------------------- /
  END SUBROUTINE TABU_ERF
  !/ ------------------------------------------------------------------- /

  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Computes the SHOWEX bottom friction with movable bed effects.
  !>
  !> @details Uses a Gaussian distribution for friction factors, and estimates
  !>  the contribution of rippled and non-rippled fractions based on the
  !>  bayesian approach of Tolman (1995).
  !>
  !> @param[in]    A        Action density spectrum.
  !> @param[in]    CG       Group velocities.
  !> @param[in]    WN       Wavenumbers.
  !> @param[in]    DEPTH    Water depth.
  !> @param[in]    D50      Median grain size.
  !> @param[in]    PSIC     Critical Shields parameter.
  !> @param[out]   TAUBBL   Components of stress leaking to the bottom.
  !> @param[inout] BEDFORM  Ripple parameters (roughness and wavelength).
  !> @param[out]   S        Source term (1-D version).
  !> @param[out]   D        Diagonal term of derivative.
  !> @param[in]    IX       Spatial grid index.
  !> @param[in]    IY       Spatial grid index.
  !>
  !> @author F. Ardhuin
  !> @author J. Lepesqueur
  !> @date   15-Mar-2012
  !>
  SUBROUTINE W3SBT4 (A, CG, WN, DEPTH, D50, PSIC, TAUBBL, BEDFORM, S, D, IX, IY )
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |            F. Ardhuin             |
    !/                  !            J. Lepesqueur          !
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         15-Mar-2012 |
    !/                  +-----------------------------------+
    !/
    !/    23-Jun-2011 : Origination.                        ( version 4.04 )
    !/    04-Jul-2011 : Adding momentum flux TAUBBL         ( version 4.05 )
    !/    15-Mar-2012 : Adding subgrid treatment for depth  ( version 4.05 )
    !/
    !  1. Purpose :
    !
    !     Computes the SHOWEX bottom friction with movable bed effects
    !
    !  2. Method :
    !     Uses a Gaussian distribution for friction factors, and estimates
    !     the contribution of rippled and non-rippled fractions based on
    !     the bayesian approach of Tolman (1995).
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       A       R.A.  I   Action density spectrum.
    !       CG      R.A.  I   Group velocities.
    !       WN      R.A.  I   Wavenumbers.
    !       DEPTH   Real  I   Water depth.
    !       D50     Real  I   Median grain size.
    !       PSIC    Real  I   Critical Shields parameter
    !       BEFORMS Real I/O  Ripple parameters (roughness and wavelength).
    !       TAUBBL  Real  O   Components of stress leaking to the bottom.
    !       S       R.A.  O   Source term (1-D version).
    !       D       R.A.  O   Diagonal term of derivative.             *)
    !       IX,IY   Int. I   Spatial grid indices
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      INSBT4    Subr. W3SBT4MD Initialization (first call).
    !     ----------------------------------------------------------------
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
    !     - The source term is computed by the C++ routine W3SBT4_CPP
    !       (BT4_ww3.cpp, declared in BT4_ww3.h).
    !     - IX and IY are no longer used, but are kept in the parameter
    !       list so that the calling routines are unchanged.
    !     - ECOS and ESIN are passed for the full spectrum, as TAUBBL
    !       uses ECOS(IS) and ESIN(IS) with the spectral index IS.
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    !     None.
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /

    USE CONSTANTS, ONLY: GRAV, SED_SG, FWTABLE, SIZEFWTABLE, ABMIN, DELAB
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, DDEN, SBTCX, ECOS, ESIN
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    !/
    INTERFACE
      SUBROUTINE W3SBT4_CPP(A, CG, WN, DEPTH, D50, PSIC, TAUBBL, BEDFORM, S, D,  &
           NTH, NK, SIG, DDEN, ECOS, ESIN, SBTCX, FWTABLE, SIZEFWTABLE,         &
           ABMIN, DELAB, ERFTABLE, SIZEERFTABLE, XERFMAX, DELXERF, GRAV,       &
           SED_SG) BIND(C, name="w3sbt4_cpp")
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(C_INT), VALUE, INTENT(IN) :: NTH, NK, SIZEFWTABLE, SIZEERFTABLE
        REAL(C_FLOAT), VALUE, INTENT(IN) :: DEPTH, D50, PSIC, ABMIN, DELAB,    &
             XERFMAX, DELXERF, GRAV, SED_SG
        REAL(C_FLOAT), INTENT(IN) :: A(NTH*NK), CG(NK), WN(NK), SIG(NK),       &
             DDEN(NK), ECOS(NTH*NK), ESIN(NTH*NK), SBTCX(7),                   &
             FWTABLE(0:SIZEFWTABLE), ERFTABLE(0:SIZEERFTABLE)
        REAL(C_FLOAT), INTENT(OUT) :: TAUBBL(2), BEDFORM(3), S(NTH*NK), D(NTH*NK)
      END SUBROUTINE W3SBT4_CPP
    END INTERFACE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    REAL, INTENT(IN)        :: CG(NK), WN(NK), DEPTH, A(NSPEC), D50
    REAL, INTENT(IN)        :: PSIC
    INTEGER, INTENT(IN)     :: IX, IY
    REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC), TAUBBL(2)
    REAL, INTENT(INOUT)     :: BEDFORM(3)
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters
    !/
    LOGICAL, SAVE           :: FIRST = .TRUE.
    !/
    !/ ------------------------------------------------------------------- /
    !/
    !
    ! 0.  Initializations ------------------------------------------------ *
    !
    IF ( FIRST ) THEN
      CALL INSBT4
      FIRST  = .FALSE.
    END IF
    !
    ! 1.  Compute source term via C++ ------------------------------------- *
    !
    ! SIG(1:NK) is passed as SIG is dimensioned 0:NK+1.
    !
    CALL W3SBT4_CPP(A, CG, WN, REAL(DEPTH, C_FLOAT), REAL(D50, C_FLOAT),      &
         REAL(PSIC, C_FLOAT), TAUBBL, BEDFORM, S, D, INT(NTH, C_INT),          &
         INT(NK, C_INT), SIG(1:NK), DDEN(1:NK), ECOS(1:NSPEC), ESIN(1:NSPEC),  &
         SBTCX(1:7), FWTABLE, INT(SIZEFWTABLE, C_INT), REAL(ABMIN, C_FLOAT),   &
         REAL(DELAB, C_FLOAT), ERFTABLE, INT(SIZEERFTABLE, C_INT),             &
         REAL(XERFMAX, C_FLOAT), REAL(DELXERF, C_FLOAT), REAL(GRAV, C_FLOAT),  &
         REAL(SED_SG, C_FLOAT))
    !
    RETURN
    !/
    !/ End of W3SBT4 ----------------------------------------------------- /
    !/
  END SUBROUTINE W3SBT4

  !/ ------------------------------------------------------------------- /


END MODULE W3SBT4MD
