!> @file
!> @brief JONSWAP bottom friction routine.
!>
!> @author H. L. Tolman
!> @date   29-May-2009
!>

#include "w3macros.h"
!/ ------------------------------------------------------------------- /
!>
!> @brief JONSWAP bottom friction routine.
!>
!> @author H. L. Tolman
!> @date   29-May-2009
!>
!> @copyright Copyright 2009-2022 National Weather Service (NWS),
!>       National Oceanic and Atmospheric Administration.  All rights
!>       reserved.  WAVEWATCH III is a trademark of the NWS.
!>       No unauthorized use without permission.
!>
MODULE W3SBT1MD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III           NOAA/NCEP |
  !/                  |           H. L. Tolman            |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         29-May-2009 |
  !/                  +-----------------------------------+
  !/
  !/    For updates see W3SBT1 documentation.
  !/
  !  1. Purpose :
  !
  !     JONSWAP bottom friction routine.
  !
  !  2. Variables and types :
  !
  !  3. Subroutines and functions :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !      W3SBT1    Subr. Public   JONSWAP source term.
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines and functions used :
  !
  !     See subroutine documentation.
  !
  !  5. Remarks :
  !
  !  6. Switches :
  !
  !     See subroutine documentation.
  !
  !  7. Source code :
  !/
  !/ ------------------------------------------------------------------- /
  !/
  PUBLIC
  !/
CONTAINS
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Bottom friction source term according to the empirical JONSWAP
  !>  formulation.
  !>
  !> @verbatim
  !>              2 GAMMA   /    CG         \      SBTC1 /     \       .
  !>       Sbt = ---------- | ------- - 0.5 | E  = ----- | ... | E    (1)
  !>             GRAV DEPTH \  SI/WN        /      DEPTH \     /
  !>
  !>     Where GAMMA = -0.038 m2/s3 (JONSWAP)
  !>                 = -0.067 m2/s3 (Bouws and Komen 1983)
  !>
  !>     In the routine, the constant 2 GAMMA / GRAV = SBTC1.
  !> @endverbatim
  !>
  !> @param[in]  A      Action density spectrum (1-D).
  !> @param[in]  CG     Group velocities.
  !> @param[in]  WN     Wavenumbers.
  !> @param[in]  DEPTH  Mean water depth.
  !> @param[out] S      Source term (1-D version).
  !> @param[out] D      Diagonal term of derivative (1-D version).
  !>
  !> @author H. L. Tolman
  !> @date   29-May-2009
  !>
  SUBROUTINE W3SBT1 (A, CG, WN, DEPTH, S, D)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           H. L. Tolman            |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         29-May-2009 |
    !/                  +-----------------------------------+
    !/
    !/    05-Dec-1996 : Final FORTRAN 77.                   ( version 1.18 )
    !/    08-Dec-1999 : Upgrade to FORTRAN 90.              ( version 2.00 )
    !/    20-Dec-2004 : Multiple model version.             ( version 3.06 )
    !/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
    !/
    !/    Copyright 2009 National Weather Service (NWS),
    !/       National Oceanic and Atmospheric Administration.  All rights
    !/       reserved.  WAVEWATCH III is a trademark of the NWS.
    !/       No unauthorized use without permission.
    !/
    !  1. Purpose :
    !
    !     Bottom friction source term according to the empirical JONSWAP
    !     formulation.
    !
    !  2. Method :
    !
    !              2 GAMMA   /    CG         \      SBTC1 /     \       .
    !       Sbt = ---------- | ------- - 0.5 | E  = ----- | ... | E    (1)
    !             GRAV DEPTH \  SI/WN        /      DEPTH \     /
    !
    !     Where GAMMA = -0.038 m2/s3 (JONSWAP)
    !                 = -0.067 m2/s3 (Bouws and Komen 1983)
    !
    !     In the routine, the constant 2 GAMMA / GRAV = SBTC1.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       A       R.A.  I   Action density spectrum (1-D)
    !       CG      R.A.  I   Group velocities.
    !       WN      R.A.  I   Wavenumbers.
    !       DEPTH   Real  I   Mean water depth.
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
    !     None.
    !
    !  7. Remarks :
    !
    !     The source term is computed by the C++ routine W3SBT1_CPP
    !     (BT1_ww3.cpp, declared in BT1_ww3.h).
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
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, SBTC1
    USE, INTRINSIC :: ISO_C_BINDING
    !
    IMPLICIT NONE
    !/
    INTERFACE
      SUBROUTINE W3SBT1_CPP(A, CG, WN, DEPTH, S, D, NTH, NK, SIG, SBTC1) BIND(C, name="w3sbt1_cpp")
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(C_INT), VALUE, INTENT(IN) :: NTH, NK
        REAL(C_FLOAT), VALUE, INTENT(IN) :: DEPTH, SBTC1
        REAL(C_FLOAT), INTENT(IN) :: A(NTH*NK), CG(NK), WN(NK), SIG(NK)
        REAL(C_FLOAT), INTENT(OUT) :: S(NTH*NK), D(NTH*NK)
      END SUBROUTINE W3SBT1_CPP
    END INTERFACE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    REAL, INTENT(IN)        :: CG(NK), WN(NK), DEPTH, A(NSPEC)
    REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
    !/
    !/ ------------------------------------------------------------------- /
    !/
    !
    ! 1.  Compute source term via C++ ------------------------------------- *
    !
    ! SIG(1:NK) is passed as SIG is dimensioned 0:NK+1.
    !
    CALL W3SBT1_CPP(A, CG, WN, REAL(DEPTH, C_FLOAT), S, D, &
                    INT(NTH, C_INT), INT(NK, C_INT), SIG(1:NK), REAL(SBTC1, C_FLOAT))
    !
    RETURN
    !/
    !/ End of W3SBT1 ----------------------------------------------------- /
    !/
  END SUBROUTINE W3SBT1
  !/
  !/ End of module W3SBT1MD -------------------------------------------- /
  !/
END MODULE W3SBT1MD
