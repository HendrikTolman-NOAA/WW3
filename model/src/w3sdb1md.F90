!> @file
!> @brief Dummy slot for bottom friction source term.
!>
!> @author J. H. Alves
!> @author H. L. Tolman
!> @date   29-May-2009
!>

#include "w3macros.h"
!/ ------------------------------------------------------------------- /
!>
!> @brief Dummy slot for bottom friction source term.
!>
!> @author J. H. Alves
!> @author H. L. Tolman
!> @date   29-May-2009
!>
!>
!> @copyright Copyright 2009-2022 National Weather Service (NWS),
!>       National Oceanic and Atmospheric Administration.  All rights
!>       reserved.  WAVEWATCH III is a trademark of the NWS.
!>       No unauthorized use without permission.
!>
MODULE W3SDB1MD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III           NOAA/NCEP |
  !/                  |           J. H. Alves             |
  !/                  |           H. L. Tolman            |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         29-May-2009 |
  !/                  +-----------------------------------+
  !/
  !/    25-Apr-2007 : Origination of module.              ( version 3.11 )
  !/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
  !/
  !/    Copyright 2009 National Weather Service (NWS),
  !/       National Oceanic and Atmospheric Administration.  All rights
  !/       reserved.  WAVEWATCH III is a trademark of the NWS.
  !/       No unauthorized use without permission.
  !/
  !  1. Purpose :
  !
  !     Dummy slot for bottom friction source term.
  !
  !  2. Variables and types :
  !
  !  3. Subroutines and functions :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !      W3SDB1    Subr. Public   Battjes and Janssen depth-induced
  !                               breaking.
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
  !> @brief Compute depth-induced breaking using Battjes and Janssen bore
  !>  model approach.
  !>
  !> @details Note that the Miche criterion can influence wave growth.
  !>
  !> @param[in]    IX      Local grid number
  !> @param[in]    A       Action density spectrum (1-D).
  !> @param[inout] DEPTH   Mean water depth.
  !> @param[inout] EMEAN   Mean wave energy.
  !> @param[inout] FMEAN   Mean wave frequency.
  !> @param[inout] WNMEAN  Mean wave number.
  !> @param[in]    CG
  !> @param[out]   LBREAK
  !> @param[out]   S       Source term (1-D version).
  !> @param[out]   D       Diagonal term of derivative (1-D version).
  !>
  !> @author J. H. Alves
  !> @author H. L. Tolman
  !> @author A. Roland
  !> @date   08-Jun-2018
  !>
  SUBROUTINE W3SDB1 (IX, A, DEPTH, EMEAN, FMEAN, WNMEAN, CG, LBREAK, S, D )
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |                        FORTRAN 90 |
    !/                  |           J. H. Alves             |
    !/                  |           H. L. Tolman            |
    !/                  !           A. Roland               |
    !/                  | Last update :         08-Jun-2018 |
    !/                  +-----------------------------------+
    !/
    !/    25-Apr-2007 : Origination of module.              ( version 3.11 )
    !/    08-Jun-2018 : Add DEBUGDB1.                       ( version 6.04 )
    !/    03-Apr-2019 : Rewrite in terms of energy density (A. Roland,version 6.07)
    !/    03-Apr-2019 : Add Thornton & Guza, 1983          (A. Roland,version 6.07)
    !/
    !  1. Purpose :
    !
    !     Compute depth-induced breaking using Battjes and Janssen bore
    !     model approach
    !
    !  2. Method : Battjes & Janssen (1978),
    !
    !           Sbr   = Dtot/Etot*WA = D * WA
    !           Dtot  = 0.25*alpha*Qb*fm*Hmax²
    !           fm    = sigma/2Pi
    !           BB    = Hrms²/Hmax² = 8Etot/Hmax²
    !           D     = Dtot/Etot = BJALFA * sigma / pi * Qb/BB = 2 * BJALFA * fm * Qb/BB
    !
    !           AR: only valid for Hrms .le. Hm, Qb .le. 1, otherwise, in the degenrative regime it is
    !           due to Qb > 1 that all wave are broken and Hrms .le. Hmax
    !           MLIM can be used to enforce this conditions, source term will smoothly converge to this limit.
    !
    !     Where CDB   = SDBC1 = BJALFA (defaults to BJALFA = 1)
    !                   modified via ww3_grid namelist parameter BJALFA
    !           HM    = GAMMA * DEP
    !           GAMMA = SDBC2 defaults to 0.73 (mean Battjes/Janssen value)
    !                   modified via ww3_grid namelist parameter BJGAM
    !
    !     And QB is estimated by iterations using the nonlinear expression
    !
    !           1 - QB = HRMS**2
    !           ------   -------
    !            ln QB    HM**2
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       A       R.A.  I   Action density spectrum (1-D)
    !       EMEAN   Real  I   Mean wave energy.
    !       FMEAN   Real  I   Mean wave frequency.
    !       WNMEAN  Real  I   Mean wave number.
    !       DEPTH   Real  I   Mean water depth.
    !       S       R.A.  O   Source term (1-D version).
    !       D       R.A.  O   Diagonal term of derivative (1-D version).
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !       None.
    !
    !  5. Called by :
    !
    !       W3SRCE   Source term integration.
    !       W3EXPO   Point output post-processor.
    !       GXEXPO   GrADS point output post-processor.
    !
    !  6. Error messages :
    !
    !       None.
    !
    !  7. Remarks :
    !
    !     - Note that the Miche criterion con influence wave growth.
    !     - The source term is computed by the C++ routine W3SDB1_CPP
    !       (DB1_ww3.cpp, declared in DB1_ww3.h).
    !     - The IWB = 2 option could not be selected (IWB was fixed to 1)
    !       and is not part of DB1. It is the separate C++ source term DB2
    !       (DB2_ww3.cpp, declared in DB2_ww3.h).
    !     - IX and CG are no longer used, but are kept in the parameter
    !       list so that the calling routines are unchanged.
    !     - If EMEAN is zero, S, D and LBREAK are not set (as before).
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
    !/
    USE W3GDATMD, ONLY: NK, NSPEC, SDBC1, SDBC2, FDONLY
    USE, INTRINSIC :: ISO_C_BINDING
    !/
    IMPLICIT NONE
    !/
    INTERFACE
      FUNCTION W3SDB1_CPP(A, DEPTH, EMEAN, FMEAN, WNMEAN, S, D, NSPEC, SDBC1, SDBC2, FDONLY) &
           BIND(C, name="w3sdb1_cpp") RESULT(IBREAK)
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(C_INT), VALUE, INTENT(IN) :: NSPEC, FDONLY
        REAL(C_FLOAT), VALUE, INTENT(IN) :: DEPTH, EMEAN, FMEAN, WNMEAN, SDBC1, SDBC2
        REAL(C_FLOAT), INTENT(IN) :: A(NSPEC)
        REAL(C_FLOAT), INTENT(OUT) :: S(NSPEC), D(NSPEC)
        INTEGER(C_INT) :: IBREAK
      END FUNCTION W3SDB1_CPP
    END INTERFACE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    INTEGER, INTENT(IN)     :: IX ! Local grid number
    REAL, INTENT(IN)        :: A(NSPEC)
    REAL, INTENT(INOUT)     :: EMEAN, FMEAN, WNMEAN, DEPTH
    REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
    REAL, INTENT(IN)        :: CG(NK)
    LOGICAL, INTENT(OUT)    :: LBREAK
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters
    !/
    INTEGER(C_INT)          :: IBREAK
    !/
    !/ ------------------------------------------------------------------- /
    !/
    !
    ! 1.  Compute source term via C++ ------------------------------------- *
    !
    ! IBREAK < 0 : EMEAN is zero, S, D and LBREAK are not set.
    !
    IBREAK = W3SDB1_CPP(A, REAL(DEPTH, C_FLOAT), REAL(EMEAN, C_FLOAT),       &
                        REAL(FMEAN, C_FLOAT), REAL(WNMEAN, C_FLOAT), S, D,  &
                        INT(NSPEC, C_INT), REAL(SDBC1, C_FLOAT),            &
                        REAL(SDBC2, C_FLOAT), MERGE(1_C_INT, 0_C_INT, FDONLY))
    IF ( IBREAK .GE. 0 ) LBREAK = IBREAK .EQ. 1
    !
    RETURN
    !/
    !/ End of W3SDB1 ----------------------------------------------------- /
    !/
  END SUBROUTINE W3SDB1
  !/
  !/
  !/ End of module W3SDB1MD -------------------------------------------- /
  !/
END MODULE W3SDB1MD
