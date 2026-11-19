#include "w3macros.h"
!/ ------------------------------------------------------------------- /
MODULE W3SLN1MD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III           NOAA/NCEP |
  !/                  |           H. L. Tolman            |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         29-May-2009 |
  !/                  +-----------------------------------+
  !/
  !/    23-Jun-2006 : Origination.                        ( version 3.09 )
  !/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
  !/
  !/    Copyright 2009 National Weather Service (NWS),
  !/       National Oceanic and Atmospheric Administration.  All rights
  !/       reserved.  WAVEWATCH III is a trademark of the NWS.
  !/       No unauthorized use without permission.
  !/
  !  1. Purpose :
  !
  !     Linear wind input according to Cavaleri and Melanotte-Rizzoli
  !     (1982) filtered for low frequencies according to Tolman (1992).
  !
  !  2. Variables and types :
  !
  !  3. Subroutines and functions :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !      W3SLN1    Subr. Public   User supplied linear input.
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines and functions used :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !      STRACE    Subr. W3SERVMD Subroutine tracing.
  !     ----------------------------------------------------------------
  !
  !  5. Remarks :
  !
  !  6. Switches :
  !
  !     !/S  Enable subroutine tracing.
  !     !/T  Test output.
  !
  !  7. Source code :
  !/
  !/ ------------------------------------------------------------------- /
  !/
  PUBLIC
  !/
CONTAINS
  !/ ------------------------------------------------------------------- /
  SUBROUTINE W3SLN1 (K, FHIGH, USTAR, USDIR, S)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           H. L. Tolman            |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         23-Jun-2006 |
    !/                  +-----------------------------------+
    !/
    !/    23-Jun-2006 : Origination.                        ( version 3.09 )
    !/
    !  1. Purpose :
    !
    !     Linear wind input according to Cavaleri and Melanotte-Rizzoli
    !     (1982) filtered for low frequencies according to Tolman (1992).
    !
    !  2. Method :
    !
    !     The expression of Cavaleri and Melanotte-Rizzoli, converted to
    !     action spectra defined in terms of wavenumber and direction
    !     becomes
    !
    !                       -1       /     /                \ \ 4
    !       Sln  = SLNC1 * k   * max | 0., | U* cos(Dtheta) | |        (1)
    !                                \     \                / /
    !
    !                             2     -2
    !              SLNC1 = 80 RHOr  GRAV   FILT                        (2)
    !
    !     Where :
    !
    !        RHOr     Density of air dev. by density of water.
    !        U*       Wind friction velocity.
    !        Dtheta   Difference in wind and wave direction.
    !        FILT     Filter based on PM and cut-off frequencies.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       K       R.A.  I   Wavenumber for entire spectrum.
    !       FHIGH   R.A.  I   Cut-off frequency in integration (rad/s)
    !       USTAR   Real  I   Friction velocity.
    !       USDIR   Real  I   Direction of USTAR.
    !       S       R.A.  O   Source term.
    !     ----------------------------------------------------------------
    !                         *) Stored as 1-D array with dimension NTH*NK
    !
    !  4. Subroutines used :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD Subroutine tracing.
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
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    !     !/S  Enable subroutine tracing.
    !     !/T  Test output.
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /
    USE CONSTANTS, ONLY: GRAV, RADE
    USE W3GDATMD, ONLY: NTH, NK, ECOS, ESIN, SIG, SLNC1, FSPM, FSHF
    USE W3ODATMD, ONLY: NDSE, NDST
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif
    USE, INTRINSIC :: ISO_C_BINDING
    !/
    IMPLICIT NONE
    !/
    INTERFACE
      SUBROUTINE W3SLN1_CPP(K, FHIGH, USTAR, USDIR, S, NTH, NK, ECOS, ESIN, SIG, SLNC1, FSPM, FSHF, GRAV) BIND(C, name="w3sln1_cpp")
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(C_INT), VALUE, INTENT(IN) :: NTH, NK
        REAL(C_FLOAT), VALUE, INTENT(IN) :: FHIGH, USTAR, USDIR, SLNC1, FSPM, FSHF, GRAV
        REAL(C_FLOAT), INTENT(IN) :: K(NK), ECOS(NTH), ESIN(NTH), SIG(NK)
        REAL(C_FLOAT), INTENT(OUT) :: S(NTH, NK)
      END SUBROUTINE W3SLN1_CPP
    END INTERFACE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    REAL, INTENT(IN)        :: K(NK), FHIGH, USTAR, USDIR
    REAL, INTENT(OUT)       :: S(NTH,NK)
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters
    !/
#ifdef W3_S
    INTEGER, SAVE           :: IENT = 0
#endif
    !/
    !/ ------------------------------------------------------------------- /
    !/
#ifdef W3_S
    CALL STRACE (IENT, 'W3SLN1')
#endif
    !
    ! 1.  Set up factors ------------------------------------------------- *
    !
#ifdef W3_T
    WRITE (NDST,900) USTAR, USDIR*RADE
#endif
    !
    ! 2.  Compute source term via C++ ------------------------------------- *
    !
    CALL W3SLN1_CPP(K, REAL(FHIGH, C_FLOAT), REAL(USTAR, C_FLOAT), REAL(USDIR, C_FLOAT), S, &
                    INT(NTH, C_INT), INT(NK, C_INT), ECOS, ESIN, SIG, &
                    REAL(SLNC1, C_FLOAT), REAL(FSPM, C_FLOAT), REAL(FSHF, C_FLOAT), REAL(GRAV, C_FLOAT))
    !
    RETURN
    !
    ! Formats
    !
#ifdef W3_T
900 FORMAT ( ' TEST W3SLN1 : USTAR, DIR :',F6.3, F6.1)
#endif
    !/
    !/ End of W3SLN1 ----------------------------------------------------- /
    !/
  END SUBROUTINE W3SLN1
  !/
  !/ End of module INSLN1MD -------------------------------------------- /
  !/
END MODULE W3SLN1MD
