C=================================================================
C     PBUTILS.F
C
C     Coefficient loader and model functions of the IBP model in
C     fixed-form FORTRAN 77. Language constraints and how they are met:
C       - no modules / derived types -> COMMON block (pbcoeff.inc)
C       - no allocatable arrays      -> fixed-size arrays
C       - no MODULO intrinsic        -> DMODLO (below)
C       - no ERF intrinsic (F2008)   -> ERF77, A&S 7.1.26 approx.
C       - no internal procedures     -> external subroutines
C
C     Coefficients are not embedded in the source. They are
C     read at run time from a plain-text data file by READCF; see
C     ibp_coefficients.txt.
C=================================================================

      SUBROUTINE READCF(FNAME, IER)
C     Reads IBP model coefficients from a plain-text data file into
C     the /PBCOEF/ common block.
      IMPLICIT NONE
      INCLUDE 'pbcoeff.inc'
      CHARACTER*(*) FNAME
      INTEGER IUCF, IOS, I, J, IER
      CHARACTER*100 LINE
      LOGICAL IOPEN
      DOUBLE PRECISION PHISUM

      IER = 1
      IBP_INITED = 0

C     Pick the first unopened unit in 90-99 so this library does not
C     clobber a unit already in use by a host model.
      IUCF = 0
      DO 5 IUCF = 90, 99
         INQUIRE(UNIT=IUCF, OPENED=IOPEN)
         IF (.NOT. IOPEN) GOTO 6
    5 CONTINUE
      RETURN
    6 CONTINUE

      OPEN(UNIT=IUCF, FILE=FNAME, STATUS='OLD', ACTION='READ',
     +     IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
         RETURN
      END IF

C     Title line (discarded)
      READ(IUCF,'(A)',IOSTAT=IOS) LINE
      IF (IOS .NE. 0) GOTO 900

C     Scalar params: C1 C2 TAU(=1/gamma, mean lifetime h) MU0 SIGMA
      READ(IUCF,'(A)',IOSTAT=IOS) LINE
      IF (IOS .NE. 0) GOTO 900
      READ(IUCF,*,IOSTAT=IOS) C1, C2, TAU, MU0, SIGMA
      IF (IOS .NE. 0) GOTO 900

C     Diurnal oscillation params: G0 GS(1) GC(1) GS(2) GC(2)
      READ(IUCF,'(A)',IOSTAT=IOS) LINE
      IF (IOS .NE. 0) GOTO 900
      READ(IUCF,*,IOSTAT=IOS) G0, GS(1), GC(1), GS(2), GC(2)
      IF (IOS .NE. 0) GOTO 900

C     Monthly LT-shift T0MON (12 values)
      READ(IUCF,'(A)',IOSTAT=IOS) LINE
      IF (IOS .NE. 0) GOTO 900
      READ(IUCF,*,IOSTAT=IOS) (T0MON(J), J=1,12)
      IF (IOS .NE. 0) GOTO 900

C     Monthly harmonic-1 sin/cos coefficients (12 values each)
      READ(IUCF,'(A)',IOSTAT=IOS) LINE
      IF (IOS .NE. 0) GOTO 900
      READ(IUCF,*,IOSTAT=IOS) (TSMON(1,J), J=1,12)
      IF (IOS .NE. 0) GOTO 900
      READ(IUCF,'(A)',IOSTAT=IOS) LINE
      IF (IOS .NE. 0) GOTO 900
      READ(IUCF,*,IOSTAT=IOS) (TCMON(1,J), J=1,12)
      IF (IOS .NE. 0) GOTO 900

C     Monthly harmonic-2 sin/cos coefficients (12 values each)
      READ(IUCF,'(A)',IOSTAT=IOS) LINE
      IF (IOS .NE. 0) GOTO 900
      READ(IUCF,*,IOSTAT=IOS) (TSMON(2,J), J=1,12)
      IF (IOS .NE. 0) GOTO 900
      READ(IUCF,'(A)',IOSTAT=IOS) LINE
      IF (IOS .NE. 0) GOTO 900
      READ(IUCF,*,IOSTAT=IOS) (TCMON(2,J), J=1,12)
      IF (IOS .NE. 0) GOTO 900

      NLONBN = 360
      LONST  = -180.0D0
      LONSTP = 1.0D0

C     Longitude x month phi table: 360 rows x 12 columns
      READ(IUCF,'(A)',IOSTAT=IOS) LINE
      IF (IOS .NE. 0) GOTO 900
      DO 10 I = 1, 360
         READ(IUCF,*,IOSTAT=IOS) (PHITAB(I,J), J=1,12)
         IF (IOS .NE. 0) GOTO 900
   10 CONTINUE

C     Reject physically impossible or non-finite coefficients before
C     marking the common block as usable.
      IF (.NOT. (TAU .GT. 0.0D0) .OR.
     +    .NOT. (SIGMA .GT. 0.0D0)) GOTO 900
      PHISUM = 0.0D0
      DO 20 I = 1, 360
         DO 15 J = 1, 12
            IF (.NOT. (PHITAB(I,J) .GE. 0.0D0)) GOTO 900
            PHISUM = PHISUM + PHITAB(I,J)
   15    CONTINUE
   20 CONTINUE
      IF (.NOT. (PHISUM .GT. 0.0D0)) GOTO 900

      IBP_INITED = 1
      CLOSE(IUCF)
      IER = 0
      RETURN
  900 CLOSE(IUCF)
      RETURN
      END


      DOUBLE PRECISION FUNCTION GOSC(IDOY)
C     Seasonal oscillation term g_osc(doy), Eq. (4) of the paper.
      IMPLICIT NONE
      INCLUDE 'pbcoeff.inc'
      INTEGER IDOY
      DOUBLE PRECISION DOYRAD, PI, TWOPI
      INTEGER I

      PI    = ACOS(-1.0D0)
      TWOPI = 2.0D0 * PI

      IF (IDOY .LT. 1 .OR. IDOY .GT. 365) THEN
         PRINT *, 'WARNING: DOY out of range: ', IDOY
      END IF

      DOYRAD = DBLE(IDOY) / 365.0D0 * TWOPI
      GOSC = G0
      DO 10 I = 1, 2
         GOSC = GOSC + GS(I) * SIN(DBLE(I)*DOYRAD)
     +               + GC(I) * COS(DBLE(I)*DOYRAD)
   10 CONTINUE

      RETURN
      END


      DOUBLE PRECISION FUNCTION TSFUNC(ALONDG, IMONTH)
C     Monthly local-time shift as a function of longitude.
C     Implements ts_month(lon), Eq. (6) of the paper.
      IMPLICIT NONE
      INCLUDE 'pbcoeff.inc'
      DOUBLE PRECISION ALONDG
      INTEGER IMONTH
      DOUBLE PRECISION LONRAD, PI, TWOPI
      INTEGER I

      IF (IMONTH .LT. 1 .OR. IMONTH .GT. 12) THEN
         PRINT *, 'ERROR: Invalid month: ', IMONTH
         STOP 1
      END IF

      PI     = ACOS(-1.0D0)
      TWOPI  = 2.0D0 * PI
      LONRAD = ALONDG / 360.0D0 * TWOPI

      TSFUNC = T0MON(IMONTH)
      DO 10 I = 1, 2
         TSFUNC = TSFUNC
     +          + TSMON(I,IMONTH) * SIN(DBLE(I)*LONRAD)
     +          + TCMON(I,IMONTH) * COS(DBLE(I)*LONRAD)
   10 CONTINUE

      RETURN
      END


      DOUBLE PRECISION FUNCTION PHIMON(ALONDG, IMONTH)
C     Longitude/month lookup into the PHITAB grid.
C     Longitudinal density lookup phi_month(lon) from PHITAB.
      IMPLICIT NONE
      INCLUDE 'pbcoeff.inc'
      DOUBLE PRECISION ALONDG
      INTEGER IMONTH
      DOUBLE PRECISION ANLON, RIDX
      INTEGER IDX

      IF (IMONTH .LT. 1 .OR. IMONTH .GT. 12) THEN
         PRINT *, 'ERROR: Invalid month: ', IMONTH
         STOP 1
      END IF

      ANLON = DMOD(ALONDG - LONST, 360.0D0)
      IF (ANLON .LT. 0.0D0) ANLON = ANLON + 360.0D0

      RIDX = ANLON / LONSTP + 1.0D0
      IDX  = INT(RIDX)
      IF (IDX .LT. 1)      IDX = 1
      IF (IDX .GT. NLONBN) IDX = NLONBN

      PHIMON = PHITAB(IDX, IMONTH)
      IF (PHIMON .LT. 0.0D0) PHIMON = 0.0D0

      RETURN
      END


      DOUBLE PRECISION FUNCTION DMODLO(A, P)
C     Emulates the FORTRAN 90 MODULO(A,P) intrinsic, which does not
C     exist in FORTRAN 77, using the specific double-precision
C     modulo intrinsic DMOD (always F77-standard).
      IMPLICIT NONE
      DOUBLE PRECISION A, P, R

      R = DMOD(A, P)
      IF (R .NE. 0.0D0 .AND.
     +    (R .LT. 0.0D0) .NEQV. (P .LT. 0.0D0)) THEN
         R = R + P
      END IF
      DMODLO = R

      RETURN
      END


      DOUBLE PRECISION FUNCTION ERF77(X)
C     Approximation to the error function ERF(X). Plain FORTRAN 77
C     has no ERF intrinsic (it was only added in Fortran 2008), so
C     this uses the Abramowitz & Stegun formula 7.1.26, which has a
C     maximum absolute error of 1.5D-7 - well within the precision
C     needed for the probability model.
      IMPLICIT NONE
      DOUBLE PRECISION X, T, Y, AX
      DOUBLE PRECISION A1, A2, A3, A4, A5, PCOF
      PARAMETER (A1   =  0.254829592D0)
      PARAMETER (A2   = -0.284496736D0)
      PARAMETER (A3   =  1.421413741D0)
      PARAMETER (A4   = -1.453152027D0)
      PARAMETER (A5   =  1.061405429D0)
      PARAMETER (PCOF =  0.3275911D0)

      AX = ABS(X)
      T  = 1.0D0 / (1.0D0 + PCOF*AX)
      Y  = 1.0D0 - ((((( A5*T + A4)*T + A3)*T + A2)*T + A1)*T)
     +      * EXP(-AX*AX)

      IF (X .LT. 0.0D0) THEN
         ERF77 = -Y
      ELSE
         ERF77 = Y
      END IF

      RETURN
      END


      INTEGER FUNCTION MONOF(IDOY)
C     Day-of-year to calendar month (non-leap calendar), as a
C     standalone external function.
      IMPLICIT NONE
      INTEGER IDOY
      INTEGER NDAYS(12)
      INTEGER ITEMP, I, IMON
      DATA NDAYS /31,28,31,30,31,30,31,31,30,31,30,31/

      IF (IDOY .LT. 1 .OR. IDOY .GT. 365) THEN
         MONOF = 0
         RETURN
      END IF

      ITEMP = IDOY
      IMON  = 0
      DO 10 I = 1, 12
         IF (ITEMP .LE. NDAYS(I)) THEN
            IMON = I
            GOTO 20
         ELSE
            ITEMP = ITEMP - NDAYS(I)
         END IF
   10 CONTINUE
   20 CONTINUE

      MONOF = IMON
      RETURN
      END


      INTEGER FUNCTION FSTNB(LINE)
C     Returns the column index of the first non-blank character in
C     LINE, or 0 if the line is entirely blank. Stands in for the
C     F90 ADJUSTL intrinsic (not available in F77): the caller only
C     ever needs the first non-blank character itself, so there is
C     no need to actually left-shift the whole string.
      IMPLICIT NONE
      CHARACTER*(*) LINE
      INTEGER I

      DO 10 I = 1, LEN(LINE)
         IF (LINE(I:I) .NE. ' ') THEN
            FSTNB = I
            RETURN
         END IF
   10 CONTINUE
      FSTNB = 0

      RETURN
      END


      INTEGER FUNCTION LSTNB(LINE)
C     Returns the column index of the last non-blank character in
C     LINE, or 0 if entirely blank. Companion to FSTNB; together
C     these let a caller print a fixed-length CHARACTER variable as
C     if TRIM(ADJUSTL(...)) had been applied, without needing the
C     F90 TRIM intrinsic.
      IMPLICIT NONE
      CHARACTER*(*) LINE
      INTEGER I

      DO 10 I = LEN(LINE), 1, -1
         IF (LINE(I:I) .NE. ' ') THEN
            LSTNB = I
            RETURN
         END IF
   10 CONTINUE
      LSTNB = 0

      RETURN
      END
