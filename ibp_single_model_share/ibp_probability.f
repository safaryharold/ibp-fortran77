C=================================================================
C     IBP_PROBABILITY.F
C
C     Callable, side-effect-free interface to the IBP (Ionospheric
C     Bubble Probability) model, for linking into a larger/global
C     model rather than running as a standalone program.
C
C     This is NOT a line-for-line port of pbsingle.f. It shares the
C     same physics functions (pbutils.f) and produces the same
C     numbers, but two things were deliberately changed because a
C     component embedded in a host model has different requirements
C     than a standalone diagnostic program:
C
C       1. NO STOP. pbsingle.f calls STOP on a bad coefficient
C          (lifetime too small); that would halt the ENTIRE host model,
C          not just this calculation. Here, that condition is
C          reported through the IER output argument instead, and
C          PROB is set to 0 so the caller gets a defined, safe value
C          either way.
C       2. NO PRINT. pbsingle.f writes prompts and DEBUG lines to
C          stdout, which a library routine should never do. Problems
C          are reported through IER/WARN_* flags instead, leaving it
C          to the caller whether to log them.
C
C     The y0 overflow term is CLAMPED before EXP(y0) rather than left
C     to overflow to Infinity for extreme inputs. pbbatch.f and
C     pbsingle.f apply the identical clamp; the only difference here
C     is that this library routine raises the WARN_Y0 flag instead of
C     printing, so a host model can log that it happened rather than
C     emit stdout. See the WARN_Y0 flag below.
C
C     A second difference, also for the same reason: when IDOY is out
C     of range, pbsingle.f defaults the *month* to January but still
C     passes the original, invalid IDOY into the seasonal-oscillation
C     function (GOSC), which both extrapolates its trig functions
C     outside their intended domain AND prints its own warning to
C     stdout - a second, hidden I/O side effect this routine cannot
C     allow. Here, an out-of-range IDOY is sanitized to 1 (January)
C     BEFORE it is used anywhere, including in GOSC, so behavior is
C     self-consistent and no output is ever printed. WARN_DOY still
C     tells the caller this happened.
C
C     Requires: pbutils.f (compiled/linked alongside), pbcoeff.inc,
C     and coefficients loaded once via READCF before the first call
C     (see IBP_INIT below).
C=================================================================

      SUBROUTINE IBP_INIT(CFFILE, IER)
C     Loads model coefficients from a text file into the shared
C     COMMON block. Call this once, before any call to IBP_PROBABILITY.
C     IER: 0 = success, 1 = could not open/read the coefficient file.
      IMPLICIT NONE
      CHARACTER*(*) CFFILE
      INTEGER IER
      EXTERNAL READCF

      IER = 0
      CALL READCF(CFFILE, IER)
      RETURN
      END


      SUBROUTINE IBP_PROBABILITY(IDOY, ALON, F107IN, TIN,
     +                           PROB, IER, WARN_DOY, WARN_F107,
     +                           WARN_Y0)
C     Computes the IBP occurrence probability for one set of inputs.
C     No I/O, no STOP - safe to call repeatedly from a host model's
C     inner loop.
C
C     Inputs:
C       IDOY    INTEGER          Day of year, nominally 1-365
C       ALON    DOUBLE PRECISION Longitude, degrees East (any range;
C                                 wraps automatically)
C       F107IN  DOUBLE PRECISION F10.7 solar flux index, sfu
C                                 (nominal valid range 80-200 sfu)
C       TIN     DOUBLE PRECISION Local time, hours (any range; wraps
C                                 automatically, e.g. 25.5 = 01:30)
C
C     Outputs:
C       PROB    DOUBLE PRECISION Occurrence probability, clamped to
C                                 [0,1]. Meaningful even when a WARN_*
C                                 flag is set; only IER=1 means PROB
C                                 could not be computed at all (0.0
C                                 is returned in that case).
C       IER     INTEGER          0 = success.
C                                 1 = coefficients were not successfully
C                                     loaded, or model coefficient TAU (1/gamma, mean lifetime)
C                                     is at or near zero - PROB set to
C                                     0.0 and not physically meaningful.
C                                 2 = a double-precision input was NaN or
C                                     non-finite - PROB set to 0.0.
C       WARN_DOY   LOGICAL       .TRUE. if IDOY was outside 1-365; the
C                                 month was defaulted to January (1)
C                                 and the calculation proceeded anyway.
C       WARN_F107  LOGICAL       .TRUE. if F107IN was outside the
C                                 80-200 sfu range the model was fit
C                                 over. Calculation proceeded anyway.
C       WARN_Y0    LOGICAL       .TRUE. if an internal exponent was
C                                 large enough that, without the
C                                 overflow clamp described above, the
C                                 unclamped calculation could have
C                                 produced Infinity/NaN. PROB is still
C                                 a safe, defined value in [0,1].
      IMPLICIT NONE
      INCLUDE 'pbcoeff.inc'

      INTEGER IDOY, IER
      DOUBLE PRECISION ALON, F107IN, TIN, PROB
      LOGICAL WARN_DOY, WARN_F107, WARN_Y0

      INTEGER IMONTH, IDOY_SAFE
      DOUBLE PRECISION TTRANS, GOSCV, TSVAL, PHIVAL
      DOUBLE PRECISION LAMBD1, LAMBD2, XMU, SIGM, TINTRN
      DOUBLE PRECISION OMULT, XIMULT, X0, Y0
      DOUBLE PRECISION IMSQRT, OMULT2
      DOUBLE PRECISION ERFARG, XINTEG, EXPARG
      DOUBLE PRECISION PI

      DOUBLE PRECISION GOSC, TSFUNC, PHIMON, ERF77, DMODLO
      INTEGER MONOF
      EXTERNAL GOSC, TSFUNC, PHIMON, ERF77, DMODLO, MONOF

      PI = ACOS(-1.0D0)
      IER = 0
      WARN_DOY  = .FALSE.
      WARN_F107 = .FALSE.
      WARN_Y0   = .FALSE.
      PROB = 0.0D0

      IF (IBP_INITED .NE. 1) THEN
         IER = 1
         RETURN
      END IF

      IF (ALON .NE. ALON .OR. F107IN .NE. F107IN .OR.
     +    TIN .NE. TIN .OR.
     +    ABS(ALON) .GT. 1.0D300 .OR.
     +    ABS(F107IN) .GT. 1.0D300 .OR.
     +    ABS(TIN) .GT. 1.0D300) THEN
         IER = 2
         RETURN
      END IF

      IDOY_SAFE = IDOY
      IMONTH = MONOF(IDOY)
      IF (IMONTH .LT. 1 .OR. IMONTH .GT. 12) THEN
         WARN_DOY = .TRUE.
         IMONTH = 1
         IDOY_SAFE = 1
      END IF

      IF (F107IN .LT. 80.0D0 .OR. F107IN .GT. 200.0D0) THEN
         WARN_F107 = .TRUE.
      END IF

      TTRANS = DMODLO(TIN + 12.0D0, 24.0D0) - 12.0D0
      GOSCV  = GOSC(IDOY_SAFE)
      PHIVAL = PHIMON(ALON, IMONTH)
      LAMBD1 = (GOSCV + C1 + C2*F107IN) * PHIVAL
      IF (LAMBD1 .LT. 0.0D0) LAMBD1 = 0.0D0

      TSVAL  = TSFUNC(ALON, IMONTH)
      XMU    = MU0 + TSVAL
      SIGM   = SIGMA
      TINTRN = TTRANS

      IF (TAU .LE. 1.0D-12) THEN
         IER  = 1
         PROB = 0.0D0
         RETURN
      END IF
      LAMBD2 = 1.0D0 / TAU

      IF (SIGM .LE. 1.0D-12) THEN
         PROB = 0.0D0
      ELSE IF (LAMBD1 .LE. 1.0D-12) THEN
         PROB = 0.0D0
      ELSE
         OMULT  = -1.0D0 / SQRT(2.0D0 * PI * SIGM**2)
         XIMULT = -1.0D0 / (2.0D0 * SIGM**2)

         IF (XIMULT .GE. -1.0D-12) THEN
            PROB = 0.0D0
         ELSE
            X0 = XMU - LAMBD2 / (2.0D0 * XIMULT)
            Y0 = LAMBD2 * (XMU - TINTRN - LAMBD2/(4.0D0*XIMULT))
            IMSQRT = SQRT(-XIMULT)

C           Clamp before EXP(y0) - see header comment. pbbatch.f and
C           pbsingle.f apply this same clamp.
            IF (Y0 .GT. 100.0D0) THEN
               WARN_Y0 = .TRUE.
               OMULT2 = -1.0D30
            ELSE IF (Y0 .LT. -100.0D0) THEN
               WARN_Y0 = .TRUE.
               OMULT2 = 0.0D0
            ELSE
               OMULT2 = SQRT(PI) * OMULT * EXP(Y0) / IMSQRT
            END IF

            ERFARG = (TINTRN - X0) * IMSQRT
            XINTEG = (0.5D0 + ERF77(ERFARG)*0.5D0) * OMULT2
            EXPARG = LAMBD1 * XINTEG

            IF (EXPARG .LT. -700.0D0) THEN
               PROB = 1.0D0
            ELSE IF (EXPARG .GT. 1.0D-9) THEN
               PROB = 0.0D0
            ELSE
               PROB = 1.0D0 - EXP(EXPARG)
            END IF
         END IF
      END IF

      PROB = MAX(0.0D0, MIN(1.0D0, PROB))

      RETURN
      END
