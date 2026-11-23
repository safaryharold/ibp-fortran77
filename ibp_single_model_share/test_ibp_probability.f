C=================================================================
C     TEST_IBP_PROBABILITY.F
C
C     Minimal example of calling the IBP model as a library routine,
C     and a self-check against reference values computed with the
C     official IBP v2.0.2 Python package (ibpcalc, exact erf) from the
C     coefficient file ..._0101.cdf. Agreement is expected to ~1e-6,
C     the accuracy of the ERF77 approximation; the test tolerance is
C     1e-5. LON=200 is evaluated as -160 and DOY=366 as DOY=1 in the
C     reference, matching the wrap-around and fallback of the code.
C
C     Build and run: see the Makefile (`make test` / `./test_ibp`).
C=================================================================
      PROGRAM TEST_IBP
      IMPLICIT NONE

      INTEGER IDOY, IER, I, NFAIL
      DOUBLE PRECISION ALON, F107IN, TIN, PROB, EXPECTED, DIFF
      LOGICAL WARN_DOY, WARN_F107, WARN_Y0

      INTEGER NCASES
      PARAMETER (NCASES = 9)
      INTEGER    T_DOY(NCASES)
      DOUBLE PRECISION T_LON(NCASES), T_F107(NCASES), T_LT(NCASES)
      DOUBLE PRECISION T_EXPECT(NCASES)

C     Reference cases and expected probabilities (IBP v2.0.2). The
C     last three cases are month-boundary regressions: DOY 91 is the
C     first of April and DOY 335 is the first of December in a
C     non-leap calendar; DOY 366 must warn and fall back to January.
      DATA T_DOY   /  80,  172,  355,  200,    1,   60,
     +               91,  335,  366 /
      DATA T_LON   /30.0D0,-75.0D0,10.0D0,200.0D0,0.0D0,-180.0D0,
     +             -75.0D0,-75.0D0,-75.0D0/
      DATA T_F107  /120.0D0,150.0D0,90.0D0,180.0D0,100.0D0,110.0D0,
     +             150.0D0,150.0D0,150.0D0/
      DATA T_LT    /20.5D0,21.0D0,19.0D0,22.25D0,0.0D0,23.9D0,
     +             21.0D0,21.0D0,21.0D0/
      DATA T_EXPECT/0.4154858D0,0.0590936D0,0.0045013D0,
     +              0.5854397D0,0.0499037D0,0.0599293D0,
     +              0.2942707D0,0.6386894D0,0.5513639D0/

      PRINT *, '=== IBP_PROBABILITY library self-test ==='

      CALL IBP_INIT('missing_ibp_coefficients.txt', IER)
      IF (IER .EQ. 0) THEN
         PRINT *, 'FAIL: missing coefficient file was accepted'
         STOP 1
      END IF
      PRINT *, 'Missing coefficient file reported through IER.'

      CALL IBP_INIT('ibp_coefficients.txt', IER)
      IF (IER .NE. 0) THEN
         PRINT *, 'FATAL: could not load coefficients, IER=', IER
         STOP 1
      END IF
      PRINT *, 'Coefficients loaded OK.'
      PRINT *, ''

      NFAIL = 0
      DO 10 I = 1, NCASES
         IDOY   = T_DOY(I)
         ALON   = T_LON(I)
         F107IN = T_F107(I)
         TIN    = T_LT(I)
         EXPECTED = T_EXPECT(I)

         CALL IBP_PROBABILITY(IDOY, ALON, F107IN, TIN, PROB, IER,
     +                         WARN_DOY, WARN_F107, WARN_Y0)

         DIFF = ABS(PROB - EXPECTED)

         PRINT *, 'Case', I, ': DOY=', IDOY, ' LON=', ALON,
     +            ' F107=', F107IN, ' LT=', TIN
         PRINT *, '   PROB=', PROB, '  expected=', EXPECTED,
     +            '  diff=', DIFF
         PRINT *, '   IER=', IER, ' WARN_DOY=', WARN_DOY,
     +            ' WARN_F107=', WARN_F107, ' WARN_Y0=', WARN_Y0

         IF (IER .NE. 0) THEN
            PRINT *, '   *** FAIL: unexpected IER ***'
            NFAIL = NFAIL + 1
         ELSE IF (DIFF .GT. 1.0D-5) THEN
            PRINT *, '   *** FAIL: probability mismatch ***'
            NFAIL = NFAIL + 1
         ELSE
            PRINT *, '   OK'
         END IF
         PRINT *, ''
   10 CONTINUE

C     --- Also demonstrate the error path: a coefficient loaded, but
C     an out-of-range DOY, should give a WARNING (not a failure) ---
      CALL IBP_PROBABILITY(400, 30.0D0, 120.0D0, 20.5D0, PROB, IER,
     +                      WARN_DOY, WARN_F107, WARN_Y0)
      PRINT *, 'Out-of-range DOY=400 case: PROB=', PROB,
     +         ' IER=', IER, ' WARN_DOY=', WARN_DOY
      IF (IER .EQ. 0 .AND. WARN_DOY) THEN
         PRINT *, '   OK (warned, did not stop the program)'
      ELSE
         PRINT *, '   *** FAIL: expected IER=0 and WARN_DOY=.TRUE. ***'
         NFAIL = NFAIL + 1
      END IF
      PRINT *, ''

      IF (NFAIL .EQ. 0) THEN
         PRINT *, '=== ALL TESTS PASSED ==='
      ELSE
         PRINT *, '=== ', NFAIL, ' TEST(S) FAILED ==='
         STOP 1
      END IF

      STOP
      END
