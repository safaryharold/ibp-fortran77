C=================================================================
C     PBSINGLE.F
C
C     Interactive single-run driver of the IBP model in fixed-form
C     FORTRAN 77 (the batch driver is pbbatch.f). Shares
C     pbutils.f / pbcoeff.inc with pbbatch.f. Uses INCLUDE
C     'pbcoeff.inc' for the coefficients, the external function MONOF
C     for DOY->month, DMODLO for MODULO and ERF77 for ERF.
C
C     NOTE: the y0 overflow term is clamped before taking EXP(y0),
C     exactly as pbbatch.f and ibp_probability.f do it (for |y0| >
C     100 the term is clamped and a warning printed). All three
C     drivers now share the identical protective clamp, so they
C     agree even for extreme synthetic inputs.
C
C     Coefficients are loaded at run time from a plain-text data
C     file (ibp_coefficients.txt) via READCF, instead of being
C     hardcoded into the source.
C=================================================================
      PROGRAM PBSINGLE
      IMPLICIT NONE
      INCLUDE 'pbcoeff.inc'

      CHARACTER*40 CFFILE

      INTEGER IDOY, IMONTH, IER
      DOUBLE PRECISION TIN, ALONIN, F107IN

      DOUBLE PRECISION TTRANS, GOSCV, TSVAL, PHIVAL
      DOUBLE PRECISION LAMBD1, LAMBD2, XMU, SIGM, TINTRN, PROB
      DOUBLE PRECISION OMULT, XIMULT, X0, Y0
      DOUBLE PRECISION IMSQRT, OMULT2
      DOUBLE PRECISION ERFARG, XINTEG, EXPARG
      DOUBLE PRECISION PI

      DOUBLE PRECISION GOSC, TSFUNC, PHIMON, ERF77, DMODLO
      INTEGER MONOF
      EXTERNAL GOSC, TSFUNC, PHIMON, ERF77, DMODLO
      EXTERNAL MONOF, READCF

      DATA CFFILE /'ibp_coefficients.txt'/

      PI = ACOS(-1.0D0)

C     --- 1. Get Input Parameters from User ---
      PRINT *, '--- IBP Model Probability Calculation (F77) ---'
      PRINT *, 'Please enter the following input parameters:'

      PRINT *, 'Enter Day of Year (DOY, 1-365):'
      READ (*,*) IDOY

      PRINT *, 'Enter Longitude (degrees East, e.g., 0.0, -75.0, ',
     +         '105.0):'
      READ (*,*) ALONIN

      PRINT *, 'Enter F10.7 solar flux index (e.g., 100.0):'
      READ (*,*) F107IN

      PRINT *, 'Enter Local Time (hours relative to midnight, ',
     +         'e.g., 22.0 for 10 PM):'
      READ (*,*) TIN

C     --- Derive month from DOY ---
      IMONTH = MONOF(IDOY)
      IF (IMONTH .LT. 1 .OR. IMONTH .GT. 12) THEN
         PRINT *, 'Error: Could not determine a valid month ',
     +            'from DOY:', IDOY
         PRINT *, 'Derived month was:', IMONTH,
     +            '. Setting DOY/month to January and continuing.'
         IMONTH = 1
         IDOY = 1
      END IF

      PRINT *, '--- Input Parameters Received ---'
      PRINT *, 'Local Time (t) [Input]:   ', TIN,
     +         ' hours (rel. to midnight)'
      PRINT *, 'Day of Year (doy):        ', IDOY
      PRINT *, 'Derived Month:            ', IMONTH
      PRINT *, 'Longitude (lon):          ', ALONIN,
     +         ' degrees East'
      PRINT *, 'F10.7:                    ', F107IN, ' s.f.u.'
      PRINT *, 'Coefficients:             loaded from data file'
      PRINT *, '-------------------------------------------------',
     +         '----------------------'

C     --- Check F10.7 range ---
      IF (F107IN .LT. 80.0D0 .OR. F107IN .GT. 200.0D0) THEN
         PRINT *, 'Warning: F10.7 =', F107IN,
     +            ' is outside the recommended range (80-200 sfu).'
      END IF

C     --- 2. Load Coefficients ---
      CALL READCF(CFFILE, IER)
      IF (IER .NE. 0) THEN
         PRINT *, 'ERROR: could not read coefficient file: ', CFFILE
         STOP 1
      END IF
      PRINT *, '--- Coefficient Sanity Check ---'
      PRINT *, ' mean lifetime tau=1/gamma (h):', TAU, ' sigma:', SIGMA,
     +         ' mu0:', MU0
      PRINT *, ' C1:', C1, ' C2:', C2
      PRINT *, ' g0:', G0, ' gs:', GS, ' gc:', GC
      PRINT *, '--------------------------------'

C     --- 3. Perform Model Calculations ---
      PRINT *, '--- Intermediate Calculations ---'

C     === Transform Local Time ===
      TTRANS = DMODLO(TIN + 12.0D0, 24.0D0) - 12.0D0
      PRINT *, 'DEBUG: Transformed Local Time (t_internal):',
     +         TTRANS

C     Calculate g_osc(doy)
      GOSCV = GOSC(IDOY)
      PRINT *, 'DEBUG: Seasonal Oscillation g_osc(doy):', GOSCV

C     Calculate phi_month(lon, month)
      PHIVAL = PHIMON(ALONIN, IMONTH)
      PRINT *, 'DEBUG: Longitudinal Density phi_month(lon,',
     +         ' month=', IMONTH, '):', PHIVAL

C     Calculate lambda
      LAMBD1 = (GOSCV + C1 + C2*F107IN) * PHIVAL
      PRINT *, 'DEBUG: Raw lambda_val before check:', LAMBD1
      IF (LAMBD1 .LT. 0.0D0) LAMBD1 = 0.0D0
      PRINT *, 'DEBUG: Final Intensity lambda (lambda_1):', LAMBD1

C     Calculate ts_month(lon, month)
      TSVAL = TSFUNC(ALONIN, IMONTH)
      PRINT *, 'DEBUG: Monthly Timeshift ts_month(lon,',
     +         ' month=', IMONTH, '):', TSVAL

C     Calculate mu = mu0 + ts_month(lon, month)
      XMU = MU0 + TSVAL
      PRINT *, 'DEBUG: Effective Mean Time mu (mu_val):', XMU

C     === Assign parameters for analytical calculation ===
      SIGM   = SIGMA
      TINTRN = TTRANS

C     --- Calculate lambda_2 = decay rate 1/tau (lifetime) ---
      IF (TAU .LE. 1.0D-12) THEN
         PRINT *, 'Error: lifetime TAU (=1/gamma) is near zero ',
     +            'or negative:', TAU
         STOP 1
      END IF
      LAMBD2 = 1.0D0 / TAU

      PRINT *, '--- Analytical Integration Parameters ---'
      PRINT *, 'DEBUG: lambda_1:', LAMBD1
      PRINT *, 'DEBUG: lambda_2 (decay rate):', LAMBD2
      PRINT *, 'DEBUG: mu:', XMU
      PRINT *, 'DEBUG: sigma:', SIGM
      PRINT *, 'DEBUG: time_internal:', TINTRN
      PRINT *, '-----------------------------------------'

C     --- 4. Calculate Final Probability using Analytical Solution ---
      PROB = 0.0D0
      IF (SIGM .LE. 1.0D-12) THEN
         PRINT *, 'Error: Sigma is non-positive or too small (',
     +            SIGM, '). Probability is zero.'
         PROB = 0.0D0
      ELSE IF (LAMBD1 .LE. 1.0D-12) THEN
         PRINT *, 'Info: Lambda_1 is effectively zero. ',
     +            'Probability is zero.'
         PROB = 0.0D0
      ELSE
         OMULT  = -1.0D0 / SQRT(2.0D0 * PI * SIGM**2)
         XIMULT = -1.0D0 / (2.0D0 * SIGM**2)

         IF (XIMULT .GE. -1.0D-12) THEN
            PRINT *, 'Error: inner_mult is non-negative or too ',
     +               'small. Cannot proceed. P=0.'
            PROB = 0.0D0
         ELSE
            X0 = XMU - LAMBD2 / (2.0D0 * XIMULT)
            Y0 = LAMBD2 * (XMU - TINTRN - LAMBD2/(4.0D0*XIMULT))

            IMSQRT = SQRT(-XIMULT)

C           Clamp the overflow term before EXP(y0), matching pbbatch.f
C           and ibp_probability.f so all three drivers behave the same.
            IF (Y0 .GT. 100.0D0) THEN
               PRINT *, 'Warning: Exponent y0 is very large (',
     +                  Y0, '). Clamping to avoid EXP(y0) overflow.'
               OMULT2 = -1.0D30
            ELSE IF (Y0 .LT. -100.0D0) THEN
               PRINT *, 'Warning: Exponent y0 is very small (',
     +                  Y0, '). Clamping OMULT2 to zero.'
               OMULT2 = 0.0D0
            ELSE
               OMULT2 = SQRT(PI) * OMULT * EXP(Y0) / IMSQRT
            END IF

            ERFARG = (TINTRN - X0) * IMSQRT
            XINTEG = (0.5D0 + ERF77(ERFARG)*0.5D0) * OMULT2
            PRINT *, "DEBUG: Analytical 'integrated' value:",
     +               XINTEG

            EXPARG = LAMBD1 * XINTEG
            PRINT *, 'DEBUG: Exponent argument ',
     +               '(lambda_1 * integrated):', EXPARG

            IF (EXPARG .LT. -700.0D0) THEN
               PROB = 1.0D0
               PRINT *, 'DEBUG: Exponent argument is very small; ',
     +                  'P approximates 1.'
            ELSE IF (EXPARG .GT. 1.0D-9) THEN
               PRINT *, 'Error: Exponent is positive (', EXPARG,
     +                  '). Expected <= 0. Setting P=0.'
               PROB = 0.0D0
            ELSE
               PROB = 1.0D0 - EXP(EXPARG)
            END IF
         END IF
      END IF

      PROB = MAX(0.0D0, MIN(1.0D0, PROB))

      PRINT *, '--- Final Result ---'
      PRINT *, 'Occurrence Probability P[X(t)=1]: ', PROB
      PRINT *, '========================================',
     +         '==============================='

      STOP
      END
