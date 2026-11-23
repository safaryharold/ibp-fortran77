C=================================================================
C     PBBATCH.F
C
C     Batch driver of the IBP model in fixed-form FORTRAN 77. Reads
C     DOY, longitude, F10.7 and local time from batch_inputs.txt and
C     writes probabilities to loopresults.txt. FORTRAN 77 has no
C     modules, ADJUSTL/TRIM, MODULO or ERF, so the code uses INCLUDE
C     'pbcoeff.inc', FSTNB/LSTNB, DMODLO and ERF77 (pbutils.f), and
C     labelled DO loops with GOTO.
C
C     Coefficients are loaded at run time from a plain-text data
C     file (ibp_coefficients.txt) via READCF.
C=================================================================
      PROGRAM PBBATCH
      IMPLICIT NONE
      INCLUDE 'pbcoeff.inc'

      CHARACTER*40 INFILE, OUTFILE, CFFILE
      INTEGER IUIN, IUOUT
      INTEGER IOS, NRUN

      INTEGER IDOY, IMONTH
      DOUBLE PRECISION TIN, ALONIN, F107IN

      DOUBLE PRECISION TTRANS, GOSCV, TSVAL, PHIVAL
      DOUBLE PRECISION LAMBD1, LAMBD2, XMU, SIGM, TINTRN, PROB
      DOUBLE PRECISION OMULT, XIMULT, X0, Y0
      DOUBLE PRECISION IMSQRT, OMULT2
      DOUBLE PRECISION ERFARG, XINTEG, EXPARG
      DOUBLE PRECISION PI

C     -- header-line detection --
      CHARACTER*256 LINE
      CHARACTER*1 FCHAR
      INTEGER IPEEK, IP, LNB

      DOUBLE PRECISION GOSC, TSFUNC, PHIMON, ERF77, DMODLO
      INTEGER MONOF, FSTNB, LSTNB
      EXTERNAL GOSC, TSFUNC, PHIMON, ERF77, DMODLO
      EXTERNAL MONOF, FSTNB, LSTNB, READCF

      DATA INFILE  /'batch_inputs.txt'/
      DATA OUTFILE /'loopresults.txt'/
      DATA CFFILE  /'ibp_coefficients.txt'/

      IUIN  = 10
      IUOUT = 11
      PI = ACOS(-1.0D0)

      PRINT *, '--- IBP Model Batch Probability Calculation ---'
      PRINT *, '    (FORTRAN 77 version)'
      PRINT *, 'Loading coefficients from: ', CFFILE

      CALL READCF(CFFILE)

      PRINT *, 'Coefficients loaded successfully.'
      PRINT *, 'Input file:  ', INFILE
      PRINT *, 'Output file: ', OUTFILE

      OPEN(UNIT=IUIN, FILE=INFILE, STATUS='OLD', ACTION='READ',
     +     IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
         PRINT *, 'Error opening input file: ', INFILE,
     +            ' IOSTAT=', IOS
         STOP 1
      END IF

      OPEN(UNIT=IUOUT, FILE=OUTFILE, STATUS='REPLACE',
     +     ACTION='WRITE', IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
         PRINT *, 'Error opening output file: ', OUTFILE,
     +            ' IOSTAT=', IOS
         CLOSE(IUIN)
         STOP 2
      END IF

C     --- Skip header line in input file if present ---
C     Peek at the first non-blank character of the first line. If
C     it looks like the start of a number (digit or minus sign),
C     rewind so the main loop reads it as the first data row.
C     Otherwise leave it consumed (it was a header) and report it.
      READ(IUIN,'(A)',IOSTAT=IPEEK) LINE
      IF (IPEEK .EQ. 0) THEN
         IP = FSTNB(LINE)
         IF (IP .EQ. 0) THEN
            FCHAR = ' '
         ELSE
            FCHAR = LINE(IP:IP)
         END IF
         IF ((FCHAR .GE. '0' .AND. FCHAR .LE. '9') .OR.
     +       FCHAR .EQ. '-') THEN
            REWIND(IUIN)
         ELSE
            LNB = LSTNB(LINE)
            IF (LNB .EQ. 0) THEN
               PRINT *, 'Header line detected and skipped: '
            ELSE
               PRINT *, 'Header line detected and skipped: ',
     +                  LINE(1:LNB)
            END IF
         END IF
      END IF

      WRITE(IUOUT,'(A)')
     +   'DOY  Lon_deg  F10.7   LT_hour  Derived_Month  Probability'
      WRITE(IUOUT,'(A)')
     +   '--- -------- ------ --------- -------------- -----------'

      NRUN = 0

  100 CONTINUE
      READ(IUIN,*,IOSTAT=IOS) IDOY, ALONIN, F107IN, TIN
      IF (IOS .LT. 0) THEN
         PRINT *, 'End of input file reached.'
         GOTO 200
      ELSE IF (IOS .GT. 0) THEN
         PRINT *, 'Error reading from input file. IOSTAT=', IOS
         GOTO 200
      END IF

      NRUN = NRUN + 1
      PRINT *, 'Processing run: ', NRUN, ' with DOY=', IDOY

      IMONTH = MONOF(IDOY)
      IF (IMONTH .LT. 1 .OR. IMONTH .GT. 12) THEN
         WRITE(IUOUT,'(I3,A,F8.2,A,F6.1,A,F7.2,A)')
     +      IDOY, ' ', ALONIN, ' ', F107IN, ' ', TIN,
     +      ' ERROR: Invalid DOY for month derivation'
         PRINT *, 'Error: Could not determine a valid month ',
     +            'from DOY:', IDOY, ' for run ', NRUN,
     +            '. Skipping.'
         GOTO 100
      END IF

      IF (F107IN .LT. 80.0D0 .OR. F107IN .GT. 200.0D0) THEN
         PRINT *, 'Warning for run ', NRUN, ': F10.7=', F107IN,
     +            ' is outside recommended range (80-200 sfu).'
      END IF

      TTRANS = DMODLO(TIN + 12.0D0, 24.0D0) - 12.0D0
      GOSCV  = GOSC(IDOY)
      PHIVAL = PHIMON(ALONIN, IMONTH)
      LAMBD1 = (GOSCV + C1 + C2*F107IN) * PHIVAL
      IF (LAMBD1 .LT. 0.0D0) LAMBD1 = 0.0D0

      TSVAL  = TSFUNC(ALONIN, IMONTH)
      XMU    = MU0 + TSVAL
      SIGM   = SIGMA
      TINTRN = TTRANS

      IF (TAU .LE. 1.0D-12) THEN
         WRITE(IUOUT,'(I3,A,F8.2,A,F6.1,A,F7.2,A,I2,A)')
     +      IDOY, ' ', ALONIN, ' ', F107IN, ' ', TIN,
     +      ' ', IMONTH, ' ERROR: lifetime TAU too small'
         PRINT *, 'Error for run ', NRUN, ': lifetime TAU ',
     +            '(lifetime) is near zero or negative. Skipping.'
         GOTO 100
      END IF
      LAMBD2 = 1.0D0 / TAU

      PROB = 0.0D0
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

            IF (Y0 .GT. 100.0D0) THEN
               OMULT2 = -1.0D30
            ELSE IF (Y0 .LT. -100.0D0) THEN
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

      WRITE(IUOUT,'(I3,F9.2,F8.1,F8.2,I15,F12.7)')
     +   IDOY, ALONIN, F107IN, TIN, IMONTH, PROB

      GOTO 100

  200 CONTINUE
      CLOSE(IUIN)
      CLOSE(IUOUT)
      PRINT *, 'Batch processing complete. Processed', NRUN, 'runs.'
      PRINT *, 'Results written to: ', OUTFILE

      STOP
      END
