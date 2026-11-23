      PROGRAM IBPCMP
C     Component-level dump for Figure 8 (make_figure8.py).
C     Reads lines "DOY LON F107 LT X" from stdin. For each line it
C     writes the month and the Fortran values of g_osc(doy),
C     ts_month(lon), phi_month(lon), lambda, mu, the final probability
C     P (from IBP_PROBABILITY) and ERF77(X), in full precision.
C     Build (from validation/):
C       gfortran -std=legacy -O2 ibp_components.f \
C         -I../ibp_single_model_share \
C         -L../ibp_single_model_share -libp -o ibp_components
      IMPLICIT NONE
      INCLUDE 'pbcoeff.inc'
      INTEGER IDOY, IER, IOS, MON
      DOUBLE PRECISION ALON, F107, T, X, PROB, G, TS, PHI, LAM, XMU
      LOGICAL W1, W2, W3
      DOUBLE PRECISION GOSC, TSFUNC, PHIMON, ERF77
      INTEGER MONOF
      EXTERNAL GOSC, TSFUNC, PHIMON, ERF77, MONOF
      CALL IBP_INIT('ibp_coefficients.txt', IER)
      IF (IER .NE. 0) STOP 2
   10 READ(*,*,IOSTAT=IOS) IDOY, ALON, F107, T, X
      IF (IOS .NE. 0) GOTO 20
      MON = MONOF(IDOY)
      G   = GOSC(IDOY)
      TS  = TSFUNC(ALON, MON)
      PHI = PHIMON(ALON, MON)
      LAM = (G + C1 + C2*F107) * PHI
      IF (LAM .LT. 0.0D0) LAM = 0.0D0
      XMU = MU0 + TS
      CALL IBP_PROBABILITY(IDOY, ALON, F107, T, PROB, IER, W1, W2, W3)
      WRITE(*,'(I3,7ES25.16E3)') MON, G, TS, PHI, LAM, XMU, PROB,
     +      ERF77(X)
      GOTO 10
   20 END
