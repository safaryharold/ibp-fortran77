      PROGRAM IBPDMP
C     Reads "DOY LON F107 LT" lines from stdin and writes the inputs and
C     the full-precision IBP_PROBABILITY result. Used by
C     compare_official.py.  Build (from validation/):
C       gfortran -std=legacy -O2 ibp_dump.f \
C         -L../ibp_single_model_share -libp -o ibp_dump
      IMPLICIT NONE
      INTEGER IDOY, IER, IOS
      DOUBLE PRECISION ALON, F107, T, PROB
      LOGICAL W1, W2, W3
      CALL IBP_INIT('ibp_coefficients.txt', IER)
      IF (IER .NE. 0) STOP 2
   10 READ(*,*,IOSTAT=IOS) IDOY, ALON, F107, T
      IF (IOS .NE. 0) GOTO 20
      CALL IBP_PROBABILITY(IDOY, ALON, F107, T, PROB, IER, W1, W2, W3)
      WRITE(*,'(I4,3F14.6,ES24.16,I3,3L2)') IDOY,ALON,F107,T,PROB,IER,
     +      W1,W2,W3
      GOTO 10
   20 END
