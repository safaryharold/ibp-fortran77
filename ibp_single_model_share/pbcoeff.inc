C-----------------------------------------------------------------
C     PBCOEF common block: IBP model coefficients.
C     Populated at run time by READCF (see pbutils.f) instead of
C     being hardcoded into the source, as in the original F90
C     version. INCLUDE this file in every routine that needs the
C     coefficients.
C-----------------------------------------------------------------
C     TAU = mean bubble lifetime 1/gamma (hours); decay rate = 1/TAU.
      DOUBLE PRECISION C1, C2, TAU, MU0, SIGMA
      DOUBLE PRECISION G0, GS(2), GC(2)
      DOUBLE PRECISION T0MON(12), TSMON(2,12), TCMON(2,12)
      DOUBLE PRECISION PHITAB(360,12)
      INTEGER NLONBN
      DOUBLE PRECISION LONST, LONSTP
      INTEGER IBP_INITED

      COMMON /PBCOEF/ C1, C2, TAU, MU0, SIGMA,
     +                G0, GS, GC,
     +                T0MON, TSMON, TCMON,
     +                LONST, LONSTP,
     +                PHITAB, NLONBN, IBP_INITED
