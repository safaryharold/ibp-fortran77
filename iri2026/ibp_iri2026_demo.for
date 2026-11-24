C     IBP_IRI2026_DEMO.FOR
C     Minimal call of IRI-2026 (IRI_SUB) returning the IBP plasma
C     bubble probability OARR(92). JF switches match the iritest.for
C     defaults; JF(38)=.TRUE. turns the IBP model on.
C     Pass VARIABLES (not constants) to IRI_SUB: it may write to them.
C     Build (in the IRI-2026 folder, with apf107.dat, ig_rz.dat and
C     ibp_emp_coeffs.dat present):
C       gfortran -std=legacy -w -O1 -o ibpdemo ibp_iri2026_demo.for
C         irisub.for irifun.for iritec.for iridreg.for igrf.for
C         cira.for iriflip.for rocdrift.for
C     Expected: Bubble probability = 0.5507 (21 Mar 2024, 20 LT,
C     lat 0, lon -76, 400 km).
      LOGICAL JF(50)
      REAL OUTF(20,1000), OARR(100)
      INTEGER I,JMAG,IYYYY,MMDD
      REAL ALAT,ALON,DHOUR,HBEG,HEND,HSTP
      CALL READ_IG_RZ
      CALL READAPF107
      DO I=1,50
         JF(I)=.TRUE.
      ENDDO
      JF(4)=.FALSE.
      JF(5)=.FALSE.
      JF(6)=.FALSE.
      JF(23)=.FALSE.
      JF(30)=.FALSE.
      JF(33)=.FALSE.
      JF(35)=.FALSE.
      JF(39)=.FALSE.
      JF(40)=.FALSE.
      JF(47)=.FALSE.
      DO I=1,100
         OARR(I)=-1.0
      ENDDO
      JMAG=0
      ALAT=0.0
      ALON=-76.0
      IYYYY=2024
      MMDD=321
      DHOUR=20.0
      HBEG=400.0
      HEND=400.0
      HSTP=1.0
      CALL IRI_SUB(JF,JMAG,ALAT,ALON,IYYYY,MMDD,DHOUR,
     &             HBEG,HEND,HSTP,OUTF,OARR)
      PRINT *, 'Bubble probability =', OARR(92)
      END
