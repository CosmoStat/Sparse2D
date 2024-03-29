*DECK DRC3JJ
      SUBROUTINE WIG3J (L2, L3, M2, M3, L1MIN, L1MAX, THRCOF, NDIM, IER)
C***BEGIN PROLOGUE  WIG3J
C***PURPOSE  Evaluate the 3j symbol f(L1) = (  L1   L2 L3)
C                                           (-M2-M3 M2 M3)
C            for all allowed values of L1, the other parameters
C            being held fixed.
C***LIBRARY   SLATEC
C***CATEGORY  C19
C***TYPE      DOUBLE PRECISION (RC3JJ-S, WIG3J-D)
C***KEYWORDS  3J COEFFICIENTS, 3J SYMBOLS, CLEBSCH-GORDAN COEFFICIENTS,
C             RACAH COEFFICIENTS, VECTOR ADDITION COEFFICIENTS,
C             WIGNER COEFFICIENTS
C***AUTHOR  Gordon, R. G., Harvard University
C           Schulten, K., Max Planck Institute
C***DESCRIPTION
C
C *Usage:
C
C        DOUBLE PRECISION L2, L3, M2, M3, L1MIN, L1MAX, THRCOF(NDIM)
C        INTEGER NDIM, IER
C
C        CALL WIG3J (L2, L3, M2, M3, L1MIN, L1MAX, THRCOF, NDIM, IER)
C
C *Arguments:
C
C     L2 :IN      Parameter in 3j symbol.
C
C     L3 :IN      Parameter in 3j symbol.
C
C     M2 :IN      Parameter in 3j symbol.
C
C     M3 :IN      Parameter in 3j symbol.
C
C     L1MIN :OUT  Smallest allowable L1 in 3j symbol.
C
C     L1MAX :OUT  Largest allowable L1 in 3j symbol.
C
C     THRCOF :OUT Set of 3j coefficients generated by evaluating the
C                 3j symbol for all allowed values of L1.  THRCOF(I)
C                 will contain f(L1MIN+I-1), I=1,2,...,L1MAX+L1MIN+1.
C
C     NDIM :IN    Declared length of THRCOF in calling program.
C
C     IER :OUT    Error flag.
C                 IER=0 No errors.
C                 IER=1 Either L2.LT.ABS(M2) or L3.LT.ABS(M3).
C                 IER=2 Either L2+ABS(M2) or L3+ABS(M3) non-integer.
C                 IER=3 L1MAX-L1MIN not an integer.
C                 IER=4 L1MAX less than L1MIN.
C                 IER=5 NDIM less than L1MAX-L1MIN+1.
C
C *Description:
C
C     Although conventionally the parameters of the vector addition
C  coefficients satisfy certain restrictions, such as being integers
C  or integers plus 1/2, the restrictions imposed on input to this
C  subroutine are somewhat weaker. See, for example, Section 27.9 of
C  Abramowitz and Stegun or Appendix C of Volume II of A. Messiah.
C  The restrictions imposed by this subroutine are
C       1. L2 .GE. ABS(M2) and L3 .GE. ABS(M3);
C       2. L2+ABS(M2) and L3+ABS(M3) must be integers;
C       3. L1MAX-L1MIN must be a non-negative integer, where
C          L1MAX=L2+L3 and L1MIN=MAX(ABS(L2-L3),ABS(M2+M3)).
C  If the conventional restrictions are satisfied, then these
C  restrictions are met.
C
C     The user should be cautious in using input parameters that do
C  not satisfy the conventional restrictions. For example, the
C  the subroutine produces values of
C       f(L1) = ( L1  2.5  5.8)
C               (-0.3 1.5 -1.2)
C  for L1=3.3,4.3,...,8.3 but none of the symmetry properties of the 3j
C  symbol, set forth on page 1056 of Messiah, is satisfied.
C
C     The subroutine generates f(L1MIN), f(L1MIN+1), ..., f(L1MAX)
C  where L1MIN and L1MAX are defined above. The sequence f(L1) is
C  generated by a three-term recurrence algorithm with scaling to
C  control overflow. Both backward and forward recurrence are used to
C  maintain numerical stability. The two recurrence sequences are
C  matched at an interior point and are normalized from the unitary
C  property of 3j coefficients and Wigner's phase convention.
C
C    The algorithm is suited to applications in which large quantum
C  numbers arise, such as in molecular dynamics.
C
C***REFERENCES  1. Abramowitz, M., and Stegun, I. A., Eds., Handbook
C                  of Mathematical Functions with Formulas, Graphs
C                  and Mathematical Tables, NBS Applied Mathematics
C                  Series 55, June 1964 and subsequent printings.
C               2. Messiah, Albert., Quantum Mechanics, Volume II,
C                  North-Holland Publishing Company, 1963.
C               3. Schulten, Klaus and Gordon, Roy G., Exact recursive
C                  evaluation of 3j and 6j coefficients for quantum-
C                  mechanical coupling of angular momenta, J Math
C                  Phys, v 16, no. 10, October 1975, pp. 1961-1970.
C               4. Schulten, Klaus and Gordon, Roy G., Semiclassical
C                  approximations to 3j  and 6j coefficients for
C                  quantum-mechanical coupling of angular momenta,
C                  J Math Phys, v 16, no. 10, October 1975,
C                  pp. 1971-1988.
C               5. Schulten, Klaus and Gordon, Roy G., Recursive
C                  evaluation of 3j and 6j coefficients, Computer
C                  Phys Comm, v 11, 1976, pp. 269-278.
C***ROUTINES CALLED 
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   880515  SLATEC prologue added by G. C. Nielson, NBS; parameters
C           HUGE and TINY revised to depend on D1MACH.
C   891229  Prologue description rewritten; other prologue sections
C           revised; LMATCH (location of match point for recurrences)
C           removed from argument list; argument IER changed to serve
C           only as an error flag (previously, in cases without error,
C           it returned the number of scalings); number of error codes
C           increased to provide more precise error information;
C           program comments revised; SLATEC error handler calls
C           introduced to enable printing of error messages to meet
C           SLATEC standards. These changes were done by D. W. Lozier,
C           M. A. McClain and J. M. Smith of the National Institute
C           of Standards and Technology, formerly NBS.
C   910415  Mixed type expressions eliminated; variable C1 initialized;
C           description of THRCOF expanded. These changes were done by
C           D. W. Lozier.
C***END PROLOGUE  WIG3J
C
      INTEGER NDIM, IER
      DOUBLE PRECISION L2, L3, M2, M3, L1MIN, L1MAX, THRCOF(NDIM)
C
      INTEGER I, INDEX, LSTEP, N, NFIN, NFINP1, NFINP2, NFINP3, NLIM,
     +        NSTEP2
      DOUBLE PRECISION A1, A1S, A2, A2S, C1, C1OLD, C2, CNORM, D1MACH,
     +                 DENOM, DV, EPS, HUGE, L1, M1, NEWFAC, OLDFAC,
     +                 ONE, RATIO, SIGN1, SIGN2, SRHUGE, SRTINY, SUM1,
     +                 SUM2, SUMBAC, SUMFOR, SUMUNI, THREE, THRESH,
     +                 TINY, TWO, X, X1, X2, X3, Y, Y1, Y2, Y3, ZERO
C
      DATA  ZERO,EPS,ONE,TWO,THREE /0.0D0,0.01D0,1.0D0,2.0D0,3.0D0/
C
C***FIRST EXECUTABLE STATEMENT  WIG3J
      IER=0
C  HUGE is the square root of one twentieth of the largest floating
C  point number, approximately.
c$$$      HUGE = SQRT(D1MACH(2)/20.0D0)
      HUGE = SQRT((10D0**307)/20.0D0)
      SRHUGE = SQRT(HUGE)
      TINY = 1.0D0/HUGE
      SRTINY = 1.0D0/SRHUGE
C     
C     LMATCH = ZERO
      M1 = - M2 - M3
C
C  Check error conditions 1 and 2.
      IF((L2-ABS(M2)+EPS.LT.ZERO).OR.
     +   (L3-ABS(M3)+EPS.LT.ZERO))THEN
         IER=1
         PRINT*, 'L2-ABS(M2) or L3-ABS(M3) less than zero.'
         RETURN
      ELSEIF((MOD(L2+ABS(M2)+EPS,ONE).GE.EPS+EPS).OR.
     +   (MOD(L3+ABS(M3)+EPS,ONE).GE.EPS+EPS))THEN
         IER=2
         PRINT*, 'L2+ABS(M2) or L3+ABS(M3) not integer.'
         RETURN
      ENDIF
C
C
C
C  Limits for L1
C
      L1MIN = MAX(ABS(L2-L3),ABS(M1))
      L1MAX = L2 + L3
C
C  Check error condition 3.
      IF(MOD(L1MAX-L1MIN+EPS,ONE).GE.EPS+EPS)THEN
         IER=3
         PRINT*, 'L1MAX-L1MIN not integer.'
         RETURN 
      ENDIF
      IF(L1MIN.LT.L1MAX-EPS)   GO TO 20
      IF(L1MIN.LT.L1MAX+EPS)   GO TO 10
C
C  Check error condition 4.
      IER=4
      PRINT*, 'L1MIN greater than L1MAX.'
      RETURN
C
C  This is reached in case that L1 can take only one value,
C  i.e. L1MIN = L1MAX
C
   10 CONTINUE
C     LSCALE = 0
      THRCOF(1) = (-ONE) ** INT(ABS(L2+M2-L3+M3)+EPS) /
     1 SQRT(L1MIN + L2 + L3 + ONE)
      RETURN
C
C  This is reached in case that L1 takes more than one value,
C  i.e. L1MIN < L1MAX.
C
   20 CONTINUE
C     LSCALE = 0
      NFIN = INT(L1MAX-L1MIN+ONE+EPS)
      IF(NDIM-NFIN)  21, 23, 23
C
C  Check error condition 5.
   21 IER = 5
      PRINT*, 'Dimension of result array for 3j coefficients too small.'
      RETURN
C
C
C  Starting forward recursion from L1MIN taking NSTEP1 steps
C
   23 L1 = L1MIN
      NEWFAC = 0.0D0
      C1 = 0.0D0
      THRCOF(1) = SRTINY
      SUM1 = (L1+L1+ONE) * TINY
C
C
      LSTEP = 1
   30 LSTEP = LSTEP + 1
      L1 = L1 + ONE
C
C
      OLDFAC = NEWFAC
      A1 = (L1+L2+L3+ONE) * (L1-L2+L3) * (L1+L2-L3) * (-L1+L2+L3+ONE)
      A2 = (L1+M1) * (L1-M1)
      NEWFAC = SQRT(A1*A2)
      IF(L1.LT.ONE+EPS)   GO TO 40
C
C
      DV = - L2*(L2+ONE) * M1 + L3*(L3+ONE) * M1 + L1*(L1-ONE) * (M3-M2)
      DENOM = (L1-ONE) * NEWFAC
C
      IF(LSTEP-2)  32, 32, 31
C
   31 C1OLD = ABS(C1)
   32 C1 = - (L1+L1-ONE) * DV / DENOM
      GO TO 50
C
C  If L1 = 1, (L1-1) has to be factored out of DV, hence
C
   40 C1 = - (L1+L1-ONE) * L1 * (M3-M2) / NEWFAC
C
   50 IF(LSTEP.GT.2)   GO TO 60
C
C
C  If L1 = L1MIN + 1, the third term in the recursion equation vanishes,
C  hence
      X = SRTINY * C1
      THRCOF(2) = X
      SUM1 = SUM1 + TINY * (L1+L1+ONE) * C1*C1
      IF(LSTEP.EQ.NFIN)   GO TO 220
      GO TO 30
C
C
   60 C2 = - L1 * OLDFAC / DENOM
C
C  Recursion to the next 3j coefficient X
C
      X = C1 * THRCOF(LSTEP-1) + C2 * THRCOF(LSTEP-2)
      THRCOF(LSTEP) = X
      SUMFOR = SUM1
      SUM1 = SUM1 + (L1+L1+ONE) * X*X
      IF(LSTEP.EQ.NFIN)   GO TO 100
C
C  See if last unnormalized 3j coefficient exceeds SRHUGE
C
      IF(ABS(X).LT.SRHUGE)   GO TO 80
C
C  This is reached if last 3j coefficient larger than SRHUGE,
C  so that the recursion series THRCOF(1), ... , THRCOF(LSTEP)
C  has to be rescaled to prevent overflow
C
C     LSCALE = LSCALE + 1
      DO 70 I=1,LSTEP
      IF(ABS(THRCOF(I)).LT.SRTINY)   THRCOF(I) = ZERO
   70 THRCOF(I) = THRCOF(I) / SRHUGE
      SUM1 = SUM1 / HUGE
      SUMFOR = SUMFOR / HUGE
      X = X / SRHUGE
C
C  As long as ABS(C1) is decreasing, the recursion proceeds towards
C  increasing 3j values and, hence, is numerically stable.  Once
C  an increase of ABS(C1) is detected, the recursion direction is
C  reversed.
C
   80 IF(C1OLD-ABS(C1))   100, 100, 30
C
C
C  Keep three 3j coefficients around LMATCH for comparison with
C  backward recursion.
C
  100 CONTINUE
C     LMATCH = L1 - 1
      X1 = X
      X2 = THRCOF(LSTEP-1)
      X3 = THRCOF(LSTEP-2)
      NSTEP2 = NFIN - LSTEP + 3
C
C
C
C
C  Starting backward recursion from L1MAX taking NSTEP2 steps, so
C  that forward and backward recursion overlap at three points
C  L1 = LMATCH+1, LMATCH, LMATCH-1.
C
      NFINP1 = NFIN + 1
      NFINP2 = NFIN + 2
      NFINP3 = NFIN + 3
      L1 = L1MAX
      THRCOF(NFIN) = SRTINY
      SUM2 = TINY * (L1+L1+ONE)
C
      L1 = L1 + TWO
      LSTEP = 1
  110 LSTEP = LSTEP + 1
      L1 = L1 - ONE
C
      OLDFAC = NEWFAC
      A1S = (L1+L2+L3)*(L1-L2+L3-ONE)*(L1+L2-L3-ONE)*(-L1+L2+L3+TWO)
      A2S = (L1+M1-ONE) * (L1-M1-ONE)
      NEWFAC = SQRT(A1S*A2S)
C
      DV = - L2*(L2+ONE) * M1 + L3*(L3+ONE) * M1 + L1*(L1-ONE) * (M3-M2)
C
      DENOM = L1 * NEWFAC
      C1 = - (L1+L1-ONE) * DV / DENOM
      IF(LSTEP.GT.2)   GO TO 120
C
C  If L1 = L1MAX + 1, the third term in the recursion formula vanishes
C
      Y = SRTINY * C1
      THRCOF(NFIN-1) = Y
      SUMBAC = SUM2
      SUM2 = SUM2 + TINY * (L1+L1-THREE) * C1*C1
C
      GO TO 110
C
C
  120 C2 = - (L1 - ONE) * OLDFAC / DENOM
C
C  Recursion to the next 3j coefficient Y
C
      Y = C1 * THRCOF(NFINP2-LSTEP) + C2 * THRCOF(NFINP3-LSTEP)
C
      IF(LSTEP.EQ.NSTEP2)   GO TO 200
C
      THRCOF(NFINP1-LSTEP) = Y
      SUMBAC = SUM2
      SUM2 = SUM2 + (L1+L1-THREE) * Y*Y
C
C  See if last unnormalized 3j coefficient exceeds SRHUGE
C
      IF(ABS(Y).LT.SRHUGE)   GO TO 110
C
C  This is reached if last 3j coefficient larger than SRHUGE,
C  so that the recursion series THRCOF(NFIN), ... ,THRCOF(NFIN-LSTEP+1)
C  has to be rescaled to prevent overflow
C
C     LSCALE = LSCALE + 1
      DO 130 I=1,LSTEP
      INDEX = NFIN - I + 1
      IF(ABS(THRCOF(INDEX)).LT.SRTINY)   THRCOF(INDEX) = ZERO
  130 THRCOF(INDEX) = THRCOF(INDEX) / SRHUGE
      SUM2 = SUM2 / HUGE
      SUMBAC = SUMBAC / HUGE
C
C
      GO TO 110
C
C
C  The forward recursion 3j coefficients X1, X2, X3 are to be matched
C  with the corresponding backward recursion values Y1, Y2, Y3.
C
  200 Y3 = Y
      Y2 = THRCOF(NFINP2-LSTEP)
      Y1 = THRCOF(NFINP3-LSTEP)
C
C
C  Determine now RATIO such that YI = RATIO * XI  (I=1,2,3) holds
C  with minimal error.
C
      RATIO = ( X1*Y1 + X2*Y2 + X3*Y3 ) / ( X1*X1 + X2*X2 + X3*X3 )
      NLIM = NFIN - NSTEP2 + 1
C
      IF(ABS(RATIO).LT.ONE)   GO TO 211
C
      DO 210 N=1,NLIM
  210 THRCOF(N) = RATIO * THRCOF(N)
      SUMUNI = RATIO * RATIO * SUMFOR + SUMBAC
      GO TO 230
C
  211 NLIM = NLIM + 1
      RATIO = ONE / RATIO
      DO 212 N=NLIM,NFIN
  212 THRCOF(N) = RATIO * THRCOF(N)
      SUMUNI = SUMFOR + RATIO*RATIO*SUMBAC
      GO TO 230
C
  220 SUMUNI = SUM1
C
C
C  Normalize 3j coefficients
C
  230 CNORM = ONE / SQRT(SUMUNI)
C
C  Sign convention for last 3j coefficient determines overall phase
C
      SIGN1 = SIGN(ONE,THRCOF(NFIN))
      SIGN2 = (-ONE) ** INT(ABS(L2+M2-L3+M3)+EPS)
      IF(SIGN1*SIGN2) 235,235,236
  235 CNORM = - CNORM
C
  236 IF(ABS(CNORM).LT.ONE)   GO TO 250
C
      DO 240 N=1,NFIN
  240 THRCOF(N) = CNORM * THRCOF(N)
      RETURN
C
  250 THRESH = TINY / ABS(CNORM)
      DO 251 N=1,NFIN
      IF(ABS(THRCOF(N)).LT.THRESH)   THRCOF(N) = ZERO
  251 THRCOF(N) = CNORM * THRCOF(N)
C
      RETURN
      END
