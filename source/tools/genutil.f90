MODULE GENUTIL
  ! generally useful utilities
  USE MT19937 ! mersenne random number generator

  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793D0
  
CONTAINS  

  INTEGER FUNCTION STRING2NUM(STRINGIN,APPENDNUM)
    ! convert the string to a unique number based on ascii characters
    ! the characters SPACE, {,},(,),[,],",`,<,> and all nonprintable characters are ignored
    ! at most the last five characters (ignoring the unacceptable characters above) at the end of the string are used
    ! any leading "!" do not affect the final number (these map to 0)
    ! if APPENDNUM is specificied, only use the last 4 characters of the string as well as the additional number modulo 84

    IMPLICIT NONE
    CHARACTER(LEN=*) :: STRINGIN
    CHARACTER(LEN=5) :: STRING
    INTEGER, OPTIONAL :: APPENDNUM
    INTEGER :: DIGARRAY(5)
    INTEGER :: ALLOWED(84)
    INTEGER :: N, I, D, COUNT
    CHARACTER*84 :: ALLOWEDSTR

    ! set the allowed characters
    ALLOWED(1:6) = (/33,35,36,37,38,39/)
    ALLOWED(7:24) = (/(I,I=42,59)/)
    ALLOWED(25:27) = (/61,63,64/)
    ALLOWED(28:53) = (/(I, I=65,90)/)
    ALLOWED(54:56) = (/92,94,95/)
    ALLOWED(57:82) = (/(I, I=97,122)/)
    ALLOWED(83:84) = (/124,126/)

    N = LEN(STRINGIN)
    IF (PRESENT(APPENDNUM)) THEN
       STRING(1:4) = STRINGIN(N-3:N)
       STRING(5:5) = ACHAR(ALLOWED(MOD(APPENDNUM,84)+1))
    ELSE
       STRING = STRINGIN(N-4:N)
    ENDIF
    N =  5


    DO I = 1,84
       ALLOWEDSTR(I:I) = ACHAR(ALLOWED(I))
    ENDDO

    DIGARRAY = 0
    COUNT = 0
    DO I = 0,N-1
       D = INDEX(ALLOWEDSTR,STRING(N-I:N-I),.FALSE.)
       IF (D.EQ.0) THEN
          print*, 'Ignoring character:', D
          CYCLE
       ENDIF

       DIGARRAY(5-COUNT) = D-1
       COUNT = COUNT + 1
       IF (COUNT.GE.5) EXIT
    ENDDO

    STRING2NUM = BASE2DEC(DIGARRAY,5,84)
  END FUNCTION STRING2NUM

  INTEGER FUNCTION BASE2DEC(DIGARRAY,N,BASE)
  ! given a number in some integer base (specified as a list of digits)
  ! convert that number to a decimal integer
  ! N is the size of the list
  ! if resulting number is too large, wrap around to negative numbers
  ! starting from the right, only use as many of the digits as 
  ! will fit into the resulting integer between -HUGE and HUGE  
  ! if any digit is greater than base-1, print error and stop

  IMPLICIT NONE
  INTEGER, DIMENSION(N) :: DIGARRAY
  INTEGER, INTENT(IN) :: N, BASE
  INTEGER :: MAXDIG, I, D

  MAXDIG = INT(LOG(2*DBLE(HUGE(BASE))+2)/LOG(DBLE(BASE)))

  BASE2DEC = 0
  DO I = 0, MIN(N-1,MAXDIG-1)
     D = DIGARRAY(N-I)
     IF (D.EQ.0) CYCLE
     IF (D.GT.BASE-1) THEN
        PRINT*, 'ERROR in BASE2DEC: digit is bigger than base.', I, D, BASE
        STOP 1
     ENDIF
     
     BASE2DEC = BASE2DEC + D*BASE**I
  ENDDO
  
  END FUNCTION BASE2DEC

  SUBROUTINE INTERPARRAY(ARRAY,NA,COL,VAL,IND,INTERP)
    ! for an 2D array with dimensions NA
    ! use the values in column COL to interpolate for the value VAL
    ! return the index IND such that ARRAY(IND,COL)<VAL<ARRAY(IND+1,COL)
    ! and the interpolation of all other columns in INTERP
    ! If val is out of bounds returns IND=0 or IND=NL
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NA(2), COL
    DOUBLE PRECISION, INTENT(IN) :: ARRAY(NA(1),NA(2))
    DOUBLE PRECISION, INTENT(IN) :: VAL
    INTEGER, INTENT(OUT) :: IND
    DOUBLE PRECISION, INTENT(OUT) :: INTERP(NA(2))
    DOUBLE PRECISION :: FRAC

    CALL INTERP1(ARRAY(:,COL),NA(1),VAL,IND)

    IF (IND.EQ.0.OR.IND.GE.NA(1)) RETURN

    FRAC = (VAL-ARRAY(IND,COL))/(ARRAY(IND+1,COL)-ARRAY(IND,COL))
    INTERP = (1-FRAC)*ARRAY(IND,:)+FRAC*ARRAY(IND+1,:)
       
  END SUBROUTINE INTERPARRAY

  SUBROUTINE INTERP1(LIST,NL,VAL,IND)
    ! for a monotonically increasing, double precision list
    ! find the index where LIST(IND) <VAL<LIST(IND+1)
    INTEGER, INTENT(IN) :: NL
    DOUBLE PRECISION, INTENT(IN) :: LIST(NL)
    DOUBLE PRECISION, INTENT(IN) :: VAL
    INTEGER, INTENT(OUT) :: IND
    DOUBLE PRECISION :: MINL, MAXL,PREVMAXL
    INTEGER :: MINI, MAXI,PREVMAXI
    LOGICAL :: VERBOSE = .FALSE.

    MINI = 1; MAXI = NL; MINL = LIST(1); MAXL = LIST(NL)
    IF (VAL.LT.MINL) THEN
       IND = 0; RETURN
    ELSEIF (VAL.EQ.MINL) THEN
       IND = 1; RETURN       
    ELSEIF (VAL.GT.MAXL) THEN
       IND = NL; RETURN
    ELSEIF (VAL.EQ.MAXL) THEN
       IND = NL-1; RETURN
    ENDIF

    DO WHILE (MAXI-MINI.GT.1.OR.MAXL.LT.VAL)
       IF (MAXL.GT.VAL) THEN
          PREVMAXL = MAXL; PREVMAXI = MAXI
          MAXI = MINI + (MAXI-MINI)/2
          MAXL = LIST(MAXI)
       ELSE
          MINI = MAXI; MAXI = PREVMAXI
          MINL = MAXL; MAXL = PREVMAXL
       ENDIF
       IF (VERBOSE) PRINT*, 'MINI, MAXI, MINL, MAXL', MINI, MAXI, MINL, MAXL,VAL
       if (maxi.eq.mini) then
          print*, 'something weird in interp1:', list(1), list(nl), val          
          stop 1
       endif
    ENDDO

    IF (.NOT.(MAXI.EQ.MINI+1.AND.LIST(MINI).LE.VAL.AND.LIST(MAXI).GE.VAL)) THEN
       PRINT*, 'SOMETHING IS WEIRD IN INTERP1', val, mini, maxi, list(mini), list(maxi)
       STOP 1
    ENDIF

    IND = MINI
  END SUBROUTINE INTERP1

  SUBROUTINE REPLACESUBSTR(INSTRING,C,REPL)
    ! replace the last instance of the substring C in INSTRING with REPL
    IMPLICIT NONE
    CHARACTER*100 :: INSTRING
    CHARACTER(LEN=*) :: C
    CHARACTER(LEN=*) :: REPL
    INTEGER :: LENC, IND

    INSTRING = ADJUSTL(INSTRING)

    LENC = LEN_TRIM(C)

    IND = INDEX(INSTRING,C,.TRUE.)
    IF (IND.GT.0) THEN! if * was found in the string

       INSTRING = INSTRING(1:IND-1) // TRIM(ADJUSTL(REPL)) // INSTRING(IND+LENC:100)   
    END IF
  END SUBROUTINE REPLACESUBSTR

  SUBROUTINE NORMALIZE(X)
    ! normalize a 3 dimensional vector

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT) :: X(3)
    DOUBLE PRECISION :: DX

    DX = SQRT(DOT_PRODUCT(X,X))
    X(:) = X(:)/DX

    RETURN
  END SUBROUTINE NORMALIZE

  DOUBLE PRECISION FUNCTION NORM(X)
    ! norm of 3D vector X

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: X(3)

    NORM = sqrt(DOT_PRODUCT(X,X))

  END FUNCTION NORM

  SUBROUTINE CROSS_PRODUCT(A, B, C)
    ! take the cross product of 3D vectors A and B; return result in C

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: A(3), B(3)
    DOUBLE PRECISION, INTENT(OUT) :: C(3)

    C(1) = A(2)*B(3) - A(3)*B(2)
    C(2) = A(3)*B(1)-A(1)*B(3)
    C(3) = A(1)*B(2) - A(2)*B(1)

    RETURN
  END SUBROUTINE CROSS_PRODUCT

  SUBROUTINE RANDOMAXIS(REFAX,CTRANGE,RANDAX)
    ! generate a random axis, within a certain range in cos(theta) 
    ! relative to the reference axis    
    DOUBLE PRECISION, INTENT(IN) :: REFAX(3), CTRANGE
    DOUBLE PRECISION, INTENT(OUT) :: RANDAX(3)
    DOUBLE PRECISION :: THETA, pHI, E1(3), E2(3), E3(3), X, Y, Z

    THETA = 1.0D0
    DO WHILE (THETA.EQ.1.0D0)
       ! get random number; MT19937 algorithm uses closed interval [0,1], 
       ! so ignore when R is exactly 1
       THETA = GRND() !get a random number
    ENDDO
    !CALL RANDOM_NUMBER(THETA)
    THETA = acos(1D0 - THETA*MAX(CTRANGE,2D0))

    PHI = 1.0D0
    DO WHILE (PHI.EQ.1.0D0)
       PHI = GRND() !get a random number
    ENDDO
    !CALL RANDOM_NUMBER(PHI)
    PHI = PHI*2*PI

    ! axis system relative to which angles are defined
    E3 = REFAX
    IF (E3(2) == 0 .AND. E3(3) == 0) THEN
       E2 = (/0D0,1D0,0D0/)
    ELSE
       CALL CROSS_PRODUCT(E3, (/1D0,0D0,0D0/),E2)
    END IF
    CALL CROSS_PRODUCT(E2, E3, E1)

    CALL NORMALIZE(E1); CALL NORMALIZE(E2); CALL NORMALIZE(E3)

    ! generate the axis around which to rotate
    X = sin(THETA)*cos(PHI)
    Y = sin(THETA)*sin(PHI)
    Z = cos(THETA)

    RANDAX = X*E1 + Y*E2 + Z*E3

  END SUBROUTINE RANDOMAXIS

  SUBROUTINE GETPERP(V1,V2)
    ! get some unit vector perpendicular to V1 and store it in V2
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: V1(3)
    DOUBLE PRECISION, INTENT(OUT) :: V2(3)

    IF (V1(2).EQ.0.AND.V1(3).EQ.0) THEN
       V2 = (/0D0,1D0,0D0/)
    ELSE
       CALL CROSS_PRODUCT(V1,(/0D0,1D0,0D0/),V2)
       CALL NORMALIZE(V2)
    ENDIF
  END SUBROUTINE GETPERP

  ! ZCS, 12/17/2021 binary search function, returns index
  INTEGER FUNCTION BINSEARCH(A,B)
    ! Given sorted array A, perform binary search algorithm to find index of element closest to B
    IMPLICIT NONE
    DOUBLE PRECISION :: A(:), B
    INTEGER :: L, R, M

    L = 1
    R = SIZE(A)

    IF (L.EQ.R) THEN
      BINSEARCH = 1
    ENDIF

    ! check if search item is less than whole array, or greater than whole array, if so return boundaries
    IF (A(L).GT.B) THEN
      BINSEARCH = L
      RETURN
    ELSEIF (A(R).LT.B) THEN
      BINSEARCH = R
      RETURN
    ENDIF

    ! while L not equal to R
    DO WHILE (L.NE.R)
      ! update midpt index
      M = INT((L+R)/2)

      ! greater, move R boundary
      IF (A(M) > B) THEN
        R = M
      ! less, move L boundary
      ELSE
        L = M + 1
      ENDIF
    ENDDO

  ! this makes sure that the index returned is greater than B
  IF ((A(M).GT.B)) THEN
    BINSEARCH = M
  ELSE
    BINSEARCH = M + 1
  ENDIF

  END FUNCTION BINSEARCH



  SUBROUTINE LLINTERSECT(X1,Y1,X2,Y2,X3,Y3,X4,Y4,CROSS,T,U)
    ! given two sets of points, (1,2) (3,4) each defining a line segment, determine whether segments cross
    ! CROSS is a logical for whether they cross. If they do cross, T tells where along (1,2) segment
    ! the crossing occurs. U is where along (3,4)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: X1, X2, X3, X4, Y1, Y2, Y3, Y4
    LOGICAL, INTENT(OUT) :: CROSS
    DOUBLE PRECISION, INTENT(OUT) :: T, U
    DOUBLE PRECISION :: DEN

    CROSS=.FALSE.
    DEN = ((X1-X2)*(Y3-Y4) - (Y1-Y2)*(X3-X4))
    IF (DEN.EQ.0D0) THEN
      RETURN
    ENDIF

    T = ((X1-X3)*(Y3-Y4) - (Y1-Y3)*(X3-X4)) / DEN
    U = ((X1-X3)*(Y1-Y2) - (Y1-Y3)*(X1-X2)) / DEN
    ! PRINT*, T, U

    IF ((T.GE.0D0).AND.(T.LE.1D0).AND.(U.GE.0D0).AND.(U.LE.1D0)) THEN
      ! crossing found
      CROSS=.TRUE.
    ENDIF


  END SUBROUTINE LLINTERSECT



  SUBROUTINE LLSLIDEINT(X1I,Y1I,X2I,Y2I,X4I,Y4I,X1F,Y1F,X2F,Y2F,X4F,Y4F,CROSS,T,TSTAR,TOUT2,TSTAR2)
    ! finds crossing of two lines changing in time assuming linear interpolation of end positions
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: X1I,Y1I,X2I,Y2I,X4I,Y4I,X1F,Y1F,X2F,Y2F,X4F,Y4F
    INTEGER, INTENT(OUT) :: CROSS
    DOUBLE PRECISION, INTENT(OUT) :: T, TSTAR, TOUT2, TSTAR2
    DOUBLE PRECISION :: DEN, TST1, TST2, SQRTARG, T1, T2

    ! cross = 0 initially, signalling no crossing event
    ! cross = 1 means 1 cross event, cross = 2 means 2 crossing events
    CROSS = 0

    SQRTARG = -4*(X4I*(-Y1I+Y2I)+X2I*(Y1I-Y4I)+X1I*(-Y2I+Y4I)) & 
            & *(-(X4F-X4I)*(Y1F-Y1I-Y2F+Y2I)+X2I*(-Y1F+Y1I+Y4F-Y4I) &
            & + X2F*(Y1F-Y1I-Y4F+Y4I)-(X1F-X1I)*(Y2F-Y2I-Y4F+Y4I)) &
            & + (X2F*Y1I-X4F*Y1I+X4I*(-Y1F+2*Y1I+Y2F-2*Y2I) &
            & + X4F*Y2I+X1I*(-Y2F+2*Y2I+Y4F-2*Y4I)-X2F*Y4I &
            & + X1F*(-Y2I+Y4I)+X2I*(Y1F-2*Y1I-Y4F+2*Y4I))**2

    ! if arg to sqrt is less than zero, result is complex, no crossing found
    IF (SQRTARG<0D0) THEN
      RETURN
    ENDIF

    ! find tstar (two different solutions)
    DEN = 2*((X4F-X4I)*(Y1F-Y1I-Y2F+Y2I)+X2F*(-Y1F+Y1I+Y4F-Y4I)+X2I*(Y1F-Y1I-Y4F+Y4I)+(X1F-X1I)*(Y2F-Y2I-Y4F+Y4I))

    IF (DEN.EQ.0D0) THEN
      RETURN
    ENDIF

    TST1 = (-X4I*Y1F+X2F*Y1I-X4F*Y1I+2*X4I*Y1I-X1I*Y2F+X4I*Y2F &
           & - X1F*Y2I+2*X1I*Y2I+X4F*Y2I-2*X4I*Y2I+X1I*Y4F+X1F*Y4I  &
           & - 2*X1I*Y4I-X2F*Y4I+X2I*(Y1F-2*Y1I-Y4F+2*Y4I) + SQRT(SQRTARG))/DEN

    TST2 = -(-X2I*Y1F+X4I*Y1F-X2F*Y1I+2*X2I*Y1I+X4F*Y1I-2*X4I*Y1I+X1I*Y2F &
           & - X4I*Y2F+X1F*Y2I-2*X1I*Y2I-X4F*Y2I+2*X4I*Y2I-X1I*Y4F &
           & + X2I*Y4F-X1F*Y4I+2*X1I*Y4I+X2F*Y4I-2*X2I*Y4I + SQRT(SQRTARG))/DEN

    IF ((TST1.LE.1D0).AND.(TST1.GE.0D0)) THEN
      ! check if tstar2 also passed
      IF ((TST2.LE.1D0).AND.(TST2.GE.0D0)) THEN

        T1 = (TST1*(X4F-X4I)*(Y1F-Y1I)-X4I*Y1I-X1I*Y4F+X4I*Y4F+X4F*(Y1I-Y4I)-TST1*(X1F-X1I)*(Y4F-Y4I)+X1I*Y4I) &
          & /((X4F-X4I)*(Y1I-Y2I)+TST1*(X4F-X4I)*(Y1F-Y1I-Y2F+Y2I)-(X1I-X2I)*(Y4F-Y4I)-TST1*(X1F-X1I-X2F+X2I)*(Y4F-Y4I))
        T2 = (TST2*(X4F-X4I)*(Y1F-Y1I)-X4I*Y1I-X1I*Y4F+X4I*Y4F+X4F*(Y1I-Y4I)-TST2*(X1F-X1I)*(Y4F-Y4I)+X1I*Y4I) &
          & /((X4F-X4I)*(Y1I-Y2I)+TST2*(X4F-X4I)*(Y1F-Y1I-Y2F+Y2I)-(X1I-X2I)*(Y4F-Y4I)-TST2*(X1F-X1I-X2F+X2I)*(Y4F-Y4I))

        IF ((T1.LE.1D0).AND.(T1.GE.0D0)) THEN
          IF ((T2.LE.1D0).AND.(T2.GE.0D0)) THEN
            IF (TST1.LE.TST2) THEN
              T=T1
              TSTAR=TST1
              TOUT2=T2
              TSTAR2=TST2
              CROSS=2
              ! this is allowed, just means crossed the line twice during time step
              ! return both crossings, with earlier time first
              RETURN

            ELSE
              T=T2
              TSTAR=TST2
              TOUT2=T1
              TSTAR2=TST1
              CROSS=2
              ! this is allowed, just means crossed the line twice during time step
              ! return both crossings, with earlier time first
              RETURN

            ENDIF
          ENDIF

          ! only have 0<T1<1
          T=T1
          TSTAR=TST1
          CROSS=1
          RETURN

        ELSEIF ((T2.LE.1D0).AND.(T2.GE.0D0)) THEN
          ! only have 0<T2<1
          T=T2
          TSTAR=TST2
          CROSS=1
          RETURN

        ENDIF

        ! if you get here, two valid TSTS, but neither T1 nor T2 were valid, no crossing
        RETURN
      ENDIF

      ! only TST1 was valid, solve for T
      ! solve for T (the fraction along first edge the intersection occurs)
      T = (TST1*(X4F-X4I)*(Y1F-Y1I)-X4I*Y1I-X1I*Y4F+X4I*Y4F+X4F*(Y1I-Y4I)-TST1*(X1F-X1I)*(Y4F-Y4I)+X1I*Y4I) &
        & /((X4F-X4I)*(Y1I-Y2I)+TST1*(X4F-X4I)*(Y1F-Y1I-Y2F+Y2I)-(X1I-X2I)*(Y4F-Y4I)-TST1*(X1F-X1I-X2F+X2I)*(Y4F-Y4I))
      TSTAR=TST1

      IF ((T.LE.1D0).AND.(T.GE.0D0)) THEN
        CROSS = 1
      ENDIF

    ELSEIF ((TST2.LE.1D0).AND.(TST2.GE.0D0)) THEN

      ! only TST2 was valid, solve for T
      ! solve for T (the fraction along first edge the intersection occurs)
      T = (TST2*(X4F-X4I)*(Y1F-Y1I)-X4I*Y1I-X1I*Y4F+X4I*Y4F+X4F*(Y1I-Y4I)-TST2*(X1F-X1I)*(Y4F-Y4I)+X1I*Y4I) &
        & /((X4F-X4I)*(Y1I-Y2I)+TST2*(X4F-X4I)*(Y1F-Y1I-Y2F+Y2I)-(X1I-X2I)*(Y4F-Y4I)-TST2*(X1F-X1I-X2F+X2I)*(Y4F-Y4I))
      TSTAR=TST2

      IF ((T.LE.1D0).AND.(T.GE.0D0)) THEN
        CROSS = 1
      ENDIF

    ! neither TST1 nor TST2 were in valid range, no crossing
    ELSE
      RETURN
    ENDIF


  END SUBROUTINE LLSLIDEINT


  SUBROUTINE LLSLIDEINT2(X1,Y1,X2,Y2,X4I,Y4I,X4F,Y4F,CROSS,T,TSTAR)
    ! finds crossing of two lines where only the last is changing in time assuming linear interpolation of end positions
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: X1,Y1,X2,Y2,X4I,Y4I,X4F,Y4F
    INTEGER, INTENT(OUT) :: CROSS
    DOUBLE PRECISION, INTENT(OUT) :: T, TSTAR
    DOUBLE PRECISION :: DEN, SQRTARG

    CROSS = 0

    DEN = -((X4F-X4I)*(Y1-Y2)) + (X1-X2)*(Y4F-Y4I)

    IF (DEN.EQ.0D0) THEN
      PRINT*, 'DEN = ZERO in LLSLIDEINT2'
      RETURN
    ENDIF

    TSTAR = (X4I*(Y1-Y2) + X1*(Y2-Y4I) + X2*(-Y1+Y4I))/DEN
  
    IF ((TSTAR.LE.1D0).AND.(TSTAR.GE.0D0)) THEN

      ! solve for T (the fraction along first edge the intersection occurs)
      T = (X4I*(Y1-Y4F) + X1*(Y4F-Y4I) + X4F*(-Y1+Y4I))/DEN

      IF ((T.LE.1D0).AND.(T.GE.0D0)) THEN
        CROSS = 1
      ENDIF

    ENDIF


  END SUBROUTINE LLSLIDEINT2


! Takes the coordinates of 2 nodes connected by an edge, 2 more nodes
! connected by another edge, the dimension of the network space (2D, 3D, etc.)
! and fills in the values of MUA and MUB, which denote the closest points
! on the first and second edges, respectively, as weights of the original coords.
  SUBROUTINE SHORTESTPATHPOINTS(P1,P2,P3,P4,DIM,MUA,MUB,MUBYR)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: DIM
    DOUBLE PRECISION, INTENT(IN) :: P1(DIM),P2(DIM),P3(DIM),P4(DIM)
    DOUBLE PRECISION, INTENT(OUT) :: MUA, MUB, MUBYR(DIM,2,4)
    DOUBLE PRECISION :: D1343, D4321, D1321, D4343, D2121, DENOM
    DOUBLE PRECISION, DIMENSION(DIM) :: R21,R13,R43
    LOGICAL :: AOUT,BOUT

    AOUT=.FALSE.; BOUT=.FALSE.

    R21 = P2-P1
    R13 = P1-P3
    R43 = P4-P3
    D4321 = DOT_PRODUCT(R43,R21)
    D4343 = DOT_PRODUCT(R43,R43)
    D2121 = DOT_PRODUCT(R21,R21)
    D1343 = DOT_PRODUCT(R13,R43)
    D1321 = DOT_PRODUCT(R13,R21)
    DENOM = D2121*D4343 - D4321*D4321

    IF(DENOM.NE.0) THEN    
       MUA = (D1343*D4321 - D1321*D4343)/DENOM
       IF(MUA.GT.1D0) THEN
         AOUT = .TRUE.; MUA = 1D0
       ELSEIF(MUA.LT.0D0) THEN
         AOUT = .TRUE.; MUA = 0D0
       ENDIF
       MUB = (D1343 + MUA*D4321)/D4343
      
       IF(MUB.GT.1D0) THEN
          BOUT = .TRUE.; MUB = 1D0
       ELSEIF(MUB.LT.0D0) THEN
          BOUT = .TRUE.; MUB = 0D0
       ENDIF

       IF(BOUT) THEN
          MUA = (MUB*D4321 - D1321)/D2121
          IF(MUA.GT.1D0) THEN
            AOUT = .TRUE.; MUA = 1D0
          ELSEIF(MUA.LT.0D0) THEN
            AOUT = .TRUE.; MUA = 0D0
          ELSE
            AOUT = .FALSE.
          ENDIF
       ENDIF

    ELSE ! parallel lines
       MUA = 0D0
       MUB = D1343/D4343
       IF(MUB.GT.1D0) THEN
          BOUT = .TRUE.; MUB = 1D0
       ELSEIF(MUB.LT.0D0) THEN
          BOUT = .TRUE.; MUB = 0D0
       ENDIF

       IF(BOUT) THEN
          MUA = (MUB*D4321-D1321)/D2121
          IF(MUA.GT.1D0) THEN
            MUA = 1D0
          ELSEIF(MUA.LT.0D0) THEN
            MUA = 0D0
          ENDIF
       ENDIF
       ! for parallel lines the derivatives are undefined so best to set them equal to zero
       AOUT = .TRUE.; BOUT = .TRUE.
    ENDIF

    ! DIST = SQRT(SUM((R21*MUA-R43*MUB+R13)**2))
    
    ! now find the derivatives
    IF(AOUT.AND.BOUT) THEN ! calculate the 4 point-point distances and 4 projection distances
       MUBYR = 0D0 ! no dependence. Both mua,mub are constants
    ELSEIF(AOUT) THEN ! Do single point Find on MUB
       !PRINT*,'MUA IS OUT'
       MUBYR(:,1,:) = 0D0
       IF(MUA.GT.5D-1) THEN
          !MUB = DOT_PRODUCT(P2-P3,P4-P3)/D4343
          MUBYR(:,2,1) = 0D0
          MUBYR(:,2,2) = R43/D4343
          MUBYR(:,2,3) = (2*P3-P2-P4 + 2*MUB*R43)/D4343
          MUBYR(:,2,4) = (P2-P3-2*MUB*R43)/D4343
       ELSE
          !MUB = DOT_PRODUCT(P1-P3,P4-P3)/D4343
          MUBYR(:,2,1) = R43/D4343
          MUBYR(:,2,2) = 0D0
          MUBYR(:,2,3) = (2*P3-P1-P4 + 2*MUB*R43)/D4343
          MUBYR(:,2,4) = (P1-P3-2*MUB*R43)/D4343
       ENDIF
    ELSEIF(BOUT) THEN ! Do single point find on MUA
       !PRINT*,'MUB IS OUT'
       MUBYR(:,2,:) = 0D0
       IF(MUB.GT.5D-1) THEN
          !MUA = DOT_PRODUCT(P4-P1,P2-P1)/D2121
          MUBYR(:,1,1) = (2*P1-P4-P2 + 2*MUA*R21)/D2121
          MUBYR(:,1,2) = (P4-P1-2*MUA*R21)/D2121
          MUBYR(:,1,3) = 0D0
          MUBYR(:,1,4) = R21/D2121
       ELSE
          !MUA = DOT_PRODUCT(P3-P1,P2-P1)/D2121
          MUBYR(:,1,1) = (2*P1-P3-P2 + 2*MUA*R21)/D2121
          MUBYR(:,1,2) = (P3-P1-2*MUA*R21)/D2121
          MUBYR(:,1,3) = R21/D2121
          MUBYR(:,1,4) = 0D0
       ENDIF
    ELSE
       ! Use the original method - no end points are used
       MUBYR(:,1,1) = (R43*(D4321-D1343) - D4343*(R21-R13) - 2*MUA*(D4321*R43 - D4343*R21)) / DENOM
       MUBYR(:,1,2) = ((D1343*R43 - D4343*R13) - 2*MUA*(D4343*R21 - D4321*R43)) / DENOM
       MUBYR(:,1,3) = (R21*(D4343-D1343) + 2*D1321*R43 - D4321*(R43+R13) - 2*MUA*(D4321*R21-D2121*R43)) / DENOM
       MUBYR(:,1,4) = (R13*D4321 + D1343*R21 - 2*D1321*R43 - 2*MUA*(D2121*R43 - D4321*R21)) / DENOM
       MUBYR(:,2,1) = ((1-MUA)*R43 + D4321*MUBYR(:,1,1)) / D4343
       MUBYR(:,2,2) = (D4321*MUBYR(:,1,2) + MUA*R43) / D4343
       MUBYR(:,2,3) = (-R43-R13 - MUA*R21 + MUBYR(:,1,3)*D4321 + 2*MUB*R43) / D4343
       MUBYR(:,2,4) = (R13 + MUA*R21 + MUBYR(:,1,4)*D4321 - 2*MUB*R43) / D4343
    ENDIF

  END SUBROUTINE SHORTESTPATHPOINTS


END MODULE GENUTIL
