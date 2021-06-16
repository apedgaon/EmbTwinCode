CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Evaluate interpolation functions, their derivatives and
C     Jacobians at Gauss point (Xi1,Xi2,Xi3) for
C     8-node brick elmt
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SHP_FUNC_3D_8N(rXiEtaPt,COORDS,iDeriv,SHP,DEV_SHP,JACOB)

C     INPUT:
C          rXiEtaPt: Natural (element) coordinates of a Gauss point
C          COORDS(3,8): Cartesian coords of eight nodal pts
C                       COORDS(1,I) is X-coord of Node I
C                       COORDS(2,I) is Y-coord of Node I
C                       COORDS(3,I) is Z-coord of Node I
C          iDeriv:      =  0, Derivatives WRT natural coords
C                       != 0, Derivatives WRT Cartesian coords

C     OUTPUT:
C          SHP(8): Eight SHP func. values at rXiEtaPt = (Xi1,Xi2,Xi3)
C          DEV_SHP(3,8): when iDeriv = 0
C                        Derivative of SHP w.r.t. Natural coords, SHP_N
C                      ( N_1,Xi1 N_2,Xi1 ... N_8,Xi1
C                        N_1,Xi2 N_2,Xi2 ... N_8,Xi2
C                        N_1,Xi3 N_2,Xi3 ... N_8,Xi3 )
C                        when iDeriv != 0
C                        Derivative of SHP w.r.t. Carterian coords, SHP_C
C                      ( N_1,X N_2,X ... N_8,X
C                        N_1,Y N_2,Y ... N_8,Y
C                        N_1,Z N_2,Z ... N_8,Z )
C          JACOB: Determinant of JACOB_MTX(3,3)

      IMPLICIT NONE

C     INPUT:
      REAL*8 rXiEtaPt(3),COORDS(3,8)
      INTEGER iDeriv

C     OUTPUT:
      REAL*8 SHP(8),DEV_SHP(3,8),JACOB

C     LOCAL VAR
      REAL*8 ZERO,EIGHTH,ONE
      PARAMETER(ZERO=0.D0,EIGHTH=0.125D0,ONE=1.D0)

      REAL*8 temp
      REAL*8 Xi1,Xi2,Xi3,JACOB_MTX(3,3),SHP_N(3,8),SHP_C(3,8)
      INTEGER I,J,K

C     Xi1,Xi2,Xi3: natural coordinates;
      Xi1 = rXiEtaPt(1)
      Xi2 = rXiEtaPt(2)
      Xi3 = rXiEtaPt(3)

C     INTERPOLATION FUNCTIONS N_A
      SHP(1) = EIGHTH*(ONE - Xi1)*(ONE - Xi2)*(ONE - Xi3)
      SHP(2) = EIGHTH*(ONE + Xi1)*(ONE - Xi2)*(ONE - Xi3)
      SHP(3) = EIGHTH*(ONE + Xi1)*(ONE + Xi2)*(ONE - Xi3)
      SHP(4) = EIGHTH*(ONE - Xi1)*(ONE + Xi2)*(ONE - Xi3)
      SHP(5) = EIGHTH*(ONE - Xi1)*(ONE - Xi2)*(ONE + Xi3)
      SHP(6) = EIGHTH*(ONE + Xi1)*(ONE - Xi2)*(ONE + Xi3)
      SHP(7) = EIGHTH*(ONE + Xi1)*(ONE + Xi2)*(ONE + Xi3)
      SHP(8) = EIGHTH*(ONE - Xi1)*(ONE + Xi2)*(ONE + Xi3)

C     DERIVATIVES WRT TO Xi1 (Partial N_A)/(Partial Xi1)
      SHP_N(1,1) = -EIGHTH*(ONE - Xi2)*(ONE - Xi3)
      SHP_N(1,2) =  EIGHTH*(ONE - Xi2)*(ONE - Xi3)
      SHP_N(1,3) =  EIGHTH*(ONE + Xi2)*(ONE - Xi3)
      SHP_N(1,4) = -EIGHTH*(ONE + Xi2)*(ONE - Xi3)
      SHP_N(1,5) = -EIGHTH*(ONE - Xi2)*(ONE + Xi3)
      SHP_N(1,6) =  EIGHTH*(ONE - Xi2)*(ONE + Xi3)
      SHP_N(1,7) =  EIGHTH*(ONE + Xi2)*(ONE + Xi3)
      SHP_N(1,8) = -EIGHTH*(ONE + Xi2)*(ONE + Xi3)

C     DERIVATIVES WRT TO Xi2 (Partial N_A)/(Partial Xi2)
      SHP_N(2,1) = -EIGHTH*(ONE - Xi1)*(ONE - Xi3)
      SHP_N(2,2) = -EIGHTH*(ONE + Xi1)*(ONE - Xi3)
      SHP_N(2,3) =  EIGHTH*(ONE + Xi1)*(ONE - Xi3)
      SHP_N(2,4) =  EIGHTH*(ONE - Xi1)*(ONE - Xi3)
      SHP_N(2,5) = -EIGHTH*(ONE - Xi1)*(ONE + Xi3)
      SHP_N(2,6) = -EIGHTH*(ONE + Xi1)*(ONE + Xi3)
      SHP_N(2,7) =  EIGHTH*(ONE + Xi1)*(ONE + Xi3)
      SHP_N(2,8) =  EIGHTH*(ONE - Xi1)*(ONE + Xi3)

C     DERIVATIVES WRT TO Xi3 (Partial N_A)/(Partial Xi3)
      SHP_N(3,1) = -EIGHTH*(ONE - Xi1)*(ONE - Xi2)
      SHP_N(3,2) = -EIGHTH*(ONE + Xi1)*(ONE - Xi2)
      SHP_N(3,3) = -EIGHTH*(ONE + Xi1)*(ONE + Xi2)
      SHP_N(3,4) = -EIGHTH*(ONE - Xi1)*(ONE + Xi2)
      SHP_N(3,5) =  EIGHTH*(ONE - Xi1)*(ONE - Xi2)
      SHP_N(3,6) =  EIGHTH*(ONE + Xi1)*(ONE - Xi2)
      SHP_N(3,7) =  EIGHTH*(ONE + Xi1)*(ONE + Xi2)
      SHP_N(3,8) =  EIGHTH*(ONE - Xi1)*(ONE + Xi2)

C     (Partial X,Y,Z)/(Partial Xi1,Xi2,Xi3)
C     JACOB_MTX(3,3): ( DXDXi1,DXDXi2,DXDXi3 
C                       DYDXi1,DYDXi2,DYDXi3
C                       DZDXi1,DZDXi2,DZDXi3 )
      DO I = 1,3
        DO J = 1,3
          temp = ZERO
          DO K = 1,8
            temp = temp + COORDS(I,K)*SHP_N(J,K)
          END DO
          JACOB_MTX(I,J) = temp
        END DO ! J
      END DO ! I

C     CALCULATION OF JACOBIAN, DET(Partial X/Partial Xi)
      JACOB = JACOB_MTX(1,1)*JACOB_MTX(2,2)*JACOB_MTX(3,3)
     1      + JACOB_MTX(1,2)*JACOB_MTX(2,3)*JACOB_MTX(3,1)
     2      + JACOB_MTX(1,3)*JACOB_MTX(2,1)*JACOB_MTX(3,2)
     3      - JACOB_MTX(1,3)*JACOB_MTX(2,2)*JACOB_MTX(3,1)
     4      - JACOB_MTX(1,1)*JACOB_MTX(2,3)*JACOB_MTX(3,2)
     5      - JACOB_MTX(1,2)*JACOB_MTX(2,1)*JACOB_MTX(3,3)

      IF (JACOB .LE. 0.d0) THEN
        CALL XPLB_ABQERR(-3,'User elmt VU1 error: '//
     $                   'Non-positive Jacobian determinant',0,0.d0,' ')
      END IF

C     DERIVATIVES WRT TO X,Y and Z, (Partial N_A/Partial X,Y,Z)
      DO I = 1,8
        SHP_C(1,I) = ( SHP_N(1,I)*(JACOB_MTX(2,2)*JACOB_MTX(3,3) 
     1                           - JACOB_MTX(2,3)*JACOB_MTX(3,2) ) 
     2             +   SHP_N(2,I)*(JACOB_MTX(2,3)*JACOB_MTX(3,1)
     3                           - JACOB_MTX(2,1)*JACOB_MTX(3,3) )
     4             +   SHP_N(3,I)*(JACOB_MTX(2,1)*JACOB_MTX(3,2)
     5                           - JACOB_MTX(2,2)*JACOB_MTX(3,1) )
     6                )/JACOB

        SHP_C(2,I) = ( SHP_N(1,I)*(JACOB_MTX(1,3)*JACOB_MTX(3,2) 
     1                           - JACOB_MTX(1,2)*JACOB_MTX(3,3) ) 
     2             +   SHP_N(2,I)*(JACOB_MTX(1,1)*JACOB_MTX(3,3)
     3                           - JACOB_MTX(1,3)*JACOB_MTX(3,1) )
     4             +   SHP_N(3,I)*(JACOB_MTX(1,2)*JACOB_MTX(3,1)
     5                           - JACOB_MTX(1,1)*JACOB_MTX(3,2) )
     6                )/JACOB

        SHP_C(3,I) = ( SHP_N(1,I)*(JACOB_MTX(1,2)*JACOB_MTX(2,3) 
     1                           - JACOB_MTX(1,3)*JACOB_MTX(2,2) ) 
     2             +   SHP_N(2,I)*(JACOB_MTX(1,3)*JACOB_MTX(2,1)
     3                           - JACOB_MTX(1,1)*JACOB_MTX(2,3) )
     4             +   SHP_N(3,I)*(JACOB_MTX(1,1)*JACOB_MTX(2,2)
     5                           - JACOB_MTX(1,2)*JACOB_MTX(2,1) )
     6                )/JACOB
      END DO  ! I

      IF (iDeriv.EQ.0) THEN
        DO I = 1,3
          DO J = 1,8
            DEV_SHP(I,J) = SHP_N(I,J) 
          END DO
        END DO
      ELSE
        DO I = 1,3
          DO J = 1,8
            DEV_SHP(I,J) = SHP_C(I,J) 
          END DO
        END DO
      END IF


      RETURN
      END


