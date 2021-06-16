CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Eight Gauss points (2x2x2) for 8-node brick elmt
C     natural coords & weights 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GSPT_3D_8N(iTIntPt, rXiEtaPt, rWeightPt)  

C     OUTPUT:
C             iTIntPt: Total number of quadrature points 
C            rXiEtaPt: Natural (element) coordinates of all Gauss points
C           rWeightPt: Integration weights for all Gauss points
      
      IMPLICIT NONE
      
      real*8 rXiEtaPt(3,8),rWeightPt(8)
      integer iTIntPt
C
      real*8 ONE,ONENEG,THREE
      PARAMETER(ONE=1.D0,ONENEG=-1.D0,THREE=3.D0)

      real*8 R,AR(2),GWE(2)
      integer I,J,K,NUMGP


      R = SQRT(ONE/THREE)
      AR(1)=ONENEG
      AR(2)=ONE
      GWE(1)=ONE
      GWE(2)=ONE
      iTIntPt = 8

      DO I = 1,2
        DO J = 1,2
          DO K = 1,2
            NUMGP = (I-1)*4+(J-1)*2+K
            rXiEtaPt(1,NUMGP) = AR(I)*R
            rXiEtaPt(2,NUMGP) = AR(J)*R
            rXiEtaPt(3,NUMGP) = AR(K)*R
            rWeightPt(NUMGP)=GWE(I)*GWE(J)*GWE(K)
          END DO ! K
        END DO   ! J
      END DO     ! I
      RETURN
      END   ! SUBROUTINE GSPT_3D_8N
