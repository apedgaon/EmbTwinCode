CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calculate the shortest length in a 3D 8-node brick elmt
C     as characteristic length
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE Clength_3D_8N(COORDS,LMIN)

      IMPLICIT NONE

      REAL*8 COORDS(3,8),LMIN
      REAL*8 L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12

      L1 = sqrt( (COORDS(1,2) - COORDS(1,1))**2
     1         + (COORDS(2,2) - COORDS(2,1))**2
     2         + (COORDS(3,2) - COORDS(3,1))**2 )
      L2 = sqrt( (COORDS(1,3) - COORDS(1,2))**2
     1         + (COORDS(2,3) - COORDS(2,2))**2
     2         + (COORDS(3,3) - COORDS(3,2))**2 )
      L3 = sqrt( (COORDS(1,4) - COORDS(1,3))**2
     1         + (COORDS(2,4) - COORDS(2,3))**2
     2         + (COORDS(3,4) - COORDS(3,3))**2 )
      L4 = sqrt( (COORDS(1,4) - COORDS(1,1))**2
     1         + (COORDS(2,4) - COORDS(2,1))**2
     2         + (COORDS(3,4) - COORDS(3,1))**2 )

      L5 = sqrt( (COORDS(1,6) - COORDS(1,5))**2
     1         + (COORDS(2,6) - COORDS(2,5))**2
     2         + (COORDS(3,6) - COORDS(3,5))**2 )
      L6 = sqrt( (COORDS(1,7) - COORDS(1,6))**2
     1         + (COORDS(2,7) - COORDS(2,6))**2
     2         + (COORDS(3,7) - COORDS(3,6))**2 )
      L7 = sqrt( (COORDS(1,8) - COORDS(1,7))**2
     1         + (COORDS(2,8) - COORDS(2,7))**2
     2         + (COORDS(3,8) - COORDS(3,7))**2 )
      L8 = sqrt( (COORDS(1,8) - COORDS(1,5))**2
     1         + (COORDS(2,8) - COORDS(2,5))**2
     2         + (COORDS(3,8) - COORDS(3,5))**2 )

      L9 = sqrt( (COORDS(1,5) - COORDS(1,1))**2
     1         + (COORDS(2,5) - COORDS(2,1))**2
     2         + (COORDS(3,5) - COORDS(3,1))**2 )
      L10 = sqrt( (COORDS(1,6) - COORDS(1,2))**2
     1         +  (COORDS(2,6) - COORDS(2,2))**2
     2         +  (COORDS(3,6) - COORDS(3,2))**2 )
      L11 = sqrt( (COORDS(1,7) - COORDS(1,3))**2
     1         +  (COORDS(2,7) - COORDS(2,3))**2
     2         +  (COORDS(3,7) - COORDS(3,3))**2 )
      L12 = sqrt( (COORDS(1,8) - COORDS(1,4))**2
     1         +  (COORDS(2,8) - COORDS(2,4))**2
     2         +  (COORDS(3,8) - COORDS(3,4))**2 )

      LMIN = min(L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12)     

      RETURN
      END

