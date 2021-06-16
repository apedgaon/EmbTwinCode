      subroutine Gstrain_3D(u,rEps,rNn,rNnCt,rdNdx,rdNdxCt,mcrd,nnode,
     $                      iVelmType,iVforType,iSpinflg,
     $                      rStrain,rSpin,rBmat,rBskw)

C     Compute strain or strain increment, and the corresponding
C     B (or B-bar) matrix, at a (Gauss) point in 8-node Brick element
C
C     Given:
C                   u: Nodal displacements 
C                rEps: Pressure mode stabilization parameter; see Ortiz et al (1987), Eq. (21).
C                 rNn: Shape functions at current quad point
C               rNnCt: Shape functions at element centroid
C               rdNdx: Shape function derivatives wrt global coords
C             rdNdxCt: Shape function derivatives wrt global coords at element centroid
C                mcrd: Number of DOFs per node (=3)
C           iVelmType: 3D element type,  1 --> 8-node brick elmt

C           iVforType: Formulation type: 0 --> displacement based
C                                        1 --> displacement based with localization behavior
C                                        2 --> B-bar incompressible
C                                        3 --> B-bar incompressible with localization behavior
C            iSpinflg: If this flag is equal to 1, compute spin axial vector
C
C
C     Return the following, evaluated at the given (Gauss) point:
C
C             rStrain: strain increment, defined such that
C        {rStrain} = < u1,1 | u2,2 | u3,3 | u2,3+u3,2 | u1,3+u3,1 | u1,2+u2,1 >^T
C               rSpin: spin axial vector (increment), defined such that
C          {rSpin} = 0.5 * < u3,2-u2,3 | u1,3-u3,1 | u2,1-u1,2 >^T
C               rBmat: Symmetric B (or B-bar) matrix, defined such that
C           [rBmat]{d} = < u1,1 | u2,2 | u3,3 | u2,3+u3,2 | u1,3+u3,1 | u1,2+u2,1 >^T
C               rBskw: Skew-symm. B (or B-bar) matrix, defined such that
C           [rBskw]{d} = < u3,2-u2,3 | u1,3-u3,1 | u2,1-u1,2 >^T
C
C
C     Notes: The localization behavior is NOT treated in this subroutine. The kinematic
C            quantities returned by the current subroutine are only the basis used
C            (in other routines) to calculate kinematic quantities for the Localized
C            and Non-localized regions.
C
C            We can compute the infinitesimal strain, if given the displacements u_(i+1),
C            and the shape function derivatives.
C
C            The most likely used of this subroutine, however, is to compute a strain
C            INCREMENT, USING THE SAME OPERATIONS AS FOR THE INFINITESIMAL STRAIN, but 
C            here u will be the INCREMENT in displacement between time steps, namely
C            u_(i+1) - u_(i), and the shape function derivatives will be WRT to the 
C            mid-step coordinates x_(i+1/2) = 1/2 * [ x_(i) + x(i+1) ]
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT NONE

      integer mcrd,nnode,iVelmType,iVforType,iSpinflg,i,n,k

      real*8 rEps
      real*8 u(mcrd,nnode), rNn(nnode), rNnCt(nnode)
      real*8 rdNdx(mcrd,nnode), rdNdxCt(mcrd,nnode)
      real*8 rStrain(6), rSpin(3),rBmat(6,mcrd*nnode),rBskw(3,mcrd*nnode)
      real*8 rBvol(6,mcrd*nnode),rB0vol(6,mcrd*nnode),rB(6,mcrd*nnode)

      real*8 r1thrd
      parameter ( r1thrd = 1.0d0 / 3.0d0 )

      DO i = 1, 6
        rStrain(i) = 0.d0
        DO n = 1, nnode
          DO k = 1, mcrd
            rBmat  (i  , (n-1)*mcrd + k ) = 0.d0
            rBvol  (i  , (n-1)*mcrd + k ) = 0.d0
            rB0vol (i  , (n-1)*mcrd + k ) = 0.d0
            rB     (i  , (n-1)*mcrd + k ) = 0.d0   
          END DO
        END DO
      END DO

      DO i = 1, 3
        rSpin(i)   = 0.d0
        DO n = 1, nnode
          DO k = 1, mcrd
            rBskw (i  , (n-1)*mcrd + k ) = 0.d0
          END DO
        END DO
      END DO

C     Displacement-based formulations
      IF ( (iVforType .EQ. 0) .OR. (iVforType .EQ. 1) ) THEN 
        DO n = 1, nnode ! Loop over nodes
          rBmat(1, (n-1)*mcrd + 1 ) = rdNdx(1,n)
          rBmat(2, (n-1)*mcrd + 2 ) = rdNdx(2,n)
          rBmat(3, (n-1)*mcrd + 3 ) = rdNdx(3,n)

          rBmat(4, (n-1)*mcrd + 2 ) = rdNdx(3,n)
          rBmat(4, (n-1)*mcrd + 3 ) = rdNdx(2,n)

          rBmat(5, (n-1)*mcrd + 1 ) = rdNdx(3,n)
          rBmat(5, (n-1)*mcrd + 3 ) = rdNdx(1,n)

          rBmat(6, (n-1)*mcrd + 1 ) = rdNdx(2,n)
          rBmat(6, (n-1)*mcrd + 2 ) = rdNdx(1,n)
        END DO

C     Incompressible (B-bar) formulations
      ELSE IF ( (iVforType .EQ. 2) .OR. (iVforType .EQ. 3) ) THEN
C       B
        DO n = 1, nnode ! Loop over nodes
          rB(1, (n-1)*mcrd + 1 ) = rdNdx(1,n)
          rB(2, (n-1)*mcrd + 2 ) = rdNdx(2,n)
          rB(3, (n-1)*mcrd + 3 ) = rdNdx(3,n)

          rB(4, (n-1)*mcrd + 2 ) = rdNdxCt(3,n)
          rB(4, (n-1)*mcrd + 3 ) = rdNdxCt(2,n)

          rB(5, (n-1)*mcrd + 1 ) = rdNdxCt(3,n)
          rB(5, (n-1)*mcrd + 3 ) = rdNdxCt(1,n)

          rB(6, (n-1)*mcrd + 1 ) = rdNdxCt(2,n)
          rB(6, (n-1)*mcrd + 2 ) = rdNdxCt(1,n)          
        END DO

C       B_vol
        DO n = 1, nnode ! Loop over nodes
          DO i = 1, mcrd
            DO k = 1, mcrd
              rBvol(i, (n-1)*mcrd+k ) = r1thrd*rdNdx(k,n)
            END DO
          END DO          
        END DO

C       B0_vol at elmt centroid
        DO n = 1, nnode ! Loop over nodes
          DO i = 1, mcrd
            DO k = 1, mcrd
              rB0vol(i, (n-1)*mcrd+k ) = r1thrd*rdNdxCt(k,n)
            END DO
          END DO          
        END DO    

C       B-bar = B + (1-eps)*B0_vol - (1-eps)*B_vol
        DO n = 1, nnode
          DO i = 1, 6
            DO k = 1, mcrd
              rBmat( i, (n-1)*mcrd + k ) 
     $      = rB( i, (n-1)*mcrd + k ) 
     $      + (1.d0 - rEps)*rB0vol( i, (n-1)*mcrd + k ) 
     $      - (1.d0 - rEps)*rBvol( i, (n-1)*mcrd + k )
            END DO
          END DO
        END DO

      END IF   ! formulation type

C     strain
      DO n = 1, nnode
        DO i = 1, 6
          DO k = 1, mcrd
            rStrain(i) = rStrain(i) + rBmat(i, (n-1)*mcrd + k )*u(k,n)
          END DO
        END DO
      END DO

C     Spin axial vector
      IF (iSpinflg.EQ.1) THEN       ! need spin tensor
        DO n = 1, nnode
          rBskw(1, (n-1)*mcrd + 2 ) = -rdNdx(3,n)
          rBskw(1, (n-1)*mcrd + 3 ) =  rdNdx(2,n)

          rBskw(2, (n-1)*mcrd + 1 ) =  rdNdx(3,n)
          rBskw(2, (n-1)*mcrd + 3 ) = -rdNdx(1,n)

          rBskw(3, (n-1)*mcrd + 1 ) = -rdNdx(2,n)
          rBskw(3, (n-1)*mcrd + 2 ) =  rdNdx(1,n)

          DO i = 1, 3
            DO k = 1, mcrd
              rSpin(i) = rSpin(i) 
     $      + 0.5d0*rBskw(i, (n-1)*mcrd + k )*u(k,n)
            END DO
          END DO      

        END DO
      END IF

      END


