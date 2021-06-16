      subroutine GdefGrad_3D(rdNdx,rCoords,mcrd,nnode,iVelmType, rF)

C     Compute deformation gradient (or its inverse) at a (Gauss) point 
C     in 3D BRICK element
C
C     Given:
C               rdNdx: Shape function derivatives wrt INITIAL (or current) global coords
C             rCoords: CURRENT (or initial) nodal coordinates
C              mcrd: Number of DOFs per node
C           iVelmType: 3D element type,  1 --> 8-node BRICK elmt
C
C     Return the following, evaluated at the given (Gauss) point:
C
C                  rF: Deformation gradient (or its inverse)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      include 'vaba_param.inc'

      dimension rdNdx(mcrd,nnode), rCoords(mcrd,nnode)
      dimension rF(3,3)

      do i = 1, 3
        do j = 1, 3
          rF(i,j) = 0.d0
        enddo
      enddo

      do n = 1, nnode    ! loop over nodes

        do i = 1, mcrd
          do j = 1, mcrd
            rF(i,j) =  rF(i,j) + rdNdx(j,n) * rCoords(i,n)
          enddo
        enddo

      enddo    ! loop over nodes

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine GnewNvec_3D(rFinv, rN)

C     Calculate up-to-date orientation vector of band
C     
C
C     Given:
C               rFinv: Inverse deformation gradient wrt configuration at bifurcation 
C               rN   : Initial orientation of localization band
C
C     Return the following, evaluated at the given (Gauss) point:
C
C               rN: Updated orientation of localization band
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      include 'vaba_param.inc'

      dimension rFinv(3,3), rN(3)
      dimension rq(3)

c     q = transpose(F_inverse) . n_initial
      rq(1)   = rFinv(1,1)*rN(1) + rFinv(2,1)*rN(2) + rFinv(3,1)*rN(3)
      rq(2)   = rFinv(1,2)*rN(1) + rFinv(2,2)*rN(2) + rFinv(3,2)*rN(3)
      rq(3)   = rFinv(1,3)*rN(1) + rFinv(2,3)*rN(2) + rFinv(3,3)*rN(3)

c     We want the unit normal
      rNvnorm = sqrt( rq(1)**2 + rq(2)**2 + rq(3)**2 )
      rN(1)   = rq(1) / rNvnorm
      rN(2)   = rq(2) / rNvnorm
      rN(3)   = rq(3) / rNvnorm
      
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
