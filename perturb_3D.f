      subroutine perturb_3D(rStress,iOutPr,rOutPr,rNvect,rMvect)

C     Perturbation analysis for the special case of dynamic adiabatic
C     conditions, with a pressure insensitive rate dependent material

      IMPLICIT NONE

      real*8 rStress(6), rSt(3,3), rNvect(3), rMvect(3)
      real*8 rEigValR(3), rEigValI(3), rEigVec(3,3)
      real*8 rOutPr(15)
      integer iOutPr(15)

      real*8 rRef(3), rtemp(3)

      real*8 rSin45,rCos45,rNorm_max,rNorm_min,rOrth
      integer i,j,iErrCode
      integer imax,imiddle,imin

      rRef(1) = 0.d0
      rRef(2) = 0.d0
      rRef(3) = 1.d0

c     cos(pi/4) = sin(pi/4) = 1/sqrt(2)

      rSin45 = 0.5d0*sqrt(2.d0)
      rCos45 = rSin45

c     Assemble stress tensor (in 2D)

      rSt(1,1) = rStress(1)
      rSt(2,2) = rStress(2)
      rSt(3,3) = rStress(3)

      rSt(1,2) = rStress(6)
      rSt(2,1) = rSt(1,2)

      rSt(1,3) = rStress(5)
      rSt(3,1) = rSt(1,3)

      rSt(2,3) = rStress(4)
      rSt(3,2) = rSt(2,3)

c     Eigenvalue decomposition of stress tensor

      do i = 1, 3
        rEigValR(i) = 0.d0
        rEigValI(i) = 0.d0
        do j = 1, 3
          rEigVec(i,j) = 0.d0
        enddo
      enddo

      call rGEigV(3,3,rSt,rEigValR,rEigValI,1,rEigVec,iErrCode)

c     Order the principal stresses and their associated directions

      do i = 1, 3
        if ((rEigValI(i).ne.0.d0).or.(iErrCode.ne.0)) then
          call XPLB_ABQERR(-3,'Problem finding principal stresses'//
     $     ' by eigenvalue decomposition',0,0.d0,' ')
        endif
      enddo

      if ( rEigValR(1).gt.rEigValR(2) ) then
        if ( rEigValR(1).gt.rEigValR(3) ) then
          imax = 1
          if ( rEigValR(2).gt.rEigValR(3) ) then
            imiddle = 2
            imin = 3
          else
            imiddle = 3
            imin = 2
          end if
        else
          imax = 3
          imiddle = 1
          imin = 2
        end if
      else
        if ( rEigValR(2).gt.rEigValR(3) ) then
          imax = 2
          if ( rEigValR(1).gt.rEigValR(3) ) then
            imiddle = 1
            imin = 3
          else
            imiddle = 3
            imin = 1
          end if
        else 
          imax = 3
          imiddle = 2
          imin = 1 
        end if       
      end if

C     Shear bands are assumed to form at an inclination of 45 degrees
C     with respect to the maximum principal stress direction in the
C     max-min principal stresses plane

C     normalize eigenvectors
      rNorm_max = sqrt( rEigVec(1,imax)**2 + rEigVec(2,imax)**2 + rEigVec(3,imax)**2 )
      do i = 1, 3
        rEigVec(i,imax) = rEigVec(i,imax)/rNorm_max
      end do

      rNorm_min = sqrt( rEigVec(1,imin)**2 + rEigVec(2,imin)**2 + rEigVec(3,imin)**2 )
      do i = 1, 3
        rEigVec(i,imin) = rEigVec(i,imin)/rNorm_min
      end do

C     check whether two eigenvectors are orthogonal
      rOrth = rEigVec(1,imax)*rEigVec(1,imin)
     &      + rEigVec(2,imax)*rEigVec(2,imin)
     &      + rEigVec(3,imax)*rEigVec(3,imin)

      if ( abs(rOrth).GT.1.d-12 ) then
        call XPLB_ABQERR(-3,'Eigenvectors are not orthoganal,'//
     $       'rOrth = %R', 0, rOrth, ' ')
      end if

C     Calculate N
      ! Does '+' or '-' make a difference ? Tao
      do i = 1, 3
        rNvect(i) = ( rEigVec(i,imax) - rEigVec(i,imin) )/sqrt(2.d0)
      end do

C     M is normal to N in the case of a shear band (by definition)
      ! M is not unique, how to decide M
      do i = 1, 3
        rMvect(i) = ( rEigVec(i,imax) + rEigVec(i,imin) )/sqrt(2.d0)
      end do

      if (abs(dot_product(rNvect, rRef)).GT.abs(dot_product(rMvect, rRef))) then
        rtemp = rMvect
        rMvect = rNvect
        rNvect = rtemp
      end if

C     Notification
      rOutPr(3) = rNvect(1)
      rOutPr(4) = rNvect(2)
      rOutPr(5) = rNvect(3)
      rOutPr(6) = rMvect(1)
      rOutPr(7) = rMvect(2)
      rOutPr(8) = rMvect(3)
      call XPLB_ABQERR(1,'PERTURBATION initiated '//
     $   'in element %I, Gauss point %I, '//
     $   'at increment number %I, step time=%R '//
     $   'total time=%R; '//
     $   'n=(%R,%R,%R), m=(%R,%R,%R) ',iOutPr,rOutPr,' ')
      do i = 1, 15
        iOutPr(i) = 0
        rOutPr(i) = 0.d0
      enddo

      end
