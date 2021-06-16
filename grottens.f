      subroutine Grottens(rSpin,rRn, rQ,rRhalf,rRone,rdQdW)

C     Compute the rotation tensors for constitutive models in Jaumann rate form
C     See Simo and Hughes, Box. 8.4
C
C     Given:
C             rSpin: The spin tensor increment
C               rRn: The (old) rotation tensor at t = t_n
C
C     Return:
C            rRhalf: The (updated) rotation tensor at t = t_n+1/2
C             rRone: The (updated) rotation tensor at t = t_n+1
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      include 'vaba_param.inc'

      dimension rSpin(3), rRn(3,3)
      dimension rQ(3,3), rW(3), rdQdW(3,3,3)
      dimension rRhalf(3,3),  rRone(3,3)

      do i = 1, 3
        rW(i) = 0.5d0*rSpin(i)                 ! Q = exp(spin/2)
      enddo
      call Expmap(rW,1, rQ,rdQdW)

      do i = 1, 3
        do j = 1, 3
          rRhalf(i,j) = 0.d0
          do k = 1, 3
            rRhalf(i,j) = rRhalf(i,j) + rQ(i,k)*rRn(k,j)
          enddo
        enddo
      enddo

      do i = 1, 3
        do j = 1, 3
          rRone(i,j) = 0.d0
          do k = 1, 3
            rRone(i,j) = rRone(i,j) + rQ(i,k)*rRhalf(k,j)
            rdQdW(i,j,k) = 0.5d0*rdQdW(i,j,k)  ! Q = exp(spin/2)
          enddo
        enddo
      enddo

      end

*************************************************************************************

      subroutine Expmap(rW,iDflg, rR,rdRdW)

C     Compute the exponential map (see Simo and Hughes Box 8.3)
C     to convert from axis-angle (rW) to rotation matrix (rR)

      include 'vaba_param.inc'

      dimension rW(3), rR(3,3), rdRdW(3,3,3)
      dimension rq(3), rdelta(3,3), rpermut(3,3,3)
      dimension rdRdq(3,3,4), rdqdw(4,3) 

      rWnrm = sqrt( rW(1)**2 + rW(2)**2 + rW(3)**2 )

      rq0   = cos(0.5d0*rWnrm)

      if (rWnrm.gt.1.d-3) then
        rqs = sin(0.5d0*rWnrm) / rWnrm
      else
        rqs = 0.5d0 * ( 1.d0 - (rWnrm**2)/24.d0 + (rWnrm**4)/1920.d0 )
      endif

      do i = 1, 3
       rq(i) = rqs*rW(i)
      enddo

c     Populate Kroenecker's delta and Levi-Civita permutation tensor
      do i = 1, 3
        do j = 1, 3
          if (i.eq.j) then
            rdelta(i,j) = 1.d0
          else
            rdelta(i,j) = 0.d0
          endif
          do k = 1, 3
            rpermut(i,j,k) = 0.5d0*(i-j)*(j-k)*(k-i)
          enddo
        enddo
      enddo

c     Rotation matrix
      do i = 1, 3
        do j = 1, 3
          rR(i,j) = ((rq0**2)-0.5d0)*rdelta(i,j) + rq(i)*rq(j)
          do k = 1, 3
            rR(i,j) = rR(i,j) - rq0*rpermut(i,j,k)*rq(k)
          enddo
          rR(i,j) = 2.d0*rR(i,j)
        enddo
      enddo

c     If we don't need derivatives, return at this point
      if (iDflg.eq.0) RETURN

c     Derivative of quaternion wrt to vector (rW)
      if (rWnrm.gt.1.d-3) then
        do i = 1, 3
          do j = 1, 3
            rdqdw(i,j) = 0.5d0*rW(i)*rW(j)/(rWnrm**2) * cos(0.5d0*rWnrm)
     $                  + ( rdelta(i,j)/rWnrm - rW(i)*rW(j)/(rWnrm**3) )
     $                                                * sin(0.5d0*rWnrm)
          enddo
        enddo
      else
        do i = 1, 3
          do j = 1, 3
            rdqdw(i,j) = rW(i)*rW(j)/24.d0 * ( (rWnrm**2)/40.d0 - 1.d0 )
     $                  + rdelta(i,j) * rqs
          enddo
        enddo
      endif

c     rdqdw(4,i) = dq0/dw  !! T.J.
      do i = 1, 3
       rdqdw(4,i) = -0.5d0*rqs*rW(i)
      enddo

c     Derivative of matrix (rR) wrt quaternion
      do i = 1, 3
        do j = 1, 3
          do k = 1, 3
            rdRdq(i,j,k) =  -rq0*rpermut(i,j,k) 
     $                     + rq(i)*rdelta(j,k) + rq(j)*rdelta(i,k)
          enddo

C         rdRdq(i,j,4) = dR/dq0
          rdRdq(i,j,4) = 2.d0*rq0*rdelta(i,j)
          do k = 1, 3
            rdRdq(i,j,4) = rdRdq(i,j,4) - rq(k)*rpermut(i,j,k)
          enddo

          do k = 1, 4
            rdRdq(i,j,k) = 2.d0*rdRdq(i,j,k)
          enddo
        enddo
      enddo

c     Derivative of matrix wrt to vector (rW) via chain rule
      do i = 1, 3
        do j = 1, 3
          do k = 1, 3
            rdRdW(i,j,k) = 0.d0
            do l = 1, 4
              rdRdW(i,j,k) = rdRdW(i,j,k) + rdRdq(i,j,l)*rdqdw(l,k)   !! Tao
            enddo
          enddo
        enddo
      enddo

      end

*************************************************************************************

      subroutine Logmap(rR, rW)

C     Compute the logarithmic map
C     to convert from rotation matrix (rR) to axis-angle (rW)

      include 'vaba_param.inc'

      dimension rR(3,3)
      dimension rW(3)

      pi = 4.d0*atan(1.d0)

      rW(1) = rR(3,2) - rR(2,3)       ! 2*x*sin(theta)
      rW(2) = rR(1,3) - rR(3,1)       ! 2*y*sin(theta)
      rW(3) = rR(2,1) - rR(1,2)       ! 2*z*sin(theta)

      rWnrm  = sqrt( rW(1)**2 + rW(2)**2 + rW(3)**2 )     ! = 2*sin(theta)
      rTrace = rR(1,1) + rR(2,2) + rR(3,3)                ! = 1+2*cos(theta)

      if (rWnrm.lt.1.d-15) then       ! 2*sin(theta) ~ 0
        do i = 1, 3
          rW(i) = 0.5d0 * rW(i)       ! theta / ( 2*sin(theta) ) --> 0.5
        enddo
        RETURN
      else
        if ( abs(rTrace-1.d0).lt.1.d-15 ) then
          rAngle = pi / 2.d0
        else
          rAngle = atan2(rWnrm,rTrace-1.d0)    ! theta > 0
        endif
        do i = 1, 3
          rW(i) = rAngle * rW(i) / rWnrm       ! theta * unit vector 
        enddo
      endif

      end

*************************************************************************************

      subroutine GrGtan(rQ,rdQdW,rRn,rRhalf,rRone,rScorot,rStress,rCmod,
     $                  rEstat,  rGtan)

      include 'vaba_param.inc'

      dimension rQ(3,3), rdQdW(3,3,3)
      dimension rRn(3,3), rRhalf(3,3), rRone(3,3)
      dimension rScorot(6),rStress(6),rCmod(6,6)
      dimension rEstat(6)
      dimension rGtan(6,3)

      dimension rSold(3,3),rSnew(3,3),rEst(3,3)
      dimension rC(3,3,3,3), rG(3,3,3)
      dimension rdelta(3,3), rZ(3,3,3,3), rX(3,3)
      dimension rER(3,3), rY(3,3,3,3), rRY(3,3,3,3)
      dimension rdSdQ(3,3,3,3), rdSdW(3,3,3)

C     Expand stresses and strains into 3x3 matrices

      k = 6
      do i = 1, 3
        rSold(i,i) = rScorot(i)
        rSnew(i,i) = rStress(i)
        rEst(i,i)  = rEstat(i)
        if (i.eq.3) exit
        do j = (i+1), 3
          rSold(i,j) = rScorot(k)
          rSold(j,i) = rScorot(k)
          rSnew(i,j) = rStress(k)
          rSnew(j,i) = rStress(k)
          rEst(i,j)  = 0.5d0 * rEstat(k)
          rEst(j,i)  = 0.5d0 * rEstat(k)
          k = k - 1
        enddo
      enddo

C     Expand moduli into 3x3x3x3 tensor

      do ij = 1, 6
        if ( ij.le.3 ) then
          i = ij
          j = ij
        elseif ( ij.eq.4 ) then
          i = 2
          j = 3
        elseif ( ij.eq.5 ) then
          i = 1
          j = 3
        else
          i = 1
          j = 2
        endif
        do kl = 1, 6
          if ( kl.le.3 ) then
            k = kl
            l = kl
          else
            if ( kl.eq.4 ) then
              k = 2
              l = 3
            elseif ( kl.eq.5 ) then
              k = 1
              l = 3
            else
              k = 1
              l = 2
            endif
          endif
          rC(i,j,k,l) = rCmod(ij,kl)
          rC(j,i,k,l) = rCmod(ij,kl)
          rC(i,j,l,k) = rCmod(ij,kl)
          rC(j,i,l,k) = rCmod(ij,kl)
        enddo   ! kl
      enddo     ! ij

C     populate delta
      do i = 1, 3
        rdelta(i,i) = 1.d0
        do j = 1, 3
          if (i.ne.j) rdelta(i,j) = 0.d0
        enddo
      enddo

C     Z: Derivative of R_n+1         !! (dR_n+1/dQ)
      do iI = 1, 3
        do iP = 1, 3
          do iQ = 1, 3
            do iC = 1, 3
              rZ(iI,iP,iQ,iC) = rdelta(iI,iP)*rRhalf(iQ,iC)
     $                         +rQ(iI,iP)*rRn(iQ,iC)
            enddo
          enddo
        enddo
      enddo

C     X: New corotational Stress * transpose(R_n+1)
      do iC = 1, 3
        do iJ = 1, 3
          rX(iC,iJ) = 0.d0
          do nD = 1, 3
            rX(iC,iJ) = rX(iC,iJ) + rSnew(iC,nD)*rRhalf(iJ,nD)
          enddo
        enddo
      enddo

C     Y: Derivative of new corotational stress
      do iP = 1, 3
        do iiF = 1, 3
          rER(iP,iiF) = 0.d0
          do nM = 1, 3
            rER(iP,iiF) = rER(iP,iiF) + rEst(iP,nM)*rRhalf(nM,iiF)
          enddo
        enddo
      enddo

      do iC = 1, 3  !! d(sigma_n=1/2)/dQ 
        do iD = 1, 3
          do iP = 1, 3
            do iQ = 1, 3
              rY(iC,iD,iP,iQ) = 0.d0
              do nE = 1, 3
                do nF = 1, 3
                  rY(iC,iD,iP,iQ) =
     $            rY(iC,iD,iP,iQ) +rC(iC,iD,nE,nF)*rRn(iQ,nE)*rER(iP,nF)
                enddo
              enddo
              rY(iC,iD,iP,iQ) = 2.d0 * rY(iC,iD,iP,iQ)
            enddo
          enddo
        enddo
      enddo

      do iI = 1, 3   !! R_n+1 * d(sigma_n=1/2)/dQ
        do iD = 1, 3
          do iP = 1, 3
            do iQ = 1, 3
              rRY(iI,iD,iP,iQ) = 0.d0
              do nC = 1, 3
                rRY(iI,iD,iP,iQ) =   rRY(iI,iD,iP,iQ)
     $                             + rRone(iI,nC)*rY(nC,iD,iP,iQ)
              enddo
            enddo
          enddo
        enddo
      enddo

C     d_Stress/d_Q
      do iI = 1, 3
        do iJ = 1, 3
          do iP = 1, 3
            do iQ = 1, 3
              rdSdQ(iI,iJ,iP,iQ) = 0.d0
              do nD = 1, 3
                rdSdQ(iI,iJ,iP,iQ) =   rdSdQ(iI,iJ,iP,iQ)
     $                               + rZ(iI,iP,iQ,nD)*rX(nD,iJ)
     $                               + rZ(iJ,iP,iQ,nD)*rX(nD,iI)
     $                               + rRY(iI,nD,iP,iQ)*rRone(iJ,nD)  !! R_n+1 * d(sigma_n=1/2)/dQ * R_n+1
              enddo
            enddo
          enddo
        enddo
      enddo

C     d_Stress/d_W
      do iI = 1, 3
        do iJ = 1, 3
          do iS = 1, 3
            rdSdW(iI,iJ,iS) = 0.d0
            do nP = 1, 3
              do nQ = 1, 3
                rdSdW(iI,iJ,iS) =   rdSdW(iI,iJ,iS)
     $                            + rdSdQ(iI,iJ,nP,nQ)*rdQdW(nP,nQ,iS)
              enddo
            enddo
          enddo
        enddo
      enddo

C     (3x3)x3 -> (Voigt-6)x3
      do ij = 1, 6
        if ( ij.le.3 ) then
          i = ij
          j = ij
        elseif ( ij.eq.4 ) then
          i = 2
          j = 3
        elseif ( ij.eq.5 ) then
          i = 1
          j = 3
        else
          i = 1
          j = 2
        endif
        do k = 1, 3
          rGtan(ij,k) = rdSdW(i,j,k)
        enddo
      enddo

      end
*************************************************************************************

      subroutine rSrt(rR,rSfactor, rT)

C     Perform forward tensor rotation R.T.R^T
C
C            rT: symmetric 3x3 tensor in Voigt form
C      rSfactor: shear factor (0.5 if rT has engineering shear strains, 1.0 otherwise)
C

      include 'vaba_param.inc'

      dimension rR(3,3), rT(6)
      dimension rTT(3,3), rTint(3,3)

C     Expand rT into 3x3 matrix

      k = 6
      do i = 1, 3
        rTT(i,i) = rT(i)
        if (i.eq.3) exit
        do j = (i+1), 3
          rTT(i,j) = rSfactor * rT(k)
          rTT(j,i) = rSfactor * rT(k)
          k = k - 1
        enddo
      enddo

C     Left multiply by R

      do i = 1, 3
        do j = 1, 3
          rTint(i,j) = 0.d0
          do k = 1, 3
            rTint(i,j) = rTint(i,j) + rR(i,k)*rTT(k,j)
          enddo
        enddo
      enddo

C     Right multiply by R^T

      do i = 1, 3
        do j = 1, 3
          rTT(i,j) = 0.d0
          do k = 1, 3
            rTT(i,j) = rTT(i,j) + rTint(i,k)*rR(j,k)
          enddo
        enddo
      enddo

C     3x3 matrix -> Voigt form

      k = 6
      do i = 1, 3
        rT(i) = rTT(i,i)
        if (i.eq.3) exit
        do j = (i+1), 3
          rT(k) = rTT(i,j) / rSfactor
          k = k - 1
        enddo
      enddo

      end

*************************************************************************************

      subroutine rtSr(rR,rSfactor, rT)

C     Perform backward tensor rotation R^T.T.R
C
C            rT: symmetric 3x3 tensor in Voigt form
C      rSfactor: shear factor (0.5 if rT has engineering shear strains, 1.0 otherwise)
C

      include 'vaba_param.inc'

      dimension rR(3,3), rT(6)
      dimension rRt(3,3)

C     Get R^T

      do i = 1, 3
        do j = 1, 3
          rRt(i,j) = rR(j,i)
        enddo
      enddo

C     Perform rotation with R replaced by R^T 

      call rSrt(rRt,rSfactor, rT)

      end
      
*************************************************************************************

      subroutine rrrrC(rRhalf, rRone, rC)

C     Perform forward rotation of the tensor of tangent moduli
C
C            rC: 4th order tensor possessing major symmetry, in Voigt form (6x6)
C
C       Note: rC is contracted with a strain tensor in Voigt form
C             containing engineering shear strains.
C

      include 'vaba_param.inc'

      dimension rRhalf(3,3),rRone(3,3), rC(6,6)
      dimension rCC(3,3,3,3), rCCint(3,3,3,3)

C     Expand into 3x3x3x3 tensor

      do ij = 1, 6
        if ( ij.le.3 ) then
          i = ij
          j = ij
        elseif ( ij.eq.4 ) then
          i = 2
          j = 3
        elseif ( ij.eq.5 ) then
          i = 1
          j = 3
        else
          i = 1
          j = 2
        endif
        do kl = 1, 6
          if ( kl.le.3 ) then
            k = kl
            l = kl
          else
            if ( kl.eq.4 ) then
              k = 2
              l = 3
            elseif ( kl.eq.5 ) then
              k = 1
              l = 3
            else
              k = 1
              l = 2
            endif
          endif
          rCC(i,j,k,l) = rC(ij,kl)
          rCC(j,i,k,l) = rC(ij,kl)
          rCC(i,j,l,k) = rC(ij,kl)
          rCC(j,i,l,k) = rC(ij,kl)
        enddo   ! kl
      enddo     ! ij

C     Perform multiplication by R

      do i = 1, 3
        do j = 1, 3
          do k = 1, 3
            do l = 1, 3
              rCCint(i,j,k,l) = 0.d0
              do m = 1, 3
                rCCint(i,j,k,l) =  rCCint(i,j,k,l)
     &                            +rRone(i,m)*rCC(m,j,k,l)
              enddo
            enddo
          enddo
        enddo
      enddo

      do i = 1, 3
        do j = 1, 3
          do k = 1, 3
            do l = 1, 3
              rCC(i,j,k,l) = 0.d0
              do m = 1, 3
                rCC(i,j,k,l) =  rCC(i,j,k,l)
     &                         +rRone(j,m)*rCCint(i,m,k,l)
              enddo
            enddo
          enddo
        enddo
      enddo

      do i = 1, 3
        do j = 1, 3
          do k = 1, 3
            do l = 1, 3
              rCCint(i,j,k,l) = 0.d0
              do m = 1, 3
                rCCint(i,j,k,l) =  rCCint(i,j,k,l)
     &                            +rRhalf(k,m)*rCC(i,j,m,l)
              enddo
            enddo
          enddo
        enddo
      enddo

      do i = 1, 3
        do j = 1, 3
          do k = 1, 3
            do l = 1, 3
              rCC(i,j,k,l) = 0.d0
              do m = 1, 3
                rCC(i,j,k,l) =  rCC(i,j,k,l)
     &                         +rRhalf(l,m)*rCCint(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      enddo

C     3x3x3x3  --> Voigt

      do ij = 1, 6
        if ( ij.le.3 ) then
          i = ij
          j = ij
        elseif ( ij.eq.4 ) then
          i = 2
          j = 3
        elseif ( ij.eq.5 ) then
          i = 1
          j = 3
        else
          i = 1
          j = 2
        endif
        do kl = 1, 6
          if ( kl.le.3 ) then
            k = kl
            l = kl
          else
            if ( kl.eq.4 ) then
              k = 2
              l = 3
            elseif ( kl.eq.5 ) then
              k = 1
              l = 3
            else
              k = 1
              l = 2
            endif
          endif
          rC(ij,kl) = rCC(i,j,k,l)
        enddo   ! kl
      enddo     ! ij

      end
