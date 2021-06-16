      subroutine hypoJ2(props,nprops,rStrain,rEPStrn,rDt,rT,rSe,rS0,rF,
     $                  iPflag,  iCvFlag,rStress,rCmod)

c     hypoelastic-plastic material model with isotropic hardening
c     
c        rStrain:   the strain increment tensor
c        rEPStrn:   the accumulated effective plastic strain
c            rDt:   Current time step size
c            rT:   local temperature (needed for MTS model)
c           rSe:   local Struct. Evol. Threshold Stress (needed for MTS model)
c         iPflag:   Plasticity flag (0 -> elastic only; 1 -> Mises; 2 -> MTS )
c        
c        rStress:   Cauchy stress tensor
c          rCmod:   Tensor of tangent (elastoplastic) moduli

      include   'vaba_param.inc'

      dimension props(nprops)
      dimension rStrain(6), rStress(6), rCmod(6,6)
      dimension rDevStrn(6), rDevStrs(6), rNormal(6)
      dimension r4IIMatr(6,6), r1x1Matr(6,6), rnxnMatr(6,6)

c     Initialize plastic multiplier

      rDGamma = 0.d0

c     Set parameters and constants

      S23   = sqrt(2.d0/3.d0)

      rE    = props(2)                               ! Young's modulus
      rNu   = props(3)                               ! Poisson's ratio
      rK    = rE / ( 3.d0*( 1.d0 - 2.d0*rNu ) )      ! Bulk modulus
      rG    = rE / ( 2.d0*( 1.d0 +      rNu ) )      ! Shear modulus

      r2G   = 2.d0*rG

      if (abs(iPflag).eq.1) then            ! Mises plasticity
        rTol  = rE*1.d-17
        rTol2 = 1.d-14
      elseif (abs(iPflag).eq.2) then        ! MTS model
        rTol  = rE*1.d-17
        rTol2 = 1.d-14
      endif

      rEPStr0 = rEPStrn                ! store value of equiv. plastic strain from prev. step 

c     Populate identity tensors, etc.

      do i = 1, 6
        do j = 1, 6
          r4IIMatr(i,j) = 0.0d0   ! 4th order symmetric identity tensor
          r1x1Matr(i,j) = 0.0d0   ! 2nd order identity tensor-multiplied by itself
          rnxnMatr(i,j) = 0.0d0   ! rNormal tensor-multiplied by itself
        enddo
      enddo

      do i = 1, 6
        if (i.le.3) then
          r4IIMatr(i,i) = 1.0d0
          do j = 1, 3
            r1x1Matr(i,j) = 1.0d0
          enddo
        else
          r4IIMatr(i,i) = 0.5d0
        endif
      enddo

c     Assemble tensor of ELASTIC moduli

      do i = 1, 6
        do j = 1, 6
          rCmod(i,j) =  rK  *  r1x1Matr(i,j)
     $                + r2G * (r4IIMatr(i,j) - 1.d0/3.d0*r1x1Matr(i,j))
        enddo
      enddo

c     Compute spherical and deviatoric stresses and strains

      rTheta =   rStrain(1) + rStrain(2) + rStrain(3)
      rHydro = ( rStress(1) + rStress(2) + rStress(3) ) / 3.d0

      do i = 1,3
        rDevStrn(i) = rStrain(i) - rTheta/3.d0
        rDevStrs(i) = rStress(i) - rHydro
      end do ! i
      do i = 4,6
        rDevStrn(i) = 0.5d0*rStrain(i)         ! engineering -> tensorial
        rDevStrs(i) =       rStress(i)
      end do ! i

c     Compute (deviatoric) trial stress

      do i = 1,6
        rDevStrs(i) = rDevStrs(i) + r2G*rDevStrn(i)
      end do ! i

c     If we're being called as a hypoelasticity (only) subroutine, skip plasticity calculations

      if (iPflag.eq.0) goto 20

c     Compute norm of trial stress

      rStrsNrm = 0.d0
      do i = 1, 3
        rStrsNrm = rStrsNrm +      rDevStrs(i)*rDevStrs(i)
      enddo
      do i = 4, 6
        rStrsNrm = rStrsNrm + 2.d0*rDevStrs(i)*rDevStrs(i)
      enddo
      rStrsNrm = sqrt(rStrsNrm)

c     Compute flow stress and trial yield function

      if (abs(iPflag).eq.1) then          ! Mises
        call vonmises(props,nprops,rEPStrn, rKflow,rKprime)
      elseif (abs(iPflag).eq.2) then      ! MTS
        rTl  = rT
        rSel = rSe
        rS0l = rS0
        rFl  = rF
        call mtsmodel(props,nprops,rStrsNrm,0.d0,rTl,rSel,rDt,
     $                rS0l,rFl,0,rKflow,rKprime,instFlg)
      endif

      rYieldF = rStrsNrm - S23*rKflow   ! trial yield function

c     Check yield condition: If met ==> Plastic increment

      if ( rYieldF .gt. rTol ) then

c       If we assumed to elastic behavior (unloading), but violated the yield condition

        if (iPflag.lt.0) then
          iCvFlag = -1
          goto 25
        endif

c       Compute normal to yield surface (vector in stress space)

        do i = 1, 6
          rNormal(i) = rDevStrs(i)/rStrsNrm
        end do ! i

c       Compute plastic consistency parameter (rDGamma)

c       A. Initialize

        iConverge = 0
        nTry      = 1
        if (iPflag.eq.1) goto 15      ! Go to Newton's method if rate-independent model

        rDGamm2  = 1.d-25                          ! Right boundary (bracketing methods)
        do nBcoun2 = 1, 10                         ! Sweep left --> right in log space
10        if (iPflag.eq.1) then       ! Mises
            rEPStr2 = rEPStr0 + S23 * rDGamm2
            call vonmises(props,nprops,rEPStr2, rKflo2,rKprim2)
          elseif (iPflag.eq.2) then   ! MTS
            rTl  = rT
            rSel = rSe
            rS0l = rS0
            rFl  = rF
            call mtsmodel(props,nprops,rStrsNrm,rDGamm2,rTl,rSel,rDt,
     &                    rS0l,rFl,0,rKflo2,rKprim2,instFlg)
          endif
          rYF2 = rStrsNrm - r2G*rDGamm2 - S23*rKflo2
          if (rYF2.ge.0.d0) then
            rDGamm1 = rDGamm2                      ! Left boundary (bracketing methods)
            rYF1    = rYF2
            rDGamm2 = rDGamm2*1.d6
            if ( (rStrsNrm - r2G*rDGamm2).lt.0.d0 ) then  ! Negative stress
              rDGamm2 = rStrsNrm/r2G
              do while ( (rStrsNrm - r2G*rDGamm2).lt.0.d0 ) ! Guard against negative stress
                rDGamm2 = rDGamm2 *( 1.d0 - 1.5d-16 )
              end do
              goto 10             ! One more iteration to update rYF2
            endif
          else
            exit
          endif
        enddo
        if (rYF1.lt.0.d0) then  ! Yield function should be positive at rDGamm1
          call XPLB_ABQERR(-2,'J2 routine: FPinit Left Br.',0,0.d0,' ')
          iCvFlag = 0
          RETURN
        endif
        if (rYF2.gt.0.d0) then  ! Yield function should be negative at rDGamm2
          call XPLB_ABQERR(-2,'J2 routine: FPinit Right Br.',0,0.d0,' ')
          iCvFlag = 0
          RETURN
        endif

c       B. Iterate

15      do while ((iConverge.eq.0) .and. (nTry.le.8))

          if (iPflag.eq.1) goto 16      ! Go to Newton's method if rate-independent model

          rBRatio = rDGamm2 / rDGamm1      ! Bisection bracket size
          nBcMax  = int(nTry**1.32)        ! Max number of bisection iterations

          ! Bisection method in log space
          do nBcount = nBcoun2+1, nBcoun2+nBcMax

            rDGamm3 = sqrt(rDGamm1*rDGamm2)
            rEPStr3 = rEPStr0 + S23 * rDGamm3

            if (iPflag.eq.1) then               ! Mises
              call vonmises(props,nprops,rEPStr3, rKflo3,rKprim3)
            elseif (iPflag.eq.2) then           ! MTS
              rTl  = rT
              rSel = rSe
              rS0l = rS0
              rFl  = rF
              call mtsmodel(props,nprops,rStrsNrm,rDGamm3,rTl,rSel,rDt,
     &                      rS0l,rFl,0,rKflo3,rKprim3,instFlg)
            endif

            rYF3 = rStrsNrm - r2G*rDGamm3 - S23*rKflo3   ! update yield function

            if ((abs(rYF3).le.rTol).or.
     &          ((rDGamm2/rDGamm1-1.d0).lt.rTol2)) then
              !write(*,*) rYF3,nBcount
              iConverge = 1
              rDGamma   = rDGamm3
              rEPStrn   = rEPStr3
              if (iPflag.eq.1) then               ! Mises
                rKflow    = rKflo3
                rKprime   = rKprim3
              elseif (iPflag.eq.2) then           ! Call MTS to update flow stress and Jacobian
                rTl  = rT
                rSel = rSe
                rS0l = rS0
                rFl  = rF
                call mtsmodel(props,nprops,rStrsNrm,rDGamm3,rTl,rSel,
     &                        rDt,rS0l,rFl,1,rKflow,rKprime,instFlg)
              endif
              EXIT
            elseif (rYF3.gt.0.d0) then
              if (     ( abs(rDGamm1-rDGamm3).lt.(rDGamm1*rTol2) )
     &            .and.( abs(rYF3).le.(rTol*10.d0) ) ) then
                !write(*,*) rYF3,nBcount
                iConverge = 1
                rDGamma   = rDGamm3
                rEPStrn   = rEPStr3
                if (iPflag.eq.1) then               ! Mises
                  rKflow    = rKflo3
                  rKprime   = rKprim3
                elseif (iPflag.eq.2) then           ! Call MTS to update flow stress and Jacobian
                  rTl  = rT
                  rSel = rSe
                  rS0l = rS0
                  rFl  = rF
                  call mtsmodel(props,nprops,rStrsNrm,rDGamma,rTl,rSel,
     &                          rDt,rS0l,rFl,1,rKflow,rKprime,instFlg)
                endif
                EXIT
              else
                rDGamm1 = rDGamm3
                rYF1    = rYF3
              endif
            elseif (rYF3.lt.0.d0) then
              if (     ( abs(rDGamm2-rDGamm3).lt.(rDGamm2*rTol2) )
     &            .and.( abs(rYF3).le.(rTol*10.d0) ) ) then
                !write(*,*) rYF3,nBcount
                iConverge = 1
                rDGamma   = rDGamm3
                rEPStrn   = rEPStr3
                if (iPflag.eq.1) then               ! Mises
                  rKflow    = rKflo3
                  rKprime   = rKprim3
                elseif (iPflag.eq.2) then           ! Call MTS to update flow stress and Jacobian
                  rTl  = rT
                  rSel = rSe
                  rS0l = rS0
                  rFl  = rF
                  call mtsmodel(props,nprops,rStrsNrm,rDGamma,rTl,rSel,
     &                          rDt,rS0l,rFl,1,rKflow,rKprime,instFlg)
                endif
                EXIT
              else
                rDGamm2 = rDGamm3
                rYF2    = rYF3
              endif
            endif

            if ( (rDGamm2/rDGamm1).le.(rBRatio*1.d-3) ) EXIT  ! Bracket got much smaller. Let's try Newton

          enddo    ! Bisection
          if (rYF1.lt.0.d0) then  ! Yield function should be positive at rDGamm1
            call XPLB_ABQERR(-2,'J2 routine: FPiter Left Br.',
     &                       0,0.d0,' ')
            iCvFlag = 0
            RETURN
          endif
          if (rYF2.gt.0.d0) then  ! Yield function should be negative at rDGamm2
            call XPLB_ABQERR(-2,'J2 routine: FPiter Right Br.',
     &                       0,0.d0,' ')
            iCvFlag = 0
            RETURN
          endif
        
          ! Newton-Raphson Method
          if (iConverge.eq.0) then
            rDGamma = (rYF2*rDGamm1-rYF1*rDGamm2)/(rYF2-rYF1)    ! initial guess (Newton)
            rEPStrn = rEPStr0 + S23 * rDGamma
          endif
16        iCount    = 0
          do while ((iConverge.eq.0) .and. (iCount.lt.7))

            if (iPflag.eq.1) then               ! Mises
              call vonmises(props,nprops,rEPStrn, rKflow,rKprime)
            elseif (iPflag.eq.2) then           ! MTS
              rTl  = rT
              rSel = rSe
              rS0l = rS0
              rFl  = rF
              call mtsmodel(props,nprops,rStrsNrm,rDGamma,rTl,rSel,rDt,
     &                      rS0l,rFl,1,rKflow,rKprime,instFlg)
            endif

            rYieldF = rStrsNrm - r2G*rDGamma - S23*rKflow   ! update yield function
            rDDGamm = rYieldF / (r2G + 2.d0/3.d0*rKprime)   ! increment of rDGamma 

            if ( (abs(rYieldF).le.rTol).or.
     &           ( (abs(rYieldF).le.(10.d0*rTol)).and.
     &             (abs(rDDGamm).le.(rDGamma*rTol2)) ) ) then
              !write(*,*) rYieldF,nBcount,(iCount+1)
              iConverge = 1
              if (rDGamma.lt.0.d0) iConverge = 0
              exit
            endif

c           update consistency parameter and equiv. plastic strain
            rDGamma = rDGamma + rDDGamm
            rEPStrn = rEPStr0 + S23 * rDGamma

            if (rDGamma.lt.0.d0) then
              exit
            endif

            iCount  = iCount + 1

          end do ! Newton iteration

          nTry = nTry + 1

        enddo ! while (no convergence) and (nTry < limit)

c       If not converged => Bail out

17      if (iConverge.eq.0) then
          call XPLB_ABQERR(-3,'J2: No Convergence',0,0.d0,' ')
          iCvFlag = 0
          RETURN
        endif

c       Add plastic correction to trial stress

        do i = 1, 6
          rDevStrs(i) = rDevStrs(i) - r2G * rDGamma * rNormal(i)
        enddo

c       Add extra terms to rCmod, to get tensor of ELASTO-PLASTIC tangent moduli

        do i = 1, 6
          do j = 1, 6
            rnxnMatr(i,j) = rNormal(i) * rNormal(j)
            rCmod(i,j) =  rCmod(i,j)
     &                  - r2G * ( 1.d0 / ( 1.d0 + rKprime/(3.d0*rG) ) )
     &                        * rnxnMatr(i,j)
     &       - r2G * rDGamma * (r2G/rStrsNrm)
     &       * (r4IIMatr(i,j) - 1.d0/3.d0*r1x1Matr(i,j) - rnxnMatr(i,j))
          end do ! j
        end do ! i

      endif    ! Plastic increment

c     We can only be here if we have achieved convergence

20    iCvFlag = 1

c     Assemble full stress tensor

25    do i = 1, 3
        rStress(i)  = rDevStrs(i) + rHydro + rK * rTheta     ! add old pressure and pressure increment
      enddo
      do i = 4, 6
        rStress(i)  = rDevStrs(i)
      enddo

c     If MTS model, update temperature and Threshold stress

      if (iPflag.eq.2) then
        rT   = rTl
        rSe  = rSel
        rS0  = rS0l
        rF   = rFl
        rEdot= S23*rDGamma/rDt

        if (instFlg.eq.1) then    ! material is unstable
          iCvFlag = 3
        endif

      endif

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine vonmises(props,nprops,rEPStrn, rKflow,rKprime)

C        i      |       PROPS(i)
C---------------|----------------------------------------------
C       34      |       Y0: Initial yield stress 
C       35      |       H0: Initial linear isotropic hardening modulus
C       36      |      Y1: Second value of yield stress
C       37      |      H1: Linear isotropic hardening modulus at intermediate strain
C       38      |      Y2: Yield stress at large strains

C Note: The Mises yield function has three hardening regimes and is given by:
C       Y = Y0 + H0* rEPStrn           if       rEPStrn <  A1
C       Y = Y1 + H1*(rEPStrn-A1)       if A1 <= rEPStrn <  A2
C       Y = Y2                         if A2 <= rEPStrn
C
C       where A1 = (Y1-Y0)/H0, A2 = A1 + (Y2-Y1)/H1

      include   'vaba_param.inc'

      dimension props(nprops)

c     Set parameters and constants

      Y0    = props(34)
      H0    = props(35)
      Y1    = props(36)
      H1    = props(37)
      Y2    = props(38)

      A1    =      (Y1-Y0)/H0
      A2    = A1 + (Y2-Y1)/H1

c     Calculate flow stress and hardening modulus

      if (rEPStrn.lt.A1) then
        rKflow = Y0 + H0 * rEPStrn          ! flow stress
        rKprime= H0                         ! isotropic hardening modulus
      elseif (rEPStrn.lt.A2) then
        rKflow = Y1 + H1 *(rEPStrn-A1)
        rKprime= H1
      else
        rKflow = Y2
        rKprime= 1.0d-14
      endif

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mtsmodel(props,nprops,rStrsNrm,rDGamma,rTn,rSen,rDt,
     &                    rS0,rF,kpflg, rKFlow,rKprime,instFlg)


      include   'vaba_param.inc'

      dimension props(nprops)

C     Parameters and material constants

      rRho     = props(1)
      rE       = props(2)                               ! Young's modulus
      rNu      = props(3)                               ! Poisson's ratio

      rTinit   = props(15)                              ! Initial temp
      rMu0     = props(41)
      rD0      = props(42)
      rT0      = props(43)                              ! Reference temp
      rK       = props(44)                              ! Boltzmann's const
      rB       = props(45)                              ! Burger's vec mag
      rEdot0i  = props(46)
      rQi      = props(47)
      rPi      = props(48)
      rG0i     = props(49)
      rEdot0e  = props(50)
      rQe      = props(51)
      rPe      = props(52)
      rG0e     = props(53)
      rSa      = props(54)
      rSi      = props(55)
      rSeinit  = props(56)
      rEdot0es = props(57)
      rG0es    = props(58)
      rSes0    = props(59)
      rA0      = props(60)
      rA1      = props(61)
      rA2      = props(62)
      rA3      = props(63)
      rA4      = props(64)
      rA5      = props(65)
      rKappa   = props(66)
      rCp0     = props(67)
      rCp1     = props(68)
      rCp2     = props(69)
      rPsi     = props(70)
      rTmelt   = props(71)
      rTdrx    = props(72)*rTmelt
      rD0drx   = props(73)
      rQRdrx   = props(74)
      rddrx    = props(75)

      S23      = sqrt(2.d0/3.d0)
      rG       = rE / ( 2.d0*( 1.d0 + rNu ) )           ! Shear modulus (elastic constant)
      r2G      = 2.d0*rG
      rEdotmax = min( rEdot0i, rEdot0e, rEdot0es )

      rMuMin   = rMu0 - rD0 / ( exp(rT0/rTmelt) - 1.d0 )
      rEdotmini= rEdot0i * exp( (-rMuMin*(rB**3)*rG0i)/(rK*rTmelt) )
      rEdotmine= rEdot0e * exp( (-rMuMin*(rB**3)*rG0e)/(rK*rTmelt) )
      rEdotmin = max(rEdotmini,rEdotmine)

C     Check for (some) errors which will lead to div. by zero

      if (rRho.le.0.d0) then
        goto 100
      endif
      if (rTinit.le.0.d0) then
        goto 100
      endif
      if (rTn.lt.rTinit) then
        goto 100
      endif
      if (rMu0.le.0.d0) then
        goto 100
      endif
      if (rB.le.0.d0) then
        goto 100
      endif
      if (rEdot0i.le.0.d0) then
        goto 100
      endif
      if (rQi.le.0.d0) then
        goto 100
      endif
      if (rPi.le.0.d0) then
        goto 100
      endif
      if (rG0i.le.0.d0) then
        goto 100
      endif
      if (rEdot0e.le.0.d0) then
        goto 100
      endif
      if (rQe.le.0.d0) then
        goto 100
      endif
      if (rPe.le.0.d0) then
        goto 100
      endif
      if (rG0e.le.0.d0) then
        goto 100
      endif
      if (rEdot0es.le.0.d0) then
        goto 100
      endif
      if (rG0es.le.0.d0) then
        goto 100
      endif
      if (rDt.eq.0.d0) then
        goto 100
      endif

C     Plastic strain rate

      rEdot    = S23*rDGamma/rDt

      meflg    = 0
      if (rEdot.lt.rEdotmin) then
        rEdot  = rEdotmin
        meflg  = 1
      elseif (rEdot.gt.rEdotmax) then
        rEdot  = rEdotmax
        meflg  = 1
      endif

      if (kpflg.eq.0) goto 10    ! skip derivatives

      ! Derivatives
      if (meflg.eq.0) then
        rdEdotdG = S23/rDt
      else
        rdEdotdG = 0.d0
      endif

C     Temperature 

10    rPDis    = (rStrsNrm-r2G*rDGamma)*rDGamma                ! Plastic dissipation
      if (rPDis.lt.0.d0) then
        goto 100                ! ERROR
      endif

      k        = 0
      rT       = rTn
      rTinc    = 100.d0
      do while ( (k.lt.15).and.(abs(rTinc).gt.1.d-15) )
        rCp    = rCp0 + rCp1*rT + rCp2/(rT**2)                 ! 1-point iteration
        if (rCp.le.0.d0) then
          goto 100              ! ERROR
        endif
        rTk    = rT
        rT     = rTn + rPsi/(rRho*rCp) * rPDis                 ! implicit update
        rTinc  = (rT - rTk)/rTk
        k      = k + 1
      enddo
c      if (abs(rTinc).gt.1.d-15) then
c        call XPLB_ABQERR(-2,'MTS: Temp. not converged',0,0.d0,' ')
c      endif

      mtflg    = 0
      if (rT.gt.rTmelt) then
        rT     = rTmelt
        mtflg  = 1
      endif

      if (kpflg.eq.0) goto 20    ! skip derivatives

      ! Derivatives
      if (mtflg.eq.0) then
        rpTpG  = rPsi/(rRho*rCp) * (rStrsNrm-4.d0*rG*rDGamma)  ! Partial deriv.
        rpTpCp =-rPsi/(rRho*rCp**2) * rPDis
        rdCpdT = rCp1 - 2.d0*rCp2/rT**3
        rdTdG  = rpTpG / ( 1.d0 - rpTpCp*rdCpdT )              ! Total deriv.
      else
        rdTdG  = 0.d0
      endif

C     Shear modulus (Temperature dependent)

20    rMu     = rMu0 - rD0 / ( exp(rT0/rT) - 1.d0 )

      if (rMu.le.0.d0) then
        goto 100   ! ERROR
      endif

      if (kpflg.eq.0) goto 30    ! skip derivatives

      ! Derivatives
      rdMudT  = -rD0*rT0*exp(rT0/rT) / ((rT**2)*(exp(rT0/rT)-1.d0)**2)
      rdMudG  = rdMudT * rdTdG

C     Rate and temperature dependent terms

c     1. Intrinsic

30    rSSi2    = rK*rT/(rMu*(rB**3)*rG0i)*log(rEdot0i/rEdot)
      rSSi1    = rSSi2**(1.d0/rQi)
      rSSi     = ( 1.d0 - rSSi1 )**(1.d0/rPi)

      if (kpflg.eq.0) goto 40    ! skip derivatives

      ! Derivatives of Si(Edot,T)
      if (rEdot.lt.rEdot0i) then
        rdSSi1   = rSSi * rSSi1 / ( rPi*rQi * (1.d0-rSSi1) )

        rdSSidEd = rdSSi1 / ( rEdot*log(rEdot0i/rEdot) )
        rdSSidT  =-rdSSi1 / rT
        rdSSidMu = rdSSi1 / rMu

        rdSSidG  = rdSSidEd*rdEdotdG + rdSSidT*rdTdG + rdSSidMu*rdMudG
      else
        rdSSidG  = 0.d0
      endif

c     2. Structure dependent

40    rSSe2    = rK*rT/(rMu*(rB**3)*rG0e)*log(rEdot0e/rEdot)
      rSSe1    = rSSe2**(1.d0/rQe)
      rSSe     = ( 1.d0 - rSSe1 )**(1.d0/rPe)

      if (kpflg.eq.0) goto 50    ! skip derivatives

      ! Derivatives of Se(Edot,T)
      if (rEdot.lt.rEdot0e) then
        rdSSe1   = rSSe * rSSe1 / ( rPe*rQe * (1.d0-rSSe1) )

        rdSSedEd = rdSSe1 / ( rEdot*log(rEdot0e/rEdot) )
        rdSSedT  =-rdSSe1 / rT
        rdSSedMu = rdSSe1 / rMu

        rdSSedG  = rdSSedEd*rdEdotdG + rdSSedT*rdTdG + rdSSedMu*rdMudG
      else
        rdSSedG  = 0.d0
      endif

      ! Dynamic recrystallization (DRX)

50    if (rTdrx.gt.0.d0) then      ! if we're considering DRX
        if (rT.ge.rTdrx) then      ! and DRX has begun

          if (rF.lt.0.d0) then       ! if we've just reached the DRX threshold
            rF  = 1.d0
            rS0 = rSen
          else
            rUn   = ( (1.d0-rF)*(rddrx**3) )**(1.d0/3.d0)
            rUdot = (rD0drx*rB**2)/(rK*rT) * exp(-rQRdrx/rT) * rS0
            rU    = rUn + rDt * rUdot

            if (rU.lt.rddrx) then           ! DRX still in progress
              rF    = 1.d0 - (rU/rddrx)**3
              rdFdU = -3.d0*(rU**2)/(rddrx**3)
              rdUdT = rDt*(rUdot/rT)*(rQRdrx/rT-1.d0)
            else                            ! DRX done
              rF    = 0.d0
              rdFdU = 0.d0
              rdUdT = 0.d0
            endif
          endif   ! Just reached DRX thershold

          rSe = rF*rS0
          if (kpflg.eq.0) then
            goto 60    ! skip derivatives and strain hardening calculations
          else
            rSeRate = 0.d0                        ! for material stability analysis
            rdFdG   = rdFdU * rdUdT * rdTdG       ! for Jacobian
            goto 55    ! skip strain hardening calculations
          endif

        endif   ! DRX has begun
      endif   ! DRX considered

      ! The evolving structure-dependent stress

      rSes     = rSes0*(rEdot/rEdot0es)**(rK*rT/(rMu*(rB**3)*rG0es))
      rh0      = rA0+rA1*log(rEdot)+rA2*sqrt(rEdot)-rA3*rT+rA4/(rT**rA5)

      rSe      = rSen
      do k = 0, 15
        if ((rSe.ge.rSes).or.(rh0.le.0.d0)) then
          rSe    = rSen                                       ! Saturated: no update
          rdSedG = 0.d0                                       ! hence no dependence on rdGamma
          rSeRate= 0.d0
          goto 55                                             ! So skip derivative calculations
        endif
        rSek = rSe
        rSeRate= rh0 * (1.d0 - rSe/rSes)**rKappa              ! Implicit update
        rSeRes = rSe - (rSen + rSeRate * S23*rDGamma)
        if (abs(rSeRes).le.1.d-10) exit
        rSeJac = 1.d0 + S23*rDGamma * rKappa*rSeRate /
     $                   (rSes * (1.d0-rSe/rSes))
        rSe    = rSek - rSeRes/rSeJac
        if (rSek.gt.0.d0) then
          rSeinc = (rSe - rSek)/rSek
        else
          rSeinc = rSe - rSek
        endif
        if (abs(rSeinc).le.1.d-15) exit
      enddo
C      if ( (abs(rSeRes).gt.1.d-10).and.(abs(rSeinc).gt.1.d-15) ) then
C        call XPLB_ABQERR(-2,'MTS: Se not converged',0,0.d0,' ')
C      endif

      if (kpflg.eq.0) goto 60    ! skip derivatives

      ! Derivatives of the evolving stress
      rdh0dEd = rA1/rEdot + 0.5d0*rA2/sqrt(rEdot)
      rdh0dT  =-rA3 - rA4*rA5/(rT*rT**rA5)
      rdh0dG  = rdh0dEd*rdEdotdG + rdh0dT*rdTdG

      rdSes1  = rSes*rK/(rB**3*rG0es)
      rdSesdEd= rdSes1*rT/(rMu*rEdot)
      rdSesdT = rdSes1*log(rEdot/rEdot0es)/rMu
      rdSesdMu=-rdSes1*rT*log(rEdot/rEdot0es)/rMu**2
      rdSesdG = rdSesdEd*rdEdotdG + rdSesdT*rdTdG + rdSesdMu*rdMudG

      rdRate13= rdh0dG*(1.d0-rSe/rSes)**rKappa                       ! Term 1: d/dh0  
     $         +rdSesdG*rKappa*rSeRate/(1.d0-rSe/rSes)*(rSe/rSes**2) ! Term 3: d/dSes 
      rdRate2 = rKappa*rSeRate / (rSes * (1.d0-rSe/rSes))            ! Term 2: d/dSe

      rdSeNUM = rdRate13*S23*rDGamma + rSeRate*S23
      rdSeDEN = 1.d0 + rdRate2*S23*rDGamma
      rdSedG  = rdSeNUM / rdSeDEN

c     Check for material instability

55    rdKdT   =  rdMudT/rMu0*( rSSi*rSi + rSSe*rSe )
     $          +   rMu/rMu0*( (rdSSidT+rdSSidMu*rdMudT)*rSi
     $                        +(rdSSedT+rdSSedMu*rdMudT)*rSe )

      rdKde   =  rMu/rMu0 * rSSe*rSeRate

      rStabnum= rPsi*rKFlow*rdKdT+rRho*rCp*rdKde

      if ( rStabnum.lt.0.d0 ) then
        instFlg = 1        ! Material has become unstable
!       call XPLB_ABQERR(1,'MTS stability number = %R ',0,rStabNum,' ')
      else
        instFlg = 0
      endif

c     Flow stress

60    rKFlow  = rSa + rMu/rMu0*( rSSi*rSi + rSSe*rSe )

      if (kpflg.eq.0) goto 70    ! skip Jacobian

c     Jacobian (derivative wrt increment of plastic strain)

      rdKdG   =  rMu/rMu0   *( rdSSidG*rSi + rdSSedG*rSe + rSSe*rdSedG)
     $          +rdMudG/rMu0*(    rSSi*rSi +    rSSe*rSe )

      rKprime = rdKdG / S23

c     Update values of temp and Se

70    rTn     = rT
      rSen    = rSe

      RETURN

100   call XPLB_ABQERR(-3,'Error in MTS model. Check input!',0,0.d0,' ')

      END

