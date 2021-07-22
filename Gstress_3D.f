      subroutine Gstress_3D(rStrain,rSpin,svars,nstatev,props,
     $                      nprops,iVelmType,
     $                      iVforType,iConsType,rTime,rDt,
     $                      rStress,rSigInc,rCmod,rGtan,rCd,iCvFlag)

C     Compute the Cauchy stress tensor at (Gauss) point in 8-node brick element

C     Finite deformation only

C     Given:
C             rStrain: Strain increment
C               rSpin: Spin increment
C               svars: State/history variables at current (Gauss) point
C             nstatev: Number of state/history variables per (Gauss) point

C           iVelmType: 3D element type:  1 --> 3D Brick (8-node) Elmt

C           iVforType: Formulation type: 0 --> displacement based
C                                        1 --> displacement based with localization behavior
C                                        2 --> B-bar incompressible
C                                        3 --> B-bar incompressible with localization behavior
C           iConsType: Constitution:     0 --> isotropic elastic
C                                        1 --> isotropic J2 plasticity
C                                        2 --> TEPLA model
C                                        3 --> Brown-Kim-Anand model
C                      in Jaumman rate (of Cauchy stress) form 
C               rTime: Step time = rTime(1), Total time = rTime(2)
C                 rDt: Current time increment
C
C
C     Return the following, evaluated at the given (Gauss) point:
C
C             rStress: Cauchy stress at given (quadrature) point
C               rCmod: Tangent moduli (d_Stress/d_Strain_increment)
C               rGtan: d_Stress/d_Spin_increment (finite def. case only)
C                 rCd: Dilatational wave speed
C
C
C     Notes:  When localization behavior is considered, the current subroutine is 
C             called twice for each point: Once for the Localized and once for
C             the Non-localized stress. In this case, the value of "nstatev" passed
C             into this subroutine should be the number of variables for the appropriate
C             (i.e. either the Localized or the Non-localized) material model, NOT the
C             total number of variables per (Gauss) point
C
C             It is assumed that the pre-existing stress is passed into this subroutine
C             in svars(i), with i varying from 1 to 6.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      include 'vaba_param.inc'

      dimension rStrain(6), rSpin(3), svars(nstatev), props(nprops)
      dimension rTime(2)
      dimension rW(3), rRn(3,3), rRhalf(3,3), rRone(3,3)
      dimension rQ(3,3), rdQdW(3,3,3)
      dimension rCmod(6,6), rGtan(6,3), rStress(6)
      dimension rEstat(6),rScorot(6)

      dimension rSigInc(6)    ! stress increment

C     local variables to be passed to constitutive update subroutines:

      character*80 cmname        ! material name string
      dimension rLcoorMp(3)      ! material point coordinates
      dimension rLprops(100)     ! re-ordered copy of PROPS array
      dimension rLstrInc(6)      ! Strain increment tensor
      dimension rLspInc(3)       ! (dummy) relative spin increment
      dimension rLUOld(6)        ! old strecth tensor
      dimension rLFOld(9)        ! old deformation gradient tensor
      dimension rLfldOld(10)     ! old field variables array
      dimension rLUNew(6)        ! new strecth tensor
      dimension rLFNew(9)        ! new deformation gradient tensor
      dimension rLfldNew(10)     ! new field variables array
      dimension rLSOld(6)        ! old stress tensor
      dimension rLstOld(nstatev) ! old state variables
      dimension rLSNew(6)        ! new stress tensor
      dimension rLstNew(nstatev) ! new state variables

C     initialize

      do i = 1, 6
        do j = 1, 6
          rCmod(i,j) = 0.d0
        enddo
        do j = 1, 3
          rGtan(i,j) = 0.d0
        enddo
      enddo

C     Finite rotations
      if ( rTime(2).eq.0.d0 ) then      ! Initialize rotation tensor
        do i = 1, 3
          rRn(i,i) = 1.d0
          if (i.eq.3) exit
          do j = (i+1), 3
            rRn(i,j) = 0.d0
            rRn(j,i) = 0.d0
          enddo
        enddo
      else
        j = nstatev - 3
        do i = 1, 3
          j = j + 1
          rW(i) = svars(j)     ! material rotation from t = t_n
        enddo
        call Expmap(rW,0, rRn,rdQdW)               ! rotation tensor R_n
      endif

      call Grottens(rSpin,rRn, rQ,rRhalf,rRone,rdQdW) ! rotation tensors R_n+1/2 and R_n+1
      call Logmap(rRone, rW)                       ! rotation "vector" at t_n+1
      j = j - 3
      do i = 1, 3
        j = j + 1
        svars(j) = rW(i)     ! write updated material rotation to history
      enddo

      rRho= props(1)          ! mass density
      rE  = props(2)          ! Young's modulus
      rNu = props(3)          ! Poisson's ratio

      if (iConsType .eq. 0) then        ! Elastic
        if ( (props(34).le.0.d0).and.(props(41).le.0.d0) ) then   ! simply hypo-elastic
          iPflag = 0
        elseif ( props(41).le.0.d0 ) then  ! Check yielding (Linear isotropic hardening)
          iPflag = -1
        else                               ! Check yielding (MTS)
          iPflag = -2
        endif
      elseif (iConsType .eq. 1) then
        if ( props(41).le.0.d0 ) then   ! Linear isotropic hardening
          iPflag = 1
        else                            ! MTS
          iPflag = 2
        endif
      endif

      if ( (iConsType .eq. 0).or.(iConsType .eq. 1) ) then
        if ( rTime(2).eq.0.d0 ) then      ! Initialize
          do i = 1, 6
            rStress(i) = 0.d0
          enddo
          rEPStrn = 0.d0 ! rEPStrn: the accumulated effective plastic strain
          if (abs(iPflag).eq.2) then
            T     = props(15)                          ! Initial temperature
            Se    = props(56)                          ! Initial value of rSe 
            S0    = 0.d0
            rF    = -1.d0
            if ( (T.le.0.d0).or.(Se.lt.0.d0) ) then
              call XPLB_ABQERR(-3,'Error in MTS initial data!',
     $                             0,0.d0,' ')
            endif
          endif
        else
          do i = 1, 6
            rStress(i) = svars(i)
          enddo
          if (abs(iPflag).eq.1) then            ! Mises with linear isotropic hardening
            rEPStrn = svars(7)
            T       = 0.d0
            Se      = 0.d0
            S0      = 0.d0
            rF      = 0.d0
          elseif (abs(iPflag).eq.2) then        ! MTS model
            rEPStrn = svars(7)     ! accumulated effective plastic strain
            T       = svars(8)     ! temperature
            Se      = svars(9)     ! Struct. Evol. Threshold Stress (MTS)
            S0      = svars(10)    ! Struct. Evol. Thres. Stress at DRX init
            rF      = svars(11)    ! Vol. Fraction of material NOT yet DRXed
          else
            rEPStrn = 0.d0
            T       = 0.d0
            Se      = 0.d0
            S0      = 0.d0
            rF      = 0.d0
          endif
        endif
      endif

C     Finite deformation case

C     1) Rotate the strain increment into the co-rotational frame

      do i = 1, 6
        rEstat(i) = rStrain(i)      ! Record strain increment in stationary frame
      enddo

C     rStrain is in corotational frame
      call rtSr(rRhalf,0.5d0, rStrain)

      do i = 1, 6
        rScorot(i) = rStress(i)     ! Record stress (at t=t_n) in corotational frame
      enddo

C     2) Compute the new stress state

      if (iConsType .eq. 0) then       ! isotropic elastic

        call hypoJ2(props,nprops,rStrain,rEPStrn,rDt,T,Se,S0,rF,
     &              iPflag, iCvFlag,rStress,rCmod)

        do i = 1, 6
          svars(i)   = rStress(i) ! stress in corotational frame
        enddo

C       Compute the current dilatational wave speed

        rCd = sqrt( rCmod(1,1) / rRho )

      elseif (iConsType .eq. 1) then   ! isotropic J2 plasticity

        do i = 1, 6
          rSigInc(i) = -rStress(i)               ! negative of old stress
        enddo

        if ( iPflag.eq.1 ) then        ! Linear isotropic hardening

          if (rTime(2).eq.0.d0) iPflag = 0   ! purely elastic calculation.
          call hypoJ2(props,nprops,rStrain,rEPStrn,rDt,T,Se,S0,rF,
     &                iPflag, iCvFlag,rStress,rCmod)

          

        elseif ( iPflag.eq.2 ) then    ! MTS

          if (rTime(2).eq.0.d0) iPflag = 0   ! purely elastic calculation.
          call hypoJ2(props,nprops,rStrain,rEPStrn,rDt,T,Se,S0,rF,
     &                iPflag, iCvFlag,rStress,rCmod)
          svars(8) = T
          svars(9) = Se
          svars(10)= S0
          svars(11)= rF
        endif

        do i = 1, 6
          svars(i)   = rStress(i)
          rSigInc(i) = rSigInc(i) + rStress(i)   ! negative of old stress + new stress
        enddo
        svars(7)     = rEPStrn

C       Compute the dilatational wave speed

        call Gwavespd(rStrain,rSigInc,rE,rNu,rRho, rCd)

        if ( props(41).gt.0.d0 ) rCd=1.3d0*rCd   ! be a bit extra conservative with MTS
            

      elseif (iConsType .eq. 2) then   ! TEPLA model
        call XPLB_ABQERR(-3,'TEPLA material model has not '//
     $                      'yet been implemented.',0,0.d0,' ')

      elseif (iConsType .eq. 3) then   ! Brown-Kim-Anand model

C       Prepare local variables to be passed to constitutive update subroutine

        nLblock     = 1             ! Process one material point in this call (Do not change!)
        nLdir       = 3             ! number of normal components in stress/strain tensors
        nLshr       = 3             ! number of shear components in stress/strain tensors
        nLstatev    = nstatev - 6   ! number of state variables
        nLfieldv    = 1             ! number of field variables (.le. size of rLfldOld, rLfldNew)
        nLprops     = 20            ! number of properties passed in via rLprops (see below)
        nLanneal    = 0             ! Annealing process? i.e. reinitialize state variables? (0=no,1=yes)

        rLstepTime  = rTime(1)      ! Time since beginning of current step
        rLtotalTime = rTime(2)      ! Total time

        rLprops(1)  = props(11)     ! PLSLMT  limit on equiv. plastic tensil strain increment
        rLprops(2)  = props(12)     ! PHI     degree of implicitness [0,1]
        rLprops(3)  = props(13)     ! OMEGA   fraction of plastic work to heating
        rLprops(4)  = props(14)     ! S0      initial value for internal variable S
        rLprops(5)  = props(15)     ! T0      initial temperature
        rLprops(6)  = rE / ( 2.d0*(1+     rNu) )  ! AMU        shear modulus
        rLprops(7)  = rE / ( 3.d0*(1-2.d0*rNu) )  ! AKAPPA     bulk modulus
        rLprops(8)  = props(16)     ! ALPHA   thermal expansion coefficient
        rLprops(9)  = props(1)      ! RHO     density
        rLprops(10) = props(17)     ! C       specific heat
        rLprops(11) = props(18)     ! A       pre-exponential factor
        rLprops(12) = props(19)     ! Q/R     activation energy divided by gas constant
        rLprops(13) = props(20)     ! XP      exponent in equation for F (1/n)
        rLprops(14) = props(21)     ! PALPHA  pressure sensitivity parameter
        rLprops(15) = props(22)     ! H0      pre-multiplier in hardening rate equation
        rLprops(16) = props(23)     ! SB      pre-multiplier for saturation value of S
        rLprops(17) = props(24)     ! AN      exponent in equation for saturation value of S
        rLprops(18) = props(25)     ! B       pre-multiplier for rate of static restoration
        rLprops(19) = props(26)     ! SA      annealed value of S
        rLprops(20) = props(27)     ! P       exponent in equation for RDOT

        rLdensity   = props(1)      ! density

        do i = 1, 3
          rLstrInc(i) = rStrain(i)            ! e_ii  (no sum)
        enddo
        rLstrInc(4) = 0.5d0 * rStrain(6)      ! e_12
        rLstrInc(5) = 0.5d0 * rStrain(4)      ! e_23
        rLstrInc(6) = 0.5d0 * rStrain(5)      ! e_13

        do i = 1, 3
          rLSOld(i) = rStress(i)              ! S_ii  (no sum)
        enddo
        rLSOld(4)   = rStress(6)              ! S_12
        rLSOld(5)   = rStress(4)              ! S_23
        rLSOld(6)   = rStress(5)              ! S_13

        rLstOld(1)  = svars(7)      ! S                plastic flow resistance, tensile
        rLstOld(2)  = svars(8)      ! TH               temperature
        rLstOld(3)  = svars(9)      ! F                plastic shear strain rate
        rLstOld(4)  = svars(10)     ! time int of F    plastic shear strain

        do i = 1, 6
          rSigInc(i) = -rStress(i)            ! negative of old stress
        enddo

C       Call constitutive update subroutine

        call vumbka(
     1    nLblock,nLdir,nLshr,nLstatev,nLfieldv,nLprops,nLanneal,
     2    rLstepTime,rLtotalTime,rDt,cmname,rLcoorMp,rLElmtSz,
     3    rLprops,rLdensity,rLstrInc,rLspInc,
     4    rLTOld,rLUOld,rLFOld,rLfldOld,
     5    rLSOld,rLstOld,rLeIntOld,rLeInlOld,
     6    rLTNew,rLUNew,rLFNew,rLfldNew,
     7    rLSNew,rLstNew,rLeIntNew,rLeInlNew)

C       Return/store stress and state variables in correct order

        do i = 1, 3
          rStress(i)  = rLSNew(i)             ! S_ii  (no sum)
        enddo
        rStress(4)  = rLSNew(5)               ! S_23
        rStress(5)  = rLSNew(6)               ! S_13
        rStress(6)  = rLSNew(4)               ! S_12

        do i = 1, 3
          svars(i)    = rLSNew(i)             ! S_ii  (no sum)
        enddo
        svars(4)    = rLSNew(5)               ! S_23
        svars(5)    = rLSNew(6)               ! S_13
        svars(6)    = rLSNew(4)               ! S_12

        svars(7)   = rLstNew(1)     ! S                plastic flow resistance, tensile
        svars(8)   = rLstNew(2)     ! TH               temperature
        svars(9)   = rLstNew(3)     ! F                plastic shear strain rate
        svars(10)  = rLstNew(4)     ! time int of F    plastic shear strain

        do i = 1, 6
          rSigInc(i) = rSigInc(i) + rStress(i)   ! negative of old stress + new stress
        enddo

C       Compute the dilatational wave speed

        call Gwavespd(rStrain,rSigInc,rE,rNu,rRho, rCd)

      endif    ! iConsType

C     3) Calculate d_Stress/d_Spin_increment (rGtan):

      call GrGtan(rQ,rdQdW,rRn,rRhalf,rRone,rScorot,rStress,rCmod,
     $            rEstat, rGtan)

C     4) Rotate the stress and moduli from co-rotational to stationary frame

      call rSrt(rRone,1.d0, rStress)
      call rrrrC(rRhalf,rRone, rCmod)

C     5) Recover the strain increment (which we rotated in step 1 above)

      do i = 1, 6
        rStrain(i) = rEstat(i)
      enddo

      end

