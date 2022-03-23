      SUBROUTINE VUEL(
c                     nblock = 136 is a parameter decided by ABAQUS
     1                nblock,
c                     below variables need to be defined
     2                rhs,amass,dtimeStable,
     3                svars,nsvars,
     4                energy,
c                     above variables need to be defined 
     5                nnode,ndofel,props,nprops,jprops,njprops,
     6                coords,mcrd,u,du,v,a,
     7                jtype,jElem,
     8                time,period,dtimeCur,dtimePrev,kstep,kinc,
     9                lflags,
     9                dMassScaleFactor,
     9                predef,npredef,
     9                jdltyp, adlmag)

C     Author: 
C          Tao Jin, T-3

C     Date: 
C          Created:        Sept 26th, 2017


C     Update:
C          1. fix the bug for writing state variables rAlphaN and rMvect;
C             (Oct 30th, 2017) 

C          2. Use the new modified quadrature rule proposed by Tao;
C             (Nov 15th, 2017)


C     Features:
C          ABAQUS/EXPLICIT
C          Finite deformation
C          3D 8-node brick element
C          Continuous shear band (CSB) propogation
C          Global tracking (solve the heat-like equation globally)
C          Displacement based or B-bar incompressible
C          Isotropic elasticity or linear J2 plasticity or MTS J2 plasticity
C          Localization flag is decided at the elmt centroid
C          New modified quadrature rule



C     Complie:
C          Define all the real variables as real*8 (double precision),
C          use "-double both" when compile

C     Voigt notation:
C          Different from FEAP or ABAQUS convection,

C          (1,1) -> 1
C          (2,2) -> 2
C          (3,3) -> 3
C          (2,3) or (3,2) -> 4
C          (1,3) or (3,1) -> 5
C          (1,2) or (2,1) -> 6


C    Contents of the PROPS array
C
C    index, i   |    PROPS(i)
C-------------------------------------------------------------------------------------
C    A)  Elastic & general properties (density, band width, numerical parameters, etc.)
C
C      1        |    Rho      Initial mass density (rho_0)
C      2        |    E        Young's modulus
C      3        |    Nu       Poisson's ratio
C      4        |    b        Localization band width
C      5        |    **UNUSED**
C      6        |    **UNUSED**
C      7        |    **UNUSED**
C      8        |    b1       coefficient of linear bulk viscosity term (default=0.06)
C      9        |    b2       coefficient of quadratic bulk viscosity term (default=1.2)
C     10        |    Epsilon  Pressure mode stabilization parameter; see Ortiz et al (1987), Eq. (21).
C
C
C    B)  Thermo-Elasto-Viscoplastic (Brown-Kim-Anand) Material Model
C
C     11        |    PLSLMT   limit on equiv. plastic tensil strain increment
C     12        |    PHI      degree of implicitness [0,1]
C     13        |    OMEGA    fraction of plastic work to heating
C     14        |    S0       initial value for internal variable S
C     15        |    T0       initial temperature (also used by MTS model)
C     16        |    ALPHA    thermal expansion coefficient
C     17        |    C        specific heat
C     18        |    A        pre-exponential factor
C     19        |    Q/R      activation energy divided by gas constant
C     20        |    XP       exponent in equation for F (1/n)
C     21        |    PALPHA   pressure sensitivity parameter
C     22        |    H0       pre-multiplier in hardening rate equation
C     23        |    SB       pre-multiplier for saturation value of S
C     24        |    AN       exponent in equation for saturation value of S
C     25        |    B        pre-multiplier for rate of static restoration
C     26        |    SA       annealed value of S
C     27        |    P        exponent in equation for RDOT
C     28        |    **UNUSED**
C     29        |    **UNUSED**
C     30        |    **UNUSED**
C     31        |    **UNUSED**
C     32        |    **UNUSED**
C     33        |    **UNUSED**
C
C    C) Isotropic Mises J2 Plasticity model
C
C     34        |    Y0       Initial yield stress
C     35        |    H0       Initial linear isotropic hardening modulus
C     36        |    Y1       Second yield stress
C     37        |    H1       Second linear isotropic hardening modulus
C     38        |    Y2       Yield stress at high strain values
C     39        |    **UNUSED**
C     40        |    **UNUSED**
C
C    D) Isotropic MTS model
C     41        |    Mu_0     Shear modulus constants
C     42        |    D_0      .
C     43        |    T_0      .
C     44        |    k        Boltzmann constant
C     45        |    b        Magnitude of Burger's vector
C     46        |    Edot_0i  Constants for the intrinsic-barriers term
C     47        |    q_i      .
C     48        |    p_i      .
C     49        |    g_0i     .
C     50        |    Edot_0e  Constants for the structure-evolution term
C     51        |    q_e      .
C     52        |    p_e      .
C     53        |    g_0e     .
C     54        |    sigma_a  Athermal MTS term
C     55        |    sigma_i  Intrinsic resistance constant
C     56        |    sigma_e  Initial structural resistance
C     57        |    Edot_0es Saturation stress constants
C     58        |    g_0es    .
C     59        |    sig_es0  .
C     60--65    |    A0--A5   Hardening premultiplier (Theta_0) constants
C     66        |    Kappa    Hardening exponent
C     67--69    |    Cp0--Cp2 Heat capacity constants
C     70        |    Psi      Fraction of plastic work converted to heat
C     71        |    T_m      Absolute melting temperature
C     72        |    T_drx    DRX homologous temperature (no DRX if T_drx <= 0)
C     73        |    D_0      DRX diffusion constant
C     74        |    Q/R      DRX activation energy divided by gas constant
C     75        |    d        DRX grain size
C     76        | iso_theta   the theta value for the isosurface


        IMPLICIT NONE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       PREDEFINED PARAMETERS & FLAGS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C       operational code keys
        integer jMassCalc,jIntForceAndDtStable,jExternForce
        parameter ( jMassCalc            = 1,
     *              jIntForceAndDtStable = 2,
     *              jExternForce         = 3)

C       flag indices
        integer iProcedure,iNlgeom,iOpCode,nFlags
        parameter (iProcedure = 1,
     *             iNlgeom    = 2,
     *             iOpCode    = 3,
     *             nFlags     = 3)

C       energy array indices
        integer iElPd,iElCd,iElIe,iElTs,iElDd,iElBv,
     *          iElDe,iElHe,iElKe,iElTh,iElDmd,iElDc,
     *          nElEnergy
        parameter ( iElPd = 1,
     *              iElCd = 2,
     *              iElIe = 3,
     *              iElTs = 4,
     *              iElDd = 5,
     *              iElBv = 6,
     *              iElDe = 7,
     *              iElHe = 8,
     *              iElKe = 9,
     *              iElTh = 10,
     *              iElDmd = 11,
     *              iElDc = 12,
     *              nElEnergy = 12)

C       predefined variables indices
        integer iPredValueNew,iPredValueOld,nPred
        parameter ( iPredValueNew = 1,
     *              iPredValueOld = 2,
     *              nPred         = 2)

C       time indices
        integer iStepTime,iTotalTime,nTime
        parameter (iStepTime  = 1,
     *             iTotalTime = 2,
     *             nTime      = 2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       PARAMETERS PASSED FROM THE INTERFACE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C       real scalar variables passed from VUEL interface
        real*8 period,dtimeCur,dtimePrev

C       integer scalar variables passed from VUEL interface
        integer nblock,nsvars,nnode,ndofel,nprops,njprops,mcrd,
     1          jtype,kstep,kinc,npredef,jdltyp

C       real Array variables passed from VUEL interface
        real*8 rhs(nblock,ndofel),amass(nblock,ndofel,ndofel),
     1         dtimeStable(nblock),
     2         svars(nblock,nsvars),energy(nblock,nElEnergy),
     3         props(nprops),time(nTime),coords(nblock,nnode,mcrd),
     4         u(nblock,ndofel), du(nblock,ndofel),
     5         v(nblock,ndofel), a(nblock, ndofel),
     6         dMassScaleFactor(nblock),  
     7         predef(nblock,nnode,npredef,nPred),adlmag(nblock)

C       integer Array variables passed from VUEL interface
        integer jprops(njprops),jElem(nblock),lflags(nFlags)

C       Twin Code Implementation Started
        integer kblock, i, j, l, n
        integer iTIntPt          ! Total # of Gauss PTS (=8) for MATRIX
        integer elmt_number      ! elmt global number
        integer iTwinFlg_Ct      ! Twinned flag at elmt centroid
        integer iumatsvst        ! umat start index for svars
        integer iumatsven        ! umat end index of svars

        real*8 rXiEtaPt(3,8)     ! Natural coords of 8 Gauss PTS
        real*8 rWeightPt(8)      ! Weights of 8 Gauss PTS
        real*8 rCoords(mcrd,nnode) ! Nodal coordinates (initial or current)
        real*8 rMSCoor(mcrd,nnode) ! Mid-step nodal coords (initial or current)
        real*8 rPrevCoords(mcrd, nnode) ! Previous coordinates
        real*8 rOrigCoords(mcrd, nnode) ! Original coordinates
        real*8 rNn(nnode),rdNndXi(mcrd,nnode) ! Eight Shape funcs and derivatives at one Gauss PT
        real*8 rJdet                          ! Jacobian at one Gauss PT
        real*8 rdOmega             ! Jacobian * Gauss_Weight
        real*8 rRho              ! Initial mass density (rho_0)
        real*8 rElmtVol            ! Elmt volume
        real*8 rElmtSize           ! Elmt size = (Elmt volume)**1/3
        real*8 rXiEtaCt(3)       ! element centroid natural coords
        real*8 rNnCt(nnode),rdNndXiCt(mcrd,nnode) ! Eight Shape funcs and derivatives at Elmt centroid
        real*8 rJdetCt             ! Jacobian at Elmt centroid
        real*8 rDfgrd0Ct(3, 3), rDfgrd1Ct(3, 3) ! Previous and Current Def grd at Elmt centroid
        real*8 rCtCoords(3)      ! Current Coordinates of element centroid
        real*8 rStressCt(6)      ! Cauchy Stress at Elmt Centroid

        !call sleep(20)
        rRho = props(78)
        call GSPT_3D_8N(iTIntPt, rXiEtaPt, rWeightPt)

C       Natural coordinates of element centroid
        rXiEtaCt(1) = 0.d0
        rXiEtaCt(2) = 0.d0
        rXiEtaCt(3) = 0.d0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                      MASS MATRIX CONSTRUCTION                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        IF ( lflags(iOpCode) .EQ. jMassCalc ) THEN
          do kblock = 1, nblock
C           Material coordinates (integrate mass terms over INITIAL config.)
            DO n = 1, nnode
              DO i = 1, mcrd
                rCoords(i,n) = coords(kblock,n,i)
              END DO
            END DO

C           Initialize the lumped mass matrix amass
            DO i = 1, ndofel
              DO j = 1, ndofel
                amass( kblock, i, j ) = 0.d0
              END DO
            END DO

            do l = 1, iTIntPt
C             Compute shape functions and Jacobian determinant
              CALL SHP_FUNC_3D_8N(rXiEtaPt(1,l),rCoords,0,rNn,rdNndXi,
     &                            rJdet)

C             Compute integration element (dOmega)
              rdOmega = rJdet * rWeightPt(l)

C             Assemble lumped mass matrix
              DO n = 1, nnode
                DO i = 1, mcrd
                  amass( kblock , (n-1)*mcrd+i , (n-1)*mcrd+i )
     $               = amass( kblock , (n-1)*mcrd+i , (n-1)*mcrd+i )
     $                 + rRho * rNn(n) * rdOmega
                END DO
              END DO
            end do
          end do

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCC            Define rhs, dtimeStable, svars, energy            CCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       Define internal force (rhs), stable time step
C       (dtimeStable), history variables (svars), and energy
        ELSE IF ( lflags(iOpCode) .EQ. jIntForceAndDtStable ) THEN
          do kblock = 1,nblock

            elmt_number = jElem(kblock)

C>>>        Get nodal coordinates, finite deformation
            DO n = 1, nnode
              DO i = 1, mcrd
                IF ( (kstep.GT.1).OR.(kinc.GT.1) ) THEN
                  rOrigCoords(i,n) = coords(kblock,n,i)
                  rCoords(i,n) = rOrigCoords(i,n)               ! current coords
     $                         + u( kblock, (n-1)*mcrd + i )
                  rMSCoor(i,n) = rCoords(i,n)                   ! mid-step coords
     $                         - 0.5d0 * du(kblock, (n-1)*mcrd + i )
                  rPrevCoords(i,n) = rCoords(i,n)
     &                            - du(kblock, (n-1)*mcrd + i )
                ELSE
                  rOrigCoords(i,n) = coords(kblock,n,i)
                  rCoords(i,n) = coords(kblock,n,i)       ! original coords
                  rMSCoor(i,n) = coords(kblock,n,i)       ! original coords
                  rPrevCoords(i,n) = coords(kblock, n, i)
                END IF
              END DO
            END DO

C           Compute element volume
            rElmtVol = 0.d0
            DO l = 1, iTIntPt
C             Jacobian determinant
              CALL SHP_FUNC_3D_8N(rXiEtaPt(1,l),rCoords,0,
     &                            rNn,rdNndXi,rJdet)
C             Integration element (dOmega)
              rdOmega = rJdet * rWeightPt(l)
C             Element surface area (even for axisymmetric elements)
              rElmtVol = rElmtVol + rdOmega
            END DO

C           Compute element size h
            rElmtSize = rElmtVol**(1.d0/3.d0)      !   h = cubic root of vol

C           If first increment of first step, initialize elmt centroid
C           twinning flag in history
            IF ( ( kstep .EQ. 1) .AND. ( kinc .EQ. 1) ) THEN
              svars(kblock, 1) = 0.d0
            END IF
            iTwinFlg_Ct = int(svars(kblock, 1))

C           Element with no twins
            if (iTwinFlg_Ct.eq.0) then

              call SHP_FUNC_3D_8N(rXiEtaCt,rCoords,1,rNnCt,
     &                            rdNndXiCt,rJdetCt)

              iumatsvst = 540 + iTIntPt * 142
              iumatsven = 539 + (iTIntPt + 1) * 142
              call calc_dfgrd(rXiEtaCt, rOrigCoords, rPrevCoords,
     &                         rCoords, mcrd, nnode, rDfgrd0Ct,
     &                         rDfgrd1Ct, rCtCoords)

              call umat_stress(rStressCt,
     &                         svars(kblock, iumatsvst:iumatsven),
     &                         rCtCoords, rDfgrd0Ct, rDfgrd1Ct,
     &                         props(77:126), 50, 142, nnode,
     &                         elmt_number, 0, mcrd, time, dtimeCur,
     &                         kstep, kinc)

            end if
          end do
        ENDIF
      END  ! Subroutine VUEL


      include 'GSPT_3D_8N.f'
      include 'SHP_FUNC_3D_8N.f'
      include 'prep_umat.f'
      include 'gdefgrad_3D.f'
      include 'fcc_umat_lbl34_test.f'