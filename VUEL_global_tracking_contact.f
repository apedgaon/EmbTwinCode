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
#include <SMAAspUserSubroutines.hdr> ! Start this line WITHOUT ANY SPACE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PREDEFINED PARAMETERS & FLAGS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     operational code keys
      integer jMassCalc,jIntForceAndDtStable,jExternForce
      parameter ( jMassCalc            = 1,
     *            jIntForceAndDtStable = 2,
     *            jExternForce         = 3)

C     flag indices
      integer iProcedure,iNlgeom,iOpCode,nFlags
      parameter (iProcedure = 1,
     *           iNlgeom    = 2,
     *           iOpCode    = 3,
     *           nFlags     = 3)

C     energy array indices
      integer iElPd,iElCd,iElIe,iElTs,iElDd,iElBv,
     *        iElDe,iElHe,iElKe,iElTh,iElDmd,iElDc,
     *        nElEnergy
      parameter ( iElPd = 1,
     *            iElCd = 2,
     *            iElIe = 3,
     *            iElTs = 4,
     *            iElDd = 5,
     *            iElBv = 6,
     *            iElDe = 7,
     *            iElHe = 8,
     *            iElKe = 9,
     *            iElTh = 10,
     *            iElDmd = 11,
     *            iElDc = 12,
     *            nElEnergy = 12)

C     predefined variables indices
      integer iPredValueNew,iPredValueOld,nPred
      parameter ( iPredValueNew = 1,
     *            iPredValueOld = 2,
     *            nPred         = 2)

C     time indices
      integer iStepTime,iTotalTime,nTime
      parameter (iStepTime  = 1,
     *           iTotalTime = 2,
     *           nTime      = 2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PARAMETERS PASSED FROM THE INTERFACE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     real scalar variables passed from VUEL interface
      real*8 period,dtimeCur,dtimePrev

C     integer scalar variables passed from VUEL interface
      integer nblock,nsvars,nnode,ndofel,nprops,njprops,mcrd,
     1        jtype,kstep,kinc,npredef,jdltyp

C     real Array variables passed from VUEL interface
      real*8 rhs(nblock,ndofel),amass(nblock,ndofel,ndofel),
     1       dtimeStable(nblock),
     2       svars(nblock,nsvars),energy(nblock,nElEnergy),
     3       props(nprops),time(nTime),coords(nblock,nnode,mcrd),
     4       u(nblock,ndofel), du(nblock,ndofel),
     5       v(nblock,ndofel), a(nblock, ndofel),
     6       dMassScaleFactor(nblock),  
     7       predef(nblock,nnode,npredef,nPred),adlmag(nblock)

C     integer Array variables passed from VUEL interface
      integer jprops(njprops),jElem(nblock),lflags(nFlags)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PARAMETERS NEEDED IN THE BODY
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     procedure flags
      integer jDynExplicit,jVUnumber
      parameter ( jDynExplicit = 17, jVUnumber = 1 )

C     Default safety factor for stable time step estimate
      real*8 rStFctr
      parameter ( rStFctr = 0.9d0 )

      integer n,i,j,l,k
      integer kblock,iVelmType,iVforType,iConsType
      integer nstatev(6)       ! # of state vars at Gauss pt
                               ! (1) - # of state vars at Gauss pt for Matrix
                               ! (2) - # of state vars at Gauss pt for Shear Band
                               ! (3) - # of state vars related with localization 
                               ! (4) - # of state vars at elmt centroid   

      integer MSflag           ! 0-> no material instability;
                               ! 1-> rate independent case, bifurcation analysis
                               ! 2-> rate dependent case, perturbation analysis


      integer iOutPr(15)
      integer iStLocF          ! 1 -> Element has localization at any Gauss point
                               !      set flag to make time step smaller 

      integer iCvFlag          ! Convergence flag: 0 -> non-convergence
                               !                   1 -> convergence
                               !                   2 -> bifurcation is detected 
                               !                        (For rate-independent material)
                               !                   3 -> material is unstable (For MTS model)
                               !                  -1 -> elastic behavior is assumed, 
                               !                        but the yield condition is violated
      integer elmt_number      ! elmt global number
      integer neigh_elmt_no,neigh_theta_counter,neigh_theta_sign
      integer iAltern,iHaveSol,iDbgFlg,iAltMax         

      real*8 rRho,rBWidth      ! Initial mass density (rho_0) & Localization band width
      real*8 rbv1,rbv2         ! linear & quadratic bulk viscosity coefficients


      real*8 pi                ! pi = 4.d0*atan(1.d0)


      real*8 rCoords(mcrd,nnode) ! Nodal coordinates (initial or current)
      real*8 rPrevCoords(mcrd, nnode) ! Previous coordinates
      real*8 rOrigCoords(mcrd, nnode) ! Original coordinates

      real*8 rdOmega             ! Jacobian * Gauss_Weight
      real*8 rMSCoor(mcrd,nnode) ! Mid-step nodal coords (initial or current)
      real*8 rElmtVol            ! Elmt volume
      real*8 rElmtSize           ! Elmt size = (Elmt volume)**1/3
      real*8 rMinElmtL           ! Size used in estimating stable time step
      real*8 rYoungsMod          ! Young's modulus from input
      real*8 rMaxWavSpd          ! Maximum wave speed

      real*8 rue(mcrd,nnode)     ! local array for passing displacements increment to strain routine

      real*8 hvars(nsvars)       ! local array for passing state vars to stress update routine

      real*8 rCmod(6,6)          ! Tangent moduli (d_Stress/d_Strain_increment)
      real*8 rGtan(6,3)          ! d_Stress/d_Spin_increment (finite def. case only)
      real*8 rCd                 ! Dilatational wave speed
      real*8 rtemp(2)
      real*8 rBVpress            ! bulk viscosity
      real*8 rSigma              ! Used in R.H.S calculation
      real*8 rOutPr(15)

      real*8 elmt_theta(8)       ! elmt nodal theta values
      real*8 theta_0             ! reference temperature for the isosurface of shear band surface

      real*8 natural_coords(12,3) ! natural coords of the intersection between shear band
                                  ! and elmt edge

      real*8 dist(8)              ! distance between gauss pt and shear band center
      real*8 dist_ave             ! average distance
      real*8 Gt_weight(8)         ! Weight-coefficient for each Gauss point based on the distance
      real*8 Gt_temp              ! Weight-coefficient at one Gauss point
      real*8 rBifCoor(mcrd,nnode), rBifF(3,3) ! X* (nodal coords at time of bifurcation) and F*=dx/dX* 
      real*8 rAlphaN,rAlphaHat     ! rAlphaHat = rTL * rAlphaN
      real*8 rTtens(6), rTtensk(3) ! symmetric and skew-symmetric parts of T tensor 
                                   ! characterizing localization mode
      real*8 rTEps                 ! [T^sym] : [grad(v)^sym]
      real*8 rTSpn                 ! [T^skw] : [grad(v)^skw]
      real*8 rLnrm                 ! norm of grad(v)
      real*8 rTL                   ! rTL = T : grad(v) = rTEps + rTSpn
      real*8 rMhat(3)              ! rMhat = rAlphaHat * rMvect
      real*8 rLn(3),rLnNrm         ! grad(v)*n and its norm
      real*8 rAlphaL               ! rAlphaL   = ( elmt_vol/SB_vol/Gt_weight - 1.d0 ) * rAlphaN
      real*8 rAlphaBar             ! rAlphaBar = ( elmt_vol/SB_vol/Gt_weight - 1.d0 ) * rAlphaHat
      real*8 rBmtL(6,mcrd*nnode)   ! B matrix in the localized band
      real*8 rBmtN(6,mcrd*nnode)   ! B matrix in the non-localized zone
      real*8 rSpinL(3)             ! vector of T^skw in shear band
      real*8 rSpinN(3)             ! vector of T^skw in non-localized zone
      real*8 rStrnL(6)             ! T^sym (Voigt notation) in the localized band
      real*8 rStrnN(6)             ! T^sym (Voigt notation) in the non-localized band
      real*8 rnCn(3,3)             ! matrix for the linearized traction continuity eq.
      real*8 rG(3),rGnrm           ! residual of the nonlinear traction continuity eq.
                                   ! and its norm
      real*8 rGnrmZ                ! initial residual norm
      real*8 rGrel                 ! relative residual norm

      real*8 rNvect(3),rMvect(3) ! n,m vectors

C     variables for elmt centroid
      integer iLocalFlg_Ct     ! Localization flag at elmt centroid
      real*8 rXiEtaCt(3)       ! element centroid natural coords
      real*8 rNnCt(nnode),rdNndXiCt(mcrd,nnode) ! Eight Shape funcs and derivatives at Elmt centroid
      real*8 rJdetCt             ! Jacobian at Elmt centroid
      real*8 rMSdNndXiCt(mcrd,nnode)   ! Shape func derivatives at Elmt centroid at Mid-step
      real*8 rMSJdetCt                 ! jacobian at Elmt centroid at Mid-step
      real*8 rEpsilon          ! Pressure mode stabilization parameter; see Ortiz et al (1987), Eq. (21).
      real*8 rStrainCt(6)          ! strain increment
      real*8 rSpinCt(3)            ! spin axial vector increment
      real*8 rMSBmCt(6,mcrd*nnode) ! B (or B-bar) operator at mid-step
      real*8 rMSBsCt(3,mcrd*nnode) ! skew-symmetric matrix operator at mid-step
      real*8 rStressCt(6),rStrIncCt(6) ! stress and stress increment

C     variables for matrix material
      integer iLocalFlg        ! Localization flag at Gauss Point
      integer iIntPt           ! # of Gauss PTS per direction for MATRIX
      integer iTIntPt          ! Total # of Gauss PTS (=8) for MATRIX
      real*8 rXiEtaPt(3,8)     ! Natural coords of 8 Gauss PTS
      real*8 rWeightPt(8)      ! Weights of 8 Gauss PTS
      real*8 rNn(nnode),rdNndXi(mcrd,nnode) ! Eight Shape funcs and derivatives at one Gauss PT
      real*8 rJdet                          ! Jacobian at one Gauss PT
      real*8 rMSdNndXi(mcrd,nnode)     ! Shape func derivatives at one Gauss PT at Mid-step 
      real*8 rMSJdet                   ! Jacobian at one Gauss PT at Mid-step
      real*8 rStrain(6)          ! strain increment
      real*8 rSpin(3)            ! spin axial vector increment
      real*8 rBmat(6,mcrd*nnode) ! current B (or B-bar) operator
      real*8 rBskw(3,mcrd*nnode) ! current skew-symmetric matrix operator
      real*8 rMSBm(6,mcrd*nnode) ! B (or B-bar) operator at mid-step
      real*8 rMSBs(3,mcrd*nnode) ! skew-symmetric matrix operator at mid-step
      real*8 rEpVolInc           ! volumetric strain increment (for bulk viscosity calculation)
      real*8 rStress(6),rStrInc(6) ! stress and stress increment 


C     variables for shear band center
      real*8 SB_ct(3)          ! natural coords of the shear band center
      real*8 SB_vol            ! shear band volume
      real*8 normalCt(3)       ! shear band surface normal in the Cartesian coordinates
      real*8 natural_normalCt(3)       ! shear band surface normal in the natural coordinates

C     Parameters for iteratively finding rAlphaN and rMvect
      integer iConverge,iIterMax,iCount,iElFlag,iErrFlg
      real*8 rDelNrm,rIterTol,rRelTol,rTfactor,rDRXfac,rGmin(5)
      real*8 rDelta(3),rCfacR,rCfacC,rCFactor

C     for output SDV
      integer iTemp
      real*8 uel_isvs(8,10)

C     Define 6 faces:
C     Face1: 1,2,3,4; Face2: 5,8,7,6; Face3: 1,5,6,2;
C     Face4: 2,6,7,3; Face5: 3,7,8,4; Face6: 4,5,8,1;
      integer FACE(6,4)
      data FACE /1,5,1,2,3,4,
     1           2,8,5,6,7,8,
     2           3,7,6,7,8,5,
     3           4,6,2,3,4,1/  

      integer elmt_conct(8,*)
      pointer(ptr_elmt_conct, elmt_conct)

      integer elmt_neighbour(6,*)
      pointer(ptr_elmt_neighbour, elmt_neighbour)

      integer elmt_status(2,*)
      pointer(ptr_elmt_status, elmt_status)

      real*8 node_theta(*)
      pointer(ptr_node_theta, node_theta)

      real*8 elmt_Nvect(3,*)
      pointer(ptr_elmt_Nvect, elmt_Nvect)

      real*8 node_coords(3,*)
      pointer(ptr_node_coords,node_coords)

      real*8 elmt_SDV(10,*)
      pointer(ptr_elmt_SDV,elmt_SDV)

      integer iumatsvst, iumatsven      ! umat_crys state variable start and stop index
      real*8 rDfgrd0(3, 3), rDfgrd1(3, 3)
      real*8 rGsCoords(3, 8)
      real*8 rDfgrd0Ct(3, 3), rDfgrd1Ct(3, 3), rCtCoords(3)
      real*8 rImat(3, 3), rdVelGrd(3, 3)
      real*8 Smat(3, 3, 12), F_tw(3, 3, 12) ! Schmid Tensor and twin dfgrd

C     READ the theta value for the iso-surface from input
      theta_0 = props(76)

      IF (kinc.GE.1) THEN
        ptr_elmt_conct = SMAIntArrayAccess(1)
        ptr_elmt_neighbour = SMAIntArrayAccess(2)
        ptr_elmt_status = SMAIntArrayAccess(3)
        ptr_node_theta = SMAFloatArrayAccess(4)
        ptr_elmt_Nvect = SMAFloatArrayAccess(8)
        ptr_node_coords = SMAFloatArrayAccess(7)
        ptr_elmt_SDV = SMAFloatArrayAccess(11)
      END IF

      IF ( (jtype.EQ.jVUnumber) .AND. 
     $     (lflags(iProcedure).EQ.jDynExplicit) ) THEN

        iVelmType = jprops(1)  ! 3D element type:  1 --> 3D Brick (8-node) Elmt


        iVforType = jprops(2)  ! Formulation type: 0 --> displacement based
                               !                   1 --> displacement based with localization behavior
                               !                   2 --> B-bar incompressible
                               !                   3 --> B-bar incompressible with localization behavior

        iConsType = jprops(3)  ! Constitution:     0 --> isotropic elastic
                               !                   1 --> isotropic J2 plasticity or MTS
                               !                   2 --> TEPLA
                               !                         (NOT IMPLEMENTED)
                               !                   3 --> Brown-Kim-Anand model
                               !                         (NOT IMPLEMENTED)
                               !                   4 --> models (1) outside and (3) inside localization band 
                               !                         (NOT IMPLEMENTED)
                               !                   5 --> models (2) outside and (3) inside localization band 
                               !                         (NOT IMPLEMENTED)

        IF (      (iConsType.EQ.2) .OR. (iConsType.EQ.3)
     $       .OR. (iConsType.EQ.4) .OR. (iConsType.EQ.5) ) THEN
          CALL XPLB_ABQERR(-3,'Constitutive model '//
     $                   'has not been implemented yet',0,0.d0,' ')
        END IF

C       Notes pertaining to STATE VARIABLES: 
C       -----------------------------------
C              (1) In formulations WITHOUT localization, the number of state variables
C                  PER GAUSS POINT is equal to nstatev(1).

C              (2) In formulations WITH localization, the number of state variables
C                  PER GAUSS POINT is equal to nstatev(1)+nstatev(2)+nstatev(3). That is,
C                  the material inside the localization band uses nstatev(2) additional
C                  state variables, and other localization-specific variables need nstatev(3)
C                  memory slots. Those are only used after the onset of localization,
C                  but the memory is reserved (via the "*USER ELEMENT, VARIABLES=" option
C                  in the input file) throughout the analysis.

C              (3) We use EIGHT (8) GAUSS points for the numerical integration.

C              (4) The onset condition of localization is checked at the element centroid. 
C                  Once the onset of localization is detected at the element centroid, all the
C                  EIGHT (8) Gauss points for the full integration have localization. At the elmt
C                  centroid, besides nstatev(4) variables, use:
C                1 for the localization flag
C             2--4 for the coordinates of the shear band center
C                5 for the shear band volume
C            6--13 for the weight coefficient based on the distance between each of 
C                  EIGHT (8) Gauss points to the shear band surface   


C              (5) The number of state variables PER ELEMENT is equal to:
C                  the number of Gauss points (8) x the number of state variables PER GAUSS POINT.
C                + the number of state variales for the element centroid
C                  which is 
C                  8 * [ nstatev(1) + nstatev(2) + nstatev(3) ]
C                  1 * [ nstatev(4) + 1 ]

C              (6) In formulations with localization, nstatev(3) = 36 state variables:
C                1 The localization flag at the current Gauss point, which is initially 0.0,
C                  and is set to 1.0 after the onset of localization.
C             2--4 Three additional variables are used to save the components of the band's 
C                  (ORIGINAL) orientation vector (n0).
C                5 Localization parameter (alphaN)
C             6--8 Three additional variables are used to save the components of the orientation
C                  vector (m), defining the strain discontinuity.
C            9--32 In the finite deformation case, twenty-four (3x8) additional variables
C                  are used to store the nodal coordinates defining the configuration (X*),
C                  where localization was initiated.
C           33--35 The natural coordinates of the shear band centroid 
C                * In formulations w/o localization, these memory slots need not be
C                  reserved via the "*USER ELEMENT, VARIABLES=" option.
C               36 The volume of the shear band 


C              (3) The number of state variables PER ELEMENT is equal to: 
C                  the number of Gauss points for MATRIX times nstatev(1) +
C                  the number of Gauss points for SHEAR BAND times nstatev(2) +
C                  nstatev(3) +
C                  nstatev(4)


C       Set number of state variables PER GAUSS POINT for matrix and band materials
C       For each model in rate/incremental form, we must keep track of 6 components
C       of stress, plus a number of model-specific history variables.
        select case (iConsType)
          case default          ! isotropic elastic
            nstatev(1) = 6
            nstatev(2) = 6
          case (1)              ! isotropic J2 plasticity / MTS
            nstatev(1) = 6+5    ! matrix
            nstatev(2) = 6+5    ! shear band
            nstatev(4) = 6+5    ! elmt centroid
          case (2)              ! TEPLA (NOT IMPLEMENTED)
            nstatev(1) = 6+5
            nstatev(2) = 6+5
          case (3)              ! Brown-Kim-Anand model (NOT IMPLEMENTED)
            nstatev(1) = 6+4
            nstatev(2) = 6+4
          case (4)              ! models (1) outside and (3) inside localization band (NOT IMPLEMENTED)
            nstatev(1) = 6+5
            nstatev(2) = 6+4
          case (5)              ! models (2) outside and (3) inside localization band (NOT IMPLEMENTED)
            nstatev(1) = 6+5
            nstatev(2) = 6+4
        end select

C       Select the kind of material stability analysis
        select case (iConsType)
          case default          ! isotropic elastic
            MSflag = 0    ! no material instability
          case (1)              ! isotropic J2 plasticity / MTS
            if (props(41).le.0.d0) then ! rate independent case
              MSflag = 1    ! bifurcation analysis
            else
              MSflag = 2    ! perturbation analysis
            endif
          case (2)              ! TEPLA
            MSflag = 2    ! perturbation analysis
          case (3)              ! Brown-Kim-Anand model
            MSflag = 2    ! perturbation analysis
          case (4)              ! models (1) outside and (3) inside localization band
            if (props(41).le.0.d0) then ! rate independent case
              MSflag = 1    ! bifurcation analysis
            else
              MSflag = 2    ! perturbation analysis
            endif
          case (5)              ! models (2) outside and (3) inside localization band
            MSflag = 2    ! perturbation analysis
        end select

C       In the finite deformation case, we need to keep track of material rotation
        nstatev(1) = nstatev(1) + 3     ! 3D rotation vector (i.e. axis-angle representation)
        nstatev(2) = nstatev(2) + 3
        nstatev(4) = nstatev(4) + 3

C       If formulation w/ localization, we need extra memory space
        IF ( (iVforType.EQ.1) .OR. (iVforType.EQ.3) ) THEN
          nstatev(3) = 1                  !  localization flag
          nstatev(3) = nstatev(3) + 3     !  n0 vector
          nstatev(3) = nstatev(3) + 1     !  alpha_N
          nstatev(3) = nstatev(3) + 3     !  m vector
          nstatev(3) = nstatev(3) + mcrd*nnode   !  nodal coordinates for X* configuration
          nstatev(3) = nstatev(3) + 3     !  natural coordinates of the shear band centroid
          nstatev(3) = nstatev(3) + 1     !  volume of the shear band
        ELSE
          nstatev(2) = 0
          nstatev(3) = 0
        END IF
        
        nstatev(5) = 106
        nstatev(6) = 9 
        iIntPt = jprops(4)  ! Number of integration points per direction for MATRIX

        IF (iIntPt.NE.2) THEN
          CALL XPLB_ABQERR(-3,'Only 2x2x2 Gauss PTS can be used for ' //
     &                        'MATRIX', 0,0.d0,' ')
        END IF

        rRho      = props(1)   ! Initial mass density (rho_0)

        rBWidth   = props(4)   ! Localization band width

        IF (props(8).LE.0.d0) THEN     ! linear bulk viscosity coefficient
          rbv1    = 0.06d0
        ELSE
          rbv1    = props(8)
        END IF

        IF (props(9).LE.0.d0) THEN     ! quadratic bulk viscosity coefficient
          rbv2    = 1.2d0
        ELSE
          rbv2    = props(9)
        END IF

        rEpsilon  = props(10)  ! Pressure mode stabilization parameter; see Ortiz et al (1987), Eq. (21).

C       Natural coordinates of element centroid
        rXiEtaCt(1) = 0.d0
        rXiEtaCt(2) = 0.d0
        rXiEtaCt(3) = 0.d0

C       Constants, etc.
        pi = 4.d0*atan(1.d0)

C       Get natural coordinates and weights for 8 Gauss points for MATRIX and SHEAR BAND
        CALL GSPT_3D_8N(iTIntPt, rXiEtaPt, rWeightPt)

C       Initialize Identity Matrix(3x3)
        do i = 1, 3
            do j = 1, 3
                if (i.eq.j) then
                    rImat(i, j) = 1.0d0;
                else
                    rImat(i, j) = 0.0d0;
                end if
            end do
        end do

C       Initialize Deformation Gradients
        rDfgrd0 = rImat
        rDfgrd1 = rImat

C       Compute Schmid Tensor and Deformation gradient for Twinning        
        call comp_twin_dfgrd(rImat, Smat, F_tw)

C       The amount of memory reserved for element state variables must be sufficient
C       The elmt centroid has one extra slot to store the localization flag
        iOutPr(1) = iTIntPt * ( nstatev(1) + nstatev(2) + nstatev(3) )
     1            +     1   * ( nstatev(4) + 13 )

        iOutPr(2) = nsvars
        IF ( nsvars .LT. iOutPr(1) ) THEN
          CALL XPLB_ABQERR(-3,'The minimum number of memory slots '//
     $                        'required for state variable storage '//
     $                        'is %I. Only %I slots have been '//
     $                        'reserved.',iOutPr,0.d0,' ')
        END IF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCC                 Define mass matrix                 CCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        IF ( lflags(iOpCode) .EQ. jMassCalc ) THEN  ! Define mass matrix

          DO kblock = 1,nblock

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

            DO l = 1, iTIntPt ! Loop over Gauss PTS, Numerical integration

C             Compute shape functions and Jacobian determinant
              CALL SHP_FUNC_3D_8N(rXiEtaPt(1,l),rCoords,0,
     1                            rNn,rdNndXi,rJdet)

C             Compute integration element (dOmega)
              rdOmega = rJdet * rWeightPt(l)

C             Assemble lumped mass matrix
              DO n = 1, nnode
                DO i = 1, mcrd
                  amass( kblock , (n-1)*mcrd+i , (n-1)*mcrd+i )
     $          = amass( kblock , (n-1)*mcrd+i , (n-1)*mcrd+i )
     $          + rRho * rNn(n) * rdOmega
                END DO
              END DO

            END DO ! l: Loop over Gauss PTS, Numerical integration

          END DO  ! kblock (to define mass matrix)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCC            Define rhs, dtimeStable, svars, energy            CCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C       Define internal force (rhs), stable time step
C       (dtimeStable), history variables (svars), and energy
        ELSE IF ( lflags(iOpCode) .EQ. jIntForceAndDtStable ) THEN

          DO kblock = 1,nblock

            elmt_number = jElem(kblock)

C           Initialize the RHS matrix
            DO n = 1, nnode
              DO i = 1, mcrd
                rhs( kblock , (n-1)*mcrd + i ) = 0.d0
              END DO
            END DO

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

            if (kinc.ge.1) then
              do n = 1, nnode
                node_coords(:, elmt_conct(n,elmt_number) )=rCoords(:,n)
              end do
            end if

C           We must pass the displacement INCREMENT "du" into GSTRAIN
            DO n = 1, nnode
              DO i = 1, mcrd
                IF ( (kstep.GT.1).OR.(kinc.GT.1) ) THEN
                  rue(i,n) = du(kblock, (n-1)*mcrd+i )
                ELSE
                  rue(i,n) = 0.d0
                END IF
              END DO
            END DO

C           Compute element volume
            rElmtVol = 0.d0
            DO l = 1, iTIntPt
C             Jacobian determinant
              CALL SHP_FUNC_3D_8N(rXiEtaPt(1,l),rCoords,0,
     1                            rNn,rdNndXi,rJdet)
C             Integration element (dOmega)
              rdOmega = rJdet * rWeightPt(l)
C             Element surface area (even for axisymmetric elements)
              rElmtVol = rElmtVol + rdOmega
            END DO

C           Prepare Deformation Gradients and Coordinates for all Gauss Pts
C           Required for UMAT type function
C           call dfgrd_prep(iTIntPt, rOrigCoords, rCoords, rPrevCoords,
C                           rXiEtaPt, mcrd, nnode, rDfgrd0, rDfgrd1,
C                           rGsCoords)

C           Compute element size h
            rElmtSize = rElmtVol**(1.d0/3.d0)      !   h = cubic root of vol
            
C           Size used in estimating stable time step
            CALL Clength_3D_8N(rCoords,rMinElmtL)
            rMinElmtL = min(rMinElmtL,rElmtSize)

C>>>        Compute the dilatational wave speed based on the elastic properties of the
C           (virgin and/or matrix) material.
 
C Note:       This is a rough estimate. A more accurate estimate will be computed below,
C             taking into consideration the constitutive response.

            rYoungsMod  = props(2)
            rMaxWavSpd  = sqrt( rYoungsMod / rRho )

C>>>        If incompressible B-bar formulation, we need to evaluate the following,
C           AT THE ELEMENT CENTROID:

            IF ( (iVforType .EQ. 1) .OR. 
     1           (iVforType .EQ. 2) .OR. 
     2           (iVforType .EQ. 3) ) THEN
C             a) shape functions and their derivatives wrt global coords
              CALL SHP_FUNC_3D_8N(rXiEtaCt,rCoords,1,rNnCt,
     1                            rdNndXiCt,rJdetCt)
C             b) shape functions and their derivatives wrt mid-step global coords
              CALL SHP_FUNC_3D_8N(rXiEtaCt,rMSCoor,1,rNnCt,
     1                            rMSdNndXiCt,rMSJdetCt)
            END IF   ! Incompressible B-bar

C>>>        Initialize localization flag for the elmt centroid
C           If formulation with localization
            IF ((iVforType.EQ.1).OR.(iVforType.EQ.3)) THEN

C             If first increment of first step, initialize elmt centroid
C             localization flag in history
              IF ( ( kstep .EQ. 1) .AND. ( kinc .EQ. 1) ) THEN
                svars(kblock, iTIntPt * nstatev(1) 
     $                      + iTIntPt * nstatev(2)
     $                      + iTIntPt * nstatev(3)
     $                      +     1   * nstatev(4)
     $                      + 1) = 0.d0
              END IF

C             Read elmt centroid localization status from history
              iLocalFlg_Ct = nint( svars(kblock, iTIntPt * nstatev(1) 
     $                                         + iTIntPt * nstatev(2)
     $                                         + iTIntPt * nstatev(3)
     $                                         +    1    * nstatev(4)
     $                                         + 1) )

            ELSE    ! form w/o localization

              iLocalFlg_Ct = 0

            END IF   ! form w/ or w/o localization

C           If formulation with localization, use the information at
C           the elmt centroid to detect localization
            IF ( (iVforType.EQ.1).OR.(iVforType.EQ.3) ) THEN

              IF (iLocalFlg_Ct.EQ.0) THEN

C               Strian and spin at elmt centroid   
                CALL Gstrain_3D(rue,rEpsilon,rNnCt,rNnCt,rMSdNndXiCt,
     $                          rMSdNndXiCt,mcrd,nnode,iVelmType,0,1,
     $                          rStrainCt,rSpinCt,rMSBmCt,rMSBsCt)

C               Select the correct state/history variables for elmt centroid
                DO i = 1, nstatev(4)
                  hvars(i) = svars(kblock, iTIntPt * nstatev(1) 
     $                                   + iTIntPt * nstatev(2) 
     $                                   + iTIntPt * nstatev(3)
     $                                   + i  )
                END DO

                iumatsvst = 540 + iTIntPt * nstatev(5)
                iumatsven = 539 + (iTIntPt + 1) * nstatev(5)
                call calc_dfgrd(rXiEtaCt, rOrigCoords, rPrevCoords,
     &                         rCoords, mcrd, nnode, rDfgrd0Ct,
     &                         rDfgrd1Ct, rCtCoords)

                call umat_stress(rStressCt,
     &                           svars(kblock, iumatsvst:iumatsven),
     &                           rCtCoords, rDfgrd0Ct, rDfgrd1Ct,
     &                           props(77:118), 42, nstatev(5), nnode,
     &                           elmt_number, 0, mcrd, time, dtimeCur,
     &                           kstep, kinc)

C               Compute the stress at elmt centroid
C               (now rStrain and rSpin are at mid-step)
                CALL Gstress_3D(rStrainCt,rSpinCt,hvars,nstatev(4),
     $                          props,nprops,iVelmType,iVforType,
     $                          iConsType,time,dtimeCur,rStressCt,
     $                          rStrIncCt,rCmod,rGtan,rCd,iCvFlag)

C               Find shear surface normal direction for each elmt and store 
C               the informaiton in the global array ID = 8
                if (kinc.GE.1) then
                  call surface_normal(rStressCt,normalCt)
                  elmt_Nvect(:,elmt_number) = normalCt 
                end if         

C               Check whether accumulated effective plastic strain increases monotonically
                IF (svars(kblock,7 + iTIntPt*nstatev(1) 
     $                             + iTIntPt*nstatev(2) 
     $                             + iTIntPt*nstatev(3)).GT.hvars(7)) THEN
                  rtemp(1) = svars( kblock,7 + iTIntPt * nstatev(1) 
     $                                       + iTIntPt * nstatev(2) 
     $                                       + iTIntPt * nstatev(3) )
                  rtemp(2) = hvars(7)
                  CALL XPLB_ABQERR(-3,'Accumulated plastic strain '//
     $                             'fails to increase monotonically,'//
     $                 'prev_step = %R, current_step = %R',0, rtemp,' ')
                END IF

C               Write back the state/history variables for elmt centroid
                DO i = 1, nstatev(4)
                  svars(kblock, iTIntPt * nstatev(1) 
     $                        + iTIntPt * nstatev(2) 
     $                        + iTIntPt * nstatev(3)
     $                        + i) = hvars(i)
                END DO

                IF (MSflag.EQ.1) THEN          ! rate independent
                  iOutPr(1) = jelem(kblock)
                  iOutPr(2) = 9  ! elmt centroid
                  iOutPr(3) = kinc
                  rOutPr(1) = time(1)
                  rOutPr(2) = time(2)
C                 call bifurc_3D(rStress,rCmod,1,
C    $                           iOutPr,rOutPr, iCvFlag,rNvect,rMvect)
                  CALL XPLB_ABQERR(-3,'Bifurcation analysis for 3D '//
     $                   'case has not been implemented yet',0,0.d0,' ')
                ELSE IF (MSflag.EQ.2) THEN      ! rate dependent
                  IF (      (iCvFlag.EQ.3)
     $                 .AND.(elmt_status(2,elmt_number).EQ.1) ! elmt has to be in the front list
     $               ) THEN       ! instability detected
                    iOutPr(1) = jelem(kblock)
                    iOutPr(2) = 9  ! elmt centroid
                    iOutPr(3) = kinc
                    rOutPr(1) = time(1)
                    rOutPr(2) = time(2)
                    CALL perturb_3D(rStressCt,iOutPr,rOutPr,
     $                              rNvect,rMvect)
                  END IF
                END IF 

                IF (    (iCvFlag.GT.1)
     $             .AND.(elmt_status(2,elmt_number).EQ.1) ! elmt has to be in the front list
     $             ) THEN         ! Localization detected just now at the elmt centroid

C                 We need to:

C                 a) Update the localization flag for the elmt centroid
                  svars(kblock,  iTIntPt*nstatev(1)
     $                         + iTIntPt*nstatev(2)
     $                         + iTIntPt*nstatev(3)
     $                         +    1   *nstatev(4) + 1 ) = 1.d0       

C                 b) Find the location of the shear band inside element
                  DO i = 1, 8
                    elmt_theta(i)= node_theta(elmt_conct(i,elmt_number))
                  END DO

                  ! Update the elmt shear and front status
                  elmt_status(1,elmt_number) = 1  ! shear band status
                  elmt_status(2,elmt_number) = 0  ! front status

C                 Calculate the natural coords of intersection on each edge,
C                 the center point of the Shear band surface and
C                 the shear band volume
                  CALL SB_ct_vol(elmt_number,theta_0,natural_coords,
     1                          rBWidth, rcoords, SB_ct, SB_vol)
                  IF (SB_vol.GE.rElmtVol) Then
                    print *, 'Elmt', elmt_number
                    print *, 'SB_vol=', SB_vol
                    print *, 'El_vol=', rElmtVol
                    CALL XPLB_ABQERR(-3,'Shear band volume is larger '//
     $                              'than elmt volume',0,0.d0,' ')
                  END IF

C                 Loop over the neighbouring elmts to update the front list
C                 status of neighboring elmts              
                  DO i = 1,6
                    neigh_elmt_no = elmt_neighbour(i,elmt_number)
                    IF (neigh_elmt_no.NE.0) THEN
                      neigh_theta_sign = 0
                      DO j = 1,4
                        IF ( node_theta( elmt_conct( FACE(i,j) ,
     1                                   elmt_number ) )
     2                     .GT.theta_0) THEN
                          neigh_theta_sign = neigh_theta_sign + 1
                        END IF
                      END DO
                      IF (   (elmt_status(1,neigh_elmt_no).NE.1)
     1                 .AND. (neigh_theta_sign.NE.4)
     2                 .AND. (neigh_theta_sign.NE.0) ) THEN
                        elmt_status(2,neigh_elmt_no) = 1
                      END IF
                    END IF
                  END DO

C                 c) save the natural coords of shear band center
C                    and the volume of the shear band
                  svars(kblock,  iTIntPt * nstatev(1)
     $                         + iTIntPt * nstatev(2)
     $                         + iTIntPt * nstatev(3)
     $                         +    1    * nstatev(4)
     $                         +    2  ) = SB_ct(1)

                  svars(kblock,  iTIntPt * nstatev(1)
     $                         + iTIntPt * nstatev(2)
     $                         + iTIntPt * nstatev(3)
     $                         +    1    * nstatev(4)
     $                         +    3  ) = SB_ct(2)

                  svars(kblock,  iTIntPt * nstatev(1)
     $                         + iTIntPt * nstatev(2)
     $                         + iTIntPt * nstatev(3)
     $                         +    1    * nstatev(4)
     $                         +    4  ) = SB_ct(3)

                  svars(kblock,  iTIntPt * nstatev(1)
     $                         + iTIntPt * nstatev(2)
     $                         + iTIntPt * nstatev(3)
     $                         +    1    * nstatev(4)
     $                         +    5  ) = SB_vol


C                 d) Calculate the dinstances between each of the EIGHT (8)
C                    Gauss points and the shear band surface, then calculate
C                    the weight-coefficient foe each Gauss point
                  call non_plane_normal(SB_ct,natural_coords,
     $                                  natural_normalCt)

                  dist_ave = 0.d0
                  DO l = 1, iTIntPt
                    dist(l) = abs( dot_product( rXiEtaPt(:,l)-SB_ct,
     $                                            natural_normalCt ) )
                    dist_ave = dist_ave + dist(l)
                  END DO
                  dist_ave = dist_ave / 8.d0
                  DO l = 1, iTIntPt
                    Gt_weight(l) = (2.d0 - dist(l)/dist_ave )
                    IF ( (SB_vol*Gt_weight(l)/rElmtVol).LE.1.d-2) THEN
                      iOutPr(1) = elmt_number
                      iOutPr(2) = l
                      rOutPr(1) = SB_vol*Gt_weight(l)/rElmtVol
                      CALL XPLB_ABQERR(1,'At Elmt %I Gauss Point %I, '//
     $                'weight coefficient %R is small, '//
     $                'it only contributes to MATRIX',
     $                 iOutPr,rOutPr,' ')
                    END IF

                    IF (rElmtVol/SB_vol/Gt_weight(l).LT.1.d0) THEN
                      Gt_weight(l) = 1.d0
                    END IF

                    svars(kblock,  iTIntPt * nstatev(1)
     $                           + iTIntPt * nstatev(2)
     $                           + iTIntPt * nstatev(3)
     $                           +    1    * nstatev(4)
     $                           +  5 + l  ) = Gt_weight(l)
                  END DO

                END IF ! iCvFlag.GT.1 & elmt is in the front list, 
                       ! localization detected just now at the elmt centroid
              END IF ! iLocalFlg_Ct.eq.0
            END IF  ! Formulation with localization, iVfortype = 1 or 3 

C>>>        Numerical integration loop

            iStLocF = 0

            DO l = 1, iTIntPt

C             Compute shape functs, their derivs (wrt global coords), and the Jacobian determinant
              CALL SHP_FUNC_3D_8N(rXiEtaPt(1,l),rCoords,1,
     $                            rNn,rdNndXi,rJdet)

C             Check Jacobian
              IF ( rJdet .LE. 0.d0) THEN
                CALL XPLB_ABQERR(-3,'Non-positive'//
     $            'Jacobian determinant in element VU1',0,0.d0,' ')
              END IF

C             Compute integration element (dOmega)
              rdOmega = rJdet * rWeightPt(l)

C>>>          Get strain or strain increment and B (or B-bar) matrix

C             Finite deformation
C             Current value of the B or B-bar operator for internal force calculation
              CALL Gstrain_3D(rue,rEpsilon,rNn,rNnCt,rdNndXi,rdNndXiCt,
     $                        mcrd,nnode,iVelmType,iVforType,1,rStrain,
     $                        rSpin,rBmat,rBskw)

C             (Mid-step) strain increment and B or B-bar operator for constitutive update
              CALL SHP_FUNC_3D_8N(rXiEtaPt(1,l),rMSCoor,1,rNn,
     $                            rMSdNndXi,rMSJdet)
              CALL Gstrain_3D(rue,rEpsilon,rNn,rNnCt,rMSdNndXi,
     $                        rMSdNndXiCt,mcrd,nnode,iVelmType,
     $                        iVforType,1, rStrain,rSpin,rMSBm,
     $                        rMSBs)
C             Until now, rStrain and rSpin are quantities wrt to mid-step    

C             Compute volumetric strain increment (for bulk viscosity calculation)
              rEpVolInc = rStrain(1) + rStrain(2) + rStrain(3)

C             If formulation w/ localization
              IF ((iVforType.EQ.1).OR.(iVforType.EQ.3)) THEN

C               If first increment of first step, initialize localization flag in history
                IF ( ( kstep .EQ. 1) .AND. ( kinc .EQ. 1) ) THEN
                  svars(kblock, iTIntPt*(nstatev(1)+nstatev(2))
     $                         +  (l-1)* nstatev(3) + 1) = 0.d0

                  iumatsvst = 1494 + (l - 1) * nstatev(6)
                  iumatsven = 1493 + l * nstatev(6)
                  svars(kblock, iumatsvst : iumatsven) = 
     &                                  reshape(rImat, (/9/))
                END IF

                call comp_vel_gr(rdVelGrd, rStrain, rSpin)
C               Retrieve previous Def Grad and form next
                if (kinc.gt.1) then
                  iumatsvst = 1494 + (l - 1) * nstatev(6)
                  iumatsven = 1493 + l * nstatev(6)
                  rDfgrd0 =
     &            reshape(svars(kblock, iumatsvst : iumatsven), (/3,3/))
                  rDfgrd1 =
     &            matmul(rImat + rdVelGrd, rDfgrd0)
                end if

                svars(kblock, iumatsvst : iumatsven) = 
     &                                  reshape(rDfgrd1, (/9/))
C               Read localization status from history
                iLocalFlg = nint( svars(kblock,
     $                            iTIntPt*(nstatev(1)+nstatev(2))
     $                           +  (l-1)* nstatev(3) + 1) )

              ELSE    ! form w/o localization

                iLocalFlg = 0

              END IF   ! form w/ or w/o localization

              IF (iLocalFlg.EQ.0) THEN  ! form w/o localization, or
                                        ! form w/  localization, but no bifurcation yet

C               Select the correct state/history variables
                DO i = 1, nstatev(1)
                  hvars(i) = svars(kblock,((l-1)*nstatev(1))+i)
                END DO
                
                iumatsvst = 540 + (l - 1) * nstatev(5) 
                iumatsven = 539 + l * nstatev(5) 
                call umat_stress(rStress,
     &                           svars(kblock, iumatsvst:iumatsven),
     &                           rGsCoords(:, l), rDfgrd0,
     &                           rDfgrd1, props(77:118),
     &                          42, nstatev(5), nnode, elmt_number,
     &                          l, mcrd, time, dtimeCur, kstep, kinc)
C               Compute the stress (now rStrain and rSpin are at mid-step)
                CALL Gstress_3D(rStrain,rSpin,hvars,nstatev(1),props,
     $                          nprops,iVelmType,iVforType,iConsType,
     $                          time,dtimeCur,rStress,rStrInc,rCmod,
     $                          rGtan,rCd,iCvFlag) 

!               DEBUGGING...
                IF (iCvFlag.EQ.0) THEN     ! Stress update has convergence problems
                  write(*,'(A,I6,A,I2,A)') "Stress update fail, Elmt",
     $                     jelem(kblock),", Gauss Pt.",l,":"
                  write(*,'("Strain increment: ",6(/,6X,EN16.6))')
     $                 ( rStrain(i), i=1,6 )
                  write(*,'("Spin increment: ",3(/,6X,EN16.6))')
     $                 ( rSpin(i), i=1,3 )
                  write(*,'("Old stress: ",6(/,6X,EN16.6))')
     $                 ( svars(kblock,((l-1)*nstatev(1))+i), i=1,6 )
                  write(*,'("History variables: ",4(/,6X,EN16.6))')
     $                 ( svars(kblock,((l-1)*nstatev(1))+i), i=7,10)
                END IF
!               ...GNIGGUBED

C               Update the maximum wave speed in the current element

                rMaxWavSpd = max(rMaxWavSpd,rCd)

C               Check whether accumulated effective plastic strain increases monotonically
                IF ( svars(kblock,7 + (l-1)*nstatev(1)  ).GT.hvars(7) ) THEN
                  rtemp(1) = svars( kblock,7 + (l-1)*nstatev(1) )
                  rtemp(2) = hvars(7)
                  CALL XPLB_ABQERR(-3,'Accumulated plastic strain '//
     $                             'fails to increase monotonically,'//
     $                'prev_step = %R, current_step = %R',0, rtemp,' ')
                END IF

C               Write back the state/history variables
                DO i = 1, nstatev(1)
                  svars(kblock,((l-1)*nstatev(1))+i) = hvars(i)
                END DO

C               Compute linear + quadratic bulk viscosity terms
                !! Tao,The quadratic bulk viscosity is applied 
                !! only if the volumetric strain rate is compressive.
                IF (dtimeCur .NE. 0.d0) THEN
                  rBVpress = -rRho *
     $                  rbv1*rCd*rElmtSize*(rEpVolInc/dtimeCur)
                  IF (rEpVolInc.LT.0.d0) THEN
                    rBVpress = rBVpress 
     $                -rRho*(rbv2*rElmtSize*rEpVolInc/dtimeCur)**2
                  END IF
                ELSE
                  rBVpress = 0.d0
                END IF

C               Compute the internal force vector 
                DO n = 1, nnode
                  DO i = 1, mcrd
                    DO j = 1, 6
                      IF (j.LE.3) THEN    ! add bulk viscosity 
                        rSigma = rStress(j) - rBVpress
                      ELSE
                        rSigma = rStress(j)
                      END IF
                      rhs( kblock , (n-1)*mcrd + i ) = 
     $                rhs( kblock , (n-1)*mcrd + i )
     $                +  rBmat( j , (n-1)*mcrd + i ) * rSigma * rdOmega
                    END DO
                  END DO
                END DO

C               use the localization flag at the elmt centroid
                IF (iLocalFlg_Ct.EQ.1) THEN

C                 Retrieve the volume of the shear band
                  SB_vol = svars(kblock,
     $                            iTIntPt * nstatev(1)
     $                          + iTIntPt * nstatev(2)
     $                          + iTIntPt * nstatev(3)
     $                          +    1    * nstatev(4) + 5)

                  Gt_temp = svars(kblock,
     $                            iTIntPt * nstatev(1)
     $                          + iTIntPt * nstatev(2)
     $                          + iTIntPt * nstatev(3)
     $                          +    1    * nstatev(4) + 5 + l)

C                 IF the contribution to the shear band is too small, 
C                 treat this Gauss point as matrix ONLY
                  IF ((SB_vol*Gt_temp/rElmtVol).GT.1.d-2) THEN

C                 We need to:

C                 a) calculate the initial band vectors (n0 and m0)
                  iOutPr(1) = jelem(kblock)
                  iOutPr(2) = l
                  iOutPr(3) = kinc
                  rOutPr(1) = time(1)
                  rOutPr(2) = time(2)
                  CALL perturb_3D(rStress,iOutPr,rOutPr,rNvect,rMvect)

C                 b) Save localization flag and initial band vectors (n0 & m0) to history
                  svars(kblock, iTIntPt * (nstatev(1)+nstatev(2))
     $                          + (l-1) * nstatev(3) + 1) = 1.d0

                  svars(kblock, iTIntPt*(nstatev(1)+nstatev(2))
     $                         +  (l-1)* nstatev(3) + 2) = rNvect(1)
                  svars(kblock, iTIntPt*(nstatev(1)+nstatev(2))
     $                         +  (l-1)* nstatev(3) + 3) = rNvect(2)
                  svars(kblock, iTIntPt*(nstatev(1)+nstatev(2))
     $                         +  (l-1)* nstatev(3) + 4) = rNvect(3)

                  svars(kblock, iTIntPt*(nstatev(1)+nstatev(2))
     $                         +  (l-1)* nstatev(3) + 6) = rMvect(1)
                  svars(kblock, iTIntPt*(nstatev(1)+nstatev(2))
     $                         +  (l-1)* nstatev(3) + 7) = rMvect(2)
                  svars(kblock, iTIntPt*(nstatev(1)+nstatev(2))
     $                         +  (l-1)* nstatev(3) + 8) = rMvect(3)

                  iumatsvst = 1494 + (l - 1) * nstatev(6)
                  iumatsven = 1493 + l * nstatev(6)
                  svars(kblock, iumatsvst : iumatsven) = 
     &                                  reshape(rDfgrd1, (/9/))

                  iumatsvst = 1566 + (l - 1) * nstatev(6)
                  iumatsven = 1565 + l * nstatev(6)
                  svars(kblock, iumatsvst : iumatsven) = 
     &                                  reshape(rDfgrd1, (/9/))

C                 c) Save current config (X*) to history (finite deformation case)
                  j = 9
                  DO n = 1, nnode
                    DO i = 1, mcrd
                      svars(kblock, iTIntPt*(nstatev(1)+nstatev(2))
     $                      +  (l-1)* nstatev(3) + j) = rCoords(i,n)
                      j = j + 1
                    END DO
                  END DO

C                 d) Initialize state variables for band material
                  SELECT CASE (iConsType)
                    CASE DEFAULT          ! same material inside and outside band
                      DO i = 1, nstatev(1)  ! in this case, nstatev(2) = nstatev(1)
                        svars(kblock,(iTIntPt*nstatev(1)
     $                               +  (l-1)*nstatev(2)) + i)
     $                  =  svars(kblock,((l-1)*nstatev(1)) + i)
                      END DO
                    CASE (4)              ! models (1) outside and (3) inside localization band
                      CONTINUE    ! insert code here
                    CASE (5)              ! models (2) outside and (3) inside localization band
                      CONTINUE    ! insert code here
                  END SELECT

C                 e) Initialize (perturb) rAlphaN and save it in history
                  svars(kblock, iTIntPt*(nstatev(1)+nstatev(2))
     $                         +  (l-1)* nstatev(3) + 5) = 1.d0
                  END IF ! SB_vol*Gt_temp/rElmt_vol > 1.d-2
                END IF  ! iLocalFlg_Ct.eq.1

C>>>          Compute stress and internal force (using modified quadrature scheme)

C             Form w/ localization, and bifurcation already took place at current point 
              ELSE IF (iLocalFlg.EQ.1) THEN

C               Element has localization at any Gauss point --> set flag to make time step smaller 
                iStLocF = 1

C               Get initial orientation of band from history
                rNvect(1) = svars(kblock,
     $                              iTIntPt*(nstatev(1)+nstatev(2))
     $                             +  (l-1)* nstatev(3) + 2)
                rNvect(2) = svars(kblock,
     $                              iTIntPt*(nstatev(1)+nstatev(2))
     $                             +  (l-1)* nstatev(3) + 3)
                rNvect(3) = svars(kblock,
     $                              iTIntPt*(nstatev(1)+nstatev(2))
     $                             +  (l-1)* nstatev(3) + 4)

C               Retrieve the volume of the shear band
                SB_vol = svars(kblock,
     $                          iTIntPt * nstatev(1)
     $                        + iTIntPt * nstatev(2)
     $                        + iTIntPt * nstatev(3)
     $                        +    1    * nstatev(4) + 5)

                Gt_temp = svars(kblock,
     $                          iTIntPt * nstatev(1)
     $                        + iTIntPt * nstatev(2)
     $                        + iTIntPt * nstatev(3)
     $                        +    1    * nstatev(4) + 5 + l)

C               Finite deformation case, calculate up-to-date orientation of band
C               read X* from history
                j = 9
                DO n = 1, nnode
                  DO i = 1, mcrd
                    rBifCoor(i,n) = svars(kblock, 
     $                     iTIntPt*(nstatev(1)+nstatev(2))
     $                    +  (l-1)* nstatev(3) + j)
                    j = j + 1
                  END DO
                END DO

C               get the INVERSE of F*
                CALL GdefGrad_3D(rdNndXi,rBifCoor,mcrd,
     $                           nnode,iVelmType,rBifF)

C               calculate up-to-date band orientation vector
                CALL GnewNvec_3D(rBifF, rNvect)

C               Get alpha_N and rMvect initial guesses (converged values from last time step)

                rAlphaN   = svars(kblock,
     $                              iTIntPt*(nstatev(1)+nstatev(2))
     $                             +  (l-1)* nstatev(3) + 5)

                IF (abs(rAlphaN).LT.1.d-2) THEN
                  rAlphaN = 1.d0
                END IF

                rMvect(1) = svars(kblock,
     $                              iTIntPt*(nstatev(1)+nstatev(2))
     $                             +  (l-1)* nstatev(3) + 6)
                rMvect(2) = svars(kblock,
     $                              iTIntPt*(nstatev(1)+nstatev(2))
     $                             +  (l-1)* nstatev(3) + 7)
                rMvect(3) = svars(kblock,
     $                              iTIntPt*(nstatev(1)+nstatev(2))
     $                             +  (l-1)* nstatev(3) + 8)

                iAltern   = 0
                iHaveSol  = 0
                iDbgFlg   = 0
                iAltMax   = 15

10              CONTINUE

C               Calculate the updated localization mode tensor (T)
                rTtens(1) =  rMvect(1)*rNvect(1)
                rTtens(2) =  rMvect(2)*rNvect(2)
                rTtens(3) =  rMvect(3)*rNvect(3)   ! T_33
                rTtens(4) = 0.5d0*( rMvect(2)*rNvect(3)
     $                            + rMvect(3)*rNvect(2) )   ! 1/2*(T_23+T_32)
                rTtens(5) = 0.5d0*( rMvect(1)*rNvect(3)
     $                            + rMvect(3)*rNvect(1) )   ! 1/2*(T_13+T_31)
                rTtens(6) = 0.5d0*( rMvect(1)*rNvect(2)
     $                            + rMvect(2)*rNvect(1) )   ! 1/2*(T_12+T_21)

                rTtensk(1)= 0.5d0*( rMvect(3)*rNvect(2)
     $                            - rMvect(2)*rNvect(3) )   ! 1/2*(T_32-T_23)
                rTtensk(2)= 0.5d0*( rMvect(1)*rNvect(3)
     $                            - rMvect(3)*rNvect(1) )   ! 1/2*(T_13-T_31)
                rTtensk(3)= 0.5d0*( rMvect(2)*rNvect(1)
     $                            - rMvect(1)*rNvect(2) )   ! 1/2*(T_21-T_12)

C               Calculate L:T and the norm of L, L = grad(v)
                rTEps   = 0.d0
                rTSpn   = 0.d0
                rLnrm   = 0.d0
                DO k = 1, 3
                  rTEps = rTEps + rTtens(k)*rStrain(k)
                  rTSpn = rTSpn + 2.d0*rTtensk(k)*rSpin(k)
                  rLnrm = rLnrm + rStrain(k)**2
                END DO
                DO k = 4, 6
                  rTEps = rTEps + rTtens(k)*rStrain(k)
                  rLnrm = rLnrm + ( 0.5d0*rStrain(k) + rSpin(k-3) )**2
     $                          + ( 0.5d0*rStrain(k) - rSpin(k-3) )**2
                END DO
                rTL   = rTEps + rTSpn
                rLnrm = sqrt(rLnrm)

                IF (iAltern.EQ.iAltMax) THEN  ! we have rAlphaHat, not rAlphaN
                  IF (abs(rTL).GT.0.d0) THEN
                    rAlphaN = rAlphaHat / rTL
                  ELSE
                    rAlphaN = 0.d0
                  END IF
                ElSE
                  IF (rTl.GT.0.d0) THEN
                    rAlphaN = abs(rAlphaN)
                  ELSE
                    rAlphaN =-abs(rAlphaN)
                  END IF
                  rAlphaHat = rAlphaN * rTL
                END IF

                DO i = 1, 3
                  rMhat(i)= rAlphaHat * rMvect(i)
                END DO

C               Calculate L.n and its norm, L = grad(v)
                rLn(1) =          rStrain(1)               * rNvect(1)
     $                  + ( 0.5d0*rStrain(6) - rSpin(3) )  * rNvect(2)
     $                  + ( 0.5d0*rStrain(5) + rSpin(2) )  * rNvect(3)

                rLn(2) =  ( 0.5d0*rStrain(6) + rSpin(3) )  * rNvect(1)
     $                  +         rStrain(2)               * rNvect(2)
     $                  + ( 0.5d0*rStrain(4) - rSpin(1) )  * rNvect(3)

                rLn(3) =  ( 0.5d0*rStrain(5) - rSpin(2) )  * rNvect(1)
     $                  + ( 0.5d0*rStrain(4) + rSpin(1) )  * rNvect(2)
     $                  +         rStrain(3)               * rNvect(3)

                rLnNrm = sqrt ( rLn(1)**2 + rLn(2)**2 + rLn(3)**2 )    ! norm

C               Iteratively find rAlphaN and rMvect

                iConverge = 0
                iIterMax  = 30
                rDelNrm   = 1.d0

                rIterTol  = max(rYoungsMod*1.d-14,1.d-15)
                rRelTol   = 1.d-12

C               If L is approximately zero or orthogonal to the band, there's no
C               physically meaningful solution for the band parameters
                IF ( rLnrm.LT.(dtimeCur*1.d-16) ) THEN
                  iConverge = 1    ! Don't iterate
                ELSEIF ( (rLnNrm/rLnrm).LT.1.d-6 ) THEN
                  iConverge = 1    ! Don't iterate
                END IF

C               If we've decided to use the best known parameters
                IF (iAltern.EQ.iAltMax) THEN
                  iConverge = 1    ! Don't iterate
                END IF

                DO iCount = 1, iIterMax

                  rtemp(1) = rElmtVol/SB_vol/Gt_temp - 1.d0
C                  if (rtemp(1).LE.0.d0) then
C                    CALL XPLB_ABQERR(-3,'The coefficient between '//
C     $             'Alpha_B and Alpha_M is negative, %R',0,rtemp(1),' ')
C                  end if
                  rAlphaBar = rtemp(1) * rAlphaHat    ! Compatibility condition

C                 If we're finished iterating and have a converged solution
                  IF ((iCount.GT.1).AND.(iConverge.EQ.1)) THEN

C                   Update T:L
                    rTEps   = 0.d0
                    rTSpn   = 0.d0
                    DO k = 1, 3
                      rTEps = rTEps + rTtens(k)*rStrain(k)
                      rTSpn = rTSpn + 2.d0*rTtensk(k)*rSpin(k)
                    END DO
                    DO k = 4, 6
                      rTEps = rTEps + rTtens(k)*rStrain(k)
                    END DO
                    rTL  = rTEps + rTSpn

C                   Calculate rAlphaN
                    IF (abs(rTL).GT.0.d0) THEN       ! avoid div. by zero
                      rAlphaN = rAlphaHat / rTL
                    ELSE                             ! L and T are orthogonal 
                      rAlphaN = 0.d0
                    END IF

                  END IF

C                 Contributions from the Localized zone:
C                 -------------------------------------

C                 Select the correct state/history variables
                  DO i = 1, nstatev(2)
                    hvars(i) = svars(kblock,(iTIntPt*nstatev(1)
     $                                       + (l-1)*nstatev(2))+i)
                  END DO

C                 Compute the modified strain and spin
                  DO j = 1, 6
                    IF (j.LE.3) THEN
                      rTfactor = 1.d0
                      rSpinL(j)= rSpin(j) + rAlphaBar * rTtensk(j)
                    ELSE
                      rTfactor = 2.d0
                    END IF
                    rStrnL(j)  = rStrain(j)                           !   continuous strain incr.
     $                         + rAlphaBar * rTfactor*rTtens(j)       ! + discontinuous incr.
                  END DO  ! j

C                 Compute the stress inside the localization band
                  CALL Gstress_3D(rStrnL,rSpinL,hvars,nstatev(2),props,
     $                         nprops,iVelmType,
     $                         iVforType,iConsType,time,dtimeCur,
     $                         rStress,rStrInc,rCmod,rGtan,rCd,iCvFlag)

                  IF (iCvFlag.EQ.0) THEN     ! Stress update has convergence problems

!                   DEBUGGING...
                    write(*,'(A,I6,A,I2,A)') "Stress update fail, Elmt",
     $                           jelem(kblock),", Gauss Pt.",l," (L):"
                    write(*,'("Strain increment: ",6(/,6X,EN16.6))')
     $                   ( rStrnL(i), i=1,6 )
                    write(*,'("Spin increment: ",3(/,6X,EN16.6))')
     $                   ( rSpinL(i), i=1,3 )
                    write(*,'("Old stress: ",6(/,6X,EN16.6))')
     $                   ( svars(kblock,(iTIntPt*nstatev(1)
     $                                   + (l-1)*nstatev(2))+i), i=1,6 )
                    write(*,'("History variables: ",4(/,6X,EN16.6))')
     $                   ( svars(kblock,(iTIntPt*nstatev(1)
     $                                   + (l-1)*nstatev(2))+i), i=7,10)
!                   ...GNIGGUBED

                    iConverge = 0  !  This initial guess has failed
                    EXIT           !  Try a different initial guess
                  END IF

C                 If we're done iterating, we need to:
                  IF (iConverge.EQ.1) THEN

C                   a) Update the maximum wave speed in the current element
                    rMaxWavSpd = max(rMaxWavSpd,rCd)

                    IF ( svars( kblock, 7 + iTIntPt*nstatev(1) +
     $                  (l-1)*nstatev(2) ).GT.hvars(7) ) THEN
                     rtemp(1) = svars( kblock, 7 + 
     $                           iTIntPt*nstatev(1) + (l-1)*nstatev(2) )
                     rtemp(2) = hvars(7)
                     CALL XPLB_ABQERR(-3,'Accumulated plastic strain '//
     $                              'fails to increase monotonically,'//
     $                 'prev_step = %R, current_step = %R',0, rtemp,' ')
                    END IF

C                   b) Write back the material state/history variables
                    DO i = 1, nstatev(2)
                      svars(kblock,(iTIntPt*nstatev(1)
     $                             +  (l-1)*nstatev(2))+i) = hvars(i)
                    END DO

C                   c) Compute linear + quadratic bulk viscosity terms
                    !! The quadratic bulk viscosity is applied 
                    !! only if the volumetric strain rate is compressive.
                    IF (dtimeCur .NE. 0.d0) THEN
                      rBVpress = -rRho *
     $                         rbv1*rCd*rElmtSize*(rEpVolInc/dtimeCur)
                      IF (rEpVolInc.LT.0.d0) THEN
                        rBVpress = rBVpress
     $                    -rRho*(rbv2*rElmtSize*rEpVolInc/dtimeCur)**2
                      END IF
                    ELSE
                      rBVpress = 0.d0
                    END IF

C                   d) Compute the internal force vector contribution
                    rAlphaL   = ( rElmtVol/SB_vol/Gt_temp - 1.d0 ) 
     $                        * rAlphaN      ! Compatibility condition
                    DO n = 1, nnode
                      DO i = 1, mcrd
                        DO j = 1, 6
                          rBmtL(j,(n-1)*mcrd+i) = rBmat(j,(n-1)*mcrd+i)
                          IF (j.LE.3) THEN
                            rTfactor  = 1.d0
                          ELSE
                            rTfactor  = 2.d0
                          END IF
                          DO k = 1, 3
                            rBmtL(j,(n-1)*mcrd+i)= rBmtL(j,(n-1)*mcrd+i)
     $                          + rAlphaL
     $                          *( rTtens(k) *rBmat(k,(n-1)*mcrd+i)
     $                            +rTtensk(k)*rBskw(k,(n-1)*mcrd+i) )
     $                          *rTfactor*rTtens(j)
                          END DO ! k
                          DO k = 4, 6
                            rBmtL(j,(n-1)*mcrd+i)= rBmtL(j,(n-1)*mcrd+i)
     $                          + rAlphaL
     $                          *rTtens(k)*rBmat(k,(n-1)*mcrd+i)
     $                          *rTfactor*rTtens(j)
                          END DO ! k
                        END DO  ! j = 1,6 
                      END DO  ! i = mcrd
                    END DO ! n = nnode

                    DO n = 1, nnode
                      DO i = 1, mcrd
                        DO j = 1, 6
                          IF (j.LE.3) THEN    ! add bulk viscosity 
                            rSigma = rStress(j) - rBVpress
                          ELSE
                            rSigma = rStress(j)
                          END IF
                          rhs( kblock , (n-1)*mcrd+i ) =
     $                       rhs( kblock , (n-1)*mcrd+i )
     $                     + SB_vol/rElmtVol*Gt_temp
     $                     * rBmtL(j,(n-1)*mcrd+i) * rSigma * rdOmega
                        END DO
                      END DO
                    END DO

C                 If we're still iterating, we need to:
                  ELSE

C                   Calculate residual for the traction continuity equations
                    rG(1) = - rNvect(1)*rStress(1)
     $                      - rNvect(2)*rStress(6)
     $                      - rNvect(3)*rStress(5)
                    rG(2) = - rNvect(1)*rStress(6)
     $                      - rNvect(2)*rStress(2)
     $                      - rNvect(3)*rStress(4)
                    rG(3) = - rNvect(1)*rStress(5)
     $                      - rNvect(2)*rStress(4)
     $                      - rNvect(3)*rStress(3)

C                   Calculate tangent matrix for the traction continuity equations

                    DO i = 1, 6
                      DO j = 1, 3
                        rCmod(i,j) = ( rElmtVol/SB_vol/Gt_temp - 1.d0 )
     $                             * rCmod(i,j)
                        rGtan(i,j) = ( rElmtVol/SB_vol/Gt_temp - 1.d0 ) 
     $                             * rGtan(i,j)
                      END DO
                      DO j = 4, 6
                        rCmod(i,j) = ( rElmtVol/SB_vol/Gt_temp - 1.d0 )
     $                             * rCmod(i,j)
                      END DO
                    END DO

C                   Whether G term has '0.5d0' or not?
                    rnCn(1,1) =   rNvect(1) * rCmod(1,1) * rNvect(1)
     $                          + rNvect(1) * rCmod(1,6) * rNvect(2)
     $                          + rNvect(1) * rCmod(1,5) * rNvect(3)
     $                          + rNvect(2) * rCmod(6,1) * rNvect(1)
     $                          + rNvect(2) * rCmod(6,6) * rNvect(2)
     $                          + rNvect(2) * rCmod(6,5) * rNvect(3)
     $                          + rNvect(3) * rCmod(5,1) * rNvect(1)
     $                          + rNvect(3) * rCmod(5,6) * rNvect(2)
     $                          + rNvect(3) * rCmod(5,5) * rNvect(3)
     $                          - rNvect(1) * rGtan(1,3) * rNvect(2) * 0.5d0
     $                          + rNvect(1) * rGtan(1,2) * rNvect(3) * 0.5d0
     $                          - rNvect(2) * rGtan(6,3) * rNvect(2) * 0.5d0
     $                          + rNvect(2) * rGtan(6,2) * rNvect(3) * 0.5d0
     $                          - rNvect(3) * rGtan(5,3) * rNvect(2) * 0.5d0
     $                          + rNvect(3) * rGtan(5,2) * rNvect(3) * 0.5d0


                    rnCn(2,2) =   rNvect(1) * rCmod(6,6) * rNvect(1)
     $                          + rNvect(1) * rCmod(6,2) * rNvect(2)
     $                          + rNvect(1) * rCmod(6,4) * rNvect(3)
     $                          + rNvect(2) * rCmod(2,6) * rNvect(1)
     $                          + rNvect(2) * rCmod(2,2) * rNvect(2)
     $                          + rNvect(2) * rCmod(2,4) * rNvect(3)
     $                          + rNvect(3) * rCmod(4,6) * rNvect(1)
     $                          + rNvect(3) * rCmod(4,2) * rNvect(2)
     $                          + rNvect(3) * rCmod(4,4) * rNvect(3)
     $                          + rNvect(1) * rGtan(6,3) * rNvect(1) * 0.5d0
     $                          - rNvect(1) * rGtan(6,1) * rNvect(3) * 0.5d0
     $                          + rNvect(2) * rGtan(2,3) * rNvect(1) * 0.5d0
     $                          - rNvect(2) * rGtan(2,1) * rNvect(3) * 0.5d0
     $                          + rNvect(3) * rGtan(4,3) * rNvect(1) * 0.5d0
     $                          - rNvect(3) * rGtan(4,1) * rNvect(3) * 0.5d0

                    rnCn(3,3) =   rNvect(1) * rCmod(5,5) * rNvect(1)
     $                          + rNvect(1) * rCmod(5,4) * rNvect(2)
     $                          + rNvect(1) * rCmod(5,3) * rNvect(3)
     $                          + rNvect(2) * rCmod(4,5) * rNvect(1)
     $                          + rNvect(2) * rCmod(4,4) * rNvect(2)
     $                          + rNvect(2) * rCmod(4,3) * rNvect(3)
     $                          + rNvect(3) * rCmod(3,5) * rNvect(1)
     $                          + rNvect(3) * rCmod(3,4) * rNvect(2)
     $                          + rNvect(3) * rCmod(3,3) * rNvect(3)
     $                          - rNvect(1) * rGtan(5,2) * rNvect(1) * 0.5d0
     $                          + rNvect(1) * rGtan(5,1) * rNvect(2) * 0.5d0
     $                          - rNvect(2) * rGtan(4,2) * rNvect(1) * 0.5d0
     $                          + rNvect(2) * rGtan(4,1) * rNvect(2) * 0.5d0
     $                          - rNvect(3) * rGtan(3,2) * rNvect(1) * 0.5d0
     $                          + rNvect(3) * rGtan(3,1) * rNvect(2) * 0.5d0

                    rnCn(1,2) =   rNvect(1) * rCmod(1,6) * rNvect(1)
     $                          + rNvect(1) * rCmod(1,2) * rNvect(2)
     $                          + rNvect(1) * rCmod(1,4) * rNvect(3)
     $                          + rNvect(2) * rCmod(6,6) * rNvect(1)
     $                          + rNvect(2) * rCmod(6,2) * rNvect(2)
     $                          + rNvect(2) * rCmod(6,4) * rNvect(3)
     $                          + rNvect(3) * rCmod(5,6) * rNvect(1)
     $                          + rNvect(3) * rCmod(5,2) * rNvect(2)
     $                          + rNvect(3) * rCmod(5,4) * rNvect(3)
     $                          + rNvect(1) * rGtan(1,3) * rNvect(1) * 0.5d0
     $                          - rNvect(1) * rGtan(1,1) * rNvect(3) * 0.5d0
     $                          + rNvect(2) * rGtan(6,3) * rNvect(1) * 0.5d0
     $                          - rNvect(2) * rGtan(6,1) * rNvect(3) * 0.5d0
     $                          + rNvect(3) * rGtan(5,3) * rNvect(1) * 0.5d0
     $                          - rNvect(3) * rGtan(5,1) * rNvect(3) * 0.5d0

                    rnCn(2,1) =   rNvect(1) * rCmod(6,1) * rNvect(1)
     $                          + rNvect(1) * rCmod(6,6) * rNvect(2)
     $                          + rNvect(1) * rCmod(6,5) * rNvect(3)
     $                          + rNvect(2) * rCmod(2,1) * rNvect(1)
     $                          + rNvect(2) * rCmod(2,6) * rNvect(2)
     $                          + rNvect(2) * rCmod(2,5) * rNvect(3)
     $                          + rNvect(3) * rCmod(4,1) * rNvect(1)
     $                          + rNvect(3) * rCmod(4,6) * rNvect(2)
     $                          + rNvect(3) * rCmod(4,5) * rNvect(3)
     $                          - rNvect(1) * rGtan(6,3) * rNvect(2) * 0.5d0
     $                          + rNvect(1) * rGtan(6,2) * rNvect(3) * 0.5d0
     $                          - rNvect(2) * rGtan(2,3) * rNvect(2) * 0.5d0
     $                          + rNvect(2) * rGtan(2,2) * rNvect(3) * 0.5d0
     $                          - rNvect(3) * rGtan(4,3) * rNvect(2) * 0.5d0
     $                          + rNvect(3) * rGtan(4,2) * rNvect(3) * 0.5d0

                    rnCn(1,3) =   rNvect(1) * rCmod(1,5) * rNvect(1)
     $                          + rNvect(1) * rCmod(1,4) * rNvect(2)
     $                          + rNvect(1) * rCmod(1,3) * rNvect(3)
     $                          + rNvect(2) * rCmod(6,5) * rNvect(1)
     $                          + rNvect(2) * rCmod(6,4) * rNvect(2)
     $                          + rNvect(2) * rCmod(6,3) * rNvect(3)
     $                          + rNvect(3) * rCmod(5,5) * rNvect(1)
     $                          + rNvect(3) * rCmod(5,4) * rNvect(2)
     $                          + rNvect(3) * rCmod(5,3) * rNvect(3)
     $                          - rNvect(1) * rGtan(1,2) * rNvect(1) * 0.5d0
     $                          + rNvect(1) * rGtan(1,1) * rNvect(2) * 0.5d0
     $                          - rNvect(2) * rGtan(6,2) * rNvect(1) * 0.5d0
     $                          + rNvect(2) * rGtan(6,1) * rNvect(2) * 0.5d0
     $                          - rNvect(3) * rGtan(5,2) * rNvect(1) * 0.5d0
     $                          + rNvect(3) * rGtan(5,1) * rNvect(2) * 0.5d0


                    rnCn(3,1) =   rNvect(1) * rCmod(5,1) * rNvect(1)
     $                          + rNvect(1) * rCmod(5,6) * rNvect(2)
     $                          + rNvect(1) * rCmod(5,5) * rNvect(3)
     $                          + rNvect(2) * rCmod(4,1) * rNvect(1)
     $                          + rNvect(2) * rCmod(4,6) * rNvect(2)
     $                          + rNvect(2) * rCmod(4,5) * rNvect(3)
     $                          + rNvect(3) * rCmod(3,1) * rNvect(1)
     $                          + rNvect(3) * rCmod(3,6) * rNvect(2)
     $                          + rNvect(3) * rCmod(3,5) * rNvect(3)
     $                          - rNvect(1) * rGtan(5,3) * rNvect(2) * 0.5d0
     $                          + rNvect(1) * rGtan(5,2) * rNvect(3) * 0.5d0
     $                          - rNvect(2) * rGtan(4,3) * rNvect(2) * 0.5d0
     $                          + rNvect(2) * rGtan(4,2) * rNvect(3) * 0.5d0
     $                          - rNvect(3) * rGtan(3,3) * rNvect(2) * 0.5d0
     $                          + rNvect(3) * rGtan(3,2) * rNvect(3) * 0.5d0

                    rnCn(2,3) =   rNvect(1) * rCmod(6,5) * rNvect(1)
     $                          + rNvect(1) * rCmod(6,4) * rNvect(2)
     $                          + rNvect(1) * rCmod(6,3) * rNvect(3)
     $                          + rNvect(2) * rCmod(2,5) * rNvect(1)
     $                          + rNvect(2) * rCmod(2,4) * rNvect(2)
     $                          + rNvect(2) * rCmod(2,3) * rNvect(3)
     $                          + rNvect(3) * rCmod(4,5) * rNvect(1)
     $                          + rNvect(3) * rCmod(4,4) * rNvect(2)
     $                          + rNvect(3) * rCmod(4,3) * rNvect(3)
     $                          - rNvect(1) * rGtan(6,2) * rNvect(1) * 0.5d0
     $                          + rNvect(1) * rGtan(6,1) * rNvect(2) * 0.5d0
     $                          - rNvect(2) * rGtan(2,2) * rNvect(1) * 0.5d0
     $                          + rNvect(2) * rGtan(2,1) * rNvect(2) * 0.5d0
     $                          - rNvect(3) * rGtan(4,2) * rNvect(1) * 0.5d0
     $                          + rNvect(3) * rGtan(4,1) * rNvect(2) * 0.5d0

                    rnCn(3,2) =   rNvect(1) * rCmod(5,6) * rNvect(1)
     $                          + rNvect(1) * rCmod(5,2) * rNvect(2)
     $                          + rNvect(1) * rCmod(5,4) * rNvect(3)
     $                          + rNvect(2) * rCmod(4,6) * rNvect(1)
     $                          + rNvect(2) * rCmod(4,2) * rNvect(2)
     $                          + rNvect(2) * rCmod(4,4) * rNvect(3)
     $                          + rNvect(3) * rCmod(3,6) * rNvect(1)
     $                          + rNvect(3) * rCmod(3,2) * rNvect(2)
     $                          + rNvect(3) * rCmod(3,4) * rNvect(3)
     $                          + rNvect(1) * rGtan(5,3) * rNvect(1) * 0.5d0
     $                          - rNvect(1) * rGtan(5,1) * rNvect(3) * 0.5d0
     $                          + rNvect(2) * rGtan(4,3) * rNvect(1) * 0.5d0
     $                          - rNvect(2) * rGtan(4,1) * rNvect(3) * 0.5d0
     $                          + rNvect(3) * rGtan(3,3) * rNvect(1) * 0.5d0
     $                          - rNvect(3) * rGtan(3,1) * rNvect(3) * 0.5d0

                  END IF   !  still iterating

C                 Contributions from the Non-localized zone:
C                 -----------------------------------------

C                 Select the correct state/history variables
                  DO i = 1, nstatev(1)
                    hvars(i) = svars(kblock,((l-1)*nstatev(1))+i)
                  END DO

C                 Compute the modified strain and spin

                  DO j = 1, 6
                    IF (j.LE.3) THEN
                      rTfactor = 1.d0
                      rSpinN(j)= rSpin(j) - rAlphaHat * rTtensk(j)
                    ELSE
                      rTfactor = 2.d0
                    END IF
                    rStrnN(j)  = rStrain(j)                           !   continuous strain incr.
     $                         - rAlphaHat * rTfactor*rTtens(j)       ! + discontinuous incr.
                  END DO  ! j

C                 If this is the first try (initial guess), only allow elastic behavior (unloading)
                  IF (iAltern.LE.1) then
                    iElFlag = 0
                  ELSE
                    iElFlag = iConsType
                  END IF

C                 Make sure DRX is suppressed outside the band
                  IF (props(72).GT.0.d0) THEN     ! If we're considering DRX
                    rDRXfac = -1.d0
                  ELSE
                    rDRXfac =  1.d0
                  END IF
                  props(72) = rDRXfac * props(72)

C                 Compute the stress outside the localization band
                  CALL Gstress_3D(rStrnN,rSpinN,hvars,nstatev(1),props,
     $                         nprops,iVelmType,
     $                         iVforType,iElFlag,time,dtimeCur,
     $                         rStress,rStrInc,rCmod,rGtan,rCd,iCvFlag)

                  props(72) = rDRXfac * props(72)  ! Re-enable DRX (inside band)

                  IF (iCvFlag.EQ.0) THEN     ! Stress update has convergence problems

!                   DEBUGGING...
                    write(*,'(A,I6,A,I2,A)') "Stress update fail, Elmt",
     $                           jelem(kblock),", Gauss Pt.",l," (N):"
                    write(*,'("Strain increment: ",6(/,6X,EN16.6))')
     $                   ( rStrnN(i), i=1,6 )
                    write(*,'("Spin increment: ",3(/,6X,EN16.6))')
     $                   ( rSpinN(i), i=1,3 )
                    write(*,'("Old stress: ",6(/,6X,EN16.6))')
     $                   ( svars(kblock,((l-1)*nstatev(1))+i), i=1,6 )
                    write(*,'("History variables: ",4(/,6X,EN16.6))')
     $                   ( svars(kblock,((l-1)*nstatev(1))+i), i=7,10)
!                   ...GNIGGUBED

                    iConverge = 0  !  This initial guess has failed
                    EXIT           !  Try a different initial guess
                  END IF

C                 If we're done iterating, we need to:
                  IF (iConverge.EQ.1) THEN

C                   a) Update the maximum wave speed in the current element
                    rMaxWavSpd = max(rMaxWavSpd,rCd)

                    IF ( svars(kblock,7 + 
     $                  (l-1)*nstatev(1) ).GT.hvars(7) ) THEN
                     rtemp(1) = svars(kblock,7 + (l-1)*nstatev(1) )
                     rtemp(2) = hvars(7)
                     CALL XPLB_ABQERR(-3,'Accumulated plastic strain '//
     $                              'fails to increase monotonically,'//
     $                 'prev_step = %R, current_step = %R',0, rtemp,' ')
                    END IF

C                   b) Write back the material state/history variables
                    DO i = 1, nstatev(1)
                      svars(kblock,((l-1)*nstatev(1))+i) = hvars(i)
                    END DO

c                   c) Compute linear + quadratic bulk viscosity terms
                    !! Tao,The quadratic bulk viscosity is applied 
                    !! only if the volumetric strain rate is compressive.
                    IF (dtimeCur .NE. 0.d0) THEN
                      rBVpress = -rRho *
     $                         rbv1*rCd*rElmtSize*(rEpVolInc/dtimeCur)
                      IF (rEpVolInc.LT.0.d0) THEN
                        rBVpress = rBVpress
     $                    -rRho*(rbv2*rElmtSize*rEpVolInc/dtimeCur)**2
                      END IF
                    ELSE
                      rBVpress = 0.d0
                    END IF

C                   d) Compute the internal force vector contribution
                    DO n = 1, nnode
                      DO i = 1, mcrd
                        DO j = 1, 6
                          rBmtN(j,(n-1)*mcrd+i) = rBmat(j,(n-1)*mcrd+i)
                          IF (j.LE.3) THEN
                            rTfactor  = 1.d0
                          ELSE
                            rTfactor  = 2.d0
                          END IF
                          DO k = 1, 3
                            rBmtN(j,(n-1)*mcrd+i)= rBmtN(j,(n-1)*mcrd+i)
     $                               - rAlphaN
     $                            * ( rTtens(k) *rBmat(k,(n-1)*mcrd+i)
     $                            +   rTtensk(k)*rBskw(k,(n-1)*mcrd+i) )
     $                               * rTfactor*rTtens(j)
                          END DO ! k
                          DO k = 4, 6
                            rBmtN(j,(n-1)*mcrd+i)= rBmtN(j,(n-1)*mcrd+i)
     $                               - rAlphaN
     $                               * rTtens(k)*rBmat(k,(n-1)*mcrd+i)
     $                               * rTfactor*rTtens(j)
                          END DO ! k
                        END DO ! j = 1, 6
                      END DO ! i = 1, mcrd
                    END DO ! n = 1, nnode

                    DO n = 1, nnode
                      DO i = 1, mcrd
                        DO j = 1, 6
                          IF (j.LE.3) THEN    ! add bulk viscosity 
                            rSigma = rStress(j) - rBVpress
                          ELSE
                            rSigma = rStress(j)
                          END IF
                          rhs( kblock , (n-1)*mcrd+i ) =
     $                        rhs( kblock , (n-1)*mcrd+i )
     $                      + (1.d0 - SB_vol/rElmtVol*Gt_temp)
     $                      * rBmtN(j,(n-1)*mcrd+i) * rSigma * rdOmega
                        END DO
                      END DO
                    END DO

                    EXIT            ! EXIT from the iteration loop

C                 If we're still iterating, we need to:
                  ELSE

C                   Calculate residual for the traction continuity equations
                    rG(1) = rG(1) + rNvect(1)*rStress(1)
     $                            + rNvect(2)*rStress(6)
     $                            + rNvect(3)*rStress(5)
                    rG(2) = rG(2) + rNvect(1)*rStress(6)
     $                            + rNvect(2)*rStress(2)
     $                            + rNvect(3)*rStress(4)
                    rG(3) = rG(3) + rNvect(1)*rStress(5)
     $                            + rNvect(2)*rStress(4)
     $                            + rNvect(3)*rStress(3)

c                   Update residual norm
                    rGnrm = sqrt( rG(1)**2 + rG(2)**2 + rG(3)**2 )

C                   Save initial residual norm
                    IF ((iCount.EQ.1).AND.(iAltern.EQ.0)) THEN
                      rGnrmZ = rGnrm
                    END IF

C                   Relative residual norm
                    rGrel = rGnrm/rGnrmZ

!                   DEBUGGING...
                    IF (iDbgFlg.EQ.1) THEN
                      iOutPr(1) = jelem(kblock)
                      iOutPr(2) = l
                      iOutPr(3) = kinc
                      iOutPr(4) = iAltern
                      iOutPr(5) = iCount

                      rOutPr(1) = rNvect(1)
                      rOutPr(2) = rNvect(2)
                      rOutPr(3) = rNvect(3)

                      rOutPr(4) = rMvect(1)
                      rOutPr(5) = rMvect(2)
                      rOutPr(6) = rMvect(3)

                      rOutPr(7) = rMHat(1)
                      rOutPr(8) = rMHat(2)
                      rOutPr(9) = rMHat(3)

                      rOutPr(10) = rAlphaHat
                      rOutPr(11) = rLnrm
                      rOutPr(12) = rGnrm
                      rOutPr(13) = rGrel

                      CALL XPLB_ABQERR(1,'Elem %I, GP %I, Inc. %I: '//
     $               'Initial Guess %I, iteration %I '//
     $               'N=(%R,%R,%R) '//
     $               'M=(%R,%R,%R) '//
     $               'M_hat=(%R,%R,%R) '//
     $               'Alpha_hat=%R, L_norm = %R '//
     $               'AbsResid=%R, RelResid=%R',iOutPr,rOutPr,' ')

                      DO i = 1, 15
                        iOutPr(i) = 0
                        rOutPr(i) = 0.d0
                      END DO

                    END IF ! iDbgFlg = 1
!                   ...GNIGGUBED

c                   Check for convergence
                    IF ((rGnrm.LE.rIterTol).OR.(rGrel.LT.rRelTol)) THEN

                      IF (iCvFlag.GT.0) THEN     ! No violation of yield condition
!                       DEBUGGING...
                        IF (iDbgFlg.EQ.1) THEN
                          CALL XPLB_ABQERR(1,'CONVERGED: RESIDUAL',
     $                                     iOutPr,rOutPr,' ')
                        END IF
!                       ...GNIGGUBED

                        iConverge = 1   ! we have a converged solution
                        CYCLE           ! skip to the final pass
                      ELSE
                        iConverge = 0  !  This initial guess has failed
                        EXIT           !  Try a different initial guess
                      END IF

                    END IF ! residual/convergence check

C                   Save minimum residual norm and corresponding parameter values
                    IF ( ( (iCount.EQ.1).AND.(iAltern.LE.1) ).OR.
     $                   ( (iCount.GT.1).AND.(rGnrm.LT.rGmin(1)) ).OR.
     $                   ( (iAltern.GE.1).AND.(rGnrm.LT.rGmin(1)) ) ) THEN
                      rGmin(1) = rGnrm
                      rGmin(2) = rAlphaHat
                      rGmin(3) = rMvect(1)
                      rGmin(4) = rMvect(2)
                      rGmin(5) = rMvect(3)
                      iHaveSol = 1

!                     Debugging...
                      IF (iDbgFlg.EQ.1) THEN
                        rOutPr(1) = rMHat(1)
                        rOutPr(2) = rMHat(2)
                        rOutPr(3) = rMHat(3)
                        rOutPr(4) = rAlphaHat
                        rOutPr(5) = rGnrm
                        CALL XPLB_ABQERR(1,'Saving minimum residual: '//
     $                                     'M=(%R,%R,%R) '//
     $                                     'Alpha_hat=%R '//
     $                                     'Resid=%R',iOutPr,rOutPr,' ')
                        DO i = 1, 15
                          iOutPr(i) = 0
                          rOutPr(i) = 0.d0
                        END DO
                      END IF ! iDbgFlg
!                     ...Gniggubed
                    END IF ! save minimum residual

C                   Calculate tangent matrix for the traction continuity equations

                    rnCn(1,1) =   rnCn(1,1)
     $                          + rNvect(1) * rCmod(1,1) * rNvect(1)
     $                          + rNvect(1) * rCmod(1,6) * rNvect(2)
     $                          + rNvect(1) * rCmod(1,5) * rNvect(3)
     $                          + rNvect(2) * rCmod(6,1) * rNvect(1)
     $                          + rNvect(2) * rCmod(6,6) * rNvect(2)
     $                          + rNvect(2) * rCmod(6,5) * rNvect(3)
     $                          + rNvect(3) * rCmod(5,1) * rNvect(1)
     $                          + rNvect(3) * rCmod(5,6) * rNvect(2)
     $                          + rNvect(3) * rCmod(5,5) * rNvect(3)
     $                          - rNvect(1) * rGtan(1,3) * rNvect(2) * 0.5d0
     $                          + rNvect(1) * rGtan(1,2) * rNvect(3) * 0.5d0
     $                          - rNvect(2) * rGtan(6,3) * rNvect(2) * 0.5d0
     $                          + rNvect(2) * rGtan(6,2) * rNvect(3) * 0.5d0
     $                          - rNvect(3) * rGtan(5,3) * rNvect(2) * 0.5d0
     $                          + rNvect(3) * rGtan(5,2) * rNvect(3) * 0.5d0

                    rnCn(2,2) =   rnCn(2,2)
     $                          + rNvect(1) * rCmod(6,6) * rNvect(1)
     $                          + rNvect(1) * rCmod(6,2) * rNvect(2)
     $                          + rNvect(1) * rCmod(6,4) * rNvect(3)
     $                          + rNvect(2) * rCmod(2,6) * rNvect(1)
     $                          + rNvect(2) * rCmod(2,2) * rNvect(2)
     $                          + rNvect(2) * rCmod(2,4) * rNvect(3)
     $                          + rNvect(3) * rCmod(4,6) * rNvect(1)
     $                          + rNvect(3) * rCmod(4,2) * rNvect(2)
     $                          + rNvect(3) * rCmod(4,4) * rNvect(3)
     $                          + rNvect(1) * rGtan(6,3) * rNvect(1) * 0.5d0
     $                          - rNvect(1) * rGtan(6,1) * rNvect(3) * 0.5d0
     $                          + rNvect(2) * rGtan(2,3) * rNvect(1) * 0.5d0
     $                          - rNvect(2) * rGtan(2,1) * rNvect(3) * 0.5d0
     $                          + rNvect(3) * rGtan(4,3) * rNvect(1) * 0.5d0
     $                          - rNvect(3) * rGtan(4,1) * rNvect(3) * 0.5d0

                    rnCn(3,3) =   rnCn(3,3)
     $                          + rNvect(1) * rCmod(5,5) * rNvect(1)
     $                          + rNvect(1) * rCmod(5,4) * rNvect(2)
     $                          + rNvect(1) * rCmod(5,3) * rNvect(3)
     $                          + rNvect(2) * rCmod(4,5) * rNvect(1)
     $                          + rNvect(2) * rCmod(4,4) * rNvect(2)
     $                          + rNvect(2) * rCmod(4,3) * rNvect(3)
     $                          + rNvect(3) * rCmod(3,5) * rNvect(1)
     $                          + rNvect(3) * rCmod(3,4) * rNvect(2)
     $                          + rNvect(3) * rCmod(3,3) * rNvect(3)
     $                          - rNvect(1) * rGtan(5,2) * rNvect(1) * 0.5d0
     $                          + rNvect(1) * rGtan(5,1) * rNvect(2) * 0.5d0
     $                          - rNvect(2) * rGtan(4,2) * rNvect(1) * 0.5d0
     $                          + rNvect(2) * rGtan(4,1) * rNvect(2) * 0.5d0
     $                          - rNvect(3) * rGtan(3,2) * rNvect(1) * 0.5d0
     $                          + rNvect(3) * rGtan(3,1) * rNvect(2) * 0.5d0

                    rnCn(1,2) =   rnCn(1,2)
     $                          + rNvect(1) * rCmod(1,6) * rNvect(1)
     $                          + rNvect(1) * rCmod(1,2) * rNvect(2)
     $                          + rNvect(1) * rCmod(1,4) * rNvect(3)
     $                          + rNvect(2) * rCmod(6,6) * rNvect(1)
     $                          + rNvect(2) * rCmod(6,2) * rNvect(2)
     $                          + rNvect(2) * rCmod(6,4) * rNvect(3)
     $                          + rNvect(3) * rCmod(5,6) * rNvect(1)
     $                          + rNvect(3) * rCmod(5,2) * rNvect(2)
     $                          + rNvect(3) * rCmod(5,4) * rNvect(3)
     $                          + rNvect(1) * rGtan(1,3) * rNvect(1) * 0.5d0
     $                          - rNvect(1) * rGtan(1,1) * rNvect(3) * 0.5d0
     $                          + rNvect(2) * rGtan(6,3) * rNvect(1) * 0.5d0
     $                          - rNvect(2) * rGtan(6,1) * rNvect(3) * 0.5d0
     $                          + rNvect(3) * rGtan(5,3) * rNvect(1) * 0.5d0
     $                          - rNvect(3) * rGtan(5,1) * rNvect(3) * 0.5d0

                    rnCn(2,1) =   rnCn(2,1)
     $                          + rNvect(1) * rCmod(6,1) * rNvect(1)
     $                          + rNvect(1) * rCmod(6,6) * rNvect(2)
     $                          + rNvect(1) * rCmod(6,5) * rNvect(3)
     $                          + rNvect(2) * rCmod(2,1) * rNvect(1)
     $                          + rNvect(2) * rCmod(2,6) * rNvect(2)
     $                          + rNvect(2) * rCmod(2,5) * rNvect(3)
     $                          + rNvect(3) * rCmod(4,1) * rNvect(1)
     $                          + rNvect(3) * rCmod(4,6) * rNvect(2)
     $                          + rNvect(3) * rCmod(4,5) * rNvect(3)
     $                          - rNvect(1) * rGtan(6,3) * rNvect(2) * 0.5d0
     $                          + rNvect(1) * rGtan(6,2) * rNvect(3) * 0.5d0
     $                          - rNvect(2) * rGtan(2,3) * rNvect(2) * 0.5d0
     $                          + rNvect(2) * rGtan(2,2) * rNvect(3) * 0.5d0
     $                          - rNvect(3) * rGtan(4,3) * rNvect(2) * 0.5d0
     $                          + rNvect(3) * rGtan(4,2) * rNvect(3) * 0.5d0

                    rnCn(1,3) =   rnCn(1,3)
     $                          + rNvect(1) * rCmod(1,5) * rNvect(1)
     $                          + rNvect(1) * rCmod(1,4) * rNvect(2)
     $                          + rNvect(1) * rCmod(1,3) * rNvect(3)
     $                          + rNvect(2) * rCmod(6,5) * rNvect(1)
     $                          + rNvect(2) * rCmod(6,4) * rNvect(2)
     $                          + rNvect(2) * rCmod(6,3) * rNvect(3)
     $                          + rNvect(3) * rCmod(5,5) * rNvect(1)
     $                          + rNvect(3) * rCmod(5,4) * rNvect(2)
     $                          + rNvect(3) * rCmod(5,3) * rNvect(3)
     $                          - rNvect(1) * rGtan(1,2) * rNvect(1) * 0.5d0
     $                          + rNvect(1) * rGtan(1,1) * rNvect(2) * 0.5d0
     $                          - rNvect(2) * rGtan(6,2) * rNvect(1) * 0.5d0
     $                          + rNvect(2) * rGtan(6,1) * rNvect(2) * 0.5d0
     $                          - rNvect(3) * rGtan(5,2) * rNvect(1) * 0.5d0
     $                          + rNvect(3) * rGtan(5,1) * rNvect(2) * 0.5d0

                    rnCn(3,1) =   rnCn(3,1)
     $                          + rNvect(1) * rCmod(5,1) * rNvect(1)
     $                          + rNvect(1) * rCmod(5,6) * rNvect(2)
     $                          + rNvect(1) * rCmod(5,5) * rNvect(3)
     $                          + rNvect(2) * rCmod(4,1) * rNvect(1)
     $                          + rNvect(2) * rCmod(4,6) * rNvect(2)
     $                          + rNvect(2) * rCmod(4,5) * rNvect(3)
     $                          + rNvect(3) * rCmod(3,1) * rNvect(1)
     $                          + rNvect(3) * rCmod(3,6) * rNvect(2)
     $                          + rNvect(3) * rCmod(3,5) * rNvect(3)
     $                          - rNvect(1) * rGtan(5,3) * rNvect(2) * 0.5d0
     $                          + rNvect(1) * rGtan(5,2) * rNvect(3) * 0.5d0
     $                          - rNvect(2) * rGtan(4,3) * rNvect(2) * 0.5d0
     $                          + rNvect(2) * rGtan(4,2) * rNvect(3) * 0.5d0
     $                          - rNvect(3) * rGtan(3,3) * rNvect(2) * 0.5d0
     $                          + rNvect(3) * rGtan(3,2) * rNvect(3) * 0.5d0

                    rnCn(2,3) =   rnCn(2,3)
     $                          + rNvect(1) * rCmod(6,5) * rNvect(1)
     $                          + rNvect(1) * rCmod(6,4) * rNvect(2)
     $                          + rNvect(1) * rCmod(6,3) * rNvect(3)
     $                          + rNvect(2) * rCmod(2,5) * rNvect(1)
     $                          + rNvect(2) * rCmod(2,4) * rNvect(2)
     $                          + rNvect(2) * rCmod(2,3) * rNvect(3)
     $                          + rNvect(3) * rCmod(4,5) * rNvect(1)
     $                          + rNvect(3) * rCmod(4,4) * rNvect(2)
     $                          + rNvect(3) * rCmod(4,3) * rNvect(3)
     $                          - rNvect(1) * rGtan(6,2) * rNvect(1) * 0.5d0
     $                          + rNvect(1) * rGtan(6,1) * rNvect(2) * 0.5d0
     $                          - rNvect(2) * rGtan(2,2) * rNvect(1) * 0.5d0
     $                          + rNvect(2) * rGtan(2,1) * rNvect(2) * 0.5d0
     $                          - rNvect(3) * rGtan(4,2) * rNvect(1) * 0.5d0
     $                          + rNvect(3) * rGtan(4,1) * rNvect(2) * 0.5d0

                    rnCn(3,2) =   rnCn(3,2)
     $                          + rNvect(1) * rCmod(5,6) * rNvect(1)
     $                          + rNvect(1) * rCmod(5,2) * rNvect(2)
     $                          + rNvect(1) * rCmod(5,4) * rNvect(3)
     $                          + rNvect(2) * rCmod(4,6) * rNvect(1)
     $                          + rNvect(2) * rCmod(4,2) * rNvect(2)
     $                          + rNvect(2) * rCmod(4,4) * rNvect(3)
     $                          + rNvect(3) * rCmod(3,6) * rNvect(1)
     $                          + rNvect(3) * rCmod(3,2) * rNvect(2)
     $                          + rNvect(3) * rCmod(3,4) * rNvect(3)
     $                          + rNvect(1) * rGtan(5,3) * rNvect(1) * 0.5d0
     $                          - rNvect(1) * rGtan(5,1) * rNvect(3) * 0.5d0
     $                          + rNvect(2) * rGtan(4,3) * rNvect(1) * 0.5d0
     $                          - rNvect(2) * rGtan(4,1) * rNvect(3) * 0.5d0
     $                          + rNvect(3) * rGtan(3,3) * rNvect(1) * 0.5d0
     $                          - rNvect(3) * rGtan(3,1) * rNvect(3) * 0.5d0

                  END IF   !  still iterating

C                 Enforce traction continuity to update alpha_N and rMvect
C                 --------------------------------------------------------

C                 If we run out of iterations, and still haven't converged
                  IF (iCount.EQ.(iIterMax-1)) THEN
                    EXIT                       ! try another initial guess

C                 If we're not out of iterations, lets ITERATE
                  ELSE
                    CALL MInvert(rnCn,3,3,iErrFlg)
                    IF (iErrFlg.NE.0) THEN     ! Problem inverting  matrix
                      iOutPr(1) = jelem(kblock)
                      iOutPr(2) = l
                      iOutPr(3) = kinc
                      rOutPr(1) = rnCn(1,1)
                      rOutPr(2) = rnCn(1,2)
                      rOutPr(3) = rnCn(1,3)
                      rOutPr(4) = rnCn(2,1)
                      rOutPr(5) = rnCn(2,2)
                      rOutPr(6) = rnCn(2,3)
                      rOutPr(7) = rnCn(3,1)
                      rOutPr(8) = rnCn(3,2)
                      rOutPr(9) = rnCn(3,3)
                      CALL XPLB_ABQERR(-2,'Problem inverting matrix: '//
     $               'in element %I, Gauss point %I, at increment %I,'//
     $               'A(1,1) = %R'//
     $               'A(1,2) = %R'//
     $               'A(1,3) = %R'//
     $               'A(2,1) = %R'//
     $               'A(2,2) = %R'//
     $               'A(2,3) = %R'//
     $               'A(3,1) = %R'//
     $               'A(3,2) = %R'//
     $               'A(3,3) = %R',iOutPr,rOutPr,' ')
                      DO i = 1, 15
                        iOutPr(i) = 0
                        rOutPr(i) = 0.d0
                      END DO
                      EXIT                   ! try another initial guess
                    END IF ! problem inverting matrix rnCn

C                   Compute solution increment
                    DO i = 1, 3
                      rDelta(i) = 0.d0
                      DO j = 1, 3
                        rDelta(i) = rDelta(i) + rnCn(i,j)*rG(j)
                      END DO
                    END DO

                    rDelNrm = sqrt( rDelta(1)**2 + rDelta(2)**2 
     $                            + rDelta(3)**2)

c                   Check for convergence
                    IF (     (rDelNrm.LE.5.d-16)
     $                  .AND.(rGnrm.LE.sqrt(rIterTol)) ) THEN

                      IF (iCvFlag.GT.0) THEN     ! No violation of yield condition
!                       DEBUGGING...
                        IF (iDbgFlg.EQ.1) THEN
                          write(*,*) "CONVERGED: SMALL CORRECTION"
                          write(*,*) "Delta = ", ( rDelta(i), i=1,3 )
                        END IF
!                       ...GNIGGUBED

                        iConverge = 1   ! we have a converged solution
                        CYCLE           ! skip to the final pass
                      ELSE
                        iConverge = 0  !  This initial guess has failed
                        EXIT           !  Try a different initial guess
                      END IF

                    END IF ! check for convergence

C                   Update solution variables

C                   First, we compute the under-relaxation factor
                    rCfacR = 1.d0

                    rCFacC = 1.d0
  
                    rCFactor = min(rCFacR, rCFacC)
 
                    DO i = 1, 3
                      rMhat(i)  =  rMhat(i) + rDelta(i) * rCFactor
                    END DO

                    rAlphaHat = sqrt( rMhat(1)**2 + rMhat(2)**2
     $                              + rMhat(3)**2 )

c                   If the discontinuous strain increment is appreciable (non-zero)
                    IF (rAlphaHat.GT.1.d-16) THEN

C                     Update the M vector 
                      DO i = 1, 3
                        rMvect(i) = rMhat(i)/rAlphaHat
                      END DO

C                     And update the T tensor
                      rTtens(1) =  rMvect(1)*rNvect(1)
                      rTtens(2) =  rMvect(2)*rNvect(2)
                      rTtens(3) =  rMvect(3)*rNvect(3)   ! T_33
                      rTtens(4) = 0.5d0*( rMvect(2)*rNvect(3)
     $                                  + rMvect(3)*rNvect(2) )   ! 1/2*(T_23+T_32)
                      rTtens(5) = 0.5d0*( rMvect(1)*rNvect(3)
     $                                  + rMvect(3)*rNvect(1) )   ! 1/2*(T_13+T_31)
                      rTtens(6) = 0.5d0*( rMvect(1)*rNvect(2)
     $                                  + rMvect(2)*rNvect(1) )   ! 1/2*(T_12+T_21)

                      rTtensk(1)= 0.5d0*( rMvect(3)*rNvect(2)
     $                                  - rMvect(2)*rNvect(3) )   ! 1/2*(T_32-T_23)
                      rTtensk(2)= 0.5d0*( rMvect(1)*rNvect(3)
     $                                  - rMvect(3)*rNvect(1) )   ! 1/2*(T_13-T_31)
                      rTtensk(3)= 0.5d0*( rMvect(2)*rNvect(1)
     $                                  - rMvect(1)*rNvect(2) )   ! 1/2*(T_21-T_12)

                    END IF ! discontinuous strain increment is appreciable

                  END IF ! not last iteration

                END DO  ! iCount

                IF (iConverge.EQ.1) THEN

                  svars(kblock, iTIntPt*(nstatev(1)+nstatev(2))
     $                         +  (l-1)* nstatev(3) + 5)   = rAlphaN
                  svars(kblock, iTIntPt*(nstatev(1)+nstatev(2))
     $                         +  (l-1)* nstatev(3) + 6) = rMvect(1)
                  svars(kblock, iTIntPt*(nstatev(1)+nstatev(2))
     $                         +  (l-1)* nstatev(3) + 7) = rMvect(2)
                  svars(kblock, iTIntPt*(nstatev(1)+nstatev(2))
     $                         +  (l-1)* nstatev(3) + 8) = rMvect(3)

                ELSEIF (iAltern.LT.iAltMax) THEN     ! try alternative initial guess

                  iAltern = iAltern + 1

                  rAlphaN   = svars(kblock,
     $                                iTIntPt*(nstatev(1)+nstatev(2))
     $                               +  (l-1)* nstatev(3) + 5)
                  rMvect(1) = svars(kblock,
     $                                iTIntPt*(nstatev(1)+nstatev(2))
     $                               +  (l-1)* nstatev(3) + 6)
                  rMvect(2) = svars(kblock,
     $                                iTIntPt*(nstatev(1)+nstatev(2))
     $                               +  (l-1)* nstatev(3) + 7)
                  rMvect(3) = svars(kblock,
     $                                iTIntPt*(nstatev(1)+nstatev(2))
     $                               +  (l-1)* nstatev(3) + 8)

                  CALL Galtinit_3D(iAltern,iAltMax,iDbgFlg,rNvect,
     $                          iHaveSol,rGmin,jElem,kblock,nblock,
     $                          l,kinc,rIterTol,time,nTime,
     $                          rGnrm,rAlphaHat,rAlphaN,rMvect)

                  GOTO 10

                END IF ! iAltern.LT.iAltMax

              ElSE
                CALL XPLB_ABQERR(-3,'Localization flag error',
     &                         0,0.d0,' ')

              END IF  ! iLocalFlg = 0 or 1

            END DO ! Numerical integration

c>>>        Compute stable time step
C
            dtimeStable(kblock) = rMinElmtL / rMaxWavSpd

            dtimeStable(kblock) = dtimeStable(kblock) * rStFctr

C           write state variables every N time steps
            IF (kinc.GT.0) THEN

              DO l = 1 , iTIntPt
C               matrix mateiral
                iTemp = (l-1)*nstatev(1)
                DO i = 1, 6
                  rStress(i) = svars(kblock,iTemp + i )
                END DO

C               Von Mises stress
                uel_isvs(l,1) = 
     $                   sqrt( 1.5d0*(  rStress(1)**2 + rStress(2)**2
     $                   + rStress(3)**2 + 2.d0*rStress(4)**2
     $                   + 2.d0*rStress(5)**2 + 2.d0*rStress(6)**2 )
     $                  -0.5d0*( rStress(1)+rStress(2)+rStress(3) )**2 )
C               Equiv. plastic strain  
                uel_isvs(l,2) = svars(kblock,iTemp + 7 )
C               Temperature
                uel_isvs(l,3) = svars(kblock,iTemp + 8 )
C               Evolving MTS
                uel_isvs(l,4) = svars(kblock,iTemp + 9 )

                ! Localization flag
                iTemp = iTIntPt*nstatev(1) + iTIntPt*nstatev(2)
     $                +   (l-1)*nstatev(3)
                uel_isvs(l,10) = svars(kblock, iTemp + 1)


C               Since all EIGHT Gauss points have localization at the same time
                IF (nint( uel_isvs(l,10) ).EQ.1) THEN    ! Localized

                  ! Band material
                  iTemp = iTIntPt*nstatev(1) + (l-1)*nstatev(2)
                  DO i = 1, 6
                    rStress(i) = svars(kblock, iTemp + i)
                  END DO
C                 Von Mises Stress
                  uel_isvs(l,5) =
     $                 sqrt( 1.5d0*(  rStress(1)**2 + rStress(2)**2
     $                     + rStress(3)**2 + 2.d0*rStress(4)**2
     $                     + 2.d0*rStress(5)**2 + 2.d0*rStress(6)**2 )
     $             -0.5d0*( rStress(1)+rStress(2)+rStress(3) )**2 )
C                 Equiv. plastic strain  
                  uel_isvs(l,6) = svars(kblock,iTemp + 7 )
C                 Temperature
                  uel_isvs(l,7) = svars(kblock,iTemp + 8 )
C                 Evolving MTS
                  uel_isvs(l,8) = svars(kblock,iTemp + 9 )

                ELSE                                    ! Not localized (yet)
c                 Von Mises stress                  
                  uel_isvs(l,5) = -1.d0
C                 Equiv. plastic strain
                  uel_isvs(l,6) = -1.d0
C                 Temperature
                  uel_isvs(l,7) = -1.d0
C                 Evolving MTS
                  uel_isvs(l,8) = -1.d0

                END IF

              END DO ! l

C             Average quantities for the elmt
              DO i = 1,4
                elmt_SDV(i,elmt_number) = 
     $            ( uel_isvs(1,i) + uel_isvs(2,i)
     $           +  uel_isvs(3,i) + uel_isvs(4,i)
     $           +  uel_isvs(5,i) + uel_isvs(6,i)
     $           +  uel_isvs(7,i) + uel_isvs(8,i))/8.d0
              END DO

              elmt_SDV(10,elmt_number) =
     $                       uel_isvs(1,10) + uel_isvs(2,10) +
     $                       uel_isvs(3,10) + uel_isvs(4,10) +
     $                       uel_isvs(5,10) + uel_isvs(6,10) +
     $                       uel_isvs(7,10) + uel_isvs(8,10)

              IF (elmt_SDV(10,elmt_number).GT.0.d0) THEN
                DO i = 5,8
                  elmt_SDV(i,elmt_number) =
     $                       ( uel_isvs(1,i)*uel_isvs(1,10) +
     $                         uel_isvs(2,i)*uel_isvs(2,10) +
     $                         uel_isvs(3,i)*uel_isvs(3,10) +
     $                         uel_isvs(4,i)*uel_isvs(4,10) +
     $                         uel_isvs(5,i)*uel_isvs(5,10) +
     $                         uel_isvs(6,i)*uel_isvs(6,10) +
     $                         uel_isvs(7,i)*uel_isvs(7,10) +
     $                         uel_isvs(8,i)*uel_isvs(8,10))/
     $                         elmt_SDV(10,elmt_number)
                END DO
              ELSE
                DO i = 5,8
                 elmt_SDV(i,elmt_number) = uel_isvs(1,i)
                END DO
              END IF

              DO l = 1, 8
                IF ( uel_isvs(l,10).GT.0.d0 ) THEN
                  elmt_SDV(10, elmt_number) = 1.d0
                END IF
              END DO 

            END IF  ! kinc.GT.0

          END DO  ! kblock (to define rhs, dtimeStable, svars, and energy)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC        Define external force from distributed loads (optional)       CCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        ELSE IF ( lflags(iOpCode) .EQ. jExternForce ) THEN
C       Define external force from distributed loads

          DO kblock = 1,nblock
C           code to define rhs
          END DO

        END IF ! lflags(iOpCode)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC        ERROR: Specified element type or analysis type incorrect      CCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ELSE   ! Write error message and exit
        CALL XPLB_ABQERR(-3,'User element VU1 usable only in '//
     $                      'Explicit Dynamic Analysis',0,0.d0,' ')
      END IF

      RETURN
      END  ! Subroutine VUEL

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)

C      include 'vaba_param.inc'
      IMPLICIT NONE
#INCLUDE <SMAAspUserSubroutines.hdr> ! Start this line WITHOUT ANY SPACE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Variables passed from interface

      integer lOp, niArray, nrArray

      integer i_Array(niArray)

      real*8 r_Array(nrArray)

C     Contents of i_Array
      integer i_int_nTotalNodes, i_int_nTotalElements,
     1        i_int_kStep, i_int_kInc, i_int_iStatus,
     2        i_int_lWriteRestart
      parameter( i_int_nTotalNodes     = 1,
     *           i_int_nTotalElements  = 2,
     *           i_int_kStep           = 3,
     *           i_int_kInc            = 4,
     *           i_int_iStatus         = 5,
     *           i_int_lWriteRestart   = 6  )

C     Possible values for the lOp argument
      integer j_int_StartAnalysis, j_int_StartStep,
     1        j_int_SetupIncrement, j_int_StartIncrement,
     2        j_int_EndIncrement, j_int_EndStep,
     3        j_int_EndAnalysis
      parameter( j_int_StartAnalysis    = 0,
     *           j_int_StartStep        = 1,
     *           j_int_SetupIncrement   = 2,
     *           j_int_StartIncrement   = 3,
     *           j_int_EndIncrement     = 4,
     *           j_int_EndStep          = 5,
     *           j_int_EndAnalysis      = 6 )

C     Possible values for i_Array(i_int_iStatus)
      integer j_int_Continue, j_int_TerminateStep,
     1        j_int_TerminateAnalysis
      parameter( j_int_Continue          = 0,
     *           j_int_TerminateStep     = 1,
     *           j_int_TerminateAnalysis = 2)

C     Contents of r_Array
      integer i_flt_TotalTime, i_flt_StepTime,
     1        i_flt_dTime
      parameter( i_flt_TotalTime   = 1,
     1           i_flt_StepTime    = 2,
     2           i_flt_dTime       = 3 )


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Variables defined locally

      integer kStep, kInc

C     elmt connectivity array for 8-node brick elmt
      integer elmt_conct(8,i_Array(i_int_nTotalElements))
      pointer(ptr_elmt_conct, elmt_conct)

C     elmt_neighbour array, each elmt has 6 neighbouring elmts at most
      integer elmt_neighbour(6,i_Array(i_int_nTotalElements))
      pointer(ptr_elmt_neighbour, elmt_neighbour)

C     elmt shear status (1) and front status (2)
      integer elmt_status(2,i_Array(i_int_nTotalElements))
      pointer(ptr_elmt_status, elmt_status)

C     elmt SDV values
      real*8 elmt_SDV(10,i_Array(i_int_nTotalElements))
      pointer(ptr_elmt_SDV, elmt_SDV)

C     node_theta, the nodal theta value
      real*8 node_theta(i_Array(i_int_nTotalNodes))
      pointer(ptr_node_theta, node_theta)

C     nodal coordinates (deformed configuration)
      real*8 node_coords(3,i_Array(i_int_nTotalNodes))
      pointer(ptr_node_coords, node_coords)

C     Relavtive length on each edge of the elmt containing a shear band
C     One elmt has 12 edges
      real*8 elmt_edge_rlength(12,i_Array(i_int_nTotalElements))
      pointer(ptr_elmt_edge_rlength, elmt_edge_rlength)

C     Shear band surface normal direction at each elmt's centroid
      real*8 elmt_Nvect(3,i_Array(i_int_nTotalElements))
      pointer(ptr_elmt_Nvect, elmt_Nvect)

C     Number of ghost elmts
      integer num_ghost_elmt
      pointer(ptr_num_ghost_elmt, num_ghost_elmt)

C     Whether a node is applied with essential B.C.
      integer BC_index(i_Array(i_int_nTotalNodes))
      pointer(ptr_BC_index, BC_index)

      integer i,ios, j,ghost_elmt_counter

      integer myThreadID, elmt_list(6), num_VUEL_elmt

      logical file_exist

      integer init_elmt, BC_node

      real*8  BC_temp

      character*20 str_kinc

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Change the directory accordingly
      character*255 cwd
      parameter(cwd='/home/pedgaonkar/abaqus_projects/NSF_Twin_Code/'//
     $              'uw_twin_code/')

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Array ID = 1, int, elmt connnection matrix;
C           ID = 2, int, elmt neighbouring elmts matrix;  
C           ID = 3, int, elmt shear and front status, initial value = 0; 
C           ID = 4, real*8, nodal theta matrix, initial value = 0.0;
C           ID = 5, real*8, relative length, initial value = -10.0;
C           ID = 7, real*8, nodal coordinates (deformed configuration)
C           ID = 8, real*8, Shear band surface normal direciton at each elmt
C                           centroid
C           ID = 9, int, whether a node is applied with essential B.C.
C                        1 - yes; 0 - no;    
C           ID = 10, number of ghost elmts in the simulation
C           ID = 11, real*8, SDV variables at the centroid of each elmt

      kStep = i_Array(i_int_kStep)
      kInc  = i_Array(i_int_kInc)


C     Start of the analysis
      IF (lOp .EQ. j_int_StartAnalysis) THEN
        print *, 'Total number of threads = ', GETNUMTHREADS()
        print *, 'Start of the analysis'

        myThreadID = get_thread_id()

        call MutexInit( 1 )      ! initialize Mutex #1


        IF (myThreadID.EQ.0) THEN

          ! create a GLOBAL array with ID=1 and SIZE = 8 * Number of elmts
          ptr_elmt_conct = SMAIntArrayCreate(1, 
     1                     8*i_Array(i_int_nTotalElements), 0)

          ! create a GLOBAL array with ID=2 and SIZE = 6 * Number of elmts
          ptr_elmt_neighbour = SMAIntArrayCreate(2, 
     1                         6*i_Array(i_int_nTotalElements), 0)

          ! create a GLOBAL array with ID=3 and SIZE = 2 * Number of elmts
          ptr_elmt_status = SMAIntArrayCreate(3,
     1                         2*i_Array(i_int_nTotalElements), 0)

          ! create a GLOBAL array with ID=4 and SIZE = 1 * Number of total nodes
          ptr_node_theta = SMAFloatArrayCreate(4,
     1                         i_Array(i_int_nTotalNodes), 0.d0)

          ! create a GLOBAL array with ID=5 and size = 12 * Number of total elmts
          ptr_elmt_edge_rlength = SMAFloatArrayCreate(5,
     1                        12*i_Array(i_int_nTotalElements), -10.d0)
         
          ! create a GLOBAL array with ID=7 and size = 3 * Number of total nodes
          ptr_node_coords = SMAFloatArrayCreate(7,
     1                       3*i_Array(i_int_nTotalNodes), 0.d0)

          ! create a GLOBAL array with ID=8 and size = 3 * Number of total elmts
          ptr_elmt_Nvect = SMAFloatArrayCreate(8,
     1                       3*i_Array(i_int_nTotalElements), 0.d0)

          ! create a GLOBAL array with ID=9 and size = 1 * Number of total nodes
          ptr_BC_index = SMAIntArrayCreate(9,
     1                         i_Array(i_int_nTotalNodes), 0)

          ! create a GLOBAL array with ID=10 and size = 1
          ptr_num_ghost_elmt = SMAIntArrayCreate(10,1,0)


          ! create a GLOBAL array with ID=11 and SIZE = 10 * Number of elmts
          ptr_elmt_SDV = SMAFloatArrayCreate(11,
     1                         10*i_Array(i_int_nTotalElements), 0.d0)

C         abaqus/Explicit, specify units 16–18 or unit numbers 
C         greater than 100 ending in 5 to 9, e.g. 105, 268, etc. 

C         Read the Input_ghost_elmt_selected.inp and get the number
C         of ghost elmts
          ghost_elmt_counter = 0
          open(unit = 105, 
     1      file = trim(cwd)//'Input_ghost_elmt_selected.inp',
     2      status = 'OLD')
          DO
            read(105,*,END=1003)
            ghost_elmt_counter = ghost_elmt_counter + 1
          END DO
1003      close(105)
          num_ghost_elmt = ghost_elmt_counter 

          num_VUEL_elmt = i_Array(i_int_nTotalElements) - num_ghost_elmt

C         CHANGE the address of the elmt connectivity matrix
          open(unit = 105, 
     1      file = trim(cwd)//'Input_elmt.inp',
     2      status = 'OLD')
          DO i = 1, num_VUEL_elmt
            read(105,*) j, elmt_conct(1:8, i)
C            print *, 'elmt', j
C            print *,  elmt_conct(1:8,i)
          END DO
          close(105)

C         Read the elmt number where the shear band initiates
          open(unit = 105, 
     1      file = trim(cwd)//'INIT_ELMT',
     2      status = 'OLD')
          DO
          read(105,*,END=1001) init_elmt
C          print *, init_elmt
          elmt_status(1, init_elmt) = 0
          elmt_status(2, init_elmt) = 1
          END DO
1001      close(105)  

C         Read the B.C. for the heat-like equation
          open(unit = 105, 
     1      file = trim(cwd)//'INIT_BC',
     2      status = 'OLD')
          DO
            read(105,*,END=1002) BC_node, BC_temp
            BC_index(BC_node) = 1
C          print *, BC_node, BC_temp
          node_theta(BC_node) = BC_temp
          END DO
1002      close(105) 


C         Form the elmt neighbour array
          DO i = 1, num_VUEL_elmt
            CALL Find_Neighbour_Elmts(i,elmt_conct(1:8,1:num_VUEL_elmt),
     1                                num_VUEL_elmt,  elmt_list)
            elmt_neighbour(1:6,i) = elmt_list
          END DO

        END IF  ! myThreadID = 0

C     Start of the step
      ELSE IF (lOp .EQ. j_int_StartStep) THEN
        print *, 'Start of the step'

C       Set up or exchange (import and export) initial values
C       with external programs.      
      
C       The initial values may need to match those at the 
C       point of restart.
        IF ( kInc .NE. 0) THEN
        END IF

C     Setup the increment       
      ELSE IF (lOp .EQ. j_int_SetupIncrement) THEN        
C       Change i_Array(i_int_lWriteRestart) and 
C              i_Array(i_int_iStatus) if desired.      
C       Change r_Array(i_flt_dTime) if desired.      

C        if (r_Array(i_flt_StepTime).GT.2.5d-5) THEN
C          r_Array(i_flt_dTime) = 1.0d-9
C        end if

C        if (r_Array(i_flt_StepTime).GT.3.0d-5) THEN
C          r_Array(i_flt_dTime) = 0.5d-9
C        end if

C     Start of the increment
      ELSE IF (lOp .EQ. j_int_StartIncrement) THEN
C       The time increment is finalized. Use r_Array(i_flt_dTime) if desired.
C       If needed, gather and export data from the configuration at the end of
C       the previous increment to external programs.      
C       Import and scatter data from external program to influence the current Abaqus increment.      

C     End of the increment
C     Solve the heat-like equation
      ELSE IF (lOp .EQ. j_int_EndIncrement) THEN
        ptr_elmt_conct = SMAIntArrayAccess(1)
        ptr_elmt_status = SMAIntArrayAccess(3)
        ptr_node_theta = SMAFloatArrayAccess(4)
        ptr_node_coords = SMAFloatArrayAccess(7)
        ptr_elmt_Nvect = SMAFloatArrayAccess(8)
        ptr_BC_index = SMAIntArrayAccess(9)
        ptr_num_ghost_elmt = SMAIntArrayAccess(10)
        ptr_elmt_SDV = SMAFloatArrayAccess(11)

        num_VUEL_elmt = i_Array(i_int_nTotalElements) - num_ghost_elmt

        myThreadID = get_thread_id()

C       update nodal theta every N time steps
        if (mod(kInc,400).eq.0) then

          call Mutexlock( 1 )

          if (myThreadID.eq.0) then
            print *, 'Start to solve heat eq. at incre', kinc
            call solve_heat_global(i_Array(i_int_nTotalNodes),
     1                         num_VUEL_elmt,
     2                         node_coords,
     3                         elmt_conct(1:8,1:num_VUEL_elmt),
     4                         node_theta,
     5                         elmt_Nvect(1:3,1:num_VUEL_elmt),
     6                         elmt_status(1:2,1:num_VUEL_elmt),
     7                         BC_index)
            print *, 'End to solve heat eq. at incre', kinc

          ! Write state variables into a txt file
          ! abaqus/Explicit, specify units 16–18 or unit numbers 
          ! greater than 100 ending in 5 to 9, e.g. 105, 268, etc.
            write(str_kinc,'(I0)') kinc
 

c            print *, 'increment=', kinc
c            print *, 'start to write elmts, ID', myThreadID
            open(unit = 105,
     1         file=trim(cwd)
     2              //'kinc_'//trim(adjustl(str_kinc))//'_elmts.dat',
     3         status = 'new')
            ! write elmt shear band status
            DO i = 1, num_VUEL_elmt
              write(105,1000) i, elmt_status(1,i)
            END DO
            close(105)
c            print *, 'end to write elmts, ID', myThreadID

c            print *, 'start to write nodes, ID', myThreadID
            open(unit = 105,
     1         file=trim(cwd)
     2              //'kinc_'//trim(adjustl(str_kinc))//'_theta.dat',
     3         status = 'new')
            DO i = 1, i_Array(i_int_nTotalNodes)
              write(105,2000) i, node_theta(i)
            END DO
            close(105)
c            print *, 'end to write nodes, ID', myThreadID
       
c            print *, 'start to write disp, ID', myThreadID
            open(unit = 105,
     1         file=trim(cwd)
     2              //'kinc_'//trim(adjustl(str_kinc))//'_coords.dat',
     3         status = 'new')
            DO i = 1, i_Array(i_int_nTotalNodes)
              write(105,3000) i, node_coords(:,i)
            END DO
            close(105)
c            print *, 'end to write disp, ID', myThreadID

c            print *, 'start to write SDV, ID', myThreadID
            open(unit = 105,
     1         file=trim(cwd)
     2              //'kinc_'//trim(adjustl(str_kinc))//'_SDV.dat',
     3         status = 'new')
            DO i = 1, num_VUEL_elmt
              write(105,4000) i, elmt_SDV(:,i)
            END DO
            close(105)
c            print *, 'end to write SDV, ID', myThreadID

          end if ! myThreadID.eq.0

          call MutexUnlock( 1 )


1000  FORMAT (2(I8))
2000  FORMAT (I8, F16.6)
3000  FORMAT (I8, 3(F16.6))
4000  FORMAT (I8, 10(F16.3))

        end if  ! mod(kinc,*).eq 0

c       Change i_Array(i_int_lWriteRestart) and 
C              i_Array(i_int_iStatus) if desired.       
C       Gather and export data from the configuration 
C       at the end of the current increment 
C       to external programs. 

C     End of the step
      ELSE IF (lOp .EQ. j_int_EndStep) THEN
        print *, 'End of the step'
C       In the case of multiple steps, prepare the transition to the next step.
C       For example, these data can serve as initial values for the next step.

C     End of the analysis
      ELSE IF (lOp .EQ. j_int_EndAnalysis) THEN
        print *, 'End of the analysis'
C       User coding to close files and disconnect any external programs, etc.

      END IF 

      RETURN
      
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Find all the elmts that share a face (4 nodes) with the 
C   give elmt 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE Find_Neighbour_Elmts(elmt_n,elmt_conct,
     1                                nTotalElmts,elmt_list)
C     Input:
C          elmt_n: element number;
C          elmt_conct: elmt connectivity
C          nTotalElmts: number of elmts
C     Output:
C          elmt_list: an Array contains the number of elmts that share
C                     a face with elmt_n;

      IMPLICIT NONE

      integer elmt_n
      integer elmt_list(6), nTotalElmts
      
      integer elmt_conct(8,nTotalElmts)

      integer i, j, k, elmt_1(8), elmt_2(8), counter, list_counter

C     Define 6 faces:
C     Face1: 1,2,3,4; Face2: 5,8,7,6; Face3: 1,5,6,2;
C     Face4: 2,6,7,3; Face5: 3,7,8,4; Face6: 4,5,8,1;
      integer FACE(6,4), temp(4)
      data FACE /1,5,1,2,3,4,
     1           2,8,5,6,7,8,
     2           3,7,6,7,8,5,
     3           4,6,2,3,4,1/             

C     Initialize elmt_list, one elmt has 6 neighbour elmts at most
      DO i = 1,6
        elmt_list(i) = 0
      END DO

C     elmt_1 is the connectivity of elmt_n
      DO i = 1, 8
        elmt_1(i) = elmt_conct(i,elmt_n)
      END DO   

C     Locate elmts that share Face 1 to 6 with elmt_n
      DO i = 1, nTotalElmts
        IF (i.EQ.elmt_n) THEN
          CYCLE
        END IF

        DO j = 1, 8
          elmt_2(j) = elmt_conct(j,i)
        END DO
      
        DO k = 1, 6  ! loop over 6 faces
          temp(1) = elmt_1(FACE(k,1))
          temp(2) = elmt_1(FACE(k,2))
          temp(3) = elmt_1(FACE(k,3))
          temp(4) = elmt_1(FACE(k,4))
          counter = 0

          DO j = 1, 4
            IF (ANY(elmt_2 == temp(j))) THEN
              counter = counter + 1
            END IF 
          END DO

          IF (counter.EQ.4) THEN
            elmt_list(k) = i 
          END IF

        END DO ! k
        
      END DO ! i

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Find the relative length parameter on each edge of the elmt
C   that contains a local shear band, the center of the shear band
C   surface and the volume of shear band 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SB_ct_vol(elmt_number, theta_0, natural_coords,
     1                     rBWidth, rcoords, SB_ct, SB_vol)
C     Calculate the relative length parameter on each elmt edge that 
C     intersect with the shear band surface and store the parameter
C     into the global array ID = 5, the length parameter is relative
C     to the 1st node of the edge

C     Input: 
C        elmt_number: element number
C        theta_0:     the reference temperature (level set)
C        rcoords:     elmt cartesian coords
C        rBWidth:     shear band width (from input)
C     Output:
C        natural_coords(12,3): the natural coordinates of the intersection
C                              point between elmt edge and the shear band if
C                              the point exists
C        SB_ct(3): the natural coords of shear band center
C        SB_vol:   shear band volume 

      IMPLICIT NONE
#include <SMAAspUserSubroutines.hdr> ! Start this line WITHOUT ANY SPACE

      integer elmt_number
      real*8  theta_0, rBWidth, rcoords(3,8), SB_vol

      real*8 natural_coords(12,3), SB_ct(3), Ct_x, Ct_y, Ct_z
      real*8 SB_ct_car(3), Ct_x_car, Ct_y_car, Ct_z_car
      real*8 car_coords(12,3), SHP(8), half_peri, area, area_1
      real*8 EIGHTH,ONE
      PARAMETER(EIGHTH=0.125D0,ONE=1.D0)
      integer counter

      real*8 n1(3), n2(3), Xi1, Xi2, Xi3, a1, a2, a3
      real*8, dimension(:), allocatable :: edge_length
      integer, dimension(:), allocatable :: vertice_order
      integer AllocateStatus

      real*8 elmt_edge_rlength(12,*)
      pointer(ptr_elmt_edge_rlength, elmt_edge_rlength)

      integer elmt_conct(8,*)
      pointer(ptr_elmt_conct, elmt_conct)

      real*8 node_theta(*)
      pointer(ptr_node_theta, node_theta)

      integer local_conct(8) 
      real*8  local_theta(8), temp, theta_1, theta_2, eta
      integer I,J

C     Define 12 edges:
C     Edge1:  1, 2; Edge2: 2, 3; Edge3: 3, 4; Edge4: 4, 1;
C     Edge5:  5, 6; Edge6: 6, 7; Edge7: 7, 8; Edge8: 8, 5;
C     Edge9:  1, 5; Edge10:2, 6; Edge11:3, 7; Edge12:4, 8;
      integer EDGE(12,2) 
      data EDGE /1,2,3,4,5,6,7,8,1,2,3,4,
     1           2,3,4,1,6,7,8,5,5,6,7,8/

      integer Natural_node(8,3)
      data Natural_node /-1, 1, 1, -1, -1, 1, 1, -1,
     1                   -1, -1, 1, 1, -1, -1, 1, 1,
     2                   -1, -1, -1, -1, 1, 1, 1, 1/

      ptr_elmt_conct = SMAIntArrayAccess(1)
      ptr_node_theta = SMAFloatArrayAccess(4)
      ptr_elmt_edge_rlength = SMAFloatArrayAccess(5)

      local_conct =  elmt_conct( 1:8, elmt_number)

      DO I = 1,8
        local_theta(I) = node_theta(local_conct(I))
      END DO

C     Initialize natural coordinates of all the possible intersection
C     points on the 12 elmt edges
      DO I = 1, 12
        DO J = 1,3
          natural_coords(I,J) = -10.d0
        END DO
      END DO

      DO I = 1, 12  ! loop over 12 edges
        theta_1 = local_theta( EDGE(I,1) )
        theta_2 = local_theta( EDGE(I,2) )

        temp = (theta_1 - theta_0)*(theta_2 - theta_0)
        IF (temp.LT.0.d0) THEN
          eta = (theta_0 - theta_1)/(theta_2 - theta_1)
          IF ((eta.GE.1.d0).OR.(eta.LE.0.d0)) THEN
            CALL XPLB_ABQERR(-3, 'Relative length eta ERROR',0,0.d0,' ')
          ELSE
            elmt_edge_rlength( I, elmt_number) = eta

            natural_coords(I,1) = Natural_node(EDGE(I,1),1) + eta*
     1         ( Natural_node(EDGE(I,2),1) - Natural_node(EDGE(I,1),1) )

            natural_coords(I,2) = Natural_node(EDGE(I,1),2) + eta*
     1         ( Natural_node(EDGE(I,2),2) - Natural_node(EDGE(I,1),2) )

            natural_coords(I,3) = Natural_node(EDGE(I,1),3) + eta*
     1         ( Natural_node(EDGE(I,2),3) - Natural_node(EDGE(I,1),3) )

          END IF
        END IF
      END DO

      counter = 0
      Ct_x = 0.d0
      Ct_y = 0.d0
      Ct_z = 0.d0
      Ct_x_car = 0.d0
      Ct_y_car = 0.d0
      Ct_z_car = 0.d0
      car_coords = 0.d0
      DO I = 1, 12
        IF (natural_coords(I,1).NE.-10.d0) THEN
C          print *, 'I=', I
          Ct_x = Ct_x + natural_coords(I,1)
          Ct_y = Ct_y + natural_coords(I,2)
          Ct_z = Ct_z + natural_coords(I,3)
          counter = counter + 1

          Xi1 = natural_coords(I,1)
          Xi2 = natural_coords(I,2)
          Xi3 = natural_coords(I,3)
          SHP(1) = EIGHTH*(ONE - Xi1)*(ONE - Xi2)*(ONE - Xi3)
          SHP(2) = EIGHTH*(ONE + Xi1)*(ONE - Xi2)*(ONE - Xi3)
          SHP(3) = EIGHTH*(ONE + Xi1)*(ONE + Xi2)*(ONE - Xi3)
          SHP(4) = EIGHTH*(ONE - Xi1)*(ONE + Xi2)*(ONE - Xi3)
          SHP(5) = EIGHTH*(ONE - Xi1)*(ONE - Xi2)*(ONE + Xi3)
          SHP(6) = EIGHTH*(ONE + Xi1)*(ONE - Xi2)*(ONE + Xi3)
          SHP(7) = EIGHTH*(ONE + Xi1)*(ONE + Xi2)*(ONE + Xi3)
          SHP(8) = EIGHTH*(ONE - Xi1)*(ONE + Xi2)*(ONE + Xi3)

          DO J = 1,8
            car_coords(I,1) = car_coords(I,1) + SHP(J)*rcoords(1,J)
            car_coords(I,2) = car_coords(I,2) + SHP(J)*rcoords(2,J)
            car_coords(I,3) = car_coords(I,3) + SHP(J)*rcoords(3,J)
          END DO
          Ct_x_car = Ct_x_car + car_coords(I,1)
          Ct_y_car = Ct_y_car + car_coords(I,2)
          Ct_z_car = Ct_z_car + car_coords(I,3)
        END IF
      END DO

      SB_ct(1) = Ct_x/counter
      SB_ct(2) = Ct_y/counter
      SB_ct(3) = Ct_z/counter

      SB_ct_car(1) = Ct_x_car/counter
      SB_ct_car(2) = Ct_y_car/counter
      SB_ct_car(3) = Ct_z_car/counter


C      print *, 'elmt', elmt_number
C      print *, 'elmt connection', local_conct
C      print *, 'elmt theta', local_theta
C      DO I = 1, 12
C        print *, natural_coords(I,:)
C      END DO
C      print *, 'SB_ct =', SB_ct
C      print *, 'SB_ct_car =', SB_ct_car
C      DO I = 1, 12
C        print *, car_coords(I,:)
C      END DO

      ALLOCATE(vertice_order(counter),STAT = AllocateStatus)
      IF (AllocateStatus.NE.0) THEN
        CALL XPLB_ABQERR(-3,'Shear band vertice ALLOCATION '//
     $                      'ERROR',0,0.d0,' ')
      END IF

      ALLOCATE(edge_length(counter),STAT = AllocateStatus)
      IF (AllocateStatus.NE.0) THEN
        CALL XPLB_ABQERR(-3,'Shear band vertice ALLOCATION '//
     $                      'ERROR',0,0.d0,' ')
      END IF

      counter = 0
      DO I = 1, 12
        IF (natural_coords(I,1).NE.-10.d0) THEN
          counter = counter + 1
          vertice_order(counter) = I
        END IF
      END DO
C      print *, 'counter=', counter
C      print *, 'vertice_order =', vertice_order

C     Calculate the shear band volume
      area = 0.d0
      IF (counter.EQ.3) Then  ! The intersected surface is triangle
        edge_length(1) 
     1   = sqrt(  ( car_coords(vertice_order(1),1) 
     2             -car_coords(vertice_order(2),1) )**2
     3          + ( car_coords(vertice_order(1),2) 
     4             -car_coords(vertice_order(2),2) )**2
     5          + ( car_coords(vertice_order(1),3) 
     6             -car_coords(vertice_order(2),3) )**2 )
        edge_length(2) 
     1   = sqrt(  ( car_coords(vertice_order(2),1) 
     2             -car_coords(vertice_order(3),1) )**2
     3          + ( car_coords(vertice_order(2),2) 
     4             -car_coords(vertice_order(3),2) )**2
     5          + ( car_coords(vertice_order(2),3) 
     6             -car_coords(vertice_order(3),3) )**2 )
        edge_length(3) 
     1   = sqrt(  ( car_coords(vertice_order(1),1) 
     2             -car_coords(vertice_order(3),1) )**2
     3          + ( car_coords(vertice_order(1),2) 
     4             -car_coords(vertice_order(3),2) )**2
     5          + ( car_coords(vertice_order(1),3) 
     6             -car_coords(vertice_order(3),3) )**2 )
        half_peri = 0.5d0*(  edge_length(1) + edge_length(2)
     1                     + edge_length(3) )
        area = sqrt(   half_peri
     1              * (half_peri - edge_length(1))
     2              * (half_peri - edge_length(2)) 
     3              * (half_peri - edge_length(3)) )
        SB_vol = area * rBWidth
      ELSE ! the intersected surface is polygon

C       order the vertices
        CALL sort_vert_2(SB_ct, natural_coords,
     1                 counter, vertice_order)     
C        print *, 'modified vertice_order =', vertice_order

        DO I = 1, (counter-1)
          edge_length(I)
     1   = sqrt(  ( car_coords(vertice_order(I+1),1) 
     2             -car_coords(vertice_order(I)  ,1) )**2
     3          + ( car_coords(vertice_order(I+1),2) 
     4             -car_coords(vertice_order(I)  ,2) )**2
     5          + ( car_coords(vertice_order(I+1),3) 
     6             -car_coords(vertice_order(I)  ,3) )**2 )
        END DO
        edge_length(counter)
     1   = sqrt(  ( car_coords(vertice_order(counter),1) 
     2             -car_coords(vertice_order(1)      ,1) )**2
     3          + ( car_coords(vertice_order(counter),2) 
     4             -car_coords(vertice_order(1)      ,2) )**2
     5          + ( car_coords(vertice_order(counter),3) 
     6             -car_coords(vertice_order(1)      ,3) )**2 )
        DO I = 1, (counter-1)
          a1 = edge_length(I)
          a2 = sqrt( ( SB_ct_car(1) 
     1               - car_coords(vertice_order(I),1) )**2
     2             + ( SB_ct_car(2) 
     3               - car_coords(vertice_order(I),2) )**2
     4             + ( SB_ct_car(3) 
     5               - car_coords(vertice_order(I),3) )**2 )
          a3 = sqrt( ( SB_ct_car(1) 
     1               - car_coords(vertice_order(I+1),1) )**2
     2             + ( SB_ct_car(2) 
     3               - car_coords(vertice_order(I+1),2) )**2
     4             + ( SB_ct_car(3) 
     5               - car_coords(vertice_order(I+1),3) )**2 )
          half_peri = 0.5d0*(a1 + a2 + a3)
          area_1 = sqrt(  half_peri
     1                  *(half_peri-a1)
     2                  *(half_peri-a2)
     3                  *(half_peri-a3)  )
          area = area + area_1
        END DO
        a1 = edge_length(counter)
        a2 = sqrt( ( SB_ct_car(1) 
     1             - car_coords(vertice_order(counter),1) )**2
     2           + ( SB_ct_car(2) 
     3             - car_coords(vertice_order(counter),2) )**2
     4           + ( SB_ct_car(3) 
     5             - car_coords(vertice_order(counter),3) )**2 )
        a3 = sqrt( ( SB_ct_car(1) 
     1             - car_coords(vertice_order(1),1) )**2
     2           + ( SB_ct_car(2) 
     3             - car_coords(vertice_order(1),2) )**2
     4           + ( SB_ct_car(3) 
     5             - car_coords(vertice_order(1),3) )**2 )
        half_peri = 0.5d0*(a1 + a2 + a3)
        area_1 = sqrt(  half_peri
     1                *(half_peri-a1)
     2                *(half_peri-a2)
     3                *(half_peri-a3)  )
        area = area + area_1
C        print *, 'area =', area
        SB_vol = area * rBWidth
      END IF 

C     Deallocate arrays
      IF (allocated(vertice_order)) DEALLOCATE(vertice_order)

      IF (allocated(edge_length)) DEALLOCATE(edge_length)

      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     INITIALIZE A MATRIX
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INIT_MTX(MATRIX,D1,D2)

      IMPLICIT NONE
      integer D1,D2,I,J
      real*8 ZERO,MATRIX(D1,D2)
      PARAMETER(ZERO=0.D0)

      DO I = 1,D1
        DO J = 1,D2
          MATRIX(I,J) = ZERO
        END DO
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     INITIALIZE A VECTOR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INIT_VEC(VECTOR,D1)

      IMPLICIT NONE

      integer D1,I
      real*8 ZERO,VECTOR(D1)
      PARAMETER(ZERO=0.D0)

      DO I = 1,D1
        VECTOR(I) = ZERO
      END DO
      
      RETURN
      END




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Sort the order of the vertices, a different scheme
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE sort_vert_2(SB_ct, natural_coords,
     1                       counter,vertice_order)

C     Output:
C       vertice_order: the sorted order

      IMPLICIT NONE
      integer counter, vertice_order(counter), pos(1)
      real*8 SB_ct(3), natural_coords(12,3)
      real*8 value(1), xmax, ymax, zmax, xdiff, ydiff, zdiff
      real*8 temp_x(counter), temp_y(counter), angle(counter)
      integer I,J,location(counter),temp_index

      xmax = 0.d0
      ymax = 0.d0
      zmax = 0.d0
      DO I = 1, counter-1
        DO J = 1, counter-I
          xdiff = abs(natural_coords(vertice_order(I),1) 
     1              - natural_coords(vertice_order(I+J),1))
          ydiff = abs(natural_coords(vertice_order(I),2) 
     1              - natural_coords(vertice_order(I+J),2))
          zdiff = abs(natural_coords(vertice_order(I),3) 
     1              - natural_coords(vertice_order(I+J),3))
          IF (xdiff.GT.xmax) THEN
            xmax = xdiff
          END IF
          IF (ydiff.GT.ymax) THEN
            ymax = ydiff
          END IF
          IF (zdiff.GT.zmax) THEN
            zmax = zdiff
          END IF
        END DO
      END DO

      pos = minloc((/xmax,ymax,zmax/))

C      print *, 'pos = ', pos(1)

      IF (pos(1).EQ.1) THEN ! yz plane
        DO I = 1, counter
          temp_x(I) = natural_coords(vertice_order(I),2) - SB_ct(2)
          temp_y(I) = natural_coords(vertice_order(I),3) - SB_ct(3)
        END DO
      ELSE IF (pos(1).EQ.2) THEN ! xz plane
        DO I = 1, counter
          temp_x(I) = natural_coords(vertice_order(I),1) - SB_ct(1)
          temp_y(I) = natural_coords(vertice_order(I),3) - SB_ct(3)
        END DO
      ELSE IF (pos(1).EQ.3) THEN ! xy plane
        DO I = 1, counter
          temp_x(I) = natural_coords(vertice_order(I),1) - SB_ct(1)
          temp_y(I) = natural_coords(vertice_order(I),2) - SB_ct(2)
        END DO
      END IF

      DO I = 1, counter
        angle(I) = atan2(temp_y(I), temp_x(I))
        IF (angle(I).LT.0.d0) THEN
          angle(I) = angle(I) + 8.d0*atan2(1.d0,1.d0)
        END IF
        location(I) = vertice_order(I)
      END DO

C      print *, 'angle = ', angle
C      print *, 'temp_x = ', temp_x
C      print *, 'temp_y = ', temp_y

      DO I = 1, (counter-1)
        pos = minloc(angle(i:counter))
        value = minval(angle(i:counter))
        temp_index = location(I)
        location(I) = location(pos(1)+I-1)
        location(pos(1)+I-1) = temp_index 
        angle(pos(1)+I-1) = angle(I)
        angle(I)  = value(1)
      END DO

      vertice_order = location
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Solve the heat equation globally
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE solve_heat_global(total_node,total_elmt,node_coords, 
     1           elmt_conct,node_theta,elmt_Nvect,elmt_status,BC_index)

C     input:
C        total_node, total number of nodes
C        total_elmt, total number of elmts
C        node_coords(3,total_node), nodal coordinates (deformed configuration)
C        elmt_conct(8,total_elmt), elmt connectivity
C        node_theta(total_node), nodal temperature
C        elmt_Nvect(3,total_elmt), shear band normal direciton at elmt
C                                  centroid
C        elmt_status(2,total_elmt), elmt shearband/frontlist status
C        BC_index(total_node), whether a node is applied with essential B.C.
C                              1 - yes; 0 - no;

C     Output:
C        node_theta(total_node), updated nodal temperature

      IMPLICIT NONE

C     Variables from the interface
      integer total_node, total_elmt
      integer elmt_conct(8,total_elmt), elmt_status(2,total_elmt)
      real*8  node_theta(total_node)
      real*8  elmt_Nvect(3,total_elmt), node_coords(3,total_node)
  
C     Variables needed in the subroutine
      real*8, dimension(:,:), allocatable ::  global_stiffness
      real*8, dimension(:), allocatable :: global_rhs
      integer, dimension(:), allocatable :: ipiv
      real*8 eyes(3,3), conduct(3,3)
      integer BC_index(total_node), neq, id(total_node)
      integer i, j, AllocateStatus, counter, e, k, k1, k2, p, q
      real*8 stable_para
      parameter(stable_para=1.d-6)
      integer iTIntPt, iErrFlg, info
      real*8  rXiEtaPt(3,8),rWeightPt(8),local_stiffness(8,8)
      real*8  SHP(8), SHP_C(3,8), JACOB, elmt_coord(3,8), WE 

      global_stiffness = 0.d0
      global_rhs = 0.d0
      eyes = 0.d0
      do i = 1,3
        eyes(i,i) = 1.d0
      end do
  
C     Find out which node has the essential B.C.    
C     If a shear band has developed inside an elmt,
C     all the nodes related with elmt will be assigned
C     as essential B.C.
      do i = 1, total_elmt
        if (elmt_status(1,i).eq.1) then
          do j = 1, 8
            BC_index( elmt_conct(j,i) ) = 1
          end do
        end if
      end do ! i
C     total number of unknowns (equations) 
      neq = total_node - sum(BC_index)
C      print *, 'total node=', total_node
C      print *, 'neq=', neq

C     Allocate the array for global stiffness matrix
      allocate ( global_stiffness(neq, neq), 
     1           STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
        CALL XPLB_ABQERR(-3,'HEAT EQ. MATRIX ALLOCATION ERROR',
     1                      0,0.d0,' ')
      end if
      global_stiffness = 0.d0

C     Allocate the RHS
      allocate ( global_rhs(neq),
     1           STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
        CALL XPLB_ABQERR(-3, 'HEAT EQ. RHS ALLOCATION ERROR',0,0.d0,' ')
      end if
      global_rhs = 0.d0

C     Allocate the RHS
      allocate ( ipiv(neq),
     1           STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
        CALL XPLB_ABQERR(-3,'HEAT EQ. IPIV ALLOCATION ERROR',0,0.d0,' ')
      end if
      ipiv = 0

C     Equation number for each node
      counter = 1
      do i = 1, total_node
        if (BC_index(i).eq.1) then
          id(i) = 0
        else
          id(i) = counter
          counter = counter + 1
        end if
      end do
      counter = counter - 1
      if (neq.ne.counter) then
        CALL XPLB_ABQERR(-3, 'HEAT EQ. number of Eqs error',0,0.d0,' ')
      end if    

C     Loop over elmt to assemble global_stiffness and global_rhs
C     DETERMINE GAUSS PTS LOCATIONS & WEIGHTS (2x2x2)
      CALL GSPT_3D_8N(iTIntPt, rXiEtaPt, rWeightPt)

      do e = 1, total_elmt
        local_stiffness = 0.d0
        do i = 1, 3
          do j = 1, 3
            conduct(i,j) = (1.d0 + stable_para)*eyes(i,j) 
     1                   - elmt_Nvect(i,e)*elmt_Nvect(j,e) 
          end do
        end do

C       deformed configuration
        do i = 1, 8
          elmt_coord(:,i) = node_coords(:, elmt_conct(i,e) )
        end do

C       Loop over eight Gauss Points
        do k = 1, 8
C         Evaluate interpolation func, their derivatives and
C         Jacobians at Gauss point
          call SHP_FUNC_3D_8N(rXiEtaPt(1,k),elmt_coord,1,
     1                        SHP,SHP_C,JACOB)
          WE = rWeightPt(k)*JACOB
          do i = 1, 8
            do j = 1, 8
              do k1 = 1,3
                do k2 = 1,3
                  local_stiffness(i,j) = local_stiffness(i,j)
     1            + WE*SHP_C(k1,i)*conduct(k1,k2)*SHP_C(k2,j)
                end do
              end do
            end do
          end do
        end do ! k

        do i = 1, 8
          do j = 1, 8
            p = id(elmt_conct(i,e))
            q = id(elmt_conct(j,e))
C           assemble global stiffness
            if (p*q.ne.0) then
              global_stiffness(p,q) = global_stiffness(p,q)
     1                              + local_stiffness(i,j)
            end if
C           assemble global rhs
            if (p.ne.0) then
              if (q.eq.0) then
                global_rhs(p) = global_rhs(p) 
     1               - local_stiffness(i,j)*node_theta(elmt_conct(j,e))
              end if
            end if
          end do
        end do
      end do ! loop over elmt

C      print *, 'Nvect = ', elmt_Nvect(:,2041)
C      print *, global_stiffness(1,1), global_stiffness(2,2)

C      print *, 'reach here'
      call dgesv( neq, 1, global_stiffness, neq, ipiv, global_rhs,
     1            neq, info )
C      print *, 'reach here'
      if (info.ne.0) then
        call XPLB_ABQERR(-3, 'HEAT EQ. Linear Eq. solve ERROR',
     1                        0,0.d0,' ')
      end if

      do i = 1, total_node
        if (id(i).ne.0) then
          node_theta(i) = global_rhs(id(i)) 
        end if
      end do

C      do i = 1, total_node
C        print *, i, node_theta(i)
C     end do

C      print *, '346:', node_theta(346)
C      print *, '5:', node_theta(5)
C      print *, '61:', node_theta(61)
C      print *, '1353:', node_theta(1353)      

      if (allocated(global_stiffness)) deallocate(global_stiffness)
      if (allocated(global_rhs))       deallocate(global_rhs)
      if (allocated(ipiv))             deallocate(ipiv)

      END ! SUBROUTINE


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Calculate the equivalent normal for the non-planar surface case 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE  non_plane_normal(SB_ct,natural_coords,normalCt)

C     Input: 
C        SB_ct: the center of the shear band
C        natural_coords: the natural coordinates of the intersections between
C                        the elmt edge and the shear band surface
C     Output:
C        normalCt:  the equivalent normal of the non-planar surface

      IMPLICIT NONE

      real*8 SB_ct(3), natural_coords(12,3), normalCt(3), rNorm
      real*8 temp(3), M(3,3), rEigValR(3), rEigValI(3), rEigVec(3,3)
      integer I,J,K,counter,iErrCode, imin(1)

      counter = 0
      M = 0.d0
      DO I = 1, 12
        IF ( natural_coords(I,1).NE.-10.d0 ) THEN
          temp = natural_coords(I,:) - SB_ct
          DO J = 1,3
            DO K = 1,3
              M(J,K) = M(J,K) + temp(J)*temp(K)
            END DO
          END DO
          counter = counter + 1
        END IF
      END DO

C     more than three intersection points
      IF (counter.Gt.3) THEN
        M = M/counter

        call rGEigV(3,3,M,rEigValR,rEigValI,1,rEigVec,iErrCode)

c       Order the principal stresses and their associated directions

        do i = 1, 3
          if ((rEigValI(i).ne.0.d0).or.(iErrCode.ne.0)) then
            call XPLB_ABQERR(-3,'Problem solving non-planar surface'//
     $       'normal by eigenvalue decomposition',0,0.d0,' ')
          endif
        enddo

        imin = minloc(abs(rEigValR))
        normalCt = rEigVec(:,imin(1))
        rNorm = sqrt(normalCt(1)**2 + normalCt(2)**2 + normalCt(3)**2)
        normalCt = normalCt/rNorm

      END IF

      END ! SUBROUTINE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Core FEM subroutines: quadrature, shape functions, etc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      include 'GSPT_3D_8N.f'
      include 'SHP_FUNC_3D_8N.f'
      include 'Clength_3D_8N.f'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Mechanics subroutines: get strains, stresses, etc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      include 'Gstrain_3D.f'
      include 'gdefgrad_3D.f'
      include 'grottens.f'

      include 'Gstress_3D.f'
      include 'hypoJ2.f'
      include 'gwavespd.f'
      include 'vumbka.f'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Subgrid method subroutines: bifurcation analysis, etc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      include 'perturb_3D.f'
      include 'surface_normal.f'
      include 'galtinit_3D.f'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   EISPACK Subroutines to compute eigenvalues and eigenvectors of a general real
C   matrix.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      include 'balanc.f'
      include 'balbak.f'
      include 'cdiv.f'
      include 'elmhes.f'
      include 'eltran.f'
      include 'hqr2.f'
      include 'hqr.f'
      include 'rg.f'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Utility subroutines: invert matrix, find roots of cubic polynomial, etc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      include 'minvert.f'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Dummy VUMAT subroutine for post-processing ISVs
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      include 'ppvumat_3D.f'
      include 'prep_umat.f'
      include 'fcc_umat_lbl34_test.f'