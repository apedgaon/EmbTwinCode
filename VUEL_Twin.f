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
#include <SMAAspUserSubroutines.hdr>
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
        integer kblock, i, j, k, l, n, m, p, iter, aidx
        integer iTIntPt          ! Total # of Gauss PTS (=8) for MATRIX
        integer elmt_number      ! elmt global number
        integer iumatsvst        ! umat start index for svars
        integer iumatsven        ! umat end index of svars
        integer iumatsvst2        ! umat start index for svars
        integer iumatsven2        ! umat end index of svars
        integer iVelmType
        integer iVforType
        integer iConsType
        integer nstatev(4)
        integer offstatev(7)
        integer iTwFlg           ! Twinning index that nucleated
        integer num_gb_elems     ! Number of Grain Boundary Elements
        integer ieulererr        ! return error from eulang routine
        integer iterErr          ! Iteration convergence error
        integer iTwEmbStp        ! Twin Embedding Step

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
        real*8 rEpsilon          ! Pressure mode stabilization parameter; see Ortiz et al (1987), Eq. (21).
        real*8 rMSBm(6,mcrd*nnode) ! B (or B-bar) operator at mid-step
        real*8 rMSBs(3,mcrd*nnode) ! skew-symmetric matrix operator at mid-step
        real*8 rMSdNndXiCt(mcrd,nnode)   ! Shape func derivatives at Elmt centroid at Mid-step
        real*8 rMSJdetCt                 ! jacobian at Elmt centroid at Mid-step
        real*8 rMSdNndXi(mcrd,nnode)     ! Shape func derivatives at one Gauss PT at Mid-step 
        real*8 rMSJdet                   ! Jacobian at one Gauss PT at Mid-step
        real*8 rStrain(6)          ! strain increment
        real*8 rSpin(3)            ! spin axial vector increment
        real*8 rue(mcrd,nnode)     ! local array for passing displacements increment to strain routine
        real*8 rdVelGrd(3, 3)     ! Velocity Gradient Increment
        real*8 del(3, 3)        ! 2nd Order Identity Tensor
        real*8 rDfgrd0(3, 3)      ! Previous Deformation Gradient
        real*8 rDfgrd1(3, 3)      ! Current Deformation Gradient
        real*8 rGsCoords(3)       ! Quadrature Point Current coordinates
        real*8 rStress(6)         ! Cauchy Stress
        real*8 rBmat(6,mcrd*nnode) ! current B (or B-bar) operator
        real*8 rBskw(3,mcrd*nnode) ! current skew-symmetric matrix operator
        real*8 Smat(3, 3, 18)      ! Twinning Schmid Tensor
        real*8 F_tw(3, 3, 18)      ! Twinning Deformation gradient
        real*8 rm11                ! C11 vs temperature slope
        real*8 c11_0               ! C11 at 0 K
        real*8 c11                 ! C11 at temperature theta
        real*8 theta               ! Temperature
        real*8 rCd                 ! Dilitational wave speed
        real*8 rdLu(3, 3)          ! Velocity gradient increment in Untwinned part
        real*8 rdLz(3, 3)          ! Velocity gradient increment in Twinned part
        real*8 rdFu(3, 3)          ! Incremental deformation gradient in the untwinned part
        real*8 rdFz(3, 3)          ! Incremental deformation gradient in the twinned part
        real*8 rFu0(3, 3)          ! Previous untwinned deformation gradient
        real*8 rFu1(3, 3)          ! Current untwinned deformation gradient
        real*8 rFz0(3, 3)          ! Previous twinned deformation gradient
        real*8 rFz0h(3, 3)         ! Fe Fp part of previous twinned deformation gradient
        real*8 rFz1(3, 3)          ! Current twinned deformation gradient
        real*8 rFz1h(3, 3)         ! Fe Fp part of current twinned deformation gradient
        real*8 rBWidth             ! Band width of twin
        real*8 F_tw_inv(3, 3)      ! Inverse of twinned deformation gradient
        real*8 rStress_U(6)     ! Stress in untwinned region
        real*8 rStress_Z(6)     ! Stress in twinned region
        real*8 rWu(3, 3, 3, 3)     ! dT/dF^U Material Jacobian for untwinned region
        real*8 rWz(3, 3, 3, 3)     ! dT/dF^Z Material Jacobian for twinned region
        real*8 rntw(3, 18)         ! Twin normal
        real*8 rmtw(3, 18)         ! Twin shear direction
        real*8 props_tw(50)        ! Properties for twinned region since euler angles have changed
        ! TODO: Remove hard coded 50
        real*8 pi                  ! PI
        real*8 g_tr(3)             ! residual for traction continuity
        real*8 z_tr(3, 3)          ! tangent matrix for traction continuity
        real*8 z_tr_inv(3, 3)      ! inverse of tangent matrix for traction continuity
        real*8 sigmaZ(3, 3)        ! Tensor Stress in twinned region
        real*8 sigmaU(3, 3)        ! Tensor Stress in untwinned region
        real*8 rdmh(3)             ! change in jump + direction
        real*8 rmh(3)              ! jump + direction
        real*8 rntwc(3)            ! Twin normal in current configuration
        real*8 rFuinv(3, 3)        ! inverse of def grad in untwinned region
        real*8 mtpn(3, 3)          ! mhat tensor product n
        real*8 rmtwc(3)            ! Twin Direction in current configuration
        real*8 rmhu(3)             ! Unit vector along mhat
        real*8 rBmatU(6,mcrd*nnode)! B-bar operator in untwinned region
        real*8 rBmatZ(6,mcrd*nnode)! B-bar operator in twinned region
        real*8 rTdis(6)            ! Discontinuity tensor (m tensor product n)
        real*8 alpha_h             ! Intermediate discontinuity magnitude
        real*8 alpha_U             ! Discontinuity magnitude in untwinned region
        real*8 alpha_Z             ! Discontinuity magnitude in twinned region
        real*8 rTdL                ! Velocity gradient increment projected on discontinous mode
        real*8 rTBmat              ! B-bar matrix projected on discontinous mode
        real*8 g_tr_err            ! Error norm in residual
        real*8 rEpVolInc           ! Volumetric strain increment
        real*8 rbv1, rbv2          ! Linear and quadratic bulk viscosity coefficent
        real*8 rBVpress            ! Bulk viscosity pressure
        real*8 rSigma              ! Adjusted stress for RHS
        real*8 rSigmaZ             ! Adjusted stress for twinned part for RHS
        real*8 rSigmaU             ! Adjusted stress for untwinned part for RHS
        real*8 tau_tw              ! RSS on twin system
        real*8 rBLength            ! Length of twin inside the element
        real*8 Smatc(3, 3)         ! Specific Twin Schmid tensor in current configuration

        integer gb_elems(*)        ! Grain boundary element array
        pointer (ptr_gb_elems, gb_elems) ! pointer for gb_elems required for SMAArrayAccess
        logical isGbElem           ! logical determining gb elem
        logical, external :: is_gb_elem

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

        rEpsilon  = props(10)  ! Pressure mode stabilization parameter; see Ortiz et al (1987), Eq. (21).

        !call sleep(20)
        rRho = props(78)
        call GSPT_3D_8N(iTIntPt, rXiEtaPt, rWeightPt)

C       Natural coordinates of element centroid
        rXiEtaCt(1) = 0.d0
        rXiEtaCt(2) = 0.d0
        rXiEtaCt(3) = 0.d0

        ! Number of state variables array and offsets array
        nstatev(1) = 148  ! svars for umat
        nstatev(2) = 9    ! svars for deformation gradient
        nstatev(3) = 2    ! svars for twinning flag
        nstatev(4) = 3    ! svars for storing discontinuity vector
        offstatev(1) = 0
        offstatev(2) = 8 * nstatev(1)
        offstatev(3) = offstatev(2) + 8 * nstatev(1)
        offstatev(4) = offstatev(3) + 8 * nstatev(2)
        offstatev(5) = offstatev(4) + 8 * nstatev(2)
        offstatev(6) = offstatev(5) + nstatev(3)
        offstatev(7) = offstatev(6) + 8 * nstatev(4)
        ! 148*8 + 148*8 + 9*8 + 9*8 + 2 + 3*8

        ! Define PI
        pi = 4.d0 * atan(1.d0)

C       Initialize Identity Matrix(3x3)
        do i = 1, 3
          do j = 1, 3
            if (i.eq.j) then
                del(i, j) = 1.0d0;
            else
                del(i, j) = 0.0d0;
            end if
          end do
        end do

C       Compute equivalent Schmid Tensor and Deformation gradient for Deformation Twinning
        call comp_twin_dfgrd(del, Smat, F_tw, rntw, rmtw)

        ! Access grain boundary elements
        if (kstep.gt.0) then
          num_gb_elems = 488  ! remove hard coding
          ptr_gb_elems = SMAIntArrayAccess(1)
        end if

        ! Compute dilitational wave speed
        rm11  = props(82)
        c11_0 = props(83)
        theta = props(77)
        c11 = c11_0 + rm11 * theta ! based on varying temperature?
        rCd = sqrt(c11 / rRho)

        ! Linear bulk viscosity coefficient
        if (props(8).le.0.d0) then
          rbv1    = 0.06d0
        else
          rbv1    = props(8)
        end if

        ! Quadratic bulk viscosity coefficient
        if (props(9).le.0.d0) then
          rbv2    = 1.2d0
        else
          rbv2    = props(9)
        end if

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
     $                 + rRho * rNn(n) * rdOmega   ! check \sum N(A)*rho*N(B)*dOmega (unlumped)
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
          ! Loop over elements in the block
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
     &                            rNn,rdNndXi,rJdet)
C             Integration element (dOmega)
              rdOmega = rJdet * rWeightPt(l)
C             Element surface area (even for axisymmetric elements)
              rElmtVol = rElmtVol + rdOmega
            END DO

C           Compute element size h
            rElmtSize = rElmtVol**(1.d0/3.d0)      !   h = cubic root of vol

C           a) shape functions and their derivatives wrt global coords
            CALL SHP_FUNC_3D_8N(rXiEtaCt,rCoords,1,rNnCt,
     1                          rdNndXiCt,rJdetCt)
C           b) shape functions and their derivatives wrt mid-step global coords
            CALL SHP_FUNC_3D_8N(rXiEtaCt,rMSCoor,1,rNnCt,
     1                          rMSdNndXiCt,rMSJdetCt)

            ! GB Element Check
            if (kstep.gt.0) then
              isGbElem = is_gb_elem(elmt_number, gb_elems, num_gb_elems)
            end if

            ! Set twinning flags 0 in the beginning
            if (kinc.eq.0) then
              svars(kblock, offstatev(5) + 1) = 0
              svars(kblock, offstatev(5) + 2) = 0
            end if

            ! Check twinning flag in history
            iTwFlg = svars(kblock, offstatev(5) + 1)
            iTwEmbStp = svars(kblock, offstatev(5) + 2)

            if (iTwFlg.eq.0) then   ! No twinning detected in previous increment
              ! Quadrature Loop
              do l = 1, iTIntPt
C               Compute shape functs, their derivs (wrt global coords), and the Jacobian determinant
                CALL SHP_FUNC_3D_8N(rXiEtaPt(1,l),rCoords,1,
     &                              rNn,rdNndXi,rJdet)

C               Compute Quadrature Point Current Coordinates
                rGsCoords = 0.d0
                do i = 1, mcrd
                  do j = 1, nnode
                    rGsCoords(i) = rGsCoords(i) + rNn(j) * rCoords(i, j)
                  end do
                end do

C               Check Jacobian
                IF ( rJdet .LE. 0.d0) THEN
                  CALL XPLB_ABQERR(-3,'Non-positive'//
     $            'Jacobian determinant in element VU1',0,0.d0,' ')
                END IF

C               Compute integration element (dOmega)
                rdOmega = rJdet * rWeightPt(l)

C               Finite deformation
C               Current value of the B or B-bar operator for internal force calculation
                CALL Gstrain_3D(rue,rEpsilon,rNn,rNnCt,rdNndXi,
     $                          rdNndXiCt,mcrd,nnode,iVelmType,
     $                          iVforType,1,rStrain,rSpin,rBmat,rBskw)

C               (Mid-step) strain increment and B or B-bar operator for constitutive update
                CALL SHP_FUNC_3D_8N(rXiEtaPt(1,l),rMSCoor,1,rNn,
     $                             rMSdNndXi,rMSJdet)
                CALL Gstrain_3D(rue,rEpsilon,rNn,rNnCt,rMSdNndXi,
     $                          rMSdNndXiCt,mcrd,nnode,iVelmType,
     $                          iVforType,1, rStrain,rSpin,rMSBm,
     $                          rMSBs)
C               Until now, rStrain and rSpin are quantities wrt to mid-step

C               Compute volumetric strain increment (for bulk viscosity calculation)
                rEpVolInc = rStrain(1) + rStrain(2) + rStrain(3)

                ! Compute Velocity Gradient
                call comp_vel_gr(rdVelGrd, rStrain, rSpin)

                ! Step 1, Inc 1, set def gradient to identity
                iumatsvst = offstatev(3) + 1 + (l - 1) * nstatev(2)
                iumatsven = offstatev(3) + l * nstatev(2)
                if ((kstep.eq.1).and.(kinc.le.1)) then
                  svars(kblock, iumatsvst : iumatsven) =
     &                                    reshape(del, (/9/))
                end if

                ! Get previous def gradients from svars and form current def grd
                rDfgrd0 =
     &            reshape(svars(kblock, iumatsvst : iumatsven), (/3,3/))
                rDfgrd1 = matmul(del + rdVelGrd, rDfgrd0)

                ! Save deformation gradients to both untwinned and twinned block
                svars(kblock, iumatsvst : iumatsven) =
     &                                  reshape(rDfgrd1, (/9/))

                iumatsvst = offstatev(4) + 1 + (l - 1) * nstatev(2)
                iumatsven = offstatev(4) + l * nstatev(2)
                svars(kblock, iumatsvst : iumatsven) =
     &                                  reshape(rDfgrd1, (/9/))

                ! Compute Stress from existing UMAT
                iumatsvst = offstatev(1) + 1 + (l - 1) * nstatev(1)
                iumatsven = offstatev(1) + l * nstatev(1)
                call umat_stress(rStress,
     &                           svars(kblock, iumatsvst:iumatsven),
     &                           rGsCoords, rDfgrd0, rDfgrd1,
     &                           props(77:126), 50, nstatev(1), nnode,
     &                           elmt_number, l, mcrd, time, dtimeCur,
     &                           kstep, kinc, .true., rWu, iterErr)

                if (iterErr.eq.1) then
                  CALL XPLB_ABQERR(-3,'UMAT Failed to Converge in '//
     &                             'unlocalized region!', 0,0.d0,' ')
                end if

                ! store same svars for twinned part which has not formed yet
                iumatsvst2 = offstatev(2) + 1 + (l - 1) * nstatev(1)
                iumatsven2 = offstatev(2) + l * nstatev(1)
                svars(kblock, iumatsvst2:iumatsven2) =
     &                                svars(kblock, iumatsvst:iumatsven)

                ! Check for twin nucleation
                if (isGbElem) then
                  if (iTwFlg.eq.0) then
                    call check_twin_nucleation(rStress, Smat, rDfgrd1,
     &                                         iTwFlg)
                  end if
                end if

C               Compute linear + quadratic bulk viscosity terms
C               Needed for damping the “ringing” in the highest element frequency
C               (Abaqus->Theory->Procedures->Nonlinear dynamics->Explicit dynamic analysis)
                !! Tao,The quadratic bulk viscosity is applied
                !! only if the volumetric strain rate is compressive.
                if (dtimeCur.ne.0.d0) then
                  rBVpress = -rRho *
     $                        rbv1*rCd*rElmtSize*(rEpVolInc/dtimeCur)
                  if (rEpVolInc.lt.0.d0) then
                    rBVpress = rBVpress
     $                      -rRho*(rbv2*rElmtSize*rEpVolInc/dtimeCur)**2
                  end if
                else
                  rBVpress = 0.d0
                end if

C               Compute the internal force vector
                DO n = 1, nnode
                  DO i = 1, mcrd
                    DO j = 1, 6
                      if (j.le.3) then    ! add bulk viscosity
                        rSigma = rStress(j) - rBVpress
                      else
                        rSigma = rStress(j)
                      end if
                      rhs( kblock , (n-1)*mcrd + i ) =
     $                rhs( kblock , (n-1)*mcrd + i )
     $                + rBmat( j , (n-1)*mcrd + i ) * rSigma * rdOmega
                    END DO
                  END DO
                END DO

              end do

              ! Store twin flag
              svars(kblock, offstatev(5) + 1) = iTwFlg

              ! If twin detected, mark next step for twin embedding
              if (iTwFlg.gt.0) then
                iTwEmbStp = 1
                svars(kblock, offstatev(5) + 2) = iTwEmbStp
              end if

            else  ! Twinning detected in the previous increment

              ! compute inverse of twinning deformation gradient
              call inverse_3x3(F_tw(:, :, iTwFlg), F_tw_inv)

              ! compute properties for twinned region due to rotation of crystal basis
              call compute_twin_props(rntw(:, iTwFlg), props(77:126),
     &                                props_tw, 50, del, pi)

              ! Quadrature Loop
              do l = 1, iTIntPt
C               Compute shape functs, their derivs (wrt global coords), and the Jacobian determinant
                CALL SHP_FUNC_3D_8N(rXiEtaPt(1,l),rCoords,1,
     &                              rNn,rdNndXi,rJdet)

C               Compute Quadrature Point Current Coordinates
                rGsCoords = 0.d0
                do i = 1, mcrd
                  do j = 1, nnode
                    rGsCoords(i) = rGsCoords(i) + rNn(j) * rCoords(i, j)
                  end do
                end do

C               Check Jacobian
                IF ( rJdet .LE. 0.d0) THEN
                  CALL XPLB_ABQERR(-3,'Non-positive'//
     $            'Jacobian determinant in element VU1',0,0.d0,' ')
                END IF

C               Compute integration element (dOmega)
                rdOmega = rJdet * rWeightPt(l)

C               Finite deformation
C               Current value of the B or B-bar operator for internal force calculation
                CALL Gstrain_3D(rue,rEpsilon,rNn,rNnCt,rdNndXi,
     $                          rdNndXiCt,mcrd,nnode,iVelmType,
     $                          iVforType,1,rStrain,rSpin,rBmat,rBskw)


C               (Mid-step) strain increment and B or B-bar operator for constitutive update
                CALL SHP_FUNC_3D_8N(rXiEtaPt(1,l),rMSCoor,1,rNn,
     $                             rMSdNndXi,rMSJdet)
                CALL Gstrain_3D(rue,rEpsilon,rNn,rNnCt,rMSdNndXi,
     $                          rMSdNndXiCt,mcrd,nnode,iVelmType,
     $                          iVforType,1, rStrain,rSpin,rMSBm,
     $                          rMSBs)
C               Until now, rStrain and rSpin are quantities wrt to mid-step

C               Compute volumetric strain increment (for bulk viscosity calculation)
                rEpVolInc = rStrain(1) + rStrain(2) + rStrain(3)

                ! Retrieve previous untwinned and twinned deformation gradient
                iumatsvst = offstatev(3) + 1 + (l - 1) * nstatev(2)
                iumatsven = offstatev(3) + l * nstatev(2)
                rFu0 =
     &            reshape(svars(kblock, iumatsvst : iumatsven), (/3,3/))
                iumatsvst = offstatev(4) + 1 + (l - 1) * nstatev(2)
                iumatsven = offstatev(4) + l * nstatev(2)
                rFz0 =
     &            reshape(svars(kblock, iumatsvst : iumatsven), (/3,3/))

                ! Compute Velocity Gradient
                call comp_vel_gr(rdVelGrd, rStrain, rSpin)

                ! Iterations to satisfy traction continuity
                rBWidth = 0.0001d0

                ! Compute twin normal in current configuration F^-T n
                call inverse_3x3(rFu0, rFuinv)
                rntwc = matmul(transpose(rFuinv), rntw(:, iTwFlg))
                rntwc = rntwc / norm2(rntwc)

                if (iTwEmbStp.eq.1) then
                  ! No twinning part in previous step
                  rFz0h = rFz0

                  ! Guess for initial rmh based on twin shear direction
                  rmtwc = matmul(rFu0, rmtw(:, iTwFlg))
                  rmtwc = rmtwc / norm2(rmtwc)
                  rmh = 1.d-4 * rmtwc
                else
                  ! Elastoplastic part of previous deformation gradient in the twinned region
                  rFz0h = matmul(rFz0, F_tw_inv)

                  ! Retrieve discontinuous mode normal and direction
                  iumatsvst = offstatev(6) + (l - 1) * nstatev(4)
                  rmh(1) = svars(kblock, iumatsvst + 1)
                  rmh(2) = svars(kblock, iumatsvst + 2)
                  rmh(3) = svars(kblock, iumatsvst + 3)
                end if

                g_tr_err = 1.d0  ! some random initial value

                do iter = 1, 100
                  ! Compute m tensor product with n
                  do j = 1, 3
                    do i = 1, 3
                      mtpn(i, j) = rmh(i) * rntwc(j)
                    end do
                  end do

                  ! Compute jump velocity gradient increments
                  rdLu = rdVelGrd - mtpn
                  rdLz = rdVelGrd + (rElmtSize/rBWidth - 1.d0) * mtpn

                  ! Compute new deformation gradient in the untwinned region
                  rdFu = matmul(rdLu, rFu0)
                  rFu1 = rFu0 + rdFu

                  ! Compute Stress (untwinned region) from existing UMAT
                  iumatsvst = offstatev(1) + 1 + (l - 1) * nstatev(1)
                  iumatsven = offstatev(1) + l * nstatev(1)
                  call umat_stress(rStress_U,
     &                             svars(kblock, iumatsvst:iumatsven),
     &                             rGsCoords, rFu0, rFu1,
     &                             props(77:126), 50, nstatev(1), nnode,
     &                             elmt_number, l, mcrd, time, dtimeCur,
     &                             kstep, kinc, .false., rWu, iterErr)

                  if (iterErr.eq.1) then
                    call random_number(rmh)
                    rmh = 1.d-3 * rmh
                    cycle
                  end if

                  ! Compute new deformation gradient in the twinned region
                  rdFz = matmul(rdLz, rFz0)
                  rFz1 = rFz0 + rdFz

                  ! Compute elastoplastic part of twinned deformation gradient
                  rFz1h = matmul(rFz1, F_tw_inv)

                  ! Compute Stress (twinned region) from existing UMAT
                  iumatsvst = offstatev(2) + 1 + (l - 1) * nstatev(1)
                  iumatsven = offstatev(2) + l * nstatev(1)
                  call umat_stress(rStress_Z,
     &                             svars(kblock, iumatsvst:iumatsven),
     &                             rGsCoords, rFz0h, rFz1h,
     &                             props_tw, 50, nstatev(1), nnode,
     &                             elmt_number, l, mcrd, time, dtimeCur,
     &                             kstep, kinc, .false., rWz, iterErr)

                  if (iterErr.eq.1) then
                    call random_number(rmh)
                    rmh = 1.d-3 * rmh
                    cycle
                  end if

                  ! vector to matrix for stress in twinned and untwinned region
                  call vec_to_sym_mat(rStress_U, sigmaU)
                  call vec_to_sym_mat(rStress_Z, sigmaZ)

                ! Compute twin normal in current configuration F^-T n
                ! TODO should this be done with F^U or homogeneous F
                  call inverse_3x3(rFu1, rFuinv)
                  rntwc = matmul(transpose(rFuinv), rntw(:, iTwFlg))
                  rntwc = rntwc / norm2(rntwc)

                  ! Form residual
                  do aidx = 1, 3
                    g_tr(aidx) = 0.d0
                    do i = 1, 3
                      g_tr(aidx) = g_tr(aidx) - rntwc(i) *
     &                            (sigmaZ(i,aidx) - sigmaU(i,aidx))
                    end do
                  end do

                  ! check for residual convergence
                  g_tr_err = norm2(g_tr)
                  if (g_tr_err.lt.1.d-8) then
                    exit
                  end if

                  ! Form residual tangent matrix
                  do j = 1, 3
                    do aidx = 1, 3
                      z_tr(aidx, j) = 0.d0
                      do i = 1, 3
                        do k = 1, 3
                          do p = 1, 3
                            z_tr(aidx, j) = z_tr(aidx, j) +
     &                      rntwc(i) * ((rElmtSize/rBWidth-1.d0)*
     &                      rWz(i,aidx,j,k)*rFz0h(p,k)+rWu(i,aidx,j,k)*
     &                      rFu0(p,k)) * rntwc(p)
                          end do
                        end do
                      end do
                    end do
                  end do

                  ! Compute change in m hat
                  call inverse_3x3(z_tr, z_tr_inv)
                  rdmh = matmul(z_tr_inv, g_tr)

                  ! Compute updated m hat
                  rmh = rmh + rdmh
                end do

                ! End analysis if traction continuity fails
                if (g_tr_err.gt.1.d-8) then
                  CALL XPLB_ABQERR(-3,'Traction Continuity failed! '//
     &                             'iTwEmbStp = %I, err = %R'  ,
     &                              iTwEmbStp, g_tr_err,' ')
                end if

                ! Store discontinuous mode normal and direction
                iumatsvst = offstatev(6) + (l - 1) * nstatev(4)
                svars(kblock, iumatsvst + 1) = rmh(1)
                svars(kblock, iumatsvst + 2) = rmh(2)
                svars(kblock, iumatsvst + 3) = rmh(3)

                ! RSS on twin system
                Smatc = matmul(rFu1, matmul(Smat(:, :, iTwFlg), rFuinv))
                tau_tw = sum(sigmaU * Smatc)
                call compute_twin_evol(props(77:132), 56, tau_tw, theta,
     &                                 dtimeCur, rBWidth, rBLength)

                ! Compute combined discontinuity magnitude and direction
                alpha_h = norm2(rmh)
                rmhu = rmh / alpha_h
                rmtwc = matmul(rFu1, rmtw(:, iTwFlg))
                rmtwc = rmtwc / norm2(rmtwc)

                ! Compute velocity gradient increment projected on discontinuity tensor
                rTdL = 0.d0
                do j = 1, 3
                  do i = 1, 3
                    rTdL = rTdL + rmhu(i) * rntwc(j) * rdVelGrd(i, j)
                  end do
                end do
                
                ! Compute discontinuity magnitude in untwinned and twinned region
                alpha_U = alpha_h / rTdL
                alpha_Z = (rElmtSize/rBWidth - 1.d0) * alpha_U


                ! Save deformation gradients to both untwinned and twinned block
                iumatsvst = offstatev(3) + 1 + (l - 1) * nstatev(2)
                iumatsven = offstatev(3) + l * nstatev(2)
                svars(kblock, iumatsvst : iumatsven) =
     &                                  reshape(rFu1, (/9/))
                iumatsvst = offstatev(4) + 1 + (l - 1) * nstatev(2)
                iumatsven = offstatev(4) + l * nstatev(2)
                svars(kblock, iumatsvst : iumatsven) =
     &                                  reshape(rFz1, (/9/))

                ! Compute Final Stress (untwinned region) from existing UMAT
                iumatsvst = offstatev(1) + 1 + (l - 1) * nstatev(1)
                iumatsven = offstatev(1) + l * nstatev(1)
                call umat_stress(rStress_U,
     &                           svars(kblock, iumatsvst:iumatsven),
     &                           rGsCoords, rFu0, rFu1,
     &                           props(77:126), 50, nstatev(1), nnode,
     &                           elmt_number, l, mcrd, time, dtimeCur,
     &                           kstep, kinc, .true., rWu, iterErr)


                ! Compute Final Stress (twinned region) from existing UMAT
                iumatsvst = offstatev(2) + 1 + (l - 1) * nstatev(1)
                iumatsven = offstatev(2) + l * nstatev(1)
                call umat_stress(rStress_Z,
     &                           svars(kblock, iumatsvst:iumatsven),
     &                           rGsCoords, rFz0h, rFz1h,
     &                           props_tw, 50, nstatev(1), nnode,
     &                           elmt_number, l, mcrd, time, dtimeCur,
     &                           kstep, kinc, .true., rWz, iterErr)

                ! Compute discontinuous mode tensor
                rTdis(1) = rmhu(1) * rntwc(1)
                rTdis(2) = rmhu(2) * rntwc(2)
                rTdis(3) = rmhu(3) * rntwc(3)
                rTdis(4) = 0.5d0 * (rmhu(2) * rntwc(3) +
     &                              rmhu(3) * rntwc(2))
                rTdis(5) = 0.5d0 * (rmhu(1) * rntwc(3) +
     &                              rmhu(3) * rntwc(1))
                rTdis(6) = 0.5d0 * (rmhu(1) * rntwc(2) +
     &                              rmhu(2) * rntwc(1))
                rTdis(4) = rmhu(2) * rntwc(3) + rmhu(3) * rntwc(2)
                rTdis(5) = rmhu(1) * rntwc(3) + rmhu(3) * rntwc(1)
                rTdis(6) = rmhu(1) * rntwc(2) + rmhu(2) * rntwc(1)

                ! Compute B-bar in twinned and untwinned region
                rBmatU = rBmat
                rBmatZ = rBmat
                do aidx = 1, mcrd*nnode
                  rTBmat = 0.d0
                  do i = 1, 6
                    rTBmat = rTBmat + rTdis(i) * rBmat(i, aidx)
                  end do

                  rBmatU(:,aidx) = rBmatU(:, aidx) -
     &                             alpha_U * rTdis(:) * rTBmat
                  rBmatZ(:, aidx) = rBmatZ(:, aidx) +
     &                             alpha_Z * rTdis(:) * rTBmat
                end do


C               Compute linear + quadratic bulk viscosity terms
C               Needed for damping the “ringing” in the highest element frequency
C               (Abaqus->Theory->Procedures->Nonlinear dynamics->Explicit dynamic analysis)
                !! Tao,The quadratic bulk viscosity is applied
                !! only if the volumetric strain rate is compressive.
                if (dtimeCur.ne.0.d0) then
                  rBVpress = -rRho *
     $                        rbv1*rCd*rElmtSize*(rEpVolInc/dtimeCur)
                  if (rEpVolInc.lt.0.d0) then
                    rBVpress = rBVpress
     $                      -rRho*(rbv2*rElmtSize*rEpVolInc/dtimeCur)**2
                  end if
                else
                  rBVpress = 0.d0
                end if

C               Compute the internal force vector
                do n = 1, nnode
                  do i = 1, mcrd
                    do j = 1, 6
                      if (j.le.3) then    ! add bulk viscosity
                        rSigmaU = rStress_U(j) - rBVpress
                        rSigmaZ = rStress_Z(j) - rBVpress
                      else
                        rSigmaU = rStress_U(j)
                        rSigmaZ = rStress_Z(j)
                      end if
                      rhs( kblock , (n-1)*mcrd + i ) =
     &                rhs( kblock , (n-1)*mcrd + i )
     &                + ((rBWidth/rElmtSize)*rBmatZ(j,(n-1)*mcrd+i)*
     &                rSigmaZ + (1.d0-rBWidth/rElmtSize)*
     &                rBmatU(j,(n-1)*mcrd+i)*rSigmaU)*rdOmega
                    end do
                  end do
                end do

              end do

              ! Twin is embedded in this step, unset the flag
              if (iTwEmbStp.eq.1) then
                iTwEmbStp = 0
                svars(kblock, offstatev(5) + 2) = iTwEmbStp
              end if


            end if

            ! Compute stable time step increment
            dtimeStable(kblock) = 0.7d0 * rElmtSize / rCd
          end do
        ENDIF
      END  ! Subroutine VUEL

      function is_gb_elem(el_num, gb_elems, num_gb_elems)
        implicit none
        logical is_gb_elem
        integer el_num, i, num_gb_elems
        integer gb_elems(num_gb_elems)

        is_gb_elem = .false.
        do i = 1, num_gb_elems
          if (el_num.eq.gb_elems(i)) then
            is_gb_elem = .true.
            exit
          end if
        end do
        return
      end function

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)
        implicit none
#include <SMAAspUserSubroutines.hdr>
C       Variables passed from interface
        integer lOp, niArray, nrArray
        integer i_Array(niArray)
        real*8 r_Array(nrArray)

C       Contents of i_Array
        integer i_int_nTotalNodes, i_int_nTotalElements,
     1          i_int_kStep, i_int_kInc, i_int_iStatus,
     2          i_int_lWriteRestart
        parameter( i_int_nTotalNodes     = 1,
     *            i_int_nTotalElements  = 2,
     *            i_int_kStep           = 3,
     *            i_int_kInc            = 4,
     *            i_int_iStatus         = 5,
     *            i_int_lWriteRestart   = 6  )

C       Possible values for the lOp argument
        integer j_int_StartAnalysis, j_int_StartStep,
     1          j_int_SetupIncrement, j_int_StartIncrement,
     2          j_int_EndIncrement, j_int_EndStep,
     3          j_int_EndAnalysis
        parameter(j_int_StartAnalysis    = 0,
     *            j_int_StartStep        = 1,
     *            j_int_SetupIncrement   = 2,
     *            j_int_StartIncrement   = 3,
     *            j_int_EndIncrement     = 4,
     *            j_int_EndStep          = 5,
     *            j_int_EndAnalysis      = 6 )

C       Possible values for i_Array(i_int_iStatus)
        integer j_int_Continue, j_int_TerminateStep,
     1          j_int_TerminateAnalysis
        parameter(j_int_Continue          = 0,
     *            j_int_TerminateStep     = 1,
     *            j_int_TerminateAnalysis = 2)

C       Contents of r_Array
        integer i_flt_TotalTime, i_flt_StepTime,
     1          i_flt_dTime
        parameter(i_flt_TotalTime   = 1,
     1            i_flt_StepTime    = 2,
     2            i_flt_dTime       = 3 )

        integer kStep
        integer kInc
        integer myThreadID

        kStep = i_Array(i_int_kStep)
        kInc  = i_Array(i_int_kInc)
        myThreadID = get_thread_id()

        if (lop.eq.j_int_StartAnalysis) then
          if (myThreadID.eq.0) then
            call read_gb_elements()
          end if
        else if (lop.eq.j_int_EndAnalysis) then
          if (myThreadID.eq.0) then
            call SMAIntArrayDelete(1)
          end if
        end if

      end subroutine

      subroutine read_gb_elements()
        implicit none
#include <SMAAspUserSubroutines.hdr>
        character(len = 50) :: gbfilename
        character (len = 1) :: pathseperator = '/'
        character (len = 256) :: outdir, gbfilepath
        integer :: lenoutdir, i, io, num_gb_elems
        integer gb_elems(*)
        pointer (ptr_gb_elems, gb_elems)

        gbfilename = 'gb_elems_1000.txt'
        call vgetoutdir(outdir, lenoutdir)
        gbfilepath = trim(outdir)//pathseperator//trim(gbfilename)

        num_gb_elems = 0
        open(120, file = gbfilepath, status = 'old', iostat = io)
        do
          read(120, *, iostat = io)
          if (io /= 0) exit
          num_gb_elems = num_gb_elems + 1
        end do
        close(120)

        ptr_gb_elems = SMAIntArrayCreate(1, num_gb_elems, -1)

        open(121, file = gbfilepath, status = 'old',action = 'read')
        do i = 1, num_gb_elems
          read(121, *) gb_elems(i)
        end do
        close(121)
      end subroutine

      include 'GSPT_3D_8N.f'
      include 'SHP_FUNC_3D_8N.f'
      include 'prep_umat.f'
      include 'gdefgrad_3D.f'
      include 'fcc_umat_lbl34_test.f'
      include 'Gstrain_3D.f'
