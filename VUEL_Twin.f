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
C          Akhilesh Pedgaonkar, TCMMG, UW-Madison

C     This code is adapted from the 3D shear band code developed by
C     Jin et. al., 2019, FEAD 155, 11.
C     This is an extension of 2D twin paper 
C     Jin et. al., 2019, JMPS 133, 103723 to three dimensions.


C     Voigt notation:
C          Different from FEAP or ABAQUS convection,

C          (1,1) -> 1
C          (2,2) -> 2
C          (3,3) -> 3
C          (2,3) or (3,2) -> 4
C          (1,3) or (3,1) -> 5
C          (1,2) or (2,1) -> 6


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
        integer num_elems        ! Number of vuel elements
        integer num_nodes        ! Number of nodes
        integer twinID           ! Twin ID
        integer tw_nxt_elmt      ! Forward element for twin propagation along length

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
        real*8 rStressCtZ(6)     ! Cauchy Stress at Elmt Centroid
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
        real*8 rStress(6, 8)         ! Cauchy Stress
        real*8 rBmat(6,mcrd*nnode) ! current B (or B-bar) operator
        real*8 rBskw(3,mcrd*nnode) ! current skew-symmetric matrix operator
        real*8 Smat(3, 3, 18)      ! Twinning Schmid Tensor
        real*8 F_tw(3, 3, 18)      ! Twinning Deformation gradient
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
        real*8 rBWidthOff          ! Band width offset
        real*8 F_tw_inv(3, 3)      ! Inverse of twinned deformation gradient
        real*8 rStress_U(6, 8)     ! Stress in untwinned region
        real*8 rStress_Z(6, 8)     ! Stress in twinned region
        real*8 rWu(3, 3, 3, 3)     ! dT/dF^U Material Jacobian for untwinned region
        real*8 rWz(3, 3, 3, 3)     ! dT/dF^Z Material Jacobian for twinned region
        real*8 rntw(3, 18)         ! Twin normal
        real*8 rmtw(3, 18)         ! Twin shear direction
        real*8 props_tw(nprops - 13)  ! Properties for twinned region since euler angles have changed
        ! TODO: Remove hard coded 50
        real*8 pi                  ! PI
        real*8 g_tr(3)             ! residual for traction continuity
        real*8 z_tr(3, 3)          ! tangent matrix for traction continuity
        real*8 z_tr_inv(3, 3)      ! inverse of tangent matrix for traction continuity
        real*8 sigmaZ(3, 3)        ! Tensor Stress in twinned region
        real*8 sigmaU(3, 3)        ! Tensor Stress in untwinned region
        real*8 rdmh(3)             ! change in jump + direction
        real*8 rmh(3)              ! jump + direction
        real*8 rmh0(3)             ! initial jump + direction
        real*8 rmhdot(3)           ! Rate of change of jump + direction
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
        real*8 Fe_m(3, 3)          ! Mean current deformation gradient
        real*8 Fe_inv(3, 3)       ! Inverse of mean current deformation gradient
        real*8 rElmtCt(3)          ! Element centroid coordinates
        real*8 tw_len_max          ! Maximum of twin length in an element
        real*8 rBLengthMax         ! Current maximum twin length
        real*8 Fe(3, 3, 8)         ! Current Elastic part of deformation gradient
        real*8 FeZ(3, 3)           ! Current Elastic part of dfgrd in twinned region

        integer gb_elems(*)        ! Grain boundary element array
        pointer (ptr_gb_elems, gb_elems) ! pointer for gb_elems required for SMAArrayAccess
        logical isGbElem           ! logical determining gb elem
        logical, external :: is_gb_elem
        integer elmt_conn(8, *)    ! Element connectivity
        pointer (ptr_elmt_conn, elmt_conn)   ! Pointer to element connectivity for SMA Access
        integer elmt_neighbour(6, *) ! List of neighbouring element
        pointer (ptr_elmt_neighbour, elmt_neighbour) ! pointer for elmt_neighbour required for SMA Access
        real*8 nodal_coords(3, *)    ! Nodal coordinates
        pointer (ptr_nodal_coords, nodal_coords)  ! Pointer to nodal coordinates for SMA Access
        integer size_list(*)    ! list storing num_elems, num_nodes, etc
        pointer (ptr_size_list, size_list)  ! pointer to size list for SMA Access
        integer elmt_stat(3, *)
        pointer (ptr_elmt_stat, elmt_stat)
        integer tw_flg(*)
        pointer (ptr_tw_flg, tw_flg)
        real*8 tw_nuc(3, *), tw_nor(3, *), tw_dir(3, *)
        pointer (ptr_tw_nuc, tw_nuc)
        pointer (ptr_tw_nor, tw_nor)
        pointer (ptr_tw_dir, tw_dir)
        real*8 tw_geom(3, *)
        pointer (ptr_tw_geom, tw_geom)
        integer tw_dec(3, *)
        pointer (ptr_tw_dec, tw_dec)
        real*8 elmt_stress(12, *)
        pointer (ptr_elmt_stress, elmt_stress)
        real*8 nuc_prob(*)
        pointer (ptr_nuc_prob, nuc_prob)
        real*8 disp(3, *)
        pointer (ptr_disp, disp)
        real*8 nuc_thr(*)
        pointer (ptr_nuc_thr, nuc_thr)
        real*8 elmt_SDV(14, *)
        pointer (ptr_elmt_SDV, elmt_SDV)

        real*8 l_dot, w_dot, b_dot, rBLatWidth
        real*8 rNRFac, rNRTol
        integer nAltInit, max_twins, num_twins
        real*8 rmhMag, bhratio, hbratio
        real*8 PtauMax, rLmin
        integer stmprop, nmprops, endmprop, eulidx(3)
        real*8 R_tw(3, 3)      ! Rotation due to twinning
        integer num_twins_gr, max_twins_gr
        real*8 P_0
        logical lNRerr
        real*8 outSDVM(7)
        real*8 outSDVMZ(7)
        real*8 temp

        iVelmType = jprops(1)  ! 3D element type:  1 --> 3D Brick (8-node) Elmt


        iVforType = jprops(2)  ! Formulation type: 0 --> displacement based
                               !                   1 --> displacement based with localization behavior
                               !                   2 --> B-bar incompressible
                               !                   3 --> B-bar incompressible with localization behavior

        rEpsilon  = props(5)  ! Pressure mode stabilization parameter; see Ortiz et al (1987), Eq. (21).

        iVforType = 0
        !call sleep(20)

        ! Set offset and nprops for umat properties
        stmprop = 14
        eulidx(1) = 31
        eulidx(2) = 32
        eulidx(3) = 33
        nmprops = nprops - stmprop + 1
        endmprop = nprops

        ! Number of state variables array and offsets array
        nstatev(1) = (nsvars - 9*8 - 9*8 - 9*8 - 1) / 16  ! svars for umat
        nstatev(2) = 9    ! svars for deformation gradient
        nstatev(3) = 9    ! svars for storing discontinuity vectors
        nstatev(4) = 1    ! svars for storing band width offset
        offstatev(1) = 0
        offstatev(2) = 8 * nstatev(1)                 ! umat svars for 8 gs pts
        offstatev(3) = offstatev(2) + 8 * nstatev(1)  ! umat svars for 8 gs pts for twinned region
        offstatev(4) = offstatev(3) + 8 * nstatev(2)  ! def grd svars for 8 gs pts
        offstatev(5) = offstatev(4) + 8 * nstatev(2)  ! def grd svars for 8 gs pts for twinned region
        offstatev(6) = offstatev(5) + 8 * nstatev(3)  ! discontinuity vector for 8 gs pts
        offstatev(7) = offstatev(6) + nstatev(4)
        ! 148*8 + 148*8 + 9*8 + 9*8 + 9*8 + 1 = 2585
        ! Set "variables=2585" for user elements in the input file

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

        call GSPT_3D_8N(iTIntPt, rXiEtaPt, rWeightPt)

C       Natural coordinates of element centroid
        rXiEtaCt(1) = 0.d0
        rXiEtaCt(2) = 0.d0
        rXiEtaCt(3) = 0.d0

C       Compute equivalent Schmid Tensor and Deformation gradient for Deformation Twinning
        call comp_twin_dfgrd(del, props(stmprop:endmprop), nmprops, pi,
     &                       Smat, F_tw, rntw, rmtw, eulidx)

        ! Twin details
        max_twins = 1
        max_twins_gr = 1

        ! Access SMA arrays created in vexternaldb
        if (kstep.gt.0) then
          ptr_size_list = SMAIntArrayAccess(1)
          ptr_gb_elems = SMAIntArrayAccess(2)
          ptr_elmt_conn = SMAIntArrayAccess(3)
          ptr_elmt_neighbour = SMAIntArrayAccess(4)
          ptr_nodal_coords = SMAFloatArrayAccess(5)
          ptr_elmt_stat = SMAIntArrayAccess(6)
          ptr_tw_flg = SMAIntArrayAccess(7)
          ptr_tw_nuc = SMAFloatArrayAccess(8)
          ptr_tw_nor = SMAFloatArrayAccess(9)
          ptr_tw_dir = SMAFloatArrayAccess(10)
          ptr_tw_geom = SMAFloatArrayAccess(11)
          ptr_tw_dec = SMAIntArrayAccess(12)
          ptr_elmt_stress = SMAFloatArrayAccess(13)
          ptr_nuc_prob = SMAFloatArrayAccess(14)
          ptr_disp = SMAFloatArrayAccess(15)
          ptr_nuc_thr = SMAFloatArrayAccess(16)
          ptr_elmt_SDV = SMAFloatArrayAccess(17)

          num_gb_elems = size_list(1)
          num_elems = size_list(2)
          num_nodes = size_list(3)
        end if

        ! Compute dilitational wave speed
        rRho = props(1)
        c11 = props(2)
        rCd = sqrt(c11 / rRho)

        ! Linear bulk viscosity coefficient
        if (props(3).le.0.d0) then
          rbv1    = 0.06d0
        else
          rbv1    = props(3)
        end if

        ! Quadratic bulk viscosity coefficient
        if (props(4).le.0.d0) then
          rbv2    = 1.2d0
        else
          rbv2    = props(4)
        end if

        ! Newton raphson correction factor and tolerance for traction continuity
        rNRFac = 0.25d0
        rNRTol = 1.d-6

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
          ! Loop over elements in the block
          do kblock = 1,nblock

            elmt_number = jElem(kblock)

            ! Set element type number in status to distinguish grains
            if ((kstep.gt.0).and.(kinc.eq.0)) then
              elmt_stat(3, elmt_number) = jtype
            end if

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

            ! Update global nodal coordinates
            if (kstep.gt.0) then
              ! do n = 1, nnode
              !   nodal_coords(:, elmt_conn(n, elmt_number))=rCoords(:, n)
              ! end do

              do n = 1, nnode
                disp(1, elmt_conn(n, elmt_number)) =
     &                                         u(kblock, (n-1)*mcrd + 1)
                disp(2, elmt_conn(n, elmt_number)) =
     &                                         u(kblock, (n-1)*mcrd + 2)
                disp(3, elmt_conn(n, elmt_number)) =
     &                                         u(kblock, (n-1)*mcrd + 3)
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

            ! GB Element detection and threshold probability
            P_0 = 0.85d0
            if (kstep.gt.0) then
              isGbElem = .false.
              do i = 1, num_gb_elems
                if (elmt_number.eq.gb_elems(i)) then
                  isGbElem = .true.
                  P_0 = nuc_thr(i)
                  exit
                end if
              end do
            end if

            ! Evaluate twinID and twinning system flag
            twinID = 0
            iTwFlg = 0
            iTwEmbStp = 0
            if (kstep.gt.0) then
              ! Check element status for twinID
              twinID = elmt_stat(1, elmt_number)
              iTwEmbStp = elmt_stat(2, elmt_number)
              if (twinID.ne.0) then
                iTwFlg = tw_flg(twinID)
              end if
            end if

            PtauMax = 0.d0

            if((twinID.eq.0).or.(iTwEmbStp.ge.2)) then   ! No twinning detected in previous increment
              !Fe_m = 0.d0                                ! or twinning detected by tracking
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
                if (kinc.lt.1) then
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
                call umat_stress(rStress(:, l),
     &                           svars(kblock, iumatsvst:iumatsven),
     &                           rGsCoords, rDfgrd0, rDfgrd1,
     &                           props(stmprop:endmprop), nmprops, 
     &                           nstatev(1), nnode, elmt_number, l,
     &                           mcrd, time, dtimeCur, kstep, kinc,
     &                           .true., rWu, iterErr, Fe(:, :, l))

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
                if (isGbElem.and.(iTwEmbStp.eq.0)) then
                  call check_twin_nucleation(rStress(:, l), Smat,
     &                                       Fe(:,:,l), iTwFlg, PtauMax,
     &                                       P_0)
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
                        rSigma = rStress(j, l) - rBVpress
                      else
                        rSigma = rStress(j, l)
                      end if
                      rhs( kblock , (n-1)*mcrd + i ) =
     $                rhs( kblock , (n-1)*mcrd + i )
     $                + rBmat( j , (n-1)*mcrd + i ) * rSigma * rdOmega
                    END DO
                  END DO
                END DO

              end do

              if (kstep.gt.0) then
                ! Initialize element average stress and SDV
                elmt_stress(:, elmt_number) = 0.d0
                elmt_SDV(:, elmt_number) = 0.d0

                ! Write avg global element stress
                rStressCt = 0.d0
                do l = 1, iTIntPt
                  rStressCt = rStressCt + rStress(:, l)
                end do

                rStressCt = rStressCt / iTIntPt
                elmt_stress(1:6, elmt_number) = rStressCt

                ! Write avg state variables
                outSDVM = 0.d0
                do l = 1, iTIntPt
                  iumatsvst = offstatev(1) + 1 + (l - 1) * nstatev(1)

                  ! Temperature
                  outSDVM(1) = outSDVM(1) + svars(kblock, iumatsvst)

                  ! Euler angles
                  outSDVM(2:4) = outSDVM(2:4) +
     &                             svars(kblock, iumatsvst + 1:3)

                  ! Cumulative dislocation density
                  temp = 0.d0
                  do i = 1, 18
                    temp = temp + svars(kblock, iumatsvst + 36 + i)
                  end do
                  outSDVM(5) = outSDVM(5) + temp / 18

                  ! Effective temperature
                  outSDVM(6) = outSDVM(6) +
     &                             svars(kblock, iumatsvst + 55)

                  ! Angular Misorientation
                  outSDVM(7) = outSDVM(7) +
     &                             svars(kblock, iumatsvst + 56)

                end do

                outSDVM = outSDVM / iTIntPt
                elmt_SDV(1:7, elmt_number) = outSDVM
              end if

              ! Save max nucleation probability
              if ((kstep.gt.0).and.isGbElem) then
                nuc_prob(elmt_number) = PtauMax
              end if

              ! If twin nucleated
              if ((iTwFlg.ne.0).and.(iTwEmbStp.eq.0)) then
                ! Check num_twins and num_twins per grain
                num_twins = size_list(4)
                if (num_twins.lt.max_twins) then
                  call MutexLock(2)
                  num_twins = size_list(4)
                  num_twins_gr = 0
                  do n = 1, num_twins
                    if (jtype.eq.elmt_stat(3, tw_dec(3, n))) then
                      num_twins_gr = num_twins_gr + 1
                    end if
                  end do
                  if (num_twins_gr.lt.max_twins_gr) then
                    size_list(4) = size_list(4) + 1
                    twinID = size_list(4)
                    elmt_stat(1, elmt_number) = twinID

                    ! Store twin type
                    tw_flg(twinID) = iTwFlg

                    ! Compute element centroid coordinates and store
                    rElmtCt = 0.d0
                    do n = 1, iTIntPt
                      rElmtCt = rElmtCt + rCoords(1:3, n)
                    end do
                    rElmtCt = 1.d0 / iTIntPt * rElmtCt
                    tw_nuc(1:3, twinID) = rElmtCt

                    ! Store global twin normal and direction (reference config.)
                    tw_nor(1:3, twinID) = rntw(:, iTwFlg)
                    tw_dir(1:3, twinID) = rmtw(:, iTwFlg)

                    ! Initialize global twin geometry parameters (length, width, thickness)
                    tw_geom(1, twinID) = 0.d0
                    tw_geom(2, twinID) = 0.d0
                    tw_geom(3, twinID) = 1.d-4

                    ! Store this element number as twin decider for all geometry paramters
                    tw_dec(:, twinID) = elmt_number
                  else
                    ! We have maxed out twins, shut of the flag
                    iTwFlg = 0
                  end if
                  call MutexUnlock(2)
                end if
              end if

              ! If twin nucleated or propagated
              if (iTwFlg.ne.0) then
                do l = 1, 8
                  ! Compute current twin normal (F^-T n / || F^-T n ||)
                  call inverse_3x3(Fe(:, :, l), Fe_inv)
                  rntwc = matmul(transpose(Fe_inv), rntw(:, iTwFlg))
                  rntwc = rntwc / norm2(rntwc)

                  ! Compute current twin direction F m / || F m ||
                  rmtwc = matmul(Fe(:, :, l), rmtw(:, iTwFlg))
                  rmtwc = rmtwc / norm2(rmtwc)

                  ! Store current normal and direction
                  iumatsvst = offstatev(5) + (l - 1) * nstatev(3)
                  svars(kblock, iumatsvst + 1) = rntwc(1)
                  svars(kblock, iumatsvst + 2) = rntwc(2)
                  svars(kblock, iumatsvst + 3) = rntwc(3)
                  svars(kblock, iumatsvst + 4) = rmtwc(1)
                  svars(kblock, iumatsvst + 5) = rmtwc(2)
                  svars(kblock, iumatsvst + 6) = rmtwc(3)
                end do

                ! Twin band width offset
                if (iTwEmbStp.eq.3) then
                  rBWidthOff = tw_geom(3, twinID)
                else
                  rBWidthOff = 0.d0
                end if
                iumatsvst = offstatev(6) + 1
                svars(kblock, iumatsvst) = rBWidthOff

                ! Set twin embedding flag to 1
                iTwEmbStp = 1
                elmt_stat(2, elmt_number) = iTwEmbStp
              end if

            else  ! Twinning detected in the previous increment

              ! compute inverse of twinning deformation gradient
              call inverse_3x3(F_tw(:, :, iTwFlg), F_tw_inv)

              ! compute properties for twinned region due to rotation of crystal basis
              call compute_twin_props(rntw(:, iTwFlg),
     &                                props(stmprop:endmprop), props_tw,
     &                                nmprops, del, pi, eulidx)

              ! Form rotation tensor due to twinning reflection
              do j = 1, 3
                do i = 1, 3
                  R_tw(i, j) = 2 * rntw(i, iTwFlg) * rntw(j, iTwFlg) -
     &                         del(i, j)
                end do
              end do

              ! Retrieve twin geometry
              rBLength = tw_geom(1, twinID)
              rBLatWidth = tw_geom(2, twinID)
              iumatsvst = offstatev(6) + 1
              rBWidthOff = svars(kblock, iumatsvst)
              rBWidth = tw_geom(3, twinID) - rBWidthOff
              if (rBWidth.lt.1.d-4) rBWidth = 1.d-4

              ! b/h ratio and handle b >= h
              bhratio = rBWidth / rElmtSize
              !if (bhratio.ge.0.9999d0) bhratio = 0.9999d0
              hbratio = 1.d0 / bhratio

              print*, "length = ", rBLength
              print*, "width = ", rBLatWidth
              print*, "thickness = ", rBWidth
              print*, "rBWidthOff = ", rBWidthOff

              ! DEBUGGING
              print*, "time = ", time(2)
              print*, "twinID = ", twinID
              print*, "elmt no = ", elmt_number
              print*, "iTwFlg = ", iTwFlg
              print*, "iTwEmbStp = ", iTwEmbStp
              print*, "rntw = ", rntw(:, iTwFlg)
              print*, "rmtw = ", rmtw(:, iTwFlg)
              print*, "twin_decider = ", tw_dec(:, twinID)
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

                ! Retrieve previous untwinned and twinned deformation gradient
                iumatsvst = offstatev(3) + 1 + (l - 1) * nstatev(2)
                iumatsven = offstatev(3) + l * nstatev(2)
                rFu0 =
     &            reshape(svars(kblock, iumatsvst : iumatsven), (/3,3/))
                iumatsvst = offstatev(4) + 1 + (l - 1) * nstatev(2)
                iumatsven = offstatev(4) + l * nstatev(2)
                rFz0 =
     &            reshape(svars(kblock, iumatsvst : iumatsven), (/3,3/))

                ! Retrieve discontinuous mode normal
                iumatsvst = offstatev(5) + (l - 1) * nstatev(3)
                rntwc(1) = svars(kblock, iumatsvst + 1)
                rntwc(2) = svars(kblock, iumatsvst + 2)
                rntwc(3) = svars(kblock, iumatsvst + 3)

                ! If fully twinned don't perform traction continuity
                if (iTwEmbStp.ne.0) then
                  ! No twinning part in previous step
                  rFz0h = rFz0
                  rFz0 = matmul(rFz0, F_tw(:, :, iTwFlg))
                  iumatsvst = offstatev(2) + 1 + (l - 1) * nstatev(1)
                  iumatsven = offstatev(2) + l * nstatev(1)
                  call rotate_twin_svars(svars(kblock,
     &                                   iumatsvst:iumatsven),
     &                                   nstatev(1), R_tw)
                else
                  ! Elastoplastic part of previous deformation gradient in the twinned region
                  rFz0h = matmul(rFz0, F_tw_inv)
                end if

                ! The element has completely twinned
                if (bhratio.gt.0.9999d0) then
                  ! Compute new deformation gradient and save
                  rFz1 = matmul(del + rdVelGrd, rFz0)
                  iumatsvst = offstatev(4) + 1 + (l - 1) * nstatev(2)
                  iumatsven = offstatev(4) + l * nstatev(2)
                  svars(kblock, iumatsvst : iumatsven) =
     &                                              reshape(rFz1, (/9/))

                  ! Compute elastoplastic part of twinned deformation gradient
                  rFz1h = matmul(rFz1, F_tw_inv)

                  ! Compute Stress (twinned region) from existing UMAT
                  iumatsvst = offstatev(2) + 1 + (l - 1) * nstatev(1)
                  iumatsven = offstatev(2) + l * nstatev(1)
                  call umat_stress(rStress_Z(:, l),
     &                             svars(kblock, iumatsvst:iumatsven),
     &                             rGsCoords, rFz0h, rFz1h,
     &                             props_tw, nmprops, nstatev(1), nnode,
     &                             elmt_number, l, mcrd, time, dtimeCur,
     &                             kstep, kinc, .true., rWz, iterErr,
     &                             FeZ)

                  if (iterErr.eq.1) then
                    print*, "kinc = ", kinc
                    call xplb_abqerr(-3,'UMAT Failed to Converge in '//
     &                             'fully twinned region!', 0,0.d0,' ')
                  end if

                  theta = svars(kblock, iumatsvst)    ! Twinned Temperature

                  ! Compute bulk viscosity terms
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

C                 Compute the internal force vector
                  do n = 1, nnode
                    do i = 1, mcrd
                      do j = 1, 6
                        if (j.le.3) then    ! add bulk viscosity
                          rSigmaZ = rStress_Z(j, l) - rBVpress
                        else
                          rSigmaZ = rStress_Z(j, l)
                        end if
                        rhs( kblock , (n-1)*mcrd + i ) =
     &                      rhs( kblock , (n-1)*mcrd + i ) +
     &                            rBmat(j,(n-1)*mcrd+i)*rSigmaZ*rdOmega
                      end do
                    end do
                  end do
                  ! Compute twin length, width and thickness evolution
                  call inverse_3x3(FeZ, Fe_inv)
                  Smatc =
     &                  matmul(FeZ, matmul(Smat(:, :, iTwFlg), Fe_inv))
                  tau_tw = sum(sigmaZ * Smatc)
                  call compute_twin_velocities(props, nprops, tau_tw,
     &                                       theta, l_dot, w_dot, b_dot)

                  ! Update length, width and thickness only if this element is a decider
                  if (tw_dec(1, twinID).eq.elmt_number) then
                    tw_geom(1, twinID) = tw_geom(1, twinID) +
     &                                 1.d0 / iTIntPt * l_dot * dtimeCur
                  end if
                  if (tw_dec(2, twinID).eq.elmt_number) then
                    tw_geom(2, twinID) = tw_geom(2, twinID) +
     &                                 1.d0 / iTIntPt * w_dot * dtimeCur
                  end if
                  if (tw_dec(3, twinID).eq.elmt_number) then
                    tw_geom(3, twinID) = tw_geom(3, twinID) +
     &                                 1.d0 / iTIntPt * b_dot * dtimeCur
                  end if
                  cycle   ! Loop, don't continue further
                end if

                ! DEBUGGING
                print*, "l = ", l

                ! Initialization of variables and indices before traction continuity
                nAltInit = 0
                iter = 0
                rmh0 = 0.d0
                rdmh = 0.d0
                g_tr_err = 1.d3
                iumatsvst = offstatev(1) + 1 + (l - 1) * nstatev(1)
                iumatsven = offstatev(1) + l * nstatev(1)
                iumatsvst2 = offstatev(2) + 1 + (l - 1) * nstatev(1)
                iumatsven2 = offstatev(2) + l * nstatev(1)

                ! Initial guess loop for traction continuity
                do while(nAltInit.le.6)
                  call init_guess(nAltInit, rmh, svars, nblock,
     &                            nsvars, kblock, l, offstatev,
     &                            nstatev, dtimeCur, rdVelGrd, rntwc,
     &                            hbratio, iTwFlg, iTwEmbStp)

                  print*, "rmh = ", rmh
                  ! Traction continuity Newton iterations
                  call traction_continuityNR(rntwc, rdVelGrd, hbratio,
     &                 rFu0, rFz0, rFz0h, del, nstatev(1),
     &                 svars(kblock, iumatsvst:iumatsven),
     &                 svars(kblock, iumatsvst2:iumatsven2), nmprops,
     &                 props(stmprop:endmprop), props_tw, Fe(:, :, l),
     &                 nnode, elmt_number, l, mcrd, time, dtimeCur,
     &                 kstep, kinc, F_tw_inv, rmh, rFu1, rFz1, rFz1h,
     &                 lNRerr)

                  if (lNRerr.eq..false.) then
                    exit
                  end if

                  nAltInit = nAltInit + 1
                  print*, "nAltInit = ", nAltInit
                end do

                ! End analysis if traction continuity fails
                if (lNRerr) then
                  ! DEBUGGING
                  print*, "twinID = ", twinID
                  print*, "elmt_no = ", elmt_number
                  print*, "itwflg = ", iTwFlg
                  print*, "rntwc = ", rntwc
                  print*, "rmtwc = ", rmtwc
                  print*, "rmh = ", rmh
                  print*, "rFu0 = ", rFu0
                  print*, "rdVelGrd = ", rdVelGrd
                  print*, "Ftwinv = ", F_tw_inv
                  print*, "l = ", l
                  iumatsvst = offstatev(1) + 1 + (l - 1) * nstatev(1)
                  iumatsven = offstatev(1) + l * nstatev(1)
                  print*, "svars = ", svars(kblock, iumatsvst:iumatsven)
                  print*, "rFu1 = ", rFu1
                  print*, "rFz1 = ", rFz1
                  print*, "rFz1h = ", rFz1h
                  print*, "kinc = ", kinc

                  ! TODO: put mutexlock?
                  call write_vuel_results(kinc, time(2))
                  CALL XPLB_ABQERR(-3,'Traction Continuity failed! '//
     &                             'iTwEmbStp = %I, err = %R'  ,
     &                              iTwEmbStp, g_tr_err,' ')
                end if

                ! Store discontinuous mode normal and direction
                call inverse_3x3(Fe(:, :, l), Fe_inv)
                rntwc = matmul(transpose(Fe_inv), rntw(:, iTwFlg))
                rntwc = rntwc / norm2(rntwc)

                rmtwc = matmul(Fe(:, :, l), rmtw(:, iTwFlg))
                rmtwc = rmtwc / norm2(rmtwc)    ! TODO: Remove, not used anymore

                iumatsvst = offstatev(5) + (l - 1) * nstatev(3)
                svars(kblock, iumatsvst + 1) = rntwc(1)
                svars(kblock, iumatsvst + 2) = rntwc(2)
                svars(kblock, iumatsvst + 3) = rntwc(3)
                if (iTwEmbStp) then ! Set rmh and rmhdot to 0 since twin
                  ! embedding has a very different rmh from rest of analysis
                  svars(kblock, iumatsvst + 4) = 0.d0
                  svars(kblock, iumatsvst + 5) = 0.d0
                  svars(kblock, iumatsvst + 6) = 0.d0

                  svars(kblock, iumatsvst + 7) = 0.d0
                  svars(kblock, iumatsvst + 8) = 0.d0
                  svars(kblock, iumatsvst + 9) = 0.d0
                else ! Store rmh and rmhdot
                  rmh0(1) = svars(kblock, iumatsvst + 4)
                  rmh0(2) = svars(kblock, iumatsvst + 5)
                  rmh0(3) = svars(kblock, iumatsvst + 6)
                  rmhdot = (rmh - rmh0) / dtimeCur

                  svars(kblock, iumatsvst + 4) = rmh(1)
                  svars(kblock, iumatsvst + 5) = rmh(2)
                  svars(kblock, iumatsvst + 6) = rmh(3)

                  svars(kblock, iumatsvst + 7) = rmhdot(1)
                  svars(kblock, iumatsvst + 8) = rmhdot(2)
                  svars(kblock, iumatsvst + 9) = rmhdot(3)
                end if


                ! Save deformation gradients to both untwinned and twinned block
                iumatsvst = offstatev(3) + 1 + (l - 1) * nstatev(2)
                iumatsven = offstatev(3) + l * nstatev(2)
                svars(kblock, iumatsvst : iumatsven) =
     &                                  reshape(rFu1, (/9/))
                iumatsvst = offstatev(4) + 1 + (l - 1) * nstatev(2)
                iumatsven = offstatev(4) + l * nstatev(2)
                svars(kblock, iumatsvst : iumatsven) =
     &                                  reshape(rFz1, (/9/))

                ! Compute combined discontinuity magnitude and direction
                alpha_h = norm2(rmh)
                rmhu = rmh / alpha_h

                ! Compute velocity gradient increment projected on discontinuity tensor
                rTdL = 0.d0
                do j = 1, 3
                  do i = 1, 3
                    rTdL = rTdL + rmhu(i) * rntwc(j) * rdVelGrd(i, j)
                  end do
                end do
                
                ! Compute discontinuity magnitude in untwinned and twinned region
                alpha_U = alpha_h / rTdL
                alpha_Z = (hbratio - 1.d0) * alpha_U

                ! Compute discontinuous mode tensor
                rTdis(1) = rmhu(1) * rntwc(1)
                rTdis(2) = rmhu(2) * rntwc(2)
                rTdis(3) = rmhu(3) * rntwc(3)
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

                theta = 0.d0

                ! Compute Final Stress (untwinned region) from existing UMAT
                iumatsvst = offstatev(1) + 1 + (l - 1) * nstatev(1)
                iumatsven = offstatev(1) + l * nstatev(1)
                call umat_stress(rStress_U(:, l),
     &                           svars(kblock, iumatsvst:iumatsven),
     &                           rGsCoords, rFu0, rFu1,
     &                           props(stmprop:endmprop), nmprops,
     &                           nstatev(1), nnode, elmt_number, l,
     &                           mcrd, time, dtimeCur, kstep, kinc,
     &                          .true., rWu, iterErr, Fe(:, :, l))

                theta = theta + svars(kblock, iumatsvst)*(1.d0-bhratio) ! Untwinned Temperature


                ! Compute Final Stress (twinned region) from existing UMAT
                iumatsvst = offstatev(2) + 1 + (l - 1) * nstatev(1)
                iumatsven = offstatev(2) + l * nstatev(1)
                call umat_stress(rStress_Z(:, l),
     &                           svars(kblock, iumatsvst:iumatsven),
     &                           rGsCoords, rFz0h, rFz1h,
     &                           props_tw, nmprops, nstatev(1), nnode,
     &                           elmt_number, l, mcrd, time, dtimeCur,
     &                           kstep, kinc, .true., rWz, iterErr,
     &                           FeZ)

                theta = theta + svars(kblock, iumatsvst)*bhratio    ! Twinned Temperature

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
                        rSigmaU = rStress_U(j, l) - rBVpress
                        rSigmaZ = rStress_Z(j, l) - rBVpress
                      else
                        rSigmaU = rStress_U(j, l)
                        rSigmaZ = rStress_Z(j, l)
                      end if
                      rhs( kblock , (n-1)*mcrd + i ) =
     &                rhs( kblock , (n-1)*mcrd + i ) +
     &                ((bhratio)*rBmatZ(j,(n-1)*mcrd+i)*rSigmaZ +
     &                (1.d0 - bhratio)*rBmatU(j,(n-1)*mcrd+i)*rSigmaU)
     &                *rdOmega
                    end do
                  end do
                end do

                ! vector to matrix for stress in twinned and untwinned region
                call vec_to_sym_mat(rStress_U(:, l), sigmaU)
                call vec_to_sym_mat(rStress_Z(:, l), sigmaZ)

                ! Compute twin length, width and thickness evolution
                Smatc =
     &             matmul(Fe(:,:,l), matmul(Smat(:, :, iTwFlg), Fe_inv))
                tau_tw = sum(sigmaU * Smatc)
                call compute_twin_velocities(props, nprops, tau_tw,
     &                                       theta, l_dot, w_dot, b_dot)

                ! Update length, width and thickness only if this element is a decider
                if (tw_dec(1, twinID).eq.elmt_number) then
                  tw_geom(1, twinID) = tw_geom(1, twinID) +
     &                                 1.d0 / iTIntPt * l_dot * dtimeCur
                end if
                if (tw_dec(2, twinID).eq.elmt_number) then
                  tw_geom(2, twinID) = tw_geom(2, twinID) +
     &                                 1.d0 / iTIntPt * w_dot * dtimeCur
                end if
                if (tw_dec(3, twinID).eq.elmt_number) then
                  tw_geom(3, twinID) = tw_geom(3, twinID) +
     &                                 1.d0 / iTIntPt * b_dot * dtimeCur
                end if

              end do

              ! Twin is embedded in this step, unset the flag
              if (iTwEmbStp.ne.0) then
                iTwEmbStp = 0
                elmt_stat(2, elmt_number) = iTwEmbStp
              end if

              if (kstep.gt.0) then
                ! Initialize element average stress and SDV
                elmt_stress(:, elmt_number) = 0.d0
                elmt_SDV(:, elmt_number) = 0.d0

                ! Write to avg global element stress for both regions
                rStressCt = 0.d0
                rStressCtZ = 0.d0
                do l = 1, iTIntPt
                  rStressCt = rStressCt + rStress_U(:, l)
                  rStressCtZ = rStressCtZ + rStress_Z(:, l)
                end do

                rStressCt = rStressCt / iTIntPt
                rStressCtZ = rStressCtZ / iTIntPt
                elmt_stress(1:6, elmt_number) = rStressCt
                elmt_stress(7:12, elmt_number) = rStressCtZ

                ! Write avg state variables
                outSDVM = 0.d0
                outSDVMZ = 0.d0
                do l = 1, iTIntPt
                  iumatsvst = offstatev(1) + 1 + (l - 1) * nstatev(1)

                  ! Temperature
                  outSDVM(1) = outSDVM(1) + svars(kblock, iumatsvst)

                  ! Euler angles
                  outSDVM(2:4) = outSDVM(2:4) +
     &                             svars(kblock, iumatsvst + 1:3)

                  ! Cumulative dislocation density
                  temp = 0.d0
                  do i = 1, 18
                    temp = temp + svars(kblock, iumatsvst + 36 + i)
                  end do
                  outSDVM(5) = outSDVM(5) + temp / 18

                  ! Effective temperature
                  outSDVM(6) = outSDVM(6) +
     &                             svars(kblock, iumatsvst + 55)

                  ! Angular Misorientation
                  outSDVM(7) = outSDVM(7) +
     &                             svars(kblock, iumatsvst + 56)

                  ! Twinned Region
                  iumatsvst = offstatev(2) + 1 + (l - 1) * nstatev(1)

                  ! Temperature
                  outSDVMZ(1) = outSDVMZ(1) + svars(kblock, iumatsvst)

                  ! Euler angles
                  outSDVMZ(2:4) = outSDVMZ(2:4) +
     &                             svars(kblock, iumatsvst + 1:3)

                  ! Cumulative dislocation density
                  temp = 0.d0
                  do i = 1, 18
                    temp = temp + svars(kblock, iumatsvst + 36 + i)
                  end do
                  outSDVMZ(5) = outSDVMZ(5) + temp / 18

                  ! Effective temperature
                  outSDVMZ(6) = outSDVMZ(6) +
     &                             svars(kblock, iumatsvst + 55)

                  ! Angular Misorientation
                  outSDVMZ(7) = outSDVMZ(7) +
     &                             svars(kblock, iumatsvst + 56)

                end do

                outSDVM = outSDVM / iTIntPt
                outSDVMZ = outSDVMZ / iTIntPt
                elmt_SDV(1:7, elmt_number) = outSDVM
                elmt_SDV(8:14, elmt_number) = outSDVMZ
              end if

            end if

            ! Compute stable time step increment
            call Clength_3D_8N(rCoords, rLmin)
            dtimeStable(kblock) = 0.7d0 * rLmin / rCd
            if (iTwEmbStp) then
              dtimeStable(kblock) = 0.07d0 * rLmin / rCd
            end if
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
        real*8 time

        kStep = i_Array(i_int_kStep)
        kInc  = i_Array(i_int_kInc)
        myThreadID = get_thread_id()
        time = r_Array(i_flt_TotalTime)

        if (lop.eq.j_int_StartAnalysis) then
          call MutexInit(1)      ! initialize Mutex #1
          call MutexInit(2)      ! initialize Mutex #2
          if (myThreadID.eq.0) then
            print *, 'Total number of threads = ', GETNUMTHREADS()
            print *, 'Start of the analysis'
            call read_gb_elements()
          end if
        else if (lop.eq.j_int_EndIncrement) then
          !if (mod(kInc,400).eq.0) then
          call Mutexlock(1)
          if (myThreadID.eq.0) then
            call solve_twin_propagation(kInc)
            if (mod(kInc,1000).eq.0) then
              call write_vuel_results(kInc, time)
            end if
          end if
          call MutexUnlock(1)
          !end if
        else if (lop.eq.j_int_EndAnalysis) then
          if (myThreadID.eq.0) then
            call SMAIntArrayDelete(1)
            call SMAIntArrayDelete(2)
            call SMAIntArrayDelete(3)
            call SMAIntArrayDelete(4)
            call SMAFloatArrayDelete(5)
            call SMAIntArrayDelete(6)
            call SMAIntArrayDelete(7)
            call SMAFloatArrayDelete(8)
            call SMAFloatArrayDelete(9)
            call SMAFloatArrayDelete(10)
            call SMAFloatArrayDelete(11)
            call SMAIntArrayDelete(12)
            call SMAFloatArrayDelete(13)
            call SMAFloatArrayDelete(14)
            call SMAFloatArrayDelete(15)
            call SMAFloatArrayDelete(16)
            call SMAFloatArrayDelete(17)
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

        character(len = 50) :: ecfilename
        character (len = 256) :: ecfilepath
        integer elmt_conn(8, *)
        pointer (ptr_elmt_conn, elmt_conn)
        integer num_elems, dummy
        integer elmt_list(6)
        integer elmt_neighbour(6, *)
        pointer (ptr_elmt_neighbour, elmt_neighbour)
        character(len = 50) :: ncfilename
        character (len = 256) :: ncfilepath
        real*8 nodal_coords(3, *)
        pointer (ptr_nodal_coords, nodal_coords)
        integer num_nodes
        integer size_list(*), elmt_stat(3, *), tw_flg(*)
        pointer (ptr_size_list, size_list)
        pointer (ptr_elmt_stat, elmt_stat)
        pointer (ptr_tw_flg, tw_flg)
        real*8 tw_nuc(3, *), tw_nor(3, *), tw_dir(3, *)
        pointer (ptr_tw_nuc, tw_nuc)
        pointer (ptr_tw_nor, tw_nor)
        pointer (ptr_tw_dir, tw_dir)
        real*8 tw_geom(3, *)
        pointer (ptr_tw_geom, tw_geom)
        integer tw_dec(3, *)
        pointer (ptr_tw_dec, tw_dec)
        real*8 elmt_stress(12, *)
        pointer (ptr_elmt_stress, elmt_stress)
        real*8 nuc_prob(*)
        pointer (ptr_nuc_prob, nuc_prob)
        real*8 disp(3, *)
        pointer (ptr_disp, disp)
        real*8 nuc_thr(*)
        pointer (ptr_nuc_thr, nuc_thr)
        real*8 elmt_SDV(14, *)
        pointer (ptr_elmt_SDV, elmt_SDV)

        ! Read and store grain boundary elements
        gbfilename = 'gb_elems_1.txt'
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

        ptr_gb_elems = SMAIntArrayCreate(2, num_gb_elems, -1)
        open(121, file = gbfilepath, status = 'old',action = 'read')
        do i = 1, num_gb_elems
          read(121, *) gb_elems(i)
        end do
        close(121)

        ! Read and store element connectivity
        ecfilename = 'input_elmt_1.inp'
        ecfilepath = trim(outdir)//pathseperator//trim(ecfilename)
        num_elems = 0
        open(120, file = ecfilepath, status = 'old', iostat = io)
        do
          read(120, *, iostat = io)
          if (io /= 0) exit
          num_elems = num_elems + 1
        end do
        close(120)

        ptr_elmt_conn = SMAIntArrayCreate(3, 8*num_elems, -1)
        open(121, file = ecfilepath, status = 'old',action = 'read')
        do i = 1, num_elems
          read(121, *) dummy, elmt_conn(1:8, i)
        end do
        close(121)

        ! Create neighbouring element list for all elements
        ptr_elmt_neighbour = SMAIntArrayCreate(4, 6*num_elems, -1)
        do i = 1, num_elems
          call Find_Neighbour_Elmts(i, elmt_conn, num_elems, elmt_list)
          elmt_neighbour(1:6, i) = elmt_list
        end do

        ! Read and store nodal coordinates
        ncfilename = 'input_node_1.inp'
        ncfilepath = trim(outdir)//pathseperator//trim(ncfilename)
        num_nodes = 0
        open(120, file = ncfilepath, status = 'old', iostat = io)
        do
          read(120, *, iostat = io)
          if (io /= 0) exit
          num_nodes = num_nodes + 1
        end do
        close(120)

        ptr_nodal_coords = SMAFloatArrayCreate(5, 3*num_nodes, 0.d0)
        open(121, file = ncfilepath, status = 'old',action = 'read')
        do i = 1, num_nodes
          read(121, *) dummy, nodal_coords(1:3, i)
        end do
        close(121)

        ! Store sizes and counters of certain quantities
        ptr_size_list = SMAIntArrayCreate(1, 4, -1)
        size_list(1) = num_gb_elems
        size_list(2) = num_elems
        size_list(3) = num_nodes
        size_list(4) = 0

        ! Create element status for twinning ID and embedding step flag
        ptr_elmt_stat = SMAIntArrayCreate(6, 3*num_elems, 0)

        ! Create arrays for twin type, nucleation coords, normal and direction
        ptr_tw_flg = SMAIntArrayCreate(7, 100, 0)
        ptr_tw_nuc = SMAFloatArrayCreate(8, 3*100, 0.d0)
        ptr_tw_nor = SMAFloatArrayCreate(9, 3*100, 0.d0)
        ptr_tw_dir = SMAFloatArrayCreate(10, 3*100, 0.d0)

        ! Create global twin size array and the decider element id array
        ptr_tw_geom = SMAFloatArrayCreate(11, 3*100, 0.d0)
        ptr_tw_dec = SMAIntArrayCreate(12, 3*100, 0)

        ! Elmt average stress components
        ptr_elmt_stress = SMAFloatArrayCreate(13, 12*num_elems, 0.d0)

        ! Nucleation probability
        ptr_nuc_prob = SMAFloatArrayCreate(14, num_elems, 0.d0)

        ! Nodal displacement
        ptr_disp = SMAFloatArrayCreate(15, 3*num_nodes, 0.d0)

        ! Read threshold probabilities for gb elems
        ! gbfilename = 'gb_elems_8gn_thr.txt'
        ! call vgetoutdir(outdir, lenoutdir)
        ! gbfilepath = trim(outdir)//pathseperator//trim(gbfilename)
        ! open(121, file = gbfilepath, status = 'old',action = 'read')
        ptr_nuc_thr = SMAFloatArrayCreate(16, num_gb_elems, 0.d0)
        do i = 1, num_gb_elems
          ! read(121, *) nuc_thr(i)
          nuc_thr(i) = 0.85d0
        end do
        ! close(121)

        ! Element state dependent variables
        ptr_elmt_SDV = SMAFloatArrayCreate(17, 14*num_elems, 0.d0)

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          Detect Twin propagation in neighbouring elements            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine solve_twin_propagation(kInc)
        implicit none
#include <SMAAspUserSubroutines.hdr>
        ! Define variables and cray pointers
        integer elmt_conn(8, *)
        pointer (ptr_elmt_conn, elmt_conn)
        integer elmt_neighbour(6, *)
        pointer (ptr_elmt_neighbour, elmt_neighbour)
        real*8 nodal_coords(3, *)
        pointer (ptr_nodal_coords, nodal_coords)
        integer num_nodes, num_elems
        integer size_list(*), elmt_stat(3, *), tw_flg(*)
        pointer (ptr_size_list, size_list)
        pointer (ptr_elmt_stat, elmt_stat)
        pointer (ptr_tw_flg, tw_flg)
        real*8 tw_nuc(3, *), tw_nor(3, *), tw_dir(3, *)
        pointer (ptr_tw_nuc, tw_nuc)
        pointer (ptr_tw_nor, tw_nor)
        pointer (ptr_tw_dir, tw_dir)
        real*8 tw_geom(3, *)
        pointer (ptr_tw_geom, tw_geom)
        integer tw_dec(3, *)
        pointer (ptr_tw_dec, tw_dec)
        real*8 elmt_stress(12, *)
        pointer (ptr_elmt_stress, elmt_stress)
        integer num_twins, itw, iel, inb, inbel
        integer iTintPt, l, ind, i
        real*8 rXiEtaPt(3, 8), rWeightPt(8), rNn(8), rdNndXi(3, 8)
        real*8 rElmtVol, rJdet, rdOmega, rElmtSize, rCoords(3, 8)
        integer, parameter :: num_dis = 10, num_dis_pts = num_dis + 1
        real*8 rad, s_tw(3), xi_dis(num_dis_pts), tw_el_len, tw_el_wid
        real*8 tw_el_thk
        logical twinned, bexp
        integer nxt_l_elmt, nxt_w_elmt, nxt_lm_elmt, nxt_wm_elmt
        integer nxt_bm_elmt, nxt_b_elmt
        character (len = 1) :: pathseperator = '/'
        character (len = 256) :: outdir, elstfilepath
        integer :: lenoutdir
        integer kInc
        character*20 str_kinc

        ! Access number of twins from Abaqus SMA array
        ptr_size_list = SMAIntArrayAccess(1)
        num_twins = size_list(4)
        if (num_twins.eq.0) return ! if no twins yet, return

        ! Access rest of the twin data
        ptr_elmt_conn = SMAIntArrayAccess(3)
        ptr_elmt_neighbour = SMAIntArrayAccess(4)
        ptr_nodal_coords = SMAFloatArrayAccess(5)
        ptr_elmt_stat = SMAIntArrayAccess(6)
        ptr_tw_flg = SMAIntArrayAccess(7)
        ptr_tw_nuc = SMAFloatArrayAccess(8)
        ptr_tw_nor = SMAFloatArrayAccess(9)
        ptr_tw_dir = SMAFloatArrayAccess(10)
        ptr_tw_geom = SMAFloatArrayAccess(11)
        ptr_tw_dec = SMAIntArrayAccess(12)
        ptr_elmt_stress = SMAFloatArrayAccess(13)
        num_elems = size_list(2)
        num_nodes = size_list(3)

        ! Num gauss pts, natural coords and weight vals required for elmt size
        call gspt_3d_8n(iTIntPt, rXiEtaPt, rWeightPt)

        ! Natural coords at discretization points (10 currently on each axis)
        do i = 1, num_dis_pts
          xi_dis(i) = -1.d0 + (i - 1) * 2.d0 / num_dis
        end do

        ! Loop over all the twins
        do itw = 1, num_twins
          ! Compute the twin lateral direction
          call cross_product(tw_nor(:, itw), tw_dir(:, itw), s_tw)

          ! Loop over all elements
          do iel = 1, num_elems
            ! Only if element has current twin status
            if (elmt_stat(1, iel).eq.itw) then
              ! Find nxt twin elmt in front (+-) and lateral (+-) direction
              call find_next_twin_elem(iel, elmt_conn, num_elems,
     &             nodal_coords, num_nodes, tw_nuc(:, itw),
     &             tw_nor(:, itw), -tw_dir(:, itw), elmt_neighbour,
     &             tw_el_len, nxt_lm_elmt)
              call find_next_twin_elem(iel, elmt_conn, num_elems,
     &             nodal_coords, num_nodes, tw_nuc(:, itw),
     &             tw_nor(:, itw), tw_dir(:, itw), elmt_neighbour,
     &             tw_el_len, nxt_l_elmt)
              call find_next_twin_elem(iel, elmt_conn, num_elems,
     &             nodal_coords, num_nodes, tw_nuc(:, itw),
     &             tw_nor(:, itw), -s_tw, elmt_neighbour,
     &             tw_el_wid, nxt_wm_elmt)
              call find_next_twin_elem(iel, elmt_conn, num_elems,
     &             nodal_coords, num_nodes, tw_nuc(:, itw),
     &             tw_nor(:, itw), s_tw, elmt_neighbour,
     &             tw_el_wid, nxt_w_elmt)
              call find_next_twin_elem(iel, elmt_conn, num_elems,
     &             nodal_coords, num_nodes, tw_nuc(:, itw),
     &             -tw_dir(:, itw), -tw_nor(:, itw), elmt_neighbour,
     &             tw_el_thk, nxt_bm_elmt)
              call find_next_twin_elem(iel, elmt_conn, num_elems,
     &             nodal_coords, num_nodes, tw_nuc(:, itw),
     &             -tw_dir(:, itw), tw_nor(:, itw), elmt_neighbour,
     &             tw_el_thk, nxt_b_elmt)

              ! Loop over neighbours
              do inb = 1, 6
                inbel = elmt_neighbour(inb, iel)
                if (inbel.eq.0) cycle ! Face is not attached to neighbouring element
                if (elmt_stat(1, inbel).ne.0) cycle ! Already twinned element
                if (elmt_stat(3, inbel).ne.elmt_stat(3, iel)) cycle ! Different grain element

                ! Retrive element nodal coords
                do ind = 1, 8
                  rCoords(:, ind) =
     &                          nodal_coords(:, elmt_conn(ind, inbel))
                end do

                ! Compute element volume
                rElmtVol = 0.d0
                do l = 1, iTIntPt
                  ! Jacobian determinant
                  call shp_func_3d_8n(rXiEtaPt(1,l),rCoords,0,
     &                            rNn,rdNndXi,rJdet)
                  ! Integration element (dOmega)
                  rdOmega = rJdet * rWeightPt(l)
                  ! Element surface area (even for axisymmetric elements)
                  rElmtVol = rElmtVol + rdOmega
                end do

                ! Compute element size h
                rElmtSize = rElmtVol**(1.d0/3.d0)      !   h = cubic root of vol

                ! Radius of search for each discretization point in detect_twinning subroutine
                rad = rElmtSize * sqrt(3.d0) / num_dis

                ! Detect Twinning in this neighbouring element
                call detect_twinning(inbel, elmt_conn, num_elems,
     &                nodal_coords, num_nodes, num_dis_pts,
     &                tw_nor(:, itw), tw_dir(:, itw), s_tw, rad,
     &                tw_nuc(:, itw), tw_geom(1, itw), tw_geom(2, itw),
     &                tw_geom(3, itw),  xi_dis, twinned, bexp)

                ! Update the decider if this is twinned and the nxt twn element in either direction
                if (twinned) then
                  elmt_stat(1, inbel) = itw
                  if (bexp) then
                    elmt_stat(2, inbel) = 3
                  else
                    elmt_stat(2, inbel) = 2   ! Set twin embedding step to 2
                  end if
                  if (inbel.eq.nxt_lm_elmt) tw_dec(1, itw) = inbel
                  if (inbel.eq.nxt_l_elmt) tw_dec(1, itw) = inbel
                  if (inbel.eq.nxt_wm_elmt) tw_dec(2, itw) = inbel
                  if (inbel.eq.nxt_w_elmt) tw_dec(2, itw) = inbel
                  if (inbel.eq.nxt_bm_elmt) tw_dec(3, itw) = inbel
                  if (inbel.eq.nxt_b_elmt) tw_dec(3, itw) = inbel
                end if
              end do
            end if
          end do
        end do
      end subroutine

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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine write_vuel_results(kInc, time)
        implicit none
#include <SMAAspUserSubroutines.hdr>
        integer kInc
        integer num_nodes, num_elems, iel, io, ind
        integer size_list(*), elmt_stat(3, *)
        pointer (ptr_size_list, size_list)
        pointer (ptr_elmt_stat, elmt_stat)
        real*8 elmt_stress(12, *)
        pointer (ptr_elmt_stress, elmt_stress)
        real*8 nuc_prob(*)
        pointer (ptr_nuc_prob, nuc_prob)
        real*8 disp(3, *)
        pointer (ptr_disp, disp)
        real*8 elmt_SDV(14, *)
        pointer (ptr_elmt_SDV, elmt_SDV)
        character (len = 1) :: pathseperator = '/'
        character (len = 256) :: outdir, filepath
        integer :: lenoutdir
        character*20 str_kinc
        real*8 time
        logical exists

        ! Retrieve data from global shared arrays
        ptr_size_list = SMAIntArrayAccess(1)
        ptr_elmt_stat = SMAIntArrayAccess(6)
        ptr_elmt_stress = SMAFloatArrayAccess(13)
        ptr_nuc_prob = SMAFloatArrayAccess(14)
        ptr_disp = SMAFloatArrayAccess(15)
        ptr_elmt_SDV = SMAFloatArrayAccess(17)
        num_elems = size_list(2)
        num_nodes = size_list(3)

        ! Write element status with twin id
        write(str_kinc, '(I0)') kInc
        call vgetoutdir(outdir, lenoutdir)
        filepath = trim(outdir)//pathseperator//'kinc_elmt_stat_'
     &                  //trim(adjustl(str_kinc))//'.dat'

        open(newunit = io, file = filepath, status = 'new')
        do iel=1,num_elems
          write(io,*) elmt_stat(1, iel)
        end do

        close(io)

        ! Write element stress for both regions
        filepath = trim(outdir)//pathseperator//'kinc_elmt_stress_'
     &                  //trim(adjustl(str_kinc))//'.dat'

        open(newunit = io, file = filepath, status = 'new')
        do iel=1,num_elems
          write(io,*) elmt_stress(:, iel)
        end do

        close(io)

        ! Write nucleation probability
        filepath = trim(outdir)//pathseperator//'kinc_nuc_prob_'
     &                  //trim(adjustl(str_kinc))//'.dat'

        open(newunit = io, file = filepath, status = 'new')
        do iel=1,num_elems
          write(io,*) nuc_prob(iel)
        end do

        close(io)

        ! Write nodal displacements
        filepath = trim(outdir)//pathseperator//'kinc_disp_'
     &                  //trim(adjustl(str_kinc))//'.dat'

        open(newunit = io, file = filepath, status = 'new')
        do ind=1,num_nodes
          write(io,*) disp(:, ind)
        end do

        close(io)

        ! Write time for this increment
        filepath = trim(outdir)//pathseperator//'kinc_time.dat'
        inquire(file=filepath, exist=exists)
        if (exists) then
          open(newunit = io, file = filepath, status = "old",
     &         position = "append")
        else
          open(newunit = io, file = filepath, status = "new")
        end if

        write(io, *) kInc, time
        close(io)

        ! Write element SDV for both regions
        filepath = trim(outdir)//pathseperator//'kinc_elmt_SDV_'
     &                  //trim(adjustl(str_kinc))//'.dat'

        open(newunit = io, file = filepath, status = 'new')
        do iel=1,num_elems
          write(io,1000) elmt_SDV(:, iel)
        end do

        close(io)

1000  FORMAT (14(E26.16))

      end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      include 'GSPT_3D_8N.f'
      include 'SHP_FUNC_3D_8N.f'
      include 'prep_umat.f'
      include 'gdefgrad_3D.f'
      include 'hcp_umat_TDT_pq_Ti.f'
      include 'Gstrain_3D.f'
      include 'calc_new_jac.f'
      include 'Clength_3D_8N.f'
      include 'trac_cont.f'
