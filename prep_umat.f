      subroutine calc_dfgrd(rXiEtaPt, rRefCoords, rPrevCoords,
     &                      rCurrCoords, mcrd, nnode, rDfgrd0,
     &                      rDfgrd1, rPtCoords)
          implicit none

C         integer inputs
          integer mcrd, nnode

C         Reference, Current and Previous Coordinates and Natural Coords
          real*8 rRefCoords(mcrd, nnode), rCurrCoords(mcrd, nnode)
          real*8 rPrevCoords(mcrd, nnode), rXiEtaPt(3)

C         Deformation Gradients(previous and current) and Point Coords
          real*8 rDfgrd0(3, 3), rDfgrd1(3, 3), rPtCoords(3)

C         Shape fcns, Derivatives and Determinant
          real*8 rN(8), rdNdXRef(3, nnode), rJdet

C         loop integers
          integer idx, jdx

C         Construct Shape Functions and Derivatives
C         // w.r.t Reference Coordinates
          call SHP_FUNC_3D_8N(rXiEtaPt, rRefCoords, 1, rN, rdNdXRef,
     &                        rJdet)

C         Deformation Gradients of previous and current coords
C         and their determinants
          call GdefGrad_3D(rdNdXRef, rPrevCoords, mcrd, nnode, 1,
     &                     rDfgrd0)
          call GdefGrad_3D(rdNdXRef, rCurrCoords, mcrd, nnode, 1,
     &                     rDfgrd1)

C         Compute Point Coordinates
          rPtCoords = 0.d0
          do idx = 1, nnode
               do jdx = 1, mcrd
                    rPtCoords(jdx) = rPtCoords(jdx) +
     &                        rN(idx) * rCurrCoords(jdx, idx)
               end do
          end do

      end subroutine

      subroutine calc_dfgrd2(rXiEtaPt, rRefCoords, rCurrCoords, mcrd,
     &                       nnode, rDfgrd1, rPtCoords)
          implicit none

C         integer inputs
          integer mcrd, nnode

C         Reference, Current and Previous Coordinates and Natural Coords
          real*8 rRefCoords(mcrd, nnode), rCurrCoords(mcrd, nnode)
          real*8 rXiEtaPt(3)

C         Deformation Gradients(previous and current) and Point Coords
          real*8 rDfgrd1(3, 3), rPtCoords(3)

C         Shape fcns, Derivatives and Determinant
          real*8 rN(8), rdNdXRef(3, nnode), rJdet

C         loop integers
          integer idx, jdx

C         Construct Shape Functions and Derivatives
C         // w.r.t Reference Coordinates
          call SHP_FUNC_3D_8N(rXiEtaPt, rRefCoords, 1, rN, rdNdXRef,
     &                        rJdet)

C         Deformation Gradients of previous and current coords
C         and their determinants
          call GdefGrad_3D(rdNdXRef, rCurrCoords, mcrd, nnode, 1,
     &                     rDfgrd1)

C         Compute Point Coordinates
          rPtCoords = 0.d0
          do idx = 1, nnode
               do jdx = 1, mcrd
                    rPtCoords(jdx) = rPtCoords(jdx) +
     &                        rN(idx) * rCurrCoords(jdx, idx)
               end do
          end do

      end subroutine

      subroutine dfgrd_prep2(iTIntPt, rRefCoords, rCurrCoords,
     &                       rXiEtaPt, mcrd, nnode, rDfgrd1, rGsCoords)
C         This subroutine prepares deformation gradient and gauss point
C         coordinates for UMAT type function.
          implicit none

C         Integer inputs
          integer iTIntPt, mcrd, nnode

C         Reference, Current and Previous Coordinates
          real*8 rRefCoords(mcrd, nnode), rCurrCoords(mcrd, nnode)
          real*8 rXiEtaPt(3, iTIntPt)

C         Initalize Deformation Gradient at each Gauss Pt
C         and Gauss Pt Coordinates
          real*8 rDfgrd1(3, 3, iTIntPt)
          real*8 rGsCoords(3, iTIntPt)

C         Additional variables for computations
          integer l
          real*8 rDetF1Bar, rSelRedIntConst
          real*8 rDetF1(iTIntPt)
          real*8, external :: determinant ! external determinant fcn

C         Compute Deformation Gradients and Gauss Pt Coordinates
C         at each guass pt
          rDetF1Bar = 0.d0
          do l = 1, iTIntPt
               call calc_dfgrd2(rXiEtaPt(:, l), rRefCoords, rCurrCoords,
     &                          mcrd, nnode, rDfgrd1(:, :, l),
     &                          rGsCoords(:, l))
               rDetF1(l) = determinant(rDfgrd1(:,:,l))
               rDetF1Bar = rDetF1Bar + rDetF1(l)
          end do

C         Scale Deformation Gradient for Selective Reduced Integration
C         See UMAT Deformation Gradient Documentation for more details
          rDetF1Bar = rDetF1Bar / iTIntPt
          do l = 1, iTIntPt
               rSelRedIntConst = (rDetF1Bar / rDetF1(l))**(1.0/3)
               rDfgrd1(:, :, l) = rDfgrd1(:, :, l) * rSelRedIntConst
          end do
      end subroutine

      subroutine comp_vel_gr(rdVelGrd, rStrainInc, rSpinInc)
          implicit none
          real*8 rdVelGrd(3, 3), rStrainInc(6), rSpinInc(3)

          rdVelGrd(1, 1) = rStrainInc(1)
          rdVelGrd(2, 1) = 0.5 * rStrainInc(6) + rSpinInc(3)
          rdVelGrd(3, 1) = 0.5 * rStrainInc(5) - rSpinInc(2)
          rdVelGrd(1, 2) = 0.5 * rStrainInc(6) - rSpinInc(3)
          rdVelGrd(2, 2) = rStrainInc(2)
          rdVelGrd(3, 2) = 0.5 * rStrainInc(4) + rSpinInc(1)
          rdVelGrd(1, 3) = 0.5 * rStrainInc(5) + rSpinInc(2)
          rdVelGrd(2, 3) = 0.5 * rStrainInc(4) - rSpinInc(1)
          rdVelGrd(3, 3) = rStrainInc(3)
      end subroutine

      subroutine comp_twin_dfgrd(rImat, props, nprops, pi, Smat, F_tw,
     &                           rnc, rmc)
          implicit none
          real*8 rn(3, 18), rm(3, 18), rmc(3, 18), rnc(3, 18)
          real*8 gamma1, gamma2, gamma3
          real*8 Smat(3, 3, 18), F_tw(3, 3, 18), rImat(3, 3)
          integer nprops
          real*8 props(nprops), th, phi, om, Q_0(3, 3), pi
          integer i, j

          ! Tensile Twins {1 0 -1 2} <-1 0 1 1>
          rn(:, 1) = (/0.585114161329908d0, 0.337815818550484d0,
     &                 0.737239371539709d0/)
          rm(:, 1) = (/-0.638468024423463d0, -0.368619685769855d0,
     &                 0.675631637100968d0/)
          rn(:, 2) = (/0.000000000000000d0, 0.675631637100968d0,
     &                 0.737239371539709d0/)
          rm(:, 2) = (/0.000000000000000d0, -0.737239371539709d0,
     &                 0.675631637100968d0/)
          rn(:, 3) = (/-0.585114161329907d0, 0.337815818550484d0,
     &                 0.737239371539709d0/)
          rm(:, 3) = (/0.638468024423463d0, -0.368619685769855d0,
     &                 0.675631637100968d0/)
          rn(:, 4) = (/-0.585114161329907d0, -0.337815818550484d0,
     &                  0.737239371539709d0/)
          rm(:, 4) = (/0.638468024423463d0, 0.368619685769855d0,
     &                 0.675631637100968d0/)
          rn(:, 5) = (/0.000000000000000d0, -0.675631637100968d0,
     &                 0.737239371539709d0/)
          rm(:, 5) = (/0.000000000000000d0, 0.737239371539709d0,
     &                 0.675631637100968d0/)
          rn(:, 6) = (/0.585114161329907d0, -0.337815818550484d0,
     &                 0.737239371539709d0/)
          rm(:, 6) = (/-0.638468024423462d0, 0.368619685769855d0,
     &                 0.675631637100968d0/)

          ! Compression Twins {1 1 -2 2} <1 1 -2 -3>
          rn(:, 7) = (/0.4230464880872873d0, 0.7327380113307633d0,
     &                 0.5330353418564765d0/)
          rm(:, 7) = (/0.26651767092823825d0, 0.4616221471626314d0,
     &                -0.8460929761745747d0/)
          rn(:, 8) = (/-0.4230464880872873d0, 0.7327380113307633d0,
     &                 0.5330353418564765d0/)
          rm(:, 8) = (/-0.26651767092823825d0, 0.4616221471626314d0,
     &                -0.8460929761745747d0/)
          rn(:, 9) = (/-0.8460929761745746d0, 0.0000000000000000d0,
     &                 0.5330353418564765d0/)
          rm(:, 9) = (/-0.5330353418564765d0, 0.00000000000000000d0,
     &                -0.8460929761745747d0/)
          rn(:, 10) = (/-0.4230464880872873d0, -0.7327380113307633d0,
     &                 0.5330353418564765d0/)
          rm(:, 10) = (/-0.26651767092823825d0, -0.4616221471626314d0,
     &                -0.8460929761745747d0/)
          rn(:, 11) = (/0.4230464880872873d0, -0.7327380113307633d0,
     &                 0.5330353418564765d0/)
          rm(:, 11) = (/0.26651767092823825d0, -0.4616221471626314d0,
     &                -0.8460929761745747d0/)
          rn(:, 12) = (/0.8460929761745746d0, 0.0000000000000000d0,
     &                 0.5330353418564765d0/)
          rm(:, 12) = (/0.5330353418564765d0, 0.00000000000000000d0,
     &                -0.8460929761745747d0/)

          ! Compression Twins {1 0 -1 1} <1 0 -1 -2>
          rn(:, 13) = (/0.7602357125290463d0, 0.4389222932762118d0,
     &                 0.47894559384193547d0/)
          rm(:, 13) = (/0.41477905129774d0, 0.23947279692096776d0,
     &                -0.8778445865524238d0/)
          rn(:, 14) = (/0.0000000000000000d0, 0.8778445865524237d0,
     &                 0.4789455938419355d0/)
          rm(:, 14) = (/0.0000000000000000d0, 0.4789455938419355d0,
     &                -0.8778445865524238d0/)
          rn(:, 15) = (/-0.7602357125290463d0, 0.4389222932762118d0,
     &                 0.47894559384193547d0/)
          rm(:, 15) = (/-0.41477905129774d0, 0.23947279692096776d0,
     &                -0.8778445865524238d0/)
          rn(:, 16) = (/-0.7602357125290463d0, -0.4389222932762118d0,
     &                 0.47894559384193547d0/)
          rm(:, 16) = (/-0.41477905129774d0, -0.23947279692096776d0,
     &                -0.8778445865524238d0/)
          rn(:, 17) = (/0.0000000000000000d0, -0.8778445865524237d0,
     &                 0.4789455938419355d0/)
          rm(:, 17) = (/0.0000000000000000d0, -0.4789455938419355d0,
     &                -0.8778445865524238d0/)
          rn(:, 18) = (/0.7602357125290463d0, -0.4389222932762118d0,
     &                 0.47894559384193547d0/)
          rm(:, 18) = (/0.41477905129774d0, -0.23947279692096776d0,
     &                -0.8778445865524238d0/)

          ! Shear in twin systems
          gamma1 = 0.17475083284085213d0     ! {1 0 -1 2} <-1 0 1 1> (tens)
          gamma2 = 0.21821257173587794d0     ! {1 1 -2 2} <1 1 -2 -3> (comp)
          gamma3 = 0.09804551470109034d0     ! {1 0 -1 1} <1 0 -1 -2> (comp)

          ! Compute rotation matrix from euler angles and rotate twin basis
          th = props(31) * pi / 180.d0
          phi = props(32) * pi / 180.d0
          om = props(33) * pi / 180.d0
          call rotmat(th, phi, om, Q_0)
          rnc = matmul(Q_0, rn)
          rmc = matmul(Q_0, rm)

          ! Form Twin Schmid Tensor
          do i = 1, 3
            do j = 1, 3
                Smat(i, j, :) = rmc(i, :) * rnc(j, :)
            end do
          end do

          ! Compute Twinning Deformation gradient
          do i = 1, 6
            F_tw(:, :, i) = rImat + gamma1 * Smat(:, :, i)
          end do
          do i = 7, 12
            F_tw(:, :, i) = rImat + gamma2 * Smat(:, :, i)
          end do
          do i = 13, 18
            F_tw(:, :, i) = rImat + gamma3 * Smat(:, :, i)
          end do

      end subroutine

      subroutine dfgrd_prep(iTIntPt, rRefCoords, rCurrCoords,
     &                      rPrevCoords, rXiEtaPt, mcrd, nnode, rDfgrd0,
     &                      rDfgrd1, rGsCoords)
C         This subroutine prepares deformation gradient and gauss point
C         coordinates for UMAT type function.
          implicit none

C         Integer inputs
          integer iTIntPt, mcrd, nnode

C         Reference, Current and Previous Coordinates
          real*8 rRefCoords(mcrd, nnode), rCurrCoords(mcrd, nnode)
          real*8 rPrevCoords(mcrd, nnode), rXiEtaPt(3, iTIntPt)

C         Initalize Deformation Gradient at each Gauss Pt
C         and Gauss Pt Coordinates
          real*8 rDfgrd0(3, 3, iTIntPt), rDfgrd1(3, 3, iTIntPt)
          real*8 rGsCoords(3, iTIntPt)

C         Additional variables for computations
          integer l
          real*8 rDetF0Bar, rDetF1Bar, rSelRedIntConst
          real*8 rDetF0(iTIntPt), rDetF1(iTIntPt)
          real*8, external :: determinant ! external determinant fcn

C         Compute Deformation Gradients and Gauss Pt Coordinates
C         at each guass pt
          rDetF0Bar = 0.d0
          rDetF1Bar = 0.d0
          do l = 1, iTIntPt
               call calc_dfgrd(rXiEtaPt(:, l), rRefCoords, rPrevCoords,
     &                         rCurrCoords, mcrd, nnode, rDfgrd0(:,:,l),
     &                         rDfgrd1(:, :, l), rGsCoords(:, l))
               rDetF0(l) = determinant(rDfgrd0(:,:,l))
               rDetF1(l) = determinant(rDfgrd1(:,:,l))
               rDetF0Bar = rDetF0Bar + rDetF0(l)
               rDetF1Bar = rDetF1Bar + rDetF1(l)
          end do

C         Scale Deformation Gradient for Selective Reduced Integration
C         See UMAT Deformation Gradient Documentation for more details
          rDetF0Bar = rDetF0Bar / iTIntPt
          rDetF1Bar = rDetF1Bar / iTIntPt
          do l = 1, iTIntPt
               rSelRedIntConst = (rDetF0Bar / rDetF0(l))**(1.0/3)
               rDfgrd0(:, :, l) = rDfgrd0(:, :, l) * rSelRedIntConst
               rSelRedIntConst = (rDetF1Bar / rDetF1(l))**(1.0/3)
               rDfgrd1(:, :, l) = rDfgrd1(:, :, l) * rSelRedIntConst
          end do

      end subroutine

      subroutine umat_stress(rStress, hvars, rGsCoords, rDfgrd0,
     &                      rDfgrd1, props, nprops, nsvars, nnode, noel,
     &                      npt, mcrd, time, dtime, kstep, kinc,
     &                      lEditSvars, rW, iterErr, Fe)
          
          implicit none
C         input arguments declaration
          integer nprops, nsvars, nnode, noel, npt, mcrd
          integer kstep, kinc
          real*8 rStress(6), hvars(nsvars), props(nprops)
          real*8 rGsCoords(3), rDfgrd0(3, 3), rDfgrd1(3, 3)
          real*8 time(2)
          real*8 dtime
          integer idx, jdx
          logical lEditSvars

C         additional arguments for umat
          real*8 ddsdde(6, 6), ddsddt(6), drplde(6), stran(6)
          real*8 dstran(6), drot(3, 3)
          real*8 sse, spd, scd, rpl, drpldt, temp, dtemp
          real*8 predef(1), dpred(1), pnewdt, celent
          character*80 cmname
          integer ndi, nshr, ntens, layer, kspt
          real*8 tmphvars(nsvars), rW(3,3,3,3), Fe(3, 3)
          integer iterErr

C         for umat stress to vuel stress
          real*8 rUmatStress(6)
          
          ndi = 3
          nshr = 3
          ntens = ndi + nshr
          rUmatStress = 0.d0
          rW = 0.d0
          if (((kstep.eq.0).or.(kstep.eq.1)).and.(kinc.eq.0)) then
              dtime = 1.d-9
          end if

          if (lEditSvars) then    ! During regular update
            call umat_crys(rUmatStress,  hvars,  ddsdde,  sse,  spd,
     &               scd,     rpl,     ddsddt,  drplde,  drpldt,
     &               stran,  dstran, time,    dtime,   temp,
     &               dtemp,   predef,  dpred,   cmname,  ndi,
     &               nshr,    ntens,   nsvars,  props,   nprops,
     &               rGsCoords,  drot,    pnewdt,  celent,  rDfgrd0,
     &               rDfgrd1,  noel,    npt,     layer,   kspt,
     &               kstep,   kinc, rW, iterErr, Fe)
          else    ! Required during iterative updates when hvars should not be edited
            tmphvars = hvars
            call umat_crys(rUmatStress,  tmphvars,  ddsdde,  sse,  spd,
     &               scd,     rpl,     ddsddt,  drplde,  drpldt,
     &               stran,  dstran, time,    dtime,   temp,
     &               dtemp,   predef,  dpred,   cmname,  ndi,
     &               nshr,    ntens,   nsvars,  props,   nprops,
     &               rGsCoords,  drot,    pnewdt,  celent,  rDfgrd0,
     &               rDfgrd1,  noel,    npt,     layer,   kspt,
     &               kstep,   kinc, rW, iterErr, Fe)
          end if

          ! if (itererr.eq.1) then
          !   CALL XPLB_ABQERR(-3,'UMAT convergence failed!',0,0.d0,' ')
          ! end if

          rStress(1) = rUmatStress(1)
          rStress(2) = rUmatStress(2)
          rStress(3) = rUmatStress(3)
          rStress(4) = rUmatStress(6)
          rStress(5) = rUmatStress(5)
          rStress(6) = rUmatStress(4)
      end subroutine
      
      ! Simple finite strain linear elastic umat
      subroutine umat_crys2(stress,  statev,  ddsdde,  sse,     spd,
     &               scd,     rpl,     ddsddt,  drplde,  drpldt,
     &               stran,  dstran, time,    dtime,   temp,
     &               dtemp,   predef,  dpred,   cmname,  ndi,
     &               nshr,    ntens,   nstatv,  props,   nprops,
     &               coords,  drot,    pnewdt,  celent,  dfgrd0,
     &               dfgrd1,  noel,    npt,     layer,   kspt,
     &               kstep,   kinc)
      
      implicit real*8(a-h,o-z)
      dimension
     &  coords(3),           ! Coordinates of Gauss point
     &  ddsdde(ntens,ntens), ! Tangent Stiffness Matrix
     &  ddsddt(ntens),       ! Derivative of stress wrt temperature
     &  dfgrd0(3,3),         ! Deformation gradient at beginning of step
     &  dfgrd1(3,3),         ! Deformation gradient at end of step
     &  dpred(1),            ! Predefined state variable array increment
     &  drplde(ntens),       ! Derivative of heat generation wrt strain
     &  drot(3,3),           ! Rotation matrix increment 
     &  dstran(ntens),       ! Strain vector increment
     &  predef(1),           ! Predefined variable array
     &  props(nprops),       ! Material properties array
     &  statev(nstatv),      ! State variables array
     &  stran(ntens),        ! Strain vector at beginning of increment
     &  stress(ntens),       ! Cauchy stress vector
     &  time(2)              ! Step time (1) and total time (2)
        character*80 cmname

        ! Local subroutine variables
        real*8 c11, c12, c44, c33, c55, c13, theta, pi
        real*8 th, phi, om, detF
        real*8 rCe(3, 3), rEe(3, 3), del(3, 3), rS(3, 3)
        real*8 C_c(3, 3, 3, 3), C(3, 3, 3, 3), rSig(3, 3)
        real*8 q_0(3, 3)

        pi = 4.d0 * atan(1.d0) ! PI

        ! Identity tensor
        del = 0.d0
        del(1, 1) = 1.d0
        del(2, 2) = 1.d0
        del(3, 3) = 1.d0

        ! Retrieve temperature, euler angels and elastic properties
        theta = props(1)
        th = props(31)*pi/180.d0     ! Euler angle theta [rad] (Canova convention)
        phi = props(32)*pi/180.d0    ! Euler angle Phi [rad]
        om = props(33)*pi/180.d0     ! Euler angle omega [rad]

        call compute_elastconst_hcp(props,nprops,theta,c11,c12,c44,
     &                              c33, c55, c13)

        ! Construct Stiffness in unrotated configuration
        C_c(:,:,:,:) = 0.d0
        do i = 1,3
          do j = 1,3
            do k = 1,3
              do l = 1,3
                C_c(i,j,k,l) = C12 * del(i,j) * del(k,l) +
     &            C44 * (del(i,k)*del(j,l)+del(i,l)*del(k,j))
              end do
            end do
          end do
        end do
        C_c(1,1,1,1) = C11
        C_c(2,2,2,2) = C11
        C_c(3,3,3,3) = C33
        C_c(1,1,3,3) = c13
        C_c(3,3,1,1) = c13
        C_c(2,2,3,3) = c13
        C_c(3,3,2,2) = c13
        C_c(1,3,1,3) = c55
        C_c(1,3,3,1) = c55
        C_c(3,1,1,3) = c55
        C_c(3,1,3,1) = c55
        C_c(2,3,2,3) = c55
        C_c(2,3,3,2) = c55
        C_c(3,2,2,3) = c55
        C_c(3,2,3,2) = c55

        ! Rotate Stiffness along crystal orientation
        call ROTMAT(TH,PHI,OM,q_0)
        call rotate_4th(q_0,C_c,C)

        ! Find Green Lagrange strain and then 2nd PK stress (linear elastic)
        rCe = matmul(transpose(dfgrd1), dfgrd1)
        rEe = 0.5d0 * (rCe - del)
        rS = reshape( (/((sum(c(i,j,:,:)*rEe),
     &          i = 1,3), j=1,3)/), (/3,3/))

        ! Push forward 2nd PK to Cauchy stress
        detF = determinant(dfgrd1)
        rSig = 1.d0 / detF *
     &    matmul(dfgrd1, matmul(rS, transpose(dfgrd1)))

        ! Vector output of Cauchy stress
        stress(1) = rSig(1, 1)
        stress(2) = rSig(2, 2)
        stress(3) = rSig(3, 3)
        stress(4) = rSig(1, 2)
        stress(5) = rSig(1, 3)
        stress(6) = rSig(2, 3)
      end subroutine

      subroutine check_twin_nucleation(stress, Smat, dfgrd, iTwinFlag,
     &                                 PtauMax)
        implicit none

        integer :: iTwinFlag, i ! Twinning nucleation flag
        real*8 stress(6)      ! Input stress vector
        real*8 Smat(3, 3, 18) ! Twinning Schmid Tensor
        real*8 SmatC(3, 3)    ! Twinning Scmid Tensor in current configuration
        real*8 sigma(3, 3)    ! Stress Tensor (3, 3)
        real*8 tau            ! Twinning Resolved shear stress
        real*8 Lambda_tau     ! Twinning Probability Distribution
        real*8 tau_c          ! Characteristic Resolved shear stress
        real*8 beta           ! Weibull shape parameter exponent
        real*8 P_tau          ! Probability of twin nucleation
        real*8 P_0            ! Threshold Probability of twin nucleation
        real*8 dfgrd(3, 3)    ! Deformation Gradient
        real*8 dfgrdinv(3, 3) ! Deformation Gradient Inverse
        real*8 PtauMax        ! Max probability in an element

        ! Set material parameters
        tau_c = 200.d0      ! TODO: Different critical stress for different types of twin
        beta = 11.1d0
        P_0 = 0.85d0

        ! Stress vector to matrix
        sigma(1, 1) = stress(1)
        sigma(2, 2) = stress(2)
        sigma(3, 3) = stress(3)
        sigma(2, 3) = stress(4)
        sigma(3, 2) = stress(4)
        sigma(1, 3) = stress(5)
        sigma(3, 1) = stress(5)
        sigma(1, 2) = stress(6)
        sigma(2, 1) = stress(6)

        call inverse_3x3(dfgrd, dfgrdinv)

        ! Resolved stress and probabilistic twin nucleation
        do i = 1, 18
          SmatC = matmul(dfgrd, matmul(Smat(:, :, i), dfgrdinv))
          tau = sum(sigma * SmatC)    ! rss
          if (tau.lt.0.d0) then       ! if rss is negative, then continue
            cycle
          end if      ! TODO check for compression twins!

          Lambda_tau = (tau / tau_c)**beta
          P_tau = 1.d0 - exp(-Lambda_tau)
          PtauMax = max(PtauMax, P_tau)
          if (P_tau.gt.P_0) then
            iTwinFlag = i
            exit
          end if

        end do

      end subroutine

      subroutine compute_twin_props(n_tw, props, props_tw, nprops,
     &                              rImat, pi)
        implicit none
        integer nprops
        real*8 n_tw(3)
        real*8 props(nprops)
        real*8 props_tw(nprops)
        real*8 rImat(3, 3)

        integer i, j, ieulererr
        real*8 R_tw(3, 3)
        real*8 Q_0(3, 3), Q_1(3, 3)
        real*8 th, th_tw
        real*8 phi, phi_tw
        real*8 om, om_tw
        real*8 pi

        ! Form rotation matrix corresponding to euler angles
        th = props(31) * pi / 180.d0
        phi = props(32) * pi / 180.d0
        om = props(33) * pi / 180.d0
        call rotmat(th, phi, om, Q_0)

        ! Form rotation matrix corresponding to twinning
        do j = 1, 3
          do i = 1, 3
            R_tw(i, j) = 2*n_tw(i)*n_tw(j) - rImat(i, j)
          end do
        end do

        ! Find new euler angles after rotation of crystal basis due to twinning
        Q_1 = matmul(R_tw, Q_0)
        call eulang(Q_1, th_tw, phi_tw, om_tw, ieulererr)

        ! Form new properties for twinning with updated euler angles
        props_tw(:) = props(:)
        props_tw(31) = th_tw * 180.d0 / pi
        props_tw(32) = phi_tw * 180.d0 / pi
        props_tw(33) = om_tw * 180.d0 / pi
      end subroutine

      ! Explicit sym mat -> vec
      subroutine sym_mat_to_vec(mat, vec)
        implicit none
        real*8 mat(3, 3), vec(6)

        vec(1) = mat(1, 1)
        vec(2) = mat(2, 2)
        vec(3) = mat(3, 3)
        vec(4) = mat(2, 3)
        vec(5) = mat(1, 3)
        vec(6) = mat(1, 2)
      end subroutine

      ! Explicit vec -> sym mat
      subroutine vec_to_sym_mat(vec, mat)
        implicit none
        real*8 vec(6), mat(3, 3)

        mat(1, 1) = vec(1)
        mat(2, 2) = vec(2)
        mat(3, 3) = vec(3)
        mat(2, 3) = vec(4)
        mat(3, 2) = vec(4)
        mat(1, 3) = vec(5)
        mat(3, 1) = vec(5)
        mat(1, 2) = vec(6)
        mat(2, 1) = vec(6)
      end subroutine

      subroutine compute_twin_evol(props, nprops, tau_tw, theta, dt,
     &                             rBWidth, rBLength)
        implicit none
        integer nprops
        real*8 props(nprops), tau_tw, theta, dt, rBWidth, rBLength
        real*8 Gtw, bk, vl, vb, sl, p, q, vldot, vbdot, ramp, fac

        ! Retrieve relevant properties
        Gtw = props(52)
        vl = props(53)
        vb = props(54)
        sl = props(55)
        p = props(56)
        q = props(57)
        bk = props(27)

        ! ramp could be positive or 0
        ramp = 1.d0 - (abs(tau_tw) / sl)**p ! TODO: Check if -ve tau_tw makes sense
        if (ramp.lt.0.d0) then
          ramp = 0.d0
        end if

        ! Compute speed of length and width propagation
        fac = exp(-Gtw / (bk * theta) * ramp**q)
        vldot = vl * fac
        vbdot = vb * fac

        ! Compute new length and width
        rBLength = rBLength + vldot * dt
        rBWidth = rBWidth + vbdot * dt
      end subroutine

      subroutine compute_twin_velocities(props, nprops, tau_tw, theta,
     &                                   l_dot, w_dot, b_dot)
        implicit none
        integer nprops
        real*8 props(nprops), tau_tw, theta, l_dot, w_dot, b_dot
        real*8 Gtw, bk, vl, vb, sl, p, q, fac
        real*8, external :: ramp

        ! Retrieve relevant properties
        Gtw = props(52)
        vl = props(53)
        vb = props(54)
        sl = props(55)
        p = props(56)
        q = props(57)
        bk = props(27)

        ! Compute speed of length and width propagation
        fac = exp(-Gtw/(bk*theta)*ramp(1.d0 - ramp(tau_tw / sl)**p)**q)
        l_dot = vl * fac
        w_dot = l_dot * 0.85d0 ! TODO: New expression for w_dot needed
        b_dot = vb * fac
      end subroutine
      subroutine find_next_twin_elem(elmt_number, elmt_conn, num_elems,
     &                        nodal_coords, num_nodes, st_coord, rn, rm,
     &                        elmt_neighbour, tw_length, nxt_elmt_num)
        implicit none
        integer elmt_number, num_elems, num_nodes, nxt_elmt_num
        integer elmt_conn(8, num_elems), elmt_neighbour(6, num_elems)
        real*8 nodal_coords(3, num_nodes), st_coord(3), rn(3), rm(3)

        ! Define 6 faces:
        ! Face1: 1,2,3,4; Face2: 5,8,7,6; Face3: 1,5,6,2;
        ! Face4: 2,6,7,3; Face5: 3,7,8,4; Face6: 4,5,8,1;
        integer FACE(6,4)
        data FACE /1,5,1,2,3,4,
     1             2,8,5,6,7,8,
     2             3,7,6,7,8,5,
     3             4,6,2,3,4,1/
        real*8 rs(3), face_node_coord(3, 5), rA, rb, t
        real*8 edge_int_coord(3, 2), rA2, rb2, t2, final_int_coord(3, 2)
        integer ifc, i, face_nodes(4), cntr1, cntr2, neigh_elmt_no
        real*8 dir_dp, tw_dist, tw_dist_max, tw_length
        integer tw_prop_face, tw_prop_face_nodes(4)

        ! Compute rs = rm X rn, rm = rn x rs
        call cross_product(rm, rn, rs)

        ! Loop over element faces
        cntr1 = 0
        tw_dist_max = 0.d0
        tw_prop_face = -1
        do ifc = 1, 6
          ! Retrieve node numbers on the nodes on the face
          face_nodes(1) = elmt_conn(Face(ifc, 1), elmt_number)
          face_nodes(2) = elmt_conn(Face(ifc, 2), elmt_number)
          face_nodes(3) = elmt_conn(Face(ifc, 3), elmt_number)
          face_nodes(4) = elmt_conn(Face(ifc, 4), elmt_number)

          ! retrieve nodal coordinates of 4 nodes on face and compute face centroid
          do i = 1, 4
            face_node_coord(1:3, i) = nodal_coords(1:3, face_nodes(i))
          end do
          face_node_coord(1:3, 5) = nodal_coords(1:3, face_nodes(1))

          ! Compute edge intersection parameter t
          cntr2 = 0
          do i = 1, 4
            rA = dot_product(face_node_coord(:, i+1) -
     &                       face_node_coord(:, i), rn)
            rb = dot_product(st_coord - face_node_coord(:, i), rn)
            if (abs(rA).lt.1d-6) then
              t = 100.d0
            else
              t = rb / rA
            end if

            ! If 0 <= t <= 1, then twin plane really intersects this edge
            if ((t.gt.0.d0).and.(t.lt.1.d0)) then
              cntr2 = cntr2 + 1
              edge_int_coord(:, cntr2) = face_node_coord(:, i) +
     &             t * (face_node_coord(:, i+1) - face_node_coord(:, i))
            end if
          end do

          ! If intersected on two edges
          !   compute intersection parameter of twin direction on line
          !     joining the intersection points of the edges
          if (cntr2.eq.2) then
            rA2 = dot_product(edge_int_coord(:, 2) -
     &                        edge_int_coord(:, 1), rs)
            rb2 = dot_product(st_coord - edge_int_coord(:, 1), rs)
            if (abs(rA2).lt.1d-6) then
              t2 = 100.d0
            else
              t2 = rb2 / rA2
            end if

            ! Compute coordinates of twin propagation direction itersection coordinates
            if ((t2.ge.0.d0).and.(t2.le.1.d0)) then
              cntr1 = cntr1 + 1
              final_int_coord(:, cntr1) = edge_int_coord(:, 1)
     &              + t2 * (edge_int_coord(:, 2) - edge_int_coord(:, 1))
              dir_dp = dot_product(final_int_coord(:, cntr1) - st_coord,
     &                             rm)
              tw_dist =sign(norm2(final_int_coord(:, cntr1) - st_coord),
     &                          dir_dp)
              if (tw_dist.gt.tw_dist_max) then
                tw_dist_max = tw_dist
                tw_prop_face = ifc
              end if
            end if
          end if
        end do

        if (tw_prop_face.eq.-1) then
          CALL XPLB_ABQERR(-3,'The Twin does not pass through '//
     &                     'the Element No. %I!',elmt_number,0.d0,' ')
          return
        end if

        tw_length = norm2(final_int_coord(:, 1) - final_int_coord(:, 2))

        tw_prop_face_nodes(1) = elmt_conn(Face(tw_prop_face, 1),
     &                                    elmt_number)
        tw_prop_face_nodes(2) = elmt_conn(Face(tw_prop_face, 2),
     &                                    elmt_number)
        tw_prop_face_nodes(3) = elmt_conn(Face(tw_prop_face, 3),
     &                                    elmt_number)
        tw_prop_face_nodes(4) = elmt_conn(Face(tw_prop_face, 4),
     &                                    elmt_number)
        call quicksort(tw_prop_face_nodes, 1, 4)

        nxt_elmt_num = 0
        do i = 1, 6
          neigh_elmt_no = elmt_neighbour(i, elmt_number)
          if (neigh_elmt_no.eq.0) then
            cycle
          end if

          do ifc = 1, 6
            face_nodes(1) = elmt_conn(Face(ifc, 1), neigh_elmt_no)
            face_nodes(2) = elmt_conn(Face(ifc, 2), neigh_elmt_no)
            face_nodes(3) = elmt_conn(Face(ifc, 3), neigh_elmt_no)
            face_nodes(4) = elmt_conn(Face(ifc, 4), neigh_elmt_no)
            call quicksort(face_nodes, 1, 4)

            if (all(face_nodes.eq.tw_prop_face_nodes)) then
              nxt_elmt_num = neigh_elmt_no
            end if
          end do
        end do

      end subroutine

      subroutine cross_product(a, b, c)
        implicit none
        real*8 a(3), b(3), c(3)

        c(1) = a(2) * b(3) - b(2) * a(3)
        c(2) = b(1) * a(3) - a(1) * b(3)
        c(3) = a(1) * b(2) - b(1) * a(2)

      end subroutine

! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
      recursive subroutine quicksort(a, first, last)
        implicit none
        integer  a(*), x, t
        integer first, last
        integer i, j

        x = a( (first+last) / 2 )
        i = first
        j = last
        do
          do while (a(i) < x)
              i=i+1
          end do
          do while (x < a(j))
              j=j-1
          end do
          if (i >= j) exit
          t = a(i);  a(i) = a(j);  a(j) = t
          i=i+1
          j=j-1
        end do
        if (first < i-1) call quicksort(a, first, i-1)
        if (j+1 < last)  call quicksort(a, j+1, last)
      end subroutine quicksort

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Detects Twinning in the neighbouring elements of iel             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine detect_twinning(iel, elmt_conn, num_elems,
     &           nodal_coords, num_nodes, num_dis_pts, n_tw, m_tw, s_tw,
     &           rad, x_0, l_tw, w_tw, b_tw, xi, twinned, bexp)
        implicit none
        ! Input variables
        integer iel, num_dis_pts, num_elems, num_nodes
        integer elmt_conn(8, num_elems)
        real*8 nodal_coords(3, num_nodes), n_tw(3), m_tw(3), s_tw(3)
        real*8 rad, x_0(3), l_tw, w_tw, b_tw, xi(num_dis_pts)
        logical twinned, bexp

        ! Subroutine used variables
        integer i, j, k, ind
        logical proj_twinned, b_twinned
        real*8 x(3), N(8), b_proj, l_proj, w_proj, proj_measure

        !! Code
        twinned = .false.
        bexp = .false.
        ! Loop over internal discrete points
        do k = 2, num_dis_pts - 1
          do j = 2, num_dis_pts - 1
            do i = 2, num_dis_pts - 1
              ! Set thickness and projection flags to false
              b_twinned = .false.
              proj_twinned = .false.

              ! Compute shape function and form position vector
              call shp_fcn_8n(xi(i), xi(j), xi(k), N)
              x = 0.d0
              do ind = 1, 8
                x = x + nodal_coords(:, elmt_conn(ind, iel)) * N(ind)
              end do

              ! Compute projections along length, width and thickness
              b_proj = dot_product(x - x_0, n_tw)
              l_proj = dot_product(x - x_0 - (b_proj * n_tw), m_tw)
              w_proj = dot_product(x - x_0 - (b_proj * n_tw), s_tw)

              ! thickness decision based on distance between internal points
              if ((abs(b_proj).lt.rad).or.(abs(b_proj).lt.b_tw)) then
                b_twinned = .true.
              end if

              ! Projection measure based on elliptical twin and decide
              ! whether the point is inside
              proj_measure = (l_proj/l_tw)**2 + (w_proj/w_tw)**2 - 1.d0
              if (proj_measure.lt.0.d0) then
                proj_twinned = .true.
              end if

              ! If both thickness and projection satisfy intersection, return true
              if (b_twinned.and.proj_twinned) then
                twinned = .true.
                if (abs(b_proj).lt.b_tw) then
                  bexp = .true.
                end if
                return
              end if
            end do
          end do
        end do
      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Computes shape functions of linear hex element                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine shp_fcn_8n(xi, eta, zeta, N)
        implicit none
        real*8 xi, eta, zeta, N(8)

        N(1) = 0.125d0 * (1 - xi) * (1 - eta) * (1 - zeta)
        N(2) = 0.125d0 * (1 + xi) * (1 - eta) * (1 - zeta)
        N(3) = 0.125d0 * (1 + xi) * (1 + eta) * (1 - zeta)
        N(4) = 0.125d0 * (1 - xi) * (1 + eta) * (1 - zeta)
        N(5) = 0.125d0 * (1 - xi) * (1 - eta) * (1 + zeta)
        N(6) = 0.125d0 * (1 + xi) * (1 - eta) * (1 + zeta)
        N(7) = 0.125d0 * (1 + xi) * (1 + eta) * (1 + zeta)
        N(8) = 0.125d0 * (1 - xi) * (1 + eta) * (1 + zeta)
      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Ramp Function                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function ramp(x)
        implicit none
        real*8 x, ramp

        if (x.lt.0.d0) then
          ramp = 0.d0
        else
          ramp = x
        end if
      end function