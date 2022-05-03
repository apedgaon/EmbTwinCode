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

      subroutine comp_twin_dfgrd(rImat, Smat, F_tw)
          implicit none
          real*8 rn(3, 18), rm(3, 18)
          real*8 gamma1, gamma2, gamma3
          real*8 Smat(3, 3, 18), F_tw(3, 3, 18), rImat(3, 3)
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

          gamma1 = 0.17475083284085197d0
          gamma2 = 0.218d0 ! Need more accuracy
          gamma3 = 0.098d0 ! Need more accuracy
          Smat = 0.0d0
          do i = 1, 3
            do j = 1, 3
                Smat(i, j, :) = rm(i, :) * rn(j, :)
            end do
          end do

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
     &                      lEditSvars, rW)
          
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
          real*8 tmphvars(nsvars), rW(3,3,3,3)

C         for umat stress to vuel stress
          real*8 rUmatStress(6)
          
          ndi = 3
          nshr = 3
          ntens = ndi + nshr
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
     &               kstep,   kinc, rW)
          else    ! Required during iterative updates when hvars should not be edited
            tmphvars = hvars
            call umat_crys(rUmatStress,  tmphvars,  ddsdde,  sse,  spd,
     &               scd,     rpl,     ddsddt,  drplde,  drpldt,
     &               stran,  dstran, time,    dtime,   temp,
     &               dtemp,   predef,  dpred,   cmname,  ndi,
     &               nshr,    ntens,   nsvars,  props,   nprops,
     &               rGsCoords,  drot,    pnewdt,  celent,  rDfgrd0,
     &               rDfgrd1,  noel,    npt,     layer,   kspt,
     &               kstep,   kinc, rW)
          end if

          rStress(1) = rUmatStress(1)
          rStress(2) = rUmatStress(2)
          rStress(3) = rUmatStress(3)
          rStress(4) = rUmatStress(6)
          rStress(5) = rUmatStress(5)
          rStress(6) = rUmatStress(4)
      end subroutine
      
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

      end subroutine

      subroutine check_twin_nucleation(stress, Smat, dfgrd, iTwinFlag)
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

        iTwinFlag = 0

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

          Lambda_tau = (abs(tau) / tau_c)**beta
          P_tau = 1.d0 - exp(-Lambda_tau)
          if (P_tau.gt.P_0) then
            iTwinFlag = i
            exit
          end if

        end do

      end subroutine
