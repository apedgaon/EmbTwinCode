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
          real*8 nt(3, 6), mt(3, 6), nc(3, 6), mc(3, 6)
          real*8 gammat, gammac
          real*8 Smat(3, 3, 12), F_tw(3, 3, 12), rImat(3, 3)
          integer i, j, k

          nt(:, 1) = (/0.5850518001199482d0, 0.3377798142891272d0,
     &                 0.7373053561688068d0/)
          mt(:, 1) = (/-0.6385251687885203d0, -0.3686526780844034d0,
     &                 0.6755596285782544d0/)
          nt(:, 2) = (/0.0000000000000000d0, 0.6755596285782544d0,
     &                 0.7373053561688068d0/)
          mt(:, 2) = (/0.0000000000000000d0, -0.7373053561688068d0,
     &                 0.6755596285782544d0/)
          nt(:, 3) = (/-0.5850518001199482d0, 0.3377798142891272d0,
     &                 0.7373053561688068d0/)
          mt(:, 3) = (/0.6385251687885203d0, -0.3686526780844034d0,
     &                 0.6755596285782544d0/)
          nt(:, 4) = (/-0.5850518001199482d0, -0.3377798142891272d0,
     &                  0.7373053561688068d0/)
          mt(:, 4) = (/0.6385251687885203d0, 0.3686526780844034d0,
     &                 0.6755596285782544d0/)
          nt(:, 5) = (/0.0000000000000000d0, -0.6755596285782544d0,
     &                 0.7373053561688068d0/)
          mt(:, 5) = (/0.0000000000000000d0, 0.7373053561688068d0,
     &                 0.6755596285782544d0/)
          nt(:, 6) = (/0.5850518001199482d0, -0.3377798142891272d0,
     &                 0.7373053561688068d0/)
          mt(:, 6) = (/-0.6385251687885203d0, 0.3686526780844034d0,
     &                 0.6755596285782544d0/)
          nc(:, 1) = (/0.4230229164297242d0, 0.7326971840222455d0,
     &                 0.5331101656329227d0/)
          mc(:, 1) = (/0.10277731996134944d0, 0.17801554003882017d0,
     &                -0.9786456406719694d0/)
          nc(:, 2) = (/-0.4230229164297242d0, 0.7326971840222455d0,
     &                 0.5331101656329227d0/)
          mc(:, 2) = (/-0.10277731996134944d0, 0.17801554003882017d0,
     &                -0.9786456406719694d0/)
          nc(:, 3) = (/-0.8460458328594485d0, 0.0000000000000000d0,
     &                 0.5331101656329227d0/)
          mc(:, 3) = (/-0.20555463992269887d0, 0.00000000000000000d0,
     &                -0.9786456406719694d0/)
          nc(:, 4) = (/-0.4230229164297242d0, -0.7326971840222455d0,
     &                 0.5331101656329227d0/)
          mc(:, 4) = (/-0.10277731996134944d0, -0.17801554003882017d0,
     &                -0.9786456406719694d0/)
          nc(:, 5) = (/0.4230229164297242d0, -0.7326971840222455d0,
     &                 0.5331101656329227d0/)
          mc(:, 5) = (/0.10277731996134944d0, -0.17801554003882017d0,
     &                -0.9786456406719694d0/)
          nc(:, 6) = (/0.8460458328594485d0, 0.0000000000000000d0,
     &                 0.5331101656329227d0/)
          mc(:, 6) = (/0.20555463992269887d0, 0.00000000000000000d0,
     &                -0.9786456406719694d0/)

          gammat = 0.216d0
          gammac = 0.177d0
          Smat = 0.0d0
          do i = 1, 3
            do j = 1, 3
                Smat(i, j, 1:6) = mt(i, :) * nt(j, :)
                Smat(i, j, 7:12) = mc(i, :) * nc(j, :)
            end do
          end do

          do i = 1, 6
            F_tw(:, :, i) = rImat + gammat * Smat(:, :, i)
          end do
          do i = 7, 12
            F_tw(:, :, i) = rImat + gammac * Smat(:, :, i)
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
     &                      npt, mcrd, time, dtime, kstep, kinc)
          
          implicit none
C         input arguments declaration
          integer nprops, nsvars, nnode, noel, npt, mcrd
          integer kstep, kinc
          real*8 rStress(6), hvars(nsvars), props(nprops)
          real*8 rGsCoords(3), rDfgrd0(3, 3), rDfgrd1(3, 3)
          real*8 time(2)
          real*8 dtime
          integer idx, jdx

C         additional arguments for umat
          real*8 ddsdde(6, 6), ddsddt(6), drplde(6), stran(6)
          real*8 dstran(6), drot(3, 3)
          real*8 sse, spd, scd, rpl, drpldt, temp, dtemp
          real*8 predef, dpred, pnewdt, celent
          character*80 cmname
          integer ndi, nshr, ntens, layer, kspt
          
          ndi = 3
          nshr = 3
          ntens = ndi + nshr
          if (((kstep.eq.0).or.(kstep.eq.1)).and.(kinc.eq.0)) then
              dtime = 1.d-9
          end if

          call umat_crys2(rStress,  hvars,  ddsdde,  sse,     spd,
     &               scd,     rpl,     ddsddt,  drplde,  drpldt,
     &               stran,  dstran, time,    dtime,   temp,
     &               dtemp,   predef,  dpred,   cmname,  ndi,
     &               nshr,    ntens,   nsvars,  props,   nprops,
     &               rGsCoords,  drot,    pnewdt,  celent,  rDfgrd0, 
     &               rDfgrd1,  noel,    npt,     layer,   kspt,
     &               kstep,   kinc)


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


      end subroutine
