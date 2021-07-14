      subroutine dfgrd_prep(iTIntPt, rRefCoords, rCurrCoords, rPrevCoords,
     &                      rXiEtaPt, mcrd, nnode, rDfgrd0, rDfgrd1,
     &                      rGsCoords)
C         This subroutine prepares deformation gradient and gauss point
C         coordinates for UMAT type function.
          implicit none

C         Integer inputs
          integer iTIntPt, mcrd, nnode

C         Reference, Current and Previous Coordinates
          real*8 rRefCoords(mcrd, nnode), rCurrCoords(mcrd, nnode)
          real*8 rPrevCoords(mcrd, nnode), rXiEtaPt(3, iTIntPt)

C         Initalize Deformation Gradient at each Gass Pt
C         and Gauss Pt Coordinates
          real*8 rDfgrd0(3, 3, iTIntPt), rDfgrd1(3, 3, iTIntPt)
          real*8 rGsCoords(3, iTIntPt)

C         Additional variables for computations
          integer l, idx, jdx
          real*8 rN(8), rdNdXRef(3, iTIntPt), rJdet
          real*8 rDetF0Bar, rDetF1Bar, rSelRedIntConst
          real*8 rDetF0(iTIntPt), rDetF1(iTIntPt)
          real*8, external :: determinant ! external determinant function

C         Compute Deformation Gradients and Gauss Pt Coordinates
          rDetF0Bar = 0.0
          rDetF1Bar = 0.0
          do l = 1, iTIntPt
C              Shape Functions and Derivatives w.r.t Reference Coordinates
               call SHP_FUNC_3D_8N(rXiEtaPt(:, l), rRefCoords, 1, rN, rdNdXRef,
     &                        rJdet)

C              Deformation Gradients of previous and current coords
C              and their determinants
               call GdefGrad_3D(rdNdXRef, rPrevCoords, mcrd, nnode,
     &                     1, rDfgrd0(:, :, l))
               call GdefGrad_3D(rdNdXRef, rCurrCoords, mcrd, nnode,
     &                     1, rDfgrd1(:, :, l))
               rDetF0(l) = determinant(rDfgrd0(:,:,l))
               rDetF1(l) = determinant(rDfgrd1(:,:,l))
               rDetF0Bar = rDetF0Bar + rDetF0(l)
               rDetF1Bar = rDetF1Bar + rDetF1(l)

C              Compute Guass Pt Coordinates
               rGsCoords(:, l) = 0.0
               do idx = 1, nnode
                   do jdx = 1, mcrd
                       rGsCoords(jdx, l) = rGsCoords(jdx, l) +
     &                        rN(idx) * rCurrCoords(jdx, idx)
                   end do
               end do
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
     &                      rDfgrd1, props, rXiEtaPt, nprops,
     &                      nsvars, nnode, noel, npt, mcrd, time,
     &                      dtime, kstep, kinc)
          
          implicit none
C         input arguments declaration
          integer nprops, nsvars, nnode, noel, npt, mcrd
          integer kstep, kinc
          real*8 rStress(6), hvars(nsvars), props(nprops)
          real*8 rGsCoords(3), rDfgrd0(3, 3), rDfgrd1(3, 3)
          real*8 time(2), rXiEtaPt(3)
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
          call umat_crys(rStress,  hvars,  ddsdde,  sse,     spd,
     &               scd,     rpl,     ddsddt,  drplde,  drpldt,
     &               stran,  dstran, time,    dtime,   temp,
     &               dtemp,   predef,  dpred,   cmname,  ndi,
     &               nshr,    ntens,   nsvars,  props,   nprops,
     &               rGsCoords,  drot,    pnewdt,  celent,  rDfgrd0, 
     &               rDfgrd1,  noel,    npt,     layer,   kspt,
     &               kstep,   kinc)


      end subroutine
      
      subroutine pumat_crys(rStrain,rSpin,svars,nstatev,props,
     &                       nprops,iVelmType,
     &                       iVforType,iConsType,rTime,rDt,
     &                       rStress,rSigInc,rCmod,rGtan,rCd,iCvFlag)
     
          dimension rStrain(6), rSpin(3), svars(nstatev), props(nprops)
          dimension rTime(2)
          dimension rCmod(6,6), rGtan(6,3), rStress(6)
          dimension rSigInc(6)    ! stress increment
          
          ntens = 6

          call umat_crys(rStress,  svars,  ddsdde,  sse,     spd,
     &               scd,     rpl,     ddsddt,  drplde,  drpldt,
     &               stran,  dstran, rTime,    rDt,   temp,
     &               dtemp,   predef,  dpred,   cmname,  ndi,
     &               nshr,    ntens,   nstatev,  props,   nprops,
     &               coords,  drot,    pnewdt,  celent,  dfgrd0,
     &               dfgrd1,  noel,    npt,     layer,   kspt,
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
      
      include 'vaba_param.inc'
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
