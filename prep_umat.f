      subroutine dfgrd_prep(iTIntPt, rRefCoords, rCurrCoords, rPrevCoords,
     &                      rXiEtaPt, mcrd, nnode, rDfgrd0, rDfgrd1,
     &                      rGsCoords)
          implicit none
          integer iTIntPt, mcrd, nnode
          real*8 rRefCoords(mcrd, nnode), rCurrCoords(mcrd, nnode)
          real*8 rPrevCoords(mcrd, nnode), rXiEtaPt(3, iTIntPt)
          real*8 rDfgrd0(3, 3, iTIntPt), rDfgrd1(3, 3, iTIntPt)
          real*8 rGsCoords(3, iTIntPt)

          integer l, idx, jdx
          real*8 rN(8), rdNdXRef(3, iTIntPt), rJdet
          real*8 rDetF0Bar, rDetF1Bar, rSelRedIntConst
          real*8 rDetF0(iTIntPt), rDetF1(iTIntPt)
          real*8, external :: determinant

          rDetF0Bar = 0.0
          rDetF1Bar = 0.0
          do l = 1, iTIntPt
               call SHP_FUNC_3D_8N(rXiEtaPt(:, l), rRefCoords, 1, rN, rdNdXRef,
     &                        rJdet)
               call GdefGrad_3D(rdNdXRef, rPrevCoords, mcrd, nnode,
     &                     1, rDfgrd0(:, :, l))
               call GdefGrad_3D(rdNdXRef, rCurrCoords, mcrd, nnode,
     &                     1, rDfgrd1(:, :, l))
               rDetF0(l) = determinant(rDfgrd0(:,:,l))
               rDetF1(l) = determinant(rDfgrd1(:,:,l))
               rDetF0Bar = rDetF0Bar + rDetF0(l)
               rDetF1Bar = rDetF1Bar + rDetF1(l)

               rGsCoords(:, l) = 0.0
               do idx = 1, nnode
                   do jdx = 1, mcrd
                       rGsCoords(jdx, l) = rGsCoords(jdx, l) +
     &                        rN(idx) * rCurrCoords(jdx, idx)
                   end do
               end do
          end do

          rDetF0Bar = rDetF0Bar / iTIntPt
          rDetF1Bar = rDetF1Bar / iTIntPt

          do l = 1, iTIntPt
               rSelRedIntConst = (rDetF0Bar / rDetF0(l))**(1.0/3)
               rDfgrd0(:, :, l) = rDfgrd0(:, :, l) * rSelRedIntConst
               rSelRedIntConst = (rDetF1Bar / rDetF1(l))**(1.0/3)
               rDfgrd1(:, :, l) = rDfgrd1(:, :, l) * rSelRedIntConst
          end do

      end subroutine

      subroutine umat_stress(rStress, hvars, rRefCoords, rCurrCoords,
     &                      rPrevCoords, props, rXiEtaPt, nprops,
     &                      nsvars, nnode, noel, npt, mcrd, time,
     &                      dtime, kstep, kinc)
          
          implicit none
          ! input arguments declaration
          integer nprops, nsvars, nnode, noel, npt, mcrd
          integer kstep, kinc
          real*8 rStress(6), hvars(nsvars), props(nprops)
          real*8 rRefCoords(mcrd, nnode), rCurrCoords(mcrd, nnode)
          real*8 rPrevCoords(mcrd, nnode), time(2), rXiEtaPt(3)
          real*8 dtime

          ! required to be computed for umat
          real*8 rGsCoords(3), rN(8), rdNdXRef(3, 8), rJdet
          real*8 rF0(3, 3), rF1(3, 3)
          integer idx, jdx

          ! additional arguments for umat
          real*8 ddsdde(6, 6), ddsddt(6), drplde(6), stran(6)
          real*8 dstran(6), drot(3, 3)
          real*8 sse, spd, scd, rpl, drpldt, temp, dtemp
          real*8 predef, dpred, pnewdt, celent
          character*80 cmname
          integer ndi, nshr, ntens, layer, kspt
          
          call SHP_FUNC_3D_8N(rXiEtaPt, rRefCoords, 1, rN, rdNdXRef,
     &                        rJdet)
          call GdefGrad_3D(rdNdXRef, rPrevCoords, mcrd, nnode, 
     &                     1, rF0)
          call GdefGrad_3D(rdNdXRef, rCurrCoords, mcrd, nnode, 
     &                     1, rF1)
          
          rGsCoords(1) = 0.0
          rGsCoords(2) = 0.0
          rGsCoords(3) = 0.0
          do idx = 1, nnode
              do jdx = 1, mcrd
                  rGsCoords(jdx) = rGsCoords(jdx) +
     &                   rN(idx) * rCurrCoords(jdx, idx)
              end do
          end do

          ndi = 3
          nshr = 3
          ntens = ndi + nshr
          call umat_crys(rStress,  hvars,  ddsdde,  sse,     spd,
     &               scd,     rpl,     ddsddt,  drplde,  drpldt,
     &               stran,  dstran, time,    dtime,   temp,
     &               dtemp,   predef,  dpred,   cmname,  ndi,
     &               nshr,    ntens,   nsvars,  props,   nprops,
     &               rGsCoords,  drot,    pnewdt,  celent,  rF0, 
     &               rF1,  noel,    npt,     layer,   kspt,
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

      subroutine umat_crys(stress,  statev,  ddsdde,  sse,     spd,
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
