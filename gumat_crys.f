      subroutine tumat_crys(rStress, svars, rOrigCoords, rCurrCoords,
     &                      props, rXiEtaPt, nprops, nstatev, nnode,
     &                      mcrd, time, dtime)
          
          include 'vaba_param.inc'

          dimension rStress(6), svars(nstatev), props(nprops)
          dimension rOrigCoords(nnode, mcrd), rCurrCoords(mcrd, nnode)
          dimension time(2), rXiEtaPt(3)
          
          real*8 rRefCoords(mcrd, nnode), rGsCoords(3)
          real*8 rN(8),rdNdXRef(3, 8), rJdet, rF(3, 3)
          
          rGsCoords(1) = 0.0
          rGsCoords(2) = 0.0
          rGsCoords(3) = 0.0
          ntens = 6
          do idx = 1, nnode
              do jdx = 1, mcrd 
                  rRefCoords(jdx, idx) = rOrigcoords(idx, jdx)
                  rGsCoords(jdx) = rGsCoords(jdx) + 
     &             rN(idx) * rCurrCoords(jdx, idx)
              end do
          end do

          call SHP_FUNC_3D_8N(rXiEtaPt,rRefCoords, 1, rN, rdNdXRef,
     &                        rJdet)
          call GdefGrad_3D(rdNdXRef, rCurrCoords, mcrd, nnode, 
     &                     iVelmType, rF)

          call umat_crys(rStress,  svars,  ddsdde,  sse,     spd,
     &               scd,     rpl,     ddsddt,  drplde,  drpldt,
     &               stran,  dstran, time,    dtime,   temp,
     &               dtemp,   predef,  dpred,   cmname,  ndi,
     &               nshr,    ntens,   nstatev,  props,   nprops,
     &               rGsCoords,  drot,    pnewdt,  celent,  dfgrd0,
     &               dfgrd1,  noel,    npt,     layer,   kspt,
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
