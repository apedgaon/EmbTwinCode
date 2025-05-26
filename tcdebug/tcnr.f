      program tcnr
        implicit none
        real*8 :: rntwc(3), rugrd(3, 3), hbratio
        real*8 :: rFu0(3, 3), rFz0(3, 3), rFz0h(3, 3), del(3, 3)
        real*8 :: hvarsU(148), hvarsZ(148), props(58)
        real*8 :: props_tw(58), Fe(3, 3), time(2), dtimeCur
        real*8 :: F_tw_inv(3, 3), rmh(3), rFu1(3, 3), rFz1(3, 3)
        real*8 :: rFz1h(3, 3), rmtwc(3), rmtw(3)
        integer :: nsvars, nmprops, nnode, elmt_number, l
        integer :: mcrd, kstep, kinc, iTwFlg, iTwEmbStp
        integer i, j
        real*8 :: rtemp, rmhMag, rElmtSize, rBWidth, rntw(3)
        real*8 R_tw(3, 3)
        logical lNRerr

        open(10, file='../debug_vars.txt', status='old', action='read')
        read(10, *) rntwc
        read(10, *) rugrd
        read(10, *) hbratio
        read(10, *) rFu0
        read(10, *) rFz0
        read(10, *) rFz0h
        read(10, *) del
        read(10, *) rtemp
        nsvars = int(rtemp)
        read(10, *) hvarsU
        read(10, *) hvarsZ
        read(10, *) rtemp
        nmprops = int(rtemp)
        read(10, *) props
        read(10, *) props_tw
        read(10, *) rtemp
        nnode = int(rtemp)
        read(10, *) rtemp
        elmt_number = int(rtemp)
        read(10, *) rtemp
        l = int(rtemp)
        read(10, *) rtemp
        mcrd = int(rtemp)
        read(10, *) time
        read(10, *) dtimeCur
        read(10, *) rtemp
        kstep = int(rtemp)
        read(10, *) rtemp
        kinc = int(rtemp)
        read(10, *) F_tw_inv
        read(10, *) rmtwc
        read(10, *) rmtw
        read(10, *) rtemp
        iTwFlg = int(rtemp)
        read(10, *) rtemp
        iTwEmbStp = int(rtemp)
        read(10, *) rElmtSize
        read(10, *) rBWidth
        read(10, *) rntw

        close(10)

        print *, 'rntwc = ', rntwc
        print *, 'rugrd = ', rugrd
        print *, 'hbratio = ', hbratio
        print *, 'rFu0 = ', rFu0
        print *, 'rFz0 = ', rFz0
        print *, 'rFz0h = ', rFz0h
        print *, 'del = ', del
        print *, 'nsvars = ', nsvars
        print *, 'hvarsU = ', hvarsU
        print *, 'hvarsZ = ', hvarsZ
        print *, 'nmprops = ', nmprops
        print *, 'props = ', props
        print *, 'props_tw = ', props_tw
        print *, 'nnode = ', nnode
        print *, 'elmt_number = ', elmt_number
        print *, 'l = ', l
        print *, 'mcrd = ', mcrd
        print *, 'time = ', time
        print *, 'dtimeCur = ', dtimeCur
        print *, 'kstep = ', kstep
        print *, 'kinc = ', kinc
        print *, 'F_tw_inv = ', F_tw_inv
        print *, 'rmtwc = ', rmtwc
        print *, 'rmtw = ', rmtw
        print *, 'iTwFlg = ', iTwFlg
        print *, 'iTwEmbStp = ', iTwEmbStp
        print *, 'rElmtSize = ', rElmtSize
        print *, 'rBWidth = ', rBWidth
        print *, 'rntw = ', rntw

    !     rntw = (/0.7327380113307633d0, 0.5330353418564765d0,
    !  &           0.4230464880872873d0/)

        do j = 1, 3
          do i = 1, 3
            R_tw(i, j) = 2 * rntw(i) * rntw(j) - del(i, j)
          end do
        end do

        hvarsZ(14) = hvarsU(14)
        hvarsZ(15) = hvarsU(15)
        hvarsZ(16) = hvarsU(16)
        hvarsZ(17) = hvarsU(17)
        hvarsZ(18) = hvarsU(18)
        hvarsZ(19) = hvarsU(19)

        call rotate_twin_svars(hvarsZ, nsvars, R_tw)
        print *, 'hvarsZ = ', hvarsZ

        rmhMag = 0.218d0 / (hbratio - 1.d0)
        rmh = rmtwc * rmhMag
        print*, 'rmh = ', rmh

        call traction_continuityNR(rntwc, rugrd, hbratio, rFu0, rFz0,
     &                             rFz0h, del, nsvars, hvarsU, hvarsZ,
     &                             nmprops, props, props_tw, Fe, nnode,
     &                             elmt_number,l, mcrd, time, dtimeCur,
     &                             kstep, kinc, F_tw_inv, rmh, rFu1,
     &                             rFz1, rFz1h, lNRerr)

        print *, 'rmh = ', rmh
        print *, 'rFu1 = ', rFu1
        print *, 'rFz1 = ', rFz1
        print *, 'rFz1h = ', rFz1h
        print *, 'lNRerr = ', lNRerr

      end program tcnr

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
          if (kinc.eq.0)  dtime = 1.d-9   ! handle 0 dtimeCur in the begining of explicit

          if (lEditSvars) then    ! During regular update
            call umat_crys(rUmatStress,  hvars,  ddsdde,  sse,  spd,
     &               scd,     rpl,     ddsddt,  drplde,  drpldt,
     &               stran,  dstran, time,    dtime,   temp,
     &               dtemp,   predef,  dpred,   cmname,  ndi,
     &               nshr,    ntens,   nsvars,  props,   nprops,
     &               rGsCoords,  drot,    pnewdt,  celent,  rDfgrd0,
     &               rDfgrd1,  noel,    npt,     layer,   kspt,
     &               kstep,   kinc, rW, iterErr, Fe, lEditSvars)
          else    ! Required during iterative updates when hvars should not be edited
            tmphvars = hvars
            call umat_crys(rUmatStress,  tmphvars,  ddsdde,  sse,  spd,
     &               scd,     rpl,     ddsddt,  drplde,  drpldt,
     &               stran,  dstran, time,    dtime,   temp,
     &               dtemp,   predef,  dpred,   cmname,  ndi,
     &               nshr,    ntens,   nsvars,  props,   nprops,
     &               rGsCoords,  drot,    pnewdt,  celent,  rDfgrd0,
     &               rDfgrd1,  noel,    npt,     layer,   kspt,
     &               kstep,   kinc, rW, iterErr, Fe, lEditSvars)

            ! Update init spk2 guess
            hvars(14) = tmphvars(14)
            hvars(15) = tmphvars(15)
            hvars(16) = tmphvars(16)
            hvars(17) = tmphvars(17)
            hvars(18) = tmphvars(18)
            hvars(19) = tmphvars(19)
          end if

          if (kinc.eq.0)  dtime = 0.d0 ! Set the vuel dtime back to 0 for initial increment

          rStress(1) = rUmatStress(1)
          rStress(2) = rUmatStress(2)
          rStress(3) = rUmatStress(3)
          rStress(4) = rUmatStress(6)
          rStress(5) = rUmatStress(5)
          rStress(6) = rUmatStress(4)
      end subroutine

      function ramp(x)
        implicit none
        real*8 x, ramp

        if (x.lt.0.d0) then
          ramp = 0.d0
        else
          ramp = x
        end if
      end function

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

      subroutine rotate_twin_svars(statev, nstatev, R_tw)
        implicit none
        ! Argument variables
        integer nstatev
        real*8 statev(nstatev), R_tw(3, 3)

        ! Local variables
        integer n, i, j
        real*8 Fpinv_n(3, 3), Spk2_n(3, 3)

        n = 4
        do j = 1, 3
          do i = 1, 3
            n = n + 1
            Fpinv_n(i, j) = statev(n)
          end do
        end do

        Spk2_n(1, 1) = statev(14)
        Spk2_n(2, 2) = statev(15)
        Spk2_n(3, 3) = statev(16)
        Spk2_n(1, 2) = statev(17)
        Spk2_n(1, 3) = statev(18)
        Spk2_n(2, 3) = statev(19)
        Spk2_n(2, 1) = statev(17)
        Spk2_n(3, 1) = statev(18)
        Spk2_n(3, 2) = statev(19)

        ! Fpinv_n = matmul(R_tw, matmul(Fpinv_n, R_tw))
        Spk2_n = matmul(R_tw, matmul(Spk2_n, R_tw))

        n = 4
        do j = 1, 3
          do i = 1, 3
            n = n + 1
            statev(n) = Fpinv_n(i, j)
          end do
        end do

        statev(14) = Spk2_n(1, 1)
        statev(15) = Spk2_n(2, 2)
        statev(16) = Spk2_n(3, 3)
        statev(17) = Spk2_n(1, 2)
        statev(18) = Spk2_n(1, 3)
        statev(19) = Spk2_n(2, 3)

      end subroutine
