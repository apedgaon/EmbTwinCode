C     Traction Continuity Newton iterations
      subroutine traction_continuityNR(rntwc, rugrd, hbratio, rFu0,
     &                                 rFz0, rFz0h, del, nsvars, hvarsU,
     &                                 hvarsZ, nmprops, props, props_tw,
     &                                 Fe, nnode, elmt_number, l, mcrd,
     &                                 time, dtimeCur, kstep, kinc,
     &                                 F_tw_inv, rdLtw, rdLtwfac, rmh,
     &                                 rFu1, rFz1, rFz1h, lNRerr)
        implicit none
        ! Global variables
        integer nsvars, nmprops, nnode, elmt_number, l, mcrd, kstep
        integer kinc
        logical lNRerr
        real*8 rmh(3), rntwc(3), rugrd(3, 3), hbratio
        real*8 rFu0(3, 3), rFz0(3, 3), rFz0h(3, 3), del(3, 3)
        real*8 hvarsU(nsvars), hvarsZ(nsvars), props(nmprops)
        real*8 props_tw(nmprops), Fe(3, 3), time(2), dtimeCur
        real*8 F_tw_inv(3, 3), rFu1(3, 3), rFz1(3, 3), rFz1h(3, 3)
        real*8 rdLtw(3, 3), rdLtwfac

        ! Local variables
        integer iter, maxiter, i, j, iterErr, aidx, k, m, n, p
        real*8 mtpn(3, 3), rdLu(3, 3), rdLz(3, 3), rstress_U(6)
        real*8 rGsCoords(3), rNRFac, rNRTol, rmh0(3), rdmh(3)
        real*8 rstress_Z(6), FeZ(3, 3), sigmaU(3, 3), sigmaZ(3, 3)
        real*8 g_tr(3), g_tr_err, z_tr(3, 3), z_tr_inv(3, 3)
        real*8 rWu(3, 3, 3, 3), rWz(3, 3, 3, 3)
        real*8 rugrdU(3, 3), rugrdZ(3, 3)
        real*8 rWut, rWzt, rdLtwU(3, 3)

        ! Set NR parameters
        maxiter = 100
        rNRFac = 0.25d0
        rNRTol = 1.d-6
        g_tr_err = 1.d4
        lNRerr = .true.
        rdLtwU = 0.d0

        ! Newton's iteration loop
        do iter = 1, maxiter
          ! DEBUGGING
          ! print*, "iter = ", iter

          ! Compute m tensor product with n
          do j = 1, 3
            do i = 1, 3
              mtpn(i, j) = rmh(i) * rntwc(j)
            end do
          end do

          ! Compute jump velocity gradient increments
          rugrdU = rugrd - mtpn
          rugrdZ = rugrd + (hbratio - 1.d0) * mtpn

          ! Compute new deformation gradient in the untwinned region
          call inverse_3x3(del - rugrdU, rFu1)

          ! Compute Stress (untwinned region)
          call umat_stress(rStress_U, hvarsU, rGsCoords, rFu0, rFu1,
     &                     props, nmprops, nsvars, nnode,
     &                     elmt_number, l, mcrd, time, dtimeCur, kstep,
     &                     kinc, .false., rWu, iterErr, Fe, rdLtwU,
     &                     0.d0)

          ! Large change suggested by NR, reduce step by correction factor
          if (iterErr.eq.1) then
            ! DEBUGGING
            print*, "Convergence umat failed in the matrix"
            print*, "rmh = ", rmh
            if (iter.eq.1) then
              return
            end if
            rdmh = rNRFac * rdmh
            rmh = rmh0 + rdmh
            cycle
          end if

          ! Compute new deformation gradient in the twinned region
          call inverse_3x3(del - rugrdZ, rFz1)

          ! Compute elastoplastic part of twinned deformation gradient
          rFz1h = matmul(rFz1, F_tw_inv)

          ! Compute Stress (twinned region)
          call umat_stress(rStress_Z, hvarsZ, rGsCoords, rFz0h, rFz1h,
     &                     props_tw, nmprops, nsvars, nnode, 
     &                     elmt_number, l, mcrd, time, dtimeCur, kstep,
     &                     kinc, .false., rWz, iterErr, FeZ, rdLtw,
     &                     rdLtwfac)

          ! Large change suggested by NR, reduce step by correction factor
          if (iterErr.eq.1) then
            ! DEBUGGING
            print*, "Convergence umat failed in the band"
            print*, "rmh = ", rmh
            if (iter.eq.1) then
              return
            end if
            rdmh = rNRFac * rdmh
            rmh = rmh0 + rdmh
            cycle
          end if

          ! vector to matrix for stress in twinned and untwinned region
          call vec_to_sym_mat(rStress_U, sigmaU)
          call vec_to_sym_mat(rStress_Z, sigmaZ)

          ! Form residual
          do aidx = 1, 3
            g_tr(aidx) = 0.d0
            do i = 1, 3
              g_tr(aidx) = g_tr(aidx) - rntwc(i) *
     &                     (sigmaZ(i,aidx) - sigmaU(i,aidx))
            end do
          end do

          ! If residual is increasing, reduce step by correction factor
          if (norm2(g_tr).gt.g_tr_err) then
            rdmh = rNRFac * rdmh
            rmh = rmh0 + rdmh
            cycle
          end if

          ! Check for residual convergence
          g_tr_err = norm2(g_tr)      ! (absolute error)
          ! DEBUGGING
          ! print*, "g_tr_err = ", g_tr_err
          if (g_tr_err.lt.rNRTol) then
            ! print*, "converged rmh = ", rmh
            lNRerr = .false.
            return
          end if

    !       ! Form residual tangent matrix
    !       do j = 1, 3
    !         do aidx = 1, 3
    !           z_tr(aidx, j) = 0.d0
    !           do i = 1, 3
    !             do k = 1, 3
    !               do p = 1, 3
    !                 z_tr(aidx, j) = z_tr(aidx, j) +
    !  &                      rntwc(i) * ((hbratio - 1.d0)*
    !  &                      rWz(i,aidx,j,k)*rFz0h(p,k)+rWu(i,aidx,j,k)*
    !  &                      rFu1sq(p,k)) * rntwc(p)
    !               end do
    !             end do
    !           end do
    !         end do
    !       end do


          ! Form residual tangent matrix
          do j = 1, 3
            do aidx = 1, 3
              z_tr(aidx, j) = 0.d0
              do i = 1, 3
                do k = 1, 3
                  rWzt = 0.d0
                  rWut = 0.d0
                  do n = 1, 3
                    do m = 1, 3
                      rWut = rWut + rWu(i,aidx,m,n)*rFu1(m,j)*rFu1(k,n)
                      rWzt = rWzt + rWz(i,aidx,m,n)*rFz1(m,j)*rFz1h(k,n)
                    end do
                  end do

                  z_tr(aidx, j) = z_tr(aidx, j) +
     &                    rntwc(i)*((hbratio-1.d0)*rWzt + rWut)*rntwc(k)
                end do
              end do
            end do
          end do

          ! Compute change in m hat
          call inverse_3x3(z_tr, z_tr_inv)
          rdmh = matmul(z_tr_inv, g_tr)

          ! Compute updated m hat
          rmh0 = rmh
          rmh = rmh0 + rdmh
        end do

      end subroutine