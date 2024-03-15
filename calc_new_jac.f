CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       Implements Material Jacobian d sigma / d F required for        C
C       traction continuity in vuel. Note this is different from       C
C       traditional abaqus material jacobian (d sigma / d Ut).         C
C       The algorithm is based on F = Fe Fp type split generally       C
C       assumed in crystal plasticity.                                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine calc_new_jac(dfgrd1, FpInv_n, props, nprops, C, nslip,
     &                     Smat, dGamdTau, Cmat, dgamma, Fe, spk2, rW)
        implicit none
        integer :: nprops, nslip, i, j, k, l, m, n, p, q, islp
        real*8 :: dfgrd1(3, 3), FpInv_n(3, 3), props(nprops)
        real*8 :: Smat(3, 3, nslip), dGamdTau(nslip)
        real*8 :: dgamma(nslip), Fe(3, 3), spk2(3, 3)
        real*8 :: rL(3, 3, 3, 3), temp, temp2, C(3, 3, 3, 3)
        real*8 :: rD(3, 3, 3, 3), rJ(3, 3, 3, 3, nslip), rB(3, 3, nslip)
        real*8 :: del(3, 3), rK(3, 3, 3, 3), Cmat(3, 3, nslip)
        real*8 :: rQ(3, 3, 3, 3), rR(3, 3, nslip), rT(3, 3, 3, 3)
        real*8 :: rW(3, 3, 3, 3), detF, determinant, FeInv(3, 3)
        real*8 :: rK2(9, 9), rK2Inv(9, 9)

        integer :: map(9, 2), indx(6)

        del = reshape( (/1.,0.,0.,0.,1.,0.,0.,0.,1. /), (/3,3/) )

        ! Compute L_ijkl
        do l = 1, 3
          do k = 1, 3
            do j = 1, 3
              do i = 1, 3
                rL(i,j,k,l) = 0.0d0
                do m = 1, 3
                  rL(i,j,k,l) = rL(i,j,k,l) +
     &              FpInv_n(l,i) * dfgrd1(k,m) * FpInv_n(m,j) +
     &              FpInv_n(m,i) * dfgrd1(k,m) * FpInv_n(l,j)
                end do
              end do
            end do
          end do
        end do

        ! Compute D_ijkl
        do l = 1, 3
          do k = 1, 3
            do j = 1, 3
              do i = 1, 3
                rD(i,j,k,l) = 0.0d0
                do n = 1, 3
                  do m = 1, 3
                    rD(i,j,k,l) = rD(i,j,k,l) +
     &                0.5 * C(i,j,m,n) * rL(m,n,k,l)
                  end do
                end do
              end do
            end do
          end do
        end do

        ! Compute J^alpha_ijkl
        do islp = 1, nslip
          do l = 1, 3
            do k = 1, 3
              do j = 1, 3
                do i = 1, 3
                  rJ(i,j,k,l,islp) = 0.0d0
                  do n = 1, 3
                    do m = 1, 3
                      temp = 0.0d0
                      do p = 1, 3
                        temp = temp + rL(m,p,k,l) * Smat(p,n,islp) +
     &                         Smat(p,m,islp) * rL(n,p,k,l)
                      end do

                      rJ(i,j,k,l,islp) = rJ(i,j,k,l,islp) +
     &                                   0.5d0 * C(i,j,m,n) * temp
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do

        ! Compute B^alpha_ij
        do islp = 1, nslip
          do j = 1, 3
            do i = 1, 3
              rB(i,j,islp) = dGamdTau(islp) * 0.5 * (Smat(i,j,islp) +
     &                   Smat(j, i, islp))
            end do
          end do
        end do


        ! Compute K_ijkl
        do l = 1, 3
          do k = 1, 3
            do j = 1, 3
              do i = 1, 3
                temp = 0.0d0
                do islp = 1, nslip
                  temp = temp + Cmat(i,j,islp) * rB(k,l,islp)
                end do

                rK(i,j,k,l) = del(i,k) * del(j,l) + temp
              end do
            end do
          end do
        end do

        ! Compute Q_ijkl
        call compute_Q(rK, rD, dgamma, nslip, rJ, rQ)

        ! Compute R^alpha_ij
        do islp = 1, nslip
          do j = 1, 3
            do i = 1, 3
              rR(i,j,islp) = 0.0d0
              do l = 1, 3
                do k = 1, 3
                  rR(i,j,islp) = rR(i,j,islp) +
     &                           rB(k,l,islp) * rQ(k,l,i,j)
                end do
              end do
            end do
          end do
        end do

        ! Compute T_ijkl
        do l = 1, 3
          do k = 1, 3
            do j = 1, 3
              do i = 1, 3
                rT(i,j,k,l) = 0.0d0
                do n = 1, 3
                  temp = 0.0d0
                  do islp = 1, nslip
                    temp = temp + dgamma(islp) * Smat(n,j,islp)
                  end do

                  rT(i,j,k,l) = rT(i,j,k,l) + del(i,k) * FpInv_n(l,n) *
     &                          (del(n,j) - temp)
                  temp = 0.0d0
                  do islp = 1, nslip
                    temp = temp + rR(k,l,islp) * Smat(n,j,islp)
                  end do

                  do m = 1, 3
                    rT(i,j,k,l) = rT(i,j,k,l) - dfgrd1(i,m) *
     &                            FpInv_n(m,n) * temp
                  end do
                end do
              end do
            end do
          end do
        end do

        ! Compute W_ijkl
        detF = determinant(Fe)
        call inverse_3x3(Fe, FeInv)
        do l = 1, 3
          do k = 1, 3
            temp2 = 0.0d0
            do q = 1, 3
              do p = 1, 3
                temp2 = temp2 + rT(p,q,k,l) * FeInv(q,p)
              end do
            end do

            do j = 1, 3
              do i = 1, 3
                temp = 0.0d0
                do n = 1, 3
                  do m = 1, 3
                    temp = temp + rT(i,m,k,l) * spk2(m,n) * Fe(j,n) +
     &                     Fe(i,m) * rQ(m,n,k,l) * Fe(j,n) +
     &                     Fe(i,m) * spk2(m,n) * rT(j,n,k,l) -
     &                     Fe(i,m) * spk2(m,n) * Fe(j,n) * temp2

                  end do
                end do

                rW(i,j,k,l) = 1.0 / detF * temp
              end do
            end do
          end do
        end do

        ! TODO Construct 6x9 W_ij = dT/dF

      end subroutine

      subroutine compute_Q(rK, rD, dgamma, nslip, rJ, rQ)
        implicit none
        integer :: nslip, map(9,2), i, j, k, l, indx(9)
        real*8 :: rK(3,3,3,3), rD(3,3,3,3), dgamma(nslip)
        real*8 :: rJ(3,3,3,3,nslip), rQ(3,3,3,3), rhs(3,3,3,3)
        real*8 :: rhs2(9,9), rK2(9,9), rk2Inv(9,9), rdval

        ! Compute rhs = D - Summation dgamma(slp) * J(slp)
        rhs = rD
        do i = 1, nslip
          rhs = rhs - dgamma(i) * rJ(:,:,:,:,i)
        end do

        ! Stacking [11, 22, 33, 23, 13, 12, 32, 31, 21]^T
        map(1, 1) = 1
        map(1, 2) = 1
        map(2, 1) = 2
        map(2, 2) = 2
        map(3, 1) = 3
        map(3, 2) = 3
        map(4, 1) = 2
        map(4, 2) = 3
        map(5, 1) = 1
        map(5, 2) = 3
        map(6, 1) = 1
        map(6, 2) = 2
        map(7, 1) = 3
        map(7, 2) = 2
        map(8, 1) = 3
        map(8, 2) = 1
        map(9, 1) = 2
        map(9, 2) = 1
        do j = 1, 9
          do i = 1, 9
            rK2(i, j) = rK(map(i,1), map(i,2), map(j,1), map(j,2))
            rhs2(i,j) = rhs(map(i,1), map(i,2), map(j,1), map(j,2))
          end do
        end do

        ! Compute Q = K^-1 rhs
        call ludcmp(rK2, 9, 9, indx, rdval)
        do j = 1, 9
          call lubksb(rK2, 9, 9, indx, rhs2(:, j))
        end do

        ! Restack rQ2 to fourth order rQ
        do j = 1, 9
          do i = 1, 9
            rQ(map(i,1), map(i,2), map(j,1), map(j,2)) = rhs2(i,j)
          end do
        end do
      end subroutine