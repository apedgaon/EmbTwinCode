      subroutine umat_crys (stress, statev, ddsdde, sse, spd, scd, rpl,
     &                 ddsddt, drplde, drpldt, stran, dstran, time,
     &                 dtime, temp, dtemp, predef, dpred, cmname, ndi,
     &                 nshr, ntens, nstatv, props, nprops, coords, drot,
     &                 pnewdt, celent, dfgrd0, dfgrd1, noel, npt, layer,
     &                 kspt, jstep, kinc, rW, iterErr, Fe)

        implicit none
!------------------------------------------------------------------
!       UMAT argument variables definition
!------------------------------------------------------------------
        integer ndi                 ! Number of direct stress components
        integer nshr                ! Number of engg shear stress components
        integer ntens               ! Size of the stress or strain component array (ndi+nshr)
        integer nstatv              ! Number of solution-dependent state variables
        integer nprops              ! Number of user defined material properties
        integer noel                ! Element number
        integer npt                 ! Integration (Gauss) point number
        integer layer               ! Layer number (for composite shells and layered solids)
        integer kspt                ! Section point number within the current layer
        integer jstep(4)            ! Step number, procedure type, nlgeom, linear perturbation
        integer kinc                ! Increment number
        character*80 cmname         ! User-defined material name
        real*8 stress(ntens)        ! Stress tensor array (Cauchy)
        real*8 statev(nstatv)       ! Solution-dependent state variables array
        real*8 ddsdde(ntens, ntens) ! Jacobian matrix (variation of stress increments w.r.t strain increments)
        real*8 sse                  ! Specific elastic strain energy
        real*8 spd                  ! Specific plastic dissipation
        real*8 scd                  ! Specific creep dissipation
        real*8 rpl                  ! Volumetric heat generation per unit time
        real*8 ddsddt(ntens)        ! Variation of the stress increments w.r.t the temperature
        real*8 drplde(ntens)        ! Variation of rpl w.r.t the strain increments
        real*8 drpldt               ! Variation of rpl w.r.t the temperature
        real*8 stran(ntens)         ! Total strain array
        real*8 dstran(ntens)        ! Strain increment array
        real*8 time(2)              ! Step time, total time
        real*8 dtime                ! Time increment
        real*8 temp                 ! Temperature
        real*8 dtemp                ! Temperature increment
        real*8 predef(1)            ! Interpolated values of predefined field variables
        real*8 dpred(1)             ! Increments of predefined field variables
        real*8 props(nprops)        ! User-specified array of material constants
        real*8 coords(3)            ! Coordinates of this material (gauss) point
        real*8 drot(3)              ! Rotation increment matrix
        real*8 pnewdt               ! Ratio of suggested new time increment
        real*8 celent               ! Characteristic element length
        real*8 dfgrd0(3, 3)         ! Deformation gradient at the beginning of the increment
        real*8 dfgrd1(3, 3)         ! Deformation gradient at the end of the increment
        real*8 rW(3, 3, 3, 3)

        ! For a mechanical analysis, atleast stress and ddsdde needs to be defined
        ! For a coupled thermo-mechanical analysis, rpl, ddsddt, drplde, drpldt in addition maybe required

!------------------------------------------------------------------
!       User variables definition
!------------------------------------------------------------------
        integer i, j, k, l, m, n, ia
        integer nslip, iterErr
        real*8 pi
        real*8 del(3, 3)
        real*8 del4th(3, 3, 3, 3)
        real*8 xmC(3, int(props(34)))
        real*8 xsC(3, int(props(34)))
        real*8 xm0(3, int(props(34)))
        real*8 xs0(3, int(props(34)))
        real*8 theta
        real*8 c11
        real*8 c12
        real*8 c13
        real*8 c33
        real*8 c44
        real*8 c55
        real*8 C_c(3, 3, 3, 3)
        real*8 C(3, 3, 3, 3)
        real*8 a11
        real*8 a33
        real*8 AThExp_c(3, 3)
        real*8 AThExp(3, 3)
        real*8 th
        real*8 phi
        real*8 om
        real*8 Qmat(3, 3)
        real*8 Fpinv_n(3, 3)
        real*8 Fpinv(3, 3)
        real*8 Fe_tr(3, 3)
        real*8 Fe(3, 3)
        real*8 Ce_tr(3, 3)
        real*8 Ee_tr(3, 3)
        real*8 Spk2_tr(3, 3)
        real*8 Spk2v_tr(6)
        real*8 Spk2_n(3, 3)
        real*8 Spk2v_n(6)
        real*8 Spk2(3, 3)
        real*8 Spk2v(6)
        real*8 Smat(3, 3, int(props(34)))
        real*8 Pmat(3, 3, int(props(34)))
        real*8 Bmat(3, 3, int(props(34)))
        real*8 Cmat(3, 3, int(props(34)))
        real*8 Pmatr(6, int(props(34)))
        real*8 Cmatr(6, int(props(34)))
        real*8 s_n(int(props(34)))
        real*8 s(int(props(34)))
        real*8 dgamma(int(props(34)))
        real*8 rdLp(3, 3)
        real*8 detFp
        real*8 detFe
        real*8 sigma(3, 3)
        real*8 dgamdtau(int(props(34)))
        real*8 rhossd_n(int(props(34)))
        real*8 rhossd(int(props(34)))
        real*8 burg(int(props(34)))
        real*8 chi_n
        real*8 chi
        real*8 amat(int(props(34)), int(props(34)))
        real*8 alphaT
        real*8 mu
        real*8 tau(int(props(34)))
        real*8 dWp
        real*8 betaT
        real*8 chi_ss
        real*8 rho
        real*8 cp
        real*8 kappax
        real*8 e_D
        real*8 theta0

        real*8, external :: determinant

!------------------------------------------------------------------
!       User code
!------------------------------------------------------------------
        theta0 = props(1)

        ! Define pi and identiity tensors
        pi = 4.d0 * atan(1.d0)
        del = reshape( (/1.,0.,0.,0.,1.,0.,0.,0.,1. /), (/3,3/) )
        do i = 1,3
          do j = 1,3
            do k = 1,3
              do l = 1,3
                del4th(i,j,k,l) = 0.5d0 * 
     &            (del(i,k) * del(j,l) + del(i,l) * del(j,k))
            end do
           end do
          end do
         end do

        ! Define slip systems for hcp Ti
        ! (1-12) -> 1st order Pyramidal <c + a> slip {1 0 -1 1} <-2 1 1 3>
        ! (12-15) -> Prismatic <a> slip {1 0 -1 0} <-1 2 -1 0>
        ! (15-18) -> Basal <a> slip {0 0 0 1} <-1 2 -1 0>
        ! xmC - Normals,        xsC - Directions
        nslip = int(props(34))
        xmC(:,1) = (/0.7602357125290463d0, 0.4389222932762118d0,
     &               0.47894559384193547d0/)
        xmC(:,2) = (/0.7602357125290463d0, 0.4389222932762118d0,
     &               0.47894559384193547d0/)
        xmC(:,3) = (/0.0000000000000000d0, 0.8778445865524237d0,
     &               0.4789455938419355d0/)
        xmC(:,4) = (/0.0000000000000000d0, 0.8778445865524237d0,
     &               0.4789455938419355d0/)
        xmC(:,5) = (/-0.7602357125290463d0, 0.4389222932762118d0,
     &               0.47894559384193547d0/)
        xmC(:,6) = (/-0.7602357125290463d0, 0.4389222932762118d0,
     &               0.47894559384193547d0/)
        xmC(:,7) = (/-0.7602357125290463d0, -0.4389222932762118d0,
     &               0.47894559384193547d0/)
        xmC(:,8) = (/-0.7602357125290463d0, -0.4389222932762118d0,
     &               0.47894559384193547d0/)
        xmC(:,9) = (/0.0000000000000000d0, -0.8778445865524237d0,
     &               0.4789455938419355d0/)
        xmC(:,10) = (/0.000000000000000d0, -0.8778445865524237d0,
     &                0.4789455938419355d0/)
        xmC(:,11) = (/0.7602357125290463d0, -0.4389222932762118d0,
     &                0.47894559384193547d0/)
        xmC(:,12) = (/0.7602357125290463d0, -0.4389222932762118d0,
     &                0.47894559384193547d0/)
        xmC(:,13) = (/0.8660254037844386d0, 0.50000000000000000d0,
     &                0.000000000000000d0/)
        xmC(:,14) = (/0.0000000000000000d0, 1.00000000000000000d0,
     &                0.0000000000000000d0/)
        xmC(:,15) = (/-0.8660254037844386d0, 0.5000000000000000d0,
     &                0.0000000000000000d0/)
        xmC(:,16) = (/0.0000000000000000d0, 0.0000000000000000d0,
     &                1.0000000000000000d0/)
        xmC(:,17) = (/0.0000000000000000d0, 0.0000000000000000d0,
     &                1.0000000000000000d0/)
        xmC(:,18) = (/0.0000000000000000d0, 0.0000000000000000d0,
     &                1.0000000000000000d0/)


        xsC(:,1)  = (/-0.26651767092823825d0, -0.4616221471626314d0,
     &                0.8460929761745747d0/)
        xsC(:,2)  = (/-0.5330353418564765d0, 0.0000000000000000d0,
     &                0.8460929761745747d0/)
        xsC(:,3)  = (/0.26651767092823825d0, -0.4616221471626314d0,
     &                0.8460929761745747d0/)
        xsC(:,4)  = (/-0.26651767092823825d0, -0.4616221471626314d0,
     &                0.8460929761745747d0/)
        xsC(:,5)  = (/0.5330353418564765d0, 0.0000000000000000d0,
     &                0.8460929761745747d0/)
        xsC(:,6)  = (/0.26651767092823825d0, -0.4616221471626314d0,
     &                0.8460929761745747d0/)
        xsC(:,7)  = (/0.26651767092823825d0, 0.4616221471626314d0,
     &                0.8460929761745747d0/)
        xsC(:,8)  = (/0.5330353418564765d0, 0.0000000000000000d0,
     &                0.8460929761745747d0/)
        xsC(:,9)  = (/-0.26651767092823825d0, 0.4616221471626314d0,
     &                0.8460929761745747d0/)
        xsC(:,10)  = (/0.26651767092823825d0, 0.4616221471626314d0,
     &                 0.8460929761745747d0/)
        xsC(:,11)  = (/-0.5330353418564765d0, 0.0000000000000000d0,
     &                 0.8460929761745747d0/)
        xsC(:,12)  = (/-0.26651767092823825d0, 0.4616221471626314d0,
     &                 0.8460929761745747d0/)
        xsC(:,13)  = (/-0.500000000000000d0, 0.8660254037844386d0,
     &                 0.000000000000000d0/)
        xsC(:,14)  = (/-1.00000000000000d0, 0.0000000000000000d0,
     &                 0.000000000000000d0/)
        xsC(:,15)  = (/-0.50000000000000d0, -0.8660254037844386d0,
     &                 0.000000000000000d0/)
        xsC(:,16)  = (/-0.500000000000000d0, 0.8660254037844386d0,
     &                 0.000000000000000d0/)
        xsC(:,17)  = (/-1.00000000000000d0, 0.0000000000000000d0,
     &                 0.000000000000000d0/)
        xsC(:,18)  = (/-0.50000000000000d0, -0.8660254037844386d0,
     &                 0.000000000000000d0/)

        ! Initiliaze variables based on props or state variables
        if (time(2) .eq. 0.0) then ! t = 0 
          theta = theta0
          Fpinv_n = del
          Spk2_n = 0.d0
          rhossd_n = props(14)
          chi_n = props(15)
        else  ! time /= 0
          theta = statev(1)

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
          Spk2_n(2, 1) = Spk2_n(1, 2)
          Spk2_n(3, 1) = Spk2_n(1, 3)
          Spk2_n(3, 2) = Spk2_n(2, 3)

          n = 19
          do i = 1, nslip
            n = n + 1
            s_n(i) = statev(n)
          end do

          do i = 1, nslip
            n = n + 1
            rhossd_n(i) = statev(n)
          end do

          chi_n = statev(56)
        end if

        ! Retrieve thermal expansion coeff and euler angles
        a11 = props(4)
        th = props(31) * pi / 180.d0
        phi = props(32) * pi / 180.d0
        om = props(33) * pi / 180.d0
        a33 = props(51)

        !  Initialize 4th rank elastic stiffness tensor in crystal basis
        call compute_elastconst_hcp(props,nprops,theta,c11,c12,c44,
     &                              c33, c55, c13)
        C_c(:,:,:,:) = 0.d0
        do i = 1,3
          do j = 1,3
            do k = 1,3
              do l = 1,3
                C_c(i,j,k,l) = c12 * del(i,j) * del(k,l) +
     &                 c44 * (del(i,k) * del(j,l) + del(i,l) * del(k,j))
              end do
            end do
          end do
        end do
        C_c(1,1,1,1) = c11
        C_c(2,2,2,2) = c11
        C_c(3,3,3,3) = c33
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

        ! Form thermal expansion tensor in crystal basis
        AThExp_c =
     &   reshape((/a11,0.d0,0.d0,0.d0,a11,0.d0,0.d0,0.d0,a33 /),(/3,3/))

        ! Rotate slip, stiffness and thermal expansion to global basis
        call rotmat(th, phi, om, Qmat)
        call rotate_4th(Qmat, C_c, C)
        xm0 = matmul(Qmat, xmC)
        xs0 = matmul(Qmat, xsC)
        AThExp = matmul(Qmat, matmul(AThExp_c, transpose(Qmat)))

        ! Compute trial quantities
        Fe_tr = matmul(dfgrd1, Fpinv_n)
        Ce_tr = matmul(transpose(Fe_tr), Fe_tr)
        Ee_tr = 0.5d0 * (Ce_tr - del)  - AThExp * (theta - theta0)
        Spk2_tr = reshape( (/((sum(C(i,j,:,:)*Ee_tr), 
     &          i = 1,3), j=1,3)/), (/3,3/))

        ! Form Schmid tensor S^alpha
        do i = 1,3
          do j = 1,3
            Smat(i,j,:) = xs0(i,:)*xm0(j,:)
          enddo
        enddo

        ! Form symmetric part of Schmid tensor sym(S^alpha)
        do i = 1,3
          do j = 1,3
            Pmat(i,j,:) = 0.5d0*(Smat(i,j,:)+Smat(j,i,:))
          enddo
        enddo

        ! Form sym(C^e trial S^alpha) required in stress update
        Bmat(:,:,:) = 0.d0
        do ia = 1,nslip
          Bmat(:,:,ia) = 0.5d0 * (matmul(Ce_tr,Smat(:,:,ia)) +
     &                            matmul(transpose(Smat(:,:,ia)),Ce_tr))
        enddo

        ! From C : sym(C^e trial S^alpha) required in stress update
        Cmat(:,:,:) = 0.d0
        do i = 1,3
          do j = 1,3
            do m = 1,3
              do n = 1,3
                Cmat(i,j,:) = Cmat(i,j,:) + C(i,j,m,n)*Bmat(m,n,:)                 
              enddo
            enddo
          enddo
        enddo

        ! Convert symmetric matrix to vectors for easier computations
        call smatvec(Spk2_n, Spk2v_n)
        call smatvec(Spk2_tr, Spk2v_tr)
        do ia = 1,nslip
          call smatvec(Cmat(:, :, ia), Cmatr(:, ia)) 
          call smatvec(Pmat(:, :, ia), Pmatr(:, ia))
        enddo
        Pmatr(4, :) = 2.d0 * Pmatr(4, :) 
        Pmatr(5, :) = 2.d0 * Pmatr(5, :) 
        Pmatr(6, :) = 2.d0 * Pmatr(6, :)

        ! Calculate dislocation interaction tensor
        call calc_amat2(nslip, props, nprops, xsC, xmC, amat)

        ! Initialize slip resistance for time 0
        mu = sqrt(c44 * c55)
        if (time(2).eq.0.0) then
          alphaT = props(12)
          call get_burg(props, nprops, nslip, burg)
          s_n = alphaT * mu * burg * sqrt(matmul(amat, rhossd_n))
        end if

        ! Implicit solve for stress, slip rate and hardening
        call solvcrys(nslip, props, nprops, theta, Spk2v_n, Spk2v_tr,
     &                Cmatr, Pmatr, s_n, rhossd_n, chi_n, amat, dtime,
     &                Spk2v, dgamma, dgamdtau, s, rhossd, tau, iterErr)

        ! Solver did not converge - tell Abaqus to cut back and repeat time step
        if (iterErr == 1) then
          pnewdt = 0.5d0
          return
        endif

        ! Compute increment in plastic velocity gradient from the converged slip increments
        rdLp = 0.d0
        do ia = 1, nslip
          rdLp = rdLp + dgamma(ia) * Smat(:, :, ia)
        end do

        ! Compute current plastic deformation gradient with volume correction
        Fpinv = matmul(Fpinv_n, del - rdLp)
        detFp = determinant(Fpinv)
        Fpinv = Fpinv / detFp**(1.d0/3.d0)

        ! Compute elastic deformation gradient and Cauchy's stress
        Fe = matmul(dfgrd1, Fpinv)
        detFe = determinant(Fe)
        call svec_to_smat(Spk2v, Spk2)
        sigma = matmul(matmul(Fe, Spk2), transpose(Fe)) / detFe

        ! Stress tensor to vector
        stress(1) = sigma(1, 1)
        stress(2) = sigma(2, 2)
        stress(3) = sigma(3, 3)
        stress(4) = sigma(1, 2)
        stress(5) = sigma(1, 3)
        stress(6) = sigma(2, 3)

        ! Compute plastic work
        dWp = dot_product(tau, dgamma)

        ! Update thermal and effective temperature
        rho = props(2)
        cp = props(3)
        kappax = props(5)
        chi_ss = props(16)
        e_D = props(23)
        betaT = chi_n / chi_ss  ! TQ coeff
        theta = theta + betaT * dWp / (rho * cp)
        chi = chi_n + kappax * e_D * (1.d0 - betaT) * dWp / mu

        ! Store state variables
        ! Store temperature
        n = 1
        statev(n) = theta

        ! Store euler angles

        ! Store FpInv
        n = 4
        do j = 1, 3
          do i = 1, 3
            n = n + 1
            statev(n) = Fpinv(i, j)
          end do
        end do

        ! Store 2nd PK stress
        statev(14) = Spk2(1, 1)
        statev(15) = Spk2(2, 2)
        statev(16) = Spk2(3, 3)
        statev(17) = Spk2(1, 2)
        statev(18) = Spk2(1, 3)
        statev(19) = Spk2(2, 3)

        ! Store hardening
        n = 19
        do i = 1, nslip
          n = n + 1
          statev(n) = s(i)
        end do

        ! Store dislocation densities
        do i = 1, nslip
          n = n + 1
          statev(n) = rhossd(i)
        end do

        ! Store effective temperature
        n = n + 1
        statev(n) = chi

        ! Calculate Stress Jacobian ddsdde
C        call calc_jac(nslip, props, nprops, C, dfgrd0, dfgrd1, Fpinv_n,
C     &                FpInv, Smat, Pmat, Cmat, dgamma, dgamdtau, del,
C     &                rdLp, sigma, ddsdde)
        call calc_new_jac(dfgrd1, FpInv_n, props, nprops, C, nslip,
     &                    Smat, dGamdTau, Cmat, dgamma, Fe, spk2, rW)
        return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          Extract independent elastic stiffness parameters from       C
C                     properties for hcp Ti                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC            
      subroutine compute_elastconst_hcp(props, nprops, theta, c11,
     &                                  c12, c44, c33, c55, c13)

        implicit none
        ! Argument Variables
        integer :: nprops
        real*8  :: theta
        real*8  :: props(nprops)
        real*8  :: c11, c12, c44, c33, c55, c13

        ! Local Variables
        real*8  :: rm11, rm12, rm44, rm33, rm55, rm13
        real*8  :: c11_0, c12_0, c44_0, c33_0, c55_0, c13_0

        rm11  = props(6)
        c11_0 = props(7)
        rm12  = props(8)
        c12_0 = props(9)
        rm44  = props(10)
        c44_0 = props(11)
        rm33  = props(45)
        c33_0 = props(46)
        rm55  = props(47)
        c55_0 = props(48)
        rm13  = props(49)
        c13_0 = props(50)

        c11 = c11_0 + rm11 * theta
        c12 = c12_0 + rm12 * theta
        !c44 = c44_0 + rm44 * theta
        c44 = 0.5d0 * (c11 - c12)
        c33 = c33_0 + rm33 * theta
        c55 = c55_0 + rm55 * theta
        c13 = c13_0 + rm13 * theta

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C            Computes rotation matrix from euler angles                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      subroutine rotmat(th,phi,om,q)
        implicit none
        real*8 th, phi, om, q(3, 3)
        real*8 sth, cth, sph, cph, som, com

        sth = dsin(th)
        cth = dcos(th)
        sph = dsin(phi)
        cph = dcos(phi)
        som = dsin(om)
        com = dcos(om)

        q(1,1) = cph*com-sph*som*cth
        q(1,2) = sph*com+cph*som*cth
        q(1,3) = som*sth

        q(2,1) = -cph*som-sph*com*cth
        q(2,2) = -sph*som+cph*com*cth
        q(2,3) =  com*sth

        q(3,1) = sph*sth
        q(3,2) = -cph*sth
        q(3,3) = cth

        return
      end

!===================================================================
!===================================================================
!
!  Rotates any 3x3x3x3 tensor by a rotation matrix.
!
!  c(i,j,k,l) = a(i,m) * a(j,n) * a(k,p) * a(l,q) * b(m,n,p,q)
!
!-------------------------------------------------------------------

      subroutine rotate_4th(a,b,c)
        implicit none
        integer i, j, k, l, m, n
        real*8 a(3,3), b(3,3,3,3), c(3,3,3,3), d(3,3,3,3)

        do m = 1,3
          do n = 1,3
            do k = 1,3
              do l = 1,3
                d(m,n,k,l) = a(k,1) * (a(l,1) * b(m,n,1,1) + 
     &                     a(l,2) * b(m,n,1,2) + a(l,3) * b(m,n,1,3)) +
     &                     a(k,2) * (a(l,1) * b(m,n,2,1) + 
     &                     a(l,2) * b(m,n,2,2) + a(l,3) * b(m,n,2,3)) +
     &                     a(k,3) * (a(l,1) * b(m,n,3,1) + 
     &                     a(l,2) * b(m,n,3,2) + a(l,3) * b(m,n,3,3))
              end do
            end do
          end do
        end do

        do i = 1,3
          do j = 1,3
            do k = 1,3
              do l = 1,3
                c(i,j,k,l) = a(i,1) * (a(j,1) * d(1,1,k,l) + 
     &                     a(j,2) * d(1,2,k,l) + a(j,3) * d(1,3,k,l)) +
     &                     a(i,2) * (a(j,1) * d(2,1,k,l) + 
     &                     a(j,2) * d(2,2,k,l) + a(j,3) * d(2,3,k,l)) +
     &                     a(i,3) * (a(j,1) * d(3,1,k,l) + 
     &                     a(j,2) * d(3,2,k,l) + a(j,3) * d(3,3,k,l))
              end do
            end do
          end do
        end do

        return
      end

C  THIS SUBROUTINE CONVERTS A SYMMETRIC MATRIX TO A VECTOR
      subroutine smatvec(a,v)
        implicit none
        real*8 a(3, 3), v(6)
        v(1) = a(1, 1)
        v(2) = a(2, 2)
        v(3) = a(3, 3)
        v(4) = a(1, 2)
        v(5) = a(1, 3)
        v(6) = a(2, 3)
        return
      end

!===================================================================
!  Convert from Voigt notation stress vector to a symmetric 
!-------------------------------------------------------------------
      subroutine svec_to_smat(stressv, stressm)

        implicit none
  
        real*8 :: stressv(6), stressm(3,3)
  
        stressm(1,1) = stressv(1)
        stressm(2,2) = stressv(2)
        stressm(3,3) = stressv(3)
        stressm(1,2) = stressv(4)
        stressm(2,1) = stressv(4)
        stressm(1,3) = stressv(5)
        stressm(3,1) = stressv(5)
        stressm(2,3) = stressv(6)
        stressm(3,2) = stressv(6)
      end subroutine

C     Implicit solve for crystal stress, slip rate and hardening (Balasubramanian 1998, Anand 2004)
      subroutine solvcrys(nslip, props, nprops, theta, Spk2v_n, 
     &                    Spk2v_tr, Cmatr, Pmatr, s_n, rhossd_n, chi_n,
     &                    amat, dtime, Spk2v, dgamma, dgamdtau, s,
     &                    rhossd, tau, iterErr)
        implicit none
        ! Argument variables definition
        integer nslip, nprops, iterErr
        real*8 props(nprops), theta, Spk2v_n(6), Spk2v_tr(6), Spk2v(6)
        real*8 Cmatr(6, nslip), Pmatr(6, nslip), dtime
        real*8 s_n(nslip), s(nslip), dgamma(nslip), dgamdtau(nslip)
        real*8 rhossd_n(nslip), rhossd(nslip), chi_n, amat(nslip, nslip)
        real*8 tau(nslip)

        ! Internal variables definition
        integer i, j, ia, indx(6), k, maxit1, l, maxit2, icorr
        real*8 gdotalpha, dgdtau, dgmax, corrbeta, dgmaxold
        real*8 beta(6), dBetadSp(6, 6), temp, dval, tol, tols
        real*8 ds(nslip), err_s, s_old(nslip), Spk2v_old(6), dgmaxtol
        real*8 drhossd(nslip), sl(nslip)
        logical kErr

        ! Initialize
        Spk2v = Spk2v_n
        Spk2v_old = Spk2v_n
        s = s_n
        rhossd = rhossd_n
        dgmax = 0.d0
        tol = 1.d-8
        tols = 1.d-8
        maxit1 = 100
        maxit2 = 20
        dgmaxtol = 0.5d0
        corrbeta = 4.d0
        iterErr = 0

        ! Retrieve lattice resistances
        do ia = 1, nslip
          if (ia.le.12) then
            sl(ia) = props(54)
          else if ((ia.ge.13).and.(ia.le.15)) then
            sl(ia) = props(57)
          else if((ia.ge.16).and.(ia.le.18)) then
            sl(ia) = props(58)
          end if
        end do

        ! Second level iterations for hardening update
        do l = 1, maxit2
          ! First level iterations for stress update
          kErr = .true.
          do k = 1, maxit1
            ! Compute slip increments
            dgmax = 0.d0
            do ia = 1, nslip
              tau(ia) = dot_product(Spk2v, Pmatr(:, ia))
              call calc_gammadot(props, nprops, theta, tau(ia),
     &                           s(ia), gdotalpha, dgdtau, sl(ia))
              dgamma(ia) = gdotalpha * dtime
              dgamdtau(ia) = dgdtau * dtime
              dgmax = max(dgmax, abs(dgamma(ia)))
            end do

            ! Check if slip increments are relevant
            if (dgmax.lt.1.d-20) then
              Spk2v = Spk2v_tr
              return
            ! if dgamma is too big, adjust newton's correction by a factor to minimize residual
            else if ((dgmax.gt.dgmaxtol).and.(k.gt.1)) then
              if (icorr.eq.0) then
                dgmaxold = 0.d0
              end if
              if (icorr.le.10) then
                if (dgmax.gt.dgmaxold) then
                  corrbeta = 1.d0 / corrbeta
                else if (dgmax.eq.dgmaxold) then
                  corrbeta = 0.25d0
                end if
                beta = corrbeta * beta
                Spk2v = Spk2v_old + beta
                dgmaxold = dgmax
                icorr = icorr + 1
                cycle
              end if
            end if
            icorr = 0

            ! Compute negative of stress residual
            beta = Spk2v_tr - Spk2v
            do ia = 1, nslip
              beta = beta - dgamma(ia) * Cmatr(:, ia)
            end do

            ! Compute Stress update Jacobian
            do j = 1, 6
              do i = 1, 6
                if (i.eq.j) then
                  dBetadSp(i, j) = 1.d0
                else
                  dBetadSp(i, j) = 0.d0
                end if

                temp = 0.d0
                do ia = 1, nslip
                  temp = temp + Cmatr(i, ia)*dgamdtau(ia)*Pmatr(j, ia)
                end do
                dBetadSp(i, j) = dBetadSp(i, j) + temp
              end do
            end do

            ! Solve Newton corrections by LU decomposition
            call ludcmp(dBetadSp, 6, 6, indx, dval)
            call lubksb(dBetadSp, 6, 6, indx, beta)

            ! If new corrections are below tolerance, exit stress loop
            if ((abs(beta(1)).lt.tol).and.
     &          (abs(beta(2)).lt.tol).and.
     &          (abs(beta(3)).lt.tol).and.
     &          (abs(beta(4)).lt.tol).and.
     &          (abs(beta(5)).lt.tol).and.
     &          (abs(beta(6)).lt.tol)) then
              kErr = .false.
              exit
            end if

            ! If not, update stress with corrections
            Spk2v_old = Spk2v
            Spk2v = Spk2v_old + beta
          end do

          ! If stress update loop failed, exit hardening loop to return error
          if (kErr) then
            exit
          end if

          ! Stress update converged, calculate hardening increment
          call calc_hardening(props, nprops, nslip, theta, dgamma,
     &                        dtime, s_n, rhossd_n, chi_n, tau, amat, 
     &                        ds, drhossd)

          ! Update hardening and check for error tolerance
          rhossd = rhossd_n + drhossd
          s_old = s
          s = s_n + ds
          err_s = maxval(abs(s - s_old))
          if (err_s.lt.tols) then
            return
          end if
        end do

        ! Two level iterative procedure failed
        iterErr = 1
      end subroutine

C     Computes slip rate on a given slip system
      subroutine calc_gammadot(props, nprops, theta, taualpha, resalpha,
     &                         gdotalpha, dgdtau, sl)
        implicit none
        ! Argument variables definition
        integer nprops
        real*8 props(nprops), theta, taualpha, resalpha
        real*8 gdotalpha, dgdtau, sl

        ! Internal variables definition
        real*8 bk, gdot0, F0, p, q
        real*8 c11, c12, c44, c33, c55, c13, c11_0, c12_0, c44_0, c55_0
        real*8 xmu, xmu0, xmu_r, val1, val2
        real*8, external :: ramp

        ! Retrieve necessary properties
        bk = props(27)
        gdot0 = props(52)
        F0 = props(53)
        p = props(55)
        q = props(56)

        ! Compute elastic constants
        call compute_elastconst_hcp(props,nprops,theta,c11,c12,c44, c33,
     &                              c55, c13)

        ! Compute elastic constants for shear at 0 K to construct mu / mu0
        c11_0 = props(7)
        c12_0 = props(9)
        c44_0 = props(11)
        c55_0 = props(48)
        xmu = sqrt(c44 * c55)
        xmu0 = sqrt(c44_0 * c55_0)
        xmu_r = xmu/xmu0

        ! PQ type flow rule (gammadot)
        gdotalpha = gdot0 * sign(1.d0, taualpha) *
     &    exp(-F0/(bk*theta)*ramp(1.d0 -
     &      ramp((abs(taualpha)-resalpha)/(xmu_r*sl))**p)**q)

        ! Derivative of gammadot w.r.t rss for pq type flow rule
        dgdtau = 0.d0
        val1 = (abs(taualpha)-resalpha)/(xmu_r*sl)
        if (val1.gt.0.d0) then
          val2 = 1.d0 - ramp(val1)**p
          if (val2.gt.0.d0) then
            dgdtau = p*q*F0 / (bk*theta*sl*xmu_r) * gdotalpha *
     &               sign(1.d0, taualpha) * val1**(p-1) * val2**(q-1) 
          end if
        end if
      end subroutine

      subroutine calc_hardening(props, nprops, nslip, theta, dgamma,
     &                          dtime, s, rhossd, chi, tau, amat, ds,
     &                          drhossd)
        implicit none
        ! Argument variables definition
        integer nprops, nslip
        real*8 props(nprops), theta, dgamma(nslip), dtime
        real*8 s(nslip), ds(nslip), rhossd(nslip), chi, tau(nslip)
        real*8 amat(nslip, nslip), drhossd(nslip)

        ! Internal variables definition
        integer ia, ja
        real*8 k1, krho, nubar, burg(nslip), a, abratio
        real*8 c11, c12, c13, c33, c44, c55, mu, e_D, rhossd_ss
        real*8 dsdrhossd(nslip, nslip)
        real*8 rhossd_i(nslip), alphaT

        ! Retrieve properties and compute shear modulus
        alphaT = props(12)
        e_D = props(23)
        k1 = props(29)
        abratio = props(37)
        call get_burg(props, nprops, nslip, burg)
        call compute_elastconst_hcp(props, nprops, theta, c11, c12, c44,
     &                              c33, c55, c13)
        mu = sqrt(c44 * c55)

        ! Compute increment in dislocation density
        do ia = 1, nslip
          nubar = tau(ia) / s(ia)
          if (nubar.eq.0.d0) then
            drhossd(ia) = 0.d0
            cycle
          end if
          krho = k1 / nubar**2
          a = abratio * burg(ia)
          rhossd_ss = 1.d0 / a**2 * exp(-e_D / chi)
          drhossd(ia) = krho / (mu * a**2) * tau(ia) * dgamma(ia) *
     &                  (1.d0 - rhossd(ia) / rhossd_ss)
        end do

        ! Compute ds^alpha / drho^alpha
        rhossd_i = matmul(amat, rhossd)
        do ja = 1, nslip
          do ia = 1, nslip
            dsdrhossd(ia, ja) = alphaT *  mu * burg(ia) * amat(ia, ja) /
     &                          (2 * sqrt(rhossd_i(ia)))
          end do
        end do

        ! Compute increment in slip resistance (taylor stress)
        ds = matmul(dsdrhossd, drhossd) 
      end subroutine


      ! function ramp(x)
      !   implicit none
      !   real*8 x, ramp

      !   if (x.lt.0.d0) then
      !     ramp = 0.d0
      !   else
      !     ramp = x
      !   end if
      ! end function

C**************************************************************
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
C
C  Given an NxN matrix [A], with physical dimension NP, this 
C  routine replaces it by the LU decomposition of a row-wise 
C  permutation of itself. [A] and N are input. [A] is output, 
C  arranged in LU form. INDX is an output vector which records
C  the row permutation effected by the partial pivoting; 
C  D is output as +1 or -1 depending on wheter the nuber of
C  row interchanges was even or odd, respectively. This routine
C  is used in combination with LUBKSB to solve linear equations 
C  or invert a matrix.
C
        IMPLICIT REAL*8 (A-H,O-Z)
        
        PARAMETER (NMAX=100,TINY=1.0E-20)
        DIMENSION A(NP,NP),INDX(N),VV(NMAX)
        D=1.
        DO I=1,N
        AAMAX=0.
        DO J=1,N
          IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
        ENDDO
!        IF (AAMAX.EQ.0.) write(*,*) 'Singular matrix.'
        VV(I)=1./AAMAX
        ENDDO
        DO J=1,N
        IF (J.GT.1) THEN
          DO I=1,J-1
            rSUM=A(I,J)
            IF (I.GT.1)THEN
              DO K=1,I-1
                rSUM=rSUM-A(I,K)*A(K,J)
              ENDDO
              A(I,J)=rSUM
            ENDIF
          ENDDO
        ENDIF
        AAMAX=0.
        DO I=J,N
          rSUM=A(I,J)
          IF (J.GT.1)THEN
            DO K=1,J-1
              rSUM=rSUM-A(I,K)*A(K,J)
            ENDDO
            A(I,J)=rSUM
          ENDIF
          DUM=VV(I)*DABS(rSUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
        ENDDO
        IF (J.NE.IMAX)THEN
          DO K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
          ENDDO
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO I=J+1,N
            A(I,J)=A(I,J)*DUM
          ENDDO
        ENDIF
        ENDDO
        IF(A(N,N).EQ.0.)A(N,N)=TINY
        RETURN
        END
C**********************************************************
       SUBROUTINE LUBKSB(A,N,NP,INDX,B)
C
C  Solves the set of N linear equations [A]{X} = {B}. 
C  Here [A] is input, not as the matrix [A], but as its LU 
C  decomposition, determined by the routine LUDCMP. INDX
C  is input as the permutation vector returned by LUDCMP. {B}
C  is input as the right-hand side vector {B}, and returns
C  with the solution vector {X}. [A], N, NP, INDX are not 
C  modified by this routine, and can be left in place
C  for succesive calls with different right-hand sides {B}.
C  This routine takes into account that {B} will begin with
C  many zero elements, so it is efficient for use in matrix 
C  inversion.
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(NP,NP),INDX(N),B(N)
       II=0
       DO I=1,N
        LL=INDX(I)
        rSUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO J=II,I-1
            rSUM=rSUM-A(I,J)*B(J)
          ENDDO
        ELSE IF (rSUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=rSUM
       ENDDO
       DO I=N,1,-1
        rSUM=B(I)
        IF(I.LT.N)THEN
          DO J=I+1,N
            rSUM=rSUM-A(I,J)*B(J)
          ENDDO
        ENDIF
        B(I)=rSUM/A(I,I)
       ENDDO
       RETURN
       END

      function determinant(a)
        implicit none
        real*8 a(3, 3), b1, b2, b3, determinant

        b1 = a(2,2) * a(3,3) - a(3,2) * a(2,3)
        b2 = a(3,1) * a(2,3) - a(2,1) * a(3,3)
        b3 = a(2,1) * a(3,2) - a(3,1) * a(2,2)

        determinant = a(1,1) * b1 + a(1,2) * b2 + a(1,3) * b3
      end

      ! Stress Jacobian for CPFE (see Balasubramanian 1998)
      subroutine calc_jac(nslip, props, nprops, C, dfgrd0, dfgrd1,
     &                    Fpinv_n, Fpinv, Smat, Pmat, Cmat, dgamma,
     &                    dgamdtau, del, rdLp, sigma, rWr)
        implicit none
        ! Argument variables definition
        integer nslip, nprops
        real*8 props(nprops), C(3, 3, 3, 3), dfgrd0(3, 3), dfgrd1(3, 3)
        real*8 Fpinv_n(3, 3), Fpinv(3, 3), Smat(3, 3, nslip)
        real*8 Pmat(3, 3, nslip), dgamma(nslip), dgamdtau(nslip)
        real*8 del(3, 3), Cmat(3, 3, nslip), rdLp(3, 3), sigma(3, 3)
        real*8 rWr(6, 6)

        ! Internal variables definition
        integer i, j, k, l, m, n, p, q, ia, indx(6)
        real*8 rL(3, 3, 3, 3), dfgrd0inv(3, 3), rD(3, 3, 3, 3)
        real*8 Frel(3, 3), Rrel(3, 3), Urel(3, 3), Fe_n(3, 3), Fe(3, 3)
        real*8 sum1, sum2, rG(3, 3, 3, 3, nslip), rJ(3, 3, 3, 3, nslip)
        real*8 rB(3, 3, nslip), rK(3, 3, 3, 3), rKr(6, 6)
        real*8 rhs(3, 3, 3, 3), rhsr(6, 6), d, rQ(3, 3, 3, 3), sum3
        real*8 sum4, rR(3, 3, nslip), rS(3, 3, 3, 3), mat1(3, 3, 3, 3)
        real*8 mat2(3, 3), mat3(3, 3), detFe, rW(3, 3, 3, 3)
        real*8 Feinv(3, 3)
        
        real*8, external :: determinant

        ! Compute Relative deformation gradient, rotation and stretch
        call inverse_3x3(dfgrd0, dfgrd0inv)
        Frel = matmul(dfgrd1, dfgrd0inv)
        call polardecomp(Frel, Rrel, Urel)

        ! Compute Elastic deformation gradients at begining and end
        Fe_n = matmul(dfgrd0, Fpinv_n)
        Fe = matmul(dfgrd1, Fpinv)
        detFe = determinant(Fe)
        call inverse_3x3(Fe, Feinv)

        ! Compute L_ijkl
        do l = 1, 3
          do k = 1, 3
            do j = 1, 3
              do i = 1, 3
                sum1 = 0.d0
                sum2 = 0.d0
                sum3 = 0.d0
                sum4 = 0.d0
                do m = 1, 3
                  sum1 = sum1 + Urel(l, m) * Fe_n(m, j)
                  sum2 = sum2 + Fe_n(m, i) * Urel(m, k)
                  sum3 = sum3 + Urel(k, m) * Fe_n(m, j)
                  sum4 = sum4 + Fe_n(m, i) * Urel(m, l)
                end do

                rL(i, j, k, l) = 0.5d0* (Fe_n(k, i) * sum1 +
     &                           sum2 * Fe_n(l, j) + Fe_n(l, i) * sum3 +
     &                           sum4 * Fe_n(k, j))
              end do
            end do
          end do
        end do

        ! Compute D_ijkl
        do l = 1, 3
          do k = 1, 3
            do j = 1, 3
              do i = 1, 3
                sum1 = 0.d0
                do n = 1, 3
                  do m = 1, 3
                    sum1 = sum1 + C(i, j, m, n) * rL(m, n, k, l)
                  end do
                end do

                rD(i, j, k, l) = 0.5d0 * sum1
              end do
            end do
          end do
        end do

        ! Compute G_ijkl^alpha
        do ia = 1, nslip
          do l = 1, 3
            do k = 1, 3
              do j = 1, 3
                do i = 1, 3
                  sum1 = 0.d0
                  do p = 1, 3
                    sum1 = sum1 + rL(i, p, k, l)  * Smat(p, j, ia) +
     &                     Smat(p, i, ia) * rL(p, j, k, l)
                  end do

                  rG(i, j, k, l, ia) = sum1
                end do
              end do
            end do
          end do
        end do

        ! Compute J_ijkl^alpha
        do ia = 1, nslip
          do l = 1, 3
            do k = 1, 3
              do j = 1, 3
                do i = 1, 3
                  sum1 = 0.d0
                  do n = 1, 3
                    do m = 1, 3
                      sum1 = sum1 + C(i, j, m, n) * rG(m, n, k, l, ia)
                    end do
                  end do

                  rJ(i, j, k, l, ia) = 0.5d0 * sum1 
                end do
              end do
            end do
          end do
        end do

        ! Compute B^alpha
        do ia = 1, nslip
          rB(:, :, ia) = dgamdtau(ia) * Pmat(:, :, ia)
        end do

        ! Compute K_ijkl
        do l = 1, 3
          do k = 1, 3
            do j = 1, 3
              do i = 1, 3
                sum1 = 0.d0
                do ia = 1, nslip
                  sum1 = sum1 + Cmat(i, j, ia) * rB(k, l, ia)
                end do

                rK(i, j, k, l) = del(i, k) * del(j, l) + sum1
              end do
            end do
          end do
        end do

        ! Reduce K_ijkl
        call ematr(rK, rKr)

        ! Compute Rhs in the reduced form required for Q
        rhs = rD
        do ia = 1, nslip
          rhs = rhs - dgamma(ia) * rJ(:, :, :, :, ia)
        end do
        call ematr(rhs, rhsr)

        ! Compute Q by solving linear equations and expand
        call ludcmp(rKr, 6, 6, indx, d)
        do i = 1, 6
          call lubksb(rKr, 6, 6, indx, rhsr(:, i))
        end do
        call ematrec(rhsr, rQ)

        ! Compute R_ij^alpha
        do ia = 1, nslip
          do j = 1, 3
            do i = 1, 3
              sum1 = 0.d0
              do l = 1, 3
                do k = 1, 3
                  sum1 = sum1 + rB(k, l, ia) * rQ(k, l, i, j)
                end do
              end do

              rR(i, j, ia) = sum1
            end do
          end do
        end do

        ! Compute Summation of R_ij^alpha * S_0_kl^alpha over slip sytem
        do l = 1, 3
          do k = 1, 3
            do j = 1, 3
              do i = 1, 3
                sum1 = 0.d0
                do ia = 1, nslip
                  sum1 = sum1 + rR(k, l, ia) * Smat(i, j, ia)
                end do

                mat1(i, j, k, l) = sum1
              end do
            end do
          end do
        end do

        ! Sub matrices required for S_ijkl
        mat2 = matmul(Fe_n, rdLp)
        mat3 = matmul(Frel, Fe_n)

        ! Compute S_ijkl
        do l = 1, 3
          do k = 1, 3
            do j = 1, 3
              do i = 1, 3
                sum1 = 0.d0
                do p = 1, 3
                  sum1 = sum1 + mat3(i, p) * mat1(p, j, k, l)
                end do

                rS(i, j, k, l) = Rrel(i, k)*(Fe_n(l, j) - mat2(l, j)) -
     &                           sum1
              end do
            end do
          end do
        end do

        ! Compute W_ijkl
        do l = 1, 3
          do k = 1, 3
            sum1 = 0.d0
            do q = 1, 3
              do p = 1, 3
                sum1 = sum1 + rS(p, q, k, l) * Feinv(q, p)
              end do
            end do

            do j = 1, 3
              do i = 1, 3
                sum2 = 0.d0
                do n = 1, 3
                  do m = 1, 3
                    sum2 = sum2 +
     &                     rS(i, m, k, l) * sigma(m, n) * Fe(j, n) +
     &                     Fe(i, m) * rQ(m, n, k, l) * Fe(j, n) +
     &                     Fe(i, m) * sigma(m, n) * rS(j, n, k, l) -
     &                     Fe(i, m) * sigma(m, n) * Fe(j, n) * sum1
                  end do
                end do

                rW(i, j, k, l) = sum2 / detFe
              end do
            end do
          end do
        end do

        ! Reduce W_ijkl
        call ematr(rW, rWr)

      end subroutine

      subroutine inverse_3x3(a,b)
        implicit none

        ! Argument variables definition
        real*8 a(3, 3), b(3, 3)

        ! Internal variables definition
        integer i, j
        real*8 det

        b(1,1) = a(2,2) * a(3,3) - a(3,2) * a(2,3)
        b(1,2) = a(3,2) * a(1,3) - a(1,2) * a(3,3)
        b(1,3) = a(1,2) * a(2,3) - a(2,2) * a(1,3)
        b(2,1) = a(3,1) * a(2,3) - a(2,1) * a(3,3)
        b(2,2) = a(1,1) * a(3,3) - a(3,1) * a(1,3)
        b(2,3) = a(2,1) * a(1,3) - a(1,1) * a(2,3)
        b(3,1) = a(2,1) * a(3,2) - a(3,1) * a(2,2)
        b(3,2) = a(3,1) * a(1,2) - a(1,1) * a(3,2)
        b(3,3) = a(1,1) * a(2,2) - a(2,1) * a(1,2)

        det = a(1,1) * b(1,1) + a(1,2) * b(2,1) + a(1,3) * b(3,1)

        do i = 1,3
          do j = 1,3
              b(i,j) = b(i,j) / det
          end do
        end do
      end

! ********************************************************************
!
!  THIS SUBROUTINE PERFORMS THE RIGHT POLAR DECOMPOSITION
!  [F]=[R][U] OF THE DEFORMATION GRADIENT [F] INTO
!  A ROTATION [R] AND THE RIGHT  STRETCH TENSOR [U].
!  THE EIGENVALUES AND EIGENVECTORS OF [U] AND
!  THE  LOGARITHMIC STRAIN [E] = LN [U]
!  ARE ALSO RETURNED.
!       
!
! *******************************************************************
      SUBROUTINE polardecomp(F,R,U)

        IMPLICIT REAL*8(A-H,O-Z)
!
! ------------------------------------------------------------------
!  VARIABLES
!
        REAL*8 F(3,3),DETF,FTRANS(3,3),
     +   C(3,3), OMEGA(3),
     +   UEIGVAL(3),EIGVEC(3,3), EIGVECT(3,3), 
     +   U(3,3),E(3,3),UINV(3,3),R(3,3),TEMP(3,3)
         
!
!       F(3,3)  -- THE DEFORMATION GRADIENT MATRIX WHOSE
!                          POLAR DECOMPOSITION IS DESIRED.
!       DETF            -- THE DETRMINANT OF [F]; DETF > 0.
!       FTRANS(3,3)     -- THE TRANSPOSE OF [F].
!       R(3,3)  -- THE ROTATION MATRIX; [R]^T [R] = [I];
!                          OUTPUT.
!       U(3,3)  -- THE RIGHT STRETCH TENSOR; SYMMETRIC
!                          AND POSITIVE DEFINITE; OUTPUT.
!       UINV(3,3)       -- THE INVERSE OF [U].
!       C(3,3)  -- THE RIGHT CAUCHY-GREEN TENSOR = [U][U];
!                           SYMMETRIC AND POSITIVE DEFINITE.
!       OMEGA(3)        -- THE SQUARES OF THE PRINCIPAL STRETCHES.
!       UEIGVAL(3)      -- THE PRINCIPAL STRETCHES; OUTPUT.
!       EIGVEC(3,3)     -- MATRIX OF EIGENVECTORS OF [U];OUTPUT.
!       EIGVECT(3,3)-- TRANSPOSE OF THE ABOVE.
!       E(3,3)  -- THE LOGARITHMIC STRAIN TENSOR, [E]=LN[U];
!                          OUTPUT.
! -----------------------------------------------------------------
!       COMPUTATION
       
!  STORE THE IDENTITY MATRIX IN  [R], [U], AND [UINV]

        r = reshape( (/1.,0.,0.,0.,1.,0.,0.,0.,1. /), (/3,3/) )
        u = reshape( (/1.,0.,0.,0.,1.,0.,0.,0.,1. /), (/3,3/) )
        uinv = reshape( (/1.,0.,0.,0.,1.,0.,0.,0.,1. /), (/3,3/) )
        
        

!  STORE THE ZERO MATRIX IN [E]

        E(:,:) = 0.d0

!  CHECK IF THE DETERMINANT OF [F] IS GREATER THAN ZERO.
!  IF NOT, THEN PRINT DIAGONOSTIC AND STOP.

!       detF = determinant(F)
!       IF (DETF .LE. 0.D0) THEN
!          WRITE(91,100)
!          RETURN
!       END IF

!  CALCULATE THE RIGHT CAUCHY GREEN STRAIN TENSOR [C]

       ftrans = transpose(f)
       c = matmul(ftrans,f)
 
!  CALCULATE THE EIGENVALUES AND EIGENVECTORS OF  [C]

         CALL SPECTRAL(C,OMEGA,EIGVEC)

!  CALCULATE THE PRINCIPAL VALUES OF [U] AND [E]

         UEIGVAL(1) = DSQRT(OMEGA(1))
         UEIGVAL(2) = DSQRT(OMEGA(2))
         UEIGVAL(3) = DSQRT(OMEGA(3))

         U(1,1) = UEIGVAL(1)
         U(2,2) = UEIGVAL(2)
         U(3,3) = UEIGVAL(3)

         E(1,1) = DLOG( UEIGVAL(1) )
         E(2,2) = DLOG( UEIGVAL(2) )
         E(3,3) = DLOG( UEIGVAL(3) )

!  CALCULATE THE COMPLETE TENSORS [U] AND [E]
        eigvecT = transpose(eigvec)
        temp = matmul(eigvec,u)
        u = matmul(temp,eigvecT)
        temp = matmul(eigvec,e)
        e = matmul(temp,eigvecT)

!  CALCULATE [UINV]
        call inverse_3x3(U,Uinv)

!  CACULATE [R]
         r = matmul(f,uinv)
! -------------------------------------------------------------------
!       FORMATS

 100  FORMAT(5X,'--ERROR IN KINEMATICS-- THE DETERMINANT OF [F]',
     +       ' IS NOT GREATER THAN 0',/,
     +       10X,'POLAR DECOMPOSITION NOT PERFORMED')
! -------------------------------------------------------------------
        RETURN
        END
! ****************************************************************
!
!  THE FOLLOWING SUBROUTINES CALCULATE THE SPECTRAL
!  DECOMPOSITION OF A SYMMETRIC THREE BY THREE MATRIX
!
! **********************************************************************
        SUBROUTINE SPECTRAL(A,D,V)
!
!  THIS SUBROUTINE CALCULATES THE EIGENVALUES AND EIGENVECTORS OF
!  A SYMMETRIC 3 BY 3 MATRIX [A]. 
!
!  THE OUTPUT CONSISTS OF A VECTOR D CONTAINING THE THREE
!  EIGENVALUES IN ASCENDING ORDER, AND
!  A MATRIX [V] WHOSE COLUMNS CONTAIN THE CORRESPONDING
!  EIGENVECTORS.
! ---------------------------------------------------------------------
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(NP=3)
        DIMENSION D(NP),V(NP,NP)
        DIMENSION A(3,3),E(NP,NP)
        DO I = 1,3
          DO J= 1,3
            E(I,J) = A(I,J)
          ENDDO
        ENDDO
        CALL JACOBI(E,3,NP,D,V,NROT)
        CALL EIGSRT(D,V,3,NP)
        RETURN
        END
! *********************************************************************
        SUBROUTINE JACOBI(A,N,NP,D,V,NROT)

!  COMPUTES ALL EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
!  MATRIX [A], WHICH IS OF SIZE N BY N, STORED IN A PHYSICAL 
!  NP BY BP ARRAY. ON OUTPUT, ELEMENTS OF [A] ABOVE THE DIAGONAL 
!  ARE DESTROYED, BUT THE DIAGONAL AND SUB-DIAGONAL ARE UNCHANGED
!  AND GIVE FULL INFORMATION ABOUT THE ORIGINAL SYMMETRIC MATRIX.
!  VECTOR D RETURNS THE EIGENVALUES OF [A] IN ITS FIRST N ELEMENTS.
!  [V] IS A MATRIX WITH THE SAME LOGICAL AND PHYSICAL DIMENSIONS AS
!  [A] WHOSE COLUMNS CONTAIN, ON OUTPUT, THE NORMALIZED
!  EIGENVECTORSOF [A]. NROT RETURNS THE NUMBER OF JACOBI ROTATIONS
!  WHICH WERE REQUIRED.
!
!  THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", PAGE 346.
 
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (NMAX =100)
        DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)

!  INITIALIZE [V] TO THE IDENTITY MATRIX

        DO IP = 1,N     
          DO IQ = 1,N
            V(IP,IQ) = 0.D0
          ENDDO
          V(IP,IP) = 1.D0
        ENDDO

!  INITIALIZE B AND D TO THE DIAGONAL OF [A], AND Z TO ZERO. THE
!  VECTOR Z WILL ACCUMULATE TERMS OF THE FORM T*A_PQ AS
!  IN EQUATION (11.1.14)

        DO IP = 1,N
          B(IP) = A(IP,IP)
          D(IP) = B(IP)
          Z(IP) = 0.D0
        ENDDO
        NROT = 0
        DO I = 1,50

!  SUM OFF-DIAGONAL ELEMENTS

        SM = 0.D0
        DO IP = 1, N-1
          DO IQ = IP + 1, N
            SM = SM + DABS ( A(IP,IQ ))
          ENDDO
        ENDDO

!  IF SUM = 0., THEN RETURN. THIS IS THE NORMAL RETURN
!  WHICH RELIES ON QUADRATIC CONVERGENCE TO MACHINE 
!  UNDERFLOW.

        IF( SM .EQ. 0.D0) RETURN

!  IF( SM .LT. 1.0D-15) RETURN

!  IN THE FIRST THREE SWEEPS CARRY OUT THE PQ ROTATION ONLY IF
!  |A_PQ| > TRESH, WHERE TRESH IS SOME THRESHOLD VALUE, C
!  SEE EQUATION (11.1.25).  THEREAFTER TRESH = 0.

        IF( I .LT. 4) THEN
           TRESH = 0.2D0*SM/N**2
        ELSE
           TRESH = 0.D0
        END IF
        DO IP = 1, N-1
          DO IQ = IP+1,N
            G = 100.D0*DABS(A(IP,IQ))

!  AFTER FOUR SWEEPS, SKIP THE ROTATION IF THE
!  OFF-DIAGONAL ELEMENT IS SMALL.

        IF( (I .GT. 4) .AND. ( DABS(D(IP))+G  .EQ.  DABS( D(IP)) )
     +               .AND. ( DABS(D(IQ))+G  .EQ.  DABS( D(IQ)) ) )THEN
                       A(IP,IQ) = 0.D0
                   ELSE IF( DABS(A(IP,IQ)) .GT. TRESH) THEN
                       H = D(IQ) - D(IP)
                       IF(DABS(H)+G .EQ. DABS(H)) THEN

!                               T= 1./(2.*THETA), EQUATION(11.1.10)

                                T =A(IP,IQ)/H
                              ELSE
                                  THETA = 0.5D0*H/A(IP,IQ)
                             T =1.D0/(DABS(THETA)+DSQRT(1.D0+THETA**2))
                                  IF(THETA .LT. 0.D0) T = -T
                              END IF
                              C = 1.D0/DSQRT(1.D0 + T**2)
                              S = T*C
                              TAU = S/(1.D0 + C)
                                H = T*A(IP,IQ)
                                Z(IP) = Z(IP) - H
                                Z(IQ) = Z(IQ) + H
                                D(IP) = D(IP) - H
                                D(IQ) = D(IQ) + H
                                A(IP,IQ) = 0.D0

!                               CASE OF ROTATIONS 1<= J < P
                               
                                DO J = 1, IP-1
                                     G = A(J,IP)
                                     H = A(J,IQ)
                                     A(J,IP) = G - S*(H + G*TAU)
                                     A(J,IQ) = H + S*(G - H*TAU)
                                ENDDO

!                               CASE OF ROTATIONS P < J < Q

                                DO J = IP+1, IQ-1
                                     G = A(IP,J)
                                     H = A(J,IQ)
                                     A(IP,J) = G - S*(H + G*TAU)
                                     A(J,IQ) = H + S*(G - H*TAU)
                                ENDDO

!                               CASE OF ROTATIONS Q < J <= N

                                DO J = IQ+1, N
                                     G = A(IP,J)
                                     H = A(IQ,J)
                                     A(IP,J) = G - S*(H + G*TAU)
                                     A(IQ,J) = H + S*(G - H*TAU)
                                ENDDO
                                DO J = 1,N
                                     G = V(J,IP)
                                     H = V(J,IQ)
                                     V(J,IP) = G - S*(H + G*TAU)
                                     V(J,IQ) = H + S*(G - H*TAU)
                                ENDDO
                            NROT = NROT + 1
                   END IF
                   ENDDO
            ENDDO

!  UPDATE D WITH THE SUM OF T*A_PQ, AND REINITIALIZE Z

          DO IP = 1, N
               B(IP) = B(IP) + Z(IP)
               D(IP) = B(IP)
               Z(IP) = 0.D0
          ENDDO
        ENDDO

!  IF THE ALGORITHM HAS REACHED THIS STAGE, THEN
!  THERE ARE TOO MANY SWEEPS, PRINT A DIAGNOSTIC
!  AND STOP

!        WRITE(*,'(/1X,A/)') '50 ITERS IN JACOBI SHOULD NEVER HAPPEN'
        RETURN
        END
! **********************************************************************
        SUBROUTINE EIGSRT(D,V,N,NP)
!
!       GIVEN THE EIGENVALUES D AND EIGENVECTORS [V] AS OUTPUT FROM
!       JACOBI, THIS ROUTINE SORTS THE EIGENVALUES INTO ASCENDING ORDER, 
!       AND REARRANGES THE COLUMNS OF [V] ACCORDINGLY.
!
!       THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", P. 348.
! ----------------------------------------------------------------------
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION D(NP),V(NP,NP)
        DO I = 1,N-1
                K = I
                P = D(I)
                DO J = I+1,N
                     IF(D(J) .GE. P) THEN
                          K = J
                          P = D(J)
                     END IF
                ENDDO
                IF(K .NE. I) THEN
                   D(K) = D(I)
                   D(I) = P
                 DO J = 1,N
                        P = V(J,I)
                        V(J,I) = V(J,K)
                        V(J,K) = P
                 ENDDO
                END IF
        ENDDO
        RETURN
        END


      ! Reduces 4th order matrix to second order (symmetry in i-j / k-l)
      subroutine ematr(a, b)
        implicit none
        real*8 a(3, 3, 3, 3), b(6, 6)

        b(1, 1) = a(1, 1, 1, 1)
        b(2, 1) = a(2, 2, 1, 1)
        b(3, 1) = a(3, 3, 1, 1)
        b(4, 1) = a(1, 2, 1, 1)
        b(5, 1) = a(1, 3, 1, 1)
        b(6, 1) = a(2, 3, 1, 1)

        b(1, 2) = a(1, 1, 2, 2)
        b(2, 2) = a(2, 2, 2, 2)
        b(3, 2) = a(3, 3, 2, 2)
        b(4, 2) = a(1, 2, 2, 2)
        b(5, 2) = a(1, 3, 2, 2)
        b(6, 2) = a(2, 3, 2, 2)

        b(1, 3) = a(1, 1, 3, 3)
        b(2, 3) = a(2, 2, 3, 3)
        b(3, 3) = a(3, 3, 3, 3)
        b(4, 3) = a(1, 2, 3, 3)
        b(5, 3) = a(1, 3, 3, 3)
        b(6, 3) = a(2, 3, 3, 3)

        b(1, 4) = a(1, 1, 1, 2)
        b(2, 4) = a(2, 2, 1, 2)
        b(3, 4) = a(3, 3, 1, 2)
        b(4, 4) = a(1, 2, 1, 2)
        b(5, 4) = a(1, 3, 1, 2)
        b(6, 4) = a(2, 3, 1, 2)

        b(1, 5) = a(1, 1, 1, 3)
        b(2, 5) = a(2, 2, 1, 3)
        b(3, 5) = a(3, 3, 1, 3)
        b(4, 5) = a(1, 2, 1, 3)
        b(5, 5) = a(1, 3, 1, 3)
        b(6, 5) = a(2, 3, 1, 3)

        b(1, 6) = a(1, 1, 2, 3)
        b(2, 6) = a(2, 2, 2, 3)
        b(3, 6) = a(3, 3, 2, 3)
        b(4, 6) = a(1, 2, 2, 3)
        b(5, 6) = a(1, 3, 2, 3)
        b(6, 6) = a(2, 3, 2, 3)
      end subroutine

      ! Expands 2nd order to 4th order matrix
      subroutine ematrec(a, b)
        implicit none
        real*8 a(6, 6), b(3, 3, 3, 3)

        b(1, 1, 1, 1) = a(1, 1)
        b(2, 1, 1, 1) = a(4, 1)
        b(3, 1, 1, 1) = a(5, 1)
        b(1, 2, 1, 1) = a(4, 1)
        b(2, 2, 1, 1) = a(2, 1)
        b(3, 2, 1, 1) = a(6, 1)
        b(1, 3, 1, 1) = a(5, 1)
        b(2, 3, 1, 1) = a(6, 1)
        b(3, 3, 1, 1) = a(3, 1)

        b(1, 1, 2, 1) = a(1, 4)
        b(2, 1, 2, 1) = a(4, 4)
        b(3, 1, 2, 1) = a(5, 4)
        b(1, 2, 2, 1) = a(4, 4)
        b(2, 2, 2, 1) = a(2, 4)
        b(3, 2, 2, 1) = a(6, 4)
        b(1, 3, 2, 1) = a(5, 4)
        b(2, 3, 2, 1) = a(6, 4)
        b(3, 3, 2, 1) = a(3, 4)

        b(1, 1, 3, 1) = a(1, 5)
        b(2, 1, 3, 1) = a(4, 5)
        b(3, 1, 3, 1) = a(5, 5)
        b(1, 2, 3, 1) = a(4, 5)
        b(2, 2, 3, 1) = a(2, 5)
        b(3, 2, 3, 1) = a(6, 5)
        b(1, 3, 3, 1) = a(5, 5)
        b(2, 3, 3, 1) = a(6, 5)
        b(3, 3, 3, 1) = a(3, 5)

        b(1, 1, 1, 2) = a(1, 4)
        b(2, 1, 1, 2) = a(4, 4)
        b(3, 1, 1, 2) = a(5, 4)
        b(1, 2, 1, 2) = a(4, 4)
        b(2, 2, 1, 2) = a(2, 4)
        b(3, 2, 1, 2) = a(6, 4)
        b(1, 3, 1, 2) = a(5, 4)
        b(2, 3, 1, 2) = a(6, 4)
        b(3, 3, 1, 2) = a(3, 4)

        b(1, 1, 2, 2) = a(1, 2)
        b(2, 1, 2, 2) = a(4, 2)
        b(3, 1, 2, 2) = a(5, 2)
        b(1, 2, 2, 2) = a(4, 2)
        b(2, 2, 2, 2) = a(2, 2)
        b(3, 2, 2, 2) = a(6, 2)
        b(1, 3, 2, 2) = a(5, 2)
        b(2, 3, 2, 2) = a(6, 2)
        b(3, 3, 2, 2) = a(3, 2)

        b(1, 1, 3, 2) = a(1, 6)
        b(2, 1, 3, 2) = a(4, 6)
        b(3, 1, 3, 2) = a(5, 6)
        b(1, 2, 3, 2) = a(4, 6)
        b(2, 2, 3, 2) = a(2, 6)
        b(3, 2, 3, 2) = a(6, 6)
        b(1, 3, 3, 2) = a(5, 6)
        b(2, 3, 3, 2) = a(6, 6)
        b(3, 3, 3, 2) = a(3, 6)

        b(1, 1, 1, 3) = a(1, 5)
        b(2, 1, 1, 3) = a(4, 5)
        b(3, 1, 1, 3) = a(5, 5)
        b(1, 2, 1, 3) = a(4, 5)
        b(2, 2, 1, 3) = a(2, 5)
        b(3, 2, 1, 3) = a(6, 5)
        b(1, 3, 1, 3) = a(5, 5)
        b(2, 3, 1, 3) = a(6, 5)
        b(3, 3, 1, 3) = a(3, 5)

        b(1, 1, 2, 3) = a(1, 6)
        b(2, 1, 2, 3) = a(4, 6)
        b(3, 1, 2, 3) = a(5, 6)
        b(1, 2, 2, 3) = a(4, 6)
        b(2, 2, 2, 3) = a(2, 6)
        b(3, 2, 2, 3) = a(6, 6)
        b(1, 3, 2, 3) = a(5, 6)
        b(2, 3, 2, 3) = a(6, 6)
        b(3, 3, 2, 3) = a(3, 6)

        b(1, 1, 3, 3) = a(1, 3)
        b(2, 1, 3, 3) = a(4, 3)
        b(3, 1, 3, 3) = a(5, 3)
        b(1, 2, 3, 3) = a(4, 3)
        b(2, 2, 3, 3) = a(2, 3)
        b(3, 2, 3, 3) = a(6, 3)
        b(1, 3, 3, 3) = a(5, 3)
        b(2, 3, 3, 3) = a(6, 3)
        b(3, 3, 3, 3) = a(3, 3)
        
      end subroutine

      subroutine get_burg(props, nprops, nslip, burg)
        implicit none
        integer nprops, nslip
        real*8 :: props(nprops), burg1, burg2, burg3
        real*8 :: burg(nslip)
        integer :: idx
        
        burg1 = props(13)
        burg2 = props(43)
        burg3 = props(44)
        do idx = 1, 12
          burg(idx) = burg1
        end do
        do idx = 13, 15
          burg(idx) = burg2
        end do
        do idx = 16, 18
          burg(idx) = burg3
        end do

      end subroutine

      subroutine calc_amat2(nslip,props,nprops,s0,n0,amat)
        implicit none
        integer nslip, nprops
        real*8 amat(nslip, nslip), props(nprops), s0(3, nslip)
        real*8 n0(3, nslip), sa_dot_sb, na_dot_nb, t0(3), na_dot_tb
        integer i, j

        ! ActMat 2011 59 7003
        amat = 0.d0
        do j = 1, nslip
          do i = 1, nslip
            !sa_dot_sb = abs(dot_product(s0(:,i),s0(:,j)))
            na_dot_nb = abs(dot_product(n0(:,i),n0(:,j)))
            if (i.eq.j) then
              amat(i, j) = 0.1d0
            else
              if (na_dot_nb.gt.0.99d0) then
                amat(i, j) = 0.12d0
              else
                t0(1) = n0(2, j) * s0(3, j) - n0(3, j) * s0(2, j)
                t0(2) = n0(3, j) * s0(1, j) - n0(1, j) * s0(3, j)
                t0(3) = n0(1, j) * s0(2, j) - n0(2, j) * s0(1, j)
                na_dot_tb = abs(dot_product(n0(:, i), t0))
                amat(i, j) = 0.12d0 * na_dot_tb
              end if
            end if
          end do
        end do
      end subroutine calc_amat2

      subroutine calc_amat3(nslip,props,nprops,s0,n0,amat)
        implicit none
        integer nslip, nprops
        real*8 amat(nslip, nslip), props(nprops), s0(3, nslip)
        real*8 n0(3, nslip), sa_dot_sb, na_dot_nb, t0(3), na_dot_tb
        integer i, j

        ! IJP 2016 80 222
        amat = 0.d0
        do j = 1, nslip
          do i = 1, nslip
            sa_dot_sb = abs(dot_product(s0(:,i),s0(:,j)))
            na_dot_nb = abs(dot_product(n0(:,i),n0(:,j)))
            if ((na_dot_nb < 0.999d0) .and. 
     &               (sa_dot_sb > 0.999d0)) then ! Collinear
              amat(i, j) = 0.7d0
            else ! Noncollinear
              amat(i, j) = 0.1d0
            end if
          end do
        end do
      end subroutine calc_amat3

      subroutine calc_amat4(nslip,props,nprops,s0,n0,amat)
        implicit none
        integer nslip, nprops
        real*8 amat(nslip, nslip), props(nprops), s0(3, nslip)
        real*8 n0(3, nslip), sa_dot_sb, na_dot_nb, t0(3), na_dot_tb
        integer i, j

        amat = 0.d0

        amat(1:12, 1:12) = 1.d0
        amat(13:15, 13:15) = 1.d0
        amat(16:18, 16:18) = 1.d0

        amat(13:15, 1:12) = 0.5d0
        amat(1:12, 13:15) = 0.5d0

        amat(16:18, 13:15) = 0.5d0
        amat(13:15, 16:18) = 0.5d0

        amat(1:12, 16:18) = 0.5d0
        amat(16:18, 1:12) = 0.5d0

      end subroutine

      ! Computes euler angles
      subroutine eulang(q,th,phi,om,ieulererr)
        implicit none
        ! Argument variables
        integer ieulererr
        real*8 q(3, 3), th, phi, om

        ! Local variables
        integer ichk, i, j
        real*8 pi, sth

        pi = 4.0*datan(1.0d0)
        ichk=0
        ieulererr = 0
        if(dabs(q(3,3))-1.0 > 1.0d-6)then
           ieulererr = 1
           return
        endif

        do i = 1,3
          do j = 1,3
            if(dabs(q(i,j)) < 1.0d-6) q(i,j) = 0.0
          enddo
        enddo

        if(dabs(dabs(q(3,3))-1.0).lt.1.0d-6)then
          call eulcheck1(q,th,phi,om,ichk)
          if(ichk /= 1)go to 20
          return
        endif

        th = dacos(q(3,3))
        sth = dsin(th)
        om = datan2(q(1,3)/sth,q(2,3)/sth)
        phi = datan2(q(3,1)/sth,-q(3,2)/sth)
        call eulcheck(q,th,phi,om,ichk)
        if(ichk == 1)return
        th = 2.0*pi-th
        sth = dsin(th)
        om = datan2(q(1,3)/sth,q(2,3)/sth)
        phi = datan2(q(3,1)/sth,-q(3,2)/sth)
        call eulcheck(q,th,phi,om,ichk)
20      if(ichk /= 1)then
            ieulererr = 1
            return
        endif

        return
      end
!===================================================================
      subroutine eulcheck(q,th,phi,om,ichk)
        implicit none
        ! Argument variables
        integer ichk
        real*8 q(3, 3), th, phi, om

        ! Local variables
        real*8 tol, a, b, c, d

        tol=1.0d-5
        ichk = 0

        a = dcos(phi)*dcos(om)-dsin(phi)*dsin(om)*dcos(th)
        b = -dsin(om)*dcos(phi)-dcos(om)*dsin(phi)*dcos(th)
        c = dcos(om)*dsin(phi)+dsin(om)*dcos(phi)*dcos(th)
        d = -dsin(phi)*dsin(om)+dcos(phi)*dcos(om)*dcos(th)

        if(dabs(a-q(1,1)).lt.tol.and.dabs(b-q(2,1)).lt.tol.
     &    and.dabs(c-q(1,2)).lt.tol.and.dabs(d-q(2,2)).lt.tol)ichk=1

        return
      end
c **********************************************************************
      subroutine eulcheck1(q,th,phi,om,ichk)
        implicit none
        ! Argument variables
        integer ichk
        real*8 q(3, 3), th, phi, om

        ! Local variables
        real*8 tol
        tol=1.0d-5
        ichk = 0

        q(3,3) = 1.0*q(3,3)/dabs(q(3,3))
        th = dacos(q(3,3))
        if(dabs(q(1,3)).gt.tol)return
        if(dabs(q(2,3)).gt.tol)return
        if(dabs(q(3,1)).gt.tol)return
        if(dabs(q(3,2)).gt.tol)return
        if(q(3,3).eq.1.0)then
           if(dabs(q(1,1)-q(2,2)).gt.tol)return
           if(dabs(q(1,2)+q(2,1)).gt.tol)return
        else
           if(dabs(q(1,1)+q(2,2)).gt.tol)return
           if(dabs(q(1,2)-q(2,1)).gt.tol)return
        endif
        phi = datan2(q(1,2),q(1,1))
        om = 0.0
        ichk = 1

        return
      end