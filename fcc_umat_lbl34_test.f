!==================================================================
!
!  This is a version of the Kalidindi/Bronkhorst/Anand UMAT that has
!  been updated to remove COMMON blocks so that it can be run in 
!  parallel. 
!
!  2/12/2019 Charles Lieou
!  Modifications have been made to incorporate the LBL theory of 
!  dislocation plasticity.
!
!  Additional code modifications/cleanup efforts have also been carried
!  out in the updating of the code, some of which has been documented
!  elsewher.
!
!  Uses total strain rate instead of plastic strain rate in the ``xnu''
!  factor.
!
!  Uses a more moderate shear modulus
!------------------------------------------------------------------

      subroutine umat_crys2 (stress,  statev,  ddsdde,  sse,  spd,
     &          scd,     rpl,     ddsddt,  drplde,  drpldt,
     &          stran,  dstran, time,    dtime,   temp,
     &          dtemp,   predef,  dpred,   cmname,  ndi,
     &          nshr,    ntens,   nstatv,  props,   nprops,
     &          coords,  drot,    pnewdt,  celent,  dfgrd0,
     &          dfgrd1,  noel,    npt,     layer,   kspt,
     &          kstep,   kinc)

      implicit real*8(a-h,o-z)
      character*80 cmname 

!------------------------------------------------------------------
!  Dimension UMAT argument arrays
!------------------------------------------------------------------

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

!------------------------------------------------------------------           
!  Dimension UMAT local arrays
!------------------------------------------------------------------           

      dimension
     &  array4  (6,6),       ! Full material Jacobian to be passed back to Abaqus 
     &  C_c  (3,3,3,3),      ! Elastic stiffness tensor in crystal basis
     &  C   (3,3,3,3),       ! Elastic stiffness tensor in reference config
     &  del (3,3),           ! Kronecker delta tensor
     &  q_0 (3,3),           ! Rotation tensor from crystal to reference config 
     &  Fe_n(3,3),           ! Elastic deformation map at beginning of step
     &  FeInv_n(3,3),        ! Inverse of Fe_n
     &  Fe    (3,3),         ! Elastic deformation map
     &  FpInv_n(3,3),        ! Inverse plastic deformation map at beginning of step
     &  FpInv(3,3),          ! Inverse plastic deformation map at end of step
     &  s_n  (props(34)),    ! slip resistance (CRSS) at beginning of step*
     &  s   (props(34)),     ! slip resistance (CRSS) at end of step*
     &  gammadot(props(34)), ! Slip strain rate array
     &  eulangs (3),         ! Euler angles
     &  sig (3,3),           ! Cauchy stress tensor
     &  Spk2    (3,3),       ! 2nd Piola-Kirchhoff stress
     &  tau (props(34)),     ! RSS
     &  xs0 (3,props(34)),   ! Slip vectors in intermediate config
     &  xm0 (3,props(34)),   ! Slip normal vectors in intermediate config
     &  xmC   (3,props(34)),   ! Slip normal vectors in crystal basis 
     &  xsC   (3,props(34)),   ! Slip vectors in crystal basis
     &  del4th(3,3,3,3),     ! Fourth rank symmetric identity tensor
     &  Fe_tr(3,3),          ! Trial elastic deformation map
     &  Ce_tr(3,3),          ! Trial elastic right stretch tensor
     &  Ee_tr(3,3),          ! Trial elastic Green strain tensor
     &  Spk2_tr(3,3),        ! Trial PK2 stress tensor
     &  Smat(3,3,props(34)), ! Schmid tensor matrix in reference config
     &  Bmat(3,3,props(34)), ! B matrix
     &  Cmat(3,3,props(34)), ! C matrix
     &  dfgrdinv0(3,3),      ! Inverse of deformation gradient at beginning of step
     &  rLBar(3,3),          ! Deformation tensor times time increment
     &  rLBartemp(3,3),      ! Temporary matrix to store intermediate product involving rLBar
     &  rLpBar(3,3),         ! Plastic velocity gradient in intermediate config
     &  smatc(props(34),3,3),! Schmid tensor matrix in crystal basis 
     &  spk2v_n(6),          ! PK2 stress vector at beginning of step
     &  spk2_n(3,3),         ! PK2 stress tensor at beginning of step
     &  spk2v_tr(6),         ! Trial PK2 stress vector
     &  pmatr(props(34),6),  ! Symmetric part of Schmid tensor in reduced form
     &  cmatr(props(34),6),  ! C matrix in reduced form
     &  spk2v(6),            ! PK2 stress vector at end of step
     &  dgamma(props(34)),   ! Plastic slip increment vector
     &  pmat(3,3,props(34)), ! Symmetric part of the Schmid tensor
     &  Re(3,3),             ! Elastic rotation tensor
     &  Ue(3,3),             ! Right elastic stretch tensor
     &  Q_n1(3,3),           ! Total lattice rotation tensor from crystal to current config

! Arrays added for dislocation density hardening law
     &  rhossd_n(props(34)), ! SSD density array at beginning of step
     &  rhossd(props(34)),   ! SSD density array at end of step
     &  drhossd(props(34)),  ! SSD density array increment over step
     &  rhossd_i(props(34)), ! SSD coupling array at beginning of step
     &  rhomf_n(props(34)),  ! MFP dislocation density array at beginning of step
     &  rhomf(props(34)),    ! MFP dislocation density array at end of step
     &  xnu(props(34)),      ! Hardening factor xnu for debugging
     &  DGT(props(34)),      ! Resolved shear strain
     &  amat(props(34),props(34)), ! Slip interaction matrix appearing in definition of CRSS
     &  dmat(props(34),props(34)) ! Slip interaction matrix appearing in definition of mean free path

      pi = 4.d0 * atan(1.d0)

!------------------------------------------------------------------
!  Assign props() array to logical variable names
!------------------------------------------------------------------

        THETA0=PROPS(1)              ! Initial temperature [K]
        RHO=PROPS(2)                 ! Initial density [Mg/mm^3]
        CP=PROPS(3)                  ! Specific Heat [mJ/(Mg.K)]
        A11=PROPS(4)                 ! Thermal expansion coefficient [1/K]
        xkappax=PROPS(5)             ! Inverse specific effective heat [-]

        rM11=PROPS(6)                ! Slope of C11 vs temp line [MPa/K]
        C11_0=PROPS(7)               ! C11 at 0 K [MPa]
        rM12=PROPS(8)                ! Slope of C12 vs temp line [MPa/K]
        C12_0=PROPS(9)               ! C12 at 0 K [MPa]
        rM44=PROPS(10)               ! Slope of C44 vs temp line [MPa/K]
        C44_0=PROPS(11)              ! C44 at 0 K [MPa]

        alphaT=PROPS(12)             ! Taylor factor [-]** Neglect this?
        burg = props(13)             ! Burgers vector [mm]
        rhossd_0 = props(14)         ! Initial SSD density [mm^-2]
        chi0 = props(15)             ! Initial effective temperature [-]*
        chi_ss = props(16)           ! Steady-state effective temperature [-]*
        a_self = props(17)           ! Dislocation interaction coefficient [-]
        a_copl = props(18)           ! Dislocation interaction coefficient [-]
        a_hirth = props(19)          ! Dislocation interaction coefficient [-]
        a_colli = props(20)          ! Dislocation interaction coefficient [-]
        a_gliss = props(21)          ! Dislocation interaction coefficient [-]
        a_lomer = props(22)          ! Dislocation interaction coefficient [-]
        e_D = props(23)              ! Dislocation formation energy [-]*
        e_copl = props(24)           ! Dislocation interaction energy [-]*

        t0=PROPS(25)                 ! Time scale [s]
        THETAP=PROPS(26)             ! Activation temperature [K] 
        BK=PROPS(27)                 ! Boltzmann's const. [J/K]
        eta = props(28)              ! Dislocation drag [MPa]*
        xkappa0=PROPS(29)            ! Dislocation storage coefficient [-]*
        xkappa2=PROPS(30)            ! Dislocation storage coefficient [-]*
        th = props(31)*pi/180.d0     ! Euler angle theta [rad] (Canova convention)
        phi = props(32)*pi/180.d0    ! Euler angle Phi [rad]
        om = props(33)*pi/180.d0     ! Euler angle omega [rad]
        nslip=int(PROPS(34))         ! Number of slip systems

        gdotr = props(35)            ! Reference strain rate for hardening [s^-1]*
        d_gr  = props(36)            ! Grain diameter [mm]
        rx = props(37)               ! Ratio a/b (length scale to Burger's vector)
        rk_inter = props(38)         ! Dislocation trapping const for latent interactions [-]
        rk_copl =  props(39)         ! Dislocation trapping const for self interactions [-]
        rs = props(40)               ! Scaling factor for Taylor stress
        xkappar = props(41)          ! Dislocation energy scaling for chi recovery [mJ]
        xmu0 = props(42)             ! Polycrystal shear modulus adjustment [MPa]
   
!------------------------------------------------------------------
!  Initialize identity matrices/tensors
!------------------------------------------------------------------
      del = reshape( (/1.,0.,0.,0.,1.,0.,0.,0.,1. /), (/3,3/) )

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          del4th(i,j,k,l) = 0.5d0 * 
     &      (del(i,k) * del(j,l) + del(i,l) * del(j,k))
         end do
        end do
       end do
      end do

!------------------------------------------------------------------
!  Assign slip system normals and slip directions for an FCC.
!           plane normal           slip direction
!------------------------------------------------------------------

        xmC(:,1)= (/-1,  1,  1/);  xsC(:,1)=(/0, -1, 1/)   !A2
        xmC(:,2)= (/-1, 1, 1/);    xsC(:,2)=(/1, 0, 1/)    !A3
        xmC(:,3)= (/-1,  1,  1/);  xsC(:,3)=(/1, 1, 0/)    !A6
        xmC(:,4)= (/1,  1,  1/);   xsC(:,4)=(/0, -1,  1/)  !B2
        xmC(:,5)= (/1,  1,  1/);   xsC(:,5)=(/-1,  0, 1/)  !B4   
        xmC(:,6)= (/1 , 1,  1/);   xsC(:,6)=(/-1,  1, 0/)  !B5
        xmC(:,7)= (/-1,  -1, 1/);  xsC(:,7)=(/0, 1,  1/)   !C1
        xmC(:,8)= (/-1,  -1, 1/);  xsC(:,8)=(/1,  0,  1/)  !C3
        xmC(:,9)= (/-1,  -1, 1/);  xsC(:,9)=(/-1,  1,  0/) !C5
        xmC(:,10)=(/1, -1,  1/);   xsC(:,10)=(/0,  1,  1/) !D1
        xmC(:,11)=(/1, -1,  1/);   xsC(:,11)=(/-1,  0, 1/) !D4
        xmC(:,12)=(/1, -1,  1/);   xsC(:,12)=(/1,  1,  0/) !D6

!  Normalize the vectors
      do i = 1,nslip
        xmCmag = sqrt(dot_product(xmC(:,i),xmC(:,i)))
        xmC(:,i) = xmC(:,i)/xmCmag
        xsCmag = sqrt(dot_product(xsC(:,i),xsC(:,i)))
        xsC(:,i) = xsC(:,i)/xsCmag
      enddo

!------------------------------------------------------------------
!  Initialize state at t = t_n
!------------------------------------------------------------------

      if (time(2) .eq. 0.0) then ! t = 0 
        
        theta = props(1)
        eulangs(1) = th 
        eulangs(2) = phi 
        eulangs(3) = om
        rhossd_n(:) = rhossd_0
        FpInv_n = del
        spk2_n = 0.d0
        chi_n = chi0

      else  ! time /= 0

        theta = statev(1)

        eulangs(1:3) = statev(2:4)

        n = 4
        do j = 1,3
          do i = 1,3
            n = n + 1
            FpInv_n(i,j) = statev(n)
          end do
        end do

        spk2_n(1,1) = statev(14)
        spk2_n(2,2) = statev(15)
        spk2_n(3,3) = statev(16)
        spk2_n(1,2) = statev(17)
        spk2_n(1,3) = statev(18)
        spk2_n(2,3) = statev(19)
        
        n = 19
        do i = 1,nslip
          n = n + 1
          rhossd_n(i) = statev(n)
        end do

        chi_n = statev(46)

      end if ! state initialization block

!  Initialize 4th rank elastic stiffness tensor
      call compute_elastConst(props,nprops,theta,c11,c12,c44)
      
      C_c(:,:,:,:) = 0.d0
      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          C_c(i,j,k,l) = C12 * del(i,j) * del(k,l) +
     &        C44 * (del(i,k)*del(j,l)+del(i,l)*del(k,j))
         end do
        end do
       end do
      end do
      C_c(1,1,1,1) = C11
      C_c(2,2,2,2) = C11
      C_c(3,3,3,3) = C11

!  Rotate stiffness tensor and slip vectors from crystal to reference configuration
      call ROTMAT(TH,PHI,OM,q_0)
      call rotate_4th(q_0,C_c,C)
      xs0 = matmul(q_0,xsC)
      xm0 = matmul(q_0,xmC)
      
!  Compute trial quantities
      Fe_tr = matmul(dfgrd1,FpInv_n)
      Ce_tr = matmul(transpose(Fe_tr),Fe_tr)
      Ee_tr = 0.5d0*(Ce_tr-del) - a11*(theta - theta0)*del
      Spk2_tr = reshape( (/((sum(C(i,j,:,:)*Ee_tr), i = 1,3), j=1,3)/), (/3,3/))

      do i = 1,3
        do j = 1,3
             Smat(i,j,:) = xs0(i,:)*xm0(j,:)
        enddo
      enddo

      do i = 1,3
        do j = 1,3
             Pmat(i,j,:) = 0.5d0*(Smat(i,j,:)+Smat(j,i,:))
        enddo
      enddo

      Bmat(:,:,:) = 0.d0
      do ia = 1,nslip
        Bmat(:,:,ia) = matmul(Ce_tr,Smat(:,:,ia)) +
     &                 matmul(transpose(Smat(:,:,ia)),Ce_tr)
      enddo

      Cmat(:,:,:) = 0.d0
      do i = 1,3
        do j = 1,3
          do m = 1,3
          do n = 1,3
            Cmat(i,j,:) = Cmat(i,j,:) + 
     &                        0.5d0*C(i,j,m,n)*Bmat(m,n,:)
          enddo
          enddo
        enddo
      enddo

!  Package arrays in a manner consistent with what the solver is
!  expecting.  Solver was adapted from Bronkhorst's UMAT.
        spk2_n(2,1) = spk2_n(1,2)
        spk2_n(3,1) = spk2_n(1,3)
        spk2_n(3,2) = spk2_n(2,3) 

        call smatvec(spk2_n,spk2v_n)
        call smatvec(spk2_tr,spk2v_tr)
        do ia = 1,nslip
          call smatvec(cmat(:,:,ia),cmatr(ia,:)) 
          call smatvec(pmat(:,:,ia),pmatr(ia,:))
        enddo
        pmatr(:,4) = 2.d0*pmatr(:,4) 
        pmatr(:,5) = 2.d0*pmatr(:,5) 
        pmatr(:,6) = 2.d0*pmatr(:,6)
        
!  Compute slip resistances from dislocation densities at t_n
      ! xmu = sqrt((c11-c12)/2.d0*c44)
      xmu = sqrt((c11-c12)/2.d0*c44) + xmu0
      call calc_amat(nslip,props,nprops,xsC,xmC,amat)
      call calc_dmat(nslip,props,nprops,xmC,amat,dmat)
      s_n = (xmu/rs)*burg*sqrt(matmul(amat,rhossd_n))
      rhomf_n = matmul(dmat,rhossd_n)
      rhossd_i = matmul(dmat,rhossd_n)

!  Compute the total strain, needed for hardening
       call inverse_3x3(dfgrd0,dfgrdinv0)
       Fe_n = matmul(dfgrd0,FpInv)
       call inverse_3x3(Fe_n,FeInv_n)
       rLBar = matmul(dfgrd1-dfgrd0,dfgrdinv0)
       rLBartemp = matmul(rLBar,Fe_n)
       rLBar = matmul(FeInv_n,rLBartemp)
       DGT = 0.d0
       do k = 1,nslip
         do i = 1,3
            do j = 1,3
               DGT(k) = DGT(k) + rLBar(i,j)*Pmat(i,j,k)
            enddo
         enddo
       enddo
       ! if (noel == 1) then
       !    WRITE(*,*) DGT/dtime
       ! endif

!  Call solver      
      dgamma = 0.d0 

      call SOLVCRYS(nslip,props,nprops,theta,dfgrd0,dfgrd1,spk2v_n,s_n,amat,dmat,
     &  rhossd_n,rhomf_n,chi_n,spk2v_tr,pmatr,cmatr,Pmat,DTIME,spk2v,s,
     &  rhossd,drhossd,rhomf,xnu,dgamma,DGT,DGMAX,tau,DOMPSC,ITERK,ITERL,ITERERR,xmC)

!  Solver did not converge - cut back and repeat time step
        if (iterErr == 1) then
          pnewdt = 0.5d0
          return
        endif

!  Unpack solver quantities needed for the following state updates
!  and Jacobian calculation
        gammadot = dgamma / dtime
        
        call svec_to_smat(spk2v,spk2)

!  Update Fp
      rLpBar = 0.d0
      do i = 1,3
        do j = 1,3
          do k = 1,nslip
            rLpBar(i,j) = rLpBar(i,j) + 
     &          Smat(i,j,k) * gammadot(k)
          end do
        end do
      end do

      FpInv = matmul(FpInv_n,del-dtime*rLpBar)
      detFp = determinant(FpInv)
      FpInv = FpInv / detFp**(1.d0/3.d0)

!  Compute Fe
      Fe = matmul(dfgrd1,FpInv)

!  Compute Cauchy stress
      detFe = determinant(Fe)
      sig = matmul(matmul(Fe,Spk2),transpose(Fe)) / detFe

!  Compute Jacobian - passed back in array4 in Voigt notation      
      do i = 1,3
      do j = 1,3
        smatc(:,i,j) = xsC(i,:)*xmC(j,:)
      enddo
      enddo

      call calc_jac(nslip,props,nprops,C,dfgrd0,dfgrd1,FpInv_n,Spk2,
     & gammadot,tau,s,rhossd,rhomf,q_0,smatc,theta,del,dtime,array4)

!  Compute the updated temperature and effective temperature
!  Include recovery term for chi here
      dWp = dot_product(tau,dgamma)
      dTheta = (chi_n/chi_ss)*dWp/(rho*cp)
      theta = theta + dTheta
      dchi = xkappax*dWp*(1.d0 - chi_n/chi_ss)/(alphaT*xmu)
      dchi2 = - dot_product(e_D+e_copl*(burg**2)*rhossd_i,drhossd)
     &          *xkappax*xkappar/(rx*burg*alphaT*xmu)
      chi = chi_n + dchi + dchi2
      ! WRITE(*,*) dchi2
      ! WRITE(*,*) time, dtime
      ! WRITE(*,*) s
      ! WRITE(*,*) xkappax, dWp, chi_ss

!  Compute Euler angles at t = t_(n+1) for texture plotting
      call polardecomp(Fe,Re,Ue)
      Q_n1 = matmul(Re,q_0)
      call eulang(Q_n1,th_n1,phi_n1,om_n1,ieulererr)
      if (ieulererr == 1) then
        eulangs(1)  = 0.d0
        eulangs(2) =  0.d0
        eulangs(3)  = 0.d0
      else
        eulangs(1)  = th_n1*180.d0/pi
        eulangs(2) =  phi_n1*180.d0/pi
        eulangs(3)  = om_n1*180.d0/pi
      endif

!  Store the stress tensor in the ABAQUS stress 'vector'
      do i = 1,ndi
         stress(i) = sig(i,i)
      end do
      if (nshr == 1) stress(ndi+1) = sig(1,2)
      if (nshr == 3) then
         stress(4) = sig(1,2)
         stress(5) = sig(1,3)
         stress(6) = sig(2,3)
      end if
      
!  Pass back the Jacobian that Abaqus is looking for
      ddsdde = 0.0

      if (ndi == 1) then          ! 1-D
         ddsdde(1,1) = array4(1,1)
      else if (ndi == 2) then         ! 2-D plane stress & axi
         ddsdde(1:2,1:2) = array4(1:2,1:2)
         ddsdde(1,3) = array4(1,4)
         ddsdde(2,3) = array4(2,4)
         ddsdde(3,1) = array4(4,1)
         ddsdde(3,2) = array4(4,2)
         ddsdde(3,3) = array4(4,4)
      else if (ndi == 3 .and. nshr == 1) then ! plane strain
         ddsdde(1:4,1:4) = array4(1:4,1:4)
      else                  ! Full 3-D
         ddsdde(1:6,1:6) = array4(1:6,1:6)
      end if
!------------------------------------------------------------------
!  Auxillary post-processing calcuations
!------------------------------------------------------------------
! Calculate the misorientation magnitude
      arg = (Re(1,1)+Re(2,2)+Re(3,3)-1.d0)/2.d0
      if (abs(arg-1.d0) < 1e-10) then
       angMis = 0.d0
      else
       angMis = dacos(arg)
      endif
      angMis = angMis*180.d0/pi

!------------------------------------------------------------------
!  Store the state variables
!------------------------------------------------------------------
         n = 0
!  Store the temperature
         n = n + 1
         statev(n) = theta

!  Store the Euler Angles
         do i = 1,3
           n = n + 1
           statev(n) = eulangs(i)
         end do
            
!  Store {F^p}^{-1}
            do j = 1,3
               do i = 1,3
                 n = n + 1
                 statev(n) = FpInv(i,j)
               end do
            end do

!  Store the PK2 stress            
        statev(14)=spk2(1,1)
        statev(15)=spk2(2,2)
        statev(16)=spk2(3,3)
        statev(17)=spk2(1,2)
        statev(18)=spk2(1,3)
        statev(19)=spk2(2,3)

!  Store s^alpha
           n = 19
           do i = 1,nslip
             n = n + 1
             statev(n) = rhossd(i)
           end do

!  Store the cumulative slips, gamma^alpha
           cumslip_all = 0.d0
           do i = 1,nslip
             n = n + 1
             statev(n) = statev(n) + abs(gammadot(i)*dtime)
             cumslip_all = cumslip_all + statev(n)
           end do

!  Store the cumulative slip on all systems
           n = n + 1
           statev(n) = cumslip_all

!  Store the misorientation angle
           n = n + 1
           statev(n) = angMis

!  Store the effective temperature chi
           n = n + 1
           statev(n) = chi

!  Store the slip resistances
           do i = 1,nslip
              n = n + 1
              statev(n) = s(i)
           enddo

!  Store the MFP dislocation densities
           do i = 1,nslip
              n = n + 1
              statev(n) = rhomf(i)
           enddo

!  Store hardening factor xnu
           do i = 1,nslip
              n = n + 1
              statev(n) = xnu(i)
           enddo

!  Store the slip strain rates
           do i = 1,nslip
              n = n + 1
              statev(n) = gammadot(i)
           enddo

!  Store the resolved strain rates
           do i = 1,nslip
              n = n + 1
              statev(n) = DGT(i)/DTIME
           enddo

       return
       end
     
!===================================================================
!===================================================================
!===================== S U B R O U T I N E S =======================
!===================================================================
!===================================================================
!  Store only select outputs from the statev array 
!-------------------------------------------------------------------
      subroutine uvarm(uvar,direct,t,time,dtime,cmname,orname,
     & nuvarm,noel,npt,layer,kspt,kstep,kinc,ndi,nshr,coord,
     & jmac,jmatyp,matlayo,laccfla)
      
      implicit real*8(a-h,o-z)
c
      character*80 cmname,orname
      character*3 flgray(106)
      dimension uvar(nuvarm),direct(3,3),t(3,3),time(2)
      dimension array(106),jarray(106),jmac(*),jmatyp(*),coord(*)

      call getvrm('SDV',array,jarray,flgray,jrcd,jmac,jmatyp,matlayo,
     &            laccfla)

      if (cmname(1:2) == 'CU') then
        uvar(1)    = array(1)
        uvar(2)    = array(44)
        uvar(3)    = array(45)
        uvar(4:6)  = array(2:4) 
        uvar(7:18) = array(20:31)
        uvar(19)   = array(46)
        uvar(20:31) = array(47:58)
        uvar(32:43) = array(59:70)
        uvar(44:55) = array(71:82)
        uvar(56:67) = array(83:94)
        uvar(68:79) = array(95:106)
      endif

      return
      end subroutine
!===================================================================
!  Convert from Voigt notation stress vector to a symmetric 
!-------------------------------------------------------------------
      subroutine svec_to_smat(stressv,stressm)

      implicit none

      real*8 :: stressv(6),stressm(3,3)

      stressm(1,1) = stressv(1)
      stressm(2,2) = stressv(2)
      stressm(3,3) = stressv(3)
      stressm(1,2) = stressv(4)
      stressm(2,1) = stressv(4)
      stressm(1,3) = stressv(5)
      stressm(3,1) = stressv(5)
      stressm(2,3) = stressv(6)
      stressm(3,2) = stressv(6)

      return
      end subroutine
!===================================================================
!===================================================================
!
!  Rotates any 3x3x3x3 tensor by a rotation matrix.
!
!  c(i,j,k,l) = a(i,m) * a(j,n) * a(k,p) * a(l,q) * b(m,n,p,q)
!
!-------------------------------------------------------------------

      subroutine rotate_4th(a,b,c)

      implicit double precision (a-h,o-z)

      dimension a(3,3), b(3,3,3,3), c(3,3,3,3), d(3,3,3,3)

      do m = 1,3
       do n = 1,3
        do k = 1,3
         do l = 1,3
          d(m,n,k,l) = a(k,1) * (a(l,1) * b(m,n,1,1) + 
     &      a(l,2) * b(m,n,1,2) + a(l,3) * b(m,n,1,3)) +
     &      a(k,2) * (a(l,1) * b(m,n,2,1) + 
     &      a(l,2) * b(m,n,2,2) + a(l,3) * b(m,n,2,3)) +
     &      a(k,3) * (a(l,1) * b(m,n,3,1) + 
     &      a(l,2) * b(m,n,3,2) + a(l,3) * b(m,n,3,3))
         end do
        end do
       end do
      end do

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          c(i,j,k,l) = a(i,1) * (a(j,1) * d(1,1,k,l) + 
     &      a(j,2) * d(1,2,k,l) + a(j,3) * d(1,3,k,l)) +
     &      a(i,2) * (a(j,1) * d(2,1,k,l) + 
     &      a(j,2) * d(2,2,k,l) + a(j,3) * d(2,3,k,l)) +
     &      a(i,3) * (a(j,1) * d(3,1,k,l) + 
     &      a(j,2) * d(3,2,k,l) + a(j,3) * d(3,3,k,l))
         end do
        end do
       end do
      end do

      return
      end

!===================================================================
!===================================================================
!
!  Calculate the inverse of a 3 x 3 matrix.
!
!-------------------------------------------------------------------

      subroutine inverse_3x3(a,b)

      implicit double precision (a-h,o-z)

      dimension a(3,3), b(3,3)

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

      return
      end
!===================================================================
!  Calculate the determinant of a 3 x 3 matrix.
!

      function determinant(a)

      implicit double precision (a-h,o-z)

      dimension a(3,3)

      b1 = a(2,2) * a(3,3) - a(3,2) * a(2,3)
      b2 = a(3,1) * a(2,3) - a(2,1) * a(3,3)
      b3 = a(2,1) * a(3,2) - a(3,1) * a(2,2)

      determinant = a(1,1) * b1 + a(1,2) * b2 + a(1,3) * b3

      return
      end



!===================================================================
!===================================================================
!======================================================================
!***********************************************************************
        SUBROUTINE ROTMAT(TH,PHI,OM,Q)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION Q(3,3)

        STH = DSIN(TH)
        CTH = DCOS(TH)
        SPH = DSIN(PHI)
        CPH = DCOS(PHI)
        SOM = DSIN(OM)
        COM = DCOS(OM)

        Q(1,1) = CPH*COM-SPH*SOM*CTH
        Q(1,2) = SPH*COM+CPH*SOM*CTH
        Q(1,3) = SOM*STH

        Q(2,1) = -CPH*SOM-SPH*COM*CTH
        Q(2,2) = -SPH*SOM+CPH*COM*CTH
        Q(2,3) =  COM*STH

        Q(3,1) = SPH*STH
        Q(3,2) = -CPH*STH
        Q(3,3) = CTH

        RETURN
        END
      
!===================================================================
        SUBROUTINE EULANG(Q,TH,PHI,OM,IEULERERR)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION Q(3,3)
        PI = 4.0*DATAN(1.0D0)
        ICHK=0
        IEULERERR = 0
        IF(DABS(Q(3,3))-1.0 > 1.0D-6)THEN
           IEULERERR = 1
           RETURN
        ENDIF
        DO I = 1,3
        DO J = 1,3
          IF(DABS(Q(I,J)) < 1.0D-6) Q(I,J) = 0.0
        ENDDO
        ENDDO
        IF(DABS(DABS(Q(3,3))-1.0).LT.1.0D-6)THEN
          CALL EULCHECK1(Q,TH,PHI,OM,ICHK)
          IF(ICHK /= 1)GO TO 20
          RETURN
        ENDIF
        TH = DACOS(Q(3,3))
        STH = DSIN(TH)
        OM = DATAN2(Q(1,3)/STH,Q(2,3)/STH)
        PHI = DATAN2(Q(3,1)/STH,-Q(3,2)/STH)
        CALL EULCHECK(Q,TH,PHI,OM,ICHK)
        IF(ICHK == 1)RETURN
        TH = 2.0*PI-TH
        STH = DSIN(TH)
        OM = DATAN2(Q(1,3)/STH,Q(2,3)/STH)
        PHI = DATAN2(Q(3,1)/STH,-Q(3,2)/STH)
        CALL EULCHECK(Q,TH,PHI,OM,ICHK)
20      IF(ICHK /= 1)THEN
            IEULERERR = 1
            RETURN
        ENDIF

        RETURN
        END
!===================================================================
        SUBROUTINE EULCHECK(Q,TH,PHI,OM,ICHK)

        IMPLICIT REAL*8(A-H,O-Z)

        DIMENSION Q(3,3)

        TOL=1.0D-3
        ICHK = 0

        A = DCOS(PHI)*DCOS(OM)-DSIN(PHI)*DSIN(OM)*DCOS(TH)
        B = -DSIN(OM)*DCOS(PHI)-DCOS(OM)*DSIN(PHI)*DCOS(TH)
        C = DCOS(OM)*DSIN(PHI)+DSIN(OM)*DCOS(PHI)*DCOS(TH)
        D = -DSIN(PHI)*DSIN(OM)+DCOS(PHI)*DCOS(OM)*DCOS(TH)

        IF(DABS(A-Q(1,1)).LT.TOL.AND.DABS(B-Q(2,1)).LT.TOL.
     &    AND.DABS(C-Q(1,2)).LT.TOL.AND.DABS(D-Q(2,2)).LT.TOL)ICHK=1

        RETURN
        END
C **********************************************************************
        SUBROUTINE EULCHECK1(Q,TH,PHI,OM,ICHK)

        IMPLICIT REAL*8(A-H,O-Z)

        DIMENSION Q(3,3)

        TOL=1.0D-3
        ICHK = 0

        Q(3,3) = 1.0*Q(3,3)/DABS(Q(3,3))
        TH = DACOS(Q(3,3))
        IF(DABS(Q(1,3)).GT.TOL)RETURN
        IF(DABS(Q(2,3)).GT.TOL)RETURN
        IF(DABS(Q(3,1)).GT.TOL)RETURN
        IF(DABS(Q(3,2)).GT.TOL)RETURN
        IF(Q(3,3).EQ.1.0)THEN
           IF(DABS(Q(1,1)-Q(2,2)).GT.TOL)RETURN
           IF(DABS(Q(1,2)+Q(2,1)).GT.TOL)RETURN
        ELSE
           IF(DABS(Q(1,1)+Q(2,2)).GT.TOL)RETURN
           IF(DABS(Q(1,2)-Q(2,1)).GT.TOL)RETURN
        ENDIF
        PHI = DATAN2(Q(1,2),Q(1,1))
        OM = 0.0
        ICHK = 1

        RETURN
        END
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

! **********************************************************************
      SUBROUTINE EMATR(ELMAT,ELMATR)
C      
C  THIS SUBROUTINE REDUCES AN ELASTIC STIFFNESS TYPE MATRIX 
C  3 X 3 X 3 X 3 MATRIX TO A 6 X 6 MATRIX
C          
      IMPLICIT REAL*8(A-H,O-Z) 
      DIMENSION ELMAT(3,3,3,3), ELMATR(6,6)
C      
        DO I = 1,6
        DO J = 1,6
          ELMATR(I,J) = 0.0
        ENDDO
        ENDDO
C 
        ELMATR(1,1)  = ELMAT(1,1,1,1)
        ELMATR(1,2)  = ELMAT(1,1,2,2)
        ELMATR(1,3)  = ELMAT(1,1,3,3) 
        ELMATR(1,4)  = 0.5 * ( ELMAT(1,1,1,2) + ELMAT(1,1,2,1) )
        ELMATR(1,5)  = 0.5 * ( ELMAT(1,1,1,3) + ELMAT(1,1,3,1) )
        ELMATR(1,6)  = 0.5 * ( ELMAT(1,1,2,3) + ELMAT(1,1,3,2) )
C
        ELMATR(2,1)  = ELMAT(2,2,1,1)
        ELMATR(2,2)  = ELMAT(2,2,2,2)
        ELMATR(2,3)  = ELMAT(2,2,3,3) 
        ELMATR(2,4)  = 0.5 * ( ELMAT(2,2,1,2) + ELMAT(2,2,2,1) )
        ELMATR(2,5)  = 0.5 * ( ELMAT(2,2,1,3) + ELMAT(2,2,3,1) )
        ELMATR(2,6)  = 0.5 * ( ELMAT(2,2,2,3) + ELMAT(2,2,3,2) )
C       
        ELMATR(3,1)  = ELMAT(3,3,1,1)
        ELMATR(3,2)  = ELMAT(3,3,2,2)
        ELMATR(3,3)  = ELMAT(3,3,3,3) 
        ELMATR(3,4)  = 0.5 * ( ELMAT(3,3,1,2) + ELMAT(3,3,2,1) )
        ELMATR(3,5)  = 0.5 * ( ELMAT(3,3,1,3) + ELMAT(3,3,3,1) )
        ELMATR(3,6)  = 0.5 * ( ELMAT(3,3,2,3) + ELMAT(3,3,3,2) )   
C
        ELMATR(4,1)  = ELMAT(1,2,1,1)
        ELMATR(4,2)  = ELMAT(1,2,2,2)
        ELMATR(4,3)  = ELMAT(1,2,3,3) 
        ELMATR(4,4)  = 0.5 * ( ELMAT(1,2,1,2) + ELMAT(1,2,2,1) )
        ELMATR(4,5)  = 0.5 * ( ELMAT(1,2,1,3) + ELMAT(1,2,3,1) )
        ELMATR(4,6)  = 0.5 * ( ELMAT(1,2,2,3) + ELMAT(1,2,3,2) )        
C
        ELMATR(5,1)  = ELMAT(1,3,1,1)
        ELMATR(5,2)  = ELMAT(1,3,2,2)
        ELMATR(5,3)  = ELMAT(1,3,3,3) 
        ELMATR(5,4)  = 0.5 * ( ELMAT(1,3,1,2) + ELMAT(1,3,2,1) )
        ELMATR(5,5)  = 0.5 * ( ELMAT(1,3,1,3) + ELMAT(1,3,3,1) )
        ELMATR(5,6)  = 0.5 * ( ELMAT(1,3,2,3) + ELMAT(1,3,3,2) ) 
C
        ELMATR(6,1)  = ELMAT(2,3,1,1)
        ELMATR(6,2)  = ELMAT(2,3,2,2)
        ELMATR(6,3)  = ELMAT(2,3,3,3) 
        ELMATR(6,4)  = 0.5 * ( ELMAT(2,3,1,2) + ELMAT(2,3,2,1) )
        ELMATR(6,5)  = 0.5 * ( ELMAT(2,3,1,3) + ELMAT(2,3,3,1) )
        ELMATR(6,6)  = 0.5 * ( ELMAT(2,3,2,3) + ELMAT(2,3,3,2) ) 
        RETURN
        END
!///////////////////////////////////////////////////////////////////////
!  Subroutine for Calculating the Const. Jacobian
!///////////////////////////////////////////////////////////////////////
       SUBROUTINE calc_jac(nslip,props,nprops,elmat,ft,ftau,fpit,tstartau,
     & gammadot,rss,stau,rhossd,rhomf,q,smatc,theta,del,DTIME,wjmatr)

        IMPLICIT REAL*8(A-H,O-Z)     
        
        DIMENSION  FTAU(3,3)
        DIMENSION  FPIT(3,3)
        DIMENSION  Q(3,3),QT(3,3),SMATG(nslip,3,3)
        DIMENSION  AMAT(3,3)
        DIMENSION  SALPHA(3,3),PALPHA(3,3),PALPHAV(6),PMAT(nslip,6)
        DIMENSION  CALPHA(3,3),CALPHAV(6),CMAT(nslip,6)
        DIMENSION  TSTARTAU(3,3),STAU(nslip)
        DIMENSION  RSS(nslip),DGAMMA(nslip)
        DIMENSION  FPITAU(3,3),FSTAU(3,3)
        DIMENSION  AUX1(3,3),AUX2(3,3)
C            arrays for analytical jacobian
        DIMENSION  FSTAUINV(3,3)
        DIMENSION  FT(3,3),FTTAU(3,3),RTTAU(3,3),UTTAU(3,3)
        DIMENSION  FSTART(3,3),CJMAT(3,3,3,3),ELMAT(3,3,3,3)
        DIMENSION  DJMAT(3,3,3,3),DJMATR(6,6)
        DIMENSION  GJMAT(nslip,3,3,3,3)
        DIMENSION  AJJMAT(nslip,3,3,3,3)
        DIMENSION  AJJMATR(nslip,6,6)
        DIMENSION  BJMAT(nslip,3,3),BJMATR(nslip,6)
        DIMENSION  AKJMATR(6,6),AKJMATRINV(6,6)
        DIMENSION  AMJMATR(6,6)
        DIMENSION  QJMATR(6,6),QJMAT(3,3,3,3)
        DIMENSION  RJMAT(nslip,3,3)
        DIMENSION  SJMAT(3,3,3,3),WJMAT(3,3,3,3),WJMATR(6,6)
        DIMENSION  INDX(6)
        DIMENSION  AUXMAT1(3,3,3,3),AUXMAT1R(6,6),AUXMAT2(3,3,3,3)
        DIMENSION  AUXV1(6)

        dimension  gammadot(nslip),rhossd(nslip),rhomf(nslip)
        dimension  smatc(nslip,3,3),del(3,3)
        dimension  props(nprops)
C
C----------------------------------------------------------------
C             CALCULATE THE JACOBIAN
C----------------------------------------------------------------

!  Intitialization
        dgamma = gammadot*dtime

C  STEP 1. CALCULATE THE MATRIX AMAT
C
        aux1 = matmul(ftau,fpit)
        aux2 = transpose(aux1)
        amat = matmul(aux2,aux1)
        
        DO ISLIP = 1,NSLIP
C
C  GET S MATRIX IN LOCAL COORDINATE SYSTEM
C
        DO I=1,3
        DO J=1,3
          SALPHA(I,J) = SMATC(ISLIP,I,J)
          SMATG(ISLIP,I,J) = 0.0
        ENDDO
        ENDDO
C
C  CONVERT S MATRIX TO GLOBAL COORDINATE SYSTEM
        salpha = matmul(q,matmul(salpha,transpose(q)))
        DO I = 1,3
        DO J = 1,3                 
          SMATG(ISLIP,I,J) =SALPHA(I,J)
        ENDDO
        ENDDO              
C
        DO I = 1,3
        DO J = 1,3
          PALPHA(I,J) = 0.5*(SALPHA(I,J)+SALPHA(J,I))
        ENDDO
        ENDDO
C
        CALL SMATVEC(PALPHA,PALPHAV)
C
C  MULTIPLY THE SHEAR COMPONENTS OF PALPHAV
C  BY A FACTOR OF 2 TO BE CONSISTENT WITH
C  THE DEFINITION OF ENGINEERING SHEAR STRAINS
C
        PALPHAV(4) = 2.0*PALPHAV(4)
        PALPHAV(5) = 2.0*PALPHAV(5)
        PALPHAV(6) = 2.0*PALPHAV(6)
        DO I=1,6
          PMAT(ISLIP,I) = PALPHAV(I)
        ENDDO
C
        aux1 = matmul(amat,salpha)
        DO I = 1,3
        DO J = 1,3
          AUX2(I,J) = 0.5*(AUX1(I,J) + AUX1(J,I))
        ENDDO
        ENDDO

        calpha = reshape( (/((sum(elmat(i,j,:,:)*aux2), i = 1,3), j=1,3)/), (/3,3/))

        CALL SMATVEC(CALPHA,CALPHAV)
        DO I=1,6
          CMAT(ISLIP,I) = CALPHAV(I)
        ENDDO

        ENDDO !ISLIP

C  STEP 5. CALCULATE FPI AT TME TAU
C          
        DO I =1,3
        DO J =1,3
          AUX1(I,J)= 0.0
          DO ISLIP=1,NSLIP
            AUX1(I,J)= AUX1(I,J) + 
     &                 DGAMMA(ISLIP)*SMATG(ISLIP,I,J)                   
          ENDDO
          AUX2(I,J) = del(I,J) - AUX1(I,J)
        ENDDO
        ENDDO    

        fpitau = matmul(fpit,aux2)

        det = determinant(fpitau)
        CDET = DET**(1./3.)
        DO I =1,3
        DO J=1,3
          FPITAU(I,J) = FPITAU(I,J)/CDET
        ENDDO
        ENDDO
C
C  STEP 6. CALCULATE FSTAU 
        fstau = matmul(ftau,fpitau)

!////////////////////////////////////////////////////////////////
C
C  STEP J1. CALCULATE THE C_IJKL MATRIX
C
        CALL inverse_3x3(FT,AUX1)           

        fttau = matmul(ftau,aux1)

        CALL polardecomp(FTTAU,RTTAU,UTTAU)

        fstart = matmul(ft,fpit)
        DO I = 1,3
        DO J = 1,3
        DO K = 1,3
        DO L = 1,3
          SUM1 = 0.0
          SUM2 = 0.0
          DO M = 1,3
            SUM1 = SUM1 + UTTAU(L,M)  * FSTART(M,J)
            SUM2 = SUM2 + FSTART(M,I) * UTTAU(M,K)
          ENDDO       
          CJMAT(I,J,K,L) = FSTART(K,I) * SUM1 +
     &                     SUM2 * FSTART(L,J)
        ENDDO
        ENDDO
        ENDDO
        ENDDO       
   
C  STEP J2. CALCULATE THE L_IJKL MATRIX IN GLOBAL
C           CORDINATES

!  ELMAT is now passed in through the calc_jac subroutine
!  interface
        
             
C  STEP J3. CALCULATE THE D_IJKL MATRIX

        DO I = 1,3
        DO J = 1,3
        DO K = 1,3
        DO L = 1,3
          DJMAT(I,J,K,L) = 0.0
          DO M  = 1,3
          DO N  = 1,3
            DJMAT(I,J,K,L) = DJMAT(I,J,K,L) +
     &                       0.5*ELMAT(I,J,M,N)*CJMAT(M,N,K,L)
          ENDDO
          ENDDO
        ENDDO
        ENDDO
        ENDDO
        ENDDO
C
        CALL EMATR(DJMAT,DJMATR)     
C
C  STEP J4. CALCULATE THE NLIP MATRICES G_MNKL^{\ALPHA}
C
        DO ISLIP = 1,NSLIP
          DO M = 1,3
          DO N = 1,3
          DO K = 1,3
          DO L = 1,3
            GJMAT(ISLIP,M,N,K,L) = 0.0
            DO IP = 1,3
              GJMAT(ISLIP,M,N,K,L) = GJMAT(ISLIP,M,N,K,L) +
     &                     CJMAT(M,IP,K,L)*SMATG(ISLIP,IP,N) +
     &                     SMATG(ISLIP,IP,M)*CJMAT(IP,N,K,L)
            ENDDO
          ENDDO
          ENDDO
          ENDDO
          ENDDO
        ENDDO
C  
C  STEP J5. CALCULATE THE NSLIP MATRICES J_IJKL^{\ALPHA}
C            
        DO ISLIP = 1,NSLIP
          DO I = 1,3
          DO J = 1,3
          DO K = 1,3
          DO L = 1,3
            AJJMAT(ISLIP,I,J,K,L) = 0.0
            DO M  = 1,3
            DO N  = 1,3
              AJJMAT(ISLIP,I,J,K,L) = AJJMAT(ISLIP,I,J,K,L) +
     &                  0.5*ELMAT(I,J,M,N)*GJMAT(ISLIP,M,N,K,L)
            ENDDO
            ENDDO
            AUXMAT1(I,J,K,L) = AJJMAT(ISLIP,I,J,K,L)
          ENDDO
          ENDDO
          ENDDO
          ENDDO
C
          CALL EMATR(AUXMAT1,AUXMAT1R)
          DO I = 1,6
          DO J = 1,6
            AJJMATR(ISLIP,I,J) = AUXMAT1R(I,J)
          ENDDO
          ENDDO            
        ENDDO
C
C  STEP J6. CALCULATE THE NSLIP MATRICES B_IJ^{\ALPHA}

        DO ISLIP = 1,NSLIP

          call calc_GAMMADOT(props,nprops,theta,rhossd(islip),
     &      rss(islip),stau(islip),rhomf(islip),gdotalpha,dgdtau)

          DO I     = 1,3
          DO J     = 1,3
            BJMAT(ISLIP,I,J) = 
     &        0.5*(SMATG(ISLIP,I,J)+SMATG(ISLIP,J,I))*
     &        DGDTAU*DTIME
            AUX1(I,J)=BJMAT(ISLIP,I,J)
          ENDDO
          ENDDO
C
          CALL SMATVEC(AUX1,AUXV1)
C
          DO I=1,6
            BJMATR(ISLIP,I) = AUXV1(I)
          ENDDO
        ENDDO
C
C  STEP J7. CALCULATE THE  REDUCED MATRIX K_IJ
C
        DO I = 1,6
        DO J = 1,6
          AKJMATR(I,J) = 0.0
          DO ISLIP = 1,NSLIP
            AKJMATR(I,J) = AKJMATR(I,J) +
     &                CMAT(ISLIP,I)*BJMATR(ISLIP,J)
          ENDDO
          IF(I.EQ.J)THEN
            AKJMATR(I,J) = 1. + AKJMATR(I,J)
          ENDIF
        ENDDO
        ENDDO
C
C  STEP J8. CALCULATE THE INVERSE OF THE MATRIX K_IJ
C
        CALL MATINV(AKJMATR,6,6,INDX,AKJMATRINV)
C
C  STEP J9. CALCULATE THE REDUCED MATRIX M_IJ
C
        DO I = 1,6
        DO J = 1,6
          AMJMATR(I,J) = 0.0
          DO ISLIP= 1, NSLIP
            AMJMATR(I,J) = AMJMATR(I,J) +
     &                     DGAMMA(ISLIP)*AJJMATR(ISLIP,I,J)
          ENDDO
          AMJMATR(I,J) =  DJMATR(I,J) - AMJMATR(I,J)
        ENDDO
        ENDDO
C
C  STEP J10. CALCULATE THE  MATRIX Q_IJKL
C     
        DO I = 1,6
        DO J = 1,6
          QJMATR(I,J) = 0.0
          DO M  = 1,6
            QJMATR(I,J) = QJMATR(I,J)+ AKJMATRINV(I,M)*
     &                    AMJMATR(M,J)
          ENDDO
        ENDDO
        ENDDO
C
          CALL EMATREC(QJMATR,QJMAT)
C
C  STEP J11. CALCULATE THE MATRICES R_IJ^{\ALPHA}
C                 
        DO  ISLIP = 1, NSLIP
        DO  I = 1,3
        DO  J = 1,3
          RJMAT(ISLIP,I,J) = 0.0
          DO K = 1,3
          DO L = 1,3
            RJMAT(ISLIP,I,J) = RJMAT(ISLIP,I,J) +
     &                         BJMAT(ISLIP,K,L)*QJMAT(K,L,I,J)
          ENDDO
          ENDDO
        ENDDO
        ENDDO
        ENDDO
C
C  STEP J12. CALCULATE THE MATRIX S_IJKL
C 
        DO IP = 1,3
        DO J = 1,3
        DO K = 1,3
        DO L = 1,3
          AUX1(IP,J)        = 0.0
          AUXMAT1(IP,J,K,L) = 0.0
          DO ISLIP = 1,NSLIP
            AUX1(IP,J) = AUX1(IP,J) 
     &                  + DGAMMA(ISLIP)*SMATG(ISLIP,IP,J)
            AUXMAT1(IP,J,K,L) = AUXMAT1(IP,J,K,L)
     &                  +SMATG(ISLIP,IP,J)*RJMAT(ISLIP,K,L)
          ENDDO           
        ENDDO
        ENDDO
        ENDDO
        ENDDO
C                   
        DO L = 1,3
        DO J = 1,3
          AUX2(L,J) = 0.0
          DO IP = 1,3
            AUX2(L,J) = AUX2(L,J) +
     &                 FSTART(L,IP)*AUX1(IP,J)        
          ENDDO               
        ENDDO
        ENDDO               
C
        DO N = 1,3
        DO J = 1,3
        DO K = 1,3
        DO L = 1,3
          AUXMAT2(N,J,K,L) = 0.0
          DO IP = 1,3                
            AUXMAT2(N,J,K,L) =  AUXMAT2(N,J,K,L) +
     &                    FSTART(N,IP)*AUXMAT1(IP,J,K,L)
          ENDDO
        ENDDO
        ENDDO
        ENDDO
        ENDDO   
C
        DO I = 1,3
        DO J = 1,3
        DO K = 1,3
        DO L = 1,3
          SJMAT(I,J,K,L) = 0.0
          AUXMAT1(I,J,K,L) = 0.0
          DO N = 1,3                 
            AUXMAT1(I,J,K,L) =  AUXMAT1(I,J,K,L) +
     &                      FTTAU(I,N)*AUXMAT2(N,J,K,L)
          ENDDO
          SJMAT(I,J,K,L) = SJMAT(I,J,K,L) 
     &           + RTTAU(I,K)*FSTART(L,J)
     &           - RTTAU(I,K)*AUX2(L,J)
     &           - AUXMAT1(I,J,K,L)
        ENDDO
        ENDDO
        ENDDO
        ENDDO   
C
C  STEP J13. CALCULATE THE JACOBIAN MATRIX FOR EACH CRYSTAL
C
        CALL inverse_3x3(FSTAU,FSTAUINV)
C
        detfstau = determinant(fstau)
C            
        DO I = 1,3
        DO J = 1,3
        DO K = 1,3
        DO L = 1,3
          WJMAT(I,J,K,L) = 0.0
          SUM1 = 0.0
          DO IP = 1,3
          DO IQ = 1,3
            SUM1 = SUM1 + SJMAT(IP,IQ,K,L)*FSTAUINV(IQ,IP)
          ENDDO
          ENDDO
          DO M = 1,3
          DO N = 1,3
            WJMAT(I,J,K,L) = WJMAT(I,J,K,L) +
     &                SJMAT(I,M,K,L)*TSTARTAU(M,N)*FSTAU(J,N) +
     &                FSTAU(I,M)*QJMAT(M,N,K,L)*FSTAU(J,N) +
     &                FSTAU(I,M)*TSTARTAU(M,N)*SJMAT(J,N,K,L) -
     &                FSTAU(I,M)*TSTARTAU(M,N)*FSTAU(J,N)*SUM1
          ENDDO
          ENDDO
          WJMAT(I,J,K,L) = WJMAT(I,J,K,L)/DETFSTAU
        ENDDO
        ENDDO
        ENDDO
        ENDDO
C
C  STEP J14. CALCULATE THE REDUCED JACOBIAN MATRIX FOR EACH CRYSTAL
C
        CALL EMATR(WJMAT,WJMATR)
C
        RETURN
        END 
C*********************************************************************
C ****************************************************************
C  THIS SUBROUTINE CONVERTS A SYMMETRIC MATRIX TO A VECTOR
        SUBROUTINE SMATVEC(A,V)
        IMPLICIT REAL*8(A-H,O-Z)       
        DIMENSION A(3,3),V(6)
        DO I = 1,3
          V(I) = A(I,I)
        ENDDO
        V(4) = A(1,2)
        V(5) = A(1,3)
        V(6) = A(2,3)
        RETURN
        END
C ******************************************************************
C ****************************************************************
      SUBROUTINE MATINV(A,N,NP,INDX,Y)
C
C  Given an NxN matrix [A], with physical dimension NP, this 
C  routine replaces it by the LU decomposition of a rowwise 
C  permutation of itself. [A] and N are input. [A] is output, 
C  arranged in LU form. INDX is an output vector which records
C  the row permutation effected by the partial pivoting; 
C  D is output as +1 or -1 depending on wheter the number of
C  row interchanges was even or odd, respectively.
C
C  Once the LU decomposition is performed, this routine 
C  calculates the inverse of [A] by using subroutine LUBKSB.
C  Note that INDX is input as the permutation vector returned by        
C  LUDCMP. {B} is input as the right-hand side vector {B}, and  
C  returns with the solution vector {X}. [A], N, NP, INDX are not 
C  modified by this routine, and are left in place
C  for succesive calls with different right-hand sides {B}.
C  This routine takes into account that {B} will begin with
C  many zero elements, so it is efficient for use in matrix 
C  inversion.
C
C  The inverse of [A] is calculated using as many unit vectors
C  {B} needed as right hand side vectors. The result is
C  returned as the matrix [Y].
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NP,NP), Y(NP,NP), INDX(NP)
C
C  Set up the identity matrix  
C
        DO I = 1,N
          DO J = 1,N
            Y(I,J) = 0.0
          ENDDO
          Y(I,I) = 1.
        ENDDO
C
C  Decompose the matrix just once
C
        CALL LUDCMP(A,N,NP,INDX,D)
C
C  Find the inverse by columns. It is necessary to recognize
C  that FORTRAN stores two dimensional matrices by column, so
C  so that Y(1,J) is the address of the Jth column of Y.
C
        DO J=1,N
          CALL LUBKSB(A,N,NP,INDX,Y(1,J))
        ENDDO
        RETURN
        END
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
C ****************************************************************
C******************************************************************************
      SUBROUTINE EMATREC(ELMATR,ELMAT)
C      
C  THIS SUBROUTINE RECONSTRUCTS A STIFFNESS MATRIX ELMATREC
C  AS A 3 X 3 X 3 X 3 MATRIX FROM ITS 
C  REDUCED 6 X 6 FORM ELMATR
C
      IMPLICIT REAL*8(A-H,O-Z) 
      DIMENSION ELMATR(6,6),ELMAT(3,3,3,3)
C
      DO I = 1,3
      DO J = 1,3
      DO K = 1,3
      DO L = 1,3
        ELMAT(I,J,K,L) = 0.0
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
C  RECONSTRUCT THE MATRIX ELMAT FROM ITS REDUCED FORM ELMATR
C
        ELMAT(1,1,1,1)  = ELMATR(1,1)
        ELMAT(1,1,2,2)  = ELMATR(1,2)
        ELMAT(1,1,3,3)  = ELMATR(1,3) 
        ELMAT(1,1,1,2)  = ELMATR(1,4)
        ELMAT(1,1,2,1)  = ELMATR(1,4)
        ELMAT(1,1,1,3)  = ELMATR(1,5)
        ELMAT(1,1,3,1)  = ELMATR(1,5)
        ELMAT(1,1,2,3)  = ELMATR(1,6)
        ELMAT(1,1,3,2)  = ELMATR(1,6)
C
        ELMAT(2,2,1,1)  = ELMATR(2,1)
        ELMAT(2,2,2,2)  = ELMATR(2,2)
        ELMAT(2,2,3,3)  = ELMATR(2,3) 
        ELMAT(2,2,1,2)  = ELMATR(2,4)
        ELMAT(2,2,2,1)  = ELMATR(2,4)
        ELMAT(2,2,1,3)  = ELMATR(2,5)
        ELMAT(2,2,3,1)  = ELMATR(2,5)
        ELMAT(2,2,2,3)  = ELMATR(2,6)
        ELMAT(2,2,3,2)  = ELMATR(2,6)
C
        ELMAT(3,3,1,1)  = ELMATR(3,1)
        ELMAT(3,3,2,2)  = ELMATR(3,2)
        ELMAT(3,3,3,3)  = ELMATR(3,3) 
        ELMAT(3,3,1,2)  = ELMATR(3,4)
        ELMAT(3,3,2,1)  = ELMATR(3,4)
        ELMAT(3,3,1,3)  = ELMATR(3,5)
        ELMAT(3,3,3,1)  = ELMATR(3,5)
        ELMAT(3,3,2,3)  = ELMATR(3,6)
        ELMAT(3,3,3,2)  = ELMATR(3,6)
C
        ELMAT(1,2,1,1)  =  ELMATR(4,1)
        ELMAT(2,1,1,1)  =  ELMATR(4,1)
        ELMAT(1,2,2,2)  =  ELMATR(4,2)
        ELMAT(2,1,2,2)  =  ELMATR(4,2)  
        ELMAT(1,2,3,3)  =  ELMATR(4,3)
        ELMAT(2,1,3,3)  =  ELMATR(4,3)
        ELMAT(1,2,1,2)  =  ELMATR(4,4) 
        ELMAT(2,1,1,2)  =  ELMATR(4,4) 
        ELMAT(1,2,2,1)  =  ELMATR(4,4)
        ELMAT(2,1,2,1)  =  ELMATR(4,4)
        ELMAT(1,2,1,3)  =  ELMATR(4,5) 
        ELMAT(2,1,1,3)  =  ELMATR(4,5) 
        ELMAT(1,2,3,1)  =  ELMATR(4,5)
        ELMAT(2,1,3,1)  =  ELMATR(4,5)
        ELMAT(1,2,2,3)  =  ELMATR(4,6) 
        ELMAT(2,1,2,3)  =  ELMATR(4,6) 
        ELMAT(1,2,3,2)  =  ELMATR(4,6)
        ELMAT(2,1,3,2)  =  ELMATR(4,6)
C
        ELMAT(1,3,1,1)  =  ELMATR(5,1)
        ELMAT(3,1,1,1)  =  ELMATR(5,1)  
        ELMAT(1,3,2,2)  =  ELMATR(5,2)
        ELMAT(3,1,2,2)  =  ELMATR(5,2)  
        ELMAT(1,3,3,3)  =  ELMATR(5,3)
        ELMAT(3,1,3,3)  =  ELMATR(5,3)
        ELMAT(1,3,1,2)  =  ELMATR(5,4) 
        ELMAT(3,1,1,2)  =  ELMATR(5,4) 
        ELMAT(1,3,2,1)  =  ELMATR(5,4)
        ELMAT(3,1,2,1)  =  ELMATR(5,4)
        ELMAT(1,3,1,3)  =  ELMATR(5,5) 
        ELMAT(3,1,1,3)  =  ELMATR(5,5) 
        ELMAT(1,3,3,1)  =  ELMATR(5,5)
        ELMAT(3,1,3,1)  =  ELMATR(5,5)
        ELMAT(1,3,2,3)  =  ELMATR(5,6) 
        ELMAT(3,1,2,3)  =  ELMATR(5,6) 
        ELMAT(1,3,3,2)  =  ELMATR(5,6)
        ELMAT(3,1,3,2)  =  ELMATR(5,6)
C
        ELMAT(2,3,1,1)  =  ELMATR(6,1)
        ELMAT(3,2,1,1)  =  ELMATR(6,1)  
        ELMAT(2,3,2,2)  =  ELMATR(6,2)
        ELMAT(3,2,2,2)  =  ELMATR(6,2)  
        ELMAT(2,3,3,3)  =  ELMATR(6,3)
        ELMAT(3,2,3,3)  =  ELMATR(6,3)
        ELMAT(2,3,1,2)  =  ELMATR(6,4) 
        ELMAT(3,2,1,2)  =  ELMATR(6,4) 
        ELMAT(2,3,2,1)  =  ELMATR(6,4)
        ELMAT(3,2,2,1)  =  ELMATR(6,4)
        ELMAT(2,3,1,3)  =  ELMATR(6,5) 
        ELMAT(3,2,1,3)  =  ELMATR(6,5) 
        ELMAT(2,3,3,1)  =  ELMATR(6,5)
        ELMAT(3,2,3,1)  =  ELMATR(6,5)
        ELMAT(2,3,2,3)  =  ELMATR(6,6) 
        ELMAT(3,2,2,3)  =  ELMATR(6,6) 
        ELMAT(2,3,3,2)  =  ELMATR(6,6)
        ELMAT(3,2,3,2)  =  ELMATR(6,6)
        RETURN 
        END  
C **************************************************************
C*************************************************************************
        SUBROUTINE SOLVCRYS(nslip,props,nprops,theta,dfgrd0,dfgrd1,TSTV,ST,amat,dmat,
     &  rhossd_n,rhomf_n,chi_n,TSTRV,PMAT,CMAT,PMAT0,DTIME,TSTAUV,STAU,
     &  rhossd,drhossd,rhomf,xnu,DGAMMA,DGT,DGMAX,RSS,DOMPSC,ITERK,ITERL,ITERERR,xmC)
C
C       THIS SUBROUTINE SOLVES FOR  (TSTARTAU,STAU)
C       HERE, WE ARE SETTING
C
C       TSTARTV   AS TSTV
C       TSTARTRV  AS TSTRV
C       TSTARTAUV AS TSTAUV
C --------------------------------------------------------------------
        IMPLICIT REAL*8(A-H,O-Z)

        PARAMETER  (MAXIT1=100,MAXIT2=20,NP=6,N=6)
        DIMENSION  ALPHAM(NP,NP),BETA(NP),INDX(NP)
        DIMENSION  TSTV(6),ST(nslip),TSTRV(6)
        DIMENSION  dfgrd0(3,3),dfgrd1(3,3)
        DIMENSION  PMAT0(3,3,nslip)
        DIMENSION  PMAT(nslip,6)
        DIMENSION  CMAT(nslip,6)
        DIMENSION  TSTAUV(6),STAU(nslip),STRESS(6)
        DIMENSION  RSS(nslip),DGAMMA(nslip),DGAMDTAU(nslip)
        DIMENSION  RESOLD(nslip),RESNEW(nslip)
        DIMENSION  TEMPSTRESS(6)
        dimension  props(nprops)
        dimension  rhossd_n(nslip),rhossd(nslip),rhossd_old(nslip)
        dimension  drhossd(nslip),ds(nslip),rhomf_n(nslip),rhomf(nslip)
        dimension  amat(nslip,nslip), dmat(nslip,nslip)
        dimension  xmC(3,nslip)
        dimension  xnu(nslip),DGT(nslip)

!  Initialize iteration counters and other parameters
        ITERK    = 0
        ITERL    = 0
        ITERLNEW = 0
        TOL      = 1.d-4
        DGMAXTOL = 0.50
        CORRBETA = 4.0
C
C  INITIAL ESTIMATES FOR THE ROOTS:
C
        DO I=1,6
           STRESS(I) = TSTV(I)
        ENDDO

        DO ISLIP=1,NSLIP
           RESOLD(ISLIP) = ST(ISLIP)
           rhossd_old(islip) = rhossd_n(islip)
        ENDDO
999     CONTINUE
        ITERERR = 0
C
C  START THE LOOP FOR SOLVING FOR THE COMPONENTS OF TSTAUV
C
        DO K=1,MAXIT1

!        if (K>5) then
!          dgmax = 0.d0
!        endif
C
C  CALCULATE:
C     1. THE RESOLVED SHEAR STRESS ON EACH
C        SLIP SYSTEM CORRESPONDING TO THE ESTIMATE OF TSTAUV.
C     2. CALCULATE THE PLASTIC SHEAR STRAIN INCREMENT FOR
C        EACH SLIP SYSTEM.
C     3. CALCULATE THE MAXIMUM SHEAR STRAIN INCREMENT 
C     4. CALCULATE THE PLASTIC WORK INCREMENT
C
           DGMAX   = 0.0
           DOMPSC = 0.0
           DO ISLIP = 1,NSLIP
             TAUALPHA = 0.0
             DO J=1,6
               TAUALPHA = TAUALPHA + STRESS(J)*PMAT(ISLIP,J)
             ENDDO
             RSS(ISLIP) = TAUALPHA
             RESALPHA   = RESOLD(ISLIP)
             RHOALPHA   = rhossd_old(ISLIP)
             RMFALPHA   = rhomf_n(ISLIP)

             call calc_GAMMADOT(props,nprops,theta,RHOALPHA,
     &         TAUALPHA,RESALPHA,RMFALPHA,GDOTALPHA,DGDTAU)

             DGAMMA(ISLIP) = DTIME*GDOTALPHA
             DGAMDTAU(ISLIP) = DTIME*DGDTAU
             DGMAX = DMAX1( DGMAX,DABS( DGAMMA(ISLIP) ) )
             TEMP =  DABS ( RSS(ISLIP)*DGAMMA(ISLIP) ) ! work increment
             DOMPSC = DOMPSC + TEMP
           ENDDO

           if (dgmax < 1d-20) then
             tstauv = tstrv
             stau = st
             rhossd = rhossd_n
             rhomf = rhomf_n
             return
           endif
             


C
C  Calculate BETA_i --- the negative of the values of
C  the functions G_i(T_k)
C
           DO I =1,6
             BETA(I) = 0.0
             DO ISLIP=1,NSLIP
               BETA(I) = BETA(I) + DGAMMA(ISLIP)*CMAT(ISLIP,I)
             ENDDO
             BETA(I) =   TSTRV(I) - BETA(I) - STRESS(I)
           ENDDO
C
C  Calculate the Jacobian for the Newton corrections
C
           DO I =1,6
           DO J =1,6
             IF (I .EQ. J) THEN      
               ALPHAM(I,J)      = 1.0D0
             ELSE
               ALPHAM(I,J)      = 0.0D0
             END IF
             DO ISLIP=1,NSLIP
               IF(RSS(ISLIP).EQ.0.0)GO TO 60
               ALPHAM(I,J)= ALPHAM(I,J) + 
     &                    DGAMDTAU(ISLIP)*CMAT(ISLIP,I)*PMAT(ISLIP,J)
60           CONTINUE
             ENDDO
           ENDDO
           ENDDO    
C
C  Calculate the Newton corrections
C

          CALL LUDCMP(ALPHAM,N,NP,INDX,D)
          CALL LUBKSB(ALPHAM,N,NP,INDX,BETA)
C
C  BETA  NOW CONTAINS THE NEWTON CORRECTIONS.
C  CHECK TO SEE IF THE NEWTON CORRECTIONS ARE
C  WITHIN SPECIFIED TOLERANCES. IF YES, THEN
C  PROCESS HAS CONVERGED. UPDATE AND GO TO LOOP FOR
C  ITERATING ON THE DEFORMATION RESISTANCES.
C  IF NOT, THEN CONTINUE ITERATING.
C
          IF (   ( DABS( BETA(1) ) .LE. TOL ) .AND.
     +           ( DABS( BETA(2) ) .LE. TOL ) .AND.
     +           ( DABS( BETA(3) ) .LE. TOL ) .AND.
     +           ( DABS( BETA(4) ) .LE. TOL ) .AND.
     +           ( DABS( BETA(5) ) .LE. TOL ) .AND.
     +           ( DABS( BETA(6) ) .LE. TOL)   )  THEN  
C
C  ITERATIVE PROCEDURE HAS CONVERGED. UPDATE 
C
              TSTAUV(1) = STRESS(1)
              TSTAUV(2) = STRESS(2)
              TSTAUV(3) = STRESS(3)
              TSTAUV(4) = STRESS(4)
              TSTAUV(5) = STRESS(5)
              TSTAUV(6) = STRESS(6)           
              ITERK    =  K
              GO TO 550
           ENDIF
C
C  TEMPORARILY ADD THE CORRECTIONS TO THE VARIABLES
C
           TEMPSTRESS(1) = STRESS(1) + BETA(1)
           TEMPSTRESS(2) = STRESS(2) + BETA(2)
           TEMPSTRESS(3) = STRESS(3) + BETA(3)
           TEMPSTRESS(4) = STRESS(4) + BETA(4)
           TEMPSTRESS(5) = STRESS(5) + BETA(5)
           TEMPSTRESS(6) = STRESS(6) + BETA(6)
C
C  CHECK IF THESE STRESS COMPONENTS GIVE RISE TO LARGE DELTA GAMMAS
C  IF THEY DO THEN REDUCE THE N-R CORRECTIONS BY A FACTOR OF 4
C  AND CONTINUE UNTIL THE DELTA GAMMAS ARE REASONABLE
C
           ICOR = 0
           dgmaxold = 0.d0
85         CONTINUE
           IF(ICOR.GT.10) THEN 
C
C  THE METHOD TO CONSTRAIN THE NEWTON-RAPHSON METHOD 
C  IS NOT WORKING! GET OUT OF THE LOOP! 
C  RETURN TO CUT BACK TIME-STEP AND RECALCULATE
C
!             WRITE(*,*) 'TOO MANY ITERATIONS IN CONSTRAINT....'
             GOTO 500
           ENDIF
           DGMAX = 0.0
           DO ISLIP = 1,NSLIP
             TAUALPHA = 0.0
             DO J=1,6
               TAUALPHA = TAUALPHA + TEMPSTRESS(J)*PMAT(ISLIP,J)
             ENDDO
             RSS(ISLIP) = TAUALPHA
             RESALPHA   = RESOLD(ISLIP)
             RHOALPHA   = rhossd_old(ISLIP)
             RMFALPHA   = rhomf_n(ISLIP)

             call calc_GAMMADOT(props,nprops,theta,RHOALPHA,
     &  TAUALPHA,RESALPHA,RMFALPHA,GDOTALPHA,DGDTAU)

             DGAMMA(ISLIP) = DTIME*GDOTALPHA
             DGAMDTAU(ISLIP) = DTIME*DGDTAU
             DGMAX = DMAX1( DGMAX,DABS( DGAMMA(ISLIP) ) )
           ENDDO
C
           IF(DGMAX.GT.DGMAXOLD)THEN
             CORRBETA=1./CORRBETA
           ENDIF
           DGMAXOLD=DGMAX
C
           IF (DGMAX .GT. DGMAXTOL) THEN
              DO I = 1,6
                 BETA(I) = CORRBETA * BETA(I)
                 TEMPSTRESS(I) = STRESS(I) + BETA(I)
              ENDDO
              ICOR = ICOR +1
              GOTO 85
           ELSE
C
C  ACCEPT THE CORRECTIONS TO THE STRESS AND CONTINUE ITERATING
C
              DO I = 1,6
                 STRESS(I) = TEMPSTRESS(I) 
              ENDDO
           ENDIF
C                    
500     CONTINUE
        ENDDO
C
C  ITERATIVE PROCEDURE FOR THE STRESS 
C  HAS NOT CONVERGED IN MAXIT1 ITERATIONS.
C  RETURN TO CUT BACK TIME-STEP, AND RECALCULATE.
C
        ITERK = MAXIT1
        ITERERR = 1
        RETURN
550   CONTINUE
C
C  ITERATIVE PROCEDURE FOR TSTAUV HAS CONVERGED
C  START ITERATING FOR THE SLIP SYSTEM DEFORMATION 
C  RESISTANCES
C
C  ESTIMATE THE THE SLIP SYSTEM DEFORMATION RESISTANCES
C  AT THE END OF THE STEP
         call calc_HARDENING(props,nprops,nslip,theta,chi_n,amat,dmat,
     &        ST,RSS,rhossd_old,xnu,DGAMMA,DGT,DTIME,ds,drhossd)

C  ESTIMATE THE SLIP SYSTEM DEFORMATION RESISTANCES AT THE
C  END OF THE INCREMENT
C             
        DO ISLIP = 1,NSLIP
          RESNEW(ISLIP) =   ds(ISLIP) + ST(ISLIP)
          rhossd(ISLIP) =   drhossd(ISLIP) + rhossd_n(ISLIP)
        ENDDO
        rhomf = matmul(dmat,rhossd)
        ERRRES = 0.0
        DO ISLIP = 1,NSLIP
         ERRRES = DMAX1(DABS(RESNEW(ISLIP)-RESOLD(ISLIP)),ERRRES)
        ENDDO
C
        TOLS = 1.D-4
        IF (ERRRES .LT. TOLS) THEN
C
C  ITERATIVE PROCEDURE HAS CONVERGED
C  UPDATE STAU AND RETURN
C
           DO ISLIP = 1,NSLIP
             STAU(ISLIP) = RESNEW(ISLIP)
           ENDDO
           RETURN
        ELSE
          DO ISLIP = 1,NSLIP
            RESOLD(ISLIP) = RESNEW(ISLIP)
            rhossd_old(islip) = rhossd(islip)
          ENDDO
          rhomf = matmul(dmat,rhossd)
          ITERLNEW = ITERLNEW + 1
          ITERL = MAX0(ITERL ,ITERLNEW)
          IF (ITERL .GT. MAXIT2) THEN
C
C  ITERATIVE PROCEDURE FOR THE DEFORMATION RESISTANCE HAS NOT
C  CONVERGED IN MAXIT2 ITERATIONS. RETURN TO CUT BACK THE TIME
C  STEP, AND RECALCULATE.
C
            ITERL = MAXIT2
            ITERERR = 1
!            WRITE (*,*) 'MAXIT2 EXCEEDED....REPEATING TIME STEP'
            RETURN
          ENDIF
          GOTO 999
        ENDIF
        END
C **********************************************************************
        SUBROUTINE calc_GAMMADOT(props,nprops,theta,RHOALPHA,
     &  TAUALPHA,RESALPHA,RMFALPHA,GDOTALPHA,DGDTAU)

!  RESALPHA is the CRSS
!  TAUALPHA is the resolved shear stress

        IMPLICIT REAL*8(A-H,O-Z)

        dimension props(nprops)

!  Initialize props needed for the calcuation
        ! C11_0=PROPS(7)
        ! C12_0=PROPS(9)
        ! C44_0=PROPS(11)

        ! call compute_elastConst(props,nprops,theta,c11,c12,c44)

        ! xmu = SQRT((c11-c12)/2.d0*c44)
        ! xmu0 = SQRT((c11_0-c12_0)/2.d0*c44_0)

        burg   = props(13)
        t0     = props(25)
        THETAP = props(26)
        bk     = props(27)
        eta    = props(28)
        rx     = props(37)

        IF (TAUALPHA.eq.0.d0) THEN
           gfunc = 0.d0
           GDOTALPHA = 0.d0
           DGDTAU = 0.d0
        ELSEIF (RHOALPHA.eq.0.d0) THEN
           gfunc = 0.d0
           GDOTALPHA = 0.d0
           DGDTAU = 0.d0
        ELSE
    !       gfunc = 2.d0*burg*SQRT(RHOALPHA)*
    ! &             DEXP(-(THETAP/THETA)*DCOSH(TAUALPHA/RESALPHA))*
    ! &             DSINH((THETAP/THETA)*DSINH(TAUALPHA/RESALPHA))
           gfunc = burg*(RHOALPHA/SQRT(RMFALPHA))*
     &             (DEXP(-(THETAP/THETA)*DEXP(-TAUALPHA/RESALPHA)) -
     &             DEXP(-(THETAP/THETA)*DEXP(TAUALPHA/RESALPHA)))
           GDOTALPHA = gfunc*rx/(t0*(1.d0 + eta*gfunc/(TAUALPHA*(burg**2)*RHOALPHA)))
    !       DgfuncDTAU = 2.d0*burg*SQRT(RHOALPHA)*
    ! &                  DEXP(-(THETAP/THETA)*DCOSH(TAUALPHA/RESALPHA))*
    ! &                  (THETAP/(THETA*RESALPHA))*
    ! &               DCOSH((THETAP/THETA)*DSINH(TAUALPHA/RESALPHA)-TAUALPHA/RESALPHA)
           DgfuncDTAU = burg*(RHOALPHA/SQRT(RMFALPHA))*(THETAP/(THETA*RESALPHA))*
     &                  (DEXP(-(THETAP/THETA)*DEXP(-TAUALPHA/RESALPHA))*
     &                  DEXP(-TAUALPHA/RESALPHA) +
     &                  DEXP(-(THETAP/THETA)*DEXP(TAUALPHA/RESALPHA))*
     &                  DEXP(TAUALPHA/RESALPHA))
           DGDTAU = rx*(DgfuncDTAU + eta*((gfunc/(TAUALPHA*burg))**2)/RHOALPHA)/
     &               (t0*((1.d0 + eta*gfunc/(TAUALPHA*(burg**2)*RHOALPHA))**2))
        ENDIF

    !    WRITE(*,*) gfunc, DGDTAU

        RETURN
        END

C **********************************************************************
! Subroutine for defining the dislocation interatction matrix that
! appears in the Taylor expression for the 
c **********************************************************************
        subroutine calc_amat(nslip,props,nprops,s0,n0,amat)

        implicit none

        real*8 :: amat(nslip,nslip),props(nprops),s0(3,nslip),
     &           n0(3,nslip),t0(3,nslip)
        integer :: i,j,nslip,nprops,ii,jj,istart,jstart,igliss1,igliss2

        real*8 :: sa_dot_sb,na_dot_nb,na_dot_tb,dp1,dp2,
     &            a_self,a_copl,a_hirth,a_colli,a_gliss,a_lomer


        a_self = props(17)
        a_copl = props(18)
        a_hirth = props(19)
        a_colli = props(20)
        a_gliss = props(21)
        a_lomer = props(22)

        amat = 0.d0
        do i = 1,nslip
          do j = 1,nslip
            sa_dot_sb = abs(dot_product(s0(:,i),s0(:,j)))
            na_dot_nb = abs(dot_product(n0(:,i),n0(:,j)))
            ! Self and/or coplanar interaction
            if (na_dot_nb > 0.999d0) then
              if (i==j) then
                amat(i,j) = a_self !0.06d0
              else 
                amat(i,j) = a_copl !0.1d0
              endif
            ! Hirth Lock
            else if ((na_dot_nb < 0.999d0) .and. 
     &               (sa_dot_sb < 0.001d0)) then
              amat(i,j) = a_hirth !0.07d0 
            ! Collinear interaction
            else if ((na_dot_nb < 0.999d0) .and. 
     &               (sa_dot_sb > 0.999d0)) then
              amat(i,j) = a_colli !0.625d0
            else
              if (i < 4) then
               istart = 1
              elseif ((i>=4) .and. (i<7)) then
               istart = 4
              elseif ((i>=7) .and. (i<10)) then
               istart = 7
              else
               istart = 10
              endif
              if (j < 4) then
               jstart = 1
              elseif ((j>=4) .and. (j<7)) then
               jstart = 4
              elseif ((j>=7) .and. (j<10)) then
               jstart = 7
              else
               jstart = 10
              endif

              igliss1 = 0
              do jj = jstart,jstart+2
                  dp1 = abs(dot_product(s0(:,i),s0(:,jj)))
                  if (dp1 > 0.999d0) then
                    igliss1 = 1
                    exit
                  endif 
              enddo
              igliss2 = 0
              do ii = istart,istart+2
                  dp2 = abs(dot_product(s0(:,ii),s0(:,j)))
                  if (dp2 > 0.999d0) then
                    igliss2 = 1
                    exit
                  endif
              enddo 

              ! Glissile Junction
              !if (igliss1 .or. igliss2) then
              if ( (igliss1==1) .or. (igliss2==1) ) then
                  amat(i,j) = a_gliss !0.137d0
              ! Lomer-Cottrell Lock
              else
                  amat(i,j) = a_lomer !0.122d0
              endif
            endif 
          enddo
        enddo
                
        return
        end subroutine calc_amat

c **********************************************************************
! Subroutine for defining the dislocation interatction matrix that
! appears in the Taylor expression for the mean-free path
        subroutine calc_dmat(nslip,props,nprops,n0,amat,dmat)

        implicit none

        real*8 :: amat(nslip,nslip),dmat(nslip,nslip),props(nprops),
     &            n0(3,nslip),na_dot_nb
        integer :: nslip,nprops,ia,ib

        real*8 :: rk_inter,rk_copl

        rk_inter = props(38)
        rk_copl =  props(39)

        do ia = 1,nslip
        do ib = 1,nslip
            na_dot_nb = abs(dot_product(n0(:,ia),n0(:,ib)))
            if (na_dot_nb > 0.99d0) then
                dmat(ia,ib) = amat(ia,ib) / rk_copl**2.d0
            else
                dmat(ia,ib) = amat(ia,ib) / rk_inter**2.d0
            endif
        enddo
        enddo

        return
        end subroutine calc_dmat

C **********************************************************************
        SUBROUTINE calc_HARDENING(props,nprops,nslip,theta,chi,amat,dmat,
     &  s,tau,rhossd,xnu,DGAMMA,DGT,DTIME,ds,drhossd)

        IMPLICIT REAL*8(A-H,O-Z)
C
        DIMENSION DGAMMA(nslip),props(nprops),ds(nslip),s(nslip),tau(nslip),
     &            rhossd(nslip),drhossd(nslip),rhossd_s(nslip),dmat(nslip,nslip),
     &            fab(nslip,nslip),rhossd_f(nslip),rhossd_d(nslip),rhossd_i(nslip),
     &            amat(nslip,nslip),xnu(nslip),DGT(nslip)

!  Initialize props needed for the calculation

        alphaT = props(12)
        burg   = props(13)
        e_D    = props(23)
        e_copl = props(24)
        t0     = props(25)
        THETAP = props(26)
        bk     = props(27)
        xkappa0 = props(29)
        xkappa2 = props(30)
        gdotr   = props(35)
        d_gr    = props(36)
        rx      = props(37)
        rs      = props(40)
        xmu0    = props(42)

        call compute_elastConst(props,nprops,theta,c11,c12,c44)
        ! xmu = sqrt((c11-c12)/2.d0*c44)
        xmu = sqrt((c11-c12)/2.d0*c44) + xmu0

!  Calculate the hardening matrix
!        call calc_dmat(nslip,props,nprops,xmC,amat,dmat)

        rhossd_f = sqrt(matmul(amat,rhossd))
        rhossd_i = matmul(dmat,rhossd)
        rhossd_d = sqrt(rhossd_i)
        do ia = 1,nslip
          do ib = 1,nslip
            fab(ia,ib) = (xmu/rs)*burg*amat(ia,ib)/2.d0/rhossd_f(ia)
          enddo
        enddo

        drhossd = 0.d0
        do ia = 1,nslip
          xkappa1 = xkappa0 + 
     &             xkappa2*(1.d0 + dabs(dgamma(ia)/dtime)/gdotr)/DSQRT(d_gr)
          if (DGT(ia) == 0) then
             xnu(ia) = 0.
          else
             xnu(ia) = DLOG(THETAP/theta) - 
     &           DLOG(DLOG(rx*burg*rhossd(ia)/(t0*rhossd_d(ia)*DABS(DGT(ia)/dtime))))
          endif
          rhossd_s(ia) = (1.d0/(rx*burg)**2)*
     &             DEXP(-(e_D + e_copl*(burg**2)*rhossd_i(ia))/chi)
          drhossd(ia) = (1.d0/rx**2)*xkappa1*tau(ia)*dgamma(ia)*
     &             (1.d0 - rhossd(ia)/rhossd_s(ia))/(alphaT*xmu*((xnu(ia))**2)*(burg**2))
        enddo
        
        DO ia = 1, NSLIP
          ds(ia) = 0.0
          DO ib = 1,NSLIP
            ds(ia) = ds(ia) + fab(ia,ib)*drhossd(ib)
          ENDDO        
        ENDDO

        RETURN 
        END      
C**************************************************************
      subroutine compute_elastConst(props,nprops,theta,c11,c12,c44)

      implicit none

! Global Variables
      integer :: nprops

      real*8  :: theta

      real*8  :: props(nprops)

      real*8  :: c11,c12,c44

! Local Variables
      real*8  :: rm11,rm12,rm44,c11_0,c12_0,c44_0

      rm11  = props(6)
      c11_0 = props(7)
      rm12  = props(8)
      c12_0 = props(9)
      rm44  = props(10)
      c44_0 = props(11)

      c11 = c11_0 + rm11*theta
      c12 = c12_0 + rm12*theta
      c44 = c44_0 + rm44*theta

      return
      end subroutine
C************************************************************************
       SUBROUTINE DISP(U,KSTEP,KINC,TIME,NODE,NOEL,JDOF,COORDS)
C
       INCLUDE 'vaba_param.inc'
C
       REAL*8 L0
       DIMENSION U(3),TIME(2),COORDS(3)
C
       L0 = 1.D0
       EDOT = -1.0D-1
C
       U(1)=L0*(DEXP(EDOT*TIME(2))-1.D0)
C
       RETURN
       END
