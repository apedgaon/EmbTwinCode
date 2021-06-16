      subroutine vumbka(
C  Read only variables
     1    nblock,ndir,nshr,nstatev,nfieldv,nprops,lanneal,
     2    stepTime,totalTime,dt,cmname,coordMp,charLength,
     3    props,density,strainInc,relSpinInc,
     4    tempOld,stretchOld,defgradOld,fieldOld,
     5    stressOld,stateOld,enerInternOld,enerInelasOld,
     6    tempNew,stretchNew,defgradNew,fieldNew,
C  Write only variables
     7    stressNew,stateNew,enerInternNew,enerInelasNew)
C
C************************************************************************

C  Isotropic Thermo-Elasto-Viscoplasticity with pressure sensitive
C  plastic flow and plastic dilatancy.  This code is based upon the 
C  work of Anand (1985) IJP 1,213 and Brown et al. (1989) IJP 5,95.
C
C  The UMAT version interfaces with the *VISCO procedure in ABAQUS.
C  Automatic timestep control is done using the CEMAX parameter.  The
C  timestep is decreased if CEMAX exceeds CETOL.
C
C  This version is not for use in plane stress or any other cases
C  where more strain terms than stress terms are used.
C
C************************************************************************
C
C  State Variables:
C    STATEV(1) = S	(plastic flow resistance, tensile)
C    STATEV(2) = TH	(temperature)
C    STATEV(3) = F	(plastic shear strain rate)
C    STATEV(4) = time integral of F    (plastic shear strain)
C
C************************************************************************
C
C Contents of *local* PROPS vector in this version:
C
C  J|	PROPS(J)
C--------------------------------------------------------------------
C  1|	PLSLMT	limit on equiv. plastic tensil strain increment
C  2|	PHI 	degree of implicitness [0,1]
C  3|	OMEGA	fraction of plastic work to heating
C  4|	S0	initial value for internal variable S
C  5|	T0	initial temperature
C  6|	AMU	shear modulus
C  7|	AKAPPA	bulk modulus
C  8|	ALPHA	thermal expansion coefficient
C  9|	RHO	density
C 10|	C	specific heat
C 11|	A	pre-exponential factor
C 12|	Q/R	activation energy divided by gas constant
C 13|	XP	exponent in equation for F (1/n)
C 14|	PALPHA	pressure sensitivity parameter
C 15|	H0	pre-multiplier in hardening rate equation
C 16|	SB	pre-multiplier for saturation value of S
C 17|	AN	exponent in equation for saturation value of S
C 18|	B	pre-multiplier for rate of static restoration
C 19|	SA	annealed value of S
C 20|	P	exponent in equation for RDOT
C 
C The above list of properties is specific to the functions chosen for 
C modeling one particular material.  For other materials a different list
C may be desirable.  For example, additional properties would be needed
C to define functions of temperature for the elastic constants and for the 
C annealed value of internal variable S.
C
C See subroutine UMPROP, GAMDOT, and SDOT for suggested units for the
C above properties.
C
C The values input for S0,A,H0,SB,B,SA should correspond to tensile
C test data.  These values are converted to equivalent shear values
C within the subroutine.
C
C The parameter PHI controls the degree of implicitness of the integration
C procedure
C	PHI = 0.0	explicit
C	PHI = 1.0	implicit
C Suggested value of PHI is >= 0.50.
C
C The parameter OMEGA controls whether the problem is isothermal or
C adiabatic
C	OMEGA = 0.0	isothermal
C	OMEGA = 1.0	adiabatic
C Suggested value for fast deformations is 0.9
C
C************************************************************************

         include 'vaba_param.inc'

        dimension props(nprops),density(nblock),coordMp(nblock,*),
     1   charLength(nblock),strainInc(nblock,ndir+nshr),
     2   relSpinInc(nblock,nshr),tempOld(nblock),
     3   stretchOld(nblock,ndir+nshr),
     4   defgradOld(nblock,ndir+nshr+nshr),
     5   fieldOld(nblock,nfieldv),stressOld(nblock,ndir+nshr),
     6   stateOld(nblock,nstatev),enerInternOld(nblock),
     7   enerInelasOld(nblock),tempNew(nblock),
     8   stretchNew(nblock,ndir+nshr),
     9   defgradNew(nblock,ndir+nshr+nshr),
     9   fieldNew(nblock,nfieldv),stress(6),
     1   stressNew(nblock,ndir+nshr),stateNew(nblock,nstatev),
     2   enerInternNew(nblock),enerInelasNew(nblock)

        character*80 cmname

!       DEBUGGING...
!       write(*,*) '>>>>>>>>>>>>>>>>>  BKA MATERIAL ROUTINE CALLED!'
!       write(*,*) 'Time = ',totalTime,' Step time = ',stepTime
!       do i = 1, (ndir+nshr)
!         write(*,*) 'dE(',i,')= ',strainInc(1,i)
!       enddo
!       ...GNIGGUBED

        do 100 km=1,nblock

         pi=4.*atan(1.0)
         sin60=sin(60.*pi/180.)
         cos60=cos(60.*pi/180.)
         zero=0.
         one=1.0
         two=2.0
         half=0.5
         sqart3=SQRT(3.)
         third=1./3.
         ndip1=ndir+1
         phidt=props(2)*dt
         ntens=ndir+nshr
         small=1.e-30
         amu=props(6)
         akappa=props(7)


C**********************************************************************
C Initialize the state variables if necessary
C*********************************************************************
        if((totalTime.lt.small).and.(stepTime.lt.small))then

!         DEBUGGING...
!         write(*,*) 'INITIALIZING STATE VARIABLES...'
!         ...GNIGGUBED

          stateNew(km,1)=props(4)/sqart3
          stateNew(km,2)=props(5)
          stateNew(km,4)=zero
          v11=akappa-two*third*amu
          dvol=strainInc(km,1)+strainInc(km,2)+strainInc(km,3)
          v4=v11*dvol   
          do k1=1,ndir
            stressNew(km,k1)=2.*amu*strainInc(km,k1)+v4
          enddo
          if(nshr.gt.0)then
            do k1=ndip1,ntens
              stressNew(km,k1)=amu*strainInc(km,k1)
            enddo
          endif
          goto 100
        endif
C**********************************************************************
C  Set the state variables
C*******************************************************************
        s=stateOld(km,1)
        th=stateOld(km,2)
C********************************************************************
C  The parameter TAUTOL is the minimum value allowed for TAUB.  This
C  avoids difficulties that could arise when dividing by TAUB.  The 
C  value used corresponds to a stress level at which F is considered
C  negligible, while 1/TAUTOL**2 can still be evaluated.
C**********************************************************
        tautol=s*1.e-6
C*********************************************************
C  Subroutine UMPROP determines AMU,AKAPPA,ALPHA,RHO, and C based
C  upon the temp. TH, using data supplied in PROPS.  For the present
C  case, the properties are assumed constant and input directly in
C  PROPS.  In other cases, additional data constants defining functions
C  of TH for each property may be input.
C***********************************************************
        call umprop(amu,akappa,alpha,rho,c,th,props,nprops)
C************************************************************
C  Pressure PB and equivalent shear stress TAUB
C***********************************************************
        do k2=1,6
          stress(k2)=stressOld(km,k2)
        enddo
        call sinv(stress,sinv1,sinv2)
        taub=max(sinv2/sqart3,tautol)
        pb=-sinv1
C****************************************************************
C  Subroutine SBETA determines the value of the plastic dilatancy
C  factor BETA.  Presently set to zero.
C***************************************************************
        call sbeta(taub,pb,th,s,beta)
C***************************************************************
        if(props(3).gt.zero)then
          con1=props(3)*(taub-beta*pb)/(rho*c)
        else
          con1=zero
        endif
C****************************************************************
C  Subroutine GAMDOT determines the equivalent plastic shear strain
C  rate F and its derivatives PDA.PDB,PDC,PDD with respect to TAUB,
C  PB,TH, and S respectively.  To make the subsequent calculations
C  more convenient, PDB, PDC, and PDD are returned as:
C	RATIOB = PDB/PDA
C	RATIOC = PDC/PDA
C	RATIOD = PDD/PDA
C  CON2 is returned for use in subroutine SDOT
C****************************************************************
        call gamdot(taub,pb,th,s,props,nprops,sqart3,
     &          f1,pda,ratiob,ratioc,ratiod,con2,pi,sin60,
     &          cos60,zero,one,two,half,third)
C*****************************************************************
C Subroutine SDOT determines the hardening rate H and the static
C restoration rate RDOT.  SA, the annealed value of S at this 
C temperature is also returned.
C****************************************************************
        call sdot(taub,pb,th,s,amu,con2,props,nprops,
     &     sqart3,f1,h,rdot,sa,pi,sin60,cos60,zero,one,two,
     &     half,third)
C*************************************************************
C  Static recovery increment
C*************************************************************
        dr=rdot*dt
C*************************************************************
C  Trace of strain increment DVOL
C***********************************************************
        dvol=zero
        do k1=1,ndir
          dvol=dvol+strainInc(km,k1)
        enddo
C************************************************************
C  Convert stress to deviatoric stress
C**********************************************************
        do k1=1,ndir
          stressOld(km,k1)=stressOld(km,k1)+pb
        enddo
C**********************************************************
C  Deviatoric stress times strain increment  SDSI
C**********************************************************
        sdsi=zero
        do k1=1,ntens
          sdsi=sdsi+stressOld(km,k1)*strainInc(km,k1)
        enddo
C***********************************************************
C  Effective plastic shear strain increment
C*********************************************************
        g=amu-(ratiob*akappa*beta+ratioc*con1+ratiod*h)
        v=phidt*pda*g
        v1=f1*dt/(one+v)
        v2=phidt*pda/(one+v)
        dgampb=v1+v2*(amu*sdsi/taub-akappa*ratiob*dvol)
C*********************************************************
C  Increments DS and DTH
C********************************************************
        ds=h*dgampb-dr
        dth=con1*dgampb
C*******************************************************
C  Constants for Jacobian and stress increment
C*******************************************************
        amub=amu*taub/(taub+amu*phidt*f1)
        v3=two*amub
        v11=akappa-two*third*amub
        v4=v11*dvol
        v5=akappa*alpha*dth/third
        v6=akappa*beta*dgampb
        v7=v4-v5-v6
        v8=amu*dgampb/taub-(amu-amub)*sdsi/taub**2
        v9=amub
C********************************************************
C Calculate complete STRESS at the end of the increment using
C the deviatoric stress at the beginning of the increment
C*******************************************************
        do k1=1,ndir
          stressNew(km,k1)=(one-v8)*stressOld(km,k1)+
     &                         v3*strainInc(km,k1)+v7-pb
        enddo
        if(nshr.gt.0)then
          do k1=ndip1,ntens
            stressNew(km,k1)=(one-v8)*stressOld(km,k1)+
     &                           v9*strainInc(km,k1)
          enddo
        endif
C***********************************************************
C Update the state variables
C************************************************************
        do k2=1,6
          stress(k2)=stressNew(km,k2)
        enddo
        call sinv(stress,sinv1,sinv2)
        taub=max(sinv2/sqart3,tautol)
        pb=-sinv1
        th=th+dth
        s=s+ds
        if(s.lt.sa)then
          s=sa
        endif 
        call gamdot(taub,pb,th,s,props,nprops,sqart3,
     &          f2,pda,ratiob,ratioc,ratiod,con2,pi,sin60,
     &          cos60,zero,one,two,half,third)
        stateNew(km,1)=s
        stateNew(km,2)=th
        stateNew(km,3)=f2
        stateNew(km,4)=stateOld(km,4)+dgampb
C
100     continue
C
!       DEBUGGING...
!       do i = 1, (ndir+nshr)
!         write(*,*) 'S(',i,')= ',stressNew(1,i)
!       enddo
!       ...GNIGGUBED
        return
        end

C************************************************************
C************************************************************
C************************************************************
C************************************************************

C************************************************************
        subroutine umprop(amu,akappa,alpha,rho,c,th,props,nprops)
C*************************************************************
C Determine the following constants for the material:
C AMU		shear modulus
C AKAPPA	bulk modulus
C ALPHA		thermal expansion coefficient
C RHO		density
C C		specific heat (J/kg/K)
C In general, these properties are functions of temperature, but in 
C the present case they are input directly assuming no temperature
C dependence.  Additional entries in PROPS could be used to define
C them as functions of temperature
C*************************************************************
        include 'vaba_param.inc'
C**************************************************************
        dimension props(nprops)
C**************************************************************
        amu=props(6)
        akappa=props(7)
        alpha=props(8)
        rho=props(9)
        c=props(10)
C***************************************************************
        return
        end
C***************************************************************
C***************************************************************
        subroutine sbeta(taub,pb,th,s,beta)
C**************************************************************
C Subroutine SBETA determines the plastic dilatancy factor
C***************************************************************
        include 'vaba_param.inc'
C****************************************************************
        beta=0.d0
C**************************************************************
        return
        end
C****************************************************************
C****************************************************************
        subroutine gamdot(taub,pb,th,s,props,nprops,sqart3,
     &          f,pda,ratiob,ratioc,ratiod,con2,pi,sin60,
     &          cos60,zero,one,two,half,third)
C****************************************************************
C Subroutine GAMDOT determines the equivalent plastic shear strain
C rate and its derivatives PDA,PDB,PDC,PDD with respect to TAUB,
C PB,TH and S repectively.  Note that the following derivative
C terms are returned:
C		RATIOB = PDB/PDA
C		RATIOC = PDC/PDA
C		RATIOD = PDD/PDA
C************************************************************
        include 'vaba_param.inc'
C***********************************************************
        dimension props(nprops)
C************************************************************
C Material parameters defining the equivalent plastic shear
C strain rate
C PROPS(11)	A	pre-exponential factor
C PROPS(12)	Q/R	activ. energy / gas constant
C PROPS(13)	1/AM	1/rate sensitivity parameter
C PROPS(14)	PALPHA	pressure sensitivity parameter
C Note:	   PALPHA=2*(DMU/DPB)/MU0  with
C		MU0   shear modulus at zero pressure
C		DMU/DPB		pressure dependency of shear modulus
C*************************************************************
        al=3.25
        con2=sqart3*props(11)*exp(-props(12)/th)
C*************************************************************
        f=con2*(sinh(taub*al/(s*(one+props(14)*pb))))**props(13)
C**************************************************************
C Note that a tolerance was set on the minimum value of TAUB allowed
C to avoid calculational difficulties
C**************************************************************
        sp=s*(one+props(14)*pb)
        pda=f*props(13)*cosh(taub*al/sp)*al/sinh(taub*al/sp)/sp
        ratiob=0.
        ratioc=props(12)*s*tanh(taub*al/s)/(al*props(13)*th**2)
        ratiod=-taub/s
C***************************************************************
        return
        end
C***************************************************************
C***************************************************************
        subroutine sdot(taub,pb,th,s,amu,con2,props,nprops,
     &     sqart3,f,h,rdot,sa,pi,sin60,cos60,zero,one,two,
     &     half,third)
C*************************************************************
C This subroutine determines the hardening rate H and the static
C restoration rate RDOT
C**************************************************************
        include 'vaba_param.inc'
C**************************************************************
        dimension props(nprops)
C***************************************************************
C Material parameters determining the rate of hardening:
C PROPS(15) = H0	pre-multiplier in hardening rate equation
C PROPS(16) = SB	pre-multiplier for saturation value of S
C PROPS(17) = AN	exponent in eqn. for saturation value of S
C PROPS(18) = B		pre-multiplier for rate of static recovery
C PROPS(19) = SA	annealed value of S
C PROPS(20) = P		exponent of ((S/SA)-1) in RDOT
C The following constant was calculated in GAMDOT
C                 CON2 = A*EXP(-Q/R*TH))
C Note that H0,SB,B, and SA should be from tensile test data because
C a conversion to pure shear data is done in this subroutine
C*************************************************************
C First calculate H
C***************************************************************
        ss=(props(16)/sqart3)*(f/con2)**props(17)
        if(s.le.ss)then
          h=third*props(15)*(one-s/ss)**1.5
        else
          h=zero
        endif
C****************************************************************
C Determine the annealed value of S.  Although in general it is a 
C function of temperature, in the present version it is input
C directly and considered constant
C************************************************************
        sa=props(19)/sqart3
C************************************************************
C Calculate RDOT
C*************************************************************
        sd=(s/sa)-one
        if(sd.gt.zero)then
          rdot=amu*third*props(18)*con2*sd**props(20)
        else
          rdot=zero
        endif
C**************************************************************
        return
        end
C**************************************************************
C**************************************************************
        subroutine sinv(stress,sinv1,sinv2)
        include 'vaba_param.inc'
        dimension stress(6)
        sinv1=(stress(1)+stress(2)+stress(3))/3.
        sdev1=stress(1)-sinv1
        sdev2=stress(2)-sinv1
        sdev3=stress(3)-sinv1
        sinv2=sdev1**2+sdev2**2+sdev3**2
        sinv2=sinv2+stress(4)**2+stress(5)**2+stress(6)**2
        sinv2=sqrt((3./2.)*sinv2)
        return
        end

