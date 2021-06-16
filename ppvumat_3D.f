      subroutine vumat (
c Read only -
     1 jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2 stepTime, totalTime, dt, cmname, coordMp, charLength,
     3 props, density, strainInc, relSpinInc,
     4 tempOld, stretchOld, defgradOld, fieldOld,
     5 stressOld, stateOld, enerInternOld, enerInelasOld,
     6 tempNew, stretchNew, defgradNew, fieldNew,
c Write only -
     7 stressNew, stateNew, enerInternNew, enerInelasNew )
c
      include 'vaba_param.inc'
c
      dimension jblock(*), props(nprops),density(*), coordMp(*),
     1 charLength(*), strainInc(*),
     2 relSpinInc(*), tempOld(*),
     3 stretchOld(*),
     4 defgradOld(*),
     5 fieldOld(*), stressOld(*),
     6 stateOld(*), enerInternOld(*),
     7 enerInelasOld(*), tempNew(*),
     8 stretchNew(*),
     9 defgradNew(*),
     1 fieldNew(*),
     2 stressNew(*), stateNew(*),
     3 enerInternNew(*), enerInelasNew(*)
c
      character*80 cmname

      parameter (
     1 i_umt_nblock = 1,
     2 i_umt_npt = 2,
     3 i_umt_layer = 3,
     4 i_umt_kspt = 4,
     5 i_umt_noel = 5 )

!     write(*,*) ' INSIDE VUMAT'

       call vumatXtrArg ( jblock(i_umt_nblock),
     1 ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2 stepTime, totalTime, dt, cmname, coordMp, charLength,
     3 props, density, strainInc, relSpinInc,
     4 tempOld, stretchOld, defgradOld, fieldOld,
     5 stressOld, stateOld, enerInternOld, enerInelasOld,
     6 tempNew, stretchNew, defgradNew, fieldNew,
     7 stressNew, stateNew, enerInternNew, enerInelasNew,
     8 jblock(i_umt_noel), jblock(i_umt_npt),
     9 jblock(i_umt_layer), jblock(i_umt_kspt))

      return
      end

c ----------------------------------------------------------------------------------

      subroutine vumatXtrArg (
c read only -
     1 nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2 stepTime, totalTime, timeinc, cmname, coordMp, charLength,
     3 props, density, strainInc, relSpinInc,
     4 tempOld, stretchOld, defgradOld, fieldOld,
     3 stressOld, stateOld, enerInternOld, enerInelasOld,
     6 tempNew, stretchNew, defgradNew, fieldNew,
c write only -
     5 stressNew, stateNew, enerInternNew, enerInelasNew,
c read only extra arguments -
     6 nElement, nMatPoint, nLayer, nSecPoint)
c
      include 'vaba_param.inc'
c
c all arrays dimensioned by (*) are not used in this algorithm
      dimension props(nprops), density(nblock),
     1 strainInc(nblock,ndir+nshr),
     2 relSpinInc(nblock,nshr), defgradOld(nblock,9),
     4 stressOld(nblock,ndir+nshr),
     5 stateOld(nblock,nstatev), enerInternOld(nblock),
     6 enerInelasOld(nblock),
     7 stretchNew(nblock,ndir+nshr), defgradNew(nblock,9),
     8 stressNew(nblock,ndir+nshr)

      dimension enerInelasNew(nblock),stateNew(nblock,nstatev),
     1 enerInternNew(nblock)

      dimension nElement(nblock)

character*80 cmname

      parameter  (maxelno=100000,maxisvs=40)
      dimension  uel_isvs(maxelno,8,maxisvs)
      common     /k_state_vars/ uel_isvs

!      write(*,*) ' INSIDE VUMATX'

      do k = 1,nblock

        if (TotalTime.eq.0.d0) then

!         write(*,*) 'Ttime=',TotalTime,'stime=',stepTime

          do i = 1, ndir+nshr
            stressNew(k,i) = 1.d0
          enddo

        else  

          do i = 1, ndir+nshr
            stressNew(k,i) = 0.d0
          enddo

          noel = nElement(k)
          print *, noel
!         npt  = nMatPoint
!         nlay = nLayer
!         nsec = nSecPoint

!         write(*,*) 'timeinc=',timeinc
!         write(*,*) 'noel=',noel,'npt=',npt
!         write(*,*) 'nlay=',nlay,'nsc=',nsec

          if (nstatev.lt.90) then
              call XPLB_ABQERR(-3,'VUMATX: need more SDVs',0,0.d0,' ')
          endif

          do i = 1, 10
            stateNew(k,   i) = uel_isvs(noel,1,i)
            stateNew(k,10+i) = uel_isvs(noel,2,i)
            stateNew(k,20+i) = uel_isvs(noel,3,i)
            stateNew(k,30+i) = uel_isvs(noel,4,i)
            stateNew(k,40+i) = uel_isvs(noel,5,i)
            stateNew(k,50+i) = uel_isvs(noel,6,i)
            stateNew(k,60+i) = uel_isvs(noel,7,i)
            stateNew(k,70+i) = uel_isvs(noel,8,i)
          enddo

C         combined / averaged quantitites for the element
          do i = 1, 4 
            stateNew(k,80+i) =( stateNew(k,   i)+stateNew(k,10+i)
     $                       +  stateNew(k,20+i)+stateNew(k,30+i)
     $                       +  stateNew(k,40+i)+stateNew(k,50+i)
     $                       +  stateNew(k,60+i)+stateNew(k,70+i) ) / 8.d0
          enddo

          stateNew(k,90) = stateNew(k,10) + stateNew(k,20)       ! localization flag
     $                   + stateNew(k,30) + stateNew(k,40)
     $                   + stateNew(k,50) + stateNew(k,60)
     $                   + stateNew(k,70) + stateNew(k,80)

          if (stateNew(k,90).gt.0.d0) then      ! element has localized points

            do i = 5, 9
              stateNew(k,80+i) =( stateNew(k,   i)*stateNew(k,10)
     $                          + stateNew(k,10+i)*stateNew(k,20)
     $                          + stateNew(k,20+i)*stateNew(k,30)
     $                          + stateNew(k,30+i)*stateNew(k,40)
     $                          + stateNew(k,40+i)*stateNew(k,50)
     $                          + stateNew(k,50+i)*stateNew(k,60)
     $                          + stateNew(k,60+i)*stateNew(k,70)
     $                          + stateNew(k,70+i)*stateNew(k,80)
     $                          )/ stateNew(k,90)
            enddo

          else

            do i = 5, 9
              stateNew(k,80+i) = stateNew(k,   i)
            enddo

          endif

!         write(*,*) 'nblock=',nblock,'nstatev=',nstatev

        endif

      enddo

!     write(*,*) 'VUMATX RETURNING'

      return
      end


