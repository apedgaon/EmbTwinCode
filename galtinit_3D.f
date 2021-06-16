      subroutine Galtinit_3D(iAltern,iAltMax,iDbgFlg,rNvect,
     $           iHaveSol,rGmin,jElem,kblock,nblock,
     $           l,kinc,rIterTol,time,nTime,
     $           rGnrm,rAlphaHat,rAlphaN,rMvect)
C
C     Return alternative initial guesses for rAlphaN and rMvect
C
C
      include 'vaba_param.inc'
C
      dimension rNvect(3), rMvect(3)
      dimension rGmin(5), jElem(nblock),time(nTime)
      dimension rOutPr(15), iOutPr(15)

      if (iAltern.eq.1) then
        continue                                            ! repeat the default initial guess,
        if (abs(rAlphaN).lt.1.d-2) rAlphaN = 1.d0
!       iDbgFlg   = 1                                       ! so that we can output debug info
      elseif (iAltern.eq.2) then
        rAlphaN   = 1.5d0*rAlphaN
        if (abs(rAlphaN).lt.1.d-2) rAlphaN = 1.5d0
      elseif (iAltern.eq.3) then
        rAlphaN   = 3.d0*rAlphaN
        if (abs(rAlphaN).lt.1.d-2) rAlphaN = 3.d0
      elseif (iAltern.eq.4) then
        rAlphaN   = 1.d0
      elseif (iAltern.eq.5) then
        rAlphaN   = 1.d1
      elseif (iAltern.eq.6) then
        rAlphaN   = 1.d2
      elseif (iAltern.eq.7) then
        rAlphaN   = 1.d3
      elseif (iAltern.eq.8) then
        rAlphaN   = 0.d0
        rMvect(1) = -rNvect(2)*sign(1.d0,rNvect(1))    ! shear
        rMvect(2) =  rNvect(1)*sign(1.d0,rNvect(1))    ! band
        rMvect(3) =  rNvect(3)*sign(1.d0,rNvect(1))
      elseif (iAltern.eq.(iAltMax-1)) then    ! last iterative try 

        call random_seed()          ! we'll use a random guess for now
        call random_number(rMvect)
        rAlphaN   = sign(10.d0*rMvect(2),rMvect(1)-rMvect(2))
        rMvect(1) = sign(rMvect(1),abs(rAlphaN)-3.d0)
        rMvect(2) = sign(rMvect(2),abs(rAlphaN)-3.d0)
        rMvect(3) = sign(rMvect(3),abs(rAlphaN)-3.d0)
        rNorm = sqrt( rMvect(1)**2 + rMvect(2)**2 + rMvect(3)**2 )
        rMvect = rMvect/rNorm

        continue ! replace the above w/ code to perform exhaustive search of parameter space?

      elseif (iAltern.eq.iAltMax) then   ! We have tried enough times

C        call XPLB_ABQERR(-3,'STOPPING TO DEBUG',0,0.d0,' ')

        if (iHaveSol.eq.1) then
          iDbgFlg = 0
          rGnrm     = rGmin(1)             ! Recall best known parameters
          rAlphaHat = rGmin(2)
          rMvect(1) = rGmin(3)
          rMvect(2) = rGmin(4)
          rMvect(3) = rGmin(5)

!         Debugging...
          write(*,*) "Retrieving min-residual parameters: ",
     $             ( rGmin(i), i=1,5)
!         ...Gniggubed

C         Print error message: Couldn't find a converged solution
          iOutPr(1) = jElem(kblock)
          iOutPr(2) = l
          iOutPr(3) = kinc
          rOutPr(1) = rGnrm
          rOutPr(2) = rIterTol
          rOutPr(3) = time(1)
          rOutPr(4) = time(2)
C          call XPLB_ABQERR(-2,'Traction continuity: '//
C     $     'residual norm = %R exceeds tolerance, (%R), '//
C     $     'in element %I, Gauss point %I, at '//
C     $     'increment number %I, step time=%R, '//
C     $     'total time=%R',iOutPr,rOutPr,' ')
          do i = 1, 15
            iOutPr(i) = 0
            rOutPr(i) = 0.d0
          enddo

        else                        ! Don't have a stored solution
          iDbgFlg = 0
          rAlphaHat = 0.d0

!         Debugging...
          write(*,*) "Something is very wrong here!"
!         ...Gniggubed

          iOutPr(1) = jelem(kblock)
          iOutPr(2) = l
          iOutPr(3) = kinc
          rOutPr(1) = time(1)
          rOutPr(2) = time(2)
C          call XPLB_ABQERR(-2,'>> TRACTION CONTINUITY: '//
C     $     'NO SOLUTION in element %I, Gauss point %I'//
C     $     'increment number %I, step time=%R, '//
C     $     'total time=%R',iOutPr,rOutPr,' ')
          do i = 1, 15
            iOutPr(i) = 0
            rOutPr(i) = 0.d0
          enddo
        endif                       ! iHaveSol

      else                        ! random initial guess
        call random_seed()
        call random_number(rMvect)
        rAlphaN   = sign(10.d0*rMvect(2),rMvect(1)-rMvect(2))
        rMvect(1) = sign(rMvect(1),abs(rAlphaN)-3.d0)
        rMvect(2) = sign(rMvect(2),abs(rAlphaN)-3.d0)
        rMvect(3) = sign(rMvect(3),abs(rAlphaN)-3.d0)
        rNorm = sqrt( rMvect(1)**2 + rMvect(2)**2 + rMvect(3)**2 )
        rMvect = rMvect/rNorm

      endif                       ! iAltern

      return
      end

