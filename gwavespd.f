      subroutine Gwavespd(rStrain,rStress,rEz,rNuz,rRho, rCd)

C     Compute the dilatational wave speed at a point. See "Explicit dynamic analysis",
C     ABAQUS (v6.8) Theory Manual, Sect. 2.4.5.
C     
C
C     Given:
C             rStrain: strain increment, defined such that
C        {rStrain} = < u1,1 | u2,2 | u3,3 | u2,3+u3,2 | u1,3+u3,1 | u1,2+u2,1 >^T
C
C             rStress: Cauchy stress increment
C                rRho: Mass density of the material
C                 rEz: Young's modulus (initial/elastic)
C                rNuz: Poisson's ratio (initial/elastic)
C
C     Return the following:
C                 rCd: Dilatational wave speed
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      include 'vaba_param.inc'

      dimension rStrain(6), rStress(6)
      dimension rDvStrn(6), rDvStrs(6)

      parameter (small = 1.d-8)
      parameter (thrd  = 1.d0/3.d0)

C     Compute volumetric strain and hydrostatic stress

      rEpVol  = 0.d0
      rHydro  = 0.d0
      do i = 1, 3
        rEpVol = rEpVol + rStrain(i)
        rHydro = rHydro + rStress(i)
      enddo
      rHydro = thrd * rHydro

C     Compute deviatoric strain and stress

      do i = 1, 3
        rDvStrn(i) = rStrain(i) - thrd*rEpVol
        rDvStrs(i) = rStress(i) -      rHydro
      enddo

      do i = 4, 6
        rDvStrn(i) = 0.5d0 * rStrain(i)
        rDvStrs(i) =         rStress(i)
      enddo

C     Compute contracted quantities

      rSEdev = 0.d0
      rEpDev = 0.d0
      do i = 1, 3
        rSEdev = rSEdev +        rDvStrs(i) * rDvStrn(i)
        rEpDev = rEpDev +        rDvStrn(i) * rDvStrn(i)
      enddo

      do i = 4, 6
        rSEdev = rSEdev + 2.d0 * rDvStrs(i) * rDvStrn(i)
        rEpDev = rEpDev + 2.d0 * rDvStrn(i) * rDvStrn(i)
      enddo

C     If insignificant volumetric strain, compute elastic (initial) Lame parameters

      if (rEpvol.lt.small) then

        rLambda = rEz*rNuz / ( (1.d0-2.d0*rNuz)*(1.d0+rNuz) ) ! elastic Lame' params.
        rMu     = rEz      / (            2.d0 *(1.d0+rNuz) )

C     If significant volumetric strain, but small deviatoric strain

      elseif (rEpDev.lt.small) then

        rLambda = rEz*rNuz / ( (1.d0-2.d0*rNuz)*(1.d0+rNuz) ) ! elastic Lame' lambda
        rMuz    = rEz      / (            2.d0 *(1.d0+rNuz) ) ! elastic shear modulus

        rKappa  = rHydro   / rEpVol                           ! effective bulk modulus

        rMu     = 0.75d0 * ( (rLambda + 2.d0*rMuz) - rKappa ) ! effective shear modulus

C     If significant volumetric and deviatoric strains

      else

        rKappa  =         rHydro / rEpVol                     ! effective bulk modulus
        rMu     = 0.5d0 * rSEdev / rEpDev                     ! effective shear modulus

        rLambda = rKappa - 2.d0 * thrd *rMu                   ! effective Lame' lambda

      endif

C     The dilatational wave speed

      rCd = sqrt( (rLambda + 2.d0*rMu) / rRho )

      end

