      subroutine MInvert(v,lv,n,iflg)

C
C      ________________________________________________________
C     |                                                        |
C     |       INVERT A GENERAL MATRIX  (NAPACK SUBROUTINE)     |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         V     --ARRAY CONTAINING MATRIX                |
C     |                                                        |
C     |         LV    --LEADING (ROW) DIMENSION OF ARRAY V     |
C     |                                                        |
C     |         N     --DIMENSION OF MATRIX STORED IN ARRAY V  |
C     |                                                        |
C     |         W     --INTEGER WORK ARRAY WITH AT LEAST N-1   |
C     |                      ELEMENTS                          |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         V     --INVERSE                                |
C     |                                                        |
C     |________________________________________________________|
C

      include 'vaba_param.inc'

      intrinsic   abs

      dimension   v(lv,*), iw(n), rvo(n,n)

!      real*8  v(lv,*),s,t
!      integer w(n),i,j,k,l,m,n,p

C     -------------------------------------------------
C     |*** INITIALIZE ERROR FLAG AND BACKUP MATRIX ***|
C     -------------------------------------------------

      iflg = 0
      do i = 1,n
        do j = 1,n
          rvo(i,j) = v(i,j)
        enddo
      enddo


      if ( n .eq. 1 ) goto 110

      l = 0
      m = 1

10    if ( l .eq. n ) goto 90

      k = l
      l = m
      m = m + 1

C     ---------------------------------------
C     |*** FIND PIVOT AND START ROW SWAP ***|
C     ---------------------------------------

      ip = l
      if ( m .gt. n ) goto 30
      s = abs(v(l,l))

      do 20 i = m,n
           t = abs(v(i,l))
           if ( t .le. s ) goto 20
           ip = i
           s  = t
20    continue

      iw(l) = ip

30    s = v(ip,l)
      v(ip,l) = v(l,l)

      if ( s .eq. 0.d0 ) goto 120

C     -----------------------------
C     |*** COMPUTE MULTIPLIERS ***|
C     -----------------------------

      v(l,l) = -1.d0
      s = 1.d0/s

      do 40 i = 1,n
40         v(i,l) = -s*v(i,l)

      j = l
50    j = j + 1

      if ( j .gt. n ) j = 1
      if ( j .eq. l ) goto 10

      t       = v(ip,j)
      v(ip,j) = v(l,j)
      v(l,j) = t

      if ( t .eq. 0.d0 ) goto 50

C     ------------------------------
C     |*** ELIMINATE BY COLUMNS ***|
C     ------------------------------

      if ( k .eq. 0 ) goto 70

      do 60 i = 1,k
60         v(i,j) = v(i,j) + t*v(i,l)

70    v(l,j) = s*t
      if ( m .gt. n ) goto 50

      do 80 i = m,n
80         v(i,j) = v(i,j) + t*v(i,l)

      GOTO 50

C     -----------------------
C     |*** PIVOT COLUMNS ***|
C     -----------------------

90    l = iw(k)

      do 100 i = 1,n
           t = v(i,l)
           v(i,l) = v(i,k)
100        v(i,k) = t

      k = k - 1

      if ( k .gt. 0 ) goto 90
      return

110   if ( v(1,1) .eq. 0. ) goto 120
      v(1,1) = 1.d0/v(1,1)
      return

C     ---------------------------------------------------
C     |*** SET ERROR FLAG AND RETURN ORIGINAL MATRIX ***|
C     ---------------------------------------------------

120   iflg = 1
      do i = 1,n
        do j = 1,n
          v(i,j) = rvo(i,j)
        enddo
      enddo

!120   call XPLB_ABQERR(-3,'MInvert error: '//
!     $  'Trying to invert singular matrix.',0,0.d0,' ')

      end

