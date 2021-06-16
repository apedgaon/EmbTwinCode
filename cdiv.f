      subroutine cdiv(ar,ai,br,bi,cr,ci)
      include 'vaba_param.inc'
c
c     complex division, (cr,ci) = (ar,ai)/(br,bi)
c
      s = abs(br) + abs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end
