      subroutine dec_1520(exmass,bmass,wmass,bapot,mespot,pxc,pyc,pzc,
     &   pabsh)

      implicit none
      include "common"

      real*8 exmass,bmass,wmass,bapot,mespot,pxc,pyc,pzc,pabsh

      bmass = rmass
      bapot = 0
      mespot =0
      pxc=0.01
      pyc=0.
      pzc=0
      pabsh=.01
      wmass = sqrt((exmass - sqrt(rmass**2+pabsh**2))**2
     &   - pabsh**2)

      return
      end
