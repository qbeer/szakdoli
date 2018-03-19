*****************************

      subroutine smear

      implicit none

      include"common"

*     storing field for smearing algorithm
      real*8 smsto(0:4,-maxx:maxx,-maxx:maxx,-maxz:maxz)
      real*8 wei(-2:2,2), fac1, fac2
      integer ism, kk, ll, mm, i, k, l, m

      wei(-2,2) = 1./16.
      wei(-1,2) = 1./4.
      wei( 0,2) = 3./8.
      wei( 1,2) = 1./4.
      wei( 2,2) = 1./16.

      wei(-2,1) = 0.
      wei(-1,1) = 1./4.
      wei( 0,1) = 1./2.
      wei( 1,1) = 1./4.
      wei( 2,1) = 0.


      fac1  = 1.0/12.0
      fac2  = 1.0/2.0



      ism = 1



*     loop over all lorentz-indices
      do i = 0, 4

*     loop over space-dimensions
        do k = -maxx+ism ,maxx-ism
          do l = -maxx+ism, maxx-ism
            do m = -maxz+ism ,maxz-ism

              smsto(i,k,l,m) = 0.0
               do kk = -ism , ism
                 do ll = -ism, ism
                   do mm = -ism, ism
                    smsto(i,k,l,m) = smsto(i,k,l,m)+ wei(kk,ism)*
     +                                wei(ll,ism)*wei(mm,ism)*
     +                               rhob_4(i,k+kk,l+ll,m+mm)
                   end do
                 end do
               end do
            end do
          end do
        end do
      end do



*     loop over all lorentz-indices
      do i = 0, 4

*     loop over space-dimensions
        do k = -maxx+ism ,maxx-ism
          do l = -maxx+ism, maxx-ism
            do m = -maxz+ism ,maxz-ism

               rhob_4(i,k,l,m) =  smsto(i,k,l,m)

            end do
          end do
        end do
      end do

******************************
******************************



      return
      end
