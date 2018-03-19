************************************************************************
*                                                                      *
      subroutine checkl
*                                                                      *
************************************************************************
      implicit none

      include"common"
      include"cominput"

      integer i

      do i = 1, maxpar
c        if(id(1,i).lt. 0) then
c          write(*,*)'check:antiparticle ',i, id(1,i), id(2,i)
c        end if
        if(id(1,i).gt. 0) then
          if(id(1,i).eq.1) then
            if(id(2,i).lt.0 .or. id(2,i).gt.1) then
              write(*,*)'id falsch ',i, id(1,i), id(2,i)
		          	write(50,*) "stop HS 4"
            stop
          end if
        end if

        if(id(1,i).ge.2 .and. id(1,i).le.(nres+1)) then
          if(resprop2(id(1,i)-1,1).eq.1) then
            if(id(2,i).lt.0 .or. id(2,i).gt.1) then
              write(*,*)'id falsch ',i, id(1,i), id(2,i)
			  	write(50,*) "stop HS 5"
              stop
            end if
          else if(resprop2(id(1,i)-1,1).eq.3) then
            if(id(2,i).lt.-1 .or. id(2,i).gt.2) then
              write(*,*)'ladg falsch ',i, id(1,i), id(2,i)
          			write(50,*) "stop HS 6"
              stop
            end if
          end if
        end if

        if(id(1,i).eq.(nres+2)) then
          if(id(2,i).ne.0) then
            write(*,*)'id falsch ',i, id(1,i), id(2,i)
          end if
        end if

        if(id(1,i).eq.(nres+3)) then
          if(id(2,i).lt.-1 .or. id(2,i).gt.1) then
            write(*,*)'id falsch ',i, id(1,i), id(2,i)
          end if
        end if

        if(e(i).le.0.8) then
          write(*,*)'problem with mass in checkl ', i, e(i),id(1,i)
        end if
        end if
      end do

      return
      end




