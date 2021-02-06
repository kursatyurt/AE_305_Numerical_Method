c------------------------------------------------------------------------
c  RK4 SOLVER for a system of ODEs                                      |
c  Course:  AE305                                                       |
c  Instructors: Ismail H. TUNCER & Nilay SEZER UZOL                     |
c  Team 8 HW2                                                           |
c------------------------------------------------------------------------
      program sysRK4
      parameter (neq=3)
      real y(neq), z1(neq), z2(neq), e(neq)
      character*40 fname,str, fname2
      real, parameter :: eallow=0.01
      real, parameter :: q=1.48e-5
      real :: dt=0.1, tmax=1e6
c..   Initialize values
      i=0
      y3i=.5
      y2i=0.
      y2temp=0.
      y3temp=0.
c..   Ask output filename to user
      print*,'  Enter the output file name [result.dat]:'
      read(*,'(a)') fname
      if( fname .eq. ' ') fname = 'result'

c..    do loop for root searching
      do while ( abs(y(2)-1) .gt. eallow)
c..   Open the output file 
c..   result- iteration number.dat is default
c..   every iteration step outputs 1 new file
      write(str,'(I3)') i
      fname2 = trim(fname)//'-'//trim(str)//'.dat'  ! results only given in dat format
      open(1,file=fname2,form='formatted')
c..   Set the initial conditions
      l=0  ! stop mechanicsm counter 
      time = 0.
      y(1) = 0.
      y(2) = 0.
      y(3) = y3i
      V=time*y(2)-y(1)
      y01=time*sqrt(q*0.1)
      y05=time*sqrt(q*0.5)
      y1=time*sqrt(q*1.)
      v01=0.5*sqrt(q/0.1)*(time*y(2)-y(1))
      v05=0.5*sqrt(q/0.5)*(time*y(2)-y(1))
      v1=0.5*sqrt(q/1.)*(time*y(2)-y(1))
      write(1,'(11E15.7)') time,(y(n),n=1,neq),V,y01,y05,y1,v01,v05,v1
c..solution loop
      DO WHILE ((time.lt.tmax))
      dttemp=dt      
c.. solution for dt=dt         
        call SRK4(time,dt,y,z1)
c.. solution for dt=dt*05        
        call SRK4(time,dt*0.5,y,z2)
        call SRK4(time,dt*0.5,z2,z2)
c..  error check for stepsize
        e(1)=abs((z1(1)-z2(1))/z2(1))
        e(2)=abs((z1(2)-z2(2))/z2(2))
        e(3)=abs((z1(3)-z2(3))/z2(3))
c..  set new dt
         dt=dt*(eallow/max(e(1),e(2),e(3)))**0.20
c.. whether continue or turn back
       if (max(e(1),e(2),e(3)) .gt. eallow) cycle
       if (max(e(1),e(2),e(3)) .lt. eallow) then 
c....   stop mechanicsm when y(2) converges       
           if ( abs(y(2)-z2(2)) .lt. (0.0000001) ) then
           l=l+1
           else
           l=0
           endif
c..    change f values to new ones 
        y(1)=z2(1)
        y(2)=z2(2)
        y(3)=z2(3)
        time = time + dttemp 
c..     Calculation of V, v @x=0.1 @x=0.5 @x=1.0 and y values   
        V=time*y(2)-y(1)
        y01=time*sqrt(q*0.1)
        y05=time*sqrt(q*0.5)
        y1=time*sqrt(q*1.)
        v01=0.5*sqrt(q/0.1)*(time*y(2)-y(1))
        v05=0.5*sqrt(q/0.5)*(time*y(2)-y(1))
        v1=0.5*sqrt(q/1.)*(time*y(2)-y(1))
        write(1,'(11E15.7)') time,(y(n),n=1,neq),V,y01,y05,y1,v01,v05,v1
        endif
         
        if (l.gt.10) exit
      ENDDO
      close(1)
c..   call iteration (root searching) subroutine
      CALL ITER(y3temp,y2temp,y(2),y3i,i)
      ENDDO 
      stop
      end
         
c    -----------------------------------------------------------------------
        subroutine  ITER(y3temp,y2temp,y2,y3i,i)
        y3new=(1-y3temp)*(y3i-y3temp)/(y2-y2temp)
        y2temp=y2
        y3temp=y3i
        y3i=y3new
        i=i+1
        end
c   ------------------------------------------------------------------------
       subroutine SRK4(t,dt,y,z)
       parameter (neq=3)
       real y(neq),ytmp(neq),k1(neq),k2(neq),k3(neq),k4(neq), z(neq)

       dt2 = 0.5*dt

       call ODES(t,y,k1)
       do n = 1,neq
         ytmp(n)  = y(n) + dt2*k1(n)
       enddo
       call ODES(t+dt2,ytmp,k2)
       do n = 1,neq
         ytmp(n)  = y(n) + dt2*k2(n)
       enddo 
       call ODES(t+dt2,ytmp,k3)
       do n = 1,neq
         ytmp(n)  = y(n) + dt*k3(n)
       enddo 
       call ODES(t+dt,ytmp,k4)
c..obtain the solution at t+dt and update y for the next step
       do n = 1,neq
          phi  = (k1(n) + 2*(k2(n)+k3(n)) + k4(n))/6.
          z(n) = y(n) + dt*phi
       enddo
       
       return
       end

c------------------------------------------------------------------------
       subroutine ODES(t,y,f)
       parameter (neq=3)
       real y(neq),f(neq)

c..define the ODE's/return the slopes in the "f" array
        f(1) =y(2)
        f(2) = y(3)
        f(3) =-(y(1)*y(3))*0.5

       return
       end



c----------------NON ADAPTIVE METHOD---------------------

c------------------------------------------------------------------------
c  RK4 SOLVER for a system of ODEs                                      |
c  Course:  AE305                                                       |
c  Instructors: Ismail H. TUNCER & Nilay SEZER UZOL                     |
c------------------------------------------------------------------------
      program sysRK4
      parameter (neq=3)
      real y(neq)
      character*40 fname
      real, parameter :: q=1.48e5
c..read the input data
      print*,' '
      print*,' Enter the step size and the end point :>'
      read(*,*) dt,tmax

c..open the output file 
      print*,'  Enter the output file name [solution.dat]:'
      read(*,'(a)') fname
      if( fname .eq. ' ') fname = 'solution.dat'
      open(1,file=fname,form='formatted')

c..set the initial conditions
      time = 0.
      y(1) = 0.
      y(2) = 0.
      y(3) = 0.33357
      V=time*y(2)-y(1)
      y01=time*sqrt(q*0.1)
      y05=time*sqrt(q*0.5)
      y1=time*sqrt(q*1.)
      v01=0.5*sqrt(q/0.1)*(time*y(2)-y(1))
      v05=0.5*sqrt(q/0.5)*(time*y(2)-y(1))
      v1=0.5*sqrt(q/1.)*(time*y(2)-y(1))
      write(1,'(11E15.7)') time,(y(n),n=1,neq),V,y01,y05,y1,v01,v05,v1
c..solution loop
      DO WHILE (time.lt.tmax)      
      call SRK4(time,dt,y)
      time = time + dt 
      V=time*y(2)-y(1)
      y01=time*sqrt(q*0.1)
      y05=time*sqrt(q*0.5)
      y1=time*sqrt(q*1.)
      v01=0.5*sqrt(q/0.1)*(time*y(2)-y(1))
      v05=0.5*sqrt(q/0.5)*(time*y(2)-y(1))
      v1=0.5*sqrt(q/1.)*(time*y(2)-y(1))
      write(1,'(11E15.7)') time,(y(n),n=1,neq),V,y01,y05,y1,v01,v05,v1
      ENDDO

      close(1)
      stop
      end

c------------------------------------------------------------------------
       subroutine SRK4(t,dt,y)
       parameter (neq=3)
       real y(neq),ytmp(neq),k1(neq),k2(neq),k3(neq),k4(neq)

       dt2 = 0.5*dt

       call ODES(t,y,k1)
       do n = 1,neq
         ytmp(n)  = y(n) + dt2*k1(n)
       enddo
       call ODES(t+dt2,ytmp,k2)
       do n = 1,neq
         ytmp(n)  = y(n) + dt2*k2(n)
       enddo 
       call ODES(t+dt2,ytmp,k3)
       do n = 1,neq
         ytmp(n)  = y(n) + dt*k3(n)
       enddo 
       call ODES(t+dt,ytmp,k4)
c..obtain the solution at t+dt and update y for the next step
       do n = 1,neq
          phi  = (k1(n) + 2*(k2(n)+k3(n)) + k4(n))/6.
          y(n) = y(n) + dt*phi
       enddo
       
       return
       end

c------------------------------------------------------------------------
       subroutine ODES(t,y,f)
       parameter (neq=3)
       real y(neq),f(neq)

c..define the ODE's/return the slopes in the "f" array
        f(1) =y(2)
        f(2) = y(3)
        f(3) =-(y(1)*y(3))*0.5

       return
       end


