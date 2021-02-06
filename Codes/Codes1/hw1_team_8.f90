c-------------------------------------------------------------------
c..AN EULER / ITERATIVE HEUNS / ITERATIVE RK2 SOLVER for 1st order ODEs
c                                   - AE305 Numerical Methods
c    Homework #1 Team #8
c    Prepared by 
c    Arda Ozuzun 
c    Yusuf Uzuntas
c    M. Kursat Yurt
c-------------------------------------------------------------------
      program solver
      character*40 fname
c.. select solution method 
40    print *, 'Select Solver: 1)Euler  2)Heuns  3)RK2'
      read (*,*) k
      if (k .lt. 1 .or. k .gt. 3) then
       print *, 'Wrong Selection Please Reselect'
       goto 40
      end if 
c..Open the output file 
      print*,'  Enter the output file name [result.dat]:'
      read(*,'(a)') fname
      if( fname .eq. ' ') fname = 'result'
      fname = trim(fname) // '.dat'  ! results only given in dat format 
      open(1,file=fname,form='formatted')
c..Set the Initial Conditions and print them out
      data time,velocity,esol,eer/4*0.0/
      write(1,'(4f12.5)') time, velocity, esol, eer
      n=0.0
      p=0.0
c...  Go to solution loop
      if (k .eq. 1) go to 10
      if (k .eq. 2) go to 20
      if (k .eq. 3) go to 30
c.. Read the solution parameters for Euler method 
10    print*, ' '
      print*, '  Enter StepSize and FinalTime :> '
      read(*,*) stepsize, finaltime
c..Solution loop for Euler method
         do while ( time .lt. finaltime )
         slope = ODE(time,velocity)
         velocity2 = velocity + stepsize*slope
         time      = time + stepsize
         esol=EXACT(time,velocity2)
         eer=abs((velocity2-esol)/esol)*100
         write(1,'(4f12.5)') time, velocity2, esol, eer
         velocity = velocity2
       enddo
       goto 50
c..Read the solution parameters for Heuns Method
20     print*, ' '
       print*, 'Enter StepSize, FinalTime and Number of Subiteration:> '
       read(*,*) stepsize, finaltime, n
c..Solution loop for heuns method
      do while ( time .lt. finaltime )
         slope1=ODE(time,velocity)
         velocity2=velocity+slope1*stepsize
            do i=1,n
            slope2=ODE(time+stepsize,velocity2)
            slopea=0.5*(slope1+slope2)
            velocity2=velocity+slopea*stepsize
            end do
      time=time+stepsize
      esol=EXACT(time,velocity2)
      eer=abs((velocity2-esol)/esol)*100
      write(1,'(4f12.5)') time,velocity2, esol, eer
      velocity=velocity2
      enddo  
      goto 50    
c..Read the solution parameters for RK2 method
30      print*, ' '
        print*, '  Enter StepSize, FinalTime, Number of Subiteration 
     &and p1 value :> '
        read(*,*) stepsize, finaltime, n, p
c.. Set parameters for RK2
          a2=1./(2*p)
          a1=1-a2
c..Solution loop for RK2
          do while ( time .lt. finaltime )
          slope1=ODE(time,velocity)
          velocity2=velocity+stepsize*p*slope1
           do i=1,n
             slope2=ODE(time+stepsize*p,velocity2)
             velocity2=velocity+slope2*stepsize*p
            end do
            time=time+stepsize
            velocityrk=velocity+(stepsize*((a1*slope1)+(a2*slope2)))
            esol=EXACT(time,velocityrk)
            eer=abs((velocityrk-esol)/esol)*100
           write(1,'(4f12.5)') time,velocityrk, esol, eer
       velocity=velocityrk
      enddo                                  
c..Close the output file
50    close(1)
      stop
      end     
c-------------------------------------------------------------------
c ..Define the ODE as a Fortran function
      function ODE(time,velocity)
      data grav/9.81/, a/3.14/, cd/1.15/, xmass/100./,rho/1.229/
      ODE  = grav - (rho*a*cd/(2*xmass))*velocity**2
      return
      end
c.. Define exact solution as a Fortran function
      function EXACT(time,velocity)
      EXACT =(21.0262*exp(0.933124*time)-21.0262)/(exp(0.933124*time)+1)
      return
      end

c-------------------------------------------------------------------

c-------------------------------------------------------------------
c..BONUS SOLVER for 1st order ODEs
c                                   - AE305 Numerical Methods
c    Homework #1 Team #8
c    Prepared by 
c    Arda Ozuzun 
c    Yusuf Uzuntas
c    M. Kursat Yurt
c-------------------------------------------------------------------
      program solver
      character*40 fname
c.. select solution method 
40    print *, 'Select Solver: 1)Euler  2)Heuns  3)RK2'
      read (*,*) k
      if (k .lt. 1 .or. k .gt. 3) then
       print *, 'Wrong Selection Please Reselect'
       goto 40
      end if 
c..Open the output file 
      print*,'  Enter the output file name [result.dat]:'
      read(*,'(a)') fname
      if( fname .eq. ' ') fname = 'result'
      fname = trim(fname) // '.dat'  ! results only given in dat format 
      open(1,file=fname,form='formatted')
c..Set the Initial Conditions and print them out
      data time,velocity,esol,eer/4*0.0/
c     write(1,'(4A12)') '#time', 'velocity', 'exact', 'error%'
      write(1,'(4f12.5)') time, velocity, esol, eer
      n=0.0
      p=0.0
c...  Go to solution loop
      if (k .eq. 1) go to 10
      if (k .eq. 2) go to 20
      if (k .eq. 3) go to 30
c.. Read the solution parameters for Euler method 
10    print*, ' '
      print*, '  Enter StepSize and FinalTime :> '
      read(*,*) stepsize, finaltime
c..Solution loop for Euler method
         do while ( time .lt. finaltime )
         slope = ODE(time,velocity)
         velocity2 = velocity + stepsize*slope
         time      = time + stepsize
         esol=EXACT(time,velocity2)
         eer=log10(abs(velocity2-esol))
         sss=log10(stepsize)
         write(1,'(3f12.5)') time,sss, eer
         velocity = velocity2
       enddo
       goto 50
c..Read the solution parameters for Heuns Method
20     print*, ' '
       print*, 'Enter StepSize, FinalTime and Number of Subiteration:> '
       read(*,*) stepsize, finaltime, n
c..Solution loop for heuns method
      do while ( time .lt. finaltime )
         slope1=ODE(time,velocity)
         velocity2=velocity+slope1*stepsize
            do i=1,n
            slope2=ODE(time+stepsize,velocity2)
            slopea=0.5*(slope1+slope2)
            velocity2=velocity+slopea*stepsize
            end do
      time=time+stepsize
      esol=EXACT(time,velocity2)
       eer=log10(abs(velocity2-esol))
          sss=log10(stepsize)
         write(1,'(3f12.5)') time,sss, eer
      velocity=velocity2
      enddo  
      goto 50    
c..Read the solution parameters for RK2 method
30      print*, ' '
        print*, '  Enter StepSize, FinalTime, Number of Subiteration 
     &and p1 value :> '
        read(*,*) stepsize, finaltime, n, p
c.. Set parameters for RK2
          a2=1./(2*p)
          a1=1-a2
c..Solution loop for RK2
          do while ( time .lt. finaltime )
          slope1=ODE(time,velocity)
          velocity2=velocity+stepsize*p*slope1
           do i=1,n
             slope2=ODE(time+stepsize*p,velocity2)
             velocity2=velocity+slope2*stepsize*p
            end do
            time=time+stepsize
            velocityrk=velocity+(stepsize*((a1*slope1)+(a2*slope2)))
            esol=EXACT(time,velocityrk)
            eer=log10(abs(velocityrk-esol))
            sss=log10(stepsize)
         write(1,'(3f12.5)') time,sss, eer
       velocity=velocityrk
      enddo                                  
c..Close the output file
50    close(1)
      stop
      end     
c-------------------------------------------------------------------
c ..Define the ODE as a Fortran function
      function ODE(time,velocity)
      data grav/9.81/, a/3.14/, cd/1.15/, xmass/100./,rho/1.229/
      ODE  = grav - (rho*a*cd/(2*xmass))*velocity**2
      return
      end
c.. Define exact solution as a Fortran function
      function EXACT(time,velocity)
      EXACT =(21.0262*exp(0.933124*time)-21.0262)/(exp(0.933124*time)+1)
      return
      end
