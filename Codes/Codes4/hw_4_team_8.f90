c-------------------------------------------------------------------------
c..AN EXPLICIT/IMPLICIT FD SOLVER FOR 1-D CONVECTION/DIFFUSION EQUATIONS |
c-------------------------------------------------------------------------
      program CONVECTION_DIFFUSION
      parameter (mxi=501)
      character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc
c...  TO SEPERATE FILES FROM DIFFERENT SOLUTIONS.      
      print *, 'enter filename to output filename-ext.curve'
      read *, abc
c..Read the input data, and initialize the problem
c...   INITIALIZATION FOR CONDUCTION/DIFFISION
c      call INIT_CD
c...   INITIALIZATION FOR CONDUCTION/DIFFISION
c      call INIT_CD_INS
c...   INITIALIZATION FOR CONVECTION/DIFFISION
        call INIT_DC
C...  INITIALIZATION FOR LINEAR CONVECTION 
c       call INIT_LC
c..Start the solution loop 
      DO nt = 2,ntmax
c..Solve for q^n+1
c................  CONDUCTION/DIFFISUON SOLUTIONS...................

c        call EXPLICIT_CD_FTBS ! CONDUCTION/DIFFISUON FTBS   !!!!! 
c        call EXPLICIT_CD_FTCS ! CONDUCTION/DIFFUSION FTCS   !!!!!
c        call EXPLICIT_CD_FTFS ! CONDUCTION/DIFFUSION FTFS   !!!!!  

c        call EXPLICIT_CD_INS  ! CONDUCTION/DIFFUSION FTCS INSULATE 1!!
                               ! BOUNDARY AT X=0

c        call IMPLICIT_CD      ! CONDUCTION/DIFFISION CODE !  

c        call IMPLICIT_CD_INS  ! CONDUCTION/DIFFISOIN CODE WITH INSULATED
c                              !  BOUNDARY AT X=0  !!!!!
C        call IMPLICIT_NICOLSON_INS ! NICOLSON METHOD FOR INSULATED
                                    ! BOUNDARY AT X=0
c       call IMPLICIT_NICOLSON    !CRANK NICOLSON SOLUTION FOR
                                  !CONDUCTION 
c..............    LINEAR CONVECTION SOLUTIONS..............

c        call EXPLICIT_LC_FTBS ! LINEAR CONVECTION FTBS 
c        call EXPLICIT_LC_FTCS ! LINEAR CONVECTION FTCS 
c        call IMPLICIT_LC      ! LINEAR CONVECTION CODE !

c....... CONVECTION/DIFFISUON SOLUTIONS....................

c        call EXPLICIT_DC      ! CONVECTION/DIFFUSION FTCS !!!
        call IMPLICIT_DC      ! CONVECTION DIFFISION CODE !
cc..Update the solution  
        do i = 2,imax-1
           qn(i) = qnp(i)
        enddo
c..Output intermediate/final solutions
         if( MOD(nt,iof) .eq. 0 .or .nt. eq. ntmax ) call IO(nt)
      ENDDO                        

      stop
      end     

c------------------------------------------------------------------------
c.... CONDUCTION DIFFISUON PROBLEM INITIALIZATION..............      
      subroutine  INIT_CD
      parameter (mxi=501)
            character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc

      print*, 'Enter Diffusion number, ntmax and iof :' !COND/DIFF
      read(*,*) dnum, ntmax, iof  !COND/DIFF. 

      imax = 201
      x(1)    = 0.
      x(imax) = 10.
      dx      = (x(imax)-x(1))/(imax-1)
cc..Initialize the qn distribution COND/DIFF. SOLUTION
      T0 = 10.
      do i = 1,imax
         x(i)  = x(1)+dx*(i-1)
         qn(i) = T0
      enddo
cc..Apply BC for CONDUCTION/DIFFISION....
      qn(1)    = 20              
      qn(imax) = 100   
      call IO(1)
      return 
      end
cc..... CONDUCTION/DIFFUSION PROBLEM INSULATED BOUNDARY COND......      
      subroutine  INIT_CD_INS
      parameter (mxi=501)
            character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc

      print*, 'Enter Diffusion number, ntmax and iof :' !COND/DIFF
      read(*,*) dnum, ntmax, iof  !COND/DIFF. 


ccc...  CONV/DIFF EQUATION CONDITIONS      
      imax = 201
      x(1)    = 0.
      x(imax) = 10.
      dx      = (x(imax)-x(1))/(imax-1)
cc..Initialize the qn distribution COND/DIFF. SOLUTION
      T0 = 10.
      do i = 1,imax
         x(i)  = x(1)+dx*(i-1)
         qn(i) = T0
      enddo
cc..Apply BC for CONDUCTION/DIFFISION....
      qn(imax) = 100   
      call IO(1)
      return 
      end
c.... INITIALIZATION FOR LINEAR CONVECTION..............
      subroutine  INIT_LC
      parameter (mxi=501)
            character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc

      print*, 'Enter Courant number, ntmax and iof :'   !LINEAR CONV
       read(*,*) cnum, ntmax, iof  !LINEAR CONV.
cc...   LINEAR CONV EQUATION CONDITIONS      
         imax = 301 
        x(1)    = -10.
         x(imax) = 20.
         dx      = (x(imax)-x(1))/(imax-1)
cc..Initialize the qn distribution LINEAR CONV. SOLUTION
cc... define PI
       PI=4*atan(1.)
      do i = 1,imax
        x(i)  = x(1)+dx*(i-1)
       if (x(i) .le. PI .and. x(i) .ge. 0) then
       qn(i)=10*cos(x(i))
       else 
       qn(i) = 0 
      endif
      enddo
      call IO(1)
      return 
      end
c.... CONVECTION DIFFISUON PROBLEM INITIALIZATION..............      
      subroutine  INIT_DC
      parameter (mxi=501)
                character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc

      print*, 'Enter Diffusion number, Courant number, ntmax and iof :' !CONV/DIFF
      read(*,*) dnum, cnum , ntmax, iof  !CONV/DIFF. 

      imax = 401
      x(1)    =-10.
      x(imax) = 30.
      dx      = (x(imax)-x(1))/(imax-1)
cc..Initialize the qn distribution CONV/DIFF. SOLUTION
      PI=4*atan(1.)
      do i = 1,imax
         x(i)  = x(1)+dx*(i-1)
       if (x(i) .le. 3  .and. x(i) .ge. -3) then
       qn(i)=abs(sin(x(i)))
       else 
       qn(i) = 0 
      endif
      enddo
      call IO(1)
      return 
      end
c------------------------------------------------------------------------

c-----------------------DIFFISION/CONDUCTION EQUATION --------------
c..... FTBS DIFFISION/CONDUCTION.......
      subroutine  EXPLICIT_CD_FTBS
      parameter (mxi=501)
            character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc
c.....  at second point FTBS formula does not work 
c      use FTCS formulation...
        qnp(2) = dnum * (qn(3) -2*qn(2) + qn(1))+ qn(2)   
        do i = 3,imax-1
        qnp(i) = dnum * (qn(i)-2*qn(i-1)+qn(i-2)) + qn(i)
        enddo
        return 
        end
c....  FTFS DIFFISION/CONDUCTION..................     
      subroutine  EXPLICIT_CD_FTFS
       parameter (mxi=501)
              character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc
c..   at imax-1 FTFS formula does not work
c..   instead of FTFS use FTCS formulation..
        qnp(imax-1)= dnum*(qn(imax)-2*qn(imax-1)+qn(imax-2))+qn(imax-1)
cc      do rest in loop  
         do i = 2,imax-2
        qnp(i) = dnum * (qn(i)-2*qn(i+1)+qn(i+2)) + qn(i)
      enddo
      return 
      end
c...  FTCS DIFFISION/CONDUCTION... c.. 
      subroutine  EXPLICIT_CD_FTCS
      parameter (mxi=501)
                character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc
          do i = 2,imax-1
         qnp(i) = dnum*(qn(i+1) - 2*qn(i) + qn(i-1)) +qn(i) 
      enddo
      return 
      end
c...  FTCS DIFFISION/CONDUCTION... INSULATED CASE    c.. 
      subroutine  EXPLICIT_CD_INS
      parameter (mxi=501)
                character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc
          do i = 2,imax-1
         qnp(i) = dnum*(qn(i+1) - 2*qn(i) + qn(i-1)) +qn(i) 
      enddo
c... insulated boundary update      
      qn(1)=qnp(2)
      return 
      end
c------------------- LINEAR CONVECTION EQUATION-------------------
c...  FTBS LINEAR CONVECTION... c
      subroutine  EXPLICIT_LC_FTBS
      parameter (mxi=501)
            character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc
      do i = 2,imax-1
       qnp(i) = -cnum* (qn(i)-qn(i-1)) + qn(i) 
      enddo
c..   Boundary Cond.       
        qn(imax)=qnp(imax-1)
      return 
      end
c...  FTCS LINEAR CONVECTION... c
      subroutine  EXPLICIT_LC_FTCS
      parameter (mxi=501)
                character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc
       do i = 2,imax-1
       qnp(i) = -0.5*cnum* (qn(i+1)-qn(i-1)) + qn(i) 
       enddo
c..    boundary cond.      
       qn(imax)=qnp(imax-1)
      return 
      end
c................. CONVECTION-DIFFISION EQUATION......................
      subroutine  EXPLICIT_DC
      parameter (mxi=501)
            character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc
          do i = 2,imax-1
         qnp(i)=(-0.5*cnum+dnum)*qn(i+1)
         qnp(i)=qnp(i)+(1-2*dnum)*qn(i)+(0.5*cnum+dnum)*qn(i-1) 
      enddo
      !!apply BC insulated
      qn(imax)=qnp(imax-1)
      return 
      end
c.............. IMPLICIT CODE FOR CONDUCTION/DIFFISION ................
      subroutine  IMPLICIT_CD
      parameter (mxi=501)
            character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc
      dimension aa(mxi),ab(mxi),ac(mxi),rhs(mxi)

c..Fill the  3 diagonal and rhs arrays
      do i = 2,imax-1
       aa(i)  =-dnum 
       ab(i)  = 1+2*dnum
       ac(i)  = -dnum
       rhs(i) = qn(i)
      enddo
c..apply BC 1... 
       rhs(2)  = rhs(2) +  dnum * qn(1)
       rhs(imax-1)  = rhs(imax-1) +  dnum * qn(imax)
c..Solve the coupled system of equations (3-diagonal)
      call THOMAS(2,imax-1, aa,ab,ac,rhs)
c..Extract the solution
      do i = 2,imax-1
       qnp(i) = rhs(i)
      enddo

      return 
      end
c.............. IMPLICIT CODE FOR CONDUCTION/DIFFISION ................
c..............      INSULATED BOUNDARY CONDITIONS     ................
      subroutine  IMPLICIT_CD_INS
      parameter (mxi=501)
            character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc
      dimension aa(mxi),ab(mxi),ac(mxi),rhs(mxi)

c..Fill the  3 diagonal and rhs arrays
      do i = 2,imax-1
       aa(i)  =-dnum 
       ab(i)  = 1+2*dnum
       ac(i)  = -dnum
       rhs(i) = qn(i)
      enddo
c..apply BC 
       ab(2)=1.+dnum
       rhs(imax-1)  = rhs(imax-1) +  dnum * qn(imax)
c..Solve the coupled system of equations (3-diagonal)
      call THOMAS(2,imax-1, aa,ab,ac,rhs)
c..Extract the solution
      do i = 1,imax-1
         qnp(i) = rhs(i)
      enddo
c....   Update Boundary..       
      qn(1)=qnp(2)
      return 
      end
c....    IMPLICIT CODE FOR LINEAR CONVECTION ......................
      subroutine  IMPLICIT_LC
      parameter (mxi=501)
            character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc
      dimension aa(mxi),ab(mxi),ac(mxi),rhs(mxi)
c..Fill the  3 diagonal and rhs arrays
       do i = 2,imax-1
       aa(i)  =-0.5*cnum
       ab(i)  =1
       ac(i)  =0.5*cnum
       rhs(i) = qn(i)
       enddo
c..apply BC
       rhs(2) = rhs(2)+ (0.5*cnum)*qn(1)  
       ab(imax-1)=1+(cnum*0.5)
c..Solve the coupled system of equations (3-diagonal)
      call THOMAS(2,imax-1,aa,ab,ac,rhs)
c..Extract the solution
      do i = 2,imax-1
         qnp(i) = rhs(i)
      enddo
c..   update boundary point 
      qn(imax)=qnp(imax-1)
      return 
      end
c....    IMPLICIT CODE FOR CONVECTION DIFFISION  ......................
      subroutine  IMPLICIT_DC
      parameter (mxi=501)
            character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc
      dimension aa(mxi),ab(mxi),ac(mxi),rhs(mxi)
c..Fill the  3 diagonal and rhs arrays
      do i = 2,imax-1
       aa(i)  =(-cnum*0.5)-dnum 
       ab(i)  =1+2*dnum
       ac(i)  =(cnum*0.5)-dnum
       rhs(i) =qn(i) 
      enddo
c..apply BC
      rhs(2)       = rhs(2) + (cnum*0.5+dnum)*qn(1)
      ab(imax-1)=1+dnum+cnum*0.5
c..Solve the coupled system of equations (3-diagonal)
      call THOMAS(2,imax-1, aa,ab,ac,rhs)
c..Extract the solution
      do i = 2,imax-1
         qnp(i) = rhs(i)
      enddo
c     update boundary point
      qn(imax)=qnp(imax-1)
      return 
      end
c........  CRANK-NICOLSON METHOD FOR CONDUCTION/DIFFISUON EQUATION......................
c........  INSULATED BOUNDARY CONDITIONS      
      subroutine  IMPLICIT_NICOLSON_INS
      parameter (mxi=501)
            character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc
      dimension aa(mxi),ab(mxi),ac(mxi),rhs(mxi)
c..Fill the  3 diagonal and rhs arrays
      do i = 2,imax-1
       aa(i)  =-dnum*0.5 
       ab(i)  =1+dnum
       ac(i)  =-dnum*0.5
       rhs(i) =dnum*0.5*qn(i-1)+(1+dnum)*qn(i)+(dnum+1)*qn(i+1)
      enddo
c..apply BC
      rhs(2) = rhs(2) + (dnum*0.5)*qn(1)
      ab(imax-1) = 1-0.5*dnum 
c..Solve the coupled system of equations (3-diagonal)
      call THOMAS(2,imax-1, aa,ab,ac,rhs)
c..Extract the solution
      do i = 2,imax-1
         qnp(i) = rhs(i)
      enddo
c     update boundary point
      qn(imax)=qnp(imax-1)
      return 
      end
c........  CRANK-NICOLSON METHOD FOR CONDUCTION/DIFFISUON EQUATION......................
      subroutine  IMPLICIT_NICOLSON
      parameter (mxi=501)
            character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc
      dimension aa(mxi),ab(mxi),ac(mxi),rhs(mxi)
c..Fill the  3 diagonal and rhs arrays
      do i = 2,imax-1
       aa(i) =-dnum*0.5 
       ab(i) =1+dnum
       ac(i) =-dnum*0.5
       rhs(i) =dnum*0.5*qn(i-1)+(1+dnum)*qn(i)+(dnum+1)*qn(i+1)
      enddo
c..apply BC
      rhs(2) = rhs(2) + (dnum*0.5)*qn(1)
      rhs(imax-1) = rhs(imax-1) + (dnum*0.5)*qn(imax)
c..Solve the coupled system of equations (3-diagonal)
      call THOMAS(2,imax-1, aa,ab,ac,rhs)
c..Extract the solution
      do i = 2,imax-1
         qnp(i) = rhs(i)
      enddo
      return 
      end
c-------------------------------------------------------------------
      subroutine  IO(nt)
      parameter (mxi=501)
            character(15) :: abc
      common /var/ imax,ntmax,iof,dnum,cnum,x(mxi),qn(mxi),qnp(mxi),abc
      character fname*32,string*7,ext*5
      write(string,'(f7.5)') float(nt)/100000
      read(string,'(2x,a5)') ext
      fname = trim(abc)//'-'//ext//'.curve' 
      open(1,file=fname,form='formatted')
      do i=1,imax
         write(1,'(2e14.6)') x(i), qn(i)
      enddo
      close(1)
      return 
      end
      
c-------------------------------------------------------------------
      SUBROUTINE THOMAS(il,iu, A,B,C,F)
c............................................................
c Solution of a tridiagonal system of equations of the form
c  A(i)*x(i-1) + B(i)*x(i) + C(i)*x(i+1) = F(i)  for k=il,iu
c  the solution X(i) is stored in F(i)
c  A(il-1) and C(iu+1) are not used.
c  A,B,C,F are arrays to be filled by the caller program
c............................................................
      parameter (mxi=501)
      dimension  a(mxi),b(mxi),c(mxi),f(mxi),tmp(mxi)
      tmp(il)=c(il)/b(il)
      f(il)=f(il)/b(il)
      ilp1 = il+1
      do i=ilp1,iu
         z=1./(b(i)-a(i)*tmp(i-1))
         tmp(i)=c(i)*Z
         f(i)=(f(i)-a(i)*f(i-1))*z
      enddo
      iupil=iu+il
      do ii=ilp1,iu
         i=iupil-ii
         f(i)=f(i)-tmp(i)*f(i+1)
      enddo
      return
      end
c------------------------------------END----------------------------
