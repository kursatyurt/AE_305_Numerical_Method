c------------------------------------------------------------------------------|
c..A ITERATIVE/DIRECT SOLVER FOR ELLIPTIC PDEs                                 |
c  Course:  AE305                                                              |
c  Instructors:  Dr. Ismail H. TUNCER - Dr. Nilay Sezer Uzol                   |
c------------------------------------------------------------------------------|
      program ELLIPTIC
      parameter (imx=201, jmx=201)
      common /grid/ imax,jmax,dx,dy,beta2,x(imx,jmx),y(imx,jmx)
      common /var/  phik(imx,jmx),phikp1(imx,jmx),f(imx,jmx)
      data kmax,iof /5000,500/ rl2tol/1.E-7/

c..Read the input data, generate the grid data and initialize the solution
      call INIT 
c..Start the iterative solution loop 
      
      open(9,file='res.txt' ,form='formatted')
      DO k=1,kmax
c..Point iterative solutions
         call ITERATIVE
c..Line iterative solution
c        call ADI(k)
c..Update phi^k array and evaluate the residual
         rl2=0.
         do j = 1,jmax
         do i = 1,imax
           rl2 = rl2 + (phikp1(i,j) - phik(i,j))**2
           phik(i,j)  = phikp1(i,j)
         enddo
         enddo
         rl2 = SQRT(rl2/(imax*jmax))
c         print*, ' L2 norm of Residual@k =',k,rl2
         write(9,*) k,rl2
         if (rl2 .lt. rl2tol) exit
c..Output intermediate solutions
         if( MOD(k,iof).eq.0 .and. k.ne.kmax) call IO(k)
         ENDDO
      call IO(k)
c... call subroutine for cross sections         
      call IIO(k)
      return 
      end

c------------------------------------------------------------------------
      subroutine  INIT
      parameter (imx=201, jmx=201)
      common /grid/ imax,jmax,dx,dy,beta2,x(imx,jmx),y(imx,jmx)
      common /var/  phik(imx,jmx), phikp1(imx,jmx), f(imx,jmx)

c..Read/set input parameters; dx,dy,imax,jmax
      dx = 0.01
      dy = 0.01
      imax = INT(1./dx + 1)
      jmax = INT(1./dy + 1)
      beta2 = (dx/dy)**2

c..Generate the grid
      do i = 1,imax
      do j = 1,jmax
         x(i,j) = (i-1)*dx
         y(i,j) = (j-1)*dy
      enddo
      enddo

c..Apply the initial guess and outer BCs
      do i = 1,imax
      do j = 1,jmax
         phik(i,j) = 0.
         f(i,j)    = 0.
      enddo
      enddo

c..Apply the hole BCs
      do i=INT(0.4/dx)+1,INT(0.6/dx)+1
      do j=INT(0.4/dy)+1,INT(0.6/dy)+1
      phikp1(i,j)=0.01
      phik(i,j)=0.01
      end do
      end do 

c..Apply the distributed forces....
      f(INT(0.2/dx)+1,INT(0.2/dy)+1)=4.
      f(INT(0.2/dx)+1,INT(0.8/dy)+1)=4.
      f(INT(0.8/dx)+1,INT(0.2/dy)+1)=4.
      f(INT(0.8/dx)+1,INT(0.8/dy)+1)=4.
      return 
      end

c-------------------------------------------------------------------
      subroutine ITERATIVE
      parameter (imx=201, jmx=201)
      common /grid/ imax,jmax,dx,dy,beta2,x(imx,jmx),y(imx,jmx)
      common /var/  phik(imx,jmx), phikp1(imx,jmx), f(imx,jmx)
      data omega /1.6/
      integer :: a,b,c,d
c..... Corner Points of inner boundary
      a=INT(0.4/dx+1)
      b=INT(0.6/dx+1)
      c=INT(0.4/dy+1)
      d=INT(0.6/dy+1)
c..Solve for phi^k+1
      do j = 2,jmax-1
      do i = 2,imax-1
c..Implement, Jacobi, Gauss-Seidel and SOR methods
c..... point jacobi      
      if((i.ge.a .and. i.le.b .and. (j.ge.c .and.j.le.d))) then
      phikp1(i,j)=0.01
      else 
      phikp1(i,j) = 0.5/(1.0+beta2)*(phik(i-1,j)+phik(i+1,j)
     &   +beta2*(phik(i,j-1)+phik(i,j+1))-f(i,j)*dx**2) 
      endif
c..... GAUSS SEIDEL      
c      if((i.ge.a .and. i.le.b .and. (j.ge.c .and.j.le.d))) then
c      phikp1(i,j)=0.01
c      else 
c      phikp1(i,j) = 0.5/(1.0+beta2)*(phikp1(i-1,j)+phik(i+1,j)
c     &   +beta2*(phikp1(i,j-1)+phik(i,j+1))-f(i,j)*dx**2) 
c      endif
c..... SOR      
c      if((i.ge.a .and. i.le.b .and. (j.ge.c .and.j.le.d))) then
c      phikp1(i,j)=0.01
c      else 
c      phikp1(i,j) = 0.5/(1.0+beta2)*(phikp1(i-1,j)+phik(i+1,j)
c     &   +beta2*(phikp1(i,j-1)+phik(i,j+1))-f(i,j)*dx**2) 
c.. apply relaxation
c       phikp1(i,j)=omega*phikp1(i,j)+phik(i,j)*(1-omega)
c      endif
      enddo
      enddo
      return 
      end]

c-------------------------------------------------------------------

c..... ALTERNATING DIRECTION IMPLICIT SOLUTION
      subroutine ADI(k)
      parameter (imx=201, jmx=201)
      common /grid/ imax,jmax,dx,dy,beta2,x(imx,jmx),y(imx,jmx)
      common /var/  phik(imx,jmx), phikp1(imx,jmx), f(imx,jmx)
c...  DO ISWEEP
      do l=2,imax-1
      call ISWEEP(l)
      enddo
      do j = 1,jmax
       do i = 1,imax
          phik(i,j)  = phikp1(i,j)
         enddo
         enddo
       do l=2,jmax-1
c.... DO JSWEEP
         call JSWEEP(l)
        enddo
      end

c.... ISWEPT ROUTINE-----------------------------------------------------
      subroutine ISWEEP(l)
      parameter (imx=201, jmx=201)
      common /grid/ imax,jmax,dx,dy,beta2,x(imx,jmx),y(imx,jmx)
      common /var/  phik(imx,jmx), phikp1(imx,jmx), f(imx,jmx)
      dimension :: aa(imx),ab(imx),ac(imx),rhs(imx)
c....  2 different case due to hole inside
      if (l .lt. int(0.4/dx+1) .or. l.gt. int(0.6/dx+1)) then
C...  FILL THE 3 DIAGONAL MATRIX AND RHS
      do i =2,imax-1
      aa(i)=1.
      ab(i)=-2-2*beta2
      ac(i)=1.
      rhs(i)=f(i,l)*(dx**2)-beta2*(phik(i,l-1)+phik(i,l+1))
c...  APPLY BCs..
      rhs(2)=rhs(2)-phik(1,l)
      rhs(imax-1)=rhs(imax-1)-phik(imax,l)
      enddo
c.... Solve the coupled system of equation (3diagonal)
      call THOMAS(imx,2,imax-1, aa,ab,ac,rhs)
c...  Extract Solution
      do i = 2, imax-1
      phikp1(i,l)=rhs(i)
      enddo
      else
C...  FILL THE 3 DIAGONAL MATRIX AND RHS for left side
      do i=2,int(0.4/dx)
      aa(i)=1.
      ab(i)=-2-2*beta2
      ac(i)=1.
      rhs(i)=f(i,l)*dx**2-beta2*(phik(i,l-1)+phik(i,l+1))
c...  APPLY BCs..
      rhs(2)=rhs(2)-phik(1,l)
      rhs(int(0.4/dx))=rhs(int(0.4/dx))-phik(int(0.4/dx+1),l)
      enddo
c...  Solve Coupled system
      call THOMAS(imx,2,int(0.4/dx),aa,ab,ac,rhs)
c...  Extract Solution
      do i=2,int(0.4/dx)
      phikp1(i,l)=rhs(i)
      enddo
C...  FILL THE 3 DIAGONAL MATRIX AND RHS for right side
c..   SHIFT FOR THOMAS ! otherwise it does not work properly
      do i=2,int(0.4/dx)
      aa(i)=1.
      ab(i)=-2-2*beta2
      ac(i)=1.
      rhs(i)=f(i+int(0.6/dx),l)*dx**2-beta2*(phik(i+int(0.6/dx),l-1)
     &              +phik(i+int(0.6/dx),l+1))
        enddo
cc...  APPLY BCs..
      rhs(2)=rhs(2)-phik(int(0.6/dx+1),l)
      rhs(int(0.4/dx))=rhs(imax-1)-phik(imax,l)
cc...  Solve Coupled system
      call THOMAS(imx,2,int(0.4/dx),aa,ab,ac,rhs)
c...  Extract Solution
      do i=2,int(0.4/dx)
      phikp1(i+int(0.6/dx),l)=rhs(i)
      enddo
       endif
      end
c------------------------------------------------------------------------

c.... JSWEEP ROUTINE------------------------------------------------------
      subroutine JSWEEP(l)
      parameter (imx=201, jmx=201)
      common /grid/ imax,jmax,dx,dy,beta2,x(imx,jmx),y(imx,jmx)
      common /var/  phik(imx,jmx), phikp1(imx,jmx), f(imx,jmx)
      dimension :: aa(imx),ab(imx),ac(imx),rhs(imx)
c....  2 different case due to hole inside
      if (l .lt. int(0.4/dx+1) .or. l.gt. int(0.6/dx+1)) then
c... no hole !
C...  FILL THE 3 DIAGONAL MATRIX AND RHS
      do j =2,jmax-1
      aa(j)=beta2
      ab(j)=-2-2*beta2
      ac(j)=beta2
      rhs(j)=f(l,j)*dx**2-(phik(l-1,j)+phik(l+1,j))
c...  APPLY BCs..
      rhs(2)=rhs(2)-beta2*phik(l,1)
      rhs(jmax-1)=rhs(jmax-1)-beta2*phik(l,jmax)
      enddo
c.... Solve the coupled system of equation (3diagonal)
      call THOMAS(imx,2,jmax-1, aa,ab,ac,rhs)
c...  Extract Solution
      do j = 2, jmax-1
      phikp1(l,j)=rhs(j)
      enddo
      else
C...  BE AWARE OF HOLE !
C...  FILL THE 3 DIAGONAL MATRIX AND RHS for left side
      do j=2,int(0.4/dy)
      aa(j)=beta2
      ab(j)=-2-2*beta2
      ac(j)=beta2
      rhs(j)=f(l,j)*dx**2-(phik(l-1,j)+phik(l+1,j))
c...  APPLY BCs..
      rhs(2)=rhs(2)-beta2*phik(l,1)
      rhs(int(0.4/dy))=rhs(int(0.4/dy))-beta2*phik(l,int(0.4/dy+1))
      enddo
c...  Solve Coupled system
      call THOMAS(imx,2,int(0.4/dy),aa,ab,ac,rhs)
c...  Extract Solution
      do i=2,int(0.4/dy)
      phikp1(l,i)=rhs(i)
      enddo
c...  FILL THE 3 DIAGONAL MATRIX AND RHS for right side
c...  SHIFT FOR THOMAS ! otherwise it does not work properly
      do j=2,int(0.4/dy)
      aa(j)=beta2
      ab(j)=-2-2*beta2
      ac(j)=beta2
      rhs(j)=f(l,j+int(0.6/dy))*dx**2-(phik(l-1,j+int(0.6/dy))
     &  +phik(l+1,j+int(0.6/dy)))
cc...  APPLY BCs..
      rhs(2)=rhs(2)-beta2*phik(l,int(0.6/dy+1))
      rhs(int(0.4/dy))=rhs(int(0.4/dy))-beta2*phik(l,jmax)
      enddo
cc...  Solve Coupled system
      call THOMAS(imx,2,int(0.4/dy),aa,ab,ac,rhs)
ccc...  Extract Solution
      do i=2,int(0.4/dy)
      phikp1(l,i+int(0.6/dy))=rhs(i)
      enddo
       endif
      end

c-------------------------------------------------------------------
      subroutine  IO(k)
c..Output solution in tecplot format
      parameter (imx=201, jmx=201)
      common /grid/ imax,jmax,dx,dy,beta2,x(imx,jmx),y(imx,jmx)
      common /var/  phik(imx,jmx), phikp1(imx,jmx), f(imx,jmx)
      character fname*32,string*8,ext*5
      write(string,'(f8.5)') float(k)/100000
      read(string,'(3x,a5)') ext
      fname = 'phi-'//ext//'.plt' 
      open(1,file=fname,form='formatted')
      write(1,*) ' variables="x","y","phi" '
      write(1,*) ' zone i=',imax, ', j=',jmax
      do j = 1,jmax
      do i = 1,imax
         write(1,*) x(i,j),y(i,j),phik(i,j)
      enddo
      enddo
      close(1)
      return
      end
c----------------------------------------------------------------------
c----- DIAGONAL VERTICAL AND HORIZONTAL DIST. TO FILE
      subroutine IIO(k)
      parameter (imx=201, jmx=201)
      common /grid/ imax,jmax,dx,dy,beta2,x(imx,jmx),y(imx,jmx)
      common /var/  phik(imx,jmx), phikp1(imx,jmx), f(imx,jmx)
c.... OPEN different files for each case. 
     open(11,file='vertical.dat',form='formatted')
      open(12,file='horizontal.dat',form='formatted')
      open(13,file='diagonal.dat',form='formatted')
      do j = 1,jmax
      do i = 1,imax
      if (i .eq. int(0.5/dx+1)) then
c...... write vertical
         write(11,*) x(i,j),y(i,j),phik(i,j)
      endif
      if (j .eq. int(0.5/dy+1)) then
c.....  write horizontal 
         write(12,*) x(i,j),y(i,j),phik(i,j)
      endif   
      if ( i .eq. j) then
c...   wirte diagonal
       write(13,*) x(i,j),y(i,j),phik(i,j)
      endif
      enddo
      enddo
      close(1)
      return
      end


c-------------------------------------------------------------------
      subroutine THOMAS(mx,il,iu,A,B,C,R)
c............................................................
c Solution of a tridiagonal system of n equations of the form
c  A(i)*x(i-1) + B(i)*x(i) + C(i)*x(i+1) = R(i)  for i=il,iu
c  the solution X(i) is stored in F(i)
c  A(il-1) and C(iu+1) are not used.
c  A,B,C,R are arrays to be provided by the user
c............................................................
      dimension  a(mx),b(mx),c(mx),r(mx),x(mx)
      x(il)=c(il)/b(il)
      r(il)=r(il)/b(il)
      do i=il+1,iu
         z=1./(b(i)-a(i)*x(i-1))
         x(i)=c(i)*z
         r(i)=(r(i)-a(i)*r(i-1))*z
      enddo
      do i=iu-1,il,-1
         r(i)=r(i)-x(i)*r(i+1)
      enddo
      return
      end

c-------------------------------------------------------------------
      subroutine GAUSS(N,A,B)
      real a(n,n),b(n)
c.. Convert to upper triangular form
      do k = 1,n-1
         if (abs(a(k,k)).gt.1.e-6) theN
         do i = k+1, n
            x = a(i,k)/a(k,k)
            do j = k+1, n
               a(i,j) = a(i,j) -a(k,j)*x
            enddo
            b(i) = b(i) - b(k)*x
         enddo
         else
           write (6,*) 'zero pivot found in line:',k
           stop
         endif
      enddo
c.. Back substitution
      do i = n,1,-1
         sum = b(i)
         if (i.lt.n) then
            do j= i+1,n
               sum = sum - a(i,j)*b(j)
            enddO
         endif
         b(i) = sum/a(i,i)
      enddo
      return
      end   
