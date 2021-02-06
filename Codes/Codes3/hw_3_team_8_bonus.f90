c------------------------------------------------------------------------------|
c..A 2-D FINITE VOLUME SOLVER FOR THE LAPLACE'S EQUATION ON UNSTRUCTURED GRIDS |
c  Course:  AE305                                                              |
c  Instructors: Ismail H. TUNCER & Nilay SEZER UZOL                            |
c------------------------------------------------------------------------------|
      program FINITEVOLUME
c..mxc : Max number of cells
c..mxn : Max number of nodes
      parameter (mxc=30000,mxn=20000)
      common/pars/ ncell,nnode,ntmax,ioper,dt,uinf,vinf
      common/grid/ node(3,mxc),neigh(3,mxc),xy(2,mxn),area(mxc)
      common/vars/ qcell(mxc),dqdx(mxc),dqdy(mxc)
      call cpu_time (t1)
c..Read the input data and initialize the solution
      call INIT
c..Start the solution loop
      DO nt=1,ntmax
c..Evaluate the gradient vector for each cell
         call GRADIENT
c..Loop over all the cells and solve for q^k+1
         dqmax = 0.
         do n = 1,ncell
           dq       = dt/area(n) * TFLUX(n)
           qcell(n) = qcell(n) - dq
           dqmax    = MAX(ABS(dq),dqmax)
         enddo
         print*, ' TimeStep, dqmax : ',nt, dqmax
c..Output the intermediate/final solutions
         if(MOD(nt,ioper) .eq. 0 .and. nt .ne. ntmax) call TECout(nt)
         if(MOD(nt,ioper) .eq. 0 .and. nt .ne. ntmax) call CPout(nt)
      ENDDO
      call TECout(nt)
      call CPout(nt)
      call cpu_time(t2)
      time=t2-t1
      print * ,'Total time spended in program:', time
      stop
      end

c------------------------------------------------------------------------
      subroutine INIT
      parameter (mxc=30000,mxn=20000)
      common/pars/ ncell,nnode,ntmax,ioper,dt,uinf,vinf
      common/grid/ node(3,mxc),neigh(3,mxc),xy(2,mxn),area(mxc)
      common/vars/ qcell(mxc),dqdx(mxc),dqdy(mxc)
      character fname*16
      logical ok

c..Read the grid data
      print*,'  Enter the grid file name [grid.dat]:'
      read(*,'(a)') fname
      if( fname .eq. ' ') fname = 'grid.dat'
      inquire(FILE=fname,EXIST=ok)
      if( .not. ok ) then
         print*, '  ', fname, ' does not exist! \n\n'
         stop
      endif
      open(1,file=fname,form='formatted')
      read(1,*) ncell,nnode
      read(1,*) (no,(xy(i,n),i=1,2),n=1,nnode)
      read(1,*) (no,(node(i,n),i=1,3),(neigh(i,n),i=1,3),n=1,ncell)
      close(1)
      print*, ' # of cells :',ncell
      print*, ' # of nodes :',nnode
      print*, ' '

c..Compute the cell areas
      do n = 1,ncell
        n1 = node(1,n)
        n2 = node(2,n)
        n3 = node(3,n)
        area(n) = 0.5*((xy(1,n2)-xy(1,n1))*(xy(2,n3)-xy(2,n1)) -
     >                 (xy(2,n2)-xy(2,n1))*(xy(1,n3)-xy(1,n1))  )
        if(area(n).le. 0.) write(*,*)'Cell area',n,' is zero/negative!'
      enddo

c..Set AoA
      print*,'  Enter AoA, dt, ntmax and ioperiod :'
      read(*,*) alphad,dt,ntmax,ioper
      alpha = alphad/180.*3.1415927
      uinf  = COS(alpha)
      vinf  = SIN(alpha)

c..Initialize the solution with freestream velocity; phi = x*u_inf + y*v_inf
      do n =1,ncell
       xc = (xy(1,node(1,n))+xy(1,node(2,n))+xy(1,node(3,n)))/3.
       yc = (xy(2,node(1,n))+xy(2,node(2,n))+xy(2,node(3,n)))/3.
       qcell(n) = xc*uinf + yc*vinf
c      qcell(n) = 0.
      enddo

      return
      end

c------------------------------------------------------------------------
      function TFLUX(n)
      parameter (mxc=30000,mxn=20000)
      common/pars/ ncell,nnode,ntmax,ioper,dt,uinf,vinf
      common/grid/ node(3,mxc),neigh(3,mxc),xy(2,mxn),area(mxc)
      common/vars/ qcell(mxc),dqdx(mxc),dqdy(mxc)

      tflux = 0.
c..Sum surface fluxes over the cell faces
      do ns = 1,3
         n1=node(ns,n)
         n2=node(ns+1,n)
         if(ns .eq. 3) n2=node(1,n)
         dx = xy(1,n2)-xy(1,n1)
         dy = xy(2,n2)-xy(2,n1)
         ne = neigh(ns,n)
c..Apply proper BCs when computing the face fluxes
         if ( ne .gt. 0 ) then        !..real neighbor...
           f =(dqdx(n)+dqdx(ne))*0.5
           g =(dqdy(n)+dqdy(ne))*0.5
         elseif ( ne .eq. -1 ) then   !..farfield, apply free stream velocities
           f =uinf
           g =vinf
         else                         !..wall, no flux!
           f =0.
           g =0.
         endif
         tflux = tflux - (f*dy - g*dx)
      enddo
      return
      end

c-------------------------------------------------------------------
      subroutine GRADIENT
      parameter (mxc=30000,mxn=20000)
      common/pars/ ncell,nnode,ntmax,ioper,dt,uinf,vinf
      common/grid/ node(3,mxc),neigh(3,mxc),xy(2,mxn),area(mxc)
      common/vars/ qcell(mxc),dqdx(mxc),dqdy(mxc)

      DO n = 1,ncell
      dqdx(n) = 0.
      dqdy(n) = 0.
      do ns = 1,3
         ne = neigh(ns,n)
         if(ne .eq. -1) then  !..farfield
           dqdx(n) = uinf
           dqdy(n) = vinf
           exit
         else
           if(ne .lt. -1)  ne = n       !.. apply q_ne = q_n  at the wall
           n1=node(ns,n)
           n2=node(ns+1,n)
           if(ns .eq. 3) n2=node(1,n)
           qns     = 0.5*(qcell(n)+qcell(ne))
           dqdx(n) = dqdx(n) + qns*(xy(2,n2)-xy(2,n1))/area(n)
           dqdy(n) = dqdy(n) - qns*(xy(1,n2)-xy(1,n1))/area(n)
         endif
      enddo
      ENDDO

      return
      end

c-------------------------------------------------------------------
      subroutine TECout(nstep)
      parameter (mxc=30000,mxn=20000)
      common/pars/ ncell,nnode,ntmax,ioper,dt,uinf,vinf
      common/grid/ node(3,mxc),neigh(3,mxc),xy(2,mxn),area(mxc)
      common/vars/ qcell(mxc),dqdx(mxc),dqdy(mxc)
      dimension qnode(mxn),Unode(mxn),Vnode(mxn)
      character fname*32,string*8,ext*5
c..Evaluate average q values at the nodes
      call QNODES(qnode,unode,vnode)
c..Output the solution and the grid in Tecplot format
      write(string,'(f8.5)') float(nstep)/100000
      read(string,'(3x,a5)') ext
      fname = 'q-'//ext//'.plt'           !...set the output file name
      open(1,file=fname, form='formatted')
      write(1,100) nnode,ncell
      write(1,101)(xy(1,n),xy(2,n),qnode(n),unode(n),vnode(n),n=1,nnode)
      write(1,102)(node(1,n),node(2,n),node(3,n),n=1,ncell)
      close(1)
  100 format (' VARIABLES= "x","y","phi","u","v"',/,
     >        ' ZONE N=', I6,' E=', I6,' F=FEPOINT ',' ET=triangle'  )
  101 format (5(1x,e12.5))
  102 format (3(1x,i6))
      return
      end
c-------------------------------------------------------------------
      subroutine QNODES(qnode,unode,vnode)
c..Evaluate variables at the nodes by averaging the cell values
      parameter (mxc=30000,mxn=20000)
      common/pars/ ncell,nnode,ntmax,ioper,dt,uinf,vinf
      common/grid/ node(3,mxc),neigh(3,mxc),xy(2,mxn),area(mxc)
      common/vars/ Qcell(mxc),dqdx(mxc),dqdy(mxc)
      dimension npass(mxn),qnode(mxn),unode(mxn),vnode(mxn)
      do n=1,nnode
         qnode(n) = 0.
         unode(n) = 0.
         vnode(n) = 0.
         npass(n) = 0
      enddo
c..Find the contribution of cells to the node values
      do n=1,ncell
      do ns=1,3
         nn = node(ns,n)
         qnode(nn)=qnode(nn)+qcell(n)
         unode(nn)=unode(nn)+dqdx(n)
         vnode(nn)=vnode(nn)+dqdy(n)
         npass(nn)=npass(nn)+1
      enddo
      enddo
c..Average the total node value with # of contributing cells
      do n=1,nnode
         qnode(n)=qnode(n)/npass(n)
         unode(n)=unode(n)/npass(n)
         vnode(n)=vnode(n)/npass(n)
      enddo
      return
      end

c-------------------------------------------------------------------
        subroutine CPOUT(nt)
        parameter (mxc=30000,mxn=20000)
        common/pars/ ncell,nnode,ntmax,ioper,dt,uinf,vinf
        common/grid/ node(3,mxc),neigh(3,mxc),xy(2,mxn),area(mxc)
        common/vars/ qcell(mxc),dqdx(mxc),dqdy(mxc)
        dimension qnode(mxn),Unode(mxn),Vnode(mxn), cp(nnode)
        character fname*32,string*8,ext*5
        write(string,'(f8.5)') float(nt)/100000
c..Evaluate average q values at the nodes
        call QNODES(qnode,unode,vnode)
c... Evaluate Cp values at the nodes
        do n=1,nnode
        cp(n)=1-(unode(n)**2+vnode(n)**2) ! calculate cp
        enddo
        read(string,'(3x,a5)') ext
        fname = 'k-'//ext//'.dat'           !...set the output file name
        open(2,file=fname, form='formatted')
        write(2,109) ((xy(i,n),i=1,2),cp(n) ,n=1,nnode)
  109   format (3(1x,e13.6))
        close(1)
        return
        end
c------------------------------END----------------------------------
