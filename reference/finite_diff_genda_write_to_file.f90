      program finite_diff
      implicit none

      integer, parameter :: n=2001, m=100,nn=2002
      integer :: i,j,nt,l,k
      double precision :: mm(nn),rr(nn,m),uu(nn,m),pp(nn,m),qq(nn,m),rho(nn,m),temp(nn,m)
      double precision :: gamma,pi,lambda,dt,Mp,GG,c0,p0,rho0,r0,Matm,cq,vesc,Mloss
      double precision :: us,lambda0,rad,mass,dr,gamma_a,rlast,tt,t0,temp0,Ma,Rgas,watermass
      double precision :: coefficient,dt0,dt_cfl
!      character(len=1024) :: filename
!      character(len=1024) :: format_string

      !parameters
      gamma=1.4d0
      gamma_a=1.4d0
      Mp=5.972d24
      GG=6.67408d-11
      pi=3.14159265d0
      Ma=0.029d0
      Rgas=8.3d0
      temp0=288d0
      p0=1.01e5!1.598e13*exp(-4.96e5/Rgas/temp0)
      rho0=p0*Ma/temp0/Rgas!1d-4 !wait how did I choose these values?
      r0=6.4d6
      cq=0.75d0



      coefficient=4d0*pi !4*pi for the sphere? I have to double check everything...

      lambda0=GG*Mp*rho0/r0/p0
      c0=sqrt(gamma_a*p0/rho0)

      vesc=sqrt(2.d0*GG*Mp/r0)
      us=0.5*vesc


      rlast=r0/(1d0-gamma_a/(gamma_a-1d0)/lambda0)
      t0=r0/c0
 

!      open(1,file='Mt.dat')
      
      !Making the initial conditions
      mm(1)=0d0
      do i=1,n
         !uu(i,1)=0d0
         rr(i,1)=r0+(rlast-r0)/(n)*(i-1)
         rho(i,1)=rho0*((gamma_a-1)/gamma_a*lambda0*(r0/rr(i,1)-1d0)+1d0)**(1d0/(gamma_a-1))
         pp(i,1)=p0*((gamma_a-1)/gamma_a*lambda0*(r0/rr(i,1)-1d0)+1d0)**(gamma_a/(gamma_a-1))
         if (i>1)then
            mm(i)=mm(i-1)+rho(i,1)*4*pi/3*(rr(i,1)**3d0-rr(i-1,1)**3d0)
         endif
         !ideal gas
         temp(i,1)=pp(i,1)/rho(i,1)*Ma/Rgas!temp0*((gamma_a-1d0)/gamma_a*lambda0*(r0/rr(i,1)-1d0)+1d0)
         !adiabat
         !temp(i,1)=temp0*(rho(i,1)/rho0)**(gamma_a-1d0)
         !equlibrium
         !temp(i,1)=-6.03*1d4/dlog(pp(i,1)/(2.23*1d13))
         uu(i,1)=0d0!sqrt(gamma_a*Rgas/Ma*temp(i,1))
      enddo
      Matm=mm(n)
      uu(1,1)=us

      !print *,Matm
      !stop

      
      !Normalizing it
      
      do i=1,n
         uu(i,1)=uu(i,1)/c0
         rr(i,1)=rr(i,1)/r0
         rho(i,1)=rho(i,1)/rho0
         pp(i,1)=pp(i,1)/p0
         mm(i)=mm(i)/(rho0*r0**3d0)
         temp(i,1)=temp(i,1)/temp0
      enddo

      pp(n+1,1)=0d0
      rho(n+1,1)=0d0

      !print *,temp(1,1)
      !stop

      
      !time-dependent part

      !dt=0.01/t0
      dt0=0.001/t0
      dt=dt0
      tt=0d0

      nt=40000*5*2
      nt=1000
      !nt=15

      do l=1,nt

!          if (l < 10) then
!            format_string = "(I1)"
!          else if (10 < l < 100) then
!            format_string = "(I2)"
!          else if (100 < l < 1000) then
!             format_string = "(I3)"
!          else
!              format_string = "(I4)"
!          endif


          open(1, file='fortran_file_test/output_'//trim(str(l))//'.dat')

         do j=1,1   
            do i=1,n-1
               if (uu(i+1,j)<uu(i,j))then
                  qq(i,j)=-cq*gamma*rho(i,j)*(uu(i+1,j)-uu(i,j))*(sqrt(pp(i,j)/rho(i,j))-(gamma+1d0)/2d0*(uu(i+1,j)-uu(i,j)))
               else
                  qq(i,j)=0d0
               endif
            enddo
            qq(n,j)=0.d0

            uu(1,j+1)=uu(1,j)-lambda0/gamma/(rr(1,j)**2d0)*dt !boundary condition
            do i=2,n-1
               uu(i,j+1)=uu(i,j)-(2d0*coefficient/gamma*rr(i,j)**2d0*(pp(i,j)-pp(i-1,j)+qq(i,j)-qq(i-1,j))/&
       &             (mm(i+1)-mm(i-1))+lambda0/gamma/rr(i,j)**2d0)*dt

            enddo
            uu(n,j+1)=uu(n-1,j+1) !free slip

            do i=1,n
               rr(i,j+1)=rr(i,j)+uu(i,j+1)*dt
            enddo
       
            do i=1,n-1
               pp(i,j+1)=pp(i,j)-coefficient*rho(i,j)*(gamma*pp(i,j)+(gamma-1d0)*qq(i,j))*(rr(i+1,j)**2d0* &
       &           uu(i+1,j)-rr(i,j)**2.d0*uu(i,j))/(mm(i+1)-mm(i))*dt
               rho(i,j+1)=3d0/coefficient*(mm(i+1)-mm(i))/(rr(i+1,j+1)**3d0-rr(i,j+1)**3d0)
            enddo

            do i=1,n
               temp(i,j+1)=pp(i,j+1)/rho(i,j+1)*Ma/Rgas
               !temp(i,j+1)=(rho(i,j+1))**(gamma_a-1d0)
               !temp(i,j+1)=-6.03*1d4/dlog(pp(i,j+1)/(2.23*1d13))/temp0
            enddo

            !outer boundary codntiion
            pp(n+1,j+1)=0d0
            rho(n+1,j+1)=0d0


            
            if(mod(l,100)==0)then
               do i=1,n
                  if  (uu(i,1)*c0/sqrt(2.d0*GG*Mp/(r0*rr(i,1)))>1d0)then
                     Mloss=1.0-mm(i)/mm(n)
!                     write(1,*)tt*t0,Mloss
                     exit
                  endif
               enddo
             
            endif
            
              
            tt=tt+dt

            do i=1,n
               pp(i,j)=pp(i,j+1)
               rho(i,j)=rho(i,j+1)
               temp(i,j)=temp(i,j+1)
               uu(i,j)=uu(i,j+1)
               rr(i,j)=rr(i,j+1)
               qq(i,j)=qq(i,j+1)
            enddo

            !CFL
            dt=dt0
            do i=1,n-1
               if  (dt > (rr(i+1,j)-rr(i,j))/uu(i,j))then
                  dt_cfl= 0.25*(rr(i+1,j)-rr(i,j))/uu(i,j)
                  dt=dt_cfl
                  
                  !print *,'CFL', dt_cfl    
                     
               endif
            enddo
            
         enddo
      enddo

      do i=1,n-1
         !if (mm(i)/mm(n)<0.80)then
         !print *,rr(i,1),pp(i,1),uu(i,1)*c0/sqrt(2.d0*GG*Mp/(r0*rr(i,1))),rho(i,1),temp(i,1),mm(i)/mm(n),tt*t0

         print *,rr(i,1),pp(i,1),uu(i,1)*c0/sqrt(2.d0*GG*Mp/r0),rho(i,1),temp(i,1),mm(i)/mm(n),tt*t0
         write(1, *) rr(i,1) * r0,pp(i,1) * p0,uu(i,1)*c0/sqrt(2.d0*GG*Mp/r0),rho(i,1) * rho0,temp(i,1) * temp0,mm(i),tt*t0
         !endif
      enddo
         
      end program finite_diff