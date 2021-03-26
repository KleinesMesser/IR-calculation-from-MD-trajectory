program IR_Combined
!This program is to calculate the IR lineshap with read-in frequency,
!coupling and transition dipole trajectories.
implicit none
integer,parameter::nsnpsht=200000,tstep=1000,f_tstep=8000,interval=100
integer,parameter::nchro=14
!tstep: maxstep to calculate average integration of m(0)F(t)m(t)
!(simplified as mFm below);
!f_tstep: max step for fourier tranfrom calculation. tstep <= f_tstep.
!When index > tstep, all mFm values are treated as zero.
!This is reflected by the command line below "mFm=zero"
!interval: step interval to calculate the average of mFm
double precision,parameter::c=2.99792458d-2,PI=acos(-1.0d0),conve=2.0*PI*c
!speed of light c in cm/ps. Thus conve in cm/ps too.
double precision,parameter::dt=0.010d0,dw=2.0d0*PI/dt/real(f_tstep)
!dt in ps, which depends on the output frequency of the MD trajectory. dw in rad/ps
double precision,parameter::minfreq=1.2d3,maxfreq=2.0d3
!maximum and minimum output frequency in cm^-1.
double precision,parameter::Tl=0.649d0*2.0d0
!Vibrational life time in ps
complex*16,parameter::one=(1.0d0,0.0d0),zero=(0.0d0,0.0d0)
complex*16,parameter::im=(0.0d0,1.0d0),it=im*dt*conve
character(100)::filename(4)
integer::s,i,j,k,l,dstep
double precision::kappa(nchro,nchro),d(nchro),P(nchro,nchro)
double precision::m(nsnpsht,nchro,3)
double precision::inten(f_tstep)
double precision::time,actualw,dintensity(0:f_tstep)
complex*16::U(nsnpsht,nchro,nchro),F(nsnpsht,nchro,nchro),mFm(0:f_tstep)
complex*16::Uni(nchro,nchro)
real::start,finish
call getarg(1,filename(1))
call getarg(2,filename(2))
call getarg(3,filename(3))
call getarg(4,filename(4))
call cpu_time(start)
open(101,file=filename(1))  !Frequency trajectory input
open(102,file=filename(2))  !Coupling trajectory input
open(103,file=filename(3))  !Transition dipole trajectory
open(11,file=filename(4))   !IR output
do i=1,nchro
   do j=1,nchro
      if(i==j)then
        Uni(i,j)=one
      else
        Uni(i,j)=zero
      endif
   enddo
enddo
do s=1,nsnpsht
   read(101,*)
   read(102,*)
   read(103,*)
   !Construct the kappa matrix:
   do i=1,nchro
      read(101,*)kappa(i,i)
      do j=i+1,nchro
         read(102,*)k,l,kappa(k,l)
         kappa(l,k)=kappa(k,l)
      enddo
      read(103,*)m(s,i,:)
   enddo
   !Calculate U matrix where U=exp(im*kappa*dt):
   call jacobi(kappa,nchro,nchro,d,P,i) !Diagonalize kappa
   do i=1,nchro
      U(s,i,:)=zero
      U(s,i,i)=exp(it*d(i))
   enddo
   U(s,:,:)=matmul(matmul(P,U(s,:,:)),transpose(P))
   F(s,:,:)=Uni
enddo
close(101)
close(102)
close(103)
!Accumulate F matrix F(n)=Uni*U1*U2*U3*...*Un and calculate the <mFm> term:
mFm=zero
do k=0,tstep
   j=0
   do s=1,nsnpsht-k,interval
      j=j+1
      if(k>0)then
        F(s,:,:)=matmul(F(s,:,:),U(s+k-1,:,:))
      endif
      do i=1,3
         mFm(k)=mFm(k)+sum(matmul(m(s,:,i),F(s,:,:))*m(s+k,:,i))
      enddo
   enddo
   mFm(k)=mFm(k)/3.0d0/real(j)
enddo
!Apply Discrete Fourier Transform to obtain IR curve:
do i=1,f_tstep
   actualw=real(i)*dw  !actualw in rad/ps
   if(actualw>=(minfreq*conve).and.actualw<=(maxfreq*conve))then
     do dstep=0,f_tstep
        time=real(dstep)*dt  !in ps
        dintensity(dstep)=real(exp(-im*actualw*time)*mFm(dstep))*exp(-time/Tl)
     enddo
     inten(i)=0.0d0
     do dstep=1,f_tstep
        inten(i)=inten(i)+dintensity(dstep-1)+dintensity(dstep)
     enddo
     inten(i)=inten(i)*dt/2.0d0
     write(11,*)actualw/conve,inten(i)
   endif
enddo
close(11)
call cpu_time(finish)
write(*,*)finish-start,"seconds"

CONTAINS
SUBROUTINE jacobi(a,n,np,d,v,nrot)
INTEGER n,np,nrot,NMAX
double precision a(np,np),d(np),v(np,np)
PARAMETER (NMAX=500)
!Adopted from "Numerical Recipes in Fortran 90".
!Computes all eigenvalues and eigenvectors of a real symmetric matrix a, which is of size n
!by n, stored in a physical np by np array. On output, elements of a above the diagonal are
!destroyed. d returns the eigenvalues of a in its first n elements. v is a matrix with the same
!logical and physical dimensions as a, whose columns contain, on output, the normalized
!eigenvectors of a. nrot returns the number of Jacobi rotations that were required.
INTEGER i,ip,iq,j
double precision c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
do ip=1,n !Initialize to the identity matrix.
   do iq=1,n
      v(ip,iq)=0.0d0
   enddo
   v(ip,ip)=1.0d0
enddo
do ip=1,n
   b(ip)=a(ip,ip) !Initialize b and d to the diagonal of a.
   d(ip)=b(ip)
   z(ip)=0.d0
enddo
nrot=0
do i=1,50
   sm=0.d0
   do ip=1,n-1 !Sum off-diagonal elements.
      do iq=ip+1,n
         sm=sm+abs(a(ip,iq))
      enddo
   enddo
   if(sm.eq.0.d0)return !The normal return, which relies on quadratic convergence to machine underflow.
   if(i.lt.4)then
     tresh=0.2d0*sm/n**2
   else
     tresh=0.0d0
   endif
   do ip=1,n-1
      do iq=ip+1,n
         g=100.0d0*abs(a(ip,iq)) !After four sweeps, skip the rotation if the off-diagonal element is small.
         if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
           a(ip,iq)=0.0d0
         else if(abs(a(ip,iq)).gt.tresh)then
           h=d(iq)-d(ip)
           if(abs(h)+g.eq.abs(h))then
             t=a(ip,iq)/h !t = 1/(2θ)
           else
             theta=0.5d0*h/a(ip,iq)
             t=1.0d0/(abs(theta)+sqrt(1.0d0+theta**2))
             if(theta.lt.0.0d0)t=-t
           endif
           c=1.0d0/sqrt(1+t**2)
           s=t*c
           tau=s/(1.0d0+c)
           h=t*a(ip,iq)
           z(ip)=z(ip)-h
           z(iq)=z(iq)+h
           d(ip)=d(ip)-h
           d(iq)=d(iq)+h
           a(ip,iq)=0.0d0
           do j=1,ip-1
              g=a(j,ip)
              h=a(j,iq)
              a(j,ip)=g-s*(h+g*tau)
              a(j,iq)=h+s*(g-h*tau)
           enddo
           do j=ip+1,iq-1
              g=a(ip,j)
              h=a(j,iq)
              a(ip,j)=g-s*(h+g*tau)
              a(j,iq)=h+s*(g-h*tau)
           enddo
           do j=iq+1,n
              g=a(ip,j)
              h=a(iq,j)
              a(ip,j)=g-s*(h+g*tau)
              a(iq,j)=h+s*(g-h*tau)
           enddo
           do j=1,n
              g=v(j,ip)
              h=v(j,iq)
              v(j,ip)=g-s*(h+g*tau)
              v(j,iq)=h+s*(g-h*tau)
           enddo
           nrot=nrot+1
         endif
      enddo
   enddo
   do ip=1,n
      b(ip)=b(ip)+z(ip)
      d(ip)=b(ip) !Update d with the sum of tapq,
      z(ip)=0.0d0 !and reinitialize z.
enddo
enddo
return
END SUBROUTINE jacobi
END PROGRAM IR_Combined
