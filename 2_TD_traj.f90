program Transition_Dipole
!This program is to calculate the transition dipole (TD) from an xyz trajactory file.
!The TD is assumed as a fixed quantity relative to the internal coordinate system.
implicit none
integer,parameter::nsolute=1
integer,parameter,dimension(5)::nbase=nsolute*[3,4,0,1,0]
!In the order of G, C, T, U, A
integer,parameter,dimension(5)::natom=[15,12,14,11,14]
!In the order of G, C, T, U, A. Only include base part.
integer,parameter::natommax=15 !Maximum of natom
integer,parameter::nbasesum=sum(nbase),nbasesumsingle=nbasesum/nsolute
integer,parameter,dimension(nbasesum)::&
seq=[1,2,4,2,2,1,1,2]
!1 means G; 2 means C; 4 means U; 5 means A
integer,parameter::nion=nbasesum-nsolute,ion=1,nsol=6900,sol=3
integer,parameter::nsnpsht=200000
integer,parameter::bb=12,ab=7 !Number of atoms befor and after base
integer,parameter::skipafter=nion*ion+nsol*sol
character::element
character(100)::filename(2)
integer::i,j,k,l,ntemp
double precision::coor(natommax,3)
double precision::TD(3)
real::start,finish
call getarg(1,filename(1))
call getarg(2,filename(2))
call cpu_time(start)
open(11,file=filename(1))  !Input coordiate xyz file
open(101,file=filename(2)) !Output transition dipole traj
99 format('# f',i8)
98 format(3es14.6)
do i=1,nsnpsht
   ntemp=0
   read(11,*)
   read(11,*)
   write(101,99)i
   do j=1,nsolute
      do k=1,nbasesumsingle
         ntemp=ntemp+1
         if(k==1)then  !5' side
           do l=1,bb-2
              read(11,*)
           enddo
         else
           do l=1,bb
              read(11,*)
           enddo
         endif
         select case(seq(ntemp))
         case(1)  !G
           do l=1,natom(1)
              read(11,*)element,coor(l,:)
           enddo
           call Tran_DipCO(coor(6,:),coor(7,:),coor(8,:),TD)
           write(101,98)TD
         case(2)  !C
           do l=1,natom(2)
              read(11,*)element,coor(l,:)
           enddo
           call Tran_DipCO(coor(11,:),coor(12,:),coor(10,:),TD)
           write(101,98)TD
           call Tran_DipCC(coor(2,:),coor(4,:),coor(1,:),TD)
           write(101,98)TD
         case(3)  !T
           do l=1,natom(3)
              read(11,*)element,coor(l,:)
           enddo
           call Tran_DipCO(coor(13,:),coor(14,:),coor(11,:),TD)
           write(101,98)TD
           call Tran_DipCO(coor(9,:),coor(10,:),coor(11,:),TD)
           write(101,98)TD
           call Tran_DipCC(coor(2,:),coor(4,:),coor(1,:),TD)
           write(101,98)TD
         case(4)  !U
           do l=1,natom(4)
              read(11,*)element,coor(l,:)
           enddo
           call Tran_DipCO(coor(10,:),coor(11,:),coor(8,:),TD)
           write(101,98)TD
           call Tran_DipCO(coor(6,:),coor(7,:),coor(8,:),TD)
           write(101,98)TD
           call Tran_DipCC(coor(2,:),coor(4,:),coor(1,:),TD)
           write(101,98)TD
         case(5)  !A
           do l=1,natom(5) 
              read(11,*)element,coor(l,:)
           enddo
           call Tran_DipCC(coor(6,:),coor(5,:),coor(7,:),TD)
           write(101,98)TD
         end select
         if(k==nbasesumsingle)then !3' side
           do l=1,ab+1
              read(11,*)
           enddo
         else
           do l=1,ab
              read(11,*)
           enddo
         endif
      enddo
   enddo
   do j=1,skipafter
      read(11,*)
   enddo
enddo
close(11)
close(101)
call cpu_time(finish)
write(*,*)finish-start,'seconds'

CONTAINS
double precision FUNCTION length(u)
implicit none
double precision,intent(in)::u(3)
length=sqrt(dot_product(u,u))
END FUNCTION length
SUBROUTINE crs_product(u,v,crsp)
implicit none
double precision,intent(in)::u(3),v(3)
double precision,intent(out)::crsp(3)
crsp(1)=u(2)*v(3)-u(3)*v(2)
crsp(2)=u(3)*v(1)-u(1)*v(3)
crsp(3)=u(1)*v(2)-u(2)*v(1)
END SUBROUTINE crs_product
SUBROUTINE Tran_DipCO(rC,rO,rN,rTD)
implicit none
double precision,parameter::Angle=3.65d0,Mag=2.57d0
double precision,parameter::PI=acos(-1.0d0),DE=PI/180.0d0
double precision,parameter::sin2A=sin(2.0d0*Angle*DE),sinA=sin(Angle*DE)
double precision,intent(in)::rC(3),rO(3),rN(3)
double precision,intent(out)::rTD(3)
double precision::vecCO(3),vecCN(3),nCO(3),nCN(3),disCO,disCN
double precision::cosT,sinT,x
vecCO=rO-rC
vecCN=rN-rC
disCO=length(vecCO)
disCN=length(vecCN)
nCO=vecCO/disCO
nCN=vecCN/disCN
cosT=dot_product(nCO,nCN)
sinT=sqrt(1.0d0-cosT**2)
x=disCO*(sinT*sin2A-2.0d0*cosT*(sinA**2))/disCN/2.0d0/(sinT**2-sinA**2)
rTD=-vecCO+x*vecCN
rTD=Mag*rTD/length(rTD)
END SUBROUTINE Tran_DipCO
SUBROUTINE Tran_DipCC(rC,rO,rN,rTD)
implicit none
double precision,parameter::Angle=0.0d0,Mag=1.44d0
double precision,parameter::PI=acos(-1.0d0),DE=PI/180.0d0
double precision,parameter::sin2A=sin(2.0d0*Angle*DE),sinA=sin(Angle*DE)
double precision,intent(in)::rC(3),rO(3),rN(3)
double precision,intent(out)::rTD(3)
double precision::vecCO(3),vecCN(3),nCO(3),nCN(3),disCO,disCN
double precision::cosT,sinT,x
vecCO=rO-rC
vecCN=rN-rC
disCO=length(vecCO)
disCN=length(vecCN)
nCO=vecCO/disCO
nCN=vecCN/disCN
cosT=dot_product(nCO,nCN)
sinT=sqrt(1.0d0-cosT**2)
x=disCO*(sinT*sin2A-2.0d0*cosT*(sinA**2))/disCN/2.0d0/(sinT**2-sinA**2)
rTD=-vecCO+x*vecCN
rTD=Mag*rTD/length(rTD)
END SUBROUTINE Tran_DipCC
END
