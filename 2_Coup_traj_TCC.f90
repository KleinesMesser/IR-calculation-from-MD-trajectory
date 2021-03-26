program Transition_Charge_Coupling
!This program is to calculate the coupling by the TCC scheme.
implicit none
integer,parameter::nsolute=1
integer,parameter,dimension(5)::nbase=nsolute*[3,4,0,1,0]
!In the order of GCTUA.
integer,parameter,dimension(5)::natom=[15,12,14,11,14]     
!In the order of GCTUA. Only include base part.
integer,parameter,dimension(5)::nmode=[1,2,3,3,1]
!In the order of GCTUA.
integer,parameter::nbasemax=nsolute*4 !Maximum of nbase
integer,parameter::natommax=15 !Maximum of natom
integer,parameter::nmodemax=3  !Maximum of nmode
integer,parameter::nbasesum=sum(nbase),nbasesumsingle=nbasesum/nsolute
integer,parameter,dimension(nbasesum)::&
seq=[1,2,4,2,2,1,1,2]
!1 means G; 2 means C; 3 means T; 4 means U; 5 means A
integer,parameter::nion=nbasesum-nsolute,ion=1,nsol=6900,sol=3
integer,parameter::nsnpsht=200000
integer,parameter::bb=12,ab=7 !Amount of atoms before and after base
integer,parameter::skipafter=nion*ion+nsol*sol
!Coupling does not concern ion and solvent. So those parts are all skipped.
double precision,parameter::C12=-10.57d0
double precision,parameter::T12=17.27d0,T13=-8.16d0,T23=6.26d0
double precision,parameter::U12=18.55d0,U13=-9.35d0,U23=16.20d0
!Intra-base coupling
character::element
character(100)::filename(2)
integer::s,a,b,i,j,k,l,cnt(5),cnt2(5),ntemp
integer::md(nbasesum) !Stores the mode numbering
double precision::ampli
double precision::coor(5,nbasemax,natommax,3) !First dimension: 1G, 2C, 3T, 4U, 5A
double precision::u(3,3)
double precision::q(5,nmodemax,natommax)  !First dimension: 1G, 2C, 3T, 4U, 5A
double precision::dq(5,nmodemax,natommax) !First dimension: 1G, 2C, 3T, 4U, 5A
double precision::v(5,nbasemax,nmodemax,natommax,3)
double precision::vtemp(5,nbasemax,nmodemax,natommax,3)
!First dimension: 1G, 2C, 3T, 4U, 5A
double precision::r(natommax,natommax,3)
real::start,finish
call getarg(1,filename(1))
call getarg(2,filename(2))
call cpu_time(start)
99 format('# f',i8)
98 format(2i3,es14.6)
open(11,file=filename(2))  !Output coupling trajactory
open(101,file=filename(1)) !Input coordiate xyz file
!########## Read in q, dq and v: ##################
!G
open(102,file="2_TCC_para/G_NMAD_q2.dat")
open(103,file="2_TCC_para/G_NMAD_dq2.dat")
open(104,file="2_TCC_para/G_NMAD_vinter2.dat")
!C
open(105,file="2_TCC_para/C_Mole2_q2.dat")
open(106,file="2_TCC_para/C_Mole2_dq2.dat")
open(107,file="2_TCC_para/C_Mole2_vinter2.dat")
!C_CC
open(108,file="2_TCC_para/C_Mole3_q.dat")
open(109,file="2_TCC_para/C_Mole3_dq.dat")
open(110,file="2_TCC_para/C_Mole3_vinter.dat")
!T_2
open(111,file="2_TCC_para/T2_Mole2_q2.dat")
open(112,file="2_TCC_para/T2_Mole2_dq2.dat")
open(113,file="2_TCC_para/T2_Mole2_vinter2.dat")
!T_4
open(114,file="2_TCC_para/T4_NMAD_q2.dat")
open(115,file="2_TCC_para/T4_NMAD_dq2.dat")
open(116,file="2_TCC_para/T4_NMAD_vinter2.dat")
!T_CC
open(117,file="2_TCC_para/T_Mole3_q.dat")
open(118,file="2_TCC_para/T_Mole3_dq.dat")
open(119,file="2_TCC_para/T_Mole3_vinter.dat")
!U_2
open(120,file="2_TCC_para/U2_Mole2_q2.dat")
open(121,file="2_TCC_para/U2_Mole2_dq2.dat")
open(122,file="2_TCC_para/U2_Mole2_vinter2.dat")
!U_4
open(123,file="2_TCC_para/U4_NMAD_q2.dat")
open(124,file="2_TCC_para/U4_NMAD_dq2.dat")
open(125,file="2_TCC_para/U4_NMAD_vinter2.dat")
!U_CC
open(126,file="2_TCC_para/U_Mole3_q.dat")
open(127,file="2_TCC_para/U_Mole3_dq.dat")
open(128,file="2_TCC_para/U_Mole3_vinter.dat")
!A_CC
open(129,file="2_TCC_para/A_Mole3_q.dat")
open(130,file="2_TCC_para/A_Mole3_dq.dat")
open(131,file="2_TCC_para/A_Mole3_vinter.dat")

ntemp=99
do i=1,5 !1G, 2C, 3T, 4U, 5A
   do j=1,nmode(i)
      ntemp=ntemp+3
      read(ntemp+2,*)ampli
      do k=1,natom(i)
         read(ntemp,*)q(i,j,k)
         read(ntemp+1,*)dq(i,j,k)
         read(ntemp+2,*)v(i,1,j,k,:)
      enddo
      !The following "if" is important: adjust the sign of CO displacement to make
      !each sign of C and O(C) the same (The most dominant component, in this case y).
      !Alternatively, the sign in the corresponding dq and v files can be changed.
      if(ntemp==105.or.ntemp==108.or.ntemp==111.or.ntemp==117.or.ntemp==120.or.ntemp==126.or.ntemp==129)then
        v(i,1,j,:,:)=-v(i,1,j,:,:)*ampli
        dq(i,j,:)=-dq(i,j,:)
      else
        v(i,1,j,:,:)=v(i,1,j,:,:)*ampli
      endif
      do k=2,nbase(i)
         v(i,k,j,:,:)=v(i,1,j,:,:)
      enddo
      close(ntemp)
      close(ntemp+1)
      close(ntemp+2)
   enddo
enddo
!########## Starting mode number of each base stores in md(): ###########
md(1)=0
do i=1,nbasesum-1
   md(i+1)=md(i)+nmode(seq(i))
enddo
!############ Read in coordinates and TCC calculations: ################
do s=1,nsnpsht
   cnt=0
   ntemp=0
   read(101,*)
   read(101,*)
   do a=1,nsolute
      do i=1,nbasesumsingle  !5' side
         ntemp=ntemp+1
         cnt(seq(ntemp))=cnt(seq(ntemp))+1
         if(i==1)then
           do j=1,bb-2
              read(101,*)
           enddo
         else
           do j=1,bb
              read(101,*)
           enddo
         endif 
         select case(seq(ntemp))
         case(1)  !G
           do j=1,natom(1)
              read(101,*)element,coor(1,cnt(1),j,:)
           enddo
           call intr_vec(coor(1,cnt(1),6,:),coor(1,cnt(1),7,:),coor(1,cnt(1),8,:),u,1.0d0)
           vtemp(1,cnt(1),1,:,:)=matmul(v(1,cnt(1),1,:,:),transpose(u))
         case(2)  !C
           do j=1,natom(2)
              read(101,*)element,coor(2,cnt(2),j,:)
           enddo
           call intr_vec(coor(2,cnt(2),11,:),coor(2,cnt(2),12,:),coor(2,cnt(2),10,:),u,1.0d0)
           vtemp(2,cnt(2),1,:,:)=matmul(v(2,cnt(2),1,:,:),transpose(u))
           call intr_vec(coor(2,cnt(2),2,:),coor(2,cnt(2),4,:),coor(2,cnt(2),1,:),u,1.0d0)
           vtemp(2,cnt(2),2,:,:)=matmul(v(2,cnt(2),2,:,:),transpose(u))
         case(3)  !T
           do j=1,natom(3)
              read(101,*)element,coor(3,cnt(3),j,:)
           enddo
           call intr_vec(coor(3,cnt(3),13,:),coor(3,cnt(3),14,:),coor(3,cnt(3),11,:),u,1.0d0)
           vtemp(3,cnt(3),1,:,:)=matmul(v(3,cnt(3),1,:,:),transpose(u))
           call intr_vec(coor(3,cnt(3),9,:),coor(3,cnt(3),10,:),coor(3,cnt(3),11,:),u,1.0d0)
           vtemp(3,cnt(3),2,:,:)=matmul(v(3,cnt(3),2,:,:),transpose(u))
           call intr_vec(coor(3,cnt(3),2,:),coor(3,cnt(3),4,:),coor(3,cnt(3),1,:),u,1.0d0)
           vtemp(3,cnt(3),3,:,:)=matmul(v(3,cnt(3),3,:,:),transpose(u))
         case(4)  !U
           do j=1,natom(4)
              read(101,*)element,coor(4,cnt(4),j,:)
           enddo
           call intr_vec(coor(4,cnt(4),10,:),coor(4,cnt(4),11,:),coor(4,cnt(4),8,:),u,1.0d0)
           vtemp(4,cnt(4),1,:,:)=matmul(v(4,cnt(4),1,:,:),transpose(u))
           call intr_vec(coor(4,cnt(4),6,:),coor(4,cnt(4),7,:),coor(4,cnt(4),8,:),u,1.0d0)
           vtemp(4,cnt(4),2,:,:)=matmul(v(4,cnt(4),2,:,:),transpose(u))
           call intr_vec(coor(4,cnt(4),2,:),coor(4,cnt(4),4,:),coor(4,cnt(4),1,:),u,1.0d0)
           vtemp(4,cnt(4),3,:,:)=matmul(v(4,cnt(4),3,:,:),transpose(u))
         case(5)  !A
           do j=1,natom(5)
              read(101,*)element,coor(5,cnt(5),j,:) 
           enddo
           call intr_vec(coor(5,cnt(5),6,:),coor(5,cnt(5),5,:),coor(5,cnt(5),7,:),u,1.0d0)
           vtemp(5,cnt(5),1,:,:)=matmul(v(5,cnt(5),1,:,:),transpose(u))
         end select
         if(i==nbasesumsingle)then !3' side
           do j=1,ab+1
              read(101,*)
           enddo
         else
           do j=1,ab
              read(101,*)
           enddo
         endif
      enddo
   enddo
   do a=1,skipafter
      read(101,*)
   enddo
   write(11,99)s
   !#################  Inter-base couplings(TCC):  ####################
   cnt=0
   do a=1,nbasesum
      cnt(seq(a))=cnt(seq(a))+1
      cnt2=0
      do b=1,nbasesum
         cnt2(seq(b))=cnt2(seq(b))+1
         if(a<b)then
           do k=1,natom(seq(a))
              do l=1,natom(seq(b))
                 r(k,l,:)=coor(seq(a),cnt(seq(a)),k,:)-coor(seq(b),cnt2(seq(b)),l,:)
              enddo
           enddo
           do i=1,nmode(seq(a))
              do j=1,nmode(seq(b))
                 write(11,98)md(a)+i,md(b)+j,&
                 TCC(natom(seq(a)),natom(seq(b)),&
                 q(seq(a),i,:natom(seq(a))),dq(seq(a),i,:natom(seq(a))),&
                 q(seq(b),j,:natom(seq(b))),dq(seq(b),j,:natom(seq(b))),&
                 r(:natom(seq(a)),:natom(seq(b)),:),&
                 vtemp(seq(a),cnt(seq(a)),i,:natom(seq(a)),:),&
                 vtemp(seq(b),cnt2(seq(b)),j,:natom(seq(b)),:))
              enddo
           enddo
         endif
      enddo
   enddo
   !#################  Intra-base couplings:  ####################
   do a=1,nbasesum
      select case(seq(a))
      !Only CTU contain intr-base couplings
      case(2) !C
        write(11,98)md(a)+1,md(a)+2,C12
      case(3) !T
        write(11,98)md(a)+1,md(a)+2,T12
        write(11,98)md(a)+1,md(a)+3,T13
        write(11,98)md(a)+2,md(a)+3,T23
      case(4) !U
        write(11,98)md(a)+1,md(a)+2,U12
        write(11,98)md(a)+1,md(a)+3,U13
        write(11,98)md(a)+2,md(a)+3,U23
      end select
   enddo
enddo
close(11)
close(101)
call cpu_time(finish)
write(*,*)finish-start,"seconds"

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
SUBROUTINE intr_vec(c1,c2,c3,u,cons)
implicit none
double precision,intent(in)::c1(3),c2(3),c3(3),cons
double precision,intent(out)::u(3,3)
u(:,2)=c2-c1
u(:,2)=cons*u(:,2)/length(u(:,2))  !y
u(:,1)=c3-c1
call crs_product(u(:,1),u(:,2),u(:,3))
u(:,3)=cons*u(:,3)/length(u(:,3))  !z
call crs_product(u(:,2),u(:,3),u(:,1))
u(:,1)=cons*u(:,1)/length(u(:,1))  !x
END SUBROUTINE
double precision FUNCTION TCC(n,m,q1,dq1,q2,dq2,r,v1,v2)
implicit none
double precision,parameter::PI=acos(-1.0d0),eps0=8.85418782d-22,h=6.626070d-34,c=2.99792458d10
double precision,parameter::uchg=1.602177d-19
!permiity in C^2/J/A; plank in Js; speed of light in cm/s
double precision,parameter::conv=4.0d0*PI*eps0*h*c/uchg**2  !Give out coupling in cm^-1
integer::n,m
integer::i,j
double precision::q1(n),dq1(n),q2(m),dq2(m) 
double precision::r(n,m,3),v1(n,3),v2(m,3)
double precision::disr
TCC=0.0d0
do i=1,n
   do j=1,m
      disr=length(r(i,j,:))
      TCC=TCC+dq1(i)*dq2(j)/disr&
          -3.0d0*q1(i)*q2(j)*dot_product(v1(i,:),r(i,j,:))*dot_product(v2(j,:),r(i,j,:))/disr**5&
          -(-dq1(i)*q2(j)*dot_product(v2(j,:),r(i,j,:))+q1(i)*dq2(j)*dot_product(v1(i,:),r(i,j,:))&
          -q1(i)*q2(j)*dot_product(v1(i,:),v2(j,:)))/disr**3
   enddo
enddo
TCC=TCC/conv
END FUNCTION TCC
END
