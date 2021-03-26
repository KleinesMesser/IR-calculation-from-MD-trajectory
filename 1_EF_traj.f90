PROGRAM ElecField
!This program is to calculate the electric field of C, O, N in intrinsic
!coordinates of RNA oligomer 3'-GCUCCGGC-5'
implicit none 
integer,parameter,dimension(4)::nbase=[3,4,1,0]
!Each dimension represents the amount of GCUA in one chain respectively.
integer,parameter,dimension(4)::natom=[34,31,30,33]
!Each dimension represents the amount of atoms in GCUA respectively.
integer,parameter,dimension(4)::nmode=[1,2,3,1]
!Each dimension represents the amount of local modes considered in GCUA respectively.
integer,parameter::nmodemax=3 !Maximum of nmode
integer,parameter::nbasesum=sum(nbase)
integer,parameter,dimension(nbasesum)::seq=[1,2,3,2,2,1,1,2]
!Sequence of each chain 1 2 3 4 represent G C U A respectively.
integer,parameter::nsolute=1,solute=dot_product(nbase,natom)-1,soluteall=nsolute*solute
!nsolute: # of oligonucleotide chains.
!solute: # of atoms in each chain. 5' side: PO2 -> H: -3+1=-2; 3' side: O3' -> O3'-HO3': +1.
!Totall # of atoms of each chain: -2+1=-1 compared to all middle residues.
integer,parameter::mode=nsolute*dot_product(nbase,nmode)
!Total # of local modes
integer,parameter::nion=(nbasesum-1)*nsolute,ion=1,nsol=6900,sol=3
!nion: total # of counter cation (asuming +1 charge). ion: # of atoms in one ion.
!nsol: total # of solvent molecules. sol: # of atoms in one solvent molecule (assuming water).
integer,parameter,dimension(1,3)::G=reshape([18,19,20],shape(G))
integer,parameter,dimension(2,3)::C=reshape([23,14,24,16,22,15],shape(C)) !C2,CC
integer,parameter,dimension(3,3)::UR=reshape([22,18,14,23,19,16,20,20,15],shape(UR)) !C2,C4,CC
integer,parameter,dimension(1,3)::A=reshape([18,17,19],shape(A)) !19 TBD
!Indicate the atom that determines the internal coordinate system.
!Each row represents one chromophore. Atom numbering is consistent with those shown in AMBER.
integer,parameter::nsnpsht=200000
double precision,parameter::chgconv=18.2223d0,disconv=0.52917721092d0**2
!convert amber intrinsic charge and distance to atomic unit (a.u).
double precision,parameter::cutoff=2.0d1 !Electric field cutoff in Angstrom.
double precision,parameter::boxtrans=2.0d0/sqrt(3.0d0)
!Dealing with the PBC of truncated octahedron
character(100)::filename(3)
integer::i,j,k,l,m,n,o,p,cnt
integer::atom(nsolute*nbasesum+1),chro(nsolute,nbasesum,nmodemax,3)
double precision::chg(soluteall+ion+sol)
double precision::rdnn(soluteall,3),rcnt(3),rsol(sol,3)
double precision::r(3),box,dis,chroctr(nsolute,nbasesum,nmodemax,3)
double precision::u(3,3) !unit vectors of the internal coordinate system
double precision::ele(nsolute,nbasesum,nmodemax,3,3),E(3,3)
!ele: electric field in original Cartesian coordinites, E: electric field in internal coordinates
real::start,finish
call getarg(1,filename(1))
call getarg(2,filename(2))
call getarg(3,filename(3))
call cpu_time(start)
open(11,file='1_Charge.dat')
open(12,file=filename(1))   !Coordinate xyz file
open(13,file=filename(2))   !Simulation box parameters (truncated octahedron)
open(101,file=filename(3))  !Output electric field trajectory
97 format(a8,9a14,4a3)
98 format(i8,9es14.6,4i3)
99 format(2x,3f8.3)

!Partial charge reading:
read(11,*)(chg(i),i=1,soluteall+ion+sol)
close(11)
chg=chg/chgconv
!Chromophore atom numbering stores in chro()
!Ending atom numbers of each residue store in atom()
cnt=0
chro=0
chro(:,1,:,:)=-2 !Dealing with 5' sitiuation where PO2 -> H
atom(1)=0
do i=1,nsolute
   do j=1,nbasesum
      cnt=cnt+1
      select case(seq(j))
      case(1) !G
        chro(i,j,1,:)=atom(cnt)+chro(i,j,1,:)+G(1,:)
      case(2) !C
        chro(i,j,1,:)=atom(cnt)+chro(i,j,1,:)+C(1,:)
        chro(i,j,2,:)=atom(cnt)+chro(i,j,2,:)+C(2,:)
      case(3) !U
        chro(i,j,1,:)=atom(cnt)+chro(i,j,1,:)+UR(1,:)
        chro(i,j,2,:)=atom(cnt)+chro(i,j,2,:)+UR(2,:)
        chro(i,j,3,:)=atom(cnt)+chro(i,j,3,:)+UR(3,:)
      case(4) !A
        chro(i,j,1,:)=atom(cnt)+chro(i,j,1,:)+A(1,:)
      end select
      select case(j)
      case(1) !5' side
        atom(cnt+1)=atom(cnt)+natom(seq(j))-2
      case(nbasesum) !3' side
        atom(cnt+1)=atom(cnt)+natom(seq(j))+1
      case default
        atom(cnt+1)=atom(cnt)+natom(seq(j))
      end select
   enddo
enddo

!Electric field calculation:
write(101,97)'#  index','Cx','Cy','Cz','O(C)x','O(C)y','O(C)z','N(H)x','N(H)y','N(H)z'&
             ,'s','b','bt','m' !s: solute #; b: base #; bt: base type (1G2C3U4A) m: mode #
do i=1,nsnpsht
   ele=0.0d0
   read(12,*)
   read(12,*)
   read(13,*)box
   box=box*boxtrans
   !nucleic acid part:
   do j=1,soluteall
      read(12,99)rdnn(j,:)
   enddo
   do j=1,nsolute
      do k=1,nbasesum
         do l=1,nmode(seq(k))
            chroctr(j,k,l,:)=&
              (rdnn(chro(j,k,l,1),:)+rdnn(chro(j,k,l,2),:)+rdnn(chro(j,k,l,3),:))/3.0d0
         enddo
      enddo
   enddo
   cnt=0
   do j=1,nsolute
      do k=1,nbasesum
         cnt=cnt+1
         do l=atom(cnt)+1,atom(cnt+1)
            do m=1,nsolute
               do n=1,nbasesum
                  do o=1,nmode(seq(n))
                     if(.not.(j==m.and.k==n))then
                       dis=length(chroctr(m,n,o,:)-rdnn(l,:))
                       if(dis<=cutoff)then
                         do p=1,3
                            r=rdnn(chro(m,n,o,p),:)-rdnn(l,:)
                            dis=length(r)
                            ele(m,n,o,p,:)=ele(m,n,o,p,:)+chg(l)*r/dis**3
                         enddo
                       endif
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !Ion part:
   do j=1,nion
      do k=1,ion
         read(12,99)rcnt
      enddo
      do k=1,nsolute
         do l=1,nbasesum
            do m=1,nmode(seq(l))
               r=chroctr(k,l,m,:)-rcnt
               call OCT_PBC(r,box)
               dis=length(r)
               if(dis<=cutoff)then
                 do n=1,3
                    do o=1,ion
                       r=rdnn(chro(k,l,m,n),:)-rcnt
                       call OCT_PBC(r,box)
                       dis=length(r)
                       ele(k,l,m,n,:)=ele(k,l,m,n,:)+chg(soluteall+o)*r/dis**3
                    enddo
                 enddo
               endif
            enddo
         enddo
      enddo
   enddo
   !solvent part:
   do j=1,nsol
      rcnt=0.0d0
      do k=1,sol
         read(12,99)rsol(k,:)
         rcnt=rcnt+rsol(k,:)
      enddo
      rcnt=rcnt/real(sol)
      do k=1,nsolute
         do l=1,nbasesum
            do m=1,nmode(seq(l))
               r=chroctr(k,l,m,:)-rcnt
               call OCT_PBC(r,box)
               dis=length(r)
               if(dis<=cutoff)then
                 do n=1,3
                    do o=1,sol
                       r=rdnn(chro(k,l,m,n),:)-rsol(o,:)
                       call OCT_PBC(r,box)
                       dis=length(r)
                       ele(k,l,m,n,:)=ele(k,l,m,n,:)+chg(soluteall+ion+o)*r/dis**3
                    enddo
                 enddo
               endif
            enddo
         enddo
      enddo
   enddo
   !Convert EF in the intrinsic coordinate:
   do j=1,nsolute
      do k=1,nbasesum
         do l=1,nmode(seq(k))
            call intr_vec(rdnn(chro(j,k,l,1),:),rdnn(chro(j,k,l,2),:),rdnn(chro(j,k,l,3),:),u,disconv)
            E=matmul(ele(j,k,l,:,:),u)
            write(101,98)i,E(1,:),E(2,:),E(3,:),j,k,seq(k),l
         enddo
      enddo
   enddo
enddo
close(12)
close(13)
close(101)
call cpu_time(finish)
write(*,*)finish-start,"seconds"

CONTAINS
SUBROUTINE OCT_PBC(u,box) !PBC of truncated octahedron
implicit none
double precision::u(3)
double precision,intent(in)::box
u=u-box*nint(u/box)
if(ABS(u(1))+ABS(u(2))+ABS(u(3))<0.75d0*box)goto 999
u=u-sign(0.5d0*box,u)
999 continue
END SUBROUTINE OCT_PBC
double precision FUNCTION length(u)
implicit none
double precision,intent(in)::u(3)
length=sqrt(dot_product(u,u))
END FUNCTION length
SUBROUTINE crs_product(u,v,crsp) !Cross product of two vectors
implicit none
double precision,intent(in)::u(3),v(3)
double precision,intent(out)::crsp(3)
crsp(1)=u(2)*v(3)-u(3)*v(2)
crsp(2)=u(3)*v(1)-u(1)*v(3)
crsp(3)=u(1)*v(2)-u(2)*v(1)
END SUBROUTINE crs_product
SUBROUTINE intr_vec(c1,c2,c3,u,cons) !Internal unit vector calculation
!c1,c2,c3 in C,O,N or C(connects to H),C,H order
implicit none
double precision,intent(in)::c1(3),c2(3),c3(3),cons
double precision,intent(out)::u(3,3)
u(:,2)=c2(:)-c1(:)
u(:,2)=cons*u(:,2)/length(u(:,2))  !y
u(:,1)=c3(:)-c1(:)
call crs_product(u(:,1),u(:,2),u(:,3))
u(:,3)=cons*u(:,3)/length(u(:,3))  !z
call crs_product(u(:,2),u(:,3),u(:,1))
u(:,1)=cons*u(:,1)/length(u(:,1))  !x
END SUBROUTINE
END PROGRAM ElecField
