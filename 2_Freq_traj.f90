program Frequency
!This program is to calculate the frequencies from the electric field
!trajactory using the frequency maps.
implicit none
integer,parameter,dimension(10)::G22=[1716,1275,5658,339,-180,-385,-13,-853,-1026,-263]
integer,parameter,dimension(10)::C44=[1691,426,5444,-31,-62,-451,21,-1257,315,-60]
integer,parameter,dimension(10)::CC=[1638,1101,1521,62,-1952,-426,-210,0,0,0]
!Frequency map parameters. The first is the intercept.
integer,parameter,dimension(4)::nbase=[3,4,1,0]
!Each dimension represents the amount of GCUA respectively
integer,parameter::nsolute=1
integer,parameter::base=nsolute*sum(nbase)
integer,parameter::nsnpsht=200000
character(100)::filename(2)
integer::i,j
integer::idx,s,b,bt,m
double precision::EF(9)
real::start,finish
call cpu_time(start)
99 format('# f',i8)
call getarg(1,filename(1))
call getarg(2,filename(2))
open(11,file=filename(1))  !Input electric field traj
open(101,file=filename(2)) !Output frequency trajactory
read(11,*)
do i=1,nsnpsht
   write(101,99)i
   do j=1,base
      read(11,*)idx,EF,s,b,bt,m
      select case(bt)
      case(1) !G
        write(101,*)G22(1)+dot_product(G22(2:),EF)
      case(2) !C
        write(101,*)C44(1)+dot_product(C44(2:),EF)
        read(11,*)idx,EF,s,b,bt,m
        write(101,*)CC(1)+dot_product(CC(2:),EF)
      case(3) !U
        write(101,*)G22(1)+dot_product(G22(2:),EF)
        read(11,*)idx,EF,s,b,bt,m
        write(101,*)C44(1)+dot_product(C44(2:),EF)
        read(11,*)idx,EF,s,b,bt,m
        write(101,*)CC(1)+dot_product(CC(2:),EF)
      case(4) !A
        write(101,*)CC(1)+dot_product(CC(2:),EF)
      endselect
   enddo
enddo
close(11)
close(101)
call cpu_time(finish)
write(*,*)finish-start,"seconds"
end
