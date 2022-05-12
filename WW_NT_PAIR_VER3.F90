!program for N(t) for WW hydrogen bond 
!considering pair of molecules hydrogen bonded
!Date: 28th March, 2015

!The code is written by Prabir Khatua

!Usage: f95 WW_NT_PAIR_VER3.F90 or gfortran WW_NT_PAIR_VER3.F90

!       ./a.out  (for interactive run and then enter the inputs interactivly 
!                 as the program will ask)
!       
!       ./a.out<wwnt.inp >wwnt.log &  
!                            (input.inp is the input file containing all the 
!                            required inputs that has to be prepared by the 
!                            user. For this, one needs to open the code or 
!                            better compile the code once in interactive mode 
!                            to know what set of inputs will be required and 
!                            prepare the input file. output.log will print 
!                            general information or error) 

!The code is written to analyse the dcd trajectories
!------------------------------------------------------------------------------
program WW_NT
integer,parameter:: mxatom=60000,maxinp=2000,mxframe=1000
integer,parameter:: mxhb=250000,mxtype=10,mxw=17000
real,dimension(mxatom):: x,y,z
real(kind=8),dimension(0:mxframe):: nt,av,avsq
real(kind=8),dimension(0:mxframe):: n
integer,dimension(mxw):: rm
integer,dimension(mxhb):: uu
integer,dimension(mxhb,mxframe):: h,hh
integer,dimension(maxinp):: atp
integer,dimension(mxtype):: site_res,num_res
integer:: funit,ounit,dummyi,res_type,total_atom
character(len=100),dimension(maxinp):: dcdfile
character(len=100):: outfile
character(len=100):: flgfile
character(len=4):: chr2,chr3
character(len=2):: chr4
logical:: there
real(kind=8):: dx,dy,dz,bx,by,bz,bar
data ounit,funit,inunit /17,18,19/
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the name of output file'
read(*,'(a)')outfile
write(*,*)outfile
open(ounit,file=outfile,status='unknown',form='formatted',iostat=ios)
!------------------------------------------------------------------------------
write(*,'(1x,a)') 'Enter number of residue types '
read(*,*) res_type
write(*,'(1x,a)')' Enter No. of resids of all types '
read(*,*) (num_res(i),i=1,res_type)
write(*,'(1x,a)')' Enter No. of sites in all resids '
read(*,*) (site_res(i),i=1,res_type)
!------------------------------------------------------------------------------
n_protein=num_res(1)*site_res(1)
n_ion=num_res(2)*site_res(2)
n_water = num_res(3)*site_res(3)
total_atom = n_protein+n_ion+n_water
write(*,*)'protein atom=====>',n_protein
write(*,*)'ions======>',n_ion
write(*,*)'water=====>',n_water
write(*,*)'total atoms=====>',total_atom

if(total_atom > mxatom)then
write(*,*)'ERROR: NO OF TOTAL ATOM EXCEEDS THE DECLARED DIMENSION'
write(*,*)'PLZ DO ADJUST THE SIZE OF MXATOM PARAMETER'
stop
endif

if(num_res(3) > mxw)then
write(*,*)'ERROR: NO OF WATER MOLECULE EXCEEDS THE DECLARED DIMENSION'
write(*,*)'PLZ DO ADJUST THE SIZE OF MXW PARAMETER'
stop
endif

nsit=site_res(3)
ni=n_protein+1
nf=total_atom-n_ion-nsit+1
write(*,*)'first water oxygen atom no=======>',ni
write(*,*)'last water oxygen atom no=======>',nf
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the skip value'
read(*,*)nskip
write(*,*)'nskip=======>',nskip
write(*,'(1x,a)')'Type first & last resd no of the selected protein surface'
read(*,*) ifirst,ilast
write(*,*)'ifirst,ilast====>',ifirst,ilast
write(*,'(1x,a)')'Enter the time gap betn two frames in ps'
read(*,*)dtstep
dtstep=dtstep*nskip
write(*,*)'dtstep=======>',dtstep
write(*,'(1x,a)')'Enter the cut-off radius for D-A bond'
read(*,*)rl
write(*,*)'cut-off radius for D-A=====>',rl
write(*,'(1x,a)')'Enter the cut-off angle'
write(*,*)'Enter the cut-off radius for A-H bond'
read(*,*)rh
write(*,*)'cut-off radius for A-H=====>',rh
read(*,*)al
write(*,*)'Cut-off angle======>',al
write(*,'(1x,a)')'Enter the inner and outer radius of selected region'
read(*,*)rin,rout
write(*,*)'rin,rout=========>',rin,rout
write(*,'(1x,a)')'Enter the no of block'
read(*,*)nset
write(*,*)'nset======>',nset
rin=rin*rin
rout=rout*rout
rh=rh*rh
rl=rl*rl
!-----------------------------------------------------------------------------
write(*,'(1x,a)')'Enter name of a reference file in PDB format'
read(*,'(a)')flgfile
write(*,*)flgfile
open(funit,file=flgfile,status='old',form='formatted',iostat=ios)

npair=0
do i=1,n_protein
read(funit,'(12x,a4,1x,a4,1x,i4,50x,a2)')chr2,chr3,nres,chr4
if(nres >= ifirst.and.nres <= ilast)then
if(chr4 /= ' H')then
npair=npair+1

if(npair > maxinp)then
write(*,*)'ERROR: NO OF HEAVY ATOM EXCEEDS THE DECLARED DIMENSION'
write(*,*)'PLZ DO RESET THE MAXINP PARAMETER'
stop
endif

atp(npair)=i
endif
endif
enddo

write(*,*)'# of protein heavy atoms in selected surface===>',npair
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the number of trajectory files you want to analyze.'
read(*,*) ninput
write(*,*)'ninput==========>',ninput

if(ninput > maxinp)then
write(*,*)'ERROR: NO OF INPUT FILE EXCEEDS THE DECLARED DIMENSION'
write(*,*)'PLZ DO ADJUST THE SIZE OF MAXINP PARAMETER'
stop
endif

do inp=1,ninput
write(*,'(1x,a)')'Enter the name of trajectory file',inp
read(*,'(a)')dcdfile(inp)
write(*,*)dcdfile(inp)
inquire(file=dcdfile(inp),exist=there)
if(.not.there)then
write(*,'(a)')'ERROR: CONFIGURATION FILE NOT FOUND'
stop
endif
enddo
!------------------------------------------------------------------------------
av(:)=0.0
avsq(:)=0.0
ng=ninput/nset
do igroup=1,ninput,ng
nu=0
nr=0
h(:,:)=0
hh(:,:)=0
nhb=0
do inp=igroup,igroup+ng-1
write(*,'(1x,a,i5,a)')'Enter name of trajectory file',inp
write(*,*)dcdfile(inp)
open(inunit,file=dcdfile(inp),status='old',form='unformatted')
read(inunit)dummyc,nframes,(dummyi,i=1,8),dummyr,(dummyi,i=1,9)
read(inunit)dummyi,dummyr
write(*,*)"dummyi,dummyr====>",dummyi,dummyr
read(inunit)natom
write(*,*)"natom====>",natom
write(*,*)"nframes===>",nframes
do ii=1,nframes
read(inunit)dx,bx,dy,by,bz,dz
read(inunit)(x(j),j=1,natom)
read(inunit)(y(j),j=1,natom)
read(inunit)(z(j),j=1,natom)
nr=nr+1
if(mod(nr,nskip) /= 0)cycle
nu=nu+1
if(nu > mxframe)exit 
!------------------------------------------------------------------------------
call CENTER(npair,atp,n_protein,ni,nf,nsit,dx,dy,dz,x,y,z)
call MIN_DIST(npair,atp,ni,nf,nsit,rin,rout,x,y,z,rm)
call HBOND(nu,nhb,ni,nf,nsit,rm,mxhb,x,y,z,al,rl,rh,uu,h,hh)
enddo
enddo
!------------------------------------------------------------------------------
nt(:)=0.d0
n(:)=0.d0

nu=nu/2
do i=1,nhb
do j=1,nu-1
if(h(i,j) == 0)cycle
do k=j+1,nu
m=k-j
nt(m)=nt(m)+h(i,j)*(1.0-h(i,k))*hh(i,k)
n(m)=n(m)+1
enddo
enddo
enddo

write(*,*)'TOTAL DIFFERENT TYPE OF HBOND=========>',nhb
!------------------------------------------------------------------------------
do i=0,nu
if(n(i) == 0)n(i)=1
if(i == 0)then
nt(i)=0.0
else
nt(i)=nt(i)/n(i)
endif
av(i)=av(i)+nt(i)
avsq(i)=avsq(i)+nt(i)*nt(i)
enddo
enddo

bar=0.0
do i=0,nu
av(i)=av(i)/real(nset)
avsq(i)=avsq(i)/real(nset)
time=i*dtstep
dif=avsq(i)-av(i)*av(i)
if(dif < 0.0)then
write(*,*)'WARNING: SOMETHING IS WRONG'
write(*,*)'STANDARD DEVIATION IS NEGATIVE'
dif=0.0
endif
xb=0.5*sqrt(dif)
write(ounit,'(3f12.6)')time,av(i),xb
xb=(xb*100.0)/av(i)
if(av(i) == 0.0)xb=0.0
bar=bar+xb
enddo

write(*,*)bar
bar=bar/real(nu+1)
write(ounit,*)'Average error bar in % for this calculation is = ',bar
write(ounit,*)'No of blocks used in this calculation = ',nset

end program WW_NT
!------------------------------------------------------------------------------
subroutine HBOND(nu,nhb,ni,nf,ns,rm,mxhb,x,y,z,al,rl,rh,uu,h,hh)
integer,dimension(mxhb,*),intent(out):: h,hh
integer,dimension(*),intent(inout):: uu
integer,intent(inout):: nhb
integer,dimension(*),intent(in):: rm
integer,intent(in):: ni,nf,ns,nu
real,dimension(*),intent(in):: x,y,z
real,intent(in):: rl,al,rh
integer:: wpair
logical:: NEW,HYD


if(nu == 1)NEW=.TRUE.
np=0
n1=0
do i=ni,nf-ns,ns
n1=n1+1
if(rm(n1) == 1)then
HYD=.TRUE.
else
HYD=.FALSE.
endif
n2=n1
jloop: do j=i+ns,nf,ns
n2=n2+1
np=np+1
if(.not.HYD.and.rm(n2) == 0)cycle jloop
rx=x(i)-x(j)
ry=y(i)-y(j)
rz=z(i)-z(j)

r=rx*rx+ry*ry+rz*rz
if(r > rl)cycle jloop
if(nu /=1)call CHECK(np,nhb,uu,wpair,NEW)
if(NEW)then
nhb=nhb+1

if(nhb > mxhb)then
write(*,*)'ERROR: TYPE OF HBOND EXCEEDS THE DECLARED DIMENSION'
write(*,*)'PLZ DO RESET MXHB PARAMETER'
stop
endif

uu(nhb)=np
hh(nhb,nu)=1
else
hh(wpair,nu)=1
endif

r=sqrt(r)
rx=rx/r
ry=ry/r
rz=rz/r

do jj=j+1,j+2     !jth water as donor
fx=x(jj)-x(i)
fy=y(jj)-y(i)
fz=z(jj)-z(i)

f=fx*fx+fy*fy+fz*fz
if(f > rh)cycle
vx=x(jj)-x(j)
vy=y(jj)-y(j)
vz=z(jj)-z(j)

v=vx*vx+vy*vy+vz*vz
v=sqrt(v)

vx=vx/v
vy=vy/v
vz=vz/v

call DOT_PROD(rx,ry,rz,vx,vy,vz,a)
if(a > al)cycle
if(nu /=1)call CHECK(np,nhb,uu,wpair,NEW)
if(NEW)then
nhb=nhb+1

if(nhb > mxhb)then
write(*,*)'ERROR: TYPE OF HBOND EXCEEDS THE DECLARED DIMENSION'
write(*,*)'PLZ DO RESET MXHB PARAMETER'
stop
endif

uu(nhb)=np
h(nhb,nu)=1
else
h(wpair,nu)=1
endif
cycle jloop
enddo

do ii=i+1,i+2     !ith water as donor
fx=x(ii)-x(j)
fy=y(ii)-y(j)
fz=z(ii)-z(j)

f=fx*fx+fy*fy+fz*fz
if(f > rh)cycle

vx=x(ii)-x(i)
vy=y(ii)-y(i)
vz=z(ii)-z(i)

v=vx*vx+vy*vy+vz*vz
v=sqrt(v)
vx=vx/v
vy=vy/v
vz=vz/v

call DOT_PROD(-rx,-ry,-rz,vx,vy,vz,a)
if(a > al)cycle
nn=nn+1
if(nu /=1)call CHECK(np,nhb,uu,wpair,NEW)
if(NEW)then
nhb=nhb+1

if(nhb > mxhb)then
write(*,*)'ERROR: TYPE OF HBOND EXCEEDS THE DECLARED DIMENSION'
write(*,*)'PLZ DO RESET MXHB PARAMETER'
stop
endif

uu(nhb)=np
h(nhb,nu)=1
else
h(wpair,nu)=1
endif
cycle jloop
enddo
enddo jloop
enddo 

end subroutine HBOND
!------------------------------------------------------------------------------
subroutine CHECK(np,nhb,uu,wpair,NEW)
integer,intent(in):: np,nhb
integer,dimension(*),intent(in):: uu
integer,intent(out):: wpair
logical,intent(out):: NEW
integer:: nn

nn=0
do i=1,nhb
if(uu(i) == np)then
wpair=i
nn=1
exit
endif
enddo

if(nn == 0)then
NEW=.TRUE.
else
NEW=.FALSE.
endif

end subroutine CHECK
!------------------------------------------------------------------------------
subroutine DOT_PROD(x1,y1,z1,x2,y2,z2,a)
real,intent(in):: x1,y1,z1,x2,y2,z2
real,intent(out):: a

pi=4.0*atan(1.0)
a=x1*x2+y1*y2+z1*z2
a=acos(a)*180.0/pi

end subroutine DOT_PROD
!------------------------------------------------------------------------------
subroutine MIN_DIST(np,atp,ni,nf,ns,rin,rout,x,y,z,rm)
real,dimension(*),intent(in):: x,y,z
real,intent(in):: rin,rout
integer,dimension(*),intent(in):: atp
integer,intent(in):: np,ni,nf,ns
integer,dimension(*),intent(out):: rm

k=0
do i=ni,nf,ns
k=k+1
rm(k)=0
do j=1,np
jj=atp(j)
rx=x(jj)-x(i)
ry=y(jj)-y(i)
rz=z(jj)-z(i)

r=rx*rx+ry*ry+rz*rz
if(r >= rin.and.r <= rout)then
rm(k)=1
exit
endif
enddo
enddo

end subroutine MIN_DIST
!------------------------------------------------------------------------------
subroutine CENTER(np,atp,npr,ni,nf,ns,dx,dy,dz,x,y,z)
real,dimension(*),intent(inout):: x,y,z
integer,dimension(*),intent(in):: atp
integer,intent(in):: np,npr,ni,nf,ns
real(kind=8),intent(in):: dx,dy,dz

cx=0.0
cy=0.0
cz=0.0
do i=1,np
ii=atp(i)
cx=cx+x(ii)
cy=cy+y(ii)
cz=cz+z(ii)
enddo

cx=cx/real(np)
cy=cy/real(np)
cz=cz/real(np)

do i=1,npr
x(i)=x(i)-cx
y(i)=y(i)-cy
z(i)=z(i)-cz
enddo

do i=ni,nf,ns
rx=x(i)-cx
ry=y(i)-cy
rz=z(i)-cz

x(i)=rx-dx*anint(rx/dx)
y(i)=ry-dy*anint(ry/dy)
z(i)=rz-dz*anint(rz/dz)
do j=i+1,i+ns-1
vx=x(j)-cx
vy=y(j)-cy
vz=z(j)-cz

x(j)=vx-dx*anint(rx/dx)
y(j)=vy-dy*anint(ry/dy)
z(j)=vz-dz*anint(rz/dz)
enddo
enddo

end subroutine CENTER
!--------------------------------end of the program----------------------------
