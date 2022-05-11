!code for reorientational correlation function
!Prabir Khatua
!date: 5th January, 2015
!modified on: 1st July, 2015

!The code is written by Prabir Khatua

!Usage: f95 CMEU_2D_FINAL.F90 or gfortran CMEU_2D_FINAL.F90

!       ./a.out  (for interactive run and then enter the inputs interactivly 
!                 as the program will ask)
!       
!       ./a.out<cmeu.inp >cmeu.log &  
!                            (input.inp is the input file containing all the 
!                            required inputs that has to be prepared by the 
!                            user. For this, one needs to open the code or 
!                            better compile the code once in interactive mode 
!                            to know what set of inputs will be required and 
!                            prepare the input file. output.log will print 
!                            general information or error) 

!The code is written to analyse the dcd trajectories
!------------------------------------------------------------------------------
program cmeu
integer,parameter:: mxatom=60000,mxwater=17000,maxinp=3000,kmax=2500,mxset=10
real,dimension(mxatom):: x,y,z
real,dimension(mxwater,kmax):: xx,yy,zz
real(kind=8),dimension(0:kmax):: sd,n,av,avsq
integer,dimension(mxset):: site_res,num_res
integer,dimension(maxinp):: atp
integer,dimension(mxwater,kmax):: rm
integer:: funit,ounit,dummyi,res_type,total_atom
character(len=100),dimension(maxinp):: dcdfile
character(len=100):: outfile
character(len=100):: flgfile
character(len=4):: chr2,chr3
character(len=2):: chr4
logical:: there
real(kind=8):: dx,dy,dz,bx,by,bz
data funit,ounit,inunit /12,14,16/
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

if(num_res(3) > mxwater)then
write(*,*)'ERROR:NO OF TOTAL WATER EXCEEDS THE DECLARED DIMENSION'
write(*,*)'PLZ DO ADJUST THE SIZE OF MXWATER PARAMETER '
stop
endif
      
nsite=site_res(3)
ni=n_protein+1
nf=total_atom-n_ion-nsite+1
write(*,*)'first water oxygen atom no=======>',ni
write(*,*)'last water oxygen atom no=======>',nf
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the skip value'
read(*,*)nskip
write(*,*)'nskip=======>',nskip
write(*,'(1x,a)')'Type first & last resd no of the selected protein surface'
read(*,*) ifirst,ilast
write(*,*)'ifirst,ilast====>',ifirst,ilast
write(*,'(1x,a)')'Type two end residue no'
read(*,*)nfirst,nlast
write(*,*)'nfirst,nlast=====>',nfirst,nlast
write(*,'(1x,a)')'Enter the time gap betn two frames in ps'
read(*,*)dtstep
dtstep=dtstep*nskip
write(*,*)'dtstep=======>',dtstep
write(*,'(1x,a)')'Enter the inner and outer radius of selected region'
read(*,*)rin,rout
write(*,*)'rin,rout=========>',rin,rout
rin=rin*rin
rout=rout*rout
write(*,'(1x,a)')'Enter the no of block'
read(*,*)nset
write(*,*)'nset======>',nset
!-----------------------------------------------------------------------------
write(*,'(1x,a)')'Enter name of a reference file in PDB format'
read(*,'(a)')flgfile
write(*,'(a)')flgfile
open(funit,file=flgfile,status='old',form='formatted',iostat=ios)

npair=0
do i=1,n_protein
read(funit,'(12x,a4,1x,a4,1x,i4,50x,a2)')chr2,chr3,nres,chr4
if(nres == nfirst.or.nres == nlast)cycle
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
nused=0
nread=0
rm(:,:)=0
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
nread=nread+1
if(mod(nread,nskip) /= 0)cycle
nused=nused+1
!------------------------------------------------------------------------------
nw=0
do i=ni,nf,nsite
nw=nw+1
call MIN_DIST(i,npair,atp,rin,rout,x,y,z,dx,dy,dz,rmin)
if(rmin < rin.or.rmin > rout)cycle
p=(x(i+1)+x(i+2))/2.0
q=(y(i+1)+y(i+2))/2.0
r=(z(i+1)+z(i+2))/2.0

p=x(i)-p
q=y(i)-q
r=z(i)-r

v=p*p+q*q+r*r
v=sqrt(v)

xx(nw,nused)=p/v
yy(nw,nused)=q/v
zz(nw,nused)=r/v
rm(nw,nused)=1
enddo
enddo
enddo
!------------------------------------------------------------------------------
sd(:)=0.0
n(:)=0
cmeunorm=0.0
nc=0.0
do i=1,nw
do j=1,nused-1
if(rm(i,j) == 0)cycle
xref=xx(i,j)
yref=yy(i,j)
zref=zz(i,j)
cmeunorm=cmeunorm+xref*xref+yref*yref+zref*zref
nc=nc+1
do k=j+1,nused
if(rm(i,k) == 0)exit
m=k-j
rx=xx(i,k)*xref
ry=yy(i,k)*yref
rz=zz(i,k)*zref
sd(m)=sd(m)+rx+ry+rz
n(m)=n(m)+1
enddo
enddo
enddo

cmeunorm=cmeunorm/real(nc)
write(*,*)'cmeunorm========>',cmeunorm

ncorr=nused/3
do i=0,ncorr
if(n(i) == 0)n(i)=1
sd(i)=sd(i)/n(i)
if(i == 0)sd(i)=1.00
av(i)=av(i)+sd(i)
avsq(i)=avsq(i)+sd(i)*sd(i)
enddo
enddo

bar=0.0
do i=0,ncorr
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

bar=bar/real(ncorr+1)
write(ounit,*)'Average error bar in % for this calculation is = ',bar
write(ounit,*)'No of blocks used in this calculation = ',nset

end program cmeu
!------------------------------------------------------------------------------
!subroutine for calculating the minimum distance of the 
!selected water from the selected protein surface
subroutine MIN_DIST(j,np,atp,rin,rout,x,y,z,dx,dy,dz,rmin)
real,dimension(*),intent(in):: x,y,z
integer,intent(in):: j,np
integer,dimension(*),intent(in):: atp
real,intent(in):: rin,rout
real,intent(out):: rmin
real(kind=8),intent(in):: dx,dy,dz

rmin=1.0e6
do i=1,np
k=atp(i)
rx=x(k)-x(j)
ry=y(k)-y(j)
rz=z(k)-z(j)

rx=rx-dx*anint(rx/dx)
ry=ry-dy*anint(ry/dy)
rz=rz-dz*anint(rz/dz)

r=rx*rx+ry*ry+rz*rz
if(r < rmin)rmin=r
if(rmin >= rin.and.rmin <= rout)exit
enddo

end subroutine MIN_DIST
!--------------------------end of the program----------------------------------
