!code for van Hove auto-correlation function 
!date: 13th October, 2014
!modified on 4th May, 2015

!The code is written by Prabir Khatua

!Usage: f95 VAN_HOVE.F90 or gfortran VAN_HOVE.F90

!       ./a.out  (for interactive run and then enter the inputs interactivly 
!                 as the program will ask)
!       
!       ./a.out<vhove.inp >vhove.log &  
!                            (input.inp is the input file containing all the 
!                            required inputs that has to be prepared by the 
!                            user. For this, one needs to open the code or 
!                            better compile the code once in interactive mode 
!                            to know what set of inputs will be required and 
!                            prepare the input file. output.log will print 
!                            general information or error) 

!The code is written to analyse the dcd trajectories
!------------------------------------------------------------------------------
program msd
integer,parameter:: mxatom=60000,mxwater=17000,maxinp=2000
integer,parameter:: mxbin=7000,kmax=5000,mxset=10
real,dimension(mxatom):: x,y,z,xc,yc,zc
real,dimension(mxwater,kmax):: rm,xx,yy,zz
real(kind=8),dimension(mxbin):: h,av,avsq
integer,dimension(mxset):: site_res,num_res
integer,dimension(maxinp):: atp
integer:: funit,ounit,dummyi,res_type,total_atom,bin
character(len=100),dimension(maxinp):: dcdfile
character(len=100):: flgfile,outfile
character(len=4):: chr2,chr3
character(len=2):: chr4
logical:: there
real(kind=8):: dx,dy,dz,bx,by,bz,ntot
data funit,ounit,inunit /12,14,16/
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the name  of output file'
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
write(*,'(1x,a)')'Enter the timestep in ps'
read(*,*)dt
write(*,*)'dt====>',dt
write(*,'(1x,a)')'Type time interval(ps) at which you want to calculate van Hove function'
read(*,*)time
write(*,*)'time====>',time
t=time/dt
dif=t-int(t)
if(dif /= 0)then
write(*,'(1x,a)')'WARNING: REQUIRED TIME IS NOT INTEGRAL MULTIPLE OF TIMESTEP'
write(*,*)'PLZ ENTER THE TIME SUCH THAT IT IS AN INTEGRAL MULTIPLE OF TIMESTEP'
stop
endif
nskip=t
write(*,*)'nskip=====>',nskip
write(*,'(1x,a)')'Enter the first and last resd no that you want to analyse'
read(*,*)ifirst,ilast
write(*,*)'ifirst,ilast=======>',ifirst,ilast
write(*,'(1x,a)')'Type two end residue no'
read(*,*)nfirst,nlast
write(*,*)'nfirst,nlast=====>',nfirst,nlast
write(*,'(1x,a)')'Enter the inner and outer radius of selected region'
read(*,*)rin,rout
write(*,*)'rin,rout=========>',rin,rout
rin=rin*rin
rout=rout*rout
write(*,'(1x,a)')'Enter the maximum distance'
read(*,*)rmax
write(*,*)'rmax===>',rmax
write(*,'(1x,a)')'Enter the bin size'
read(*,*)dr
write(*,*)'dr====>',dr
write(*,'(1x,a)')'Enter the no of block'
read(*,*)nset
write(*,*)'nset======>',nset

maxbin=int(rmax/dr)+1
write(*,*)'maxbin====>',maxbin
if(maxbin > mxbin)then
write(*,*)'ERROR: THE NO OF BIN EXCEEDS THE DECLARED DIMENSION'
write(*,*)'PLZ DO RESET MXBIN PARAMETER'
stop
endif
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
ng=ninput/nset
av(:)=0.0
avsq(:)=0.0
do igroup=1,ninput,ng
nused=0
nread=0
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
!if(mod(nread,nskip) /= 0)cycle
nused=nused+1
!------------------------------------------------------------------------------
call CENTER(npair,atp,ni,nf,nsite,x,y,z,dx,dy,dz,xc,yc,zc)
nw=0
do i=ni,nf,nsite
nw=nw+1
call MIN_DIST(nw,npair,rin,rout,xc,yc,zc,rmin)
xx(nw,nused)=x(i)
yy(nw,nused)=y(i)
zz(nw,nused)=z(i)
if(rmin >= rin.and.rmin<= rout)then
rm(nw,nused)=1
else
rm(nw,nused)=0
endif
enddo
enddo
enddo
!------------------------------------------------------------------------------
h(:)=0.0
ntot=0.0
do j=1,nused-nskip
do i=1,nw
if(rm(i,j) == 0)cycle
xref=xx(i,j)
yref=yy(i,j)
zref=zz(i,j)
k=j+nskip
rx=xx(i,k)-xref
ry=yy(i,k)-yref
rz=zz(i,k)-zref
r=rx*rx+ry*ry+rz*rz
if(r <= rmax*rmax)then
r=sqrt(r)
bin=int(r/dr)+1
h(bin)=h(bin)+1
ntot=ntot+1
endif
enddo
enddo

do i=1,maxbin
h(i)=h(i)/ntot
av(i)=av(i)+h(i)
avsq(i)=avsq(i)+h(i)*h(i)
enddo
enddo

bar=0.0
do i=1,maxbin
av(i)=av(i)/real(nset)
avsq(i)=avsq(i)/real(nset)
r=(i-0.5)*dr
dif=avsq(i)-av(i)*av(i)
if(dif < 0.0)then
write(*,*)'WARNING: SOMETHING IS WRONG'
write(*,*)'STANDARD DEVIATION IS NEGATIVE'
dif=0.0
endif
xb=0.5*sqrt(dif)
write(ounit,'(3f12.6)')r,av(i),xb
xb=(xb*100.0)/av(i)
if(av(i) == 0.0)xb=0.0
bar=bar+xb
enddo

bar=bar/real(maxbin)
write(ounit,*)'Average error bar in % for this calculation is = ',bar


end program msd
!------------------------------------------------------------------------------
!subroutine for calculating the minimum distance of the 
!selected water from the selected protein surface
subroutine MIN_DIST(i,np,rin,rout,x,y,z,rmin)
real,dimension(*),intent(in):: x,y,z
integer,intent(in):: i,np
real,intent(in):: rin,rout
real,intent(out):: rmin

k=i+np
rmin=1.0e6
do j=1,np
rx=x(k)-x(j)
ry=y(k)-y(j)
rz=z(k)-z(j)

r=rx*rx+ry*ry+rz*rz
if(r < rmin)rmin=r
if(rmin >= rin.and.rmin <= rout)exit
enddo
end subroutine MIN_DIST
!------------------------------------------------------------------------------
subroutine CENTER(np,atp,ni,nf,ns,x,y,z,dx,dy,dz,xc,yc,zc)
real,dimension(*),intent(in):: x,y,z
real,dimension(*),intent(out):: xc,yc,zc
integer,dimension(*),intent(in):: atp
integer,intent(in):: ni,nf,ns
real(kind=8),intent(in):: dx,dy,dz

cx=0.0
cy=0.0
cz=0.0
do i=1,np
j=atp(i)
cx=cx+x(j)
cy=cy+y(j)
cz=cz+z(j)
enddo

cx=cx/real(np)
cy=cy/real(np)
cz=cz/real(np)

do i=1,np
j=atp(i)
xc(i)=x(j)-cx
yc(i)=y(j)-cy
zc(i)=z(j)-cz
enddo

k=np
do i=ni,nf,ns
k=k+1
rx=x(i)-cx
ry=y(i)-cy
rz=z(i)-cz

xc(k)=rx-dx*anint(rx/dx)
yc(k)=ry-dy*anint(ry/dy)
zc(k)=rz-dz*anint(rz/dz)
enddo

end subroutine CENTER
!--------------------------end of the program----------------------------------
