!------------------------------------------------------------------------------
!program for water entropy around protein using 2PT model. First one should
!need to check first the VACF. Then decide the no of data points to be
!considered for fourier transform depending upon convergence of VACF to zero.
!Here decision of no of data points is very important since we will be using
!Filon's method for fourier transform which is a discrete fourier transform of
!finite no of data points and highly depends upon no of data points.
!______________________________________________________________________________
!References: I.  Lin et al. JCP 119,11792, 2003
!            II. Lin et al. JPCB 114, 24, 2010
!______________________________________________________________________________
!The program is written by Prabir Khatua
!date: 27th October, 2014
!modified on 14th April, 2015

!Usage: f95 HYD_2PT_ENTROPY_VER3.F90 or gfortran HYD_2PT_ENTROPY_VER3.F90

!       ./a.out  (for interactive run and then enter the inputs interactivly 
!                 as the program will ask)
!       
!       ./a.out<2pt.inp >2pt.log &  
!                            (input.inp is the input file containing all the 
!                            required inputs that has to be prepared by the 
!                            user. For this, one needs to open the code or 
!                            better compile the code once in interactive mode 
!                            to know what set of inputs will be required and 
!                            prepare the input file. output.log will print 
!                            general information or error) 

!The code is written to analyse the dcd trajectories
!------------------------------------------------------------------------------
program ENTROPY
integer,parameter:: mxatom=50000,maxinp=5000
integer,parameter:: mxw=18033,kmax=4900,mxtype=10
real,dimension(mxatom):: x,y,z,vx,vy,vz,xc,yc,zc
real(kind=8),dimension(0:kmax):: vt,vr,n,vrn,vtn
real,dimension(0:kmax):: s,sg,ss,sws,ws
real,dimension(3):: dd
real,dimension(3,3):: xi,v
real,dimension(mxw,kmax):: wx,wy,wz,vvx,vvy,vvz
integer,dimension(mxw,kmax):: h
integer,dimension(mxtype):: site_res,num_res
integer,dimension(maxinp):: atp
integer:: ounit,funit,vunit,dummyi,res_type
integer:: punit,tunit,cunit,total_atom
character(len=100),dimension(maxinp):: dcdfile,vdcdfile
character(len=100):: vfile1,vfile2,pdbfile
character(len=100):: efile1,efile2,outfile
character(len=4):: chr2,chr3
character(len=2):: chr4
logical:: there
real(kind=8):: dx,dy,dz,bx,by,bz,rho
data inunit,vunit,ounit,tunit /10,11,12,13/
data funit,cunit,punit,kunit /14,15,18,19/
!------------------------------------------------------------------------------
xmas=18.0154     !xmax = mass of water in g/mol unit
cf=20.45482706   !cf = conversion factor to convert namd velocity in A/ps unit
sigma=2.0        !Rotational symmetry for water
!------------------------------------------------------------------------------
write(*,'(1x,a)') 'Enter two output file for trans and rot VACF'
read(*,'(a)')vfile1
read(*,'(a)')vfile2
write(*,'(1x,a)')'Enter two output file for trans and rot spectra'
read(*,'(a)')efile1
read(*,'(a)')efile2
write(*,'(1x,a)')'Enter output file for entropy'
read(*,'(a)')outfile

open(ounit,file=vfile1,status='unknown',form='formatted',iostat=ios)
open(tunit,file=vfile2,status='unknown',form='formatted',iostat=ios)
open(funit,file=efile1,status='unknown',form='formatted',iostat=ios)
open(punit,file=efile2,status='unknown',form='formatted',iostat=ios)
open(cunit,file=outfile,status='unknown',form='formatted',iostat=ios)
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
write(*,*)'no of atoms per water molecule===>',nsit
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the skip value'
read(*,*)nskip
write(*,*)'nskip=======>',nskip
write(*,'(1x,a)')'Type first & last resd no that to be analysed'
read(*,*) ifirst,ilast
write(*,*)'ifirst,ilast========>',ifirst,ilast
write(*,'(1x,a)')'Type two end residue no'
read(*,*)nfirst,nlast
write(*,*)'nfirst,nlast=====>',nfirst,nlast
write(*,'(1x,a)')'Enter the inner and outer radius of selected region'
read(*,*)rin,rout
write(*,*)'rin,rout=========>',rin,rout
rin=rin*rin
rout=rout*rout
write(*,'(1x,a)')'Enter the time gap betn two frames in ps'
read(*,*)dt
dt=dt*nskip
write(*,*)'dtstep=======>',dt
write(*,'(1x,a)')'Enter the temperature in Kelvin'
read(*,*)temp
write(*,*)'temp===>',temp
write(*,'(1x,a)')'Type # of data reqd for fourier transfrm of trans, rot motion'
read(*,*)ndt,ndr
write(*,*)'ndt,ndr====>',ndt,ndr
write(*,'(1x,a)')'Enter the number density in angs^3'
read(*,*)rho
write(*,*)'rho=====>',rho
write(*,'(1x,a)')'Enter the no of block'
read(*,*)nset
write(*,*)'nset======>',nset

if(mod(ndt,2).ne.0)ndt=ndt-1
if(mod(ndr,2).ne.0)ndr=ndr-1


vt(:)=0.0
vr(:)=0.0
vrn(:)=0.0
vtn(:)=0.0
n(:)=0.0
!------------------------------------------------------------------------------
write(*,'(1x,a)')'Enter the name of pdbfile'
read(*,'(a)')pdbfile
write(*,*)pdbfile
open(kunit,file=pdbfile,status='old',form='formatted',iostat=ios)

np=0
do i=1,n_protein
read(kunit,'(12x,a4,1x,a4,1x,i4,50x,a2)')chr2,chr3,nres,chr4
if(chr4 == ' H')cycle
if(nres == nfirst.or.nres == nlast)cycle
if(nres >= ifirst.and.nres <= ilast)then
np=np+1
if(np > maxinp)then
write(*,*)'ERROR: THE NO OF SELECTED ATOM EXCEEDS THE DIMENSION'
write(*,*)'PLZ DO RESET MAXINP PARAMETER'
stop
endif

atp(np)=i
endif
enddo

write(*,*)'npair======>',np
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
write(*,'(a)')'ERROR: DCDFILE NOT FOUND'
stop
endif

write(*,'(1x,a)')'Enter the name of veldcd file',inp
read(*,'(a)')vdcdfile(inp)
write(*,*)vdcdfile(inp)
inquire(file=vdcdfile(inp),exist=there)
if(.not.there)then
write(*,'(a)')'ERROR: VELDCD FILE NOT FOUND'
stop
endif
enddo
!------------------------------------------------------------------------------
ng=ninput/nset
do igroup=1,ninput,ng
nu=0
nr=0
h(:,:)=0
inploop: do inp=igroup,igroup+ng-1
write(*,'(1x,a,i5,a)')'Enter name of trajectory file',inp
write(*,*)dcdfile(inp)
open(inunit,file=dcdfile(inp),status='old',form='unformatted')
read(inunit)dummyc,nframes,(dummyi,i=1,8),dummyr,(dummyi,i=1,9)
read(inunit)dummyi,dummyr
write(*,*)"dummyi,dummyr====>",dummyi,dummyr
read(inunit)natom
write(*,*)"natom====>",natom
write(*,*)"nframes===>",nframes

write(*,'(1x,a,i5,a)')'Enter name of veldcd file',inp
write(*,*)vdcdfile(inp)
open(vunit,file=vdcdfile(inp),status='old',form='unformatted')
read(vunit)dummyc,nframes,(dummyi,i=1,8),dummyr,(dummyi,i=1,9)
read(vunit)dummyi,dummyr
write(*,*)"dummyi,dummyr====>",dummyi,dummyr
read(vunit)natom
write(*,*)"natom====>",natom
write(*,*)"nframes===>",nframes
      
do ii=1,nframes
read(inunit)dx,bx,dy,by,bz,dz
read(inunit)(x(j),j=1,natom)
read(inunit)(y(j),j=1,natom)
read(inunit)(z(j),j=1,natom)

read(vunit)(vx(j),j=1,natom)
read(vunit)(vy(j),j=1,natom)
read(vunit)(vz(j),j=1,natom)
nr=nr+1
if(mod(nr,nskip) /= 0)cycle
nu=nu+1
if(nu > kmax)exit inploop
!------------------------------------------------------------------------------
call CENTER(np,atp,ni,nf,nsit,x,y,z,dx,dy,dz,xc,yc,zc)
nw=0
nc=0
do  i=ni,nf,nsit
nw=nw+1

!calculation of center of mass and center of mass weigheted velocity
!of water molecule for the calculation of translational VAC

call MIN_DIST(nw,np,xc,yc,zc,rin,rout,rm)
if(rm < rin.or.rm > rout)cycle
nc=nc+1
h(nw,nu)=1

cx=0.0
cy=0.0
cz=0.0
cvx=0.0
cvy=0.0
cvz=0.0
tm=0.0
k=0
do j=i,i+nsit-1
k=k+1
if(k == 1)then
xm=15.99940
else
if(k == 2.or.k == 3)then
xm=1.0080
else
xm=0.0
endif
endif

cx=cx+xm*x(j)
cy=cy+xm*y(j)
cz=cz+xm*z(j)

cvx=cvx+xm*vx(j)
cvy=cvy+xm*vy(j)
cvz=cvz+xm*vz(j)
tm=tm+xm
enddo

vvx(nw,nu)=cf*cvx/tm
vvy(nw,nu)=cf*cvy/tm
vvz(nw,nu)=cf*cvz/tm

cx=cx/tm
cy=cy/tm
cz=cz/tm

!calculation of angular velocity of water molecule

xi(:,:)=0.0 !xi(i,j)=ith row and jth column element of moment of inertia tensor
xl=0.0      !xl= X component of total angular momentum of water molecule
yl=0.0      !yl= Y component of total angular momentum of water molecule
zl=0.0      !zl= Z component of total angular momentum of water molecule
k=0
do j=i,i+nsit-1
k=k+1
if(k ==1)then
xm=15.99940
else
if(k == 2.or.k == 3)then
xm=1.0080
else
xm=0.0
endif
endif

vx(j)=cf*vx(j)
vy(j)=cf*vy(j)
vz(j)=cf*vz(j)

rx=x(j)-cx
ry=y(j)-cy
rz=z(j)-cz

xl=xl+xm*(ry*vz(j)-rz*vy(j))
yl=yl+xm*(rz*vx(j)-rx*vz(j))
zl=zl+xm*(rx*vy(j)-ry*vx(j))

xi(1,1)=xi(1,1)+xm*(ry*ry+rz*rz)
xi(2,2)=xi(2,2)+xm*(rx*rx+rz*rz)
xi(3,3)=xi(3,3)+xm*(ry*ry+rx*rx)
xi(1,2)=xi(1,2)-xm*rx*ry
xi(2,3)=xi(2,3)-xm*ry*rz
xi(1,3)=xi(1,3)-xm*rx*rz
enddo

xi(2,1)=xi(1,2)
xi(3,2)=xi(2,3)
xi(3,1)=xi(1,3)

if(nu == 1.and.nc == 1)then
!diagonalization of moment of inertia tensor xi will give the principal
!component of moment of inertia 
call jacobi(xi,3,3,dd,v,nrot)
write(*,*)(dd(l),l=1,3)
endif

call invt33(xi)

!inverse matrix of moment of inertia tensor multiplied with total angular 
!momentum vector will give the angular velocity

wx(nw,nu)=xi(1,1)*xl+xi(1,2)*yl+xi(1,3)*zl
wy(nw,nu)=xi(2,1)*xl+xi(2,2)*yl+xi(2,3)*zl
wz(nw,nu)=xi(3,1)*xl+xi(3,2)*yl+xi(3,3)*zl
enddo
enddo
enddo inploop
!------------------------------------------------------------------------------
if(nu > kmax)nu=nu-1
do i=1,nw
do j=1,nu
if(h(i,j) == 0)cycle
vxr=vvx(i,j)
vyr=vvy(i,j)
vzr=vvz(i,j)
wxr=wx(i,j)
wyr=wy(i,j)
wzr=wz(i,j)
do k=j,nu
if(h(i,k) == 0)exit
m=k-j
vxt=vvx(i,k)
vyt=vvy(i,k)
vzt=vvz(i,k)
wxt=wx(i,k)
wyt=wy(i,k)
wzt=wz(i,k)

vt(m)=vt(m)+vxr*vxt+vyr*vyt+vzr*vzt
vr(m)=vr(m)+dd(1)*wxr*wxt+dd(2)*wyr*wyt+dd(3)*wzr*wzt
vtn(m)=vtn(m)+vxr*vxr+vyr*vyr+vzr*vzr
vrn(m)=vrn(m)+dd(1)*wxr*wxr+dd(2)*wyr*wyr+dd(3)*wzr*wzr
n(m)=n(m)+1.0
enddo
enddo
enddo
enddo
!------------------------------------------------------------------------------
nu=nu/2
do i=0,nu
tv=vt(i)/vtn(i)
rv=vr(i)/vrn(i)
vt(i)=xmas*vt(i)/n(i)
vr(i)=vr(i)/n(i)
t=i*dt
write(ounit,'(3f12.6)')t,vt(i),tv
write(tunit,'(3f12.6)')t,vr(i),rv
enddo
!---------------calculation of translational entropy---------------------------
dw=1.0/(dt*real(ndt))
call FTRANS(ndt,kmax,dt,dw,vt,s,temp)

s0=s(0)
call DELTA(s0,temp,xmas,rho,d)
call NEW_RAP(d,f)
call DOS(ndt,kmax,dw,f,s,s0,sg,ss)
call WFRAC(ndt,kmax,dw,temp,ws)

write(funit,*)'# frequency (cm^-1), total DOS (cm), solid(cm), gas(cm)'
do i=0,ndt
sws(i)=ss(i)*ws(i)
w=dw*i*5.3078 
     
!5.3078 converts angular frequency w (in ps^-1) to to wave number;
!neu_bar in (cm^-1). neu_bar = w/(2xpixc); c=velocity of light
!Here DOS is in ps unit and calculated from fourier transform of VACF by
!integration over c(t)cos(neuxt). 0.03 is multiplied with the DOS to change
!its unit to cm; neu_bar=neu/c

write(funit,'(4f10.5)')w,s(i)*0.03,ss(i)*0.03,sg(i)*0.03
enddo

call SIMPSON(sws,dw/6.283185,kmax,ndt,es)
call SIMPSON(sg,dw/6.283185,kmax,ndt,eg)
call HST(f,d,xmas,temp,rho,yy,zy,shs,wg)

es=es*8.314
eg=wg*eg*8.314
ts1=es+eg
difc=(s0*0.8314*temp)/(12.0*xmas)
write(cunit,*)'TRANSLATIONAL MOTION:'
write(cunit,*)'Diffusion constant =',difc
write(cunit,*)'delta=',d
write(cunit,*)'fluidicity=',f
write(cunit,*)'y=',yy
write(cunit,*)'compressibility factor z(y)=',zy
write(cunit,*)'Hard-sphere entropy SHS=',shs*8.314
write(cunit,*)'entropy_solid=',es
write(cunit,*)'entropy_gas=',eg
write(cunit,*)'entropy_solid+entropy_gas=',ts1
!---------------calculation of rotational entropy---------------------------
dw=1.0/(dt*real(ndr))
call FTRANS(ndr,kmax,dt,dw,vr,s,temp)
s0=s(0)
call DELTA(s0,temp,xmas,rho,d)
call NEW_RAP(d,f)
call DOS(ndr,kmax,dw,f,s,s0,sg,ss)
call WFRAC(ndr,kmax,dw,temp,ws)

write(punit,*)'# frequency (cm^-1), total DOS (cm), solid(cm), gas(cm)'
do i=0,ndr
sws(i)=ss(i)*ws(i)

!5.3078 converts angular frequency w (in ps^-1) to to wave number;
!neu_bar in (cm^-1). neu_bar = w/(2xpixc); c=velocity of light
!Here DOS is in ps unit and calculated from fourier transform of VACF by
!integration over c(t)cos(neu x t). 0.03 is multiplied with the DOS to change
!its unit to cm; neu_bar=neu/c

w=dw*i*5.3078
write(punit,'(4f10.5)')w,s(i)*0.03,ss(i)*0.03,sg(i)*0.03
enddo

call SIMPSON(sws,dw/6.283185,kmax,ndr,es)
call SIMPSON(sg,dw/6.283185,kmax,ndr,eg)
call HSR(temp,sigma,dd(1),dd(2),dd(3),shs,wg)

es=es*8.314
eg=wg*eg*8.314
ts2=es+eg

difc=(s0*0.8314*temp)/(12.0*xmas)
write(cunit,*)
write(cunit,*)'ROTATIONAL MOTION:'
write(cunit,*)'Diffusion constant =',difc
write(cunit,*)'delta=',d
write(cunit,*)'fluidicity=',f
write(cunit,*)'Hard-sphere entropy SHS=',shs*8.314
write(cunit,*)'entropy_solid=',es
write(cunit,*)'entropy_gas=',eg
write(cunit,*)'entropy_solid+entropy_gas=',ts2

write(cunit,*)
write(cunit,*)'Total entropy =',ts1+ts2
write(cunit,*)'Entropy is given in J mol^-1 K^-1 unit'
write(cunit,*)'Diffusion coefficient in Angs^2 ps^-1 unit'


end program ENTROPY
!------------------------------------------------------------------------------
!This routine calculates different properties of hard sphere gas
!for the rotational motion
subroutine HSR(t,s,ai,bi,ci,shs,wg)
real,intent(in):: t,ai,bi,ci,s
real,intent(out):: shs,wg

pi=3.1415925
th=(6.625**2*100)/(8.0*pi**2*1.38)
tha=th/ai
thb=th/bi
thc=th/ci

xx=sqrt(t**3/(tha*thb*thc))
shs=log(sqrt(pi)/s)+log(xx)+(3.0/2.0)

wg=shs/3.0

end subroutine HSR
!------------------------------------------------------------------------------
!This routine calculates different properties of hard sphere gas
!for the translational motion
!------------------------------------------------------------------------------
subroutine HST(f,d,xm,t,rho,y,zy,shs,wg)
real,intent(in):: f,d,xm,t
real(kind=8),intent(in):: rho
real,intent(out):: y,zy,shs,wg

pi=3.1415925
y=f**(5.0/2.0)/d**(3.0/2.0)

zy=(1.0+y+y*y-y**3)/(1.0-y)**3

f1=(3.0/2.0)*log((2.0*pi*1.38*1.66*xm*t)/(6.625**2*100))
f2=log(zy/(rho*f))
f3n=y*(3.0*y-4.0)
f3d=(1.0-y)**2
f3=f3n/f3d

shs=(5.0/2.0)+f1+f2+f3
wg=shs/3.0

end subroutine HST
!------------------------------------------------------------------------------
!This routine does numerical integration by simpson's 1/3 rule
subroutine SIMPSON(f,h,nmax,n,xs)
integer,intent(in):: n,nmax
real,dimension(0:nmax),intent(in):: f(0:nmax)
real,intent(in):: h
real,intent(out):: xs

xs=0.0
do i=1,n-1
if(mod(i,2).eq.0)then
a=2.0
else
a=4.0
endif

xs=xs+a*f(i)
enddo

xs=xs+f(0)+f(n)
xs=xs*h/3.0

end subroutine SIMPSON
!------------------------------------------------------------------------------
!This routine calculats weighting fraction for solid component
!------------------------------------------------------------------------------
subroutine WFRAC(ndata,nmax,dw,t,ws)
integer,intent(in):: ndata,nmax
real,intent(in):: dw,t
real,dimension(0:nmax),intent(out):: ws

xk=66.25/(1.38*t)         !xk=h/kt    h = Planck's constant
                          !k = Boltzman constant  t = temperature
                          !unit of xk is in ps
do i=0,ndata
if(i == 0)then
ws(i)=0.0
else
w=real(i)*dw/6.283185
xp=xk*w
f1=xp/(exp(xp)-1.0)
f2=log(1.0-exp(-xp))
ws(i)=f1-f2
endif
enddo

end subroutine WFRAC
!------------------------------------------------------------------------------
!This routine calculates of density of state (DOS) for the gas and solid
!component separately. Gas component is calculated considering hard sphere.
!Solid component is calculated subtracing gas component from total DOS which 
!has been calculated from fourier tranform of VACF
!------------------------------------------------------------------------------
subroutine DOS(ndata,nmax,dw,f,s,s0,sg,ss)
real,dimension(0:nmax),intent(in):: s
real,dimension(0:nmax):: sg,ss

do i=0,ndata
w=real(i)*dw/6.283185
xd=(3.1415925*s0*w)/(6.0*f)
xd=1.0+xd*xd
sg(i)=s0/xd

ss(i)=s(i)-sg(i)
enddo

end subroutine DOS
!------------------------------------------------------------------------------
!This routine calculates fluidicity parameter (f) from the dimensionless
!diffusivity constant value by Newton-Rhapson's method
!------------------------------------------------------------------------------
subroutine NEW_RAP(d,f)
real,intent(in):: d
real,intent(out):: f
real,dimension(1000):: x

xl=0.00005
x(1)=0.0
l=0
45    k=1
46    f=x(k)
y1=2.0*d**(-9.0/2.0)*f**(15.0/2.0)-6.0*d**(-3.0)*f**5.0
y2=-d**(-3.0/2.0)*f**(7.0/2.0)+6.0*d**(-3.0/2.0)*f**(5.0/2.0)
y3=2.0*f-2.0

y=y1+y2+y3

dy1=15.0*d**(-9.0/2.0)*f**(13.0/2.0)-30.0*d**(-3.0)*f**4.0
dy2=-(7.0/2.0)*d**(-3.0/2.0)*f**(5.0/2.0)
dy3=15.0*d**(-3.0/2.0)*f**(3.0/2.0)+2.0

dy=dy1+dy2+dy3

x(k+1)=x(k)-(y/dy)
err=x(k+1)-x(k)

if(abs(err).le.xl)goto 47
k=k+1
if(k.le.500)then
goto 46
else
l=l+1
x(1)=x(1)+0.02
goto 45
endif

47    f=x(k+1)

end subroutine NEW_RAP      
!------------------------------------------------------------------------------
!This routine calculates dimensionless diffusivity (delta) from  
!zero frequency DOS value corresponding to one atom/molecule
!------------------------------------------------------------------------------
subroutine DELTA(s0,t,xm,rho,d)
real,intent(in):: s0,t,xm
real(kind=8),intent(in):: rho
real,intent(out):: d

pi=3.1415925

f1=(2.0*s0)/9.0
f2=sqrt((pi*0.8314*t)/xm)        !unit of R= gA^2ps^-2mol^-1K^-1
f3=rho**(1.0/3.0)
f4=(6.0/pi)**(2.0/3.0)

d=f1*f2*f3*f4

end subroutine DELTA
!------------------------------------------------------------------------------
!This routine calculates fourier transform of VACF by Filon's method to
!get density of state
!------------------------------------------------------------------------------
subroutine FTRANS(ndata,nmax,dt,dw,vac,sp,temp)
real(kind=8),dimension(0:nmax),intent(in):: vac
real,dimension(0:nmax),intent(out):: sp

sps=0.0
tmax=dt*real(ndata)
do i=0,ndata
w=real(i)*dw
th=w*dt
!----------------------calculate the Filon parameters)-------------------------
sint=sin(th)
cost=cos(th)
sinsq=sint*sint
cossq=cost*cost
thsq=th*th
thcub=thsq*th

if(th.eq.0)then
a=0.0
b=2.0/3.0
g=4.0/3.0
else
a=(1.0/thcub)*(thsq+th*sint*cost-2.0*sinsq)
b=(2.0/thcub)*(th*(1.0+cossq)-2.0*sint*cost)
g=(4.0/thcub)*(sint-th*cost)
endif
!----------------------(do sum over even ordinates)----------------------------
s2=0.0
do j=0,ndata,2
s2=s2+vac(j)*cos(th*real(j))
enddo

s2=s2-0.5*(vac(0)+vac(ndata)*cos(w*tmax))
!----------------------(do sum over odd ordinates)-----------------------------
s1=0.0
do j=1,ndata-1,2
s1=s1+vac(j)*cos(th*real(j))
enddo

sp(i)=4.0*(a*vac(ndata)*sin(w*tmax)+b*s2+g*s1)*dt/(0.8314*temp)
sps=sps+dw*sp(i)/6.283185              !6.283185 convert omega to neu
enddo

cf=3.0/sps             !this factor is multiplied to make sum of 
                             !power spectrum 3N (in this case N=1)
do i=0,ndata
sp(i)=sp(i)*cf
enddo
write(*,*)'sps====>',sps

end subroutine FTRANS
!------------------------------------------------------------------------------
!this subroutine calculates the inverse matrix of a 3x3 matrix
subroutine invt33(a)
real,dimension(3,3),intent(inout):: a
real,dimension(3,3):: aa

aa(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
aa(1,2)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
aa(1,3)=a(2,1)*a(3,2)-a(2,2)*a(3,1)
aa(2,1)=a(1,3)*a(3,2)-a(1,2)*a(3,3)
aa(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1)
aa(2,3)=a(1,2)*a(3,1)-a(1,1)*a(3,2)
aa(3,1)=a(1,2)*a(2,3)-a(2,2)*a(1,3)
aa(3,2)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
aa(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      
det=a(1,1)*aa(1,1)+a(1,2)*aa(1,2)+a(1,3)*aa(1,3)
      
do i=1,3
do j=1,3
a(j,i)=aa(i,j)/det
enddo
enddo

end subroutine invt33                   
!------------------------------------------------------------------------------
!this routine diagonalizes a symmetric matrix to get eigen vector and eigen
!value by jacobi rotation method
!------------------------------------------------------------------------------
subroutine jacobi(a,n,np,d,v,nrot)
integer,parameter:: NMAX=500
integer,intent(in):: n,np
real,dimension(np,np),intent(inout):: a
real,dimension(np,np),intent(out):: v
real,dimension(np),intent(out):: d
real,dimension(NMAX):: b,z

do ip=1,n
do iq=1,n
v(ip,iq)=0
enddo
v(ip,ip)=1
enddo

do ip=1,n
b(ip)=a(ip,ip)
d(ip)=b(ip)
z(ip)=0
enddo

nrot=0
do i=1,50
sm=0
do ip=1,n-1
do iq=ip+1,n
sm=sm+abs(a(ip,iq))
enddo
enddo

if(sm == 0)return
if(i < 4)then
tresh=0.2*sm/n**2
else
tresh=0
endif

do ip=1,n-1
do iq=ip+1,n
g=100.*abs(a(ip,iq))
if((i > 4).and.(abs(d(ip))+g == abs(d(ip))).and.(abs(d(iq))+g == abs(d(iq))))then
a(ip,iq)=0.
else if(abs(a(ip,iq)) > tresh)then
h=d(iq)-d(ip)
if(abs(h)+g == abs(h))then
t=a(ip,iq)/h
else
theta=0.5*h/a(ip,iq)
t=1./(abs(theta)+sqrt(1.+theta**2))
if(theta < 0.)t=-t
endif

c=1./sqrt(1+t**2)
s=t*c
tau=s/(1.+c)
h=t*a(ip,iq)

z(ip)=z(ip)-h
z(iq)=z(iq)+h
d(ip)=d(ip)-h
d(iq)=d(iq)+h
a(ip,iq)=0.
do j=1,ip-1 !Case of rotations 1  j < p.
g=a(j,ip)
h=a(j,iq)
a(j,ip)=g-s*(h+g*tau)
a(j,iq)=h+s*(g-h*tau)
enddo

do j=ip+1,iq-1 !Case of rotations p < j < q.
g=a(ip,j)
h=a(j,iq)
a(ip,j)=g-s*(h+g*tau)
a(j,iq)=h+s*(g-h*tau)
enddo

do j=iq+1,n !Case of rotations q < j  n.
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
z(ip)=0.    !and reinitialize z.
enddo
enddo
print*, 'too many iterations in jacobi'

END subroutine jacobi
!------------------------------------------------------------------------------
subroutine MIN_DIST(i,np,x,y,z,rin,rout,rmin)
real,dimension(*),intent(in):: x,y,z
real,intent(in):: rin,rout
real,intent(out):: rmin
integer,intent(in):: i,np

k=i+np
rmin=1.0e6
do  j=1,np
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
integer,intent(in):: ni,nf,ns,np
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

nh=np
do i=ni,nf,ns
nh=nh+1
rx=x(i)-cx
ry=y(i)-cy
rz=z(i)-cz

xc(nh)=rx-dx*anint(rx/dx)
yc(nh)=ry-dy*anint(ry/dy)
zc(nh)=rz-dz*anint(rz/dz)
enddo

end subroutine CENTER
!------------------------------------------------------------------------------
