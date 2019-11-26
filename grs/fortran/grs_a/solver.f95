module solver

use common

!declaration of shared variables
integer :: lband
integer,parameter :: Nconstrain=1
real(rtype) :: alpha !streching coef. for retrieval of aerosol mixing ratio
real(rtype),dimension(:,:),allocatable :: s_af, s_ac
real(rtype),dimension(:),allocatable :: s_cextf, s_cextc ! extinction coefficients from lut (fine and coarse aerosol modes)
real(rtype),dimension(:),allocatable :: s_rot ! rayleigh optical thickness
real(rtype),dimension(:),allocatable :: s_m   ! air mass
real(rtype),dimension(:),allocatable :: s_rg_ratio ! sunglint spectral variation (referred to as epsilon in Harmel et al., RSE; 201?)
real(rtype),dimension(:),allocatable :: rsim, rmes ! reflectance from measurements
real(rtype),dimension(:),allocatable :: sigma ! uncertainty applied to each band
real(rtype) :: s_cextf550, s_cextc550  ! extinction coefficients from lut at 550 nm
real(rtype) :: aotguess


contains

subroutine aero_glint(aot550,beta,brdf)
!--------------------------------------------------------------------
!aot550: aerosol optical thickness at 550 nm
!beta: mixing ratio of the fine aerosol mode
!brdf: brdf in the glint direction for the SWIRest band
!--------------------------------------------------------------------

implicit none

real(rtype),intent(inout) :: aot550,beta,brdf

integer, parameter :: n=3 !number of unknows in optimization
integer :: info
real(rtype) :: norm
! epsfcn is an input variable used in determining a suitable
!       step length for the forward-difference approximation. this
!       approximation assumes that the relative errors in the
!       functions are of the order of epsfcn.
! factor is a positive input variable used in determining the
!       initial step bound. this bound is set to the product of
!       factor and the euclidean norm of diag*x if nonzero, or else
!       to factor itself. in most cases factor should lie in the
!       interval (.1,100.). 100. is a generally recommended value.
real(rtype) :: epsfcn,factor
integer :: mode
integer :: l,m,maxfev,nprint,nfev,ldfjac,ix,iy,iexist
integer,dimension(n):: ipvt
real(rtype) :: ftol,xtol,gtol,cri,eps,R2
real(rtype),dimension(n):: x,diag,qtf,wa1,wa2,wa3,fjnorm
real(rtype) :: enorm
real(rtype),dimension(:),allocatable :: fvec,wa4
real(rtype),dimension(:,:),allocatable :: fjac
integer :: i

m=lband+Nconstrain
!----------------------------------------------

!----------------------------------------------
!      LEVENBERG-MARQUARDT SETTINGS
!-- optimization parameters
factor=1d2
epsfcn=1d-6 !epsilon(epsfcn)
mode=1
ftol=sqrt(epsilon(1._rtype)) !1d-12
xtol=sqrt(epsilon(1._rtype))
gtol=sqrt(epsilon(1._rtype))
maxfev=100
nprint=0
ldfjac=m
alpha=1d-3
sigma=1d-3
sigma(1)=1d-2
aotguess=aot550
!----------------------------------------------

allocate(fvec(m),fjac(ldfjac,n),wa4(m))

!--- first guess

x(1)=dsqrt(aot550)
x(2)=-1d0/alpha * dlog(1d0/beta-1d0)
x(3)=dsqrt(brdf)

!--------Levenberg-Marquardt
call lmdif(cost_func,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,&
&                     diag,mode,factor,nprint,info,nfev,fjac,ldfjac,&
&                     ipvt,qtf,wa1,wa2,wa3,wa4)
!print*,'info ',info,'m ',m,' nfev',nfev ,'  fvec ',sqrt(sum(fvec**2d0))/m
!write(*,'(7(f8.4,x))')exp(x)

!print*,'fvec'
!do i=1,m
 !write(*,'(f15.4)')fvec(i)
   ! write(*,'(2(f12.5,x))')rmes(i),rsim(i)
!enddo
norm=dsqrt(sum(fvec**2d0))


!--- Compute upper bounds
eps=0.05
do i=1,n
 l=ipvt(i)
 fjnorm(l)=enorm(i,fjac(1,i))
enddo

!do i=1,n
!sigma(i)=sqrt(eps*(norm/fjnorm(i))**2d0)
!print*,x(i),' +-',sigma(i)
!enddo

aot550=(x(1))**2
beta=1._rtype / ( 1._rtype + dexp(-alpha*x(2)))
brdf=x(3)**2

print*,' aot, beta, brdf ',aot550,beta,brdf

deallocate(fvec,fjac,wa4)
return
end subroutine aero_glint


!-----------------------------------------------------------
! Cost function
!-----------------------------------------------------------

subroutine cost_func(m,n,x,fvec,iflag)
!Define cost function for Levenberg-Marquardt algo

implicit none

integer,intent(in) :: m,n
integer,intent(out) :: iflag
real(rtype),dimension(m),intent(out) :: fvec
real(rtype),dimension(n),intent(inout) :: x
integer :: i,ib

real(rtype) :: rsimf, rsimc, rg, aot, tud, aot550, beta, brdf

! x(1): aot550
! x(2): beta (mixing ratio of fine mode aerosols)
! x(3): brdf(2200nm)

aot550=(x(1))**2
beta=1._rtype / ( 1._rtype + dexp(-alpha*x(2)))
brdf=x(3)**2

do ib=1,lband
    rsimf = s_af(ib,0)
    rsimc = s_ac(ib,0)
    do i=1,norder
       rsimf = rsimf + s_af(ib,i)*aot550**i
       rsimc = rsimc + s_ac(ib,i)*aot550**i
    enddo

    aot = beta* s_cextf(ib)/s_cextf550 * aot550 + (1.-beta)* s_cextc(ib)/s_cextc550 * aot550

    tud = exp(- (s_rot(ib) + aot) / s_m(ib) )
    rg = tud * s_rg_ratio(ib) * brdf

    rsim(ib)=beta*rsimf + (1.-beta)*rsimc + rg
enddo


!--------------------------------------
! COST FUNCTION
!--------------------------------------
do i=1,m-Nconstrain
  if(i .ge. lband-2)then
      fvec(i)= ( rmes(i) - rsim(i) ) / sigma(i)
  else
      fvec(i)=0
      if((rmes(i) - rsim(i))<0.)fvec(i)= 1d-1 !( rmes(i) - rsim(i) ) / sigma(i)
  end if
enddo

!--------------------------------------
!  CONSTRAINTS
!--------------------------------------

fvec(m-Nconstrain+1) = 0. !( aot550 - aotguess ) !/ aotguess

return
end subroutine cost_func

end module