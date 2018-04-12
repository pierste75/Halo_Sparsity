PROGRAM HALO_SPARSITY

  USE COSMO
  USE POWER_SPECTRUM
  USE MF
  
  IMPLICIT NONE
  REAL(8) :: REDSHIFT, SPARSITY
  
  WRITE(*,*) 'Insert Omega_m,Omegabh2,h,sigma_8,ns'
  READ(*,*) OMEGA_M
  READ(*,*) OMEGABH2
  READ(*,*) h
  READ(*,*) SIGMA_8
  READ(*,*) ANS
  
  WRITE(*,*) 'Insert w_0,w_a'
  READ(*,*) w_0
  READ(*,*) w_a

  OMEGA_B=OMEGABH2/h**2
 
  
  OMEGA_R=2.48d0*1.d-5/h**2
  RHO_M=OMEGA_M*2.78*1.d+11   ! M_sun h^-1 / (h^-1 Mpc)^3
  
  R_to_M=3.d0/4.d0/PI/RHO_M

  ! Compute Growth Factor

  CALL GROWTH_FACTOR

  ! Compute Normalized Linear Matter Power Spectrum
  
  CALL MATTER_VARIANCE

  WRITE(*,*) 'Insert Redshift'
  READ(*,*) REDSHIFT

  write(*,'(a2,f5.3,a3,f9.5)') 's(',REDSHIFT,') =',SPARSITY(REDSHIFT)
        
END PROGRAM HALO_SPARSITY


FUNCTION SPARSITY(REDSHIFT)
  USE COSMO
  USE MF
    IMPLICIT NONE
    REAL(8) :: SPARSITY_TEMP,SPARSITY,REDSHIFT
    REAL(8) :: SD0_OLD,SD1_OLD,HS0_OLD,HS1_OLD
    REAL(8) :: CHANGE
    REAL(8), PARAMETER :: M1=1.d+13, M2=1.d+16
    INTEGER :: ITER
    INTEGER, PARAMETER :: NITER=300
    


    SD0_OLD=1.2d0
    SD1_OLD=1.3d0
    HS0_OLD=SD0_OLD*MF_M500c_INTEGRAL(SD0_OLD*M1,SD0_OLD*M2,REDSHIFT)- &
         & MF_M1000c_INTEGRAL(M1,M2,REDSHIFT)
    HS1_OLD=SD1_OLD*MF_M500c_INTEGRAL(SD1_OLD*M1,SD1_OLD*M2,REDSHIFT)- &
         & MF_M1000c_INTEGRAL(M1,M2,REDSHIFT)

    DO 100 ITER=1,NITER
       SPARSITY_TEMP=SD1_OLD-HS1_OLD*(SD1_OLD-SD0_OLD)/ &
            & (HS1_OLD-HS0_OLD)
       CHANGE=2.d0*ABS((SPARSITY_TEMP-SD1_OLD)/(SPARSITY_TEMP+SD1_OLD))
       IF(CHANGE.LE.1.d-5) GOTO 120
       SD0_OLD=SD1_OLD
       SD1_OLD=SPARSITY_TEMP
       HS0_OLD=SD0_OLD*MF_M500c_INTEGRAL(SD0_OLD*M1,SD0_OLD*M2,REDSHIFT)- &
            & MF_M1000c_INTEGRAL(M1,M2,REDSHIFT)
       HS1_OLD=SD1_OLD*MF_M500c_INTEGRAL(SD1_OLD*M1,SD1_OLD*M2,REDSHIFT)- &
            & MF_M1000c_INTEGRAL(M1,M2,REDSHIFT)
100    CONTINUE
       print*, 'Divergent Secant Method'
       print*,ITER,HS0_OLD,HS1_OLD
       GOTO 200
120    CONTINUE
200    CONTINUE

       SPARSITY=SPARSITY_TEMP
       
     END FUNCTION SPARSITY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                       NUMERICAL RECIPES SUBRTOUINES     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE NSPLINT(xa, ya, y2a, n, x, y)
  !   USE nrtype
  !
  ! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function
  ! (with the xa(i) in order), and given the array y2a(1:n), which is the output
  ! from the subroutine spline, and given a value of x, this routine returns a
  ! cubic spline interpolated value y.
  ! (adopted from Numerical Recipes in FORTRAN 77)
  !
  INTEGER:: n
  REAL(8):: x, y, xa(n), y2a(n), ya(n)
  INTEGER:: k, khi, klo
  REAL(8):: a, b, h
  klo=1
  khi=n
1 if (khi-klo.gt.1) then
     k=(khi+klo)/2
     if (xa(k).gt.x) then
        khi=k
     else
        klo=k
     endif
     goto 1
  endif
  h=xa(khi)-xa(klo)
  if (h.eq.0.) stop 'bad xa input in splint'
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
END SUBROUTINE NSPLINT

SUBROUTINE spline(x, y, n, yp1, ypn, y2)
  !   use nrtype
  !
  ! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e.
  ! y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1 and ypn for
  ! the first derivative of the interpolating function at points 1 and n,
  ! respectively, this routine returns an array y2(1:n) of length n which
  ! contains the second derivatives of the interpolating function at the
  ! tabulated points x(i).  If yp1 and/or ypn are equal to 1.e30 or larger,
  ! the routine is signaled to set the corresponding boundary condition for a
  ! natural spline with zero second derivative on that boundary.
  ! Parameter: nmax is the largest anticipiated value of n
  ! (adopted from Numerical Recipes in FORTRAN 77)
  !
  INTEGER:: n
  INTEGER, PARAMETER :: nmax=1000000
  REAL(8):: yp1, ypn, x(n), y(n), y2(n)
  INTEGER:: i, k
  REAL(8):: p, qn, sig, un, u(nmax)

  if (yp1.ge.1.d+30) then
     y2(1)=0.
     u(1)=0.
  else
     y2(1)=-0.5
     u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  endif
  do i=2, n-1
     sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
     p=sig*y2(i-1)+2.d0
     y2(i)=(sig-1.)/p
     u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/&
          & (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  enddo
  if (ypn.ge.1.d+30) then
     qn=0.
     un=0.
  else
     qn=0.5d0
     un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  endif

  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)

  do k=n-1, 1, -1
     y2(k)=y2(k)*y2(k+1)+u(k)
  enddo
END SUBROUTINE SPLINE
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine odeint(ystart,nvar,alfa1,alfa2,EPS,hh1,hmin,nok,nbad, &
     & derivs,rkqs)

  implicit double precision (a-h,o-z)

  integer nbad,nok,nvar,KMAXX,MAXSTP,NMAX,kmax,kount
  real*8 EPS,hh1,hmin,ystart(nvar),TINY,alfa1,alfa2
  external derivs,rkqs
  parameter (MAXSTP=1000000,NMAX=100,KMAXX=200,TINY=1.d-30)
  integer i,nstp
  real*8 dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX), &
       & yp(NMAX,KMAXX),yscal(NMAX)
  common /path/ kmax,kount,dxsav,xp,yp


  x=alfa1
  h=sign(hh1,alfa2-alfa1)
  nok=0
  nbad=0
  kount=0

  do i=1,nvar
     y(i)=ystart(i)
  enddo
  if (kmax.gt.0) xsav=x-2.*dxsav
  do nstp=1,MAXSTP
     call derivs(nvar,x,y,dydx)
     do i=1,nvar
        yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
     enddo
     if(kmax.gt.0)then
        if(abs(x-xsav).gt.abs(dxsav)) then
           if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do i=1,nvar
                 yp(i,kount)=y(i)
              enddo
              xsav=x
           endif
        endif
     endif
     if((x+h-alfa2)*(x+h-alfa1).gt.0.) h=alfa2-x
     call rkqs(y,dydx,nvar,x,h,EPS,yscal,hdid,hnext,derivs)
     if(hdid.eq.h)then
        nok=nok+1
     else
        nbad=nbad+1
     endif
     if((x-alfa2)*(alfa2-alfa1).ge.0.)then
        do i=1,nvar
           ystart(i)=y(i)
        enddo

        if(kmax.ne.0)then
           kount=kount+1
           xp(kount)=x
           do i=1,nvar
              yp(i,kount)=y(i)
           enddo
        endif
        return
     endif
     h=hnext
     if(abs(hnext).lt.hmin) then
        h=hmin
     end if
  enddo
  write(*,*) 'too many steps in odeint; x,hmin',x,alfa1,alfa2
  return
end subroutine odeint

subroutine rkqs(y,dydx,n,x,htry,EPS,yscal,hdid,hnext,derivs)

  implicit double precision (a-h,o-z)

  integer n,NMAX
  real*8 EPS,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
  external derivs
  parameter (NMAX=100)
  integer i
  real*8 errmax,h,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW
  real*8   ERRCON,PSHRNK
  parameter (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89d-4)     
  h=htry
1 call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
  errmax=0.
  do i=1,n
     errmax=max(errmax,abs(yerr(i)/yscal(i)))
  enddo
  errmax=errmax/EPS
  if(errmax.gt.1.)then
     h=SAFETY*h*(errmax**PSHRNK)
     if(h.lt.0.1*h)then
        h=.1*h
     endif
     xnew=x+h
     goto 1
  else
     if(errmax.gt.ERRCON)then
        hnext=SAFETY*h*(errmax**PGROW)
     else
        hnext=5.*h
     endif
     hdid=h
     x=x+h

     do i=1,n
        y(i)=ytemp(i)
     enddo
     return
  endif
END subroutine rkqs

subroutine rkck(y,dydx,n,x,h,yout,yerr,derivs)

  implicit double precision (a-h,o-z)

  integer n,NMAX
  real*8 h,x,dydx(n),y(n),yerr(n),yout(n)
  external derivs
  parameter (NMAX=100)
  integer i
  real*8 ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX), &
       & ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53, &
       & B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
  parameter (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40., &
       & B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5, &
       & B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512., &
       & B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378., &
       & C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648., &
       & DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336., &
       & DC6=C6-.25)

  do i=1,n
     ytemp(i)=y(i)+B21*h*dydx(i)
  enddo
  call derivs(n,x+A2*h,ytemp,ak2)
  do i=1,n
     ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
  enddo
  call derivs(n,x+A3*h,ytemp,ak3)
  do i=1,n
     ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
  enddo
  call derivs(n,x+A4*h,ytemp,ak4)
  do i=1,n
     ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53* &
          & ak3(i)+B54*ak4(i))
  enddo
  call derivs(n,x+A5*h,ytemp,ak5)
  do i=1,n
     ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+ &
          & B64*ak4(i)+B65*ak5(i))
  enddo
  call derivs(n,x+A6*h,ytemp,ak6)
  do i=1,n
     yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
  enddo
  do i=1,n
     yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6* &
          & ak6(i))
  enddo
  return
END subroutine rkck
