MODULE COSMO

  REAL(8), PARAMETER :: PI=3.141592653589793
  REAL(8), PARAMETER :: c=299790.d0

  REAL(8) :: OMEGA_M,h,SIGMA_8,w_0,w_a
  REAL(8) :: OMEGA_R,OMEGA_B,OMEGABH2,ANS,RHO_M
  REAL(8) :: R_to_M
  
  INTEGER, PARAMETER :: NZ=400
  REAL(8), DIMENSION(NZ) :: AFIN,DPLUS,D2DPLUS


CONTAINS

  !--- COMPUTE H(z) for given cosmology ------------------------------------

  FUNCTION f_DE(z)
    IMPLICIT NONE
    REAL(8) :: f_DE,z
    REAL(8) :: a
    a=1.d0/(1.d0+z)
    f_DE=exp(3.d0*w_a*(a-1.d0))*(1.d0+z)**(3.d0*(1.d0+w_0+w_a))
  END FUNCTION f_DE

  FUNCTION inv_hubblez(z)
    IMPLICIT NONE
    REAL(8) :: inv_hubblez,z,hz
    hz=OMEGA_M*(1.d0+z)**3+OMEGA_R*(1.d0+z)**4+(1.d0-OMEGA_M-OMEGA_R)*f_DE(z)
    inv_hubblez=1.d0/sqrt(hz)
  END FUNCTION inv_hubblez


  FUNCTION D_PLUS(z)
    IMPLICIT NONE
    REAL(8) :: D_PLUS,z
    REAL(8) :: a
    a=1.d0/(1.d0+z)
    CALL NSPLINT(AFIN,DPLUS,D2DPLUS,NZ,a,D_PLUS)
  END FUNCTION D_PLUS

  !--- COMPUTE D+(z) for given cosmology ------------------------------------

  SUBROUTINE FDERIVS(n,x,y,dydx)
    IMPLICIT NONE
    INTEGER n
    REAL(8) :: x
    REAL(8), DIMENSION(2) :: y, dydx
    REAL(8) :: OMEGA_DE,w_de,z,f_OM,f_OR,f_Q,f

    z=exp(-x)-1.d0
    OMEGA_DE=1.d0-OMEGA_R-OMEGA_M
    f_OM=OMEGA_M*exp(-3.d0*x)*inv_hubblez(z)**2
    f_OR=OMEGA_R*exp(-4.d0*x)*inv_hubblez(z)**2
    f_Q=Omega_DE*f_DE(z)*inv_hubblez(z)**2
    w_de=w_0+w_a*(exp(x)-1.d0)
    f=y(2)     
    dydx(1)=f
    dydx(2)=-f**2-f/2.d0*(1.d0-f_OR-3.d0*w_de*f_Q)+3.D0/2.d0*f_OM    
  END SUBROUTINE FDERIVS
  
  SUBROUTINE GROWTH_FACTOR
    IMPLICIT NONE
    REAL(8), PARAMETER :: EPS=1.d-7
    INTEGER, PARAMETER :: NMAX=100, KMAXX=100
    INTEGER :: nok,nbad,nvar
    INTEGER :: IA
    REAL(8) :: h1,hmin
    REAL(8) :: XINI,AINI,XFIN,ZINI
    REAL(8), DIMENSION(2) :: y,dydx
    REAL(8) :: ALFAFIN,ALFASTART,DIN,DEND
    EXTERNAL rkqs
    REAL(8), DIMENSION(NZ) :: DPLUS_TEMP
    nvar=2
    h1=1.d-7
    hmin=0.d0
    ZINI=100.d0
    AINI=1.d0/(ZINI+1.d0)
    XINI=log(1.D0/(1.D0+ZINI))    
    XFIN=0.d0
    y(1)=log(3.d0/5.d0*AINI)
    y(2)=0.6*AINI/exp(y(1))
    DO IA=1,NZ
       ALFAFIN  = XINI + DBLE(IA)*(XFIN-XINI)/DBLE(NZ)
       ALFASTART = XINI + DBLE(IA-1)*(XFIN-XINI)/DBLE(NZ)
       call ODEINT(y,nvar,ALFASTART,ALFAFIN,EPS,h1,hmin,nok,nbad,FDERIVS,rkqs)
       DPLUS_TEMP(IA)=exp(y(1))
       AFIN(IA)=exp(ALFAFIN)
    END DO
    DO IA=1,NZ
       DPLUS(IA)=DPLUS_TEMP(IA)/DPLUS_TEMP(NZ)
    END DO
    DIN=1d+30
    DEND=1d+30
    CALL SPLINE(AFIN,DPLUS,NZ,DIN,DEND,D2DPLUS)
  END SUBROUTINE GROWTH_FACTOR


  FUNCTION DELTAC(OMEGAM,zred)
    IMPLICIT NONE
    REAL(8), PARAMETER :: pi=3.141592653589793
    REAL(8) :: DELTAC,OMEGAM,zred,OMEGAMZ

    OMEGAMZ=OMEGAM*(1.d0+zred)**3*inv_hubblez(zred)**2
    DELTAC=3.d0/20.d0*(12.d0*pi)**(2./3.)*(1.d0+0.0123*log10(OMEGAMZ))

  END FUNCTION DELTAC
  
END MODULE COSMO
