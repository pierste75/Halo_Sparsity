MODULE POWER_SPECTRUM

  USE COSMO

  INTEGER, PARAMETER :: NK=800,NM=800,NMRED=NM-2
  REAL(8), DIMENSION(NM) :: M,R  
  REAL(8), DIMENSION(NK) :: K_WAV,DELTA2K,D2DELTA2K,TK_MATTER
  REAL(8), DIMENSION(NMRED) :: SIGMA,DSIGMADR,RADIUS,D2SIGMA,D2DSIGMADR

CONTAINS

  FUNCTION SIGMA_R(R_Filter)
    IMPLICIT NONE
    REAL(8) :: SIGMA_R,R_Filter
    CALL NSPLINT(RADIUS,SIGMA,D2SIGMA,NMRED,R_Filter,SIGMA_R)
  END FUNCTION SIGMA_R


  FUNCTION DSIGMADR_R(R_Filter)
    IMPLICIT NONE
    REAL(8) :: DSIGMADR_R,R_Filter
    CALL NSPLINT(RADIUS,DSIGMADR,D2DSIGMADR,NMRED,R_Filter,DSIGMADR_R)
  END FUNCTION DSIGMADR_R


  SUBROUTINE MATTER_VARIANCE
    IMPLICIT NONE
    INTEGER :: IK,IM
    REAL(8) :: PKINI,PK
    REAL(8) :: DIN,DEND
    REAL(8) :: LG10M_MIN,LG10M_MAX

    CALL  TFfit(NK,OMEGA_M,OMEGA_B,h,K_WAV,TK_MATTER)

    DO IK=1,NK
       PKINI=(K_WAV(IK)*h)**(ANS-1.d0)
       PK=2*PI**2*PKINI*(K_WAV(IK)*h)*TK_MATTER(IK)**2*h**3
       DELTA2K(IK)=1./2.d0/PI**2*K_WAV(IK)**3*PK
    END DO

    DIN=1d+30
    DEND=1d+30
    CALL SPLINE(K_WAV,DELTA2K,NK,DIN,DEND,D2DELTA2K)

    LG10M_MIN=8.d0       ! M_MIN = 10^8 M_sun/h
    LG10M_MAX=19.d0      ! M_MAX = 10^19 M_sun/h

    DO IM=1,NM
       M(IM)=10**(LG10M_MIN+(LG10M_MAX-LG10M_MIN)*dble(IM-1)/dble(NM-1))
       R(IM)=(3.d0/4.d0/PI*M(IM)/RHO_M)**(1.d0/3.d0)   ! h^-1 Mpc
    END DO

    CALL VARIANCE(K_WAV(1),K_WAV(NK),NM,R,SIGMA_8,SIGMA,DSIGMADR)

    DO IM=2,NMRED
       RADIUS(IM-1)=R(IM)
    END DO

    DIN=1d+30
    DEND=1d+30
    CALL SPLINE(RADIUS,SIGMA,NMRED,DIN,DEND,D2SIGMA)
    CALL SPLINE(RADIUS,DSIGMADR,NMRED,DIN,DEND,D2DSIGMADR)

  END SUBROUTINE MATTER_VARIANCE

  SUBROUTINE VARIANCE(KMIN,KMAX,NRAD,RAD,SIGMA8,SIG,DSIGDR)
    IMPLICIT NONE
    INTEGER :: NRAD
    REAL(8), DIMENSION(NRAD) :: RAD,SIG2
    REAL(8), DIMENSION(NRAD-1) :: SIG,DSIGDR
    INTEGER, PARAMETER :: NKSTEP=1000
    REAL(8) :: ANORM,ANORM_MIN,ANORM_MAX,SIGMA8
    REAL(8) :: LN_KMIN,LN_KMAX,XMIN,XMAX,RADIUS,X,K
    REAL(8) :: DSIG2_MIN,DSIG2_MAX,DSIG2,KMIN,KMAX,DELTALNK
    REAL(8) :: WINDOW,DELTA2OFK
    INTEGER :: IK,IR

    LN_KMIN=log(KMIN)
    LN_KMAX=log(KMAX)
    DELTALNK=(LN_KMAX-LN_KMIN)/dble(NKSTEP)

    RADIUS=8.D0
    XMIN=RADIUS*KMIN
    XMAX=RADIUS*KMAX

    ANORM_MIN=WINDOW(XMIN)**2*DELTA2OFK(KMIN)
    ANORM_MAX=WINDOW(XMAX)**2*DELTA2OFK(KMAX)       
    ANORM=0.D0
    
    DO IK=1,NKSTEP
       K=exp(LN_KMIN+IK*DELTALNK)
       X=RADIUS*K
       ANORM=ANORM+WINDOW(X)**2*DELTA2OFK(K)
    END DO
    ANORM=DELTALNK*(ANORM+(ANORM_MIN+ANORM_MAX)/2.D0)

    ANORM=SIGMA8**2/ANORM
    
    DO IR=1,NRAD
       RADIUS=RAD(IR)
       XMIN=RADIUS*KMIN
       XMAX=RADIUS*KMAX
       DSIG2_MIN=WINDOW(XMIN)**2*ANORM*DELTA2OFK(KMIN)
       DSIG2_MAX=WINDOW(XMAX)**2*ANORM*DELTA2OFK(KMAX)       
       DSIG2=0.D0       
       DO IK=1,NKSTEP
          K=exp(LN_KMIN+IK*DELTALNK)
          X=RADIUS*K
          DSIG2=DSIG2+WINDOW(X)**2*ANORM*DELTA2OFK(K)
       END DO
       SIG2(IR)=DELTALNK*(DSIG2+(DSIG2_MIN+DSIG2_MAX)/2.D0)
    END DO
    
    DO IR=2,NRAD-1
       SIG(IR-1)=sqrt(SIG2(IR))
       DSIGDR(IR-1)=(sqrt(SIG2(IR+1))-sqrt(SIG2(IR-1)))/ &
            & (RAD(IR+1)-RAD(IR-1))/2.d0
    END DO

  END SUBROUTINE VARIANCE

END MODULE POWER_SPECTRUM


FUNCTION DELTA2OFK(K)
  USE POWER_SPECTRUM
  IMPLICIT NONE
  REAL(8) :: DELTA2OFK,K
  CALL NSPLINT(K_WAV,DELTA2K,D2DELTA2K,NK,K,DELTA2OFK)
END FUNCTION DELTA2OFK

FUNCTION WINDOW(X)
  IMPLICIT NONE
  REAL(8) WINDOW,X
  WINDOW=3.D0*(sin(X)-X*cos(X))/X**3
END FUNCTION WINDOW

