!@sum This file contains the radiation subroutines which don''t use
!@+   the module RADPAR. They are used by RADIATION and/or ALBEDO.
#include "rundeck_opts.h"
 
      MODULE GTAU_STATE_MOD
      SAVE 
      REAL*8  ::  GTAU(51,11,143)
      REAL*8  ::  TAUGSA(1001,14), SALBTG(768,14), TAUTGS(768),         &
                  TAUTGD(122)
      END MODULE GTAU_STATE_MOD

      SUBROUTINE RXSNOW(RBSNO,XCOSZ,GGSNO,RXSNO)
!@sum RXSNOW calculate zenith angle dependence for snow/ice albedo
!@auth A. Lacis (modified by G. Schmidt)
  !    USE RADPAR, only : gtsalb,sgpgxg
      IMPLICIT NONE
!@var RBSNO diffuse albedo
      REAL*8, INTENT(IN)  ::  RBSNO
!@var XCOSZ zenith angle
      REAL*8, INTENT(IN)  ::  XCOSZ
!@var GGSNO Asymmetry parameter for snow
      REAL*8, INTENT(IN)  ::  GGSNO
!@var RXSNO direct albedo
      REAL*8, INTENT(OUT)  ::  RXSNO
      INTEGER NDBLS, NN
      REAL*8 XXG, XXT, GGSN, RBSN, FRTOP, TAU, TAUSN, GPFF, PR, PT,     &
             DBLS, SECZ, XANB, XANX, TANB, TANX, RASB, RASX, BNORM,     &
             XNORM, RARB, RARX, XATB, DENOM, DB, DX, UB, UX, DRBRAT,    &
             RBBOUT
 
      IF ( RBSNO<0.05D0 ) THEN
         RXSNO = RBSNO
         RETURN
      ENDIF
      XXG = 0.D0
      XXT = 0.D0
      GGSN = GGSNO
      IF ( GGSNO>0.9D0 ) GGSN = 0.9D0
      RBSN = RBSNO
      FRTOP = 1.D0
      IF ( RBSNO>0.5D0 ) THEN
         RBSN = 0.5D0
         FRTOP = ((1.D0-RBSNO)/0.5D0)**2
      ENDIF
 
      CALL GTSALB(XXG,XXT,RBBOUT,RBSN,GGSN,TAUSN,2)
      CALL SGPGXG(XCOSZ,TAUSN,GGSN,GPFF)
      PR = 1.D0 - GPFF
      PT = 1.D0 + GPFF
      DBLS = 10.D0 + 1.44269D0*LOG(TAUSN)
      NDBLS = DBLS
      TAU = TAUSN/2**NDBLS
!     Set optically thin limit values of R,T,X using PI0 renormalization
!     ------------------------------------------------------------------
!
      SECZ = 1.D0/XCOSZ
      XANB = EXP(-TAU-TAU)
      XANX = EXP(-TAU*SECZ)
      TANB = PT*XANB
      XXT = (SECZ-2.D0)*TAU
      TANX = PT*SECZ*                                                   &
             (.5D0+XXT*(.25D0+XXT*(.0833333D0+XXT*(.0208333D0+XXT))))   &
             *XANX
      RASB = PR*(1.D0-TAU*(2.D0-2.66667D0*TAU*(1.D0-TAU)))
      XXT = (SECZ+2.D0)*TAU
      RASX = PR*SECZ*                                                   &
             (.5D0-XXT*(.25D0-XXT*(.0833333D0-XXT*(.0208333D0-XXT))))
      BNORM = (1.D0-XANB)/(RASB+TANB)
      XNORM = (1.D0-XANX)/(RASX+TANX)
      RASB = RASB*BNORM
      RASX = RASX*XNORM
      TANB = TANB*BNORM
      TANX = TANX*XNORM
      DO NN = 1, NDBLS
         RARB = RASB*RASB
         RARX = XANX*RASX
         XATB = XANB + TANB
         DENOM = 1.D0 - RARB
         DB = (TANB+XANB*RARB)/DENOM
         DX = (TANX+RARX*RASB)/DENOM
         UB = RASB*(XANB+DB)
         UX = RARX + RASB*DX
         RASB = RASB + XATB*UB
         RASX = RASX + XATB*UX
         TANB = XANB*TANB + XATB*DB
         TANX = XANX*TANX + XATB*DX
         XANB = XANB*XANB
         XANX = XANX*XANX
      ENDDO
      DRBRAT = RASX/RBSN - 1.D0
      RXSNO = RBSNO*(1.D0+DRBRAT*FRTOP)
      END SUBROUTINE RXSNOW 
 
      SUBROUTINE SETGTS(tgdata_in)
      USE GTAU_STATE_MOD
      IMPLICIT NONE
      REAL*8, INTENT(IN)  ::  tgdata_in(122,13)
      REAL*8 CWM, CWE, TIJ, RBB, RBBI, BTAU
      INTEGER I, J
 
      REAL*8  ::  tgdata(122,13)
                               ! why cant we just
      TGDATA = tgdata_in       ! pass tgdata_in to spline
 
      DO I = 1, 122
         TAUTGD(I) = (I-1)*0.1D0
         IF ( I>24 ) TAUTGD(I) = (I-24)*0.2D0 + 2.2D0
         IF ( I>48 ) TAUTGD(I) = (I-48)*0.5D0 + 7.0D0
         IF ( I>72 ) TAUTGD(I) = (I-72) + 19.0D0
         IF ( I>96 ) TAUTGD(I) = (I-96)*5.0D0 + 40.0D0
         IF ( I>112 ) TAUTGD(I) = (I-112)*100.0D0 + 100.0D0
         IF ( I==121 ) TAUTGD(I) = 9999.99D0
         IF ( I==122 ) TAUTGD(I) = 12000.0D0
      ENDDO
 
      DO I = 1, 768
         IF ( I<602 ) TAUTGS(I) = (I-1)*0.05D0
         IF ( I>601 ) TAUTGS(I) = (I-601)*0.50D0 + 30.0D0
         IF ( I>741 ) TAUTGS(I) = (I-741)*50.0D0 + 100.D0
         IF ( I>758 ) TAUTGS(I) = (I-758)*1000.D0
      ENDDO
 
      DO J = 1, 13
         DO I = 1, 768
            CWM = 0.5
            CWE = 0.5
            IF ( I>759 ) CWM = 0.0
            IF ( I>759 ) CWE = 0.0
            TIJ = TAUTGS(I)
            CALL SPLINE(TAUTGD,TGDATA(1,J),122,TIJ,RBBI,CWM,CWE,0)
            SALBTG(I,J) = RBBI
         ENDDO
      ENDDO
      DO J = 1, 13
         DO I = 2, 1000
            RBB = (I-1)*0.001D0
            CWM = 0.5
            CWE = 0.5
            CALL SPLINE(SALBTG(1,J),TAUTGS,768,RBB,BTAU,CWM,CWE,0)
            TAUGSA(I,J) = BTAU
         ENDDO
      ENDDO
      SALBTG(1,:) = 0                                ! 1:14
      TAUGSA(1,:) = 0
      TAUGSA(1001,:) = 10000
 
      SALBTG(:,14) = SALBTG(:,13)*2 - SALBTG(:,12)   ! 1:768
      TAUGSA(:,14) = TAUGSA(:,13)*2 - TAUGSA(:,12)   ! 1:1001
 
      END SUBROUTINE SETGTS
 
      SUBROUTINE GTSALB(GIN,TAUIN,RBBOUT,RBBIN,EGIN,TAUOUT,KGTAUR)
      USE GTAU_STATE_MOD
      IMPLICIT NONE
      REAL*8, INTENT(IN)  ::  GIN, TAUIN, RBBIN, EGIN
      INTEGER, INTENT(IN)  ::  KGTAUR
      REAL*8, INTENT(OUT)  ::  RBBOUT, TAUOUT
 
      REAL*8 FFKG(4,3), RBBK(3)
      REAL*8, PARAMETER, DIMENSION(14)                                  &
                                    ::  GVALUE = (/.0,.25,.45,.50,.55,  &
                                   .60,.65,.70,.75,.80,.85,.90,.95,1./)
      REAL*8 RBB, G, TAU, EG, DELTAU, TI, WTJ, WTI, GI, WGI, WGJ, F1,   &
             F2, F3, F4, A, B, C, RB2, RB3, TBB, TB2, TB3, XG, XM, XP,  &
             RBBB, RI, WRJ, WRI, EI, WEI, WEJ, DELALB, X1, X2, X3, X4,  &
             XX, BB, DTAU
      INTEGER K, KTERPL, IT, JT, IG, JG, ITERPL, IGM, JGP, KG, IR, JR,  &
              IE, JE, IEM, JEP
      REAL*8, EXTERNAL  ::  COMPUTE
      KTERPL = 0
 
      G = GIN
      TAU = TAUIN
      RBB = RBBIN
      EG = EGIN
 
      RBBOUT = 0.0
      TAUOUT = 0.0
!                                           ---------------------------
!                                           OPTICAL DEPTH INTERPOLATION
!                                           0.05 ON (0.00 < TAU < 30.0)
!                                           0.50 ON (30.0 < TAU < 100.)
!                                           50.0 ON (100. < TAU < 1000)
!                                           ---------------------------
 
      IF ( KGTAUR==2 ) GOTO 300
 
 
 200  DELTAU = 0.05D0
      TI = TAU/DELTAU
      IT = TI
      IF ( IT>599 ) THEN
         DELTAU = 0.50D0
         TI = TAU/DELTAU
         IT = TI
         IF ( IT>199 ) THEN
            DELTAU = 50.0D0
            TI = TAU/DELTAU
            IT = TI
            IF ( IT>19 ) THEN
               DELTAU = 1000.0D0
               TI = TAU/DELTAU
               IT = TI
               WTJ = TI - IT
               WTI = 1.0 - WTJ
               IT = IT + 758
            ELSE
               WTJ = TI - IT
               WTI = 1.0 - WTJ
               IT = IT + 649
            ENDIF
         ELSE
            WTJ = TI - IT
            WTI = 1.0 - WTJ
            IT = IT + 541
         ENDIF
      ELSE
         WTJ = TI - IT
         WTI = 1.D0 - WTJ
         IT = IT + 1
      ENDIF
      JT = IT + 1
 
!                                    ---------------------------------
!                                    ASYMMETRY PARAMETER INTERPOLATION
!                                    0.05 CUBIC SPLINE (0.5 < G < 0.9)
!                                    0.25 QUADRATIC ON (0.0 < G < 0.5)
!                                    LINEAR EXTRAP FOR (.95 < G < 1.0)
!                                    ---------------------------------
 
      GI = G*20.D0
      IF ( GI>10.0 ) THEN
 
         ITERPL = 4
         IG = GI
         WGJ = GI - IG
         WGI = 1.D0 - WGJ
         IG = IG - 6
         IF ( IG>12 ) THEN
            ITERPL = 2
            IG = 12
         ENDIF
         JG = IG + 1
      ELSE
         IG = 2
         JG = 3
         ITERPL = 1
      ENDIF
 
 
      IGM = IG - 1
      JGP = JG + 1
 
      K = 0
      DO KG = IGM, JGP
         K = K + 1
         F1 = SALBTG(IT-1,KG)
         F2 = SALBTG(IT,KG)
         F3 = SALBTG(JT,KG)
         F4 = SALBTG(JT+1,KG)
         IF ( IT==1 ) F1 = -F3
         FFKG(K,1) = COMPUTE(F1,F2,F3,F4,WTJ)
         FFKG(K,2) = F2
         FFKG(K,3) = F3
      ENDDO
 
      IF ( ITERPL<4 ) THEN
 
         XG = G*2.D0 - 0.5D0
         IF ( ITERPL==2 ) XG = G*10.D0 - 9.D0
         XM = 1.D0 - XG - XG
         XP = 1.D0 + XG + XG
         RBB = XM*XP*FFKG(ITERPL+1,1) - XG*XM*FFKG(ITERPL,1)            &
               + XG*XP*FFKG(4,1)
         RB2 = XM*XP*FFKG(ITERPL+1,2) - XG*XM*FFKG(ITERPL,2)            &
               + XG*XP*FFKG(4,2)
         RB3 = XM*XP*FFKG(ITERPL+1,3) - XG*XM*FFKG(ITERPL,3)            &
               + XG*XP*FFKG(4,3)
 
         IF ( KGTAUR==1 ) RETURN
         IF ( KTERPL==1 ) GOTO 400
      ELSE
 
         DO K = 1, 3
            F1 = FFKG(1,K)
            F2 = FFKG(2,K)
            F3 = FFKG(3,K)
            F4 = FFKG(4,K)
            RBBK(K) = COMPUTE(F1,F2,F3,F4,WGJ)
         ENDDO
         RBB = RBBK(1)
         RB2 = RBBK(2)
         RB3 = RBBK(3)
         TBB = TAU
         TB2 = TAUTGS(IT)
         TB3 = TAUTGS(JT)
         IF ( KGTAUR==1 ) RETURN
         IF ( KTERPL==1 ) GOTO 400
      ENDIF
 
 300  RBBB = RBB
 
      RI = RBB*1000.D0
      IR = RI
      WRJ = RI - IR
      WRI = 1.D0 - WRJ
      IR = IR + 1
      JR = IR + 1
 
      EI = EG*20.D0
      IF ( EI>10.0 ) THEN
 
         ITERPL = 4
         IE = EI
         WEJ = EI - IE
         WEI = 1.D0 - WEJ
         IE = IE - 6
         IF ( IE>12 ) THEN
            ITERPL = 2
            IE = 12
         ENDIF
         JE = IE + 1
      ELSE
         IE = 2
         JE = 3
         ITERPL = 1
      ENDIF
 
      DELALB = 0.001D0
      IEM = IE - 1
      JEP = JE + 1
      K = 0
      DO KG = IEM, JEP
         K = K + 1
         F1 = TAUGSA(IR-1,KG)
         F2 = TAUGSA(IR,KG)
         F3 = TAUGSA(JR,KG)
         F4 = TAUGSA(JR+1,KG)
         IF ( IR==1 ) F1 = -F3
         FFKG(K,1) = COMPUTE(F1,F2,F3,F4,WRJ)
         FFKG(K,2) = F2
         FFKG(K,3) = F3
      ENDDO
      X1 = GVALUE(IE-1)
      X2 = GVALUE(IE)
      X3 = GVALUE(JE)
      X4 = GVALUE(JE+1)
      XX = WEJ
      IF ( ITERPL<4 ) THEN
 
         XG = EG*2.D0 - 0.5D0
         IF ( ITERPL==2 ) XG = G*10.D0 - 9.D0
         XM = 1.D0 - XG - XG
         XP = 1.D0 + XG + XG
         TBB = XM*XP*FFKG(ITERPL+1,1) - XG*XM*FFKG(ITERPL,1)            &
               + XG*XP*FFKG(4,1)
         TB2 = XM*XP*FFKG(ITERPL+1,2) - XG*XM*FFKG(ITERPL,2)            &
               + XG*XP*FFKG(4,2)
         TB3 = XM*XP*FFKG(ITERPL+1,3) - XG*XM*FFKG(ITERPL,3)            &
               + XG*XP*FFKG(4,3)
         IF ( KTERPL==1 ) GOTO 400
      ELSE
 
         DO K = 1, 3
            F1 = FFKG(1,K)
            F2 = FFKG(2,K)
            F3 = FFKG(3,K)
            F4 = FFKG(4,K)
            RBBK(K) = COMPUTE(F1,F2,F3,F4,WEJ)
         ENDDO
         TBB = RBBK(1)
         TB2 = RBBK(2)
         TB3 = RBBK(3)
 
         IF ( KTERPL==1 ) GOTO 400
      ENDIF
      KTERPL = 1
      TAU = TBB
      G = EGIN
      GOTO 200
 400  IF ( ABS(WTI*WTJ)<0.1D0 ) DTAU = (RBBB-RB2)/(RB3-RB2)
      IF ( ABS(WTI*WTJ)>=0.1D0 ) THEN
         C = (RB3-RBB)/WTI - (RBB-RB2)/WTJ
         B = (RBB-RB2)/WTJ - WTJ*C
         A = RB2
         BB = B*B + 4.D0*C*(RBBB-A)
         IF ( BB>0.D0 ) DTAU = (SQRT(BB)-B)/(C+C)
      ENDIF
      TAUOUT = (IT-1+DTAU)*DELTAU
      RBBOUT = RBBB
 
 
      END SUBROUTINE GTSALB
 
      SUBROUTINE SGPGXG(XMU,TAU,G,GG)
      USE GTAU_STATE_MOD
      IMPLICIT NONE
!     ----------------------------------------------------------------
!     COSBAR ADJUSTMENT TO REPRODUCE THE SOLAR ZENITH ANGLE DEPENDENCE
!     FOR AEROSOL ALBEDO FOR OPTICAL THICKNESSES [0.0 < TAU < 10000.0]
!     ----------------------------------------------------------------
      REAL*8, INTENT(IN)  ::  XMU, TAU, G
      REAL*8, INTENT(OUT)  ::  GG
      REAL*8 XI, WXI, WXJ, GI, WGI, WGJ, TI, WTJ, WTI
      INTEGER IX, JX, IG, JG, IT, IT0, JT
!                          -------------------------------------------
!                          XMU (COSZ) SOLAR ZENITH ANGLE INTERPOLATION
!                          DATA INTERVAL:  0.02  ON  [0.0 < XMU < 1.0]
!                          -------------------------------------------
 
      XI = XMU*50.D0 + 0.999999D0
                               ! >1 since XMU=COSZ>.001
      IX = XI
      JX = IX + 1
      WXJ = XI - IX
      WXI = 1.D0 - WXJ
 
!                                      -------------------------------
!                                      COSBAR DEPENDENCE INTERPOLATION
!                                         0.10 ON [0.0 < COSBAR < 1.0]
!                                      -------------------------------
 
      GI = G*10.D0
      IG = GI
      WGJ = GI - IG
      WGI = 1.D0 - WGJ
      IG = IG + 1
      JG = IG + 1
 
!                            -----------------------------------------
!                               AEROSOL TAU INTERPOLATION INTERVALS
!                            -----------------------------------------
!                            dTau      1      1  (Lin Int)  61     62
!                            0.10 ON [0.00 , 0.00 < TAU < 6.00 , 6.10]
!                                     63     64            92     93
!                            0.50 ON [5.50 , 6.00 < TAU < 20.0 , 20.5]
!                                     94     95           111    112
!                            5.00 ON [15.0 , 20.0 < TAU < 100. , 105.]
!                                     113    114          132    133
!                            50.0 ON [50.0 , 100. < TAU < 1000 , 1050]
!                                            134          143
!                            1000 ON [     , 1000 < TAU < 10000,     ]
!                            -----------------------------------------
 
      IF ( TAU<6.D0 ) THEN
         TI = TAU*10.D0 + 1.
         IT = TI
         WTJ = TI - IT
         IT0 = 0
 
      ELSEIF ( TAU<20.D0 ) THEN
         TI = (TAU-6.D0)*2.00D0 + 2.D0
         IT = TI
         WTJ = TI - IT
         IT0 = 62
 
      ELSEIF ( TAU<100.D0 ) THEN
         TI = (TAU-20.D0)*0.20D0 + 2.D0
         IT = TI
         WTJ = TI - IT
         IT0 = 93
 
      ELSEIF ( TAU<1000.D0 ) THEN
         TI = (TAU-100.D0)*0.02D0 + 2.D0
         IT = TI
         WTJ = TI - IT
         IT0 = 112
 
      ELSE
         TI = TAU*0.001D0 + 1.D-6
         IT = TI
         WTJ = TI - IT
         IF ( IT>9 ) IT = 9
         IT0 = 133
      ENDIF
 
      WTI = 1.D0 - WTJ
      IT = IT + IT0
      JT = IT + 1
      GG = WGI*(WTI*(WXI*GTAU(IX,IG,IT)+WXJ*GTAU(JX,IG,IT))             &
           +WTJ*(WXI*GTAU(IX,IG,JT)+WXJ*GTAU(JX,IG,JT)))                &
           + WGJ*(WTI*(WXI*GTAU(IX,JG,IT)+WXJ*GTAU(JX,JG,IT))           &
           +WTJ*(WXI*GTAU(IX,JG,JT)+WXJ*GTAU(JX,JG,JT)))
 
      RETURN
 
      END SUBROUTINE SGPGXG
 
      SUBROUTINE SET_SGPGXG(GTAU_IN)
      USE GTAU_STATE_MOD
      REAL*8, INTENT(IN)  ::  GTAU_IN(51,11,143)
      GTAU = GTAU_IN
      END SUBROUTINE SET_SGPGXG
 
      SUBROUTINE SPLINE(X,F,NXF,XX,FF,CUSPWM,CUSPWE,KXTRAP)
      IMPLICIT NONE
 
      INTEGER, INTENT(IN)  ::  NXF, KXTRAP
      REAL*8, INTENT(IN)  ::  X(NXF), F(NXF), XX, CUSPWM, CUSPWE
      REAL*8, INTENT(OUT)  ::  FF
      REAL*8  ::  FFVEC(1)
 
!---------------------------------------------------------------------
!
!    SPLINE locates XX between points (F2,X2)(F3,X3) on 4-point spread
!       and returns 4-point Cubic Spline interpolated value FF = F(XX)
!
!    Quadratic Derivatives of Spline are continuous at (F2,X2),(F3,X3)
!    (X-Coordinate may be specified in increasing or decreasing order)
!
!---------------------------------------------------------------------
!
!    Custom Control Parameters:  CUSPWM,CUSPWE,KXTRAP
!------------------------------
!
!    In cases where data points are unevenly spaced and/or data points
!    exhibit abrupt changes in value, Spline Interpolation may produce
!    undesirable bulging of interpolated values. In more extreme cases
!    Linear Interpolation may be less problematic to use.
!
!    Interpolation can be weighted between: Cubic Spline and Linear by
!    adjusting weights CUSPWM and CUSPWE to values between 1.0 and 0.0
!
!    CUSPWM = Cubic Spline Weight at the (X2-X3) Interval Mid-point
!    CUSPWE = Cubic Spline Weight at the (X2-X3) Interval End-points
!
!    For example, with:
!
!    CUSPWM=1.0,CUSPWE=1.0  FF returns Cubic Spline interpolated value
!    CUSPWM=0.0,CUSPWE=0.0  FF returns   Linearly   interpolated value
!
!---------------------------------------------------------------------
!
!     Extrapolation for XX outside of defined interval:  X(1)<->X(NXF)
!
!               KXTRAP = 0    No Extrapolation  (i.e., sets F(XX)=0.0)
!                        1    Fixed Extrapolation (F(XX) = edge value)
!                        2    Linear Extrapolation using 2 edge points
!
!---------------------------------------------------------------------
 
      FFVEC(1) = FF
      CALL SPLINEVECTOR(X,F,1,NXF,XX,FFVEC,CUSPWM,CUSPWE,KXTRAP)
      FF = FFVEC(1)
 
      END SUBROUTINE SPLINE
 
      SUBROUTINE SPLINEVECTOR(X,F,NVEC,NXF,XX,FF,CUSPWM,CUSPWE,KXTRAP)
      IMPLICIT NONE
 
      INTEGER, INTENT(IN)  ::  NVEC, NXF, KXTRAP
      REAL*8, INTENT(IN)  ::  X(NXF), F(NVEC,NXF), XX, CUSPWM, CUSPWE
      REAL*8, INTENT(OUT)  ::  FF(NVEC)
 
!---------------------------------------------------------------------
!
!    SPLINEVector is identical to SPLINE, except operates on vector functions
!    rather than scalar functions.  More efficient than calling in a loop.
!
!---------------------------------------------------------------------
 
      REAL*8 x1, x2, x3, x4, x21, x32, x43, x31, x42, betw, CUSPWT
      REAL*8, DIMENSION(NVEC)  ::  f1, f2, f3, f4
      REAL*8, DIMENSION(NVEC)  ::  f21, f32, f43, f3221, f4332
      REAL*8, DIMENSION(NVEC)  ::  A, B, C, D, FFCUSP, FFLINR
      REAL*8 xf, xe, xexm
      INTEGER K
 
      K = 2
      X2 = X(K)
      X3 = X(NXF-1)
      BETW = (XX-X2)*(X3-XX)
      IF ( BETW<=0.D0 ) THEN
 
!                Edge Point Interval Interpolation and/or Extrapolation
!                ------------------------------------------------------
         BETW = (X2-XX)*(X3-X2)
         IF ( BETW<0.D0 ) THEN
 
!                    X(NXF-1),X(NXF)  Edge Point Interval Interpolation
!                    --------------------------------------------------
            F3 = F(:,NXF)
            X3 = X(NXF)
            F2 = F(:,NXF-1)
            X2 = X(NXF-1)
            X32 = X3 - X2
            F32 = (F3-F2)/X32
            XF = XX - X3
            BETW = (X2-XX)*(XX-X3)
            IF ( BETW<0.D0 ) THEN
 
!              Extrapolation for X Outside of Interval  X(NXF-1)-X(NXF)
!              --------------------------------------------------------
!                  IF(KXTRAP == 0)  (No Extrapolation:  sets F(XX)=0.0)
!                  IF(KXTRAP == 1)  (Extrapolation at Fixed Edge Value)
!                  IF(KXTRAP == 2)  (2 Edge Point Linear Extrapolation)
 
               IF ( KXTRAP==0 ) FF = 0.D0
               IF ( KXTRAP==1 ) FF = F3
               IF ( KXTRAP==2 ) FF = F3 + XF*(F3-F2)/(X3-X2)
            ELSE
               F1 = F(:,NXF-2)
               X1 = X(NXF-2)
               X21 = X2 - X1
               X31 = X3 - X1
               F21 = (F2-F1)/X21
               XF = XX - X2
 
!                    3-Point Quadratic Interpolation for Edge Intervals
!                    --------------------------------------------------
!
!      (Edge Option)     ----------------------------------------------
!                        For Linear Interpolation within Edge Intervals
!                        between X(1),X(2), and between X(NXF-1),X(NXF)
!                        set the value of coefficient C below, to C=0.0
!                        ----------------------------------------------
 
               C = (F32-F21)/X31
               B = F21 + X21*C
               A = F2
               FFCUSP = A + XF*(B+XF*C)
               FFLINR = A + XF*F32
               XE = 1.D0 - 2.D0*XF/X32
               IF ( XE<0.D0 ) XE = -XE
               XEXM = XE**2
               CUSPWT = (1.D0-XEXM)*CUSPWM + XEXM*CUSPWE
               FF = FFCUSP*CUSPWT + FFLINR*(1.D0-CUSPWT)
            ENDIF
         ELSE
 
!                          X(1),X(2)  Edge Point Interval Interpolation
!                          --------------------------------------------
            X1 = X(1)
            F1 = F(:,1)
            F2 = F(:,2)
            X21 = X2 - X1
            F21 = (F2-F1)/X21
            XF = XX - X1
            BETW = (X2-XX)*XF
            IF ( BETW<0.D0 ) THEN
 
!                  Extrapolation for XX Outside of Interval X(1) - X(2)
!                  ----------------------------------------------------
!                  IF(KXTRAP == 0)  (No Extrapolation:  sets F(XX)=0.0)
!                  IF(KXTRAP == 1)  (Extrapolation at Fixed Edge Value)
!                  IF(KXTRAP == 2)  (2 Edge Point Linear Extrapolation)
 
               IF ( KXTRAP==0 ) FF = 0.D0
               IF ( KXTRAP==1 ) FF = F1
               IF ( KXTRAP==2 ) FF = F1 + XF*F21
            ELSE
               F3 = F(:,3)
               X3 = X(3)
               X32 = X3 - X2
               X31 = X3 - X1
               C = ((F3-F2)/X32-F21)/X31
               B = F21 - X21*C
               A = F1
               FFCUSP = A + XF*(B+XF*C)
               FFLINR = A + XF*F21
               XE = 1.D0 - 2.D0*XF/X21
               IF ( XE<0.D0 ) XE = -XE
               XEXM = XE**2
               CUSPWT = (1.D0-XEXM)*CUSPWM + XEXM*CUSPWE
               FF = FFCUSP*CUSPWT + FFLINR*(1.D0-CUSPWT)
            ENDIF
         ENDIF
      ELSE
         DO
 
            K = K + 1
            X3 = X(K)
            BETW = (XX-X2)*(X3-XX)
            IF ( BETW>=0.D0 ) THEN
 
               F3(:) = F(:,K)
               F4(:) = F(:,K+1)
               X4 = X(K+1)
               F2(:) = F(:,K-1)
               X2 = X(K-1)
               F1(:) = F(:,K-2)
               X1 = X(K-2)
               X21 = X2 - X1
               X31 = X3 - X1
               X32 = X3 - X2
               X43 = X4 - X3
               X42 = X4 - X2
               F21(:) = (F2(:)-F1(:))/(X21*X21)
               F32(:) = (F3(:)-F2(:))/(X32*X32)
               F43(:) = (F4(:)-F3(:))/(X43*X43)
               F3221(:) = (F32(:)+F21(:))/X31*X21
               F4332(:) = (F43(:)+F32(:))/X42*X43
               A = F2
               B = X32*F3221
               C = 3.D0*F32 - F3221 - F3221 - F4332
               D = (F3221+F4332-F32-F32)/X32
               XF = XX - X2
 
!                             FFCUSP= Cubic Spline Interpolation Result
!                             -----------------------------------------
 
               FFCUSP = A + XF*(B+XF*(C+XF*D))
               XE = (X3+X2-XX-XX)/X32
               IF ( XE<0.D0 ) XE = -XE
               XEXM = XE**2
               CUSPWT = (1.D0-XEXM)*CUSPWM + XEXM*CUSPWE
 
!                                   FFLINR= Linear Interpolation Result
!                                   -----------------------------------
               FFLINR = A + XF*F32*X32
               FF = FFCUSP*CUSPWT + FFLINR*(1.D0-CUSPWT)
               EXIT
            ELSE
               X2 = X3
            ENDIF
         ENDDO
      ENDIF
 
      END SUBROUTINE SPLINEVECTOR 
!cc the following subroutines were just moved from RADIATION.f to
!cc reduce its size. Only MODULE RADPAR subroutines or those that
!cc USE RADPAR module were left.
 
      SUBROUTINE BOXAV1(DEGLAT,TAULAT,NLAT,JALIM,JBLIM,TAU)
      IMPLICIT NONE
!
!--------------------------------------------------------------------
!     BOXAV1  Performs:
!                       Latitudinal average (area-weighted) of TAULAT
!
!              DEGLAT   Center latitude of grid-box variable (TAULAT)
!                       of the form:  DEGLAT = -90+(J-1)*180/(NLAT-1)
!
!              TAULAT   Zonal average value is constant over grid-bos
!
!        JALIM, JBLIM   Latitude boxes for which (TAULAT) is averaged
!
!                 TAU   Area-weighted (TAULAT) latitude average value
!--------------------------------------------------------------------
!
      INTEGER, INTENT(IN)  ::  NLAT, JALIM, JBLIM
      REAL*8, DIMENSION(NLAT), INTENT(IN)  ::  DEGLAT, TAULAT
      REAL*8, INTENT(OUT)  ::  TAU
      REAL*8  ::  ASUM, TSUM
      REAL*8  ::  ONES(NLAT)
 
      ONES = 1.0D0
      CALL BOXAV(DEGLAT,ONES,TAULAT,NLAT,JALIM,JBLIM,TSUM,ASUM)
      TAU = TSUM/ASUM
 
      END SUBROUTINE BOXAV1
 
      SUBROUTINE BOXAV2(DEGLAT,TAULAT,SIZLAT,NLAT,JALIM,JBLIM,SIZ)
      IMPLICIT NONE
!
!--------------------------------------------------------------------
!     BOXAV2  Performs:
!                       TAULAT-weighted latitudinal average of SIZLAT
!
!              DEGLAT   Center latitude of grid-box variable (TAULAT)
!                       of the form:  DEGLAT = -90+(J-1)*180/(NLAT-1)
!
!              TAULAT   Zonal average value is constant over grid-box
!              SIZLAT   Zonal average value is constant over grid-box
!
!        JALIM, JBLIM   Latitude boxes for which variable is averaged
!
!                 SIZ   TAULAT-weighted latitudinal average of SIZLAT
!--------------------------------------------------------------------
!
      INTEGER, INTENT(IN)  ::  NLAT, JALIM, JBLIM
      REAL*8, DIMENSION(NLAT), INTENT(IN)  ::  DEGLAT, TAULAT, SIZLAT
      REAL*8, INTENT(OUT)  ::  SIZ
      REAL*8 ASUM, TSUM
 
      CALL BOXAV(DEGLAT,TAULAT,SIZLAT,NLAT,JALIM,JBLIM,TSUM,ASUM)
      SIZ = (1.D-20+TSUM)/(1.D-10+ASUM)
 
      END SUBROUTINE BOXAV2
 
      SUBROUTINE BOXAV(DEGLAT,W1,ARR,NLAT,JALIM,JBLIM,TSUM,ASUM)
      IMPLICIT NONE
!
!--------------------------------------------------------------------
!     BOXAV  Performs:
!                       W1 weighted sums of ARR
!
!              DEGLAT   Center latitude of grid-box variable (W1)
!                       of the form:  DEGLAT = -90+(J-1)*180/(NLAT-1)
!
!              W1       Zonal average value is constant over grid-box
!              SIZLAT   Zonal average value is constant over grid-box
!
!        JALIM, JBLIM   Latitude boxes for which variable is averaged
!
!--------------------------------------------------------------------
!
      INTEGER, INTENT(IN)  ::  NLAT, JALIM, JBLIM
      REAL*8, DIMENSION(NLAT), INTENT(IN)  ::  DEGLAT, W1, ARR
      REAL*8, INTENT(OUT)  ::  TSUM, ASUM
      REAL*8 PI, RADIAN, RLAT1, RLAT2, ALAT1, ALAT2, ALATJ
      INTEGER J, J1, J2
 
      ASUM = 0.D0
      TSUM = 0.D0
      PI = ACOS(-1.D0)
      RADIAN = 180.D0/PI
      J1 = JALIM - 1
      IF ( J1<1 ) J1 = 1
      RLAT1 = (0.5D0*(DEGLAT(J1)+DEGLAT(JALIM))+90.D0)/RADIAN
      ALAT1 = SIN(RLAT1)
      DO J = JALIM, JBLIM
         J2 = J + 1
         IF ( J2>NLAT ) J2 = NLAT
         RLAT2 = (0.5D0*(DEGLAT(J)+DEGLAT(J2))+90.D0)/RADIAN
         ALAT2 = SIN(RLAT2)
         ALATJ = 0.5D0*(ALAT1+ALAT2)/(RLAT2-RLAT1)
         ASUM = ASUM + ALATJ*W1(J)
         TSUM = TSUM + ALATJ*W1(J)*ARR(J)
         RLAT1 = RLAT2
         ALAT1 = ALAT2
      ENDDO
      END SUBROUTINE BOXAV
 
      SUBROUTINE PHATMO(P,H,D,T,O,Q,S,OCM,WCM,NPHD,NATM)
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     -------------     MCCLATCHY (1972) ATMOSPHERE DATA     -----------
!     ------------------------------------------------------------------
!
!        INPUT DATA
!------------------
!                  NATM=0  GIVES ABREVIATED DATA FOR  STANDARD ATMOSPHER
!                                (INPUT: P OR H) (RETURNS: H OR P   D,T)
!
!                  NATM=1  GIVES ATMOSPHERE DATA FOR  TROPICAL LATITUDES
!                  NATM=2  GIVES ATMOSPHERE DATA FOR  MIDLATITUDE SUMMER
!                  NATM=3  GIVES ATMOSPHERE DATA FOR  MIDLATITUDE WINTER
!                  NATM=4  GIVES ATMOSPHERE DATA FOR  SUBARCTIC SUMMER
!                  NATM=5  GIVES ATMOSPHERE DATA FOR  SUBARCTIC WINTER
!                  NATM=6  GIVES ATMOSPHERE DATA FOR  STANDARD ATMOSPHER
!
!                  NPHD=1  RETURNS H,D,T,O,Q,S DATA FOR GIVEN PRESSURE P
!                  NPHD=2  RETURNS P,D,T,O,Q,S DATA FOR GIVEN   HEIGHT H
!                  NPHD=3  RETURNS P,H,T,O,Q,S DATA FOR GIVEN  DENSITY D
!
!       OUTPUT DATA
!------------------
!                  P = PRESSURE IN MILLIBARS
!                  H = HEIGHT IN KILOMETERS
!                  D = DENSITY IN GRAMS/METER**3
!                  T = TEMPERATURE (ABSOLUTE)
!                  O = OZONE MIXING RATIO (GRAMS OZONE)/(GRAMS AIR)
!                  Q = SPECIFIC HUMIDITY (GRAMS WATER VAPOR)/(GRAMS AIR)
!                  S = SATURATION RATIO (GRAMS WATER VAPOR)/(GRAMS AIR)
!                  OCM = OZONE (CM-STP) ABOVE GIVEN HEIGHT
!                  WCM = WATER VAPOR (CM-STP) ABOVE GIVEN HEIGHT
!
!           REMARKS
!------------------
!                  INPUT P,H,D PARAMETERS ARE NOT ALTERED
!                  P,D INTERPOLATION IS EXPONENTIAL WITH HEIGHT
!                  NO EXTRAPOLATION IS MADE OUTSIDE 0-100 KM INTERVAL
!                  S  IS NOT COMPUTED ABOVE 40 KM (FORMULA NOT ACCURATE)
!
!                  R = Q/S          GIVES RELATIVE HUMIDITY
!                  W = Q/(1-Q)      GIVES WATER VAPOR MIXING RATIO
!                  N = D*2.079E 16  GIVES NUMBER DENSITY PER CM**3
!
      REAL*8, DIMENSION(33)  ::  PRS1, PRS2, PRS3, PRS4, PRS5, PRS6,    &
                                 DNS1, DNS2, DNS3, DNS4, DNS5, DNS6,    &
                                 TMP1, TMP2, TMP3, TMP4, TMP5, TMP6,    &
                                 WVP1, WVP2, WVP3, WVP4, WVP5, WVP6,    &
                                 OZO1, OZO2, OZO3, OZO4, OZO5, OZO6
      REAL*8, DIMENSION(33,6)  ::  PRES, DENS, TEMP, WVAP, OZON
 
      EQUIVALENCE (PRES(1,1),PRS1(1))
      EQUIVALENCE (DENS(1,1),DNS1(1))
      EQUIVALENCE (TEMP(1,1),TMP1(1))
      EQUIVALENCE (PRES(1,2),PRS2(1))
      EQUIVALENCE (DENS(1,2),DNS2(1))
      EQUIVALENCE (TEMP(1,2),TMP2(1))
      EQUIVALENCE (PRES(1,3),PRS3(1))
      EQUIVALENCE (DENS(1,3),DNS3(1))
      EQUIVALENCE (TEMP(1,3),TMP3(1))
      EQUIVALENCE (PRES(1,4),PRS4(1))
      EQUIVALENCE (DENS(1,4),DNS4(1))
      EQUIVALENCE (TEMP(1,4),TMP4(1))
      EQUIVALENCE (PRES(1,5),PRS5(1))
      EQUIVALENCE (DENS(1,5),DNS5(1))
      EQUIVALENCE (TEMP(1,5),TMP5(1))
      EQUIVALENCE (PRES(1,6),PRS6(1))
      EQUIVALENCE (DENS(1,6),DNS6(1))
      EQUIVALENCE (TEMP(1,6),TMP6(1))
      EQUIVALENCE (WVAP(1,1),WVP1(1))
      EQUIVALENCE (OZON(1,1),OZO1(1))
      EQUIVALENCE (WVAP(1,2),WVP2(1))
      EQUIVALENCE (OZON(1,2),OZO2(1))
      EQUIVALENCE (WVAP(1,3),WVP3(1))
      EQUIVALENCE (OZON(1,3),OZO3(1))
      EQUIVALENCE (WVAP(1,4),WVP4(1))
      EQUIVALENCE (OZON(1,4),OZO4(1))
      EQUIVALENCE (WVAP(1,5),WVP5(1))
      EQUIVALENCE (OZON(1,5),OZO5(1))
      EQUIVALENCE (WVAP(1,6),WVP6(1))
      EQUIVALENCE (OZON(1,6),OZO6(1))
 
      REAL*8, PARAMETER, DIMENSION(33)                                  &
                                    ::  HTKM = (/1D-9,1D0,2D0,3D0,4D0,  &
                                   5D0,6D0,7D0,8D0,9D0,10D0,11D0,12D0,  &
                                   13D0,14D0,15D0,16D0,17D0,18D0,19D0,  &
                                   20D0,21D0,22D0,23D0,24D0,25D0,30D0,  &
                                   35D0,40D0,45D0,50D0,70D0,99.9D0/)
 
!----------------------------------------------------------------------
!0000 GLOBAL   U.S. (1976) STANDARD ATMOSPHERE   P, T, GEO H  PARAMETERS
!----------------------------------------------------------------------
 
      REAL*8, PARAMETER, DIMENSION(8)  ::  SPLB = (/1013.25D0,226.32D0, &
                                   54.748D0,8.6801D0,1.109D0,.66938D0,  &
                                   .039564D0,3.7338D-03/),              &
                                   STLB = (/288.15D0,216.65D0,216.65D0, &
                                   228.65D0,270.65D0,270.65D0,214.65D0, &
                                   186.87D0/),                          &
                                   SHLB = (/0D0,11D0,20D0,32D0,47D0,    &
                                   51D0,71D0,84.852D0/),                &
                                   SDLB = (/-6.5D0,0D0,1D0,2.8D0,0D0,   &
                                   -2.8D0,-2D0,0D0/)
      REAL*8, PARAMETER  ::  HPCON = 34.16319D0
 
!-----------------------------------------------------------------------
!1111 TROPICAL LATITUDES      MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT
!-----------------------------------------------------------------------
 
      DATA PRS1/1.013D03, 9.040D02, 8.050D02, 7.150D02, 6.330D02,       &
           5.590D02, 4.920D02, 4.320D02, 3.780D02, 3.290D02, 2.860D02,  &
           2.470D02, 2.130D02, 1.820D02, 1.560D02, 1.320D02, 1.110D02,  &
           9.370D01, 7.890D01, 6.660D01, 5.650D01, 4.800D01, 4.090D01,  &
           3.500D01, 3.000D01, 2.570D01, 1.220D01, 6.000D00, 3.050D00,  &
           1.590D00, 8.540D-01, 5.790D-02, 3.000D-04/
      DATA DNS1/1.167D03, 1.064D03, 9.689D02, 8.756D02, 7.951D02,       &
           7.199D02, 6.501D02, 5.855D02, 5.258D02, 4.708D02, 4.202D02,  &
           3.740D02, 3.316D02, 2.929D02, 2.578D02, 2.260D02, 1.972D02,  &
           1.676D02, 1.382D02, 1.145D02, 9.515D01, 7.938D01, 6.645D01,  &
           5.618D01, 4.763D01, 4.045D01, 1.831D01, 8.600D00, 4.181D00,  &
           2.097D00, 1.101D00, 9.210D-02, 5.000D-04/
      DATA TMP1/300.0, 294.0, 288.0, 284.0, 277.0, 270.0, 264.0, 257.0, &
           250.0, 244.0, 237.0, 230.0, 224.0, 217.0, 210.0, 204.0,      &
           197.0, 195.0, 199.0, 203., 207.0, 211.0, 215.0, 217.0, 219.0,&
           221.0, 232.0, 243.0, 254.0, 265.0, 270., 219.0, 210.0/
      DATA WVP1/1.9D01, 1.3D01, 9.3D00, 4.7D00, 2.2D00, 1.5D00, 8.5D-01,&
           4.7D-01, 2.5D-01, 1.2D-01, 5.0D-02, 1.7D-02, 6.0D-03,        &
           1.8D-03, 1.0D-03, 7.6D-04, 6.4D-04, 5.6D-04, 5.0D-04,        &
           4.9D-04, 4.5D-04, 5.1D-04, 5.1D-04, 5.4D-04, 6.0D-04,        &
           6.7D-04, 3.6D-04, 1.1D-04, 4.3D-05, 1.9D-05, 6.3D-06,        &
           1.4D-07, 1.0D-09/
      DATA OZO1/5.6D-05, 5.6D-05, 5.4D-05, 5.1D-05, 4.7D-05, 4.5D-05,   &
           4.3D-05, 4.1D-05, 3.9D-05, 3.9D-05, 3.9D-05, 4.1D-05,        &
           4.3D-05, 4.5D-05, 4.5D-05, 4.7D-05, 4.7D-05, 6.9D-05,        &
           9.0D-05, 1.4D-04, 1.9D-04, 2.4D-04, 2.8D-04, 3.2D-04,        &
           3.4D-04, 3.4D-04, 2.4D-04, 9.2D-05, 4.1D-05, 1.3D-05,        &
           4.3D-06, 8.6D-08, 4.3D-11/
 
!-----------------------------------------------------------------------
!2222 MIDLATITUDE SUMMER      MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT
!-----------------------------------------------------------------------
 
      DATA PRS2/1.013D03, 9.020D02, 8.020D02, 7.100D02, 6.280D02,       &
           5.540D02, 4.870D02, 4.260D02, 3.720D02, 3.240D02, 2.810D02,  &
           2.430D02, 2.090D02, 1.790D02, 1.530D02, 1.300D02, 1.110D02,  &
           9.500D01, 8.120D01, 6.950D01, 5.950D01, 5.100D01, 4.370D01,  &
           3.760D01, 3.220D01, 2.770D01, 1.320D01, 6.520D00, 3.330D00,  &
           1.760D00, 9.510D-01, 6.710D-02, 3.000D-04/
      DATA DNS2/1.191D03, 1.080D03, 9.757D02, 8.846D02, 7.998D02,       &
           7.211D02, 6.487D02, 5.830D02, 5.225D02, 4.669D02, 4.159D02,  &
           3.693D02, 3.269D02, 2.882D02, 2.464D02, 2.104D02, 1.797D02,  &
           1.535D02, 1.305D02, 1.110D02, 9.453D01, 8.056D01, 6.872D01,  &
           5.867D01, 5.014D01, 4.288D01, 1.322D01, 6.519D00, 3.330D00,  &
           1.757D00, 9.512D-01, 6.706D-02, 5.000D-04/
      DATA TMP2/294.0, 290.0, 285.0, 279.0, 273.0, 267.0, 261.0, 255.0, &
           248.0, 242.0, 235.0, 229.0, 222.0, 216.0, 216.0, 216.0,      &
           216.0, 216.0, 216.0, 217., 218.0, 219.0, 220.0, 222.0, 223.0,&
           224.0, 234.0, 245.0, 258.0, 270.0, 276., 218.0, 210.0/
      DATA WVP2/1.4D01, 9.3D00, 5.9D00, 3.3D00, 1.9D00, 1.0D00, 6.1D-01,&
           3.7D-01, 2.1D-01, 1.2D-01, 6.4D-02, 2.2D-02, 6.0D-03,        &
           1.8D-03, 1.0D-03, 7.6D-04, 6.4D-04, 5.6D-04, 5.0D-04,        &
           4.9D-04, 4.5D-04, 5.1D-04, 5.1D-04, 5.4D-04, 6.0D-04,        &
           6.7D-04, 3.6D-04, 1.1D-04, 4.3D-05, 1.9D-05, 6.3D-06,        &
           1.4D-07, 1.0D-09/
      DATA OZO2/6.0D-05, 6.0D-05, 6.0D-05, 6.2D-05, 6.4D-05, 6.6D-05,   &
           6.9D-05, 7.5D-05, 7.9D-05, 8.6D-05, 9.0D-05, 1.1D-04,        &
           1.2D-04, 1.5D-04, 1.8D-04, 1.9D-04, 2.1D-04, 2.4D-04,        &
           2.8D-04, 3.2D-04, 3.4D-04, 3.6D-04, 3.6D-04, 3.4D-04,        &
           3.2D-04, 3.0D-04, 2.0D-04, 9.2D-05, 4.1D-05, 1.3D-05,        &
           4.3D-06, 8.6D-08, 4.3D-11/
 
!-----------------------------------------------------------------------
!3333 MIDLATITUDE WINTER      MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT
!-----------------------------------------------------------------------
 
      DATA PRS3/1.018D03, 8.973D02, 7.897D02, 6.938D02, 6.081D02,       &
           5.313D02, 4.627D02, 4.016D02, 3.473D02, 2.992D02, 2.568D02,  &
           2.199D02, 1.882D02, 1.610D02, 1.378D02, 1.178D02, 1.007D02,  &
           8.610D01, 7.350D01, 6.280D01, 5.370D01, 4.580D01, 3.910D01,  &
           3.340D01, 2.860D01, 2.430D01, 1.110D01, 5.180D00, 2.530D00,  &
           1.290D00, 6.820D-01, 4.670D-02, 3.000D-04/
      DATA DNS3/1.301D03, 1.162D03, 1.037D03, 9.230D02, 8.282D02,       &
           7.411D02, 6.614D02, 5.886D02, 5.222D02, 4.619D02, 4.072D02,  &
           3.496D02, 2.999D02, 2.572D02, 2.206D02, 1.890D02, 1.620D02,  &
           1.388D02, 1.188D02, 1.017D02, 8.690D01, 7.421D01, 6.338D01,  &
           5.415D01, 4.624D01, 3.950D01, 1.783D01, 7.924D00, 3.625D00,  &
           1.741D00, 8.954D-01, 7.051D-02, 5.000D-04/
      DATA TMP3/272.2, 268.7, 265.2, 261.7, 255.7, 249.7, 243.7, 237.7, &
           231.7, 225.7, 219.7, 219.2, 218.7, 218.2, 217.7, 217.2,      &
           216.7, 216.2, 215.7, 215.2, 215.2, 215.2, 215.2, 215.2,      &
           215.2, 215.2, 217.4, 227.8, 243.2, 258.5, 265.7, 230.7,      &
           210.2/
      DATA WVP3/3.5D00, 2.5D00, 1.8D00, 1.2D00, 6.6D-01, 3.8D-01,       &
           2.1D-01, 8.5D-02, 3.5D-02, 1.6D-02, 7.5D-03, 6.9D-03,        &
           6.0D-03, 1.8D-03, 1.0D-03, 7.6D-04, 6.4D-04, 5.6D-04,        &
           5.0D-04, 4.9D-04, 4.5D-04, 5.1D-04, 5.1D-04, 5.4D-04,        &
           6.0D-04, 6.7D-04, 3.6D-04, 1.1D-04, 4.3D-05, 1.9D-05,        &
           6.3D-06, 1.4D-07, 1.0D-09/
      DATA OZO3/6.0D-05, 5.4D-05, 4.9D-05, 4.9D-05, 4.9D-05, 5.8D-05,   &
           6.4D-05, 7.7D-05, 9.0D-05, 1.2D-04, 1.6D-04, 2.1D-04,        &
           2.6D-04, 3.0D-04, 3.2D-04, 3.4D-04, 3.6D-04, 3.9D-04,        &
           4.1D-04, 4.3D-04, 4.5D-04, 4.3D-04, 4.3D-04, 3.9D-04,        &
           3.6D-04, 3.4D-04, 1.9D-04, 9.2D-05, 4.1D-05, 1.3D-05,        &
           4.3D-06, 8.6D-08, 4.3D-11/
 
!-----------------------------------------------------------------------
!4444 SUBARCTIC SUMMER        MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT
!-----------------------------------------------------------------------
 
      DATA PRS4/1.010D03, 8.960D02, 7.929D02, 7.000D02, 6.160D02,       &
           5.410D02, 4.730D02, 4.130D02, 3.590D02, 3.107D02, 2.677D02,  &
           2.300D02, 1.977D02, 1.700D02, 1.460D02, 1.250D02, 1.080D02,  &
           9.280D01, 7.980D01, 6.860D01, 5.890D01, 5.070D01, 4.360D01,  &
           3.750D01, 3.227D01, 2.780D01, 1.340D01, 6.610D00, 3.400D00,  &
           1.810D00, 9.870D-01, 7.070D-02, 3.000D-04/
      DATA DNS4/1.220D03, 1.110D03, 9.971D02, 8.985D02, 8.077D02,       &
           7.244D02, 6.519D02, 5.849D02, 5.231D02, 4.663D02, 4.142D02,  &
           3.559D02, 3.059D02, 2.630D02, 2.260D02, 1.943D02, 1.671D02,  &
           1.436D02, 1.235D02, 1.062D02, 9.128D01, 7.849D01, 6.750D01,  &
           5.805D01, 4.963D01, 4.247D01, 1.338D01, 6.614D00, 3.404D00,  &
           1.817D00, 9.868D-01, 7.071D-02, 5.000D-04/
      DATA TMP4/287.0, 282.0, 276.0, 271.0, 266.0, 260.0, 253.0, 246.0, &
           239.0, 232.0, 225.0, 225.0, 225.0, 225.0, 225.0, 225.0,      &
           225.0, 225.0, 225.0, 225., 225.0, 225.0, 225.0, 225.0, 226.0,&
           228.0, 235.0, 247.0, 262.0, 274.0, 277., 216.0, 210.0/
      DATA WVP4/9.1D00, 6.0D00, 4.2D00, 2.7D00, 1.7D00, 1.0D00, 5.4D-01,&
           2.9D-01, 1.3D-02, 4.2D-02, 1.5D-02, 9.4D-03, 6.0D-03,        &
           1.8D-03, 1.0D-03, 7.6D-04, 6.4D-04, 5.6D-04, 5.0D-04,        &
           4.9D-04, 4.5D-04, 5.1D-04, 5.1D-04, 5.4D-04, 6.0D-04,        &
           6.7D-04, 3.6D-04, 1.1D-04, 4.3D-05, 1.9D-05, 6.3D-06,        &
           1.4D-07, 1.0D-09/
      DATA OZO4/4.9D-05, 5.4D-05, 5.6D-05, 5.8D-05, 6.0D-05, 6.4D-05,   &
           7.1D-05, 7.5D-05, 7.9D-05, 1.1D-04, 1.3D-04, 1.8D-04,        &
           2.1D-04, 2.6D-04, 2.8D-04, 3.2D-04, 3.4D-04, 3.9D-04,        &
           4.1D-04, 4.1D-04, 3.9D-04, 3.6D-04, 3.2D-04, 3.0D-04,        &
           2.8D-04, 2.6D-04, 1.4D-04, 9.2D-05, 4.1D-05, 1.3D-05,        &
           4.3D-06, 8.6D-08, 4.3D-11/
 
!-----------------------------------------------------------------------
!5555 SUBARCTIC WINTER        MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT
!-----------------------------------------------------------------------
 
      DATA PRS5/1.013D03, 8.878D02, 7.775D02, 6.798D02, 5.932D02,       &
           5.158D02, 4.467D02, 3.853D02, 3.308D02, 2.829D02, 2.418D02,  &
           2.067D02, 1.766D02, 1.510D02, 1.291D02, 1.103D02, 9.431D01,  &
           8.058D01, 6.882D01, 5.875D01, 5.014D01, 4.277D01, 3.647D01,  &
           3.109D01, 2.649D01, 2.256D01, 1.020D01, 4.701D00, 2.243D00,  &
           1.113D00, 5.719D-01, 4.016D-02, 3.000D-04/
      DATA DNS5/1.372D03, 1.193D03, 1.058D03, 9.366D02, 8.339D02,       &
           7.457D02, 6.646D02, 5.904D02, 5.226D02, 4.538D02, 3.879D02,  &
           3.315D02, 2.834D02, 2.422D02, 2.071D02, 1.770D02, 1.517D02,  &
           1.300D02, 1.113D02, 9.529D01, 8.155D01, 6.976D01, 5.966D01,  &
           5.100D01, 4.358D01, 3.722D01, 1.645D01, 7.368D00, 3.330D00,  &
           1.569D00, 7.682D-01, 5.695D-02, 5.000D-04/
      DATA TMP5/257.1, 259.1, 255.9, 252.7, 247.7, 240.9, 234.1, 227.3, &
           220.6, 217.2, 217.2, 217.2, 217.2, 217.2, 217.2, 217.2,      &
           216.6, 216., 215.4, 214.8, 214.1, 213.6, 213.0, 212.4, 211.8,&
           211.2, 216.0, 222.2, 234.7, 247., 259.3, 245.7, 210.0/
      DATA WVP5/1.2D00, 1.2D00, 9.4D-01, 6.8D-01, 4.1D-01, 2.0D-01,     &
           9.8D-02, 5.4D-02, 1.1D-02, 8.4D-03, 5.5D-03, 3.8D-03,        &
           2.6D-03, 1.8D-03, 1.0D-03, 7.6D-04, 6.4D-04, 5.6D-04,        &
           5.0D-04, 4.9D-04, 4.5D-04, 5.1D-04, 5.1D-04, 5.4D-04,        &
           6.0D-04, 6.7D-04, 3.6D-04, 1.1D-04, 4.3D-05, 1.9D-05,        &
           6.3D-06, 1.4D-07, 1.0D-09/
      DATA OZO5/4.1D-05, 4.1D-05, 4.1D-05, 4.3D-05, 4.5D-05, 4.7D-05,   &
           4.9D-05, 7.1D-05, 9.0D-05, 1.6D-04, 2.4D-04, 3.2D-04,        &
           4.3D-04, 4.7D-04, 4.9D-04, 5.6D-04, 6.2D-04, 6.2D-04,        &
           6.2D-04, 6.0D-04, 5.6D-04, 5.1D-04, 4.7D-04, 4.3D-04,        &
           3.6D-04, 3.2D-04, 1.5D-04, 9.2D-05, 4.1D-05, 1.3D-05,        &
           4.3D-06, 8.6D-08, 4.3D-11/
 
!----------------------------------------------------------------------
!6666 GLOBAL   U.S. (1976) STANDARD ATMOSPHERE   P, T, GEO H  PARAMETERS
!----------------------------------------------------------------------
 
      DATA PRS6/1.01325D+03, 8.987D+02, 7.950D+02, 7.011D+02, 6.164D+02,&
           5.402D+02, 4.718D+02, 4.106D+02, 3.560D+02, 3.074D+02,       &
           2.644D+02, 2.263D+02, 1.933D+02, 1.651D+02, 1.410D+02,       &
           1.204D+02, 1.029D+02, 8.787D+01, 7.505D+01, 6.410D+01,       &
           5.475D+01, 4.678D+01, 4.000D+01, 3.422D+01, 2.931D+01,       &
           2.511D+01, 1.172D+01, 5.589D+00, 2.775D+00, 1.431D+00,       &
           7.594D-01, 4.634D-02, 2.384D-04/
      DATA DNS6/1.225D+03, 1.112D+03, 1.006D+03, 9.091D+02, 8.191D+02,  &
           7.361D+02, 6.597D+02, 5.895D+02, 5.252D+02, 4.663D+02,       &
           4.127D+02, 3.639D+02, 3.108D+02, 2.655D+02, 2.268D+02,       &
           1.937D+02, 1.654D+02, 1.413D+02, 1.207D+02, 1.031D+02,       &
           8.803D+01, 7.487D+01, 6.373D+01, 5.428D+01, 4.627D+01,       &
           3.947D+01, 1.801D+01, 8.214D+00, 3.851D+00, 1.881D+00,       &
           9.775D-01, 7.424D-02, 4.445D-04/
      DATA TMP6/288.150, 281.650, 275.150, 268.650, 262.150, 255.650,   &
           249.150, 242.650, 236.150, 229.650, 223.150, 216.650,        &
           216.650, 216.650, 216.650, 216.650, 216.650, 216.650,        &
           216.650, 216.650, 216.650, 217.650, 218.650, 219.650,        &
           220.650, 221.650, 226.650, 237.050, 251.050, 265.050,        &
           270.650, 217.450, 186.870/
      DATA WVP6/1.083D+01, 6.323D+00, 3.612D+00, 2.015D+00, 1.095D+00,  &
           5.786D-01, 2.965D-01, 1.469D-01, 7.021D-02, 3.226D-02,       &
           1.419D-02, 5.956D-03, 5.002D-03, 4.186D-03, 3.490D-03,       &
           2.896D-03, 2.388D-03, 1.954D-03, 1.583D-03, 1.267D-03,       &
           9.967D-04, 8.557D-04, 7.104D-04, 5.600D-04, 4.037D-04,       &
           2.406D-04, 5.404D-05, 2.464D-05, 1.155D-05, 5.644D-06,       &
           2.932D-06, 2.227D-07, 1.334D-09/
      DATA OZO6/7.526D-05, 3.781D-05, 6.203D-05, 3.417D-05, 5.694D-05,  &
           3.759D-05, 5.970D-05, 4.841D-05, 7.102D-05, 6.784D-05,       &
           9.237D-05, 9.768D-05, 1.251D-04, 1.399D-04, 1.715D-04,       &
           1.946D-04, 2.300D-04, 2.585D-04, 2.943D-04, 3.224D-04,       &
           3.519D-04, 3.714D-04, 3.868D-04, 3.904D-04, 3.872D-04,       &
           3.728D-04, 2.344D-04, 9.932D-05, 3.677D-05, 1.227D-05,       &
           4.324D-06, 5.294D-08, 1.262D-10/
 
      REAL*8, INTENT(INOUT)  ::  H, P, D
      INTEGER, INTENT(IN)  ::  NATM, NPHD
      REAL*8, INTENT(OUT)  ::  O, Q, S, OCM, WCM, T
      REAL*8  ::  XX, XI, XJ, DELTA, RAT, PI, PJ, DI, DJ, DP, ES, RS,   &
                  OI, OJ, QI, QJ
      INTEGER  ::  I, J, K, N
 
      IF ( NATM>0 ) THEN
 
         IF ( NPHD/=1 ) THEN
            IF ( NPHD/=2 ) THEN
               XX = D
               XI = DENS(1,NATM)
               IF ( D>XI ) XX = XI
               IF ( D<5.0E-04 ) GOTO 280
               DO J = 2, 33
                  XJ = DENS(J,NATM)
                  IF ( XX>XJ ) GOTO 260
                  XI = XJ
               ENDDO
            ENDIF
            XX = H
            XI = HTKM(1)
            IF ( H<XI ) XX = XI
            IF ( H>99.9 ) GOTO 280
            DO J = 2, 33
               XJ = HTKM(J)
               IF ( XX<XJ ) GOTO 260
               XI = XJ
            ENDDO
         ENDIF
         XX = P
         XI = PRES(1,NATM)
         IF ( P>XI ) XX = XI
         IF ( P<3.0E-04 ) GOTO 280
         DO J = 2, 33
            XJ = PRES(J,NATM)
            IF ( XX>XJ ) EXIT
            XI = XJ
         ENDDO
      ELSE
         O = 1.E-10
         Q = 1.E-10
         S = 1.E-10
         OCM = 1.E-10
         WCM = 1.E-10
         IF ( NPHD<2 ) THEN
 
            DO N = 2, 8
               IF ( P>SPLB(N) ) GOTO 170
            ENDDO
            N = 9
 170        N = N - 1
            IF ( ABS(SDLB(N))<1.E-04 ) THEN
               H = SHLB(N) + STLB(N)/HPCON*LOG(SPLB(N)/P)
            ELSE
               H = SHLB(N) + STLB(N)/SDLB(N)                            &
                   *((SPLB(N)/P)**(SDLB(N)/HPCON)-1.)
            ENDIF
            T = STLB(N) + SDLB(N)*(H-SHLB(N))
            D = P/T*28.9644E05/8.31432E03
            RETURN
         ELSE
            DO N = 2, 8
               IF ( H<SHLB(N) ) GOTO 120
            ENDDO
            N = 9
 120        N = N - 1
            IF ( ABS(SDLB(N))<1.E-04 ) THEN
               P = SPLB(N)*EXP(-HPCON/STLB(N)*(H-SHLB(N)))
            ELSE
               P = SPLB(N)*(1.+SDLB(N)/STLB(N)*(H-SHLB(N)))             &
                   **(-HPCON/SDLB(N))
            ENDIF
            T = STLB(N) + SDLB(N)*(H-SHLB(N))
            D = P/T*28.9644E05/8.31432E03
                                       ! P/T*mair/gasc
            RETURN
         ENDIF
      ENDIF
 260  DELTA = (XX-XI)/(XJ-XI)
      I = J - 1
      IF ( NPHD/=2 ) THEN
         RAT = LOG(XX/XI)/LOG(XJ/XI)
         H = HTKM(I) + (HTKM(J)-HTKM(I))*RAT
         T = TEMP(I,NATM) + (TEMP(J,NATM)-TEMP(I,NATM))*RAT
      ENDIF
      PI = PRES(I,NATM)
      PJ = PRES(J,NATM)
      DI = DENS(I,NATM)
      DJ = DENS(J,NATM)
      IF ( NPHD==1 ) D = DI + DELTA*(DJ-DI)
      IF ( NPHD==3 ) P = PI + DELTA*(PJ-PI)
      IF ( NPHD==2 ) THEN
         P = PI*(PJ/PI)**DELTA
         D = DI*(DJ/DI)**DELTA
         T = TEMP(I,NATM) + DELTA*(TEMP(J,NATM)-TEMP(I,NATM))
      ENDIF
      O = OZON(I,NATM)/DI + DELTA*(OZON(J,NATM)/DJ-OZON(I,NATM)/DI)
      Q = WVAP(I,NATM)/DI + DELTA*(WVAP(J,NATM)/DJ-WVAP(I,NATM)/DI)
      ES = 10.D0**(9.4051D0-2353.D0/T)
      IF ( P<PI ) PI = P
      S = 1.D+06
      RS = (PI-ES+0.622*ES)/(0.622*ES)
      IF ( RS>1.E-06 ) S = 1./RS
      OI = O
      QI = Q
      OCM = 0.D0
      WCM = 0.D0
      DO K = J, 33
         PJ = PRES(K,NATM)
         DJ = DENS(K,NATM)
         OJ = OZON(K,NATM)/DJ
         QJ = WVAP(K,NATM)/DJ
         DP = PI - PJ
         OCM = OCM + 0.5D0*(OI+OJ)*DP
         WCM = WCM + 0.5D0*(QI+QJ)*DP
         OI = OJ
         QI = QJ
         PI = PJ
      ENDDO
      WCM = WCM/0.980D0*22420.7D0/18.D0
      OCM = OCM/0.980D0*22420.7D0/48.D0
      RETURN
 280  T = 210.D0
      IF ( NATM==6 ) T = 186.87
      O = 1.D-10
      Q = 1.D-10
      S = 1.D-10
      OCM = 1.D-10
      WCM = 1.D-10
      IF ( NPHD/=1 ) P = 1.D-05
      IF ( NPHD/=2 ) H = 99.99
      IF ( NPHD/=3 ) D = 2.D-05
      END SUBROUTINE PHATMO
 
      REAL*8 FUNCTION PFOFTK(WAVNA,WAVNB,TK)
!     ------------------------------------------------------------------
!
!        INPUT DATA
!                  WAVNA,WAVNB  SPECTRAL INTERVAL IN WAVENUMBERS
!                               (ORDER OF WAVNA,WAVNB NOT IMPORTANT)
!
!                  TK           ABSOLUTE TEMPERATURE IN DEGREES KELVIN
!
!       OUTPUT DATA
!                  PFofTK       PLANCK FLUX (W/m^2)
!
!
!           REMARKS
!                   PLANCK INTENSITY (W/m^2*STER) IS GIVEN BY PFofTK/PI
!
!     ------------------------------------------------------------------
      USE CONSTANT, ONLY:stbo  ! (W m-2 K-4) Stefan-Boltzmann
      IMPLICIT NONE
      REAL*8, PARAMETER, DIMENSION(21)                                  &
                                    ::  BN = (/1D0,-1D0,1D0,-1D0,1D0,   &
                                   -1D0,5D0,-691D0,7D0,-3617D0,43867D0, &
                                   -174611D0,854513D0,-236364091D0,     &
                                   8553103D0,-23749461029D0,            &
                                   8615841276005D0,-7709321041217D0,    &
                                   2577687858367D0,-2631527155305348D4, &
                                   2929993913841559D0/),                &
                                   BD = (/1D0,2D0,6D0,30D0,42D0,30D0,   &
                                   66D0,2730D0,6D0,510D0,798D0,330D0,   &
                                   138D0,2730D0,6D0,870D0,14322D0,510D0,&
                                   6D0,1919190D0,6D0/)
      REAL*8, PARAMETER  ::  PI4 = 97.40909103400244D0
!     REAL*8, PARAMETER :: PI =3.141592653589793D0
      REAL*8, PARAMETER  ::  HCK = 1.43879D0
      REAL*8, PARAMETER  ::  DGXLIM = 1D-06
 
      REAL*8, INTENT(IN)  ::  WAVNA, WAVNB, TK
      REAL*8 GSUM, B, DG, DGB, GX, PNORM, GTERM, GXA, GXB, X, XX, XN,   &
             XN3, XNN, XNM, XNF, XNX
      INTEGER II, NB, NNB, N
 
      PFOFTK = 0D0
      IF ( TK<1D-06 ) RETURN
      DO II = 1, 2
         IF ( II==1 ) X = HCK*WAVNA/TK
         IF ( II==2 ) X = HCK*WAVNB/TK
         IF ( X>2.3D0 ) THEN
            GSUM = PI4/15.D0
            DO N = 1, 20
               NNB = N
               XN = N
               XNN = XN*XN
               XNX = XN*X
               IF ( XNX>100.D0 ) EXIT
               GTERM = (X*X*(3.D0+XNX)+6.D0*(1.D0+XNX)/XNN)/XNN
               DG = GTERM*EXP(-XNX)
               GSUM = GSUM - DG
               DGB = DG
               IF ( DG<DGXLIM ) EXIT
            ENDDO
            GX = GSUM
         ELSE
            XX = X*X
            GSUM = 1D0/3D0 - X/8D0 + XX/60D0
            NB = 3
            XNF = XX/2D0
            DO N = 4, 38, 2
               NB = NB + 1
               NNB = NB
               B = BN(NB)/BD(NB)
               XN3 = N + 3
               XNM = N*(N-1)
               XNF = XNF*(XX/XNM)
               DG = B/XN3*XNF
               GSUM = GSUM + DG
               DGB = DG
               IF ( ABS(DG)<DGXLIM ) EXIT
            ENDDO
            GX = GSUM*XX*X
         ENDIF
         IF ( II==1 ) GXA = GX
         IF ( II==2 ) GXB = GX
      ENDDO
      PNORM = 15.D0/PI4
      PFOFTK = ABS(GXB-GXA)*PNORM
      PFOFTK = PFOFTK*stbo*TK**4
      END FUNCTION PFOFTK
 
      REAL*8 FUNCTION TKOFPF(WAVNA,WAVNB,FLUXAB)
!     ------------------------------------------------------------------
!
!        INPUT DATA
!------------------
!                  WAVNA,WAVNB  SPECTRAL INTERVAL IN WAVENUMBERS
!                               (ORDER OF WAVNA,WAVNB NOT IMPORTANT)
!                  FLUXAB       PLANCK FLUX (W/m^2) IN INTERVAL
!                                                       (WAVNA,WAVNB)
!
!       OUTPUT DATA
!------------------
!                  TK           BRIGHTNESS TEMPERATURE IN DEGREES KELVIN
!
!
!           REMARKS
!------------------
!                   TKofPF IS INVERSE FUNCTION OF PFofTK(WAVNA,WAVNB,TK)
!                   THE OUTPUT OF TKofPF SATISFIES THE IDENTITY
!                                 FLUXAB=PFofTK(WAVNA,WAVNB,TK)
!                   (UNITS FOR FLUXAB AND PFofTK MUST BE IDENTICAL)
!
!     ------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL LOGFIT
      REAL*8, PARAMETER  ::  DELFIT = 1.D-06
      INTEGER, PARAMETER  ::  NMAX = 20
      REAL*8, INTENT(IN)  ::  WAVNA, WAVNB, FLUXAB
      REAL*8 XA, YA, XB, YB, XX, YY, XC, YC, XBA, XCA, XBC, YBA, YCA,   &
             YBC, YXBA, YXCA, A, B, C, ROOT, PF, PFOFTK
      INTEGER NFIT
 
      IF ( FLUXAB<=0.D0 ) RETURN
      LOGFIT = .FALSE.
      NFIT = 0
      PF = FLUXAB
      XA = 0.D0
      YA = 0.D0
      XB = 250.D0
      YB = PFOFTK(WAVNA,WAVNB,XB)
      XX = PF*XB/YB
      YY = PFOFTK(WAVNA,WAVNB,XX)
      IF ( ABS(YY-PF)<DELFIT ) THEN
         TKOFPF = XX
         GOTO 99999
      ELSEIF ( (YY/PF)<0.5D0 ) THEN
         XA = XX
         YA = YY
         DO
            XC = XB
            YC = YB
            XB = XB/2.D0
            YB = PFOFTK(WAVNA,WAVNB,XB)
            IF ( YB<YA ) GOTO 190
            IF ( YB<=PF ) THEN
               XA = XB
               YA = YB
               GOTO 190
            ENDIF
         ENDDO
      ELSEIF ( (YY/PF)>2.0D0 ) THEN
         XC = XX
         YC = YY
         DO
            XA = XB
            YA = YB
            XB = XB*2.D0
            YB = PFOFTK(WAVNA,WAVNB,XB)
            IF ( YB>YC ) GOTO 190
            IF ( YB>=PF ) THEN
               XC = XB
               YC = YB
               GOTO 190
            ENDIF
         ENDDO
      ELSEIF ( XX>XB ) THEN
         XC = XX
         YC = YY
      ELSE
         XC = XB
         YC = YB
         XB = XX
         YB = YY
      ENDIF
 120  DO
         XBA = XB - XA
         XCA = XC - XA
         XBC = XB - XC
         YBA = YB - YA
         YCA = YC - YA
         YBC = YB - YC
         NFIT = NFIT + 1
         IF ( NFIT>NMAX ) THEN
            TKOFPF = XX
            GOTO 99999
         ELSE
            YXBA = YBA/XBA
            YXCA = YCA/XCA
            C = (YXBA-YXCA)/XBC
            B = YXBA - (XB+XA)*C
            A = YA - XA*(B+XA*C)
            ROOT = SQRT(B*B+4.D0*C*(PF-A))
            XX = 0.5D0*(ROOT-B)/C
            IF ( XX<XA .OR. XX>XC ) XX = -0.5D0*(ROOT+B)/C
            YY = PFOFTK(WAVNA,WAVNB,XX)
            IF ( LOGFIT ) YY = LOG(YY)
            IF ( ABS(YY-PF)<DELFIT ) THEN
               TKOFPF = XX
               GOTO 99999
            ELSE
               IF ( XX>XB ) THEN
                  XA = XB
                  YA = YB
               ELSE
                  XC = XB
                  YC = YB
               ENDIF
               XB = XX
               YB = YY
            ENDIF
         ENDIF
      ENDDO
 190  XB = XA + (PF-YA)*(XC-XA)/(YC-YA)
      YB = PFOFTK(WAVNA,WAVNB,XB)
      XX = XB
      IF ( ABS(YB-PF)<DELFIT ) THEN
         TKOFPF = XX
      ELSE
         PF = LOG(PF)
         YA = LOG(YA)
         YB = LOG(YB)
         YC = LOG(YC)
         LOGFIT = .TRUE.
         GOTO 120
      ENDIF
99999 END FUNCTION TKOFPF
 
      SUBROUTINE REPART(FXL,XLB,NXB,GYL,YLB,NYB)
      IMPLICIT NONE
 
!     ------------------------------------------------------------------
!     See procedure REPARTINT below
!     ------------------------------------------------------------------
      INTEGER, INTENT(IN)  ::  NXB, NYB
      REAL*8, INTENT(IN)  ::  FXL(NXB-1), XLB(NXB), YLB(NYB)
      REAL*8, INTENT(OUT)  ::  GYL(NYB-1)
 
      CALL REPARTINT(FXL,XLB,NXB,GYL,YLB,NYB,'repart')
      END SUBROUTINE REPART
 
      SUBROUTINE RETERP(FXL,XLB,NXB,GYL,YLB,NYB)
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     See procedure REPARTINT below
!     ------------------------------------------------------------------
      INTEGER, INTENT(IN)  ::  NXB, NYB
      REAL*8, INTENT(IN)  ::  FXL(NXB-1), XLB(NXB), YLB(NYB)
      REAL*8, INTENT(OUT)  ::  GYL(NYB-1)
 
      CALL REPARTINT(FXL,XLB,NXB,GYL,YLB,NYB,'interp')
 
      END SUBROUTINE RETERP
 
      SUBROUTINE REPARTINT(FXL,XLB,NXB,GYL,YLB,NYB,oper)
      IMPLICIT NONE
 
!     ------------------------------------------------------------------
!
!     REPART/RETERP
!              Repartitions or Interpolates FXL (a histogram-type
!              distribution function) where XLB depicts the NXB
!              partitions that define FXL data. FXL is assumed to be
!              constant between XLB(N) AND XLB(N+1)
!
!              GYL(N) is the new histogram distribution function defined
!              by the (input) NYB partitions YLB(N).  The YLB partitions
!              can differ from XLB in number and/or spacing, or in upper
!              and lower limits.  XLB and YLB coordinates are assumed to
!              be linear both increasing or decreasing in the same sense
!
!     RETURNS  GYL as a histogram-type distribution function, with NYB-1
!              values assumed to be constant between YLB(N) and YLB(N+1)
!
!       NOTE:  The column amount of a vertically distributed quantity is
!              conserved, within the Repartition Interval that is common
!              to to XLB and YLB (layer bottom edge and top edge) limits
!
!     ------------------------------------------------------------------
      INTEGER, INTENT(IN)  ::  nxb, nyb
      REAL*8, INTENT(IN)  ::  fxl(nxb-1), xlb(nxb), ylb(nyb)
      REAL*8, INTENT(OUT)  ::  gyl(nyb-1)
      CHARACTER(LEN=*), INTENT(IN)  ::  oper
                                           ! either 'repart' or 'interp'
      INTEGER NXF, NYG
      REAL*8 SUMG, SUMY, XA, YA, XB, YB, XAYA, PART
      INTEGER I, J
 
      NXF = NXB - 1
      NYG = NYB - 1
      SUMG = 0.D0
      DO I = 1, NYG
         GYL(I) = 0.D0
      ENDDO
      SUMY = 0.D0
      I = 1
      XA = XLB(I)
      J = 1
      YA = YLB(J)
      XB = XLB(I+1)
      IF ( XB<XA ) THEN
         DO
 
            YB = YLB(J+1)
            IF ( YB<XA ) THEN
               DO
                  XB = XLB(I+1)
                  IF ( XB<YA ) THEN
                     DO
                        XAYA = XA
                        IF ( YA<XA ) XAYA = YA
                        IF ( YB<XB ) THEN
                           PART = (XB-XAYA)/(XB-XA)
                           SUMG = SUMG + PART*FXL(I)
                           SUMY = SUMY + PART
                           I = I + 1
                           IF ( I>NXF ) THEN
                              SELECT CASE (oper)
                              CASE ('repart')
                                 GYL(J) = SUMG
                              CASE ('interp')
                                 GYL(J) = SUMG/SUMY
                              ENDSELECT
                              DO
                                 J = J + 1
                                 IF ( J>NYG ) GOTO 160
                                 GYL(J) = 0.D0
                              ENDDO
                           ELSE
                              XA = XB
                              XB = XLB(I+1)
                           ENDIF
                        ELSE
                           PART = (YB-XAYA)/(XB-XA)
                           SUMG = SUMG + PART*FXL(I)
                           SUMY = SUMY + PART
                           SELECT CASE (oper)
                           CASE ('repart')
                              GYL(J) = SUMG
                           CASE ('interp')
                              GYL(J) = SUMG/SUMY
                           ENDSELECT
                           J = J + 1
                           IF ( J>NYG ) GOTO 160
                           SUMG = 0.D0
                           SUMY = 0.D0
                           YA = YB
                           YB = YLB(J+1)
                        ENDIF
                     ENDDO
                  ELSE
                     I = I + 1
                     IF ( I>NXF ) GOTO 160
                     XA = XB
                  ENDIF
               ENDDO
            ELSE
               GYL(J) = 0.D0
               J = J + 1
               IF ( J>NYG ) EXIT
               YA = YB
            ENDIF
         ENDDO
      ELSE
         DO
            YB = YLB(J+1)
            IF ( YB>XA ) THEN
               DO
                  XB = XLB(I+1)
                  IF ( XB>YA ) THEN
                     DO
                        XAYA = XA
                        IF ( YA>XA ) XAYA = YA
                        IF ( YB>XB ) THEN
                           PART = (XB-XAYA)/(XB-XA)
                           SUMG = SUMG + PART*FXL(I)
                           SUMY = SUMY + PART
                           I = I + 1
                           IF ( I>NXF ) THEN
                              SELECT CASE (oper)
                              CASE ('repart')
                                 GYL(J) = SUMG
                              CASE ('interp')
                                 GYL(J) = SUMG/SUMY
                              ENDSELECT
                              DO
                                 J = J + 1
                                 IF ( J>NYG ) GOTO 160
                                 GYL(J) = 0.D0
                              ENDDO
                           ELSE
                              XA = XB
                              XB = XLB(I+1)
                           ENDIF
                        ELSE
                           PART = (YB-XAYA)/(XB-XA)
                           SUMG = SUMG + PART*FXL(I)
                           SUMY = SUMY + PART
                           SELECT CASE (oper)
                           CASE ('repart')
                              GYL(J) = SUMG
                           CASE ('interp')
                              GYL(J) = SUMG/SUMY
                           ENDSELECT
                           J = J + 1
                           IF ( J>NYG ) GOTO 160
                           SUMG = 0.D0
                           SUMY = 0.D0
                           YA = YB
                           YB = YLB(J+1)
                        ENDIF
                     ENDDO
                  ELSE
                     I = I + 1
                     IF ( I>NXF ) GOTO 160
                     XA = XB
                  ENDIF
               ENDDO
            ELSE
               GYL(J) = 0.D0
               J = J + 1
               IF ( J>NYG ) EXIT
               YA = YB
            ENDIF
         ENDDO
      ENDIF
 
 160  END SUBROUTINE REPARTINT
 
      SUBROUTINE FABINT(F,X,NX,ALIM,BLIM,ABINT)
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     FABINT  PERFORMS NUMERICAL INTEGRATION (AREA UNDER CURVE) OF F(X)
!             BETWEEN THE LIMITS X=ALIM AND X=BLIM  (WITH BLIM GT ALIM)
!
!       F(X)  IS DEFINED BY CONNECTING SUCCESSIVE F(X) DATA POINTS USING
!             STRAIGHT-LINE SEGMENTS, I.E. F(X) IS PIECE-WISE CONTINUOUS
!             THE  X  COORDINATE CAN BE IN ASCENDING OR DESCENDING ORDER
!
!             (F(X) IS ZERO OUTSIDE THE INTERVAL BETWEEN X(1) AND X(NX))
!     ------------------------------------------------------------------
      INTEGER, INTENT(IN)  ::  NX
      REAL*8, INTENT(IN)  ::  F(NX), X(NX), ALIM, BLIM
      REAL*8, INTENT(OUT)  ::  ABINT
      REAL*8, PARAMETER  ::  DELTA = 1.D-07
      REAL*8 XA, XB, XX, XMIN, XMAX, XJ, XI, FI, FJ, BF, AF, DINT, X2,  &
             X1
      INTEGER JX, KX, IX
 
      ABINT = 0.D0
      JX = 1
      KX = 1
      XA = X(JX)
      XB = X(NX)
      XX = XA
      IF ( XB<=XA ) THEN
         XA = XB
         XB = XX
         JX = NX
         KX = -1
      ENDIF
      XMIN = XA
      XMAX = XB
      IF ( XMIN>=BLIM ) RETURN
      IF ( XMAX<=ALIM ) RETURN
      IF ( XMIN<ALIM ) XMIN = ALIM
      IF ( XMAX>BLIM ) XMAX = BLIM
      DO
         JX = JX + KX
         XJ = X(JX)
         IF ( XJ>XMIN ) THEN
            IX = JX - KX
            XI = X(IX)
            IF ( (XJ-XI)>=DELTA ) THEN
               FI = F(IX)
               FJ = F(JX)
               BF = (FJ-FI)/(XJ-XI)
               AF = FJ - BF*XJ
               X2 = XMIN
               DO
                  X1 = X2
                  X2 = XJ
                  IF ( X2>XMAX ) X2 = XMAX
                  DINT = AF*(X2-X1) + BF*(X2**2-X1**2)/2.D0
                  ABINT = ABINT + DINT
                  IF ( DABS(X2-XMAX)<DELTA ) RETURN
                  IF ( (XJ-X2)<=DELTA ) THEN
                     DO
                        XI = XJ
                        FI = FJ
                        IX = JX
                        JX = JX + KX
                        XJ = X(JX)
                        FJ = F(JX)
                        IF ( DABS(XJ-XI)>=DELTA ) THEN
                           BF = (FJ-FI)/(XJ-XI)
                           AF = FJ - BF*XJ
                           EXIT
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDDO
      END SUBROUTINE FABINT
 
      SUBROUTINE FXGINT(F,X,NX,G,Y,NY,ALIM,BLIM,ABINT)
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     FXGINT  PERFORMS NUMERICAL INTEGRATION (AREA UNDER CURVE) OF  F*G
!             BETWEEN THE LIMITS X=ALIM AND X=BLIM  (WITH BLIM GT ALIM)
!
!       F(X)  IS DEFINED BY CONNECTING SUCCESSIVE F(X) DATA POINTS USING
!             STRAIGHT-LINE SEGMENTS, I.E. F(X) IS PIECE-WISE CONTINUOUS
!             THE  X  COORDINATE CAN BE IN ASCENDING OR DESCENDING ORDER
!
!       G(Y)  IS DEFINED BY CONNECTING SUCCESSIVE G(Y) DATA POINTS USING
!             STRAIGHT-LINE SEGMENTS, I.E. G(Y) IS PIECE-WISE CONTINUOUS
!             THE  Y  COORDINATE CAN BE IN ASCENDING OR DESCENDING ORDER
!
!             (X,Y ARE THE SAME LINEAR COORDINATE INDEPENDENTLY DEFINED)
!
!             (F(X) IS ZERO OUTSIDE THE INTERVAL BETWEEN X(1) AND X(NX))
!             (G(Y) IS ZERO OUTSIDE THE INTERVAL BETWEEN Y(1) AND Y(NY))
!     ------------------------------------------------------------------
      INTEGER, INTENT(IN)  ::  NX, NY
      REAL*8, INTENT(IN)  ::  F(NX), X(NX), G(NY), Y(NY), ALIM, BLIM
      REAL*8, INTENT(OUT)  ::  ABINT
      REAL*8, PARAMETER  ::  DELTA = 1.D-07
      REAL*8 XA, YA, XB, YB, XX, XMIN, XMAX, XJ, XI, FI, FJ, BF, AF, YI,&
             YJ, GI, GJ, AG, BG, DINT, X2, X1
      INTEGER JX, JY, KX, KY, IX, IY
 
      ABINT = 0.D0
      JX = 1
      JY = 1
      KX = 1
      KY = 1
      XA = X(JX)
      YA = Y(JY)
      XB = X(NX)
      YB = Y(NY)
      XX = XA
      IF ( XB<=XA ) THEN
         XA = XB
         XB = XX
         JX = NX
         KX = -1
      ENDIF
      XX = YA
      IF ( YB<=YA ) THEN
         YA = YB
         YB = XX
         JY = NY
         KY = -1
      ENDIF
      XMIN = MAX(XA,YA)
      XMAX = MIN(XB,YB)
      IF ( XMIN>=BLIM ) RETURN
      IF ( XMAX<=ALIM ) RETURN
      IF ( XMIN<ALIM ) XMIN = ALIM
      IF ( XMAX>BLIM ) XMAX = BLIM
      DO
         JX = JX + KX
         XJ = X(JX)
         IF ( XJ>XMIN ) THEN
            IX = JX - KX
            XI = X(IX)
            IF ( (XJ-XI)>=DELTA ) THEN
               FI = F(IX)
               FJ = F(JX)
               BF = (FJ-FI)/(XJ-XI)
               AF = FJ - BF*XJ
               DO
                  JY = JY + KY
                  YJ = Y(JY)
                  IF ( YJ>XMIN ) THEN
                     IY = JY - KY
                     YI = Y(IY)
                     IF ( (YJ-YI)>=DELTA ) THEN
                        GI = G(IY)
                        GJ = G(JY)
                        BG = (GJ-GI)/(YJ-YI)
                        AG = GJ - BG*YJ
                        X2 = XMIN
                        GOTO 160
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDDO
 160  X1 = X2
      X2 = MIN(XJ,YJ)
      IF ( X2>XMAX ) X2 = XMAX
      DINT = (AF*AG)*(X2-X1) + (AF*BG+BF*AG)*(X2**2-X1**2)              &
             /2.D0 + (BF*BG)*(X2**3-X1**3)/3.D0
      ABINT = ABINT + DINT
      IF ( DABS(X2-XMAX)<DELTA ) RETURN
      IF ( (XJ-X2)<=DELTA ) THEN
         DO
            XI = XJ
            FI = FJ
            IX = JX
            JX = JX + KX
            XJ = X(JX)
            FJ = F(JX)
            IF ( DABS(XJ-XI)>=DELTA ) THEN
               BF = (FJ-FI)/(XJ-XI)
               AF = FJ - BF*XJ
               EXIT
            ENDIF
         ENDDO
      ENDIF
      DO WHILE ( YJ<=X2 )
         YI = YJ
         GI = GJ
         IY = JY
         JY = JY + KY
         YJ = Y(JY)
         GJ = G(JY)
         IF ( DABS(YJ-YI)>=DELTA ) THEN
            BG = (GJ-GI)/(YJ-YI)
            AG = GJ - BG*YJ
            EXIT
         ENDIF
      ENDDO
      GOTO 160
      END SUBROUTINE FXGINT
 
      SUBROUTINE CTREND(JYEAR,IDEC,JDEC,CWTI,CWTJ)
      IMPLICIT NONE
 
!-------------------------------------------------------------------
!     Black Carbon interdecadal TAU interpolation is based on linear
!     TAU trend (between decadal global TAUmaps) with a superimposed
!     intra-decadal time dependence scaled to the Black Carbon Total
!     emission rate.
!
!        INPUT: JYEAR   (Julian year)
!
!                 CTREND coefficients refer to sep2003_OCI_Koch maps
!                 CTREND coefficients refer to sep2003_BCI_Koch maps
!                 --------------------------------------------------
!
!                 Map=  1850 1875 1900 1925 1950 1960 1970 1980 1990
!       OUTPUT:  IDEC=   (0)   1    2    3    4    5    6    7    8
!                JDEC=  IDEC + 1    (returned IDEC,JDEC are (1 to 8)
!
!                CWTI=   (Multiplicative Weight for BC DataMap IDEC)
!                CWTJ=   (Multiplicative Weight for BC DataMap JDEC)
!
!        NOTE:  Time dependence is linear before 1950. Industrial BC
!               is assumed 0 in 1850 so CWTI=0, and IDEC is set to 1
!-------------------------------------------------------------------
 
      INTEGER, INTENT(IN)  ::  JYEAR
      INTEGER, INTENT(OUT)  ::  IDEC, JDEC
      REAL*8, INTENT(OUT)  ::  CWTI, CWTJ
 
!    Global Annual Emissions of BC U   Emission (Mt/yr)
 
!      Year    Hard_Coal    Brown_Coal      Diesel        Total
      REAL*8, PARAMETER, DIMENSION(5,45)  ::  BCE = RESHAPE((/50.0,     &
                                   2.280581713,0.4449132979,            &
                                   0.1599090248,2.885536671,51.0,       &
                                   2.443193913,0.4855868816,            &
                                   0.1884280443,3.117194653,52.0,       &
                                   2.473641872,0.5115299225,            &
                                   0.2027695477,3.187930107,53.0,       &
                                   2.481340885,0.5448409319,            &
                                   0.2149295360,3.241089582,54.0,       &
                                   2.505670071,0.5780177116,            &
                                   0.2343477309,3.317960978,55.0,       &
                                   2.698692560,0.6238067150,            &
                                   0.2733324766,3.595800638,56.0,       &
                                   2.855226278,0.6531309485,            &
                                   0.3043369055,3.812692404,57.0,       &
                                   2.975781679,0.6821750998,            &
                                   0.3207367063,3.978575468,58.0,       &
                                   3.341105223,0.7035279870,            &
                                   0.3370627165,4.381746292,59.0,       &
                                   3.638528824,0.7075053453,            &
                                   0.3695519567,4.715488434,60.0,       &
                                   3.770926714,0.7416650057,            &
                                   0.3832504749,4.896034241,61.0,       &
                                   3.392980337,0.7805693150,            &
                                   0.4217525721,4.595387459,62.0,       &
                                   3.288835049,0.8179932237,            &
                                   0.4603823125,4.567360401,63.0,       &
                                   3.359177589,0.8604368567,            &
                                   0.5090782642,4.728550911,64.0,       &
                                   3.432664871,0.8952696323,            &
                                   0.5388473868,4.866865158,65.0,       &
                                   3.529418945,0.8819132447,            &
                                   0.5785927773,4.989773750,66.0,       &
                                   3.577459812,0.8817394972,            &
                                   0.6323299408,5.091631413,67.0,       &
                                   3.418204546,0.8635972142,            &
                                   0.6592246890,4.941041946,68.0,       &
                                   3.452457905,0.8943673372,            &
                                   0.7338049412,5.080585003,69.0,       &
                                   3.626069546,0.9298774004,            &
                                   0.7889106274,5.344810009,70.0,       &
                                   3.264039755,0.9229136109,            &
                                   0.8880128860,5.074741840,71.0,       &
                                   3.437611580,0.9374827743,            &
                                   0.9531223178,5.328329086,72.0,       &
                                   3.473345757,0.7836616039,            &
                                   1.0180075170,5.274850368,73.0,       &
                                   3.495583296,0.8056778908,            &
                                   1.1174367670,5.418928623,74.0,       &
                                   3.506143808,0.8251076341,            &
                                   1.0828053950,5.413989067,75.0,       &
                                   3.906814098,0.8527192473,            &
                                   1.0454736950,5.804963112,76.0,       &
                                   4.005736828,0.8900613785,            &
                                   1.1400985720,6.035901546,77.0,       &
                                   4.236912251,0.9103702307,            &
                                   1.2190728190,6.366260529,78.0,       &
                                   4.459666252,0.9303293228,            &
                                   1.2408012150,6.630728722,79.0,       &
                                   4.697422504,0.9856286645,            &
                                   1.3019220830,6.984815121,80.0,       &
                                   4.796229839,0.9959300756,            &
                                   1.2336660620,7.026207924,81.0,       &
                                   4.789204121,1.0459070210,            &
                                   1.1664049630,7.001126766,82.0,       &
                                   4.872739315,1.0975246430,            &
                                   1.1601715090,7.130136490,83.0,       &
                                   4.983223438,1.1424025300,            &
                                   1.1732926370,7.298912525,84.0,       &
                                   5.265352249,1.2178678510,            &
                                   1.2251536850,7.708741188,85.0,       &
                                   5.763637543,1.2965050940,            &
                                   1.2428865430,8.303324699,86.0,       &
                                   5.924767494,1.3386499880,            &
                                   1.2930148840,8.556744576,87.0,       &
                                   6.155550480,1.3738890890,            &
                                   1.3162037130,8.845513344,88.0,       &
                                   6.379704475,1.3670797350,            &
                                   1.3813229800,9.127896309,89.0,       &
                                   6.594299316,1.4169263840,            &
                                   1.4029121400,9.414231300,90.0,       &
                                   6.566919804,1.4685817960,            &
                                   1.4224120380,9.458042145,91.0,       &
                                   6.661097050,1.2067918780,            &
                                   1.4163945910,9.284657478,92.0,       &
                                   7.737902641,1.3509917260,            &
                                   1.4471185210,10.53625107,93.0,       &
                                   7.393332005,1.2448183300,            &
                                   1.4543261530,10.09271908,94.0,       &
                                   7.515841007,1.2333894970,            &
                                   1.4780857560,10.22745800/),(/5,45/))
 
      REAL*8 XDEC
      INTEGER IBCDEC, JBCDEC, IJYEAR
 
      IF ( JYEAR<1876 ) THEN
         CWTJ = (JYEAR-1850)/25.D0
         IF ( CWTJ<0.D0 ) CWTJ = 0.D0
         CWTI = 0.D0
         IDEC = 1
         JDEC = 1
         GOTO 100
      ENDIF
 
      IF ( JYEAR<1950 ) THEN
         XDEC = (JYEAR-1850)/25.D0
         IDEC = XDEC
         JDEC = IDEC + 1
         CWTJ = XDEC - IDEC
         CWTI = 1.D0 - CWTJ
         GOTO 100
      ENDIF
 
      IF ( JYEAR<1990 ) THEN
         IDEC = (JYEAR-1910)/10
         JDEC = IDEC + 1
         IBCDEC = 1 + (IDEC-4)*10
         JBCDEC = IBCDEC + 10
         IJYEAR = JYEAR - 1949
         CWTJ = (BCE(5,IJYEAR)-BCE(5,IBCDEC))                           &
                /(BCE(5,JBCDEC)-BCE(5,IBCDEC))
         CWTI = 1.D0 - CWTJ
         GOTO 100
      ENDIF
 
      IF ( JYEAR>1989 ) THEN
         IDEC = 7
         JDEC = 8
         IJYEAR = JYEAR - 1949
         IF ( IJYEAR>45 ) IJYEAR = 45
         CWTJ = BCE(5,IJYEAR)/BCE(5,41)
         CWTI = 0.D0
      ENDIF
 
 100  END SUBROUTINE CTREND
 
      SUBROUTINE STREND(JYEAR,IDEC,JDEC,SWTI,SWTJ)
      IMPLICIT NONE
 
!-------------------------------------------------------------------
!     Anthropogenic Sulfate inter-decadal TAU interpolation is based
!     on a linear TAU trend (between decadal global TAU-maps) with a
!     superimposed intradecadal time dependence scaled in proportion
!     to the Anthropogenic Sulfate global emission rate.
!
!        INPUT: JYEAR   (Julian year)
!
!                 CTREND coefficients refer to sep2003_SUI_Koch maps
!                 --------------------------------------------------
!
!                 Map=  1850 1875 1900 1925 1950 1960 1970 1980 1990
!       OUTPUT:  IDEC=   (0)   1    2    3    4    5    6    7    8
!                JDEC=  IDEC + 1    (returned IDEC,JDEC are (1 to 8)
!
!                SWTI=  (Multiplicative Weight for SUI DataMap IDEC)
!                SWTJ=  (Multiplicative Weight for SUI DataMap JDEC)
!
!        NOTE:  Time dependence linear before 1950.   Industrial SUI
!               is assumed 0 in 1850 so SWTI=0, and IDEC is set to 1
!-------------------------------------------------------------------
 
      INTEGER, INTENT(IN)  ::  JYEAR
      INTEGER, INTENT(OUT)  ::  IDEC, JDEC
      REAL*8, INTENT(OUT)  ::  SWTI, SWTJ
 
!     Global Emission of Sulfate
 
!     Emission (Mt/yr)
!               year      Anthropogenic_Sulfate Natural_Sulfate
      REAL*8, PARAMETER, DIMENSION(3,41)  ::  SUE = RESHAPE((/1950.0,   &
                                   30.46669769,14.4,1951.0,32.38347244, &
                                   14.4,1952.0,32.18632889,14.4,1953.0, &
                                   32.83379745,14.4,1954.0,32.79270935, &
                                   14.4,1955.0,35.79611969,14.4,1956.0, &
                                   39.93603897,14.4,1957.0,38.68806839, &
                                   14.4,1958.0,39.35904312,14.4,1959.0, &
                                   41.06065369,14.4,1960.0,42.67050934, &
                                   14.4,1961.0,41.32410431,14.4,1962.0, &
                                   41.80470276,14.4,1963.0,43.26312637, &
                                   14.4,1964.0,44.68368530,14.4,1965.0, &
                                   45.81701660,14.4,1966.0,46.61584091, &
                                   14.4,1967.0,46.42276001,14.4,1968.0, &
                                   47.77438354,14.4,1969.0,49.30817032, &
                                   14.4,1970.0,52.81050873,14.4,1971.0, &
                                   52.95043945,14.4,1972.0,54.10167694, &
                                   14.4,1973.0,55.93037415,14.4,1974.0, &
                                   57.31056213,14.4,1975.0,58.52788162, &
                                   14.4,1976.0,59.71361542,14.4,1977.0, &
                                   62.59599304,14.4,1978.0,61.98198318, &
                                   14.4,1979.0,64.71042633,14.4,1980.0, &
                                   65.28986359,14.4,1981.0,63.23768234, &
                                   14.4,1982.0,62.88000488,14.4,1983.0, &
                                   61.45023346,14.4,1984.0,63.85008621, &
                                   14.4,1985.0,66.47412872,14.4,1986.0, &
                                   68.00902557,14.4,1987.0,69.87956238, &
                                   14.4,1988.0,70.52937317,14.4,1989.0, &
                                   72.06355286,14.4,1990.0,71.29174805, &
                                   14.4/),(/3,41/))
 
      REAL*8 xdec
      INTEGER ISUDEC, JSUDEC, IJYEAR
 
      IF ( JYEAR<1876 ) THEN
         SWTJ = (JYEAR-1850)/25.D0
         IF ( SWTJ<0.D0 ) SWTJ = 0.D0
         SWTI = 0.D0
         IDEC = 1
         JDEC = 1
         GOTO 100
      ENDIF
 
      IF ( JYEAR<1950 ) THEN
         XDEC = (JYEAR-1850)/25.D0
         IDEC = XDEC
         JDEC = IDEC + 1
         SWTJ = XDEC - IDEC
         SWTI = 1.D0 - SWTJ
         GOTO 100
      ENDIF
 
      IF ( JYEAR<1990 ) THEN
         IDEC = (JYEAR-1910)/10
         JDEC = IDEC + 1
         ISUDEC = 1 + (IDEC-4)*10
         JSUDEC = ISUDEC + 10
         IJYEAR = JYEAR - 1949
         SWTJ = (SUE(2,IJYEAR)-SUE(2,ISUDEC))                           &
                /(SUE(2,JSUDEC)-SUE(2,ISUDEC))
         SWTI = 1.D0 - SWTJ
         GOTO 100
      ENDIF
 
      IF ( JYEAR>1989 ) THEN
         IDEC = 7
         JDEC = 8
         IJYEAR = JYEAR - 1949
         IF ( IJYEAR>41 ) IJYEAR = 41
         SWTJ = SUE(2,IJYEAR)/SUE(2,41)
         SWTI = 0.D0
      ENDIF
 
 100  END SUBROUTINE STREND
 
      SUBROUTINE SPLINV(X,F,NXF,XX,FF,CUSPWM,CUSPWE,KXTRAP)
      IMPLICIT NONE
 
      INTEGER, INTENT(IN)  ::  NXF, KXTRAP
      REAL*8, INTENT(IN)  ::  X(NXF), F(NXF), FF, CUSPWM, CUSPWE
      REAL*8, INTENT(OUT)  ::  XX
 
!---------------------------------------------------------------------
!    Inverse spline:
!    SPLINV locates FF between points (F2,X2)(F3,X3) on 4-point spread
!    and returns 4-point Cubic Spline value of XX such that FF = F(XX)
!
!    Quadratic Derivatives of Spline are continuous at (F2,X2),(F3,X3)
!    (X-Coordinate may be specified in increasing or decreasing order)
!
!---------------------------------------------------------------------
!
!    Custom Control Parameters:  CUSPWM,CUSPWE
!------------------------------
!
!    In cases where data points are unevenly spaced and/or data points
!    exhibit abrupt changes in value, Spline Interpolation may produce
!    undesirable bulging of interpolated values. In more extreme cases
!    Linear Interpolation may be less problematic to use.
!
!    Interpolation can be weighted between: Cubic Spline and Linear by
!    adjusting weights CUSPWM and CUSPWE to values between 1.0 and 0.0
!
!    CUSPWM = Cubic Spline Weight at the (X2-X3) Interval Mid-point
!    CUSPWE = Cubic Spline Weight at the (X2-X3) Interval End-points
!
!    For example, with:
!
!    CUSPWM=1.0,CUSPWE=1.0  FF returns Cubic Spline interpolated value
!    CUSPWM=0.0,CUSPWE=0.0  FF returns   Linearly   interpolated value
!
!---------------------------------------------------------------------
!
!     Extrapolation for XX outside of defined interval:  X(1)<->X(NXF)
!
!               KXTRAP = 0    No Extrapolation   (i.e., sets XX = 0.0)
!                        1    Fixed Extrapolation (sets XX=edge value)
!                        2    Linear Extrapolation using 2 edge points
!
!---------------------------------------------------------------------
!
!
!    NOTE:  F(X) is assumed to be monotonic between F(1) and F(NXF)
!
!------------------------------------------------------------------
 
      REAL*8 x1, x2, x3, x4, x21, x32, x43, x31, x42, BETW, FFCUSP,     &
             FFLINR, CUSPWT
      REAL*8 f1, f2, f3, f4, f21, f32, f43, f3221, f4332, a, b, c, d,   &
             xf, xe, xexm
      REAL*8 DX, gg, xg, xy, deltx, slopec, slopel, slopes
      INTEGER k, kk
 
      BETW = (F(2)-FF)*(F(NXF)-F(1))
      IF ( BETW>0.D0 ) THEN
 
!                Edge Point Interval Interpolation and/or Extrapolation
!                ------------------------------------------------------
         BETW = (F(1)-FF)*(F(NXF)-F(1))
         IF ( BETW>0.D0 ) THEN
 
!                  Extrapolation for FF Outside of Interval F(1) - F(2)
!                  ----------------------------------------------------
!                  IF(KXTRAP == 0)  (No Extrapolation:   sets XX = 0.0)
!                  IF(KXTRAP == 1)  (Extrapolation at Fixed Edge Value)
!                  IF(KXTRAP == 2)  (2 Edge Point Linear Extrapolation)
 
            IF ( KXTRAP==0 ) XX = 0.D0
            IF ( KXTRAP==1 ) XX = X(1)
            IF ( KXTRAP==2 ) XX = X(1) - (F(1)-FF)/(F(2)-F(1))          &
                                  *(X(2)-X(1))
         ELSE
 
!                          F(1),F(2)  Edge Point Interval Interpolation
!                          --------------------------------------------
            DO KK = 2, 6
               X1 = X(1)
               X2 = X(2)
               X3 = X(3)
               F1 = F(1)
               F2 = F(2)
               F3 = F(3)
               XX = X1 + (FF-F(1))/(F(2)-F(1))*(X2-X1)
               XF = XX - X1
               X21 = X2 - X1
               F21 = (F2-F1)/X21
               X32 = X3 - X2
               X31 = X3 - X1
               C = ((F3-F2)/X32-F21)/X31
               B = F21 - X21*C
               A = F1
               FFCUSP = A + XF*(B+XF*C)
               FFLINR = A + XF*F21
               XE = 1.D0 - 2.D0*XF/X21
               XEXM = XE**2
               CUSPWT = (1.D0-XEXM)*CUSPWM + XEXM*CUSPWE
               GG = FFCUSP*CUSPWT + FFLINR*(1.D0-CUSPWT)
               SLOPEC = B + 2.D0*C*XF
               SLOPEL = F21
               SLOPES = SLOPEC*CUSPWT + SLOPEL*(1.D0-CUSPWT)
               XG = XF
               DELTX = (GG-FF)/SLOPES
               XX = XF - (GG-FF)/SLOPES + X1
            ENDDO
         ENDIF
      ELSE
         BETW = (FF-F(NXF-1))*(F(NXF)-F(1))
         IF ( BETW>0.D0 ) THEN
 
            BETW = (FF-F(NXF))*(F(NXF)-F(1))
            IF ( BETW>0.D0 ) THEN
 
!              Extrapolation for F Outside of Interval  F(NXF-1)-F(NXF)
!              --------------------------------------------------------
!                  IF(KXTRAP == 0)  (No Extrapolation:   sets XX = 0.0)
!                  IF(KXTRAP == 1)  (Extrapolation at Fixed Edge Value)
!                  IF(KXTRAP == 2)  (2 Edge Point Linear Extrapolation)
 
               IF ( KXTRAP==0 ) XX = 0.D0
               IF ( KXTRAP==1 ) XX = X(NXF)
               IF ( KXTRAP==2 ) XX = X(NXF) - (F(NXF)-FF)               &
                                     /(F(NXF-1)-F(NXF))                 &
                                     *(X(NXF-1)-X(NXF))
            ELSE
 
!                    F(NXF-1),F(NXF)  Edge Point Interval Interpolation
!                    --------------------------------------------------
               DO KK = 3, 7
                  X1 = X(NXF-2)
                  X2 = X(NXF-1)
                  X3 = X(NXF)
                  F1 = F(NXF-2)
                  F2 = F(NXF-1)
                  F3 = F(NXF)
                  XX = X2 + (FF-F2)/(F3-F2)*(X3-X2)
                  XF = XX - X2
                  X32 = X3 - X2
                  F32 = (F3-F2)/X32
                  X21 = X2 - X1
                  X31 = X3 - X1
                  F21 = (F2-F1)/X21
 
!                    3-Point Quadratic Interpolation for Edge Intervals
!                    --------------------------------------------------
!
!      (Edge Option)     ----------------------------------------------
!                        For Linear Interpolation within Edge Intervals
!                        between F(1),F(2), and between F(NXF-1),F(NXF)
!                        set the value of coefficient C below, to C=0.0
!                        ----------------------------------------------
 
                  C = (F32-F21)/X31
                  B = F21 + X21*C
                  A = F2
                  FFCUSP = A + XF*(B+XF*C)
                  FFLINR = A + XF*F32
                  XE = 1.D0 - 2.D0*XF/X32
                  IF ( XE<0.D0 ) XE = -XE
                  XEXM = XE**2
                  CUSPWT = (1.D0-XEXM)*CUSPWM + XEXM*CUSPWE
                  GG = FFCUSP*CUSPWT + FFLINR*(1.D0-CUSPWT)
                  SLOPEC = B + 2.D0*C*XF
                  SLOPEL = F21
                  SLOPES = SLOPEC*CUSPWT + SLOPEL*(1.D0-CUSPWT)
                  XG = XF
                  DELTX = (GG-FF)/SLOPES
                  XX = XF - (GG-FF)/SLOPES + X2
               ENDDO
            ENDIF
         ELSE
 
            DO K = 3, NXF - 1
               BETW = (FF-F(K-1))*(F(K)-FF)
               DX = (FF-F(K-1))/(F(K)-F(K-1))
               XX = X(K-1) + DX*(X(K)-X(K-1))
               IF ( BETW>=0.D0 ) EXIT
            ENDDO
 
            DO KK = 1, 5
               X1 = X(K-2)
               X2 = X(K-1)
               X3 = X(K)
               X4 = X(K+1)
               F1 = F(K-2)
               F2 = F(K-1)
               F3 = F(K)
               F4 = F(K+1)
               X21 = X2 - X1
               X31 = X3 - X1
               X32 = X3 - X2
               X43 = X4 - X3
               X42 = X4 - X2
               F21 = (F2-F1)/(X21*X21)
               F32 = (F3-F2)/(X32*X32)
               F43 = (F4-F3)/(X43*X43)
               F3221 = (F32+F21)/X31*X21
               F4332 = (F43+F32)/X42*X43
               A = F2
               B = X32*F3221
               C = 3.D0*F32 - F3221 - F3221 - F4332
               D = (F3221+F4332-F32-F32)/X32
               XF = XX - X2
 
!                             FFCUSP= Cubic Spline Interpolation Result
!                             -----------------------------------------
 
               FFCUSP = A + XF*(B+XF*(C+XF*D))
               XE = (X3+X2-XX-XX)/X32
               IF ( XE<0.D0 ) XE = -XE
               XEXM = XE**2
               CUSPWT = (1.D0-XEXM)*CUSPWM + XEXM*CUSPWE
 
!                                   FFLINR= Linear Interpolation Result
!                                   -----------------------------------
               FFLINR = A + XF*F32*X32
               GG = FFCUSP*CUSPWT + FFLINR*(1.D0-CUSPWT)
               SLOPEC = B + 2.D0*C*XF + 3.D0*D*XF**2
               SLOPEL = F32*X32
               SLOPES = SLOPEC*CUSPWT + SLOPEL*(1.D0-CUSPWT)
               XG = XF
               XY = XX
               DELTX = (GG-FF)/SLOPES
               XX = XF - (GG-FF)/SLOPES + X2
            ENDDO
         ENDIF
      ENDIF
 
      END SUBROUTINE SPLINV
 
      SUBROUTINE THREEPTQUADINTERPOLATION(NVEC,X21,X31,X32,XF,CUSPWM,   &
         CUSPWE,F21,F32,F2,FF)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  ::  NVEC
      REAL*8, INTENT(IN)  ::  X21, X31, X32, XF, CUSPWM, CUSPWE
      REAL*8, DIMENSION(NVEC), INTENT(IN)  ::  F32, F21, F2
      REAL*8, INTENT(OUT)  ::  FF(NVEC)
 
!                    3-Point Quadratic Interpolation for Edge Intervals
!                    --------------------------------------------------
!
!      (Edge Option)     ----------------------------------------------
!                        For Linear Interpolation within Edge Intervals
!                        between X(1),X(2), and between X(NXF-1),X(NXF)
!                        set the value of coefficient C below, to C=0.0
!                        ----------------------------------------------
      INTEGER  ::  K
      REAL*8  ::  A, B, C, FFCUSP, XE, XEXM, CUSPWT, FFLINR
 
      DO K = 1, NVEC
         C = (F32(K)-F21(K))/X31
         B = F21(K) + X21*C
         A = F2(K)
         FFCUSP = A + XF*(B+XF*C)
         FFLINR = A + XF*F32(K)
         XE = 1.D0 - 2.D0*XF/X32
         IF ( XE<0.D0 ) XE = -XE
         XEXM = XE**2
         CUSPWT = (1.D0-XEXM)*CUSPWM + XEXM*CUSPWE
         FF(K) = FFCUSP*CUSPWT + FFLINR*(1.D0-CUSPWT)
      ENDDO
      END SUBROUTINE THREEPTQUADINTERPOLATION
 
      SUBROUTINE SPLN44(Q,NI,NJ,IR,DR,JN,DN,QQ)
      IMPLICIT NONE
 
      INTEGER, INTENT(IN)  ::  NI, NJ, IR, JN
      REAL*8, INTENT(IN)  ::  Q(NI,NJ), DR, DN
      REAL*8, INTENT(OUT)  ::  QQ
 
!nu   REAL*8,save :: CUSPWM=1., CUSPWE=1. ,CUSPWT,fflinr
      REAL*8 QK(4)
      REAL*8 f1, f2, f3, f4
      INTEGER k, kr, irm, irp
      REAL*8, EXTERNAL  ::  COMPUTE2
 
      K = 0
      IRM = IR - 1
      IRP = IR + 2
      DO KR = IRM, IRP
         K = K + 1
         F1 = Q(KR,JN-1)
         F2 = Q(KR,JN)
         F3 = Q(KR,JN+1)
         F4 = Q(KR,JN+2)
         QK(K) = COMPUTE2(F1,F2,F3,F4,DN)
      ENDDO
      F1 = QK(1)
      F2 = QK(2)
      F3 = QK(3)
      F4 = QK(4)
      QQ = COMPUTE2(F1,F2,F3,F4,DR)
      END SUBROUTINE SPLN44
 
      SUBROUTINE SPLNI4(Q,NI,NJ,IR,JN,DN,QQ)
      IMPLICIT NONE
 
      INTEGER, INTENT(IN)  ::  NI, NJ, IR, JN
      REAL*8, INTENT(IN)  ::  Q(NI,NJ), DN
      REAL*8, INTENT(OUT)  ::  QQ
 
!nu   REAL*8,save :: CUSPWM=1., CUSPWE=1. ,CUSPWT,fflinr
      REAL*8 f1, f2, f3, f4
      REAL*8, EXTERNAL  ::  COMPUTE2
 
      F1 = Q(IR,JN-1)
      F2 = Q(IR,JN)
      F3 = Q(IR,JN+1)
      F4 = Q(IR,JN+2)
      QQ = COMPUTE2(F1,F2,F3,F4,DN)
      END SUBROUTINE SPLNI4
  
                                                    !  input parameters
                                                    !  SW input ( 6,110)
                                                    !  LW input (33,110)
                                                    !     SETREL  output
                                                    !  RQ input    (110)
                                                    !  3(6,190),(33,190)
      SUBROUTINE SETREL(REFF0,NAER,KDREAD,SRUQEX,SRUQSC,SRUQCB,TRUQEX,  &
                        TRUQSC,TRUQCB,REFU22,Q55U22,FRSULF,SRHQEX,      &
                        SRHQSC,SRHQCB,TRHQAB,RHDATA)
                                                    !  RH info   (190,9)
 
      USE FILEMANAGER, ONLY:OPENUNIT, CLOSEUNIT
      USE DOMAIN_DECOMP_ATM, ONLY:AM_I_ROOT
      IMPLICIT NONE
 
      INTEGER NAER, KDREAD
      REAL*8 REFF0, SRUQEX(6,110), SRUQSC(6,110), SRUQCB(6,110),        &
             TRUQEX(33,110), TRUQSC(33,110), TRUQCB(33,110), REFU22(110)&
             , Q55U22(110), FRSULF(8)
      REAL*8 SRHQEX(6,190), SRHQSC(6,190), SRHQCB(6,190), TRHQAB(33,190)&
             , RHDATA(190,15)
 
!     ------------------------------------------------------------------
!     REFF0  = Effective radius for dry aerosol seed size (in microns)
!     NAER   = Aerosol composition index
!     KDREAD = IO READ unit number for Q(m,r),g(m,r) data used by SETREL
!     ------------------------------------------------------------------
!     Aerosol index = NAER    Composition         Input data order = NNA
!                      1      SO4  Sulfate                            1
!                      2      SEA  Sea Salt                           2
!                      3      NO3  Nitrate                            3
!                                  Pure Water                         4
!                      4      ORG  Organic                            5
!     ------------------------------------------------------------------
 
      CHARACTER*40, SAVE  ::  dtfile =                                  &
                             &'oct2003.relhum.nr.Q633G633.table'
      LOGICAL qexist
      INTEGER i, j, j1, k, k1, in1, ir1, jdry, jwet, jhimax, khimax,    &
              maxdry, maxwet
      INTEGER n, n0, n1, nn, np, nrhn1
      REAL*8 x, xx, xi, xn0, xn1, xr1, ff, fi, gi, gd1, gd2, gw1, gw2,  &
             grh, qrh, rrh
      REAL*8 rh, rhi, rr0, rd1, rd2, rw1, rw2, dwr, qd1, qd2, qw1, qw2, &
             xdry, sdry
      REAL*8 xwet, swet, qqdmax, qqwmax, rqdmax, rqwmax, q55dry, q63dry,&
             dwn
      REAL*8 aermas, ddry, dwet, reffi, rhrhi, sum, sumw, vd1, vd2, vw1,&
             vw2
      REAL*8 w1, w2, w3, w4, wd1, wd2, ww1, ww2, wtx, wty, wtz, wts,    &
             wta, xfdry
      REAL*8 q55rh1, q55rh2, q55rh3, q55rh4, q550, q633, qgaerx, qscqcb
 
!     Output variables (RHDATA/RHINFO)
 
      REAL*8 RHRHRH(190), RHTAUF(190), RHREFF(190), RHWGM2(190),        &
             RHDGM2(190), RHTGM2(190), RHXMFX(190), RHDENS(190),        &
             RHQ550(190), TAUM2G(190), XNRRHX(190), ANBCM2(190),        &
             COSBAR(190), PIZERO(190), ANGSTR(190), RHINFO(190,15)
 
      EQUIVALENCE (RHINFO(1,1),RHRHRH(1))
      EQUIVALENCE (RHINFO(1,2),RHTAUF(1))
      EQUIVALENCE (RHINFO(1,3),RHREFF(1))
      EQUIVALENCE (RHINFO(1,4),RHWGM2(1))
      EQUIVALENCE (RHINFO(1,5),RHDGM2(1))
      EQUIVALENCE (RHINFO(1,6),RHTGM2(1))
      EQUIVALENCE (RHINFO(1,7),RHXMFX(1))
      EQUIVALENCE (RHINFO(1,8),RHDENS(1))
      EQUIVALENCE (RHINFO(1,9),RHQ550(1))
      EQUIVALENCE (RHINFO(1,10),TAUM2G(1))
      EQUIVALENCE (RHINFO(1,11),XNRRHX(1))
      EQUIVALENCE (RHINFO(1,12),ANBCM2(1))
      EQUIVALENCE (RHINFO(1,13),COSBAR(1))
      EQUIVALENCE (RHINFO(1,14),PIZERO(1))
      EQUIVALENCE (RHINFO(1,15),ANGSTR(1))
 
!     ------------------------------------------------------------------
!     RHDATA/    Local
!     RHINFO   Variable                   Description
!     ------   --------   ----------------------------------------------
!        1      RHRHRH    Relative humidity index RH (DO 110  0.0-0.999)
!        2      RHTAUF    Dry TAU multiplication factor due to RH effect
!        3      RHREFF    RH dependent effective radius
!        4      RHWGM2    Liquid water content (g/m2) per unit (dry) TAU
!        5      RHDGM2    Dry mass density     (g/m2) per unit (dry) TAU
!        6      RHTGM2    Total mass density   (g/m2) per unit (dry) TAU
!        7      RHXMFX    Dry mass fraction X of total aerosol mass
!        8      RHDENS    RH dependent density (g/cm3)
!        9      RHQ550    RH dependent Mie extinction efficiency (550nm)
!       10      TAUM2G    RH dependent TAU factor  (m2/g) of dry aerosol
!       11      XNRRHX    RH dependent real refractive index
!       12      ANBCM2    Aerosol Number density (Billion)/cm2
!       13      COSBAR    RH dependent Mie asymmetry parameter (visible)
!       14      PIZERO    RH dependent single scattering albedo(visible)
!       15      ANGSTR    Angstrom exponent = -(1-SRHQEX(5)/(0.55-0.815)
!     ------------------------------------------------------------------
 
 
!     Local variables
 
      REAL*8 R633NR(890), XNR(31), Q633NR(890,31), G633NR(890,31)
      REAL*8 Q880M1(890), G880M1(890), Q880M0(890), G880M0(890)
      REAL*8 Q880N1(890), Q880N0(890), R550NR(890), SMOOTH(890)
      REAL*8 RR0RHX(190), QRH633(190), GRH633(190), DNRX(190)
 
 
      REAL*8 QXAERN(33), QSAERN(33), QGAERN(33), SR1QEX(6), SR1QSC(6),  &
             SR1QCB(6), SR2QEX(6), SR2QSC(6), SR2QCB(6), SR3QEX(6),     &
             SR3QSC(6), SR3QCB(6), SR4QEX(6), SR4QSC(6), SR4QCB(6),     &
             TR1QEX(33), TR1QSC(33), TR1QCB(33), TR2QEX(33), TR2QSC(33),&
             TR2QCB(33), TR3QEX(33), TR3QSC(33), TR3QCB(33), TR4QEX(33),&
             TR4QSC(33), TR4QCB(33), TRHQEX(33), TRHQSC(33), TRHQCB(33)
 
      INTEGER, PARAMETER, DIMENSION(4)  ::  NRHCRY = (/38,47,28,38/)
 
      CHARACTER*8 AERTYP(4)
      DATA AERTYP/'Sulfate ', 'SeaSalt ', 'Nitrate ', 'Organic '/
 
!     ------------------------------------------------------------------
!     Hygroscopic aerosols (Sulfate,SeaSalt,Nitrate) physical properties
!     formulas from Tang and Munkelwitz (1994, 1996) in JGR 99, JGR 101.
!
!     AW=water activity RO=density  BX=growth factor RX=refractive index
!     SO4 = ammonium sulfate;   SEA = sea salt;   NO3 = ammonium nitrate
!     ------------------------------------------------------------------
 
!     functions
 
      REAL*8 AWSO4, DWSO4, ROSO4, BXSO4, RXSO4, DRWSO4, DRDSO4
      REAL*8 AWSEA, DWSEA, ROSEA, BXSEA, RXSEA, DRWSEA, DRDSEA
      REAL*8 RRSEA, VVSEA, GXSEA
      REAL*8 AWNO3, DWNO3, RONO3, BXNO3, R1NO3, R2NO3, DRXNO3
      REAL*8 AWOCX, DWOCX, ROOCX, BXOCX, RXOCX, DRWOCX, DRDOCX
 
      !        Sulfate parametric formulas from Tang   Munkelwitz(94,96)
      AWSO4(X) = 1.D0 - 0.2715*X + 0.3113*X**2 - 2.336*X**3 + 1.412*X**4
                                                                  ! TM94
      DWSO4(X) = -0.2715D0 + 0.6226*X - 7.008*X**2 + 5.648*X**3
      ROSO4(X) = 0.9971D0 + 5.92D-01*X - 5.036D-02*X**2 + 1.024D-02*X**3
                                                                  ! TM94
      BXSO4(X) = (1.D0/X*1.760D0/ROSO4(X))**(1.D0/3.D0)           ! TM96
      RXSO4(X) = 1.3330 + 0.16730*X - 0.0395*X**2                 ! TM91
      DRWSO4(RH) = 1.002146 - 0.00149*RH + 0.001*RH/(1.0+0.911*RH**10)
      DRDSO4(RH) = 1.002503   ! ratio of wet   dry nr(0.550) / nr(0.633)
 
      !        SeaSalt parametric formulas from Tang   Munkelwitz(94,96)
      AWSEA(X) = 1.0D0 - 0.6366*X + 0.8624*X**2 - 11.58*X**3 +          &
                 15.18*X**4                                       ! TM96
      DWSEA(X) = -0.6366D0 + 1.7248*X - 34.74*X**2 + 60.72*X**3
      ROSEA(X) = 0.9971 + 0.741*X - 0.3741*X**2 + 2.252*X**3 -          &
                 2.060*X**4                                       ! TM96
      BXSEA(X) = (1.D0/X*2.165D0/ROSEA(X))**(1.D0/3.D0)
      RRSEA(X) = 3.70958 + (8.95-3.70958)/(1.D0+(1.0-X)/X*58.448/18.0)
      VVSEA(X) = (18.0+(58.448-18.0)/(1.0+(1.0-X)/X*58.448/18.0))       &
                 /ROSEA(X)
      GXSEA(X) = SQRT((2.D0*RRSEA(X)+VVSEA(X))/(VVSEA(X)-RRSEA(X)))
                                                                  ! TM96
      RXSEA(X) = 1.333 + (GXSEA(X)-1.333)*(1.490-1.333)/(1.544-1.333)
      DRWSEA(RH) = 1.00212 - 0.001625*RH + 0.00131*RH/(1.0+0.928*RH**3)
      DRDSEA(RH) = 1.003007   ! ratio of wet   dry nr(0.550) / nr(0.633)
 
      !        Nitrate parametric formulas from Tang   Munkelwitz(94,96)
      AWNO3(X) = 1.D0 - 3.65D-01*X - 9.155D-02*X**2 - 2.826D-01*X**3
                                                                  ! TM96
      DWNO3(X) = -3.65D-01 - 18.31D-02*X - 8.478D-01*X**3
      RONO3(X) = 0.9971D0 + 4.05D-01*X + 9.0D-02*X**2             ! TM96
      BXNO3(X) = (1.D0/X*1.725D0/RONO3(X))**(1.D0/3.D0)           ! TM96
      R1NO3(X) = 1.3330 + 0.119D0*X      !  (X<0.205)              TWM81
      R2NO3(X) = 1.3285 + 0.145D0*X      !  (X>0.205)              TWM81
      DRXNO3(RH) = 1.001179   ! ratio of wet   dry nr(0.550) / nr(0.633)
 
      !        Organic Carbon - adapted from Sulfate parametric formulas
      !        yields growth factor G=1.1 at RH=0.84 Virkkula et al 1999
      AWOCX(X) = 1D0 - X**8D0
      DWOCX(X) = -8D0*X**7D0
      ROOCX(X) = 1D0 + .5D0*X
      BXOCX(X) = (1.5D0/(X*ROOCX(X)))**(1D0/3D0)
      RXOCX(X) = 1.3330D0 + .193D0*X
      DRWOCX(RH) = 1.00253 - 0.00198*RH + 0.00184*RH/(1.0+0.656*RH**1.1)
      DRDOCX(RH) = 1.00253
 
!     ------------------------------------------------------------------
!     Q,G Mie data (879x31) at 0.633 microns, use 31 points to cover the
!     refractive index from 1.30 to 1.60 with equal spacing of 0.01
!
!     Q,G data effective radius spans the range from 0.0 to 20.4 microns
!     in (3) segments of equally spaced data for optimized 4-point Cubic
!     Spline interpolation.  The equally spaced segments are as follows:
!
!     Index:    1 - 303   304 - 603   604 -  879   881 - 885   886 - 890
!     Reff:   0.00-3.02   3.04-9.02   9.04-20.04   2.98-3.04   8.96-9.08
!     Delta:     0.01        0.02        0.04         0.02        0.04
!
!     The last two intervals are constructed to accommodate transitions
!     between the (3) segments using 4-point Cubic Spline interpolation
!     ------------------------------------------------------------------
 
 
      INQUIRE (FILE=dtfile,EXIST=qexist)
      IF ( .NOT.qexist ) dtfile = 'RH_QG_Mie  '
                                           ! generic name used by GCM
      INQUIRE (FILE=dtfile,EXIST=qexist)
      IF ( .NOT.qexist ) CALL STOP_MODEL('setrel: no RH_QG files',255)
      CALL OPENUNIT(dtfile,kdread,.FALSE.,.TRUE.) ! formatted, old
 
      READ (KDREAD,7000) (XNR(J),J=1,31)
      DO I = 1, 880
         READ (KDREAD,7001) R633NR(I), (Q633NR(I,J),J=1,31)
      ENDDO
      READ (KDREAD,7000) (XNR(J),J=1,31)
      DO I = 1, 880
         READ (KDREAD,7001) R633NR(I), (G633NR(I,J),J=1,31)
      ENDDO
      CALL CLOSEUNIT(KDREAD)
 
      J = 880
      DO K = 299, 305
         IF ( K/=300 ) THEN
            IF ( K/=302 ) THEN
               J = J + 1
               R633NR(J) = R633NR(K)
               DO I = 1, 31
                  Q633NR(J,I) = Q633NR(K,I)
                  G633NR(J,I) = G633NR(K,I)
               ENDDO
            ENDIF
         ENDIF
      ENDDO
      DO K = 600, 606
         IF ( K/=601 ) THEN
            IF ( K/=603 ) THEN
               J = J + 1
               R633NR(J) = R633NR(K)
               DO I = 1, 31
                  Q633NR(J,I) = Q633NR(K,I)
                  G633NR(J,I) = G633NR(K,I)
               ENDDO
            ENDIF
         ENDIF
      ENDDO
 
!     Apply 13-point quadratic least-squares smoothing to large particle
!     portion of Mie Qx data to eliminate low-amplitude ripple in Q633NR
!     (Monotonic size dependence is needed for inverse Qx interpolation)
!     (Smoothing affects 4th decimal of Q633NR for large particle sizes)
!     ------------------------------------------------------------------
      DO I = 1, 31
         DO J = 1, 880
            SMOOTH(J) = Q633NR(J,I)
         ENDDO
         DO J = 881, 886
            SMOOTH(J) = SMOOTH(880)
         ENDDO
         DO J = 250, 880
            J1 = J - 2
            IF ( SMOOTH(J)>=SMOOTH(J-1) ) EXIT
         ENDDO
         DO J = J1, 880
            SUM = 4550.D0/13.D0*SMOOTH(J)
            DO K = 1, 6
               SUM = SUM + (4550.D0/13.D0-14*K*K)                       &
                     *(SMOOTH(J-K)+SMOOTH(J+K))
            ENDDO
            Q633NR(J,I) = SUM/2002.D0
         ENDDO
      ENDDO
 
!                                     Set relative humidity RHRHRH scale
!                                     ----------------------------------
      DO I = 1, 190
         RHRHRH(I) = (I-1)/100.D0
         IF ( I>91 ) RHRHRH(I) = 0.90D0 + (I-91)/1000.D0
      ENDDO
 
!         Define RH (=AW), RO, BX, RX as functions of X for NAER aerosol
!         --------------------------------------------------------------
      NRHN1 = NRHCRY(NAER) + 1
      DO I = 1, 190
         RHI = RHRHRH(I)
         RR0RHX(I) = 1.D0
         RHXMFX(I) = 1.D0
         IF ( NAER==1 ) THEN   !    Dry Sulfate refrac index and density
            XNRRHX(I) = 1.526
            RHDENS(I) = 1.760
            IF ( I<NRHN1 ) DNRX(I) = DRDSO4(RHI)
            IF ( I>=NRHN1 ) DNRX(I) = DRWSO4(RHI)
         ENDIF
         IF ( NAER==2 ) THEN   !    Dry SeaSalt refrac index and density
            XNRRHX(I) = 1.490
            RHDENS(I) = 2.165
            IF ( I<NRHN1 ) DNRX(I) = DRDSEA(RHI)
            IF ( I>=NRHN1 ) DNRX(I) = DRWSEA(RHI)
         ENDIF
         IF ( NAER==3 ) THEN   !    Dry Nitrate refrac index and density
            XNRRHX(I) = 1.554
            RHDENS(I) = 1.725
            DNRX(I) = DRXNO3(RHRHRH(I))
         ENDIF
         IF ( NAER==4 ) THEN   !    Dry Organic refrac index and density
            XNRRHX(I) = 1.526  !                  (representative value)
            RHDENS(I) = 1.5    !                  (representative value)
            IF ( I<NRHN1 ) DNRX(I) = DRDOCX(RHI)
            IF ( I>=NRHN1 ) DNRX(I) = DRWOCX(RHI)
         ENDIF
      ENDDO
 
!            Invert X, RO, BX, RX functions of (X) to be functions of RH
!            -----------------------------------------------------------
      I = 191
      FF = 1.D0
      XX = 0.D0
      IF ( NAER==1 ) GI = DWSO4(XX)
      IF ( NAER==2 ) GI = DWSEA(XX)
      IF ( NAER==3 ) GI = DWNO3(XX)
      IF ( NAER==4 ) THEN
         FF = .9995D0
         XX = (1D0-FF)**.125D0
         GI = DWOCX(XX)
      ENDIF
      DO
         I = I - 1
         FI = RHRHRH(I)
         DO K = 1, 5
            XI = XX - (FF-FI)/GI
            IF ( NAER==1 ) FF = AWSO4(XI)
            IF ( NAER==2 ) FF = AWSEA(XI)
            IF ( NAER==3 ) FF = AWNO3(XI)
            IF ( NAER==4 ) FF = AWOCX(XI)
            IF ( I>0 ) THEN
            ENDIF
            XX = XI
            IF ( NAER==1 ) GI = DWSO4(XX)
            IF ( NAER==2 ) GI = DWSEA(XX)
            IF ( NAER==3 ) GI = DWNO3(XX)
            IF ( NAER==4 ) GI = DWOCX(XX)
         ENDDO
         RHXMFX(I) = XX
         IF ( NAER==1 ) THEN      !       RH dependent Sulfate X,R,NR,RO
            RHDENS(I) = ROSO4(XX)
            RR0RHX(I) = BXSO4(XX)
            XNRRHX(I) = RXSO4(XX)
         ENDIF
         IF ( NAER==2 ) THEN      !       RH dependent SeaSalt X,R,NR,RO
            RHDENS(I) = ROSEA(XX)
            RR0RHX(I) = BXSEA(XX)
            XNRRHX(I) = RXSEA(XX)
         ENDIF
         IF ( NAER==3 ) THEN      !       RH dependent Nitrate X,R,NR,RO
            RHDENS(I) = RONO3(XX)
            RR0RHX(I) = BXNO3(XX)
            XNRRHX(I) = R1NO3(XX)
            IF ( XX>0.205D0 ) XNRRHX(I) = R2NO3(XX)
         ENDIF
         IF ( NAER==4 ) THEN      !       RH dependent Organic X,R,NR,RO
            RHDENS(I) = ROOCX(XX)
            RR0RHX(I) = BXOCX(XX)
            XNRRHX(I) = RXOCX(XX)
         ENDIF
         IF ( I<=NRHN1 ) THEN
 
!     ------------------------------------------------------------------
!     Find Qdry(r),gdry(r) from Q(m,r),g(m,r) maps for each aerosol type
!     Find Qwet(r),gwet(r) from Q(m,r),g(m,r) maps for each aerosol type
!     also locate MAXDRY,MAXWET pts where Qdry(r),Qwet(r) are at maximum
!          (M1 refers to mass fraction X of 1.0, i.e., "dry" aerosol)
!          (M0 refers to mass fraction X of 0.0, i.e., "wet" aerosol)
!     ------------------------------------------------------------------
            MAXDRY = 1
            MAXWET = 1
            QQDMAX = 0.D0
            QQWMAX = 0.D0
            XDRY = XNRRHX(1)
!     IF(MCRYON == 1) XDRY=XNRRHX(NRHN1) ! If "dry" = RHC reference line
            SDRY = XDRY*100.D0 - 129
            JDRY = SDRY
            DDRY = SDRY - JDRY
            XWET = 1.3330D0             !  Pure water Nr = "wet" aerosol
            SWET = XWET*100.D0 - 129
            JWET = SWET
            DWET = SWET - JWET
            DO I = 1, 880
               CALL SPLNI4(Q633NR,890,31,I,JDRY,DDRY,Q880M1(I))
               CALL SPLNI4(G633NR,890,31,I,JDRY,DDRY,G880M1(I))
               CALL SPLNI4(Q633NR,890,31,I,JWET,DWET,Q880M0(I))
               CALL SPLNI4(G633NR,890,31,I,JWET,DWET,G880M0(I))
               IF ( Q880M1(I)>QQDMAX ) THEN
                  QQDMAX = Q880M1(I)
                  MAXDRY = I
               ENDIF
               IF ( Q880M0(I)>QQWMAX ) THEN
                  QQWMAX = Q880M0(I)
                  MAXWET = I
               ENDIF
            ENDDO
            RQDMAX = R633NR(MAXDRY)
            RQWMAX = R633NR(MAXWET)
 
!     Define:  Qdry(r) and Qwet(r) at the reference wavelength of 550 nm
!              using refractive index off-set and size parameter scaling
!     ------------------------------------------------------------------
            XDRY = XNRRHX(1)*DNRX(1)     !      Dry aerosol Nr at 550 nm
!     IF(MCRYON == 1) XDRY=XNRRHX(NRHN1) ! If "dry" = RHC reference line
            SDRY = XDRY*100.D0 - 129
            JDRY = SDRY
            DDRY = SDRY - JDRY
            XWET = 1.3330D0*1.001179
                                 !       Pure water aerosol Nr at 550 nm
            SWET = XWET*100.D0 - 129
            JWET = SWET
            DWET = SWET - JWET
            DO I = 1, 880
               CALL SPLNI4(Q633NR,890,31,I,JDRY,DDRY,Q880N1(I))
               CALL SPLNI4(Q633NR,890,31,I,JWET,DWET,Q880N0(I))
               R550NR(I) = R633NR(I)*(0.550/0.633)
                                        !  Size shift refers Q to 550 nm
            ENDDO
            CALL SPLINE(R550NR,Q880N1,880,REFF0,Q55DRY,1.D0,1.D0,1)
            CALL SPLINE(R633NR,Q880M1,880,REFF0,Q63DRY,1.D0,1.D0,1)
 
!     Find Q(RH),g(RH) paths in Q(m,r),g(m,r) maps for seed size = REFF0
!     2-coordinate paths defined via XN0=XNRRHX(I)   RR0=REFF0*RR0RHX(I)
!     ------------------------------------------------------------------
            DO I = 1, 190
               XN0 = XNRRHX(I)
               XN1 = XN0*100.D0 - 129
               IN1 = XN1
               DWN = XN1 - IN1
               RR0 = REFF0*RR0RHX(I)
               IF ( RR0<0.01 ) RR0 = 0.01
               IF ( RR0<=3.00D0 ) XR1 = RR0*100.D0 + 1
               IF ( RR0>3.00D0 .AND. RR0<3.04D0 ) XR1 = RR0*50.0D0 + 732
               IF ( RR0>=3.04D0 .AND. RR0<=9.00D0 ) XR1 = RR0*50.0D0 +  &
                    152
               IF ( RR0>9.00D0 .AND. RR0<9.08D0 ) XR1 = RR0*25.0D0 + 662
               IF ( RR0>=9.08D0 ) THEN
                  XR1 = RR0*25.0D0 + 378
                  IF ( XR1>877.9999D0 ) XR1 = 877.9999D0
               ENDIF
               IR1 = XR1
               DWR = XR1 - IR1
               CALL SPLN44(Q633NR,890,31,IR1,DWR,IN1,DWN,QRH633(I))
               CALL SPLN44(G633NR,890,31,IR1,DWR,IN1,DWN,GRH633(I))
            ENDDO
 
!     Define Q55(RH) by tracing path in Q(m,r) map for RH dependent size
!     via 2-coordinate path XN0=XNRRHX(I)*DNRX(I), RR0=RRH(I)*(.633/.55)
!     ------------------------------------------------------------------
            DO I = 1, 190
               XN0 = XNRRHX(I)*DNRX(I)
               XN1 = XN0*100.D0 - 129
               IN1 = XN1
               DWN = XN1 - IN1
               RR0 = REFF0*RR0RHX(I)*(0.633D0/0.550D0)
               IF ( RR0<0.01 ) RR0 = 0.01
               IF ( RR0<=3.00D0 ) XR1 = RR0*100.D0 + 1
               IF ( RR0>3.00D0 .AND. RR0<3.04D0 ) XR1 = RR0*50.0D0 + 732
               IF ( RR0>=3.04D0 .AND. RR0<=9.00D0 ) XR1 = RR0*50.0D0 +  &
                    152
               IF ( RR0>9.00D0 .AND. RR0<9.08D0 ) XR1 = RR0*25.0D0 + 662
               IF ( RR0>=9.08D0 ) THEN
                  XR1 = RR0*25.0D0 + 378
                  IF ( XR1>877.9999D0 ) XR1 = 877.9999D0
               ENDIF
               IR1 = XR1
               DWR = XR1 - IR1
               CALL SPLN44(Q633NR,890,31,IR1,DWR,IN1,DWN,RHQ550(I))
               RHREFF(I) = RR0RHX(I)*REFF0
            ENDDO
 
!        Aerosol liquid water content is in kg/m2 per unit optical depth
!     of dry aerosol with aerosol effective radius expressed in microns.
!     ------------------------------------------------------------------
 
            DO I = 1, 190
               RHTAUF(I) = (RHQ550(I)/Q55DRY)*RR0RHX(I)**2
               AERMAS = 1.33333333D0*RHREFF(I)*RHDENS(I)/RHQ550(I)      &
                        *RHTAUF(I)
               RHTGM2(I) = AERMAS
               RHDGM2(I) = AERMAS*RHXMFX(I)
               RHWGM2(I) = RHTGM2(I) - RHDGM2(I)
               TAUM2G(I) = 0.75D0/RHDENS(1)/RHREFF(1)*RHQ550(1)         &
                           *RHTAUF(I)
               ANBCM2(I) = TAUM2G(I)/(1.5080*RHQ550(I)*RHREFF(I)**2)
            ENDDO
 
!     Determination of RH dependent Mie scattering tables for GCM input.
!     Find equivalent aersol dry sizes (RD1,RD2) and wet sizes (RW1,RW2)
!     and corresponding weights to match the RH dependent Q(r) and g(r).
!     Fits made to form: QRH=X*[Y*QD1+(1-Y)*QD2]+(1-X)*[Z*WD1+(1-Z)*WD2]
!     ------------------------------------------------------------------
            J1 = MAXWET
            JHIMAX = 881 - MAXWET
            K1 = MAXDRY
            KHIMAX = 881 - MAXDRY
            NP = 190 - NRHN1 + 1
            DO I = 1, 190
               RHRHI = RHRHRH(I)
               XFDRY = RHXMFX(I)
               REFFI = RHREFF(I)
               RRH = RR0RHX(I)*REFF0
               GRH = GRH633(I)
               QRH = QRH633(I)
               QD1 = QRH
               QD2 = QRH
               QW1 = QRH
               QW2 = QRH
               IF ( QW1>QQWMAX ) QW1 = QQWMAX
               IF ( QW2>QQWMAX ) QW2 = QQWMAX
               CALL SPLINV(R633NR,Q880M0,MAXWET,RW1,QW1,1.D0,1.D0,1)
               CALL SPLINV(R633NR(J1),Q880M0(J1),JHIMAX,RW2,QW2,1.D0,   &
                           1.D0,1)
               CALL SPLINE(R633NR,G880M0,880,RW1,GW1,1.D0,1.D0,1)
               CALL SPLINE(R633NR,G880M0,880,RW2,GW2,1.D0,1.D0,1)
               IF ( I>=NRHN1 .AND. QRH>QQWMAX ) THEN
                  QD1 = QQWMAX + (QRH-QQWMAX)/XFDRY
                                    ! QD1 such that  QRH=X*QD1+(1-X)*QW1
                  QD2 = 2.3D0       ! 2 dry sizes are used if QD1>QQWMAX
               ENDIF
               CALL SPLINV(R633NR,Q880M1,MAXDRY,RD1,QD1,1.D0,1.D0,1)
               CALL SPLINV(R633NR(K1),Q880M1(K1),KHIMAX,RD2,QD2,1.D0,   &
                           1.D0,1)
               CALL SPLINE(R633NR,G880M1,880,RD1,GD1,1.D0,1.D0,1)
               CALL SPLINE(R633NR,G880M1,880,RD2,GD2,1.D0,1.D0,1)
 
               IF ( I<NRHN1 ) THEN
                                 !          Pure dry aerosol region (1)
                  WTX = 1.D0
                  WTY = 1.D0
                  IF ( REFF0>RQDMAX ) WTY = 0.D0
                  WTZ = 1.D0
               ELSE   !         Dry/wet weighted average regions (2)-(4)
                  IF ( QRH<=QQWMAX .AND. REFFI<RW1 ) THEN
                                              !   Small-size region (2)
                     WTZ = 1.D0
                     WTY = 1.D0
                     WTX = (GRH-GW1)/(GD1-GW1)
                  ENDIF
                  IF ( QRH>QQWMAX ) THEN      !  Medium-size region (3)
!     Fit form: QRH=X*(Y*QD1+(1-Y)*QD2)+(1-X)*QWmax   QRH=/QD1=/QD2=/QW1
                     WTZ = 1.D0
                     WTY = ((GRH-GD2)*QRH*QD2+(GD2-GW1)                 &
                           *QD2*QW1+(GW1-GRH)*QRH*QW1)                  &
                           /((GD1-GRH)*QRH*QD1+(GRH-GD2)                &
                           *QRH*QD2+(GW1-GD1)*QD1*QW1+(GD2-GW1)*QD2*QW1)
                     WTX = (QRH-QW1)/(WTY*(QD1-QD2)+(QD2-QW1))
                  ENDIF
                  IF ( QRH<=QQWMAX .AND. REFFI>RW1 ) THEN
                                               !  Large size region (4)
                     WTY = 0.D0
                     WTZ = 0.D0
                     WTX = (GRH-GW2)/(GD2-GW2)
                  ENDIF
               ENDIF
               IF ( REFFI>RQWMAX .AND. RHRHI>0.995 ) THEN
                                                  ! High RH region (5)
                  WTY = 0.D0
                  WTX = XFDRY
                  WTZ = ((GRH-GW2)-(GD2-GW2)*WTX)/((1.D0-WTX)*(GW1-GW2))
               ENDIF
 
               VD1 = WTX*WTY
               VD2 = WTX*(1.D0-WTY)
               VW1 = WTZ*(1.D0-WTX)
               VW2 = (1.D0-WTZ)*(1.D0-WTX)
               RD1 = MIN(RD1,10.D0)
               RD2 = MIN(RD2,10.D0)
               RW1 = MIN(RW1,10.D0)
               RW2 = MIN(RW2,10.D0)
 
!     Computed weight factors are for Lab reference wavelength of 633nm.
!     Rescale spectral extinction to 550 nm   renormalize weight factors
!     ------------------------------------------------------------------
               CALL SPLINE(R550NR,Q880N1,880,RD1,Q550,1.D0,1.D0,1)
               CALL SPLINE(R633NR,Q880M1,880,RD1,Q633,1.D0,1.D0,1)
               WD1 = VD1*(Q550/Q633)
               CALL SPLINE(R550NR,Q880N1,880,RD2,Q550,1.D0,1.D0,1)
               CALL SPLINE(R633NR,Q880M1,880,RD2,Q633,1.D0,1.D0,1)
               WD2 = VD2*(Q550/Q633)
               CALL SPLINE(R550NR,Q880N0,880,RW1,Q550,1.D0,1.D0,1)
               CALL SPLINE(R633NR,Q880M0,880,RW1,Q633,1.D0,1.D0,1)
               WW1 = VW1*(Q550/Q633)
               CALL SPLINE(R550NR,Q880N0,880,RW2,Q550,1.D0,1.D0,1)
               CALL SPLINE(R633NR,Q880M0,880,RW2,Q633,1.D0,1.D0,1)
               WW2 = VW2*(Q550/Q633)
               SUMW = WD1 + WD2 + WW1 + WW2
               W1 = WD1/SUMW
               W2 = WD2/SUMW
               W3 = WW1/SUMW
               W4 = WW2/SUMW
 
!     ------------------------------------------------------------------
!     Tabulate relative humidity dependent solar, thermal Mie scattering
!     parameters SRHQEX,SRHQSC,SRHQCS, TRHQAB for each aerosol type NAER
!     These are mass weighted averages of equivalent dry and wet aerosol
!     parameters for sizes matching the relative humidity dependent Q(r)
!     ------------------------------------------------------------------
 
               N0 = 0           !      Select Mie parameters for Sulfate
               IF ( NAER==2 ) N0 = 22
                                !      Select Mie parameters for SeaSalt
               IF ( NAER==3 ) N0 = 44
                                !      Select Mie parameters for Nitrate
               IF ( NAER==4 ) N0 = 88
                                !      Select Mie parameters for Organic
               N1 = N0 + 1
               DO K = 1, 6                    !   SW dry sizes RD1   RD2
                  DO N = 1, 22
                     NN = N0 + N
                     WTS = FRSULF(NAER)
                     WTA = 1.D0 - WTS
                     QXAERN(N) = SRUQEX(K,NN)*WTA + SRUQEX(K,N)*WTS
                     QSAERN(N) = SRUQSC(K,NN)*WTA + SRUQSC(K,N)*WTS
                     QGAERX = SRUQCB(K,NN)*SRUQSC(K,NN)                 &
                              *WTA + SRUQCB(K,N)*SRUQSC(K,N)*WTS
                     QGAERN(N) = QGAERX/QSAERN(N)
                  ENDDO
                  CALL SPLINE(REFU22,QXAERN,22,RD1,SR1QEX(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QSAERN,22,RD1,SR1QSC(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QGAERN,22,RD1,SR1QCB(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QXAERN,22,RD2,SR2QEX(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QSAERN,22,RD2,SR2QSC(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QGAERN,22,RD2,SR2QCB(K),1.D0,1.D0, &
                              1)
               ENDDO
 
               DO K = 1, 33                   !   LW dry sizes RD1   RD2
                  DO N = 1, 22
                     NN = N0 + N
                     WTS = FRSULF(NAER)
                     WTA = 1.D0 - WTS
                     QXAERN(N) = TRUQEX(K,NN)*WTA + TRUQEX(K,N)*WTS
                     QSAERN(N) = TRUQSC(K,NN)*WTA + TRUQSC(K,N)*WTS
                     QGAERX = TRUQCB(K,NN)*TRUQSC(K,NN)                 &
                              *WTA + TRUQCB(K,N)*TRUQSC(K,N)*WTS
                     QGAERN(N) = QGAERX/(QSAERN(N)+1D-10)
                  ENDDO
                  CALL SPLINE(REFU22,QXAERN,22,RD1,TR1QEX(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QSAERN,22,RD1,TR1QSC(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QGAERN,22,RD1,TR1QCB(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QXAERN,22,RD2,TR2QEX(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QSAERN,22,RD2,TR2QSC(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QGAERN,22,RD2,TR2QCB(K),1.D0,1.D0, &
                              1)
               ENDDO
               CALL SPLINE(REFU22,Q55U22(N1),22,RD1,Q55RH1,1.D0,1.D0,1)
               CALL SPLINE(REFU22,Q55U22(N1),22,RD2,Q55RH2,1.D0,1.D0,1)
 
               N0 = 66         !   Select Mie parameters for pure water
               N1 = N0 + 1
               DO K = 1, 6                   !   SW wet sizes RW1   RW2
                  DO N = 1, 22
                     NN = N0 + N
                     QXAERN(N) = SRUQEX(K,NN)
                     QSAERN(N) = SRUQSC(K,NN)
                     QGAERN(N) = SRUQCB(K,NN)
                  ENDDO
                  CALL SPLINE(REFU22,QXAERN,22,RW1,SR3QEX(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QSAERN,22,RW1,SR3QSC(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QGAERN,22,RW1,SR3QCB(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QXAERN,22,RW2,SR4QEX(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QSAERN,22,RW2,SR4QSC(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QGAERN,22,RW2,SR4QCB(K),1.D0,1.D0, &
                              1)
               ENDDO
 
               DO K = 1, 33                  !   LW wet sizes RW1   RW2
                  DO N = 1, 22
                     NN = N0 + N
                     QXAERN(N) = TRUQEX(K,NN)
                     QSAERN(N) = TRUQSC(K,NN)
                     QGAERN(N) = TRUQCB(K,NN)
                  ENDDO
                  CALL SPLINE(REFU22,QXAERN,22,RW1,TR3QEX(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QSAERN,22,RW1,TR3QSC(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QGAERN,22,RW1,TR3QCB(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QXAERN,22,RW2,TR4QEX(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QSAERN,22,RW2,TR4QSC(K),1.D0,1.D0, &
                              1)
                  CALL SPLINE(REFU22,QGAERN,22,RW2,TR4QCB(K),1.D0,1.D0, &
                              1)
               ENDDO
               CALL SPLINE(REFU22,Q55U22(N1),22,RW1,Q55RH3,1.D0,1.D0,1)
               CALL SPLINE(REFU22,Q55U22(N1),22,RW2,Q55RH4,1.D0,1.D0,1)
 
                        !      Weighted GCM SW Mie scattering parameters
               DO K = 1, 6
                  SRHQEX(K,I) = W1*SR1QEX(K) + W2*SR2QEX(K)             &
                                + W3*SR3QEX(K) + W4*SR4QEX(K)
                  SRHQSC(K,I) = W1*SR1QSC(K) + W2*SR2QSC(K)             &
                                + W3*SR3QSC(K) + W4*SR4QSC(K)
                  QSCQCB = W1*SR1QCB(K)*SR1QSC(K) + W2*SR2QCB(K)        &
                           *SR2QSC(K) + W3*SR3QCB(K)*SR3QSC(K)          &
                           + W4*SR4QCB(K)*SR4QSC(K)
                  SRHQCB(K,I) = QSCQCB/SRHQSC(K,I)
               ENDDO
                        !      Weighted GCM LW Mie scattering parameters
               DO K = 1, 33
                  TRHQEX(K) = W1*TR1QEX(K) + W2*TR2QEX(K) + W3*TR3QEX(K)&
                              + W4*TR4QEX(K)
                  TRHQSC(K) = W1*TR1QSC(K) + W2*TR2QSC(K) + W3*TR3QSC(K)&
                              + W4*TR4QSC(K)
                  QSCQCB = W1*TR1QCB(K)*TR1QSC(K) + W2*TR2QCB(K)        &
                           *TR2QSC(K) + W3*TR3QCB(K)*TR3QSC(K)          &
                           + W4*TR4QCB(K)*TR4QSC(K)
                  TRHQCB(K) = QSCQCB/TRHQSC(K)
                  TRHQAB(K,I) = TRHQEX(K) - TRHQSC(K)
               ENDDO
 
               COSBAR(I) = SRHQCB(6,I)
               PIZERO(I) = SRHQSC(6,I)/SRHQEX(6,I)
               ANGSTR(I) = -(1.D0-SRHQEX(5,I))/(0.550D0-0.815D0)
 
!              Transfer EQUIVALENCEd SETREL output information to RHDATA
               DO J = 1, 15
                  RHDATA(I,J) = RHINFO(I,J)
               ENDDO
            ENDDO
 
!                                                      Diagnostic output
            IF ( AM_I_ROOT() ) THEN
               DO I = 1, 190
                  IF ( I==1 ) WRITE (99,6000) AERTYP(NAER), NAER, REFF0
                  IF ( I==82 ) WRITE (99,6000) AERTYP(NAER), NAER, REFF0
                  IF ( I==137 ) WRITE (99,6000) AERTYP(NAER), NAER,     &
                                       REFF0
                  IF ( I>=27 ) THEN
                     WRITE (99,6100) I, (RHINFO(I,N),N=1,15),           &
                                     SRHQEX(6,I), SRHQEX(5,I),          &
                                     SRHQEX(1,I), TRHQAB(1,I)
 6100                FORMAT (I3,F5.3,18F8.4)
                  ENDIF
               ENDDO
            ENDIF
            EXIT
         ENDIF
      ENDDO
 7000 FORMAT (12X,F5.3,30F8.3)
 7001 FORMAT (3X,F6.2,31F8.5)
 6000 FORMAT (T90,A8,'  NAER=',I2,'   REFF0=',F5.2//                    &
             &'      RH  RHTAUF  RHREFF  RHWGM2  RHDGM2  RHTGM2  RHXMFX'&
              ,                                                         &
             &'  RHDENS  RHQ550  TAUM2G  XNRRHX  ANBCM2  COSBAR  PIZERO'&
              ,'  ANGSTR SRHQEX6 SRHQEX5 SRHQEX1 TRHQAB1')
 
      END SUBROUTINE SETREL
 
      REAL*8 FUNCTION RHDTNA(TK,NA)
      IMPLICIT NONE
 
      INTEGER, INTENT(IN)  ::  NA
      REAL*8, INTENT(IN)  ::  TK
!     functions
      REAL*8 RHDSO4, RHDSEA, RHDNO3, RHDOCX
 
      RHDSO4(TK) = MIN(1.D0,0.80D0*EXP(25.D0*(298.0D0-TK)/(298.0D0*TK)))
      RHDSEA(TK) = MIN(1.D0,0.75D0*EXP(80.D0*(298.0D0-TK)/(298.0D0*TK)))
      RHDNO3(TK) = MIN(1.D0,0.62D0*EXP(852.D0*(298.0D0-TK)/(298.0D0*TK))&
                   )
      RHDOCX(TK) = MIN(1.D0,0.80D0*EXP(25.D0*(298.0D0-TK)/(298.0D0*TK)))
 
      IF ( NA==1 ) RHDTNA = RHDSO4(TK)
      IF ( NA==2 ) RHDTNA = RHDSEA(TK)
      IF ( NA==3 ) RHDTNA = RHDNO3(TK)
      IF ( NA==4 ) RHDTNA = RHDOCX(TK)
      END FUNCTION RHDTNA
 
      REAL*8 FUNCTION COMPUTE(f1,f2,f3,f4,var)RESULT(RESULT)
      IMPLICIT NONE
      REAL*8, INTENT(IN)  ::  f1, f2, f3, f4
      REAL*8, INTENT(IN)  ::  var
      REAL*8  ::  F21, F32, F43, F3221, F4332
      REAL*8  ::  A, B, C, D, XF, FFCUSP, XEXM, CUSPWT
      REAL*8  ::  FFLINR
      REAL*8, PARAMETER  ::  CUSPWM = 0.5
      REAL*8, PARAMETER  ::  CUSPWE = 0.5
      F21 = (F2-F1)
      F32 = (F3-F2)
      F43 = (F4-F3)
      F3221 = (F32+F21)*0.5D0
      F4332 = (F43+F32)*0.5D0
      A = F2
      B = F3221
      C = 3.D0*F32 - F3221 - F3221 - F4332
      D = (F3221+F4332-F32-F32)
      XF = var
      FFCUSP = A + XF*(B+XF*(C+XF*D))
      XEXM = (1-2*XF)**2        ! = XE**2
      CUSPWT = (1.0-XEXM)*CUSPWM + XEXM*CUSPWE
      FFLINR = A + XF*F32
      RESULT = FFCUSP*CUSPWT + FFLINR*(1.D0-CUSPWT)
      END FUNCTION COMPUTE
 
      REAL*8 FUNCTION COMPUTE2(f1,f2,f3,f4,var)RESULT(RESULT)
      IMPLICIT NONE
      REAL*8, INTENT(IN)  ::  f1, f2, f3, f4
      REAL*8, INTENT(IN)  ::  var
      REAL*8  ::  F21, F32, F43, F3221, F4332
      REAL*8  ::  A, B, C, D, XF, FFCUSP, XEXM, XE
      F21 = (F2-F1)
      F32 = (F3-F2)
      F43 = (F4-F3)
      F3221 = (F32+F21)*0.5D0
      F4332 = (F43+F32)*0.5D0
      A = F2
      B = F3221
      C = 3.D0*F32 - F3221 - F3221 - F4332
      D = (F3221+F4332-F32-F32)
      XF = var
      FFCUSP = A + XF*(B+XF*(C+XF*D))
      XE = 1.D0 - XF - XF
      IF ( XE<0.0 ) XE = -XE
      XEXM = XE**2
!=1   CUSPWT=(1.D0-XEXM)*CUSPWM+XEXM*CUSPWE
!nu   FFLINR=A+XF*F32
      RESULT = FFCUSP
      END FUNCTION COMPUTE2

      MODULE O3MOD
!@sum O3mod administers reading of ozone files
!@auth M. Kelley and original development team
      USE TIMESTREAM_MOD, ONLY:TIMESTREAM
      IMPLICIT NONE
      SAVE 
!@var O3stream interface for reading and time-interpolating O3 files
!@+   See usage notes in timestream_mod
      TYPE (TIMESTREAM) :: O3stream, delta_O3stream
#ifdef HIGH_FREQUENCY_O3_INPUT
      TYPE (TIMESTREAM) :: OxHFstream, PSFforO3stream
#endif
#ifdef GCAP
      REAL*8, ALLOCATABLE  ::  save_to3(:,:)
#endif
!@dbparam use_sol_Ox_cycle if =1, a cycle of ozone is appled to
!@+ o3year, as a function of the solar constant cycle.
      INTEGER  ::  use_sol_Ox_cycle = 0
      REAL*8  ::  S0min, S0max
 
!@dbparam ozone_use_ppm_interp = 1 uses ppm interpolation in the
!@+ timestream. Otherwise uses linm2m.
      INTEGER  ::  ozone_use_ppm_interp = 1
 
!@var have_o3_file whether an O3file was specified in the rundeck
      LOGICAL  ::  have_o3_file
 
!@param NLO3_traditional assumed number of layers in ozone data files.
      INTEGER, PARAMETER  ::  NLO3_TRADITIONAL = 49
!@var NLO3 number of layers in ozone data files, as read from file.
      INTEGER  ::  NLO3 = 0
!@var PLBO3_traditional assumed edge pressures in O3 input file.
      REAL*8  ::  PLBO3_traditional(NLO3_TRADITIONAL+1)                 &
                & = (/984D0,934D0,854D0,720D0,550D0,390D0,285D0,210D0,  &
                & 150D0,125D0,100D0,80D0,60D0,55D0,50D0,45D0,40D0,35D0, &
                & 30D0,25D0,20D0,15D0,10.D0,7.D0,5.D0,4.D0,3.D0,2.D0,   &
                & 1.5D0,1.D0,7D-1,5D-1,4D-1,3D-1,2D-1,1.5D-1,1D-1,7D-2, &
                & 5D-2,4D-2,3D-2,2D-2,1.5D-2,1D-2,7D-3,5D-3,4D-3,3D-3,  &
                & 1D-3,1D-7/)
!@var PLBO3 edge pressures in O3 input file, as read from file
      REAL*8, ALLOCATABLE  ::  PLBO3(:)
 
 
      CONTAINS
 
      SUBROUTINE UPDO3D(JYEARO,JJDAYO,O3JDAY,O3JREF)
      USE DICTIONARY_MOD
      USE RESOLUTION, ONLY:psf
      USE DOMAIN_DECOMP_ATM, ONLY:grid, GETDOMAINBOUNDS
      USE TIMESTREAM_MOD, ONLY:INIT_STREAM, READ_STREAM,                &
        & GETNAME_FIRSTFILE
      USE PARIO, ONLY:PAR_OPEN, PAR_CLOSE, READ_DIST_DATA,              &
        & VARIABLE_EXISTS, GET_DIMLEN, READ_DATA
      USE FILEMANAGER, ONLY:FILE_EXISTS
      IMPLICIT NONE
      INTEGER, INTENT(IN)  ::  JYEARO, JJDAYO
      REAL*8, DIMENSION(:,:,:), POINTER  ::  o3jday, o3jref
 
      INTEGER  ::  i, j, l, jyearx, fid
      LOGICAL, SAVE  ::  init = .FALSE.
      LOGICAL  ::  cyclic, exists
      REAL*8, ALLOCATABLE  ::  o3arr(:,:,:)
      CHARACTER(LEN=6)  ::  method
      CHARACTER(LEN=32)  ::  fname1st
 
      INTEGER  ::  j_0, j_1, i_0, i_1
 
      CALL GETDOMAINBOUNDS(grid,J_STRT=j_0,J_STOP=j_1,I_STRT=i_0,       &
                         & I_STOP=i_1)
 
      jyearx = ABS(jyearo)
 
      IF ( .NOT.init ) THEN
         init = .TRUE.
 
         have_o3_file = FILE_EXISTS('O3file')
 
         IF ( have_o3_file ) THEN
 
          ! Initialize the timestream for the O3 data file:
            cyclic = jyearo<0
 
            CALL SYNC_PARAM("ozone_use_ppm_interp",ozone_use_ppm_interp)
            IF ( ozone_use_ppm_interp==1 ) THEN
               method = 'ppm'
            ELSE
               method = 'linm2m'
            ENDIF
            CALL INIT_STREAM(grid,O3stream,'O3file','O3',0D0,1D30,      &
                           & TRIM(method),jyearx,jjdayo,cyclic=cyclic)
          ! query the layering
            CALL GETNAME_FIRSTFILE(O3stream,fname1st)
            fid = PAR_OPEN(grid,TRIM(fname1st),'read')
            IF ( VARIABLE_EXISTS(grid,fid,'ple') ) THEN
               nlo3 = GET_DIMLEN(grid,fid,'ple') - 1
                                                ! coord var but one less
               IF ( nlo3/=GET_DIMLEN(grid,fid,'plm') )                  &
                   & CALL STOP_MODEL('ple/plm dim problem in O3file',   &
                  & 255)
               ALLOCATE (plbo3(nlo3+1))
               CALL READ_DATA(grid,fid,'ple',plbo3,BCAST_ALL=.TRUE.)
            ELSE
               CALL STOP_MODEL('missing ple info in o3file',255)
            ENDIF
            CALL PAR_CLOSE(grid,fid)
         ELSE
            nlo3 = NLO3_TRADITIONAL
            ALLOCATE (plbo3(nlo3+1))
            plbo3(:) = plbo3_traditional(:)
         ENDIF
 
#ifdef GCAP
         ALLOCATE (save_to3(grid%I_STRT:grid%I_STOP,grid%J_STRT:grid%   &
                 & J_STOP))
         save_to3 = 0.
#endif
 
         ALLOCATE (o3jday(nlo3,grid%I_STRT:grid%I_STOP,grid%J_STRT:grid%&
                 & J_STOP))
         o3jday = 0.
 
         ALLOCATE (o3jref(NLO3_TRADITIONAL,grid%I_STRT:grid%I_STOP,     &
                 & grid%J_STRT:grid%J_STOP))
         o3jref = 0.
 
! The next line is brought over from the original UPDO3D. I think
! it is to prevent "losing" some ozone in the REPART interpolation
! if the (fixed) lowest O3 level pressure is at lower pressure than
! the the (fixed) lowest nominal model pressure:
         IF ( plbo3(1)<psf ) plbo3(1) = psf
         IF ( plbo3_traditional(1)<psf ) plbo3_traditional(1) = psf
 
! Read the 3D field for O3 RCOMPX reference calls.
! (There is no need to allow for this one on flexible # of levels)
         INQUIRE (FILE='Ox_ref',EXIST=exists)
         IF ( exists ) THEN
            ALLOCATE (o3arr(grid%I_STRT_HALO:grid%I_STOP_HALO,          &
                    & grid%J_STRT_HALO:grid%J_STOP_HALO,                &
                    & NLO3_TRADITIONAL))
            fid = PAR_OPEN(grid,'Ox_ref','read')
            CALL READ_DIST_DATA(grid,fid,'O3',o3arr)
            CALL PAR_CLOSE(grid,fid)
            DO j = j_0, j_1
               DO i = i_0, i_1
                  O3JREF(:,I,J) = O3ARR(I,J,:)
               ENDDO
            ENDDO
            DEALLOCATE (o3arr)
                            ! note quick deallocation as will be resized below
         ENDIF
 
      ENDIF  ! end init
 
      IF ( have_o3_file ) THEN
         ALLOCATE (o3arr(grid%I_STRT_HALO:grid%I_STOP_HALO,             &
                 & grid%J_STRT_HALO:grid%J_STOP_HALO,nlo3))
                                                          ! resizing
 
         CALL READ_STREAM(grid,O3stream,jyearx,jjdayo,o3arr)
         DO j = j_0, j_1
            DO i = i_0, i_1
               O3JDAY(:,I,J) = O3ARR(I,J,:)
            ENDDO
         ENDDO
         DEALLOCATE (o3arr)
      ENDIF ! there is an o3file
 
      END SUBROUTINE UPDO3D
 
 
#ifdef HIGH_FREQUENCY_O3_INPUT
      SUBROUTINE UPDO3D_HIGHFREQUENCY(JYEARO,JJDAYO,                    &
                                    & o3jday_HF_modelLevels)
      USE RESOLUTION, ONLY:LM, mtop, MFIX, MFRAC, mfixs
      USE DOMAIN_DECOMP_ATM, ONLY:grid, GETDOMAINBOUNDS
      USE TIMESTREAM_MOD, ONLY:INIT_STREAM, READ_STREAM
      USE PARIO, ONLY:PAR_OPEN, PAR_CLOSE, READ_DIST_DATA
      USE FILEMANAGER, ONLY:FILE_EXISTS
      USE ATM_COM, ONLY:pedn
      USE CONSTANT, ONLY:BYGRAV, tf, rgas, mb2kg, KG2MB
      IMPLICIT NONE
      INTEGER, INTENT(IN)  ::  JYEARO, JJDAYO
      REAL*8, DIMENSION(:,:,:), POINTER  ::  o3jday_HF_modelLevels
 
      INTEGER  ::  i, j, l, jyearx, fid
      LOGICAL, SAVE  ::  init = .FALSE.
      LOGICAL  ::  cyclic, exists
      REAL*8, ALLOCATABLE  ::  OxHFarr(:,:,:)
      REAL*8, ALLOCATABLE  ::  psf4o3arr(:,:)
      REAL*8, DIMENSION(LM)  ::  OxHFarr_Interpolated,                  &
                               & OxHFarr_Converted, airmass
      REAL*8, DIMENSION(LM+1)  ::  modelPressureBottoms,                &
                                 & filePressureBottoms
      REAL*8  ::  numerator, denominator, mvar
 
      INTEGER  ::  j_0, j_1, i_0, i_1
 
      CALL GETDOMAINBOUNDS(grid,J_STRT=j_0,J_STOP=j_1,I_STRT=i_0,       &
                         & I_STOP=i_1)
 
      ALLOCATE (OxHFarr(grid%I_STRT_HALO:grid%I_STOP_HALO,              &
              & grid%J_STRT_HALO:grid%J_STOP_HALO,LM))
      ALLOCATE (psf4o3arr(grid%I_STRT_HALO:grid%I_STOP_HALO,            &
              & grid%J_STRT_HALO:grid%J_STOP_HALO))
 
      jyearx = ABS(jyearo)
 
      IF ( .NOT.init ) THEN
         init = .TRUE.
         ALLOCATE (o3jday_HF_modelLevels(LM,grid%I_STRT:grid%I_STOP,    &
                 & grid%J_STRT:grid%J_STOP))
         o3jday_HF_modelLevels = 0.
 
         cyclic = jyearo<0
 
         CALL INIT_STREAM(grid,OxHFstream,'OxHFfile','Ox',0D0,1D30,     &
                         &'none',jyearx,jjdayo,cyclic=cyclic)
         CALL INIT_STREAM(grid,PSFforO3stream,'OxHFfile','p_surf',0D0,  &
                        & 1D30,'none',jyearx,jjdayo,cyclic=cyclic)
 
      ENDIF  ! end init
 
      CALL READ_STREAM(grid,OxHFstream,jyearx,jjdayo,OxHFarr)
      CALL READ_STREAM(grid,PSFforO3stream,jyearx,jjdayo,psf4o3arr)
 
      DO j = j_0, j_1
         DO i = i_0, i_1
            modelPressureBottoms(:) = pedn(:,i,j)
        ! approximate the air mass concurrent with ozone input:
            mvar = psf4o3arr(i,j)*mb2kg - mfixs - mtop
            filePressureBottoms(LM+1) = mtop
            DO L = LM, 1, -1
               airmass(L) = MFIX(L) + mvar*MFRAC(L)
               filePressureBottoms(L) = filePressureBottoms(L+1) + 
 
            ENDDO
!, adjust
        ! bottom level edge (similar to how routine UPDO3D does:
        ! if(plbo3(1) < psf) plbo3(1) = psf ) Though we are not altering
        ! the airmass at the same time. Should we?
            filePressureBottoms(1) = MAX(filePressureBottoms(1),        &
                                   & modelPressureBottoms(1))
 
        ! Convert units from pppv (volume mixing ratio input) to atm-cm.
        ! For now using the approx. input file air mass, but this could
        ! be changed such that the vmr is the fundamental quantity and
        ! the mass is based on current model air mass...
        !
        ! Explanation of conversion:
        ! Get "numerator" which is the amount of ozone we are inputting
        ! in a given layer in units of kg(O3)/m2, where m2 is horizontal
        ! surface area. Then get a "denominator" that is the density of
        ! ozone at standard atmospheric pressure and 0 deg C in units of
        ! kg(O3)/m3. This is a constant.
        ! Then the num/den ratio has units of {kg/m2} / {kg/m3} = m and
        ! represents "how much" ozone you would have if it were at those
        ! pressure and temperature conditions, expressed as a thickness.
        ! Call that an "atm-m". Then in the end conversion to atm-cm or
        ! Dobson Unit are just powers of 10.
        !
        ! The numerator is obtained starting with read-in mole fraction
        ! and denote a mole as n, our starting units are: n(O3)/n(air).
        ! n(O3)/n(air) * [ratio of molecular weights, ozone to air] is:
        ! n(O3)/n(air) * [48. g(O3)/n(O3)  /  28.9655d g(air)/n(air)]
        !  --> g(O3)/g(air) = kg(O3)/kg(air). So now we have a mass
        ! mixing ratio. Then multiply by the air mass in kg/m2:
        ! kg(O3)/kg(air) * [kg(air)/m2] --> kg(O3)/m2.
        !
        ! The demoninator is obtained starting with the density of ozone
        ! at 1 atmosphere and 0 deg C. At those conditions, air density
        ! is p/RT. I.e. p, R, T are constants here and reference Earth's
        ! atmosphere: p=101325 Pa, R=rgas in J kg-1 K-1, T=tf in K,
        ! so units work out to kg(air)/m3. Convert from air to ozone
        ! again using the ratio of molecular weights, e.g. on Earth:
        ! 1.2922 kg(air)/m3 * [48. g(O3)/n(O3) / 28.9655d g(air)/n(air)]
        ! --> 2.1415 kg(O3)/m3. Note that this ratio of molecular weights
        ! appears in the numerator and denominator so is skipped below.
        !
        ! Now do numerator/denominator and obtain atm-m units, and
        ! multiply by 100 to get the desired atm-cm units. This is the
        ! cm thickness of O3 one would have under those those specific
        ! atmoserpheric conditions. All that results in just:
 
            DO L = 1, LM
               numerator = OxHFarr(i,j,L)*airmass(L)
                                              ! kg O3 / m2 we have
               denominator = 101325.D0/(rgas*tf)
                                              ! kg O3 / m3 @ 1 atm and 0 deg C
               OxHFarr_Converted(L) = 1.D2*numerator/denominator
            ENDDO
 
        ! Now, interpolate vertically, but this interpolation is not onto
        ! the rad code O3 levels, it is just an adjustment over the same
        ! LM levels but allowing for different surface pressure than was
        ! concurrent when this model input was saved from a previous run:
                                                                     ! IN
            CALL REPART(OxHFarr_Converted,filePressureBottoms,LM+1,     &
                      & OxHFarr_Interpolated,modelPressureBottoms,LM+1)
                                                                     ! OUT
        ! save for use in rad code proper:
            o3jday_HF_modelLevels(:,i,j) = OxHFarr_Interpolated(:)
         ENDDO
      ENDDO
 
      DEALLOCATE (OxHFarr,psf4o3arr)
 
      END SUBROUTINE UPDO3D_HIGHFREQUENCY
#endif /* HIGH_FREQUENCY_O3_INPUT */
 
 
      SUBROUTINE UPDO3D_SOLAR(jjdayo,S0,o3jday)
!@sum UPDO3D_solar adds solar cycle variability to O3JDAY
      USE DICTIONARY_MOD
      USE DOMAIN_DECOMP_ATM, ONLY:grid, GETDOMAINBOUNDS, AM_I_ROOT
      USE TIMESTREAM_MOD, ONLY:INIT_STREAM, READ_STREAM
      USE PARIO, ONLY:PAR_OPEN, PAR_CLOSE, READ_DATA
      IMPLICIT NONE
      INTEGER  ::  jjdayo
      REAL*8  ::  S0
      REAL*8, DIMENSION(:,:,:), POINTER  ::  o3jday
!@var delta_o3_now the difference in O3 between solar max and solar min,
!@+   interpolated to the current day
      REAL*8, ALLOCATABLE  ::  delta_o3_now(:,:,:)
!@var add_sol is [S00WM2(now)-1/2(S00WM2min+S00WM2max)]/
!@+ [S00WM2max-S00WM2min] so that O3(altered) = O3(default) +
!@+ add_sol*delta_O3_now
      REAL*8  ::  add_sol
      LOGICAL, SAVE  ::  init = .FALSE.
      INTEGER  ::  i, j, l, fid, jyearx
 
      INTEGER  ::  j_0, j_1, i_0, i_1
 
      jyearx = 2000 ! nominal year
 
      IF ( .NOT.init ) THEN
         init = .TRUE.
 
         CALL SYNC_PARAM("use_sol_Ox_cycle",use_sol_Ox_cycle)
 
         IF ( use_sol_Ox_cycle/=1 ) RETURN
 
         fid = PAR_OPEN(grid,'delta_O3','read')
         CALL READ_DATA(grid,fid,'S0min',S0min,BCAST_ALL=.TRUE.)
         CALL READ_DATA(grid,fid,'S0max',S0max,BCAST_ALL=.TRUE.)
         CALL PAR_CLOSE(grid,fid)
 
         CALL INIT_STREAM(grid,delta_O3stream,'delta_O3','O3',-1D30,    &
                        & 1D30,'linm2m',jyearx,jjdayo,CYCLIC=.TRUE.)
 
      ENDIF
 
      IF ( use_sol_Ox_cycle/=1 ) RETURN
 
      CALL GETDOMAINBOUNDS(grid,J_STRT=j_0,J_STOP=j_1,I_STRT=i_0,       &
                         & I_STOP=i_1)
 
      add_sol = (S0-0.5D0*(S0min+S0max))/(S0max-S0min)
      IF ( AM_I_ROOT() ) THEN
         WRITE (6,661) JJDAYO, S0, S0min, S0max, add_sol
 
 661     FORMAT ('JJDAYO,S0,S0min,S0max,frac=',I4,3F9.2,F7.3)
      ENDIF
 
      ALLOCATE (delta_o3_now(grid%I_STRT_HALO:grid%I_STOP_HALO,         &
              & grid%J_STRT_HALO:grid%J_STOP_HALO,nlo3))
 
      CALL READ_STREAM(grid,delta_O3stream,jyearx,jjdayo,delta_o3_now)
      DO j = j_0, j_1
         DO i = i_0, i_1
            O3JDAY(:,I,J) = O3JDAY(:,I,J) + add_sol*delta_O3_now(i,j,:)
         ENDDO
      ENDDO
 
      DEALLOCATE (delta_o3_now)
      END SUBROUTINE UPDO3D_SOLAR
 
      END MODULE O3MOD

      SUBROUTINE SET_FPXCO2(PL,FPXCO2,NL,KFPCO2)
      USE FILEMANAGER, ONLY:FILE_EXISTS, OPENUNIT, CLOSEUNIT
      USE DICTIONARY_MOD, ONLY:SYNC_PARAM, SET_PARAM
      IMPLICIT NONE
      INTEGER J, N, NL, iu, np, NCOL
      REAL*8 PL(NL), FPXCO2(NL)
      INTEGER, PARAMETER  ::  NCOLS = 4
      REAL*8 FPI, FPJ, PFI, PFJ, pf(NCOLS)
      REAL*8, ALLOCATABLE  ::  FPX(:), PFP(:)
      CHARACTER*80 title
!@dbparam KFPCO2 selects CO2 profile absorber scaling (if >0 )
      INTEGER  ::  KFPCO2
                        ! KFPCO2 will be set from NL or from rundeck
!
! FPXCO2 scaling factors: 1.0 for P > 50mb, linear in P for P < 50mb
! PFP Pressure scale inflection points: continuous linear line segments
! PL=layerL mean pressure, FPXCO2=CO2 absorber scaling factor
! NL=total number of radiation layers incl. the top 3 rad. only layers
!
! CO2 profile absorber scaling: KFPCO2=0  FPXCO2=1, no scaling
!                               KFPCO2=1  FPXCO2:  43-layer scaling
!                               KFPCO2=2  FPXCO2:  99-layer scaling
!                               KFPCO2=3  FPXCO2: 105-layer scaling (2017)
!                               KFPCO2=4  FPXCO2: 105-layer scaling plus
!                               adjust up-flux corr. factors in top 10 layers
!                               KFPCO2<0  NL determines scaling
!                               KFPCO2>4  FPXCO2=1, no scaling
 
      CALL SYNC_PARAM("KFPCO2",KFPCO2)
 
      FPXCO2 = 1. ! default
 
! KFPCO2>2 is reserved for a specific 102-layer modelE (used in year 2017)
! but may also work if the layering is the same above 50 mb; the criterion
! below only checks the number of layers above 50 mb - it may be necessary
! to specify KFPCO2 rather than rely on the automatic selection
      IF ( NL>40 ) THEN
         IF ( KFPCO2<0 .AND. ABS(PL(NL-38)-45.)<5. ) THEN
            KFPCO2 = 4
            CALL SET_PARAM("KFPCO2",KFPCO2,'o')
         ENDIF
      ENDIF
      IF ( KFPCO2>2 .OR. KFPCO2==0 ) RETURN
 
      IF ( .NOT.FILE_EXISTS('CO2profile') ) THEN
         KFPCO2 = 0
         CALL SET_PARAM("KFPCO2",KFPCO2,'o')
         RETURN
      ENDIF
 
      CALL OPENUNIT('CO2profile',iu,.FALSE.,.TRUE.)
 
      READ (iu,'(a)') title
      READ (title,*) np
      READ (iu,'(a)') title
      READ (iu,'(a)') title
      READ (iu,'(a)') title
 
! Find appropriate column for current layering
      IF ( KFPCO2<0 ) THEN
         IF ( nl<30 ) THEN
            KFPCO2 = 0
         ELSEIF ( nl<80 ) THEN
            KFPCO2 = 1
         ELSE
            KFPCO2 = 2
         ENDIF
         CALL SET_PARAM("KFPCO2",KFPCO2,'o')
      ENDIF
 
      IF ( KFPCO2>2 .OR. KFPCO2<1 ) RETURN
 
      ncol = 2*KFPCO2 - 1
 
      ALLOCATE (FPX(np),PFP(np))
      DO n = 1, np
         READ (iu,*) pf
         pfp(n) = pf(ncol)
         fpx(n) = pf(ncol+1)
      ENDDO
 
      CALL CLOSEUNIT(iu)
 
! FPX CO2 scaling profile: (1.0 for P > 50mb) (linear in P for P < 50mb)
      j = 1
      FPj = FPX(j)
      PFj = PFP(j)
      N = 1
      DO WHILE ( PL(N)>=PFj )
         FPXCO2(N) = FPj
         N = N + 1
         IF ( N>NL ) GOTO 100
      ENDDO
 
      DO j = 2, np
         FPI = FPj
         PFI = PFj
         FPj = FPX(j)
         PFj = PFP(j)
         DO WHILE ( PL(N)>=PFj )
            FPXCO2(N) = FPI - (FPI-FPj)*(PFI-PL(N))/(PFI-PFj)
            N = N + 1
            IF ( N>NL ) GOTO 100
         ENDDO
      ENDDO
 
 100  DEALLOCATE (FPX,PFP)
      END SUBROUTINE SET_FPXCO2
 
      SUBROUTINE GET_FPXCO2_105(FPZCO2,JLAT,MLAT46,JDAY)
      IMPLICIT NONE
      INTEGER JLAT, MLAT46, JDAY
 
      INTENT (IN)JLAT, MLAT46, JDAY
      INTENT (OUT)FPZCO2                     ! FPZCO2 <==> FPXCO2
 
      REAL*8 FPZCO2(39)
      REAL*8, DIMENSION(39)  ::  FPZ_JAN, FPZ_JUL
      REAL*8 REFLAT, WTJLAT, REFDAY, WTJDAY, WT1, WT2, WT3
      REAL*8, PARAMETER  ::  FPX_SPEQNP_JAN(39,3)                       &
                           & = RESHAPE((/0.106139D+01,0.106783D+01,     &
                           & 0.106256D+01,0.106812D+01,0.106650D+01,    &
                           & 0.105212D+01,0.100292D+01,0.978585D+00,    &
                           & 0.973002D+00,0.104304D+01,0.103306D+01,    &
                           & 0.969076D+00,0.958428D+00,0.101544D+01,    &
                           & 0.100945D+01,0.993654D+00,0.991752D+00,    &
                           & 0.979623D+00,0.941461D+00,0.937047D+00,    &
                           & 0.928036D+00,0.909102D+00,0.900246D+00,    &
                           & 0.901978D+00,0.874685D+00,0.890840D+00,    &
                           & 0.947327D+00,0.980658D+00,0.101200D+01,    &
                           & 0.102130D+01,0.727029D+00,0.750237D+00,    &
                           & 0.856221D+00,0.852355D+00,0.991601D+00,    &
                           & 0.100000D+01,0.100000D+01,0.100000D+01,    &
                           & 0.100000D+01,0.795649D+00,0.704146D+00,    &
                           & 0.727177D+00,0.792645D+00,0.748645D+00,    &
                           & 0.872383D+00,0.778125D+00,0.761280D+00,    &
                           & 0.769040D+00,0.762770D+00,0.787812D+00,    &
                           & 0.795996D+00,0.812415D+00,0.806050D+00,    &
                           & 0.931991D+00,0.879633D+00,0.936397D+00,    &
                           & 0.942009D+00,0.914533D+00,0.932472D+00,    &
                           & 0.905136D+00,0.880848D+00,0.810632D+00,    &
                           & 0.878675D+00,0.901477D+00,0.951314D+00,    &
                           & 0.101205D+01,0.109697D+01,0.105704D+01,    &
                           & 0.120352D+01,0.128532D+01,0.152551D+01,    &
                           & 0.107507D+01,0.993298D+00,0.993298D+00,    &
                           & 0.100000D+01,0.100000D+01,0.100000D+01,    &
                           & 0.100000D+01,0.105912D+01,0.104105D+01,    &
                           & 0.103011D+01,0.101845D+01,0.991954D+00,    &
                           & 0.980846D+00,0.948530D+00,0.898252D+00,    &
                           & 0.900887D+00,0.931042D+00,0.913461D+00,    &
                           & 0.885751D+00,0.882327D+00,0.962325D+00,    &
                           & 0.978728D+00,0.984428D+00,0.997767D+00,    &
                           & 0.974346D+00,0.962522D+00,0.951995D+00,    &
                           & 0.935819D+00,0.925180D+00,0.923809D+00,    &
                           & 0.912133D+00,0.915321D+00,0.915390D+00,    &
                           & 0.905317D+00,0.931876D+00,0.976460D+00,    &
                           & 0.974906D+00,0.100000D+01,0.100000D+01,    &
                           & 0.100000D+01,0.100000D+01,0.100000D+01,    &
                           & 0.100000D+01,0.100000D+01,0.100000D+01,    &
                           & 0.100000D+01/),(/39,3/))
      REAL*8, PARAMETER  ::  FPX_SPEQNP_JUL(39,3)                       &
                           & = RESHAPE((/0.101719D+01,0.988671D+00,     &
                           & 0.986317D+00,0.992046D+00,0.979385D+00,    &
                           & 0.962768D+00,0.917461D+00,0.884814D+00,    &
                           & 0.904524D+00,0.954613D+00,0.967657D+00,    &
                           & 0.923239D+00,0.926085D+00,0.942071D+00,    &
                           & 0.948953D+00,0.923623D+00,0.947999D+00,    &
                           & 0.894875D+00,0.928884D+00,0.938342D+00,    &
                           & 0.904316D+00,0.906018D+00,0.893131D+00,    &
                           & 0.876373D+00,0.869938D+00,0.849842D+00,    &
                           & 0.890331D+00,0.913173D+00,0.933819D+00,    &
                           & 0.880337D+00,0.853570D+00,0.710237D+00,    &
                           & 0.826221D+00,0.772355D+00,0.991601D+00,    &
                           & 0.100000D+01,0.100000D+01,0.100000D+01,    &
                           & 0.100000D+01,0.813710D+00,0.718946D+00,    &
                           & 0.774912D+00,0.764653D+00,0.790177D+00,    &
                           & 0.781660D+00,0.768685D+00,0.778010D+00,    &
                           & 0.742164D+00,0.823045D+00,0.810695D+00,    &
                           & 0.809154D+00,0.839186D+00,0.890938D+00,    &
                           & 0.912096D+00,0.957669D+00,0.940655D+00,    &
                           & 0.970092D+00,0.949381D+00,0.925820D+00,    &
                           & 0.904259D+00,0.919764D+00,0.842020D+00,    &
                           & 0.892613D+00,0.930514D+00,0.978044D+00,    &
                           & 0.972728D+00,0.108045D+01,0.115574D+01,    &
                           & 0.121910D+01,0.135634D+01,0.165630D+01,    &
                           & 0.109376D+01,0.991629D+00,0.991629D+00,    &
                           & 0.100000D+01,0.100000D+01,0.100000D+01,    &
                           & 0.100000D+01,0.109553D+01,0.103307D+01,    &
                           & 0.101880D+01,0.104636D+01,0.105279D+01,    &
                           & 0.103332D+01,0.974511D+00,0.948682D+00,    &
                           & 0.946840D+00,0.100858D+01,0.997002D+00,    &
                           & 0.924039D+00,0.896250D+00,0.963271D+00,    &
                           & 0.977617D+00,0.996577D+00,0.992301D+00,    &
                           & 0.974108D+00,0.947932D+00,0.929054D+00,    &
                           & 0.928114D+00,0.912474D+00,0.912335D+00,    &
                           & 0.915498D+00,0.901183D+00,0.925813D+00,    &
                           & 0.971726D+00,0.104123D+01,0.110262D+01,    &
                           & 0.120227D+01,0.100000D+01,0.100000D+01,    &
                           & 0.100000D+01,0.100000D+01,0.100000D+01,    &
                           & 0.100000D+01,0.100000D+01,0.100000D+01,    &
                           & 0.100000D+01/),(/39,3/))
!
! FPXCO2 scaling factors: 1.D0 for P>50mb, for layers 1-66 for NL=105
! FPXCO2 scaling is applied only to the topmost 39 layers (67-105)
! PLZ=layerN layer-mean pressure, FPZCO2(N)=CO2 absorber scaling factor
! PLZ= input GCM variable PL(N),  FPZCO2= output GCM variable FPXCO2(N)
! NL =total number of radiation layers (includes top 3 rad-only layers)
!
!                           #### operation details ####
!----------------------------------------------------------------------
!x    IF(KFPCO2.GE.3) CALL GET_FPXCO2_105(PL,FPXCO2,JLAT,MLAT46,JDAY)
!x
!x    above CALL should be placed in RADIA inside the C**** MAIN J LOOP
!x    CO2 absorber scaling in top 39 layers is JLAT and JDAY dependent.
!x    Applicable for NL=105, scaling is invoked by KFPCO2=3 or KFPCO2=4
!x
!x    KFPCO2=4 invokes additional cooling rate control in top 10 layers
!x    via CALL GET_DXTRU3_CORR in TAUGAS to adjust XTRU(96:105,3) coeff
!x
!x    TAPER option is being utilized in TAUGAS
!x        xtru(l,2:nrcf+1) = wt_one*1d0 + (1d0-wt_one)*xtru(l,2:nrcf+1)
!x        xtrd(l,2:nrcf+1) = wt_one*1d0 + (1d0-wt_one)*xtrd(l,2:nrcf+1)
!x    with: xtrd(l,2:nrcf+1) also now included (small smoothing effect)
!----------------------------------------------------------------------
 
      FPZCO2 = 1.D0
 
      REFLAT = MLAT46/2
      WTJLAT = JLAT/REFLAT
      REFDAY = 183.D0
      WTJDAY = JDAY/REFDAY
 
      IF ( WTJLAT<1.D0 ) THEN
         WT1 = 1.D0 - WTJLAT
         IF ( WT1>0.9D0 ) WT1 = 1.D0
         WT2 = 1.D0 - WT1
         FPZ_JAN(:) = FPX_SPEQNP_JAN(:,1)*WT1 + FPX_SPEQNP_JAN(:,2)*WT2
         FPZ_JUL(:) = FPX_SPEQNP_JUL(:,1)*WT1 + FPX_SPEQNP_JUL(:,2)*WT2
      ELSE
         WT2 = 2.D0 - WTJLAT
         IF ( WT2<0.1D0 ) WT2 = 0.D0
         WT3 = 1.D0 - WT2
         FPZ_JAN(:) = FPX_SPEQNP_JAN(:,2)*WT2 + FPX_SPEQNP_JAN(:,3)*WT3
         FPZ_JUL(:) = FPX_SPEQNP_JUL(:,2)*WT2 + FPX_SPEQNP_JUL(:,3)*WT3
      ENDIF
 
      WT1 = ABS(1.D0-WTJDAY)
      WT2 = 1.D0 - WT1
      FPZCO2(:) = FPZ_JAN(:)*WT1 + FPZ_JUL(:)*WT2
 
      END SUBROUTINE GET_FPXCO2_105
 
      SUBROUTINE GET_DXTRU3_CORR(DXTRU3_10,JLAT,MLAT46,JDAY)
      IMPLICIT NONE
      INTEGER JLAT, MLAT46, JDAY
 
      INTENT (IN)JLAT, MLAT46, JDAY
      INTENT (OUT)DXTRU3_10
 
      REAL*8 DXTRU3_10(10)
      REAL*8, DIMENSION(10)  ::  DX3_JAN(10), DX3_JUL(10)
      REAL*8 REFLAT, WTJLAT, REFDAY, WTJDAY, WT1, WT2, WT3
 
      REAL*8, PARAMETER  ::  DXTRU3_SPEQNP_JAN(10,3)                    &
                           & = RESHAPE((/0.000000D+00,0.111138D-03,     &
                           & 0.594676D-04,0.770460D-04,0.694530D-04,    &
                           & 0.645742D-04,0.245626D-04,0.248727D-04,    &
                           & -.520744D-05,0.165988D-04,0.000000D+00,    &
                           & -.891406D-05,-.743531D-04,-.160635D-04,    &
                           & -.119857D-04,0.873209D-05,0.473174D-05,    &
                           & 0.158620D-04,0.548129D-05,0.121151D-04,    &
                           & 0.000000D+00,-.176254D-04,-.395480D-04,    &
                           & 0.200113D-04,0.172199D-04,-.661705D-05,    &
                           & 0.399818D-04,0.130381D-04,0.198000D-04,    &
                           & 0.263105D-04/),(/10,3/))
      REAL*8, PARAMETER  ::  DXTRU3_SPEQNP_JUL(10,3)                    &
                           & = RESHAPE((/0.000000D+00,0.120614D-04,     &
                           & 0.809026D-05,0.473324D-05,0.482039D-05,    &
                           & -.863633D-05,-.856459D-06,0.335398D-04,    &
                           & 0.422371D-04,0.512790D-04,0.000000D+00,    &
                           & 0.128452D-04,0.811048D-05,0.153064D-04,    &
                           & 0.799234D-05,0.214280D-04,0.146487D-04,    &
                           & 0.225772D-05,0.327750D-05,0.835161D-05,    &
                           & 0.000000D-04,0.652637D-04,0.651595D-04,    &
                           & 0.637362D-04,0.538371D-04,0.429989D-04,    &
                           & 0.948327D-05,0.194748D-04,-.200184D-06,    &
                           & 0.110353D-04/),(/10,3/))
 
! XTRU(L,3)= CO2 LW up-flux correction factor: layers 96-105 for NL=105
! DXTRU3_SPEQNP(L,1)= SP region, DXTRU3(L,2)= EQ region, DXTRU3(L,3)=NP
! MLAT46=total number of latitude points, interpolation utilizes JLAT
! JDAY interpolation in time: (JDAY=1 =>JAN data) (JDAY=183 =>JUL data)
!
!                           #### operation details ####
!----------------------------------------------------------------------
!x     IF(KFPCO2.EQ.4) THEN
!x     CALL GET_DXTRU3_CORR(DXTRU3_10,JLAT,MLAT46,JDAY)
!x     XTRU(96:105,3)=1.D0+DXTRU3_10
!x     ENDIF
!x     above CALL sequence should appear just before RETURN from TAUGAS
!----------------------------------------------------------------------
 
      REFLAT = MLAT46/2
      WTJLAT = JLAT/REFLAT
      REFDAY = 183.D0
      WTJDAY = JDAY/REFDAY
 
      IF ( WTJLAT<1.D0 ) THEN
         WT1 = 1.D0 - WTJLAT
         IF ( WT1>0.9D0 ) WT1 = 1.D0
         WT2 = 1.D0 - WT1
         DX3_JAN(:) = DXTRU3_SPEQNP_JAN(:,1)                            &
                    & *WT1 + DXTRU3_SPEQNP_JAN(:,2)*WT2
         DX3_JUL(:) = DXTRU3_SPEQNP_JUL(:,1)                            &
                    & *WT1 + DXTRU3_SPEQNP_JUL(:,2)*WT2
      ELSE
         WT2 = 2.D0 - WTJLAT
         IF ( WT2<0.1D0 ) WT2 = 0.D0
         WT3 = 1.D0 - WT2
         DX3_JAN(:) = DXTRU3_SPEQNP_JAN(:,2)                            &
                    & *WT2 + DXTRU3_SPEQNP_JAN(:,3)*WT3
         DX3_JUL(:) = DXTRU3_SPEQNP_JUL(:,2)                            &
                    & *WT2 + DXTRU3_SPEQNP_JUL(:,3)*WT3
      ENDIF
 
      WT1 = ABS(1.D0-WTJDAY)
      WT2 = 1.D0 - WT1
      DXTRU3_10(:) = DX3_JAN(:)*WT1 + DX3_JUL(:)*WT2
 
      END SUBROUTINE GET_DXTRU3_CORR						
