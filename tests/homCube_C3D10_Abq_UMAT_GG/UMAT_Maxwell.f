CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C UMAT FOR VISCOELASTIC VISCOPLASTIC DAMAGE CONSTITUTIVE MODEL         C
C                                                                      C
C CREATED BY JAKOB SCHWIEDRZIK - ISTB - 2012                           C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C                
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C       FEATURES                                                       C
C                                                                      C
C       * ANISOTROPIC ELASTICITY AND VISCOELASTICITY                   C
C       * QUADRIC YIELD SURFACE AND PLASTIC POTENTIAL                  C
C       * SCALAR DAMAGE                                                C
C       * VISCOPLASTICITY AND ISOTROPIC HARDENING                      C
C       * LINE SEARCH AND PRIMAL CPPA FOR GLOBAL CONVERGENCE           C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                                                                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C       MODIFICATIONS                                                  C
C                                                                      C
C       GABRIELA GERBER - MSB - 2024                                   C
C                                                                      C
C       * ADDED VISCOELASTICITY (MAXWELL LAYERS)                       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  * D E P V A R                                                       C
C                                                                      C
C    SDV 1:       CUMULATED PLASTIC STRAIN IN ELASTOPLASTIC LAYER      C
C    SDV 2-11:    CUMULATED PLASTIC STRAIN IN MAXWELL LAYERS           C
C    SDV 12:      SCALAR DAMAGE VARIABLE IN ELASTOPLASTIC LAYER        C
C    SDV 13-22:   SCALAR DAMAGE VARIABLE IN MAXWELL LAYERS             C
C    SDV 23:      ORIGINAL BVTV                                        C
C    SDV 24-29:   PLASTIC STRAIN IN ELASTOPLASTIC LAYER                C
C    SDV 30-89:   PLASTIC STRAIN IN MAXWELL LAYERS                     C
C    SDV 90-149:  VISCOUS STRAIN IN MAXWELL LAYERS                     C
C    SDV 150:     STRAIN RATE                                          C
C                                                                      C
C  * U S E R   M A T E R I A L                                         C
C                                                                      C
C     PROPS(1)    VISCOSITY FLAG: 0 - NO VISCOELASTICITY               C
C                                 1 - VISCOELASTICITY                  C
C     PROPS(2)    POSTYIELD FLAG: 0 - PERFECT PLASTICITY               C
C                                 1 - LINEAR HARDENING                 C
C                                 2 - EXPONENTIAL HARDENING            C
C                                 3 - SIMPLE SOFTENING                 C
C                                 4 - PIECEWISE SOFTENING              C
C                                                                      C
C     PROPS(3-11)  STIFFNESS TENSOR SSSP                               C
C     PROPS(12-20) YIELD TENSOR FFFF                                   C
C     PROPS(21-23) YIELD TENSOR FF                                     C
C                                                                      C
C     PROPS(24-33) SCALE FACTORS FOR MAXWELL LAYERS                    C
C     PROPS(34-43) CHARACTERISTIC TIMES FOR MAXWELL LAYERS             C
C     PROPS(44)    NUMBER OF MAXWELL LAYERS (MAX 10)                   C
C                                                                      C
C     PROPS(45)   DAMAGE PARAMETER KONSTD                              C
C     PROPS(46)   DAMAGE PARAMETER CRITD                               C
C     PROPS(47)   HARDENING PARAMETER RDY                              C
C     PROPS(48)   HARDENING PARAMETER KSLOPE                           C
C     PROPS(49)   HARDENING PARAMETER KMAX                             C
C     PROPS(50)   HARDNEING PARAMETER KWIDTH                           C 

C     PROPS(51)   TOLERANCE OF NUMERICAL ERROR                         C
C     PROPS(52)   MAX NUMBER OF ITERATIONS                             C
C     PROPS(53)   LINE SEARCH PARAMETER BETA                           C
C     PROPS(54)   LINE SEARCH PARAMETER ETA                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     ABAQUS INPUT FILE TEMPLATE                                       C
C                                                                      C
C     ** XXX to be replaced by user                                    C
C     ** CardName BVTV, EigVal.., replace by medtool                   C
C     *MATERIAL, NAME=CardName                                         C
C     *DEPVAR                                                          C
C     24                                                               C
C     *USER MATERIAL, CONSTANTS=9, UNSYMM                              C
C     XXX, BVTV, EigVal3, EigVal2,  EigVal1, XXX, XXX, XXX, XXX        C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1                DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,
     2                PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,
     3                NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,
     4                NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      IMPLICIT NONE
C
C     INPUT/OUTPUT VARIABLES:
      INTEGER NTENS, NDI, NSHR
      INTEGER NSTATV
      INTEGER NPROPS
C
C     INPUT/OUTPUT VARIABLES THAT HAVE TO BE DEFINED BY THE UMAT
C     IN ALL SITUATIONS:
      DOUBLE PRECISION STRESS(NTENS), STATEV(NSTATV)
      DOUBLE PRECISION DDSDDE(NTENS,NTENS)
      DOUBLE PRECISION SSE,SPD,SCD
C
C     ONLY IN FULLY COUPLED TEMPERATURE-DISPLACEMENT ANALYSIS:
      DOUBLE PRECISION RPL,DDSDDT(NTENS),DRPLDE(NTENS),DRPLDT
C
C     OUTPUT VARIABLES THAT CAN BE UPDATED WITHIN THE UMAT (TIME INC):
      DOUBLE PRECISION PNEWDT
C
C     INPUT VARIABLES: SUPPLIED BY ABAQUS:
      DOUBLE PRECISION STRAN(NTENS), DSTRAN(NTENS)
      DOUBLE PRECISION TIME(2)     , DTIME
      DOUBLE PRECISION TEMP        , DTEMP
      DOUBLE PRECISION PREDEF(1)   , DPRED(1)
      CHARACTER*8      CMNAME
      DOUBLE PRECISION PROPS(NPROPS)
      DOUBLE PRECISION COORDS(3)
      DOUBLE PRECISION DROT(3,3)
      DOUBLE PRECISION CELENT
      DOUBLE PRECISION DFGRD0(3,3) , DFGRD1(3,3)
      INTEGER          NOEL,NPT,LAYER,KSPT,KSTEP,KINC
C
C     EXTERNAL FUNCTIONS
      DOUBLE PRECISION RADK,DRADK,DAM,DDAM,VECNORM 
C
C     VISCOELASTIC-VISCOPLASTIC DAMAGE UMAT VARIABLES
      INTEGER ITER,MAXITER,K1,K2,VELFL
      INTEGER NLAYER,I,L
      DOUBLE PRECISION TAUV(10),LAMBDAV(10),TOL
      DOUBLE PRECISION VVVI(6,6,10),SSSV(6,6,10)
      DOUBLE PRECISION CCCV(6,6,10),VVVV(6,6,10)
      DOUBLE PRECISION SSSP(6,6),CCCP(6,6),TANMP(6,6)
      DOUBLE PRECISION VVVI_MX(6,6),SSSV_MX(6,6),CCCC_MX(6,6)
      DOUBLE PRECISION SSP1(6),STR(6),SSI(6),SS1(6),DSS1(6)
      DOUBLE PRECISION SSMX(6,10),TANMMX(6,6,10),SSL(6)
      DOUBLE PRECISION SSSA_MX(6,6),DDK1_MX,DSS1_MX(6),SSI_MX(6)
      DOUBLE PRECISION ETOT1(6),ETOT0(6)
      DOUBLE PRECISION EVISC_MXI(6,10),EPLAS_MXI(6,10)
      DOUBLE PRECISION EPLAS0I(6)
      DOUBLE PRECISION EPLAS0_MX(6),EPLAS1_MX(6),EVISC0_MX(6)
      DOUBLE PRECISION EVISC_TR(6),STR_MX(6),FFS_MX(6),SFFS_MX,YMX
      DOUBLE PRECISION EPLAS_MX(6,10),EVISC_MX(6,10),KAPPA1_MX
      DOUBLE PRECISION KAPPA0,DKAPPAI,DKAPPA1,DDK1,KAPPA1
      DOUBLE PRECISION KAPPA0_MX,DMG0_MX,KAPPA_MX(10),DKAPPA1_MX
      DOUBLE PRECISION DMG,DMG0,DDMG,RAD,DRAD,DMG_MX,DDMG_MX
      DOUBLE PRECISION DKAPPAI_MX,EPLAS0(6),EPLAS1(6)
      DOUBLE PRECISION FFFF(6,6),FF(6),GGGG(6,6),GG(6)
      DOUBLE PRECISION FFFF_MX(6,6),FF_MX(6),RAD_MX,GGGG_MX(6,6)
      DOUBLE PRECISION DDSG(6,6),DSG(6)
      DOUBLE PRECISION DRAD_MX,GGS_MX(6),SGGS_MX,GG_MX(6)
      DOUBLE PRECISION DDSG_MX(6,6),DSG_MX(6),HIMX,NPMX(6)
      DOUBLE PRECISION SFFS,FFS(6),YSTR,YS,DYDS(6),DYDK,ABSY
      DOUBLE PRECISION SGGS,GGS(6),NP(6),DNPDS(6,6),HI,DHDS(6)
      DOUBLE PRECISION DHDS_MX(6),DNPDS_MX(6,6)
      DOUBLE PRECISION DYDS_MX(6),DYDK_MX,DRRDS_MX1(6,6),DRRDK_MX1(6)
      DOUBLE PRECISION RR(6),DRRDS(6,6),DRRDK(6),NORMRR
      DOUBLE PRECISION RR_MX1(6),ABSY_MX,NORMRR_MX1
      DOUBLE PRECISION EVISCI(6)
      DOUBLE PRECISION IIII(6,6)
      DOUBLE PRECISION SSSA(6,6),TANM(6,6)
      DOUBLE PRECISION RR_MX2(6),NORMRR_MX2,DRRDS_MX2(6,6)
      DOUBLE PRECISION DRRDEV_MX1(6,6)
      DOUBLE PRECISION EVISC1_MX(6),DEV(6),NONCONV
      DOUBLE PRECISION KAPPA_MX_BUFF(10), EPLAS_MX_BUFF(6,10)
      DOUBLE PRECISION EVISC_MX_BUFF(6,10), SSMX_BUFF(6,10)
C
C     PRIMAL CPP ALGORITHM AND LINE SEARCH VARIABLES
C
      INTEGER ITERL,FLAG
      DOUBLE PRECISION J11(6,6),J22,J12(6),J21(6)
      DOUBLE PRECISION RRR(7),RRI(7),XXI(7),XX1(7),DDD(7),JJJ(7,7)
      DOUBLE PRECISION INVJ(7,7),DDJ(7,7),AUX,DMK,UPLIM
      DOUBLE PRECISION ETAL,BETA,MKI,MK1,ALPHA1,ALPHA2,ALPHA
C
C     INTERFACE FOR ARRAY VALUED EXTERNAL FUNCTION VECDYAD
C
      INTERFACE
      FUNCTION VECDYAD(BB,CC)
      DOUBLE PRECISION VECDYAD(6,6)
      DOUBLE PRECISION BB(6),CC(6)
      END FUNCTION VECDYAD
      END INTERFACE
C
C     COMMON BLOCK
C
      INTEGER PYFL
      DOUBLE PRECISION KSLOPE,KMAX,KWIDTH,RDY,KMIN,GMIN,ND,EXPS
      DOUBLE PRECISION OFFS,KONSTD,CRITD
      COMMON KSLOPE,KMAX,KWIDTH,RDY,KMIN,GMIN,ND,EXPS,KONSTD
      COMMON CRITD,OFFS,PYFL,VELFL
C
C     INTERNAL FUNCTIONS USED
C
C     DOT_PRODUCT(A,B), MATMUL(A,B), DEXP(A), DSQRT(A), DLOG(A)
C     
C     TOLERANCE OF NUMERICAL ERROR  
C 
      TOL=PROPS(51)
C
C     MAXIMUM ALLOWED ITERATIONS  
C 
      MAXITER=PROPS(52)
C
C     LINE SEARCH PARAMETERS
C
      BETA=PROPS(53)
      ETAL=PROPS(54)
C
C     MATERIAL FLAGS
C      
      PYFL=PROPS(2)
      VELFL=PROPS(1)
C
C     VISCOELASTICITY PARAMETERS
C 
      NLAYER=PROPS(44)
C   
      LAMBDAV(1)=PROPS(24)
      LAMBDAV(2)=PROPS(25)
      LAMBDAV(3)=PROPS(26)
      LAMBDAV(4)=PROPS(27)
      LAMBDAV(5)=PROPS(28)
      LAMBDAV(6)=PROPS(29)
      LAMBDAV(7)=PROPS(30)
      LAMBDAV(8)=PROPS(31)
      LAMBDAV(9)=PROPS(32)
      LAMBDAV(10)=PROPS(33)
C   
      TAUV(1)=PROPS(34)
      TAUV(2)=PROPS(35)
      TAUV(3)=PROPS(36)
      TAUV(4)=PROPS(37)
      TAUV(5)=PROPS(38)
      TAUV(6)=PROPS(39)
      TAUV(7)=PROPS(40)
      TAUV(8)=PROPS(41)
      TAUV(9)=PROPS(42)
      TAUV(10)=PROPS(43)
C
C     PLASTIC AND DAMAGE PROPERTIES
C
C     HARDENING PARAMETERS
C 
      RDY=PROPS(47)
      KSLOPE=PROPS(48)
      KMAX=PROPS(49)
      KWIDTH=PROPS(50)
C
C     DAMAGE PARAMETERS (WOLFRAM J BIOMECH 2011)
C 
      KONSTD=PROPS(45)
      CRITD=PROPS(46)
C
C     FABRIC-BASED ELASTOPLASTIC STIFFNESS TENSOR SSSP
C
      SSSP=0.0D0
      SSSP(1,1)=PROPS(3)
      SSSP(2,2)=PROPS(4)
      SSSP(3,3)=PROPS(5)
      SSSP(4,4)=PROPS(6)
      SSSP(5,5)=PROPS(7)
      SSSP(6,6)=PROPS(8)
      SSSP(2,1)=PROPS(9)
      SSSP(3,1)=PROPS(10)
      SSSP(3,2)=PROPS(11)
      SSSP(1,2)=SSSP(2,1)
      SSSP(1,3)=SSSP(3,1)
      SSSP(2,3)=SSSP(3,2)
C
C     FABRIC-BASED ELASTOPLASTIC COMPLIANCE TENSOR CCCP
C     
      CALL MIGS(SSSP,6,CCCP)
C
C     FABRIC-BASED QUADRIC FOURTH ORDER TENSOR FFFF
C 
      FFFF=0.0D0
      FFFF(1,1)=PROPS(12)
      FFFF(2,2)=PROPS(13)
      FFFF(3,3)=PROPS(14)
      FFFF(4,4)=PROPS(15)
      FFFF(5,5)=PROPS(16)
      FFFF(6,6)=PROPS(17)
      FFFF(1,2)=PROPS(18)
      FFFF(1,3)=PROPS(19)
      FFFF(2,3)=PROPS(20)
      FFFF(2,1)=FFFF(1,2)
      FFFF(3,1)=FFFF(1,3)
      FFFF(3,2)=FFFF(2,3)
C 
C     FABRIC-BASED QUADRIC SECOND ORDER TENSOR FF
C
      FF=0.0D0
      FF(1)=PROPS(21)
      FF(2)=PROPS(22)
      FF(3)=PROPS(23)
      FF(4)=0.0D0
      FF(5)=0.0D0
      FF(6)=0.0D0
C
C     PLASTIC POTENTIAL TENSORS GGGG AND GG
C
      GGGG=FFFF
      GG=FF
C
C     VISCOELASTIC STIFFNESS TENSOR SSSV
C
      SSSV=0.0D0
      CCCV=0.0D0
C
      IF(VELFL.GT.0) THEN
        DO I=1,NLAYER      
            SSSV(1:6,1:6,I)=LAMBDAV(I)*SSSP
            CCCV(1:6,1:6,I)=CCCP/LAMBDAV(I)
        END DO
      ENDIF
C
C     INVERSE OF VISCOSITY TENSOR VVVI AND VVVV
C
      VVVI=0.0D0
      VVVV=0.0D0
C
      IF(VELFL.GT.0) THEN
        DO I=1,NLAYER
            VVVI(1:6,1:6,I)=CCCV(1:6,1:6,I)/TAUV(I)
            VVVV(1:6,1:6,I)=SSSV(1:6,1:6,I)*TAUV(I)
        END DO
      ENDIF
C
C     6X6 IDENTITY MATRIX
      IIII =0.0D0
      IIII(1,1)=1.0D0
      IIII(2,2)=1.0D0
      IIII(3,3)=1.0D0
      IIII(4,4)=1.0D0
      IIII(5,5)=1.0D0
      IIII(6,6)=1.0D0
      
      NONCONV=0
C     __________________________________________________________________
C
C     RECOVER STATE VARIABLES
C
C     INPUT STRAIN
C
      ETOT1=STRAN+DSTRAN
      ETOT0=STRAN
C
C     CUMULATED PLASTIC STRAIN AT THE BEGINNING OF THE CURRENT INCREMENT
C     
      KAPPA0=STATEV(1)
C
      IF (KAPPA0<0.0D0) THEN
        KAPPA0=0.0D0
      ENDIF 
C
      DO I=1,NLAYER
            KAPPA_MX(I)=STATEV(1+I)
            IF (KAPPA_MX(I)<0.0D0) THEN
                  KAPPA_MX(I)=0.0D0
            ENDIF
      END DO
C
C     RECOVER PLASTIC STRAIN IN ELASTOPLASTIC LAYER
C       
      EPLAS0I(1)=STATEV(24)
      EPLAS0I(2)=STATEV(25)
      EPLAS0I(3)=STATEV(26)
      EPLAS0I(4)=STATEV(27)
      EPLAS0I(5)=STATEV(28)
      EPLAS0I(6)=STATEV(29)
C
C     RECOVER VISCOUS AND PLASTIC STRAIN IN MAXWELL LAYERS
C 
      DO L=1,NLAYER
            DO I=1,6
                  EPLAS_MXI(I,L)=STATEV(29+(L-1)*6+I)
                  EVISC_MXI(I,L)=STATEV(89+(L-1)*6+I)
            END DO
      END DO
C
C     ROTATE RECOVERED STRAINS
C
      CALL ROTSIG(EPLAS0I,DROT,EPLAS0,2,NDI,NSHR)
C
      DO L=1,NLAYER
            CALL ROTSIG(EPLAS_MXI(:,L),DROT,EPLAS_MX(:,L),2,NDI,NSHR)
            CALL ROTSIG(EVISC_MXI(:,L),DROT,EVISC_MX(:,L),2,NDI,NSHR)
      END DO

C
C     CONVERT INITIAL STRAINS TO MANDEL NOTATION
C 
      DO K1=4,6
         ETOT1(K1)=ETOT1(K1)/DSQRT(2.0D0)
         ETOT0(K1)=ETOT0(K1)/DSQRT(2.0D0)
         EPLAS0(K1)=EPLAS0(K1)/DSQRT(2.0D0)
         DO L=1,NLAYER
            EPLAS_MX(K1,L)=EPLAS_MX(K1,L)/DSQRT(2.0D0)
            EVISC_MX(K1,L)=EVISC_MX(K1,L)/DSQRT(2.0D0)
         END DO
      END DO
C
C     __________________________________________________________________
C
C     SOLUTION OF THE PRANDTEL LAYER
C
C     PLASTIC AND DAMAGE STATE AT THE BEGINNING OF THE CURRENT INCREMENT
C
      DMG=DAM(KAPPA0)
      DMG0=DAM(KAPPA0)
      DDMG=DDAM(KAPPA0)
      RAD=RADK(KAPPA0)
      DRAD=DRADK(KAPPA0)
C
C
C     ELASTIC TRIAL STRESS
C 
      STR=(1.0D0-DMG)*MATMUL(SSSP,ETOT1-EPLAS0)
C
C     YIELD CRITERION WITH TRIAL STRESS
C
      FFS=MATMUL(FFFF,STR)
      SFFS=DABS(DOT_PRODUCT(STR,FFS))
C 
      YSTR=DSQRT(SFFS)+DOT_PRODUCT(FF,STR)-RAD
C     __________________________________________________________________
C
C     ELASTIC CASE
C     __________________________________________________________________
C 
      IF ((YSTR.LE.TOL)) THEN    
C
C       UPDATE PLASTIC STATE VARIABLES
C    
        KAPPA1=KAPPA0
        EPLAS1=EPLAS0 
C
C       PRANDTEL LAYER STRESS UPDATE
C 
        SSP1 = STR
C
C       ELASTIC JACOBIAN
C
        TANMP=(1.0D0-DMG)*SSSP
C
C     __________________________________________________________________
C
C     INELASTIC CASE
C     __________________________________________________________________
C
C
      ELSE IF (YSTR>TOL) THEN
C
C       REINITIALIZE ELASTIC TRIAL STATE
C    
        ITER=0
        EPLAS1=EPLAS0
        DKAPPA1=0.0D0
        KAPPA1=KAPPA0+DKAPPA1
        SSP1=STR
C
        DMG=DAM(KAPPA1)
        DDMG=DDAM(KAPPA1)
        RAD=RADK(KAPPA1)
        DRAD=DRADK(KAPPA1)
C 
        GGS=MATMUL(GGGG,SSP1)
        SGGS=DOT_PRODUCT(SSP1,GGS)
        DSG=1.0D0/DSQRT(SGGS)*GGS+GG
        DDSG=-VECDYAD(GGS,GGS)*(SGGS)**(-1.5D0)+(SGGS)**(-0.5D0)*GGGG
        HI=VECNORM(DSG,6)
        NP=DSG/HI
C
        FFS=MATMUL(FFFF,SSP1)
        SFFS=DOT_PRODUCT(SSP1,FFS)
C 
        YS=DSQRT(SFFS)+DOT_PRODUCT(FF,SSP1)-RAD
        RR=MATMUL(CCCP,STR-SSP1)-(DMG-DMG0)*(ETOT1-EPLAS0)
     1     -(1.0D0-DMG)*DKAPPA1*NP
C   
        ABSY=DABS(YS)
        NORMRR=VECNORM(RR,6)
C
C     __________________________________________________________________
C
C         NEWTON CLOSEST POINT PROJECTION W/O LINE SEARCH
C     __________________________________________________________________ 
C
C   
        DO WHILE ((NORMRR>TOL.OR.ABSY>TOL).AND.ITER<=10)
C
          ITER=ITER+1 
C 
          SSI=SSP1
          DKAPPAI=DKAPPA1
C
C         PARTIAL DERIVATIVES FOR RECURSIVE FORMULAE
C 
          DHDS=1.0D0/DSQRT(DOT_PRODUCT(DSG,DSG))*MATMUL(DDSG,DSG)
          DNPDS=(DDSG*HI-VECDYAD(DSG,DHDS))/HI**2
C
          DYDS=1.0D0/DSQRT(SFFS)*FFS+FF
          DYDK=-DRAD
          DRRDS=-CCCP-(1.0D0-DMG)*DKAPPAI*DNPDS
          DRRDK=-DDMG*(ETOT1-EPLAS0-DKAPPAI*NP)-
     1           (1.0D0-DMG)*NP
C
C         INVERSION OF JACOBIAN
C
          CALL MIGS(-DRRDS,6,SSSA)
C
C         RECURSIVE FORMULAE
C 
          DDK1=-(YS+DOT_PRODUCT(DYDS,MATMUL(SSSA,RR)))
     1                      /(DOT_PRODUCT(DYDS,MATMUL(SSSA,DRRDK))+DYDK)
          DSS1=MATMUL(SSSA,(RR+DRRDK*DDK1))
C
C         VARIABLE UPDATE
C 
          SSP1=SSI+DSS1
C 
          DKAPPA1=DKAPPAI+DDK1         
C 
C         IF DKAPPA1 BECOMES SMALLER THAN 0 AND THEREFORE INADMISSIBLE,
C         RESTART THE NEWTON SCHEME WITH A DIFFERENT STARTING POINT
C 
          IF (DKAPPA1<0.0D0) THEN
            SSP1=STR
            DKAPPA1=1.D-6
          END IF
C 
          KAPPA1=KAPPA0+DKAPPA1
C
C         CHECK RESIDUAL AND YIELD CRITERION FOR CONVERGENCE
C 
          DMG=DAM(KAPPA1)
          DDMG=DDAM(KAPPA1)
          RAD=RADK(KAPPA1)
          DRAD=DRADK(KAPPA1)
C 
          GGS=MATMUL(GGGG,SSP1)
          SGGS=DOT_PRODUCT(SSP1,GGS)
          DSG=1.0D0/DSQRT(SGGS)*GGS+GG
          DDSG=-VECDYAD(GGS,GGS)*(SGGS)**(-1.5D0)+(SGGS)**(-0.5D0)*GGGG
          HI=VECNORM(DSG,6)
          NP=DSG/HI
C                  
          FFS=MATMUL(FFFF,SSP1)
          SFFS=DOT_PRODUCT(SSP1,FFS)
C 
          YS=DSQRT(SFFS)+DOT_PRODUCT(FF,SSP1)-RAD
          RR=MATMUL(CCCP,STR-SSP1)-(DMG-DMG0)*(ETOT1-EPLAS0)
     1        -(1.0D0-DMG)*DKAPPA1*NP
C   
          ABSY=DABS(YS)
          NORMRR=VECNORM(RR,6)
C 
        END DO 
C     __________________________________________________________________
C
C         PRIMAL CPPA FOR GLOBAL CONVERGENCE
C     __________________________________________________________________ 
C
C       ONLY ACTIVATED IN CASE OF NONCONVERGENCE OF NEWTON CPPA
C
        IF (ITER.GT.10.AND.(NORMRR.GT.TOL.OR.ABSY.GT.TOL)) THEN 
C
C         REINITIALIZATION OF ELASTIC TRIAL STATE
C 
          ITER=0
          EPLAS1=EPLAS0
          DKAPPA1=0.0D0
          KAPPA1=KAPPA0+DKAPPA1
          SSP1=STR
C
          DMG=DAM(KAPPA0)
          DDMG=DDAM(KAPPA0)
          RAD=RADK(KAPPA0)
          DRAD=DRADK(KAPPA0)
C 
          GGS=MATMUL(GGGG,SSP1)
          SGGS=DOT_PRODUCT(SSP1,GGS)
          DSG=1.0D0/DSQRT(SGGS)*GGS+GG
          DDSG=-VECDYAD(GGS,GGS)*(SGGS)**(-1.5D0)+(SGGS)**(-0.5D0)*GGGG
          HI=VECNORM(DSG,6)
          NP=DSG/HI
C                
          FFS=MATMUL(FFFF,SSP1)
          SFFS=DOT_PRODUCT(SSP1,FFS)
C 
          YS=DSQRT(SFFS)+DOT_PRODUCT(FF,SSP1)-RAD
          RR=MATMUL(CCCP,STR-SSP1)-(DMG-DMG0)*(ETOT1-EPLAS0)
     1        -(1.0D0-DMG)*DKAPPA1*NP
C   
          ABSY=DABS(YS)
          NORMRR=VECNORM(RR,6)
C
          DO K1 = 1,6
            RRI(K1) = RR(K1)
          ENDDO
            RRI(7)=YS
C
          DO WHILE ((NORMRR>TOL.OR.ABSY>TOL).AND.ITER<=MAXITER)
C 
            ITER=ITER+1 
C 
            SSI=SSP1
            DKAPPAI=DKAPPA1
C
C           PARTIAL DERIVATIVES FOR SYSTEM JACOBIAN
C 
            DHDS=1.0D0/DSQRT(DOT_PRODUCT(DSG,DSG))*MATMUL(DDSG,DSG)
            DNPDS=(DDSG*HI-VECDYAD(DSG,DHDS))/HI**2
C
            DYDS=1.0D0/DSQRT(SFFS)*FFS+FF
            DYDK=-DRAD
            DRRDS=-CCCP-(1.0D0-DMG)*DKAPPAI*DNPDS
            DRRDK=-DDMG*(ETOT1-EPLAS0-DKAPPAI*NP)-
     1             (1.0D0-DMG)*NP
C 
C           BUILD SYSTEM RESIDUAL VECTOR
C
            DO K1 = 1,6
              RRR(K1) = RR(K1)
            ENDDO
            RRR(7)=YS
C
            J11=DRRDS
            J22=DYDK
            J12=DRRDK
            J21=DYDS
C
C           BUILD SYSTEM JACOBIAN
C
            DO K1 = 1,6
              DO K2 = 1,6
                JJJ(K1,K2) = J11(K1,K2)
              ENDDO
              JJJ(7,K1)=J21(K1)
              JJJ(K1,7)=J12(K1)
            ENDDO
            JJJ(7,7)=J22
C
C           INVERT SYSTEM JACOBIAN
C
            CALL MIGS(JJJ,7,INVJ)
C
            AUX=0.0D0
            IF (DKAPPA1.EQ.0.0D0) THEN
              AUX=DOT_PRODUCT(NP,(SSP1-STR))
            END IF
C
C           DETERMINE SYSTEM UPDATES
C
            IF (AUX.LE.0.0D0) THEN
C
              FLAG=1
C
              DDD=MATMUL(-INVJ,RRR)
C
            ELSE
C
              FLAG=0
C
              DDJ=MATMUL(INVJ,TRANSPOSE(INVJ))
              DO K1 = 1,6
                DDJ(7,K1)=0.0D0
                DDJ(K1,7)=0.0D0
              ENDDO
C
              DDD=MATMUL(-MATMUL(DDJ,TRANSPOSE(JJJ)),RRR)
C
            END IF  
C
C           RECOVER STRESS AND KAPPA UPDATE
C
            DO K1 = 1,6
              DSS1(K1)=DDD(K1)
            ENDDO
            DDK1=DDD(7)            
C 
            SSP1=SSI+DSS1
            DKAPPA1=DKAPPAI+DDK1
C
            DO K1 = 1,6
              XXI(K1)=SSP1(K1)
            ENDDO
            XXI(7)=DKAPPA1
C 
C           ENFORCE CONSTRAIN DKAPPA >= 0.0D0
C
            IF (DKAPPA1<0.0D0) THEN
              DKAPPA1=0.0D0
            END IF
C 
            KAPPA1=KAPPA0+DKAPPA1
C
C           CHECK RESIDUAL AND YIELD CRITERION FOR CONVERGENCE
C 
            DMG=DAM(KAPPA1)
            DDMG=DDAM(KAPPA1)
            RAD=RADK(KAPPA1)
            DRAD=DRADK(KAPPA1)

            GGS=MATMUL(GGGG,SSP1)
            SGGS=DOT_PRODUCT(SSP1,GGS)
            DSG=1.0D0/DSQRT(SGGS)*GGS+GG
            DDSG=-VECDYAD(GGS,GGS)*(SGGS)**(-1.5D0)
     1           +(SGGS)**(-0.5D0)*GGGG
            HI=VECNORM(DSG,6)
            NP=DSG/HI
C
            FFS=MATMUL(FFFF,SSP1)
            SFFS=DOT_PRODUCT(SSP1,FFS)
C 
C           LINE SEARCH MERIT FUNCTION AND DERIVATIVE OF ITERATION I
C
            MKI=0.5D0*(DOT_PRODUCT(RR,RR)+YS**2)
C
            IF (FLAG.EQ.1) THEN
              DMK=-2.0D0*MKI
            ELSE
              DMK=DOT_PRODUCT(RRR,MATMUL(JJJ,DDD))
            END IF
C           
C           RESIDUAL AND LINE SEARCH MERIT FUNCTION OF ITERATION I+1
C     
            YS=DSQRT(SFFS)+DOT_PRODUCT(FF,SSP1)-RAD
            RR=MATMUL(CCCP,STR-SSP1)-(DMG-DMG0)*(ETOT1-EPLAS0)
     1            -(1.0D0-DMG)*DKAPPA1*NP
C
            DO K1 = 1,6
              RRR(K1) = RR(K1)
            ENDDO
            RRR(7)=YS
C
            MK1=0.5D0*(DOT_PRODUCT(RR,RR)+YS**2)
C
C           LINE SEARCH ALGORITHM FOR CONSTRAINED PROBLEMS (ARMERO 2002)
C           ACTIVATED IN CASE OF NONCONVERGENCE (MERIT FUNCTION > UPLIM)
C
            IF (FLAG.EQ.1.AND.DKAPPA1.GE.0.0D0) THEN
              UPLIM=(1.0D0-2.0D0*BETA*ALPHA)*MKI
            ELSE
              UPLIM=MKI+BETA*DOT_PRODUCT(RRR,MATMUL(JJJ,XX1-XXI))
            END IF
C
            IF (MK1.GT.UPLIM) THEN
C
              ITERL=0
              ALPHA=1.0D0
              XX1=XXI
C
              DO WHILE (MK1.GT.UPLIM.AND.ALPHA.GE.BETA.AND.ITERL.LT.20)
C
                ITERL=ITERL+1
C
                ALPHA1=ETAL*ALPHA
                ALPHA2=-ALPHA**2*DMK/2.0D0/(MK1-MKI-ALPHA*DMK)
C
                IF (ALPHA1.GE.ALPHA2) THEN
                  ALPHA=ALPHA1
                ELSE
                  ALPHA=ALPHA2
                END IF
C
                DKAPPA1=DKAPPAI+ALPHA*DDK1
C
                IF (DKAPPA1<0.0D0) THEN
                  DKAPPA1=0.0D0
                END IF
C
C               PULL BACK STATE VARIABLE UPDATE ACCORDING TO LS CRITERIA
C
                KAPPA1=KAPPA0+DKAPPA1
                SSP1=SSI+ALPHA*DSS1
C
                DO K1 = 1,6
                  XX1(K1)=SSP1(K1)
                ENDDO
                XX1(7)=DKAPPA1
C
C               CONVERGENCE CHECK IN LINE SEARCH
C
                DMG=DAM(KAPPA1)
                RAD=RADK(KAPPA1)

                GGS=MATMUL(GGGG,SSP1)
                SGGS=DOT_PRODUCT(SSP1,GGS)
                DSG=1.0D0/DSQRT(SGGS)*GGS+GG
                DDSG=-VECDYAD(GGS,GGS)*(SGGS)**(-1.5D0)
     1                +(SGGS)**(-0.5D0)*GGGG
                HI=VECNORM(DSG,6)
                NP=DSG/HI 
C                
                FFS=MATMUL(FFFF,SSP1)
                SFFS=DOT_PRODUCT(SSP1,FFS)
C   
                YS=DSQRT(SFFS)+DOT_PRODUCT(FF,SSP1)-RAD
C
                RR=MATMUL(CCCP,STR-SSP1)-(DMG-DMG0)*(ETOT1-EPLAS0)
     1             -(1.0D0-DMG)*DKAPPA1*NP
C
                DO K1 = 1,6
                  RRR(K1) = RR(K1)
                ENDDO
                RRR(7)=YS
C
                MK1=0.5D0*(DOT_PRODUCT(RR,RR)+YS**2)
C
                IF (FLAG.EQ.1.AND.DKAPPA1.GE.0.0D0) THEN
                  UPLIM=(1.0D0-2.0D0*BETA*ALPHA)*MKI
                ELSE
                  UPLIM=MKI+BETA*DOT_PRODUCT(RRI,MATMUL(JJJ,XX1-XXI))
                END IF
C       
              END DO
C
C             CONVERGENCE CHECK OUT OF LINE SEARCH
C
              DMG=DAM(KAPPA1)
              DDMG=DDAM(KAPPA1)
              RAD=RADK(KAPPA1)
              DRAD=DRADK(KAPPA1)
C 
              GGS=MATMUL(GGGG,SSP1)
              SGGS=DOT_PRODUCT(SSP1,GGS)
              DSG=1.0D0/DSQRT(SGGS)*GGS+GG
              DDSG=-VECDYAD(GGS,GGS)*(SGGS)**(-1.5D0)
     1                       +(SGGS)**(-0.5D0)*GGGG
              HI=VECNORM(DSG,6)
              NP=DSG/HI
C
              FFS=MATMUL(FFFF,SSP1)
              SFFS=DOT_PRODUCT(SSP1,FFS)
C 
              YS=DSQRT(SFFS)+DOT_PRODUCT(FF,SSP1)-RAD
C
              RR=MATMUL(CCCP,STR-SSP1)-(DMG-DMG0)*(ETOT1-EPLAS0)
     1            -(1.0D0-DMG)*DKAPPA1*NP
C
              DO K1 = 1,6
                RRI(K1) = RR(K1)
              ENDDO
              RRI(7)=YS
C
            END IF
C 
            NORMRR=VECNORM(RR,6)
            ABSY=DABS(YS)
C 
          END DO
        END IF
C
C       TANGENT STIFFNESS MATRIX
C
        DMG=DAM(KAPPA1)
        DDMG=DDAM(KAPPA1)
        RAD=RADK(KAPPA1)
        DRAD=DRADK(KAPPA1)
C 
        GGS=MATMUL(GGGG,SSP1)
        SGGS=DOT_PRODUCT(SSP1,GGS)
        DSG=1.0D0/DSQRT(SGGS)*GGS+GG
        DDSG=-VECDYAD(GGS,GGS)*(SGGS)**(-1.5D0)+(SGGS)**(-0.5D0)*GGGG
        HI=VECNORM(DSG,6)
        NP=DSG/HI
C
        FFS=MATMUL(FFFF,SSP1)
        SFFS=DOT_PRODUCT(SSP1,FFS)
C 
        YS=DSQRT(SFFS)+DOT_PRODUCT(FF,SSP1)-RAD
        RR=MATMUL(CCCP,STR-SSP1)-(DMG-DMG0)*(ETOT1-EPLAS0)
     1     -(1.0D0-DMG)*DKAPPA1*NP
C 
        DHDS=1.0D0/DSQRT(DOT_PRODUCT(DSG,DSG))*MATMUL(DDSG,DSG)
        DNPDS=(DDSG*HI-VECDYAD(DSG,DHDS))/HI**2
C
        DYDS=1.0D0/DSQRT(SFFS)*FFS+FF
        DYDK=-DRAD
        DRRDS=-CCCP-(1.0D0-DMG)*DKAPPA1*DNPDS
        DRRDK=-DDMG*(ETOT1-EPLAS0-DKAPPA1*NP)-
     1         (1.0D0-DMG)*NP 
C 
        CALL MIGS(-DRRDS,6,SSSA)
C 
        TANMP=(1.0D0-DMG)*(SSSA-MATMUL(MATMUL(SSSA,VECDYAD(DRRDK,NP)),
     1      SSSA)/(DOT_PRODUCT(NP,MATMUL(SSSA,DRRDK))+DYDK
     2      /VECNORM(DYDS,6)))
C 
C       FINAL PLASTIC STRAIN
C 
        EPLAS1=EPLAS0+DKAPPA1*NP
C 
      END IF
C
C     CHECK FOR LOCAL CONVERGENCE OF BACK PROJECTION
C 
      IF ((YSTR>=0.0D0.AND.ITER>MAXITER).OR.NONCONV>0) THEN
        WRITE(*,*) 'HIGH NUMBER OF ITERATIONS NEEDED'
        WRITE(*,*) 'Step : ',KSTEP
        WRITE(*,*) 'Inc : ',KINC
        WRITE(*,*) 'Elem : ',NOEL
        WRITE(*,*) 'Point : ',NPT
        WRITE(*,*) 'THE INCREMENT SIZE WILL BE REDUCED BY 50%'
        PNEWDT=0.5
        KAPPA1=KAPPA0
        EPLAS1=EPLAS0
        SS1=STR
        TANM=(1.0D0-DAM(KAPPA0))*SSSP
        NONCONV=1
      END IF
C
C     __________________________________________________________________
C
C     SOLUTION OF THE MAXWELL LAYER
C
C     __________________________________________________________________
C
      IF(VELFL.GT.0) THEN
C
C       BUFFER VARIABLES
        KAPPA_MX_BUFF=KAPPA_MX
        EPLAS_MX_BUFF=EPLAS_MX
        EVISC_MX_BUFF=EVISC_MX
        SSMX_BUFF=0.0D0
C
C       LOOP OVER ALL MAXWELL LAYERS
        DO L=1,NLAYER
C
C           LAYER SPECIFIC VARIABLES
            KAPPA0_MX=KAPPA_MX(L)
            DMG0_MX=DAM(KAPPA0_MX)
C
            EPLAS0_MX=EPLAS_MX(1:6,L)
            EVISC0_MX=EVISC_MX(1:6,L)
C
            VVVI_MX=VVVI(1:6,1:6,L)
            SSSV_MX=SSSV(1:6,1:6,L)
            CCCC_MX=CCCV(1:6,1:6,L)
C
            FFFF_MX = FFFF/LAMBDAV(L)/LAMBDAV(L)
            FF_MX = FF/LAMBDAV(L)
            RAD_MX = RADK(KAPPA0_MX)
            GGGG_MX = FFFF_MX
            GG_MX = FF_MX

C           DETERMINE VISCOUS TRIAL STRAIN   
            EVISC_TR=((1-DMG0_MX)*DTIME/TAUV(L)*(ETOT1-EPLAS0_MX)
     1            +EVISC0_MX)/(1+(1-DMG0_MX)*DTIME/TAUV(L))
C
C           DETERMINE TRIAL STRESS
            STR_MX = (1-DMG0_MX)*MATMUL(SSSV_MX,ETOT1-EVISC_TR-
     1            EPLAS0_MX)
            SSMX_BUFF(1:6,L)=STR_MX
C
C           TEST YIELD CRITERION
            FFS_MX=MATMUL(FFFF_MX,STR_MX)
            SFFS_MX=DABS(DOT_PRODUCT(STR_MX,FFS_MX))
C 
            YMX=DSQRT(SFFS_MX)+DOT_PRODUCT(FF_MX,STR_MX)-RAD_MX

C           ELASTIC CASE
            IF ((YMX.LE.TOL)) THEN 
C
C               UPDATE PLASTIC STATE VARIABLES    
                KAPPA_MX(L)=KAPPA0_MX
                EPLAS_MX(1:6,L)=EPLAS0_MX
C
C               FINAL VISCOUS STRAIN AND EQUIVALENT VISCOUS STRAIN
                EVISC_MX(1:6,L)=EVISC_TR
C
C               MAXWELL LAYER STRESS UPDATE 
                SSMX(1:6,L) = STR_MX
C
C               ELASTIC JACOBIAN
                TANMMX(1:6,1:6,L)=MATMUL((1.0D0-DMG0_MX)*SSSV_MX,
     1                  IIII-(1-DMG0_MX)*DTIME/(TAUV(L)+
     2                  (1.0D0-DMG0_MX)*DTIME)*IIII)
C
C           PLASTIC CASE
            ELSE IF (YMX>TOL) THEN
C
C           REINITIALIZE ELASTIC TRIAL STATE
C    
                ITER=0
                EPLAS1_MX=EPLAS0_MX
                DKAPPA1_MX=0.0D0
                KAPPA1_MX=KAPPA0_MX+DKAPPA1_MX
                SSL=STR_MX
                EVISC1_MX = EVISC0_MX
C
                DMG_MX=DAM(KAPPA1_MX)
                DDMG_MX=DDAM(KAPPA1_MX)
                RAD_MX=RADK(KAPPA1_MX)
                DRAD_MX=DRADK(KAPPA1_MX)               
C 
                GGS_MX=MATMUL(GGGG_MX,SSL)
                SGGS_MX=DOT_PRODUCT(SSL,GGS_MX)
                DSG_MX=1.0D0/DSQRT(SGGS_MX)*GGS_MX+GG_MX
                DDSG_MX=-VECDYAD(GGS_MX,GGS_MX)*(SGGS_MX)**(-1.5D0)
     1           +(SGGS_MX)**(-0.5D0)*GGGG_MX
                HIMX=VECNORM(DSG_MX,6)
                NPMX=DSG_MX/HIMX
C
                FFS_MX=MATMUL(FFFF_MX,SSL)
                SFFS_MX=DOT_PRODUCT(SSL,FFS_MX)
C
                EVISC_TR=((1-DMG_MX)*DTIME/TAUV(L)*(ETOT1-EPLAS1_MX)
     1            +EVISC0_MX)/(1+(1-DMG_MX)*DTIME/TAUV(L))
                EVISCI=EVISC_TR

                STR_MX = (1-DMG0_MX)*MATMUL(SSSV_MX,ETOT1-EVISC_TR-
     1            EPLAS0_MX)
C 
                YMX=DSQRT(SFFS_MX)+DOT_PRODUCT(FF_MX,SSL)-RAD_MX
                RR_MX1=MATMUL(CCCC_MX,STR_MX-SSL)-(DMG_MX-DMG0_MX)*
     1              (ETOT1-EVISC_TR-EPLAS0_MX)-(1.0D0-DMG_MX)
     2                  *(EVISCI-EVISC_TR+DKAPPA1_MX*NPMX)
                RR_MX2=EVISCI-MATMUL(CCCC_MX,SSL)*DTIME/TAUV(L)
     1            -EVISC0_MX
C   
                ABSY_MX=DABS(YMX)
                NORMRR_MX1=VECNORM(RR_MX1,6)
                NORMRR_MX2=VECNORM(RR_MX2,6)
C
C               __________________________________________________________________
C
C               NEWTON CLOSEST POINT PROJECTION W/O LINE SEARCH
C               __________________________________________________________________ 
C
C   
                DO WHILE ((NORMRR_MX1>TOL.OR.ABSY_MX>TOL.OR.
     1                  NORMRR_MX2>TOL).AND.ITER<=100)
C
                    ITER=ITER+1 
C 
                    SSI_MX=SSL
                    DKAPPAI_MX=DKAPPA1_MX
                    EVISCI=EVISC1_MX
C
C                   PARTIAL DERIVATIVES FOR RECURSIVE FORMULAE
C 
                    DHDS_MX=1.0D0/DSQRT(DOT_PRODUCT(DSG_MX,DSG_MX))
     1               *MATMUL(DDSG_MX,DSG_MX)
                    DNPDS_MX=(DDSG_MX*HIMX-VECDYAD(DSG_MX,DHDS_MX))/
     1                       HIMX**2
C
                    DYDS_MX=1.0D0/DSQRT(SFFS_MX)*FFS_MX+FF_MX
                    DYDK_MX=-DRAD_MX
C
                    DRRDS_MX1=-CCCC_MX-(1-DMG_MX)*DKAPPAI_MX*DNPDS_MX
                    DRRDK_MX1=-DDMG_MX*(ETOT1-EVISC_TR-EPLAS0_MX)+
     1                  DDMG_MX*(EVISCI-EVISC_TR+DKAPPAI_MX*NPMX)
     2                  -(1-DMG_MX)*NPMX
                    DRRDEV_MX1=-(1-DMG_MX)*IIII
C
                    DRRDS_MX2=-CCCC_MX*DTIME/TAUV(L)
C
C                   INVERSION OF JACOBIAN
C
                    CALL MIGS(-DRRDS_MX1+MATMUL(DRRDEV_MX1,DRRDS_MX2)
     1                  ,6,SSSA_MX)
C
C                   RECURSIVE FORMULAE
C 
                    DDK1_MX=-(YMX+DOT_PRODUCT(DYDS_MX,MATMUL(SSSA_MX,
     1                      RR_MX1-MATMUL(DRRDEV_MX1,RR_MX2))))/
     2               (DOT_PRODUCT(DYDS_MX,MATMUL(SSSA_MX,DRRDK_MX1))
     3                  +DYDK_MX)
C
                    DSS1_MX=MATMUL(SSSA_MX,(RR_MX1+DRRDK_MX1*DDK1_MX-
     1                  MATMUL(DRRDEV_MX1,RR_MX2)))
C
                    DEV=-RR_MX2-MATMUL(DRRDS_MX2,DSS1_MX)
C
C                   VARIABLE UPDATE
C 
                    SSL=SSI_MX+DSS1_MX
C 
                    DKAPPA1_MX=DKAPPAI_MX+DDK1_MX   
                    EVISC1_MX=EVISCI+DEV 
C 
C                   IF DKAPPA1 BECOMES SMALLER THAN 0 AND THEREFORE INADMISSIBLE,
C                   RESTART THE NEWTON SCHEME WITH A DIFFERENT STARTING POINT
C 
                    IF (DKAPPA1_MX<0.0D0) THEN
                        SSL=STR_MX
                        DKAPPA1_MX=1.D-6
                        EVISC1_MX=EVISC_TR
                    END IF
C 
                    KAPPA1_MX=KAPPA0_MX+DKAPPA1_MX
C
C                   RECALCULATE PLASTIC STRAIN
                    EPLAS1_MX=EPLAS0_MX+DKAPPA1_MX*NPMX
C
C                   CHECK RESIDUAL AND YIELD CRITERION FOR CONVERGENCE 
                    DMG_MX=DAM(KAPPA1_MX)
                    DDMG_MX=DDAM(KAPPA1_MX)
                    RAD_MX=RADK(KAPPA1_MX)
                    DRAD_MX=DRADK(KAPPA1_MX)
                    GGS_MX=MATMUL(GGGG_MX,SSL)
                    SGGS_MX=DOT_PRODUCT(SSL,GGS_MX)
                    DSG_MX=1.0D0/DSQRT(SGGS_MX)*GGS_MX+GG_MX
                    DDSG_MX=-VECDYAD(GGS_MX,GGS_MX)*(SGGS_MX)**(-1.5D0)+
     1               (SGGS_MX)**(-0.5D0)*GGGG_MX
                    HIMX=VECNORM(DSG_MX,6)
                    NPMX=DSG_MX/HIMX
C                  
                    FFS_MX=MATMUL(FFFF_MX,SSL)
                    SFFS_MX=DOT_PRODUCT(SSL,FFS_MX)
C 
                    YMX=DSQRT(SFFS_MX)+DOT_PRODUCT(FF_MX,SSL)-RAD_MX
C   
                    RR_MX1=MATMUL(CCCC_MX,STR_MX-SSL)-(DMG_MX-DMG0_MX)*
     1                  (ETOT1-EVISC_TR-EPLAS0_MX)-(1.0D0-DMG_MX)
     2                  *(EVISC1_MX-EVISC_TR+DKAPPA1_MX*NPMX)
                    RR_MX2=EVISC1_MX-MATMUL(CCCC_MX,SSL)*DTIME/TAUV(L)
     1                  -EVISC0_MX
C   
                    ABSY_MX=DABS(YMX)
                    NORMRR_MX1=VECNORM(RR_MX1,6)
                    NORMRR_MX2=VECNORM(RR_MX2,6)
C 
                END DO 
C
C               UPDATE PLASTIC STATE VARIABLES    
                KAPPA_MX(L)=KAPPA1_MX
C
C               MAXWELL LAYER STRESS UPDATE 
                SSMX(1:6,L) = SSL
C
C               RECALCUALTE PLASTIC STRAIN
                EPLAS1_MX=EPLAS0_MX+DKAPPA1_MX*NPMX
                EPLAS_MX(:,L)=EPLAS1_MX
C
C               UPDATE VISCOUS STRAIN
                EVISC_MX(:,L)=EVISC1_MX
C
C               PLASTIC JACOBIAN 
                DMG_MX=DAM(KAPPA1_MX)
                DDMG_MX=DDAM(KAPPA1_MX)
                RAD_MX=RADK(KAPPA1_MX)
                DRAD_MX=DRADK(KAPPA1_MX)
                GGS_MX=MATMUL(GGGG_MX,SSL)
                SGGS_MX=DOT_PRODUCT(SSL,GGS_MX)
                DSG_MX=1.0D0/DSQRT(SGGS_MX)*GGS_MX+GG_MX
                DDSG_MX=-VECDYAD(GGS_MX,GGS_MX)*(SGGS_MX)**(-1.5D0)+
     1            (SGGS_MX)**(-0.5D0)*GGGG_MX
                HIMX=VECNORM(DSG_MX,6)
                NPMX=DSG_MX/HIMX
C
                FFS_MX=MATMUL(FFFF_MX,SSL)
                SFFS_MX=DOT_PRODUCT(SSL,FFS_MX)
C 
                YMX=DSQRT(SFFS_MX)+DOT_PRODUCT(FF_MX,SSL)-RAD_MX
                RR_MX1=MATMUL(CCCC_MX,STR_MX-SSL)-(DMG_MX-DMG0_MX)*
     1                  (ETOT1-EVISC_TR-EPLAS0_MX)-(1.0D0-DMG_MX)
     2                  *(EVISC1_MX-EVISC_TR+DKAPPA1_MX*NPMX)
                RR_MX2=EVISC1_MX-MATMUL(CCCC_MX,SSL)*DTIME/TAUV(L)
     1                  -EVISC0_MX
C
                DHDS_MX=1.0D0/DSQRT(DOT_PRODUCT(DSG_MX,DSG_MX))
     1               *MATMUL(DDSG_MX,DSG_MX)
                DNPDS_MX=(DDSG_MX*HIMX-VECDYAD(DSG_MX,DHDS_MX))/HIMX**2
C
                DYDS_MX=1.0D0/DSQRT(SFFS_MX)*FFS_MX+FF_MX
                DYDK_MX=-DRAD_MX
                DRRDS_MX1=-CCCC_MX-(1-DMG_MX)*DKAPPA1_MX*DNPDS_MX
                DRRDK_MX1=-DDMG_MX*(ETOT1-EVISC_TR-EPLAS0_MX)+
     1                  DDMG_MX*(EVISC1_MX-EVISC_TR+DKAPPA1_MX*NPMX)
     2                  -(1-DMG_MX)*NPMX
                DRRDEV_MX1=-(1-DMG_MX)*IIII
                DRRDS_MX2=-CCCC_MX*DTIME/TAUV(L)
C 
                CALL MIGS(-DRRDS_MX1+MATMUL(DRRDEV_MX1,DRRDS_MX2)
     1                  ,6,SSSA_MX)
C 
                TANMMX(1:6,1:6,L)=(1.0D0-DMG_MX)*(SSSA_MX-MATMUL(MATMUL(
     1            SSSA_MX,VECDYAD(DRRDK_MX1,NPMX)),SSSA_MX)/(
     2            DOT_PRODUCT(NPMX,MATMUL(SSSA_MX,DRRDK_MX1))+DYDK_MX
     3            /VECNORM(DYDS_MX,6)))
C           CHECK FOR LOCAL CONVERGENCE OF BACK PROJECTION
C 
                 IF ((YMX>=0.0D0.AND.ITER>MAXITER).OR.NONCONV>0) THEN
                 WRITE(*,*) 'HIGH NUMBER OF ITERATIONS NEEDED'
                 WRITE(*,*) 'Step : ',KSTEP
                 WRITE(*,*) 'Inc : ',KINC
                 WRITE(*,*) 'Elem : ',NOEL
                 WRITE(*,*) 'Point : ',NPT
                 WRITE(*,*) 'THE INCREMENT SIZE WILL BE REDUCED BY 50%'
                 PNEWDT=0.5
                 KAPPA1=KAPPA0
                 EPLAS1=EPLAS0
                 SS1=STR
                 TANM=(1.0D0-DAM(KAPPA0))*SSSP
                 KAPPA_MX=KAPPA_MX_BUFF
                 EPLAS_MX=EPLAS_MX_BUFF
                 EVISC_MX=EVISC_MX_BUFF
                 SSMX=SSMX_BUFF
                 NONCONV=1
                 END IF
            END IF
C
        END DO
C
C     
      ELSE
            SSMX=0.0D0
            EVISC_MX=0.0D0
            EPLAS_MX=0.0D0
            TANMMX=0.0D0
      END IF      
C     __________________________________________________________________
C
C     COMBINATION OF THE TWO LAYERS
C
C     __________________________________________________________________
C
C     TOTAL STRESS
      SS1 = SSP1
      DO I=1,NLAYER
            SS1 = SS1 + SSMX(1:6,I)
      END DO
C
C     TOTAL JACOBIAN
      TANM = TANMP
      DO I=1,NLAYER
            TANM = TANM + TANMMX(1:6,1:6,I)
      END DO
C
C    __________________________________________________________________
C
C     RETURN VARIABLES TO ABAQUS 
C     __________________________________________________________________
C
C     VISCOELASTIC DISSIPATION (NOT CALCULATED)
C
      SCD = 0
C
C     PLASTIC AND DAMAGE DISSIPATION (NOT CALCULATED)
C     
      SPD=0 
C    
C     CONVERT TANGENT STIFFNESS OPERATOR MATRIX TO ABAQUS CONVENTION 
C 
      DO K1=4,6
        DO K2=4,6
          TANM(K1,K2)=TANM(K1,K2)/2.0D0
          TANM(K1-3,K2)=TANM(K1-3,K2)/DSQRT(2.0D0)
          TANM(K1,K2-3)=TANM(K1,K2-3)/DSQRT(2.0D0)
        END DO
      END DO
C     __________________________________________________________________
C
C     CONVERT STRESSES AND STRAINS TO ABAQUS CONVENTION 
C     
      DO K1=4,6
        SS1(K1)=SS1(K1)/DSQRT(2.0D0)
        EPLAS1(K1)=EPLAS1(K1)*DSQRT(2.0D0)
        EPLAS_MX(K1,:)=EPLAS_MX(K1,:)*DSQRT(2.0D0)
        EVISC_MX(K1,:)=EVISC_MX(K1,:)*DSQRT(2.0D0)
      END DO 
C     __________________________________________________________________     
C
C     UPDATE OF FIELD VARIABLES
C 
C     STRESS TENSOR
C 
      STRESS=SS1
C    
C     TANGENT STIFFNESS OPERATOR 
C 
      DDSDDE=TANM
C
C     STRAIN ENERGY
C 
      SSE = 0.0D0
C 
C     PLASTIC AND VISCOUS STRAIN TENSOR
      DO I=1,6
        STATEV(23+I)=EPLAS1(I)
      END DO
C
C     CUMULATED PLASTIC STRAIN AND DAMAGE IN ELASTOPLASTIC LAYER     
      STATEV(1)=KAPPA1
      STATEV(12)=DAM(KAPPA1)
C
C     VISCOUS AND PLASTIC STRAIN IN MAXWELL ELEMENT
      DO I=1,NLAYER
            DO L=1,6
                  STATEV(29+(I-1)*6+L)=EPLAS_MX(L,I)
                  STATEV(89+(I-1)*6+L)=EVISC_MX(L,I)
            END DO
C
C           CUMULATED PLASTIC STRAIN AND DAMAGE IN MAXWELL LAYERS
            STATEV(1+I)=KAPPA_MX(I)
            STATEV(12+I)=DAM(KAPPA_MX(I))
      END DO     
C
C     Strain Rate
      STATEV(150) = VECNORM(ETOT1-ETOT0,6)/DTIME*1000
C
C     BTVT
      STATEV(23)=PROPS(55)
C
      RETURN
      END
C     __________________________________________________________________
C
C     MATERIAL MODEL FUNCTIONS
C     __________________________________________________________________
C
C     DAMAGE FUNCTION DAM(k)
C 
      DOUBLE PRECISION FUNCTION DAM(KAPPA)
      IMPLICIT NONE
      DOUBLE PRECISION KAPPA
      INTEGER PYFL
      DOUBLE PRECISION KSLOPE,KMAX,KWIDTH,RDY,KMIN,GMIN,ND,EXPS
      DOUBLE PRECISION OFFS,KONSTD,CRITD
      COMMON KSLOPE,KMAX,KWIDTH,RDY,KMIN,GMIN,ND,EXPS,KONSTD
      COMMON CRITD,OFFS,PYFL
C     MEAN KONSTK OF WOLFRAM (JBIOMECH 2011) FOR COMPRESSION AND 
C     TENSION, AXIAL AND TRANSVERSE: 8.0075, CRITD=0.85
      DAM=CRITD*(1.0D0-DEXP(-KONSTD*KAPPA))
      RETURN
      END
C
C     DERIVATIVE OF DAMAGE FUNCTION dDAM/dk
C 
      DOUBLE PRECISION FUNCTION DDAM(KAPPA)
      IMPLICIT NONE
      DOUBLE PRECISION KAPPA
      INTEGER PYFL
      DOUBLE PRECISION KSLOPE,KMAX,KWIDTH,RDY,KMIN,GMIN,ND,EXPS
      DOUBLE PRECISION OFFS,KONSTD,CRITD
      COMMON KSLOPE,KMAX,KWIDTH,RDY,KMIN,GMIN,ND,EXPS,KONSTD
      COMMON CRITD,OFFS,PYFL
C     MEAN KONSTK OF WOLFRAM (JBIOMECH 2011) FOR COMPRESSION AND 
C     TENSION, AXIAL AND TRANSVERSE: 8.0075, CRITD=0.85
      DDAM=CRITD*KONSTD*DEXP(-KONSTD*KAPPA)
      RETURN
      END
C
C     POSTYIELD FUNCTION r(k)
C 
      DOUBLE PRECISION FUNCTION RADK(KAPPA)
      IMPLICIT NONE
      DOUBLE PRECISION KAPPA
      INTEGER PYFL
      DOUBLE PRECISION KSLOPE,KMAX,KWIDTH,RDY,KMIN,GMIN,ND,EXPS
      DOUBLE PRECISION OFFS,KONSTD,CRITD
      COMMON KSLOPE,KMAX,KWIDTH,RDY,KMIN,GMIN,ND,EXPS,KONSTD
      COMMON CRITD,OFFS,PYFL
C 
C     PERFECT PLASTICITY
      IF (PYFL.EQ.0) THEN
        RADK=1.0D0
C     LINEAR HARDENING
      ELSEIF (PYFL.EQ.1) THEN
        RADK=1.0D0+KSLOPE*KAPPA
C     EXPONENTIAL HARDENING
      ELSEIF (PYFL.EQ.2) THEN
        RADK=1.0D0+(RDY-1.0D0)*(1.0D0-DEXP(-KSLOPE*KAPPA))
C     SIMPLE SOFTENING FUNCTION
      ELSEIF (PYFL.EQ.3) THEN
        OFFS=1.0D0/KWIDTH 
        RADK=RDY+(1.0D0-RDY)*(DEXP(-((KAPPA-KMAX)**2)
     1                      /(KWIDTH*KMAX**2))-DEXP(-OFFS-KSLOPE*KAPPA))
C     PIECEWISE SOFTENING FUNCTION
      ELSEIF (PYFL.EQ.4) THEN
        IF (KAPPA.LT.KMAX) THEN
          RADK=1.0D0-((KMAX-KAPPA)/KMAX)**(EXPS*KMAX)
        ELSEIF ((KAPPA.GE.KMAX).AND.(KAPPA.LT.(KMIN+KMAX)/2.0D0)) THEN
          RADK=1.0D0-((1.0D0-GMIN)/2.0D0)*((2.0D0*(KAPPA-KMAX))
     1             /(KMIN-KMAX))**ND
        ELSEIF ((KAPPA.GE.(KMIN+KMAX)/2.0D0).AND.(KAPPA.LT.KMIN)) THEN
          RADK=GMIN+((1.0D0-GMIN)/2.0D0)*((2.0D0*(KMIN-KAPPA))
     1             /(KMIN-KMAX))**ND
        ELSE 
          RADK=GMIN
        ENDIF
        RADK=RDY+(1.0D0-RDY)*RADK
      ENDIF
C
      RETURN
      END
C
C     DERIVATIVE OF POSTYIELD FUNCTION dr/dk
C 
      DOUBLE PRECISION FUNCTION DRADK(KAPPA)
      IMPLICIT NONE
      DOUBLE PRECISION KAPPA
      INTEGER PYFL
      DOUBLE PRECISION KSLOPE,KMAX,KWIDTH,RDY,KMIN,GMIN,ND,EXPS
      DOUBLE PRECISION OFFS,KONSTD,CRITD
      COMMON KSLOPE,KMAX,KWIDTH,RDY,KMIN,GMIN,ND,EXPS,KONSTD
      COMMON CRITD,OFFS,PYFL
C
C     PERFECT PLASTICITY
      IF (PYFL.EQ.0) THEN
        DRADK=0.0D0
C     LINEAR HARDENING
      ELSEIF (PYFL.EQ.1) THEN
        DRADK=KSLOPE
C     EXPONENTIAL HARDENING
      ELSEIF (PYFL.EQ.2) THEN
        DRADK=(RDY-1.0D0)*KSLOPE*DEXP(-KSLOPE*KAPPA)  
C     SIMPLE SOFTENING FUNCTION
      ELSEIF (PYFL.EQ.3) THEN
        OFFS=1.0D0/KWIDTH 
        DRADK=(1.0D0-RDY)*((-KAPPA+KMAX)*DEXP(-(KAPPA-KMAX)**2
     1                        /(KWIDTH*KMAX**2))/(0.5D0*KWIDTH*KMAX**2)
     2                                 +KSLOPE*DEXP(-OFFS-KSLOPE*KAPPA))
C     PIECEWISE SOFTENING FUNCTION
      ELSEIF (PYFL.EQ.4) THEN
        IF (KAPPA.LT.KMAX) THEN
          DRADK=EXPS*((KMAX-KAPPA)/KMAX)**(EXPS*KMAX-1.0D0)
        ELSEIF ((KAPPA.GE.KMAX).AND.(KAPPA.LT.(KMIN+KMAX)/2.0D0))THEN
          DRADK=ND*(GMIN-1.0D0)/(KMIN-KMAX)*((2.0D0*(KAPPA-KMAX))
     1          /(KMIN-KMAX))**(ND-1.0D0)
        ELSEIF ((KAPPA.GE.(KMIN+KMAX)/2.0D0).AND.(KAPPA.LT.KMIN))THEN
          DRADK=ND*(1.0D0-GMIN)/(KMIN-KMAX)*((2.0D0*(KMIN-KAPPA))
     1          /(KMIN-KMAX))**(ND-1.0D0)
        ELSE 
          DRADK=0.0D0
        ENDIF
        DRADK=(1.0D0-RDY)*DRADK
      ENDIF
C 
      RETURN
      END
C
C     __________________________________________________________________
C 
C     MATHEMATICAL FUNCTIONS
C     __________________________________________________________________
C     
C     INTRINSIC FUNCTIONS:
C
C     MATMUL(AAAA,BBBB)
C     MULTIPLICATION MATRIX WITH MATRIX
C
C     MATMUL(AAAA,BB)
C     MULTIPLICATION MATRIX WITH VECTOR
C
C     DOT_PRODUCT(AA,BB)
C     SCALAR PRODUCT OF TWO VECTORS
C 
C     EXTERNAL FUNCTIONS:
C     _________________________________________________________________
C 
      FUNCTION VECDYAD(BB,CC)
      IMPLICIT NONE
C     DYADIC PRODUCT AB OF VECTORS BB AND CC WITH DIMENSION 6x1  
      DOUBLE PRECISION BB(6),CC(6)
      DOUBLE PRECISION VECDYAD(6,6)
      INTEGER I,J
C 
      VECDYAD=0.0D0
C 
      VECDYAD=SPREAD(BB,2,SIZE(CC))*SPREAD(CC,1,SIZE(BB))
C 
      RETURN
      END
C     _________________________________________________________________
C 
      DOUBLE PRECISION FUNCTION VECNORM(FF,N)
      IMPLICIT NONE
C     NORM OF VECTOR AA WITH DIMENSION Nx1
      INTEGER N
      DOUBLE PRECISION FF(N)
C 
      VECNORM=DSQRT(DOT_PRODUCT(FF,FF))
C 
      RETURN
      END
C     __________________________________________________________________ 
C 
C     SUBROUTINES FOR MATRIX INVERSION
C     __________________________________________________________________
C 
      SUBROUTINE MIGS(A,N,X)
C
C     Subroutine to invert matrix A(N,N) with the inverse stored
C     in X(N,N) in the output. A is stored in STOA and is resituted
C     as outpout
C 
      IMPLICIT NONE
C 
      DOUBLE PRECISION A(N,N),STOA(N,N),B(N,N),X(N,N)
      INTEGER INDX(N),N,I,J,K
C 
      DO 140 I=1,N
        DO 130 J=1,N
          STOA(I,J)=A(I,J)
          B(I,J)=0.0D0
  130   CONTINUE
  140 CONTINUE
      DO 150 I=1,N
        B(I,I)=1.0D0
  150 CONTINUE
C 
      CALL ELGS(A,N,INDX)
C 
      DO 180 I=1,N-1
        DO 170 J=I+1,N
          DO 160 K=1,N
            B(INDX(J),K)=B(INDX(J),K)
     1                    -A(INDX(J),I)*B(INDX(I),K)
  160     CONTINUE
  170   CONTINUE
  180 CONTINUE
C 
      DO 210 I=1,N
        X(N,I)=B(INDX(N),I)/A(INDX(N),N)
        DO 200 J=N-1,1,-1
          X(J,I)=B(INDX(J),I)
          DO 190 K=J+1,N
            X(J,I)=X(J,I)-A(INDX(J),K)*X(K,I)
  190     CONTINUE
          X(J,I)= X(J,I)/A(INDX(J),J)
  200   CONTINUE
  210 CONTINUE
C
C Restitution of A
C 
      DO 230 I=1,N
        DO 220 J=1,N
          A(I,J)=STOA(I,J)
  220   CONTINUE
  230 CONTINUE
C 
      RETURN
      END
C     __________________________________________________________________
C 
      SUBROUTINE ELGS(A,N,INDX)
C
C     Subroutine to perform the partial-pivoting Gaussian elimination.
C     A(N,N) is the original matrix in the input and transformed
C     matrix plus the pivoting element ratios below the diagonal in
C     the output.  INDX(N) records the pivoting order.
C 
      IMPLICIT NONE
C 
      DOUBLE PRECISION A(N,N),C(N)
      INTEGER N, INDX(N)
      DOUBLE PRECISION C1,PI1,PI,PJ
      INTEGER K,ITMP,I,J
C
C     Initialize the index
C 
      DO 240 I=1,N
        INDX(I)=I
  240 CONTINUE
C
C     Find the rescaling factors, one from each row
C 
        DO 260 I=1,N
          C1= 0.0
          DO 250 J=1,N
            C1=DMAX1(C1,DABS(A(I,J)))
  250     CONTINUE
          C(I)=C1
  260   CONTINUE
C
C     Search the pivoting (largest) element from each column
C 
      DO 300 J=1,N-1
        PI1=0.0
        DO 270 I=J,N
          PI=DABS(A(INDX(I),J))/C(INDX(I))
          IF (PI.GT.PI1) THEN
            PI1=PI
            K=I
          ELSE
          END IF
  270   CONTINUE
C
C     Interchange the rows via INDX(N) to record pivoting order
C 
        ITMP=INDX(J)
        INDX(J)=INDX(K)
        INDX(K)=ITMP
        DO 290 I=J+1,N
          PJ=A(INDX(I),J)/A(INDX(J),J)
C
C     Record pivoting ratios below the diagonal
C 
          A(INDX(I),J)=PJ
C
C     Modify other elements accordingly
C 
          DO 280 K=J+1,N
            A(INDX(I),K)=A(INDX(I),K)-PJ*A(INDX(J),K)
  280     CONTINUE
  290   CONTINUE
  300 CONTINUE
C 
      RETURN
      END
C     __________________________________________________________________
C
      character*(*) FUNCTION dmkname(fname,dname,exten)
C     COMPOSE A FILENAME AS directory/jobname.exten
C
      character*(*) fname,dname,exten
C     fname  I   jobname
C     dname  I   directory
C     exten  I   extension
C     dmkname O directory/jobname.exten

      ltot = len(fname)
      lf = 0
      DO k1 = ltot,2,-1
        if (lf.EQ.0.AND.fname(k1:k1).NE.' ')  lf = k1
      END DO

      ltot = len(dname)
      ld = 0
      DO k1 = ltot,2,-1
        IF (ld.EQ.0.AND.dname(k1:k1).NE.' ')  ld = k1
      END DO

      ltot = len(exten)
      le = 0
      DO k1 = ltot,2,-1
        IF (le.EQ.0.AND.exten(k1:k1).NE.' ')  le = k1
      END DO

      IF ((lf + ld + le) .LE. len(dmkname)) THEN
        dmkname = dname(1:ld)//'/'//fname(1:lf)
        ltot = ld + lf + 1
        IF ( le.GT.0) THEN
           dmkname = dmkname(1:ltot)//exten(1:le)
        END IF
      END IF
C
      RETURN
      END
C
