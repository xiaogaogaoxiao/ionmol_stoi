      IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION V(14),XS(5,14),NV(10)

  	 ZP=2.D0

      u0=2.8002*10.D0**(-17)/10.D0**(-16) !/zp**2   ! units

      OPEN(UNIT=15, FILE="NV.dat", STATUS='unknown')
	   
            INIT15=0



      NV(1)=1
      NV(2)=4
      NV(3)=5
      NV(4)=6
	DO 100 IROUTE=1,2

      IF(IROUTE.EQ.1) THEN
      OPEN(UNIT=4, FILE="P_Brn.dat", STATUS='unknown')
      OPEN(UNIT=7, FILE="He2_Brn_adn.dat", STATUS='unknown')
      OPEN(UNIT=8, FILE="He2_Brn_CHn.dat", STATUS='unknown')
      OPEN(UNIT=9, FILE="He2_Brn_NOH.dat", STATUS='unknown')
      OPEN(UNIT=10,FILE="He2_Brn_PYR.dat", STATUS='unknown')

      OPEN(UNIT=11,FILE="He2_Brn_adn_nv.dat", STATUS='unknown')
      OPEN(UNIT=12,FILE="He2_Brn_CHn_NV.dat", STATUS='unknown')
      OPEN(UNIT=13,FILE="He2_Brn_NOH_NV.dat", STATUS='unknown')
      OPEN(UNIT=14,FILE="He2_Brn_PYR_NV.dat", STATUS='unknown')

                      ENDIF

      IF(IROUTE.EQ.2) THEN
      OPEN(UNIT=4, FILE="He2_CDW.dat", STATUS='old')
      OPEN(UNIT=7, FILE="He2_CDW_adn.dat", STATUS='unknown')
      OPEN(UNIT=8, FILE="He2_CDW_CHn.dat", STATUS='unknown')
      OPEN(UNIT=9, FILE="He2_CDW_NOH.dat", STATUS='unknown')
      OPEN(UNIT=10,FILE="He2_CDW_PYR.dat", STATUS='unknown')
 
      OPEN(UNIT=11,FILE="He2_CDW_adn_nv.dat", STATUS='unknown')
      OPEN(UNIT=12,FILE="He2_CDW_CHn_NV.dat", STATUS='unknown')
      OPEN(UNIT=13,FILE="He2_CDW_NOH_NV.dat", STATUS='unknown')
      OPEN(UNIT=14,FILE="He2_CDW_PYR_NV.dat", STATUS='unknown')
                      ENDIF
C-----IF cdw

      do 1 IA=1,5

	read(4,*) iv0
	do 2 iv=1,14
      
	read(4,400) v(iv), xs(ia,iv) 
	IF(iv .GT. iv0)  xs(ia,iv)=xs(ia,iv)*ZP**2
  400 format(3x,1pd9.2,19x,1pd10.3)
    2 continue
c      read(4,*)
    1 continue


      WRITE( 7,701)   !"He2_CDW_adn.dat"
      WRITE(11,702)   !"He2_CDW_adn_nv.dat"
 
      WRITE( 8,801)   !"He2_CDW_CHn.dat"
	WRITE(12,802)   !"He2_CDW_CHn_NV.dat"


      WRITE(9 ,901)   !"He2_CDW_NOH.dat"
      WRITE(13,902)   !"He2_CDW_NOH_NV.dat"

      WRITE(10,1001)  !"He2_CDW_PYR.dat"
      WRITE(14,1002)  !"He2_CDW_PYR_NV.dat"
	!

      DO 3 IV=1,14
c     1=H
c     2=C
c     3=N
c     4=O

C----------------------------------------------------------
c     U=Uracilo

c     A=adenina
c     T=Timina

c     C=citosina
c     G=guanina

c     C

      U=4.D0*XS(1,IV)+4.D0*XS(2,IV)+2.D0*XS(3,IV)+2.D0*XS(4,IV)
      A=5.D0*XS(1,IV)+5.D0*XS(2,IV)+5.D0*XS(3,IV)
      T=6.D0*XS(1,IV)+5.D0*XS(2,IV)+2.D0*XS(3,IV)+2.D0*XS(4,IV)
      C=5.D0*XS(1,IV)+4.D0*XS(2,IV)+3.D0*XS(3,IV)+1.D0*XS(4,IV)
      G=5.D0*XS(1,IV)+5.D0*XS(2,IV)+5.D0*XS(3,IV)+1.D0*XS(4,IV)
      
 	WRITE(6,600) V(IV)**2*25.D0, U*u0, A*u0, T*u0, C*u0, G*u0
	WRITE(7,600) V(IV)**2*25.D0, U*u0, A*u0, T*u0, C*u0, G*u0

c..........................................................	
      NVU=4.D0*NV(1)+4.D0*NV(2)+2.D0*NV(3)+2.D0*NV(4)
      NVA=5.D0*NV(1)+5.D0*NV(2)+5.D0*NV(3)
      NVT=6.D0*NV(1)+5.D0*NV(2)+2.D0*NV(3)+2.D0*NV(4)
      NVC=5.D0*NV(1)+4.D0*NV(2)+3.D0*NV(3)+1.D0*NV(4)
      NVG=5.D0*NV(1)+5.D0*NV(2)+5.D0*NV(3)+1.D0*NV(4)

      IF(INIT15.EQ.0) 
     1   WRITE (15,1500) NVU, NVA, NVT, NVC, NVG
 1500 FORMAT(" NVU, NVA, NVT, NVC, NVG= ",10(I3,1X))

	X=0.d0
 	WRITE(6, 600) V(IV)**2*25.D0, U/NVU*V(IV)**X, A/NVA*V(IV)**X,
     |                              T/NVT*V(IV)**X,
     |                              C/NVC*V(IV)**X, G/NVG*V(IV)**X
 	WRITE(11,600) V(IV)**2*25.D0, U/NVU*V(IV)**X, A/NVA*V(IV)**X,
     |                              T/NVT*V(IV)**X,
     |                              C/NVC*V(IV)**X, G/NVG*V(IV)**X

C----------------------------------------------------------
      CH4 =4.D0*XS(1,IV)+1.D0*XS(2,IV)
      C2H2=2.D0*XS(1,IV)+2.D0*XS(2,IV)
      C2H4=4.D0*XS(1,IV)+2.D0*XS(2,IV)
      C2H6=6.D0*XS(1,IV)+2.D0*XS(2,IV)
      C6H6=6.D0*XS(1,IV)+6.D0*XS(2,IV)
      C3H8=8.D0*XS(1,IV)+3.D0*XS(2,IV)

	WRITE(6,600) V(IV)**2*25.D0, CH4*u0,  C2H2*u0, C2H4*u0,
     |                             C2H6*u0, C6H6*u0, C3H8*u0
	WRITE(8,600) V(IV)**2*25.D0, CH4*u0,  C2H2*u0, C2H4*u0,
     |                             C2H6*u0, C6H6*u0, C3H8*u0

c..........................................................	
      NVCH4   =4.D0*NV(1)+1.D0*NV(2)
      NVC2H2  =2.D0*NV(1)+2.D0*NV(2)
      NVC2H4  =4.D0*NV(1)+2.D0*NV(2)
      NVC2H6  =6.D0*NV(1)+2.D0*NV(2)
      NVC6H6  =6.D0*NV(1)+6.D0*NV(2)
      NVC3H8  =8.D0*NV(1)+3.D0*NV(2)


      IF(INIT15.EQ.0) 
     |   WRITE (15,1501) NVCH4, NVC2H2, NVC2H4, NVC2H6, NVC6H6,NVC3H8
 1501 FORMAT(" NVCH4, NVC2H2, NC2H4, NVC2H6, NVC6H6, NVC3H8=",10(I3,1X))

 	WRITE(6,600) V(IV)**2*25.D0,  CH4/NVCH4,  C2H2/NVC2H2, 
     |                              C2H4/NVC2H4,C2H6/NVC2H6,
     |                              C6H6/NVC6H6,C3H8/NVC3H8
 	WRITE(12,600) V(IV)**2*25.D0, CH4/NVCH4,  C2H2/NVC2H2, 
     |                              C2H4/NVC2H4,C2H6/NVC2H6,
     |                              C6H6/NVC6H6,C3H8/NVC3H8


C----------------------------------------------------------
       H2 =2.D0*XS(1,IV)
      XN2 =              2.D0*XS(3,IV)
       O2 =                             2.D0*XS(4,IV)
      XNH3=3.D0*XS(1,IV)+1.D0*XS(3,IV)
       OH2=2.D0*XS(1,IV)              + 1.D0*XS(4,IV)

	WRITE(6,600) V(IV)**2*25.D0, H2*u0,XN2*u0, O2*u0, 
     |                             XNH3*u0,OH2*u0
	WRITE(9,600) V(IV)**2*25.D0, H2*u0,XN2*u0, O2*u0, 
     |                             XNH3*u0,OH2*u0
c..........................................................	
       NVH2 =2.D0*NV(1)
      NVXN2 =           2.D0*NV(3)
       NVO2 =                         2.D0*NV(4)
      NVNH3 =3.D0*NV(1)+1.D0*NV(3)
      NVOH2 =2.D0*NV(1)             + 1.D0*NV(4)

      IF(INIT15.EQ.0) 
     |   WRITE (15,1502) NVH2, NVXN2, NVO2, NVNH3, NVOH2
 1502 FORMAT(" NVH2, NVN2, NVO2, NVNH3, NVOH2= ",10(I3,1X))

	WRITE(13,600) V(IV)**2*25.D0, H2/NVH2,XN2/NVXN2, O2/NVO2, 
     |                              XNH3/NVNH3     ,OH2/NVOH2
	WRITE(6 ,600) V(IV)**2*25.D0, H2/NVH2,XN2/NVXN2, O2/NVO2, 
     |                              XNH3/NVNH3     ,OH2/NVOH2


C----------------------------------------------------------
C  PYRAMIDINAS
      C4H4N2 =4.D0*XS(1,IV)+4.D0*XS(2,IV)+2.D0*XS(3,IV)
	C2H7N1 =7.D0*XS(1,IV)+2.D0*XS(2,IV)+1.D0*XS(3,IV)
	C1H5N1 =5.D0*XS(1,IV)+1.D0*XS(2,IV)+1.D0*XS(3,IV)
      C5H5N1 =5.D0*XS(1,IV)+5.D0*XS(2,IV)+1.D0*XS(3,IV)

	WRITE(6,600)  V(IV)**2*25.D0, C4H4N2*u0, C2H7N1*u0,
     |                              C1H5N1*u0, C5H5N1*u0
	WRITE(10,600) V(IV)**2*25.D0, C4H4N2*u0, C2H7N1*u0,
     |                              C1H5N1*u0, C5H5N1*u0
c..........................................................	
      NVC4H4N2 =4.D0*NV(1)+4.D0*NV(2)+2.D0*NV(3)
	NVC2H7N1 =7.D0*NV(1)+2.D0*NV(2)+1.D0*NV(3)
	NVC1H5N1 =5.D0*NV(1)+1.D0*NV(2)+1.D0*NV(3)
	NVC5H5N1 =5.D0*NV(1)+5.D0*NV(2)+1.D0*NV(3)

      IF(INIT15.EQ.0) 
     |   WRITE (15,1503) NVC4H4N2, NVC2H7N1, NVC1H5N1,NVC5H5N1
 1503 FORMAT(" NVC4H4N2, NVC2H7N1, NVC1H5N1, NVC5H5N1= ",10(I3,1X))

	WRITE(6,600)  V(IV)**2*25.D0, C4H4N2/NVC4H4N2,  C2H7N1/NVC2H7N1,
     |                              C1H5N1/NVC1H5N1,  C5H5N1/NVC5H5N1
	WRITE(14,600) V(IV)**2*25.D0, C4H4N2/NVC4H4N2,  C2H7N1/NVC2H7N1,
     |                              C1H5N1/NVC1H5N1,  C5H5N1/NVC5H5N1
      INIT15=1
                                
    3 CONTINUE

      CLOSE(UNIT=4)

      CLOSE(UNIT=7)
      CLOSE(UNIT=8)
      CLOSE(UNIT=9)
      CLOSE(UNIT=10)

      CLOSE(UNIT=11)
      CLOSE(UNIT=12)
      CLOSE(UNIT=13)
      CLOSE(UNIT=14)
  100 CONTINUE

  600 format(1pd10.3,2x,6(1x,1pd10.3))
 
  701 FORMAT(" ENERGY       URACILO    ADENINA    TIMINA    ",
     |       " CITOSINA   GUANINA")
  702 FORMAT(" ENERGY     SLURACILO  SLADENINA  SLTIMINA  ",
     |       " SLCITOSINA   SLGUANINA")

  801 FORMAT(" ENERGY       CH4       C2H2       C2H4      ",
     |       " C2H6        C6H6       C3H8")
  802 FORMAT(" ENERGY       SLCH4      SLC2H2     SLC2H4     ",
     |       " SLC2H6       SLC6H6     SLC3H8")

  901 FORMAT(" ENERGY       H2         N2         O2       ",
     |       " NH3          OH2  ")
  902 FORMAT(" ENERGY       SLH2        SLN2       SLO2    ",
     |       " SLNH3        SLOH2  ")

 1001 FORMAT(" ENERGY       C4H4N2     C2H7N1    ",
     |       " C1H5N1     C5H5N1"     )
 1002 FORMAT(" ENERGY       SLC4H4N2   SLC2H7N1 ",
     |       " SLC1H5N1   SLC5H5N1" )

      STOP
      END


      