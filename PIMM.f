        PROGRAM PIMM
C     VERSION OCTOBER 2003
C     !!!
C     INPUT FORMATES FOR INTERNAL COORDINATES AND
C     ATOM SORT NUMBERS DIFFER FROM EARLIER VERSIONS
C     !!!
C
C     CLOSED SHELL PI-SCF-LCAO-MO-MOLECULAR MECHANICS PROGRAM
C
C     HANS J. LINDNER,
C     INSTITUT FUER ORGANISCHE CHEMIE DER TECHNISCHEN UNIVERSITAET
C     DARMSTADT,
C     PETERSENSTRASSE 22
C     D-64287 DARMSTADT
C
C     SCF CALCULATION DERIVED FROM A PROGRAM OF U. MUELLER WESTERHOFF,
C     PRIVATE COMMUNICATION 1968.
C     SUBROUTINE DEWGEO USES PARTS OF THE MINDO/2 PROGRAM (M.J.S. DEWAR
C     AND COWORKERS, PRIVATE COMMUNICATION FROM P. BISCHOF, 1973).
C     THE MUTAGH-SARGENT-OPTIMIZATION ROUTINE OF J. PANCIR IS USED IN
C     SUBROUTINE GEO.
C
C     PARAMETRIZATION OF THE PI-SCF-LCAO-MO-CALCULATION:
C     M.J.S. DEWAR AND COWORKERS, PROC. ROY. SOC. LONDON SER A 315,
C     443, 457 (1970), J. AMER. CHEM. SOC. 92, 1453 (1970).
C
C     PARAMETRIZATION OF THE MOLECULAR-MECHANICS-CALCULATION:
C     H. J. LINDNER, TETRAHEDRON 30, 1127 (1974).
C     THIS VERSION CONTAINS PARAMETERS FOR C,N,O,H,BIVALENT S,SI,CL,
C     PRELIMINARY PARAMETERS FOR SULFONIC ACIDS AND DERIVATIVES,
C     PHOSPHORIC ACID, AND PARAMETERS FOR HYDROGEN BONDS. SIMULATION
C     OF LI(+),K(+),CA(2+), MG(2+), O(-), S(-) AND N(+) IS POSSIBLE.
C
C     SIGMA CHARGES ARE ESTIMATED BY THE METHOD OF MARSILLI AND
C     GASTEIGER. THEY INFLUENCE THE PI-SCF-CALCULATION BY
C     SIMPLIFIED CNDO-EQUATIONS.
C
C     This version includes test versions for
C
c                           - calc of crystal packing
c                           - frequency analysis
c                           - Born-Still continuum approximation
c                             for solvation 
C
C-----------------------------------------------------------------------
C
C
C     THIS PROGRAM USES ANGSTROM UNITS - 1 A = 1.E-10 M = 100 PM
C
C
C     INPUT
C     -----
C
C     (1) TITLE
C     FORMAT (I3,I1,19A4/20A4)
C     1. CARD
C     COL 1-3    ISEQ = 0  COMPLETE CALCULATION.
C                     = 1  CALCULATION USES GEOMETRY AND BONDORDER
C                          MATRIX OF THE PREVIOUS CALCULATION;
C                          INPUT ONLY CARDS (1) AND (4).
C                     = 2  PERMUTATION OF TORSIONAL ANGLES
C                          SELECTED BY ITOR;
C                          INPUT ONLY CARDS (1) AND (5).
C                     = 3  RESTART PREVIOUS CALCULATION SAVED ON UNIT 9
C                          INPUT ONLY CARDS (1)
C                     =-1  END OF CALCULATION.
C     COL 4      INFORM    INDICATOR FOR INPUT FORMAT:
C                     = 0  FORMATTED INPUT WITH FORMATES GIVEN
C                          BELOW.
C                     = 1  LIST-DIRECTED READ STATEMENTS.
C                          CONSTANTS MUST BE SEPARATED BY BLANKS
C                          OR COMMAS.
C                          CLOSE INCOMPLETE INPUT RECORDS WITH
C                          A SLASH.
C     COL 5-80   TITLE
C     2. CARD
C     COL 1-80   TITLE
C
C     (2) OPTIONS
C     FOR LIST-DIRECTED INPUT SEE SEQUENCE OF OPTIONS BELOW.
C
C     FORMAT (4I3,F6.2,10I3)
C     COL 1-3    NA   NUMBER OF ATOMS CONTRIBUTING TO THE PI-SYSTEM.
C     COL 4-6    NB   NUMBER OF OCCUPIED PI-ORBITALS.
C     COL 7-9    NC   = 0  NO OUTPUT OF PI-SCF-DATA.
C                          IF ISEQ EQ 2,
C                          ONLY OUTPUT OF PERMUTATED TORSION ANGLES
C                          AND FINAL HEAT OF FORMATION.
C                     = 1  OUTPUT OF PI-SCF-DATA.
C     COL 10-12  ND   OPTION FOR OUTPUT FILES IPLO (INPUT FOR PLOT
C                     PROGRAM PLUTO) AND IPU (OUTPUT OF REFINED
C                     COORDINATES AS NEW PIMM-INPUT):
C                              IPLO        IPU
C                     = 0        NO         NO
C                     = 1        NO        YES
C                     = 2       YES        YES
C                     = 3       YES         NO
c                     if crystal dynamic calculation: ND=0 means 
c                      output of the coordinates of the asymmetric
c                      units only on pdb file.
C     COL 13-18  ZBB  BLOWUP FACTOR, DEFAULT ZBB=2.
C     COL 19-21  NG   = 0  PREPARE PLOT FILE WITH INPUT GEOMETRY.
C                     = 1  NO MOLECULAR-MECHANICS-CALCULATION;
C                          CALCULATION OF HEAT OF FORMATION FOR
C                          INPUT GEOMETRY.
C                     = 2  MOLECULAR-MECHANICS-CALCULATION,
C                          DEFAULT NG=2
C                     = 3  SYMMETRY RESTRICTIONS IN PI-BOND ORDERS,
C                          INPUT: SEE (6).
C                     = 4  MOLECULAR DYNAMICS CALCULATION
C                          INPUT: SEE (7)
C     COL 22-24  NH   =-1  MD only: Monitoring of selected torsional
c                          angles (degrees)--> MD input
c                     =-2  MD only: Monitoring of selected torsional
c                          angles (g+,g-,t) --> MD input          
c                     = 0  UNRESTRICTED OPTIMIZATION
C                     = 1  INDICATES CALCULATION WITH RESTRICTED
C                          OPTIMIZATION. IN PERMUTATIONS, THE INITIAL
C                          STARTING GEOMETRY IS USED FOR EACH STEP
C                     = 2  RESTRICTED OPTIMIZATION WITH PERMUTATIONS RELATIVE 
C                          TO THE OPTIMIZED STRUCTURE FROM THE PREVIOUS STEP 
C     COL 25-27  NI   NUMBER OF SCF-CYCLES, DEFAULT NI=40
C     COL 28-30  NK   GEOMETRY OPTIMIZATION EVERY NK-TH SCF-CYCLE,
C                     DEFAULT NK=3 (NK=1 FOR MD SIMULATIONS)
C     COL 31-33  NL   SCF-TOLERANCE = 5.*(20)**(-NL), DEFAULT NL=5
C     COL 34-36  NM   OPTION FOR INPUT GEOMETRY:
C                     = 1  INPUT OF INTERNAL COORDINATES.
C                          INPUT OF CARTESIAN COORDINATES WITH NM.LE.0:
C                     = 0  INPUT OF ONE SET OF COORDINATES.
C                     =-1  INPUT OF TWO SETS OF COORDINATES, THE SECOND
C                          OF WHICH CAN BE TRANSLATED AND ROTATED WITH
C                          RESPECT TO THE FIRST.
C                     =-2  INPUT OF TWO SETS OF COORDINATES, THE SECOND
C                          BEING TRANSFORMED ACCORDING TO SPECIFIED
C                          POINTS, TO WICH THREE GIVEN ATOMS ARE
C                          APPROXIMATED.
C                     =-3  INPUT OF CARTESIAN COORDINATES WITH
C                          DEFINITION OF BOND MATRIX.
c                     =-4  free.
c                     =-5  MM or MD starting with crystal coordinates.
c                     =-6  MM of MD crystal packing, cell constants fixed.
c                     =-7  mm of MD crystal packing, cell constants optimized.
c
c     
c    
C     COL 37-39  IGEO NUMBER OF CYCLES FOR GEOMETRY OPTIMIZATION,
C                     DEFAULT IGEO=10
C     COL 40-42  NE   NUMBER OF ITERATIONS IN ONE CYCLE OF GEOMETRY
C                     OPTIMIZATION, DEFAULT NE = 50
C     COL 43-45  NEPS = 0  Vacuum calculation
c                     = 1  Solvent calculation,
c                          dielectric constant must be given
c                     = 2  solute calculation,
c                          solvent calculation must precede                          
C     COL 46-48  IHB  = 0  NORMAL TREATMENT OF HYDROGEN BONDS
C                     = 1  IGNORE HYDROGEN BONDS
C     col 49-51  nf   = 0  calc. as usual
C                nf   = 1  analysis of molecular vibrations
c                       For molecules with less than IX/3 atoms only!
c                                                     
C
C     SEQUENCE OF OPTIONS IN LIST-DIRECTED INPUT:
C
C     NA, NB, NC, ND, NM, NG, NH, NI, IGEO, NE
C
C
C     (3) GEOMETRY INPUT
C     DEFINITION OF ATOM SORTS:
C     ATOMS BELONGING TO THE PI-SYSTEM
C         ATOM                SORT
C         C(SP2)               1
C         N(PYRROL TYPE)       2
C         N(PYRIDINE TYPE)     3
C         O(ETHER TYPE)        4
C         O(CARBONYL)          5
C         S(THIOETHER)         6
C         CL                   7
C         F                    8
C         BR                   9 
C         C(SP)               10
C         N(SP)               11
C
C     ATOMS NOT BELONGING TO THE PI-SYSTEM
C         C(SP3)               0
C         MG(2+)              12
C         CA(2+)              13
C         SI                  14
C         P(H3 P O4)          15
C         S(H2 S O4)          16
C         LI(+)               17
C         K(+)                18
C         H                   19
C
C     HETEROATOMS ISOLATED FROM THE PI-SYSTEM
C         N(AMINE TYPE)      102
C         N(+)(AMMONIUM)     202  FOUR BONDS TO ATOMS OF TYPE 0
C                                 OR 9
C         O(ETHER TYPE)      104
C         O(-)(ALCOHOLATE)   204  ONLY ONE NEIGHBOUR
C         S(THIOETHER)       106
C         S(-)(THIOLATE)     106  ONLY ONE NEIGHBOUR
C         CL                 107
C         F                  108
C         BR                 109 
C
C       S-O, P-O, S-N, P-N BONDS OF SULFATES, PHOSPHATES AND THEIR
C       DERIVATIVES ARE TREATED WITHOUT PI-BOND FORMALISM. FOR O AND
C       N ATOM SORTS 102 AND 104 ARE EXPECTED. DOUBLE BONDS ARE
C       SIMULATED BY SHORTER BOND LENGTHS AND ADDITIONAL CHARGES
C       OF +/- 0.333 ELECTRONS ON S/P AND O/N.
C
C
C     (3A) INPUT OF CARTESIAN COORDINATES:
C
C     FORMAT (4I3,F8.4,2F10.4,1X,3I3)
C     COL 1-3    IAT    ATOMIC NUMBER.
C     COL 4-6    ISORT  ATOM SORT NUMBER.
C     COL 7-9    IFIX  =0  ALL COORDINATES ARE OPTIMIZED.
C                      =1  ALL COORDINATES ARE FIXED.
C                      =2  Z-COORDINATE IS FIXED.
C     COL 10-12  IMOL  MOLECULE NUMBER (ATOMS WITH DIFFERENT
C                      IMOL-VALUES ARE NOT BONDED). 
C                      AN IMOL NUMBER OF 99 DENOTES A BOUNDARY MOLECULE
C                      WITH ZERO INTERNAL ENERGY AND ZERO INTERACTION
C                      ENERGY WITH OTHER IMOL=99 MOLECULES, SUITABLE 
C                      FOR SUPRAMOLECULAR CRYSTAL PACKING CALCULATIONS 
C     COL 13-20  X-COORDINATE.
C     COL 21-30  Y-COORDINATE.
C     COL 31-40  Z-COORDINATE.
C
C     COL 42-50  NA(I),NB(I),NC(I). THE THREE ATOMS WICH DEFINE INTERNAL
C       (3I3)    COORDINATES OF (I) (SEE BELOW).
C                IF NA(I) IS NOT SPECIFIED FOR THE FIRST TWO ATOMS, THE
C                PROGRAM CREATES A NEW LIST FOR DEFINITION OF INTERNAL
C                COORDINATES.
C     IF(NM.EQ.-3), THE ATOMS BONDED TO (I) ARE SPECIFIED HERE.
C
C     END INPUT WITH IAT.EQ.0
C
C     IF (NM.EQ.-1.OR.NM.EQ.-2), A SECOND SET OF CARTESIAN COORDINATES
C                                IS READ IN:
C     IF (NM.EQ.-1), THE FIRST CARD OF THIS SET GIVES TRANSLATIONS AND
C                    ROTATIONS:
C     FORMAT (6F10.4)
C                   TRX, TRY, TRZ FOR TRANSLATIONS,
C                   ROX, ROY, ROZ FOR ROTATIONS.
C                   IF NG.EQ.0, THESE VALUES ARE PROMPTED ON TAPE 16
C                   FOR FORMAT FREE INPUT ON TAPE 15.
C
C     IF (NM.EQ.-2), THREE CARDS FOLLOW WITH THE NUMBERS OF THE ATOMS
C                    TO BE APPROXIMATED TO POINTS IN SPACE WICH ARE
C                    SPECIFIED BY THEIR COORDINATES.
C     FORMAT (I3,7X,3F10.4)
C
C     INPUT OF THE SECOND SET FOLLOWS THE DESCRIPTION GIVEN ABOVE.
C
c
c     Subroutine incrys reads crystal data produced by GSTAT with
c     OUTPUT COORDS FRAC
c     with some additional integers
c
c     (1) Title
c     (2) Options
c     as decribed above
c
c     (3) CELL as put out by GSTAT
c     (4) INSERT: ksym, number of SYMM lines
c                 iscl, crystal system indicator
c                       =1 cubic
c                       =2 tetragonal
c                       =3 orthorhombic
c                       =4 monoclinic
c                       =5 triclinic
c     (5) ksym SYMM lines as put out by GSTAT.
c     (6) atomic coordinates - change element symbol to atomic number
c         and PIMM atom sort number, format as given above in sect.(3A)
c
c     (7) two empty lines
c
c     Example:
c
c   0     benzene orthorhombic modification         
c         test
c  6  3  0  2  2.00  2  0 80  3  5 -7 20 20  0
cCELL     7.440   9.550  6.920  90.000  90.000  90.000
c   4  3  1
cSYMM      1.  0   0.  .00000   0.  1.  0.  .00000   0.  0.  1.  .00000 
cSYMM      1.  0.  0.  .50000   0. -1.  0.  .50000   0.  0. -1.  .00000  
cSYMM     -1.  0.  0.  .00000   0.  1.  0.  .50000   0.  0. -1.  .50000
cSYMM     -1.  0.  0.  .50000   0. -1.  0.  .00000   0.  0.  1.  .50000
c  6  1       -.05690    .13870   -.00540
c  6  1       -.13350    .04600    .12640
c  6  1        .07740    .09250   -.12950
c  6  1       -.07740   -.09250    .12950
c  6  1        .05690   -.13870    .00540
c  6  1        .13350   -.04600   -.12640
c  1 19       -.09760    .24470   -.01770
c  1 19       -.24090    .07940    .22180
c  1 19        .13710    .16310   -.23120
c  1 19       -.13710   -.16310    .23120
c  1 19        .24090   -.07940   -.22180
c  1 19        .09760   -.24470    .01770
c
c
c -1
c
C
C     (3B) INPUT INTERNAL COORDINATES
C
C     TRIAL GEOMETRY    ONE CARD PER ATOM
C      1 -  2  NAT(I) I2       ATOMIC NUMBER OF ATOM (I).
C                              NAT(I).LE.0, TO END INPUT OF GEOMETRY.
C      3 -  5 NST(I)  I3       ATOM SORT NUMBER.
C      6 -  8 IFIX(I) I3       IF IFIX(I).EQ.1, THE CARTESIAN
C                              COORDINATES OF THE ATOM (I) ARE FIXED;
C                              IF IFIX(I).EQ.2, THE Z-COORDINATE
C                              IS FIXED.
C      9 - 11 IMOL(I) I3       ATOMS WITH DIFFERENT IMOL-VALUES ARE
C                              NOT BONDED.
C     12 - 20 A(1,I) F10.5     INTERATOMIC SEPARATION (IN ANGSTROMS)
C                              BETWEEN ATOMS (I) AND NA(I).
C     31 - 40 A(2,I) F10.5     THE ANGLE NB(I)-NA(I)-(I) IN DEGREES.
C     51 - 60 A(3,I) F10.5     THE ANGLE BETWEEN THE VECTORS
C                              NC(I)-NB(I) AND NA(I)-(I) IN DEGREES,
C                              MEASURED CLOCKWISE FROM THE DIRECTION
C                              B TO A.
C     61 - 62 ITOR(I) I2       INDICATOR FOR PERMUTATION OF DIHEDRAL
C                              ANGLE:
C                              ITOR(I)=0 NO PERMUTATION.
C                                     =1 0. AND 180. DEGREES.
C                                     =2 60., 180., 300. DEGREES.
C                                     =3 0., 90., 180., 270. DEGREES.
C                                     =4 0., 45., 90., 135.,
C                                        180. DEGREES.
C                                     =5 0., 60., 120., 180.,
C                                        240., 360. DEGREES.
C                              THERE MUST NOT BE MORE THAN
C                              10 ITOR(I).NE.0.
C     71 - 73 NA(I)  I3        THE NUMBER OF THE ATOM, A(1,I) ANG-
C                              STROMS FROM (I). NA(I) MUST BE LESS THAN
C                              (I).
C     74 - 76 NB(I)  I3        NB(I) MUST BE LESS THAN (I). USED TO
C                              DEFINE A(2,I).
C     77 - 79 NC(I)  I3        NC(I) MUST BE LESS THAN (I). USED TO
C                              DEFINE A(3,I).
C
C
C     (4) GEOMETRY VARIATIONS FOR SEQUENTIAL CALCULATION (ISEQ=1):
C     FORMAT (4I3,F8.2)
C        IF NH.GT.0, THE INTERNAL COORDINATE DEFINED BY M1,
C        M2, M3, AND M4 IS FIXED ON VALUE ARSTA.
C     COL  1- 3     M1  FIRST ATOM TO DEFINE INTERNAL COORDINATE.
C     COL  4- 6     M2  SECOND ATOM TO DEFINE BOND LENGTH M1-M2.
C     COL  7- 9     M3  THIRD ATOM TO DEFINE BOND ANGLE M1-M2-M3.
C     COL 10-12     M4  FOURTH ATOM TO DEFINE TORS. ANGLE M1-M2-M3-M4.
C     COL 13-20  ARSTA  NEW VALUE FOR THE INTERNAL COORDINATE
C                       DEFINED BY M1, M2, M3, AND M4.
C     M1.EQ.0 ENDS INPUT.
C
C     (5) VALUES FOR TORSIONAL ANGLE PERMUTATION (ONLY IF ISEQ = 2)
C     1. (I3) PERMUTATION CARD:
C     COL. 1 - 3   KTOR  NUMBER OF PERMUTATED TORSIONAL ANGLES.
C                        IF (KTOR. EQ. 0) INTERNAL LIST OF VALUES
C                        FOR TORSIONAL ANGLE PERMUTATION AS DESCRIBED
C                        FOR ITOR(I) IS USED.
C                        IF (KTOR. GT. 0), THE COMPLETE SET OF VALUES
C                        FOR TORSION ANGLE PERMUTATION MUST BE INPUT.
C     2. (10X,6F10.4) ANGLE SET CARDS:
C     COL 11 -20   BTOR(1,I) FIRST TORSIONAL ANGLE IN THE ANGLE SET
C                  FOR TORSIONAL ANGLE PERMUTATION INDICATED BY I-TH
C                  ITOR(J).NE.0.
C     COL 21 -30   BTOR(2,I) SECOND TORSIONAL ANGLE.
C     ETC.         SPECIFY ITOR(I)+1 ANGLES.
C     NEXT CARD    INPUT OF SECOND ANGLE SET AS DESCRIBED ABOVE.
C     ETC.
C     THERE MUST BE KTOR ANGLE SET CARDS.
C
C     (6) INPUT OF SYMMETRY RESTRICTIONS FOR BOND ORDERS (IF NG.EQ.3):
C     FORMAT (I4) ISBO NUMBER OF BOND ORDER GROUPS WITH EQUAL BOND
C                      ORDERS.
C     FORMAT (I4,20I3) JSBO NUMBER OF BONDS IN FIRST BOND ORDER GROUP
C                      IPA,IPB  CENTRES DEFINING THE BONDS OF THIS
C                      GROUP. THERE MUST BE JSBO PAIRS OF CENTRES.
C            ISBO CARDS OF THIS TYPE ARE NEEDED FOR INPUT.
C
C     (7) INPUT OF MD DATA (IF NG.EQ.4):
c
c     if nh.lt.0: 
c     Format (f10.5) dystar
c     torsion angle monitoring starts after dystar ps.
c     followed by ndyn (max 15) lines defining the torsional angles to be
c     monitored as
c     md1,md2,md3,md4
c     stop with empty line
c
C     FORMAT (2F6.2,3I6)
C     COL  1 - 6   TEMP  STARTING TEMPERATURE
C                        NEGATIVE TEMP CAUSES SYMMETRICAL MOTION
C     COL  7 -12   TIME  TIME STEP IN FEMTOSECONDS, DEFAULT=0.01
C     COL  13-18   MDCYC NUMBER OF STEPS, DEFAULT = 1000
C     COL  18-23   MDUMP WRITE COORDINATE FILE EVERY MDUMP STEPS,
C                  DEFAULT MDUMP=1000
C     COL  24-30   MCOOL COOLING TO ABSOLUTE ZERO TO START AFTER 
C                        STEP MCOOL, DEFAULT MCOOL=0, I.E. NO COOLING
c
C     Format (I6,5F6.1,I6)
C     COL 1-6      IAV   COORDINATES AVERAGED AFTER IAV STEPS
C     COL 7-12     TAUTE COUPLING PARAMETER FOR TEMPERATURE
C     COL 13-18    TAUPE COUPLING PARAMETER FOR PRESSURE
C     COL 19-24    PSTART STARTING PRESSURE
C     COL 25-30    PEND   END PRESSURE
C     COL 31-36    MPRESS PRESSURE STEPS
C
C
c     (8) if neps.eq.1: Solvent calculation
c     col 1-4     Dielectric constant of solvent
c
C     ------------------------------------------------------------------
C
C     USED TAPES:
C
C          INPUT:       IIN = 5
C          OUTPUT:      IUT = 6
C          PLUTO INPUT: IPLO = 1
C          COORDINATES OUTPUT: IPU = 7
C          COORDINATES OUTPUT (BROOKHAVEN FORMAT): IPDB = 2
C          OUTPUT OF CHARGES (MOLCAD FORMAT): IESP=3
C          PERMUTED COORDINATE VS. ENERGY: IARC=8
C          RESTART:     IRES=9 
C
C
C
C-----------------------------------------------------------------------
C
C
C     COMMON AND DIMENSION STATEMENTS
C
C     MAXIMAL NUMBER OF ATOMS IS <IX> (<IX>.LE.999).
C
C     INSERT ACTUAL VALUE FOR <IX> IN ALL PARAMETER-STATEMENTS
C     ( IX=301 --> IX=VALUE ) BEFORE COMPILATION OF THE PROGRAM.
C
C
      CHARACTER*2 ATA,ATN
      double precision xkc,ykc,zkc,vircx,vircy,vircz,viryx,virzx,virzy
      LOGICAL ISTDA
      CHARACTER*4 RESNA,RESNAM
      CHARACTER*5 RESTY,RESTYP
      CHARACTER*30 INFILE,FILE2,FILE3,FILE6,FILE7,FILE8,FILE9,dfil
      PARAMETER(IX=301,IY=IX,IX4=IX*4,IX6=IX*6)
      COMMON /LABEL/ ATA(50)
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      COMMON /DIAG/ F(IY,IY),E(IY),V(IY,IY),ISENT,DBETA(IY)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM2/NUMAT,IATD(IX),NFIX(IX),ILES,NNAT(IX),NAD
     1(IX),NBD(IX),NCD(IX),NDD(IX),NCTOR,MCTOR,NPERM(20),ITOR(IX),
     2NST(IX),NMOL(IX),KTOR
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /BETA/ H(IY,IY)
      COMMON /DYN/ TIME,TEMP,GMASS,MDCYC,TEMITT,EMITT,TMASS,MDUMP,MCOOL,
     1SCHW(3),VELO(3,IX),ACC(3,IX),X3(3,IX),TSOLL,IRAND,TE,iav
     2,pnull,taupe,taupi,pmitt,taute,pstart,pdiff,pend,mpress
      COMMON /GAMMMA/ RR(IX,IX)
      COMMON /STRAIN/ SPANB(IX),SPANA(IX),SPANT(IX),SPANK(IX),SPANE(IX),
     1SPANV(IX)
      COMMON /PDBNAM/ RESNA(IX),NURES(IX),RESNAM(IX),NUMRES(IX),
     1RESTY(IX),RESTYP(IX)
cp
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1 cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2 imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3 ,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
      common /pacplo/ npac,ispac(2000),cpac(3,2000),ihhyd(100),
     1 iahyd(100),dhyd(100),nhyd,jhyd(ix),xo(3,ix),xsum(3,ix),
     2 ibhyd(100),hhyd(100),abhyd(100),anghyd(100)   
cp
cs
      common /solvat/ eps,dielec,dsgl,fav(ix),fam(ix),rxy(ix,ix),
     1 rvdw(50),ion(50)
      common /dyna/ idz,idd(4,40),irdz(ix6),idyn,dystar
cs
      DIMENSION JSBO(20),IPA(20,20),IPB(20,20)
C
C    TAPE ASSIGNMENT
C
      IIN=5
      IUT=6
      IPLO=1
      IPU=7
      IPDB=2
      IESP=3
      IARC=8
      IRES=9
      idyn=20
      CALL GETARG (1,INFILE)
 1003 INQUIRE (FILE=INFILE,EXIST=ISTDA)
      IF (ISTDA) GOTO 1004
      WRITE(*,1005) INFILE
 1005 FORMAT(' Datei ',A15,'existiert nicht - Dateiname ?')
      READ(*,1006) INFILE
 1006 FORMAT(A)
      GOTO 1003
 1004 OPEN(IIN,FILE=INFILE)
      DO 1001 I=1,15
      IF (INFILE(I:I).EQ.'.') GOTO 1002
      IF (INFILE(I:I).EQ.' ') GOTO 1002
 1001 CONTINUE
 1002 I=I-1
      dfil=infile(1:i)
      FILE6=INFILE(1:I)//'P.log'
      FILE7=INFILE(1:I)//'P.opt'
      FILE2=INFILE(1:I)//'P.pdb'
      FILE3=INFILE(1:I)//'P.esp'
      FILE8=INFILE(1:I)//'P.arc'
      FILE9=INFILE(1:I)//'P.rst'
      OPEN(IUT,FILE=FILE6,STATUS='NEW')
      OPEN(IPU,FILE=FILE7,STATUS='NEW')
      OPEN(IPDB,FILE=FILE2,STATUS='NEW')
      OPEN(IESP,FILE=FILE3,STATUS='NEW')
      OPEN(IARC,FILE=FILE8,STATUS='NEW')
      OPEN(IRES,FILE=FILE9,STATUS='UNKNOWN',FORM='UNFORMATTED')
C
C    CALL PARAMETER SUBROUTINE
C
      CALL PARSUB
C
C    CLEAR SOME ARRAYS
C
      ISEQ=0
      ISER=0
      NCTOR=0
      MCTOR=1
  20  DO 31 I=1,IX
         ARST(I)=0.
cs
         fav(i)=0.
         fam(i)=0.
         dsgl=0.
cs
         NUMM(I)=0
         NUMRES(I)=0
         RESTYP(I)='     '
         RESNAM(I)='    '
         QSEFF(I)=0.
cp
         xpac(i)=0.
         jhyd(i)=0
cp
         DO 31 J=1,4
   31 IRST(J,I)=0
      ISER=ISER+1
      IFLAG1=1
      ISYM=0
      EPA=0.
      ISTOP=0
      DO 311 I=1,IY
      E(I)=0.
      DBETA(I)=0.
      DO 311 J=1,IY
      F(I,J)=0.
      V(I,J)=0.
  311 CONTINUE
      IF (ISEQ.EQ.2) GOTO 160
C
C    READ AND WRITE TITLE AND OPTIONS, SET DEFAULT VALUES
C
      READ(IIN,560) ISEQ,INFORM,TITLE
C
C RESTART FILE
C
      IF (ISEQ.NE.3) GOTO 101
      WRITE(IUT,990)
  990 FORMAT(' RESTART FILE READ, CALCULATION RESTARTED')
      READ(IRES) NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,IGEO,ZBB,NOSCCI,IHB,IRZ
      IFLAG1=1
      IFLAG2=1
      NOSCCI=0
      WRITE(IUT,610) NA,NB,NC,ND,ZBB,NG,NH,NI,NK,NL,NM,IGEO,NE
      DO 950 I=1,NAA
      READ(IRES) IAT(I),ISORT(I),IFIX(I),IMOL(I),NNAT(I),QADD(I),
     1CS(I),X(1,I),X(2,I),X(3,I),numm(i)
  950 CONTINUE
      if (nm.lt.-5) then
      read (ires) (ck(i),i=1,3),cal,cbe,cga,ksym,imix,imax,imiy,
     1imay,imiz,imaz,iscl
      do 991 i=1,3
      read (ires) (x(jj,(naa+i)),jj=1,3)
      read (ires) (cl(i,jj),jj=1,3),(ascl(i,jj),jj=1,3)
      read (ires) (tsym(i,jj),jj=1,ksym)
      do 991 k=1,3
  991 read (ires) (rsym(i,k,jj),jj=1,ksym)
      do 992 j=1,naa
  992 read (ires) (xo(ij,j),ij=1,3),jhyd(j)
      end if  
C GEOMETRY RESTRICTIONS
      DO 955 I=1,IRZ
      READ(IRES) IRST(1,I),IRST(2,I),IRST(3,I),IRST(4,I),ARST(I)
  955 CONTINUE
      IF (NA.EQ.0) GOTO 970
      DO 961 I=1,NA
  961 READ (IRES) (H(I,J),P(I,J),j=i,na)
      DO 960 I=1,NA
      DO 960 J=I,NA
      H(J,I)=H(I,J)
      P(J,I)=P(I,J)
  960 CONTINUE
  970 IF (NG.NE.4) GOTO 50
      READ (IRES) TIME,TEMP,GMASS,MDCYC,TEMITT,EMITT,TMASS,MDUMP,MCOOL,
     1SCHW(1),SCHW(2),SCHW(3),TSOLL,te,iav
      DO 980 I=1,NAA
      READ (IRES) VELO(1,I),VELO(2,I),VELO(3,I),ACC(1,I),ACC(2,I),
     1ACC(3,I),X3(1,I),X3(2,I),X3(3,I)
  980 CONTINUE
      GOTO 50
C
C END RESTART SEQUENCE
C
  101 IF (ISEQ.EQ.2) GOTO 160
      IF (ISEQ.EQ.1) GOTO 909
      DO 10 I=1,IX
         X(1,I)=0.
         X(2,I)=0.
         X(3,I)=0.
         IFIX(I)=0
         IMOL(I)=0
         ITOR(I)=0
         QADD(I)=0.
         QSEFF(I)=0.
   10    ISORT(I)=0
      DO 11 I=1,IY
         E(I)=0.
      DO 11 J=1,IY
         F(I,J)=0.
         P(I,J)=0.
   11 H(I,J)=0.
      DO 30 I=1,IX
         SWI(I)=0.
         CS(I)=0.
      DO 30 J=1,IX
         VV(I,J)=0.
   30 VS(I,J)=0.
  909 WRITE (IUT,760) ISER
      WRITE(IUT,600) TITLE,ISEQ
      IF (ISEQ.LT.0) STOP
      IF (ISEQ.EQ.1) GOTO 110
      WRITE (IUT,601) INFORM
      IF (INFORM.EQ.0) THEN
         READ(IIN,570) NA,NB,NC,ND,ZBB,NG,NH,NI,NK,NL,NM,IGEO,NE,NEPS
     1  ,IHB,nf
      ELSE
         READ(IIN,*) NA,NB,NC,ND,NM,NG,NH,NI,IGEO,NE,NEPS,nf
      END IF
cs
      if (neps.gt.1.and.dielec.eq.0.) then
      write (iut,830)
  830 format ('ERROR - MISSING DIELECTRIC CONSTANT')
      stop
      end if
cs  
      IF (NI.LT.0) THEN
      NI=10**(-NI)
      IGEO=NI/NL
      END IF
      IF (NE.EQ.0) NE=50
      IF (NL.EQ.0) NL=5
      ZL=10.**(-NL)*5.0
      IF (NI.EQ.0) NI=40
      IF (NK.EQ.0) NK=3
      IF (ZBB.LE.0.) ZBB=2.
      IF (IGEO.EQ.0) IGEO=10
      IF (NG.LT.2) IGEO=1
      if (ng.eq.4) nk=1
      WRITE(IUT,610) NA,NB,NC,ND,ZBB,NG,NH,NI
      write(iut,611) NK,NL,NM,IGEO,NE,NEPS,IHB,nf
C
C    READ COORDINATES
C
      ILES=0
      CALL DEWGEO
      NAA=NUMAT
cs
cs read dielectric constant and calculate eps=(e-1)/e
cs
      if (neps.eq.0) eps=0.
      if (neps.eq.1) then
      read (iin,fmt='(f12.2)') eps
      write (iut,841) eps
  841 format ('       DIELECTRIC CONSTANT', f10.2)
      if (eps.eq.0.) then
      write (iut,840)
  840 format (' ERROR - NO DIELECTRIC CONSTANT GIVEN')
      stop
      end if
      dielec=eps
      eps=(eps-1.)/eps
      end if
cs
C
C    WRITE CARTESIAN COORDINATES
C
      WRITE(IUT,620)
      DO 100 I=1,NAA
         L=NNAT(I)
         NUMM(L)=I
         ISORT(L)=MOD(ISORT(L),100)
         ISO=ISORT(L)+1
         ATN=ATA(ISO)
  100 IF (NC.GE.0)
     1WRITE(IUT,630)ATN,I,IAT(L),ISORT(L),IMOL(L),IFIX(L),(X(J,L),J=1,3)
C
C    READ AND STORE NEW VALUES FOR INTERNAL COORDINATES DEFINED BY M1
C    TO M4 FOR A NEW CALCULATION. IF NH.GT.0, THIS VALUE WILL BE FIXED
C    DURING GEOMETRY OPTIMIZATION.
C
  110 IRZ=1
  120 IF (INFORM.EQ.0) THEN
         READ(IIN,590) M1,M2,M3,M4,ARSTA
      ELSE
         READ(IIN,*) M1,M2,M3,M4,ARSTA
      END IF
      IF (M1.GT.0) THEN
         MM1=NNAT(M1)
         MM2=NNAT(M2)
         IRST(1,IRZ)=MIN0(MM1,MM2)
         IRST(2,IRZ)=MAX0(MM1,MM2)
         IF (M3.EQ.0) GOTO 130
         ARSTA=ARSTA*0.0174532925
         MM3=NNAT(M3)
         IRST(1,IRZ)=MM2
         IRST(2,IRZ)=MIN0(MM1,MM3)
         IRST(3,IRZ)=MAX0(MM1,MM3)
         IF (M4.EQ.0) GOTO 130
         MM4=NNAT(M4)
         IRST(1,IRZ)=MIN0(MM2,MM3)
         IRST(2,IRZ)=MAX0(MM2,MM3)
         IRST(3,IRZ)=MM1
         IRST(4,IRZ)=MM4
         IF (MM3.GT.MM2) GO TO 130
         IRST(3,IRZ)=MM4
         IRST(4,IRZ)=MM1
  130    ARST(IRZ)=ARSTA
         IF(M4.NE.0) ARST(IRZ)=-ARST(IRZ)
         IF(M3.NE.0) ARSTA=ARSTA*57.29577951
C
C    OUTPUT OF THE NEW INTERNAL COORDINATE
C    STORE NEW INTERNAL COORDINATE IN THE INTERNAL COORDIATE STORAGE
C
         WRITE(IUT,640) M1,M2,M3,M4,ARSTA
         DO 140 I=1,NAA
         DO 140 J=1,3
         IF (IRST(1,IRZ).NE.NRA(1,J,I)) GO TO 140
         IF (IRST(2,IRZ).NE.NRA(2,J,I)) GO TO 140
         IF (IRST(3,IRZ).NE.NRA(3,J,I)) GO TO 140
         IF (IRST(4,IRZ).NE.NRA(4,J,I)) GO TO 140
         A(J,I)=ARST(IRZ)
  140    CONTINUE
         IRZ=IRZ+1
         GOTO 120
      END IF
cmm
cm
      idz=0
      if (nh.lt.0) then
         read(iin,5590) dystar
 5590    format (f10.5)
         write (iut,8141) dystar
 8141    format (///,'Collecting starts after',f6.1,' ps',///)
 8120 IF (INFORM.EQ.0) THEN
         READ(IIN,590) Md1,Md2,Md3,Md4
      ELSE
         READ(IIN,*) Md1,Md2,Md3,Md4
      END IF
      IF (Md1.GT.0) THEN
         idz=idz+1
         MMd1=NNAT(Md1)
         MMd2=NNAT(Md2)
         MMd3=NNAT(Md3)
         MMd4=NNAT(Md4)
         Idd(1,IdZ)=MIN0(MMd2,MMd3)
         Idd(2,IdZ)=MAX0(MMd2,MMd3)
         Idd(3,IdZ)=MMd1
         Idd(4,IdZ)=MMd4
         IF (MMd3.GT.MMd2) GO TO 8130
         Idd(3,IdZ)=MMd4
         Idd(4,Idz)=MMd1
 8130    continue
         write (iut,8140) md1,md2,md3,md4,idz
 8140    format(' Torsional angle',4i3,' monitored as dy',i2)
         goto 8120
         end if
      dfil=dfil(1:6)//'P.dyn'
      open(idyn,file=dfil,status='new')
      end if
cm
  160 IF (ISEQ.EQ.0) GOTO 170
      ILES=1
      WRITE (IUT,770) ISER
      IF (ISEQ.EQ.2.AND.NC.GT.0) WRITE(IUT,600) TITLE,ISEQ
      IF (ISEQ.EQ.2) CALL TORVAR(*20)
      CALL DEWGEO
C
C    OUTPUT OF NEW CARTESIAN COORDINATES
C
      IF (ISEQ.NE.2.AND.NC.GT.0) WRITE(IUT,620)
      DO 161 I=1,NAA
         L=NNAT(I)
         NUMM(L)=I
         ISO=ISORT(L)+1
         IF (ISEQ.NE.2.AND.NC.GT.0) THEN
            ATN=ATA(ISO)
            WRITE(IUT,630)
     1      ATN,I,IAT(L),ISORT(L),IMOL(L),IFIX(L),(X(J,L),J=1,3)
         END IF
  161 CONTINUE
  170 IF (NG.EQ.3) THEN
C
C    INPUT OF SETS OF BOND ORDERS THAT WILL BE KEPT EQUAL, IF
C    NG.EQ.3
C
      ISYM=1
      IF (INFORM.EQ.0) READ (IIN,800) ISBO
      IF (INFORM.EQ.1) READ (IIN,*) ISBO
      DO 172 I=1,ISBO
         IF (INFORM.EQ.0)
     1   READ (IIN,810) JS,(IPA(I,J),IPB(I,J),J=1,JS)
         IF (INFORM.EQ.1)
     1   READ (IIN,*) JS,(IPA(I,J),IPB(I,J),J=1,JS)
         WRITE(IUT,820) JS,(IPA(I,J),IPB(I,J),J=1,JS)
  172 JSBO(I)=JS
  800 FORMAT (I4)
  810 FORMAT (I4,20I3)
  820 FORMAT (1H0,10X,'EQUAL BOND ORDERS',I5,10(I5,I3))
      END IF
C
C    CALCULATE TRIAL ATOM-ATOM DISTANCES, SELECT BONDS
C    PERMUTATIONS REUSE THE OLD BOND MATRIX
C
   50 IF (ISEQ.EQ.1) GOTO 182
      DO 179 I=1,NAA
      DO 179 J=I+1,NAA
         VV(I,J)=SQRT((X(1,I)-X(1,J))**2+(X(2,I)-X(2,J))**2+(X(3,I)-X(3
     1   ,J))**2)
         VV(J,I)=VV(I,J)
  179 CONTINUE
      IF (NM.NE.-3) THEN
      DO 180 I=1,NAA
      DO 180 J=1,NAA
         ISO1=ISORT(I)+1
         ISO2=ISORT(J)+1
         RMAX=1.1*(RBOND(ISO1)+RBOND(ISO2))
         VS(I,J)=0.
         IF (VV(I,J).GT.0..AND.VV(I,J).LT.RMAX) VS(I,J)=VV(I,J)
         IF (IMOL(I).NE.IMOL(J)) VS(I,J)=0.
  180 IF (ISOrt(i).eq.19.and.isort(j).eq.19) VS(I,J)=0.
      ELSE
         DO 186 I=1,NAA
            I1=NUMM(I)
            IF(NAD(I1).EQ.0) GOTO 186
            J1=NAD(I1)
            NAD(I1)=0
            J=NNAT(J1)
            IS1=ISORT(I)+1
            IS2=ISORT(J)+1
            VS(I,J)=AL(IS1,IS2)
            VS(J,I)=VS(I,J)
            IF(NBD(I1).EQ.0) GOTO 186
            J1=NBD(I1)
            J=NNAT(J1)
            IS1=ISORT(I)+1
            IS2=ISORT(J)+1
            VS(I,J)=AL(IS1,IS2)
            VS(J,I)=VS(I,J)
            IF(NCD(I1).EQ.0) GOTO 186
            J1=NCD(I1)
            J=NNAT(J1)
            IS1=ISORT(I)+1
            IS2=ISORT(J)+1
            VS(I,J)=AL(IS1,IS2)
            VS(J,I)=VS(I,J)
            IF(NDD(I1).EQ.0) GOTO 186
            J1=NDD(I1)
            J=NNAT(J1)
            IS1=ISORT(I)+1
            IS2=ISORT(J)+1
            VS(I,J)=AL(IS1,IS2)
            VS(J,I)=VS(I,J)
  186   CONTINUE
      END IF
C
C    TEST FOR ERRORS IN GEOMETRY INPUT
C
      CALL GEOTES
      WRITE (IUT,4711) NAA,IX
 4711 FORMAT(//,' TOTAL NUMBER OF ATOMS IS ',I4,' (MAXIMUM ',I4,')',/)
C
C    IF NG.LT.0, PLOT INPUT AND START NEXT CALCULATION
C
      IF (NG.LE.0) CALL PLUPLO
      IF (ISTOP.GT.0.OR.NG.EQ.0) GOTO 20
      IF (ISEQ.EQ.2.AND.NC.EQ.0) GOTO 182
      IF (NG.EQ.4) GOTO 182
      WRITE (IUT,740)
  182 NAAM=NAA-1
      NAM=NA-1
C
C    ENTER SETGEO TO FIND AND STORE ALL BOND LENGTHS, BOND ANGLES,
C    TORSION ANGLES, AND BEND ANGLES.
C
C    ENTER CHARGE FOR FIRST COMPUTATION OF SIGMA-CHARGES
C
C    ENTER SUBROUTINE GEO FOR FIRST GEOMETRY CALCULATION
C    IF NO PI-SYSTEM PRESENT, COMPLETE OPTIMIZATION
C
      CALL SETGEO
c opt
      anaa=naa*2
      if (neps.gt.0) anaa=anaa*0.5
      convmm=amin1(0.5e-3,4.e-3*sqrt(1./anaa))
      coufac=1.
      IF (NG.EQ.4.AND.ISEQ.NE.3) CALL MDINIT
      IF (NA.GT.0.AND.ISEQ.NE.3) GOTO 190
      DO 181 I=1,IGEO
         NN=NE
         IF (I.EQ.IGEO) IFLAG1=0
         IF (I.EQ.IGEO) IFLAG2=IGEO
         CALL GAMMA(NAA,ISORT,RV)
         CALL CHARGE
      nn=ne
      CALL GEO(ifin,convmm)
      IF (NG.EQ.4) GOTO 185
  181 if (ifin.eq.1) goto 420
      GOTO 420
  185 CALL CHARGE
      DO 187 I=1,MDCYC
      IF (TEMP.LE.0.5) iflag2=mdcyc
      IFLAG2=I
      CALL MDYN
      IF (MOD(I,MDUMP).EQ.0) THEN
      CALL COROUT
      CALL PLUPLO
      END IF
  187 CONTINUE
      GOTO 420
  190 NN=ne
      NGT=NG
      NG=1
      CALL GAMMA(NAA,ISORT,RV)
      CALL CHARGE
      CALL GEO(ifin,convmm)
      IF (NGT.EQ.4) THEN
      IFLAG2=1
      CALL MDYN
      END IF
      NG=NGT
c
C SKIP HMO CALCULATION WHEN RESTARTING
C
      IF (ISEQ.EQ.3) GOTO 240
C
      DO 200 I=1,NA
         ISO1=ISORT(I)
         ISP1=ISO1+1
         f(i,i)=0.
         H(I,I)=HV(ISO1)
C
C    SET UP HUCKEL MATRIX FOR INITIAL HMO-CALCULATION
C    SET UP FIRST H-MATRIX
C
C
      DO 2001 J=1,I-1
         IF (VS(I,J).LE.0.) THEN
         F(I,J)=0.
         ELSE
         ISO2=ISORT(J)
         ISP2=ISO2+1
         f(i,j)=1.
         f(j,i)=f(i,j)
         END IF
 2001 continue
  200 CONTINUE
C
C    SUBROUTINE BOSYM SETS SELECTED MATRIX ELEMENTS TO THEIR MEAN
C    VALUES (HERE ELEMENTS F(I,J)).
C
      IF (ISYM.EQ.1) CALL BOSYM(ISBO,JSBO,IPA,IPB,F)
C
C    DIAGONALIZE HUCKEL MATRIX
C
      ISENT=0
      CALL THOQR
C
C    CALCULATE HMO BOND ORDER MATRIX
C
      DO 230 I=1,NA
      DO 220 J=1,NA
         SI=0.
         DO 210 K=1,NB
  210    SI=SI+V(K,I)*V(K,J)
  220 P(I,J)=2.*SI
  230 CS(I)=ZV(ISORT(I))-P(I,I)
C
C    MEAN VALUES FOR BOND ORDERS
C
      IF (ISYM.EQ.1) CALL BOSYM(ISBO,JSBO,IPA,IPB,P)
c
      IFLAG2=1
      IFLAG1=1
      NOSCCI=-10
C
C    BEGIN OF THE SCF-CYCLE
C
  240 NOSCCI=NOSCCI+1
c
      if (ng.lt.4) then
      call gamma(NAA,ISORT,RV)
      call charge
      call betov
      end if
      if (ng.eq.4) then
      if (mod(iflag2,10).ne.0) goto 341
      call gamma(NAA,ISORT,RV)
      call charge
      call betov
      end if 
C
C    SET UP F-MATRIX
C
      DO 270 I=1,NA
         ISO1=ISORT(I)
         ISP1=ISO1+1
cs
         rrii=0.25*eps*rxy(i,i)*cs(i)
         F(I,I)=-HV(ISO1)-rrii+(P(I,I)*0.5-QSIG(I)+QQ(ISP1))*RR(I,I)
         DO 270 J=NA+1,NAA
cs
         rrij=rr(i,j)
         if (vs(i,j).le.0.) rrij=rr(i,j)-0.25*eps*rxy(i,j)
  270    F(I,I)=F(I,I)-QSIG(J)*RRIJ
      DO 280 I=2,NA
      DO 280 J=1,I-1
cs
        rrij=rr(i,j)
        F(I,I)=F(I,I)-(CS(J)+QSIG(J))*RRIJ
        F(J,J)=F(J,J)-(CS(I)+QSIG(I))*RRIJ
        F(I,J)=H(I,J)-(0.5*P(I,J)-QSIG(I)*QSIG(J))*RRIJ
  280 F(J,I)=F(I,J)
cp
      if (nm.lt.-5) then
      do 281 i=1,na
  281 f(i,i)=f(i,i)-xpac(i)
      end if
cp
C
C    MEAN VALUES FOR FOCK-MATRIX ELEMENTS
C
      IF (ISYM.EQ.1) CALL BOSYM(ISBO,JSBO,IPA,IPB,F)
C
C    DIAGONALIZE F-MATRIX
C
      ISENT=1
      CALL THOQR
      IFLAG1=0
C
C    CALCULATE BOND ORDER MATRIX AND TEST FOR SELFCONSISTENCY
C
      DO 345 I=1,NA
      DO 340 J=1,NA
         SI=0.
         DO 310 K=1,NB
  310    SI=SI+V(K,I)*V(K,J)
         PP=2.*SI
         DPJ=PP-P(I,J)
         IF (ABS(DPJ).GT.5.E-2) THEN
         IFLAG3=1
         IFLAG1=1
         P(I,J)=P(I,J)+SIGN(5.E-2,DPJ)
         GOTO 340
         END IF
         IF (ABS(DPJ)-ZL) 330,330,320
  320    IFLAG1=1
  330    P(I,J)=PP
  340    CONTINUE
  345 CS(I)=ZV(ISORT(I))-P(I,I)
        IF (IFLAG3.NE.0) THEN
        WRITE(IUT,343)
  343   FORMAT (' WARNING - SCF CALCULATION MAY FAIL')
        IFLAG3=0
        END IF
C
C    MEAN VALUES FOR BOND ORDERS
C    SWITCH OFF SCF TEST,IF NG.EQ.3
C
      IF (ISYM.EQ.1) CALL BOSYM(ISBO,JSBO,IPA,IPB,P)
      IF (ISYM.EQ.1.AND.IFLAG2.GE.IGEO) IFLAG1=0
C
C    CHECK, IF GEOMETRY OPTIMIZATION MUST FOLLOW
C
  341 IF (noscci.lt.0.or.MOD(NOSCCI,NK).NE.0) GOTO 240
C
C    CALCULATE EQUILIBRIUM BOND LENGTHS FOR PI-BONDS FROM PI-BOND
C    ORDERS
C
      IF (ISYM.EQ.1) CALL QSYM(ISBO,JSBO,IPA,IPB,QSIG)
c
      DO 370 I=1,NAM
         ISO1=ISORT(I)
         ISP1=ISO1+1
         I1=I+1
      DO 370 J=I1,NA
         IF (VS(I,J).LE.0.) GOTO 370
         ISO2=ISORT(J)
         ISP2=ISO2+1
         VS(I,J)=AL(ISP1,ISP2)-BL(ISO1,ISO2)*P(I,J)
         IB=0
         DO 360 K=1,NA
            IF (K.EQ.I.OR.K.EQ.J) GOTO 360
            IF (VS(I,K).GT.0.) IB=IB+1
            IF (VS(J,K).GT.0.) IB=IB+1
  360    CONTINUE
         IF (IB.EQ.4) THEN
            VS(I,J)=VS(I,J)+0.015
         END IF
         VS(J,I)=VS(I,J)
  370 CONTINUE
      NN=NE
      IF (NG.EQ.4) GOTO 380
      IF (IFLAG1.NE.0.AND.IFLAG2.GT.IGEO.AND.NG.NE.1) GOTO 390
C
C   ENTER SUBROUTINE GEO TO OPTIMIZE GEOMETRY
C
      CALL GEO(ifin,convmm)
  380 IF (NG.EQ.4) CALL MDYN
      IFLAG2=IFLAG2+1
  390 continue
c
      IF (ISYM.EQ.1) CALL BOSYM(ISBO,JSBO,IPA,IPB,H)
C
C    CALCULATE SIGMA-BOND ENERGY
C
  420 SUMSIG=0.
      DO 430 I=1,NAA
        IF (IMOL(I).EQ.99) GOTO 430
         ISO1=ISORT(I)+1
         IF (ISO1.EQ.1) THEN
            IB=0
            DO 432 M=1,NAA
  432       IF (VS(I,M).GT.0..AND.ISORT(M).EQ.0) IB=IB+1
            IF (IB.EQ.3) SUMSIG=SUMSIG+0.10
            IF (IB.EQ.4) SUMSIG=SUMSIG+0.30
         END IF
      DO 431 J=I,NAA
         IF (VS(I,J).LE.0.) GOTO 431
         ISO2=ISORT(J)+1
         EBSIG=(DE(ISO1,ISO2)+DF(ISO1,ISO2)*(QSIG(I)*QSIG(J)+
     1   QSIG(I)*CS(J)+CS(I)*QSIG(J))*RR(I,J))
     2   *(1.-(1.-EXP(AME(ISO1,ISO2)*(AL(ISO1,ISO2)-VV(I,J))))**2)
         SUMSIG=SUMSIG+EBSIG
  431 CONTINUE
  430 CONTINUE
C
C    CALCULATE PI-BOND ENERGY
C
      EPIB=0.
      IF (NA.EQ.0) GOTO 460
      DO 450 I=1,NA
        IF (IMOL(I).EQ.99) GOTO 450
         ISO1=ISORT(I)
         ISP1=ISO1+1
cs
         rrii=0.25*eps*rxy(i,i)*cs(i)
         EPIB=EPIB+(0.25*P(I,I)*P(I,I)+0.5*(QSIG(I)-QQ(ISP1))*CS(I))*
     1   RR(I,I)-(hv(iso1)+rrii)*p(i,i)+h(i,i)*zv(iso1)
     2   -(ZV(ISO1)-1.)*rr(i,i)
         DO 440 J=1,I-1
         rrij=rr(i,j)
            EPIB=EPIB+2.*P(I,J)*(H(I,J)+QSIG(I)*QSIG(J)*RRIJ)
     1      -(0.5*P(I,J)*P(I,J)-CS(I)*CS(J))*RRIJ
  440    CONTINUE
  450 CONTINUE
C
C    CALCULATE HEAT OF FORMATION FROM SIGMA-, PI-ENERGIES AND
C    PARTS OF THE STRAIN ENERGY, CALCULATED IN SUBROUTINE GEO
C
  460 EGES=EPIB-SUMSIG+EP
      DHFM=EGES*23.061
      DO 470 I=1,NAA
c       IF (IMOL(I).EQ.99) GOTO 470
         NS=ISORT(I)+1
C
         IF (NS.EQ.5.AND.QADD(I).LE.-1.) DHFM=DHFM+18.2
         IF (NS.EQ.7.AND.QADD(I).LE.-1.) DHFM=DHFM-12.1
         IF (NS.EQ.3.AND.QADD(I).GE.1.)  DHFM=DHFM+79.5
C
C    CORRECTION OF HEAT OF FORMATION FOR O(-) ESTIMATED FROM
C    ELECTRON AFFINITIES OF METHOXY RADICALS
C    AND FOR S(-) FROM ELECTRON AFFINITIES OF MERCAPTO RADICALS.
C    LIT: B.K. JANOUSEK AND J.I. BRAUMANN, GAS PHASE ION CHEMISTRY,
C    VOL 2, ED. M.T. BOWERS, ACADEMIC PRESS, 1979, P. 53
C
         DHFM=DHFM+DHFA(NS)
  470 CONTINUE
c
c    add solvatation energy
c
      dhfm=dhfm+dsgl    
      DHFM=DHFM*4.184
      IF (ISEQ.EQ.2.AND.NC.LE.0) GOTO 475
      IF (NG.EQ.4) GOTO 475
      WRITE(IUT,650) EGES,EPIB,EP,DHFM
      IF (ISEQ.NE.1) GOTO 475
      IF (NA.GT.0.AND.IFLAG2.LE.IGEO) GOTO 475
      IF (IFLAG1.EQ.1.AND.NOSCCI.LE.NI) GOTO 475
      ARST1=-57.29577951*ARST(1)
      IF (ARST1.GT.180.) ARST1=ARST1-360.
      ARST2=-57.29577951*ARST(2)
      IF (ARST2.GT.180.) ARST2=ARST2-360.
      WRITE (IARC,474)ARST1,ARST2,DHFM
  474 FORMAT (3F8.1)
  475 IF (NA.EQ.0) GOTO 530
  480 IF (NG.EQ.4.AND.IFLAG2.LE.MDCYC) THEN
      IF (MOD(IFLAG2,MDUMP).EQ.0.and.nh.ge.0) THEN
      CALL COROUT
      CALL PLUPLO
      END IF
      IF (TEMP.LE.0.) GOTO 490
      GOTO 240
      END IF
      IF (NOSCCI.GT.NI) GOTO 490
      if (iflag1.eq.0.and.ifin.eq.1) goto 500
      IF (IFLAG2.LE.IGEO.OR.IFLAG1.EQ.1)  GOTO 240
C
C    END OF SCF-CYCLE
C
C    OUTPUT OF SCF-DATA AND MATRIX OF ATOM-ATOM DISTANCES
C
      GOTO 500
  490 IF (ISEQ.EQ.2.AND.NC.LE.0) GOTO 530
      WRITE(IUT,660)
  500 IF (ISEQ.EQ.2.AND.NC.LE.0) GOTO 530
      WRITE(IUT,670) NOSCCI
C
C    OUTPUT OF GEOMETRY DATA
C
      CALL OUTGEO
C     IF (ISEQ.NE.1) GOTO 505
C     ARST1=-57.29577951*ARST(1)
C     IF (ARST1.GT.180.) ARST1=ARST1-360.
C     ARST2=-57.29577951*ARST(2)
C     IF (ARST2.GT.180.) ARST2=ARST2-360.
C     WRITE (IARC,474)ARST1,ARST2,DHFM
C 474 FORMAT (3F8.1)
  505 IF (NC.LE.0) GOTO 530
      WRITE(IUT,680)
      DO 510 I=1,NA
  510 WRITE(IUT,690) (NUMM(I),NUMM(J),P(I,J),J=I,NA)
      WRITE(IUT,700)
      WRITE(IUT,710) (I,E(I),I=1,NA)
      WRITE(IUT,720)
      DO 520 I=1,NA
  520 WRITE(IUT,690) (I,NUMM(J),V(I,J),J=1,NA)
  530 IF (NA.EQ.0) CALL OUTGEO
      IF (NC.LE.0) GOTO 541
      WRITE(IUT,730)
      DO 540 I=1,NAAM
         II=I+1
  540 WRITE(IUT,690) (NUMM(I),NUMM(J),VV(I,J),J=II,NAA)
C
C    SET INTERNAL COORDINATES FOR COORDINATE OUTPUT
C
  541 IF (ND.GT.0.AND.NH.NE.1) CALL CORSET
C
C    PREPARE PLOT FILE
C
      IF (ND.GT.1) CALL PLUPLO
C
C    PREPARE COORDINATE OUTPUT FILE
C
      IF (ND.GT.0.AND.ND.LT.3) CALL COROUT
      IF (ISEQ.EQ.2.AND.NC.LE.0) GOTO 545
cs 
cs    write some solvatation data
cs
      if (neps.gt.0) call SOLVut(NAA,ISORT,NUMM,cs,qseff)
cs
      if (iseq.eq.0) then
      write(iarc,559) dhfm,(title(i),i=1,9)
  559 format (f12.1,9a4,2x,2f8.1)
      end if
c
C WRITE RESTART FILE
C
      if (nc.le.0) goto 550
      REWIND(IRES)
      WRITE(IRES) NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,
     1IFLAG2,NAA,IGEO,ZBB,NOSCCI,IHB,IRZ
      DO 900 I=1,NAA
      WRITE (IRES) IAT(I),ISORT(I),IFIX(I),IMOL(I),NNAT(I),QADD(I),
     1CS(I),X(1,I),X(2,I),X(3,I),numm(i)
  900 CONTINUE
      if (nm.lt.-5) then
      write (ires) (ck(i),i=1,3),cal,cbe,cga,ksym,imix,imax,imiy,
     1imay,imiz,imaz,iscl
      do 901 i=1,3
      write (ires) (x(jj,(naa+i)),jj=1,3)
      write (ires) (cl(i,jj),jj=1,3),(ascl(i,jj),jj=1,3)
      write (ires) (tsym(i,jj),jj=1,ksym)
      do 901 k=1,3
  901 write (ires) (rsym(i,k,jj),jj=1,ksym)
      do 902 j=1,naa
  902 write (ires) (xo(ij,j),ij=1,3),jhyd(j)
      end if  
C GEOMETRY RESTRICTIONS
      DO 905 I=1,IRZ
      WRITE(IRES) IRST(1,I),IRST(2,I),IRST(3,I),IRST(4,I),ARST(I)
  905 CONTINUE
      IF (NA.EQ.0) GOTO 920
      DO 910 I=1,NA
      WRITE (IRES) (H(I,J),P(I,J),j=i,na)
  910 CONTINUE
  920 IF (NG.NE.4) GOTO 550
      WRITE (IRES) TIME,TEMP,GMASS,MDCYC,TEMITT,EMITT,TMASS,MDUMP,MCOOL,
     1SCHW(1),SCHW(2),SCHW(3),TSOLL,te,iav
      DO 930 I=1,NAA
      WRITE (IRES) VELO(1,I),VELO(2,I),VELO(3,I),ACC(1,I),ACC(2,I),
     1ACC(3,I),X3(1,I),X3(2,I),X3(3,I)
  930 CONTINUE
C
C END RESTART FILE
C
C
C    CALCULATE DIPOLE MOMENT
C
  550 SCS=0.
      DMX=0.0
      DMY=0.0
      DMZ=0.0
      DO 535 I=1,NAA
         DMX=DMX+(CS(I)+QSIG(I))*X(1,I)
         DMY=DMY+(CS(I)+QSIG(I))*X(2,I)
         DMZ=DMZ+(CS(I)+QSIG(I))*X(3,I)
  535 SCS=SCS+(CS(I)+QSIG(I))
      IF (ABS(SCS).GT.0.1) THEN
      WRITE(IUT,790)SCS
      GOTO 2000
      END IF
      DMX=4.8023*DMX
      DMY=4.8023*DMY
      DMZ=4.8023*DMZ
      DM=SQRT(DMX**2+DMY**2+DMZ**2)
      WRITE(IUT,780) DM,DMX,DMY,DMZ
C
C    START NEXT CALCULATION
C
      GOTO 2000
C
  545 WRITE(IUT,750) DHFM
 2000 if (nf.eq.1) call freqan
      go to 20
C
  560 FORMAT (I3,I1,19A4/20A4)
  570 FORMAT (4I3,F6.2,12I3)
  580 FORMAT (3I3,1X,3F10.4)
  590 FORMAT (4I3,F8.2)
  600 FORMAT (10X,34HPI-SCF-MOLECULAR-MECHANICS PROGRAM/10X,19A4
     1/10X,20A4/10X,5HISEQ=,I5)
  601 FORMAT (10X,'INPUT FORMAT INDICATOR INFORM IS', I2)
  610 FORMAT (12X,40H   NA   NB   NC   ND   ZB   NG   NH   NI,
     1/12x,4i5,f5.2,3i5)
  611 format (12x,40H   NK   NL   NM IGEO   NE NEPS  IHB  NF ,
     1/12X,8I5)
  620 FORMAT (10X,23HCOORDINATES (ANGSTROMS)/15X,1HI,7X,6HATOMIC
     1,4X,4HIMOL,3X,4HIFIX,5X,1HX,9X,1HY,9X,1HZ/1H ,23X,2HNR,2X,4HSORT/)
  630 FORMAT (13X,A2,I3,I7,3I6,3F10.4)
  640 FORMAT (10X,'NEW INTERNAL COORDINATE ',4I5,' =',F8.2)
  650 FORMAT (1H ,3(F8.2),F10.1)
  660 FORMAT (10X,33HUNABLE TO REACH SELFCONSISTENCY  )
  670 FORMAT (10X,10HSCF CYCLES,I5)
  680 FORMAT (1X,41H FINAL BOND ORDER - CHARGE DENSITY MATRIX)
  690 FORMAT (5(2I4,F9.4,7X))
  700 FORMAT (2X,28H FINAL EIGENVALUES ( IN EV ))
  710 FORMAT (6(I3,F12.6,4X))
  720 FORMAT (2X,19H FINAL EIGENVECTORS)
  730 FORMAT (10X,21HINTERATOMIC DISTANCES)
C  740 FORMAT ('  BOND ENERGY PI-BOND ENERGY STRAIN ENERGY    HEAT OF',
C     1'   NUMBER NO CONVERGENCE   AVERAGE OF    AVERAGE OF     ENERGY OF
C     2'/46X,'FORMATION    OF      (EP-EPA)    COORD.CHANGES',
C     3'   GRADIENTS    DEFORMATION'/8X,'EV',12X,'EV',12X,'EV',9X,
C     4'KJ/MOLE   CYCLES ',6X,'EV',9X,'ANGSTROM     EV/ANGSTROM',8X,
C     5'EV'/)
  740 FORMAT('  BOND  PI-BOND   STRAIN   HEAT OF  NUMBER  AVERAGE  AVERA
     1GE  ENERGY',/,' ENERGY  ENERGY   ENERGY  FORMATION   OF     COORD.
     2  GRADIENT   OF',/,'  (EV)    (EV)     (EV)   (KJ/MOL)  CYCLES  CH
     3ANGES  (EV/A)   DEFORM.')
  750 FORMAT(50X,'FINAL HEAT OF FORMATION',F8.1 ,' KJ/MOL')
  760 FORMAT(' CALCULATION NO.',I4)
  770 FORMAT(1H ,'CALCULATION NO.',I4)
  780 FORMAT(10X,'DIPOLE MOMENT',F9.3,' WITH COMPONENTS',3F8.3)
  790 FORMAT(1H ,'NET CHARGE ',F7.2)
C
      END
      SUBROUTINE PARSUB
C
C    ENDE DER PARAMETRISIERUNG:  OKTOBER 1991
C
C    SET ALL PARAMETERS USED FOR CALCULATION
C
      CHARACTER*2 ATA
      COMMON /LABEL/ ATA(50)
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      COMMON /POINT/ XPP(3,288)
      DIMENSION XP1(3,48),XP2(3,48),XP3(3,48),XP4(3,48),XP5(3,48),
     1XP6(3,48)
      DATA AMASS /12.,12.,14.,14.,16.,16.,32.,35.,19.,81.,12.,14.,
     124.,40.,28.,31.,32.,7.,39.,1.,63.,112.,56.,59.,59.,65.,140.,
     2115.,140.,91.,232.,59.,56.,59.,52.,12.,4*0.,11.,9*0./
      DATA XN/7.98,8.79,11.54,12.87,14.18,17.07,10.14,11.,14.66,
     110.08,8.79,15.68,1.,1.,7.3,8.79,10.14,1.,1.,7.17,15*1.,8.79,
     24*1.,5.98,9*1./
      DATA B/9.18,9.32,10.82,11.15,12.92,13.79,9.13,9.69,13.85,
     18.47,9.32,11.70,0.,0.,6.57,9.32,9.13,0.,0.,6.24,15*0.,9.32,
     24*0.,6.82,9*0./
      DATA C/1.88,1.51,1.36,0.85,1.39,0.47,1.38,1.35,2.31,
     11.16,1.51,-0.27,0.,0.,0.66,0.51,1.38,0.,0.,-0.56,15*0.,1.51,
     24*0.,1.61,9*0./
      DATA XP /19.04,19.62,23.72,24.87,28.94,31.33,20.65,22.04,30.82,
     119.71,19.62,27.11,0.,0.,14.53,18.73,20.65,0.,0.,15.00,15*0.,
     219.62,4*0.,14.41,9*0./
C
C
      DATA XP1/
     +  0.199, 0.000, 0.981, 0.101, 0.172, 0.981,-0.099, 0.173, 0.981,
     + -0.199, 0.000, 0.981,-0.101,-0.172, 0.981, 0.099,-0.173, 0.981,
     +  0.400, 0.000, 0.918, 0.347, 0.199, 0.918, 0.203, 0.345, 0.918,
     +  0.000, 0.400, 0.918,-0.199, 0.347, 0.918,-0.345, 0.203, 0.918,
     + -0.400, 0.000, 0.918,-0.347,-0.199, 0.918,-0.203,-0.345, 0.918,
     +  0.000,-0.400, 0.918, 0.199,-0.347, 0.918, 0.345,-0.203, 0.918,
     +  0.582,-0.001, 0.815, 0.549, 0.195, 0.815, 0.449, 0.371, 0.815,
     +  0.295, 0.502, 0.815, 0.105, 0.573, 0.815,-0.099, 0.574, 0.815,
     + -0.290, 0.505, 0.815,-0.445, 0.376, 0.815,-0.547, 0.200, 0.815,
     + -0.582, 0.001, 0.815,-0.547,-0.200, 0.815,-0.449,-0.371, 0.815,
     + -0.295,-0.502, 0.815,-0.105,-0.573, 0.815, 0.099,-0.574, 0.815,
     +  0.290,-0.505, 0.815, 0.445,-0.376, 0.815, 0.547,-0.200, 0.815,
     +  0.739,-0.001, 0.675, 0.711, 0.205, 0.675, 0.622, 0.400, 0.675,
     +  0.488, 0.556, 0.675, 0.309, 0.672, 0.675, 0.111, 0.731, 0.675,
     + -0.104, 0.732, 0.675,-0.302, 0.675, 0.675,-0.482, 0.561, 0.675,
     + -0.622, 0.400, 0.675,-0.708, 0.212, 0.675,-0.739, 0.001, 0.675/
       DATA XP2/
     + -0.711,-0.205, 0.675,-0.622,-0.400, 0.675,-0.488,-0.556, 0.675,
     + -0.309,-0.672, 0.675,-0.111,-0.731, 0.675, 0.104,-0.732, 0.675,
     +  0.309,-0.672, 0.675, 0.482,-0.561, 0.675, 0.622,-0.400, 0.675,
     +  0.708,-0.212, 0.675, 0.864,-0.001, 0.507, 0.839, 0.206, 0.507,
     +  0.766, 0.400, 0.507, 0.649, 0.570, 0.507, 0.495, 0.708, 0.507,
     +  0.313, 0.805, 0.507, 0.112, 0.856, 0.507,-0.104, 0.858, 0.507,
     + -0.305, 0.809, 0.507,-0.488, 0.713, 0.507,-0.644, 0.577, 0.507,
     + -0.762, 0.407, 0.507,-0.837, 0.214, 0.507,-0.864, 0.001, 0.507,
     + -0.839,-0.206, 0.507,-0.766,-0.400, 0.507,-0.649,-0.570, 0.507,
     + -0.495,-0.708, 0.507,-0.313,-0.805, 0.507,-0.104,-0.858, 0.507,
     +  0.104,-0.858, 0.507, 0.305,-0.809, 0.507, 0.488,-0.713, 0.507,
     +  0.644,-0.577, 0.507, 0.762,-0.407, 0.507, 0.839,-0.206, 0.507,
     +  0.950,-0.001, 0.316, 0.929, 0.199, 0.316, 0.864, 0.397, 0.316,
     +  0.762, 0.568, 0.316, 0.620, 0.721, 0.316, 0.448, 0.839, 0.316,
     +  0.263, 0.913, 0.316, 0.058, 0.949, 0.316,-0.152, 0.938, 0.316,
     + -0.344, 0.885, 0.316,-0.529, 0.789, 0.316,-0.688, 0.655, 0.316/
       DATA XP3 /
     + -0.810, 0.497, 0.316,-0.899, 0.309, 0.316,-0.945, 0.105, 0.316,
     + -0.946,-0.095, 0.316,-0.902,-0.300, 0.316,-0.815,-0.489, 0.316,
     + -0.695,-0.648, 0.316,-0.537,-0.784, 0.316,-0.353,-0.882, 0.316,
     + -0.162,-0.936, 0.316, 0.049,-0.949, 0.316, 0.254,-0.916, 0.316,
     +  0.440,-0.843, 0.316, 0.613,-0.727, 0.316, 0.756,-0.576, 0.316,
     +  0.859,-0.406, 0.316, 0.927,-0.208, 0.316, 0.995,-0.001, 0.111,
     +  0.975, 0.199, 0.111, 0.912, 0.397, 0.111, 0.810, 0.579, 0.111,
     +  0.671, 0.735, 0.111, 0.504, 0.859, 0.111, 0.313, 0.945, 0.111,
     +  0.110, 0.989, 0.111,-0.100, 0.990, 0.111,-0.304, 0.948, 0.111,
     + -0.495, 0.863, 0.111,-0.664, 0.742, 0.111,-0.804, 0.587, 0.111,
     + -0.909, 0.407, 0.111,-0.973, 0.208, 0.111,-0.995, 0.001, 0.111,
     + -0.973,-0.208, 0.111,-0.909,-0.407, 0.111,-0.810,-0.579, 0.111,
     + -0.671,-0.735, 0.111,-0.504,-0.859, 0.111,-0.313,-0.945, 0.111,
     + -0.110,-0.989, 0.111, 0.100,-0.990, 0.111, 0.304,-0.948, 0.111,
     +  0.495,-0.863, 0.111, 0.664,-0.742, 0.111, 0.804,-0.587, 0.111,
     +  0.909,-0.407, 0.111, 0.973,-0.208, 0.111, 0.996,-0.001,-0.101/
       DATA XP4 /
     +  0.976, 0.199,-0.101, 0.913, 0.398,-0.101, 0.811, 0.579,-0.101,
     +  0.672, 0.735,-0.101, 0.504, 0.859,-0.101, 0.314, 0.946,-0.101,
     +  0.110, 0.990,-0.101,-0.100, 0.991,-0.101,-0.305, 0.949,-0.101,
     + -0.496, 0.865,-0.101,-0.665, 0.742,-0.101,-0.805, 0.588,-0.101,
     + -0.909, 0.407,-0.101,-0.974, 0.209,-0.101,-0.996, 0.001,-0.101,
     + -0.974,-0.209,-0.101,-0.909,-0.407,-0.101,-0.811,-0.579,-0.101,
     + -0.672,-0.735,-0.101,-0.504,-0.859,-0.101,-0.314,-0.946,-0.101,
     + -0.110,-0.990,-0.101, 0.100,-0.991,-0.101, 0.305,-0.949,-0.101,
     +  0.496,-0.865,-0.101, 0.665,-0.742,-0.101, 0.805,-0.588,-0.101,
     +  0.909,-0.407,-0.101, 0.974,-0.209,-0.101, 0.953,-0.001,-0.306,
     +  0.932, 0.199,-0.306, 0.866, 0.398,-0.306, 0.765, 0.570,-0.306,
     +  0.622, 0.723,-0.306, 0.450, 0.841,-0.306, 0.264, 0.916,-0.306,
     +  0.058, 0.952,-0.306,-0.152, 0.941,-0.306,-0.346, 0.889,-0.306,
     + -0.531, 0.792,-0.306,-0.691, 0.658,-0.306,-0.812, 0.499,-0.306,
     + -0.902, 0.310,-0.306,-0.948, 0.105,-0.306,-0.949,-0.096,-0.306,
     + -0.905,-0.300,-0.306,-0.818,-0.490,-0.306,-0.698,-0.651,-0.306/
       DATA XP5 /
     + -0.539,-0.787,-0.306,-0.354,-0.885,-0.306,-0.162,-0.939,-0.306,
     +  0.049,-0.952,-0.306, 0.255,-0.919,-0.306, 0.441,-0.845,-0.306,
     +  0.614,-0.729,-0.306, 0.759,-0.577,-0.306, 0.862,-0.407,-0.306,
     +  0.930,-0.209,-0.306, 0.869,-0.001,-0.498, 0.844, 0.207,-0.498,
     +  0.771, 0.402,-0.498, 0.653, 0.574,-0.498, 0.498, 0.712,-0.498,
     +  0.315, 0.810,-0.498, 0.113, 0.862,-0.498,-0.105, 0.862,-0.498,
     + -0.306, 0.813,-0.498,-0.491, 0.717,-0.498,-0.648, 0.580,-0.498,
     + -0.766, 0.410,-0.498,-0.842, 0.216,-0.498,-0.869, 0.001,-0.498,
     + -0.844,-0.207,-0.498,-0.771,-0.402,-0.498,-0.653,-0.574,-0.498,
     + -0.498,-0.712,-0.498,-0.315,-0.810,-0.498,-0.105,-0.862,-0.498,
     +  0.105,-0.862,-0.498, 0.306,-0.813,-0.498, 0.491,-0.717,-0.498,
     +  0.648,-0.580,-0.498, 0.766,-0.410,-0.498, 0.844,-0.207,-0.498,
     +  0.746,-0.001,-0.668, 0.719, 0.199,-0.668, 0.640, 0.384,-0.668,
     +  0.514, 0.541,-0.668, 0.345, 0.661,-0.668, 0.156, 0.729,-0.668,
     + -0.045, 0.745,-0.668,-0.249, 0.703,-0.668,-0.428, 0.611,-0.668,
     + -0.575, 0.475,-0.668,-0.684, 0.298,-0.668,-0.738, 0.105,-0.668/
       DATA XP6 /
     + -0.740,-0.097,-0.668,-0.684,-0.298,-0.668,-0.580,-0.470,-0.668,
     + -0.434,-0.607,-0.668,-0.249,-0.703,-0.668,-0.053,-0.744,-0.668,
     +  0.149,-0.731,-0.668, 0.345,-0.661,-0.668, 0.509,-0.546,-0.668,
     +  0.636,-0.390,-0.668, 0.717,-0.206,-0.668, 0.591,-0.001,-0.809,
     +  0.557, 0.197,-0.809, 0.455, 0.377,-0.809, 0.299, 0.510,-0.809,
     +  0.106, 0.581,-0.809,-0.100, 0.582,-0.809,-0.293, 0.512,-0.809,
     + -0.451, 0.380,-0.809,-0.554, 0.203,-0.809,-0.591, 0.001,-0.809,
     + -0.554,-0.203,-0.809,-0.455,-0.377,-0.809,-0.299,-0.510,-0.809,
     + -0.106,-0.581,-0.809, 0.100,-0.582,-0.809, 0.293,-0.512,-0.809,
     +  0.451,-0.380,-0.809, 0.554,-0.203,-0.809, 0.409, 0.000,-0.914,
     +  0.355, 0.203,-0.914, 0.207, 0.353,-0.914, 0.000, 0.409,-0.914,
     + -0.203, 0.355,-0.914,-0.353, 0.207,-0.914,-0.409, 0.000,-0.914,
     + -0.355,-0.203,-0.914,-0.207,-0.353,-0.914, 0.000,-0.409,-0.914,
     +  0.203,-0.355,-0.914, 0.353,-0.207,-0.914, 0.209, 0.000,-0.979,
     +  0.106, 0.181,-0.979,-0.104, 0.182,-0.979,-0.209, 0.000,-0.979,
     + -0.106,-0.181,-0.979, 0.104,-0.182,-0.979, 0.001, 0.000,-1.001/
C
      DO 1001 I=1,48
      II=I
      XPP(1,II)=XP1(3,I)
      XPP(2,II)=XP1(1,I)
      XPP(3,II)=XP1(2,I)
 1001 CONTINUE
      DO 1002 I=1,48
      II=I+48
      XPP(1,II)=XP2(3,I)
      XPP(2,II)=XP2(1,I)
      XPP(3,II)=XP2(2,I)
 1002 CONTINUE
      DO 1003 I=1,48
      II=I+96
      XPP(1,II)=XP3(3,I)
      XPP(2,II)=XP3(1,I)
      XPP(3,II)=XP3(2,I)
 1003 CONTINUE
      DO 1004 I=1,48
      II=I+144
      XPP(1,II)=XP4(3,I)
      XPP(2,II)=XP4(1,I)
      XPP(3,II)=XP4(2,I)
 1004 CONTINUE
      DO 1005 I=1,48
      II=I+192
      XPP(1,II)=XP5(3,I)
      XPP(2,II)=XP5(1,I)
      XPP(3,II)=XP5(2,I)
 1005 CONTINUE
      DO 1006 I=1,48
      II=I+240
      XPP(1,II)=XP6(3,I)
      XPP(2,II)=XP6(1,I)
      XPP(3,II)=XP6(2,I)
 1006 CONTINUE
C
      DO 1 I=1,50
      RV(I)=0.
      RBOND(I)=0.
      DHFA(I)=0.
      DBS(I)=0.
      GANG(I)=0.
      DSW(I)=0.
      DECUB(I)=0.
      QQ(I)=0.
      ISPI(I)=0.
      ATA(I)='  '
      HV(I)=0.
      ZS(I)=0.
      ZV(I)=0.
      DO 1 J=1,50
      BL(I,J)=0.
      BK(I,J)=0.
      PTT(I,J)=0.
      AL(I,J)=0.
      AME(I,J)=0.
      DE(I,J)=0.
      DF(I,J)=0.
      PK2(I,J)=0.
      PK4(I,J)=0.
      PK6(I,J)=0.
      PPIT(I,J)=0.
      PPIZ(I,J)=0.
      ppix(i,j)=0.
      AVDW(I,J)=0.
      BVDW(I,J)=0.
      CVDW(I,J)=0.
      IVDW(I,J)=0
    1 CONTINUE
CX      DO 2 I=1,11
CX    2 CONTINUE
C
C    PARAMETERS
C    XP,B,C,XN PARAMETERS OF THE GASTEIGER SIGMA CHARGE ALGORITHM
C    AMASS ATOMIC MASSES
C    ZV  NUMBER OF PI-ELECTRONS CONTRIBUTED BY ATOM MU
C    HV  SCF-PARAMETER H(MU,MU)
C    RV  SCF-PARAMETER R(MU,MU)
C    ZS  EFFECTIVE Z OF ATOM MU
C    BK  SCF-PARAMETER TO CALCULATE RESONANCE INTEGRAL BETA(MU,NU)
C    DE  BOND ENERGY OF THE SINGLE BOND MU,NU
C    DF  BOND ENERGY CORRECTION FOR CHARGE EFFECTS IN POLAR BONDS
C    QQ  STANDARD SIGMA-CHARGE OF AN ATOM IN THE PI-SYSTEM
C    AME MORSE CONSTANT OF THE SINGLE BOND MU,NU
C    DHFA HEATS OF FORMATION FOR FREE ATOMS
C    RBOND MAXIMAL BOND RADII
C    PARAMETERS AL AND BL TO CALCULATE PI-BOND LENGTHS
C    GANG  EQUILIBRIUM BOND ANGLES
C    DSW   HARMONIC FORCE CONSTANTS FOR BOND ANGLE DISTORTION
C    DBS   HARMONIC FORCE CONSTANTS FOR BENDS OF (SP2)-CENTRES
C    PK2, PK4, PK6  PARAMETERS FOR BOND LENGTH - FORCE CONSTANT
C          RELATION
C    DECUB ADDITIONAL CONSTANT FOR BOND ANGLE POTENTIAL
C    PPIZ  CONSTANTS FOR TORSION POTENTIALS OF SINGLE BONDS
C    PPIT CONSTANTS FOR ADDITIONAL TORSIONAL POTENTIAL OF PI BONDS
C    AVDW, BVDW, CVDW, IVDW.. PARAMETERS FOR VAN DER WAALS
C    INTERACTIONS (E.GIGLIO, NATURE 222, 339 (1969))
C
      ATA(1)=' C'
      ATA(2)=' C'
      ATA(3)=' N'
      ATA(4)=' N'
      ATA(5)=' O'
      ATA(6)=' O'
      ATA(7)=' S'
      ATA(8)='CL'
      ATA(9)=' F'
      ATA(10)='BR'
      ATA(11)=' C'
      ATA(12)=' N'
      ATA(13)='MG'
      ATA(14)='CA'
      ATA(15)='SI'
      ATA(16)=' P'
      ATA(17)=' S'
      ATA(18)='LI'
      ATA(19)=' K'
      ATA(20)=' H'
      ATA(21)='CU'
      ATA(22)='CD'
      ATA(23)='FE'
      ATA(24)='CO'
      ATA(25)='NI'
      ATA(26)='ZN'
      ATA(27)='CE'
      ATA(28)='IN'
      ATA(29)='CE'
      ATA(30)='ZR'
      ATA(31)='TH'
      ATA(32)='NI'
      ATA(33)='CO'
      ATA(34)='CR'
      ATA(35)='FE'
      ATA(36)=' C'
      ATA(41)=' B'
      ZV(1)=1.
      ZV(2)=2.
      ZV(3)=1.
      ZV(4)=2.
      ZV(5)=1.
      ZV(6)=2.
      ZV(7)=2.
      ZV(8)=2.
      ZV(9)=2.
      ZV(10)=1.
      ZV(11)=1.
      ZV(35)=1.
      ISPI(2)=1
      ISPI(3)=1
      ISPI(4)=1
      ISPI(5)=1
      ISPI(6)=1
      ISPI(7)=1
      ISPI(8)=1
      ISPI(9)=1
      ISPI(10)=1
      ISPI(11)=1
      ISPI(12)=1
      ISPI(36)=1
      HV(1)=11.16
      HV(2)=28.586
      HV(3)=14.12
      HV(4)=33.901
      HV(5)=17.697
      HV(6)=22.88
      HV(7)=25.32
      HV(8)=39.82
      HV(9)=25.32
      HV(10)=11.19
      HV(11)=14.18
      HV(35)=11.19
      RV(1)=12.27
      RV(2)=11.13
      RV(3)=16.628
      RV(4)=12.341
      RV(5)=18.600
      RV(6)=15.227
      RV(7)=11.90
      RV(8)=10.81
      RV(9)=21.025
      RV(10)=10.
      RV(11)=11.09
      RV(12)=12.52
      RV(13)=12.85
      RV(14)=12.85
      RV(15)=9.04
      RV(16)=9.73
      RV(17)=15.
      RV(18)=12.85
      RV(19)=12.85
      RV(20)=12.85
      RV(21)=12.85
      RV(22)=12.85
      RV(23)=12.85*0.5
      RV(24)=12.85
      RV(25)=12.85
      RV(26)=12.85
      RV(27)=12.85
      RV(28)=12.85
      RV(29)=12.85
      RV(30)=12.85
      RV(31)=12.85
      RV(32)=12.85
      RV(33)=12.85
      RV(34)=12.85
      RV(35)=12.85
      RV(36)=12.27
      RV(41)=11.96
      ZS(1)=3.18
      ZS(2)=4.92
      ZS(3)=3.55
      ZS(4)=5.43
      ZS(5)=4.18
      ZS(6)=4.55
      ZS(7)=4.50
      ZS(8)=6.01
      ZS(9)=4.50
      ZS(10)=3.25
      ZS(11)=3.75
      ZS(35)=3.25
      BL(1,1)=0.174
      BL(1,2)=0.178
      BL(1,3)=0.178
      BL(1,4)=0.185
      BL(1,5)=0.185
      BL(1,6)=0.229
      BL(1,7)=0.230
      BL(1,8)=0.230
      BL(1,9)=0.230
      BL(1,10)=0.174
      BL(1,11)=0.000
      BL(2,2)=0.177
      BL(2,3)=0.177
      BL(2,4)=0.190
      BL(2,5)=0.190
      BL(2,6)=0.000
      BL(2,7)=0.000
      BL(2,8)=0.000
      BL(2,9)=0.000
      BL(2,10)=0.178
      BL(2,11)=0.000
      BL(3,3)=0.177
      BL(3,4)=0.190
      BL(3,5)=0.190
      BL(3,6)=0.000
      BL(3,7)=0.000
      BL(3,8)=0.000
      BL(3,9)=0.000
      BL(3,10)=0.178
      BL(3,11)=0.000
      BL(3,41)=0.01
      BL(4,4)=0.000
      BL(4,5)=0.000
      BL(4,6)=0.000
      BL(4,7)=0.000
      BL(4,8)=0.000
      BL(4,9)=0.000
      BL(4,10)=0.000
      BL(4,11)=0.000
      BL(5,5)=0.000
      BL(5,6)=0.000
      BL(5,7)=0.000
      BL(5,8)=0.000
      BL(5,9)=0.000
      BL(5,11)=0.000
      BL(5,35)=0.2
      BL(6,6)=0.230
      BL(6,7)=0.000
      BL(6,8)=0.000
      BL(6,9)=0.000
      BL(6,10)=0.000
      BL(6,11)=0.000
      BL(7,7)=0.000
      BL(7,8)=0.000
      BL(7,9)=0.000
      BL(7,10)=0.000
      BL(7,11)=0.000
      BL(8,8)=0.000
      BL(8,9)=0.000
      BL(8,10)=0.000
      BL(8,11)=0.000
      BL(9,9)=0.000
      BL(9,10)=0.000
      BL(9,11)=0.000
      BL(10,10)=0.220
      BL(10,11)=0.205
      BL(11,11)=0.000
      AL(1,1)=1.533
      AL(1,2)=1.508
      AL(1,3)=1.467
      AL(1,4)=1.472
      AL(1,5)=1.420
      AL(1,6)=1.426
      AL(1,7)=1.820
      AL(1,8)=1.772
      AL(1,9)=1.330
      AL(1,10)=1.933
      AL(1,11)=1.462
      AL(1,12)=0.000
      AL(1,13)=0.000
      AL(1,14)=0.000
      AL(1,15)=1.850
      AL(1,16)=1.830
      AL(1,17)=1.807
      AL(1,18)=0.000
      AL(1,19)=0.000
      AL(1,20)=1.1
      AL(1,41)=1.6
      AL(2,2)=1.512
      AL(2,3)=1.442
      AL(2,4)=1.454
      AL(2,5)=1.424
      AL(2,6)=1.398
      AL(2,7)=1.807
      AL(2,8)=1.755
      AL(2,9)=1.380
      AL(2,10)=1.892
      AL(2,11)=1.380
      AL(2,12)=1.467
      AL(2,13)=0.000
      AL(2,14)=0.000
      AL(2,15)=1.872
      AL(2,16)=1.797
      AL(2,17)=1.787
      AL(2,18)=0.000
      AL(2,19)=0.000
      AL(2,20)=1.08
      AL(3,3)=1.446
      AL(3,4)=1.417
      AL(3,5)=1.405
      AL(3,6)=1.380
      AL(3,7)=0.000
      AL(3,8)=0.000
      AL(3,9)=0.000
      AL(3,10)=0.000
      AL(3,11)=0.000
      AL(3,12)=0.000
      AL(3,13)=2.200
      AL(3,14)=2.750
      AL(3,15)=1.750
      AL(3,16)=1.600
      AL(3,17)=1.650
      AL(3,18)=0.000
      AL(3,19)=0.000
      AL(3,20)=1.03
      AL(3,41)=1.54
      AL(4,4)=1.430
      AL(4,5)=1.436
      AL(4,6)=1.401
      AL(4,7)=0.000
      AL(4,8)=0.000
      AL(4,9)=0.000
      AL(4,10)=0.000
      AL(4,11)=0.000
      AL(4,12)=0.000
      AL(4,13)=0.000
      AL(4,14)=0.000
      AL(4,15)=0.000
      AL(4,16)=0.000
      AL(4,17)=0.000
      AL(4,18)=0.000
      AL(4,19)=0.000
      AL(4,20)=1.01
      AL(5,5)=0.000
      AL(5,6)=0.000
      AL(5,7)=0.000
      AL(5,8)=0.000
      AL(5,9)=0.000
      AL(5,10)=0.000
      AL(5,11)=0.000
      AL(5,12)=0.000
      AL(5,13)=0.000
      AL(5,14)=0.000
      AL(5,15)=1.650
      AL(5,16)=1.60
      AL(5,17)=1.592
      AL(5,18)=0.000
      AL(5,19)=0.000
      AL(5,20)=0.972
      AL(6,6)=0.000
      AL(6,7)=0.000
      AL(6,8)=0.000
      AL(6,9)=0.000
      AL(6,10)=0.000
      AL(6,11)=0.00
      AL(6,12)=0.000
      AL(6,13)=0.000
      AL(6,14)=0.000
      AL(6,15)=0.000
      AL(6,16)=0.000
      AL(6,17)=0.000
      AL(6,18)=0.000
      AL(6,19)=0.000
      AL(6,20)=1.00
      AL(6,36)=1.33
      AL(7,7)=2.040
      AL(7,8)=0.000
      AL(7,9)=0.000
      AL(7,10)=0.000
      AL(7,12)=0.000
      AL(7,13)=0.000
      AL(7,14)=0.000
      AL(7,15)=0.000
      AL(7,16)=0.000
      AL(7,17)=2.040
      AL(7,18)=0.000
      AL(7,19)=0.000
      AL(7,20)=1.335
      AL(8,8)=0.000
      AL(8,9)=0.000
      AL(8,10)=0.000
      AL(8,11)=0.000
      AL(8,12)=0.000
      AL(8,13)=0.000
      AL(8,14)=0.000
      AL(8,15)=0.000
      AL(8,16)=0.000
      AL(8,17)=0.000
      AL(8,18)=0.000
      AL(8,19)=0.000
      AL(8,20)=0.000
      AL(9,9)=0.000
      AL(9,10)=0.000
      AL(9,11)=0.000
      AL(9,12)=0.000
      AL(9,13)=0.000
      AL(9,14)=0.000
      AL(9,15)=0.000
      AL(9,16)=0.000
      AL(9,17)=0.000
      AL(9,18)=0.000
      AL(9,19)=0.000
      AL(9,20)=0.000
      AL(10,10)=0.000
      AL(10,11)=0.000
      AL(10,12)=0.000
      AL(10,13)=0.000
      AL(10,14)=0.000
      AL(10,15)=0.000
      AL(10,16)=0.000
      AL(10,17)=0.000
      AL(10,18)=0.000
      AL(10,19)=0.000
      AL(10,20)=0.000
      AL(11,11)=1.425
      AL(11,12)=1.360
      AL(11,13)=1.370
      AL(11,14)=0.000
      AL(11,15)=1.79
      AL(11,16)=0.000
      AL(11,17)=0.000
      AL(11,18)=0.000
      AL(11,19)=0.000
      AL(11,20)=1.06
      AL(15,15)=2.331
      AL(15,20)=1.484
      DHFA(1)=170.89
      DHFA(2)=170.89
      DHFA(3)=113.0
      DHFA(4)=113.0
      DHFA(5)=59.559
      DHFA(6)=59.559
      DHFA(7)=65.65
      DHFA(8)=28.95
      DHFA(9)=18.86
      DHFA(10)=26.73
      DHFA(11)=170.89
      DHFA(12)=113.0
      DHFA(13)=502.
      DHFA(14)=460.29
      DHFA(15)=88.04
      DHFA(16)=75.18
      DHFA(17)=65.65
      DHFA(18)=162.42
      DHFA(19)=122.92
      DHFA(20)=52.102
      DHFA(36)=170.89
      RBOND(1)=0.8
      RBOND(2)=0.8
      RBOND(3)=0.8
      RBOND(4)=0.8
      RBOND(5)=0.8
      RBOND(6)=0.8
      RBOND(7)=1.1
      RBOND(8)=1.2
      RBOND(9)=0.7
      RBOND(10)=1.2
      RBOND(11)=0.8
      RBOND(12)=0.6
      RBOND(13)=0.0
      RBOND(14)=0.0
      RBOND(15)=1.2
      RBOND(16)=1.0
      RBOND(17)=1.1
      RBOND(18)=0.0
      RBOND(19)=0.0
      RBOND(20)=0.4
      RBOND(21)=0.4
      RBOND(22)=0.4
      RBOND(23)=0.4
      RBOND(24)=0.4
      RBOND(25)=0.4
      RBOND(26)=0.4
      RBOND(27)=0.4
      RBOND(28)=0.4
      RBOND(29)=0.4
      RBOND(30)=0.4
      RBOND(31)=0.4
      RBOND(33)=0.1
      RBOND(34)=0.1
      RBOND(35)=0.1
      RBOND(36)=0.6
      RBOND(41)=1.0
      QQ(1)=-0.0
      QQ(2)=-0.22
      QQ(3)=-0.235
      QQ(4)=-0.08
      QQ(5)=-0.24
      QQ(6)=-0.058
      QQ(7)=-0.158
      QQ(8)=-0.0900
      QQ(9)=-0.1600
      QQ(10)=-0.0900
      QQ(11)=-0.20
      QQ(12)=-0.270
      QQ(13)=0.0922
      QQ(14)=0.0922
      QQ(15)=-0.1224
      QQ(16)=-0.2000
      QQ(17)=-0.2200
      QQ(18)=0.0922
      QQ(19)=0.0922
      QQ(20)=0.0
      QQ(36)=0.24
      BK(1,1)=7.743
      BK(1,2)=15.600
      BK(1,3)=7.200
      BK(1,4)=16.000
      BK(1,5)=9.622
      BK(1,6)=18.000
      BK(1,7)=18.50
      BK(1,8)=15.00
      BK(1,9)=18.5
      BK(1,10)=6.927
      BK(1,11)=0.00
      BK(2,2)=17.876
      BK(2,3)=14.559
      BK(2,4)=23.698
      BK(2,5)=20.20
      BK(2,6)=0.00
      BK(2,7)=0.00
      BK(2,8)=0.00
      BK(2,9)=0.00
      BK(2,10)=17.876
      BK(2,11)=0.00
      BK(3,3)=6.820
      BK(3,4)=14.371
      BK(3,5)=9.66
      BK(3,6)=0.00
      BK(3,7)=0.00
      BK(3,8)=0.00
      BK(3,9)=0.00
      BK(3,10)=0.00
      BK(3,11)=0.00
      BK(4,4)=0.00
      BK(4,5)=0.00
      BK(4,6)=0.00
      BK(4,7)=0.00
      BK(4,8)=0.00
      BK(4,9)=0.00
      BK(4,10)=0.00
      BK(4,11)=0.00
      BK(5,5)=0.00
      BK(5,6)=0.00
      BK(5,7)=0.00
      BK(5,8)=0.00
      BK(5,9)=0.00
      BK(5,10)=0.00
      BK(5,35)=9.622
      BK(6,7)=0.00
      BK(6,8)=0.00
      BK(6,9)=0.00
      BK(6,10)=0.00
      BK(6,11)=0.00
      BK(7,7)=0.00
      BK(7,8)=0.00
      BK(7,9)=0.00
      BK(7,10)=0.00
      BK(7,11)=0.00
      BK(8,8)=0.00
      BK(8,9)=0.00
      BK(8,10)=0.00
      BK(8,11)=0.00
      BK(9,9)=0.00
      BK(9,10)=0.00
      BK(9,11)=0.00
      BK(10,10)=10.600
      BK(10,11)= 11.830
      BK(11,11)=0.000
      DF(1,1)=-0.5
      DF(1,2)=-0.5
      DF(1,3)=0.1
      DF(1,4)=0.1
      DF(1,5)=-0.5
      DF(1,7)=0.5
      DF(1,8)=0.55
      DF(1,9)=-0.85
      DF(1,10)=-0.4
      DF(1,20)=-0.5
      DF(2,2)=-0.25
      DF(2,3)=0.1
      DF(2,4)=0.1
      DF(2,5)=0.4
      DF(2,6)=0.600
      DF(2,7)=-0.500
      DF(2,8)=0.500
      DF(2,9)=0.3
      DF(2,10)=0.5
      DF(2,11)=0.700
      DF(2,20)=-0.5
      DF(3,3)=0.250
      DF(3,4)=0.250
      DF(3,5)=0.250
      DF(3,6)=-.7
      DF(3,13)=0.65
      DF(3,14)=0.65
      DF(3,20)=-0.3682
      DF(3,41)=0.1
      DF(4,4)=0.200
      DF(4,5)=-.75
      DF(4,6)=0.250
      DF(4,20)=-0.25
      DF(5,16)=-0.3
      DF(5,20)=-0.22
      DF(7,20)=-0.5
      DF(11,11)=0.5
      DF(11,12)=.8
      DE(1,1)=3.5900
      DE(1,2)=3.8800
      DE(1,3)=3.0450
      DE(1,4)=3.5288
      DE(1,5)=3.60
      DE(1,6)=3.5780
      DE(1,7)=2.9224
      DE(1,8)=3.41
      DE(1,9)=4.5
      DE(1,10)=2.88
      DE(1,11)=4.4
      DE(1,15)=3.00
      DE(1,16)=2.8000
      DE(1,17)=3.0824
      DE(1,20)=4.2681
      DE(2,2)=4.152
      DE(2,3)=3.0950
      DE(2,4)=3.4778
      DE(2,5)=4.1
      DE(2,6)=4.6
      DE(2,7)=3.5494
      DE(2,8)=3.7236
      DE(2,9)=5.
      DE(2,10)=3.1
      DE(2,11)=4.5
      DE(2,15)=4.20
      DE(2,16)=3.1
      DE(2,17)=3.45
      DE(2,20)=4.3891
      DE(3,3)=1.9168
      DE(3,4)=3.1660
      DE(3,5)=2.8300
      DE(3,6)=3.8
      DE(3,13)=5.
      DE(3,14)=3.
      DE(3,15)=2.9
      DE(3,16)=3.0000
      DE(3,17)=3.6300
      DE(3,20)=3.6070
      DE(3,41)=2.
      DE(4,4)=3.1900
      DE(4,5)=3.3500
      DE(4,6)=3.4092
      DE(4,20)=3.8128
      DE(5,15)=3.8244
      DE(5,16)=3.6500
      DE(5,17)=3.0000
      DE(5,20)=4.0
      DE(6,20)=4.0202
      DE(6,36)=1.07
      DE(7,7)=2.8
      DE(7,15)=2.3527
      DE(7,20)=3.5802
      DE(11,11)=4.625
      DE(11,12)=5.000
      DE(11,20)=4.8
      DE(15,15)=1.8345
      DE(15,20)=3.15
      AME(1,1)=1.8407
      AME(1,2)=1.8641
      AME(1,3)=1.2959
      AME(1,4)=1.2959
      AME(1,5)=1.0838
      AME(1,6)=1.0838
      AME(1,7)=1.9663
      AME(1,8)=1.8000
      AME(1,9)=1.8000
      AME(1,10)=1.8000
      AME(1,11)=1.8641
      AME(1,15)=2.92
      AME(1,16)=1.9000
      AME(1,17)=1.9663
      AME(1,20)=1.3800
      AME(2,2)=2.0022
      AME(2,3)=1.9209
      AME(2,4)=1.9209
      AME(2,5)=1.7870
      AME(2,6)=1.7870
      AME(2,7)=1.9663
      AME(2,8)=1.8000
      AME(2,9)=1.8000
      AME(2,10)=1.8000
      AME(2,11)=2.0022
      AME(2,15)=2.97
      AME(2,16)=1.9000
      AME(2,17)=1.9663
      AME(2,20)=1.3320
      AME(3,3)=2.5290
      AME(3,4)=2.5290
      AME(3,5)=2.0275
      AME(3,6)=2.0275
      AME(3,13)=0.1
      AME(3,14)=1.
      AME(3,16)=1.9000
      AME(3,17)=1.9663
      AME(3,20)=1.4690
      AME(3,41)=1.3
      AME(4,4)=2.5290
      AME(4,5)=2.0275
      AME(4,6)=2.0275
      AME(4,20)=1.4690
      AME(5,15)=2.9000
      AME(5,16)=1.9000
      AME(5,17)=1.9663
      AME(5,20)=1.2400
      AME(6,20)=1.2400
      AME(6,36)=1.0838
      AME(7,7)=2.
      AME(7,20)=1.5000
      AME(11,11)=2.3177
      AME(11,12)=2.3500
      AME(11,20)=1.3320
      AME(15,15)=2.95
      AME(15,20)=2.72
      GANG(1)=1.91986
      GANG(2)=2.09440
      GANG(3)=2.00713
      GANG(4)=2.09440
      GANG(5)=1.92
      GANG(6)=2.09440
      GANG(7)=1.83260
      GANG(11)=3.14159
      GANG(15)=2.0944
      GANG(16)=1.72788
      GANG(17)=1.72988
      GANG(36)=3.14159
      GANG(41)=2.00713
      DSW(1)=0.572E-3
      DSW(2)=0.442E-3
      dsw(3)=0.54e-3
      DSW(4)=0.798E-3
      DSW(5)=0.668E-3
      DSW(6)=0.752E-3
      DSW(7)=0.800E-3
      DSW(11)=0.221E-3
      DSW(15)=0.30E-3
      DSW(16)=0.280E-3
      DSW(17)=0.280E-3
      DSW(36)=0.221E-3
      DSW(41)=0.486E-3
      DBS(2)=0.411E-3
      pk2(1,1)=-24.
      pk2(1,2)=-27.
      pk2(1,3)=-25.
      PK2(1,4)=-22.292
      pk2(1,5)=-22.0
      PK2(1,6)=-20.067
      PK2(1,7)=-66.726
      PK2(1,8)=10.770
      PK2(1,9)=11.380
      PK2(1,10)=10.770
      PK2(1,11)=-31.514
      PK2(1,15)=10.
      PK2(1,16)=20.
      PK2(1,17)=-66.726
      pk2(1,20)=7.25
      pk2(2,2)=-33.5
      PK2(2,3)=-22.292
      PK2(2,4)=-22.292
      PK2(2,5)=-20.067
      pk2(2,6)=-17.5
      PK2(2,7)=-66.726
      PK2(2,8)=10.770
      PK2(2,9)=11.380
      PK2(2,10)=10.770
      PK2(2,11)=-31.514
      PK2(2,15)=10.
      PK2(2,16)=20.
      PK2(2,17)=-66.726
      PK2(2,20)=8.16
      PK2(3,3)=-56.632
      PK2(3,4)=-56.632
      PK2(3,5)=18.267
      pk2(3,6)=23.5
      PK2(3,11)=-22.292
      PK2(3,13)=7.5
      PK2(3,14)=10.
      PK2(3,15)=10.
      PK2(3,17)=-66.726
      PK2(3,16)=20.
      pk2(3,20)=8.5 
      PK2(3,41)=10.
      PK2(4,4)=-56.632
      PK2(4,5)=18.267
      PK2(4,6)=18.267
      PK2(4,11)=-22.292
      PK2(4,20)=7.14
      PK2(5,15)=20.
      PK2(5,16)=22.
      PK2(5,17)=-66.726
      PK2(5,11)=-20.067
      pk2(5,20)=8.20
      PK2(6,20)=7.00
      PK2(6,36)=-20.067
      PK2(7,7)=-55.826
      PK2(7,11)=-66.726
      pk2(7,20)=9.0
      pk2(11,11)=-21.5
      PK2(11,15)=10.
      pk2(11,12)=-16.0
      pk2(11,20)=8.0
      PK2(15,15)=10.0
      PK2(15,20)=9.00
      PK4(1,1)=130.072
      PK4(1,2)=130.072
      PK4(1,3)=67.003
      PK4(1,4)=67.003
      PK4(1,5)=52.681
      PK4(1,6)=52.681
      PK4(1,7)=261.520
      PK4(1,11)=130.072
      PK4(1,17)=261.520
      PK4(2,2)=130.072
      PK4(2,3)=67.003
      PK4(2,4)=67.003
      PK4(2,5)=52.681
      PK4(2,6)=52.681
      PK4(2,7)=261.520
      PK4(2,11)=130.072
      PK4(2,17)=261.520
      PK4(3,3)=179.565
      PK4(3,4)=179.565
      PK4(3,5)=-61.330
      PK4(3,6)=-61.330
      PK4(3,11)=67.003
      PK4(3,17)=261.520
      PK4(4,4)=179.565
      PK4(4,5)=-61.330
      PK4(4,6)=-61.330
      PK4(4,11)=67.003
      PK4(5,11)=52.681
      PK4(5,17)=261.520
      PK4(6,36)=52.681
      PK4(7,7)=408.020
      PK4(7,11)=261.520
      PK4(11,11)=130.072
      PK4(11,12)=67.003
      PK6(1,1)=-70.269
      PK6(1,2)=-70.269
      PK6(1,3)=-6.018
      PK6(1,4)=-6.018
      PK6(1,5)=3.629
      PK6(1,6)=3.629
      PK6(1,11)=-70.269
      PK6(2,2)=-70.269
      PK6(2,3)=-6.018
      PK6(2,4)=-6.018
      PK6(2,5)=3.629
      PK6(2,6)=3.629
      PK6(2,11)=-70.269
      PK6(3,3)=-94.043
      PK6(3,4)=-94.043
      PK6(3,5)=75.710
      PK6(3,6)=75.710
      PK6(3,11)=-6.018
      PK6(4,4)=-94.043
      PK6(4,5)=75.710
      PK6(4,6)=75.710
      PK6(4,11)=-6.018
      PK6(5,11)=3.629
      PK6(6,36)=3.629
      PK6(11,11)=-70.269
      PK6(11,12)=-6.018
      if (ng.lt.4) DECUB(1)=-0.1246E-2
      if (ng.lt.4) decub(5)=-0.125e-2
      ppix(2,20)=2.0e-5
      ppix(1,2)=3.0e-5
      ppix(6,20)=1.5e-5
      ppix(1,6)=2.0e-5
      ppix(5,20)=-0.5e-5
      ppix(1,5)=-0.5e-5
      PPIZ(1,1)=1.67E-5
      PPIZ(1,2)=1.2E-5
      PPIZ(1,3)=1.3E-5
      PPIZ(1,4)=1.3E-5
      PPIZ(1,5)=1.2E-5
      PPIZ(1,6)=1.2E-5
      PPIZ(1,7)=1.6E-5
      PPIZ(1,8)=1.6E-5
      PPIZ(1,9)=1.6E-5
      PPIZ(1,10)=1.6E-5
      PPIZ(1,11)=1.6E-5
      PPIZ(1,15)=2.1E-5
      PPIZ(1,16)=1.6E-5
      PPIZ(1,17)=1.6E-5
      PPIZ(1,20)=1.27E-5
      PPIZ(2,2)=1.6E-5
      PPIZ(2,3)=1.3E-5
      PPIZ(2,4)=1.3E-5
      PPIZ(2,5)=1.2E-5
      PPIZ(2,6)=1.2E-5
      PPIZ(2,7)=1.6E-5
      PPIZ(2,8)=1.6E-5
      PPIZ(2,9)=1.6E-5
      PPIZ(2,10)=1.6E-5
      PPIZ(2,11)=1.6E-5
      PPIZ(2,15)=1.6E-5
      PPIZ(2,16)=1.6E-5
      PPIZ(2,17)=1.6E-5
      PPIZ(2,20)=1.3E-5
      PPIZ(3,3)=1.07E-5
      PPIZ(3,4)=1.07E-5
      PPIZ(3,5)=1.07E-5
      PPIZ(3,6)=1.07E-5
      PPIZ(3,7)=1.6E-5
      PPIZ(3,8)=1.6E-5
      PPIZ(3,9)=1.6E-5
      PPIZ(3,10)=1.6E-5
      PPIZ(3,11)=1.3E-5
      PPIZ(3,15)=1.3E-5
      PPIZ(3,16)=1.6E-5
      PPIZ(3,17)=1.6E-5
      PPIZ(3,20)=1.3E-5
      PPIZ(3,41)=1.3E-5
      PPIZ(4,4)=1.07E-5
      PPIZ(4,5)=1.07E-5
      PPIZ(4,6)=1.07E-5
      PPIZ(4,7)=1.6E-5
      PPIZ(4,8)=1.6E-5
      PPIZ(4,9)=1.6E-5
      PPIZ(4,10)=1.6E-5
      PPIZ(4,11)=1.3E-5
      PPIZ(4,16)=1.6E-5
      PPIZ(4,17)=1.6E-5
      PPIZ(4,20)=1.3E-5
      PPIZ(5,5)=1.07E-5
      PPIZ(5,6)=1.07E-5
      PPIZ(5,7)=1.6E-5
      PPIZ(5,8)=1.6E-5
      PPIZ(5,9)=1.6E-5
      PPIZ(5,10)=1.6E-5
      PPIZ(5,11)=1.2E-5
      PPIZ(5,15)=1.3E-5
      PPIZ(5,16)=1.07E-5
      PPIZ(5,17)=1.6E-5
      PPIZ(5,20)=1.2E-5
      PPIZ(6,6)=1.07E-5
      PPIZ(6,7)=1.6E-5
      PPIZ(6,8)=1.6E-5
      PPIZ(6,9)=1.6E-5
      PPIZ(6,10)=1.6E-5
      PPIZ(6,16)=1.6E-5
      PPIZ(6,17)=1.6E-5
      PPIZ(6,20)=1.2E-5
      PPIZ(7,7)=1.07E-5
      PPIZ(7,8)=1.6E-5
      PPIZ(7,9)=1.6E-5
      PPIZ(7,10)=1.6E-5
      PPIZ(7,16)=1.6E-5
      PPIZ(7,17)=1.6E-5
      PPIZ(7,20)=1.3E-5
      PPIZ(8,8)=1.6E-5
      PPIZ(8,9)=1.6E-5
      PPIZ(8,10)=1.6E-5
      PPIZ(8,11)=1.6E-5
      PPIZ(8,16)=1.6E-5
      PPIZ(8,17)=1.6E-5
      PPIZ(8,20)=1.3E-5
      PPIZ(9,9)=1.6E-5
      PPIZ(9,10)=1.6E-5
      PPIZ(9,11)=1.6E-5
      PPIZ(9,16)=1.6E-5
      PPIZ(9,17)=1.6E-5
      PPIZ(9,20)=1.3E-5
      PPIZ(10,10)=1.6E-5
      PPIZ(10,11)=1.6E-5
      PPIZ(10,16)=1.6E-5
      PPIZ(10,17)=1.6E-5
      PPIZ(10,20)=1.3E-5
      PPIZ(11,11)=1.6E-5
      PPIZ(11,16)=1.6E-5
      PPIZ(11,17)=1.6E-5
      PPIZ(15,15)=2.0E-5
      PPIZ(15,20)=1.3E-5
      PPIZ(16,16)=1.6E-5
      PPIZ(16,17)=1.6E-5
      PPIZ(16,20)=1.6E-5
      PPIZ(17,17)=1.6E-5
      PPIZ(11,20)=1.3E-5
      PPIZ(17,20)=1.3E-5
      PPIZ(20,20)=0.858E-5
      PPIT(1,1)=1.67E-5
      PPIT(1,2)=1.67E-5
      PPIT(1,3)=1.0E-5
      PPIT(1,4)=1.0E-5
      PPIT(1,5)=0.9E-5
      PPIT(1,6)=0.9E-5
      PPIT(1,7)=1.2E-5
      PPIT(1,8)=1.2E-5
      PPIT(1,9)=1.2E-5
      PPIT(1,10)=1.2E-5
      PPIT(1,11)=1.67E-5
      PPIT(1,15)=1.6E-5
      PPIT(1,16)=1.2E-5
      PPIT(1,17)=1.2E-5
      PPIT(1,20)=1.1E-5
      PPIT(2,2)=1.67E-5
      PPIT(2,3)=1.0E-5
      PPIT(2,4)=1.0E-5
      PPIT(2,5)=0.9E-5
      PPIT(2,6)=0.9E-5
      PPIT(2,7)=1.2E-5
      PPIT(2,8)=1.2E-5
      PPIT(2,9)=1.2E-5
      PPIT(2,10)=1.2E-5
      PPIT(2,11)=1.67E-5
      PPIT(2,15)=1.6E-5
      PPIT(2,16)=1.2E-5
      PPIT(2,17)=1.2E-5
      PPIT(2,20)=1.1E-5
      PPIT(3,3)=0.8E-5
      PPIT(3,4)=0.8E-5
      PPIT(3,5)=0.8E-5
      PPIT(3,6)=0.8E-5
      PPIT(3,7)=1.2E-5
      PPIT(3,8)=1.2E-5
      PPIT(3,9)=1.2E-5
      PPIT(3,10)=1.2E-5
      PPIT(3,11)=1.0E-5
      PPIT(3,15)=1.0E-5
      PPIT(3,16)=1.2E-5
      PPIT(3,17)=1.2E-5
      PPIT(3,20)=1.0E-5
      PPIT(3,41)=1.0E-5
      PPIT(4,4)=0.8E-5
      PPIT(4,5)=0.8E-5
      PPIT(4,6)=0.8E-5
      PPIT(4,7)=1.2E-5
      PPIT(4,8)=1.2E-5
      PPIT(4,9)=1.2E-5
      PPIT(4,10)=1.2E-5
      PPIT(4,11)=1.0E-5
      PPIT(4,16)=1.2E-5
      PPIT(4,17)=1.2E-5
      PPIT(4,20)=1.0E-5
      PPIT(5,5)=0.8E-5
      PPIT(5,6)=0.8E-5
      PPIT(5,7)=1.2E-5
      PPIT(5,8)=1.2E-5
      PPIT(5,9)=1.2E-5
      PPIT(5,10)=1.2E-5
      PPIT(5,11)=0.9E-5
      PPIT(5,15)=1.0E-5
      PPIT(5,16)=1.2E-5
      PPIT(5,17)=1.2E-5
      PPIT(5,20)=0.9E-5
      PPIT(6,6)=0.8E-5
      PPIT(6,7)=1.2E-5
      PPIT(6,8)=1.2E-5
      PPIT(6,9)=1.2E-5
      PPIT(6,10)=1.2E-5
      PPIT(6,11)=0.9E-5
      PPIT(6,16)=1.2E-5
      PPIT(6,17)=1.2E-5
      PPIT(6,20)=0.9E-5
      PPIT(7,7)=0.8E-5
      PPIT(7,8)=1.2E-5
      PPIT(7,9)=1.2E-5
      PPIT(7,10)=1.2E-5
      PPIT(7,11)=1.2E-5
      PPIT(7,16)=1.2E-5
      PPIT(7,17)=1.2E-5
      PPIT(7,20)=1.0E-5
      PPIT(8,8)=1.2E-5
      PPIT(8,9)=1.2E-5
      PPIT(8,10)=1.2E-5
      PPIT(8,11)=1.2E-5
      PPIT(8,16)=1.2E-5
      PPIT(8,17)=1.2E-5
      PPIT(8,20)=1.0E-5
      PPIT(9,9)=1.2E-5
      PPIT(9,10)=1.2E-5
      PPIT(9,11)=1.2E-5
      PPIT(9,16)=1.2E-5
      PPIT(9,17)=1.2E-5
      PPIT(9,20)=1.0E-5
      PPIT(10,10)=1.2E-5
      PPIT(10,11)=1.2E-5
      PPIT(10,16)=1.2E-5
      PPIT(10,17)=1.2E-5
      PPIT(10,20)=1.0E-5
      PPIT(11,11)=1.67E-5
      PPIT(11,16)=1.2E-5
      PPIT(11,17)=1.2E-5
      PPIT(11,20)=1.1E-5
      PPIT(15,15)=1.6E-5
      PPIT(15,20)=1.0E-5
      PPIT(16,16)=1.2E-5
      PPIT(16,17)=1.2E-5
      PPIT(16,20)=1.2E-5
      PPIT(17,17)=1.2E-5
      PPIT(17,20)=1.2E-5
      PPIT(20,20)=0.58E-5
      AVDW(1,1)=301.20E3
      AVDW(1,2)=301.20E3
      AVDW(1,3)=340.00E3
      AVDW(1,4)=340.00E3
      AVDW(1,5)=278.70E3
      AVDW(1,6)=278.70E3
      AVDW(1,7)=255.40E3
      AVDW(1,8)=255.40E3
      AVDW(1,9)=112.80E3
      AVDW(1,10)=291.10E3
      AVDW(1,11)=301.20E3
      AVDW(1,12)=340.00E3
      AVDW(1,13)=25.E3
      AVDW(1,14)=100.E3
      AVDW(1,15)=291.1E3
      AVDW(1,16)=255.40E3
      AVDW(1,17)=255.40E3
      AVDW(1,18)=25.E3
      AVDW(1,19)=96.47E3
      AVDW(1,20)=96.47E3
      AVDW(1,21)=96.47E3
      AVDW(1,22)=96.47E3
      AVDW(1,23)=96.47E3
      AVDW(1,24)=96.47E3
      AVDW(1,25)=80.E3
      AVDW(1,26)=96.47E3
      AVDW(1,27)=240.E3
      AVDW(1,28)=96.47E3
      AVDW(1,29)=96.47E3
      AVDW(1,30)=96.47E3
      AVDW(1,31)=96.47E3
      AVDW(1,32)=35.E3
      AVDW(1,33)=25.E3
      AVDW(1,34)=96.47E3
      AVDW(1,35)=96.47E3
      AVDW(1,36)=301.20E3/3.
      AVDW(1,41)=AVDW(1,1)*0.8
      AVDW(2,2)=609.70E3
      AVDW(2,3)=340.00E3
      AVDW(2,4)=340.00E3
      AVDW(2,5)=278.70E3
      AVDW(2,6)=278.70E3
      AVDW(2,7)=255.40E3
      AVDW(2,8)=255.40E3
      AVDW(2,9)=112.80E3
      AVDW(2,10)=291.1E3
      AVDW(2,11)=301.20E3
      AVDW(2,12)=340.00E3
      AVDW(2,13)=52.0E3
      AVDW(2,14)=120.0E3
      AVDW(2,15)=291.1E3
      AVDW(2,16)=255.40E3
      AVDW(2,17)=255.40E3
      AVDW(2,18)=35.00E3
      AVDW(2,19)=44.80E3
      AVDW(2,20)=44.80E3
      AVDW(2,21)=44.80E3
      AVDW(2,22)=44.80E3
      AVDW(2,23)=255.E3
      AVDW(2,24)=44.80E3
      AVDW(2,25)=95.0E3
      AVDW(2,26)=64.80E3
      AVDW(2,27)=240.0E3
      AVDW(2,28)=44.80E3
      AVDW(2,29)=44.80E3
      AVDW(2,30)=44.80E3
      AVDW(2,31)=44.80E3
      AVDW(2,32)=680.E3
      AVDW(2,33)=35.E3
      AVDW(2,34)=38.E4
      AVDW(2,35)=25.E4
      AVDW(2,36)=301.20E3/10.
      AVDW(2,41)=AVDW(1,2)*0.8
      AVDW(3,3)=387.00E3
      AVDW(3,4)=387.00E3
      AVDW(3,5)=316.20E3
      AVDW(3,6)=316.20E3
      AVDW(3,7)=288.16E3
      AVDW(3,8)=288.16E3
      AVDW(3,9)=122.40E3
      AVDW(3,10)=325.90E3
      AVDW(3,11)=340.00E3
      AVDW(3,12)=387.00E3
      AVDW(3,13)=54.E3
      AVDW(3,14)=25.E3
      AVDW(3,15)=325.9E3
      AVDW(3,16)=288.16E3
      AVDW(3,17)=288.16E3
      AVDW(3,18)=84.0E3
      AVDW(3,19)=950.E3
      AVDW(3,20)=52.10E3
      AVDW(3,21)=52.10E3
      AVDW(3,22)=52.10E3
      AVDW(3,23)=52.10E3
      AVDW(3,24)=52.10E3
      AVDW(3,25)=45.E3
      AVDW(3,26)=52.10E3
      AVDW(3,27)=52.10E3
      AVDW(3,28)=52.10E3
      AVDW(3,29)=52.10E3
      AVDW(3,30)=52.10E3
      AVDW(3,31)=52.10E3
      AVDW(3,33)=56.E3
      AVDW(3,34)=52.1E3
      AVDW(3,35)=52.1E3
      AVDW(3,36)=340.00E3
      AVDW(3,41)=AVDW(1,3)*0.8
      AVDW(4,4)=387.00E3
      AVDW(4,5)=316.20E3
      AVDW(4,6)=316.20E3
      AVDW(4,7)=288.16E3
      AVDW(4,8)=288.16E3
      AVDW(4,9)=122.40E3
      AVDW(4,10)=325.90E3
      AVDW(4,11)=340.00E3
      AVDW(4,12)=387.00E3
      AVDW(4,13)=152.E3
      AVDW(4,14)=100.0E3
      AVDW(4,15)=325.9E3
      AVDW(4,16)=288.16E3
      AVDW(4,17)=288.16E3
      AVDW(4,18)=54.0E3
      AVDW(4,19)=52.10E3*1.1
      AVDW(4,20)=52.10E3
      AVDW(4,21)=20.E3
      AVDW(4,22)=98.0E3
      AVDW(4,23)=23.0E3
      AVDW(4,24)=29.0E3
      AVDW(4,25)=20.0E3
      AVDW(4,26)=2600.0E3*10.
      AVDW(4,27)=220.0E3
      AVDW(4,28)=65.0E3
      AVDW(4,29)=420.0E3
      AVDW(4,30)=300.0E3
      AVDW(4,31)=420.0E3
      AVDW(4,33)=54.E3
      AVDW(4,34)=23.E3
      AVDW(4,35)=23.E3
      AVDW(4,36)=340.00E3
      AVDW(4,41)=AVDW(1,4)/10.
      AVDW(5,5)=179.70E3
      AVDW(5,6)=259.00E3
      AVDW(5,7)=239.20E3
      AVDW(5,8)=239.20E3
      AVDW(5,9)=97.00E3
      AVDW(5,10)=272.70E3
      AVDW(5,11)=278.70E3
      AVDW(5,12)=316.20E3
      AVDW(5,13)=25.E3
      avdw(5,14)=60.0e3
      AVDW(5,15)=272.7E3
      AVDW(5,16)=239.20E3
      AVDW(5,17)=239.20E3
      AVDW(5,18)=15.00E3
      AVDW(5,19)=490.00E3
      AVDW(5,20)=42.00E3
      AVDW(5,21)=42.00E3
      AVDW(5,22)=42.00E3
      AVDW(5,23)=42.00E3
      AVDW(5,24)=42.00E3
      AVDW(5,25)=42.00E3
      AVDW(5,26)=50.00E3
      AVDW(5,27)=130.00E3
      AVDW(5,28)=42.00E3
      AVDW(5,29)=42.00E3
      AVDW(5,30)=42.00E3
      AVDW(5,31)=42.00E3
      AVDW(5,33)=15.E3
      AVDW(5,34)=42.E3
      AVDW(5,35)=42.E3
      AVDW(5,36)=278.70E3
      AVDW(5,41)=AVDW(1,5)*0.8
      AVDW(6,6)=259.00E3
      AVDW(6,7)=239.20E3
      AVDW(6,8)=239.20E3
      AVDW(6,9)=97.00E3
      AVDW(6,10)=272.70E3
      AVDW(6,11)=278.70E3
      AVDW(6,12)=316.20E3
      AVDW(6,13)=25.00E3
      AVDW(6,14)=100.0E3
      AVDW(6,15)=272.7E3
      AVDW(6,16)=239.20E3
      AVDW(6,17)=239.20E3
      AVDW(6,18)=18.E3
      AVDW(6,19)=42.00E3
      AVDW(6,20)=42.00E3
      AVDW(6,21)=42.00E3
      AVDW(6,22)=42.00E3
      AVDW(6,23)=82.00E3
      AVDW(6,24)=42.00E3
      AVDW(6,25)=28.00E3
      AVDW(6,26)=40.00E3*5.
      AVDW(6,27)=280.00E3
      AVDW(6,28)=42.00E3
      AVDW(6,29)=42.00E3
      AVDW(6,30)=42.00E3
      AVDW(6,31)=420.00E3
      AVDW(6,32)=140.E3
      AVDW(6,33)=18.E3
      AVDW(6,34)=200.E3
      AVDW(6,35)=380.E3
      AVDW(6,36)=27.E3
      AVDW(6,41)=AVDW(1,6)*0.8
      AVDW(7,7)=220.80E3
      AVDW(7,8)=220.80E3
      AVDW(7,9)=86.40E3
      AVDW(7,10)=251.60E3
      AVDW(7,11)=255.40E3
      AVDW(7,12)=288.16E3
      AVDW(7,13)=40.50E3
      AVDW(7,14)=40.50E3
      AVDW(7,15)=251.60E3
      AVDW(7,16)=220.80E3
      AVDW(7,17)=220.80E3
      AVDW(7,18)=40.50E3
      AVDW(7,19)=40.50E3
      AVDW(7,20)=40.50E3
      AVDW(7,21)=40.50E3
      AVDW(7,22)=40.50E3
      AVDW(7,23)=40.50E3
      AVDW(7,24)=40.50E3
      AVDW(7,25)=40.50E3
      AVDW(7,26)=4.10E2
      AVDW(7,27)=40.50E3
      AVDW(7,28)=40.50E3
      AVDW(7,29)=40.50E3
      AVDW(7,30)=40.50E3
      AVDW(7,31)=40.50E3
      AVDW(7,33)=40.5E3
      AVDW(7,34)=40.5E3
      AVDW(7,35)=40.50E3
      AVDW(7,36)=255.40E3
      AVDW(7,41)=AVDW(1,7)*0.8
      AVDW(8,8)=220.80E3
      AVDW(8,9)=86.40E3
      AVDW(8,10)=251.60E3
      AVDW(8,11)=255.40E3
      AVDW(8,12)=288.16E3
      AVDW(8,13)=40.50E3
      AVDW(8,14)=40.50E3
      AVDW(8,15)=251.6E3
      AVDW(8,16)=220.80E3
      AVDW(8,17)=220.80E3
      AVDW(8,18)=40.50E3
      AVDW(8,19)=40.50E3
      AVDW(8,20)=40.50E3
      AVDW(8,21)=40.50E3
      AVDW(8,22)=40.50E3
      AVDW(8,23)=40.50E3
      AVDW(8,24)=40.50E3
      AVDW(8,25)=40.50E3
      AVDW(8,26)=40.50E3
      AVDW(8,27)=40.50E3
      AVDW(8,28)=40.50E3
      AVDW(8,29)=40.50E3
      AVDW(8,30)=40.50E3
      AVDW(8,31)=40.50E3
      AVDW(8,33)=40.50E3
      AVDW(8,34)=40.50E3
      AVDW(8,35)=40.50E3
      AVDW(8,36)=255.40E3
      AVDW(8,41)=AVDW(1,8)*0.8
      AVDW(9,9)=54.20E3
      AVDW(9,10)=251.60E3
      AVDW(9,11)=112.80E3
      AVDW(9,12)=122.40E3
      AVDW(9,13)=12.50E3
      AVDW(9,14)=12.50E3
      AVDW(9,15)=49.10E3
      AVDW(9,16)=86.40E3
      AVDW(9,17)=86.40E3
      AVDW(9,18)=12.50E3
      AVDW(9,19)=12.50E3
      AVDW(9,20)=12.50E3
      AVDW(9,21)=12.50E3
      AVDW(9,22)=12.50E3
      AVDW(9,23)=12.50E3
      AVDW(9,24)=12.50E3
      AVDW(9,25)=12.50E3
      AVDW(9,26)=12.50E3
      AVDW(9,27)=12.50E3
      AVDW(9,28)=12.50E3
      AVDW(9,29)=12.50E3
      AVDW(9,30)=12.50E3
      AVDW(9,31)=12.50E3
      AVDW(9,33)=12.50E3
      AVDW(9,34)=12.50E3
      AVDW(9,35)=12.50E3
      AVDW(9,36)=112.80E3
      AVDW(9,41)=AVDW(1,9)*0.8
      AVDW(10,10)=273.90E3
      AVDW(10,11)=291.10E3
      AVDW(10,12)=325.90E3
      AVDW(10,13)=49.10E3
      AVDW(10,14)=49.10E3
      AVDW(10,15)=273.90E3
      AVDW(10,16)=251.60E3
      AVDW(10,17)=251.60E3
      AVDW(10,18)=49.10E3
      AVDW(10,19)=49.10E3
      AVDW(10,20)=49.10E3
      AVDW(10,21)=49.10E3
      AVDW(10,22)=49.10E3
      AVDW(10,23)=49.10E3
      AVDW(10,24)=49.10E3
      AVDW(10,25)=49.10E3
      AVDW(10,26)=49.10E3
      AVDW(10,27)=49.10E3
      AVDW(10,28)=49.10E3
      AVDW(10,29)=49.10E3
      AVDW(10,30)=49.10E3
      AVDW(10,31)=49.10E3
      AVDW(10,33)=49.10E3
      AVDW(10,34)=49.10E3
      AVDW(10,35)=49.10E3
      AVDW(10,36)=291.10E3
      AVDW(10,41)=AVDW(1,10)*0.8
      AVDW(11,11)=301.20E3
      AVDW(11,12)=340.00E3
      AVDW(11,13)=44.80E3
      AVDW(11,14)=44.80E3
      AVDW(11,15)=291.1E3
      AVDW(11,16)=255.40E3
      AVDW(11,17)=255.40E3
      AVDW(11,18)=44.80E3
      AVDW(11,19)=44.80E3
      AVDW(11,20)=44.80E3
      AVDW(11,21)=44.80E3
      AVDW(11,22)=44.80E3
      AVDW(11,23)=20.E3
      AVDW(11,24)=44.80E3
      AVDW(11,25)=44.80E3
      AVDW(11,26)=44.80E3
      AVDW(11,27)=44.80E3
      AVDW(11,28)=44.80E3
      AVDW(11,29)=44.80E3
      AVDW(11,30)=44.80E3
      AVDW(11,31)=44.80E3
      AVDW(11,32)=20.E3
      AVDW(11,33)=44.8E3
      AVDW(11,34)=20.E3
      AVDW(11,35)=20.E3
      AVDW(11,36)=301.20E3/10.
      AVDW(11,41)=AVDW(1,11)*0.8
      AVDW(12,12)=387.0E3
      AVDW(12,13)=52.10E3
      AVDW(12,14)=52.10E3
      AVDW(12,15)=325.90E3
      AVDW(12,16)=288.16E3
      AVDW(12,17)=288.16E3
      AVDW(12,18)=52.10E3
      AVDW(12,19)=52.10E3
      AVDW(12,20)=52.10E3
      AVDW(12,21)=52.10E3
      AVDW(12,22)=52.10E3
      AVDW(12,23)=52.10E3
      AVDW(12,24)=52.10E3
      AVDW(12,25)=52.10E3
      AVDW(12,26)=52.10E3
      AVDW(12,27)=52.10E3
      AVDW(12,28)=52.10E3
      AVDW(12,29)=52.10E3
      AVDW(12,30)=52.10E3
      AVDW(12,31)=52.10E3
      AVDW(12,33)=52.10E3
      AVDW(12,34)=52.10E3
      AVDW(12,35)=52.10E3
      AVDW(12,36)=340.00E3
      AVDW(12,41)=AVDW(1,12)*0.8
      AVDW(13,13)=52.1E3
      AVDW(13,16)=52.1E3
      AVDW(13,17)=288.16E3
      AVDW(13,20)=6.60E3
      AVDW(13,21)=6.60E3
      AVDW(13,22)=6.60E3
      AVDW(13,23)=6.60E3
      AVDW(13,24)=6.60E3
      AVDW(13,25)=6.60E3
      AVDW(13,26)=6.60E3
      AVDW(13,27)=6.60E3
      AVDW(13,28)=6.60E3
      AVDW(13,29)=6.60E3
      AVDW(13,30)=6.60E3
      AVDW(13,31)=6.60E3
      AVDW(13,34)=6.60E3
      AVDW(13,35)=6.60E3
      AVDW(13,41)=13.E3
      AVDW(14,15)=100.0E3
      AVDW(14,16)=40.50E3
      AVDW(14,17)=40.50E3
      AVDW(14,20)=6.60E3
      AVDW(14,21)=6.60E3
      AVDW(14,22)=6.60E3
      AVDW(14,23)=6.60E3
      AVDW(14,24)=6.60E3
      AVDW(14,25)=6.60E3
      AVDW(14,26)=6.60E3
      AVDW(14,27)=6.60E3
      AVDW(14,28)=6.60E3
      AVDW(14,29)=6.60E3
      AVDW(14,30)=6.60E3
      AVDW(14,31)=6.60E3
      AVDW(14,34)=6.60E3
      AVDW(14,35)=6.60E3
      AVDW(14,41)=AVDW(1,14)*0.8
      AVDW(15,15)=273.90E3
      AVDW(15,16)=273.90E3
      AVDW(15,17)=251.60E3
      AVDW(15,20)=49.10E3
      AVDW(15,21)=49.10E3
      AVDW(15,22)=49.10E3
      AVDW(15,23)=49.10E3
      AVDW(15,24)=49.10E3
      AVDW(15,25)=49.10E3
      AVDW(15,26)=49.10E3
      AVDW(15,27)=380.E3
      AVDW(15,28)=49.10E3
      AVDW(15,29)=49.10E3
      AVDW(15,30)=49.10E3
      AVDW(15,31)=49.10E3
      AVDW(15,34)=4.E3
      AVDW(15,35)=4.E3
      AVDW(15,36)=290.E3
      AVDW(15,41)=AVDW(1,15)*0.8
      AVDW(16,16)=220.80E3
      AVDW(16,17)=220.80E3
      AVDW(16,18)=40.50E3
      AVDW(16,19)=40.50E3
      AVDW(16,20)=40.50E3
      AVDW(16,21)=40.50E3
      AVDW(16,22)=40.50E3
      AVDW(16,23)=40.50E3
      AVDW(16,24)=40.50E3
      AVDW(16,25)=19.0E3
      AVDW(16,26)=40.50E3
      AVDW(16,27)=40.50E3
      AVDW(16,28)=40.50E3
      AVDW(16,29)=40.50E3
      AVDW(16,30)=40.50E3
      AVDW(16,31)=40.50E3
      AVDW(16,33)=40.50E3
      AVDW(16,34)=40.50E3
      AVDW(16,35)=40.50E3
      AVDW(16,41)=AVDW(1,16)*0.8
      AVDW(17,17)=220.80E3
      AVDW(17,18)=40.50E3
      AVDW(17,19)=40.50E3
      AVDW(17,20)=40.50E3
      AVDW(17,21)=40.50E3
      AVDW(17,22)=40.50E3
      AVDW(17,23)=40.50E3
      AVDW(17,24)=40.50E3
      AVDW(17,25)=40.50E3
      AVDW(17,26)=40.50E3
      AVDW(17,27)=40.50E3
      AVDW(17,28)=40.50E3
      AVDW(17,29)=40.50E3
      AVDW(17,30)=40.50E3
      AVDW(17,31)=40.50E3
      AVDW(17,33)=40.50E3
      AVDW(17,34)=40.50E3
      AVDW(17,35)=40.50E3
      AVDW(17,41)=AVDW(1,17)*0.8
      AVDW(18,20)=6.60E3
      AVDW(18,21)=6.60E3
      AVDW(18,22)=6.60E3
      AVDW(18,23)=6.60E3
      AVDW(18,24)=6.60E3
      AVDW(18,25)=6.60E3
      AVDW(18,26)=6.60E3
      AVDW(18,27)=6.60E3
      AVDW(18,28)=6.60E3
      AVDW(18,29)=6.60E3
      AVDW(18,30)=6.60E3
      AVDW(18,31)=6.60E3
      AVDW(18,34)=6.60E3
      AVDW(18,35)=6.60E3
      AVDW(18,41)=AVDW(1,18)*0.8
      AVDW(19,20)=6.60E3
      AVDW(19,21)=6.60E3
      AVDW(19,22)=6.60E3
      AVDW(19,23)=6.60E3
      AVDW(19,24)=6.60E3
      AVDW(19,25)=6.60E3
      AVDW(19,26)=6.60E3
      AVDW(19,27)=6.60E3
      AVDW(19,28)=6.60E3
      AVDW(19,29)=6.60E3
      AVDW(19,30)=6.60E3
      AVDW(19,31)=6.60E3
      AVDW(19,34)=6.60E3
      AVDW(19,35)=6.60E3
      AVDW(19,41)=AVDW(1,19)*0.8
      AVDW(20,20)=6.60E3
      AVDW(20,21)=6.60E3
      AVDW(20,22)=6.60E3
      AVDW(20,23)=6.60E3
      AVDW(20,24)=6.60E3
      AVDW(20,25)=50.E3
      AVDW(20,26)=6.60E3
      AVDW(20,27)=100.0E3
      AVDW(20,28)=6.60E3
      AVDW(20,29)=6.60E3
      AVDW(20,30)=6.60E3
      AVDW(20,31)=6.60E3
      AVDW(20,32)=80.E3
      AVDW(20,34)=60.E3
      AVDW(20,35)=60.E3
      AVDW(20,36)=44.80E3/2.
      AVDW(20,41)=AVDW(1,20)*0.8
      AVDW(23,36)=5.E3
      AVDW(34,34)=750.E3
      AVDW(34,36)=16.E3
      AVDW(35,35)=750.E3
      AVDW(35,36)=15.E3
      AVDW(36,36)=31.20E3
      AVDW(41,41)=AVDW(1,1)*0.8
      BVDW(1,7)=1.811
      BVDW(1,8)=1.811
      BVDW(1,9)=1.940
      BVDW(1,10)=1.665
      BVDW(1,13)=1.665
      BVDW(1,14)=1.665
      BVDW(1,15)=1.665
      BVDW(1,16)=1.811
      BVDW(1,17)=1.811
      BVDW(1,18)=1.665
      BVDW(1,19)=1.795
      BVDW(1,20)=1.795
      BVDW(1,21)=1.795
      BVDW(1,22)=1.795
      BVDW(1,23)=1.795
      BVDW(1,24)=1.795
      BVDW(1,25)=1.665
      BVDW(1,26)=1.795
      BVDW(1,27)=1.665
      BVDW(1,28)=1.795
      BVDW(1,29)=1.795
      BVDW(1,30)=1.795
      BVDW(1,31)=1.795
      BVDW(1,32)=1.665
      BVDW(1,33)=1.665
      BVDW(1,34)=1.665
      BVDW(1,35)=1.665
      BVDW(2,7)=1.811
      BVDW(2,8)=1.811
      BVDW(2,9)=1.940
      BVDW(2,10)=1.665
      BVDW(2,13)=1.665
      BVDW(2,14)=1.665
      BVDW(2,15)=1.665
      BVDW(2,16)=1.811
      BVDW(2,17)=1.811
      BVDW(2,18)=1.665
      BVDW(2,19)=2.040
      BVDW(2,20)=2.040
      BVDW(2,21)=2.040
      BVDW(2,22)=2.040
      BVDW(2,23)=2.040
      BVDW(2,24)=2.040
      BVDW(2,25)=2.040
      BVDW(2,26)=2.040
      BVDW(2,27)=2.040
      BVDW(2,28)=2.040
      BVDW(2,29)=2.040
      BVDW(2,30)=2.040
      BVDW(2,31)=2.040
      BVDW(2,32)=2.040
      BVDW(2,33)=1.665
      BVDW(2,34)=1.
      BVDW(2,35)=1.
      BVDW(2,36)=1.
      BVDW(3,7)=1.811
      BVDW(3,8)=1.811
      BVDW(3,9)=1.940
      BVDW(3,10)=1.665
      BVDW(3,13)=1.665
      BVDW(3,14)=1.665
      BVDW(3,15)=1.665
      BVDW(3,16)=1.811
      BVDW(3,17)=1.811
      BVDW(3,18)=1.665
      BVDW(3,19)=1.665
      BVDW(3,20)=2.040
      BVDW(3,21)=2.040
      BVDW(3,22)=2.040
      BVDW(3,23)=2.040
      BVDW(3,24)=2.040
      BVDW(3,25)=2.040
      BVDW(3,26)=2.040
      BVDW(3,27)=2.040
      BVDW(3,28)=2.040
      BVDW(3,29)=2.040
      BVDW(3,30)=2.040
      BVDW(3,31)=2.040
      BVDW(3,33)=1.665
      BVDW(3,34)=2.040
      BVDW(3,35)=2.040
      BVDW(4,7)=1.811
      BVDW(4,8)=1.811
      BVDW(4,9)=1.940
      BVDW(4,10)=1.665
      BVDW(4,13)=1.665
      BVDW(4,14)=1.665
      BVDW(4,15)=1.665
      BVDW(4,16)=1.811
      BVDW(4,17)=1.811
      BVDW(4,18)=1.665
      BVDW(4,19)=2.040
      BVDW(4,20)=2.040
      BVDW(4,21)=1.665
      BVDW(4,22)=2.040
      BVDW(4,23)=1.665
      BVDW(4,24)=1.665
      BVDW(4,25)=1.665
      BVDW(4,26)=1.665
      BVDW(4,27)=1.665
      BVDW(4,28)=1.665
      BVDW(4,29)=1.665
      BVDW(4,30)=1.665
      BVDW(4,31)=1.665
      BVDW(4,33)=1.665
      BVDW(4,34)=1.665
      BVDW(4,35)=1.665
      BVDW(5,7)=1.811
      BVDW(5,8)=1.811
      BVDW(5,9)=1.940
      BVDW(5,10)=1.665
      BVDW(5,13)=1.665
      bvdw(5,14)=1.665
      BVDW(5,15)=1.665
      BVDW(5,16)=1.811
      BVDW(5,17)=1.811
      BVDW(5,18)=1.665
      BVDW(5,19)=1.665
      BVDW(5,20)=2.040
      BVDW(5,21)=2.040
      BVDW(5,22)=2.040
      BVDW(5,23)=2.040
      BVDW(5,24)=2.040
      BVDW(5,25)=2.040
      BVDW(5,26)=2.040
      BVDW(5,27)=1.665
      BVDW(5,28)=2.040
      BVDW(5,29)=2.040
      BVDW(5,30)=2.040
      BVDW(5,31)=2.040
      BVDW(5,33)=1.665
      BVDW(5,34)=2.040
      BVDW(5,35)=2.040
      BVDW(6,7)=1.811
      BVDW(6,8)=1.811
      BVDW(6,9)=1.940
      BVDW(6,10)=1.665
      BVDW(6,13)=1.665
      BVDW(6,14)=1.665
      BVDW(6,15)=1.665
      BVDW(6,16)=1.811
      BVDW(6,17)=1.811
      BVDW(6,18)=1.665
      BVDW(6,19)=2.040
      BVDW(6,20)=2.040
      BVDW(6,21)=2.040
      BVDW(6,22)=2.040
      BVDW(6,23)=2.040
      BVDW(6,24)=2.040
      BVDW(6,25)=2.040
      BVDW(6,26)=1.665
      BVDW(6,27)=2.040
      BVDW(6,28)=2.040
      BVDW(6,29)=2.040
      BVDW(6,30)=2.040
      BVDW(6,31)=2.040
      BVDW(6,32)=1.665
      BVDW(6,33)=1.665
      BVDW(6,34)=1.665
      BVDW(6,35)=1.665
      BVDW(6,36)=1.665
      BVDW(7,7)=3.621
      BVDW(7,8)=3.621
      BVDW(7,9)=3.621
      BVDW(7,10)=3.475
      BVDW(7,11)=1.811
      BVDW(7,12)=1.811
      BVDW(7,13)=3.851
      BVDW(7,14)=3.851
      BVDW(7,15)=3.475
      BVDW(7,16)=3.851
      BVDW(7,17)=3.621
      BVDW(7,18)=3.851
      BVDW(7,19)=3.851
      BVDW(7,20)=3.851
      BVDW(7,21)=3.851
      BVDW(7,22)=3.851
      BVDW(7,23)=3.851
      BVDW(7,24)=3.851
      BVDW(7,25)=3.851
      BVDW(7,26)=0.665
      BVDW(7,27)=3.851
      BVDW(7,28)=3.851
      BVDW(7,29)=3.851
      BVDW(7,30)=3.851
      BVDW(7,31)=3.851
      BVDW(7,33)=3.851
      BVDW(7,34)=3.851
      BVDW(7,35)=3.851
      BVDW(7,36)=1.811
      BVDW(7,41)=1.811
      BVDW(8,8)=3.621
      BVDW(8,9)=3.621
      BVDW(8,10)=3.475
      BVDW(8,11)=1.811
      BVDW(8,12)=1.811
      BVDW(8,13)=3.851
      BVDW(8,14)=3.851
      BVDW(8,15)=3.475
      BVDW(8,16)=3.621
      BVDW(8,17)=3.621
      BVDW(8,18)=3.851
      BVDW(8,19)=3.851
      BVDW(8,20)=3.851
      BVDW(8,21)=3.851
      BVDW(8,22)=3.851
      BVDW(8,23)=3.851
      BVDW(8,24)=3.851
      BVDW(8,25)=3.851
      BVDW(8,26)=3.851
      BVDW(8,27)=3.851
      BVDW(8,28)=3.851
      BVDW(8,29)=3.851
      BVDW(8,30)=3.851
      BVDW(8,31)=3.851
      BVDW(8,33)=3.851
      BVDW(8,34)=3.851
      BVDW(8,35)=3.851
      BVDW(8,36)=1.811
      BVDW(8,41)=1.811
      BVDW(9,9)=3.621
      BVDW(9,10)=3.475
      BVDW(9,11)=1.940
      BVDW(9,12)=1.940
      BVDW(9,13)=4.080
      BVDW(9,14)=4.080
      BVDW(9,15)=3.475
      BVDW(9,16)=3.621
      BVDW(9,17)=3.621
      BVDW(9,18)=4.080
      BVDW(9,19)=4.080
      BVDW(9,20)=4.080
      BVDW(9,21)=4.080
      BVDW(9,22)=4.080
      BVDW(9,23)=4.080
      BVDW(9,24)=4.080
      BVDW(9,25)=4.080
      BVDW(9,26)=4.080
      BVDW(9,27)=4.080
      BVDW(9,28)=4.080
      BVDW(9,29)=4.080
      BVDW(9,30)=4.080
      BVDW(9,31)=4.080
      BVDW(9,33)=4.080
      BVDW(9,34)=4.080
      BVDW(9,35)=4.080
      BVDW(9,36)=1.940
      BVDW(9,41)=1.940
      BVDW(10,10)=3.329
      BVDW(10,11)=1.665
      BVDW(10,12)=1.665
      BVDW(10,13)=3.705
      BVDW(10,14)=3.705
      BVDW(10,15)=3.329
      BVDW(10,16)=3.475
      BVDW(10,17)=3.475
      BVDW(10,18)=3.705
      BVDW(10,19)=3.705
      BVDW(10,20)=3.705
      BVDW(10,21)=3.705
      BVDW(10,22)=3.705
      BVDW(10,23)=3.705
      BVDW(10,24)=3.705
      BVDW(10,25)=3.705
      BVDW(10,26)=3.705
      BVDW(10,27)=3.705
      BVDW(10,28)=3.705
      BVDW(10,29)=3.705
      BVDW(10,30)=3.705
      BVDW(10,31)=3.705
      BVDW(10,33)=3.705
      BVDW(10,34)=3.705
      BVDW(10,35)=3.705
      BVDW(10,36)=1.665
      BVDW(10,41)=1.665
      BVDW(11,13)=1.665
      BVDW(11,14)=1.665
      BVDW(11,15)=1.665
      BVDW(11,16)=1.811
      BVDW(11,17)=1.811
      BVDW(11,18)=2.040
      BVDW(11,19)=2.040
      BVDW(11,20)=2.040
      BVDW(11,21)=3.705
      BVDW(11,22)=3.705
      BVDW(11,23)=1.665
      BVDW(11,24)=3.705
      BVDW(11,25)=3.705
      BVDW(11,26)=3.705
      BVDW(11,27)=3.705
      BVDW(11,28)=3.705
      BVDW(11,29)=3.705
      BVDW(11,30)=3.705
      BVDW(11,31)=3.705
      BVDW(11,32)=1.665
      BVDW(11,33)=2.040
      BVDW(11,34)=1.665
      BVDW(11,35)=1.665
      BVDW(11,36)=1.
      BVDW(12,13)=1.665
      BVDW(12,14)=1.665
      BVDW(12,15)=1.665
      BVDW(12,16)=1.811
      BVDW(12,17)=1.811
      BVDW(12,18)=2.040
      BVDW(12,19)=2.040
      BVDW(12,20)=2.040
      BVDW(12,21)=2.040
      BVDW(12,22)=2.040
      BVDW(12,23)=2.040
      BVDW(12,24)=2.040
      BVDW(12,25)=2.040
      BVDW(12,26)=2.040
      BVDW(12,27)=2.040
      BVDW(12,28)=2.040
      BVDW(12,29)=2.040
      BVDW(12,30)=2.040
      BVDW(12,31)=2.040
      BVDW(12,33)=2.040
      BVDW(12,34)=2.040
      BVDW(12,35)=2.040
      BVDW(13,13)=1.665
      BVDW(13,16)=1.665
      BVDW(13,17)=3.851
      BVDW(13,20)=4.896
      BVDW(13,20)=4.896
      BVDW(13,21)=4.896
      BVDW(13,22)=4.896
      BVDW(13,23)=4.896
      BVDW(13,24)=4.896
      BVDW(13,25)=4.896
      BVDW(13,26)=4.896
      BVDW(13,27)=4.896
      BVDW(13,28)=4.896
      BVDW(13,29)=4.896
      BVDW(13,30)=4.896
      BVDW(13,31)=4.896
      BVDW(13,34)=4.896
      BVDW(13,35)=4.896
      BVDW(14,15)=1.665
      BVDW(14,16)=3.851
      BVDW(14,17)=3.851
      BVDW(14,20)=4.896
      BVDW(14,20)=4.896
      BVDW(14,21)=4.896
      BVDW(14,22)=4.896
      BVDW(14,23)=4.896
      BVDW(14,24)=4.896
      BVDW(14,25)=4.896
      BVDW(14,26)=4.896
      BVDW(14,27)=4.896
      BVDW(14,28)=4.896
      BVDW(14,29)=4.896
      BVDW(14,30)=4.896
      BVDW(14,31)=4.896
      BVDW(14,34)=4.896
      BVDW(14,35)=4.896
      BVDW(15,15)=3.329
      BVDW(15,16)=3.475
      BVDW(15,20)=3.705
      BVDW(15,27)=3.329
      BVDW(15,34)=1.665
      BVDW(15,35)=1.665
      BVDW(15,36)=1.665
      BVDW(16,16)=3.621
      BVDW(16,17)=3.621
      BVDW(16,18)=3.851
      BVDW(16,19)=3.851
      BVDW(16,20)=3.851
      BVDW(16,25)=1.665
      BVDW(17,17)=3.621
      BVDW(17,18)=3.851
      BVDW(17,19)=3.851
      BVDW(17,20)=3.851
      BVDW(18,20)=4.896
      BVDW(18,21)=4.896
      BVDW(18,22)=4.896
      BVDW(18,23)=4.896
      BVDW(18,24)=4.896
      BVDW(18,25)=4.896
      BVDW(18,26)=4.896
      BVDW(18,27)=4.896
      BVDW(18,28)=4.896
      BVDW(18,29)=4.896
      BVDW(18,30)=4.896
      BVDW(18,31)=4.896
      BVDW(18,41)=1.665
      BVDW(19,20)=4.896
      BVDW(19,21)=4.896
      BVDW(19,22)=4.896
      BVDW(19,23)=4.896
      BVDW(19,24)=4.896
      BVDW(19,25)=4.896
      BVDW(19,26)=4.896
      BVDW(19,27)=4.896
      BVDW(19,28)=4.896
      BVDW(19,29)=4.896
      BVDW(19,30)=4.896
      BVDW(19,31)=4.896
      BVDW(19,41)=1.795
      BVDW(20,20)=4.896
      BVDW(20,21)=4.896
      BVDW(20,22)=4.896
      BVDW(20,23)=4.896
      BVDW(20,24)=4.896
      BVDW(20,25)=4.896
      BVDW(20,26)=4.896
      BVDW(20,27)=1.665
      BVDW(20,28)=4.896
      BVDW(20,29)=4.896
      BVDW(20,30)=4.896
      BVDW(20,31)=4.896
      BVDW(20,32)=1.665
      BVDW(20,33)=4.896
      BVDW(20,34)=4.896
      BVDW(20,35)=4.896
      BVDW(20,36)=2.040
      BVDW(20,41)=1.795
      BVDW(23,36)=1.665
      BVDW(34,34)=1.
      BVDW(34,36)=1.665
      BVDW(35,35)=1.
      BVDW(35,36)=1.665
      CVDW(1,1)=327.2
      CVDW(1,2)=327.2
      CVDW(1,3)=340.0
      CVDW(1,4)=340.0
      CVDW(1,5)=342.3
      CVDW(1,6)=342.3
      CVDW(1,7)=684.0
      CVDW(1,8)=684.0
      CVDW(1,9)=154.1
      CVDW(1,10)=981.1
      CVDW(1,11)=327.2
      CVDW(1,12)=340.0
      CVDW(1,13)=1026.3*5.
      CVDW(1,14)=1026.3
      CVDW(1,15)=981.1
      CVDW(1,16)=684.0
      CVDW(1,17)=684.0
      CVDW(1,18)=1026.3
      CVDW(1,19)=269.0
      CVDW(1,20)=269.0
      CVDW(1,21)=269.0
      CVDW(1,22)=269.0
      CVDW(1,23)=269.0
      CVDW(1,24)=269.0
      CVDW(1,25)=2200.
      CVDW(1,26)=269.0
      CVDW(1,27)=269.0
      CVDW(1,28)=269.0
      CVDW(1,29)=269.0
      CVDW(1,30)=269.0
      CVDW(1,31)=269.0
      CVDW(1,32)=1500.
      CVDW(1,33)=1026.3
      CVDW(1,34)=269.
      CVDW(1,35)=269.0
      CVDW(1,36)=327.2
      CVDW(1,41)=CVDW(1,1)*0.8
      CVDW(2,2)=465.4
      CVDW(2,3)=340.0
      CVDW(2,4)=340.0
      CVDW(2,5)=342.3
      CVDW(2,6)=342.3
      CVDW(2,7)=684.0
      CVDW(2,8)=684.0
      CVDW(2,9)=154.1
      CVDW(2,10)=981.1
      CVDW(2,11)=327.2
      CVDW(2,12)=340.0
      CVDW(2,13)=1026.3
      CVDW(2,14)=1026.3
      CVDW(2,15)=981.1
      CVDW(2,16)=684.0
      CVDW(2,17)=684.0
      CVDW(2,18)=1026.3
      CVDW(2,19)=125.0
      CVDW(2,20)=125.0
      CVDW(2,21)=125.0
      CVDW(2,22)=125.0
      CVDW(2,23)=125.0
      CVDW(2,24)=125.0
      CVDW(2,25)=125.0
      CVDW(2,26)=125.0
      CVDW(2,27)=125.0
      CVDW(2,28)=125.0
      CVDW(2,29)=125.0
      CVDW(2,30)=125.0
      CVDW(2,31)=125.0
      CVDW(2,32)=300.
      CVDW(2,33)=1600.
      CVDW(2,34)=1000.
      CVDW(2,35)=1000.
      CVDW(2,36)=327.2
      CVDW(2,41)=CVDW(1,2)*0.8
      CVDW(3,3)=354.0
      CVDW(3,4)=354.0
      CVDW(3,5)=356.0
      CVDW(3,6)=356.0
      CVDW(3,7)=711.5
      CVDW(3,8)=711.5
      CVDW(3,9)=154.1
      CVDW(3,10)=1020.5
      CVDW(3,11)=340.0
      CVDW(3,12)=354.0
      CVDW(3,13)=1020.5
      CVDW(3,14)=1026.3
      CVDW(3,15)=1020.5
      CVDW(3,16)=711.5
      CVDW(3,17)=711.5
      CVDW(3,18)=3500.
      CVDW(3,19)=15000.
      CVDW(3,20)=132.0
      CVDW(3,21)=132.0
      CVDW(3,22)=132.0
      CVDW(3,23)=132.0
      CVDW(3,24)=132.0
      CVDW(3,25)=1320.0
      CVDW(3,26)=132.0
      CVDW(3,27)=132.0
      CVDW(3,28)=132.0
      CVDW(3,29)=132.0
      CVDW(3,30)=132.0
      CVDW(3,31)=132.0
      CVDW(3,33)=132.0
      CVDW(3,34)=132.0
      CVDW(3,36)=340.0
      CVDW(3,41)=CVDW(1,3)*0.8
      CVDW(4,4)=354.0
      CVDW(4,5)=356.0
      CVDW(4,6)=356.0
      CVDW(4,7)=711.5
      CVDW(4,8)=711.5
      CVDW(4,9)=154.1
      CVDW(4,10)=1020.5
      CVDW(4,11)=340.0
      CVDW(4,12)=354.0
      CVDW(4,13)=1020.5*5.
      CVDW(4,14)=1026.3
      CVDW(4,15)=1020.5
      CVDW(4,16)=711.5
      CVDW(4,17)=711.5
      CVDW(4,18)=1020.5
      CVDW(4,19)=132.0*7.
      CVDW(4,20)=132.0
      CVDW(4,21)=132.0
      CVDW(4,22)=132.0
      CVDW(4,23)=1000.
      CVDW(4,24)=1400.
      CVDW(4,25)=1000.
      CVDW(4,26)=200.0*250.
      CVDW(4,27)=132.0
      CVDW(4,28)=900.
      CVDW(4,29)=7200.
      CVDW(4,30)=5800.
      CVDW(4,31)=5000.
      CVDW(4,33)=1020.5
      CVDW(4,34)=132.0
      CVDW(4,35)=132.0
      CVDW(4,36)=340.0
      CVDW(4,41)=CVDW(1,4)*0.8
      CVDW(5,5)=298.2
      CVDW(5,6)=358.0
      CVDW(5,7)=715.5
      CVDW(5,8)=715.5
      CVDW(5,9)=496.7
      CVDW(5,10)=1026.3
      CVDW(5,11)=342.3
      CVDW(5,12)=356.0
      CVDW(5,13)=1026.3
      cvdw(5,14)=1026.3
      CVDW(5,15)=1026.3
      CVDW(5,16)=715.5
      CVDW(5,17)=715.5
      CVDW(5,18)=1026.3
      CVDW(5,19)=7200.
      CVDW(5,20)=132.7
      CVDW(5,21)=132.7
      CVDW(5,22)=132.7
      CVDW(5,23)=132.7
      CVDW(5,24)=132.7
      CVDW(5,25)=132.7
      CVDW(5,26)=132.7
      CVDW(5,27)=2000.
      CVDW(5,28)=132.7
      CVDW(5,29)=132.7
      CVDW(5,30)=132.7
      CVDW(5,31)=132.7
      CVDW(5,33)=1026.3
      CVDW(5,34)=132.7
      CVDW(5,35)=132.7
      CVDW(5,36)=342.3
      CVDW(5,41)=CVDW(1,5)*0.8
      CVDW(6,6)=358.0
      CVDW(6,7)=715.5
      CVDW(6,8)=715.5
      CVDW(6,9)=496.7
      CVDW(6,10)=1026.3
      CVDW(6,11)=342.3
      CVDW(6,12)=356.0
      CVDW(6,13)=1026.3
      CVDW(6,14)=1026.3
      CVDW(6,15)=1026.3
      CVDW(6,16)=715.5
      CVDW(6,17)=715.5
      CVDW(6,18)=1026.3
      CVDW(6,19)=132.7
      CVDW(6,20)=132.7
      CVDW(6,21)=132.7
      CVDW(6,22)=132.7
      CVDW(6,23)=132.7
      CVDW(6,24)=132.7
      CVDW(6,25)=132.7
      CVDW(6,26)=132.7*5.
      CVDW(6,27)=132.7
      CVDW(6,28)=132.7
      CVDW(6,29)=132.7
      CVDW(6,30)=132.7
      CVDW(6,31)=132.7
      CVDW(6,33)=1026.3
      CVDW(6,34)=2000.
      CVDW(6,35)=150.
      CVDW(6,36)=342.3
      CVDW(6,41)=CVDW(1,6)*0.8
      CVDW(7,7)=1430.0
      CVDW(7,8)=1430.0
      CVDW(7,9)=633.5
      CVDW(7,10)=2051.1
      CVDW(7,11)=684.0
      CVDW(7,12)=711.5
      CVDW(7,13)=265.2
      CVDW(7,14)=265.2
      CVDW(7,15)=2051.1
      CVDW(7,16)=1430.0
      CVDW(7,17)=1430.0
      CVDW(7,18)=265.2
      CVDW(7,19)=265.2
      CVDW(7,20)=265.2
      CVDW(7,21)=265.2
      CVDW(7,22)=265.2
      CVDW(7,23)=265.2
      CVDW(7,24)=265.2
      CVDW(7,25)=265.2
      CVDW(7,26)=26.2
      CVDW(7,27)=265.2
      CVDW(7,28)=265.2
      CVDW(7,29)=265.2
      CVDW(7,30)=265.2
      CVDW(7,31)=265.2
      CVDW(7,33)=265.2
      CVDW(7,34)=265.2
      CVDW(7,35)=265.2
      CVDW(7,36)=684.0
      CVDW(7,41)=CVDW(1,7)*0.8
      CVDW(8,8)=1430.0
      CVDW(8,9)=633.5
      CVDW(8,10)=2051.1
      CVDW(8,11)=684.0
      CVDW(8,12)=711.5
      CVDW(8,13)=265.2
      CVDW(8,14)=265.2
      CVDW(8,15)=2051.1
      CVDW(8,16)=1430.0
      CVDW(8,17)=1430.0
      CVDW(8,18)=265.2
      CVDW(8,19)=265.2
      CVDW(8,20)=265.2
      CVDW(8,21)=265.2
      CVDW(8,22)=265.2
      CVDW(8,23)=265.2
      CVDW(8,24)=265.2
      CVDW(8,25)=265.2
      CVDW(8,26)=265.2
      CVDW(8,27)=265.2
      CVDW(8,28)=265.2
      CVDW(8,29)=265.2
      CVDW(8,30)=265.2
      CVDW(8,31)=265.2
      CVDW(8,33)=265.2
      CVDW(8,34)=265.2
      CVDW(8,35)=265.2
      CVDW(8,36)=684.0
      CVDW(8,41)=CVDW(1,8)*0.8
      CVDW(9,9)=135.2
      CVDW(9,10)=2051.1
      CVDW(9,11)=154.1
      CVDW(9,12)=154.1
      CVDW(9,13)=118.8
      CVDW(9,14)=118.8
      CVDW(9,15)=2051.1
      CVDW(9,16)=633.5
      CVDW(9,17)=633.5
      CVDW(9,18)=118.8
      CVDW(9,19)=118.8
      CVDW(9,20)=118.8
      CVDW(9,21)=118.8
      CVDW(9,22)=118.8
      CVDW(9,23)=118.8
      CVDW(9,24)=118.8
      CVDW(9,25)=118.8
      CVDW(9,26)=118.8
      CVDW(9,27)=118.8
      CVDW(9,28)=118.8
      CVDW(9,29)=118.8
      CVDW(9,30)=118.8
      CVDW(9,31)=118.8
      CVDW(9,33)=118.8
      CVDW(9,34)=118.8
      CVDW(9,35)=118.8
      CVDW(9,36)=154.1
      CVDW(9,41)=CVDW(1,9)*0.8
      CVDW(10,10)=2942.0
      CVDW(10,11)=981.1
      CVDW(10,12)=1020.5
      CVDW(10,13)=380.5
      CVDW(10,14)=380.5
      CVDW(10,15)=2942.0
      CVDW(10,16)=2051.1
      CVDW(10,17)=2051.1
      CVDW(10,18)=380.5
      CVDW(10,19)=380.5
      CVDW(10,20)=380.5
      CVDW(10,21)=380.5
      CVDW(10,22)=380.5
      CVDW(10,23)=380.5
      CVDW(10,24)=380.5
      CVDW(10,25)=380.5
      CVDW(10,26)=380.5
      CVDW(10,27)=380.5
      CVDW(10,28)=380.5
      CVDW(10,29)=380.5
      CVDW(10,30)=380.5
      CVDW(10,31)=380.5
      CVDW(10,33)=380.5
      CVDW(10,34)=380.5
      CVDW(10,35)=380.5
      CVDW(10,36)=981.1
      CVDW(10,41)=CVDW(1,10)*0.8
      CVDW(11,11)=327.2
      CVDW(11,12)=340.0
      CVDW(11,13)=1026.3
      CVDW(11,14)=1026.3
      CVDW(11,15)=981.1
      CVDW(11,17)=684.0
      CVDW(11,18)=125.0
      CVDW(11,19)=125.0
      CVDW(11,20)=125.0
      CVDW(11,21)=125.0
      CVDW(11,22)=125.0
      CVDW(11,23)=1600.
      CVDW(11,24)=125.0
      CVDW(11,25)=125.0
      CVDW(11,26)=125.0
      CVDW(11,27)=125.0
      CVDW(11,28)=125.0
      CVDW(11,29)=125.0
      CVDW(11,30)=125.0
      CVDW(11,31)=125.0
      CVDW(11,32)=1600.
      CVDW(11,33)=125.
      CVDW(11,34)=225.
      CVDW(11,35)=225.
      CVDW(11,36)=327.2
      CVDW(12,12)=354.0
      CVDW(12,13)=1026.3
      CVDW(12,14)=1026.3
      CVDW(12,15)=1020.5
      CVDW(12,16)=711.5
      CVDW(12,17)=711.5
      CVDW(12,18)=132.0
      CVDW(12,19)=132.0
      CVDW(12,20)=132.0
      CVDW(12,21)=132.0
      CVDW(12,22)=132.0
      CVDW(12,23)=132.0
      CVDW(12,24)=132.0
      CVDW(12,25)=132.0
      CVDW(12,26)=132.0
      CVDW(12,27)=132.0
      CVDW(12,28)=132.0
      CVDW(12,29)=132.0
      CVDW(12,30)=132.0
      CVDW(12,31)=132.0
      CVDW(12,33)=132.0
      CVDW(12,34)=132.0
      CVDW(12,35)=132.0
      CVDW(13,13)=1026.3
      CVDW(13,16)=1026.3
      CVDW(13,17)=265.2
      CVDW(13,20)=16.5
      CVDW(13,21)=16.5
      CVDW(13,22)=16.5
      CVDW(13,23)=16.5
      CVDW(13,24)=16.5
      CVDW(13,25)=16.5
      CVDW(13,26)=16.5
      CVDW(13,27)=16.5
      CVDW(13,28)=16.5
      CVDW(13,29)=16.5
      CVDW(13,30)=16.5
      CVDW(13,31)=16.5
      CVDW(13,34)=16.5
      CVDW(13,35)=16.5
      CVDW(13,41)=16.5
      CVDW(14,15)=981.1
      CVDW(14,16)=265.2
      CVDW(14,17)=265.2
      CVDW(14,20)=16.5
      CVDW(14,21)=16.5
      CVDW(14,22)=16.5
      CVDW(14,23)=16.5
      CVDW(14,24)=16.5
      CVDW(14,25)=16.5
      CVDW(14,26)=16.5
      CVDW(14,27)=16.5
      CVDW(14,28)=16.5
      CVDW(14,29)=16.5
      CVDW(14,30)=16.5
      CVDW(14,31)=16.5
      CVDW(14,34)=16.5
      CVDW(14,35)=16.5
      CVDW(15,15)=2942.0
      CVDW(15,20)=380.5
      CVDW(15,21)=380.5
      CVDW(15,22)=380.5
      CVDW(15,23)=380.5
      CVDW(15,24)=380.5
      CVDW(15,25)=380.5
      CVDW(15,26)=380.5
      CVDW(15,27)=1000.
      CVDW(15,28)=380.5
      CVDW(15,29)=380.5
      CVDW(15,30)=380.5
      CVDW(15,31)=380.5
      CVDW(15,34)=380.5
      CVDW(15,35)=380.5
      CVDW(15,36)=380.5
      CVDW(16,16)=1430.0
      CVDW(16,17)=1430.0
      CVDW(16,18)=265.2
      CVDW(16,19)=265.2
      CVDW(16,20)=265.2
      CVDW(16,21)=265.2
      CVDW(16,22)=265.2
      CVDW(16,23)=265.2
      CVDW(16,24)=265.2
      CVDW(16,25)=1220.
      CVDW(16,26)=265.2
      CVDW(16,27)=265.2
      CVDW(16,28)=265.2
      CVDW(16,29)=265.2
      CVDW(16,30)=265.2
      CVDW(16,31)=265.2
      CVDW(16,33)=265.2
      CVDW(16,34)=265.2
      CVDW(16,35)=265.2
      CVDW(17,17)=1430.0
      CVDW(17,18)=265.2
      CVDW(17,19)=265.2
      CVDW(17,20)=265.2
      CVDW(17,21)=265.2
      CVDW(17,22)=265.2
      CVDW(17,23)=265.2
      CVDW(17,24)=265.2
      CVDW(17,25)=265.2
      CVDW(17,26)=265.2
      CVDW(17,27)=265.2
      CVDW(17,28)=265.2
      CVDW(17,29)=265.2
      CVDW(17,30)=265.2
      CVDW(17,31)=265.2
      CVDW(17,33)=265.2
      CVDW(17,34)=265.2
      CVDW(17,35)=265.2
      CVDW(17,41)=CVDW(1,17)*0.8
      CVDW(18,20)=16.5
      CVDW(18,21)=16.5
      CVDW(18,22)=16.5
      CVDW(18,23)=16.5
      CVDW(18,24)=16.5
      CVDW(18,25)=16.5
      CVDW(18,26)=16.5
      CVDW(18,27)=16.5
      CVDW(18,28)=16.5
      CVDW(18,29)=16.5
      CVDW(18,30)=16.5
      CVDW(18,31)=16.5
      CVDW(18,34)=16.5
      CVDW(18,35)=16.5
      CVDW(18,41)=CVDW(1,18)*0.8
      CVDW(19,20)=16.5
      CVDW(19,21)=16.5
      CVDW(19,22)=16.5
      CVDW(19,23)=16.5
      CVDW(19,24)=16.5
      CVDW(19,25)=16.5
      CVDW(19,26)=16.5
      CVDW(19,27)=16.5
      CVDW(19,28)=16.5
      CVDW(19,29)=16.5
      CVDW(19,30)=16.5
      CVDW(19,31)=16.5
      CVDW(19,34)=16.5
      CVDW(19,35)=16.5
      CVDW(19,41)=CVDW(1,19)*0.8
      CVDW(20,20)=16.5
      CVDW(20,21)=16.5
      CVDW(20,22)=16.5
      CVDW(20,23)=16.5
      CVDW(20,24)=16.5
      CVDW(20,25)=300.
      CVDW(20,26)=16.5
      CVDW(20,27)=16.5
      CVDW(20,28)=16.5
      CVDW(20,29)=16.5
      CVDW(20,30)=16.5
      CVDW(20,31)=16.5
      CVDW(20,32)=3800.
      CVDW(20,33)=16.5
      CVDW(20,34)=16.5
      CVDW(20,35)=16.5
      CVDW(20,36)=125.0
      CVDW(20,41)=CVDW(1,20)*0.8
      CVDW(23,36)=350.
      CVDW(34,34)=350.
      CVDW(34,36)=1250.
      CVDW(35,35)=350.
      CVDW(35,36)=1250.
      CVDW(36,36)=327.2
      IVDW(1,1)=12
      IVDW(1,2)=12
      IVDW(1,3)=12
      IVDW(1,4)=12
      IVDW(1,5)=12
      IVDW(1,6)=12
      IVDW(1,7)=6
      IVDW(1,8)=6
      IVDW(1,9)=6
      IVDW(1,10)=6
      IVDW(1,11)=12
      IVDW(1,12)=12
      IVDW(1,13)=6
      IVDW(1,14)=6
      IVDW(1,15)=6
      IVDW(1,16)=6
      IVDW(1,17)=6
      IVDW(1,18)=6
      IVDW(1,19)=6
      IVDW(1,20)=6
      IVDW(1,21)=6
      IVDW(1,22)=6
      IVDW(1,23)=6
      IVDW(1,24)=6
      IVDW(1,25)=6
      IVDW(1,26)=6
      IVDW(1,27)=6
      IVDW(1,28)=6
      IVDW(1,29)=6
      IVDW(1,30)=6
      IVDW(1,31)=6
      IVDW(1,32)=6
      IVDW(1,33)=6
      IVDW(1,34)=6
      IVDW(1,35)=6
      IVDW(1,36)=12
      IVDW(1,41)=12
      IVDW(2,2)=12
      IVDW(2,3)=12
      IVDW(2,4)=12
      IVDW(2,5)=12
      IVDW(2,6)=12
      IVDW(2,7)=6
      IVDW(2,8)=6
      IVDW(2,9)=6
      IVDW(2,10)=6
      IVDW(2,11)=12
      IVDW(2,12)=12
      IVDW(2,13)=6
      IVDW(2,14)=6
      IVDW(2,15)=6
      IVDW(2,16)=6
      IVDW(2,17)=6
      IVDW(2,18)=6
      IVDW(2,19)=6
      IVDW(2,20)=6
      IVDW(2,21)=6
      IVDW(2,22)=6
      IVDW(2,23)=6
      IVDW(2,24)=6
      IVDW(2,25)=6
      IVDW(2,26)=6
      IVDW(2,27)=6
      IVDW(2,28)=6
      IVDW(2,29)=6
      IVDW(2,30)=6
      IVDW(2,31)=6
      IVDW(2,32)=12
      IVDW(2,33)=6
      IVDW(2,34)=12
      IVDW(2,35)=12
      IVDW(2,36)=6
      IVDW(2,41)=12
      IVDW(3,3)=12
      IVDW(3,4)=12
      IVDW(3,5)=12
      IVDW(3,6)=12
      IVDW(3,7)=6
      IVDW(3,8)=6
      IVDW(3,9)=6
      IVDW(3,10)=6
      IVDW(3,11)=12
      IVDW(3,12)=12
      IVDW(3,13)=6
      IVDW(3,14)=6
      IVDW(3,15)=6
      IVDW(3,16)=6
      IVDW(3,17)=6
      IVDW(3,18)=6
      IVDW(3,19)=6
      IVDW(3,20)=6
      IVDW(3,21)=6
      IVDW(3,22)=6
      IVDW(3,23)=6
      IVDW(3,24)=6
      IVDW(3,25)=6
      IVDW(3,26)=6
      IVDW(3,27)=6
      IVDW(3,28)=6
      IVDW(3,29)=6
      IVDW(3,30)=6
      IVDW(3,31)=6
      IVDW(3,33)=6
      IVDW(3,34)=6
      IVDW(3,35)=6
      IVDW(3,36)=12
      IVDW(3,41)=12
      IVDW(4,4)=12
      IVDW(4,5)=12
      IVDW(4,6)=12
      IVDW(4,7)=6
      IVDW(4,8)=6
      IVDW(4,9)=6
      IVDW(4,10)=6
      IVDW(4,11)=12
      IVDW(4,12)=12
      IVDW(4,13)=6
      IVDW(4,14)=6
      IVDW(4,15)=6
      IVDW(4,16)=6
      IVDW(4,17)=6
      IVDW(4,18)=6
      IVDW(4,19)=6
      IVDW(4,20)=6
      IVDW(4,21)=6
      IVDW(4,22)=6
      IVDW(4,23)=6
      IVDW(4,24)=6
      IVDW(4,25)=6
      IVDW(4,26)=12
      IVDW(4,27)=6
      IVDW(4,28)=6
      IVDW(4,29)=6
      IVDW(4,30)=6
      IVDW(4,31)=6
      IVDW(4,33)=6
      IVDW(4,34)=6
      IVDW(4,35)=6
      IVDW(4,36)=12
      IVDW(4,41)=12
      IVDW(5,5)=12
      IVDW(5,6)=12
      IVDW(5,7)=6
      IVDW(5,8)=6
      IVDW(5,9)=6
      IVDW(5,10)=6
      IVDW(5,11)=12
      IVDW(5,12)=12
      IVDW(5,13)=6
      ivdw(5,14)=6
      IVDW(5,15)=6
      IVDW(5,16)=6
      IVDW(5,17)=6
      IVDW(5,18)=6
      IVDW(5,19)=6
      IVDW(5,20)=6
      IVDW(5,21)=6
      IVDW(5,22)=6
      IVDW(5,23)=6
      IVDW(5,24)=6
      IVDW(5,25)=6
      IVDW(5,26)=6
      IVDW(5,27)=6
      IVDW(5,28)=6
      IVDW(5,29)=6
      IVDW(5,30)=6
      IVDW(5,31)=6
      IVDW(5,33)=6
      IVDW(5,34)=6
      IVDW(5,35)=6
      IVDW(5,36)=12
      IVDW(5,41)=12
      IVDW(6,6)=12
      IVDW(6,7)=6
      IVDW(6,8)=6
      IVDW(6,9)=6
      IVDW(6,10)=6
      IVDW(6,11)=12
      IVDW(6,12)=12
      IVDW(6,13)=6
      IVDW(6,14)=6
      IVDW(6,15)=6
      IVDW(6,16)=6
      IVDW(6,17)=6
      IVDW(6,18)=6
      IVDW(6,19)=6
      IVDW(6,20)=6
      IVDW(6,21)=6
      IVDW(6,22)=6
      IVDW(6,23)=6
      IVDW(6,24)=6
      IVDW(6,25)=6
      IVDW(6,26)=6
      IVDW(6,27)=6
      IVDW(6,28)=6
      IVDW(6,29)=6
      IVDW(6,30)=6
      IVDW(6,31)=6
      IVDW(6,32)=6
      IVDW(6,34)=6
      IVDW(6,35)=6
      IVDW(6,35)=6
      IVDW(6,36)=6
      IVDW(6,41)=12
      IVDW(7,11)=6
      IVDW(7,12)=6
      IVDW(7,36)=6
      IVDW(7,41)=6
      IVDW(8,11)=6
      IVDW(8,12)=6
      IVDW(8,36)=6
      IVDW(8,41)=6
      IVDW(9,11)=6
      IVDW(9,12)=6
      IVDW(9,36)=6
      IVDW(9,41)=6
      IVDW(10,11)=6
      IVDW(10,12)=6
      IVDW(10,36)=6
      IVDW(10,41)=6
      IVDW(11,11)=12
      IVDW(11,12)=12
      IVDW(11,13)=6
      IVDW(11,14)=6
      IVDW(11,15)=6
      IVDW(11,16)=6
      IVDW(11,17)=6
      IVDW(11,18)=6
      IVDW(11,19)=6
      IVDW(11,20)=6
      IVDW(11,21)=6
      IVDW(11,22)=6
      IVDW(11,23)=6
      IVDW(11,24)=6
      IVDW(11,25)=6
      IVDW(11,26)=6
      IVDW(11,27)=6
      IVDW(11,28)=6
      IVDW(11,29)=6
      IVDW(11,30)=6
      IVDW(11,31)=6
      IVDW(11,32)=6
      IVDW(11,33)=6
      IVDW(11,34)=6
      IVDW(11,35)=6
      IVDW(11,36)=12
      IVDW(12,12)=12
      IVDW(12,13)=6
      IVDW(12,14)=6
      IVDW(12,16)=6
      IVDW(12,17)=6
      IVDW(12,18)=6
      IVDW(12,19)=6
      IVDW(12,20)=6
      IVDW(12,21)=6
      IVDW(12,22)=6
      IVDW(12,23)=6
      IVDW(12,24)=6
      IVDW(12,25)=6
      IVDW(12,26)=6
      IVDW(12,27)=6
      IVDW(12,28)=6
      IVDW(12,29)=6
      IVDW(12,30)=6
      IVDW(12,31)=6
      IVDW(12,33)=6
      IVDW(12,34)=6
      IVDW(12,35)=6
      IVDW(12,36)=12
      IVDW(13,13)=6
      IVDW(13,15)=6
      IVDW(13,16)=6
      IVDW(13,41)=6
      IVDW(14,15)=6
      IVDW(15,27)=6
      IVDW(15,34)=6
      IVDW(15,35)=6
      IVDW(15,36)=6
      IVDW(16,25)=6
      IVDW(20,27)=6
      IVDW(20,32)=6
      IVDW(20,34)=6
      IVDW(20,35)=6
      IVDW(20,36)=6
      IVDW(20,41)=6
      IVDW(23,36)=6
      IVDW(34,34)=6
      IVDW(34,36)=6
      IVDW(35,36)=6
      IVDW(35,35)=6
      IVDW(36,36)=12
      IVDW(41,41)=12
      PTT(1,1)=0.14
      PTT(1,2)=0.14
      PTT(1,3)=0.14
      PTT(1,4)=0.02
      PTT(1,11)=0.14
      PTT(7,7)=0.02
C
      DO 3 I=1,50
      DO 3 J=I,50
      AL(J,I)=AL(I,J)
      AME(J,I)=AME(I,J)
      DE(J,I)=DE(I,J)
      DF(J,I)=DF(I,J)
      PK2(J,I)=PK2(I,J)
      PK4(J,I)=PK4(I,J)
      PK6(J,I)=PK6(I,J)
      PPIT(J,I)=PPIT(I,J)
      PPIZ(J,I)=PPIZ(I,J)
      PPIx(J,I)=PPIx(I,J)
      AVDW(J,I)=AVDW(I,J)
      BVDW(J,I)=BVDW(I,J)
      CVDW(J,I)=CVDW(I,J)
      IVDW(J,I)=IVDW(I,J)
      BL(J,I)=BL(I,J)
      BK(J,I)=BK(I,J)
      PTT(J,I)=PTT(I,J)
    3 CONTINUE
      RETURN
      END
      SUBROUTINE GEOTES
C
C    THIS SUBROUTINE CHECKS INPUT GEOMETRY
C
      double precision xkc,ykc,zkc,vircx,vircy,vircz,viryx,virzx,virzy
      CHARACTER*2 ATA,ATN,W1,W2
      COMMON /LABEL/ ATA(50)
      COMMON /DIRECT / NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2
     1,NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      PARAMETER(IX=301,IY=IX,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1    IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM2/NUMAT,IATD(IX),NFIX(IX),ILES,NNAT(IX),NAD
     1(IX),NBD(IX),NCD(IX),NDD(IX),NCTOR,MCTOR,NPERM(20),ITOR(IX),
     2NST(IX),NMOL(IX),KTOR
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
cp
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1 cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2 imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3 ,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
      common /pacplo/ npac,ispac(2000),cpac(3,2000),ihhyd(100),
     1 iahyd(100),dhyd(100),nhyd,jhyd(ix),xo(3,ix),xsum(3,ix),
     2 ibhyd(100),hhyd(100),abhyd(100),anghyd(100)   
cp
      DIMENSION IBGR(50),JJ(1000),W1(4),N1(4),W2(4),N2(4),R(4)
C
C    IBGR STORES THE CORRECT NUMBER OF BONDS FOR EACH ATOM SORT
C
      DATA IBGR /4,3,3,2,2,1,2,1,1,1,2,1,2,2,4,3,2,2,2,1,19*1,3,10*1/
      IWARN=0
      ISTOP=0
      naaold=naa
      DO 40 I=1,NAA
      IF(ISORT(I).EQ.12.OR.ISORT(I).EQ.13) THEN
      QADD(I)=2.
C TREAT NA AS CA 1+
      IF (IAT(I).EQ.11) QADD(I)=1.
      GOTO 40
      END IF
C TREAT AL AS SI 1-
      IF (IAT(I).EQ.13) QADD(I)=-1.
C VARIOUS METAL IONS RESIDE BETWEEN 20 AND 30
      IF (ISORT(I).GT.39) GOTO 4711
      IF (ISORT(I).EQ.32) THEN
      QADD(I)=1.
      GOTO 40
      END IF
      IF (ISORT(I).GT.30) GOTO 40 
      IF (ISORT(I).GT.27) THEN
      QADD(I)=4.
      GOTO 40
      END IF
      IF (ISORT(I).GT.25) THEN
      QADD(I)=3.
      GOTO 40
      END IF
      IF (ISORT(I).GT.19) THEN
      QADD(I)=2.
      GOTO 40
      END IF
 4711 IBO=0
      DO 10 J=1,NAA
      IF(VS(I,J).LE.0.) GOTO 10
      IF(ISORT(J).EQ.12.OR.ISORT(J).EQ.13) GOTO 10
      IBO=IBO+1
      IF (IBO.GT.4) GOTO 11
      JJ(IBO)=J
   10 CONTINUE
C
C    IF AN ATOM HAS LESS THAN IBGR BONDS A WARNING IS ISSUED,
C    AND CALCULATION CONTINUES
C    IF AN ATOM HAS MORE THAN IBGR BONDS, THE CALCULATION IS
C    TERMINATED
C
   11 ISO=ISORT(I)+1
      IF(IBO-IBGR(ISO)) 20,40,30
   20 IF (ISO.EQ.18.OR.ISO.EQ.19) THEN
      QADD(I)=1.
      ELSE IF (ISO.EQ.7.OR.NST(NUMM(I)).GE.200) THEN
      QADD(I)=-1.
      DO 31 J=1,NAA
      IF (VS(I,J).LE.0.) GOTO 31
      IF (ISORT(J).EQ.16) QADD(I)=QADD(I)*0.3333
      IF (ISORT(J).EQ.15) QADD(I)=QADD(I)*0.5000
   31 CONTINUE
      ELSE
      IWARN=IWARN+1
      ATN=ATA(ISO)
      WRITE(IUT,60) ATN,NUMM(I),(X(K,I),K=1,3)
      IF (ISO.GT.5.AND.ISO.NE.11.AND.ISO.NE.15) GOTO 32
      IF (ISO.EQ.4) GOTO 32
      CALL SETHYD(JJ,IBGR,I,IBO,ISO,IWARN)
   32 CONTINUE
      END IF
      GOTO 40
   30 IF (ISO.EQ.3) THEN
      QADD(I)=1.
      ELSE IF (ISO.EQ.41) THEN
      QADD(I)=-1.
      ELSE IF (ISO.EQ.17) THEN
      QADD(I)=FLOAT(IBO-IBGR(ISO))*0.3333
      ELSE IF (ISO.EQ.16) THEN
      QADD(I)=FLOAT(IBO-IBGR(ISO))*0.5000
      ELSE
      IWARN=IWARN+1
      ISTOP=ISTOP+1
      ATN=ATA(ISO)
      WRITE(IUT,70) ATN,NUMM(I),(X(K,I),K=1,3)
      END IF
   40 CONTINUE
      DO 41 I=1,NAA
      IF (ISORT(I).NE.16.AND.ISORT(I).NE.15) GOTO 41
      QA=0.
      QSUM=0.
      DO 42 J=1,NAA
      IF (VS(I,J).LE.0.) GOTO 42
      IF (QADD(J).EQ.0.) GO TO 42
      QA=QA+1.
      QSUM=QSUM+QADD(J)
   42 CONTINUE
      DO 43 J=1,NAA
      IF (VS(I,J).LE.0.) GOTO 43
      IF (QADD(J).EQ.0.) GO TO 43
      IF (ISORT(I).EQ.16) QADD(J)=QADD(J)+(QADD(I)+QSUM)*2./QA
      IF (ISORT(I).EQ.15) QADD(J)=QADD(J)+(QADD(I)+QSUM)/QA
   43 CONTINUE
   41 CONTINUE
      if (naa.gt.naaold.and.nm.lt.-4) then
      x(1,naa+1)=cl(1,1)
      x(2,naa+1)=0.
      x(3,naa+1)=0.
      x(1,naa+2)=cl(1,2)
      x(2,naa+2)=cl(2,2)
      x(3,naa+2)=0.
      x(1,naa+3)=cl(1,3)
      x(2,naa+3)=cl(2,3)
      x(3,naa+3)=cl(3,3)
      call crtcry
      end if
      IF(IWARN.EQ.0) RETURN
      IF (ISEQ.EQ.2.AND.NC.LE.0) RETURN
      WRITE(IUT,80)
      NPRINT=0
      DO 50 I=1,NAA
      DO 50 J=I,NAA
      IF(VS(I,J).LE.0.) GOTO 50
      NPRINT=NPRINT+1
      ISO1=ISORT(I)+1
      ISO2=ISORT(J)+1
      W1(NPRINT)=ATA(ISO1)
      W2(NPRINT)=ATA(ISO2)
      N1(NPRINT)=NUMM(I)
      N2(NPRINT)=NUMM(J)
      R(NPRINT)=VS(I,J)
      IF (NPRINT.LT.2) GOTO 50
      WRITE(IUT,90) (W1(M),N1(M),W2(M),N2(M),R(M),M=1,2)
      NPRINT=0
   50 CONTINUE
      WRITE(IUT,90) (W1(M),N1(M),W2(M),N2(M),R(M),M=1,NPRINT)
      RETURN
   60 FORMAT(1H ,'WARNING',7X,A2,I3,3F10.4,' ATOM HAS NOT ENOUGH BONDS'
     1)
   70 FORMAT (1H0,'ERROR',10X,A2,I3,3F10.4,' ATOM HAS TOO MANY BONDS')
   80 FORMAT (1H0,'BOND LENGTHS')
   90 FORMAT(' ',4(A2,I3,A2,I3,F8.3,2X))
      END
      SUBROUTINE PLUPLO 
C
C    THIS SUBROUTINE CREATES THE INPUT FOR THE PLOT PROGRAM PLUTO
C
      DOUBLE PRECISION XK,YK,ZK,virx,viry,virz
      CHARACTER*2 ATA,AS
      CHARACTER*5 RESTY,RESTYP,restypc
      CHARACTER*4 RESNA,RESNAM,resnamc
C      CHARACTER*14 FMT
      COMMON /LABEL/ATA(50)
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      PARAMETER(IX=301,IY=IX,IX3=IX*3,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1    IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /PDBNAM/ RESNA(IX),NURES(IX),RESNAM(IX),NUMRES(IX)
     1,RESTY(IX),RESTYP(IX)
      COMMON /OPTIM/virx,viry,virz,XD(IX3),XDD(IX3),QZ(IX3),DX(IX3),
     1XK(IX),YK(IX),ZK(IX),NAM,EPP,EPW
cp    
      common /pacplo/ npac,ispac(2000),cpac(3,2000),ihhyd(100),
     1 iahyd(100),dhyd(100),nhyd,jhyd(ix),xo(3,ix),xsum(3,ix),
     2 ibhyd(100),hhyd(100),abhyd(100),anghyd(100)   
cp
cs
      common /solvat/ eps,dielec,dsgl,fav(ix),fam(ix),rxy(ix,ix),
     1 rvdw(50),ion(50)
      common /dyna/ idz,idd(4,40),irdz(ix6),idyn,dystar
cs
      iset=iser
      if (ng.eq.4) then
      iset=iflag2    
      EPOT=EPP+EPW
      dhfm=96.49*EPOT+dsgl*4.184
      end if
      WRITE(IPDB,14) ISET,DHFM,(TITLE(I),I=1,8)
   14 FORMAT (7HHEADER ,I6,F8.1,2X,8A4)
      DO 100 I=1,NAA
      ISO=ISORT(I)+1
      AS=ATA(ISO)
       IF (RESTYP(I).EQ.'     ')RESTYP(I)=ATA(ISO)
       IF (RESNAM(I).EQ.'    ')RESNAM(I)='PIM '
       IF (NUMRES(I).EQ.0)NUMRES(I)=IMOL(I)+1
       WRITE(IPDB,5)I,RESTYP(I),RESNAM(I),NUMRES(I),(X(J,I),J=1,3)
    5 FORMAT(4HATOM,4X,I3,1X,A5,A4,2X,I3,4X,3F8.3,2X,
     110H0.00 00.00)   
  100 CONTINUE
cp
      if (nm.lt.-5.and.nd.gt.0) then
      call outpac
      if (npac.gt.0) then
      do 1000 ipac=1,npac
      nii=naa+ipac
      as=ata(ispac(ipac))        
      RESTYPc=as
      RESNAMc='PIM '
      NUMRESc=1
      WRITE(IPDB,5)nIi,RESTYPc,RESNAMc,NUMRESc,(cpac(j,ipac),j=1,3)
 1000 continue
      end if
      restypc='Au'
      anull=0.
      WRITE(IPDB,5)nIi,RESTYPc,RESNAMc,NUMRESc,anull,anull,anull
      do 1001 jj=1,3
      nii=nii+jj
 1001 WRITE(IPDB,5)nIi,RESTYPc,RESNAMc,NUMRESc,x(1,naa+jj),x(2,naa+jj),
     1x(3,naa+jj) 
      end if
      WRITE(IPDB,12)
   12 FORMAT ('END')
      RETURN
      END
      SUBROUTINE TORVAR(*)
C
C    THIS SUBROUTINE PERMUTATES TORSIONAL ANGLES SELECTED BY ITOR
C
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      PARAMETER(IX=301)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1    IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM2/NUMAT,IATD(IX),NFIX(IX),ILES,NNAT(IX),NAD
     1(IX),NBD(IX),NCD(IX),NDD(IX),NCTOR,MCTOR,NPERM(20),ITOR(IX),
     2NST(IX),NMOL(IX),KTOR
      COMMON /DYN/ TIME,TEMP,GMASS,MDCYC,TEMITT,EMITT,TMASS,MDUMP,MCOOL,
     1SCHW(3),VELO(3,IX),ACC(3,IX),X3(3,IX),TSOLL,IRAND,TE,iav
     2,pnull,taupe,taupi,pmitt,taute,pstart,pdiff,pend,mpress
      DIMENSION ATOR(6,5),BTOR(6,10)
      DATA ATOR / 0.   , 3.141592654, 0. , 0. ,0. , 0. ,
     1            1.047197551, 3.141592654, 5.235987757, 0. , 0. , 0. ,
     2            0. , 1.570796327, 3.141592654, 4.712388981, 0. , 0. ,
     3            0. ,0.785398164, 1.570796327, 2.35619449,
     4            3.141592654, 0. ,
     5            0. ,1.047197551, 2.094395102, 3.141592654,
     6            4.188790205, 5.235987757/
      NCTOR=NCTOR+1
C
C    SET VALUES OF THE TORSIONAL ANGLES TO INITIAL VALUES AND
C    CALCULATE THE NUMBER OF PERMUTATIONS
C
      IF (NCTOR.GT.1) GOTO 20
      IF (INFORM.EQ.0) READ(IIN,70) KTOR
      IF (INFORM.EQ.1) READ(IIN,*) KTOR
      IF (KTOR.EQ.0) GOTO 12
C
C    INPUT OF ANGLE SET FOR TORSION ANGLE PERMUTATION.
C
      DO 11 I=1,KTOR
      IF (INFORM.EQ.0) READ(IIN,80) (BTOR(J,I),J=1,6)
      IF (INFORM.EQ.1) READ(IIN,*) (BTOR(J,I),J=1,6)
   11 CONTINUE
   12 IPERM=1
      DO 10 I=1,NAA
         IF (ITOR(I).EQ.0) GOTO 10
         ITO=ITOR(I)
         IF(ITO.EQ.6) GOTO 15
         IF (KTOR.EQ.0) GOTO 13
         A(3,I)=-BTOR(1,IPERM)*0.0174532925
         GOTO 14
   13    A(3,I)=-ATOR(1,ITO)
         GOTO 14
   15    A(3,I)=AMOD(RANDOM(),180.)*0.0174532925
   14    NPERM(IPERM)=1
         IPERM=IPERM+1
         MCTOR=MCTOR*(ITO+1)
   10 CONTINUE
      GOTO 50
C
C   SET TORSIONAL ANGLES FOR EACH PERMUTATION
C
   20 IF (NCTOR.GT.MCTOR) GOTO 60
      IPERM=1
      DO 40 I=1,NAA
         IF (ITOR(I).EQ.0) GOTO 40
         ITO=ITOR(I)
         NTO=NPERM(IPERM)+1
         IF (NTO.GT.(ITO+1)) GOTO 30
         IF(ITO.EQ.6) GOTO 33
         IF (KTOR.EQ.0) GOTO 21
         A(3,I)=-BTOR(NTO,IPERM)*0.0174532925
         GOTO 22
   21    A(3,I)=-ATOR(NTO,ITO)
         GOTO 22
   22    NPERM(IPERM)=NPERM(IPERM)+1
         GOTO 50
   30    NPERM(IPERM)=1
         IF(ITO.EQ.6) GOTO 33
         IF (KTOR.EQ.0) GOTO 31
         A(3,I)=-BTOR(1,IPERM)*0.0174532925
         GOTO 32
   31    A(3,I)=-ATOR(1,ITO)
         GOTO 32
   33    A(3,I)=AMOD(RANDOM(),180.)*0.0174532925
   32    IPERM=IPERM+1
   40 CONTINUE
   50 WRITE (IUT,90)
      DO 51 I=1,NAA
         IF (ITOR(I).EQ.0) GOTO 51
         WW=-A(3,I)*57.295
         WRITE(IUT,95) NCD(I),NBD(I),NAD(I),I,WW,ITOR(I)
   51 CONTINUE
      RETURN
C
C    RESET IF ALL PERMUTATIONS ARE CALCULATED
C
   60 WRITE(IUT,100)
   70 FORMAT (I3)
   80 FORMAT (10X,6F10.4)
   90 FORMAT (/)
   95 FORMAT (1H ,'PERMUTATED ANGLE',4I5,F8.2,' FOR ITOR=',I3)
  100 FORMAT(1H0,10X,'END OF CALCULATION')
      ISEQ=0
      NCTOR=0
      MCTOR=1
      DO 200 I=1,IPERM
      NPERM(I)=0
  200 CONTINUE
      IPERM=0
      RETURN 1
      END
      SUBROUTINE BETOV
C
C    THIS SUBROUTINE CALCULATES RESONANCE INTEGRALS BETA(I,J)
C
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      PARAMETER(IX=301,IY=IX,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /BETA/ H(IY,IY)
      DIMENSION E(3),NP(6),TA(3),TB(3),TR(3,3),CA(3),CB(3)
C
      NAM=NA-1
      DO 120 I=1,NAM
      K=ISORT(I)
      I1=I+1
      DO 120 J=I1,NA
C
C    FIND ALL NEIGHBOURS OF THE ATOMS I AND J
C
      DO 1 N=1,6
    1 NP(N)=0
      MQ=1
      NQ=4
      DO 20 N=1,NAA
      IF (VS(I,N).LE.0.) GOTO 10
      NP(MQ)=N
      MQ=MQ+1
 10   IF (VS(J,N).LE.0.) GOTO 20
      NP(NQ)=N
      NQ=NQ+1
 20   CONTINUE
   21 N1=NP(1)
      N2=NP(2)
      N3=NP(3)
      IF (N3.NE.0) GOTO 25
C
C    IF AN ATOM HAS ONLY ONE NEIGHBOUR, LOOK FOR A NEXT NEAREST ATOM
C
      IF (N2) 22,22,24
   22 KK=N1
      DO 220 N=1,3
  220 NP(N)=0
      MQ=1
      DO 23 N=1,NAA
      IF (VS(KK,N).LE.0.) GOTO 23
      NP(MQ)=N
      MQ=MQ+1
   23 CONTINUE
      IF (NP(3).EQ.0) NP(3)=KK
      GOTO 21
   24 N3=I
   25 N4=NP(4)
      N5=NP(5)
      N6=NP(6)
CMK FUER CO      
      IF (N2.EQ.0) GOTO 70
      IF (N6.NE.0) GOTO 29
      IF (N5) 26,26,28
   26 KK=N4
      DO 260 N=4,6
  260 NP(N)=0
      NQ=4
      DO 27 N=1,NAA
      IF (VS(KK,N).LE.0.) GOTO 27
      NP(NQ)=N
      NQ=NQ+1
   27 CONTINUE
      IF (NP(6).EQ.0) NP(6)=KK
      GOTO 25
   28 N6=J
C
C    CALCULATE NORMAL VECTORS OF THE PLANES FORMED BY THE ATOMS
C    SURROUNDING ATOMS I AND J
C
   29 TA(1)=(X(2,N1)-X(2,N3))*(X(3,N1)-X(3,N2))-(X(3,N1)-X(3,N3))*
     1(X(2,N1)-X(2,N2))
      TA(2)=(X(3,N1)-X(3,N3))*(X(1,N1)-X(1,N2))-(X(1,N1)-X(1,N3))*
     1(X(3,N1)-X(3,N2))
      TA(3)=(X(1,N1)-X(1,N3))*(X(2,N1)-X(2,N2))-(X(2,N1)-X(2,N3))*
     1(X(1,N1)-X(1,N2))
      TB(1)=(X(2,N4)-X(2,N6))*(X(3,N4)-X(3,N5))-(X(3,N4)-X(3,N6))*
     1(X(2,N4)-X(2,N5))
      TB(2)=(X(3,N4)-X(3,N6))*(X(1,N4)-X(1,N5))-(X(1,N4)-X(1,N6))*
     1(X(3,N4)-X(3,N5))
      TB(3)=(X(1,N4)-X(1,N6))*(X(2,N4)-X(2,N5))-(X(2,N4)-X(2,N6))*
     1(X(1,N4)-X(1,N5))
      TN=SQRT(TA(1)*TA(1)+TA(2)*TA(2)+TA(3)*TA(3))
      IF(TN.EQ.0.)TN=1.E-2
      DO 30 N=1,3
 30   TA(N)=TA(N)/TN
      TN=SQRT(TB(1)*TB(1)+TB(2)*TB(2)+TB(3)*TB(3))
      IF(TN.EQ.0.)TN=1.E-2
      DO 40 N=1,3
 40   TB(N)=TB(N)/TN
      DO 50 N=1,3
 50   E(N)=(X(N,J)-X(N,I))/VV(I,J)
C
C    CALCULATE ANGLE BETWEEN NORMAL VECTORS
C
      COST=E(3)
      IF ((1.-COST**2).GT.0.01) GOTO 60
      SINT=0.
      GOTO 70
 60   SINT=SQRT(1.-COST**2)
 70   CONTINUE
      IF (SINT.GT.0.01) GOTO 80
C
C    SEPARATE OVERLAP BETWEEN ATOMS I AND J INTO COMPONENTS
C    (2P SIGMA/2P SIGMA) AND (2P PI/2P PI)
C
      COSP=1.
      SINP=0.
      GOTO 90
 80   COSP=E(1)/SINT
      SINP=E(2)/SINT
 90   TR(1,1)=SINT*COSP
      TR(1,2)=COST*COSP
      TR(1,3)=-SINP
      TR(2,1)=SINT*SINP
      TR(2,2)=COST*SINP
      TR(2,3)=COSP
      TR(3,1)=COST
      TR(3,2)=-SINT
      TR(3,3)=0.
      DO 100 N=1,3
      CA(N)=0.
      CB(N)=0.
      DO 100 M=1,3
      CA(N)=CA(N)+TA(M)*TR(M,N)
 100  CB(N)=CB(N)+TB(M)*TR(M,N)
      TCA=SQRT(CA(1)*CA(1)+CA(2)*CA(2)+CA(3)*CA(3))
      TCB=SQRT(CB(1)*CB(1)+CB(2)*CB(2)+CB(3)*CB(3))
      IF(TCA.EQ.0.)TCA=1.
      IF(TCB.EQ.0.)TCB=1.
      DO 101 N=1,3
      CA(N)=CA(N)/TCA
 101  CB(N)=CB(N)/TCB
      L=ISORT(J)
C
C    EVALUATE OVERLAP INTEGRAL (2P PI/2P PI)
C
      K01=K+1
      L01=L+1
      Z1=ZS(K)*SQRT(1.+(QSIG(I)-QQ(K01))*RV(K01)/HV(K))
      Z2=ZS(L)*SQRT(1.+(QSIG(J)-QQ(L01))*RV(L01)/HV(L))
      CALL OVLAP(VV(I,J),Z1,Z2,OVLP,OVLPS)
      SPI=(CA(2)*CB(2)+CA(3)*CB(3))*OVLP
      IF (ISORT(I).EQ.10.AND.ISORT(J).EQ.10)
     1SPI=OVLP
      IF (ISORT(I).EQ.10.AND.ISORT(J).EQ.11)
     1SPI=OVLP
      IF (ISORT(I).EQ.11.AND.ISORT(J).EQ.10)
     1SPI=OVLP
      IF (ISORT(I).EQ.35.OR.ISORT(J).EQ.35)
     1SPI=OVLP
      IF (SPI.GE.0.) GOTO 110
      SPI=-SPI
      CB(1)=-CB(1)
C
C    EVALUATE OVERLAP INTEGRAL (2P SIGMA/2P SIGMA)
C
 110  SSIG=-CA(1)*CB(1)*OVLPS
C
C    CALCULATE BETA(I,J) FROM OVERLAP INTEGRALS AND PARAMETER BK.
C    IF A PI BOND EXISTS BETWEEN I AND J BETA IS CALCULATED BY USING
C    BOTH OVERLAP INTEGRALS.
C    IF NO PI BOND EXISTS BETWEEN I AND J BETA IS CALCULATED BY USING
C    (2P SIGMA/2P SIGMA) ONLY TO ACCOUNT FOR HOMOKONJUGATIVE
C    EFFECTS.
C
      BETA=-BK(K,L)*2.*SSIG
      IF (K.EQ.35.OR.L.EQ.35) BETA=0.
      IF (VS(I,J).GT.0.) BETA=-BK(K,L)*(SPI+2.*SSIG)
      H(I,J)=BETA
      H(J,I)=BETA
  120 CONTINUE
      RETURN
      END
      SUBROUTINE OVLAP(RR,ZA,ZB,OVLP,OVLPS)
C
C    THIS SUBROUTINE EVALUATES OVERLAP INTEGRALS (2P SIGMA/2P SIGMA)
C    AND (2P PI/2P PI) USING FUNCTIONS OA AND OB
C
      DIMENSION A(5),B(5)
      SA=0.5*ZA
      SB=0.5*ZB
      RAB=RR/0.529167
      ALPHA=0.5*RAB*(SA+SB)
      BETA=0.5*RAB*(SB-SA)
      DO 3 I=1,5
      N=I-1
      A(I)=OA(ALPHA,N)
  3   B(I)=OB(BETA,N)
      W=-SQRT((SA*SB)**5)*(RAB**5)*0.0625
      OVLP=0.5*W*(A(5)*(B(1)-B(3))-B(5)*(A(1)-A(3))-A(3)*B(1)+B(3)*A(1))
      OVLPS=W*(B(3)*(A(5)+A(1))-A(3)*(B(5)+B(1)))
      RETURN
      END
      FUNCTION OA(A,K)
C
C     THIS FUNCTION CALCULATES THE  A  INTEGRALS ASSOCIATED WITH THE
C     OVERLAP ROUTINE LISTED BELOW.
C
      B=1.0/A
      S=1.0
      OA=1.0
      DO 1 M=1,K
      L=K-M+1
      S=L*S*B
    1 OA=OA+S
      OA=OA*B*EXP(-A)
      RETURN
C
      END
      FUNCTION OB(BETA,N)
C
C     THIS FUNCTION CALCULATES THE B INTEGRALS ASSOCIATED WITH THE
C     OVERLAP ROUTINE LISTED BELOW.
C
      B=BETA**2
      FN=N
      IF (MOD(N,2)) 1,2,1
    1 FNUMER=FN+2.0
      SUM=BETA/FNUMER
      FACTOR=-2.0
      FI=3.0
      GOTO 3
    2 FNUMER=FN+1.
      SUM=1./FNUMER
      FACTOR=2.0
      FI=2.0
    3 T=SUM
    4 DENOM=FNUMER+2.0
      T=T/FI*B/(FI-1.0)*FNUMER/DENOM
      TSUM=SUM+T
      IF (SUM-TSUM) 5,6,5
    6 OB=FACTOR*SUM
      RETURN
    5 SUM=TSUM
      FI=FI+2.0
      FNUMER=DENOM
      GOTO 4
      END
      SUBROUTINE THOQR
C
C    THIS IS THE SUBROUTINE HOQR, WHICH FIRST TRANSFORMS MATRIX A
C    INTO A TRIDIAGONAL MATRIX AND DIAGONALIZES A BY RUTISHAUSER
C    Q/R - TRANSFORMATION USING THE ADDITIONAL SUBROUTINE DCQR.
C    SOME CHANGES HAVE BEEN MADE TO FIT THIS ROUTINE TO THE PURPOSE
C    OF THIS PROGRAM.
C
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      PARAMETER(IX=301,IY=IX)
      COMMON /DIAG/ F(IY,IY),E(IY),V(IY,IY),ISENT,DBETA(IY)
      DIMENSION A(IY,IY)
C
C     DIAGONALIZATION IS PERFORMED ON A(I,J) WHILE THE 
C     ORIGINAL MATRIX IS RETAINED IN F(I,J).
C     ON RETURN E HOLDS THE EIGENVALUES AND V THE EIGENVECTORS OF A
C
      DO 10 I=1,NA
      DO 10 J=I,NA
         A(I,J)=F(I,J)
         A(J,I)=A(I,J)
   10 CONTINUE
C
      V(1,1)=1.
      JCOL=NA
      DO 140 J=1,JCOL
  140 E(J)=A(J,J)
  150 IF (JCOL-2) 290,290,160
  160 SF=ABS(A(1,JCOL))
      JCOL1=JCOL-1
      DO 170 II=2,JCOL1
         X1=AMAX1(SF,ABS(A(II,JCOL)))
         SF=X1
  170 CONTINUE
      IF (ABS(SF).LE.1.E-8) GOTO 280
  180 XSUM=0.
      DO 190 I2=3,JCOL
         T=A(I2-2,JCOL)/SF
         A(I2-2,JCOL)=T
  190 XSUM=XSUM+T*T
      T=A(JCOL-1,JCOL)/SF
      XTOT=SQRT(XSUM+T*T)
      DBETA(JCOL)=SIGN(XTOT*SF,T)
      T=T+SIGN(XTOT,T)
      A(JCOL-1,JCOL)=T
      T=SQRT(0.5*(XSUM+T*T))
      DO 200 I1=2,JCOL
  200 A(I1-1,JCOL)=A(I1-1,JCOL)/T
      XSUM=0.
      DO 210 I1=2,JCOL
  210 XSUM=XSUM+A(I1-1,JCOL)*A(1,I1-1)
      DBETA(1)=XSUM
      DO 240 J1=3,JCOL
         XSUM=0.
         DO 220 I2=3,J1
  220    XSUM=XSUM+A(I2-2,JCOL)*A(I2-2,J1-1)
         DO 230 I1=J1,JCOL
  230    XSUM=XSUM+A(I1-1,JCOL)*A(J1-1,I1-1)
  240 DBETA(J1-1)=XSUM
      XSUM=0.
      DO 250 I1=2,JCOL
  250 XSUM=XSUM+DBETA(I1-1)*A(I1-1,JCOL)
      WAW=XSUM
      DO 260 J1=2,JCOL
      DO 260 I1=2,J1
         WI=A(I1-1,JCOL)
         WJ=A(J1-1,JCOL)
  260 A(I1-1,J1-1)=WI*WAW*WJ-DBETA(I1-1)*WJ-DBETA(J1-1)*WI+A(I1-1,J1-1)       
  270 JCOL=JCOL-1
      GOTO 150
  280 DBETA(JCOL)=0.
      GOTO 270
  290 DBETA(2)=A(1,2)
      DO 300 I=1,NA
         T=A(I,I)
         A(I,I)=E(I)
         IF (I.EQ.NA) GOTO 300
         DBETA(I)=DBETA(I+1)
  300 E(I)=T
      CALL DCQR
      IF (NA-2) 340,340,310
  310 DO 330 JCOL=3,NA
      DO 330 J=1,NA
         XSUM=0.
         DO 320 I1=2,JCOL
  320    XSUM=A(I1-1,JCOL)*V(I1-1,J)+XSUM
         WAW=-XSUM
      DO 330 I1=2,JCOL
  330 V(I1-1,J)=-(V(I1-1,J)+WAW*A(I1-1,JCOL))
  340 CONTINUE
C
C    THE EIGENVALUES AND EIGENVECTORS ARE SORTED ACCORDING TO
C     THE VALUE OF (ISENT) UPWARDS OR DOWNWARDS.
C
   20 NM1=NA-1
      DO 80 I=1,NM1
         IP1=I+1
         DO 70 J=IP1,NA
            IF (ISENT) 40,30,40
   30       IF (E(I)-E(J)) 50,50,70
   40       IF (E(I)-E(J)) 70,70,50
   50       ETEMP=E(I)
            E(I)=E(J)
            E(J)=ETEMP
            DO 60 K=1,NA
               VTEMP=V(K,I)
               V(K,I)=V(K,J)
               V(K,J)=VTEMP
   60       CONTINUE
   70    CONTINUE
   80 CONTINUE
      IF (V(1,1)) 90,110,110
   90 DO 100 I=1,NA
      DO 100 J=I,NA
         TEMP=-V(I,J)
         V(I,J)=-V(J,I)
         V(J,I)=TEMP
  100 CONTINUE
      RETURN
  110 DO 120 I=1,NM1
         IP1=I+1
      DO 120 J=IP1,NA
         TEMP=V(I,J)
         V(I,J)=V(J,I)
         V(J,I)=TEMP
  120 CONTINUE
      RETURN
C
      END
      SUBROUTINE DCQR
C
C    THIS SUBROUTINE IS USED FOR MATRIX DIAGONALIZATION BY RUTISHAUSER
C    Q/R-TRANSFORMATION
C
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      PARAMETER(IX=301,IY=IX)
      COMMON /DIAG/ F(IY,IY),E(IY),V(IY,IY),ISENT,DBETA(IY)
      M=NA
      icount=0
      DBETA(M)=0.
      DO 20 J=1,M
         DO 10 I=1,M
   10    V(I,J)=0.
   20 V(J,J)=1.
      SC=AMAX1(ABS(E(1)),ABS(DBETA(1)))
      DO 30 II=2,M
         X1=AMAX1(SC,ABS(E(II)))
         Y1=AMAX1(X1,ABS(DBETA(II)))
         SC=Y1
   30 CONTINUE
      DO 40 I=1,M
         E(I)=E(I)/SC
   40 DBETA(I)=DBETA(I)/SC
   50 IF (M) 60,60,70
   60 RETURN
   70 DBETA(M)=0.
      M1=M-1
      IF (M1) 90,80,90
   80 E(M)=E(M)*SC
      ICOUNT=0
      M=M1
      GOTO 50
   90 K=M
  100 K=K-1
      IF (K) 120,120,110
  110 IF (ABS(DBETA(K))-1.E-8) 130,130,100
  120 K=0
      ICOUNT=ICOUNT+1
      IF (ICOUNT.GT.10000) THEN
      WRITE (6,999)
  999 FORMAT(' WARNING - FOCK MATRIX APPEARS TO BE SINGULAR')
      GOTO 80
      END IF
      GOTO 150
  130 IF (K-M1) 140,80,140
  140 DBETA(K)=0.
  150 K=K+1
      C=1.
      S=E(M)-E(M-1)
      IF (ABS((ABS(DBETA(M-1))+ABS(S))-ABS(DBETA(M-1))).LE.1.E-6) 
     1   GOTO 170
      C=(2.*DBETA(M-1)/S)/(1.+SQRT(1.+(2.*DBETA(M-1)/S)**2))
  170 H=E(K)-(E(M)+C*DBETA(M-1))
      W=DBETA(K)
      DO 210 I=K,M1
         T=SQRT(H*H+W*W)
         IF (I-K) 190,190,180
  180    DBETA(I-1)=T
  190    C=H/T
         S=W/T
         FF=C*E(I)+S*DBETA(I)
         Q=C*DBETA(I)+S*E(I+1)
         W=S*DBETA(I+1)
         DBETA(I+1)=C*DBETA(I+1)
         EINEW=FF*C+Q*S
         H=-FF*S+Q*C
         E(I+1)=E(I)+E(I+1)-EINEW
         E(I)=EINEW
         DO 200 J=1,NA
            P=V(J,I)
            Q=V(J,I+1)
            V(J,I)=C*P+S*Q
  200    V(J,I+1)=-S*P+C*Q
  210 CONTINUE
      DBETA(M-1)=H
      GOTO 90
C
      END
      SUBROUTINE DEWGEO
C
C    THIS SUBROUTINE TRANSFORMS INTERNAL COORDINATES TO CARTESIAN
C    COORDINATES.
C    IT USES PARTS OF THE MINDO/2' PROGRAM (M.J.S DEWAR AT AL.)
C
      double precision xkc,ykc,zkc,vircx,vircy,vircz,viryx,virzx,virzy
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      PARAMETER(IX=301)
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM2/NUMAT,IATD(IX),NFIX(IX),ILES,NNAT(IX),NAD
     1(IX),NBD(IX),NCD(IX),NDD(IX),NCTOR,MCTOR,NPERM(20),ITOR(IX),
     2NST(IX),NMOL(IX),KTOR
      COMMON /PDBNAM/RESNA(IX),NURES(IX),RESNAM(IX),NUMRES(IX)   
     1,RESTY(IX),RESTYP(IX)    
cp
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1 cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2 imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3 ,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
      common /pacplo/ npac,ispac(2000),cpac(3,2000),ihhyd(100),
     1 iahyd(100),dhyd(100),nhyd,jhyd(ix),xo(3,ix),xsum(3,ix),
     2 ibhyd(100),hhyd(100),abhyd(100),anghyd(100)   
cp
C ARRAYS X1,Y1,Z1 HOLD THE INITIAL COORDINATES OF THE UNSORTED ATOMS     
      DIMENSION X1(IX),Y1(IX),Z1(IX)
      LOGICAL LINEAR
      CHARACTER*4 RESNA,RESNAM
      CHARACTER*5 RESTY,RESTYP
C
C    CLEAR SOME ARRAYS
C
      DO 10 I=1,IX
      NNAT(I)=0
      NFIX(I)=0
      NMOL(I)=0
      NURES(I)=0
      NUMRES(I)=0
      RESNA(I)='    '
      RESNAM(I)='    '
      RESTY(I)='     '
      RESTYP(I)='     '
      DO 10 J=1,3
      DO 10 K=1,4
   10 NRA(K,J,I)=0
      IF (NM.GT.0) GOTO 19
cp
      if (nm.lt.-4) call incrys
cp
C
C     READ CARTESIAN COORDINATES
C
      J=0
      I=0
   11 I=I+1
      IF (INFORM.EQ.0) THEN
       READ(IIN,260)
     1 IATD(I),NST(I),NFIX(I),NMOL(I),X1(I),Y1(I),Z1(I),NAD(I),NBD(I),       
     2 NCD(I),NDD(I),RESTY(I),RESNA(I),NURES(I)
         ELSE
       READ(IIN,*)
     1 IATD(I),NST(I),NFIX(I),NMOL(I),X1(I),Y1(I),Z1(I),NAD(I),NBD(I),
     2 NCD(I),NDD(I)
      END IF
      IF (IATD(I).GT.0) GOTO 11
      NUMAT=I-1
      IF (NM.EQ.0.OR.NM.le.-3) GOTO 190
      I=I-1
C
C    ENTER ROUTINE TRANSF FOR COORDINATE-TRANSFORMATION OF SECOND SET
C
      CALL TRANSF(INFORM,IIN,J,NG,NM,NMOL,NUMAT,X,A,X1,Y1,Z1,*11,*190)
   19 I=0
   20 I=I+1
C
C    READ INTERNAL COORDINATES
C
      IF (ILES.EQ.0) THEN
         IF (INFORM.EQ.0) READ(IIN,250)IATD(I),NST(I),NFIX(I),NMOL(I),
     1                    (A(J,I),J=1,3),ITOR(I),NAD(I),NBD(I),NCD(I)
         IF (INFORM.EQ.1) READ(IIN,*)  IATD(I),NST(I),NFIX(I),NMOL(I),
     1                    (A(J,I),J=1,3),ITOR(I),NAD(I),NBD(I),NCD(I)
      END IF
C
C    CHECK THE INPUT DATA
C
      IF (IATD(I).LE.0) GOTO 40
      IF (NAD(I).GE.I.OR.NBD(I).GE.I.OR.NCD(I).GE.I) GOTO 100
C
C    CONVERT THE ANGLES TO RADIANS
C
      IF (ILES.EQ.0) THEN
         A(2,I)=A(2,I)*0.01745329252
         A(3,I)=-A(3,I)*0.01745329252
      END IF
      AB(1,I)=A(1,I)
      AB(2,I)=A(2,I)
      AB(3,I)=-A(3,I)
      GOTO 20
c   40 continue
   40 NAD(2)=1
      NAD(3)=2
      NBD(3)=1
      NUMAT=I-1
      IF (ISEQ.EQ.2.AND.NC.LE.0) GOTO 110
      IF (NC.LT.0) GOTO 110
      WRITE(IUT,310)
C
C    OUTPUT OF ALL INTERNAL COODINATES
C
      WRITE(IUT,330)
      WRITE(IUT,340) IATD(1)
      WRITE(IUT,350) IATD(2),A(1,2),NAD(2)
      W=A(2,3)*57.29577951
      WRITE(IUT,360) IATD(3),A(1,3),W,NAD(3),NBD(3)
      IF (NUMAT.LT.4) GOTO 110
      DO 90 I=4,NUMAT
         W=A(2,I)*57.29577951
         WW=-A(3,I)*57.29577951
   90 WRITE(IUT,370) I,IATD(I),A(1,I),W,WW,ITOR(I),NAD(I),NBD(I),NCD(I)
  110 LINEAR=.FALSE.
C
C    FIND FIRST ATOM (IF ANY) WHICH IS NOT COLLINEAR WITH ATOMS 1 AND 2
C
      DO 120 I=3,NUMAT
         CCOS=COS(A(2,I))
         IF (ABS(CCOS).LT.0.99999) GOTO 130
  120 CONTINUE
C
C    THE MOLECULE IS LINEAR
C
      LINEAR=.TRUE.
      I=NUMAT+1
  130 X1(I-1)=A(1,I-1)
      Y1(I-1)=0.
      Z1(I-1)=0.
      II=I-2
      DO 140 J=1,II
         K=II-J+1
         X1(K)=X1(K+1)-A(1,K+1)
         Y1(K)=0.
         Z1(K)=0.
  140 CONTINUE
      IF (LINEAR) GOTO 190
      SSIN=SIN(A(2,I))
      X1(I)=X1(I-1)-A(1,I)*CCOS
      Y1(I)=A(1,I)*SSIN
      Z1(I)=0.
      IF (I.EQ.NUMAT) GOTO 190
      II=I+1
      DO 180 I=II,NUMAT
         COSA=COS(A(2,I))
         MB=NBD(I)
         MC=NAD(I)
         XB=X1(MB)-X1(MC)
         YB=Y1(MB)-Y1(MC)
         ZB=Z1(MB)-Z1(MC)
         RBC=1./SQRT(XB*XB+YB*YB+ZB*ZB)
         IF (ABS(COSA).LT.0.99999) GOTO 150
C
C    ATOMS MC, MB, AND (I) ARE COLLINEAR
C
         RBC=A(1,I)*RBC
         X1(I)=X1(MC)-XB*RBC
         Y1(I)=Y1(MC)-YB*RBC
         Z1(I)=Z1(MC)-ZB*RBC
         GOTO 180
C
C    THE ATOMS ARE NOT COLLINEAR
C
  150    MA=NCD(I)
         XA=X1(MA)-X1(MC)
         YA=Y1(MA)-Y1(MC)
         ZA=Z1(MA)-Z1(MC)
C
C    ROTATE ABOUT THE Z-AXIS TO MAKE YB=0, AND XB POSITIVE.
C    IF XYB IS TOO SMALL, FIRST ROTATE THE Y-AXIS BY 90 DEGREES.
C
         XYB=SQRT(XB*XB+YB*YB)
         K=-1
         IF (XYB.GT.0.1) GOTO 160
         XPA=ZA
         ZA=-XA
         XA=XPA
         XPB=ZB
         ZB=-XB
         XB=XPB
         XYB=SQRT(XB*XB+YB*YB)
         K=1
C
C    ROTATE ABOUT THE Y-AXIS TO MAKE ZB VANISH
C
  160    COSTH=XB/XYB
         SINTH=YB/XYB
         XPA=XA*COSTH+YA*SINTH
         YPA=YA*COSTH-XA*SINTH
         SINPH=ZB*RBC
         COSPH=SQRT(ABS(1.-SINPH*SINPH))
         ZQA=ZA*COSPH-XPA*SINPH
C
C    ROTATE ABOUT THE X-AXIS TO MAKE ZA=0, AND YA POSITIVE
C
         YZA=SQRT(YPA*YPA+ZQA*ZQA)
         COSKH=YPA/YZA
         SINKH=ZQA/YZA
C
C    COORDINATES:  A=(XQA,YZA,0), B=(RBC,0,0), C(0,0,0)
C    NONE ARE NEGATIVE.
C    THE COORDINATES OF I ARE EVALUATED IN THE NEW FRAME.
C
         SINA=SIN(A(2,I))
         SIND=SIN(A(3,I))
         COSD=COS(A(3,I))
         XD=A(1,I)*COSA
         YD=A(1,I)*SINA*COSD
         ZD=A(1,I)*SINA*SIND
C
C    TRANSFORM THE COORDINATES BACK TO THE ORIGINAL SYSTEM
C
         YPD=YD*COSKH-ZD*SINKH
         ZPD=ZD*COSKH+YD*SINKH
         XPD=XD*COSPH-ZPD*SINPH
         ZQD=ZPD*COSPH+XD*SINPH
         XQD=XPD*COSTH-YPD*SINTH
         YQD=YPD*COSTH+XPD*SINTH
         IF (K.LT.1) GOTO 170
         XRD=-ZQD
         ZQD=XQD
         XQD=XRD
  170    X1(I)=XQD+X1(MC)
         Y1(I)=YQD+Y1(MC)
         Z1(I)=ZQD+Z1(MC)
  180 CONTINUE
C
C    STORE THE COORDINATES ETC OF ALL ATOMS THAT CONTRIBUTE TO THE
C    PI-SYSTEM INTO THE ARRAYS TRANSFERRED TO THE MAIN PROGRAM
C
  190 J=0
      DO 200 I=1,NUMAT
CXX      IF (NST(I).LE.0.OR.NST(I).GE.12) GOTO 200
         IF(NST(I).GT.100) GOTO 200
         IF (ISPI(NST(I)+1).NE.1) GOTO 200
         J=J+1
         X(1,J)=X1(I)
         X(2,J)=Y1(I)
         X(3,J)=Z1(I)
         IAT(J)=IATD(I)
         ISORT(J)=NST(I)
         IFIX(J)=NFIX(I)
         IMOL(J)=NMOL(I)
         NNAT(I)=J
         RESTYP(J)=RESTY(I)
         RESNAM(J)=RESNA(I)
         NUMRES(J)=NURES(I)
  200 CONTINUE
      IF (J.NE.NA) WRITE(IUT,380)NA,J
C
C    NOW STORE THE COORDINATES OF ALL C(SP3)- AND H-ATOMS
C
      DO 210 I=1,NUMAT
CXX      IF (NST(I).NE.0.AND.NST(I).LT.12) GOTO 210
         IF (NST(I).GT.100) GOTO 211 
         IF (ISPI(NST(I)+1).NE.0) GOTO 210
  211    J=J+1
         X(1,J)=X1(I)
         X(2,J)=Y1(I)
         X(3,J)=Z1(I)
         IAT(J)=IATD(I)
         ISORT(J)=MOD(NST(I),100)
         IFIX(J)=NFIX(I)
         IMOL(J)=NMOL(I)
         NNAT(I)=J
         RESTYP(J)=RESTY(I)
         RESNAM(J)=RESNA(I)
         NUMRES(J)=NURES(I)
  210 CONTINUE
cp
      if (nm.lt.-4) then
      epce=0.
      epcv=0.
      naa=numat
      call crycrt
      do 212 i=1,numat
      do 212 j=1,3
 212  xo(j,i)=x(j,i)
      call outcry
      end if
      xdif=0.
      ydif=0.
      zdif=0.
      DO 2121 I=1,Numat
      xdif=xdif+x(1,i)
      ydif=ydif+x(2,i)
      zdif=zdif+x(3,i)
 2121 CONTINUE
      xdif=xdif/numat
      ydif=ydif/numat
      zdif=zdif/numat
      if (nm.gt.-6) then
      do 3212 i=1,numat
      x(1,i)=x(1,i)-xdif
      x(2,i)=x(2,i)-ydif
 3212 x(3,i)=x(3,i)-zdif
      end if
      if (nm.lt.-5) then
      xpos(1)=xdif
      xpos(2)=ydif
      xpos(3)=zdif
      end if
C
C    STORE DEFINITIONS OF INTERNAL COORDINATES
C
      DO 240 I=1,NUMAT
      ID=NNAT(I)
      JA=NAD(I)
      JB=NBD(I)
      JC=NCD(I)
      IF (JA.EQ.0) GOTO 230
      IA=NNAT(JA)
      NRA(1,1,I)=MIN0(IA,ID)
      NRA(2,1,I)=MAX0(IA,ID)
      IF(JB.EQ.0) GOTO 230
      IB=NNAT(JB)
      NRA(1,2,I)=IA
      NRA(2,2,I)=MIN0(IB,ID)
      NRA(3,2,I)=MAX0(IB,ID)
      IF (JC.EQ.0) GOTO 230
      IC=NNAT(JC)
      NRA(1,3,I)=MIN0(IA,IB)
      NRA(2,3,I)=MAX0(IA,IB)
      NRA(3,3,I)=IC
      NRA(4,3,I)=ID
      IF (IA.GT.IB) GOTO 231
      NRA(3,3,I)=ID
      NRA(4,3,I)=IC
  231 CONTINUE
  230 CONTINUE
  240 CONTINUE
      RETURN
C
C    STOP IF GEOMETRY DEFINITION IS INCORRECT
C    I.E. THE INTERNAL COORDINATES OF AN ATOM ARE DEFINED USING AN
C    ATOM THAT NOT YET OCCURED IN THE INPUT
C
  100 WRITE(IUT,320) I
      STOP
C
  250 FORMAT (I2,3I3,F9.5,2(10X,F10.5),I2,8X,3I3)
  260 FORMAT (4I3,F8.4,2F10.4,1X,4I3,4X,A5,1X,A4,I4)
  310 FORMAT (//5X,25HTRIAL GEOMETRY PARAMETERS/)
  320 FORMAT (///5X,38H***** GEOMETRY IS ILL DEFINED - ATOM #,I3)
  330 FORMAT (/4HATOM,2X,6HATOMIC,5X,11HBOND LENGTH,2X,10HBOND ANGLE
     1,6X,11HTWIST ANGLE/6HNUMBER,6HNUMBER,5X,11H(ANGSTROMS),4
     2X,9H(DEGREES),7X,9H(DEGREES)/3H(I),17X,4HNA I,10X,7HNB NA I,6
     3X,10HNC NB NA I,8X,2HNA,4X,2HNB,4X,2HNC/)
  340 FORMAT (2X,1H1,5X,I2)
  350 FORMAT (2X,1H2,5X,I2,F13.2,40X,I3)
  360 FORMAT (2X,1H3,5X,I2,2(F13.2,2X),23X,2(I3,3X))
  370 FORMAT (I3,5X,I2,2(F13.2,2X),F13.2,I2,5X,3(3X,I3))
  380 FORMAT (1X,'WARNING - EXPECTED',I5,'PI ATOMS, FOUND',I5)
C
      END
      SUBROUTINE TRANSF(I1,I2,I3,NG,NM,NMOL,NUMAT,XX,YY,X1,X2,X3,*,*)
C
C    THIS SUBROUTINE TRANSFORMS THE COORDINATES OF THE SECOND SET OF
C    ATOMS BY GIVEN TRANSLATION AND ROTATION INFORMATION (IF NM.EQ.-2)
C    OR ACCOMODATION OF THREE ATOMS NEAR GIVEN POINTS (IF NM.EQ.-3)
C
      PARAMETER(IX=301)
      DIMENSION NMOL(IX),X1(IX),X2(IX),X3(IX),XX(3,IX),YY(3,IX),A(4,4),B
     1(4,4),O(3,3),P(3,3),N(3),X(3),Y(3),D(3)
      SAVE N,NFIRST,D,O,X
      I3=I3+1
      IF (I3.GT.1) GOTO 20
      NFIRST=NUMAT
      IF (NM.EQ.-2) THEN
C
C    INPUT THE ATOMS WICH ARE TO BE POSITIONED AND RESPECTIVE POINTS
C
         DO 10 I=1,3
            D(I)=0.
   10    READ(I2,510) N(I),(O(J,I),J=1,3)
      ELSE
C
C    INPUT TRANSLATION AND ROTATION DATA
C
         IF (NG.EQ.0) THEN
            WRITE(16,600)
            READ(15,*) D,X
         ELSE
            IF (I1.EQ.0) READ(I2,500) D,X
            IF (I1.EQ.1) READ(I2,*) D,X
         END IF
      END IF
C
C    RETURN TO MAIN AND READ IN SECOND COORDINATE SET
C
      RETURN 1
   20 IF (NM.EQ.-1) THEN
C
C    PREPARE TRANSFORMATION MATRIX
C
      DO 25 I=1,4
      DO 25 J=1,4
   25 A(I,J)=0.
      DO 26 I=1,4
   26 A(I,I)=1.
         DO 30 I=1,3
            X(I)=X(I)*0.0174533
   30    CALL MATROT(A,SIN(X(I)),COS(X(I)),I,1)
      ELSE
         DO 40 I=1,3
            P(1,I)=X1(N(I))
            P(2,I)=X2(N(I))
   40    P(3,I)=X3(N(I))
         CALL SETMAT(A,P)
         CALL SETMAT(B,O)
         CALL INVMAT(B,4)
         CALL MATMLT(B,A)
         DO 60 J=1,3
         DO 50 I=1,3
         Y(I)=P(I,J)
   50    X(I)=O(I,J)
         CALL MATVEC(A,Y)
         DO 60 K=1,3
   60    D(K)=D(K)+X(K)-Y(K)
         DO 70 I=1,3
   70    D(I)=D(I)/3.
      END IF
C
C    INCLUDE LAST TRANSLATION VECTOR
C
      CALL TRXYZ(A,D)
C
C    GENERATE NEW COORDINATES BY APLICATION OF THE
C    TRANSFORMATIONMATRIX TO THE POSITION VECTORS
C
      CALL GEOUT(NFIRST,NUMAT,XX,YY,NMOL,A,X1,X2,X3)
      RETURN 2
  500 FORMAT(6F10.4)
  510 FORMAT(I3,7X,3F10.4)
  600 FORMAT ('      INPUT:  TRX,  TRY,  TRZ,  ROX,  ROY,  ROZ')
      END
      SUBROUTINE SETMAT(A,X)
C
C    TRANSFORMATION OF P1,P2,P3 TO P1(0,0,0),P2(X2,0,0),P3(X3,Y3,0)
C
      DIMENSION A(4,4),X(3,3),Y(3)
      DO 10 I=1,3
   10 Y(I)=-X(I,1)
      DO 15 I=1,4
      DO 15 J=1,4
   15 A(I,J)=0.
      DO 16 I=1,4
   16 A(I,I)=1.
      CALL TRXYZ(A,Y)
      DO 20 I=1,3
   20 Y(I)=X(I,2)
      CALL MATVEC(A,Y)
      CALL MATROT(A,Y(1),Y(2),3,0)
      DO 30 I=1,3
   30 Y(I)=X(I,2)
      CALL MATVEC(A,Y)
      CALL MATROT(A,Y(1),-Y(3),2,0)
      DO 40 I=1,3
   40 Y(I)=X(I,3)
      CALL MATVEC(A,Y)
      CALL MATROT(A,Y(2),Y(3),1,0)
      RETURN
      END
      SUBROUTINE TRXYZ(A,X)
      DIMENSION A(4,4),B(4,4),X(3),JR(3,2)
C
C    SETUP TRANSLATION AND ROTATION MATRIX
C
      DO 5 I=1,4
      DO 5 J=1,4
    5 B(I,J)=0.
      DO 6 I=1,4
    6 B(I,I)=1.
      DO 10 I=1,3
   10 B(I,4)=X(I)
      CALL MATMLT(B,A)
      RETURN
      ENTRY MATROT(A,X1,X2,IR,IT)
      DATA JR/2,3,1,3,1,2/
      IF (IT.EQ.0) THEN
         R=SQRT(X1*X1+X2*X2)
         CCOS=X1/R
         SSIN=X2/R
      ELSE
         SSIN=X1
         CCOS=X2
      END IF
      N1=JR(IR,1)
      N2=JR(IR,2)
      DO 15 I=1,4
      DO 15 J=1,4
   15 B(I,J)=0.
      DO 16 I=1,4
   16 B(I,I)=1.
      B(N1,N1)=CCOS
      B(N2,N2)=CCOS
      B(N1,N2)=SSIN
      B(N2,N1)=-SSIN
      CALL MATMLT(B,A)
      RETURN
      END
      SUBROUTINE MATMLT(B,C)
C
C    SOME MATRIX OPERATIONS
C
      DIMENSION A(4,4),B(4,4),C(4,4),X(3),Y(3)
      DO 10 I=1,4
      DO 10 J=1,4
         A(I,J)=C(I,J)
   10 C(I,J)=0.
      DO 20 I=1,4
      DO 20 J=1,4
      DO 20 K=1,4
   20 C(I,J)=B(I,K)*A(K,J)+C(I,J)
      RETURN
      ENTRY MATVEC(B,Y)
      DO 30 I=1,3
         X(I)=Y(I)
   30 Y(I)=B(I,4)
      DO 40 I=1,3
      DO 40 J=1,3
   40 Y(I)=B(I,J)*X(J)+Y(I)
      RETURN
      END
      SUBROUTINE INVMAT(A,IMAT)
C
C    MATRIX-INVERSION ROUTINE BY THE PIVOT STRATEGY
C
      DIMENSION A(IMAT,IMAT),IROW(4),JCOL(4),Y(4)
      DO 18 K=1,IMAT
         KM1=K-1
         PIVOT=0.
        DO 11 I=1,IMAT
        DO 11 J=1,IMAT
            IF(K.EQ.1) GOTO 9
            DO 8 ISCAN=1,KM1
            DO 8 JSCAN=1,KM1
               IF(I.EQ.IROW(ISCAN)) GOTO 11
               IF(J.EQ.JCOL(JSCAN)) GOTO 11
    8       CONTINUE
    9       IF(ABS(A(I,J)).LE.ABS(PIVOT)) GOTO 11
            PIVOT=A(I,J)
            IROW(K)=I
            JCOL(K)=J
   11    CONTINUE
         IF(PIVOT.EQ.0.) GOTO 33
         IROWK=IROW(K)
         JCOLK=JCOL(K)
         DO 14 J=1,IMAT
   14    A(IROWK,J)=A(IROWK,J)/PIVOT
         A(IROWK,JCOLK)=1./PIVOT
         DO 18 I=1,IMAT
            AIJCK=A(I,JCOLK)
            IF(I.EQ.IROWK) GOTO 18
            A(I,JCOLK)=-AIJCK/PIVOT
         DO 17 J=1,IMAT
   17       IF(J.NE.JCOLK) A(I,J)=A(I,J)-AIJCK*A(IROWK,J)
   18 CONTINUE
      DO 28 J=1,IMAT
      DO 27 I=1,IMAT
            IROWI=IROW(I)
            JCOLI=JCOL(I)
   27    Y(JCOLI)=A(IROWI,J)
      DO 28 I=1,IMAT
   28 A(I,J)=Y(I)
      DO 30 I=1,IMAT
      DO 29 J=1,IMAT
            IROWJ=IROW(J)
            JCOLJ=JCOL(J)
   29    Y(IROWJ)=A(I,JCOLJ)
      DO 30 J=1,IMAT
   30 A(I,J)=Y(J)
      RETURN
C
C    STOP IF THE MATRIX TO BE INVERTED IS SINGULAR
C
   33 WRITE(6,*) ' ***** TRANSFORMATION-MATRIX IS SINGULAR'
      STOP
      END
      SUBROUTINE GEOUT(NA,NE,QA,Q,NMOL,A,X,Y,Z)
C
C    COORDINATE TRANSFORMATION    P'  <----  P
C
      PARAMETER(IX=301)
      DIMENSION NMOL(IX),Q(3,IX),QA(3,IX),A(4,4),X(IX),Y(IX),Z(IX)
      ITRANS=NE-NA
      DO 30 K=1,ITRANS
      K1=K+NA
      QA(1,K)=X(K1)
      QA(2,K)=Y(K1)
      QA(3,K)=Z(K1)
      DO 10 M=1,3
   10 Q(M,K)=A(M,4)
      DO 20 I=1,3
      DO 20 J=1,3
   20 Q(I,K)=A(I,J)*QA(J,K)+Q(I,K)
      IF (NMOL(K1).EQ.0) NMOL(K1)=1
      X(K1)=Q(1,K)
      Y(K1)=Q(2,K)
   30 Z(K1)=Q(3,K)
      RETURN
      END
      SUBROUTINE SETGEO
C
C    THIS SUBROUTINE SEARCHES AND STORES DEFINITIONS FOR ALL GEOMETRY
C    PARAMETERS (BOND LENGTHS, BOND ANGLES, TORSION ANGLES, BEND ANGLES)
C    BASED ON INPUT TOPOLOGY
C
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      PARAMETER(IX=301,IY=IX,IX3=IX*3,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      common /dyna/ idz,idd(4,40),irdz(ix6),idyn,dystar
      NAM=NAA-1
C
C    SET EQUILIBRIUM BOND LENGTHS FOR SOME BONDS
C
      DO 10 I=1,NAA
      ISO1=ISORT(I)+1
      DO 10 J=1,NAA
      IF (J.LE.NA) GOTO 10
      ISO2=ISORT(J)+1
      IF (VS(I,J).LE.0.) GOTO 10
      VS(I,J)=AL(ISO1,ISO2)
      IF (QADD(I)*QADD(J).LT.-0.05) VS(I,J)=VS(I,J)*0.92
      VS(J,I)=VS(I,J)
   10 CONTINUE
C
C    FIND ALL BONDS
C
      IBOND=0
      DO 20 I=1,NAA
      DO 20 J=I,NAA
      IF (VS(I,J).LE.0.) GOTO 20
      IBOND=IBOND+1
      NBOND(1,IBOND)=I
      NBOND(2,IBOND)=J
      RBON(IBOND)=0.
   20 CONTINUE
C
C    FIND ALL BOND ANGLES
C
      IANG=0
      DO 40 I=1,NAA
      IF (ISORT(I).EQ.12.OR.ISORT(I).EQ.13) GOTO 40
      DO 39 J=1,NAM
      IF (ISORT(J).EQ.12.OR.ISORT(J).EQ.13) GOTO 39
      IF (VS(I,J).LE.0.) GOTO 39
      J1=J+1
      DO 30 K=J1,NAA
      IF (ISORT(K).EQ.12.OR.ISORT(K).EQ.13) GOTO 30
      IF (VS(I,K).LE.0.) GOTO 30
      IANG=IANG+1
      NANG(1,IANG)=I
      NANG(2,IANG)=J
      NANG(3,IANG)=K
      RANG(IANG)=0.
      if (vs(j,k).eq.0..and.vs(k,j).eq.0.) then
      vs(j,k)=-2.
      vs(k,j)=-2.
      end if
      DO 26 M=1,NAA
      IF (M.EQ.I.OR.VS(M,J).LE.0.) GOTO 26
      DO 24 N=1,NAA
      IF (N.EQ.I.OR.VS(N,K).LE.0.) GOTO 24
      IF (VS(J,K).GT.0..OR.M.EQ.N.OR.VS(M,N).GT.0.)
     1   NANG(1,IANG)=-IABS(NANG(1,IANG))
      IF (VS(N,J).GT.0..OR.VS(M,K).GT.0.)
     1   NANG(1,IANG)=-IABS(NANG(1,IANG))
   24 CONTINUE
   26 CONTINUE
   30 CONTINUE
   39 CONTINUE
   40 CONTINUE
C
C    FIND ALL TORSION ANGLES AND BEND ANGLES
C
      ITORS=0
      IBEND=0
      DO 70 I=1,NAA
      ISP1=ISORT(I)
      IF (ISP1.EQ.12.OR.ISP1.EQ.13) GOTO 70
      IBFLA=0
      DO 69 J=1,NAA
      IF (VS(I,J).LE.0.) GOTO 69
      ISP2=ISORT(J)
      IF (ISP2.EQ.12.OR.ISP2.EQ.13) GOTO 69
      DO 60 K=1,NAA
      IF (K.EQ.I.OR.K.EQ.J) GOTO 60
      IF (VS(I,K).LE.0.) GOTO 60
      IF (ISORT(K).EQ.12.OR.ISORT(K).EQ.13) GOTO 60
      IF (IBFLA.EQ.1.OR.ISP1.LT.1) GOTO 85
      DO 90 M=1,NAA
      IF (VS(I,M).LE.0.) GOTO 90
      IF (M.EQ.J.OR.M.EQ.K) GOTO 90
      IF (ISORT(M).EQ.12.OR.ISORT(M).EQ.13) GOTO 90
      IBFLA=1
      IBEND=IBEND+1
      NBEND(1,IBEND)=I
      NBEND(2,IBEND)=J
      NBEND(3,IBEND)=K
      NBEND(4,IBEND)=M
      RBEND(IBEND)=0.
      GOTO 85
  90  CONTINUE
  85  DO 50 L=1,NAA
      IF (L.EQ.I.OR.L.EQ.J.OR.L.EQ.K) GOTO 50
      IF (VS(J,L).LE.0.) GOTO 50
      IF (ISORT(L).EQ.12.OR.ISORT(L).EQ.13) GOTO 50
      IF (J.LT.I) GOTO 50
      ITORS=ITORS+1
      NTORS(1,ITORS)=I
      NTORS(2,ITORS)=J
      NTORS(3,ITORS)=K
      NTORS(4,ITORS)=L
      RTORS(ITORS)=0.
      if (vs(l,k).eq.0..and.vs(k,l).eq.0.) then
      vs(k,l)=-1.
      vs(l,k)=-1.
      end if
      ctors(itors)=1.
cm
      irdz(itors)=0.
      if (nh.lt.0) then
      do 8100 idy=1,idz
      if (i.ne.idd(1,idy)) goto 8100
      if (j.ne.idd(2,idy)) goto 8100
      if (k.ne.idd(3,idy)) goto 8100
      if (l.ne.idd(4,idy)) goto 8100
      irdz(itors)=idy
 8100 continue
      endif
cm
  50  CONTINUE
  60  CONTINUE
  69  CONTINUE
  70  CONTINUE
      do 71 i=1,itors
      do 71 j=i+1,itors
      if (ntors(3,I).eq.ntors(3,j).or.ntors(3,i).eq.ntors(4,j)) then
      if (ntors(4,I).eq.ntors(3,j).or.ntors(4,i).eq.ntors(4,j)) then
      ctors(i)=0.5
      ctors(j)=0.5
      do 72 k=1,itors
      if (ntors(1,k).eq.ntors(1,i).and.ntors(2,k).eq.ntors(2,i)) then
      ctors(k)=0.5 
      goto 72
      end if
      if (ntors(1,k).eq.ntors(1,j).and.ntors(2,k).eq.ntors(2,j)) then
      ctors(k)=0.5 
      goto 72
      end if
      if (ntors(1,k).eq.ntors(2,i).and.ntors(2,k).eq.ntors(1,i)) then
      ctors(k)=0.5 
      goto 72
      end if
      if (ntors(1,k).eq.ntors(2,j).and.ntors(2,k).eq.ntors(1,j)) then
      ctors(k)=0.5 
      goto 72
      end if
 72   continue
      end if
      end if
 71   continue

C
C    FIND ALL TORSION ANGLES OF TYPE DELTA AND TAU
C
      IDELT=0
      DO 140 I=1,NAM
      IF (ISORT(I).EQ.12.OR.ISORT(I).EQ.13) GOTO 140
      II=I+1
      DO 139 J=II,NAA
      IF (VS(I,J).LE.0.) GOTO 139
      IF (ISORT(J).EQ.12.OR.ISORT(J).EQ.13) GOTO 139
      DO 130 K=1,NAM
      IF (K.EQ.I.OR.K.EQ.J) GOTO 130
      IF (VS(I,K).LE.0..AND.VS(J,K).LE.0.) GOTO 130
      IF (ISORT(K).EQ.12.OR.ISORT(K).EQ.13) GOTO 130
      KK=K+1
      DO 120 L=KK,NAA
      IF (L.EQ.I.OR.L.EQ.J.OR.L.EQ.K) GOTO 120
      IF (ISORT(L).EQ.12.OR.ISORT(L).EQ.13) GOTO 120
      IF (VS(I,K).GT.0..AND.VS(I,L).GT.0.) GOTO 110
      IF (VS(J,K).GT.0..AND.VS(J,L).GT.0.) GOTO 110
      GOTO 120
  110 IDELT=IDELT+1
      NDELT(1,IDELT)=I
      NDELT(2,IDELT)=J
      NDELT(3,IDELT)=K
      NDELT(4,IDELT)=L
      RDELT(IDELT)=0.
  120 CONTINUE
  130 CONTINUE
  139 CONTINUE
  140 CONTINUE
C
c
c
      do 5010 i=1,naa
      if (isort(i).lt.2.or.isort(i).gt.9) goto 5010
      do 5009 j=1,naa
      if (isort(j).eq.19.and.vs(i,j).le.0.) then
      do 5008 k=1,naa
      if (k.eq.i) goto 5008
      if (isort(k).lt.2.or.isort(k).gt.9) goto 5008
      if (vs(j,k).gt.0..and.vs(i,k).eq.0.) then
      vs(i,j)=-3.
      vs(j,i)=-3.
      goto 5009
      endif
 5008 continue
      endif
 5009 continue
 5010 continue     
      do 5011 i=1,naa
      do 5011 j=i,naa
      if (vs(i,j).gt.0.) goto 5011
      vs(i,j)=0.
      vs(j,i)=0.
 5011 continue

      RETURN
      END
      SUBROUTINE GEO(ifin,convmm)
C
C   THIS SUBROUTINE DOES THE FORCE FIELD CALCULATION.
C   THE ENERGY OF DEFORMATION AND THE GRADIENTS ARE
C   EVALUATED FOR DEFORMATIONS OF BOND LENGTHS, BOND
C   ANGLES, TORSIONAL ANGLES, BEND ANGLES OF (SP2)-CENTRES
C   AND NONBONDING INTERACTIONS.
C
      DOUBLE PRECISION XK,YK,ZK,xkc,ykc,zkc,gmax,virx,viry,virz,
     1vircx,vircy,vircz,viryx,virzx,virzy
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      PARAMETER(IX=301,IY=IX,IX3=IX*3,IX4=IX*4,IX6=IX*6,
     1IXA=(IX3*(IX3+1))/2)
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /OPTIM/ virx,viry,virz,XD(IX3),XDD(IX3),QZ(IX3),DX(IX3),
     1XK(IX),YK(IX),ZK(IX),NAM,EPP,EPW
      COMMON /OPTIM1/ AMAT(IXA)
cp
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1 cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2 imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3 ,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
      common /pacplo/ npac,ispac(2000),cpac(3,2000),ihhyd(100),
     1 iahyd(100),dhyd(100),nhyd,jhyd(ix),xo(3,ix),xsum(3,ix),
     2 ibhyd(100),hhyd(100),abhyd(100),anghyd(100)   
cp
C
C    SET SOME FLAGS AND CLEAR ARRAYS
C
      naax=naa
      if (nm.lt.-6) naax=naa+3
      NAM=NAAx-1
      NXX=NN-10
      SLOPE=0.
      if (nn.eq.ne) then
      do 9 i=1,naa
    9 jhyd(i)=0
      end if
      DO 10 I=1,NAAx
         XK(I)=0.D0
         YK(I)=0.D0
   10 ZK(I)=0.D0
      IBLO=0
      SGA=0.
      EPP=0.
      EPW=0.
      EPA=0.
      virx=0.
      viry=0.
      virz=0.
      KA=3*NAAx
C
C    SET FIRST TRIAL OF HESSIAN MATRIX
C
      K=0
      DO 20 I=1,KA
         K=K+1
         DX(I)=0.
         XDD(I)=0.
         AMAT(K)=60.
      DO 20 J=I+1,KA
         K=K+1
         AMAT(K)=0.
   20 CONTINUE
   40 DO 30 I=1,NAA
      DO 30 J=1,NAA
         IF (VS(I,J).GT.0.) GOTO 30
         VS(I,J)=0.
   30 CONTINUE
C
C    BOND LENGTH SECTION
C
      CALL BOND
C
C    BOND ANGLE SECTION
C
      CALL ANGLE
C
C    TORSIONAL ANGLE SECTION
C
      CALL TORS
      CALL BEND
C
C    NONBONDING INTERACTION SECTION
cp
      if (nm.lt.-5) then
      call crypac
      do 99 i=1,naa
      xk(i)=xk(i)+xkc(i)
      yk(i)=yk(i)+ykc(i)
      zk(i)=zk(i)+zkc(i)
   99 continue
      if (nm.lt.-6) then
      if (iscl.eq.1) then
      xkcm=(xkc(naa+1)+ykc(naa+2)+zkc(naa+3))/3.
      xkc(naa+1)=xkcm
      ykc(naa+2)=xkcm
      zkc(naa+3)=xkcm
      end if
      if (iscl.eq.2) then
      xkcm=(xkc(naa+1)+ykc(naa+2))/2.
      xkc(naa+1)=xkcm
      ykc(naa+2)=xkcm
      end if
      xk(naa+1)=xkc(naa+1)*ascl(1,1)
      yk(naa+1)=ykc(naa+1)*ascl(1,2)
      zk(naa+1)=zkc(naa+1)*ascl(1,3)
      xk(naa+2)=xkc(naa+2)*ascl(2,1)
      yk(naa+2)=ykc(naa+2)*ascl(2,2)
      zk(naa+2)=zkc(naa+2)*ascl(2,3)
      xk(naa+3)=xkc(naa+3)*ascl(3,1)
      yk(naa+3)=ykc(naa+3)*ascl(3,2)
      zk(naa+3)=zkc(naa+3)*ascl(3,3)
      end if
      end if
C
      CALL NONBON
C
C    SET GRADIENTS TO ZERO, IF ATOM IS FIXED TO ITS POSITION
C
      DO 100 I=1,NAA
      IF (IFIX(I).EQ.0) GOTO 100
         ZK(I)=0.D0
      IF (IFIX(I).EQ.2) GOTO 100
         XK(I)=0.D0
         YK(I)=0.D0
  100 CONTINUE
C
C    OPTIMIZATION
C    MURTAGH-SARGENT PROCEDURE ADOPTED FROM A PROCEDURE
C    WRITTEN BY J. PANCIR
C
      SG=0.
      gmax=0.
      ifin=0
      DO 410 I=1,NAAx
         IF (IBLO.EQ.0) GOTO 800
         XK(I)=XK(I)*ZBB
         YK(I)=YK(I)*ZBB
         ZK(I)=ZK(I)*ZBB
         IF (DABS(XK(I)).LT.0.1D-9) XK(I)=0.D0
         IF (DABS(YK(I)).LT.0.1D-9) YK(I)=0.D0
         IF (DABS(ZK(I)).LT.0.1D-9) ZK(I)=0.D0
  800    gmax=dmax1(dabs(xk(i)),dabs(yk(i)),dabs(zk(i)),gmax)
         II=NAAx+I
         III=2*NAAx+I
         XD(I)=-XK(I)
         SG=SG+XK(I)**2
         XD(II)=-YK(I)
         SG=SG+YK(I)**2
         XD(III)=-ZK(I)
  410 SG=SG+ZK(I)**2
      AKA=KA
      EP=EPW
      EPP=EPP+EPW
      IF (NG.LT.2) return
      SGT=1.E-18*AKA*AKA
      IF (SG.LE.SGT) GOTO 581
      IF (NN.GE.NXX) GOTO 510
      IF (IBLO.EQ.1) GOTO 510
      IF (SGA.GT.SG) GOTO 440
      A1=0.
      A2=0.
      DO 420 I=1,KA
         A1=XD(I)*DX(I)+A1
  420 A2=(XDD(I)-XD(I))*DX(I)+A2
      A1=A1/A2
      IF (A1.LT.(-0.5)) A1=A1*0.5
      DO 430 I=1,KA
  430 DX(I)=DX(I)*A1
      GOTO 540
  440 DO 450 I=1,KA
         QZ(I)=DX(I)
         K=I
      DO 451 J=1,I
      QZ(I)=QZ(I)-AMAT(K)*(XD(J)-XDD(J))
  451 K=K+KA-J
      K=K-KA+I
      DO 450 J=I+1,KA
      K=K+1
  450 QZ(I)=QZ(I)-AMAT(K)*(XD(J)-XDD(J))
      CN=0.
      C1=0.
      C2=0.
      DO 460 I=1,KA
         CN=CN+(XD(I)-XDD(I))*QZ(I)
         C2=C2-QZ(I)*XD(I)
  460 C1=C1+QZ(I)*QZ(I)
      A1=ABS(CN)/C1
      IF (A1.LT.0.1E-4) GOTO 470
      A1=C2/CN
      IF (A1.LT.0.1E-5) GOTO 470
      GOTO 490
  470 CN=C1
      K=0
      DO 480 I=1,KA
         K=K+1
         AMAT(K)=60.
      DO 480 J=I+1,KA
         K=K+1
         AMAT(K)=0.
  480 CONTINUE
      IF (SLOPE.GT.0.5E-3) SLOPE=SLOPE*0.5
      GOTO 520
  490 K=0
      DO 500 I=1,KA
      DO 500 J=I,KA
      K=K+1
  500 AMAT(K)=AMAT(K)+QZ(I)*QZ(J)/CN
      GOTO 520
  510 A1=0.0
      SLOPE=1.
  520 DO 530 I=1,KA
         K=I
         DX(I)=0.
      DO 531 J=1,I
         DX(I)=DX(I)-SLOPE*AMAT(K)*XD(J)
  531    K=K+KA-J
      K=K-KA+I
      DO 530 J=I+1,KA
         K=K+1
  530 DX(I)=DX(I)-SLOPE*AMAT(K)*XD(J)
  540 A1=0.
      DO 550 I=1,KA
  550 A1=A1+DX(I)*DX(I)
      IF (A1.LE.0.02) GOTO 570
      FAD=0.02/A1
      DO 560 I=1,KA
  560 DX(I)=DX(I)*FAD
  570 DO 580 I=1,KA
  580 XDD(I)=XD(I)
      GOTO 583
  581 DO 582 I=1,KA
  582 DX(I)=0.
  583 nnn=nxx-nn+10 
      SGA=SG
      A1=SQRT(A1)/AKA
      SG=SQRT(SG)/AKA*6242.
C
C   CHECK FOR CONVERGENCE OF OPTIMIZATION
C   TEST IF OPTIMIZATION MUST BE RESTARTED
C   SET FLAGS FOR OUTPUT IN THE FINAL CYCLE OF OPTIMIZATION
C
  595 IF (gmax.lt.0.5d-5.and.SG.LT.convmm.AND.NN.LT.NE) GOTO 600
      DO 590 I=1,NAAx
         II=NAAx+I
         III=2*NAAx+I
         X(1,I)=X(1,I)+DX(I)
         X(2,I)=X(2,I)+DX(II)
         X(3,I)=X(3,I)+DX(III)
         xo(1,i)=x(1,i)
         xo(2,i)=x(2,i)
         xo(3,i)=x(3,i)
         XK(I)=0.D0
         YK(I)=0.D0
  590 ZK(I)=0.D0
         virx=0.
         viry=0.
         virz=0.
      if (nm.lt.-5) call crtcry
      IF (A1.GE.0.1E-7) GOTO 871
      SLOPE=1.
      IBLO=1
      K=0
      DO 872 I=1,KA
         K=K+1
         AMAT(K)=60.
      DO 872 J=I+1,KA
         K=K+1
  872 AMAT(K)=0.
      GOTO 873
  871 IBLO=0
  873 EPA=EPP
      EPP=0.
      EPW=0.
      NN=NN-1
      IF (NNN.GE.100) GOTO 600
      IF (NN.GE.0) GOTO 40
C
C   DO NOT END OPTIMIZATION IN AN INSTABLE STATE.
C
      IF (SG.GE.0.1E-1) GOTO 40
C
C   OUTPUT OF ENERGY OF DEFORMATION AND SOME CHARACTERISTIC VALUES
C   TO CHECK OPTIMIZATION
C
  600 IF (ISEQ.EQ.2.AND.NC.LE.0) goto 6000
      WRITE(IUT,660) NNN,A1,SG,EPA
 6000 if (nnn.le.1) ifin=1
      RETURN
C
  660 FORMAT (1H ,34X,I5,2X,2E9.2,F9.1)
C
      END
      SUBROUTINE BOND
C
C    THIS SUBROUTINE CALCULATES ENERGIES AND FORCES FOR BOND LENGTH
C    DISTORTIONS
C
      DOUBLE PRECISION XK,YK,ZK,virx,viry,virz
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      PARAMETER(IX=301,IY=IX,IX3=IX*3,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /OPTIM/virx,viry,virz,XD(IX3),XDD(IX3),QZ(IX3),DX(IX3),
     1XK(IX),YK(IX),ZK(IX),NAM,EPP,EPW
      COMMON /STRAIN/ SPANB(IX),SPANA(IX),SPANT(IX),SPANK(IX),SPANE(IX),
     1SPANV(IX)
C*      SPANB=0.
      NAM=NAA-1
C***************
      DO 1 I=1,NAA
      SPANB(I)=0.
    1 CONTINUE
C
C
C    CALCULATE ALL INTERATOMIC DISTANCES AND FIND THE BONDS
C
      DO 10 I=1,NAM
         ISO1=ISORT(I)+1
         I1=I+1
      DO 10 J=I1,NAA
        VV(I,J)=SQRT((X(1,I)-X(1,J))**2+(X(2,I)-X(2,J))**2+(X(3,I)-X(3,J
     1  ))**2)
   10   VV(J,I)=VV(I,J)
      DO 100 N=1,IBOND
         I=NBOND(1,N)
         J=NBOND(2,N)
         ISO1=ISORT(I)+1
         ISO2=ISORT(J)+1
C
C    CALCULATE FORCE CONSTANT OF THE BOND I-J
C
         PKON=(PK2(ISO1,ISO2)/VS(I,J)**2+PK4(ISO1,ISO2)/VS(I,J)**4+
     1   PK6(ISO1,ISO2)/VS(I,J)**6)*0.001
C
C    EVALUATE THE STRAIN ENERGY OF THE BOND I-J AND THE COMPONENTS
C    OF THE FORCES ACTING ON THE ATOMS I AND J
C
         DELV=VV(I,J)-VS(I,J)
         SPAN=ABS(DELV)**2*PKON*3121.
         IF (IMOL(I).EQ.99) SPAN=0.
         EPP=EPP+SPAN
         SPANB(I)=SPANB(I)+SPAN/2.
         SPANB(J)=SPANB(J)+SPAN/2.
         DELF=PKON*(1.-VS(I,J)/VV(I,J))
         IF (NH.EQ.0) GOTO 50
         DO 40 JRST=1,IRZ
            IF (NBOND(1,N).NE.IRST(1,JRST)) GOTO 40
            IF (NBOND(2,N).NE.IRST(2,JRST)) GOTO 40
            IF (IRST(3,JRST).NE.0) GOTO 40
            DELF=0.01*(1.-ARST(JRST)/VV(I,J))
            GOTO 50
   40    CONTINUE
   50    XK(I)=XK(I)-DELF*(X(1,I)-X(1,J))
         YK(I)=YK(I)-DELF*(X(2,I)-X(2,J))
         ZK(I)=ZK(I)-DELF*(X(3,I)-X(3,J))
         XK(J)=XK(J)+DELF*(X(1,I)-X(1,J))
         YK(J)=YK(J)+DELF*(X(2,I)-X(2,J))
         ZK(J)=ZK(J)+DELF*(X(3,I)-X(3,J))
         RBON(N)=VV(I,J)
  100 CONTINUE
      RETURN
      END
      SUBROUTINE ANGLE
C
C    THIS SUBROUTINE CALCULATES ENERGIES AND FORCES FOR BOND ANGLE
C    DISTORTIONS
C
      DOUBLE PRECISION XK,YK,ZK,virx,viry,virz
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      PARAMETER(IX=301,IY=IX,IX3=IX*3,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /OPTIM/virx,viry,virz,XD(IX3),XDD(IX3),QZ(IX3),DX(IX3),
     1XK(IX),YK(IX),ZK(IX),NAM,EPP,EPW
      COMMON /STRAIN/ SPANB(IX),SPANA(IX),SPANT(IX),SPANK(IX),SPANE(IX),
     1SPANV(IX)
C*      SPANA=0.
C
      DO 5 I=1,NAA
C**************
      SPANA(I)=0.
C**************
    5 SWI(I)=0.
C
C    FIND BOND ANGLE J -- I -- K
C
      DO 100 N=1,IANG
         I=IABS(NANG(1,N))
         J=NANG(2,N)
         K=NANG(3,N)
         IANGA=NANG(1,N)/I
         ISO1=ISORT(I)+1
C
C    CALCULATE BOND ANGLE
C
            VR=(X(1,I)-X(1,J))*(X(1,I)-X(1,K))+(X(2,I)-X(2,J))*(X(2,I)-X
     1      (2,K))+(X(3,I)-X(3,J))*(X(3,I)-X(3,K))
            COSG=VR/(VV(I,J)*VV(I,K))
            ANG=ACOSF(COSG)
c
            angina=amax1(ang,1.57)
C
C    FIND DELTA(ALFA) AND THE FORCE CONSTANT
C
            SANG=GANG(ISO1)
            DST=DSW(ISO1)
            if (isort(j).eq.19.and.isort(k).eq.19) dst=dst*1.2
            if (isort(i).eq.4.and.isort(j).eq.19) dst=dsw(iso1)*0.5
            if (isort(i).eq.4.and.isort(k).eq.19) dst=dsw(iso1)*0.5
            IF (ISO1.EQ.1.AND.ISORT(J).NE.19.AND.ISORT(K).NE.19.AND.
     1      IANGA.GT.0) SANG=1.975
            IF (ISO1.EQ.5.AND.ISORT(J).EQ.14.AND.ISORT(K).EQ.14)
     1      SANG=2.53
            IF (I.GT.NA) GOTO 10
            IF (P(I,I).LE.0.) GOTO 10
            IF (ISO1.NE.3) GOTO 10
            SANG=2.0944-1.7453*(P(I,I)-1.90)
            IF (SANG.GT.2.0944) SANG=2.0944
  10        DELTA=(SANG-ANGina)
            IF (NH.EQ.0) GOTO 30
            DO 20 JRST=1,IRZ
               IF (I.NE.IRST(1,JRST)) GOTO 20
               IF (J.NE.IRST(2,JRST)) GOTO 20
               IF (K.NE.IRST(3,JRST)) GOTO 20
               IF (IRST(4,JRST).NE.0) GOTO 20
               DST=0.5
               DELTA=ARST(JRST)-ANG
   20       CONTINUE
   30       SWI(I)=SWI(I)+ABS(DELTA)
C
C    CALCULATE STRAIN ENERGY DUE TO BOND ANGLE DEFORMATION
C    AND THE COMPONENTS OF THE RESULTING FORCES ON THE ATOMS
C    I, J, K
C
            APK=DELTA*DST+DELTA**3*DECUB(ISO1)*0.5
            SPAN=DST*DELTA**2*3121.+DECUB(ISO1)*DELTA**4*1561.
            IF (IMOL(I).EQ.99) SPAN=0.
            EPW=EPW+SPAN
            SPANA(I)=SPANA(I)+SPAN
            LJ=J
            LK=K
   40      FN=SQRT((VR*(X(1,I)-X(1,LJ))-VV(I,LJ)**2*(X(1,I)-X(1,LK)))**2
     1     +(VR*(X(2,I)-X(2,LJ))-VV(I,LJ)**2*(X(2,I)-X(2,LK)))**2+(VR*(X
     2     (3,I)-X(3,LJ))-VV(I,LJ)**2*(X(3,I)-X(3,LK)))**2)*VV(I,LJ)
           IF(FN.EQ.0.)FN=1.E-2
           PAX=APK*(VR*(X(1,I)-X(1,LJ))-VV(I,LJ)**2*(X(1,I)-X(1,LK)))/FN
           PAY=APK*(VR*(X(2,I)-X(2,LJ))-VV(I,LJ)**2*(X(2,I)-X(2,LK)))/FN
           PAZ=APK*(VR*(X(3,I)-X(3,LJ))-VV(I,LJ)**2*(X(3,I)-X(3,LK)))/FN
           XK(I)=XK(I)+PAX
           YK(I)=YK(I)+PAY
           ZK(I)=ZK(I)+PAZ
           XK(LJ)=XK(LJ)-PAX
           YK(LJ)=YK(LJ)-PAY
           ZK(LJ)=ZK(LJ)-PAZ
           IF (LJ.EQ.K) GOTO 50
           LJ=K
           LK=J
           GOTO 40
C
C    SET FLAG TO EXCLUDE NONBONDING INTERACTION J -- K
C
   50      IF (VS(K,J).LE.0.) VS(K,J)=-2.
           IF (VS(J,K).LE.0.) VS(J,K)=-2.
           RANG(N)=ANG
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE TORS
C
C    THIS SUBROUTINE CALCULATES ENERGIES AND FORCES FOR TORSIONAL
C    ANGLES
C
      DOUBLE PRECISION XK,YK,ZK,virx,viry,virz
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      PARAMETER(IX=301,IY=IX,IX3=IX*3,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /OPTIM/virx,viry,virz,XD(IX3),XDD(IX3),QZ(IX3),DX(IX3),
     1XK(IX),YK(IX),ZK(IX),NAM,EPP,EPW
      COMMON /STRAIN/ SPANB(IX),SPANA(IX),SPANT(IX),SPANK(IX),SPANE(IX),
     1SPANV(IX)
cs
      common /solvat/ eps,dielec,dsgl,fav(ix),fam(ix),rxy(ix,ix),
     1 rvdw(50),ion(50)
      common /dyna/ idz,idd(4,40),irdz(ix6),idyn,dystar
cs
C*      SPANT=0.
      IF (ITORS.EQ.0) RETURN
C***********
      DO 1 I=1,NAA
      SPANT(I)=0.
    1 CONTINUE
C
      DO 100 N=1,ITORS
C
C
C    FIND THE ATOMS THAT DEFINE A TORSINAL ANGLE K -- I -- J -- L
C
         I=NTORS(1,N)
         J=NTORS(2,N)
         K=NTORS(3,N)
         L=NTORS(4,N)
         ISO1=ISORT(I)+1
         ISP1=ISORT(I)
         ISO2=ISORT(J)+1
         ISP2=ISORT(J)
            ISO3=ISORT(K)+1
            ISO4=ISORT(L)+1
C
C    CALCULATE NORMAL VECTOR OF THE PLANE K - I - J
C
               TAX=(X(2,I)-X(2,J))*(X(3,I)-X(3,K))-(X(3,I)-X(3,J))*(X(2,
     1         I)-X(2,K))
               TAY=(X(3,I)-X(3,J))*(X(1,I)-X(1,K))-(X(1,I)-X(1,J))*(X(3,
     1         I)-X(3,K))
               TAZ=(X(1,I)-X(1,J))*(X(2,I)-X(2,K))-(X(2,I)-X(2,J))*(X(1,
     1         I)-X(1,K))
               TNA=SQRT(TAX*TAX+TAY*TAY+TAZ*TAZ)
               IF(TNA.EQ.0.)TNA=1.E-2
               TAX=TAX/TNA
               TAY=TAY/TNA
               TAZ=TAZ/TNA
C
C    CALCULATE NORMAL VECTOR OF THE PLANE I - J - L
C
               TBX=(X(2,I)-X(2,J))*(X(3,J)-X(3,L))-(X(3,I)-X(3,J))*(X(2,
     1         J)-X(2,L))
               TBY=(X(3,I)-X(3,J))*(X(1,J)-X(1,L))-(X(1,I)-X(1,J))*(X(3,
     1         J)-X(3,L))
               TBZ=(X(1,I)-X(1,J))*(X(2,J)-X(2,L))-(X(2,I)-X(2,J))*(X(1,
     1         J)-X(1,L))
               TNB=SQRT(TBX*TBX+TBY*TBY+TBZ*TBZ)
               IF(TNB.EQ.0.)TNB=1.E-2
               TBX=TBX/TNB
               TBY=TBY/TNB
               TBZ=TBZ/TNB
C
C    FIND THE VECTOR BISECTING THE TORSIONAL ANGLE TO DEFINE
C    SIGN OF THE TORSIONAL ANGLE
C
               ZAX=TAY*TBZ-TAZ*TBY
               ZAY=TAZ*TBX-TAX*TBZ
               ZAZ=TAX*TBY-TAY*TBX
               SZ=ZAX*(X(1,I)-X(1,J))+ZAY*(X(2,I)-X(2,J))+ZAZ*(X(3,I)-X(
     1         3,J))
               SIGTAU=SIGN(1.,SZ)
C
C    EVALUATE TORSIONAL ANGLE
C
               COSTAU=TAX*TBX+TAY*TBY+TAZ*TBZ
               IF (COSTAU.LT.-1.) COSTAU=-1.
               IF (COSTAU.GT.1.) COSTAU=1.
               TAU=-ACOSF(COSTAU)*SIGTAU
               VOR=0.
C
C    IF I--J IS A PI-BOND, DERIVE THE CONSTANT FOR THE TORSIONAL
C    POTENTIAL FROM THE BOND LENGTH I - J
C
               IF (I.GT.NA.OR.J.GT.NA) GOTO 10
               PIBO=(AL(ISO1,ISO2)-VS(I,J))/BL(ISP1,ISP2)
               VOR=(PIBO+0.15*PIBO**2
     1         +PTT(ISP1,ISP2)*SIN(TAU)**2)*0.75E-4
  10           TTAU=-TAU
               PTOR=0.
               IF (ISO1.EQ.11.AND.ISO2.EQ.11) GOTO 60
               TFIX=0.
               IF (NH.EQ.0) GOTO 30
               DO 20 JRST=1,IRZ
                  IF (I.NE.IRST(1,JRST)) GOTO 20
                  IF (J.NE.IRST(2,JRST)) GOTO 20
                  IF (K.NE.IRST(3,JRST)) GOTO 20
                  IF (L.NE.IRST(4,JRST)) GOTO 20
                  TFIX=1.
                  FTAU=-ARST(JRST)
                  IF (FTAU.GT.3.1416) FTAU=FTAU-6.283
                 GOTO 30
  20           CONTINUE
  30           IF (VOR.LE.0.) VOR=0.
C
C    CALCULATE TORSIONAL STRAIN ENERGY FOR TORSIONS OF PI-BONDS
C    AND THE FIRST DERIVATIVE
C
               IF (I.GT.NA.OR.J.GT.NA) GOTO 40
               PPITA=PPIT(ISO3,ISO4)
               if (iso1.eq.4.or.iso2.eq.4) vor=vor*2.
               if (iso1.eq.5.or.iso2.eq.5) vor=vor*2.
               if (iso1.eq.7.or.iso2.eq.7) vor=vor*2.
               PTOR=SIN(TAU)*COS(TAU)*VOR*2.+SIN(TAU)**3*COS(TAU)
     1         *PTT(ISP1,ISP2)*1.50E-4
               EEP=VOR*SIN(TAU)**2*6242.
               IF (IMOL(I).EQ.99) EEP=0.
               EPP=EPP+EEP
               PTOR=PTOR-0.5*PPITA*SIN(2.*(-TTAU))
               EEPW=PPITA*(1.-COS(2.*(-TTAU)))*1561.
               IF (IMOL(I).EQ.99) EEPW=0.
               EPW=EPW-EEPW
               SPAN=VOR*SIN(TAU)**2*6242.
     1         -PPITA*(1.-COS(2.*(-TTAU)))*1561.
               IF (IMOL(I).EQ.99) SPAN=0.
               SPANT(I)=SPANT(I)+SPAN/2.
               SPANT(J)=SPANT(J)+SPAN/2.
               GOTO 50
C
C    CALCULATE TORSIONAL STRAIN ENERGY FOR TORSIONS OF SIGMA-BONDS
C    AND THE FIRST DERIVATIVE
C
  40           if (i.le.na.or.j.le.na) goto 140
               IF (VS(K,L).GT.0.OR.VS(K,L).LT.-1.) GOTO 50
               QTOR=0.25
               IF (ISO3.EQ.16.OR.ISO4.EQ.16) GOTO 45
               IF (ISO1.GT.16.OR.ISO2.GT.16) THEN
               QTOR=0.05
               GOTO 90
               END IF
               IF (ISO1.NE.5.AND.ISO2.NE.5) GOTO 45
               IF (ISO3.EQ.2.OR.ISO4.EQ.2) THEN
               QTOR=2.
               GOTO 90
               END IF
  45           PPIZA=PPIZ(ISO3,ISO4)*ctors(n)
           if (iso1.eq.1.and.iso2.eq.2) ppiza=ppiza*0.1
           if (iso2.eq.1.and.iso1.eq.2) ppiza=ppiza*0.1
           if (iso1.eq.1.and.iso2.eq.3.and.j.le.na) ppiza=ppiza*0.1
           if (iso2.eq.1.and.iso1.eq.3.and.i.le.na) ppiza=ppiza*0.1
               IF (ISO1.GT.2.OR.ISO2.GT.2) PPIZA=PPIZA*1.5
               IF (ISO1.EQ.5.OR.ISO2.EQ.5) PPIZA=PPIZA*2.0
               if (iso1.eq.5.and.iso3.eq.20) ppiza=ppiza*0.67
               if (iso2.eq.5.and.iso4.eq.20) ppiza=ppiza*0.67
               PTOR=-0.3333*SIN(3.*TAU)*PPIZA+PTOR
               SPAN=PPIZA*(1.+COS(3.*TAU))*694.6
               IF (IMOL(I).EQ.99) SPAN=0.
               EPW=EPW+SPAN
               SPANT(I)=SPANT(I)+SPAN/2.
               SPANT(J)=SPANT(J)+SPAN/2.
C
C ORBITAL OVERLAP EFFECTS (ANOMERIC, GAUCHE,...)
C
  49           IF (ISO1.LE.3.AND.ISO2.LE.3) GOTO 50
               IF (ISO3.EQ.16.OR.ISO4.EQ.16) GOTO 50
               IF (ISO1.EQ.16.OR.ISO2.EQ.16) QTOR=0.5
               IF (QADD(K).LT.0.OR.QADD(L).LT.0) GOTO 50
  90           QI=QSEFF(I)+CS(I)
               QJ=QSEFF(J)+CS(J)
               S=-SIGN(QTOR,QI*QJ)
               DQIK=QSEFF(I)+CS(I)-QSEFF(K)-CS(K)
               DQLJ=QSEFF(L)+CS(L)-QSEFF(J)-CS(J)
               DDQ=DQIK*DQLJ*(1.-eps)
               EPW=EPW+S*DDQ*(1.+COS(2.*TAU))
               SPANT(I)=SPANT(I)+0.5*S*DDQ*(1.+COS(2.*TAU))
               SPANT(J)=SPANT(J)+0.5*S*DDQ*(1.+COS(2.*TAU))
               PTOR=PTOR-S*DDQ*0.5*SIN(2.*TAU)/1561.
        goto 50    
c
  140          PPIxA=-PPIx(ISO3,ISO4)*ctors(n)*(p(i,k)**2+p(j,l)**2)
               IF (ISO1.GT.2.OR.ISO2.GT.2) PPIxA=PPIxA*1.5
               IF (ISO1.EQ.5.OR.ISO2.EQ.5) PPIxA=PPIxA*2.0
               if (iso1.eq.5.and.iso3.eq.20) ppixa=ppixa*0.67
               if (iso2.eq.5.and.iso4.eq.20) ppixa=ppixa*0.67
               PTOR=-0.3333*SIN(3.*TAU)*PPIxA+PTOR
               SPAN=PPIxA*(1.+COS(3.*TAU))*694.6
               IF (IMOL(I).EQ.99) SPAN=0.
               EPW=EPW+SPAN
               SPANT(I)=SPANT(I)+SPAN/2.
               SPANT(J)=SPANT(J)+SPAN/2.
C
C    ADD COMPONENTS OF THE FORCES DUE TO TORSIONAL STRAIN TO THE
C    RESULTING FORCES
C
   50          IF (TFIX.GT.0.) THEN
               DELTAU=TAU-FTAU
C IF DELTAU IS LARGE, A CHANGE OF SIGN HAS OCURRED AND WE ARE HEADING
C IN THE WRONG DIRECTION (GOING ALMOST FULL-CIRCLE)
               IF (ABS(DELTAU).GT.1.)DELTAU=-SIGN(0.05,DELTAU)
               PTOR=DELTAU*.01
               END IF
               PFA=PTOR/VV(I,K)
               PFB=PTOR/VV(J,L)
               XK(I)=XK(I)-PFA*TAX
               YK(I)=YK(I)-PFA*TAY
               ZK(I)=ZK(I)-PFA*TAZ
               XK(J)=XK(J)+PFB*TBX
               YK(J)=YK(J)+PFB*TBY
               ZK(J)=ZK(J)+PFB*TBZ
               XK(K)=XK(K)+PFA*TAX
               YK(K)=YK(K)+PFA*TAY
               ZK(K)=ZK(K)+PFA*TAZ
               XK(L)=XK(L)-PFB*TBX
               YK(L)=YK(L)-PFB*TBY
               ZK(L)=ZK(L)-PFB*TBZ
  60           IF (VS(K,L).GT.0.) GOTO 70
               IF (VS(K,L).LT.-1.) GOTO 70
C
C   SET FLAGS TO NEGLECT THE NONBONDING INTERACTION K -- L
C
               VS(K,L)=-1.
               VS(L,K)=-1.
  70           RTORS(N)=-TTAU
  100 CONTINUE
      RETURN
      END
      SUBROUTINE BEND
C
C    THIS SUBROUTINE CALCULATES ENERGIES AND FORCES FOR BEND
C    ANGLES
C
      DOUBLE PRECISION XK,YK,ZK,virx,viry,virz
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      PARAMETER(IX=301,IY=IX,IX3=IX*3,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /OPTIM/virx,viry,virz,XD(IX3),XDD(IX3),QZ(IX3),DX(IX3),
     1XK(IX),YK(IX),ZK(IX),NAM,EPP,EPW
      COMMON /STRAIN/ SPANB(IX),SPANA(IX),SPANT(IX),SPANK(IX),SPANE(IX),
     1SPANV(IX)
      IF (IBEND.EQ.0) RETURN
      DO 1 I=1,NAA
      SPANK(I)=0.
    1 CONTINUE
C
      DO 100 N=1,IBEND
C
C
C    FIND THE ATOMS THAT DEFINE A TORSIONAL ANGLE I -- J -- K -- M
C
         I=NBEND(1,N)
         J=NBEND(2,N)
         K=NBEND(3,N)
         M=NBEND(4,N)
         ISO1=ISORT(I)+1
C
C    CALCULATE NORMAL VECTOR OF THE PLANE K - I - J
C
               TAX=(X(2,I)-X(2,J))*(X(3,I)-X(3,K))-(X(3,I)-X(3,J))*(X(2,
     1         I)-X(2,K))
               TAY=(X(3,I)-X(3,J))*(X(1,I)-X(1,K))-(X(1,I)-X(1,J))*(X(3,
     1         I)-X(3,K))
               TAZ=(X(1,I)-X(1,J))*(X(2,I)-X(2,K))-(X(2,I)-X(2,J))*(X(1,
     1         I)-X(1,K))
               TNA=SQRT(TAX*TAX+TAY*TAY+TAZ*TAZ)
               IF(TNA.EQ.0.)TNA=1.E-2
               TAX=TAX/TNA
               TAY=TAY/TNA
               TAZ=TAZ/TNA
                  COSBET=(TAX*(X(1,M)-X(1,I))+TAY*(X(2,M)-X(2,I))+TAZ*(X
     1            (3,M)-X(3,I)))/VV(I,M)
                  BETA=1.570796327-ACOSF(ABS(COSBET))
                  IF (ABS(BETA).LT.0.001) GOTO 80
                  DBEND=DBS(ISO1)
                  IF (I.GT.NA) GOTO 80
                  IF (ISO1.NE.3) GOTO 10
                  IF (P(I,I).LE.0.) GOTO 10
                  DBEND=0.1E-2*(2.-P(I,I))
C
C    CALCULATE STRAIN ENERGY DUE TO BENDING AND THE RESULTING
C    FORCES
C
   10             SPAN=DBEND*BETA**2*3121.
                 IF(IMOL(I).EQ.99) SPAN=0.
                  EPW=EPW+SPAN
                  SPANK(I)=SPANK(I)+SPAN
                  PBJ=BETA*DBEND/VV(I,J)*SIGN(0.66667,COSBET)
                  PBK=BETA*DBEND/VV(I,K)*SIGN(0.66667,COSBET)
                  PBM=BETA*DBEND/VV(I,M)*SIGN(0.66667,COSBET)
               TAX=(X(2,M)-X(2,J))*(X(3,M)-X(3,K))-(X(3,M)-X(3,J))*(X(2,
     1         M)-X(2,K))
               TAY=(X(3,M)-X(3,J))*(X(1,M)-X(1,K))-(X(1,M)-X(1,J))*(X(3,
     1         M)-X(3,K))
               TAZ=(X(1,M)-X(1,J))*(X(2,M)-X(2,K))-(X(2,M)-X(2,J))*(X(1,
     1         M)-X(1,K))
               TNA=SQRT(TAX*TAX+TAY*TAY+TAZ*TAZ)
               IF(TNA.EQ.0.)TNA=1.E-2
               TAX=TAX/TNA
               TAY=TAY/TNA
               TAZ=TAZ/TNA
                  XK(I)=XK(I)+(PBJ+PBK+PBM)*TAX
                  YK(I)=YK(I)+(PBJ+PBK+PBM)*TAY
                  ZK(I)=ZK(I)+(PBJ+PBK+PBM)*TAZ
                  XK(J)=XK(J)-PBJ*TAX
                  YK(J)=YK(J)-PBJ*TAY
                  ZK(J)=ZK(J)-PBJ*TAZ
                  XK(K)=XK(K)-PBK*TAX
                  YK(K)=YK(K)-PBK*TAY
                  ZK(K)=ZK(K)-PBK*TAZ
                  XK(M)=XK(M)-PBM*TAX
                  YK(M)=YK(M)-PBM*TAY
                  ZK(M)=ZK(M)-PBM*TAZ
C
C    BEND ANGLE DEFINED AS ANGLE BETWEEN PLANE I - J - K
C    AND VECTOR I - M
C
   80 RBEND(N)=BETA*57.29577951
  100 CONTINUE
      RETURN
      END
      SUBROUTINE NONBON
C
C    THIS SUBROUTINE CALCULATES ENERGIES AND FORCES FOR HYDROGEN
C    BONDS, COULOMB INTERACTIONS AND VAN DER WAALS INTERACTIONS
C
      DOUBLE PRECISION XK,YK,ZK,virx,viry,virz
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      PARAMETER(IX=301,IY=IX,IX3=IX*3,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /OPTIM/virx,viry,virz,XD(IX3),XDD(IX3),QZ(IX3),DX(IX3),
     1XK(IX),YK(IX),ZK(IX),NAM,EPP,EPW
      COMMON /STRAIN/ SPANB(IX),SPANA(IX),SPANT(IX),SPANK(IX),SPANE(IX),
     1SPANV(IX)
      COMMON /GAMMMA/ RR(IX,IX)
cs
      common /solvat/ eps,dielec,dsgl,fav(ix),fam(ix),rxy(ix,ix),
     1 rvdw(50),ion(50)
cs
C
C     SET EFFECTIVE ATOM CHARGE EQUAL SIGMA-CHARGE
C
      DO 10 I=1,NAA
      SPANE(I)=0.
      SPANV(I)=0.
      if (nm.gt.-6) QSEFF(I)=QSIG(I)
   10 continue 
C
C    FIND ELECTRONEGATIVE ATOM X
C
      DO 3230 I=1,NAA
         IF (IHB.EQ.1) GOTO 3230
         ISO1=ISORT(I)
         IF(ISO1.EQ.11) ISO1=2
         IF (ISO1.LT.2.OR.ISO1.GT.9) GOTO 3230
C
C    LOOK FOR A HYDROGEN ATOM IN SUITABLE DISTANCE FROM X
C    TO FORM X...H OF A HYDROGEN BOND OR A POLAR NONBONDING
C    INTERACTION
C
         DO 323 J=1,NAA
            IF (ISORT(J).NE.19) GOTO 323
            ISO2=ISORT(J)
            IF (VS(I,J).Lt.-1..OR.VS(I,J).GT.0.1) GOTO 323
 3235       DIST=RBOND(ISO1+1)+RBOND(ISO2+1)+1.6
            IF (VV(I,J).GT.DIST) GOTO 323
            IHBO=0
C
C    FIND ATOM Y BONDED TO H AND ACCEPT HYDROGEN BOND X...H-Y
C
            DO 321 K=1,NAA
               IF(VS(K,J).LE.0.) GOTO 321
               IF(VS(I,K).LT.0..OR.VS(I,K).GT.0.) GOTO 321
               ISO3=ISORT(K)
               IF (ISO3.LT.2.OR.ISO3.GT.8) GOTO 321
C
C CHECK ANGLE X - H - Y
C
            AKR=(VV(I,J)**2+VV(J,K)**2-VV(I,K)**2)/(2*VV(I,J)*VV(J,K))
            IF (AKR.GT.-0.25) GOTO 321
            IHBO=1
            GOTO 3211
  321       CONTINUE
            IF (IHBO.EQ.0) GOTO 323
 3211       zq=2.8/vv(i,j)*(-akr)*(1.-eps)
            CHARG=0.
            IF (QADD(K).NE.0.) CHARG=1.
C
C    INCREASE EFFECTIVE CHARGE TO ACCOUNT FOR POLAR NONBONDING
C    INTERACTION
C
            QSEFF(I)=QSEFF(I)-0.019*ZQ
            QSEFF(J)=QSEFF(J)+0.019*ZQ
C
C    SET FLAGS TO NEGLECT NONBONDING INTERACTIONS OF J,
C    IF I...J-K IS A HYDROGEN BOND
C
            VS(I,J)=-3.-CHARG
            VS(J,I)=-3.-CHARG
 3231 continue
  323    CONTINUE
 3230 CONTINUE
cs
c
c
c    electrostatics for solvation
c
c
      if (neps.eq.0) goto 840
      do 800 i=1,naa
      do 800 j=i,naa
      ccc=(qseff(i)+cs(i))*(qseff(j)+cs(j))
      span=-ccc*rxy(i,j)*eps
      epw=epw+span
      spane(i)=spane(i)+span*0.5
      spane(j)=spane(j)+span*0.5
c$$$$$$$$$
      if (i.ne.j) then
      ph=ccc*eps*(7.73E-7)*rxy(I,J)**3
     1*(1.-0.50/(fam(i)*fam(j))*
     1exp(-vv(i,j)**2*0.50/(fam(i)*fam(j))))
      XK(I)=XK(I)-PH*(X(1,I)-X(1,J))
      YK(I)=YK(I)-PH*(X(2,I)-X(2,J))
      ZK(I)=ZK(I)-PH*(X(3,I)-X(3,J))
      XK(J)=XK(J)+PH*(X(1,I)-X(1,J))
      YK(J)=YK(J)+PH*(X(2,I)-X(2,J))
      ZK(J)=ZK(J)+PH*(X(3,I)-X(3,J))
      end if
 800  continue
 840  continue
c
C
C    COULOMB INTERACTIONS
C
      DO 324 I=1,NAA
      DO 325 J=I+1,NAA
      CCC=QSEFF(I)*CS(J)+CS(I)*QSEFF(J)+QSEFF(I)*QSEFF(J)
      PH=(CCC+CS(I)*CS(J))*(-7.73E-7)*RR(I,J)**3
      SPAN=CCC*RR(I,J)
      EPW=EPW+SPAN
      SPANE(I)=SPANE(I)+SPAN/2.
      SPANE(J)=SPANE(J)+SPAN/2.
      XK(I)=XK(I)-PH*(X(1,I)-X(1,J))
      YK(I)=YK(I)-PH*(X(2,I)-X(2,J))
      ZK(I)=ZK(I)-PH*(X(3,I)-X(3,J))
      XK(J)=XK(J)+PH*(X(1,I)-X(1,J))
      YK(J)=YK(J)+PH*(X(2,I)-X(2,J))
      ZK(J)=ZK(J)+PH*(X(3,I)-X(3,J))
  325 CONTINUE
  324 CONTINUE
C
C    NONBONDING INTERACTIONS
C
      DO 401 I=1,NAM
      DO 400 J=I+1,NAA
      IF (VS(I,J).NE.0..AND.VS(I,J).NE.-4.) GOTO 400
         ISO1=ISORT(I)+1
         ISO2=ISORT(J)+1
         VINT=amax1(1.5,VV(I,J))
         IXP=IVDW(ISO1,ISO2)
         DXP=IXP
         EXX=(-BVDW(ISO1,ISO2)*VINT)
         AEXPB=AVDW(ISO1,ISO2)*EXP(EXX)
         CEXP=CVDW(ISO1,ISO2)*VINT**(-6)
         VEXP=VINT**(-IXP)
         IF (VS(I,J).EQ.-4..AND.NG.NE.4) THEN
         AEXPB=AEXPB*0.1
         CEXP=CEXP*0.1
         END IF
         SPAN=0.04334*(AEXPB*VEXP-CEXP)
         EPW=EPW+SPAN
         SPANV(I)=SPANV(I)+SPAN/2.
         SPANV(J)=SPANV(J)+SPAN/2.
         IF (NH.EQ.0) GOTO 410
         DO 420 JRST=1,IRZ
         IF (I.NE.IRST(1,JRST)) GOTO 420
         IF (J.NE.IRST(2,JRST)) GOTO 420
         IF (IRST(3,JRST).NE.0) GOTO 420
         PH=0.0033*(1.-ARST(JRST)/VINT)
         GOTO 450
  420    CONTINUE
  410   PH=0.6944E-5/VV(I,J)*(-AEXPB*BVDW(ISO1,ISO2)*VEXP+(-AEXPB*DXP
     1  *VEXP+6.*CEXP)/VINT)
  450    XK(I)=XK(I)-PH*(X(1,I)-X(1,J))
         YK(I)=YK(I)-PH*(X(2,I)-X(2,J))
         ZK(I)=ZK(I)-PH*(X(3,I)-X(3,J))
         XK(J)=XK(J)+PH*(X(1,I)-X(1,J))
         YK(J)=YK(J)+PH*(X(2,I)-X(2,J))
         ZK(J)=ZK(J)+PH*(X(3,I)-X(3,J))
  400 CONTINUE
  401 CONTINUE
      RETURN
C
      END
      FUNCTION ACOSF(X)
      IF (X.NE.0.) GOTO 2
      ACOSF=1.570796327
      GOTO 20
    2 IF (ABS(X).LT.1.) GOTO 5
      ACOSF=0.
      GOTO 10
    5 ACOSF=ATAN(SQRT(1.-X*X)/X)
   10 IF (X.LT.0.) ACOSF=ACOSF+3.141592654
   20 RETURN
      END
      SUBROUTINE SECTOR
C
C    CALCULATION OF TORSION ANGLES OF TYPE DELTA (I,J,K,L - I,L BOTH
C    BONDED TO J) AND TAU (I,J,K,L - I,L BOTH BONDED TO K)
C    SEE F.H. ALLEN AND D.ROGERS, ACTA CRYST B25, 1326 (1969).
C
      PARAMETER(IX=301,IY=IX,IX3=IX*3,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
C
C
C
C    FIND THE ATOMS THAT DEFINE A TORSION ANGLE OF TYPE DELTA OR TAU
C
      DO 100 N=1,IDELT
C
C
C
         I=NDELT(1,N)
         J=NDELT(2,N)
         K=NDELT(3,N)
         L=NDELT(4,N)
C
C    CALCULATE NORMAL VECTOR OF THE PLANE K - I - J
C
               TAX=(X(2,I)-X(2,J))*(X(3,I)-X(3,K))-(X(3,I)-X(3,J))*(X(2,
     1         I)-X(2,K))
               TAY=(X(3,I)-X(3,J))*(X(1,I)-X(1,K))-(X(1,I)-X(1,J))*(X(3,
     1         I)-X(3,K))
               TAZ=(X(1,I)-X(1,J))*(X(2,I)-X(2,K))-(X(2,I)-X(2,J))*(X(1,
     1         I)-X(1,K))
          TNA=SQRT(TAX*TAX+TAY*TAY+TAZ*TAZ)
          IF(TNA.EQ.0.)TNA=1.E-2
               TAX=TAX/TNA
               TAY=TAY/TNA
               TAZ=TAZ/TNA
C
C    CALCULATE NORMAL VECTOR OF THE PLANE I - J - L
C
               TBX=(X(2,I)-X(2,J))*(X(3,J)-X(3,L))-(X(3,I)-X(3,J))*(X(2,
     1         J)-X(2,L))
               TBY=(X(3,I)-X(3,J))*(X(1,J)-X(1,L))-(X(1,I)-X(1,J))*(X(3,
     1         J)-X(3,L))
               TBZ=(X(1,I)-X(1,J))*(X(2,J)-X(2,L))-(X(2,I)-X(2,J))*(X(1,
     1         J)-X(1,L))
        TNB=SQRT(TBX*TBX+TBY*TBY+TBZ*TBZ)
        IF(TNB.EQ.0.)TNB=1.E-2
               TBX=TBX/TNB
               TBY=TBY/TNB
               TBZ=TBZ/TNB
C
C    FIND THE VECTOR BISECTING THE TORSIONAL ANGLE TO DEFINE
C    SIGN OF THE TORSIONAL ANGLE
C
               ZAX=TAY*TBZ-TAZ*TBY
               ZAY=TAZ*TBX-TAX*TBZ
               ZAZ=TAX*TBY-TAY*TBX
               SZ=ZAX*(X(1,I)-X(1,J))+ZAY*(X(2,I)-X(2,J))+ZAZ*(X(3,I)-X(
     1         3,J))
               SIGTAU=SIGN(1.,SZ)
C
C    EVALUATE TORSIONAL ANGLE DELTA OR TAU
C
               COSTAU=TAX*TBX+TAY*TBY+TAZ*TBZ
               IF (COSTAU.LT.-1.) COSTAU=-1.
               IF (COSTAU.GT.1.) COSTAU=1.
               TAU=-ACOSF(COSTAU)*SIGTAU
  100 RDELT(N)=TAU
      RETURN
      END
      SUBROUTINE OUTGEO
C
C     OUTPUT OF MOLECULAR GEOMETRY
C
      CHARACTER*2 ATA,ATN,W1(4),W2(4),W3(4),W4(4)
      COMMON /LABEL/ ATA(50)
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      PARAMETER(IX=301,IY=IX,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /STRAIN/ SPANB(IX),SPANA(IX),SPANT(IX),SPANK(IX),SPANE(IX),
     1SPANV(IX)
      DIMENSION N1(4),N2(4),N3(4),N4(4),R(4)
      IF (NC.LT.0) GOTO 999
C
C     CALCULATE TORSION ANGLES OF TYPE DELTA AND TAU
C
      CALL SECTOR
C
C     OUTPUT OF CARTESIAN COORDINATES
C
      WRITE (IUT,100)
      DO 10 I=1,NAA
      IS=ISORT(I)+1
      CSUM=CS(I)+QSEFF(I)
      ATN=ATA(IS)
      WRITE (IUT,110) ATN,NUMM(I),ISORT(I),(X(J,I),J=1,3),
     1CSUM,QSEFF(I),CS(I)
   10 CONTINUE
C
C     OUTPUT OF BOND LENGTHS
C
      WRITE (IUT,120)
      NPRINT=0
      DO 20 I=1,NAA
      DO 20 J=I,NAA
      IF(VS(I,J).LE.0.) GOTO 20
      NPRINT=NPRINT+1
      ISO1=ISORT(I)+1
      ISO2=ISORT(J)+1
      W1(NPRINT)=ATA(ISO1)
      W2(NPRINT)=ATA(ISO2)
      N1(NPRINT)=NUMM(I)
      N2(NPRINT)=NUMM(J)
      R(NPRINT)=VV(I,J)
      IF (NPRINT.LT.4) GOTO 20
      WRITE(IUT,130) (W1(M),N1(M),W2(M),N2(M),R(M),M=1,4)
      NPRINT=0
   20 CONTINUE
      WRITE(IUT,130) (W1(M),N1(M),W2(M),N2(M),R(M),M=1,NPRINT)
C
C     OUTPUT OF BOND ANGLES
C
      WRITE (IUT,140)
      KP=0
      DO 30 N=1,IANG
      I=IABS(NANG(1,N))
      J=NANG(2,N)
      K=NANG(3,N)
      I1=ISORT(I)+1
      I2=ISORT(J)+1
      I3=ISORT(K)+1
      KP=KP+1
      W2(KP)=ATA(I1)
      W1(KP)=ATA(I2)
      W3(KP)=ATA(I3)
      N2(KP)=NUMM(I)
      N1(KP)=NUMM(J)
      N3(KP)=NUMM(K)
      R(KP)=RANG(N)*57.29577951
      IF (KP.LT.3) GOTO 30
      WRITE(IUT,150) (W1(M),N1(M),W2(M),N2(M),W3(M),N3(M),R(M),M=1,3)
      KP=0
   30 CONTINUE
      WRITE(IUT,150) (W1(M),N1(M),W2(M),N2(M),W3(M),N3(M),R(M),M=1,KP)
C
C     OUTPUT OF TORSION ANGLES
C
      WRITE (IUT,160)
      KP=0
      DO 40 N=1,ITORS
         I=NTORS(1,N)
         J=NTORS(2,N)
         K=NTORS(3,N)
         L=NTORS(4,N)
      I1=ISORT(I)+1
      I2=ISORT(J)+1
      I3=ISORT(K)+1
      I4=ISORT(L)+1
      KP=KP+1
      W2(KP)=ATA(I1)
      W3(KP)=ATA(I2)
      W1(KP)=ATA(I3)
      W4(KP)=ATA(I4)
      N2(KP)=NUMM(I)
      N3(KP)=NUMM(J)
      N1(KP)=NUMM(K)
      N4(KP)=NUMM(L)
      R(KP)=RTORS(N)*57.29577951
      IF (KP.LT.2) GOTO 40
      WRITE(IUT,170) (W1(M),N1(M),W2(M),N2(M),W3(M),N3(M),W4(M),N4(M),
     1R(M),M=1,2)
      KP=0
   40 CONTINUE
      WRITE(IUT,170) (W1(M),N1(M),W2(M),N2(M),W3(M),N3(M),W4(M),N4(M),
     1R(M),M=1,KP)
C
C     OUTPUT OF TORSION ANGLES OF TYPE DELTA AND TAU
C
      WRITE (IUT,180)
      KP=0
      DO 50 N=1,IDELT
         I=NDELT(1,N)
         J=NDELT(2,N)
         K=NDELT(3,N)
         L=NDELT(4,N)
      I1=ISORT(I)+1
      I2=ISORT(J)+1
      I3=ISORT(K)+1
      I4=ISORT(L)+1
      KP=KP+1
      W2(KP)=ATA(I1)
      W3(KP)=ATA(I2)
      W1(KP)=ATA(I3)
      W4(KP)=ATA(I4)
      N2(KP)=NUMM(I)
      N3(KP)=NUMM(J)
      N1(KP)=NUMM(K)
      N4(KP)=NUMM(L)
      R(KP)=RDELT(N)*57.29577951
      IF (KP.LT.2) GOTO 50
      WRITE(IUT,190) (W1(M),N1(M),W2(M),N2(M),W3(M),N3(M),W4(M),N4(M),
     1R(M),M=1,2)
      KP=0
   50 CONTINUE
      WRITE(IUT,190) (W1(M),N1(M),W2(M),N2(M),W3(M),N3(M),W4(M),N4(M),
     1R(M),M=1,KP)
C
C     SEARCH FOR HYDROGEN BONDS AND OUTPUT
C
C
C    FIND ATOM X
C
  999 DO 90 I=1,NAA
         I1=ISORT(I)+1
         IF (I1.EQ.12) I1=3
         IF (I1.LT.3.OR.I1.GT.9) GOTO 90
C
C    LOOK FOR HYDROGEN ATOM TO FORM A HYDROGEN BOND X...H
C
        DO 80 J=1,NAA
            IHBO=0
            IF (VS(I,J).GT.-3.) GOTO 80
            I2=ISORT(J)+1
            DIST=RBOND(I1)+RBOND(I2)+1.6
            IF (VV(I,J).GT.DIST) GOTO 80
C
C    FIND ATOM Y BONDED TO H AND ACCEPT HYDROGEN BOND X...H-Y
C
            DO 70 K=1,NAA
               IF(VS(K,J).LE.0.) GOTO 70
               IHBO=1
               I3=ISORT(K)+1
            AKR=(VV(I,J)**2+VV(J,K)**2-VV(I,K)**2)/(2*VV(I,J)*VV(J,K))
            ANG=ACOSF(AKR)/3.1415926*180.
            if (ang.lt.110.) goto 80
               GOTO 71
   70       CONTINUE
   71       IF (IHBO.EQ.0) GOTO 80
      WRITE (IUT,200)
      WRITE (IUT,210) ATA(I1),NUMM(I),ATA(I2),NUMM(J),ATA(I3),
     1NUMM(K),VV(I,J),VV(J,K),VV(I,K),ANG
   80 CONTINUE
   90 CONTINUE
C
C SEARCH FOR METAL ATOM AND OUTPUT OF COORDINATION SPHERE
C
      DO 300 I=1,NAA
      IS=ISORT(I)+1
      IF(IS.LT.21.AND.IS.NE.13.AND.IS.NE.14) THEN
      IF(IS.NE.18.AND.IS.NE.19) GOTO 300
      END IF
      IF (IS.GE.36) GOTO 300
      WRITE(IUT,310) ATA(IS),NUMM(I)
      KP=0
      DO 320 J=1,NAA
      IF (VV(I,J).GT.4.0) GOTO 320
      IF (I.EQ.J) GOTO 320
      ISJ=ISORT(J)+1
      KP=KP+1
      W1(KP)=ATA(IS)
      W2(KP)=ATA(ISJ)
      N1(KP)=NUMM(I)
      N2(KP)=NUMM(J)
      R(KP)=VV(I,J)
      IF(KP.LT.4) GOTO 320
      WRITE (IUT,330) (W1(M),N1(M),W2(M),N2(M),R(M),M=1,4)
      KP=0
  320 CONTINUE
      WRITE(IUT,330) (W1(M),N1(M),W2(M),N2(M),R(M),M=1,KP)
  300 CONTINUE
C
C    OUTPUT TERMES OF STRAIN ENERGY
C
      SSPANB=0.
      SSPANA=0.
      SSPANT=0.
      SSPANK=0.
      SSPANE=0.
      SSPANV=0.
      DO 4711 I=1,NAA
      SSPANB=SSPANB+SPANB(I)
      SSPANA=SSPANA+SPANA(I)
      SSPANT=SSPANT+SPANT(I)
      SSPANK=SSPANK+SPANK(I)
      SSPANE=SSPANE+SPANE(I)
      SSPANV=SSPANV+SPANV(I)
 4711 CONTINUE
      SPGES=SSPANB+SSPANA+SSPANT+SSPANK+SSPANE+SSPANV
      SSPANB=SSPANB*96.533
      SSPANA=SSPANA*96.533
      SSPANT=SSPANT*96.533
      SSPANK=SSPANK*96.533
      SSPANE=SSPANE*96.533
      SSPANV=SSPANV*96.533
      WRITE(IUT,220) SPGES,SSPANB,SSPANA,SSPANT,SSPANK,SSPANE,SSPANV
cp
      if (nm.lt.-5) then
      call crtcry
      call outcry
      end if
      IF (NC.LT.1) RETURN
      WRITE (IUT,800)
      DO 410 I=1,NAA
      ATN=ATA(ISORT(I)+1)
      WRITE (IUT,810) ATN,NUMM(I),ISORT(I),SPANB(I),SPANA(I),
     1SPANT(I),SPANK(I),SPANE(I),SPANV(I)
  410 CONTINUE

C
C
C    FORMATES
C
  100 FORMAT (1H ,'NEW COORDINATES AND CHARGES'/1H ,'    I    ISORT(I)'
     1,'    X         Y         Z       CHARGE   SIGMA-  PI-CHARGE')
  110 FORMAT (2X,A2,I3,I8,7F10.4)
  120 FORMAT (1H ,'      BOND LENGTHS')
  130 FORMAT(' ',4(A2,I3,A2,I3,F8.3,2X))
  140 FORMAT (1H ,'      BOND ANGLES')
  150 FORMAT (1H ,3(A2,I3,1H-,A2,I3,1H-,A2,I3,F6.1,5X))
  160 FORMAT (1H ,'      TORSION ANGLES')
  170 FORMAT (1H ,2(A2,I3,1H-,A2,I3,1H-,A2,I3,1H-,A2,I3,F7.1,5X))
  180 FORMAT (1H ,'      TORSION ANGLES TYPES DELTA AND TAU')
  190 FORMAT (1H ,2(A2,I3,1H-,A2,I3,1H-,A2,I3,1H-,A2,I3,F7.1,5X))
  200 FORMAT (1H ,10X,'HYDROGEN BOND X...H-Y'/1H ,14X,1HX,4X,1HH,4X,1HY,
     1        3X,'R(X..H)',4X,'R(H-Y)',3X,'R(X..Y)',4X,'ANGLE X..H-Y')
  210 FORMAT (1H ,6X,3(1X,A2,I3),4F10.3)
  310 FORMAT(24HCOORDINATION SPHERE FOR ,A2,I3)
  330 FORMAT(4(1H ,A2,I3,1H-,A2,I3,F7.3,1X))
  220 FORMAT ('  SPLIT TERMS OF STRAIN ENERGY (ES=',F8.4,' EV)',
     1' IN KJ/MOL:',/,'   E(BONDS)  =',F8.1,/,'   E(ANGLES) =',F8.1,/
     2,'   E(TWIST)  =',F8.1,/,'   E(BEND)   =',F8.1,/,'   E(COULOMB)=',
     3 F8.1,/,'   E(VDW)    =',F8.1,/)
  800 FORMAT (' REMAINING STRAIN ON INDIVIDUAL ATOMS',/,
     1'ATOM    BOND   ANGLE  TORSION  BENDING  COULOMB  VDW')
  810 FORMAT (A2,I3,I8,7F10.5)
      RETURN
      END
      SUBROUTINE CORSET
C
C     SAVE NEW INTERNAL COORDINATES FOR NEXT CALCULATION OR OUTPUT
C
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      PARAMETER(IX=301,IY=IX,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM2/NUMAT,IATD(IX),NFIX(IX),ILES,NNAT(IX),NAD
     1(IX),NBD(IX),NCD(IX),NDD(IX),NCTOR,MCTOR,NPERM(20),ITOR(IX),
     2NST(IX),NMOL(IX),KTOR
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
C
      IF(NAD(1).EQ.0.AND.NAD(2).EQ.0)CALL NABC(NAA,NAD,NBD,NCD,VV,NNAT)
C
      DO 10 I1=1,NAA
         J=NAD(I1)
         IF (J.EQ.0) GOTO 10
         I=NNAT(I1)
         J=NNAT(J)
         AB(1,I1)=VV(I,J)
   10 CONTINUE
C
      DO 20 J1=1,NAA
         I=NAD(J1)
         K=NBD(J1)
         IF (I.EQ.0.OR.K.EQ.0) GOTO 20
         I=NNAT(I)
         J=NNAT(J1)
         K=NNAT(K)
         VR=(X(1,I)-X(1,J))*(X(1,I)-X(1,K))+(X(2,I)-X(2,J))*(X(2,I)-X(2,
     1   K))+(X(3,I)-X(3,J))*(X(3,I)-X(3,K))
         COSG=VR/(VV(I,J)*VV(I,K))
         AB(2,J1)=ACOSF(COSG)
   20 CONTINUE
C
      DO 30 K1=1,NAA
         I=NAD(K1)
         J=NBD(K1)
         L=NCD(K1)
         IF (I.EQ.0.OR.J.EQ.0.OR.L.EQ.0) GOTO 30
         I=NNAT(I)
         J=NNAT(J)
         K=NNAT(K1)
         L=NNAT(L)
         TAX=(X(2,I)-X(2,J))*(X(3,I)-X(3,K))-(X(3,I)-X(3,J))*(X(2,I)-X(2
     1   ,K))
         TAY=(X(3,I)-X(3,J))*(X(1,I)-X(1,K))-(X(1,I)-X(1,J))*(X(3,I)-X(3
     1   ,K))
         TAZ=(X(1,I)-X(1,J))*(X(2,I)-X(2,K))-(X(2,I)-X(2,J))*(X(1,I)-X(1
     1   ,K))
          TNA=SQRT(TAX*TAX+TAY*TAY+TAZ*TAZ)
          IF(TNA.EQ.0.)TNA=1.E-2
         TAX=TAX/TNA
         TAY=TAY/TNA
         TAZ=TAZ/TNA
         TBX=(X(2,I)-X(2,J))*(X(3,J)-X(3,L))-(X(3,I)-X(3,J))*(X(2,J)-X(2
     1   ,L))
         TBY=(X(3,I)-X(3,J))*(X(1,J)-X(1,L))-(X(1,I)-X(1,J))*(X(3,J)-X(3
     1   ,L))
         TBZ=(X(1,I)-X(1,J))*(X(2,J)-X(2,L))-(X(2,I)-X(2,J))*(X(1,J)-X(1
     1   ,L))
         TNB=SQRT(TBX*TBX+TBY*TBY+TBZ*TBZ)
         IF(TNB.EQ.0.)TNB=1.E-2
         TBX=TBX/TNB
         TBY=TBY/TNB
         TBZ=TBZ/TNB
         ZAX=TAY*TBZ-TAZ*TBY
         ZAY=TAZ*TBX-TAX*TBZ
         ZAZ=TAX*TBY-TAY*TBX
         SZ=ZAX*(X(1,I)-X(1,J))+ZAY*(X(2,I)-X(2,J))+ZAZ*(X(3,I)-X(3,J))
         SIGTAU=SIGN(1.,SZ)
         COSTAU=TAX*TBX+TAY*TBY+TAZ*TBZ
         IF (COSTAU.LT.-1.) COSTAU=-1.
         IF (COSTAU.GT.1.) COSTAU=1.
         AB(3,K1)=-ACOSF(COSTAU)*SIGTAU
   30 CONTINUE
C
C     IF (NH.le.0) RETURN
C
      DO 40 I=1,NAA
      IF (A(1,I).NE.0..AND.NH.le.0.) GOTO 40
         A(1,I)=AB(1,I)
         A(2,I)=AB(2,I)
         A(3,I)=-AB(3,I)
   40 CONTINUE
      RETURN
      END
      SUBROUTINE NABC(NAT,NA,NB,NC,V,N)
C
C    THIS ROUTINE SETS THE DEFINITION OF INTERNAL GEOMETRY
C
      PARAMETER(IX=301)
      DIMENSION NA(IX),NB(IX),NC(IX),V(IX,IX),N(IX)
      DO 10 I=1,NAT
         NA(I)=0
         NB(I)=0
   10 NC(I)=0
         NA(2)=1
         NA(3)=2
         NB(3)=1
      DO 20 I=4,NAT
         NI=I
         NA(I)=NEXT(0,0,NI,V,N)
         IF (NA(I).EQ.0) GOTO 20
         NAI=NA(I)
         NB(I)=NA(NAI)
         IF (NB(I).NE.0) GOTO 12
         NB(I)=NEXT(NAI,0,NI,V,N)
         IF (NB(I).EQ.0) GOTO 20
   12    NBI=NB(I)
         NC(I)=NB(NAI)
         IF (NC(I).NE.0) GOTO 20
         NC(I)=NEXT(NAI,NBI,NI,V,N)
   20 CONTINUE
      RETURN
      END
      FUNCTION NEXT(NA,NB,K,R,N)
C
C    FINDS ATOM NEXT TO K
C
      PARAMETER (IX=301)
      DIMENSION R(IX,IX),N(IX)
      KM=K-1
      DO 10 J=1,KM
         IF (J.NE.NA.AND.J.NE.NB) GOTO 15
   10 CONTINUE
   15 DO 20 I=1,KM
      IF (I.EQ.NA.OR.I.EQ.NB) GOTO 20
      IR=N(I)
      JR=N(J)
      KR=N(K)
      IF (R(IR,KR).LT.R(JR,KR).AND.K.NE.I) J=I
   20 CONTINUE
      NEXT=J
      END
      SUBROUTINE COROUT
C
C    THIS SUBROUTINE PREPARES PIMM-INPUT WITH REFINED PARAMETERS ON FILE
C    IPU.
C
      double precision xkc,ykc,zkc,vircx,vircy,vircz,viryx,virzx,virzy
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      PARAMETER(IX=301,IY=IX,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM2/NUMAT,IATD(IX),NFIX(IX),ILES,NNAT(IX),NAD
     1(IX),NBD(IX),NCD(IX),NDD(IX),NCTOR,MCTOR,NPERM(20),ITOR(IX),
     2NST(IX),NMOL(IX),KTOR
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /PDBNAM/RESNA(IX),NURES(IX),RESNAM(IX),NUMRES(IX),
     1RESTY(IX),RESTYP(IX)
cp
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1 cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2 imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3 ,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
cp
      CHARACTER*1 FIX(3)
      CHARACTER*5 RESTY,RESTYP
      CHARACTER*4 RESNA,RESNAM
      DATA FIX /' ','1','2'/
      IS=0
C
C    CALCULATION NUMBER IS OUTPUT ON FIRST TITLE CARD
C
      WRITE (IPU,100) IS,ISER,DHFM,(TITLE(I),I=1,36)
      WRITE (IPU,200) NA,NB,NC,ND,ZBB,NG,NH,NI,NK,NL,NM,IGEO,NE,NEPS
      if (nm.gt.-6) then
      DO 10 I1=1,NAA
      I=NNAT(I1)
      II=IFIX(I)+1
   10 WRITE (IPU,300) IATD(I1),NST(I1),FIX(II),(X(J,I),J=1,3)
      WRITE (IPU,400)
      else
      calp=cal*57.2957795
      cbet=cbe*57.2957795
      cgam=cga*57.2957795
      write(ipu,9000) (ck(i),i=1,3),calp,cbet,cgam
 9000 format('CELL',6f10.4)
      write(ipu,9001) ksym,iscl
 9001 format(2i4)
      do 9002 k=1,ksym
      write(ipu,9003) (rsym(1,j,k),j=1,3),tsym(1,k),
     1(rsym(2,j,k),j=1,3),tsym(2,k),(rsym(3,j,k),j=1,3),tsym(3,k)
 9002 continue
 9003 Format('SYMM', 3(3f4.0,f8.4))
      DO 9010 I1=1,NAA
      I=NNAT(I1)
      II=IFIX(I)+1
 9010 WRITE (IPU,300) IATD(I1),NST(I1),FIX(II),(CX(J,I),J=1,3)
      WRITE (IPU,400)
      end if
      IF (NG.EQ.4) RETURN
      IF (NAD(1).EQ.0.AND.NAD(2).EQ.0) RETURN
      WRITE (IPU,500) IATD(1),NST(1),FIX(IFIX(NNAT(1))+1)
      WRITE (IPU,500) IATD(2),NST(2),FIX(IFIX(NNAT(2))+1),AB(1,2)
      W=AB(2,3)*57.29577951
      WRITE (IPU,500) IATD(3),NST(3),FIX(IFIX(NNAT(3))+1),AB(1,3),W
      IF (NAA.LT.4) GOTO 40
      DO 30 I=4,NUMAT
      W=AB(2,I)*57.29577951
      WW=AB(3,I)*57.29577951
   30 WRITE (IPU,500) IATD(I),NST(I),FIX(IFIX(NNAT(I))+1),AB(1,I),W,WW,
     1                NAD(I),NBD(I),NCD(I)
   40 WRITE (IPU,400)
      RETURN
C
  100 FORMAT(I3,1X,I4,F8.1,16A4/20A4)
  200 FORMAT(4I3,F6.2,10I3)
  300 FORMAT(2I3,2X,A1,1X,3F10.4,1X,4I3,4X,A5,1X,A4,I4)
  400 FORMAT(/)
  500 FORMAT(I2,I3,2X,A1,2X,3(F10.5,10X),3I3)
C
      END
      SUBROUTINE BOSYM(ISBO,JSBO,IPA,IPB,P)
C
C    THIS SUBROUTINE SETS MEAN VALUES FOR SELECTED MATRIX
C    ELEMENTS P(I,J)
C
      PARAMETER(IX=301,IY=IX)
      DIMENSION JSBO(20),IPA(20,20),IPB(20,20),P(IY,IY)
      DO 30 I=1,ISBO
      PS=0.
      I1=JSBO(I)
      DO 10 J=1,I1
      I2=IPA(I,J)
      I3=IPB(I,J)
  10  PS=PS+P(I2,I3)
      PS=PS/FLOAT(I1)
      DO 20 J=1,I1
      I2=IPA(I,J)
      I3=IPB(I,J)
      P(I2,I3)=PS
  20  P(I3,I2)=PS
  30  CONTINUE
      RETURN
      END
      SUBROUTINE QSYM(ISBO,JSBO,IPA,IPB,C)
C
C    THIS SUBROUTINE SETS MEAN VALUES FOR SELECTED SIGMA
C    CHARGES C(I)
C
      PARAMETER(IX=301)
      DIMENSION JSBO(20),IPA(20,20),IPB(20,20),C(IX)
      DO 30 I=1,ISBO
      PS=0.
      I1=JSBO(I)
      DO 10 J=1,I1
      I2=IPA(I,J)
      I3=IPB(I,J)
      IF (I2.NE.I3) GOTO 10
      PS=PS+C(I2)
  10  CONTINUE
      PS=PS/FLOAT(I1)
      DO 20 J=1,I1
      I2=IPA(I,J)
      I3=IPB(I,J)
      IF (I2.NE.I3) GOTO 20
      C(I2)=PS
  20  CONTINUE
  30  CONTINUE
      RETURN
      END
      SUBROUTINE CHARGE
C
C    EMPIRICAL EVALUATION OF SIGMA-CHARGES USING THE
C    FORMALISM OF M.MARSILLI AND J. GASTEIGER.
C    TETRAHEDRON 36,3219-3228 (1980)
C
      double precision xkc,ykc,zkc,vircx,vircy,vircz,viryx,virzx,virzy
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      PARAMETER(IX=301,IY=IX,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /GAMMMA/ RR(IX,IX)
cp
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1 cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2 imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3 ,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
cp
cs
      common /solvat/ eps,dielec,dsgl,fav(ix),fam(ix),rxy(ix,ix),
     1 rvdw(50),ion(50)
cs
      DIMENSION XE(IX),QK(IX),qsa(ix)
C
C   SET CORE CHARGES AS INITIAL VALUES
C
      DO 10 I=1,NAA
      qsa(i)=qsig(i)
   10 QSIG(I)=QADD(I)
      DO 40 K=1,10
C
C    CALCULATE ELECTRONEGATIVITIES FOR NEXT ITERATION STEP
C
         DO 20 I=1,NAA
            QK(I)=0.
            ISO1=ISORT(I)+1
            XE(I)=(XN(ISO1)+B(ISO1)*QSIG(I)+C(ISO1)*QSIG(I)*QSIG(I))
     1      *(1.+SWI(I)*0.03333/XN(ISO1))
            xe(i)=xe(i)+(cs(i)+qsig(i))*rr(i,i)
     1      -0.5*qsig(i)*eps*rxy(i,i)
         DO 20 J=1,NAA
            if (i.eq.j) goto 20
            rrij=rr(i,j)
            if (isort(j).eq.13) rrij=rr(i,j)*(1.-fav(j))
            XE(I)=XE(I)+(CS(J)+QSIG(J))*RRIJ
   20     CONTINUE
cp
      if (nm.lt.-5) then
      do 21 i=1,naa
   21 xe(i)=xe(i)+xpac(i)
      end if
cp
C
C    CALCULATE CHARGE TRANSFER IN THE SIGMA-SKELETON
C
      DO 40 I=1,NAA
      nq=0.
         DO 30 J=I,NAA
         IF (VS(I,J).LE.0.) GOTO 30
         IF (ISORT(I).EQ.12.OR.ISORT(J).EQ.12) GOTO 30
         IF (ISORT(I).EQ.13.OR.ISORT(J).EQ.13) GOTO 30
 4711    ISO2=ISORT(I)+1
         IF (XE(I).GT.XE(J)) ISO2=ISORT(J)+1
            QIJ=(XE(J)-XE(I))/XP(ISO2)*0.5**K
            QK(I)=QK(I)+QIJ
            QK(J)=QK(J)-QIJ
   30    CONTINUE
         QSIG(I)=QSIG(I)+QK(I)
         IF (K.EQ.1) QSIG(I)=QK(I)
   40 CONTINUE
C
C   ADD CORE CHARGES TO FINAL CHARGE-DISTRIBUTION
C
 555  DO 50 I=1,NAA
      QSIG(I)=QSIG(I)+QADD(I)
   50 continue 
      RETURN
      END
      SUBROUTINE SETHYD(JJ,IBGR,I,IBO,ISO,IWARN)
C
C    THIS SUBROUTINE ADDS HYDROGEN ATOMS
C
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      PARAMETER(IX=301,IY=IX,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM2/NUMAT,IATD(IX),NFIX(IX),ILES,NNAT(IX),NAD
     1(IX),NBD(IX),NCD(IX),NDD(IX),NCTOR,MCTOR,NPERM(20),ITOR(IX),
     2NST(IX),NMOL(IX),KTOR
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /PDBNAM/RESNA(IX),NURES(IX),RESNAM(IX),NUMRES(IX)
     1,RESTY(IX),RESTYP(IX)
      CHARACTER*4 RESNAM,RESNA
      CHARACTER*5 RESTYP,RESTY    
      DIMENSION VECPRO(3),VECSUM(3),IBGR(50),JJ(1000)
C
      WRITE(IUT,61) IBGR(ISO)-IBO
   61 FORMAT(1H ,I1,' HYDROGEN ATOM(S) WILL BE GENERATED')
      IF (IBGR(ISO)-IBO-2) 100,110,120
C
C     ONE BOND MISSING ON SP,SP2 OR SP3 CARBON
C     OR ON PYRROLE OR AMINE TYPE NITROGEN
C
  100 IF (ISO.EQ.5) GOTO 120
      IF (ISO.EQ.11.OR.NST(NUMM(I)).EQ.200) RETURN
      NAA=NAA+1
      NUMM(NAA)=NAA
      NNAT(NAA)=NAA
      NUMAT=NUMAT+1
      IBOND=IBOND+1
C
C     CALCULATE AND MIRROR THE SUM OF THE BOND VECTORS
C
      VLANG=0.
      DO 101 K=1,3
      IF (ISO.GT.1) X(K,NAA)=-(X(K,JJ(1))+X(K,JJ(2))-2*X(K,I))
      IF (ISO.EQ.1) X(K,NAA)=-(X(K,JJ(1))+X(K,JJ(2))+X(K,JJ(3))
     1-3*X(K,I))
      IF (ISO.EQ.11) X(K,NAA)=-(X(K,JJ(1))-X(K,I))
      VLANG=VLANG+X(K,NAA)**2
  101 CONTINUE
C
C     REDUCE LENGTH OF NEW BOND TO 1.1 ANGSTROMS
C
      VCORR=1.1/SQRT(VLANG)
      DO 102 K=1,3
      X(K,NAA)=X(K,NAA)*VCORR+X(K,I)
  102 CONTINUE
      ISORT(NAA)=19
      VS(I,NAA)=1.1
      VS(NAA,I)=1.1
      NST(NAA)=19
      IAT(NAA)=1
      IATD(NAA)=1
      IMOL(NAA)=IMOL(I)
      RESNAM(NAA)=RESNAM(I)
      RESTYP(NAA)=' H   '
      NUMRES(NAA)=NUMRES(I)
      NAD(NAA)=NUMM(I)
      NBD(NAA)=NUMM(JJ(1))
      NCD(NAA)=NUMM(JJ(2))
      DO 103 K=1,NAA-1
      VV(K,NAA)=SQRT((X(1,K)-X(1,NAA))**2+(X(2,K)-X(2,NAA))**2
     1+(X(3,K)-X(3,NAA))**2)
      VV(NAA,K)=VV(K,NAA)
  103 CONTINUE
      GO TO 39
C
C     TWO BONDS MISSING ON SP3 CARBON
C
  110 IF (ISO.GT.1.AND.ISO.NE.15) GO TO 120
  104 DO 111 K=1,3
C
C     CALCULATE SUM AND PRODUCT OF THE TWO BOND VECTORS
C
      VECSUM(K)=-(X(K,JJ(1))+X(K,JJ(2))-2*X(K,I))
  111 CONTINUE
      VECPRO(1)=(X(2,I)-X(2,JJ(1)))*(X(3,I)-X(3,JJ(2)))
      VECPRO(1)=VECPRO(1)-(X(2,I)-X(2,JJ(2)))*(X(3,I)-X(3,JJ(1)))
      VECPRO(2)=(X(3,I)-X(3,JJ(1)))*(X(1,I)-X(1,JJ(2)))
      VECPRO(2)=VECPRO(2)-(X(3,I)-X(3,JJ(2)))*(X(1,I)-X(1,JJ(1)))
      VECPRO(3)=(X(1,I)-X(1,JJ(1)))*(X(2,I)-X(2,JJ(2)))
      VECPRO(3)=VECPRO(3)-(X(1,I)-X(1,JJ(2)))*(X(2,I)-X(2,JJ(1)))
      NAA=NAA+2
      IBOND=IBOND+2
      NUMAT=NUMAT+2
      VLANG=0.
      DO 112 K=1,3
C
C     THE MISSING BONDS LIE APPROXIMATELY IN THE DIRECTION OF
C     THE SUM AND DIFFERENCE OF THESE TWO VECTORS
C
      X(K,NAA-1)=VECSUM(K)+VECPRO(K)
      X(K,NAA)=VECSUM(K)-VECPRO(K)
      VLANG=VLANG+X(K,NAA)**2
  112 CONTINUE
      VCORR=1.1/SQRT(VLANG)
      DO 113 K=1,3
      X(K,NAA-1)=X(K,NAA-1)*VCORR+X(K,I)
      X(K,NAA)=X(K,NAA)*VCORR+X(K,I)
  113 CONTINUE
      NNAT(NAA-1)=NAA-1
      NNAT(NAA)=NAA
      ISORT(NAA-1)=19
      ISORT(NAA)=19
      IAT(NAA-1)=1
      IAT(NAA)=1
      IATD(NAA)=1
      IATD(NAA-1)=1
      IMOL(NAA)=IMOL(I)
      IMOL(NAA-1)=IMOL(I)
      RESNAM(NAA-1)=RESNAM(I)
      RESNAM(NAA)=RESNAM(I)
      NUMRES(NAA-1)=NUMRES(I)
      NUMRES(NAA)=NUMRES(I)
      RESTYP(NAA-1)=' H   '
      RESTYP(NAA)=' H  '
      VS(I,NAA-1)=1.1
      VS(NAA-1,I)=1.1
      DO 114 K=1,NAA-2
      VV(K,NAA-1)=SQRT((X(1,K)-X(1,NAA-1))**2+(X(2,K)
     1-X(2,NAA-1))**2+(X(3,K)-X(3,NAA-1))**2)
      VV(NAA-1,K)=VV(K,NAA-1)
  114 CONTINUE
      VS(I,NAA)=1.1
      VS(NAA,I)=1.1
      DO 115 K=1,NAA-1
      VV(K,NAA)=SQRT((X(1,K)-X(1,NAA))**2+(X(2,K)-X(2,NAA))**2
     1+(X(3,K)-X(3,NAA))**2)
      VV(NAA,K)=VV(K,NAA)
  115 CONTINUE
      NST(NAA-1)=19
      NST(NAA)=19
      NAD(NAA-1)=NUMM(I)
      NAD(NAA)=NUMM(I)
      NBD(NAA-1)=NUMM(JJ(1))
      NBD(NAA)=NUMM(JJ(1))
      NCD(NAA-1)=NUMM(JJ(2))
      NCD(NAA)=NUMM(JJ(2))
      NUMM(NAA-1)=NAA-1
      NUMM(NAA)=NAA
      GO TO 39
C
C     THREE BONDS MISSING ON SP3 CARBON
C     OR TWO BONDS MISSING ON SP2 CARBON OR AMINE NITROGEN
C     OR ONE BOND MISSING ON OXYGEN
C
  120 NAA=NAA+1
      NUMAT=NUMAT+1
      IBOND=IBOND+1
C
C FIND ANOTHER ATOM BONDED TO THE BOND PARTNER OF I
C
      IK=JJ(1)
C
C NEAREST NEIGHBOR IS SP CARBON - NO GOOD
C
      IF (ISORT(IK).EQ.10) GOTO 132
      DO 130 L=1,NAA-1
      IF (L.EQ.I) GOTO 130
      IF (VS(IK,L).GT.0.) GOTO 131
  130 CONTINUE
C
C NO ATOMS BONDED TO IK
C
      GOTO 132
  131 VLANG=0.
      DO 121 K=1,3
      X(K,NAA)=X(K,L)+X(K,I)-2.*X(K,IK)
      VLANG=VLANG+X(K,NAA)**2
  121 CONTINUE
C
C     REDUCE LENGTH OF NEW BOND TO 1.1 ANGSTROMS
C
      VCORR=1.1/SQRT(VLANG)
      DO 122 K=1,3
      X(K,NAA)=X(K,NAA)*VCORR+X(K,I)
  122 CONTINUE
      GOTO 134
C
C USE A GUESSING ALGORITHM AS LAST RESORT :
C THE FIRST NEW BOND LIES APPROXIMATELY IN THE DIRECTION OF
C THE PRODUCT OF THE SINGLE BOND VECTOR AND AN ARBITRARY SECOND VECTOR
C
  132 VLANG=0.
      VECPRO(1)=(X(2,I)-X(2,JJ(1)))*(X(3,I)+X(3,JJ(1))+.1)
      VECPRO(1)=VECPRO(1)-(X(2,I)+X(2,JJ(1)))*(X(3,I)-X(3,JJ(1)))
      VECPRO(2)=(X(3,I)-X(3,JJ(1)))*(X(1,I)+X(1,JJ(1)))
      VECPRO(2)=VECPRO(2)-(X(3,I)+X(3,JJ(1))+.1)*(X(1,I)-X(1,JJ(1)))
      VECPRO(3)=(X(1,I)-X(1,JJ(1)))*(X(2,I)+X(2,JJ(1)))
      VECPRO(3)=VECPRO(3)-(X(1,I)+X(1,JJ(1)))*(X(2,I)-X(2,JJ(1)))
      VLANG=SQRT(VECPRO(1)**2+VECPRO(2)**2+VECPRO(3)**2)
      VCORR=1.1/VLANG
      DO 133 K=1,3
      X(K,NAA)=VECPRO(K)*VCORR+X(K,I)
  133 CONTINUE
  134 ISORT(NAA)=19
      IAT(NAA)=1
      IATD(NAA)=1
      IMOL(NAA)=IMOL(I)
      RESNAM(NAA)=RESNAM(I)
      NUMRES(NAA)=NUMRES(I)
      RESTYP(NAA)=' H   '
      NNAT(NAA)=NAA
      NAD(NAA)=NUMM(I)
      NBD(NAA)=NUMM(JJ(1))
      IF (JJ(1).GT.1 .AND. JJ(1)-1.NE.I)THEN
      NCD(NAA)=NUMM(JJ(1)-1)
      ELSE
      NCD(NAA)=NUMM(JJ(1)+1)
      ENDIF
      VS(I,NAA)=1.1
      VS(NAA,I)=1.1
      DO 123 K=1,NAA-1
      VV(K,NAA)=SQRT((X(1,K)-X(1,NAA))**2+(X(2,K)-X(2,NAA))**2
     1+(X(3,K)-X(3,NAA))**2)
      VV(NAA,K)=VV(K,NAA)
  123 CONTINUE
      NST(NAA)=19
      NUMM(NAA)=NAA
      JJ(2)=NAA
      IF (ISO.EQ.5) GOTO 39
C
C SECOND HYDROGEN
C
      NAA=NAA+1
      NUMM(NAA)=NAA
      NNAT(NAA)=NAA
      NUMAT=NUMAT+1
      IBOND=IBOND+1
C
C     CALCULATE AND MIRROR THE SUM OF THE BOND VECTORS
C
      VLANG=0.
      DO 141 K=1,3
      X(K,NAA)=-(X(K,IK)+X(K,NAA-1)-2*X(K,I))
      VLANG=VLANG+X(K,NAA)**2
  141 CONTINUE
C
C     REDUCE LENGTH OF NEW BOND TO 1.1 ANGSTROMS
C
      VCORR=1.1/SQRT(VLANG)
      DO 142 K=1,3
      X(K,NAA)=X(K,NAA)*VCORR+X(K,I)
  142 CONTINUE
      ISORT(NAA)=19
      VS(I,NAA)=1.1
      VS(NAA,I)=1.1
      NST(NAA)=19
      IAT(NAA)=1
      IATD(NAA)=1
      IMOL(NAA)=IMOL(I)
      RESNAM(NAA)=RESNAM(I)
      NUMRES(NAA)=NUMRES(I)
      RESTYP(NAA)=' H   '
      NAD(NAA)=NUMM(I)
      NBD(NAA)=NUMM(JJ(1))
      NCD(NAA)=NUMM(JJ(2))
      DO 143 K=1,NAA-1
      VV(K,NAA)=SQRT((X(1,K)-X(1,NAA))**2+(X(2,K)-X(2,NAA))**2
     1+(X(3,K)-X(3,NAA))**2)
      VV(NAA,K)=VV(K,NAA)
  143 CONTINUE
C
C COMPLETED TREATMENT OF SP2 C AND N HERE
C
      IF (I.LE.NA) GOTO 39
C
C PROTONATE ALIPHATIC AMINES (I.E. TREAT THEM AS METHYLS) IF TYPE 202
C
      INS=NST(NUMM(I))
      IF (INS.EQ.102.OR.INS.EQ.2) GOTO 39
      IF (INS.EQ.202) QADD(I)=1.
C
C THE FIRST HYDROGEN ECLIPSES ATOM IK. DELETE IT AND GENERATE TWO NEW
C HYDROGENS
C
      X(1,NAA-1)=X(1,NAA)
      X(2,NAA-1)=X(2,NAA)
      X(3,NAA-1)=X(3,NAA)
      NAA=NAA-1
      NUMAT=NUMAT-1
      GOTO 104
   39 IWARN=IWARN-1
      RETURN
      END
      SUBROUTINE GAMMA(NAA,ISORT,RV)
C
C    ELECTRONIC INTERACTION ENERGIES ARE DETERMINED BY A
C    OHNO/KLOPMAN-LIKE RELATION
C
      PARAMETER(IX=301,IY=IX,IX3=IX*3,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM3/VV(IX,IX),VS(IX,IX),EPP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /GAMMMA/ RR(IX,IX)
cs
      common /solvat/ eps,dielec,dsgl,fav(ix),fam(ix),rxy(ix,ix),
     1 rvdw(50),ion(50)
      dimension isort(ix),rv(50)
cs
      R(Q1,Q2,Q3)=14.397/SQRT(Q1**2+0.25*(14.397/Q2+14.397/Q3)**2)
      DO 10 I=1,NAA
      RR(I,I)=RV(ISORT(I)+1)
      DO 10 J=1,I-1
      RR(I,J)=R(VV(I,J),RV(ISORT(I)+1),RV(ISORT(J)+1))
   10 RR(J,I)=RR(I,J)
      if (neps.eq.0) RETURN
C
C    STRUCT:  CALCULATE STRUCTURAL FACTORS (%) OF SOLVATION SHELL
C
      CALL STRUCT(NAa)
C
C     ATOMPAAR-STRUKURFAKTOREN
C
      DO 9100 I=1,NAa
      rxyz=fam(i)**2
      RXY(I,i)=7.199/sqrt(rxyz)
      DO 9100 J=I+1,NAa
      RXYZ=FAm(I)*FAm(J)
      RXY(I,J)=14.397/sqrt(vv(i,j)**2+rxyz*exp(-vv(i,j)**2*0.5/rxyz))
      RXY(J,I)=RXY(I,J)
 9100 CONTINUE
C
      RETURN
      END
      SUBROUTINE STRUCT(NAa)
C
C    THIS SUBROUTINE IS ENTERED ONCE TO EVALUATE THE ATOMIC
C    STRUCTURAL FACTORS USED TO QUANTIFY THE EFFECTS OF SOLVATION
C
      PARAMETER(IX=301,IY=IX,IX3=IX*3,IX4=IX*4,IX6=IX*6)
      logical*2 ctest(ix)
      COMMON /GEOM1/ X(3,IX),isort(ix),a(3,ix),nra(4,3,ix),irst(4,ix),
     1irz,arst(ix),numm(ix),ab(3,ix),ifix(ix),imol(ix),iat(ix)
      COMMON /GEOM3/VV(IX,IX),VS(IX,IX),EPP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /GAMMMA/ RR(IX,IX)
      COMMON /POINT/ XPP(3,288)
cs
      common /solvat/ eps,dielec,dsgl,fav(ix),fam(ix),rxy(ix,ix),
     1 rvdw(50),ion(50)
cs
      DIMENSION XYZNP(3),QVDW(IX),
     1qborn(ix),qsvdw(ix),qsborn(ix),FFACT(50)
      do 1 i=1,50
      rvdw(i)=1.8
      ion(i)=0
   1  continue
      ion(13)=1
      ion(14)=1
      ion(18)=1
      ion(19)=1
      do 2 i=21,50
   2  ion(i)=1
      RVDW(1)=1.7
      RVDW(2)=1.9
      RVDW(3)=1.6
      RVDW(4)=1.6
      RVDW(5)=1.5
      RVDW(6)=1.6
      RVDW(7)=1.8
      RVDW(8)=1.45
      RVDW(9)=1.47
      RVDW(10)=2.
      RVDW(11)=1.7
      RVDW(12)=1.55
      rvdw(13)=1.35
      rvdw(14)=1.5
      RVDW(15)=2.
      RVDW(16)=1.9
      RVDW(17)=1.8
      rvdw(18)=0.7
      RVDW(19)=1.8
      RVDW(20)=1.2
      do 11 i=1,50
      ffact(i)=0.7
   11 continue      
      ffact(3)=0.8
      ffact(5)=0.85
      ffact(20)=0.667
      DO 12 I=1,NAa
      QVDW(I)=RVDW(ISORT(I)+1)+0.3
      qsvdw(i)=qvdw(i)**2
   12 CONTINUE
      SURF=0.
      DO 40 I=1,NAa
      FAV(i)=0.
      do 21 l=1,naa
      ctest(l)=.false.
      IF (l.EQ.i) ctest(l)=.true.
      if (vv(i,l).gt.(qvdw(i)+qvdw(l))) ctest(l)=.true.
   21 continue
      DO 30 J=1,288,6
C
C    CALCULATION FAV() WITH 288 POINTS >J< ARE PLACED AROUND ATOM (I)
C    but only every sixth point will be used to speed up the program
C
      do 31 m=1,3
   31 XYZNP(m)=X(m,I)+xpp(m,j)*QVDW(I)
      ZPA=1.
C
C    THE CONTRIBUTION OF >J< TO THE STRUCTURAL FACTOR OF (I),
C    FAV(I), IS CHARACTERIZED BY "ZP", WHICH DEPENDS ON MOLECULAR
C    GEOMETRY AND ON THE DIMENSIONS OF THE SOLVENT
C
      DO 20 K=1,NAa
      if (ctest(k)) goto 20
      DISTA=(XYZNP(1)-X(1,K))**2+(XYZNP(2)-X(2,K))**2+
     &(XYZNP(3)-X(3,K))**2
      IF (DISTA.gt.qsvdw(k)) goto 20
      zpa=0.
      goto 30
   20 CONTINUE
   30 FAV(I)=FAV(i)+ZPA
      FAV(i)=FAV(i)/288*6
      SURF=SURF+FAV(I)*4.*3.1415926*(QVDW(I))**2
   40 CONTINUE
      DSGL=-0.044*surf
c           0.046122   
c
c versuch born-radien nach still
c
      do 110 i=1,naa
      qborn(i)=rvdw(isort(i)+1)
  110 qsborn(i)=qborn(i)**2
c
      do 160 i=1,naa
      akon=0.
      delta=0.1
      rtest=qborn(i)-0.5*delta
      fam(i)=0.
      do 150 n=1,20
      favi=0.
      do 111 l=1,naa
      ctest(l)=.false.
      if (l.eq.i) ctest(l)=.true.
      if (vv(i,l).gt.(qborn(l)+rtest)) ctest(l)=.true.
      if (rtest.gt.(vv(i,l)+qborn(l))) ctest(l)=.true.
  111 continue
      do 140 j=1,288,6
      do 120 m=1,3
  120 xyznp(m)=x(m,i)+xpp(m,j)*rtest
      za=1.
      do 130 k=1,naa
      if (ctest(k)) goto 130
      dista=(xyznp(1)-x(1,k))**2+(xyznp(2)-x(2,k))**2
     &+(xyznp(3)-x(3,k))**2
      if (dista.gt.qsborn(k)) goto 130
      za=0.
      goto 140
  130 continue
  140 favi=favi+za
      if (favi.eq.48.) then
      akon=akon+1./(rtest+delta)
      fam(i)=1./akon*ffact(isort(i)+1)
      goto 160
      end if
      akon=akon+favi/48.*(1./rtest-1./(rtest+delta))
      rtest=rtest+delta
      delta=delta*1.5
  150 continue
  160 continue
c
c
c
      RETURN
      END
      SUBROUTINE SOLVut(NAA,ISORT,NUMM,cs,qseff)
c
      PARAMETER(IX=301,IY=IX*2/3,IX4=IX*4,IX6=IX*6)
      dimension isort(ix),numm(ix),cs(ix),qseff(ix)
c
cs
      common /solvat/ eps,dielec,dsgl,fav(ix),fam(ix),rxy(ix,ix),
     1 rvdw(50),ion(50)
cs
C
      WRITE(6,FMT='('' ===> CHARACTERISTICS OF SOLVATION:''/
     18X,''I SORT(I)  FAV(I)       FAM(I)        Q(I)   '')')
      DO 777 I=1,NAA
  777 WRITE(6,66) NUMM(I),ISORT(I),fav(I),FAM(I),QSeff(I)+CS(I)
   66 FORMAT(I9,I4,2X,3F12.6)
      END
      subroutine incrys
c
c     This subroutine reads crystal data produced by GSTAT with
c     OUTPUT COORDS FRAC
c     with some additional integers
c
c     (1) Title
c     (2) Options
c     as decribed above
c
c     (3) CELL as put out by GSTAT
c     (4) INSERT: ksym, number of SYMM lines
c                 iscl, crystal system indicator
c                       =1 cubic
c                       =2 tetragonal
c                       =3 orthorhombic
c                       =4 monoclinic
c                       =5 triclinic
c     (5) ksym SYMM lines as put out by GSTAT.
c     (6) atomic coordinates - change element symbol to atomic number
c         and PIMM atom sort number, format as given above in sect.(3A)
c
c     (7) two empty lines
c
c     Example:
c
c   0     benzene orthorhombic modification         
c         test
c  6  3  0  2  2.00  2  0 80  3  5 -7 20 20  0
cCELL     7.440   9.550  6.920  90.000  90.000  90.000
c   4  3  1
cSYMM      1.  0   0.  .00000   0.  1.  0.  .00000   0.  0.  1.  .00000 
cSYMM      1.  0.  0.  .50000   0. -1.  0.  .50000   0.  0. -1.  .00000  
cSYMM     -1.  0.  0.  .00000   0.  1.  0.  .50000   0.  0. -1.  .50000
cSYMM     -1.  0.  0.  .50000   0. -1.  0.  .00000   0.  0.  1.  .50000
c  6  1       -.05690    .13870   -.00540
c  6  1       -.13350    .04600    .12640
c  6  1        .07740    .09250   -.12950
c  6  1       -.07740   -.09250    .12950
c  6  1        .05690   -.13870    .00540
c  6  1        .13350   -.04600   -.12640
c  1 19       -.09760    .24470   -.01770
c  1 19       -.24090    .07940    .22180
c  1 19        .13710    .16310   -.23120
c  1 19       -.13710   -.16310    .23120
c  1 19        .24090   -.07940   -.22180
c  1 19        .09760   -.24470    .01770
c
c
c -1
c
      character*4 text
      double precision xkc,ykc,zkc,vircx,vircy,vircz,viryx,virzx,virzy
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      parameter(IX=301)
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1 cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2 imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3 ,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
      do 10 i=1,192
      do 10 j=1,3
         tsym(j,i)=0.
      do 10 k=1,3
   10 rsym(j,k,i)=0.
      do 20 i=1,3
         rsym(i,i,1)=1.
      do 20 j=1,3
      ascl(i,j)=0.
   20 cl(i,j)=0.
c
      read(iin,*) text,(ck(i),i=1,3),calp,cbet,cgam
      cal=calp/57.2957795
      cbe=cbet/57.2957795
      cga=cgam/57.2957795
      read(iin,*) ksym,iscl 
      do 30 k=1,ksym
      read(iin,*) text,(rsym(1,j,k),j=1,3),tsym(1,k),
     1(rsym(2,j,k),j=1,3),tsym(2,k),(rsym(3,j,k),j=1,3),tsym(3,k)
   30 continue
      if (iscl.le.0) goto 90
      ascl(1,1)=1.
      ascl(3,3)=1.
      ascl(2,2)=1.
      if (iscl.le.3) goto 90
      ascl(1,3)=1.
      if (iscl.eq.4) goto 90
      ascl(1,2)=1.
      ascl(2,3)=1.
   90 return
      end
      subroutine crycrt
      double precision xkc,ykc,zkc,vircx,vircy,vircz,viryx,virzx,virzy    
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      parameter(IX=301)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM2/NUMAT,IATD(IX),NFIX(IX),ILES,NNAT(IX),NAD
     1(IX),NBD(IX),NCD(IX),NDD(IX),NCTOR,MCTOR,NPERM(20),ITOR(IX),
     2NST(IX),NMOL(IX),KTOR
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1 cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2 imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3 ,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
      dimension xmax(3),xmin(3)
      do 10 i=1,3
      do 10 j=1,3
   10 cl(i,j)=0.
      do 20 i=1,naa
      do 20 j=1,3
   20 cx(j,i)=x(j,i)
      ca=cos(cal)
      sa=sin(cal)
      cb=cos(cbe)
      sb=sin(cbe)
      cg=cos(cga)
      sg=sin(cga)
      cl(1,1)=ck(1)
      cl(1,2)=ck(2)*cg
      cl(2,2)=ck(2)*sg
      cl(1,3)=ck(3)*cb
      cl(2,3)=ck(3)*(ca-cb*cg)/sg
      cl(3,3)=ck(3)*sqrt(sa**2-cb**2-cg**2+2.*ca*cb*cg)/sg
      do 30 i=1,naa
      x(1,i)=cx(1,i)*cl(1,1)+cx(2,i)*cl(1,2)+cx(3,i)*cl(1,3)
      x(2,i)=cx(2,i)*cl(2,2)+cx(3,i)*cl(2,3)
   30 x(3,i)=cx(3,i)*cl(3,3)
      x(1,naa+1)=cl(1,1)
      x(2,naa+1)=0.
      x(3,naa+1)=0.
      x(1,naa+2)=cl(1,2)
      x(2,naa+2)=cl(2,2)
      x(3,naa+2)=0.
      x(1,naa+3)=cl(1,3)
      x(2,naa+3)=cl(2,3)
      x(3,naa+3)=cl(3,3)
      cutoff=10.
      do 50 j=1,3
      xmin(j)=1000.
      xmax(j)=-1000.
      do 40 i=1,naa
      xmax(j)=amax1(xmax(j),x(j,i))
      xmin(j)=amin1(xmin(j),x(j,i))
   40 continue
      xmax(j)=xmax(j)+cutoff
      xmin(j)=xmin(j)-cutoff
   50 continue
      imix=int(xmin(1)/cl(1,1))-1
      imiy=int(xmin(2)/cl(2,2))-1
      imiz=int(xmin(3)/cl(3,3))-1
      imax=int(xmax(1)/cl(1,1))+1
      imay=int(xmax(2)/cl(2,2))+1
      imaz=int(xmax(3)/cl(3,3))+1
      vvol=float((imax-imix+1)*(imay-imiy+1)*(imaz-imiz+1))
      return
      end
      subroutine crtcry
      double precision xkc,ykc,zkc,vircx,vircy,vircz,viryx,virzx,virzy
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      parameter(IX=301)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM2/NUMAT,IATD(IX),NFIX(IX),ILES,NNAT(IX),NAD
     1(IX),NBD(IX),NCD(IX),NDD(IX),NCTOR,MCTOR,NPERM(20),ITOR(IX),
     2NST(IX),NMOL(IX),KTOR
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
      dimension abc(3)
      if (nm.lt.-5) then
      do 10 i=1,3
      abc(i)=0.
      do 10 j=1,3
   10 cl(i,j)=x(i,naa+j)
      do 30 i=1,3
      do 20 j=1,3
   20 abc(i)=abc(i)+cl(j,i)*cl(j,i)
   30 abc(i)=sqrt(abc(i))
      gam=(cl(1,1)*cl(1,2)+cl(2,1)*cl(2,2)+cl(3,1)*cl(3,2))/
     1(abc(1)*abc(2))
      bet=(cl(1,1)*cl(1,3)+cl(2,1)*cl(2,3)+cl(3,1)*cl(3,3))/
     1(abc(1)*abc(3))
      alp=(cl(1,2)*cl(1,3)+cl(2,2)*cl(2,3)+cl(3,2)*cl(3,3))/
     1(abc(2)*abc(3))
      cal=acos(alp)
      cbe=acos(bet)
      cga=acos(gam)
      do 40 i=1,3
   40 ck(i)=abc(i)
      end if
      rf1=1./cl(1,1)
      rf2=-cl(1,2)/(cl(1,1)*cl(2,2))
      rf3=1./cl(2,2)
      rf4=-(cl(2,2)*cl(1,3)-cl(1,2)*cl(2,3))/(cl(1,1)*cl(2,2)*cl(3,3))
      rf5=-cl(2,3)/(cl(2,2)*cl(3,3))
      rf6=1./cl(3,3)
      do 60 i=1,naa
      cx(1,i)=rf1*x(1,i)+rf2*x(2,i)+rf4*x(3,i)
      cx(2,i)=rf3*x(2,i)+rf5*x(3,i)
   60 cx(3,i)=rf6*x(3,i)
      return
      end
      subroutine crypac            
      DOUBLE PRECISION XK,YK,ZK,xkc,ykc,zkc,virx,viry,virz,
     1vircx,vircy,vircz,viryx,virzx,virzy
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      PARAMETER(IX=301,IY=IX,IX3=IX*3,IX4=IX*4,IX6=IX*6,
     1IXA=(IX3*(IX3+1))/2)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /OPTIM/virx,viry,virz,XD(IX3),XDD(IX3),QZ(IX3),DX(IX3),
     1XK(IX),YK(IX),ZK(IX),NAM,EPP,EPW
      COMMON /OPTIM1/ AMAT(IXA)
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
      common /pacplo/ npac,spac(2000),cpac(3,2000),ihhyd(100),
     1 iahyd(100),dhyd(100),nhyd,jhyd(ix),xo(3,ix),xsum(3,ix),
     2 ibhyd(100),hhyd(100),abhyd(100),anghyd(100)   
      dimension tr(3,16),ace(3),trans(3),xi(3)
      nhyd=0
      epce=0.
      epcv=0.
      do 1 i=1,naa
    1 qseff(i)=qsig(i)
      if (ihb.eq.0) call cryhyd
      do 10 i=1,ksym
      do 10 j=1,3
      tr(j,i)=0.
      do 10 k=1,3
   10 tr(j,i)=tr(j,i)+cl(j,k)*tsym(k,i)
      do 20 i=1,3
      xkc(naa+i)=0.
      ykc(naa+i)=0.
      zkc(naa+i)=0.
   20 continue
      vircx=0.
      vircy=0.
      vircz=0.
      viryx=0.
      virzx=0.
      virzy=0.
      do 120 i=1,naa
      is1=isort(i)+1
      xkc(i)=0.
      ykc(i)=0.
      zkc(i)=0.
      xpac(i)=0.
      do 110 ia=imix,imax
      do 110 ib=imiy,imay
      do 110 ic=imiz,imaz
      ace(1)=ia
      ace(2)=ib
      ace(3)=ic
      do 30 n=1,3
      trans(n)=0.
      do 30 m=1,3
   30 trans(n)=trans(n)+cl(n,m)*ace(m)
      do 100 isym=1,ksym
      if (isym.eq.1.and.ia.eq.0.and.ib.eq.0.and.ic.eq.0) goto 100
      do 90 j=1,naa
      do 40 m=1,3
   40 xi(m)=0.
      do 60 m=1,3
      do 50 n=1,3
      xi(m)=xo(n,j)*rsym(m,n,isym)+xi(m)
   50 continue
      xi(m)=xi(m)+tr(m,isym)+trans(m)
   60 continue
      dxs=0.
      do 70 m=1,3
   70 dxs=dxs+(x(m,i)-xi(m))**2
      dxr=sqrt(dxs)
      is2=isort(j)+1
      dxrr=dxr
      ph=0.
      rrr=14.397/sqrt(dxs+0.25*(14.397/RV(IS1)+14.397/RV(IS2))**2)
      ccc=(qseff(i)+cs(i))*(qseff(j)+cs(j))
      ph=ccc*(-7.73e-7)*rrr**3*1.5
      epce=epce+ccc*rrr*0.5
      xpac(i)=xpac(i)+0.5*(qsig(j)+cs(j))*rrr
ctest
      if ((jhyd(i)+jhyd(j)).eq.1) goto 80
      IXP=IVDW(IS1,IS2)
      DXP=IXP
      EXX=-BVDW(IS1,IS2)*dxrr
      AEXPB=AVDW(IS1,IS2)*EXP(EXX)
      CEXP=CVDW(IS1,IS2)*dxrr**(-6)*0.9
      VEXP=dxrr**(-IXP)
      epcvv=0.02167*(AEXPB*VEXP-CEXP)
      epcv=epcv+epcvv
      PH=ph+0.6944E-5/dxr*(-AEXPB*BVDW(IS1,IS2)*VEXP+(-AEXPB*DXP
     1  *VEXP+6.*CEXP)/dxrr)
  80  phx=PH*(X(1,I)-xi(1))
      phy=PH*(X(2,I)-xi(2))
      phz=PH*(X(3,I)-xi(3))
      XKc(I)=XKc(I)-phx
      YKc(I)=YKc(I)-phy
      ZKc(I)=ZKc(I)-phz
      if (nm.gt.-6) goto 90
      phx=PH*(X(1,I)-xi(1))
      phy=PH*(X(2,I)-xi(2))
      phz=PH*(X(3,I)-xi(3))
      xkc(naa+1)=xkc(naa+1)+phx*(tsym(1,isym)+ace(1))
      xkc(naa+2)=xkc(naa+2)+phx*(tsym(2,isym)+ace(2))
      ykc(naa+2)=ykc(naa+2)+phy*(tsym(2,isym)+ace(2))
      xkc(naa+3)=xkc(naa+3)+phx*(tsym(3,isym)+ace(3))
      ykc(naa+3)=ykc(naa+3)+phy*(tsym(3,isym)+ace(3))
      zkc(naa+3)=zkc(naa+3)+phz*(tsym(3,isym)+ace(3))
   90 continue
  100 continue
  110 continue
  120 continue
 9999 format (3e15.6)
      vircx=xkc(naa+1)*x(1,naa+1)
      vircy=ykc(naa+2)*x(2,naa+2)
      vircz=zkc(naa+3)*x(3,naa+3)
      viryx=xkc(naa+2)*x(2,naa+2)
      virzx=xkc(naa+3)*x(3,naa+3)
      virzy=ykc(naa+3)*x(3,naa+3)
      epw=epw+epce+epcv
      return 
      end
      subroutine cryhyd            
      DOUBLE PRECISION XK,YK,ZK,xkc,ykc,zkc,virx,viry,virz,
     1vircx,vircy,vircz,viryx,virzx,virzy
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      PARAMETER(IX=301,IY=IX,IX3=IX*3,IX4=IX*4,IX6=IX*6,
     1IXA=(IX3*(IX3+1))/2)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /OPTIM/virx,viry,virz,XD(IX3),XDD(IX3),QZ(IX3),DX(IX3),
     1XK(IX),YK(IX),ZK(IX),NAM,EPP,EPW
      COMMON /OPTIM1/ AMAT(IXA)
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
      common /pacplo/ npac,spac(2000),cpac(3,2000),ihhyd(100),
     1 iahyd(100),dhyd(100),nhyd,jhyd(ix),xo(3,ix),xsum(3,ix),
     2 ibhyd(100),hhyd(100),abhyd(100),anghyd(100)   
      dimension tr(3,16),ace(3),trans(3),xi(3)
      nhyd=0
      epce=0.
      epcv=0.
      do 10 i=1,ksym
      do 10 j=1,3
      tr(j,i)=0.
      do 10 k=1,3
   10 tr(j,i)=tr(j,i)+cl(j,k)*tsym(k,i)
      do 120 i=1,naa
      is1=isort(i)+1
      if (is1.ne.20) goto 120
      do 11 ii=1,naa
      if (vs(i,ii).gt.0.) then
      if (isort(ii).lt.2.or.isort(ii).gt.9) goto 120
      nhy=ii
      goto 111
      endif
   11 continue
  111 xkc(i)=0.
      ykc(i)=0.
      zkc(i)=0.
      xpac(i)=0.
      do 110 ia=imix,imax
      do 110 ib=imiy,imay
      do 110 ic=imiz,imaz
      ace(1)=ia
      ace(2)=ib
      ace(3)=ic
      do 30 n=1,3
      trans(n)=0.
      do 30 m=1,3
   30 trans(n)=trans(n)+cl(n,m)*ace(m)
      do 100 isym=1,ksym
      if (isym.eq.1.and.ia.eq.0.and.ib.eq.0.and.ic.eq.0) goto 100
      do 90 j=1,naa
      do 40 m=1,3
   40 xi(m)=0.
      do 60 m=1,3
      do 50 n=1,3
      xi(m)=xo(n,j)*rsym(m,n,isym)+xi(m)
   50 continue
   60 xi(m)=xi(m)+tr(m,isym)+trans(m)
      dxs=0.
      do 70 m=1,3
   70 dxs=dxs+(x(m,i)-xi(m))**2
      dxr=sqrt(dxs)
      is2=isort(j)+1
      if (dxr.gt.2.8) goto 75
      if (is2.lt.3.or.is2.ge.11) goto 75
      if (is2.eq.3.and.cs(j).ne.0.) goto 75
      iadd=0
      if (qadd(nhy).ne.0.) iadd=1
      dys=0.
      do 701 m=1,3
  701 dys=dys+(x(m,nhy)-xi(m))**2
      dyr=sqrt(dys)
      dzs=vv(i,nhy)**2
      dzr=vv(i,nhy)
      AKR=(dxs+dzs-dys)/(2*dxr*dzr)
      IF (AKR.GT.-0.5) GOTO 75
      ZQ=2.8/dxr*(-AKR)
      nhyd=nhyd+1
      ihhyd(nhyd)=i
      iahyd(nhyd)=j
      ibhyd(nhyd)=nhy
      abhyd(nhyd)=dyr
      anghyd(nhyd)=acosf(akr)*57.296
      hhyd(nhyd)=dzr
      dhyd(nhyd)=dxr
      jhyd(i)=1
      qseff(i)=qseff(i)+0.019*zq
      qseff(j)=qseff(j)-0.019*zq
      if (vs(i,j).eq.0.) vs(i,j)=-1.
   75 continue
   90 continue
  100 continue
  110 continue
  120 continue
      return 
      end
      subroutine outcry
      CHARACTER*2 ATA,asym,leer,xplus,xminus,yplus,yminus,zplus,zminus,
     1ato 
      double precision xkc,ykc,zkc,vircx,vircy,vircz,viryx,virzx,virzy
      COMMON /LABEL/ ATA(50)
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      PARAMETER(IX=301,IY=IX,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /STRAIN/ SPANB(IX),SPANA(IX),SPANT(IX),SPANK(IX),SPANE(IX),
     1SPANV(IX)
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1 cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2 imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3 ,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
      common /pacplo/ npac,ispac(2000),cpac(3,2000),ihhyd(100),
     1 iahyd(100),dhyd(100),nhyd,jhyd(ix),xo(3,ix),xsum(3,ix),
     2 ibhyd(100),hhyd(100),abhyd(100),anghyd(100)   
      dimension ifac(50),asym(3,3,16),ato(50)
      data ifac /1,1,3,3,4,4,5,5,5,5,1,3,5,5,5,5,5,5,5,2,30*5/
      data ato /'C','C','N','N','O','O','S','F','F','F','C','N',
     17*'X','H',30*'X'/
      leer='  '
      xplus='+X'
      yplus='+Y'
      zplus='+Z'
      xminus='-X'
      yminus='-Y'
      zminus='-Z'
      calp=cal*57.2957795
      cbet=cbe*57.2957795
      cgam=cga*57.2957795
      write(iarc,1100) dhfm,(title(i),i=1,8)
 1100 format (5x,f8.1,2x,8a4)     
      write(iarc,100)(ck(i),i=1,3),calp,cbet,cgam
      write(iarc,110)
      do 10 i=1,ksym
      write(iarc,120)(rsym(1,j,i),j=1,3),tsym(1,i),
     1(rsym(2,j,i),j=1,3),tsym(2,i),(rsym(3,j,i),j=1,3),tsym(3,i)
   10 continue
      write(iarc,130)
      do 20 i=1,naa
      is=isort(i)+1
      write(iarc,140)ata(is),numm(i),(cx(n,i),n=1,3)
   20 continue
      if (nhyd.gt.0.and.numm(1).gt.0) then
      write(iarc,160) 
      do 30 i=1,nhyd
      is=isort(iahyd(i))+1
      it=isort(ihhyd(i))+1
      iv=isort(ibhyd(i))+1
      write(iarc,170)ata(is),numm(iahyd(i)),ata(it),numm(ihhyd(i)),
     1ata(iv),numm(ibhyd(i)),dhyd(i),hhyd(i),abhyd(i),anghyd(i)
   30 continue
      end if
      epcea=epce*96.533
      epcva=epcv*96.533
      epaca=epcea+epcva
      write(iarc,150)epaca,epcea,epcva
      if (numm(1).gt.0) then
      write(iesp,1000) dhfm,(title(i),i=1,8)
 1000 format ('TITL',f8.1,2x,8a4)     
      write(iesp,1001)(ck(i),i=1,3),calp,cbet,cgam
 1001 format ('CELL  1.5418 ',6f10.4)
      write(iesp,1002)
 1002 format ('LATT  -1')
      if (ksym.gt.1) then
      do 1005 i=2,ksym
      if (rsym(1,1,i).eq.0.) asym(1,1,i)=leer
      if (rsym(1,1,i).eq.1.) asym(1,1,i)=xplus
      if (rsym(1,1,i).eq.-1.) asym(1,1,i)=xminus
      if (rsym(2,1,i).eq.0.) asym(2,1,i)=leer
      if (rsym(2,1,i).eq.1.) asym(2,1,i)=xplus
      if (rsym(2,1,i).eq.-1.) asym(2,1,i)=xminus
      if (rsym(3,1,i).eq.0.) asym(3,1,i)=leer
      if (rsym(3,1,i).eq.1.) asym(3,1,i)=xplus
      if (rsym(3,1,i).eq.-1.) asym(3,1,i)=xminus
      if (rsym(1,2,i).eq.0.)  asym(1,2,i)=leer
      if (rsym(1,2,i).eq.1.)  asym(1,2,i)=yplus
      if (rsym(1,2,i).eq.-1.) asym(1,2,i)=yminus
      if (rsym(2,2,i).eq.0.)  asym(2,2,i)=leer
      if (rsym(2,2,i).eq.1.)  asym(2,2,i)=yplus
      if (rsym(2,2,i).eq.-1.) asym(2,2,i)=yminus
      if (rsym(3,2,i).eq.0.)  asym(3,2,i)=leer
      if (rsym(3,2,i).eq.1.)  asym(3,2,i)=yplus
      if (rsym(3,2,i).eq.-1.) asym(3,2,i)=yminus
      if (rsym(1,3,i).eq.0.)  asym(1,3,i)=leer
      if (rsym(1,3,i).eq.1.)  asym(1,3,i)=zplus
      if (rsym(1,3,i).eq.-1.) asym(1,3,i)=zminus
      if (rsym(2,3,i).eq.0.)  asym(2,3,i)=leer
      if (rsym(2,3,i).eq.1.)  asym(2,3,i)=zplus
      if (rsym(2,3,i).eq.-1.) asym(2,3,i)=zminus
      if (rsym(3,3,i).eq.0.)  asym(3,3,i)=leer
      if (rsym(3,3,i).eq.1.)  asym(3,3,i)=zplus
      if (rsym(3,3,i).eq.-1.) asym(3,3,i)=zminus
      write(iesp,1003)tsym(1,i),(asym(1,j,i),j=1,3),
     1tsym(2,i),(asym(2,j,i),j=1,3),tsym(3,i),(asym(3,j,i),j=1,3)
 1003 format('SYMM',f10.3,3a2,',',f10.3,3a2,',',f10.3,3a2)
 1005 continue
      end if
      write(iesp,1009)
 1009 format ('SFAC  C H N O X')  
      do 1004 i=1,naa
      is=isort(i)+1
      write(iesp,1006)ato(is),numm(i),ifac(is),(cx(n,i),n=1,3)
 1006 format(a1,i3,i10,3f10.5,'   11.0000   0.05') 
 1004 continue
      write(iesp,1007)
 1007 format('HKLF 4'/'END') 
      end if
  100 format(10x,'lattice constants',/6f10.4)
  110 format(10x,'crystal symmetry')
  120 format(3(3f5.0,f10.3))
  130 format(10x,'fractional atomic coordinates')
  140 format(10x,a2,i3,10x,3f10.4)
  150 format(10x,'lattice energy , -electrostatic, -van-der-Waals ',/,
     110x,3f10.2,'  Kj/mol')
  160 format(10x,'Intermolecular hydrogen bonds'/
     1' atom(a) .. H(b) -  atom(c)  r(a,b)   r(b,c)   r(a,c) a(a,b,c)')
  170 format(2x,a2,i3,5x,a2,i3,5x,a2,i3,4f9.3)
      return
      end
      subroutine outpac            
      CHARACTER*2 ATA
      double precision xkc,ykc,zkc,vircx,vircy,vircz,viryx,virzx,virzy
      COMMON /LABEL/ ATA(50)
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      PARAMETER(IX=301,IY=IX,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /STRAIN/ SPANB(IX),SPANA(IX),SPANT(IX),SPANK(IX),SPANE(IX),
     1SPANV(IX)
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1 cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2 imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3 ,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
      common /pacplo/ npac,ispac(2000),cpac(3,2000),ihhyd(100),
     1 iahyd(100),dhyd(100),nhyd,jhyd(ix),xo(3,ix),xsum(3,ix),
     2 ibhyd(100),hhyd(100),abhyd(100),anghyd(100)   
      dimension tr(3,16),ice(3),trans(3),xi(3),xs(3)
      ana=naa
      do 5 i=1,3
    5 xs(i)=0.
      do 6 i=1,naa
      do 6 j=1,3
    6 xs(j)=xs(j)+x(j,i)
      do 7 i=1,3
    7 xs(i)=xs(i)/ana
      do 10 i=1,ksym
      do 10 j=1,3
      tr(j,i)=0.
      do 10 k=1,3
   10 tr(j,i)=tr(j,i)+cl(j,k)*tsym(k,i)
      npac=0
      do 110 ia=imix,imax
      do 110 ib=imiy,imay
      do 110 ic=imiz,imaz        
      ice(1)=ia
      ice(2)=ib
      ice(3)=ic
      do 30 n=1,3
      trans(n)=0.
      do 30 m=1,3
   30 trans(n)=trans(n)+cl(n,m)*ice(m)
      do 100 isym=1,ksym
      if (isym.eq.1.and.ia.eq.0.and.ib.eq.0.and.ic.eq.0) goto 100
      do 90 j=1,naa
      do 40 m=1,3
   40 xi(m)=0.
      do 60 m=1,3
      do 50 n=1,3
   50 xi(m)=xo(n,j)*rsym(m,n,isym)+xi(m)
   60 xi(m)=xi(m)+tr(m,isym)+trans(m)
      dxs=0.
      do 70 m=1,3
   70 dxs=dxs+(xi(m)-xs(m))**2
      dxr=sqrt(dxs)
      if (dxr.gt.12.) goto 90
      is=isort(j)+1
      if (is.ne.20) then
      npac=npac+1
      ispac(npac)=is
      do 80 nn=1,3
   80 cpac(nn,npac)=xi(nn)
      end if
   90 continue
  100 continue
  110 continue
      write(iarc,120) npac
  120 format (10x,'Number of non hydrogen atoms in the',
     1' plot file',i5)
      return
      end
      subroutine freqan
      DOUBLE PRECISION XK,YK,ZK,xkc,ykc,zkc,virx,viry,virz,
     1vircx,vircy,vircz,viryx,virzx,virzy
      CHARACTER*4 RESNA,RESNAM
      CHARACTER*5 RESTY,RESTYP
      PARAMETER(IX=301,IY=IX,ix3=ix*3,IX4=IX*4,IX6=IX*6)
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM2/NUMAT,IATD(IX),NFIX(IX),ILES,NNAT(IX),NAD
     1(IX),NBD(IX),NCD(IX),NDD(IX),NCTOR,MCTOR,NPERM(20),ITOR(IX),
     2NST(IX),NMOL(IX),KTOR
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /BETA/ H(IY,IY)
      COMMON /DYN/ TIME,TEMP,GMASS,MDCYC,TEMITT,EMITT,TMASS,MDUMP,MCOOL,
     1SCHW(3),VELO(3,IX),ACC(3,IX),X3(3,IX),TSOLL,IRAND,TE,iav
     2,pnull,taupe,taupi,pmitt,taute,pstart,pdiff,pend,mpress
      COMMON /GAMMMA/ RR(IX,IX)
      COMMON /STRAIN/ SPANB(IX),SPANA(IX),SPANT(IX),SPANK(IX),SPANE(IX),
     1SPANV(IX)
      COMMON /OPTIM/virx,viry,virz,XD(IX3),XDD(IX3),QZ(IX3),DX(IX3),
     1XK(IX),YK(IX),ZK(IX),NAM,EPP,EPW
      COMMON /PDBNAM/ RESNA(IX),NURES(IX),RESNAM(IX),NUMRES(IX),
     1RESTY(IX),RESTYP(IX)
cp
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1 cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2 imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3 ,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
      common /pacplo/ npac,ispac(2000),cpac(3,2000),ihhyd(100),
     1 iahyd(100),dhyd(100),nhyd,jhyd(ix),xo(3,ix),xsum(3,ix),
     2 ibhyd(100),hhyd(100),abhyd(100),anghyd(100)   
cp
cf
      COMMON /DIAG/ hess(IY,IY),xkl(IY),V(IY,IY),ISENT,xkr(IY)
cf
      dimension w(ix)
      nas=naa
      ka=naa*3
      del=0.025
      do 10 i=1,naa
      do 10 j=1,3
      x(j,i)=x(j,i)-del
      ii=3*(i-1)+j
      call forgeo
      do 20 k=1,ka,3
      nk=(k+2)/3
      kk=k+1
      kkk=k+2
      xkr(k)=xk(nk)
      xkr(kk)=yk(nk)
   20 xkr(kkk)=zk(nk)
      x(j,i)=x(j,i)+2.*del
      call forgeo
      do 30 k=1,ka,3
      nk=(k+2)/3
      kk=k+1
      kkk=k+2
      xkl(k)=xk(nk)
      xkl(kk)=yk(nk)
   30 xkl(kkk)=zk(nk)
      x(j,i)=x(j,i)-del
      do 40 k=1,ka
   40 hess(ii,k)=(xkr(k)-xkl(k))*500./del
   10 continue
      isent=1
      na=ka
      call thoqr
      write(iut,601)
      write(iut,602) (i,xkl(i),i=1,ka)
      do 60 i=1,ka,3
      ik=(i+2)/3
      ii=i+1
      iii=i+2
      iw=isort(ik)+1
      w(i)=amass(iw)
      w(ii)=amass(iw)
   60 w(iii)=amass(iw)
      do 70 i=1,ka
   70 w(i)=1./sqrt(w(i))
      do 80 i=1,ka
      do 80 j=1,ka
   80 hess(i,j)=w(i)*hess(i,j)*w(j)
      isent=1
      call thoqr
      do 90 i=1,ka
      xkl(i)=sign(1.,xkl(i))*sqrt(abs(1.2888607e6*xkl(i)))
      do 90 j=1,ka
   90 hess(i,j)=v(i,j)*w(j)
      do 100 i=1,ka
      dn=0.
      do 110 j=1,ka
  110 dn=dn+hess(i,j)**2
      dn=sqrt(dn)
      do 100 j=1,ka
  100 hess(i,j)=hess(i,j)/dn
      write (iut,603)
      write (iut,604) (i,xkl(i),i=1,ka)
      write (iut,605)
      naa=nas
      return
  601 format (/,10x,'eigenvalues of the force constant matrix',/)
  602 format (5(i3,f10.5,3x))
  603 format (/,10x,'eigenfrequencies in cm**(-1)',/)
  604 format (5(i3,f10.1,3x))
  605 format (/,10x,'normal modes',/)
  606 format (5(2i3,f10.5,3x))
  607 format (/)
      end
      SUBROUTINE forgeo
C
      DOUBLE PRECISION XK,YK,ZK,xkc,ykc,zkc,virx,viry,virz,
     1vircx,vircy,vircz,viryx,virzx,virzy
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      PARAMETER(IX=301,IY=IX,IX3=IX*3,IX4=IX*4,IX6=IX*6,
     1IXA=(IX3*(IX3+1))/2)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /OPTIM/virx,viry,virz,XD(IX3),XDD(IX3),QZ(IX3),DX(IX3),
     1XK(IX),YK(IX),ZK(IX),NAM,EPP,EPW
      COMMON /OPTIM1/ AMAT(IXA)
cp
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1 cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2 imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3 ,vvol,xops(3),vircx,vircy,vircz,viryx,virzx,virzy
      common /pacplo/ npac,ispac(2000),cpac(3,2000),ihhyd(100),
     1 iahyd(100),dhyd(100),nhyd,jhyd(ix),xo(3,ix),xsum(3,ix),
     2 ibhyd(100),hhyd(100),abhyd(100),anghyd(100)   
cp
C
C    SET SOME FLAGS AND CLEAR ARRAYS
C
      naax=naa
      if (nm.lt.-6) naax=naa+3
      NAM=NAAx-1
      NXX=NN-10
      SLOPE=0.
      DO 10 I=1,NAAx
         XK(I)=0.D0
         YK(I)=0.D0
   10 ZK(I)=0.D0
      IBLO=0
      SGA=0.
      EPP=0.
      EPW=0.
      EPA=0.
      virx=0.
      viry=0.
      virz=0.
      KA=3*NAAx
      do 30 i=1,naax
      do 30 j=1,naax
      if (vs(i,j).gt.0) goto 30
      vs(i,j)=0.
   30 continue
C
C    BOND LENGTH SECTION
C
      CALL BOND
C
C    BOND ANGLE SECTION
C
      CALL ANGLE
C
C    TORSIONAL ANGLE SECTION
C
      CALL TORS
      CALL BEND
C
C    NONBONDING INTERACTION SECTION
cp
      if (nm.lt.-5)then
      call crypac
      do 99 i=1,naa
      xk(i)=xk(i)+xkc(i)
      yk(i)=yk(i)+ykc(i)
      zk(i)=zk(i)+zkc(i)
   99 continue
      if (nm.lt.-6)then
      if (iscl.eq.1) then
      xkcm=(xkc(naa+1)+ykc(naa+2)+zkc(naa+3))/3.
      xkc(naa+1)=xkcm
      ykc(naa+2)=xkcm
      zkc(naa+3)=xkcm
      end if
      if (iscl.eq.2) then
      xkcm=(xkc(naa+1)+ykc(naa+2))/2.
      xkc(naa+1)=xkcm
      ykc(naa+2)=xkcm
      end if
      xk(naa+1)=xkc(naa+1)*ascl(1,1)
      yk(naa+1)=ykc(naa+1)*ascl(2,1)
      zk(naa+1)=zkc(naa+1)*ascl(3,1)
      xk(naa+2)=xkc(naa+2)*ascl(1,2)
      yk(naa+2)=ykc(naa+2)*ascl(2,2)
      zk(naa+2)=zkc(naa+2)*ascl(3,2)
      xk(naa+3)=xkc(naa+3)*ascl(1,3)
      yk(naa+3)=ykc(naa+3)*ascl(2,3)
      zk(naa+3)=zkc(naa+3)*ascl(3,3)
      end if
      end if
cp
C
      CALL NONBON
C
      return
      end
      SUBROUTINE MDINIT
      DOUBLE PRECISION xkc,ykc,zkc,vircx,vircy,vircz,viryx,virzx,virzy
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      PARAMETER(IX=301)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /DYN/ TIME,TEMP,GMASS,MDCYC,TEMITT,EMITT,TMASS,MDUMP,MCOOL,
     1SCHW(3),VELO(3,IX),ACC(3,IX),X3(3,IX),TSOLL,IRAND,TE,iav
     2,pnull,taupe,taupi,pmitt,taute,pstart,pdiff,pend,mpress
cp
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1 cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2 imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3 ,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
      common /pacplo/ npac,ispac(2000),cpac(3,2000),ihhyd(100),
     1 iahyd(100),dhyd(100),nhyd,jhyd(ix),xo(3,ix),xsum(3,ix),
     2 ibhyd(100),hhyd(100),abhyd(100),anghyd(100)   
cp
      dimension angi(3,3)
C
      EKIN1=0.
      MDCYC=0
      TIME=0.
      IF(INFORM.EQ.0) THEN
      READ (IIN,100) TSOLL,TIME,MDCYC,MDUMP,MCOOL,IRAND
      read (iin,101) iav,taute,taupe,pstart,pdiff,pend,mpress,coufac
      ELSE
      READ (IIN,*) TSOLL,TIME,MDCYC,MDUMP,MCOOL,IRAND
      read (iin,*) iav,taute,taupe,pstart,pdiff,pend,mpress,coufac
      END IF
      IF (TIME.LE.0.) TIME=0.01
      IF (MDCYC.LE.0) MDCYC=1000
      IF (MDUMP.LE.0) MDUMP=1000
      if (iav.eq.0) iav=1
      temp=tsoll
      WRITE(IUT,150) TSOLL,TIME,MDCYC,MDUMP,iav
      if (taupe.eq.0.) taupe=1000.*time
      if (taute.eq.0.) taute=time
      if (pstart.eq.0.) pstart=1.035
      if (pend.eq.0.) pend=1.035
      if (coufac.eq.0.) coufac=1.
      write(iut,203) pstart,pend,pdiff,mpress,coufac
      if (mpress.eq.0) mpress=9999999
      pnull=pstart
      taupi=taupe
      write (iut,201) taute
      write (iut,202) taupe
      IF (MCOOL.GT.0) THEN
      WRITE (IUT,160) MCOOL
      ELSE
      MCOOL=99999999
      END IF
C
      GMASS=0.
      SCHW(1)=0.
      SCHW(2)=0.
      SCHW(3)=0.
      DO 10 I=1,NAA
      DO 11 J=1,3
      schw(j)=schw(j)+x(j,i)
   11 CONTINUE
   10 CONTINUE
      SCHW(1)=SCHW(1)/naa
      SCHW(2)=SCHW(2)/naa
      SCHW(3)=SCHW(3)/naa
      TEMITT=0.
      EMITT=0.
      pmitt=0.
      DO 20 I=1,NAA
      IF (IFIX(I).GT.0) GOTO 20
      ISO1=ISORT(I)+1
      IF (TSOLL.GE.0.) THEN
      VEL=0.0015373*SQRT(TEMP/AMASS(ISO1))
      XIS=RANDOM()
      YIS=RANDOM()
      ZIS=RANDOM()
      ELSE
      VEL=0.0015373*SQRT(TEMP*NAA/GMASS)
      XIS=SCHW(1)-X(1,I)
      YIS=SCHW(2)-X(2,I)
      ZIS=SCHW(3)-X(3,I)
      END IF
      RIS=SQRT(XIS*XIS+YIS*YIS+ZIS*ZIS)
      VELO(1,I)=VEL*XIS/RIS
      VELO(2,I)=VEL*YIS/RIS
      VELO(3,I)=VEL*ZIS/RIS
      ACC(1,I)=0.
      ACC(2,I)=0.
      ACC(3,I)=0.
      X3(1,I)=0.
      X3(2,I)=0.
      X3(3,I)=0.
      VELX=VELO(1,I)**2
      VELY=VELO(2,I)**2
      VELZ=VELO(3,I)**2
      EKIN1=EKIN1+AMASS(ISO1)*(VELX+VELY+VELZ)*51.82
   20 CONTINUE
      TE=EKIN1/NAA*8165.5
      WRITE(IUT,4711) TE
c
c
c
c      call mdrotr
c
 4711 FORMAT(//,' INITIAL TEMPERATURE IS ',F8.3,/)
      WRITE(IUT,200)
      IF (TSOLL.LT.0.) TSOLL=-TSOLL
      if (nm.gt.-8) coufac=1.
      RETURN
  100 FORMAT(F6.0,F6.3,5I6)
  101 format (i6,f6.1,e6.0,3f6.1,i6,f6.1)
  150 FORMAT (//,
     1' MD SIMULATION',/,'  Temperature ',F6.0,' K, time step '
     2,F6.3,' femtoseconds,'/,I6,' cycles, writing coords every ',I6,
     3' cycles',//' in crystal dynamics: coords averaged after'
     4,I6,' step(s)' )
  160 FORMAT(/'The system will be cooled to 0 K after',I8,'
     1steps')
  200 FORMAT (//1H ,'  TIME    TEMP    PRESS   TOTAL.E   POT.E '
     1,'  KIN.E    MEAN.T    MEAN.E   MEAN.P',/,
     2             '   (FS)    (K)     (bar)   (KJ/mol) (KJ/mol)',
     3 '(KJ/mol)    (K)     (KJ/mol)   (bar) ')
  201 format(/' Temperature coupling parameter tau(t)            '
     1,f6.1)
  202 format(/' Pressure coupling parameter tau(p)               '
     1,f10.1)
  203 format(/' Pressure will be increased from',f7.3,' bar to',f7.3,
     1' bar'/' in steps of',f7.4,' bar after',i5,' steps.'//
     2' Electrostatic packing forces are increased for the factor',
     3f6.1,'.')
      END
      SUBROUTINE MDrotr
      DOUBLE PRECISION xkc,ykc,zkc,vircx,vircy,vircz,viryx,virzx,virzy
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      PARAMETER(IX=301)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /DYN/ TIME,TEMP,GMASS,MDCYC,TEMITT,EMITT,TMASS,MDUMP,MCOOL,
     1SCHW(3),VELO(3,IX),ACC(3,IX),X3(3,IX),TSOLL,IRAND,TE,iav
     2,pnull,taupe,taupi,pmitt,taute,pstart,pdiff,pend,mpress
cp
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1 cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2 imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3 ,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
      common /pacplo/ npac,ispac(2000),cpac(3,2000),ihhyd(100),
     1 iahyd(100),dhyd(100),nhyd,jhyd(ix),xo(3,ix),xsum(3,ix),
     2 ibhyd(100),hhyd(100),abhyd(100),anghyd(100)   
cp
      dimension angi(3,3)
C
c      entferne Transation und Rotation
c
      naay=naa
      xdif=0.
      ydif=0.
      zdif=0.
      velox=0.
      veloy=0.
      veloz=0.
      DO 121 I=1,NAAy
      iso=isort(i)+1
      xdif=xdif+x(1,i)
      ydif=ydif+x(2,i)
      zdif=zdif+x(3,i)
      velox=velox+velo(1,i)
      veloy=veloy+velo(2,i)
      veloz=veloz+velo(3,i)
  121 CONTINUE
      xdif=xdif/naay
      ydif=ydif/naay
      zdif=zdif/naay
      velox=velox/naay
      veloy=veloy/naay
      veloz=veloz/naay
      do 125 i=1,naay
      velo(1,i)=velo(1,i)-velox
      velo(2,i)=velo(2,i)-veloy
      velo(3,i)=velo(3,i)-veloz
      x(1,i)=x(1,i)-xdif
      x(2,i)=x(2,i)-ydif
      x(3,i)=x(3,i)-zdif
  125 continue
c      if (nm.le.-6) goto 999
      angmx=0.
      angmy=0.
      angmz=0.
      do 126 m=1,3
      do 126 n=1,3
      angi(m,n)=0.
  126 continue
      do 130 i=1,naay
      dr=(x(1,i)**2+x(2,i)**2+x(3,i)**2)
      angmx=angmx+(velo(3,i)*x(2,i)-x(3,i)*velo(2,i))/dr
      angmy=angmy+(velo(1,i)*x(3,i)-x(1,i)*velo(3,i))/dr
      angmz=angmz+(velo(2,i)*x(1,i)-x(2,i)*velo(1,i))/dr
      angi(1,1)=angi(1,1)+(x(2,i)**2+x(3,i)**2)/dr
      angi(2,2)=angi(2,2)+(x(1,i)**2+x(3,i)**2)/dr
      angi(3,3)=angi(3,3)+(x(1,i)**2+x(2,i)**2)/dr
      angi(1,2)=angi(1,2)-x(1,i)*x(2,i)/dr
      angi(1,3)=angi(1,3)-x(1,i)*x(3,i)/dr
      angi(2,3)=angi(2,3)-x(2,i)*x(3,i)/dr      
 130  continue
      angi(2,1)=angi(1,2)
      angi(3,1)=angi(1,3)
      angi(3,2)=angi(2,3)
      do 132 m=1,3
      do 132 n=1,3
      angi(m,n)=angi(m,n)/naa
 132  continue
      call invmat(angi,3)
      angmx=angmx/naa
      angmy=angmy/naa
      angmz=angmz/naa
      angx=angmx*angi(1,1)+angmy*angi(1,2)+angmz*angi(1,3)
      angy=angmy*angi(2,1)+angmy*angi(2,2)+angmz*angi(2,3)
      angz=angmz*angi(3,1)+angmy*angi(3,2)+angmz*angi(3,3)
      do 124 i=1,naay
      xj=(angy*x(3,i)-x(2,i)*angz)
      yj=(angz*x(1,i)-x(3,i)*angx)
      zj=(angx*x(2,i)-x(1,i)*angy)
      velo(1,i)=velo(1,i)-xj
      velo(2,i)=velo(2,i)-yj
      velo(3,i)=velo(3,i)-zj
 124  continue
 999  return
      end
      SUBROUTINE MDYN
C THIS SUBROUTINE DOES THE MOLECULAR DYNAMICS CALCULATION
C
      DOUBLE PRECISION XK,YK,ZK,xkc,ykc,zkc,virx,viry,virz,
     1vircx,vircy,vircz,viryx,virzx,virzy
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      PARAMETER(IX=301,IY=IX,IX3=IX*3,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),rbend(ix),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /OPTIM/virx,viry,virz,XD(IX3),XDD(IX3),QZ(IX3),DX(IX3),
     1XK(IX),YK(IX),ZK(IX),NAM,EPP,EPW
      COMMON /DYN/ TIME,TEMP,GMASS,MDCYC,TEMITT,EMITT,TMASS,MDUMP,MCOOL,
     1SCHW(3),VELO(3,IX),ACC(3,IX),X3(3,IX),TSOLL,IRAND,TE,iav
     2,pnull,taupe,taupi,pmitt,taute,pstart,pdiff,pend,mpress
C
cp
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1 cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2 imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3 ,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
      common /pacplo/ npac,ispac(2000),cpac(3,2000),ihhyd(100),
     1 iahyd(100),dhyd(100),nhyd,jhyd(ix),xo(3,ix),xsum(3,ix),
     2 ibhyd(100),hhyd(100),abhyd(100),anghyd(100)   
      common /dim/ xmax(3),xmin(3),vela(3,ix),xa(3,ix),
     1rneu(ix),thneu(ix),phineu(ix),thalt(ix),phialt(ix)
cp
cs
      common /solvat/ eps,dielec,dsgl,fav(ix),fam(ix),rxy(ix,ix),
     1 rvdw(50),ion(50)
      common /dyna/ idz,idd(4,40),irdz(ix6),idyn,dystar
cs
      save idump
      naax=naa
      NAM=NAAx-1
      KA=3*NAAx
C
C    SET FIRST TRIAL OF HESSIAN MATRIX
C
      EPP=0.
      EPW=0.
      DO 1 I=1,NAAx
      XK(I)=0.D0
      YK(I)=0.D0
      ZK(I)=0.D0
    1 CONTINUE
      virx=0.
      viry=0.
      virz=0.
      DO 2 I=1,NAA
      DO 2 J=1,NAA
      IF (VS(I,J).GT.0.) GOTO 2
      VS(I,J)=0.
    2 CONTINUE
      t2=time/2.
C
C recalculate the forces
C
      CALL BOND
      CALL ANGLE
      if (mod(iflag2,20).eq.0) then 
      CALL GAMMA(NAA,ISORT,RV)
      CALL CHARGE
      end if
      CALL TORS
      CALL BEND
cp
      if (nm.lt.-5) then
      call crypac
      do 99 i=1,naa
      xk(i)=xk(i)+xkc(i)
      yk(i)=yk(i)+ykc(i)
      zk(i)=zk(i)+zkc(i)
   99 continue
      end if
C
      CALL NONBON
      EP=epp+EPW
C
C calculate unconstrained accelerations at time t 
C
c
c      entferne Transation und Rotation
c        
c      if (nm.gt.-6)
      call mdrotr
c
      ekinx=0.
      ekiny=0.
      ekinz=0.
      DO 10 I=1,NAAx
      ISO1=ISORT(I)+1
      acc(1,i)=xk(i)/amass(iso1)*120.44
      acc(2,i)=yk(i)/amass(iso1)*120.44
      acc(3,i)=zk(i)/amass(iso1)*120.44
      ekinx=ekinx+amass(iso1)*velo(1,i)**2
      ekiny=ekiny+amass(iso1)*velo(2,i)**2
      ekinz=ekinz+amass(iso1)*velo(3,i)**2
      do 1010 m=1,3
      vela(m,i)=velo(m,i)
 1010 xa(m,i)=x(m,i)
   10 continue
      ekinx=ekinx*51.82
      ekiny=ekiny*51.82
      ekinz=ekinz*51.82
      ekin1=ekinx+ekiny+ekinz
      ekin=ekin1
      TE=EKIN1/(NAAx)*8165.5
      tex=ekinx/(naax)*24496.5
      tey=ekiny/(naax)*24496.5
      tez=ekinz/(naax)*24496.5
c
c  pressure
c
      prm=0.
      if (nm.lt.-5) then
      IF (IFLAG2.GE.Mpress) THEN
      if (pnull.lt.pend) pnull=pnull+pdiff
      end if
      taupi=taupe
      if (iflag2.ge.1000) taupi=taupe
      zuv=1066000.*ksym/(cl(1,1)*cl(2,2)*cl(3,3))/vvol
      vvirx=3121.*vircx
      vviry=3121.*vircy
      vvirz=3121.*vircz
      prx=zuv*(ekinx+vvirx)
      pry=zuv*(ekiny+vviry)
      prz=zuv*(ekinz+vvirz)
      prxy=zuv*(3121.*viryx)
      prxz=zuv*(3121.*virzx) 
      pryz=zuv*(3121.*virzy) 
      prm=(prx+pry+prz)/3.
      if (iscl.eq.1) then
      prx=prm
      pry=prm
      prz=prm
      end if
      if (iscl.eq.2) then
      prx=(prx+pry)/2
      pry=prx
      end if
      xmu=1.-time/(3.*taupi)*(pnull-prx)
      ymu=1.-time/(3.*taupi)*(pnull-pry)
      zmu=1.-time/(3.*taupi)*(pnull-prz)
      xymu=-time/(3.*taupi)*(-prxy)
      xzmu=-time/(3.*taupi)*(-prxz)
      yzmu=-time/(3.*taupi)*(-pryz)
      end if
C
C  ADJUST REQUIRED TEMPERATURE IF NECESSARY
C
      IF (TEMP.LT.TSOLL) TEMP=TEMP+0.1
      IF (IFLAG2.GE.MCOOL) THEN
      temp=temp-0.01
      TSOLL=temp
      TEMP=amax1(TEMP,0.)
      END IF
C
C CALCULATE CORRECTION FACTOR FOR TEMPERATURE
C
       tcor=sqrt(1.+time/taute*(temp/te-1.))
       tcorx=sqrt(1.+time/taute*(temp/tex-1.))
       tcory=sqrt(1.+time/taute*(temp/tey-1.))
       tcorz=sqrt(1.+time/taute*(temp/tez-1.))
C
c calculate constrained velocities at time t+dt/2
c                                             
c
      DO 20 I=1,NAAx   
      ISO=ISORT(I)+1   
      velo(1,i)=tcorx*(velo(1,i)+time*acc(1,i))
      velo(2,i)=tcory*(velo(2,i)+time*acc(2,i))
      velo(3,i)=tcorz*(velo(3,i)+time*acc(3,i))
   20 continue
      do 15  i=1,naax
      iso=isort(i)+1
      x(1,i)=x(1,i)+time*velo(1,i)
      x(2,i)=x(2,i)+time*velo(2,i)
      x(3,i)=x(3,i)+time*velo(3,i)
   15 continue
c
      if (nm.lt.-6) then
      x(1,naa+1)=x(1,naa+1)*xmu
      x(2,naa+2)=x(2,naa+2)*ymu
      x(3,naa+3)=x(3,naa+3)*zmu
      x(1,naa+2)=x(1,naa+2)+xymu*ascl(1,2)
      x(1,naa+3)=x(1,naa+3)+xzmu*ascl(1,3)
      x(2,naa+3)=x(2,naa+3)+yzmu*ascl(2,3)
      end if
c
      EPOT=EPP+EPW
      EEPOT=96.49*EPOT+dsgl*4.184
      EEKIN=EKIN*96.49
      EGES=EEKIN+EEPOT
      if (temp.ge.tsoll) then
      TEMITT=(TEMITT*idump+TE)/(Idump+1)
      if (nm.lt.-5) pMITT=(pMITT*idump+prm)/(idump+1)
      EMITT=(EMITT*Idump+EGES)/(idump+1)
      end if
      idump=mod(iflag2,mdump)
      TIMES=TIME*IFLAG2
      if (nm.lt.-5) call crtcry
      IF (idump.EQ.0) then 
      WRITE(IUT,100) TIMES,TE,prm,EGES,EEPOT,EEKIN,TEMITT,EMITT,
     1pmitt
c
      if (nh.lt.0) call dynout
      if (nm.lt.-5) then 
      call outcry
      do 1050 j=1,3
      xmin(j)=1000.
      xmax(j)=-1000.
      do 1040 i=1,naa
      xmax(j)=amax1(xmax(j),xo(j,i))
      xmin(j)=amin1(xmin(j),xo(j,i))
 1040 continue
      xmax(j)=xmax(j)+cutoff
      xmin(j)=xmin(j)-cutoff
 1050 continue
      imix=int(xmin(1)/cl(1,1))-1
      imiy=int(xmin(2)/cl(2,2))-1
      imiz=int(xmin(3)/cl(3,3))-1
      imax=int(xmax(1)/cl(1,1))+1
      imay=int(xmax(2)/cl(2,2))+1
      imaz=int(xmax(3)/cl(3,3))+1
      end if
c
  100 FORMAT(2F8.1,f9.2,6F9.2)
      end if
      if (nm.lt.-5) then
      do 201 njj=1,naax
      do 201 nii=1,3
  201 xsum(nii,njj)=xsum(nii,njj)+x(nii,njj)
      if (mod(iflag2,iav).eq.0.and.iflag2.gt.20) then
      do 1001 njj=1,naax
      do 1001 nii=1,3
      xo(nii,njj)=x(nii,njj)
 1001 xsum(nii,njj)=0.
      end if
      end if     
      RETURN
      END
      FUNCTION RANDOM()
      PARAMETER(IX=301)
      COMMON /DYN/ TIME,TEMP,GMASS,MDCYC,TEMITT,EMITT,TMASS,MDUMP,MCOOL,
     1SCHW(3),VELO(3,IX),ACC(3,IX),X3(3,IX),TSOLL,IRAND,TE,iav
     2,pnull,taupe,taupi,pmitt,taute,pstart,pdiff,pend,mpress
      IRAND=IRAND*17+289
      IRAND=IRAND-INT(FLOAT(IRAND)/34679.)*34679
      RANDOM=FLOAT(17340-IRAND)
      RETURN
      END
      subroutine dynout
      DOUBLE PRECISION XK,YK,ZK,xkc,ykc,zkc,virx,viry,virz,
     1vircx,vircy,vircz,viryx,virzx,virzy
      character*4 cc,ct
      COMMON /DIRECT/ NA,NB,NC,ND,NE,NF,NG,NH,NI,NK,NL,NM,IFLAG1,IFLAG2,
     1NAA,ZBB,IIN,IUT,IPLO,IPU,IGEO,TITLE(39),ISTOP,ISEQ,ISER,DHFM,
     2INFORM,IPDB,IESP,IHB,IRES,IARC
      COMMON /PARAM/ AL(50,50),BL(50,50),HV(50),RV(50),ZS(50),
     1AME(50,50),DE(50,50),RBOND(50),ZV(50),DHFA(50),DF(50,50),DBS(50),
     2GANG(50),PK2(50,50),PK4(50,50),PK6(50,50),DSW(50),DECUB(50),
     3PPIZ(50,50),QQ(50),AVDW(50,50),BVDW(50,50),CVDW(50,50),IVDW(50,
     450),PTT(50,50),BK(50,50),PPIT(50,50),AMASS(50),XN(50),B(50),C(50),
     5XP(50),ISPI(50),ppix(50,50)
      PARAMETER(IX=301,IY=IX,IX3=IX*3,IX4=IX*4,IX6=IX*6)
      COMMON /GEOM1/ X(3,IX),ISORT(IX),A(3,IX),NRA(4,3,IX),IRST(4,IX),
     1IRZ,ARST(IX),NUMM(IX),AB(3,IX),IFIX(IX),IMOL(IX),IAT(IX)
      COMMON /GEOM3/ VV(IX,IX),VS(IX,IX),EP,SWI(IX),NN,P(IY,IY),CS(IX),
     1QSIG(IX),IBOND,NBOND(2,IX4),RBON(IX4),IANG,NANG(3,IX6),RANG(IX6),
     2ITORS,NTORS(4,IX6),RTORS(IX6),IBEND,NBEND(4,IX),RBEND(IX),
     3IDELT,NDELT(4,IX6),RDELT(IX6),QSEFF(IX),QADD(IX),NEPS,ctors(ix6)
      COMMON /OPTIM/virx,viry,virz,XD(IX3),XDD(IX3),QZ(IX3),DX(IX3),
     1XK(IX),YK(IX),ZK(IX),NAM,EPP,EPW
      COMMON /DYN/ TIME,TEMP,GMASS,MDCYC,TEMITT,EMITT,TMASS,MDUMP,MCOOL,
     1SCHW(3),VELO(3,IX),ACC(3,IX),X3(3,IX),TSOLL,IRAND,TE,iav
     2,pnull,taupe,taupi,pmitt,taute,pstart,pdiff,pend,mpress
C
cp
      common /crys/ ck(3),cal,cbe,cga,ksym,tsym(3,192),rsym(3,3,192),
     1 cl(3,3),epce,epcv,xpac(ix),cx(3,ix),xkc(ix),ykc(ix),zkc(ix),
     2 imix,imax,imiy,imay,imiz,imaz,iscl,ascl(3,3),cutoff,coufac
     3 ,vvol,xpos(3),vircx,vircy,vircz,viryx,virzx,virzy
      common /pacplo/ npac,ispac(2000),cpac(3,2000),ihhyd(100),
     1 iahyd(100),dhyd(100),nhyd,jhyd(ix),xo(3,ix),xsum(3,ix),
     2 ibhyd(100),hhyd(100),abhyd(100),anghyd(100)   
      common /dyna/ idz,idd(4,40),irdz(ix6),idyn,dystar
      dimension at(40),cc(8),ct(40),jt(40)
      data cc /'  c ','  c ',' g+ ',' g+t','  t ','  t ',' g-t',
     1         ' g- '/ 
      times=time*iflag2*0.001
      if (times.lt.dystar) return
      do 10 n=1,itors
c     FIND  THE ATOMS THAT DEFINE A TORSINAL ANGLE K -- I -- J -- L
C
         I=NTORS(1,N)
         J=NTORS(2,N)
         K=NTORS(3,N)
         L=NTORS(4,N)
         ators=rtors(n)*57.296
         do 20 nn=1,idz
         if (irdz(n).eq.nn) then
         at(nn)=ators
         end if
   20 continue 
   10 continue
      if (nh.eq.-1) then
      write(idyn,200) times
      write(idyn,201) (at(j),j=1,idz)
      end if
  200 format (f5.1)
  201 format (40f5.0)
      if (nh.le.-2) then
      do 100 i=1,idz
      if (at(i).ge.-30..and.at(i).lt.0.) then
      ct(i)=cc(1)
      jt(i)=1
      goto 100
      endif
      if (at(i).ge.0..and.at(i).lt.30.) then
      ct(i)=cc(2)
      jt(i)=1
      goto 100
      end if
      if (at(i).ge.30..and.at(i).lt.90.) then
      ct(i)=cc(3)
      jt(i)=2
      goto 100
      end if
      if (at(i).ge.90..and.at(i).lt.150) then
      ct(i)=cc(4)
      jt(i)=3
      goto 100
      end if
      if (at(i).ge.150..and.at(i).le.181.) then
      ct(i)=cc(5)
      jt(i)=3
      goto 100
      end if
      if (at(i).ge.-90..and.at(i).lt.-30.) then
      ct(i)=cc(8)
      jt(i)=5
      goto 100
      end if
      if (at(i).ge.-150..and.at(i).lt.-90.) then
      ct(i)=cc(7)
      jt(i)=6
      goto 100
      end if
      if (at(i).ge.-181..and.at(i).lt.-150.) then
      ct(i)=cc(6)
      jt(i)=4
      goto 100
      end if
      ct(i)=' ?  '
      jt(i)=10
  100 continue
      write (idyn,400) times,(jt(j),j=1,idz)
  400 format (f5.1,10x,40i1)
      end if
      return
      end
