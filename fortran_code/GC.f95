MODULE DATA
  !This is the data module used with PROGRAM GC.
  IMPLICIT NONE

  !Declare the fundamental constants that are named constants.
  !KELVIN is the conversion factor for dividing the temperature in K to
  !  get the energy in hartree.
  !IDEALGAS is the molar gas constant in J/mole-K.
  !AVOGADRO is the number of items per mole.
  !HARTREE is the number of Joules in one hartree.
  !BOHR is the bohr radius of a 1s electron in a H atom, in m.
  !FARADAY is the amount of charge in one mole of protons, in C/mole.
  !TOWNSEND is the unit of E/N, in V-m^2.
  !LOSCHMIDT is the standard number density of an ideal gas, in 1/m**3.
  REAL(SELECTED_REAL_KIND(18)),PARAMETER::KELVIN=3.1577504D5,IDEALGAS=&
    8.3144621D0,AVOGADRO=6.02214129D23,HARTREE=4.35974434D-18,&
    BOHR=5.2917721092D-11,FARADAY=9.64853365D4,TOWNSEND=1.0D-21,&
    LOSCHMIDT=2.6867805D25

  !Declare the global real quantities that are named constants.
  !ZERO, ONE, TWO, THREE, FOUR and FIVE are commonly used numbers.
  REAL(SELECTED_REAL_KIND(18)),PARAMETER::ZERO=0,ONE=1,TWO=2,THREE=3,&
    FOUR=4,FIVE=5

  !Declare the global integers that are named constants.
  !MAXE is the maximum number of energies in each of the 3 regions.
  !  It should be one larger than 2 to some power.
  !MAXL is the number of transport cross sections to be calculated.
  INTEGER,PARAMETER::MAXE=65,MAXL=30

  !Declare the global integers.
  !NMAXCH is the maximum order of approximation in using the Curtiss-
  !  Clenshaw method of integration.  It is important that the value
  !  of NMAXCH here be the same as in CCgen.f90.
  !NPTSCH is the size of the arrays CHEPTS and CHEWTS.
  !ERR allows the program to recover from input/output errors.
  !N1, N2 and N3 are the number of cross sections in each energy region.
  INTEGER::NMAXCH,NPTSCH,ERR,N1,N2,N3

  !Declare the global real variables and arrays.
  !PI is the usual 3.14159....
  !EMIN1, EMIN2, EMIN3, and EMIN4 are the boundaries of the three
  !  energy regions for the cross sections.
  !Array H holds factorials such that H(N)=GAMMA(N)
  !Array G holds gamma functions such G(N)=GAMMA(N-0.5D0)/GAMMA(0.5D0)
  !The quadrature points and weights are placed in CHEPTS and CHEWTS.
  !FILE2 holds the cross sections computed previously using program PC.
  REAL(SELECTED_REAL_KIND(18))::PI,EMIN1,EMIN2,EMIN3,EMIN4
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(31)::H,G
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(:),ALLOCATABLE::CHEPTS,CHEWTS
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(3*MAXE,MAXL)::FILE2

  !Declare the global character variables.
  !  CCFILE is the fully qualified name of the input file.
  !  FNUM indicates the extension of output files, i.e. the temperature.
  CHARACTER(LEN=78)::CCFILE
  CHARACTER(LEN=3)::FNUM
END MODULE DATA

PROGRAM GC
  !This is a Fortran-95 version of the Gram-Charlier program based on a
  !paper in Chemical Physics 179 (1994) 71.  This version was written
  !in October, 2018, by Dr. Larry A. Viehland, Science Department,
  !Chatham University, Pittsburgh, Pennsylvania 15232, USA.  He may be
  !contacted by phone at 1-412-365-2752, by fax at 1-412-365-1505, or
  !by email at viehland@chatham.edu.
  USE DATA
  IMPLICIT NONE

  !Declare local variables.
  !  I is a loop variable.
  INTEGER::I

  !Open the input file.
  OPEN(UNIT=10,FILE='CCgen.in',STATUS='OLD',IOSTAT=ERR)
  IF(ERR/=0)THEN
    PRINT *,'Could not open the old file named CCgen.in.'
    STOP
  ENDIF

  !Read NMAXCH and CCFILE.
  READ(10,*)NMAXCH
  READ(10,'(A)')CCFILE
  CLOSE(10)

  !Assign the value to PI.
  PI=FOUR*ATAN(ONE)

  !Set up arrays H and G.
  H(1)=ONE
  G(1)=ONE
  DO I=2,31
    H(I)=DBLE(I-1)*H(I-1)
    G(I)=DBLE(2*I-3)*G(I-1)/TWO
  END DO

  !Allocate space for the quadrature points and weights.
  NPTSCH=0
  DO I=1,NMAXCH
    NPTSCH=NPTSCH+1+2**I
  ENDDO
  ALLOCATE(CHEPTS(NPTSCH),CHEWTS(NPTSCH))

  !Read the quadrature points and weights.
  OPEN(10,FILE=CCFILE,IOSTAT=ERR,STATUS='OLD')
  IF(ERR/=0)THEN
    PRINT *,'Could not open the old file named ',TRIM(CCFILE),'.'
    PRINT *,'Be sure that you have permission to read this file.'
    STOP
  ENDIF
  DO I=1,NPTSCH
    READ(10,'(18X,2F23.20)',IOSTAT=ERR)CHEPTS(I),CHEWTS(I)
    IF(ERR/=0)THEN
      PRINT *,'File ',TRIM(CCFILE),' is incompatible with this &
        &program.'
      PRINT *,'Be sure that it was created with program CCgen.f95.'
      STOP
    ENDIF
  ENDDO

  !If everything is consistent, the next READ statement should fail.
  !Note that EMIN1 and EMIN2 are dummy variables at this point.
  READ(10,'(18X,2F23.20)',IOSTAT=ERR)EMIN1,EMIN2
  IF(ERR==0)THEN
    PRINT *,'File ',TRIM(CCFILE),' is incompatible with this &
      &program.'
    PRINT *,'Be sure that it was created with program CCgen.f95.'
    STOP
  ENDIF
  CLOSE(10)
  
  !Open the files with the parameters describing the cross sections.
  OPEN(UNIT=16,FILE='crosec',STATUS='OLD',ACCESS='DIRECT',RECL=496,&
    &IOSTAT=ERR)
  IF(ERR/=0)THEN
    PRINT *,'Error opening the old file named crosec.'
    STOP
  END IF
  OPEN(UNIT=17,FILE='coeffs',STATUS='OLD',ACCESS='DIRECT',RECL=480,&
    &IOSTAT=ERR)
  IF(ERR/=0)THEN
    PRINT *,'Error opening the old file named coeffs.'
    STOP
  END IF
  OPEN(18,FILE='ninval',IOSTAT=ERR)
  IF(ERR/=0)THEN
    PRINT *,'Error opening the old file named ninval.'
    STOP
  ENDIF

  !Read the number of cross sections in each energy region.
  READ(18,*,IOSTAT=ERR)N1
    IF(ERR/=0)THEN
      PRINT *,'Error reading N1 from file ninval.'
      STOP
    END IF
    READ(18,*,IOSTAT=ERR)N2
    IF(ERR/=0)THEN
      PRINT *,'Error reading N2 from file ninval.'
      STOP
    END IF
    READ(18,*,IOSTAT=ERR)N3
    IF(ERR/=0)THEN
      PRINT *,'Error reading N3 from file ninval.'
      STOP
    END IF
  CLOSE(18)

  EMIN1=ZERO
  EMIN2=ZERO
  EMIN3=ZERO
  EMIN4=ZERO
  IF(N1>0)THEN
    READ(16,REC=1,IOSTAT=ERR)EMIN1
    IF(ERR/=0)THEN
      PRINT *,'Error reading EMIN1 from line 1 of file crosec.'
      STOP
    END IF
    READ(16,REC=N1,IOSTAT=ERR)EMIN2
    IF(ERR/=0)THEN
      PRINT *,'Error reading EMIN2 from line ',N1,' of file crosec.'
      STOP
    END IF
  END IF
  IF(N1==0.AND.N2>0)THEN
    READ(16,REC=MAXE+1,IOSTAT=ERR)EMIN2
    IF(ERR/=0)THEN
      PRINT *,'Error reading EMIN2 from line ',MAXE+1,' of file crosec.'
      STOP
    END IF
  END IF
  IF(N2>0)THEN
    READ(16,REC=(MAXE+N2),IOSTAT=ERR)EMIN3
    IF(ERR/=0)THEN
      PRINT *,'Error reading EMIN3 from line ',N2,' of file crosec.'
      STOP
    END IF
  END IF
  IF(N2==0.AND.N3>0)THEN
    READ(16,REC=2*MAXE+1,IOSTAT=ERR)EMIN3
    IF(ERR/=0)THEN
      PRINT *,'Error reading EMIN3 from line ',2*MAXE+1,' of file &
        &crosec.'
      STOP
    END IF
  END IF
  IF(N3>0)THEN
    READ(16,REC=(2*MAXE+N3),IOSTAT=ERR)EMIN4
    IF(ERR/=0)THEN
      PRINT *,'Error reading EMIN4 from line ',2*MAXE+N3,' of file &
        &crosec.'
      STOP
    END IF
  END IF

  !Read the cross section parameters.
  DO I=1,3*MAXE
    READ(17,REC=I,IOSTAT=ERR)FILE2(I,:)
    IF(ERR/=0)THEN
      PRINT *,'Error reading FILE2 from file coeffs.'
      STOP
    END IF
  END DO
  CLOSE(16)
  CLOSE(17)

  !Open the file containing the temperatures of interest.
  OPEN(UNIT=26,FILE='temps',IOSTAT=ERR)
  IF(ERR/=0)THEN
    PRINT *,'Error opening file temps.'
    STOP
  ENDIF

  !In a loop, look at each temperature of interest.
  DO
    READ(26,*,IOSTAT=ERR)FNUM
    IF(ERR/=0)EXIT
    CALL GRAMCH
  ENDDO
END PROGRAM GC

SUBROUTINE GRAMCH
  USE DATA
  IMPLICIT NONE

  !Declare the common variables.
  !ACCT, ACCT1, ACCT2 and ACCT3 are the accuracies of the cross
  !  sections, mobilities, diffusion coefficients and other transport
  !  properties, respectively.
  !AM1 and AM0 are the ion mass and neutral mass in g/mole.
  !TEMP is the gas temperature in Kelvin, later converted to hartree.
  !RAM is the ratio of the ion mass to the neutral mass.
  !TPARA and TPERP are the ion temperatures parallel and perpendicular
  !  to the electric field.
  !ALL is the skewness of the distribution.
  !BEL and BET are the kurtosis along and perpendicular to the field.
  !CL1 and CL2 are the two correlation coefficients.
  !TT1 is the value of e in the equations for the matrix elements.
  !TT2 and TT3 are the values of f_perp and f_para.
  !T10 and T11 are the values of d_perp and d_para.
  !SAVE1 is used to save the matrix elements.
  !SAVE2 is used to save (part of) the reducible collision integrals.
  !The dimensions of SAVE1 and SAVE2 depend on the maximum values
  !  allowed for IAPPM, which is presently 9.
  !XTRANS is the argument of the first integral of the pair involved
  !  in eq. A.15.
  !SAVEANS is used to save the irreducible collisions integrals.
  !IPTRAN, IQTRAN and ILTRAN are the indices for the irreducible
  !  collision integrals defined by eq. A.15.
  REAL(SELECTED_REAL_KIND(18))::ACCT,ACCT1,ACCT2,ACCT3,AM0,AM1,TEMP,&
    RAM,TPARA,TPERP,VDBAR,ALL,BEL,BET,CL1,CL2,TT1,TT2,TT3,TT10,TT11,&
    XTRANS
  REAL(SELECTED_REAL_KIND(18)),&
    DIMENSION(0:11,0:11,0:13,0:23,0:23,0:23)::SAVE1,SAVE2
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(0:MAXL,0:MAXL,1:MAXL)::&
    SAVEANS
  INTEGER::IPTRAN,IQTRAN,ILTRAN
  COMMON/REAL0/ACCT,ACCT1,ACCT2,ACCT3
  COMMON/REAL1/AM0,AM1,TEMP,RAM,TPARA,TPERP,VDBAR
  COMMON/REAL2/ALL,BEL,BET,CL1,CL2
  COMMON/REAL3/TT1,TT2,TT3,TT10,TT11
  COMMON/REAL4/SAVE1,SAVE2
  COMMON/REAL5/XTRANS
  COMMON/REAL6/SAVEANS
  COMMON/INT1/ILTRAN,IPTRAN,IQTRAN

  !Declare the local character variables.
  !SYSTEM is a short (15 characters or less) name for the system.
  !COMM is a long (78 characters or less) desription of the potential.
  !FORMAT is used to adjust the printing format when numbers get large.
  CHARACTER(LEN=15)::SYSTEM
  CHARACTER(LEN=76)::COMM
  CHARACTER(LEN=34)::FORMAT

  !Declare the local real variables and arrays.
  !VDS and VDEND are the initial and final values of the dimensionless
  !  drift velocity VDSN=SQRT(M0/(2*k_B*TEMP)*Vd, where k_B is
  !  Boltzmann's constant, vd is the drift velocity, M0 is the neutral
  !  mass, and TEMP is the gas temperature--all in basic SI units.
  !ALPHMU is the neutral polarizability in cubic Angstrom.  It is later
  !  changed to the product of the neutral polarizability in cubic
  !  Angstrom and the reduced mass in g/mole.  It has little use,
  !  since the program no longer computes the Langevin mobility, but is
  !  retained so that the input file GC.in is still like it was.
  !ELEM0, ELEM1, ELEM3, XN and XM are external functions.
  !W is the normalized difference between VDSN, the next estimate of
  !  VDS, and the present value.
  !AC is the value of E/N in Townsend
  !RMOB is the drift velocity in 10**4 cm/s.
  !ROK is the reduced mobility, Ko, in cm**2/Vs.
  !V1 and V2 are estimates of the ion temperatures parallel and
  !  perpendicular, respectively, to the electrostatic field.
  !EINS1 and EINS2 are the ratios qD/Kk_BT.
  !DIFF1 and DIFF2 are the values of NDpara and NDperp, in 10^20 /ms.
  !ROKOLD, V1OLD and V2OLD are previous value of ROK, V1 and V2.
  !ALLOLD, BELOLD, BETOLD, CL1OLD, CL2OLD, D1OLD and D2OLD are
  !  similarly related to BEL, BET, CL1, CL2, DIFF1 and DIFF2.
  !ACCN1 is the fractional difference to be compared to ACCT1.
  !ACCN2 is the fractional difference to be compared to ACCT2.
  !ACCN3 is the fractional difference to be compared to ACCT3.
  !AMAT contains the matrix elements in the matrix equations.
  !BMAT and CMAT are the inhomogenous terms in the matrix equations.
  !ALLNEW is the value of ALL obtained either from the input file (for
  !  the first value of VDS) or as the final result from the previous
  !  value of VDS.
  !BELNEW, BETNEW, CL1NEW, and CL2NEW are similarly related to BEL, BET,
  ! CL1 and CL2.
  !VAL1-VAL4 are used to print breakpoints when STANDERR is true.
  REAL(SELECTED_REAL_KIND(18))::VDS,VDEND,ALLNEW,BELNEW,BETNEW,CL1NEW,&
    CL2NEW,BB
  REAL(SELECTED_REAL_KIND(18))::ALPHMU,ELEM0,ELEM1,ELEM3,XN,XM,W,VDSN,&
    AC,RMOB,ROK,V1,V2,EINS1,EINS2,DIFF1,DIFF2,ROKOLD,V1OLD,V2OLD,&
    ALLOLD,BELOLD,BETOLD,CL1OLD,CL2OLD,D1OLD,D2OLD,ACCN1,ACCN2,ACCN3,&
    VAL1,VAL2,VAL3,VAL4
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(:),ALLOCATABLE::BMAT,CMAT
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(:,:),ALLOCATABLE::AMAT

  !Declare the local integer variables and arrays.
  !ICHARG is the ion charge in units of the charge of a proton.  Note:
  !  the mobility is independent of the sign of ICHARG, except to the
  !  extent that there is a difference in the interaction potential.
  !I and K are loop indices.
  !II1, JJ1 and JJ2 label the matrix elements and moments.
  !IP, IQ, IR, IS, IT and IU are the indices p, q, r, s, t, and u.
  !KK is a flag for errors in subroutine SIMQE.
  !MAXEQ is the maximum number of simultaneous equations that subroutine
  !  SIMQE will be expected to handle.
  !IAPPM is the maximum number of approximations.
  !IAPPN is the minimum number of approximations.
  !IAPP is the current order of approximation.
  INTEGER::ICHARG,I,K,II1,JJ1,JJ2,IP,IQ,IR,IS,IT,IU,KK,MAXEQ,IAPP,&
    &IAPPM,IAPPN

  !Declare the local logical variables.
  !HALT is true if the program should stop when it fails to converge
  !  at any value of E/N
  !PRINTALL is true if best estimate should be printed even when the
  !  program fails to converge.  It is false if only converged values
  !  should be printed.
  !STANDERR indicates whether the standard sequence is used for
  !  successive values of ACCT1, ACCT2 and ACCT3.
  !MOBERR is true when only convergence of the mobility will be checked.
  !RESTART indicates that a calculation is being restarted at a high
  !  value of E/N.
  LOGICAL::HALT,PRINTALL,STANDERR,MOBERR,RESTART

  !Open the input and output files.
  OPEN(UNIT=13,FILE='GC.in.'//FNUM,STATUS='OLD',IOSTAT=ERR)
  IF(ERR/=0)THEN
    PRINT *,'Error opening the old file named GC.in.'//FNUM
    RETURN
  END IF
  OPEN(UNIT=14,FILE='output.'//FNUM,IOSTAT=ERR)
  IF(ERR/=0)THEN
    PRINT *,'Error opening the output file named output.'//FNUM
    RETURN
  END IF
  OPEN(UNIT=15,FILE='junk.'//FNUM,IOSTAT=ERR)
  IF(ERR/=0)THEN
    PRINT *,'Error opening the output file named junk.'//FNUM
    RETURN
  END IF

  !Read the system name, a comment line, the ion mass, the neutral mass,
  !  the gas temperature, the ion charge, and the polarizability of the
  !  neutral gas.  Write the same quantities to the file named output 
  !  and at the end of summary file named transport.
  READ(13,'(A)',IOSTAT=ERR)SYSTEM
  IF(ERR/=0)THEN
    PRINT *,'Error reading line 1 of file GC.in.'//FNUM
    RETURN
  END IF
  WRITE(14,'(25X,A)')TRIM(SYSTEM)
  WRITE(15,'(/,25X,A)')TRIM(SYSTEM)
  READ(13,'(A)',IOSTAT=ERR)COMM
  IF(ERR/=0)THEN
    PRINT *,'Error reading line 2 of file GC.in.'//FNUM
    RETURN
  END IF
  WRITE(14,'(5X,A,/)')TRIM(COMM)
  WRITE(15,'(5X,A,/)')TRIM(COMM)
  READ(13,*,IOSTAT=ERR)AM1
  IF(ERR/=0)THEN
    PRINT *,'Error reading line 3 of file GC.in.'//FNUM
    RETURN
  END IF
  WRITE(14,'(5X,A,F12.6,A)')'The ion mass is ',AM1,' g/mole.'
  READ(13,*,IOSTAT=ERR)AM0
  IF(ERR/=0)THEN
    PRINT *,'Error reading line 4 of file GC.in.'//FNUM
    RETURN
  END IF
  WRITE(14,'(5X,A,F12.6,A)')'The neutral mass is ',AM0,' g/mole.'
  READ(13,*,IOSTAT=ERR)TEMP
  IF(ERR/=0)THEN
    PRINT *,'Error reading line 5 of file GC.in.'//FNUM
    RETURN
  END IF
  WRITE(14,'(5X,A,F12.4,A)')'The gas temperature is ',TEMP,' K.'
  WRITE(15,'(5X,A,F12.2,A)')'The gas temperature is ',TEMP,' K.'
  READ(13,*,IOSTAT=ERR)ICHARG
  IF(ERR/=0)THEN
    PRINT *,'Error reading line 6 of file GC.in.'//FNUM
    RETURN
  END IF
  WRITE(14,'(5X,A,I2,A)')'The ion charge is ',ICHARG,&
    &' times the proton charge.'
  ICHARG=ABS(ICHARG)
  READ(13,*,IOSTAT=ERR)ALPHMU
  IF(ERR/=0)THEN
    PRINT *,'Error reading line 7 of file GC.in.'//FNUM
    RETURN
  END IF
  WRITE(14,'(5X,A,F10.6,A)')'The neutral polarizability is ',ALPHMU,&
    &' cubic Angstroms.'

  !Change TEMP to atomic units of energy and ALPHMU to be the product of
  !  the polarizability in cubic Angstrom and the ion-neutral reduced
  !  mass in g/mole.  Print the polarization mobility on the long 
  !  output file.
  TEMP=TEMP/KELVIN
  ALPHMU=ALPHMU*AM1*AM0/(AM1+AM0)
  WRITE(14,'(5X,A,F7.4,A)')'The polarization limit of the mobility is &
    &',13.853D0/SQRT(ALPHMU),' cm^2/Vs'

  !Calculate RAM.  Then read the initial value of vd*, the final value
  !  of vd*, TPARA, TPERP and the maximum order of approximation
  RAM=AM1/AM0
  READ(13,*,IOSTAT=ERR)VDS
  IF(ERR/=0)THEN
    PRINT *,'Error reading line 8 of file GC.in.'//FNUM
    RETURN
  END IF
  READ(13,*,IOSTAT=ERR)VDEND
  IF(ERR/=0)THEN
    PRINT *,'Error reading line 9 of file GC.in.'//FNUM
    RETURN
  END IF
  WRITE(14,'(5X,A,F12.6,A,F12.6,A)')'The range of vd* is ',VDS,' to ',&
    &VDEND,'.'
  READ(13,*,IOSTAT=ERR)IAPPN,IAPPM
  IF(ERR/=0)THEN
    PRINT *,'Error reading line 10 of file GC.in.'//FNUM
    RETURN
  END IF

  !Don't use anything less than 3 for the minimum order of approximation.
  IF(IAPPN<3)THEN
    PRINT *,'The program requires at least three approximations,'
    PRINT *,'so the value of the first integer in line 10 of file'
    PRINT *,'GC.in.'//FNUM//' has been changed to 3.'
    IAPPN=3
  ENDIF

  !The program can probably go to 28th order by changing the dimensions
  !  of various arrays, but so far it has been checked only through 10th
  !  approximation; this version is limited to 9th approximation.
  IF(IAPPM>9)THEN
    PRINT *,'The program is presently limited to no more than nine'
    PRINT *,'approximations, so the value of the second integer in'
    PRINT *,'line 10 of file C.in.'//FNUM//' has been changed to 9.'
    IAPPM=9
  END IF

  !Read the accuracies.
  READ(13,*,IOSTAT=ERR)ACCT
  IF(ERR/=0)THEN
    PRINT *,'Error reading line 11 of file GC.in.'//FNUM
    RETURN
  END IF
  READ(13,*,IOSTAT=ERR)ACCT1
  IF(ERR/=0)THEN
    PRINT *,'Error reading line 12 of file GC.in.'//FNUM
    RETURN
  END IF
  READ(13,*,IOSTAT=ERR)ACCT2
  IF(ERR/=0)THEN
    PRINT *,'Error reading line 13 of file GC.in.'//FNUM
    RETURN
  END IF
  READ(13,*,IOSTAT=ERR)ACCT3
  IF(ERR/=0)THEN
    PRINT *,'Error reading line 14 of file GC.in.'//FNUM
    RETURN
  END IF
  IF(ACCT>MAX(ACCT1,ACCT2,ACCT3))THEN
    PRINT *,'Since the accuracy of the transport coefficients may'
    PRINT *,'not exceed the accuracy of the cross sections, they'
    PRINT *,'have been increased above the requested values.'
    ACCT1=ACCT
    ACCT2=MAX(ACCT,ACCT2)
    ACCT3=MAX(ACCT,ACCT3)
  END IF

  !Read the logical variables controlling the calculations.
  READ(13,*,IOSTAT=ERR)HALT
  IF(ERR/=0)THEN
    HALT=.TRUE.
  ENDIF
  READ(13,*,IOSTAT=ERR)PRINTALL
  IF(ERR/=0)THEN
    PRINTALL=.FALSE.
  ENDIF
  READ(13,*,IOSTAT=ERR)STANDERR
  IF(ERR/=0)THEN
    STANDERR=.FALSE.
  ENDIF
  READ(13,*,IOSTAT=ERR)MOBERR
  IF(ERR/=0)THEN
    STANDERR=.FALSE.
  ENDIF
  READ(13,*,IOSTAT=ERR)RESTART
  IF(ERR/=0)THEN
    RESTART=.FALSE.
  ENDIF

  !If RESTART is true, read the parameter values needed to set up 
  !the first, non-negligible, E/N value.
  IF(RESTART)THEN
    READ(13,*,IOSTAT=ERR)ALLNEW
    IF(ERR/=0)THEN
      PRINT *,'Error reading line 14 of file GC.in.'//FNUM
      RETURN
    END IF
    READ(13,*,IOSTAT=ERR)BELNEW
    IF(ERR/=0)THEN
      PRINT *,'Error reading line 15 of file GC.in.'//FNUM
      RETURN
    END IF
    READ(13,*,IOSTAT=ERR)BETNEW
    IF(ERR/=0)THEN
      PRINT *,'Error reading line 16 of file GC.in.'//FNUM
      RETURN
    END IF
    READ(13,*,IOSTAT=ERR)CL1NEW
    IF(ERR/=0)THEN
      PRINT *,'Error reading line 17 of file GC.in.'//FNUM
      RETURN
    END IF
    READ(13,*,IOSTAT=ERR)CL2NEW
    IF(ERR/=0)THEN
      PRINT *,'Error reading line 18 of file GC.in.'//FNUM
      RETURN
    END IF
    READ(13,*,IOSTAT=ERR)V1
    IF(ERR/=0.OR.V1/KELVIN<=TEMP)THEN
      !In cases where V1 can't be read from the input file or is too low,
      !  use the Maxwell model with A*=1 to make initial guesses; see eq.
      !  68 of Viehland and Mason, Ann. Phys. (NY) 110 (1978) 287.
      V1=(ONE+TWO*VDS*VDS*(ONE+THREE*RAM)/(THREE+FIVE*RAM))*TEMP
    ELSE
      V1=V1/KELVIN
    END IF
    READ(13,*,IOSTAT=ERR)V2
    IF(ERR/=0.OR.V2/KELVIN<=TEMP)THEN
      !In cases where V2 can't be read from the input file or is too low,
      !  use the Maxwell model with A*=1 to make initial guesses; see eq.
      !  69 of Viehland and Mason, Ann. Phys. (NY) 110 (1978) 287.
      V2=(ONE+TWO*VDS*VDS*(ONE+RAM)/(THREE+FIVE*RAM))*TEMP
    ELSE
      V2=V2/KELVIN
    END IF
  ELSE
    ALLNEW=ZERO
    BELNEW=THREE
    BETNEW=THREE
    CL1NEW=ZERO
    CL2NEW=ONE
    V1=TEMP
    V2=TEMP
  ENDIF
  CLOSE(13)

  !Write the appropriate comments to the large output file.
  WRITE(14,'(5X,A)')&
    &'Successive approximations of the kinetic theory are'
  WRITE(14,'(5X,A,I2,A)')'continued through at least the ',IAPPN,&
    &'-th approximation, and'
  WRITE(14,'(5X,A,I2,A)')'stopped after either the ',IAPPM,&
    &'-th approximation, or'
  IF(MOBERR)THEN
    WRITE(14,'(5X,A,F9.6,A)')'the mobility Ko has stabilized within ',&
      &100.0*ACCT1,'%.'
    WRITE(14,'(5X,A)')'However, convergence is not checked for the'
    WRITE(14,'(5X,A)')'diffusion coefficients or other quantities.'
  ELSE
    WRITE(14,'(5X,A,F9.6,A)')'the mobility Ko has stabilized within ',&
      &100.0*ACCT1,'%,'
    WRITE(14,'(5X,A,F9.6,A)')'the diffusion coefficients within ',&
      &100.0D0*ACCT2,'%,'
    WRITE(14,'(5X,A,F9.6,A)')'and the other quantities within ',&
      &100.0D0*ACCT3,'%.'
  ENDIF
  WRITE(14,'(5X,A,F9.6,A)')'The collision integrals have an accuracy &
    &of ',100.0D0*ACCT,'%'
  CLOSE(14)

  !Write the appropriate comments to the small output file.
  IF(MOBERR)THEN
    WRITE(15,'(5X,A,F4.2,A)')&
      &'Fractional accuracy for Ko is ',100.0*ACCT1,'%,'
    WRITE(15,'(5X,A)')'However, convergence is not checked for the'
    WRITE(15,'(5X,A)')'diffusion coefficients or other quantities'
    IF(IAPPN<IAPPM)THEN
      WRITE(15,*)
    ELSE
      WRITE(15,'(5X,A,I2)')'Minimum and maximum approximations were &
        &both ',IAPPM
    ENDIF
  ELSE
    WRITE(15,'(5X,A,F4.2,A)')&
      &'Fractional accuracy for Ko is ',100.0*ACCT1,'%,'
    WRITE(15,'(5X,A,F11.4,A)')'for NDx and NDz is ',100.0*ACCT2,'%,'
    IF(IAPPN<IAPPM)THEN
      WRITE(15,'(5X,A,F11.4,A,/)')'and for the other quantities is ',&
        &100.0*ACCT3,'%.'
    ELSE
      WRITE(15,'(5X,A,F11.4,A)')'and for the other quantities is ',&
        &100.0*ACCT3,'%.'
      WRITE(15,'(5X,A,I2)')'Minimum and maximum approximations were &
        &both ',IAPPM
    ENDIF
  ENDIF
  WRITE(15,'(2A,/)')'   E/N    Ko      Tz      Tx     NDz     NDx  ',&
    &' 1+SKEW KURz  KURx  1+C1  C2'
  CLOSE(15)

  !Allocate AMAT, BMAT and CMAT.
  MAXEQ=0
  DO IP=0,IAPPM+1,2
    DO IQ=0,IAPPM+1-IP,2
      DO IR=0,IAPPM+1-IP-IQ
        IF(IP+IQ+IR>0)MAXEQ=MAXEQ+1
      ENDDO
    ENDDO
  ENDDO
  ALLOCATE(AMAT(MAXEQ,MAXEQ),BMAT(MAXEQ),CMAT(MAXEQ))

  !This is the top of the loop for different values of VDS (indirectly,
  !  for different values of E/N).
  VDSloop: DO
    !Initialize the present approximations to the Gram-Charlier
    !  parameters at this value of VDS as equal to the values obtained
    !  from the previous E/N results.  Make sure, however, that these
    !  values are reasonable.
    ALL=MAX(ALLNEW,ZERO)
    BEL=MAX(BELNEW,THREE)
    BET=MAX(BETNEW,THREE)
    CL1=MAX(CL1NEW,ZERO)
    CL2=MAX(CL2NEW,ONE)
    TPARA=MAX(V1,(ONE+TWO*VDS*VDS*(ONE+THREE*RAM)/(THREE+FIVE*RAM))*TEMP)
    TPERP=MAX(V2,(ONE+TWO*VDS*VDS*(ONE+RAM)/(THREE+FIVE*RAM))*TEMP)
    VDBAR=VDS*SQRT(RAM/(RAM+TPARA/TEMP))

    !Establish the values used to compute the matrix elements.
    TT1=ONE/(ONE+RAM)
    TT2=RAM*(TPERP/TEMP-ONE)/(ONE+RAM)/(RAM+TPERP/TEMP)
    TT3=RAM*(TPARA/TEMP-ONE)/(ONE+RAM)/(RAM+TPARA/TEMP)
    TT10=TPERP/(RAM*TEMP+TPERP)
    TT11=TPARA/(RAM*TEMP+TPARA)

    !Before starting the successive approximations at this E/N value,
    !initialize SAVEANS and SAVE2.  Then set BB=ELEM1(0,0,1,0,0,0) as
    !the matrix element used to normalize the moment equations.
    !Finally, write VDS to file 14.
    SAVEANS=ZERO
    SAVE1=ZERO
    SAVE2=ZERO
    BB=ELEM1(0,0,1,0,0,0)
    OPEN(UNIT=14,FILE='output.'//FNUM,POSITION='APPEND')
    WRITE(14,*)
    WRITE(14,*)'               VDS=',VDS
    WRITE(14,'(/,3A)')'   E/N      Ko        Tz         Tx   1+SKEW',&
      &'  KURz   KURx 1+CORR1 CORR2     NDz        Dz/K       NDx',&
      &'       Dx/K  APP'

    !Now start the successive approximations.
    IAPPloop: DO IAPP=1,IAPPM
      !Close and reopen the long output file, so the tail can be read.
      CLOSE(14)
      OPEN(UNIT=14,FILE='output.'//FNUM,POSITION='APPEND')

      !We first set up the density-gradient-independent moment eq. (71)
      !  of Viehland, Chem. Phys. 179 (1994) 71.  We use BMAT for the
      !  inhomogenous terms in the matrix equations, and we label these
      !  terms with the index II1.
      II1=0
      IPloop1: DO IP=0,IAPP+1,2
        IF(IAPP==1)EXIT IPloop1
      IQloop1: DO IQ=0,IAPP+1-IP,2
      IRloop1: DO IR=0,IAPP+1-IP-IQ
        IF(IP+IQ+IR==0)CYCLE IRloop1
        II1=II1+1
        !Note that we assume that the quantity given by eq. 74 is
        !  equal to 1/2, which is why there is an extra TWO in the
        !  denominator of the next line.  Note also that XN is an
        !  external function to compute the quantity in eq. (47).
        BMAT(II1)=XN(IP,IQ,IR-1,0,0,0)*TWO**(IP+IQ+IR)/TWO*H(IP+1)*&
          &H(IQ+1)*H(IR+1)-ELEM0(IP,IQ,IR,0,0,0)/BB
        !We use AMAT for the matrix elements in the matrix equations,
        !  and we label these terms with the index JJ1.  Note that XM
        !  is an external function to compute the quantity in eq. (58).
        JJ1=0
        ISloop1: DO IS=0,IAPP+1,2
        ITloop1: DO IT=0,IAPP+1-IS,2
        IUloop1: DO IU=0,IAPP+1-IS-IT
          IF(IS+IT+IU==0)CYCLE IUloop1
          JJ1=JJ1+1
          AMAT(II1,JJ1)=ELEM0(IP,IQ,IR,IS,IT,IU)/BB-&
            &ELEM0(IP,IQ,IR,0,0,0)/BB*XN(0,0,0,IS,IT,IU)/&
            &TWO**(IS+IT+IU)/H(IS+1)/H(IT+1)/H(IU+1)-&
            &XM(IP,IQ,IR-1,IS,IT,IU)*TWO**(IP+IQ+IR)/&
            &TWO**(1+IS+IT+IU)*H(IP+1)/H(IS+1)*H(IQ+1)/H(IT+1)*&
            &H(IR+1)/H(IU+1)
        END DO IUloop1
        END DO ITloop1
        END DO ISloop1
      END DO IRloop1
      END DO IQloop1
      END DO IPloop1

      !Call SIMQE to solve for the moments needed for the mobility.
      IF(IAPP>1)THEN
        KK=0
        CALL SIMQE(AMAT,MAXEQ,BMAT,II1,KK)
        IF(KK==1)THEN
          WRITE(14,*)' Error in SIMQE for mobility at T='//FNUM
          RETURN
        END IF
      END IF

      !From the moments which subroutine SIMQE left in BMAT, calculate
      !  the next approximation to W by using eq. (54) of Viehland, 
      !  Chem. Phys. 179 (1994) 71.
      W=ZERO
      II1=0
      ISloop2: DO IS=0,IAPP+1,2
        IF(IAPP==1)EXIT ISloop2
      ITloop2: DO IT=0,IAPP+1-IS,2
      IUloop2: DO IU=0,IAPP+1-IS-IT
        IF(IS+IT+IU==0)CYCLE IUloop2
        II1=II1+1
        IF(BMAT(II1)==ZERO)CYCLE IUloop2

        !Change BMAT from the moment form to the expansion coefficient
        !  form, by using eq. 69 of Viehland, Chem. Phys. 179 (1994) 71.
        BMAT(II1)=BMAT(II1)/TWO**(IS+IT+IU)/H(IS+1)/H(IT+1)/H(IU+1)

        !Now we can use eq. (54).
        W=W+BMAT(II1)*XN(0,0,1,IS,IT,IU)
      END DO IUloop2
      END DO ITloop2
      END DO ISloop2
      
      !Change W into VDSN, a better estimate of VDS.
      VDSN=VDS+SQRT(TPARA/TEMP/RAM)*W

      !Calculate AC from the Ideal Gas Constant, the number of Kelvin
      !  per hartree, the number of Joules per hartree, Avogadro's
      !  number, the number of meters per bohr, and Faraday's number.
      AC=SQRT(TWO*IDEALGAS*AM0*TEMP*KELVIN*HARTREE*AVOGADRO)/&
        &TOWNSEND*BOHR**2/FARADAY/DBLE(ICHARG)*AM1*VDS/(AM0+AM1)*&
        &(ELEM3(0,0,1)+ELEM3(0,1,1)/VDBAR)

      !Calculate ROK and RMOB from the same fundamental constants and
      !  from Loschmidt's constant.
      ROK=SQRT(TWO*IDEALGAS*KELVIN*1.0D-5)/LOSCHMIDT*1.0D29
      ROK=ROK*VDSN/AC*SQRT(TEMP/AM0)
      RMOB=LOSCHMIDT*ROK*AC*1.0D-23

      !Now calculate the values of the ion temperatures, the skewness
      !  coefficient, the kurtosis coefficients, and the correlations.
      V1=ONE+TWO*W*W
      V2=ONE
      ALLNEW=TWO*THREE*XN(0,0,3,0,0,0)-THREE/TWO*W-W**3
      BELNEW=THREE/FOUR+TWO*THREE*FOUR*(XN(0,0,4,0,0,0)-W*&
        &XN(0,0,3,0,0,0))+THREE*W*W+W**4
      BETNEW=THREE/FOUR+TWO*THREE*FOUR*XN(4,0,0,0,0,0)
      CL1NEW=TWO*XN(2,0,1,0,0,0)-W/TWO
      CL2NEW=ONE/FOUR+FOUR*(XN(2,0,2,0,0,0)-W*XN(2,0,1,0,0,0))+&
        &W*W/TWO
      II1=0
      ISloop3: DO IS=0,IAPP+1,2
      ITloop3: DO IT=0,IAPP+1-IS,2
      IUloop3: DO IU=0,IAPP+1-IS-IT
        IF(IS+IT+IU==0)CYCLE IUloop3
        II1=II1+1
        IF(BMAT(II1)==ZERO)CYCLE IUloop3
        V1=V1+FOUR*BMAT(II1)*(XN(0,0,2,IS,IT,IU)-W*XN(0,0,1,IS,IT,IU))
        V2=V2+FOUR*BMAT(II1)*XN(2,0,0,IS,IT,IU)
        ALLNEW=ALLNEW+BMAT(II1)*(THREE/TWO*XN(0,0,1,IS,IT,IU)+TWO*&
          &THREE*(XM(0,0,3,IS,IT,IU)-W*XN(0,0,2,IS,IT,IU))+THREE*&
          &W*W*XN(0,0,1,IS,IT,IU))
        BELNEW=BELNEW+BMAT(II1)*(TWO*THREE*(XN(0,0,2,IS,IT,IU)+FOUR*&
          &XM(0,0,4,IS,IT,IU)-W*XN(0,0,1,IS,IT,IU)-FOUR*W*&
          &XM(0,0,3,IS,IT,IU)+TWO*W*W*XN(0,0,2,IS,IT,IU))-FOUR*&
          &W**3*XN(0,0,1,IS,IT,IU))
        BETNEW=BETNEW+BMAT(II1)*TWO*THREE*(XN(2,0,0,IS,IT,IU)+FOUR*&
          &XM(4,0,0,IS,IT,IU))
        CL1NEW=CL1NEW+BMAT(II1)*(XN(0,0,1,IS,IT,IU)/TWO+TWO*&
          &(XM(2,0,1,IS,IT,IU)-W*XM(2,0,0,IS,IT,IU)))
        CL2NEW=CL2NEW+BMAT(II1)*(XN(2,0,0,IS,IT,IU)+&
          &XN(0,0,2,IS,IT,IU)+FOUR*XM(2,0,2,IS,IT,IU)-W*&
          &XN(0,0,1,IS,IT,IU)-FOUR*W*XM(2,0,1,IS,IT,IU)+TWO*W*W*&
          &XN(2,0,0,IS,IT,IU))
      END DO IUloop3
      END DO ITloop3
      END DO ISloop3

      !The absolute values in the following section should not be
      !  necessary, but are used to prevent terrible convergence from
      !  giving negative values of V1 that stop the program.
      ALLNEW=ALLNEW*SQRT(ABS(TWO/V1))**3
      BELNEW=BELNEW*(TWO/V1)**2
      BETNEW=BETNEW*(TWO/V2)**2
      CL1NEW=CL1NEW*SQRT(ABS(TWO/V1))*(TWO/V2)
      CL2NEW=CL2NEW*(TWO/V1)*(TWO/V2)
      V1=V1*TPARA
      V2=V2*TPERP

      !Skip all of the diffusion coefficient calculations in first
      !  approximation.
      IF(IAPP==1)THEN
        FORMAT(1:6)='(F9.4,'
        IF(ROK>FIVE*TWO)THEN
          FORMAT(7:11)='F8.4,'
        ELSE
          FORMAT(7:11)='F8.5,'
        ENDIF
        IF(V1*KELVIN<1.D4.AND.V2*KELVIN<1.D4)THEN
          FORMAT(12:)='2F10.3,5F7.4,44X,I2)'
        ELSE
          FORMAT(12:)='2F10.2,5F7.4,44X,I2)'
        ENDIF
        WRITE(14,FORMAT)AC,ROK,V1*KELVIN,V2*&
          &KELVIN,1+ALLNEW,BELNEW,BETNEW,1+CL1NEW,CL2NEW,IAPP
        CYCLE IAPPloop
      ENDIF

      !We next set up the density-gradient-dependent moment equations
      !  parallel to the electric field.  We again use AMAT for the
      !  matrix elements, but use CMAT for the inhomogenous terms in
      !  the matrix equations.  Then call SIMQE to solve for the moments
      !  needed for the diffusion coefficient parallel to the field.
      II1=0
      IPloop4: DO IP=0,IAPP,2
      IQloop4: DO IQ=0,IAPP-IP,2
      IRloop4: DO IR=0,IAPP-IP-IQ
        IF(IP+IQ+IR==0)CYCLE IRloop4
        II1=II1+1
        CMAT(II1)=DBLE(2*IR+2)*XN(IP,IQ,IR+1,0,0,0)+&
          &XN(IP,IQ,IR-1,0,0,0)-TWO*W*XN(IP,IQ,IR,0,0,0)
        JJ1=0
        JJ2=0
        ISloop4: DO IS=0,IAPP+1,2
        ITloop4: DO IT=0,IAPP+1-IS,2
        IUloop4: DO IU=0,IAPP+1-IS-IT
          IF(IS+IT+IU==0)CYCLE IUloop4
          JJ1=JJ1+1
          CMAT(II1)=CMAT(II1)+BMAT(JJ1)*(DBLE(2*IR+2)*&
            &XM(IP,IQ,IR+1,IS,IT,IU)+XM(IP,IQ,IR-1,IS,IT,IU)-TWO*W*&
            &XM(IP,IQ,IR,IS,IT,IU))
          IF(IS<=IAPP.AND.IT<=IAPP-IS.AND.IU<=IAPP-IS-IT)THEN
            JJ2=JJ2+1
            AMAT(II1,JJ2)=ELEM0(IP,IQ,IR,IS,IT,IU)/BB-&
              &ELEM0(IP,IQ,IR,0,0,0)/BB*XN(0,0,0,IS,IT,IU)/&
              &TWO**(IS+IT+IU)/H(IS+1)/H(IT+1)/H(IU+1)-&
              &XM(IP,IQ,IR-1,IS,IT,IU)*TWO**(IP+IQ+IR)/&
              &TWO**(1+IS+IT+IU)*H(IP+1)/H(IS+1)*H(IQ+1)/H(IT+1)*&
              &H(IR+1)/H(IU+1)
          ENDIF
        END DO IUloop4
        END DO ITloop4
        END DO ISloop4
        CMAT(II1)=CMAT(II1)*TWO**(IP+IQ+IR)*H(IP+1)*H(IQ+1)*H(IR+1)
      END DO IRloop4
      END DO IQloop4
      END DO IPloop4
      CALL SIMQE(AMAT,MAXEQ,CMAT,II1,KK)
      IF(KK==1)THEN
         WRITE(14,*)' Error in SIMQE for parallel diffusion &
           &for T='//FNUM
         RETURN
      END IF

      !Parallel to the electric field, calculate EINS1 and DIFF1.
      EINS1=ZERO
      II1=0
      ISloop5: DO IS=0,IAPP,2
      ITloop5: DO IT=0,IAPP-IS,2
      IUloop5: DO IU=0,IAPP-IS-IT
        IF(IS+IT+IU==0)CYCLE IUloop5
        II1=II1+1
        IF(CMAT(II1)==ZERO)CYCLE IUloop5
        CMAT(II1)=CMAT(II1)/TWO**(IS+IT+IU)/H(IS+1)/H(IT+1)/H(IU+1)
        EINS1=EINS1+CMAT(II1)*XN(0,0,1,IS,IT,IU)
      END DO IUloop5
      END DO ITloop5
      END DO ISloop5
      EINS1=EINS1/TWO/SQRT(RAM*TEMP/TPARA)*TPARA/TEMP/VDSN
      DIFF1=EINS1*RMOB*TEMP/AC/ICHARG*KELVIN*&
        &IDEALGAS/FARADAY/FIVE/TWO

      !We next set up the density-gradient-dependent moment equations
      !  perpendicular to the electric field.  We again use AMAT for the 
      !  matrix elements and CMAT for the inhomogenous terms in the
      !  matrix equations.  Then we call SIMQE to solve for the moments
      !  needed for the diffusion coefficient perpendicular to the field.
      II1=0
      IPloop6: DO IP=1,IAPP,2
      IQloop6: DO IQ=0,IAPP-IP,2
      IRloop6: DO IR=0,IAPP-IP-IQ
        II1=II1+1
        CMAT(II1)=DBLE(2*IP+2)*XN(IP+1,IQ,IR,0,0,0)+&
          &XN(IP-1,IQ,IR,0,0,0)
      JJ1=0
      ISloop6: DO IS=0,IAPP+1,2
      ITloop6: DO IT=0,IAPP+1-IS,2
      IUloop6: DO IU=0,IAPP+1-IS-IT
        IF(IS+IT+IU==0)CYCLE IUloop6
        JJ1=JJ1+1
        CMAT(II1)=CMAT(II1)+BMAT(JJ1)*(DBLE(2*IP+2)*&
          &XM(IP+1,IQ,IR,IS,IT,IU)+XM(IP-1,IQ,IR,IS,IT,IU))
      END DO IUloop6
      END DO ITloop6
      END DO ISloop6
      CMAT(II1)=CMAT(II1)*TWO**(IP+IQ+IR)*H(IP+1)*H(IQ+1)*H(IR+1)
      JJ2=0
      ISloop7: DO IS=1,IAPP,2
      ITloop7: DO IT=0,IAPP-IS,2
      IUloop7: DO IU=0,IAPP-IS-IT
        JJ2=JJ2+1
        AMAT(II1,JJ2)=ELEM0(IP,IQ,IR,IS,IT,IU)/BB-&
          &ELEM0(IP,IQ,IR,0,0,0)/BB*XN(0,0,0,IS,IT,IU)/&
          &TWO**(IS+IT+IU)/H(IS+1)/H(IT+1)/H(IU+1)-&
          &XM(IP,IQ,IR-1,IS,IT,IU)*TWO**(IP+IQ+IR)/&
          &TWO**(1+IS+IT+IU)*H(IP+1)/H(IS+1)*H(IQ+1)/H(IT+1)*&
          &H(IR+1)/H(IU+1)
      END DO IUloop7
      END DO ITloop7
      END DO ISloop7
      END DO IRloop6
      END DO IQloop6
      END DO IPloop6
      CALL SIMQE(AMAT,MAXEQ,CMAT,II1,KK)
      IF(KK==1)THEN
        WRITE(14,*)' Error in SIMQE for perpendicular diffusion &
          &for T='//FNUM
        RETURN
      END IF

      !Perpendicular to the electric field, calculate EINS2 and DIFF2.
      EINS2=ZERO
      II1=0
      ISloop8: DO IS=1,IAPP,2
      ITloop8: DO IT=0,IAPP-IS,2
      IUloop8: DO IU=0,IAPP-IS-IT
        II1=II1+1
        IF(CMAT(II1)==ZERO)CYCLE IUloop8
        CMAT(II1)=CMAT(II1)/TWO**(IS+IT+IU)/H(IS+1)/H(IT+1)/H(IU+1)
        EINS2=EINS2+CMAT(II1)*XN(1,0,0,IS,IT,IU)
      END DO IUloop8
      END DO ITloop8
      END DO ISloop8
      EINS2=EINS2/TWO/SQRT(RAM*TEMP/TPARA)*TPERP/TEMP/VDSN
      DIFF2=EINS2*RMOB*TEMP/AC/ICHARG*KELVIN*&
        &IDEALGAS/FARADAY/FIVE/TWO

      !Write the results to the long output file.
      FORMAT(1:6)='(F9.4,'
      IF(ROK>FIVE*TWO)THEN
        FORMAT(7:11)='F8.4,'
      ELSE
        FORMAT(7:11)='F8.5,'
      ENDIF
      IF(V1*KELVIN<1.D5.AND.V2*KELVIN<1.D5)THEN
        FORMAT(12:24)='2F10.3,5F7.4,'
      ELSE
        FORMAT(12:24)='2F10.2,5F7.4,'
      ENDIF
      IF(DIFF1<1.D5.AND.EINS1<1.D5.AND.DIFF2<1.D5.AND.EINS2<1.D5)THEN
        FORMAT(25: )='4F11.5,I2)'
      ELSE
        FORMAT(25: )='4F11.4,I2)'
      ENDIF
      WRITE(14,FORMAT)AC,ROK,V1*KELVIN,V2*KELVIN,&
        &1+ALLNEW,BELNEW,BETNEW,1+CL1NEW,CL2NEW,DIFF1,EINS1,DIFF2,&
        &EINS2,IAPP

      !IF IAPP equals 2, establish the old values to which higher
      !  approximations will be compared.
      IF(IAPP==2)THEN
        ROKOLD=ROK
        V1OLD=V1
        V2OLD=V2
        ALLOLD=ALLNEW
        BELOLD=BELNEW
        BETOLD=BETNEW
        CL1OLD=CL1NEW
        CL2OLD=CL2NEW
        D1OLD=DIFF1
        D2OLD=DIFF2
      !If IAPP is greater than 2, compute the fractional changes
      !  between this approximation and the previous approximation.
      ELSEIF(IAPP>2)THEN
        ACCN1=ABS(ONE-ROKOLD/ROK)
        ACCN2=MAX(ABS(ONE-D1OLD/DIFF1),ABS(ONE-D2OLD/DIFF2))
        ACCN3=MAX(ABS(ONE-V1OLD/V1),ABS(ONE-V2OLD/V2),&
          &ABS(ONE-(ONE+ALLOLD)/(ONE+ALLNEW)),ABS(ONE-BELOLD/BELNEW),&
          &ABS(ONE-BETOLD/BETNEW),ABS(ONE-(ONE+CL1OLD)/(ONE+CL1NEW)),&
          &ABS(ONE-CL2OLD/CL2NEW))
        !Check whether convergence has been achieved.  If so, write the
        !  results to the summary file 'transport'.  If not, save the
        !  present values and go the next higher approximation.
        IF((MOBERR.AND.ACCN1<=ACCT1).OR.(ACCN1<=ACCT1.AND.ACCN2<=ACCT2&
          &.AND.ACCN3<=ACCT3.AND.IAPP>=IAPPN))THEN
          FORMAT(1:6)='(F7.2,'
          IF(ROK>FIVE*TWO)THEN
            FORMAT(7:11)='F7.3,'
          ELSE
            FORMAT(7:11)='F7.4,'
          ENDIF
          IF(V1*KELVIN<1.D5.AND.V2*KELVIN<1.D5)THEN
            FORMAT(12:17)='2F8.1,'
          ELSE
            FORMAT(12:17)='2F8.0,'
          ENDIF
          IF(DIFF1<1.D3.AND.DIFF2<1.D3)THEN
            FORMAT(18: )='2F8.3,5F6.3)'
          ELSE
            FORMAT(18: )='2F8.2,5F6.3)'
          ENDIF
          OPEN(UNIT=15,FILE='junk.'//FNUM,POSITION='APPEND')
          WRITE(15,FORMAT)AC,ROK,V1*KELVIN,V2*&
            &KELVIN,DIFF1,DIFF2,ONE+ALLNEW,BELNEW,BETNEW,&
            &ONE+CL1NEW,CL2NEW
          CLOSE(15)
          EXIT IAPPloop
        ELSE
          ROKOLD=ROK
          V1OLD=V1
          V2OLD=V2
          ALLOLD=ALLNEW
          BELOLD=BELNEW
          BETOLD=BETNEW
          CL1OLD=CL1NEW
          CL2OLD=CL2NEW
          D1OLD=DIFF1
          D2OLD=DIFF2
        END IF
      END IF
    END DO IAPPloop

    IF(IAPP>IAPPM)THEN
      !If convergence is not achieved, print best estimate if desired.
      IF(PRINTALL)THEN
        FORMAT(1:6)='(F7.2,'
        IF(ROK>FIVE*TWO)THEN
          FORMAT(7:11)='F7.3,'
        ELSE
          FORMAT(7:11)='F7.4,'
        ENDIF
          IF(V1*KELVIN<1.D5.AND.V2*KELVIN<1.D5)THEN
          FORMAT(12:17)='2F8.1,'
        ELSE
          FORMAT(12:17)='2F8.0,'
        ENDIF
        IF(DIFF1<1.D3.AND.DIFF2<1.D3)THEN
          FORMAT(18: )='2F8.3,5F6.3,A1)'
        ELSE
          FORMAT(18: )='2F8.2,5F6.3,A1)'
        ENDIF
        OPEN(UNIT=15,FILE='junk.'//FNUM,POSITION='APPEND')
        WRITE(15,FORMAT)AC,ROK,V1*KELVIN,V2*KELVIN,&
          &DIFF1,DIFF2,ONE+ALLNEW,BELNEW,BETNEW,ONE+CL1NEW,CL2NEW,'*'
        CLOSE(15)
      ENDIF
      !If convergence is not achieved, check STANDERR to see if
      !  we want to change ACCT1, ACCT2 and ACCT3
      IF(STANDERR)THEN
        IF(ACCT1<0.00495D0)THEN
          ACCT1=0.005D0
          ACCT2=0.08D0
          ACCT3=0.05D0
        ELSEIF(ACCT1<0.010D0)THEN
          ACCT1=0.01D0
          ACCT2=0.15D0
          ACCT3=0.10D0
        ELSEIF(ACCT1<0.015D0)THEN
          ACCT1=0.015D0
          ACCT2=0.25D0
          ACCT3=0.15D0
        ELSEIF(ACCT1<0.020D0)THEN
          ACCT1=0.020D0
          ACCT2=0.50D0
          ACCT3=0.20D0
          OPEN(UNIT=15,FILE='junk.'//FNUM,POSITION='APPEND')
          WRITE(15,'(A,3(F5.3,1X))')'*****Accuracies changed to ',&
            &REAL(ACCT1),REAL(ACCT2),REAL(ACCT3)
          CLOSE(15)
        ELSEIF(HALT)THEN
          EXIT
        ENDIF
        IF(ACCT1<0.020D0)THEN
          OPEN(UNIT=15,FILE='junk.'//FNUM,POSITION='APPEND')
          WRITE(15,'(A,3(F5.3,1X))')'*****Accuracies changed to ',&
            &REAL(ACCT1),REAL(ACCT2),REAL(ACCT3)
          CLOSE(15)
          CYCLE VDSloop
        ENDIF
      !If convergence is not achieved, halt if desired.
      ELSEIF(HALT)THEN
        EXIT
      ENDIF
    ENDIF
    IF(VDS<0.01D0)THEN
      VDS=VDS*4.00D0
    ELSEIF(VDS<0.04D0)THEN
      VDS=VDS*1.75D0
    ELSEIF(VDS<0.07D0)THEN
      VDS=VDS*1.40D0
    ELSEIF(VDS<0.13D0)THEN
      VDS=VDS*1.20D0
    ELSEIF(VDS<0.22D0)THEN
      VDS=VDS*1.15D0
    ELSEIF(VDS<0.29D0)THEN
      VDS=VDS*1.10D0
    ELSEIF(VDS<0.5D0)THEN
      VDS=VDS*1.05D0
    ELSE
      VDS=VDS*1.02D0
    END IF

    IF(VDS>VDEND)EXIT VDSloop
  ENDDO VDSloop

  IF(STANDERR)THEN
    OPEN(UNIT=15,FILE='junk.'//FNUM)
    VAL1=0.0D0
    VAL2=0.0D0
    VAL3=0.0D0
    VAL4=0.0D0
    DO I=1,11
      READ(15,'(A)')COMM
    END DO
    K=0
    DO 
      READ(15,'(A)',IOSTAT=ERR)COMM
      IF(ERR/=0)EXIT
      IF(COMM(1:1)=='*')THEN
        BACKSPACE(15)
        BACKSPACE(15)
        IF(K==0)THEN
          READ(15,*)VAL1
          K=1
        ELSEIF(K==1)THEN
          READ(15,*,IOSTAT=ERR)VAL2
          IF(ERR/=0)VAL2=0.0D0
          K=2
        ELSEIF(K==2)THEN
          READ(15,*,IOSTAT=ERR)VAL3
          IF(ERR/=0)VAL3=0.0D0
          K=3
        ELSEIF(K==3)THEN
          READ(15,*,IOSTAT=ERR)VAL4
          IF(ERR/=0)VAL4=0.0D0
          K=4
          EXIT
        ENDIF
        READ(15,'(A)')COMM
      END IF
    END DO
    REWIND(15)
    OPEN(UNIT=25,FILE='transport.'//FNUM)
    DO I=1,5
      READ(15,'(A)')COMM
      WRITE(25,'(A)')TRIM(COMM)
    ENDDO
    IF(K==0)THEN
      DO I=1,4
        READ(15,'(A)')COMM
        WRITE(25,'(A)')TRIM(COMM)
      ENDDO
    ELSEIF(VAL1<=0.0D0)THEN
      DO I=1,3
        READ(15,'(A)')COMM
        WRITE(25,'(A)')TRIM(COMM)
      ENDDO
      READ(15,'(A)')COMM
      WRITE(25,'(A)')'Warning: VAL1<=0.  Check file junk.'//FNUM
    ELSEIF(K==1)THEN
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A,I3,A)')'Fractional accuracy for Ko is ',&
        &100.0*ACCT1,'% below and 0.5% above ',INT(VAL1+0.99),' Td.'
      READ(15,'(A29,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For NDx and NDz, it is ',100.0*ACCT1,&
        &'% and 8%.'
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For the other quantities, it is ',&
        &100.0*ACCT1,'% and 5%.'
      READ(15,'(A)')COMM
      WRITE(25,'(A)')' '
    ELSEIF(K==2.AND.VAL2<=0.0D0)THEN
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A,I3,A)')'Fractional accuracy for Ko is ',&
        &100.0*ACCT1,'% below and 1.0% above ',INT(VAL1+0.99),' Td.'
      READ(15,'(A29,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For NDx and NDz, it is ',100.0*ACCT1,&
        &'% and 15%.'
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For the other quantities, it is ',&
        &100.0*ACCT1,'% and 10%.'
      READ(15,'(A)')COMM
      WRITE(25,'(A)')' '
    ELSEIF(K==2)THEN
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,3(A,I3),A)')'Fractional accuracy for Ko &
        &is ',100.0*ACCT1,'% below ',INT(VAL1+0.99),' Td, 0.5% for ',&
        &INT(VAL1+0.99),'-',INT(VAL2+0.99),' Td,'
      WRITE(25,'(5X,A,I3,A)')'and 1.0% above ',INT(VAL2+0.99),' Td.'
      READ(15,'(A29,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For NDx and NDz, it is ',100.0*ACCT1,&
        &'%, 8% and 15%.'
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For the other quantities, it is ',&
        &100.0*ACCT1,'%, 5% and 10%.'
      READ(15,'(A)')COMM
    ELSEIF(K==3.AND.VAL2<=0.0D0.AND.VAL3<=0.0D0)THEN
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A,I3,A)')'Fractional accuracy for Ko is ',&
        &100.0*ACCT1,'% below and 1.5% above ',INT(VAL1+0.99),' Td.'
      READ(15,'(A29,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For NDx and NDz, it is ',100.0*ACCT1,&
        &'% and 25%.'
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For the other quantities, it is ',&
        &100.0*ACCT1,'% and 15%.'
      READ(15,'(A)')COMM
      WRITE(25,'(A)')' '
    ELSEIF(K==3.AND.VAL2<=0.0D0.AND.VAL3>0.0D0)THEN
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,3(A,I3),A)')'Fractional accuracy for Ko &
        &is ',100.0*ACCT1,'% below ',INT(VAL1+0.99),' Td, 1.0% for ',&
        &INT(VAL1+0.99),'-',INT(VAL3+0.99),' Td,'
      WRITE(25,'(5X,A,I3,A)')'and 1.5% above ',INT(VAL3+0.99),' Td.'
      READ(15,'(A29,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For NDx and NDz, it is ',100.0*ACCT1,&
        &'%, 15% and 25%.'
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For the other quantities, it is ',&
        &100.0*ACCT1,'%, 10% and 15%.'
      READ(15,'(A)')COMM
    ELSEIF(K==3.AND.VAL2>0.0D0.AND.VAL3<=0.0D0)THEN
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,3(A,I3),A)')'Fractional accuracy for Ko &
        &is ',100.0*ACCT1,'% below ',INT(VAL1+0.99),' Td, 0.5% for ',&
        &INT(VAL1+0.99),'-',INT(VAL2+0.99),' Td,'
      WRITE(25,'(5X,A,I3,A)')'and 1.5% above ',INT(VAL2+0.99),' Td.'
      READ(15,'(A29,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For NDx and NDz, it is ',100.0*ACCT1,&
        &'%, 8% and 25%.'
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For the other quantities, it is ',&
        &100.0*ACCT1,'%, 5% and 15%.'
      READ(15,'(A)')COMM
    ELSEIF(K==3)THEN
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,3(A,I3),A)')'Fractional accuracy for Ko &
        &is ',100.0*ACCT1,'% below ',INT(VAL1+0.99),' Td, 0.5% for ',&
        &INT(VAL1+0.99),'-',INT(VAL2+0.99),' Td,'
      WRITE(25,'(5X,3(A,I3),A)')'1.0% for ',INT(VAL2+0.99),'-',&
        &INT(VAL3+0.99),' Td, and 1.5% above ',INT(VAL3+0.99),' Td.'
      READ(15,'(A29,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For NDx and NDz, it is ',100.0*ACCT1,&
        &'%, 8%, 15% and 25%.'
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For the other quantities, it is ',&
        &100.0*ACCT1,'%, 5%, 10% and 15%.'
    ELSEIF(K==4.AND.VAL2<=0.0D0.AND.VAL3<=0.0D0.AND.VAL4<=0.0D0)THEN
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A,I3,A)')'Fractional accuracy for Ko is ',&
        &100.0*ACCT1,'% below and 2.0% above ',INT(VAL1+0.99),' Td.'
      READ(15,'(A29,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For NDx and NDz, it is ',100.0*ACCT1,&
        &'% and 50%.'
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For the other quantities, it is ',&
        &100.0*ACCT1,'% and 20%.'
      READ(15,'(A)')COMM
      WRITE(25,'(A)')' '
    ELSEIF(K==4.AND.VAL2<=0.0D0.AND.VAL3<=0.0D0.AND.VAL4>0.0D0)THEN
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,3(A,I3),A)')'Fractional accuracy for Ko &
        &is ',100.0*ACCT1,'% below ',INT(VAL1+0.99),' Td, 1.5% for ',&
        &INT(VAL1+0.99),'-',INT(VAL4+0.99),' Td,'
      WRITE(25,'(5X,A,I3,A)')'and 2.0% above ',INT(VAL4+0.99),' Td.'
      READ(15,'(A29,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For NDx and NDz, it is ',100.0*ACCT1,&
        &'%, 25% and 50%.'
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For the other quantities, it is ',&
        &100.0*ACCT1,'%, 15% and 20%.'
      READ(15,'(A)')COMM
    ELSEIF(K==4.AND.VAL2<=0.0D0.AND.VAL3>0.0D0.AND.VAL4<=0.0D0)THEN
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,3(A,I3),A)')'Fractional accuracy for Ko &
        &is ',100.0*ACCT1,'% below ',INT(VAL1+0.99),' Td, 1.0% for ',&
        &INT(VAL1+0.99),'-',INT(VAL3+0.99),' Td,'
      WRITE(25,'(5X,A,I3,A)')'and 2.0% above ',INT(VAL3+0.99),' Td.'
      READ(15,'(A29,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For NDx and NDz, it is ',100.0*ACCT1,&
        &'%, 15% and 50%%.'
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For the other quantities, it is ',&
        &100.0*ACCT1,'%, 10% and 20%.'
      READ(15,'(A)')COMM
    ELSEIF(K==4.AND.VAL2<=0.0D0.AND.VAL3>0.0D0.AND.VAL4>0.0D0)THEN
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,3(A,I3),A)')'Fractional accuracy for Ko &
        &is ',100.0*ACCT1,'% below ',INT(VAL1+0.99),' Td, 1.0% for ',&
        &INT(VAL1+0.99),'-',INT(VAL3+0.99),' Td,'
      WRITE(25,'(5X,3(A,I3),A)')'1.5% for ',INT(VAL3+0.99),'-',&
        &INT(VAL4+0.99),' Td and 2.0% above ',INT(VAL4+0.99),' Td.'
      READ(15,'(A29,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For NDx and NDz, it is ',100.0*ACCT1,&
        &'%, 15%, 25% and 50%%.'
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For the other quantities, it is ',&
        &100.0*ACCT1,'%, 10%, 15% and 20%.'
      READ(15,'(A)')COMM
    ELSEIF(K==4.AND.VAL2>0.0D0.AND.VAL3<=0.0D0.AND.VAL4<=0.0D0)THEN
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,3(A,I3),A)')'Fractional accuracy for Ko &
        &is ',100.0*ACCT1,'% below ',INT(VAL1+0.99),' Td, 0.5% for ',&
        &INT(VAL1+0.99),'-',INT(VAL2+0.99),' Td,'
      WRITE(25,'(5X,A,I3,A)')'and 2.0% above ',INT(VAL2+0.99),' Td.'
      READ(15,'(A29,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For NDx and NDz, it is ',100.0*ACCT1,&
        &'%, 8% and 50%.'
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For the other quantities, it is ',&
        &100.0*ACCT1,'%, 5% and 20%.'
      READ(15,'(A)')COMM
    ELSEIF(K==4.AND.VAL2>0.0D0.AND.VAL3<=0.0D0.AND.VAL4>0.0D0)THEN
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,3(A,I3),A)')'Fractional accuracy for Ko &
        &is ',100.0*ACCT1,'% below ',INT(VAL1+0.99),' Td, 0.5% for ',&
        &INT(VAL1+0.99),'-',INT(VAL2+0.99),' Td,'
      WRITE(25,'(5X,3(A,I3),A)')'1.5% for ',INT(VAL2+0.99),'-',&
        &INT(VAL4+0.99),' Td and 2.0% above ',INT(VAL4+0.99),' Td.'
      READ(15,'(A29,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For NDx and NDz, it is ',100.0*ACCT1,&
        &'%, 8%, 25% and 50%.'
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For the other quantities, it is ',&
        &100.0*ACCT1,'%, 5%, 15% and 20%.'
      READ(15,'(A)')COMM
    ELSEIF(K==4.AND.VAL2>0.0D0.AND.VAL3>0.0D0.AND.VAL4<=0.0D0)THEN
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,3(A,I3),A)')'Fractional accuracy for Ko &
        &is ',100.0*ACCT1,'% below ',INT(VAL1+0.99),' Td, 0.5% for ',&
        &INT(VAL1+0.99),'-',INT(VAL2+0.99),' Td,'
      WRITE(25,'(5X,3(A,I3),A)')'1.0% for ',INT(VAL2+0.99),'-',&
        &INT(VAL3+0.99),' Td and 2.0% above ',INT(VAL3+0.99),' Td.'
      READ(15,'(A29,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For NDx and NDz, it is ',100.0*ACCT1,&
        &'%, 8%, 15% and 50%.'
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For the other quantities, it is ',&
        &100.0*ACCT1,'%, 5%, 10% and 20%.'
      READ(15,'(A)')COMM
    ELSEIF(K==4)THEN
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,3(A,I3),A)')'Fractional accuracy for Ko &
        &is ',100.0*ACCT1,'% below ',INT(VAL1+0.99),' Td, 0.5% for ',&
        &INT(VAL1+0.99),'-',INT(VAL2+0.99),' Td,'
      WRITE(25,'(5X,5(A,I3),A)')'1.0% for ',INT(VAL2+0.99),'-',&
        &INT(VAL3+0.99),' Td, 1.5% for ',INT(VAL3+0.99),'-',&
        &INT(VAL4+0.99),' Td, and 2.0% above ',INT(VAL4+0.99),' Td.'
      READ(15,'(A29,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For NDx and NDz, it is ',100.0*ACCT1,&
        &'%, 8%, 15%, 25% and 50%.'
      READ(15,'(A42,F6.4)')COMM,ACCT1
      WRITE(25,'(5X,A,F4.2,A)')'For the other quantities, it is ',&
        &100.0*ACCT1,'%, 5%, 10%, 15% and 20%.'
      READ(15,'(A)')COMM
    ENDIF
    READ(15,'(A)')COMM
    WRITE(25,'(A)')TRIM(COMM)
    READ(15,'(A)')COMM
    WRITE(25,'(A)')TRIM(COMM)
  ELSE
    OPEN(UNIT=15,FILE='junk.'//FNUM)
    CLOSE(25)
    OPEN(UNIT=25,FILE='transport.'//FNUM)
    DO I=1,11
      READ(15,'(A)')COMM
      WRITE(25,'(A)')TRIM(COMM)
    ENDDO
  END IF

  CLOSE(25)
  OPEN(UNIT=25,FILE='transport.'//FNUM)
  OPEN(UNIT=35,FILE='t.'//FNUM)
  DO I=1,11
    READ(25,'(A)')COMM
    WRITE(35,'(A)')TRIM(COMM)
  ENDDO
  CLOSE(25)
  OPEN(UNIT=25,FILE='transport.'//FNUM,POSITION='APPEND')

  DO 
    READ(15,'(A)',IOSTAT=ERR)COMM
    IF(ERR/=0)EXIT
    IF(COMM(1:1)/='*')THEN
      WRITE(25,'(A)')TRIM(COMM)
      WRITE(35,'(A)')TRIM(ADJUSTL(COMM(1:7)))//','//&
        &TRIM(ADJUSTL(COMM(8:14)))//','//TRIM(ADJUSTL(COMM(15:22)))&
        &//','//TRIM(ADJUSTL(COMM(23:30)))//','//&
        &TRIM(ADJUSTL(COMM(31:38)))//','//&
        &TRIM(ADJUSTL(COMM(39:46)))//','//TRIM(ADJUSTL(COMM(47:52)))//&
        &','//TRIM(ADJUSTL(COMM(53:58)))//','//&
        &TRIM(ADJUSTL(COMM(59:64)))//','//TRIM(ADJUSTL(COMM(65:70)))//&
        &','//TRIM(ADJUSTL(COMM(71:76)))
    ENDIF
  ENDDO
  CLOSE(25)
  CLOSE(35)
END SUBROUTINE GRAMCH

SUBROUTINE TERP(ENER,Q,IERR)
  !This program interpolates in array FILE2 to obtain the values of Q
  !  appropriate to energy ENER. IERR=1 if ENER is below the minimum
  !  and IERR=2 if above the maximum.
  USE DATA
  IMPLICIT NONE

  !Declare the arguments.
  !ENER is the energy in hartree.
  !Q is the cross section in square bohr and other units (see
  !  function F3).  All of the units are taken care of when the
  !  transport coefficients are determined just prior to being
  !  printed.
  !IERR is a flag to indicate when results are extrapolated either to
  !  energies below the minimum (IERR=1) or above the maximum (IERR=2)
  !  used in program PROGCROS.
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::ENER
  REAL(SELECTED_REAL_KIND(18)),INTENT(OUT)::Q
  INTEGER,INTENT(OUT)::IERR

  !Declare the common variables.
  INTEGER::IPTRAN,IQTRAN,ILTRAN
  COMMON/INT1/ILTRAN,IPTRAN,IQTRAN

  !Declare the local variables.
  !EMIN and EMAX are local copies of the minimum and maximum energies.
  !EA, EB and EL are used to interpolate energies logarithmically.
  !A, B and C are simply used for temporary storage.
  !N, J, JJ and JK are used to interpolate energies and cross sections.
  REAL(SELECTED_REAL_KIND(18))::EMIN,EMAX,EA,EB,EL,A,B,C
  INTEGER::N,J,JJ,JK

  !Initialize IERRR and locate the energy region appropriate to ENER.
  IF(ENER<=ZERO)THEN
     PRINT *,'Error: subroutine TERP was called with ENER=',ENER
     STOP
  ENDIF
  IERR=0
  DO
    IF(N1/=0)THEN
      N=N1
      J=1
      EMIN=EMIN1
      IF(ENER<EMIN)IERR=1
      EMAX=EMIN2
      IF(ENER<=EMAX)EXIT
    END IF
    IF(N2/=0)THEN
      N=N2
      J=2
      EMIN=EMIN2
      EMAX=EMIN3
      IF(ENER<=EMAX)EXIT
    END IF
    IF(N3==0)THEN
      IERR=2
      EXIT
    ELSE
      EMIN=EMIN3
      EMAX=EMIN4
      IF(ENER>EMAX)IERR=2
      N=N3
      J=3
      EXIT
    END IF
  END DO

  !Interpolate to find Q(ENER) by means of the Chebychev polynomial
  !  fitting values tabulated in array FILE2.  Note that this is a fit 
  !  of ln(Q) versus ln(ENER).
  EA=LOG(EMIN)
  EB=LOG(EMAX)
  EL=LOG(ENER)
  EL=TWO*(EL-EA)/(EB-EA)-ONE
  A=ZERO
  JJ=MAXE*(J-1)+N
  B=FILE2(JJ,ILTRAN)/TWO
  DO JK=2,N
    JJ=MAXE*(J-1)+N+1-JK
    C=TWO*EL*B-A+FILE2(JJ,ILTRAN)
    A=B
    B=C
  END DO
  Q=EXP(B-EL*A-FILE2(JJ,ILTRAN)/TWO)
END SUBROUTINE TERP

SUBROUTINE F1(X,ANS)
  !This subroutine guides the first integral in eq. A.15 (not A.18) of
  !  L. A. Viehland, Chem. Phys. 179 (1994) 71.
  USE DATA
  IMPLICIT NONE

  !Declare the arguments.
  !X is the square of gamma_r in eq. A.15.
  !ANS is the value of the integrand returned to the calling program.
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::X
  REAL(SELECTED_REAL_KIND(18)),INTENT(OUT)::ANS

  !Declare the common variables.
  REAL(SELECTED_REAL_KIND(18))::ACCT,ACCT1,ACCT2,ACCT3,XTRANS
  COMMON/REAL0/ACCT,ACCT1,ACCT2,ACCT3
  COMMON/REAL5/XTRANS

  !Declare the local variables.
  !A and B are the endpoints for the second integral.
  !EA and EB are the values of the integrands at the endpoints.
  !BNS is the value obtained by performing the second integral.
  !ACC is the accuracy of the integral.
  !IPR is a flag for convergence failure, which we ignore.
  REAL(SELECTED_REAL_KIND(18))::A,B,EA,EB,BNS,ACC
  INTEGER::IPR
  EXTERNAL F2

  !Put X into XTRANS so that it can be transferred to F2.
  XTRANS=X

  !Start the second integral.
  A=ZERO
  B=ONE
  CALL F2(A,EA)
  CALL F2(B,EB)
  ACC=1.1D0*ACCT
  CALL INTEG(F2,A,B,EA,EB,BNS,ACC,IPR)
  ANS=BNS
END SUBROUTINE F1

SUBROUTINE F2(Y,ANS)
  !This subroutine guides the second integral in eq. A.15 (not A.18) of
  !  L. A. Viehland, Chem. Phys. 179 (1994) 71.
  USE DATA
  IMPLICIT NONE

  !Declare the arguments.
  !Y is the integration variable v_d+gamma_z in eq. A.16.
  !ANS is the value of the integrand returned to the calling program.
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::Y
  REAL(SELECTED_REAL_KIND(18)),INTENT(OUT)::ANS

  !Declare the common variables.
  REAL(SELECTED_REAL_KIND(18))::XTRANS
  COMMON/REAL5/XTRANS

  !Declare the local variables.
  !F3 is the function that evaluates the integrand for both integrals.
  REAL(SELECTED_REAL_KIND(18))::F3

  !The function is never called with XTRANS=0, but it is called with
  !  Y=0 so we must be careful.
  IF(ABS(Y)<1.0D-20)THEN
    ANS=F3(XTRANS,Y)+F3(ONE/XTRANS,Y)/XTRANS**2
  ELSE
    ANS=F3(XTRANS,Y)+F3(ONE/XTRANS,Y)/XTRANS**2+F3(XTRANS,ONE/Y)/Y**2+&
      &F3(ONE/XTRANS,ONE/Y)/(XTRANS*Y)**2
  END IF
END SUBROUTINE F2

FUNCTION F3(X,Y)
  !This is the function that evaluates the integrand for the double
  !  integral in eq. A.15 of L. A. Viehland, Chem. Phys. 179 (1994) 71.
  USE DATA
  IMPLICIT NONE

  !Declare the arguments.
  !X and Y are the current values of the integration variables.
  !GG is the collision energy given by the bracketed term in eq. A.16.
  !Q is the cross section.
  !T1 and T2 are used to avoid recalculating terms.
  !IERR is a flag for extrapolation to get the cross section.
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::X,Y
  REAL(SELECTED_REAL_KIND(18))::F3,GG,Q,T1,T2
  INTEGER::IERR

  !Declare the common variables.
  REAL(SELECTED_REAL_KIND(18))::AM0,AM1,TEMP,RAM,TPARA,TPERP,VDBAR
  INTEGER::IPTRAN,IQTRAN,ILTRAN
  COMMON/REAL1/AM0,AM1,TEMP,RAM,TPARA,TPERP,VDBAR
  COMMON/INT1/ILTRAN,IPTRAN,IQTRAN

  !Determine T1, T2 and GG.
  T1=Y-VDBAR
  T2=Y+VDBAR
  GG=((RAM*TEMP+TPERP)*X+(RAM*TEMP+TPARA)*Y**2)/(RAM+ONE)

  !We must avoid extrapolating in the cross sections, since this can
  !  easily lead to divisions by zero or similar problems.  Here we
  !  simply assume constant cross sections beyond the energy limits,
  !  but this section of code could be modified to assume linear
  !  behavior for ln(Q) vs ln(E) in extrapolations.
  IF(GG<EMIN1)THEN
    CALL TERP(EMIN1,Q,IERR)
  ELSEIF(EMIN4<=ZERO)THEN
    PRINT *,'EMAX is too low.  See line 1006 of vary.f95'
    STOP
  ELSEIF(GG>EMIN4)THEN
    CALL TERP(EMIN4,Q,IERR)
  ELSE
    CALL TERP(GG,Q,IERR)
  ENDIF

  !Now multiply the cross section by the relative velocity.  A factor
  !  of 1000 is needed inside the square root if we want to be in basic
  !  SI units, but this is taken into account later, when the
  !  dimensional quantities are computed.
  Q=Q*SQRT(TWO*GG*(AM0+AM1)/AM0/AM1)

  !For the rest of the integrand, Be sure to avoid underflow (which is
  !  assumed to occur when T1 or T2 exceeds 26) or raising things to the
  !  zero power.
  F3=ZERO
  IF(IQTRAN==0)THEN
    IF(T1<DBLE(26))F3=EXP(-T1**2)
    IF(T2<DBLE(26))F3=F3+EXP(-T2**2)
  ELSE
    IF(T1<DBLE(26))F3=T1**IQTRAN*EXP(-T1**2)
    IF(T2<DBLE(26))THEN
      IF(MOD(IQTRAN,2)==0)THEN
        F3=F3+T2**IQTRAN*EXP(-T2**2)
      ELSE
        F3=F3-T2**IQTRAN*EXP(-T2**2)
      ENDIF
    ENDIF
  ENDIF
  IF(IPTRAN/=0)F3=F3*SQRT(X)**IPTRAN
  IF(X<26*26)THEN
    F3=F3*EXP(-X)/SQRT(PI)
  ELSE
    F3=ZERO
  ENDIF
  F3=F3*Q
END FUNCTION F3

RECURSIVE SUBROUTINE INTEG(F,A,B,EA,EB,ANS,ACC,IPR)
  !This subroutine uses the Clenshaw-Curtiss procedure to integrate
  !  the function given by external procedure F.
  !The lower and upper limits on the integral are A and B.  
  !In order to allow this procedure to be used when there are integrable
  !  singularities, the value of F at A must be supplied in EA and the
  !  value of F at B must be supplied in EB.
  !The relative accuracy of the calculated result in ANS is ACC.
  !If convergence is achieved with ANS nonzero, the subroutine returns
  !  with IPR=0.  If convergence is not achieved but the value of ANS is
  !  equal to 0 or very small, the subroutine returns with IPR=1.  Only
  !  if convergence is not achieved otherwise does it return with IPR=2.
  !This subroutine was modified in March, 2014, to incorporate extrapo-
  !  lations.
  USE DATA
  IMPLICIT NONE

  !Declare the dummy arguments explained above.
  EXTERNAL F
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::A,B,EA,EB,ACC
  REAL(SELECTED_REAL_KIND(18)),INTENT(OUT)::ANS
  INTEGER,INTENT(OUT)::IPR

  !Declare the local real arrays.
  !  FILE3 contains the function values used to calculate ANS.
  !  ERRMAX contains the maxima of the error estimates.
  !  PRANS holds previous estimates for ANS.
  !  EXTRAP holds estimates for ANS based on extrapolating from PRANS.
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(NPTSCH)::FILE3
  REAL(SELECTED_REAL_KIND(18))::ERRMAX
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(2)::PRANS,EXTRAP

  !Declare the local real variables.
  !  By a change of variables, ANS equals HALF times the integral from
  !    Y=-1 to Y=1 of F(AMID+HALF*Y).
  REAL(SELECTED_REAL_KIND(18))::AMID,HALF

  !Declare the local integer variables.
  !  JSTART locates the position of the start of the points and weights
  !    for the N-th order of approximation.
  !  I is a loop index.
  INTEGER::JSTART,N,I

  !Declare the local logical variables.
  !  CONVERGE is true if convergence is achieved.
  LOGICAL::CONVERGE

  !Determine AMID and HALF
  AMID=(B+A)/TWO
  HALF=(B-A)/TWO

  !Initialize FILE3, PRANS and EXTRAP.
  FILE3=ZERO
  PRANS=ZERO
  EXTRAP=ZERO

  !For integration of order 1, place the values at the endpoints in
  !positions 1 and 3 of FILE3.  Then determine the value at the center
  !of the integration range (where Y=0 and the integration variable is
  !equal to AMID) and put it in position 2.
  FILE3(1)=EA
  FILE3(3)=EB
  CALL F(AMID,FILE3(2))

  !Start the main calculation.  Setting JSTART=1 means we are starting
  !  at the first values in CHEPTS and CHEWTS, which corresponds to N=1
  !  in the loop.
  JSTART=1
  DO N=1,NMAXCH
    !FILE3 is already established for N=1, so the next loop is used
    !only for N>1.  Note that the loop increments by 2 and that we use
    !JSTART+I.  This is because previous passes have already computed
    !the values at every other position (see rearrangement of FILE3
    !below).
    IF(N>1)THEN
      DO I=1,2**N-1,2
        CALL F(AMID+HALF*CHEPTS(JSTART+I),FILE3(JSTART+I))
      ENDDO
    ENDIF

    !Calculate ANS in order N.
    ANS=ZERO
    DO I=JSTART,JSTART+2**N
      ANS=ANS+CHEWTS(I)*FILE3(I)
    ENDDO
    ANS=ANS*HALF

    !A test for convergence is appropriate when N>1.
    IF(N>1)THEN
      !Set IPR=1 in the special case where all of the integrals are 0.
      IF(ABS(ANS)<=1.0D-14)THEN
        IPR=1
        RETURN
      ENDIF

      !A conservative estimate of convergence requires comparing the
      !most recent answer in PRANS(2) with the present answers in ANS.  
      CONVERGE=.TRUE.
      ERRMAX=ZERO
      IF(ABS(ANS)>=1.0D-14)ERRMAX=ABS(ANS-PRANS(2))/ABS(ANS)
      IF(ERRMAX>ACC)CONVERGE=.FALSE.

      !If we have converged, set IPR=0 and return.
      IF(CONVERGE)THEN
        IPR=0
        RETURN
      ENDIF
    ENDIF

    !Now consider successive extrapolations to infinite order.
    IF(N>1)THEN
      EXTRAP(2)=ANS-(ANS-PRANS(2))**2/(ANS-2*PRANS(2)+PRANS(1))
    ENDIF
    IF(N>2.AND.EXTRAP(2)/=ZERO)THEN
      ERRMAX=ABS(EXTRAP(2)-EXTRAP(1))/ABS(EXTRAP(2))
      IF(ERRMAX<ACC)THEN
        IPR=0
        ANS=EXTRAP(2)
        RETURN
      ENDIF
    ENDIF

    !If we have reached the maximum order of quadrature without
    !converging in a relative sense, check for convergence to a very
    !small number.  If it is achieved, set IPR=1 and return.
    !If it is not achieved either, set IPR=2 and return.
    IF(N==NMAXCH)THEN
      ERRMAX=ABS(ANS-PRANS(2))
      IF(ERRMAX<ACC)THEN
        IPR=1
      ELSE
        IPR=2
      ENDIF
      RETURN
    ENDIF

    !Place ANS into PRANS in preparation for the next N value.
    PRANS(1)=PRANS(2)
    PRANS(2)=ANS
    EXTRAP(1)=EXTRAP(2)

    !Copy the values in FILE3 to new positions appropriate for the
    !next higher approximation.
    DO I=JSTART,JSTART+2**N
      FILE3(2*I+3-N)=FILE3(I)
    ENDDO

    !Change JSTART so that it points to the first values in CHEPTS and
    !CHEWTS that apply to the next order of approximation.
    JSTART=JSTART+2**N+1
  ENDDO
END SUBROUTINE INTEG

FUNCTION ELEM0(P,Q,R,S,T,U)
!This function evaluates eq. A.1 of L. A. Viehland, Chem. Phys. 179
!  (1994) 71.
  USE DATA
  IMPLICIT NONE

  !Declare the arguments.
  !P, Q, R, S, T and U are the indices on the left of eq. A.1.
  !ELEM0 and ELEM1 are function subprograms.
  !TERM is an intermediate term in calculating ELEM0.
  INTEGER,INTENT(IN)::P,Q,R,S,T,U
  REAL(SELECTED_REAL_KIND(18))::ELEM0,ELEM1,TERM

  !Declare the common variables.
  REAL(SELECTED_REAL_KIND(18))::ALL,BEL,BET,CL1,CL2
  COMMON/REAL2/ALL,BEL,BET,CL1,CL2

  !Start with the diagonal term.
  ELEM0=ELEM1(P,Q,R,S,T,U)

  !Add in the terms that are multiplied by the skewness.
  TERM=TWO*DBLE((U+3)*(U+2)*(U+1))*ELEM1(P,Q,R,S,T,U+3)
  IF(U>=1)TERM=TERM+THREE*DBLE((U+1)*U)*ELEM1(P,Q,R,S,T,U+1)
  IF(U>=2)TERM=TERM+THREE/TWO*DBLE(U-1)*ELEM1(P,Q,R,S,T,U-1)
  IF(U>=3)TERM=TERM+ELEM1(P,Q,R,S,T,U-3)/FOUR
  ELEM0=ELEM0+ALL*TERM/SQRT(TWO)/THREE

  !Add in the terms that are multiplied by the parallel kurtosis.
  TERM=FOUR*DBLE((U+4)*(U+3)*(U+2)*(U+1))*ELEM1(P,Q,R,S,T,U+4)
  IF(U>=1)TERM=TERM+TWO*FOUR*DBLE((U+2)*(U+1)*U)*ELEM1(P,Q,R,S,T,U+2)
  IF(U>=2)TERM=TERM+THREE*TWO*DBLE((U-1)*U)*ELEM1(P,Q,R,S,T,U)
  IF(U>=3)TERM=TERM+TWO*DBLE(U-2)*ELEM1(P,Q,R,S,T,U-2)
  IF(U>=4)TERM=TERM+ELEM1(P,Q,R,S,T,U-4)/FOUR
  ELEM0=ELEM0+(BEL-THREE)/(TWO*THREE*FOUR)*TERM

  !Add in the terms that are multiplied by the perpendicular kurtosis.
  TERM=FOUR*DBLE((S+4)*(S+3)*(S+2)*(S+1))*ELEM1(P,Q,R,S+4,T,U)
  IF(S>=1)TERM=TERM+FOUR*TWO*DBLE((S+2)*(S+1)*S)*ELEM1(P,Q,R,S+2,T,U)
  IF(S>=2)TERM=TERM+THREE*TWO*DBLE((S-1)*S)*ELEM1(P,Q,R,S,T,U)
  IF(S>=3)TERM=TERM+TWO*DBLE(S-2)*ELEM1(P,Q,R,S-2,T,U)
  IF(S>=4)TERM=TERM+ELEM1(P,Q,R,S-4,T,U)/FOUR
  TERM=TERM+FOUR*DBLE((T+4)*(T+3)*(T+2)*(T+1))*ELEM1(P,Q,R,S,T+4,U)
  IF(T>=1)TERM=TERM+FOUR*TWO*DBLE((T+2)*(T+1)*T)*ELEM1(P,Q,R,S,T+2,U)
  IF(T>=2)TERM=TERM+THREE*TWO*DBLE((T-1)*T)*ELEM1(P,Q,R,S,T,U)
  IF(T>=3)TERM=TERM+TWO*DBLE(T-2)*ELEM1(P,Q,R,S,T-2,U)
  IF(T>=4)TERM=TERM+ELEM1(P,Q,R,S,T-4,U)/FOUR
  ELEM0=ELEM0+(BET-THREE)/(TWO*THREE*FOUR)*TERM

  !Add in the terms that are multiplied by the first correlation.
  TERM=DBLE((S+2)*(S+1)*(U+1))*ELEM1(P,Q,R,S+2,T,U+1)+DBLE((T+2)*(T+1)*&
    &(U+1))*ELEM1(P,Q,R,S,T+2,U+1)
  IF(S>=1)TERM=TERM+DBLE(S*(U+1))*ELEM1(P,Q,R,S,T,U+1)
  IF(T>=1)TERM=TERM+DBLE(T*(U+1))*ELEM1(P,Q,R,S,T,U+1)
  IF(S>=2)TERM=TERM+DBLE(U+1)*ELEM1(P,Q,R,S-2,T,U+1)/FOUR
  IF(T>=2)TERM=TERM+DBLE(U+1)*ELEM1(P,Q,R,S,T-2,U+1)/FOUR
  IF(U>=1)THEN
    TERM=TERM+DBLE((S+2)*(S+1))*ELEM1(P,Q,R,S+2,T,U-1)/TWO+&
      &DBLE((T+2)*(T+1))*ELEM1(P,Q,R,S,T+2,U-1)/TWO
    IF(S>=1)TERM=TERM+DBLE(S)*ELEM1(P,Q,R,S,T,U-1)/TWO
    IF(T>=1)TERM=TERM+DBLE(T)*ELEM1(P,Q,R,S,T,U-1)/TWO
    IF(S>=2)TERM=TERM+ELEM1(P,Q,R,S-2,T,U-1)/(TWO*FOUR)
    IF(T>=2)TERM=TERM+ELEM1(P,Q,R,S,T-2,U-1)/(TWO*FOUR)
  END IF
  ELEM0=ELEM0+CL1*SQRT(TWO)*TERM

  !Add in the terms that are multiplied by the second correlation.
  TERM=DBLE((S+2)*(S+1)*(U+2)*(U+1))*ELEM1(P,Q,R,S+2,T,U+2)+&
    &DBLE((S+T)*(U+2)*(U+1))*ELEM1(P,Q,R,S,T,U+2)+&
    &DBLE(U*(S+2)*(S+1))*ELEM1(P,Q,R,S+2,T,U)+DBLE((S+T)*U)*&
    &ELEM1(P,Q,R,S,T,U)+DBLE((T+2)*(T+1)*(U+2)*(U+1))*&
    &ELEM1(P,Q,R,S,T+2,U+2)+DBLE(U*(T+2)*(T+1))*ELEM1(P,Q,R,S,T+2,U)
  IF(S>=2)TERM=TERM+DBLE((U+2)*(U+1))*ELEM1(P,Q,R,S-2,T,U+2)/FOUR+&
    &DBLE(U)*ELEM1(P,Q,R,S-2,T,U)/&
    &FOUR
  IF(U>=2)TERM=TERM+DBLE((S+2)*(S+1))*ELEM1(P,Q,R,S+2,T,U-2)/FOUR+&
    &DBLE(S+T)*ELEM1(P,Q,R,S,T,U-2)/FOUR+DBLE((T+2)*(T+1))*&
    &ELEM1(P,Q,R,S,T+2,U-2)/FOUR
  IF(S>=2.AND.U>=2)TERM=TERM+ELEM1(P,Q,R,S-2,T,U-2)/FOUR/FOUR
  IF(T>=2)TERM=TERM+DBLE((U+2)*(U+1))*ELEM1(P,Q,R,S,T-2,U+2)/FOUR+&
    &DBLE(U)*ELEM1(P,Q,R,S,T-2,U)/FOUR
  IF(T>=2.AND.U>=2)TERM=TERM+ELEM1(P,Q,R,S,T-2,U-2)/FOUR/FOUR
  ELEM0=ELEM0+(CL2-ONE)*TERM
END FUNCTION ELEM0

FUNCTION ELEM1(P,Q,R,S,T,U)
!Rather than use eq. A.4 of L. A. Viehland, Chem. Phys. 179 (1994) 71,
!  this function evaluates the 3T matrix elements by using a pair of
!  equations equivalent to eq. A.31 of S. L. Lin et al., Chem. Phys. 37
!  (1979) 411.  This function gives the matrix elements in terms of
!  reducible collision integrals.
  USE DATA
  IMPLICIT NONE

  !Declare the arguments.
  !P, Q, R, S, T and U are the indices on the left of eq. A.4.
  INTEGER,INTENT(IN)::P,Q,R,S,T,U
  REAL(SELECTED_REAL_KIND(18))::ELEM1,ELEM2

  !Declare the common variables.
  REAL(SELECTED_REAL_KIND(18))::AM0,AM1,TEMP,RAM,TPARA,TPERP,VDBAR,&
    TT1,TT2,TT3,TT10,TT11
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(0:11,0:11,0:13,0:23,0:23,&
    0:23)::SAVE1,SAVE2
  COMMON/REAL1/AM0,AM1,TEMP,RAM,TPARA,TPERP,VDBAR
  COMMON/REAL3/TT1,TT2,TT3,TT10,TT11
  COMMON/REAL4/SAVE1,SAVE2

  !Declare the local variables.
  !P1, P2, P3, Q1, Q2, Q3, R1, R2, R3, S1, S2, T1, T2, U1 and U2
  !  are the 15 summation indices.
  !EXPON is an exponent
  !TERM is an intermediate term in the calculations of ELEM1
  INTEGER::P1,P2,P3,Q1,Q2,Q3,R1,R2,R3,S1,S2,T1,T2,U1,U2,EXPON
  INTEGER::P1MAX,R1MAX,S1MAX,T1MAX,U1MAX
  REAL(SELECTED_REAL_KIND(18))::TERM

  !Check whether this has already been calculated.
  ELEM1=SAVE1(P,Q,R,S,T,U)
  IF(ELEM1/=ZERO)RETURN

  !Start the 15 sums.
  P1MAX=P/2
  DO P1=0,P1MAX
  DO P2=0,P-2*P1
  DO P3=0,P2
  DO Q1=0,Q/2
  DO Q2=0,Q-2*Q1
  DO Q3=0,Q2
  R1MAX=R/2
  DO R1=0,R1MAX
  DO R2=0,R-2*R1
  DO R3=0,R2
    IF(P3+Q3+R3==0)CYCLE
  S1MAX=S/2
  DO S1=0,S1MAX
  DO S2=0,S-2*S1
    IF(MOD(P+P2+S+S2,2)/=0)CYCLE
  T1MAX=T/2
  DO T1=0,T1MAx
  DO T2=0,T-2*T1
    IF(MOD(Q+Q2+T+T2,2)/=0)CYCLE
  U1MAX=U/2
  DO U1=0,U1MAX
  DO U2=0,U-2*U1
    IF(MOD(R+R2+U+U2,2)/=0)CYCLE

    !Start the next term by getting the value of the appropriate
    !  reducible collision integral.
    TERM=ELEM2(P3,Q3,R3,P2-P3+S2,Q2-Q3+T2,R2-R3+U2)

    !Include the factors that depend on the effective temperatures.
    EXPON=P3+Q3+R3
    IF(EXPON/=0)TERM=TERM*TT1**EXPON
    EXPON=S2+T2-P2-Q2
    IF(EXPON/=0)TERM=TERM*SQRT(TT10)**EXPON
    EXPON=U2-R2
    IF(EXPON/=0)TERM=TERM*SQRT(TT11)**EXPON
    EXPON=P2+Q2-P3-Q3
    IF(EXPON/=0)TERM=TERM*TT2**EXPON
    EXPON=R2-R3
    IF(EXPON/=0)TERM=TERM*TT3**EXPON
    EXPON=P+Q+S+T-2*(P1+Q1+S1+T1)-P2-Q2-S2-T2
    IF(EXPON/=0)TERM=TERM*SQRT(ONE-TT10)**EXPON
    EXPON=R+U-2*(R1+U1)-R2-U2
    IF(EXPON/=0)TERM=TERM*SQRT(ONE-TT11)**EXPON

    !Now include the factors that are constants.  Note, however, that
    !  the three factors of G(...) incorporate only PI**3/2 in the
    !  denominator.  The other factor of PI**3/2 is in ELEM2.
    IF(MOD(P1+Q1+R1+S1+T1+U1,2)/=0)TERM=-TERM
    EXPON=P+Q+R-2*(P1+Q1+R1+S1+T1+U1)
    IF(EXPON/=0)TERM=TERM*TWO**EXPON
    TERM=TERM*H(P+1)*H(Q+1)*H(R+1)/H(P1+1)/H(Q1+1)/H(R1+1)/H(S1+1)/&
      &H(T1+1)/H(U1+1)
    TERM=TERM/H(P-2*P1-P2+1)/H(Q-2*Q1-Q2+1)/H(R-2*R1-R2+1)
    TERM=TERM/H(S-2*S1-S2+1)/H(T-2*T1-T2+1)/H(U-2*U1-U2+1)
    TERM=TERM/H(S2+1)/H(T2+1)/H(U2+1)
    TERM=TERM/H(P3+1)/H(P2-P3+1)/H(Q3+1)/H(Q2-Q3+1)/H(R3+1)/H(R2-R3+1)
    TERM=TERM*G((P+S-2*P1-2*S1-P2-S2)/2+1)
    TERM=TERM*G((Q+T-2*Q1-2*T1-Q2-T2)/2+1)
    TERM=TERM*G((R+U-2*R1-2*U1-R2-U2)/2+1)
    ELEM1=ELEM1+TERM
  END DO !U2
  END DO !U1
  END DO !T2
  END DO !T1
  END DO !S2
  END DO !S1
  END DO !R3
  END DO !R2
  END DO !R1
  END DO !Q3
  END DO !Q2
  END DO !Q1
  END DO !P3
  END DO !P2
  END DO !P1
  SAVE1(P,Q,R,S,T,U)=ELEM1

  !Make use of symmetry.
  SAVE1(Q,P,R,T,S,U)=ELEM1
END FUNCTION ELEM1

FUNCTION ELEM2(P3,Q3,R3,S3,T3,U3)
!This function gives the reducible collision integrals in terms of the
!  irreducible collision integrals of function ELEM3.
  USE DATA
  IMPLICIT NONE

  !Declare the arguments.
  !P3, Q3, R3, S3, T3 and U3 are the indices of ELEM2
  INTEGER,INTENT(IN)::P3,Q3,R3,S3,T3,U3
  REAL(SELECTED_REAL_KIND(18))::ELEM2,ELEM3

  !Declare the common variables.
  REAL(SELECTED_REAL_KIND(18))::AM0,AM1,TEMP,RAM,TPARA,TPERP,VDBAR,&
    TT1,TT2,TT3,TT10,TT11
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(0:11,0:11,0:13,0:23,0:23,&
    0:23)::SAVE1,SAVE2
  COMMON/REAL1/AM0,AM1,TEMP,RAM,TPARA,TPERP,VDBAR
  COMMON/REAL3/TT1,TT2,TT3,TT10,TT11
  COMMON/REAL4/SAVE1,SAVE2

  !Declare the local variables.
  !P4, P5, Q4, Q5, R4, R5, L, L1, K and K1 are the 10 summation indices.
  !EXPON is an exponent
  !TERM is an intermediate term in the calculations of ELEM2
  INTEGER::P4,P5,Q4,Q5,R4,R5,L,L1,K,K1,EXPON
  INTEGER::LMAX,KMAX
  REAL(SELECTED_REAL_KIND(18))::TERM

  !Check whether this has already been calculated.
  ELEM2=SAVE2(P3,Q3,R3,S3,T3,U3)
  IF(ELEM2/=ZERO)RETURN

  !Start the 10 sums.
  DO P4=0,P3
  DO P5=0,P4
  DO Q4=0,Q3
  DO Q5=0,Q4
    IF(MOD(P5+Q5,2)/=0)CYCLE
  DO R4=0,R3
  DO R5=0,R4
    IF(MOD(P4+Q4+R5,2)/=0)CYCLE
  LMAX=(P5+Q5)/2
  DO L=0,LMAX
  DO L1=0,P4-P5+Q4-Q5+R4-R5+2*L
  KMAX=(P4+Q4+R5)/2
  DO K=0,KMAX
  DO K1=0,R4-R5
    IF(P3-P4+Q3-Q4+2*K+K1<1)CYCLE

    !Start the next term by getting the value of the appropriate
    !  irreducible collision integral.
    TERM=ELEM3(P3+Q3+S3+T3-P4+P5-Q4+Q5+R5-2*L,P4-P5+Q4-Q5+R3-R5+U3+&
      &2*L-L1,P3-P4+Q3-Q4+2*K+K1)

    !Include the factors that depend on the effective temperatures.
    EXPON=(P4-P5+Q4-Q5-R5)/2+L
    IF(EXPON/=0)TERM=TERM*((ONE-TT10)/(ONE-TT11))**EXPON
    IF(L1/=0)TERM=TERM*VDBAR**L1

    !Now include the factors that are constants.  Note that the three
    !  factors of G(...) incorporate the factor of PI**3/2 that was
    !  not put into ELEM1.
    IF(MOD(Q5+R4+K+K1,2)/=0)TERM=-TERM
    TERM=TERM*H(P3+1)*H(Q3+1)*H(R3+1)/H(P5+1)/H(Q5+1)/H(R5+1)
    TERM=TERM/H(P3-P4+1)/H(P4-P5+1)/H(Q3-Q4+1)/H(Q4-Q5+1)/H(R3-R4+1)
    TERM=TERM*H((P5+Q5)/2+1)/H(L+1)/H((P5+Q5)/2-L+1)
    TERM=TERM*H(P4-P5+Q4-Q5+R4-R5+2*L+1)/H(L1+1)/&
      &H(P4-P5+Q4-Q5+R4-R5+2*L-L1+1)
    TERM=TERM/H(K+1)/H((P4+Q4+R5)/2-K+1)/H(K1+1)/H(R4-R5-K1+1)
    TERM=TERM*G((P4-P5+Q4-Q5+R5)/2+1)*G((P5+Q5)/2+1)
    TERM=TERM*G((P3-P5+Q5+S3)/2+1)*G((P5+Q3-Q5+T3)/2+1)
    TERM=TERM/H((P3+Q3+S3+T3)/2+1)
    ELEM2=ELEM2+TERM
  END DO !K1
  END DO !K
  END DO !L1
  END DO !L
  END DO !R5
  END DO !R4
  END DO !Q5
  END DO !Q4
  END DO !P5
  END DO !P4
  SAVE2(P3,Q3,R3,S3,T3,U3)=ELEM2

  !Make use of symmetry.
  SAVE2(Q3,P3,R3,T3,S3,U3)=ELEM2
END FUNCTION ELEM2

FUNCTION ELEM3(P,Q,L)
!This function evaluates eq. A.11 of L. A. Viehland, Chem. Phys. 179
!  (1994) 71, which is the same as eq. A.34 of Lin et al., Chem. Phys.
!  37 (1979) 411.  Here we do not use the change of variables that leads
!  to eq. A.18 of Viehland.  Note from eq. A.35 of Lin et al., however,
!  that the left-hand side of eq. A.16 of Viehland should be g^2, not g.
  USE DATA
  IMPLICIT NONE

  !Declare the arguments.
  !P, Q and L are the indices on the left of eq. A.11.
  REAL(SELECTED_REAL_KIND(18))::ELEM3
  INTEGER,INTENT(IN)::P,Q,L

  !Declare the common variables.
  REAL(SELECTED_REAL_KIND(18))::ACCT,ACCT1,ACCT2,ACCT3,TT1,TT2,TT3,&
    TT10,TT11
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(0:MAXL,0:MAXL,1:MAXL)::&
    SAVEANS
  INTEGER::IPTRAN,IQTRAN,ILTRAN
  COMMON/REAL0/ACCT,ACCT1,ACCT2,ACCT3
  COMMON/REAL3/TT1,TT2,TT3,TT10,TT11
  COMMON/REAL6/SAVEANS
  COMMON/INT1/ILTRAN,IPTRAN,IQTRAN

  !Declare the local variables.
  !IPR is a flag for lack of convergence.
  !A and B are the endpoints for the integration of F1.
  !EA and EB are the values of the integrands at A and B.
  !ACC is the accuracy desired to computed ANS, the integral.
  INTEGER::IPR
  REAL(SELECTED_REAL_KIND(18))::A,B,EA,EB,ACC,ANS
  EXTERNAL F1

  IPTRAN=P
  IQTRAN=Q
  ILTRAN=L

  !Check to see if the irreducible collision integral is available.
  !  If not, determine it.
  ELEM3=SAVEANS(IPTRAN,IQTRAN,ILTRAN)
  IF(ELEM3/=ZERO)RETURN
  A=ZERO
  B=ONE
  EA=ZERO
  CALL F1(B,EB)
  ACC=ACCT
  IPR=1
  CALL INTEG(F1,A,B,EA,EB,ANS,ACC,IPR)
  SAVEANS(IPTRAN,IQTRAN,ILTRAN)=ANS
  ELEM3=ANS
END FUNCTION ELEM3

FUNCTION XN(IP,IQ,IR,IS,IT,IU) 
  !This function evaluates eq. (47) of L. A. Viehland, Chem. Phys. 179
  !  (1994) 71.
  USE DATA
  IMPLICIT NONE

  !Declare the arguments.
  INTEGER,INTENT(IN)::IP,IQ,IR,IS,IT,IU

  !Declare the common variables.
  REAL(SELECTED_REAL_KIND(18))::ALL,BEL,BET,CL1,CL2
  COMMON/REAL2/ALL,BEL,BET,CL1,CL2

  !Declare the variables.
  REAL(SELECTED_REAL_KIND(18))::XN,DELTA

  XN=ZERO
  IF(IP<0.OR.IQ<0.OR.IR<0.OR.IS<0.OR.IT<0.OR.IU<0)RETURN
  IF(IP==IS.AND.IQ==IT)THEN
    IF(IR==IU)XN=ONE
    XN=XN+ALL*SQRT(TWO)/(TWO*THREE*FOUR)*DELTA(3,IR,IU)
    XN=XN+(BEL-THREE)/(TWO*THREE*FOUR*FOUR)*DELTA(4,IR,IU)
  END IF
  IF(IR==IU)THEN
    IF(IQ==IT)XN=XN+(BET-THREE)/(TWO*THREE*FOUR*FOUR)*DELTA(4,IP,IS)
    IF(IP==IS)XN=XN+(BET-THREE)/(TWO*THREE*FOUR*FOUR)*DELTA(4,IQ,IT)
  END IF
  IF(IQ==IT)XN=XN+DELTA(2,IP,IS)*(CL1*SQRT(TWO)/(TWO*FOUR)*&
    &DELTA(1,IR,IU)+(CL2-ONE)/(FOUR*FOUR)*DELTA(2,IR,IU))
  IF(IP==IS)XN=XN+DELTA(2,IQ,IT)*(CL1*SQRT(TWO)/(TWO*FOUR)*&
    &DELTA(1,IR,IU)+(CL2-ONE)/(FOUR*FOUR)*DELTA(2,IR,IU))
END FUNCTION XN

FUNCTION DELTA(I,J,K) 
  !This function evaluates eq. (49) of L. A. Viehland, Chem. Phys. 179
  !  (1994) 71.
  USE DATA
  IMPLICIT NONE

  !Declare the arguments.
  INTEGER,INTENT(IN)::I,J,K

  !Declare the local variables.
  REAL(SELECTED_REAL_KIND(18))::DELTA
  INTEGER::IS

  DELTA=ZERO
  IF(MOD(I+J+K,2)/=0)RETURN
  IF(ABS(I-J)>K.OR.K>I+J)RETURN 
  IS=(I+J+K)/2 
  IF(IS==J)THEN
    DELTA=ONE
  ELSE
    DELTA=TWO**(IS-J) 
  END IF
  DELTA=DELTA/H(IS-I+1)*H(I+1)/H(IS-J+1)*H(K+1)/H(IS-K+1) 
END FUNCTION DELTA

FUNCTION XM(IP,IQ,IR,IS,IT,IU) 
  !This function evaluates eq. (58) of L. A. Viehland, Chem. Phys. 179
  !  (1994) 71.
  USE DATA
  IMPLICIT NONE

  !Declare the arguments.
  INTEGER,INTENT(IN)::IP,IQ,IR,IS,IT,IU

  !Declare the local variables.
  REAL(SELECTED_REAL_KIND(18))::XM,XN

  XM=ZERO
  IF(IP<0.OR.IQ<0.OR.IR<0.OR.IS<0.OR.IT<0.OR.IU<0)RETURN
  XM=XN(IP,IQ,IR,IS,IT,IU)-XN(IP,IQ,IR,0,0,0)*XN(0,0,0,IS,IT,IU)
END FUNCTION XM

SUBROUTINE SIMQE(D,NN,B,N,KS)
  !This subroutine solves linear, simultaneous algebraic equations.
  USE DATA
  IMPLICIT NONE

  !Declare the arguments.
  INTEGER,INTENT(IN)::NN,N
  INTEGER,INTENT(OUT)::KS
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(NN*NN),INTENT(IN)::D
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(NN),INTENT(INOUT)::B

  !Declare the local variables.
  REAL(SELECTED_REAL_KIND(18))::TOL,BIGA,SAVE
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(NN*NN)::A
  INTEGER::NI,IJ,NM,K,L,JJ,J,JY,IT,I,UMAX,I1,I2,IQS,IX,IXJ,JX,IXJX,JJX,&
    &NY,IA,IB,IC
  NI=NN-N
  IJ=0
  NM=0
  DO K=1,N
    DO L=1,N
      IJ=IJ+1
      NM=NM+1
      A(IJ)=D(NM)
    END DO
    NM=NM+NI
  END DO

  !Forward solution.
  TOL=ZERO
  KS=0
  JJ=-N
  Jloop: DO J=1,N
    JY=J+1
    JJ=JJ+N+1
    BIGA=ZERO
    IT=JJ-J
    
    !Search for the maximum coefficient in the column.
    DO I=J,N
      IJ=IT+I
      IF(ABS(BIGA)<ABS(A(IJ)))THEN
        BIGA=A(IJ)
        UMAX=I
      END IF
    END DO

    !Test for the pivot less than the tolerance, which indicates that
    !  the matrix is singular.
    IF(ABS(BIGA)<=TOL)THEN
      KS=1
      RETURN
    END IF

    !Interchange the rows if necessary.  Then divide the rows by the
    !  leading coefficient.
    I1=J+N*(J-2)
    IT=UMAX-J
    DO K=J,N
      I1=I1+N
      I2=I1+IT
      SAVE=A(I1)
      A(I1)=A(I2)
      A(I2)=SAVE
      A(I1)=A(I1)/BIGA
    END DO
    SAVE=B(UMAX)
    B(UMAX)=B(J)
    B(J)=SAVE/BIGA

    !Eliminate the next variable.
    IF(J/=N)THEN
      IQS=N*(J-1)
      DO IX=JY,N
        IXJ=IQS+IX
        IT=J-IX
        DO JX=JY,N
          IXJX=N*(JX-1)+IX
          JJX=IXJX+IT
          A(IXJX)=A(IXJX)-(A(IXJ)*A(JJX))
        END DO
        B(IX)=B(IX)-(B(J)*A(IXJ))
      END DO
    END IF
  END DO Jloop

  !Back solution.
  NY=N-1
  IT=N*N
  DO J=1,NY
    IA=IT-J
    IB=N-J
    IC=N
    DO K=1,J
      B(IB)=B(IB)-A(IA)*B(IC)
      IA=IA-N
      IC=IC-1
    END DO
  END DO
END SUBROUTINE SIMQE
