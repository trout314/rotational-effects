MODULE GLOBAL
  !This is global module for program COMBINE.
  IMPLICIT NONE

  !Declare the global named constants.
  !  MAXE is the maximum number of energies in each of the 3 regions.
  !    It should be one larger than 2 to some power.
  !  MAXL is the number of transport cross sections to be calculated.
  !  MAXORBITS is the maximum number of sets of orbiting parameters.
  !  ZERO, ONE, TWO, THREE, FOUR and TEN are high-precision numbers.
  !  INOOMAX is the maximum number of regions of multiple orbiting.
  INTEGER,PARAMETER::MAXE=65,MAXL=30,MAXORBITS=600,INOOMAX=40
  REAL(SELECTED_REAL_KIND(18)),PARAMETER::ZERO=0,ONE=1,TWO=2,THREE=3,&
    FOUR=4,TEN=10

  !Declare the global integer variables.
  !  MAXANG is 90 for homonuclear diatomics and 180 for heteronuclears.
  !  NANGLES is the number of angles to be studied.
  !  ERR allows for recovery from input/output errors.
  !  NMAXCH is the maximum order of approximation in using the Curtiss-
  !    Clenshaw method of integration.  It is important that the value
  !    of NMAXCH here be the same as in CCgen.f90.
  !  NPTSCH is the size of the arrays CHEPTS and CHEWTS.
  !  NV is the number of tabulated potential pairs (R,V).
  !  NLONG is the exponent of the long-range potential.
  !  NO is the number of valid sets of orbiting parameters.
  !  INOO is the number of orbiting intervals at any energy.
  !  NAITK is the number of orbiting values at a particular E, but it is
  !    reset to indicate the number of valid entries in BO and RO.
  INTEGER::MAXANG,NANGLES,ERR,NMAXCH,NPTSCH,NV,NLONG,NO,INOO,NAITK

  !Declare the global integer arrays.
  !  N1, N2 and N3 are the number of energies in each region.
  !  NOO stores the starting and ending indices for orbiting parameters.
  INTEGER,DIMENSION(:),ALLOCATABLE::N1,N2,N3
  INTEGER,DIMENSION(2,INOOMAX)::NOO=0

  !Declare the global real variables.
  !  ORBIT stores the sets of orbiting parameters.
  !  ECC stores the energy ranges for orbiting.
  !  B0 and RO contain the orbiting impact parameters and orbiting
  !    separations corresponding to orbiting energy EE.
  !  ENERGY gives the individual energies.
  !  POTEN contains the separations and the potential values.
  !  FILES3 gives the individual cross sections.
  !  FILES4 gives the Chebychev coefficients.
  !  The quadrature points and weights are placed in CHEPTS and CHEWTS.
  !  CPOTEN contains the parameters for a spline fit of the potential.
  !  EMIN1-EMIN4 are the energies delimiting each region.
  !  EPS gives the well depths at the various angles.
  !  HPI is PI/2, where PI is the usual parameter for circles.
  !  EMIN and EMAX are the minimum and maximum energies (in hartree).
  !  ACCURACY is the desired accuracy of the transport cross sections.
  !  EXPON is the inverse power of the R dependence at short range.
  !  CSHORT is the coefficient of the short-range inverse power.
  !  CLONG is the coefficient of the long-range inverse power.
  !  ACC1 is used to compute the integrands with greater accuracy
  !     than is required (via ACCURACY) for the integrals.
  !  ED and EC are the minimum and maximum orbiting energies.
  !  EC2 is always equal to ten times global variable EC.
  !  EE is an energy of interest.
  !  ROMAX and BOMAX are the maximum values of RO and BO.
  !  EING and BING are used to transfer E and B to functions GA and GB.
  !  R1ING and R2ING are used to similarly transfer separations.
  !  EAING and EBING are used to similarly transfer values at endpoints.
  !  T1 and T2 are temporary values used to determine EPS.
  !  WEIGHT is the sum of all the separate weights.
  !  TERM is a part of WEIGHT.
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(3,MAXORBITS)::ORBIT
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(2,INOOMAX)::ECC
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(INOOMAX)::BO,RO
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(:,:),ALLOCATABLE::ENERGY,POTEN
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(:,:,:),ALLOCATABLE::FILES3,&
    FILES4
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(:),ALLOCATABLE::CHEPTS,CHEWTS,&
    CPOTEN,EMIN1,EMIN2,EMIN3,EMIN4,EPS
  REAL(SELECTED_REAL_KIND(18))::HPI,PI,EMIN,EMAX,ACCURACY,EXPON,CSHORT,&
    CLONG,ACC1,ED,EC,EC2,EE,ROMAX,BOMAX,EING,BING,R1ING,R2ING,EAING,&
    EBING

  !Declare the global character variables.
  !  CCFILE is the fully qualified name of the quadrature file .
  !  FILENAME is FILENAME1//FILENAME2//FILENAME3//FILENAME4
  !  COMMENT is a line of commentary describing the system.
  CHARACTER(LEN=78)::CCFILE
  CHARACTER(LEN=80)::FILENAME,COMMENT
  CHARACTER(LEN=7)::FILENAME1
  CHARACTER(LEN=1)::FILENAME2,FILENAME3,FILENAME4
END MODULE GLOBAL

PROGRAM PCandBMM
  !This program was written in July, 2014 by Larry A. Viehland.  It
  !computes transport cross sections at different angles and then
  !combines them as prescribed by the BMM approximation.  The cross
  !sections are computed using the techniques in L. A. Viehland and
  !Y. Chang, Comp. Phys. Comm. 181 (2010) 1687-1696.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the local integer variables.
  !  I is a loop index for the different energies.
  !  NUMBER indicates which data set (angle) is being considered.
  !  NUMBER2 is an angle equivalent to NUMBER.
  INTEGER::I,NUMBER,NUMBER2

  !Ask the user for information
  DO
    PRINT *,'Enter 90 if the diatom is homonuclear, 180 if not.'
    READ *,MAXANG
    IF(MAXANG==90)THEN
      PRINT *,'What is N if 90/N is the angular spacing?'
      PRINT *,'Warning: 90 must be evenly divisible by N.'
      READ *,NANGLES
      EXIT
    ELSEIF(MAXANG==180)THEN
      PRINT *,'What is N if 180/N is the angular spacing?'
      PRINT *,'Warning: 180 must be evenly divisible by N.'
      READ *,NANGLES
      EXIT
    ELSE
      PRINT *,'You made an incorrect response.'
    ENDIF
  ENDDO

  !Allocate space for the arrays.
  NANGLES=NANGLES+1
  ALLOCATE(N1(NANGLES),N2(NANGLES),N3(NANGLES),ENERGY(NANGLES,3*MAXE),&
    &EMIN1(NANGLES),EMIN2(NANGLES),EMIN3(NANGLES),EMIN4(NANGLES),&
    &EPS(NANGLES),FILES3(0:NANGLES,3*MAXE,MAXL),&
    &FILES4(0:NANGLES,3*MAXE,MAXL))

  !Open the input file for the Curtiss-Clenshaw quadrature points and
  !weights and read NMAXCH and CCFILE.  Finally, close the file.
  OPEN(UNIT=10,FILE='CCgen.in',STATUS='OLD',IOSTAT=ERR)
  IF(ERR/=0)THEN
    PRINT *,'Could not open the old file named CCgen.in.'
    STOP
  ENDIF
  READ(10,*)NMAXCH
  READ(10,'(A)')CCFILE
  CLOSE(10)

  !Assign values to HPI and PI.
  HPI=TWO*ATAN(ONE)
  PI=TWO*HPI

  !Allocate space for the quadrature points and weights.
  NPTSCH=0
  DO I=1,NMAXCH
    NPTSCH=NPTSCH+1+2**I
  ENDDO
  ALLOCATE(CHEPTS(NPTSCH),CHEWTS(NPTSCH))

  !Read the quadrature points and weights.
  OPEN(10,FILE=CCFILE,IOSTAT=ERR,STATUS='OLD')
  IF(ERR/=0)THEN
    PRINT *,'Could not open the old file named ',TRIM(CCFILE)
    PRINT *,'Be sure that you have permission to read this file.'
    STOP
  ENDIF
  DO I=1,NPTSCH
    READ(10,'(18X,2F23.20)',IOSTAT=ERR)CHEPTS(I),CHEWTS(I)
    IF(ERR/=0)THEN
      PRINT *,'File ',TRIM(CCFILE),' does not work with this program.'
      PRINT *,'Be sure that it was created with program CCgen.f95.'
      STOP
    ENDIF
  ENDDO

  !If everything is consistent, the next READ statement should fail.
  !Note that EMIN and EMAX are dummy variables at this point.
  READ(10,'(18X,2F23.20)',IOSTAT=ERR)EMIN,EMAX
  IF(ERR==0)THEN
    PRINT *,'File ',TRIM(CCFILE),' does not work with this program.'
    PRINT *,'Be sure that it was created with program CCgen.f95.'
    STOP
  ENDIF
  CLOSE(10)

  !Determine COMMENT, ACCURACY, EMIN and EMAX.
  PRINt *,'Enter the line of comment that will be added to each &
    &output file.'
  READ '(A)',COMMENT
  PRINT *,'It is assumed that the input files are named PC.in.xxx'
  PRINT *,'where xxx is 000 for 0 degrees, 005 for 5 degrees, etc.'
  PRINT *,'This program must be revised if you want fractional degrees.'
  PRINT *
  PRINT *,'Enter the fractional accuracy or 0 to abort.'
  READ *,ACCURACY
  IF(ACCURACY<=0)STOP
  PRINT *,'Enter the minimum and maximum energies in hartree.'
  PRINT *,'Typically, EMIN=1E-9.'
  PRINT *,'You will need to look at each of the input files to find'
  PRINT *,'the smallest value from the various first lines.'
  READ *,EMIN,EMAX

  !Initialize global variable ACC1, thus requiring the integrals to be
  !calculated a bit more accurately than the overall results.
  ACC1=0.8D0*ACCURACY

  !Determine the long-range coefficient.
  PRINT *,'The long-range potential should vary as 1/r^3 (NLONG=3)'
  PRINT *,'if the system is an atomic ion with a diatomic neutral.'
  PRINT *,'For a diatomic ion with an atomic neutral, it should vary'
  PRINT *,'as 1/r^4 (NLONG=4).  Enter NLONG now.'
  READ *,NLONG

  !Determine the angle at which to start.
  PRINT *,'Enter the angle at which you want the calculations to start.'
  PRINT *,'For homonucelar diatoms, enter 1 or 2 digits from 0-90.'
  PRINT *,'For heteronucelar diatoms, enter 1-3 digits from 0-180.'
  PRINT *,'This allows you to skip angles for which may already have'
  PRINT *,'been calculated.'
  READ *,I
  IF(MAXANG==90)THEN
    I=I/(90/(NANGLES-1))+1
  ELSE
    I=I/(180/(NANGLES-1))+1
  ENDIF

  !Cycle through the input files.  For each angle, call subroutine
  !INPUT to read the potential information, subroutine ORBITS to
  !determine the orbiting parameters, subroutine REGIONS to determine
  !the (possibly multiple) orbiting regions, subroutine EFIT to
  !calculate and write the cross sections, and subroutine GEN to
  !calculate and write the momentum-transfer collision frequency.
  DO NUMBER=I,NANGLES
    !Determine FILENAME for the present angle and open the input file.
    NUMBER2=(MAXANG/(NANGLES-1))*(NUMBER-1)
    FILENAME1='PC.in.'
    FILENAME2=CHAR((NUMBER2)/100+48)
    FILENAME3=CHAR((NUMBER2)/10-(NUMBER2)/100*10+48)
    FILENAME4=CHAR((NUMBER2)-(NUMBER2)/10*10+48)
    FILENAME=TRIM(FILENAME1)//FILENAME2//FILENAME3//FILENAME4
    OPEN(UNIT=11,FILE=FILENAME,IOSTAT=ERR)
    PRINT *,'OPENED FILE ',FILENAME
    IF(ERR/=0)THEN
      PRINT *,'Could not open a file named ',FILENAME
      STOP
    ENDIF

    !Open the output files.
    FILENAME1='Q1.'
    FILENAME=TRIM(FILENAME1)//FILENAME2//FILENAME3//FILENAME4
    OPEN(UNIT=12,FILE=FILENAME,IOSTAT=ERR)
    IF(ERR/=0)THEN
      PRINT *,'Could not open a file named ',FILENAME
      STOP
    ENDIF
    FILENAME1='PC.out.'
    FILENAME=FILENAME1//FILENAME2//FILENAME3//FILENAME4
    OPEN(UNIT=13,FILE=FILENAME,IOSTAT=ERR)
    WRITE(13,'(A)')COMMENT
    WRITE(13,*)ACCURACY,EMIN,EMAX
    FILENAME1='crsect.'
    FILENAME=FILENAME1//FILENAME2//FILENAME3//FILENAME4
    OPEN(15,FILE=FILENAME,IOSTAT=ERR)
    IF(ERR/=0)THEN
      PRINT *,'Could not open a file named ',FILENAME
      STOP
    ENDIF
    FILENAME1='crosec.'
    FILENAME=FILENAME1//FILENAME2//FILENAME3//FILENAME4
    OPEN(16,FILE=FILENAME,ACCESS='DIRECT',IOSTAT=ERR,RECL=496)
    IF(ERR/=0)THEN
      PRINT *,'Could not open a file named ',FILENAME
      STOP
    ENDIF
    FILENAME1='coeffs.'
    FILENAME=FILENAME1//FILENAME2//FILENAME3//FILENAME4
    OPEN(17,FILE=FILENAME,ACCESS='DIRECT',IOSTAT=ERR,RECL=480)
    IF(ERR/=0)THEN
      PRINT *,'Could not open a file named ',FILENAME
      STOP
    ENDIF
    FILENAME1='ninval.'
    FILENAME=FILENAME1//FILENAME2//FILENAME3//FILENAME4
    OPEN(18,FILE=FILENAME,IOSTAT=ERR)
    IF(ERR/=0)THEN
      PRINT *,'Could not open a file named ',FILENAME
      STOP
    ENDIF

    CALL INPUT
    CALL ORBITS
    CALL REGIONS
    CALL EFIT
    CALL GEN
    CLOSE(11)
    CLOSE(12)
    CLOSE(13)
    CLOSE(15)
    CLOSE(16)
    CLOSE(17)
    CLOSE(18)

    !Deallocate the arrays that are allocated in subroutine INPUT.
    DEALLOCATE(POTEN,CPOTEN)
  ENDDO

END PROGRAM PCandBMM

SUBROUTINE INPUT
  !This subroutine reads information about the interaction potential.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the local integer variables.
  !  I is a loop index.
  INTEGER::I

  !Read the value of NV from the input file and check if is too low.
  READ(11,*)NV
  IF(NV<=3)THEN
    PRINT *,'The NV parameter must be 4 or more on line 1 of file '
    PRINT *,TRIM(FILENAME)
    STOP
  ENDIF

  !Allocate space for the potential values.
  ALLOCATE(CPOTEN(NV),POTEN(2,NV))

  !Read in the first tabulated point on the potential energy curve.
  READ(11,*,IOSTAT=ERR)POTEN(1,1),POTEN(2,1)
  IF(ERR/=0)THEN
    PRINT *,'Error reading the first tabulated potential point &
      &from file '
    PRINT *,TRIM(FILENAME)
    STOP
  ENDIF
  WRITE(13,*)POTEN(1,1),POTEN(2,1)

  !Read in the rest of the points, making sure that R increases.
  DO I=2,NV
    READ(11,*,IOSTAT=ERR)POTEN(1,I),POTEN(2,I)
    IF(ERR/=0)THEN
      PRINT '(A,I3,A)','Error reading the ',I,'-th tabulated &
        &potential from file '
      PRINT *,TRIM(FILENAME)
      STOP
    ELSEIF(POTEN(1,I)<=POTEN(1,I-1))THEN
      PRINT *,'The separations do not increase monotonically at &
        &separation ',POTEN(1,I),' in file'
      PRINT *,TRIM(FILENAME)
      STOP
    ENDIF
    WRITE(13,*)POTEN(1,I),POTEN(2,I)
  ENDDO

  !Make sure that the values are positive at the 3 smallest separations.
  DO I=1,3
    IF(POTEN(2,I)<ZERO)THEN
      PRINT *,'In input file '//TRIM(FILENAME)
      PRINT '(A,I1,A)','the tabulated potential value ',I,&
        &' is negative.'
      PRINT *,'This is not allowed, so add more points at smaller'
      PRINT *,'separation before running this program again.'
      STOP
    ENDIF
  ENDDO

  !Make sure that the values are negative at the 3 largest separations
  !when NLONG is 4 or 6.
  IF(NLONG==4.OR.NLONG==6)THEN
    DO I=NV-2,NV
      IF(POTEN(2,I)>ZERO)THEN
        PRINT *,'In input file '//TRIM(FILENAME)
        PRINT *,'the tabulated potential value, ',I,', is &
          &positive while NLONG=',NLONG
        PRINT *,'This is not allowed, so add more points at larger'
        PRINT *,'separation before running this program again.'
        STOP
      ENDIF
    ENDDO
  ENDIF

  !Determine EXPON, CSHORT and CLONG.  Note the asymmetry in the way
  !  that the short- and long-range potentials are evaluated in VF.
  EXPON=LOG(POTEN(2,2)/POTEN(2,1))/LOG(POTEN(1,1)/POTEN(1,2))
  CSHORT=POTEN(2,1)
  CLONG=-POTEN(2,NV)*POTEN(1,NV)**NLONG

  !Write information to the output file.
  WRITE(13,'(A,/,A)')'Transport cross sections for ',TRIM(COMMENT)
  WRITE(13,'(A)')'The long-range potential is assumed to be:'
  IF(CLONG>ZERO)THEN
    WRITE(13,'(A,G15.8,A,I1,A)')'V(R)=-',CLONG,'*R**(-',NLONG,')'
  ELSE
    WRITE(13,'(A,G15.8,A,I1,A)')'V(R)=',-CLONG,'*R**(-',NLONG,')'
  ENDIF
  WRITE(13,'(A,I4,A)')'There are ',NV,' pairs of potential values &
    &entered.'
  WRITE(13,'(A)')'Listed below are the input values of R and V'
  WRITE(13,'(A)')'that are fit by a clamped cubic spline.'
  DO I=1,NV
    WRITE(13,'(2(1X,G22.14))')POTEN(1,I),POTEN(2,I)
  ENDDO

  WRITE(13,'(A,G15.8,A)')'Minimum energy is EMIN=',EMIN,' hartree.'
  WRITE(13,'(A,G15.8,A)')'Maximum energy is EMAX=',EMAX,' hartree.'
  WRITE(13,'(A,I2,A)')'The maximum L value is ',MAXL,'.'
  WRITE(13,'(A,G12.6,A)')'The relative accuracy is ',ACCURACY,'.'

  !Spline fit the potential.
  CALL SPLINE(NV,POTEN(1,1:NV),POTEN(2,1:NV),-POTEN(2,1)*EXPON/&
    &POTEN(1,1),-POTEN(2,NV)*NLONG/POTEN(1,NV),CPOTEN)
END SUBROUTINE INPUT

SUBROUTINE ORBITS
  !This subroutine computes the orbiting parameters.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the local real variables.
  !  R is the separation of interest.
  !  C1-C3 are used to check the orbiting conditions.
  !  VF is the value of the potential at R.
  !  VD is the first derivative of VF(R) at R.
  !  VD2 is the second derivative of VF(R) at R.
  !  VEFF is the effective potential at orbiting.
  REAL(SELECTED_REAL_KIND(18))::R,C1,C2,C3,VF,VD,VD2,VEFF

  !Declare the local integer variables.
  !  I is a loop index.
  !  J labels a particular set of orbiting parameters.
  !  K is used to retain the value of J while sets are added.
  INTEGER::I,J,K

  !Initialize global variables NO and ORBIT.  Note that NO is designed
  !to work down at first; later it will work up to its final value.
  NO=MAXORBITS+1
  ORBIT=ZERO

  !There are three criteria for orbiting.  The first is that VD(ro)>0.
  !The second is that VD2(ro)+3*VD(ro)/ro < 0.  The third is that
  !Veff>0.  If these are satisfied, then Eo=V(ro)+ro*VD(ro)/2 and
  !bo=SQRT(ro**3*VD(ro)/2/Eo).  Start the loop looking for orbiting
  !conditions at the first tabulated separation on the potential.
  R=POTEN(1,1)
  DO
    !At this point we have a designated value of r0, called R.  It was
    !either set before the loop or was reset at the end of the loop.
    C1=VD(R)
    C2=VD2(R)+THREE*C1/R
    C3=VF(R)
    VEFF=C3+R*C1/TWO

    !Check the orbiting conditions.
    IF(C1>ZERO.AND.C2<ZERO.AND.VEFF>ZERO)THEN
      !There is at least one orbiting region.  We'll look for more
      !later, after we store the set we have found.
      NO=NO-1
      IF(NO<=0)THEN
        PRINT *,'More room is needed for orbiting parameters'
        PRINT *,'for angle ',FILENAME2
        STOP
      ENDIF
      ORBIT(1,NO)=VEFF
      ORBIT(2,NO)=SQRT(R**3*C1/VEFF/TWO)
      ORBIT(3,NO)=R

      !Improve the guess by many changes in R.  When NO=MAXORBITS, we
      !must find the maximum orbiting energy to a very high degree of
      !accuracy.  For other values of NO, we use a resonably high degree
      !of accuraccy, checking to see if this is an energy where one
      !region of multiple orbiting ends.
      DO I=1,7
        DO
          R=(ONE-0.001**I)*ORBIT(3,NO)
          C1=VD(R)
          C2=VD2(R)+THREE*C1/R
          C3=VF(R)
          VEFF=C3+R*C1/TWO
          IF(C1>ZERO.AND.C2<ZERO.AND.VEFF>ORBIT(1,NO))THEN
            ORBIT(1,NO)=VEFF
            ORBIT(2,NO)=SQRT(R**3*C1/VEFF/TWO)
            ORBIT(3,NO)=R
            IF(NO==MAXORBITS)CYCLE
          ENDIF
          EXIT
        ENDDO
      ENDDO
    ENDIF

    !If there is at least one orbiting region, we exit when the
    !orbiting energy is less than EMIN but greater than 0.  The latter
    !condition arises because we initially set ORBIT=0, so a value of
    !0 indicates that no orbiting was found.
    IF(NO<=MAXORBITS)THEN
      IF(ORBIT(1,NO)<=EMIN.AND.ORBIT(1,NO)>ZERO)EXIT
    ENDIF

    !Move to a larger separation by making rather big steps (5%) beyond
    !the last tabulated separation but smaller steps (1%) otherwise.
    IF(R>POTEN(1,NV))THEN
      R=1.05*R
    ELSE
      R=1.01*R
    ENDIF

    IF(R>1.0D3)THEN
      IF(NO>MAXORBITS)THEN
        PRINT *,'Even at 1000 bohr, the potential did not support &
          &orbiting in file '
        PRINT *,FILENAME
        STOP
      ENDIF
      EXIT
    ENDIF
  ENDDO

  !Determine ED for potentials with NLONG=3 and CLONG>0.  We must use
  !very high degree of accuracy in order to avoid problems with
  !accidential singularities when performing integrals.
  IF(NLONG==3.AND.CLONG<ZERO)THEN
    ORBIT(:,1)=ORBIT(:,NO)
    DO I=1,7
      DO
        R=(ONE+0.001**I)*ORBIT(3,1)
        C1=VD(R)
        C2=VD2(R)+THREE*C1/R
        C3=VF(R)
        VEFF=C3+R*C1/TWO
        IF(C1>ZERO.AND.C2<ZERO.AND.VEFF<ORBIT(1,1))THEN
          ORBIT(1,1)=VEFF
          ORBIT(2,1)=SQRT(R**3*C1/VEFF/TWO)
          ORBIT(3,1)=R
        ENDIF
        EXIT
      ENDDO
    ENDDO
    !Since we have already finished the first orbiting triplet, set J=1
    !and increment NO.
    J=1
    NO=NO+1
    ED=ORBIT(1,1)
  ELSE
    !Since we have not considered the first orbiting triplet, set J=0
    !and ED=ZERO, but do not change NO.
    J=0
    ED=ZERO
  ENDIF

  !We now make a pass through the orbiting parameters in order to add
  !sets in regions where the energy is decreasing as the separation
  !decreases.
  DO
    J=J+1
    !Move the next set down from position NO to position J in ORBIT.
    ORBIT(:,J)=ORBIT(:,NO)

    !Reset NO so that it will be ready for the next pass.
    NO=NO+1
    IF(NO>MAXORBITS)EXIT

    !If J ever becomes larger than NO-19, we are in trouble since there
    !won't be enough room to add more sets.
    IF(J>NO-19)THEN
      PRINT *,'There is not enough room for the orbiting parameters &
        &in file'
      PRINT *,FILENAME
      STOP
    ENDIF

    !If the energies decrease, add up to 19 sets.
    IF(NO<MAXORBITS)THEN
      IF(ORBIT(1,J)>ORBIT(1,NO))THEN
        !Use K in order to retain the value of J before the additional
        !orbiting sets are added in a loop.  This means that we will
        !test separations that are linearly spaced between K and NO.
        K=J
        DO I=1,19
          R=(DBLE(20-I)*ORBIT(3,K)+DBLE(I)*ORBIT(3,NO))/TEN/TWO
          C1=VD(R)
          C2=VD2(R)+THREE*C1/R
          C3=VF(R)
          VEFF=C3+R*C1/TWO
          IF(C1>ZERO.AND.C2<ZERO.AND.VEFF>ORBIT(1,J))THEN
            J=J+1
            ORBIT(1,J)=VEFF
            ORBIT(2,J)=SQRT(R**3*C1/VEFF/TWO)
            ORBIT(3,J)=R
          ENDIF
        ENDDO
      ENDIF
    ELSE
      IF(ORBIT(1,J)<ORBIT(1,J-1))THEN
        DO I=1,19
          R=ORBIT(3,J)+DBLE(I)/TEN/TWO*(ORBIT(3,J)-ORBIT(3,J-1))
          C1=VD(R)
          C2=VD2(R)+THREE*C1/R
          C3=VF(R)
          VEFF=C3+R*C1/TWO
          IF(C1>ZERO.AND.C2<ZERO.AND.VEFF>ORBIT(1,J))THEN
            J=J+1
            ORBIT(1,J)=VEFF
            ORBIT(2,J)=SQRT(R**3*C1/VEFF/TWO)
            ORBIT(3,J)=R
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ENDDO  

  !Chop off any points where the energies have yet to start going up.
  DO
    IF(ORBIT(1,J)<ORBIT(1,J-1))THEN
      J=J-1
    ELSE
      EXIT
    ENDIF
  ENDDO

  !There are J sets of orbiting parameters, so put J into NO and
  !write the parameters to the output file.
  NO=J
  WRITE(13,'(/,A,I3,A)')' The ',NO,' sets of orbiting parameters are:'
  WRITE(13,'(3A)')'          E                     B&
    &                     RM'
  WRITE(13,'(3(1X,G22.14E4))')((ORBIT(J,I),J=1,3),I=1,NO)
END SUBROUTINE ORBITS

REAL(SELECTED_REAL_KIND(18)) FUNCTION VF(R)
  !This function evaluates the interaction potential, which is assumed
  !to be repulsive at short separations R.  It is also assumed that VF
  !smoothly connects to the proper long-range (asymptotic) tail.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the dummy variable.
  !  R is the ion-neutral separation, in bohr.
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::R

  !Declare the local real variables.
  !  FIT is an external function for the spline fit.
  REAL(SELECTED_REAL_KIND(18))::FIT

  !First, a possible short-range extrapolation must be considered.
  IF(R<POTEN(1,1))THEN
    VF=CSHORT*(POTEN(1,1)/R)**EXPON
    RETURN

  !Second, we consider a possible long-range extrapolation.
  ELSEIF(R>POTEN(1,NV))THEN
    VF=-CLONG/R**NLONG
    RETURN

  !Get the potential from the spline fit of the interpolation function
  !and change from the interpolation function to the potential.
  ELSE
    VF=FIT(NV,POTEN(1,1:NV),POTEN(2,1:NV),CPOTEN(1:NV),R)
  ENDIF
END FUNCTION VF

REAL(SELECTED_REAL_KIND(18)) FUNCTION VD(R)
  !This function evaluates the derivative of the potential energy
  !with respect to separation.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the dummy variable.
  !  R is the ion-neutral separation, in bohr.
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::R

  !Declare the local real variables.
  !  FIT2 is an external function for the spline fit.
  REAL(SELECTED_REAL_KIND(18))::FIT2

  !First, a possible short-range extrapolation must be considered.
  IF(R<POTEN(1,1))THEN
    VD=-CSHORT*EXPON/POTEN(1,1)*(POTEN(1,1)/R)**(EXPON+ONE)
    RETURN

  !Second, we consider a possible long-range extrapolation.
  ELSEIF(R>POTEN(1,NV))THEN
    VD=CLONG*NLONG/R**(NLONG+ONE)
    RETURN

  !What remains is interpolation in tabulated values.
  ELSE
    VD=FIT2(NV,POTEN(1,1:NV),POTEN(2,1:NV),CPOTEN(1:NV),R)
  ENDIF
END FUNCTION VD

REAL(SELECTED_REAL_KIND(18)) FUNCTION VD2(R)
  !This function evaluates the second derivative of the potential
  !energy with respect to separation.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the dummy variable.
  !  R is the ion-neutral separation, in bohr.
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::R

  !Declare the local real variables.
  !  FIT3 is an external function for the spline fit.
  REAL(SELECTED_REAL_KIND(18))::FIT3

  !First, a possible short-range extrapolation must be considered.
  IF(R<POTEN(1,1))THEN
    VD2=CSHORT*EXPON*(EXPON+ONE)/POTEN(1,1)**2*(POTEN(1,1)/R)**&
      &(EXPON+TWO)
    RETURN

  !Second, we consider a possible long-range extrapolation.
  ELSEIF(R>POTEN(1,NV))THEN
    VD2=-CLONG*NLONG*(NLONG+ONE)/R**(NLONG+TWO)
    RETURN

  !What remains is interpolation in tabulated values.
  ELSE
    VD2=FIT3(NV,POTEN(1,1:NV),POTEN(2,1:NV),CPOTEN(1:NV),R)
  ENDIF
END FUNCTION VD2

SUBROUTINE SPLINE(N,X,Y,YP1,YPN,Y2)
  !This subroutine is based on Sec. 3.3 of "Numerical Recipes" (FORTRAN
  !VERSION) by W. H. Press, B. P. Flannery, S. A. Teukolsky and W. T.
  !Vetterling (Cambridge University Press, Cambridge, England, 1989).
  !It computes the second derivative matrix Y2 that is needed for a
  !spline fit using function FIT of the N points in arrays X and Y,
  !with Y as a function of X.  Since the fit will generally be a
  !clamped spline, the derivatives at the endpoints must be supplied in
  !YP1 and YPN; if the values in YP1 and YPN exceed 0.99D30 then a
  !natural spline fit will be used.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the dummy variables explained above.
  INTEGER,INTENT(IN)::N
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(N),INTENT(IN)::X,Y
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::YP1,YPN
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(N),INTENT(OUT)::Y2

  !Declare the local variables.
  !  U, SIG, P, QN and UN are used for temporary storage.
  !  I is a loop index.
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(N)::U
  REAL(SELECTED_REAL_KIND(18))::SIG,P,QN,UN
  INTEGER::I

  !Iniitialize the lower end values.
  IF(YP1>0.99D30)THEN
    Y2(1)=ZERO
    U(1)=ZERO
  ELSE
    Y2(1)=-ONE/TWO
    U(1)=THREE/(X(2)-X(1))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
  ENDIF

  !Decompose the tridiagonal matrix.
  DO I=2,N-1
    SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
    P=SIG*Y2(I-1)+TWO
    Y2(I)=(SIG-ONE)/P
    U(I)=(THREE*TWO*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))/&
      &(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
  ENDDO

  !Initialize the upper end values.
  IF(YPN>0.99D30)THEN
    QN=ZERO
    UN=ZERO
  ELSE
    QN=ONE/TWO
    UN=(THREE/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
  ENDIF

  !Compute the derivative matrix by back-substitution.
  Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+ONE)
  DO I=N-1,1,-1
    Y2(I)=Y2(I)*Y2(I+1)+U(I)
  ENDDO
END SUBROUTINE SPLINE

REAL(SELECTED_REAL_KIND(18)) FUNCTION FIT(N,X,Y,Y2,XBAR)
  !This subroutine is based on Sec. 3.3 of "Numerical Recipes" (FORTRAN
  !VERSION) by W. H. Press, B. P. Flannery, S. A. Teukolsky and W. T.
  !Vetterling (Cambridge University Press, Cambridge, England, 1989).
  !It computes the spline-fit function (See subroutine SPLINE).
  !See the comments in subroutine SPLINE.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the dummy variables explained in subroutine SPLINE.
  INTEGER,INTENT(IN)::N
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(N),INTENT(IN)::X,Y,Y2
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::XBAR

  !Declare the local variables.
  INTEGER::KLOW,KHI,KMID
  REAL(SELECTED_REAL_KIND(18))::XH,XA,XB

  !Check for the case N=1.
  IF(N<=1)THEN
    PRINT *,'Subroutine FIT was called with N=',N,' for angle'&
      &//FILENAME2
    STOP
  ENDIF

  !Search for the correct subinterval containing XBAR.
  KLOW=1
  KHI=N
  DO
    IF(KHI-KLOW<=1)EXIT
    KMID=(KHI+KLOW)/TWO
    IF(X(KMID)>XBAR)THEN
      KHI=KMID
    ELSE
      KLOW=KMID
    ENDIF
  ENDDO

  !Evaluate the appropriate piece of Y(X).
  XH=X(KHI)-X(KLOW)
  XA=(X(KHI)-XBAR)/XH
  XB=(XBAR-X(KLOW))/XH
  FIT=XA*Y(KLOW)+XB*Y(KHI)+(XA*(XA*XA-ONE)*Y2(KLOW)+XB*(XB*XB-ONE)*&
    &Y2(KHI))*XH*XH/THREE/TWO
END FUNCTION FIT

REAL(SELECTED_REAL_KIND(18)) FUNCTION FIT2(N,X,Y,Y2,XBAR)
  !This subroutine is based on Sec. 3.3 of "Numerical Recipes" (FORTRAN
  !VERSION) by W. H. Press, B. P. Flannery, S. A. Teukolsky and W. T.
  !Vetterling (Cambridge University Press, Cambridge, England, 1989).
  !It computes the derivative of the spline-fit function.
  !See the comments in subroutine SPLINE.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the dummy variables explained in subroutine SPLINE.
  INTEGER,INTENT(IN)::N
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(N),INTENT(IN)::X,Y,Y2
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::XBAR

  !Declare the local variables.
  INTEGER::KLOW,KHI,KMID
  REAL(SELECTED_REAL_KIND(18))::XH,XA,XB

  !Check for the case N=1.
  IF(N<=1)THEN
    PRINT *,'Subroutine FIT2 was called with N=',N,' for angle ',&
      &FILENAME2
    STOP
  ENDIF

  !Search for the correct subinterval containing XBAR.
  KLOW=1
  KHI=N
  DO
    IF(KHI-KLOW<=1)EXIT
    KMID=(KHI+KLOW)/TWO
    IF(X(KMID)>XBAR)THEN
      KHI=KMID
    ELSE
      KLOW=KMID
    ENDIF
  ENDDO

  !Evaluate the appropriate piece of Y'(X).
  XH=X(KHI)-X(KLOW)
  XA=(X(KHI)-XBAR)/XH
  XB=(XBAR-X(KLOW))/XH
  FIT2=(Y(KHI)-Y(KLOW))/XH+XH/THREE/TWO*&
    &(-(THREE*XA*XA-ONE)*Y2(KLOW)+(THREE*XB*XB-ONE)*Y2(KHI))
END FUNCTION FIT2

REAL(SELECTED_REAL_KIND(18)) FUNCTION FIT3(N,X,Y,Y2,XBAR)
  !This subroutine is based on Sec. 3.3 of "Numerical Recipes" (FORTRAN
  !VERSION) by W. H. Press, B. P. Flannery, S. A. Teukolsky and W. T.
  !Vetterling (Cambridge University Press, Cambridge, England, 1989).
  !It computes the second derivative of the spline-fit function.
  !See the comments in subroutine SPLINE.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the dummy variables explained in subroutine SPLINE.
  INTEGER,INTENT(IN)::N
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(N),INTENT(IN)::X,Y,Y2
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::XBAR

  !Declare the local variables.
  INTEGER::KLOW,KHI,KMID
  REAL(SELECTED_REAL_KIND(18))::XH,XA,XB

  !Check for the case N=1.
  IF(N<=1)THEN
    PRINT *,'Subroutine FIT3 was called with N=',N,' for angle ',&
      &FILENAME2
    STOP
  ENDIF

  !Search for the correct subinterval containing XBAR.
  KLOW=1
  KHI=N
  DO
    IF(KHI-KLOW<=1)EXIT
    KMID=(KHI+KLOW)/TWO
    IF(X(KMID)>XBAR)THEN
      KHI=KMID
    ELSE
      KLOW=KMID
    ENDIF
  ENDDO

  !The next line does nothing, but it keeps the compiler from sending
  !a warning message about Y not being used.  Y is used as a dummy
  !argument simply to keep FIT3 parallel to FIT2.
  XH=Y(1)

  !Evaluate the appropriate piece of Y''(X).
  XH=X(KHI)-X(KLOW)
  XA=(X(KHI)-XBAR)/XH
  XB=(XBAR-X(KLOW))/XH
  FIT3=XA*Y2(KLOW)+XB*Y2(KHI)
END FUNCTION FIT3

SUBROUTINE REGIONS
  !This subroutine determine the orbiting regions, the collections of
  !orbiting parameters in which the orbiting energy always increases.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the local integer variable.
  !  J is a loop index.
  INTEGER::J

  !Declare the local logical variable.
  !  INCREASING is true if energies are rising.
  LOGICAL::INCREASING

  !Write the line that heads this output region.
  WRITE(13,'(A)')' The orbiting energy regions are '

  !Set the first energy in the first region, and the pointer that
  !indicates the corresponding number of the orbiting set.
  INOO=1
  ECC(1,1)=ORBIT(1,1)
  NOO(1,1)=1

  !Initially, the energies must be increasing.
  INCREASING=.TRUE.

  !Use a loop to look for the end of each region and the start of the
  !next.
  DO J=2,NO
    !Treat the situation where the energies quit increasing.
    IF(INCREASING.AND.ORBIT(1,J)<ORBIT(1,J-1))THEN
      !Place the ending energy into ECC(2,INOO) and the number of the
      !corresponding orbiting set into NOO(2,INOO)
      ECC(2,INOO)=ORBIT(1,J-1)
      NOO(2,INOO)=J-1

      !Write information about the region.
      WRITE(13,*)ECC(1,INOO),'-',ECC(2,INOO),' Region ',INOO

      !Get ready for the next region.
      INOO=INOO+1
      IF(INOO>INOOMAX)THEN
        PRINT *,'INOOMAX is too small for angle ',FILENAME
        PRINT *,'You must reset this named constant, recompile the '
        PRINT *,'program and try again.'
        STOP
      ENDIF
      INCREASING=.FALSE.

    !Treat the situation where the energies resume increasing.
    ELSEIF(.NOT.INCREASING.AND.ORBIT(1,J)>ORBIT(1,J-1))THEN
      ECC(1,INOO)=ORBIT(1,J-1)
      NOO(1,INOO)=J-1
      INCREASING=.TRUE.
    ENDIF
  ENDDO

  !Place the ending energy into ECC(2,INOO) and the number of the
  !corresponding orbiting set into NOO(2,INOO)
  ECC(2,INOO)=ORBIT(1,NO)
  NOO(2,INOO)=NO

  !Write information about the final region.
  WRITE(13,*)ECC(1,INOO),'-',ECC(2,INOO),' Region ',INOO
  WRITE(13,*)
END SUBROUTINE REGIONS

SUBROUTINE EFIT
  !This subroutine orders the calculation of the cross sections in three
  !regions of energy: less than EC, the maximum orbiting energy; from
  !EC to 10*EC; and above 10*EC.  It also calculates the Chebyshev
  !coefficients that describe them.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the local integer variables.
  !  A non-zero value of NQEC signals that endpoints must be matched
  !    between regions.  
  !  INDEX is a loop index indicating one of the three energy regions.
  !  NSTAR=1 if we must calculate values at the first energy in the
  !    region, but NSTAR=2 if we have already calculated them as the
  !    last energy in the previous region.
  !  INM is 1 when we must calculate cross sections at every energy.
  !    INM is 2 when we can avoid recalculating values.
  !  NM is the number of energies in the region (5, 9, 17 ... MAXE)
  !  I, J and L are loop indices.
  !  NPRINT is the number of sets of cross sections to be printed.
  !  JSTRT is the starting value for a loop indexed by J.
  INTEGER::NQEC,INDEX,NSTAR,INM,NM,I,J,L,NPRINT,JSTRT

  !Declare the local integer array.
  !  NIN holds the number of Chebyshev coefficients for each energy
  !     region.
  INTEGER,DIMENSION(3)::NIN

  !Declare the local real variables.
  !  E1 and E2 are the minimum and maximum energies in the region.
  !  E contains the energy being considered.
  !  E2P1 and E2M1 are used to find E, on a logarithmic scale,
  !    between E1 and E2.
  !  XCHEB, TJ, ACHEB and CCHEB are used to determine Chebyshev
  !    coefficients.
  REAL(SELECTED_REAL_KIND(18))::E1,E2,E,E2P1,E2M1,XCHEB,TJ,ACHEB,&
    CCHEB

  !Declare the local real arrays.
  !  Array QEC contains up to MAXL cross sections at an energy that
  !     is at the bottom of one region and the top of another.
  !  Array COEFF contains preliminary values for QEC.
  !  Array EST contains up to MAXL estimates of the accuracies at one
  !     energy.
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(MAXL)::QEC,COEFF,EST

  !Determine EC and EC2.
  EC=MAXVAL(ORBIT(1,:))
  EC2=TEN*EC

  !Initialize the variables.
  NQEC=0
  NIN=0

  !Start a loop over the three energy regions.
  DO INDEX=1,3
    !Set E1 and E2 to the minimum and maximum energies of the region.
    !We must be careful, in case EMIN and/or EMAX prevent us from
    !using three energy regions.
    IF(INDEX==1)THEN
      IF(EMIN>=EC)CYCLE
      E1=EMIN
      E2=MIN(EC,EMAX)
    ELSEIF(INDEX==2)THEN
      IF(EMAX<=EC)EXIT
      IF(EMIN>=EC2)CYCLE
      E1=MAX(EC,EMIN)
      E2=MIN(EC2,EMAX)
    ELSE
      IF(EMAX<=EC2)EXIT
      E1=MAX(EC2,EMIN)
      E2=EMAX
    ENDIF

    !NQEC will be equal to 0 for the first valid energy region, as
    !there are no previous calculations to match endpoints with.
    !This means we must start at the first possible energy (NSTAR=1)
    !and make single steps (INM=1).
    IF(NQEC==0)THEN
      NSTAR=1
      INM=1

    !If the endpoints must match, write the stored values E and QEC
    !from the bottom of this loop.  This means that we can start at the
    !second possible energy (NSTAR=2) but must use single steps (INM=1).
    ELSE
      WRITE(16,REC=MAXE*(INDEX-1)+1)E,QEC
      NSTAR=2
      INM=1
    ENDIF

    !Establish the parameters for the mimimum number (NM=5) of Chebyshev
    !coefficients.
    NM=5
    E1=LOG10(E1)
    E2=LOG10(E2)
    E2P1=E2+E1
    E2M1=E2-E1
    
    !Start a loop to use 5, 9, 17, ... MAXE values.
    DO
      !Use a loop to add cross sections at the missing energies.
      DO I=NSTAR,NM,INM
        !Determine the next energy, first on a logarthmic scale and then
        !converted to the actual energy in hartree.  Since NM-1 is
        !even, it doesn't matter what sign precedes E2M1 below; using
        !the negative sign makes the energies go up.
        E=(E2P1-E2M1*COS((I-1)*PI/(NM-1)))/TWO
        E=TEN**E

        !Call CRSCT to place the cross sections in COEFF.
        CALL CRSCT(E,COEFF)

        !Write out the cross sections as they are calculated.
        !Note that the cross sections are not being saved,
        !so they will need to be read again when needed below.
        WRITE(13,'(1X,G25.17)')E
        WRITE(13,'(1X,5G15.7)')COEFF
        WRITE(16,REC=MAXE*(INDEX-1)+I)E,COEFF
      ENDDO

      !Begin to compute the Chebychev coefficients for file 17
      !by putting zero into the unused parts of file 17.
      COEFF=ZERO
      DO I=MAXE*(INDEX-1)+1,195
        WRITE(17,REC=I)COEFF
      END DO

      !For each term in file 17, first read the cross section.
      DO J=1,NM
        READ(17,REC=MAXE*(INDEX-1)+J)COEFF

        !Read and use the cross sections in a loop.
        DO I=1,NM
          READ(16,REC=MAXE*(INDEX-1)+I)E,QEC
          XCHEB=-COS((I-1)*PI/(NM-1))
          TJ=XCHEB
          IF(J==1)THEN
            TJ=ONE
          ELSEIF(J>2)THEN
            ACHEB=ONE
            DO L=3,J
              CCHEB=TWO*XCHEB*TJ-ACHEB
              ACHEB=TJ
              TJ=CCHEB
            END DO
          END IF
          IF(I==1.OR.I==NM)THEN
            COEFF=COEFF+TJ*LOG(QEC)/(NM-1)
          ELSE
            COEFF=COEFF+TJ*LOG(QEC)/(NM-1)*TWO
          END IF
          WRITE(17,REC=MAXE*(INDEX-1)+J)COEFF
        END DO
      END DO

      !Interpolations using file 17 will be involved in integrations
      !to form the collision integrals.  Estimate the error now by
      !using the "empirical" estimator of H. O'Hara and F. J. Smith,
      !J. Computat. Phys. 5 (1970) 328.
      READ(17,REC=MAXE*(INDEX-1)+NM)EST
      READ(17,REC=MAXE*(INDEX-1)+NM-1)COEFF
      EST=(ABS(EST)+ABS(COEFF))*TWO/(NM-1)

      !If convergence has been achieved, or if the maximum number of
      !energies have been sampled, write the appropriate quantities.
      IF(MAXVAL(EST)<=ACCURACY.OR.NM==MAXE)THEN
        IF(MAXVAL(EST)>ACCURACY/0.8D0)THEN
          PRINT *,' Convergence failure in region ',INDEX
          PRINT *,' restricts the accuracy to ',MAXVAL(EST)
          STOP
        ENDIF
        NIN(INDEX)=NM
        EXIT

      !Otherwise, double the number of energies to sample and continue
      !the calculation in the do loop, making double steps (INM=2).
      ELSE
        DO I=NM,2,-1
          READ(16,REC=MAXE*(INDEX-1)+I)E,QEC
          WRITE(16,REC=MAXE*(INDEX-1)+2*I-1)E,QEC
        ENDDO
        NM=2*NM-1
        NSTAR=2
        INM=2
      ENDIF
    ENDDO
    NQEC=1
    READ(16,REC=MAXE*(INDEX-1)+NM)E,QEC
  ENDDO

  !Write the cross sections into file 15 in such manner that the
  !energies increase monotonically.
  NPRINT=NIN(1)+NIN(2)+NIN(3)
  IF(NIN(1)/=0.AND.NIN(2)/=0)NPRINT=NPRINT-1
  IF(NIN(2)/=0.AND.NIN(3)/=0)NPRINT=NPRINT-1
  WRITE(15,*)NPRINT
  DO I=1,3
    IF(NIN(I)==0)CYCLE
    IF((I==2.AND.NIN(1)/=0).OR.(I==3.AND.NIN(2)/=0))THEN
      JSTRT=2
    ELSE
      JSTRT=1
    ENDIF
    DO J=JSTRT,NIN(I)
      READ(16,REC=J+MAXE*(I-1))E,QEC
      WRITE(15,'(D15.7)')E
      WRITE(15,'(5F15.6)')QEC
    ENDDO
  ENDDO

  !Write NIN.
  DO I=1,3
    WRITE(18,*)NIN(I)
  ENDDO
END SUBROUTINE EFIT

SUBROUTINE CRSCT(E,X1)
  !This subroutine calculates the transport cross sections at energy
  !E and puts them in array X1.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the dummy arguments.
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::E
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(MAXL),INTENT(OUT)::X1

  !Declare the local real arrys.
  !  EXX, BXX and RXX are the arrays used for interpolations to find
  !    global arrays BO and RO.
  !  CEB contains the coefficients to spline fit orbiting impact
  !    parameters as a function of orbiting energies.
  !  CBR contains the coefficients to spline fit orbiting separations
  !    as a function of orbiting impact parameters.
  !  EA and EB hold the integrand values at the endpoints A and B.
  !  X2 and X3 temporarily hold the value of X1.
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(MAXORBITS)::EXX,BXX,RXX,CEB,CBR
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(MAXL)::EA,EB,X2,X3

  !Declare the local real variables.
  !  BCC is a local copy of global variable ACC1.
  !  FA and FB are the endpoints of the integral.
  !  FLAG is used to signal a natural cubic spline
  !  FIT is the external function used in spline fits.
  !  BC is the orbiting impact parameter.
  !  RM is the orbiting separation.
  !  VF is the potential function.
  REAL(SELECTED_REAL_KIND(18))::BCC,FA,FB,FLAG,FIT,BC,RM,VF

  !Declare the external subroutine.
  EXTERNAL::QINT

  !Declare the local integers.
  !  LMAX is a local copy of global variable MAXL.
  !  IPR is a flag for lack of convergence.
  !  I and J are loop indices.
  INTEGER::LMAX,IPR,I,J

  !We cannot move dummy variables from one subroutine to another, so we
  !must initialize the global variable EE that can be used instead of E
  !in the subroutines.
  EE=E

  !Subroutine INTEG was written so that variables are used to tell it
  !how many functions are being integrated, the desired fractional
  !accuracy, the endpoints, and what to do if there is a convergence
  !error.  We need to initialize variables to be used this way.
  LMAX=MAXL
  BCC=ACC1
  FA=-ONE
  FB=ONE
  IPR=0

  !Initialize X1 and FLAG.
  X1=ZERO
  FLAG=1.0D30
  NAITK=0

  !Cycle through the potential to find the impact parameters and
  !orbiting separations when E is low.  Note that round-off error may
  !make E slightly smaller or slightly larger than EC even when we think
  !it is exactly equal to EC.  The different parts of the next IF
  !statement will sort the calculations correctly, as will QINT.
  IF(E>=ED.AND.E<=EC)THEN
    !Find the (possibly multiple) orbiting impact parameters and
    !separations corresponding to this energy. 
    DO I=1,INOO
      !If the energy does not fall into this particular range of
      !orbiting energies, set BO and RO to zero.
      IF(E<ECC(1,I).OR.E>ECC(2,I))THEN
        BO(I)=ZERO
        RO(I)=ZERO

      !IF E is exactly equal to the lowest orbiting energy in this 
      !region, then we know the value of BO(I) and RO(I).
      ELSEIF(E<=ECC(1,I))THEN
        BO(I)=ORBIT(2,NOO(1,I))
        RO(I)=ORBIT(3,NOO(1,I))

      !IF E is exactly equal to the highest orbiting energy in this
      !region, then we know the value of BO(I) and RO(I).
      ELSEIF(E>=ECC(2,I))THEN
        BO(I)=ORBIT(2,NOO(2,I))
        RO(I)=ORBIT(3,NOO(2,I))

      !When E is within this region, then we must interpolate.
      !However, we must not use too many points from ORBIT because
      !subroutine AITKE can give results that oscillate badly.
      ELSE
        !Add the values to the arrays that will be used for
        !spline fits.
        DO J=0,NOO(2,I)-NOO(1,I)
          EXX(J+1)=ORBIT(1,NOO(1,I)+J)
          BXX(J+1)=ORBIT(2,NOO(1,I)+J)
          RXX(J+1)=ORBIT(3,NOO(1,I)+J)
        ENDDO

        !Spline fit BXX as a function of EXX and determine BO(I).
        CALL SPLINE(NOO(2,I)-NOO(1,I)+1,EXX,BXX,FLAG,FLAG,CEB)
        BO(I)=FIT(NOO(2,I)-NOO(1,I)+1,EXX,BXX,CEB,E)

        !Spline fit RXX as a function of EXX and determine RO(I).
        CALL SPLINE(NOO(2,I)-NOO(1,I)+1,EXX,RXX,FLAG,FLAG,CBR)
        RO(I)=FIT(NOO(2,I)-NOO(1,I)+1,EXX,BXX,CBR,E)
      ENDIF
      IF(RO(I)<ZERO)THEN
        PRINT *,'Unusual behavior at point 1 in CRSECT for angle',&
          &FILENAME2
        PRINT *,'You need to change the potential at large '
        PRINT *,'separations, increase the value for the minimum'
        PRINT *,'energy, or lower the accuracy requested.'
        STOP
      ENDIF
    ENDDO

    !Set NAITK to the number of valid entries in arrays BO and RO.
    NAITK=INOO

    !Eliminate any empty sets from the arrays B0 and RO.  At the same
    !time, adjust BO to exactly match RO; this avoids having a small
    !error in the orbiting separation make a large error in the
    !orbiting impact parameter and thus in the scattering angle.
    I=1
    DO
      IF(RO(I)>ZERO)THEN
        BC=ONE-VF(RO(I))/E
        IF(BC>=ZERO)BO(I)=RO(I)*SQRT(BC)
        I=I+1
        IF(I>NAITK)EXIT
      ELSEIF(I==NAITK)THEN
        NAITK=NAITK-1
        EXIT
      ELSE
        BO(I)=BO(NAITK)
        RO(I)=RO(NAITK)
        NAITK=NAITK-1
      ENDIF
    ENDDO

    !Sort the sets of orbiting parameters so that the impact parameters
    !always increase.  Note that BC and RM are used here only for 
    !temporary storage.
    DO I=1,NAITK
      DO J=2,NAITK
        IF(BO(I)>BO(J))THEN
          BC=BO(I)
          BO(I)=BO(J)
          BO(J)=BC
          RM=RO(I)
          RO(I)=RO(J)
          RO(J)=RM
        ENDIF
      ENDDO
    ENDDO

    !Set the values for ROMAX and BOMAX that will be used when E>EC.
    ROMAX=MAXVAL(RO)
    BOMAX=MAXVAL(BO)
  ELSE
    ROMAX=ORBIT(3,NO)
    BOMAX=ORBIT(2,NO)
  ENDIF

  !Evaluate the integrand at the endpoints.
  CALL QINT(FA,EA,LMAX)
  CALL QINT(FB,EB,LMAX)

  !Evaluate the integral.
  CALL INTEG(LMAX,QINT,FA,FB,EA,EB,X2,BCC,IPR,X3)
  IF(IPR/=0)THEN
    PRINT *,' Convergence failure with energy ',EE,' and angle'
    PRINT *,FILENAME2//FILENAME3//FILENAME4,&
      ' restricts the accuracy to ',MAXVAL(ABS(ONE-X3/X2))
  ENDIF
  X1=X2
END SUBROUTINE CRSCT

SUBROUTINE QINT(Y,FUN,LMAX)
  !This subroutine provides the integrand for CRSCT.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the dummy variables.
  !  LMAX is the number of functions.
  !  Y is the value of integration variable.
  !  FUN is the value of the integrands.
  INTEGER,INTENT(IN)::LMAX
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::Y
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(LMAX),INTENT(OUT)::FUN

  !Declare the local real variables.
  !  B2 is an orbiting impact parameter.
  !  R1, R2, R3, R4 and R5 are orbiting separations.
  !  RBAR is a separation at which the integral is split.
  !  TEMP is used for the parts of the cross sections that are
  !    independent of the scattering angle.
  !  X11 and X12 are values of X1 and X2 at endpoints.
  !  X1, X2 and X3 result from changes of variable for Y.
  !  VF is the potential function.
  !  VD is the first derivative of VF(R) at R.
  !  ANS stores part of FUN.
  REAL(SELECTED_REAL_KIND(18))::B2,R1,R2,R3,R4,R5,RBAR,TEMP,X11,&
    X12,X1,X2,X3,VF,VD
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(LMAX)::ANS

  !Declare the local integer variables.
  !  I and L are loop indices.
  INTEGER::I,L
  
  !First, examine the integrand if EE is such that orbiting can occur.
  !Suboutine CRSCT organized the orbiting parameters for this energy.
  IF(EE>=ED.AND.EE<=EC)THEN
    !The first terms are for the smallest orbiting impact parameter.
    IF(Y<=-ONE.OR.Y>=ONE)THEN
      FUN=ZERO
    ELSE
      B2=BO(1)*COS(PI*(Y+ONE)/FOUR)
      TEMP=(HPI*BO(1))**2*SIN(HPI*(Y+ONE))
      CALL EofB(B2,TEMP,FUN)
    ENDIF
  
    !The next terms are for the intermediate orbiting impact parameters.
    !There are two pieces because we split the integral at RBAR.  Note
    !the use of MIN and MAX, since in rare cases the values in RO do
    !not increase monotonically.
    IF(Y<ONE)THEN
      DO I=1,NAITK-1
        RBAR=(RO(I)+RO(I+1))/TWO
        X11=FOUR/PI*ASIN(MIN(RO(I),RO(I+1))/RBAR)-ONE
        X1=((ONE-X11)*Y+X11+ONE)/TWO
        R1=RBAR*SIN(PI*(X1+ONE)/FOUR)
        B2=ONE-VF(R1)/EE
        IF(B2<ZERO)THEN
          PRINT *,'Error at position 1 in QINT for angle'//FILENAME2
          PRINT *,'You need to change the number of potential values'
          PRINT *,'at large separations, increase the value for the'
          PRINT *,'minimum energy, or lower the accuracy in PC.in.'
          STOP
        ENDIF
        B2=R1*SQRT(B2)
        TEMP=(HPI*RBAR)**2/TWO*(ONE-X11)*SIN(HPI*(X1+ONE))*&
          &(ONE-VF(R1)/EE-R1/TWO*VD(R1)/EE)
        CALL EofB(B2,TEMP,ANS)
        !Whether ANS is added to or subtracted from FUN depends on the
        !relative sizes of RO(I) and RO(I+1).
        IF(RO(I)<RO(I+1))THEN
          FUN=FUN+ANS
        ELSE
          FUN=FUN-ANS
        ENDIF
      ENDDO
    ENDIF
    IF(Y>-ONE)THEN
      DO I=1,NAITK-1
        RBAR=(RO(I)+RO(I+1))/TWO
        IF(RBAR/MAX(RO(I),RO(I+1))<-ONE.OR.RBAR/MAX(RO(I),RO(I+1))>ONE)THEN
          PRINT *,'AARON ERROR: BAD ANGLE IN ACOS (line #1461)'
          PRINT *,'AARON ERROE: Tried to take ACOS of ', RBAR/MAX(RO(I),RO(I+1))
          STOP
        ENDIF
        X12=FOUR/PI*ACOS(RBAR/MAX(RO(I),RO(I+1)))-ONE
        X2=((ONE+X12)*Y+X12-ONE)/TWO
        R2=MAX(RO(I),RO(I+1))*COS(PI*(X2+ONE)/FOUR)
        B2=ONE-VF(R2)/EE
        IF(B2<ZERO)THEN
          PRINT *,'Error at position 2 in QINT and angle',FILENAME2
          PRINT *,'You need to change the number of potential values'
          PRINT *,'at large separations, increase the value for the'
          PRINT *,'minimum energy, or lower the accuracy required.'
          STOP
        ENDIF
        B2=R2*SQRT(B2)
        TEMP=(HPI*MAX(RO(I),RO(I+1)))**2/TWO*(ONE+X12)*SIN(HPI*(X2+&
          &ONE))*(ONE-VF(R2)/EE-R2/TWO*VD(R2)/EE)
        CALL EofB(B2,TEMP,ANS)
        !Whether ANS is added to or subtracted from FUN depends on the
        !relative sizes of RO(I) and RO(I+1).
        IF(RO(I)<RO(I+1))THEN
          FUN=FUN+ANS
        ELSE
          FUN=FUN-ANS
        ENDIF
      ENDDO
    ENDIF

    !The final terms are for large orbiting impact parameters.
    IF(Y>-ONE.AND.Y<ONE)THEN
      X3=(ONE+Y)/TWO
      R3=SIN(HPI*X3)
      IF(R3>ZERO)THEN
        R3=RO(NAITK)/R3
        B2=ONE-VF(R3)/EE
        IF(B2>ZERO)THEN
          B2=R3*SQRT(B2)
          TEMP=(HPI*RO(NAITK))**2*TWO*(R3/RO(NAITK))**3*COS(HPI*X3)*&
            &(ONE-VF(R3)/EE-R3/TWO*VD(R3)/EE)
          CALL EofB(B2,TEMP,ANS)
          FUN=FUN+ANS
        ENDIF
      ENDIF
    ENDIF

  !Second, examine the integrand if EC<EE<EC2.
  ELSEIF(EE>ED.AND.EE<=EC2)THEN
    !Estimate the root of eq. (48) in the text; call it R1 here,
    !although it was designated ~Rc in the manuscript.  A first estimate
    !of R1 is ROMAX, but ROMAX may not have a value if EMIN>EC.  Hence
    !we make a wild guess when ROMAX<=0.
    IF(ROMAX>ZERO)THEN
      R1=ROMAX
    ELSE
      R1=POTEN(1,NV)
    ENDIF
  
    !Call MINIM to precisely determine this root.
    CALL MINIM(EE,ROMAX,R1)

    !Determine the root of eq. (49) in the text; call it R2 here,
    !although it was designated ~Ro in the manuscript.
    B2=ZERO
    R2=R1
    CALL MINIM(EE,B2,R2)

    !The expressions are different if Y=-1, Y+1, or -1<Y<1.
    IF(Y<=-ONE)THEN
      FUN=ZERO
      DO L=2,LMAX,2
        FUN(L)=-PI*R2**2*(R1-R2)*VD(R2)/EE
      ENDDO
    ELSEIF(Y>=ONE)THEN
      TEMP=PI*R1*(R1-R2)*(ONE-VF(R1)/EE-R1/TWO*VD(R1)/EE)
      CALL EofB(ROMAX,TEMP,FUN)
    ELSE
      !Start piece one.
      R4=((R1-R2)*Y+R1+R2)/TWO
      B2=ONE-VF(R4)/EE
      IF(B2<ZERO)THEN
        PRINT *,'Error at position 4 in QINT and angle',FILENAME2
        PRINT *,'You need to change the number of potential values'
        PRINT *,'at large separations, increase the value for the'
        PRINT *,'minimum energy, or lower the accuracy requested.'
        STOP
      ENDIF
      B2=R4*SQRT(B2)
      TEMP=PI*(R1-R2)*R4*(ONE-VF(R4)/EE-R4/TWO*VD(R4)/EE)
      CALL EofB(B2,TEMP,FUN)

      !Start piece two.
      X3=(ONE+Y)/TWO
      R5=R1/SIN(HPI*X3)
      B2=ONE-VF(R5)/EE
      IF(B2<ZERO)THEN
        PRINT *,'Error at position 5 in QINT and angle',FILENAME2
        PRINT *,'You need to change the number of potential values'
        PRINT *,'at large separations, increase the value for the'
        PRINT *,'minimum energy, or lower the accuracy requested.'
        STOP
      ENDIF
      B2=R5*SQRT(B2)
      TEMP=(HPI*R1)**2*TWO*(R5/R1)**3*COS(HPI*X3)*(ONE-VF(R5)/EE-&
        &R5/TWO*VD(R5)/EE)
      CALL EofB(B2,TEMP,ANS)
      FUN=FUN+ANS
    ENDIF

  !Finally, examine the integrand if EE<ED or EE>EC2.
  ELSE
    IF(Y<=-ONE)THEN
      FUN=ZERO
    ELSEIF(Y>=ONE)THEN
      IF(EE<ED)THEN
        B2=ORBIT(2,1)
      ELSE
        B2=ORBIT(2,NO)
      ENDIF
      TEMP=TWO*PI*B2**2
      CALL EofB(B2,TEMP,FUN)
    ELSE
      !Start piece one.
      IF(EE<ED)THEN
        B2=ORBIT(2,1)/TWO*(Y+ONE)
        TEMP=HPI*ORBIT(2,1)**2*(Y+ONE)
      ELSE
        B2=ORBIT(2,NO)/TWO*(Y+ONE)
        TEMP=HPI*ORBIT(2,NO)**2*(Y+ONE)
      ENDIF
      CALL EofB(B2,TEMP,FUN)

      !Start piece two.
      IF(EE<ED)THEN
        B2=ORBIT(2,1)*TWO/(Y+ONE)
        TEMP=PI*ORBIT(2,1)**2*(TWO/(Y+ONE))**3
      ELSE
        B2=ORBIT(2,NO)*TWO/(Y+ONE)
        TEMP=PI*ORBIT(2,NO)**2*(TWO/(Y+ONE))**3
      ENDIF
      CALL EofB(B2,TEMP,ANS)
      FUN=FUN+ANS
    ENDIF
  ENDIF
END SUBROUTINE QINT

RECURSIVE SUBROUTINE INTEG(LMAX,F,A,B,EA,EB,ANS,ACC,IPR,PRANS)
  !This subroutine uses the Clenshaw-Curtiss procedure to simultaneously
  !  integrate LAMX functions given by external procedure F.
  !The lower and upper limits on the integral are A and B.  
  !In order to allow this procedure to be used when there are integrable
  !  singularities, the values of F at A must be supplied in EA and the
  !  values of F at B must be supplied in EB.
  !The relative accuracy of the calculated results in ANS is ACC.
  !If convergence is achieved with at least one value of ANS nonzero,
  !  the subroutine returns with IPR=0.  If convergence is not achieved
  !  but all of the values of ANS are equal to 0 or are very small, the
  !  subroutine returns with IPR=1.  Only if convergence is not
  !  achieved otherwise does it return with IPR=2.  PRANS is the
  !  previous estimate of ANS.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the dummy arguments explained above.
  INTEGER,INTENT(IN)::LMAX
  EXTERNAL F
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::A,B,ACC
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN),DIMENSION(LMAX)::EA,EB
  REAL(SELECTED_REAL_KIND(18)),INTENT(OUT),DIMENSION(LMAX)::ANS,PRANS
  INTEGER,INTENT(OUT)::IPR

  !Declare the local real arrays.
  !  FILE3 contains the function values used to calculate ANS.
  !  ERRMAX contains the maxima of the error estimates.
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(NPTSCH,LMAX)::FILE3
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(LMAX)::ERRMAX

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

  !LMAX must be positive.
  IF(LMAX<1)THEN
    PRINT *,'ERROR!  INTEG was called with LMAX=',LMAX
    PRINT *,'for the potential at angle ',&
      FILENAME2//FILENAME3//FILENAME4
    STOP
  ENDIF

  !Determine AMID and HALF
  AMID=(B+A)/TWO
  HALF=(B-A)/TWO

  !Initialize FILE3 and PRANS
  FILE3=ZERO
  PRANS=ZERO

  !For integration of order 1, place the values at the endpoints in
  !positions 1 and 3 of FILE3.  Then determine the value at the center
  !of the integration range (where Y=0 and the integration variable is
  !equal to AMID) and put it in position 2.
  FILE3(1,:)=EA
  FILE3(3,:)=EB
  CALL F(AMID,FILE3(2,:),LMAX)

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
        CALL F(AMID+HALF*CHEPTS(JSTART+I),FILE3(JSTART+I,:),LMAX)
      ENDDO
    ENDIF

    !Calculate ANS in order N.
    ANS=ZERO
    DO I=JSTART,JSTART+2**N
      ANS=ANS+CHEWTS(I)*FILE3(I,:)
    ENDDO
    ANS=ANS*HALF

    !A test for convergence is appropriate when N>1.
    IF(N>1)THEN
      !Set IPR=1 in the special case where all of the integrals are 0.
      IF(MAXVAL(ABS(ANS))<=1.0D-14)THEN
        IPR=1
        RETURN
      ENDIF

      !A conservative estimate of convergence requires comparing the
      !previous answers in PRANS with the present answers in ANS.  
      CONVERGE=.TRUE.
      ERRMAX=ZERO
      WHERE(ABS(ANS)>=1.0D-14)ERRMAX=ABS(ANS-PRANS)/ABS(ANS)
      IF(MAXVAL(ERRMAX)>ACC)CONVERGE=.FALSE.

      !If we have converged, set IPR=0 and return.
      IF(CONVERGE)THEN
        IPR=0
        RETURN
      ENDIF
    ENDIF

    !If we have reached the maximum order of quadrature without
    !converging in a relative sense, check for convergence to a very
    !small number.  If it is achieved, set IPR=1 and return.
    !If it is not achieved either, set IPR=2 and return.
    IF(N==NMAXCH)THEN
      ERRMAX=ABS(ANS-PRANS)
      IF(MAXVAL(ERRMAX)<ACC)THEN
        IPR=1
      ELSE
        IPR=2
      ENDIF
      RETURN
    ENDIF

    !Place ANS into PRANS in preparation for the next N value.
    PRANS=ANS

    !Copy the values in FILE3 to new positions appropriate for the
    !next higher approximation.
    DO I=JSTART,JSTART+2**N
      FILE3(2*I+3-N,:)=FILE3(I,:)
    ENDDO

    !Change JSTART so that it points to the first values in CHEPTS and
    !CHEWTS that apply to the next order of approximation.
    JSTART=JSTART+2**N+1
  ENDDO
END SUBROUTINE INTEG

SUBROUTINE MINIM(E,B,R)
  !This subroutine determines the root R of the equation Y supplied
  !as an external function.  It is important to note that an estimate
  !of R must be passed in initially.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the dummy variables.
  !  E is the energy.
  !  B is the impact parameter
  !  R is the separation.
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::E
  REAL(SELECTED_REAL_KIND(18)),INTENT(INOUT)::B,R

  !Declare the local real arrays.
  !  RC and YC are used to search for R by Lagrangian interpolations.
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(20)::RC,YC

  !Declare the local real variables.
  !  Y is the value of the external function whose root we want.
  !  VF is the potential energy given by an external function.
  !  VD is the potential derivative given by an external function.
  !  RSMAL is a benchmark for testing convergence.
  !  SCALE is a multiplicative factor used to find RA and RB.
  !  RA and RB bracket R.
  !  YA and YB are the values of external function Y at RA and RB.
  !  RT and YT are estimates better than RA and YA or RB and YB.
  !  SIGN is used to determine BING more accurately.
  REAL(SELECTED_REAL_KIND(18))::Y,VF,VD,RSMAL,SCALE,RA,YA,RB,YB,RT,&
    YT,SIGN

  !Declare the local integer variables.
  !  J is a loop index.
  !  N is the order of Lagrangian interpolation.
  INTEGER::J,N

  !We can solve analytically if NLONG=4 and R is large enough that
  !the long-range potential is the only thing to consider.
  IF(NLONG==4)THEN
    IF(E*B**4>=FOUR*CLONG)THEN
      RA=SQRT(B*B/TWO*(ONE+SQRT(ONE-FOUR*CLONG/E/B**4)))
      IF(RA>=POTEN(1,NV))THEN
        R=RA
        RETURN
      ENDIF
    ENDIF

  !We can solve numerically for any value of NLONG if CLONG<0.
  ELSEIF(CLONG<ZERO)THEN
    RA=(-CLONG/E)**(ONE/DBLE(NLONG))
    RB=MAX(B,(-CLONG/E+B*B*RA**(NLONG-2))**(ONE/DBLE(NLONG)))
    IF(RB>POTEN(1,NV))THEN
      DO J=1,20
        IF(ABS(ONE-RA/RB)<=1.0D-14)THEN
          R=RB
          RETURN
        ENDIF
        RA=RB
        IF(RA>POTEN(1,NV))EXIT
        RB=(-CLONG/E+B*B*RA**(NLONG-2))**(ONE/DBLE(NLONG))
      ENDDO
    ENDIF
    R=RB
  ENDIF

  !Initialize RSMAL and SCALE.  Use a large value of SCALE if E
  !is close to EC, in order to be very careful.
  RSMAL=MAX(1.0D-14,ACC1**6)
  SCALE=0.95D0
  IF(E>0.99*EC.AND.E<1.01*EC)THEN
    SCALE=0.999D0
  ENDIF
  
  !Be sure that the initial guess for R is small enough.
  R=MAX(R,POTEN(1,1))
  DO
    IF(Y(R,B,E)<ZERO)EXIT
    R=R*SCALE
    IF(R<0.5*POTEN(1,1))THEN
      PRINT *,'R is below ',POTEN(1,1),' in subroutine MINIM'
      PRINT *,'for the potential at angle ',&
        FILENAME2//FILENAME3//FILENAME4
      STOP
    ENDIF
  ENDDO

  !Be sure that the initial guess for R is large enough.
  DO
    IF(Y(R,B,E)<=ZERO)THEN
      R=R/SCALE
    ELSE
      EXIT
    ENDIF
    IF(R>1.0D12)THEN
      PRINT *,'R exceeded ',1.0D12,' in subroutine MINIM'
      PRINT *,'for the potential at angle ',&
        FILENAME2//FILENAME3//FILENAME4
      STOP
    ENDIF
  ENDDO
  !Bracket the turning point R between RA and RB.
  !Starting at POTEN(1,1) allows this subroutine to be used with
  !B=0, i.e. to determine the root of E=V(r).
  IF(B<POTEN(1,1))THEN
    DO J=1,NV
      IF(POTEN(2,J)<ZERO)THEN
        RA=POTEN(1,J)
        EXIT
      ENDIF
    ENDDO
  ELSE
    RA=R
    DO
      IF(Y(RA,B,E)>ZERO)EXIT
      RA=RA/SCALE
    ENDDO
  ENDIF

  !We have a value for RA but must still initialize RB.
  DO J=1,20
    YA=Y(RA,B,E)
    RB=SCALE*RA
    YB=Y(RB,B,E)
    IF(YA*YB<=ZERO)THEN
      EXIT
    ELSE
      R=MIN(RA,RB)
      RB=MAX(RA,RB)
      RA=R
    ENDIF
    YA=Y(RA,B,E)
    RB=RA/SQRT(SCALE)
    IF(RB>POTEN(1,NV))RB=POTEN(1,NV)
    YB=Y(RB,B,E)
    IF(YA*YB<=ZERO)THEN
      EXIT
    ELSE
      R=MIN(RA,RB)
      RB=MAX(RA,RB)
      RA=R
    ENDIF
    IF(J>=20)THEN
      RA=R
      YA=Y(RA,B,E)
      RB=POTEN(1,NV)
      YB=Y(RB,B,E)
    ENDIF 
    SCALE=0.95D0*SCALE
  ENDDO  
  !Make sure RA gives a value YA that is negative and that RB gives a
  !value YB that is positive.
  IF(YA>ZERO)THEN
    IF(YB>ZERO)THEN
      PRINT *,'There is a logic error in subroutine MINIM.'
      PRINT *,'for the potential at angle ',&
        FILENAME2//FILENAME3//FILENAME4
      STOP
    ENDIF
    RT=RA
    YT=YA
    RA=RB
    YA=YB
    RB=RT
    YB=YT
  ELSE
    IF(YB<ZERO)THEN
      PRINT *,'There is a logic error in subroutine MINIM.'
      PRINT *,'for the potential at angle ',&
        FILENAME2//FILENAME3//FILENAME4
      STOP
    ENDIF
  ENDIF

  !Narrow the range by a factor of 32 by looping 6 times.
  DO J=1,6
    RT=(RA+RB)/TWO
    YT=Y(RT,B,E)
    IF(YT>ZERO)THEN
      YB=YT
      RB=RT
    ELSEIF(YT<ZERO)THEN
      YA=YT
      RA=RT
    ELSE
      R=RT
      RETURN
    ENDIF
  ENDDO

  !Set up the arrays for Lagrangian interpolation.
  RC(1)=RA
  YC(1)=YA
  RC(2)=RB
  YC(2)=YB

  DO N=2,9
    !Calculate R accurately by interpolation.
    CALL AITKE(YC,RC,ZERO,R,N)
    RC(N+1)=R
    YC(N+1)=Y(RC(N+1),B,E)
    IF(ABS(YC(N+1))<RSMAL)EXIT
    IF(ABS(RC(N+1)-RC(N))<RSMAL)EXIT
    IF(ABS(RC(N+1)-RC(N-1))<RSMAL)EXIT
  ENDDO

  !A small error in R makes a large error in B.  Compensate by
  !adjusting BING so that any slight errors give Y(R,BING,E) positive. 
  DO
    BING=ONE-VF(R)/E
    IF(BING>ZERO)THEN
      BING=R*SQRT(BING)
      B=BING
      IF(Y(R,BING,E)>=ZERO)EXIT
    ELSE
      !The next line is crude, but the rest of this loop is devoted
      !to improving BING, so don't worry.
      BING=R*SQRT(ABS(BING))
    ENDIF
    SIGN=TWO*E*BING*BING/R/R/R-VD(R)
    IF(SIGN<ZERO.OR.SIGN>ZERO)THEN
      SIGN=SIGN/ABS(SIGN)
    ELSE
      RETURN
    ENDIF
    R=R+SIGN*RSMAL
    RSMAL=TWO*RSMAL
    IF(RSMAL>1.0D-12)EXIT
  ENDDO
END SUBROUTINE MINIM

SUBROUTINE EofB(B2,TEMP,ANS)
  !The subroutine completes the calculation of the integrands
  !QINT by evaluating the parts that depend upon the scattering
  !angle.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the dummy variables.
  !  B2 is the impact parameter.
  !  TEMP contains the constants that are multipliers.
  !  ANS contains the final answers.
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::B2,TEMP
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(MAXL),INTENT(OUT)::ANS

  !Declare the local variables.
  !  E is the local version of EE.
  !  B is the local version of B2.
  !  THETA is the scattering cross section.
  !  L is a loop index.
  REAL(SELECTED_REAL_KIND(18))::E,B,THETA
  INTEGER::L

  !Establish the values of E and B.
  E=EE
  B=B2

  !Determine the scattering angle.
  CALL ANGLE(E,B,THETA)

  !Determine the integrands.
  DO L=1,MAXL
    ANS(L)=TEMP*(ONE-COS(THETA)**L)
  ENDDO
END SUBROUTINE EofB

REAL(SELECTED_REAL_KIND(18)) FUNCTION Y(R,B,E)
  USE GLOBAL
  IMPLICIT NONE

  !Declare the dummy variable.
  !  R is the separation.
  !  B is the impact parameter
  !  E is the energy.
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::R,B,E

  !Declare the local real variables.
  !  VF is the potential energy given by an external function.
  REAL(SELECTED_REAL_KIND(18))::VF

  Y=E-VF(R)-E*B*B/R/R
END FUNCTION Y

SUBROUTINE AITKE(X,Y,XX,YY,N)
  !This subroutine interpolates in the lists X and Y to find the value
  !YY appropriate to the X value XX, using N-point Lagrangian
  !interpolation.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the dummy arguments explained above.
  INTEGER,INTENT(IN)::N
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(N),INTENT(IN)::X,Y
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::XX
  REAL(SELECTED_REAL_KIND(18)),INTENT(OUT)::YY

  !Declare the local real variables.
  !  TERM is temporary storage for a portion of YY.
  REAL(SELECTED_REAL_KIND(18))::TERM

  !Declare the local integer variables.
  !  I and J are loop indices.
  INTEGER::I,J

  IF(N<=1)THEN
    PRINT *,' Error!  AITKE called with N=',N,' at angle ',FILENAME2
    STOP
  ENDIF

  YY=ZERO
  DO I=1,N
    TERM=ONE
    DO J=1,N
      IF(J==I)CYCLE
      !Watch out for round-off errors that make X(I)=X(J).
      IF(X(I)<X(J).OR.X(I)>X(J))THEN
        TERM=TERM*(XX-X(J))/(X(I)-X(J))
      ELSE
        CYCLE
      ENDIF
    ENDDO
    YY=YY+TERM*Y(I)
  ENDDO 
END SUBROUTINE AITKE

SUBROUTINE ANGLE(E,B,THETA)
  !This subYoutine calculates the scattering angle, THETA, and the
  !phase shift, SHIFT, at energy E and impact parameter B.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the dummy variables.
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::E,B
  REAL(SELECTED_REAL_KIND(18)),INTENT(OUT)::THETA

  !Declare the local real arrays.
  !  EA and EB are the values of the integrand at the endpoints.
  !  ANS and BNS are the values of the integral.
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(1)::EA,EB,ANS,BNS

  !Declare the local real variables.
  !  FA and FB are the endpoints.
  !  R is the turning point.
  !  RBAR is the separation used to split the integral at low E or B.
  !  VD is the derivative of the potential energy.
  REAL(SELECTED_REAL_KIND(18))::FA,FB,R,RBAR,VD, ACOSARG

  !Declare the local integer variables.
  !  LMAX indicates the number of functions to be integrated.
  !  IPR is a flag for lack of convergence.
  !  I is a loop index.
  INTEGER::LMAX,IPR,I

  !Declare the local logical variable.
  !  CASE1 is true when Sec. 6.3.1 of the text is used.
  LOGICAL::CASE1

  !Declare the external functions.
  EXTERNAL GA,GB

  !Trap errors and print an expanation. 
  IF(B<=ZERO)THEN
    PRINT *,'Error at top of ANGLE, where B=',B
    PRINT *,'for the potential at angle ',&
      FILENAME2//FILENAME3//FILENAME4
    PRINT *,'You probably need to increase the minimum energy.'
    STOP
  ENDIF

  !Initialize variables that are common in both parts below.
  LMAX=1
  IPR=0
  EING=E
  BING=B
  FA=-ONE
  FB=ONE

  !Determine the value for CASE1.
  IF(E<ED.OR.E>EC.OR.NAITK<2)THEN
    CASE1=.TRUE.
  ELSE
    CASE1=.TRUE.
    DO I=1,NAITK-1
      IF(B>=BO(I).AND.B<BO(I+1))THEN
        CASE1=.FALSE.
        EXIT
      ENDIF
    ENDDO
  ENDIF

  !Calculate the scattering angle using Sec. 6.3.1 when appropriate.
  IF(CASE1)THEN
    R=ROMAX    
    CALL MINIM(EING,BING,R)

    !Trap errors and print an explanation.
    EA=ONE-R**3*VD(R)/TWO/B**2/E
    IF(EA(1)<=ZERO)THEN
      PRINT *,'Error at point 1 in ANGLE, with R=',R
      PRINT *,'for the potential at angle ',&
        FILENAME2//FILENAME3//FILENAME4
      PRINT *,'You need to change the potential at large '
      PRINT *,'separations, increase the value for the minimum'
      PRINT *,'energy, or lower the accuracy requested.'
      STOP
    ENDIF
    EA=ONE-ONE/SQRT(EA)
    EB=ONE-B/R

    !Pass R into subroutine GA.
    R1ING=R
    CALL INTEG(LMAX,GA,FA,FB,EA,EB,ANS,ACC1,IPR,BNS)

    !Watch out for lack of convergence for 1-cos(THETA), not THETA.
    IF(IPR==1)THEN
      ANS(1)=HPI*ANS(1)
      IF((ONE-ANS(1)**2/TWO)<-ONE.OR.(ONE-ANS(1)**2/TWO)>ONE)THEN
	PRINT *,'AARON ERROR: BAD ANGLE IN ACOS (line #2167)'
        PRINT *,'AARON ERROE: Tried to take ACOS of ', (ONE-ANS(1)**2/TWO)
        STOP
      ENDIF
      THETA=ACOS(ONE-ANS(1)**2/TWO)
      IF(ANS(1)<ZERO)THETA=-THETA

    !Otherwise, accept the value whether it is converged or not.  Errors
    !will be inconsequential unless they occur so many times that INTEG
    !cannot obtain converged values for the cross sections.
    ELSE
      THETA=HPI*ANS(1)
    ENDIF

  !Otherwise, calculate the scattering angle using Sec. 6.3.2.
  ELSE
    !Determine RBAR and R=Rm.
    RBAR=MAXVAL(RO)
    R=0.5*RBAR
    CALL MINIM(EING,BING,R)
    IF(R>=RBAR)THEN
      EA(1)=R
      R=RBAR
      RBAR=EA(1)
    ENDIF
    CALL MINIM(EING,BING,R)

    !Establish the integrand values at the endpoints.
    EA=ONE
    IF(B<=ZERO)THEN
      PRINT *,'Error at point 2 in ANGLE for the potential at angle ',&
        FILENAME2//FILENAME3//FILENAME4
      PRINT *,'You need to change the potential at large '
      PRINT *,'separations, increase the value for the minimum'
      PRINT *,'energy, or lower the accuracy requested.'
      STOP
    ELSE
      !In principle, EB in the next line should be positive.
      EB=ONE-R**3*VD(R)/TWO/E/B**2

      !If EB<0, call MINIM again to determine R more accurately.
      IF(EB(1)<ZERO)THEN
        CALL MINIM(EING,BING,R)
        EB=ONE-R**3*VD(R)/TWO/E/B**2
        !Give up if EB is still negative.
        IF(EB(1)<=ZERO)THEN
          PRINT *,'Error at point 2 in ANGLE for the potential at &
            &angle ',FILENAME2//FILENAME3//FILENAME4
          PRINT *,'You need to change the potential at large '
          PRINT *,'separations, increase the value for the minimum'
          PRINT *,'energy, or lower the accuracy requested.'
          STOP
        ENDIF
      ENDIF
      ACOSARG=R/RBAR
      IF(ACOSARG<-ONE.OR.ACOSARG>ONE)THEN
	PRINT *,'AARON ERROR: BAD ANGLE IN ACOS (line #2222)'
        PRINT *,'AARON ERROE: Tried to take ACOS of ', ACOSARG
	!STOP
	IF(ACOSARG<-ONE)THEN
	  ACOSARG=-ONE
	ELSE
	  ACOSARG=ONE
	ENDIF
      ENDIF
      EB=ONE-B/RBAR-ACOS(ACOSARG)/SQRT(EB)

    ENDIF

    !Pass R, RBAR, EA and EB into subroutine GB.
    R1ING=R
    R2ING=RBAR
    EAING=EA(1)
    EBING=EB(1)
    CALL INTEG(LMAX,GB,FA,FB,EA,EB,ANS,ACC1,IPR,BNS)

    !Watch out for lack of convergence for 1-cos(THETA), not THETA.
    IF(IPR==1)THEN
      ANS(1)=HPI*ANS(1)
      IF((ONE-ANS(1)**2/TWO)<-ONE.OR.(ONE-ANS(1)**2/TWO)>ONE)THEN
	PRINT *,'AARON ERROR: BAD ANGLE IN ACOS (line #2240)'
        PRINT *,'AARON ERROE: Tried to take ACOS of ', (ONE-ANS(1)**2/TWO)
        STOP
      ENDIF
      THETA=ACOS(ONE-ANS(1)**2/TWO)
      IF(ANS(1)<ZERO)THETA=-THETA
      RETURN
    !Otherwise, accept the value whether it is converged or not.  Errors
    !will be inconsequential unless they occur so many times that INTEG
    !cannot obtain converged values for the cross sections.
    ELSE
      THETA=HPI*ANS(1)
    ENDIF
  ENDIF
END SUBROUTINE ANGLE

SUBROUTINE GA(Y,FUNC,LMAX)
  !Subroutines GA and GB compute the integrands for calculating the
  !scattering angle.  They differ because of changes of variable.
  USE GLOBAL
  IMPLICIT NONE

  !Declare the dummy variables.
  !   Y is the integration variable.
  !   LMAX is the number of integrands to work with.
  !   FUNC is the computed value of the integral.
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::Y
  INTEGER,INTENT(IN)::LMAX
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(LMAX),INTENT(OUT)::FUNC

  !Declare the local variables.
  !  VF is the potential value returned by subprogram VF.
  !  VD is the potential derivative returned by subprogram VD.
  !  R is the separation at which VF is required.  
  REAL(SELECTED_REAL_KIND(18))::VF,R

  !Calculate R=R6 using R1ING [called Rm in the manuscript].
  R=R1ING/COS(PI*(Y+ONE)/FOUR)
  !In principle, FUNC(1) in the next line should be positive.
  !If not, evaluate the next line by a Taylor series.
  FUNC(1)=ONE-(BING/R)**2-VF(R)/EING
  IF(FUNC(1)<=ZERO)THEN
    FUNC(1)=(R-R1ING)*(TWO*BING**2/R1ING**3-VF(R1ING)/EING)
  ENDIF

  !If it is still not positive, give up.
  IF(FUNC(1)<=ZERO)THEN
    FUNC(1)=ZERO
  ELSE
    FUNC(1)=ONE-BING/R1ING*SIN(PI*(Y+ONE)/FOUR)/SQRT(FUNC(1))
  ENDIF
END SUBROUTINE GA

SUBROUTINE GB(Y,FUNC,LMAX)
  !Subroutines GA and GB compute the integrands for calculating the
  !scattering angle.  They differ because of changes of vari:e PC
  USE GLOBAL
  IMPLICIT NONE

  !Declare the dummy variables.
  !  Y is the integration variable.
  !  LMAX=1 is the number of integrands to work with.
  !  FUNC is the computed value of the integral.
  REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::Y
  INTEGER,INTENT(IN)::LMAX
  REAL(SELECTED_REAL_KIND(18)),DIMENSION(LMAX),INTENT(OUT)::FUNC

  !Declare the local variables.
  !  VF is the potential value returned by subprogram VF.
  !  R is the separation at which VF is required.
  !  Z1, Z2 and Z3 are used for temporary storage.
  REAL(SELECTED_REAL_KIND(18))::VF,R,Z1,Z2,Z3,ACOSARG

  !Calculate R7 using R2ING (r bar).
  R=R2ING/COS(PI*(Y+ONE)/FOUR)

  !Watch out for roundoff errors.
  Z3=ONE-(BING/R)**2-VF(R)/EING
  IF(Z3>ZERO.OR.Z3<ZERO)THEN
    Z3=ABS(Z3)
    FUNC(1)=ONE-BING/R2ING*SIN(PI*(Y+ONE)/FOUR)/SQRT(Z3)
  ELSE
    PRINT *,'Zero value at point 1 in subroutine GB for the potential &
      &at angle',FILENAME2//FILENAME3//FILENAME4
    STOP
  ENDIF

  !Calculate R8 using R1ING (called Rm) and R2ING.
  !Again, watch out for roundoff errors.
  ACOSARG=R1ING/R2ING
  IF(ACOSARG<-ONE.OR.ACOSARG>ONE)THEN
    PRINT *,'AARON ERROR: BAD ANGLE IN ACOS (line #2330)'
    PRINT *,'AARON ERROE: Tried to take ACOS of ', ACOSARG
    !STOP
    IF(ACOSARG<-ONE)THEN
      ACOSARG=-ONE
    ELSE
      ACOSARG=ONE
    ENDIF
  ENDIF
  Z1=ACOS(ACOSARG)
  Z2=Z1*COS(PI*(Y+ONE)/FOUR)
  R=R1ING/COS(Z2)

  !In principle, Z3 in the next line should be positive.
  Z3=ONE-(BING/R)**2-VF(R)/EING
  IF(Z3>ZERO)THEN
    FUNC(1)=FUNC(1)-BING/R1ING*Z1*SIN(Z2)*SIN(PI*(ONE+Y)/FOUR)/SQRT(Z3)
  !If not, set FUNC(1) equal to EA or EB.
  ELSE
    IF(ABS(Y)<ZERO)THEN
      FUNC(1)=EAING
    ELSE
      FUNC(1)=EBING
    ENDIF
  ENDIF
END SUBROUTINE GB

SUBROUTINE GEN
  USE GLOBAL
  IMPLICIT NONE

  !Declare the real arrays.
  !  The energy and first 30 transport cross sections are in DATA.
  REAL(SELECTED_REAL_KIND(18))::DATA(31)

  !Declare the integer variables.
  !  COL indicates a particular column in a line of data.
  INTEGER::COL

  !Extract the data and put it into the output file.
  REWIND (15)
  READ(15,*)
  DO
    READ(15,*,IOSTAT=ERR)(DATA(COL),COL=1,31)
    IF(ERR/=0)EXIT
    WRITE(12,'(2(1X,F30.14))')DATA(1),DATA(2)*SQRT(DATA(1))
  ENDDO
END SUBROUTINE GEN
