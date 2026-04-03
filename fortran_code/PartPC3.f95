MODULE GLOBAL
 !This is the data module for program PC3.
 IMPLICIT NONE

 !Declare the global named constants.
 !  INOOMAX is the maximum number of regions of multiple orbiting.
 !  MAXL is the number of transport cross sections to be calculated.
 !  MAXORBITS is the maximum number of sets of orbiting parameters.
 !  MAXE is the maximum number of energies in each of the 3 regions.
 !    It should be one larger than 2 to some power.
 !  ZERO, ONE, TWO, THREE, FOUR and TEN are high-precision numbers.
 INTEGER,PARAMETER::INOOMAX=40,MAXL=30,MAXORBITS=600,MAXE=65
 REAL(SELECTED_REAL_KIND(18)),PARAMETER::ZERO=0,ONE=1,TWO=2,THREE=3,&
   FOUR=4,TEN=10

 !Declare the global integers.
 !  ERROR allows for recovery from input/output errors.
 !  NMAXCH is the maximum order of approximation in using the Curtiss-
 !    Clenshaw method of integration.  It is important that the value
 !    of NMAXCH here be the same as in CCgen.in.
 !  NPTSCH is the size of the arrays CHEPTS and CHEWTS.
 !  NV is the number of equally-spaced separations.
 !  NANG is the number of equally-spaced angles.
 !  SYM is 90 if the diatomic ion is homonuclear, but 180 otherwise.
 !  NLONG is the exponent of the long-range potential.
 INTEGER::ERROR,NMAXCH,NPTSCH,NV,NANG,SYM,NLONG

 !Declare the global integer arrays.
 !  NAITK is the number of orbiting values at a particular E, but it is
 !    reset to indicate the number of valid entries in BO and RO.
 !  INOO is the number of orbiting intervals at any energy.
 !  NO is the number of valid sets of orbiting parameters.
 !  NOO stores the starting and ending indices for orbiting parameters.
 INTEGER,DIMENSION(:),ALLOCATABLE::NAITK,INOO,NO
 INTEGER,DIMENSION(2,INOOMAX)::NOO=0

 !Declare the global real variables.
 !  HPI is half of the usual pi for circular functions.
 !  PI is the usual pi for circular functions.
 !  EMIN and EMAX are the minimum and maximum energies (in hartree).
 !  ACCURACY is the desired accuracy of the transport cross sections.
 !  ACC1 is used to compute the integrands with greater accuracy
 !     than is required (via ACCURACY) for the integrals.
 !  CLONG is the coefficient of the long-range inverse power at each angle.
 !  ED and EC are the minimum and maximum orbiting energies.
 !  EC2 is always equal to ten times global variable EC.
 !  EE is an energy of interest.
 !  ROMAX and BOMAX are the maximum values of RO and BO.
 !  RM1 is used to transfer orbiting separations between subroutines.
 !  THO is the scattering angle transferred to functions GA and GB.
 !  EING and BING are used to transfer E and B to functions GA and GB.
 !  R1ING and R2ING are used to similarly transfer separations.
 !  EAING and EBING are used to similarly transfer values at endpoints.
 REAL(SELECTED_REAL_KIND(18))::HPI,PI,EMIN,EMAX,ACCURACY,ACC1,CLONG
   !ED,EC,EC2,EE,ROMAX,BOMAX,RM1,THO,EING,BING,R1ING,R2ING,EAING,EBING

 !Declare the global real arrays.
 !  The quadrature points and weights are placed in CHEPTS and CHEWTS.
 !  POTENR contains the separations used for the potential values.
 !  EXPON is the power of the 1/R dependence of the short range potential.
 !  CSHORT is the coefficient of the short-range inverse power.
 !  POTENV contains the potential values.
 !  CPOTEN contains the parameters for a spline fit of the potential.
 !  ECC stores the energy ranges for orbiting.
 !  ORBIT stores the sets of orbiting parameters.
 !  B0 and RO contain the orbiting impact parameters and orbiting
 !    separations corresponding to orbiting energy EE.
 REAL(SELECTED_REAL_KIND(18)),DIMENSION(:),ALLOCATABLE::CHEPTS,CHEWTS,&
   POTENR,EXPON,CSHORT
 REAL(SELECTED_REAL_KIND(18)),DIMENSION(:,:),ALLOCATABLE::POTENV,CPOTEN
 REAL(SELECTED_REAL_KIND(18)),DIMENSION(:,:,:),ALLOCATABLE::ECC,ORBIT
 REAL(SELECTED_REAL_KIND(18)),DIMENSION(INOOMAX)::BO,RO

 !Declare the global character variables.
 !  CCFILE is the fully qualified name of an input file.
 CHARACTER(LEN=78)::CCFILE
END MODULE GLOBAL

PROGRAM PC3
 !This program was written in 2024, by Larry A. Viehland, Department of
 !Science, Chatham University, Pittsburgh, Pennsylvania 15232 USA.  It
 !computes the transport cross sections for diatomic ions interacting with
 !atomic neutrals, using the techniques of L. A. Viehland et al. (2025).
 USE GLOBAL
 IMPLICIT NONE

 !Declare the local integer variables.
 !  I and J are loop indices.
 INTEGER::I,J

 !Open the input file.
 OPEN(UNIT=10,FILE="CCgen.in",STATUS="OLD",IOSTAT=ERROR)
 IF(ERROR/=0)THEN
   PRINT *,"Could not open the old file named CCgen.in."
   STOP
 ENDIF

 !Read NMAXCH and CCFILE.
 READ(10,*)NMAXCH
 READ(10,"(A)")CCFILE
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
 OPEN(10,FILE=CCFILE,IOSTAT=ERROR,STATUS="OLD")
 IF(ERROR/=0)THEN
   PRINT *,"Could not open the old file named ",TRIM(CCFILE),"."
   PRINT *,"Be sure that you have permission to read this file."
   STOP
 ENDIF
 DO I=1,NPTSCH
   READ(10,"(18X,2F23.20)",IOSTAT=ERROR)CHEPTS(I),CHEWTS(I)
   IF(ERROR/=0)THEN
     PRINT *,"File ",TRIM(CCFILE)," is incompatible with this program."
     PRINT *,"Be sure that it was created with program CCgen.f95."
     STOP
   ENDIF
 ENDDO

 !If everything is consistent, the next READ statement should fail
 !since EMIN and EMAX are dummy variables at this point.
 READ(10,"(18X,2F23.20)",IOSTAT=ERROR)EMIN,EMAX
 IF(ERROR==0)THEN
   PRINT *,"File ",TRIM(CCFILE)," is incompatible with this program."
   PRINT *,"Be sure that it was created with program CCgen.f95."
   STOP
 ENDIF
 CLOSE(10)
 
 !READ the potential information.
 CALL INPUT

 !Spline fit the potentials.
 !CALL SPLINEFIT

 !Determine the orbiting parameters and the orbiting regions.
 ALLOCATE(NO(NANG),ORBIT(3,MAXORBITS,NANG))
 DO J=1,NANG
   CALL ORBITS(J)  
   CALL REGIONS(J)
 ENDDO

 !Calculate and print the cross sections.
 !CALL EFIT
END PROGRAM PC3

SUBROUTINE INPUT
 !This subroutine reads information about the interaction potentials.
 USE GLOBAL
 IMPLICIT NONE

 !Declare the local character variables.
 !  COMMENT has lines of commentary describing the system.
 CHARACTER(LEN=72)::COMMENT(4)

 !Declare the local integer variables.
 !  I and J are loop indices.
 INTEGER::I,J

 !Open the main input and output files.
 OPEN(UNIT=12,FILE="PC3.in",STATUS="OLD",IOSTAT=ERROR)
 IF(ERROR/=0)THEN
   PRINT *,"Could not open the old file named PC3.in."
   STOP
 ENDIF
 OPEN(UNIT=13,FILE="PC3.out",IOSTAT=ERROR)
 IF(ERROR/=0)THEN
   PRINT *,"Could not open an output file named PC3.out."
   STOP
 ENDIF
 WRITE(13,"(A)")"Transport cross sections for the"

 !Read in and write out the four comment lines.
 DO I=1,4
   READ(12,"(A)",IOSTAT=ERROR)COMMENT(I)
   IF(ERROR/=0)THEN
     PRINT *,"Something is wrong with line ",I," of file PC3.in."
     STOP
   ENDIF
   WRITE(13,"(A)")COMMENT(I)
 ENDDO

 !Read in ACCURACY, EMIN, and EMAX.
 READ(12,*,IOSTAT=ERROR)ACCURACY,EMIN,EMAX
 IF(ERROR/=0)THEN
   PRINT *,"Something is wrong with line 5 of file PC3.in."
   STOP
 ENDIF

 !Initialize global variable ACC1, requiring the integrals to be
 !calculated a bit more accurately than the overall results.
 ACC1=0.8D0*ACCURACY

 !Read in NV, NANG and SYM.
 READ(12,*,IOSTAT=ERROR)NV,NANG
 IF(ERROR/=0)THEN
   PRINT *,"Something is wrong with line 6 of file PC3.in."
   STOP
 ENDIF
 READ(12,*,IOSTAT=ERROR)SYM
 IF(ERROR/=0.OR.(SYM/=90.AND.SYM/=180))THEN
   PRINT *,"Something is wrong with line 7 of file PC3.in."
   STOP
 ENDIF

 !Check NV for a value that is too low.
   IF(NV<=3)THEN
     PRINT *,"The NV value on line 6 of file PC3.in must be 4 or more."
     STOP
   ENDIF
 
 !Allocate the arrays
 !TO DO LARRY: Rank mismatch NOO below
 ALLOCATE(POTENR(NV),EXPON(NANG),CSHORT(NANG),POTENV(NV,NANG),NAITK(NANG),&
   CPOTEN(NV,NANG),NOO(NV,2,INOOMAX),INOO(NANG),ECC(NV,2,INOOMAX))

 !Set the value of NLONG for ion-neutral potentials.
 NLONG=4

 !Read the potentials. '
 DO J=1,NANG
   READ(12,*,IOSTAT=ERROR)
   IF(ERROR/=0)THEN
     PRINT *,"There is an error on line ",J+7," of file PC3.in"
     STOP
   ENDIF
   OPEN(11,FILE=TRIM(CCFILE),IOSTAT=ERROR)
   IF(ERROR/=0)THEN
     PRINT *,"Could not open file ",TRIM(CCFILE)
     PRINT *,"Check line ",J+7," of file PC3.in"
     STOP
   ENDIF
   DO I=1,NV
     READ(11,*,IOSTAT=ERROR)POTENR(I),POTENV(I,J)
     IF(ERROR/=0)THEN
       PRINT *,' There is an error on line ',I,' of file PC3.in'
       STOP
     ENDIF
   ENDDO
   CLOSE(11)
 ENDDO
 CLOSE(12)

 !Make sure that EMAX is not greater than any first potential value.
 EMAX=MIN(EMAX,MINVAL(POTENV(1,:)))

 !Make sure that the POTENR values increase monotonically.
 DO I=2,NV
   IF(POTENR(I)<=POTENR(I-1))THEN
     PRINT *,"The separations do not increase monotonically."
     PRINT *,"Check the poten files at or near separation ",POTENR(I)
     STOP
   ENDIF
 ENDDO

 !Make sure that the first 3 POTENV values are positive.
 DO I=1,3
   DO J=1,NANG
     IF(POTENV(I,J)<ZERO)THEN
       PRINT *,"Tabulated potential value ",POTENV(I,J)," is negative."
       PRINT *,"This is not allowed, so add more points at small R"
       PRINT *,"for angle ",J," before running this program again."
       STOP
     ENDIF
   ENDDO
 ENDDO

 !Make sure that the last 3 POTENV values are negative.
 DO I=NV-2,NV
   DO J=1,NANG
     IF(POTENV(I,J)>ZERO)THEN
       PRINT *,"Tabulated potential value ",POTENV(I,J)," is positive."
       PRINT *,"This is not allowed, so add more points at large R"
       PRINT *,"for angle ",J," before running this program again."
       STOP
     ENDIF
   ENDDO
 ENDDO

 !Determine EXPON, CSHORT and CLONG.  Note the asymmetry in the way 
 !that the short- and long-range potentials are evaluated in VF.
 !Note also that CLONG is the same for all potentials.
 DO J=1,NANG
   EXPON(J)=LOG(POTENV(2,J)/POTENV(1,J))/LOG(POTENR(1)/POTENR(2))
   CSHORT(J)=POTENV(1,J)*POTENR(1)**EXPON(J)
 END DO
 CLONG=-POTENV(NV,NANG)*POTENR(NV)**NLONG

 !Write information to the output file.
 WRITE(13,*)"At every angle, the long-range potential is assumed to be:"
 IF(CLONG>ZERO)THEN
   WRITE(13,"(A,G15.8,A,I1,A)")"V(R)=-",CLONG,"*R**(-",NLONG,")"
 ELSE
   WRITE(13,"(A,G15.8,A,I1,A)")"V(R)=",-CLONG,"*R**(-",NLONG,")"
 ENDIF
 WRITE(13,"(A,I4,A)")"and there are ",NV," pairs of potential values."
 WRITE(13,"(A)")"Listed below are the input values of R and V."
 DO J=1,NANG
   WRITE(13,*) DBLE(J-1)*DBLE(SYM)/DBLE(NANG-1)
   DO I=1,NV
     WRITE(13,"(2(1X,G22.14))")POTENR(I),POTENV(I,J)
   ENDDO
 ENDDO
 WRITE(13,*)
 WRITE(13,"(A,G15.8,A)")"Minimum energy is EMIN=",EMIN," hartree."
 WRITE(13,"(A,G15.8,A)")"Maximum energy is EMAX=",EMAX," hartree."
 WRITE(13,"(A,G12.6,A)")"The relative accuracy is ",ACCURACY,"."
END SUBROUTINE INPUT

SUBROUTINE SPLINEFIT
 !This subroutine guides the interpolations for the spline of V(R,theta).
 USE GLOBAL
 IMPLICIT NONE

 !Declare the local integer variable.
 !  J is a loop index.
 INTEGER::J

 !Cycle through the angles.
 DO J=1,NANG
   CALL SPLINE(NV,POTENR,POTENV(1:NV,J),-POTENV(1,J)*EXPON(J)/POTENR(1),&
     -POTENV(NV,J)*NLONG/POTENR(NV),CPOTEN(1:NV,J))
 ENDDO
END SUBROUTINE SPLINEFIT

SUBROUTINE SPLINE(N,X,Y,YP1,YPN,Y2)
 !This subroutine is based on Sec. 3.3 of "Numerical Recipes" (FORTRAN
 !VERSION) by W. H. Press, B. P. Flannery, S. A. Teukolsky and W. T.
 !Vetterling (Cambridge University Press, Cambridge, England, 1989).
 !It computes the second derivative matrix Y2 that is needed for a
 !spline fit, using function FIT, of the N points in arrays X and Y,
 !with Y as a function of X.  Since the fit will generally be a
 !clamped spline, the derivatives at the endpoints must be supplied in
 !YP1 and YPN; if the values in YP1 and YPN exceed 0.99D30, then a
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

SUBROUTINE ORBITS(JANGLE)
 !This subroutine computes the orbiting parameters for angle JANGLE.
 USE GLOBAL
 IMPLICIT NONE

 !Declare the dummy arguments and local real variables.
 !  JANGLE is the fixed angle.
 !  THETA is the angle and R is the separation of interest.
 !  VF is the value of the potential at R.
 !  VD is the first derivative of VF(R) at R.
 !  VD2 is the second derivative of VF(R) at R.
 !  C1-C3 are used to check the orbiting conditions.
 !  VEFF is the effective potential at orbiting.
 INTEGER,INTENT(IN)::JANGLE
 REAL(SELECTED_REAL_KIND(18))::THETA
 REAL(SELECTED_REAL_KIND(18))::R,VF,VD,VD2,C1,C2,C3,VEFF

 !Declare the local integer variables.
 !  I is a loop index.
 !  J labels a particular set of orbiting parameters.
 !  K is used to retain the value of J while sets are added.
 INTEGER::I,J,K

 !Initialize THETA and global variables NO, ORBIT and J.  Note that NO is
 !designed to work down at first; later it will work up to its final value.
 THETA=DBLE(JANGLE-1)/DBLE(NANG-1)*DBLE(SYM)
 NO(JANGLE)=MAXORBITS+1
 ORBIT(:,:,JANGLE)=ZERO
 J=0

 !There are three criteria for orbiting at fixed THETA.  The first is that
 !VD(ro)>0. The second is that VD2(ro)+3*VD(ro)/ro < 0.  The third is that 
 !Veff>0.  If these are satisfied, then Eo=V(ro)+ro*VD(ro)/2 and 
 !bo=SQRT(ro**3*VD(ro)/2/Eo).  Start the loop looking for orbiting 
 !conditions at the first tabulated separation on the potential. 
 R=POTENR(1) 
 DO
   !At this point we have a designated value of r0, called R.  It was
   !either set before the loop or was reset at the end of the loop.
   !C1=VD(R,THETA)
   !C2=VD2(R,THETA)+THREE*C1/R
   C3=VF(R,THETA)
   VEFF=C3+R*C1/TWO
   STOP

   !Check the orbiting conditions.
   IF(C1>ZERO.AND.C2<ZERO.AND.VEFF>ZERO)THEN
     !There is at least one orbiting region.  We"ll look for more
     !later, after we store the set we have found.
     NO(JANGLE)=NO(JANGLE)-1
     IF(NO(JANGLE)<=0)THEN
       PRINT *," More room is needed for orbiting parameters!"
       STOP
     ENDIF
     ORBIT(1,NO(JANGLE),JANGLE)=VEFF
     ORBIT(2,NO(JANGLE),JANGLE)=SQRT(R**3*C1/VEFF/TWO)
     ORBIT(3,NO(JANGLE),JANGLE)=R

     !Improve the guess by many changes in R.  When NO(JANGLE)=MAXORBITS,
     !we must find the maximum orbiting energy to a very high degree of
     !accuracy.  For other values of NO(JANGLE), we use a reasonably high
     !degree of accuracy, checking to see if this is an energy where one
     !region of multiple orbiting ends.
     DO I=1,7
       DO
         R=(ONE-0.001**I)*ORBIT(3,NO(JANGLE),JANGLE)
         C1=VD(R,THETA)
         C2=VD2(R,THETA)+THREE*C1/R
         C3=VF(R,THETA)
         VEFF=C3+R*C1/TWO
         IF(C1>ZERO.AND.C2<ZERO.AND.VEFF>ORBIT(1,NO(JANGLE),JANGLE))THEN
           ORBIT(1,NO(JANGLE),JANGLE)=VEFF
           ORBIT(2,NO(JANGLE),JANGLE)=SQRT(R**3*C1/VEFF/TWO)
           ORBIT(3,NO(JANGLE),JANGLE)=R
           IF(NO(JANGLE)==MAXORBITS)CYCLE
         ENDIF
         EXIT
       ENDDO
     ENDDO
   ENDIF

   !If there is at least one orbiting region, we exit the DO loop when the
   !orbiting energy is less than EMIN but greater than 0.  The latter
   !condition arises because we initially set ORBIT=0, so a value of 0
   !indicates that no orbiting was found.
   IF(NO(JANGLE)<=MAXORBITS)THEN
     IF(ORBIT(1,NO(JANGLE),JANGLE)<=EMIN.AND.ORBIT(1,NO(JANGLE),JANGLE)&
       >ZERO)EXIT
   ENDIF

   !Move to a larger separation by making rather big steps (5%) beyond 
   !the last tabulated separation but smaller steps (1%) otherwise. 
   IF(R>POTENR(NV))THEN
     R=1.05*R
   ELSE
     R=1.01*R
   ENDIF

   IF(R>1.0D3)THEN
     IF(NO(JANGLE)>MAXORBITS)THEN
       PRINT *,"Even at 1000 bohr, the potential did not support orbiting."
       STOP
     ENDIF
     EXIT
   ENDIF
 ENDDO

 !We now make a pass through the orbiting parameters in order to add
 !sets in regions where the energy is decreasing as the separation
 !decreases.
 DO
   J=J+1
   !Move a set down from position NO to position J in ORBIT.
   ORBIT(:,J,JANGLE)=ORBIT(:,NO(JANGLE),JANGLE)

   !Reset NO so that it will be ready for the next pass.
   NO(JANGLE)=NO(JANGLE)+1
   IF(NO(JANGLE)>MAXORBITS)EXIT

   !If J ever becomes larger than NO-19, we are in trouble since there 
   !won"t be enough room to add more sets.
   IF(J>NO(JANGLE)-19)THEN
     PRINT "(A)","There is not enough room for the orbiting parameters."
     STOP
   ENDIF

   !If the energies decrease, add up to 19 sets.
   IF(NO(JANGLE)<MAXORBITS)THEN
     IF(ORBIT(1,J,JANGLE)>ORBIT(1,NO(JANGLE),JANGLE))THEN
       !Use K in order to retain the value of J before the additional
       !orbiting sets are added in a loop.  This means that we will
       !test separations that are linearly spaced between K and NO.
       K=J
       DO I=1,19
         R=(DBLE(20-I)*ORBIT(3,K,JANGLE)+DBLE(I)*&
           ORBIT(3,NO(JANGLE),JANGLE))/TEN/TWO
         C1=VD(R,THETA)
         C2=VD2(R,THETA)+THREE*C1/R
         C3=VF(R,THETA)
         VEFF=C3+R*C1/TWO
         IF(C1>ZERO.AND.C2<ZERO.AND.VEFF>ORBIT(1,J,JANGLE))THEN
           J=J+1
           ORBIT(1,J,JANGLE)=VEFF
           ORBIT(2,J,JANGLE)=SQRT(R**3*C1/VEFF/TWO)
           ORBIT(3,J,JANGLE)=R
         ENDIF
       ENDDO
     ENDIF
   ELSE
     IF(ORBIT(1,J,JANGLE)<ORBIT(1,J-1,JANGLE))THEN
       DO I=1,19
         R=ORBIT(3,J,JANGLE)+DBLE(I)/TEN/TWO*(ORBIT(3,J,JANGLE)-&
           ORBIT(3,J-1,JANGLE))
         C1=VD(R,THETA)
         C2=VD2(R,THETA)+THREE*C1/R
         C3=VF(R,THETA)
         VEFF=C3+R*C1/TWO
         IF(C1>ZERO.AND.C2<ZERO.AND.VEFF>ORBIT(1,J,JANGLE))THEN
           J=J+1
           ORBIT(1,J,JANGLE)=VEFF
           ORBIT(2,J,JANGLE)=SQRT(R**3*C1/VEFF/TWO)
           ORBIT(3,J,JANGLE)=R
         ENDIF
       ENDDO
     ENDIF
   ENDIF
 ENDDO 

 !Chop off any points where the energies have yet to start going up.
 DO
   IF(ORBIT(1,J,JANGLE)<ORBIT(1,J-1,JANGLE))THEN
     J=J-1
   ELSE
     EXIT
   ENDIF
 ENDDO

 !There are J sets of orbiting parameters, so put J into NO(JANGLE) and
 !write the parameters to the output file.
 NO(JANGLE)=J
 WRITE(13,"(/,A,I3,A)")" The ",NO(JANGLE)," sets of orbiting parameters are:"
 WRITE(13,"(3A)")"          E                     B&
   &                     RM"
 WRITE(13,"(3(1X,G22.14E4))")((ORBIT(J,I,JANGLE),J=1,3),I=1,NO(JANGLE))
END SUBROUTINE ORBITS

SUBROUTINE REGIONS(JANGLE)
 !This subroutine determine the orbiting regions, the collections of
 !orbiting parameters in which the orbiting energy always increases.
 USE GLOBAL
 IMPLICIT NONE

 !Declare the dummy arguments and local real variables.
 !  JANGLE labels the angle of interest.
 !  J is a loop index. 
 INTEGER,INTENT(IN)::JANGLE
 INTEGER::J

 !Declare the local logical variable.
 !  INCREASING is true if energies are rising.
 LOGICAL::INCREASING

 INOO(JANGLE)=1
 !Set the first energy in the first region, and the pointer that
 !indicates the corresponding number of the orbiting set.
 ECC(JANGLE,1,1)=ORBIT(1,JANGLE,1)
 !TO DO LARRY: Fix rank mismatch for NOO error below and everywhere else!
 NOO(JANGLE,1,1)=1

 !Initially, the energies must be increasing.
 INCREASING=.TRUE.

 !Use a loop to look for the end of one region and the start of the next.
 DO J=2,NO(JANGLE)
   !Treat the situation where the energies quit increasing.
   IF(INCREASING.AND.ORBIT(1,J,JANGLE)<ORBIT(1,J-1,JANGLE))THEN
     !Place the ending energy into ECC and the number of the
     !corresponding orbiting set into NOO.
     ECC(JANGLE,2,INOO(JANGLE))=ORBIT(1,J-1,JANGLE)
     NOO(JANGLE,2,INOO(JANGLE))=J-1

     !Write information about the region.
     WRITE(13,*)ECC(JANGLE,1,INOO(JANGLE)),"-",ECC(JANGLE,2,INOO(JANGLE)),&
       " Region ",INOO(JANGLE)

     !Get ready for the next region.
     INOO(JANGLE)=INOO(JANGLE)+1
     IF(INOO(JANGLE)>INOOMAX)THEN
       PRINT *,"INOOMAX is too small.  You must reset this named"
       PRINT *,"constant, recompile the program and try again."
       STOP
     ENDIF
     INCREASING=.FALSE.

   !Treat the situation where the energies resume increasing.
   ELSEIF(.NOT.INCREASING.AND.ORBIT(1,J,JANGLE)>ORBIT(1,J-1,JANGLE))THEN
     ECC(JANGLE,1,INOO(JANGLE))=ORBIT(1,J-1,JANGLE)
     ! LARRY TO DO: The array access needs three indexes in next line. Fix this.
     !NOO(1,INOO(JANGLE))=J-1
     INCREASING=.TRUE.
   ENDIF
 ENDDO

 !Place the ending energy into ECC and the number of the
 !corresponding orbiting set into NOO.
 !TO DO LARRY: Fix incompatible ranks in assignment below
 ECC(JANGLE,2,INOO(JANGLE))=ORBIT(1,NO,JANGLE)
 ! LARRY TO DO: The array access needs three indexes in next line. Fix this.
 ! NOO(2,INOO(JANGLE))=NO(JANGLE)

 !Write information about the final region.
 WRITE(13,*)ECC(JANGLE,1,INOO(JANGLE)),"-",ECC(JANGLE,2,INOO(JANGLE)),&
   " Region ",INOO(JANGLE)
 WRITE(13,*)
END SUBROUTINE REGIONS

REAL(SELECTED_REAL_KIND(18)) FUNCTION VF(R,PHI)
 !This function evaluates the interaction potential, which is assumed
 !to be repulsive at short separations.  It is also assumed to smoothly
 !connect to the proper long-range (asymptotic) tail.
 USE GLOBAL
 IMPLICIT NONE

 !Declare the dummy variables.
 !  R is the ion-neutral separation in bohr.
 !  PHI is the potential angle, from 0 to SYM degrees
 REAL(SELECTED_REAL_KIND(18)),INTENT(IN)::R,PHI

 !Declare the local real variable.
 !  FIT is an external function for the spline fit.
 REAL(SELECTED_REAL_KIND(18))::FIT

 !Declare the local integer variable.
 !  JANGLE labels which particular angle being considered.
 INTEGER::JANGLE

 !Determine JANGLE
 ! TO DO: LARRY Check the formula below, should contain #Angles
 JANGLE=NINT(ABS(PHI)/SYM+1)
 IF(JANGLE>180.0)JANGLE=JANGLE-180.0D0

 !First, a possible short-range extrapolation must be considered.
 ! TO DO, LARRY: Consider if we should average the extrapolations for the
 ! two closest angles. (Same as long range extrapolation below.)
 IF(R<POTENR(1))THEN
   VF=CSHORT(JANGLE)*(POTENR(1)/R)**EXPON(JANGLE)
   RETURN

 !Second, we consider a possible long-range extrapolation.
 ! TO DO LARRY: Should this be CLONG(JANG) instead?
 ! Also, maybe should average adjacent angle values.
 ELSEIF(R>POTENR(NV))THEN
   VF=-CLONG/R**NLONG
   RETURN

 !Finally, get the potential from the spline fit of the interpolation
 !function of the distances and a linear fit of the angles.
 ELSE
   VF=(FIT(NV,POTENR(1:NV),POTENV(1:NV,JANGLE+1),CPOTEN(1:NV,JANGLE+1),R) &
      - FIT(NV,POTENR(1:NV),POTENV(1:NV,JANGLE),CPOTEN(1:NV,JANGLE),R)) &
      /ABS(PHI)*180.0D0/DBLE(NANG-1)
   ! TO DO LARRY: Check the formula above.
 ENDIF
 PRINT *,R,PHI,VF
 STOP
END FUNCTION VF

