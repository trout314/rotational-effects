PROGRAM CCGEN
  !This program calculates the Clenshaw-Curtiss points and weights for
  !evaluating the integral I = int {from -1 to 1} F[x]dx.  This method
  !was developed by C. W. Clenshaw and A. R. Curtis, Numer. Math.  2
  !(1960) 197-205, and is discussed further by F. J. Smith, Numerische
  !Mathematik 7 (1965) 406-411 and H. O'Hara and F. J. Smith, J.
  !Computat. Phys. 5 (1970) 328-344.
  !
  !According to the Clenshaw-Curtis quadrature method, the N-th
  !approximation to I is
  !             I^N = sum {from s=0 to N} w_s^N F[x_s^N].
  !The Clenshaw-Curtis quadrature points and weights are
  !             x_s^N = cos(pi*s/N).
  !             w_0^N = h_N^N = 1/(N*N-1)            and
  !h_s^N = [2(-1)^s]/[N*N-1] + [4/N]* sin (pi*s/N) * sum {from i=1 to
  !  N/2} [sin((2i-1)*pi*s/N/(2i-1)] for  1<=s<=N-1
  !
  !The error estimates are given in terms of the a_N-2r quantities in
  !eq. (10) of the paper by O'Hara and Smith.  I have rewritten these
  !so that a_0 = sum {from s=0 to N } w0_s^N F[x_s^N]
  !        a_2 = sum {from s=0 to N } w2_s^N F[x_s^N]
  !        a_4 = sum {from s=0 to N } w4_s^N F[x_s^N]
  !where
  !   w0_s^N = (-1)^s (1/N)  if s=0 or s=N
  !   w0_s^N = (-1)^s (2/N)  if 0<s<N
  !   w2_s^N = (-1)^s (1/N) cos (2*pi*s/N) if s=0 or s=N
  !   w2_s^N = (-1)^s (2/N) cos (2*pi*s/N) if 0<s<N
  !   w4_s^N = (-1)^s (1/N) cos (4*pi*s/N) if s=0 or s=N
  !   w4_s^N = (-1)^s (2/N) cos (4*pi*s/N) if 0<s<N
  !The points and weights are written to file CCPTS&WTS.  However, the
  !a-N are not written, as they have not proved useful for integrals
  !that arise in the kinetic theory of gaseous ions.
  IMPLICIT NONE

  !Declare the integer variables.
  !  NNMAXCH is the maximum order of approximation in using the Curtiss-
  !    Clenshaw method of integration with NPTSCH pivots.
  !  LENGTH is the number of characters in CCFILE.
  !  ERROR allows for recovery from I/O errors.
  !  I, J, N and S are do-loop indices.
  INTEGER::NMAXCH,LENGTH,ERROR,I,J,N,S

  !Declare the real variables.
  !  PI is the usual 3.141592...
  !  ZERO throuigh FOUR are double precision versions of 0-4.
  !  CHEPTS is a quadrature point.
  !  CHEWTS is a quadrature weight.
  !  WT0 is a weight for the error estimate of a_0.
  !  WT2 is a weight for the error estimate of a_2.
  !  WT4 is a weight for the error estimate of a_4.
  !  Y is temporary storage for an intermediate quantity.
  REAL(SELECTED_REAL_KIND(18))::PI,ZERO=0,ONE=1,TWO=2,FOUR=4,CHEPTS,CHEWTS,&
    WT0,WT2,WT4,Y

  !Declare the character array.
  !  CCFILE is the fully qualified name of the output file.
  CHARACTER(LEN=78)::CCFILE

  !Open the input file.
  OPEN(UNIT=10,FILE='CCgen.in',STATUS='OLD',IOSTAT=ERROR)
  IF(ERROR/=0)THEN
    PRINT *,'Could not open the old file named CCgen.in.'
    STOP
  ENDIF

  !Read NMAXCH and CCFILE.
  READ(10,*)NMAXCH
  READ(10,'(A)')CCFILE
  CLOSE(10)

  !Determine LENGTH
  LENGTH=INDEX(CCFILE//'  ','  ')-1

  !Establish the value for PI.
  PI=FOUR*ATAN(ONE)

  !Open the output file.
  OPEN(10,FILE=CCFILE,IOSTAT=ERROR)
  IF(ERROR/=0)THEN
    PRINT *,'Could not open file '//CCFILE(1:LENGTH)
    PRINT *,'Do you have permission to write to this file?'
    STOP
  ENDIF

  !Calculate the points and weights for orders 1, 2 .. NMAXCH,
  !corresponding to N=2**1, 2**2, ... 2**NMAXCH points.  The loop
  !for S runs backwards in order to get the points in ascending
  !order
  DO J=1,NMAXCH
    N=2**J
    DO S=N,0,-1
      CHEPTS=COS(PI*S/N)
      IF(S==0.OR.S==N)THEN
        CHEWTS=ONE/(N*N-1)
        WT0=ONE/N
        WT2=ONE/N
        WT4=ONE/N
      ELSE
        CHEWTS=TWO/(N*N-1)
        WT0=TWO/N
        WT2=TWO/N
        WT4=TWO/N
        IF(MOD(S,2)/=0)THEn
          CHEWTS=-CHEWTS
          WT0=-WT0
          WT2=-WT2
          WT4=-WT4
        ENDIf
        Y=ZERO
        DO I=1,N/2
          Y=Y+SIN(DBLE(2*I-1)*PI*S/N)/(2*I-1)
        ENDDO
        CHEWTS=CHEWTS+FOUR/N*SIN(PI*S/N)*Y
        WT2=WT2*COS(TWO*PI*S/N)
        WT4=WT4*COS(FOUR*PI*S/N)
      ENDIF
      !Choose one of the next two lines, by adding ! in front of the
      !other.  Both lines print values for the quadrature points
      !and weights, but the first also prints the WTn values.
      WRITE(10,'(3(I5,1X),2F23.20)')J,N,S,CHEPTS,CHEWTS
    ENDDO
  ENDDO
END PROGRAM CCGen
