program split
implicit none
character(25)                           :: fName1, fName2, fName3, fff1, fff2
integer                                 :: i, j, k, fStat
integer, parameter                      :: inUn = 100, outUn1 = 201, outUn2 = 202
real(8), allocatable, dimension(:)      :: psiR, psiI, psi2, hR, hI, V
real(8), parameter                      :: h = 1.0D0
real(8)                                 :: L, m, sig, k0, dx, dt, tmax, N, norm, x0, x
real(8)                                 :: iniBarr, finBarr, valBarr, normBef, normAft
real(8)                                 :: h1, h2
integer                                 :: Nx, Nt, pstep, x_lim1, x_lim2
complex(8), allocatable, dimension(:)   :: psi, pot, kin


call get_command_argument(1, fName1, status=fStat)
if (fStat /= 0) then
        print*, 'Any file given ---> Exitting program (1)'
        call exit()
end if

call get_command_argument(2, fName2, status=fStat)
if (fStat /= 0) then
        print*, 'Any file given ---> Exitting program (2)'
        call exit()
end if

call get_command_argument(3, fName3, status=fStat)
if (fStat /= 0) then
        print*, 'Any file given ---> Exitting program (3)'
        call exit()
end if


fff1 = '(5F14.4)'
fff2 = '(I8, 2F10.4)'
open(unit=inUn  , file=trim(fName1))
open(unit=outUn1, file=trim(fName2))
open(unit=outUn2, file=trim(fName3))

read(inUn,*) L
read(inUn,*) m
read(inUn,*) sig
read(inUn,*) k0
read(inUn,*) x0
read(inUn,*) dx
read(inUn,*) dt
read(inUn,*) tmax
read(inUn,*) pstep
read(inUn,*) iniBarr
read(inUn,*) finBarr
read(inUn,*) valbarr

Nx = L/dx
Nt = tmax/dt

allocate(psi(Nx), V(Nx))
allocate(pot(Nx), kin(Nx))

! INITIALIZATION OF THE WAVE-PACKET:
call init_wp(sig, k0, Nx, L, psi)
call operators(Nx, dt, pot, kin, L, iniBarr, finBarr, valBarr, dx, x_lim1, x_lim2, V)
call fourier(0, Nx, psi)
!call fourier(-1, Nx, psi)


! INITIALIZE THE MAIN LOOP


do i = 1, Nt, 1
        psi(:) = psi(:)*pot(:)
        call fourier(1, Nx, psi)
        psi(:) = psi(:)*kin(:)
        call fourier(-1, Nx, psi)



        if (mod(i,pstep) == 0) then
                do j = 1, Nx, 1
                        x = dble(j)*dx
                        write(outUn1,fff1) x, V(j), abs(psi(j)), real(psi(j)), aimag(psi(j))
                end do
                norm = sum(dx*abs(psi)**2)
                !normBef = sum(dx*psi(1:x_lim1))
                !normAft = sum(dx*psi(x_lim2:))
                write(outUn2,fff2) i, dt*i, norm
                write(*,*) i, dt*i, norm
        end if
end do


contains

subroutine init_wp(sig, k0, Nx, L, psi)
implicit none
real(8), parameter                              :: pi = 3.141592653589793D0
real(8), intent(in)                             :: sig, k0, L
integer, intent(in)                             :: Nx
complex(8), dimension(Nx), intent(out)          :: psi
integer                                         :: i, j
real(8)                                         :: k, dk, nr, dx, x, N
real(8), dimension(Nx)                          :: psiR, psiI

N = 1.0D0/dsqrt(1.77156762398*sig)
dx = L/Nx
do j = 1, Nx, 1
        x = dx*dble(j)
        psiI(j) = N*dexp(-(x-x0)**2/(2*sig**2))*dsin(k0*x)
        psiR(j) = N*dexp(-(x-x0)**2/(2*sig**2))*dcos(k0*x)
end do

psiR(1) = 0.0D0; psiI(1) = 0.0D0
psiR(Nx) = 0.0D0; psiI(Nx) = 0.0D0
do j = 1, Nx, 1
        psi(j) = dcmplx(psiR(j), psiI(j))
end do


end subroutine init_wp

subroutine fourier(dir, Nx, psi)
implicit none
real(8), parameter                              :: pi = 3.141592653589793D0
integer, intent(in)                             :: dir, Nx
complex(8), dimension(Nx), intent(inout)        :: psi
integer                                         :: i
real(8), dimension(:), allocatable, save        :: wsave
real(8)                                         :: nr

if (dir == 1) then
        call dcfftf(Nx,psi,wsave)
        nr=1.d0/float(Nx)
        psi(:) = psi(:)*nr
else if (dir == -1) then
        call dcfftb(Nx,psi,wsave)
else if (dir == 0) then
        allocate(wsave(4*Nx+20))
        call dcffti(Nx,wsave)
end if
end subroutine fourier

subroutine operators(Nx, dt, pot, kin, L, iniBarr, finBarr, valBarr, dx, x_lim1, x_lim2, V)
implicit none
real(8), parameter                              :: pi = 3.141592653589793D0
integer, intent(in)                             :: Nx
real(8), intent(in)                             :: iniBarr, finBarr, valBarr
real(8), intent(in)                             :: dt, dx, L
complex(8), dimension(Nx), intent(out)          :: pot, kin 
real(8), dimension(Nx), intent(out)             :: V
integer                                         :: x_lim1, x_lim2
real(8)                                         :: dk, k, h1, h2, x, k0
integer                                         :: i, j

pot(:) = 1.0D0; kin(:) = 0.0D0
V(:) = 0.0D0
dk=2.d0*pi/L
k0 = 20
do j = 1, Nx, 1
        x = dble(j)*dx
        k = dble(j)*dk
        if ((x <= finBarr).and.(x >= iniBarr)) then
                V(j) = valBarr
                pot(j) = exp(-dt*(0,1)*(valBarr))
        end if
        if (abs(x - iniBarr) < h1) then
                h1 = abs(x - iniBarr)
                x_lim1 = j
        end if
        if (abs(x - finBarr) < h2) then
                h2 = abs(x - finBarr)
                x_lim2 = j
        end if
        kin(j)=exp(-dt*0.5d0*(0.0D0,1.0D0)*(k)**2)!**2/2)
end do
pot(1) = exp(-dt*(0,1)*(9999999))
pot(Nx) = exp(-dt*(0,1)*(9999999))

end subroutine operators

end program split
































!------------------------!

      SUBROUTINE DCFFTF(N,C,WSAVE)
!***BEGIN PROLOGUE  DCFFTF
!***DATE WRITTEN   790601   (YYMMDD)
!***REVISION DATE  860115   (YYMMDD)
!***CATEGORY NO.  J1A2
!***KEYWORDS  FOURIER TRANSFORM
!***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
!***PURPOSE  Forward transform of a complex, periodic sequence.
!***DESCRIPTION
!           From the book, "Numerical Methods and Software" by
!                D. Kahaner, C. Moler, S. Nash
!                Prentice Hall, 1988
!
!  Subroutine DCFFTF computes the forward complex discrete Fourier
!  transform (the Fourier analysis).  Equivalently, DCFFTF computes
!  the Fourier coefficients of a complex periodic sequence.
!  The transform is defined below at output parameter C.
!
!  The transform is not normalized.  To obtain a normalized transform
!  the output must be divided by N.  Otherwise a call of DCFFTF
!  followed by a call of DCFFTB will multiply the sequence by N.
!
!  The array WSAVE which is used by subroutine DCFFTF must be
!  initialized by calling subroutine DCFFTI(N,WSAVE).
!
!  Input Parameters
!
!
!  N      the length of the complex sequence C.  The method is
!         more efficient when N is the product of small primes.
!
!  C      a complex array of length N which contains the sequence
!
!  WSAVE   a d.p. work array which must be dimensioned at least 4*N+15
!          in the program that calls DCFFTF.  The WSAVE array must be
!          initialized by calling subroutine DCFFTI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.  Thus subsequent
!          transforms can be obtained faster than the first.
!          The same WSAVE array can be used by DCFFTF and DCFFTB.
!
!  Output Parameters
!
!  C      for J=1,...,N
!
!             C(J)=the sum from K=1,...,N of
!
!                   C(K)*EXP(-I*J*K*2*PI/N)
!
!                         where I=SQRT(-1)
!
!  WSAVE   contains initialization calculations which must not be
!          destroyed between calls of subroutine DCFFTF or DCFFTB
!
!  *   References                                                      *
!  *                                                                   *
!  *   1. P.N. Swarztrauber, Vectorizing the FFTs, in Parallel         *
!  *      Computations (G. Rodrigue, ed.), Academic Press, 1982,       *
!  *      pp. 51-83.                                                   *
!  *   2. B.L. Buzbee, The SLATEC Common Math Library, in Sources      *
!  *      and Development of Mathematical Software (W. Cowell, ed.),   *
!  *      Prentice-Hall, 1984, pp. 302-318.                            *
!  *                                                                   *
!  *********************************************************************
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DCFTF1
!***END PROLOGUE  DCFFTF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  C(10000) ,WSAVE(10000)
!***FIRST EXECUTABLE STATEMENT  DCFFTF
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL DCFTF1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END


      SUBROUTINE DCFFTB(N,C,WSAVE)
!***BEGIN PROLOGUE  DCFFTB
!***DATE WRITTEN   790601   (YYMMDD)
!***REVISION DATE  860115   (YYMMDD)
!***CATEGORY NO.  J1A2
!***KEYWORDS  FOURIER TRANSFORM
!***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
!***PURPOSE  Unnormalized inverse of DCFFTF.
!***DESCRIPTION
!           From the book, "Numerical Methods and Software" by
!                D. Kahaner, C. Moler, S. Nash
!                Prentice Hall, 1988
!
!  Subroutine DCFFTB computes the backward complex discrete Fourier
!  transform (the Fourier synthesis).  Equivalently, DCFFTB computes
!  a complex periodic sequence from its Fourier coefficients.
!  The transform is defined below at output parameter C.
!
!  A call of DCFFTF followed by a call of DCFFTB will multiply the
!  sequence by N.
!
!  The array WSAVE which is used by subroutine DCFFTB must be
!  initialized by calling subroutine DCFFTI(N,WSAVE).
!
!  Input Parameters
!
!
!  N      the length of the complex sequence C.  The method is
!         more efficient when N is the product of small primes.
!
!  C      a complex array of length N which contains the sequence
!
!  WSAVE   a d.p. work array which must be dimensioned at least 4*N+15
!          in the program that calls DCFFTB.  The WSAVE array must be
!          initialized by calling subroutine DCFFTI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.  Thus subsequent
!          transforms can be obtained faster than the first.
!          The same WSAVE array can be used by DCFFTF and DCFFTB.
!
!  Output Parameters
!
!  C      For J=1,...,N
!
!             C(J)=the sum from K=1,...,N of
!
!                   C(K)*EXP(I*J*K*2*PI/N)
!
!                         where I=SQRT(-1)
!
!  WSAVE   contains initialization calculations which must not be
!          destroyed between calls of subroutine DCFFTF or DCFFTB
!
!  *   References                                                      *
!  *                                                                   *
!  *   1. P.N. Swarztrauber, Vectorizing the FFTs, in Parallel         *
!  *      Computations (G. Rodrigue, ed.), Academic Press, 1982,       *
!  *      pp. 51-83.                                                   *
!  *   2. B.L. Buzbee, The SLATEC Common Math Library, in Sources      *
!  *      and Development of Mathematical Software (W. Cowell, ed.),   *
!  *      Prentice-Hall, 1984, pp. 302-318.                            *
!  *                                                                   *
!  *********************************************************************
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DCFTB1
!***END PROLOGUE  DCFFTB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(10000),WSAVE(10000)
!***FIRST EXECUTABLE STATEMENT  DCFFTB
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL DCFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END

      SUBROUTINE DCFFTI(N,WSAVE)
!***BEGIN PROLOGUE  DCFFTI
!***DATE WRITTEN   790601   (YYMMDD)
!***REVISION DATE  860115   (YYMMDD)
!***CATEGORY NO.  J1A2
!***KEYWORDS  FOURIER TRANSFORM
!***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
!***PURPOSE  Initialize for DCFFTF and DCFFTB.
!***DESCRIPTION
!           From the book, "Numerical Methods and Software" by
!                D. Kahaner, C. Moler, S. Nash
!                Prentice Hall, 1988
!
!  Subroutine DCFFTI initializes the array WSAVE which is used in
!  both DCFFTF and DCFFTB.  The prime factorization of N together with
!  a tabulation of the trigonometric functions are computed and
!  stored in WSAVE.
!
!  Input Parameter
!
!  N       the length of the sequence to be transformed
!
!  Output Parameter
!
!  WSAVE   a work array which must be dimensioned at least 4*N+15.
!          The same work array can be used for both DCFFTF and DCFFTB
!          as long as N remains unchanged.  Different WSAVE arrays
!          are required for different values of N.  The contents of
!          WSAVE must not be changed between calls of DCFFTF or DCFFTB.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DCFTI1
!***END PROLOGUE  DCFFTI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  WSAVE(10000)
!***FIRST EXECUTABLE STATEMENT  DCFFTI
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL DCFTI1 (N,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
      SUBROUTINE DCFTF1(N,C,CH,WA,IFAC)
!***BEGIN PROLOGUE  DCFTF1
!***REFER TO  DCFFTF
!***ROUTINES CALLED  DPASSF,DPASF2,DPASF3,DPASF4,DPASF5
!***END PROLOGUE  DCFTF1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CH(10000),C(10000),WA(10000),IFAC(10000)
!***FIRST EXECUTABLE STATEMENT  DCFTF1
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL DPASF4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL DPASF4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL DPASF2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL DPASF2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL DPASF3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL DPASF3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL DPASF5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL DPASF5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL DPASSF (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL DPASSF (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END

      SUBROUTINE DCFTI1(N,WA,IFAC)
!***BEGIN PROLOGUE  DCFTI1
!***REFER TO  DCFFTI
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DCFTI1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WA(10000),IFAC(10000),NTRYH(10000)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/3,4,2,5/
!***FIRST EXECUTABLE STATEMENT  DCFTI1
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 8.0D0*ATAN(1.0D0)
      ARGH = TPI/DBLE(N)
      I = 2
      L1 = 1
      DO 110 K1=1,NF
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO+IDO+2
         IPM = IP-1
         DO 109 J=1,IPM
            I1 = I
            WA(I-1) = 1.0D0
            WA(I) = 0.0D0
            LD = LD+L1
            FI = 0.0D0
            ARGLD = DBLE(LD)*ARGH
            DO 108 II=4,IDOT,2
               I = I+2
               FI = FI+1.0D0
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IF (IP .LE. 5) GO TO 109
            WA(I1-1) = WA(I-1)
            WA(I1) = WA(I)
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END

      SUBROUTINE DPASF2(IDO,L1,CC,CH,WA1)
!***BEGIN PROLOGUE  DPASF2
!***REFER TO  DCFFTF
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DPASF2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CC(IDO,2,L1),CH(IDO,L1,2)
      DIMENSION  WA1(10000)
!***FIRST EXECUTABLE STATEMENT  DPASF2
      IF (IDO .GT. 2) GO TO 102
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
      RETURN
  102 IF(IDO/2.LT.L1) GO TO 105
      DO 104 K=1,L1
!DIR$ IVDEP
         DO 103 I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
      RETURN
  105 DO 107 I=2,IDO,2
!DIR$ IVDEP
      DO 106 K=1,L1
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
  106    CONTINUE
  107 CONTINUE
      RETURN
      END
      SUBROUTINE DPASF3(IDO,L1,CC,CH,WA1,WA2)
!***BEGIN PROLOGUE  DPASF3
!***REFER TO  DCFFTF
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DPASF3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CC(IDO,3,L1),CH(IDO,L1,3)
      DIMENSION  WA1(10000),WA2(10000)
!***FIRST EXECUTABLE STATEMENT  DPASF3
      TAUR = -.5D0
      TAUI = -.5D0*SQRT(3.0D0)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TR2 = CC(1,2,K)+CC(1,3,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         TI2 = CC(2,2,K)+CC(2,3,K)
         CI2 = CC(2,1,K)+TAUR*TI2
         CH(2,K,1) = CC(2,1,K)+TI2
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
         CH(2,K,2) = CI2+CR3
         CH(2,K,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 IF(IDO/2.LT.L1) GO TO 105
      DO 104 K=1,L1
!DIR$ IVDEP
         DO 103 I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
  103    CONTINUE
  104 CONTINUE
      RETURN
  105 DO 107 I=2,IDO,2
!DIR$ IVDEP
         DO 106 K=1,L1
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
  106    CONTINUE
  107 CONTINUE
      RETURN
      END
      SUBROUTINE DPASF4(IDO,L1,CC,CH,WA1,WA2,WA3)
!***BEGIN PROLOGUE  DPASF4
!***REFER TO  DCFFTF
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DPASF4
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CC(IDO,4,L1),CH(IDO,L1,4)
      DIMENSION WA1(10000),WA2(10000),WA3(10000)
!***FIRST EXECUTABLE STATEMENT  DPASF4
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,2,K)-CC(2,4,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,4,K)-CC(1,2,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 IF(IDO/2.LT.L1) GO TO 105
      DO 104 K=1,L1
!DIR$ IVDEP
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,2,K)-CC(I,4,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,4,K)-CC(I-1,2,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
  105 DO 107 I=2,IDO,2
!DIR$ IVDEP
         DO 106 K=1,L1
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,2,K)-CC(I,4,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,4,K)-CC(I-1,2,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
  106    CONTINUE
  107 CONTINUE
      RETURN
      END
      SUBROUTINE DPASF5(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
!***BEGIN PROLOGUE  DPASF5
!***REFER TO  DCFFTF
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DPASF5
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CC(IDO,5,L1),CH(IDO,L1,5)
      DIMENSION  WA1(10000),WA2(10000),WA3(10000),WA4(10000)
!***FIRST EXECUTABLE STATEMENT  DPASF5
      PI = 4.0D0*ATAN(1.0D0)
      TR11 = SIN(.1D0*PI)
      TI11 = -SIN(.4D0*PI)
      TR12 = -SIN(.3D0*PI)
      TI12 = -SIN(.2D0*PI)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 IF(IDO/2.LT.L1) GO TO 105
      DO 104 K=1,L1
!DIR$ IVDEP
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
  105 DO 107 I=2,IDO,2
!DIR$ IVDEP
         DO 106 K=1,L1
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
  106    CONTINUE
  107 CONTINUE
      RETURN
      END
      SUBROUTINE DPASSF(NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
!***BEGIN PROLOGUE  DPASSF
!***REFER TO  DCFFTF
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DPASSF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  CH(IDO,L1,IP),CC(IDO,IP,L1)
      DIMENSION C1(IDO,L1,IP),WA(10000) ,C2(IDL1,IP)
      DIMENSION CH2(IDL1,IP)
!***FIRST EXECUTABLE STATEMENT  DPASSF
      IDOT = IDO/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
!
      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
!DIR$ IVDEP
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
!DIR$ IVDEP
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
!DIR$ IVDEP
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
!DIR$ IVDEP
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
!DIR$ IVDEP
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = -WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
!DIR$ IVDEP
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
!DIR$ IVDEP
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
!DIR$ IVDEP
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
!DIR$ IVDEP
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
!DIR$ IVDEP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
!DIR$ IVDEP
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
!DIR$ IVDEP
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END


      SUBROUTINE DCFTB1(N,C,CH,WA,IFAC)
!***BEGIN PROLOGUE  DCFTB1
!***REFER TO  DCFFTB
!***ROUTINES CALLED  DPASSB,DPASB2,DPASB3,DPASB4,DPASB5
!***END PROLOGUE  DCFTB1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CH(10000),C(10000),WA(10000),IFAC(10000)
!***FIRST EXECUTABLE STATEMENT  DCFTB1
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL DPASB4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL DPASB4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL DPASB2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL DPASB2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL DPASB3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL DPASB3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL DPASB5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL DPASB5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL DPASSB (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL DPASSB (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END
      SUBROUTINE DPASB4(IDO,L1,CC,CH,WA1,WA2,WA3)
!***BEGIN PROLOGUE  DPASB4
!***REFER TO  DCFFTB
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DPASB4
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CC(IDO,4,L1),CH(IDO,L1,4)
      DIMENSION WA1(10000) ,WA2(10000),WA3(10000)
!***FIRST EXECUTABLE STATEMENT  DPASB4
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,4,K)-CC(2,2,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,2,K)-CC(1,4,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 IF(IDO/2.LT.L1) GO TO 105
      DO 104 K=1,L1
!DIR$ IVDEP
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,4,K)-CC(I,2,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,2,K)-CC(I-1,4,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
  105 DO 107 I=2,IDO,2
!DIR$ IVDEP
         DO 106 K=1,L1
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,4,K)-CC(I,2,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,2,K)-CC(I-1,4,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
  106    CONTINUE
  107 CONTINUE
      RETURN
      END
      SUBROUTINE DPASB5(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
!***BEGIN PROLOGUE  DPASB5
!***REFER TO  DCFFTB
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DPASB5
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CC(IDO,5,L1),CH(IDO,L1,5)
      DIMENSION WA1(10000),WA2(10000),WA3(10000),WA4(10000)
!***FIRST EXECUTABLE STATEMENT  DPASB5
      PI = 4.0D0*ATAN(1.0D0)
      TR11 = SIN(.1D0*PI)
      TI11 = SIN(.4D0*PI)
      TR12 = -SIN(.3D0*PI)
      TI12 = SIN(.2D0*PI)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 IF(IDO/2.LT.L1) GO TO 105
      DO 104 K=1,L1
!DIR$ IVDEP
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
  105 DO 107 I=2,IDO,2
!DIR$ IVDEP
         DO 106 K=1,L1
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
  106    CONTINUE
  107 CONTINUE
      RETURN
      END

      SUBROUTINE DPASSB(NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
!***BEGIN PROLOGUE  DPASSB
!***REFER TO  DCFFTB
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DPASSB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION   CH(IDO,L1,IP),CC(IDO,IP,L1)
      DIMENSION  C1(IDO,L1,IP),WA(10000),C2(IDL1,IP)
      DIMENSION  CH2(IDL1,IP)
!***FIRST EXECUTABLE STATEMENT  DPASSB
      IDOT = IDO/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
!
      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
!DIR$ IVDEP
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
!DIR$ IVDEP
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
!DIR$ IVDEP
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
!DIR$ IVDEP
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
!DIR$ IVDEP
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
!DIR$ IVDEP
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
!DIR$ IVDEP
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
!DIR$ IVDEP
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
!DIR$ IVDEP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
!DIR$ IVDEP
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
!DIR$ IVDEP
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END
      SUBROUTINE DPASB2(IDO,L1,CC,CH,WA1)
!***BEGIN PROLOGUE  DPASB2
!***REFER TO  DCFFTB
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DPASB2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CC(IDO,2,L1),CH(IDO,L1,2)
      DIMENSION WA1(10000)
!***FIRST EXECUTABLE STATEMENT  DPASB2
      IF (IDO .GT. 2) GO TO 102
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
      RETURN
  102 IF(IDO/2.LT.L1) GO TO 105
      DO 104 K=1,L1
!DIR$ IVDEP
         DO 103 I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
      RETURN
  105 DO 107 I=2,IDO,2
!DIR$ IVDEP
         DO 106 K=1,L1
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
  106    CONTINUE
  107 CONTINUE
      RETURN
      END
      SUBROUTINE DPASB3(IDO,L1,CC,CH,WA1,WA2)
!***BEGIN PROLOGUE  DPASB3
!***REFER TO  DCFFTB
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DPASB3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  CC(IDO,3,L1),CH(IDO,L1,3)
      DIMENSION WA1(10000),WA2(10000)
!***FIRST EXECUTABLE STATEMENT  DPASB3
      TAUR = -.5D0
      TAUI = .5D0*SQRT(3.0D0)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TR2 = CC(1,2,K)+CC(1,3,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         TI2 = CC(2,2,K)+CC(2,3,K)
         CI2 = CC(2,1,K)+TAUR*TI2
         CH(2,K,1) = CC(2,1,K)+TI2
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
         CH(2,K,2) = CI2+CR3
         CH(2,K,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 IF(IDO/2.LT.L1) GO TO 105
      DO 104 K=1,L1
!DIR$ IVDEP
         DO 103 I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
  103    CONTINUE
  104 CONTINUE
      RETURN
  105 DO 107 I=2,IDO,2
!DIR$ IVDEP
         DO 106 K=1,L1
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
  106    CONTINUE
  107 CONTINUE
      RETURN
      END

