!INCLUDE 'MATHIO.F90'
!INCLUDE 'TENSOR.F90'
! ############### CONST ###################
MODULE CONST
	COMPLEX, PARAMETER :: Z0 = (0.,0.), Z1 = (1.,0.), ZI = (0.,1.)
	REAL, PARAMETER :: PI = 4*ATAN(1.)
END MODULE CONST
! ############### MODEL ###################
MODULE MODEL
	INTEGER :: LX = 32
	INTEGER :: LY = 8
	INTEGER :: LA = 0
	REAL :: SPEED = 4.
END MODULE MODEL
! ############# PHYSICS ###################
MODULE PHYSICS
CONTAINS
! ---------- bilayer Hamiltonian ------------
! set Hamiltonian
FUNCTION SET_H(L1, L2) RESULT(H)
! in: L1, L2 - layer index
! out: H - Hamiltonian matrix (4LX x 4LX)
! H =
! [ HA, HB]
! [ HC, HD]
	USE CONST
	USE MODEL
	INTEGER, INTENT(IN) :: L1, L2
	COMPLEX, TARGET :: H(4*LX,4*LX)
	! local variable
	INTEGER :: M, M1
	COMPLEX, POINTER :: HA(:,:), HB(:,:), HC(:,:), HD(:,:)
	
	H = Z0 ! initialization
	! set pointer to blocks
	HA => H(:2*LX,:2*LX)
	HB => H(:2*LX,2*LX+1:)
	HC => H(2*LX+1:,:2*LX)
	HD => H(2*LX+1:,2*LX+1:)
	! set diagonal block
	CALL SET_H0(HA, L1)
	CALL SET_H0(HD, L2)
	! set off-diagonal block (coupling between layers)
	CALL SET_H1(HB, L1, L2)
	HC = TRANSPOSE(CONJG(HB)) ! Hermitian conj
END FUNCTION SET_H
! set diagonal block Hamiltonian
SUBROUTINE SET_H0(H0, L)
! inout: H0 - diagonal Hamiltonian to be set
! in: L - layer index
	USE MODEL
	COMPLEX, INTENT(INOUT) :: H0(:,:)
	INTEGER, INTENT(IN) :: L
	! local variables
	INTEGER :: M, SIGN
	
	! set layer velocity sign, according to the pairy
	SIGN = (-1)**L
	IF (L <= LA) THEN
		! long ring layer
		FORALL (M = 1:2*LX)
			H0(M, M) = SIGN*EH(M)
		END FORALL
	ELSE
		! short rings layer
		FORALL (M = 2:2*LX:2)
			H0(M/2, M/2) = SIGN*EH(M)	
			H0(LX+M/2, LX+M/2) = H0(M/2, M/2)	
		END FORALL
	END IF
END SUBROUTINE SET_H0
! set off-diagonal block Hamiltonian
SUBROUTINE SET_H1(H1, L1, L2)
! inout: H1- off-diagonali block to be set
! in: L1, L2 - layer indices
	USE CONST
	USE MODEL
	COMPLEX, INTENT(INOUT) :: H1(:,:)
	INTEGER, INTENT(IN) :: L1, L2
	! local variables
	INTEGER :: M, M1
	
	IF (L1 <= LA .EQV. L2 <= LA) THEN
		! both layers are of the same kind, coupled trivially
		FORALL (M = 1:2*LX)
			H1(M,M) = Z1
		END FORALL
	ELSE ! layers are not of the same kind
		IF (L1 <= LA) THEN
			! L1 long, L2 short, from short to long
			! overlap <k1,k> transfer |k1> to |k>
			FORALL (M = 1:2*LX, M1 = 2:2*LX:2)
				H1(M, M1/2) = KOVERLAP(M, M1) ! cuff 1
				H1(M, LX+M1/2) = (-1)**(M-M1)*H1(M, M1/2) ! cuff 2
			END FORALL
		ELSE
			! L1 short, L2 long, from long to short
			! overlap <k1,k>* transfer |k> to |k1>
			FORALL (M = 1:2*LX, M1 = 2:2*LX:2)
				H1(M1/2, M) = CONJG(KOVERLAP(M, M1)) ! cuff 1
				H1(LX+M1/2, M) = (-1)**(M-M1)*H1(M1/2, M) ! cuff 2
			END FORALL		
		END IF
	END IF
	! use negative coupling to stabilize the phase
	H1 = - H1
END SUBROUTINE SET_H1
! energy function
PURE FUNCTION EH(M) RESULT(Z)
! input: momentum index M
! output: energy
	USE CONST
	USE MODEL
	INTEGER, INTENT(IN) :: M
	COMPLEX :: Z
	
	Z = CMPLX((SPEED/LX)*(MODULO(M+LX,2*LX)-LX))
END FUNCTION EH
! momentum overlap
PURE FUNCTION KOVERLAP(M, M1) RESULT(Z)
! input: momentum index M, M1
! output: overlap as a complex number
	USE CONST
	USE MODEL
	INTEGER, INTENT(IN) :: M, M1
	COMPLEX :: Z
	! local variables
	INTEGER :: N
	
	N = M - M1
	IF (MODULO(N,2) == 0) THEN ! even
		IF (N == 0) THEN
			Z = Z1/SQRT(2.)
		ELSE
			Z = Z0
		END IF
	ELSE ! odd
		Z = EXP(ZI*PI*N/LX)
		Z = Z*(1-Z**LX)/(1-Z)/LX/SQRT(2.)
	END IF
END FUNCTION KOVERLAP
! diagonalization of Hamiltonian
FUNCTION GET_OCC(H) RESULT(OCC)
! input: H - Hamiltonian
! output: OCC - occupied states
	USE MODEL
	COMPLEX, INTENT(IN) :: H(4*LX,4*LX)
	COMPLEX :: OCC(4*LX,2*LX)
	! local variable
	INTEGER :: N, LWORK,INFO,I
	COMPLEX, ALLOCATABLE :: A(:,:), WORK(:)
	REAL, ALLOCATABLE :: W(:), RWORK(:)
	INTEGER, ALLOCATABLE :: INDS(:)
	
	N = 4*LX ! get size
	A = H ! move data
	! prepare for diagonalization
	LWORK = 65*N
	ALLOCATE(W(N),WORK(LWORK),RWORK(3*N))
	CALL ZHEEV('V','U',N,A,N,W,WORK,LWORK,RWORK,INFO)
	IF (INFO>0) PRINT *,'GET_OCC::diag: diagonalization failed to converge.'
	! pick occupied states
	INDS = PACK([(I,I=1,N)],W < 0)
	! check the spectrum is gapped from half
	IF (SIZE(INDS) /= N/2) THEN
		PRINT '(A,2I5)', 'GET_OCC::pick: number of occupied state is not half of the Hilbert space dim, check if there are gapless modes.', SIZE(INDS), N/2
		STOP
	END IF
	! take occupied states
	OCC = A(:,INDS)
END FUNCTION GET_OCC
! ------------ overlap matrix --------------
! construct overlap matrix 
FUNCTION GET_A() RESULT(A)
! out: A - overlap matrix
	USE CONST
	USE MODEL
	COMPLEX :: A(LX*LY,LX*LY)
	! local variables
	COMPLEX :: OCC0(4*LX,2*LX), OCC1(4*LX,2*LX)
	INTEGER :: I0, I1, L, D
	
	D = 2*LX
	A = Z0 ! clear
	! initialize OCC1
	I0 = -1
	I1 = 0
	OCC1 = GET_OCC(SET_H(LY,1))
	! begin filling A
	DO L = 1, LY
		IF (MODULO(L,2) == 1) THEN ! L odd
			I0 = MODULO(I0 + 1, LY/2)
			OCC0 = GET_OCC(SET_H(L,MODULO(L,LY)+1))
			A(I0*D+1:(I0+1)*D,I1*D+1:(I1+1)*D) = MATMUL(TRANSPOSE(CONJG(OCC0(:D,:))),OCC1(D+1:,:))
		ELSE ! L even
			I1 = MODULO(I1 + 1, LY/2)
			OCC1 = GET_OCC(SET_H(L,MODULO(L,LY)+1))
			A(I0*D+1:(I0+1)*D,I1*D+1:(I1+1)*D) = MATMUL(TRANSPOSE(CONJG(OCC0(D+1:,:))),OCC1(:D,:))
		END IF
	END DO
END FUNCTION GET_A
! cal Renyi entropy from overlap matrix by SVD
FUNCTION GET_H2() RESULT(H2)
! out: S - Renyi entropy
	USE MODEL
	REAL :: H2
	! local variable
	INTEGER :: N, LWORK, INFO
	COMPLEX, ALLOCATABLE :: A(:,:), U(:,:), VT(:,:), WORK(:)
	REAL, ALLOCATABLE :: S(:), RWORK(:)
	INTEGER, ALLOCATABLE :: IWORK(:)
	REAL, PARAMETER :: TOL = 1.E-8
	
	A = GET_A() ! get the overlap matrix
	N = LX*LY ! array dim
	! cal SVD
	LWORK = 9*N
	! solver 1 ++++++++++++++
!	ALLOCATE(S(N),U(0,0),VT(0,0),WORK(LWORK),RWORK(5*N))
!	CALL ZGESVD('N','N',N,N,A,N,S,U,N,VT,N,WORK,LWORK,RWORK,INFO)
	! solver 2 ++++++++++++++
	ALLOCATE(S(N),U(0,0),VT(0,0),WORK(LWORK),RWORK(5*N),IWORK(8*N))
	CALL ZGESDD('N',N,N,A,N,S,U,N,VT,N,WORK,LWORK,RWORK,IWORK,INFO)
	IF (INFO /= 0) PRINT *,'GET_H2::svd: SVD can not converge, result is not reliable.'
	PRINT '(1000F10.6)', S
	H2 = SUM(LOG(PACK(S,S>TOL)))
END FUNCTION GET_H2
! end of module PHYSICS
END MODULE PHYSICS
! ############### TASK ####################
MODULE TASK
CONTAINS
! ----------- data -----------------
SUBROUTINE COLLECT()
	USE MODEL
	USE PHYSICS
	USE MATHIO
	! local variables
	REAL :: S(0:LY)
	
	DO LA = 0, LY
		PRINT *, LA
		S(LA) = GET_H2()
	END DO
	CALL EXPORT('S',S)
END SUBROUTINE COLLECT
! ----------- tests -----------------
! test
SUBROUTINE TEST()
	USE MODEL
	USE PHYSICS
	USE MATHIO
	REAL :: H2
	
	H2 = GET_H2()
	PRINT *, H2
END SUBROUTINE TEST
! test occupation space
SUBROUTINE TEST_OCC()
	USE MODEL
	USE PHYSICS
	USE MATHIO
	COMPLEX, ALLOCATABLE :: H(:,:), OCC(:,:)
	
	H = SET_H(8,1)
	CALL EXPORT('H',H)
	OCC = GET_OCC(H)
	CALL EXPORT('OCC',OCC)
END SUBROUTINE TEST_OCC
! test overlap matrix
SUBROUTINE TEST_A()
	USE MODEL
	USE PHYSICS
	USE MATHIO
	COMPLEX, ALLOCATABLE :: A(:,:)
	
	A = GET_A()
	CALL EXPORT('A',A)
END SUBROUTINE TEST_A
! end of module task
END MODULE TASK
! ############## PROGRAM ##################
PROGRAM MAIN
	USE TASK
	PRINT *, '------------ QSH -------------'
	
	CALL COLLECT()
!	CALL TEST()
!	CALL TEST_OCC()
!	CALL TEST_A()
END PROGRAM MAIN