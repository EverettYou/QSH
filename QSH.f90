!INCLUDE 'MATHIO.F90'
!INCLUDE 'TENSOR.F90'
! ############### CONST ###################
MODULE CONST
	COMPLEX, PARAMETER :: Z0 = (0.,0.), Z1 = (1.,0.), ZI = (0.,1.)
	REAL, PARAMETER :: PI = 4*ATAN(1.)
END MODULE CONST
! ############### MODEL ###################
MODULE MODEL
	INTEGER :: LX = 16
END MODULE MODEL
! ############### TASK ####################
MODULE TASK
CONTAINS
! test
SUBROUTINE TEST()

END SUBROUTINE TEST
! end of module task
END MODULE TASK
! ############## PROGRAM ##################
PROGRAM MAIN
	USE TASK
	PRINT *, '------------ QSH -------------'
	
	CALL TEST()
END PROGRAM MAIN