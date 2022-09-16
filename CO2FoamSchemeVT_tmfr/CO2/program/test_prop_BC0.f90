PROGRAM test_pressure

! Short test program for the function: 'Helmholtz.f90'.
!
! Pressure corresponds to: -(df / dv)_T.
!
! Let us check if pressure values from tables at the end of the paper
! (from page 1562 on) are obtained using 'Helmholtz.f90'.
!
! The derivative is approximated by a 2nd order finite difference.
!
      USE def_constants
      USE Grid
      USE non_linear_solvers
      USE properties
      USE solver_eos
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: Ntab = 27
      INTEGER :: Niter, exitflag, i   
      REAL(pr) :: p_in, v_in, out2, out3, resnorm, vguess, e_in  

      REAL(pr), ALLOCATABLE :: tabT(:)

      ALLOCATE(tabT(Ntab))
     
     DO i = 1, Ntab
        tabT(i) = 336.d0 - i + 1
     END DO

     p_in = 8.d6
     vguess = 1.d0/200 ! Initial guess

     PRINT*, "  rho (kg/m3)", "               e (J/kg)"
     PRINT*, " "

OPEN (UNIT = 22,             FILE = 'e_SW.txt', &
&   FORM = 'formatted',     ACTION = 'write',   &
&   STATUS = 'replace')

      DO i = 1, Ntab
         
         CALL New_Rap1D(3, v_in, out2, resnorm, Niter,&
         &                             exitflag, p_in, vguess, tabT(i), out3)    ! Here output is v_in as a function of ranged T at fixed p 
         CALL inter_energy(tabT(i),v_in,e_in)    ! Here output is e_in=f(T,v)
         PRINT*, 1/v_in, e_in + 5.d5
         WRITE(22,*) e_in + 5.d5
      
      ENDDO
CLOSE(UNIT = 22)

END PROGRAM test_pressure
