PROGRAM write_read_write
!      
      USE def_constants
      USE def_variables
      USE non_linear_solvers
      USE properties
      
      IMPLICIT NONE

      INTEGER, PARAMETER :: Ncells=4557012
      INTEGER :: i, Niter, exitflag
      REAL(pr) :: v_in, p_in, T_in, vguess, out2, out3, resnorm
     !REAL(pr) :: cp
     !REAL(pr) :: p_out, T_out, c_out, x_out, vp_out, xxx
      
      REAL(pr), DIMENSION(Ncells) ::  pcells, Tcells, rhocells_CP, rhocells_SW
      
 !   CALL MAKE_GRID()

      DO i = 1, Ncells
         OPEN(UNIT=22,FILE='pcells.txt',STATUS='old')
         READ(22,*) pcells(i)
      ENDDO
      DO i = 1, Ncells
         OPEN(UNIT=22,FILE='Tcells.txt',STATUS='old')
         READ(22,*) Tcells(i)
      ENDDO
      DO i = 1, Ncells
         OPEN(UNIT=22,FILE='rhocells_CP.txt',STATUS='old')
         READ(22,*) rhocells_CP(i)
      ENDDO
!     DO i = 1, Ncells
!        PRINT*, pcells(i), Tcells(i), rhocells_CP(i)
!     ENDDO

DO i = 1, Ncells

   p_in = pcells(i)
   T_in = Tcells(i)
   vguess = 1/rhocells_CP(i)
   CALL New_Rap1D(3, v_in, out2, resnorm, Niter,&
     &                             exitflag, p_in, vguess, T_in, out3)     
   rhocells_SW(i) = 1/v_in 
!  PRINT*, (rhocells_CP(i) - rhocells_SW(i))/rhocells_CP(i)

ENDDO

! Save cells
 
OPEN (UNIT = 22,             FILE = 'rhocells_SW.txt', &
&   FORM = 'formatted',     ACTION = 'write',   &
&   STATUS = 'replace')

DO i = 1, Ncells
      WRITE(22,*) rhocells_SW(i)
ENDDO

 CLOSE(UNIT = 22)

END PROGRAM write_read_write
