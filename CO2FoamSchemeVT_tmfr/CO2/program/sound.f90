PROGRAM sound
     
     USE non_linear_solvers
     USE properties
     USE Grid
     IMPLICIT NONE

! Without table

!      INTEGER::Niter,flag
!      REAL(8)::vliq,vvap,Tsat,res,p,in_2,vliq_guess,vvap_guess,T_guess

!p = 6.d6
!vliq_guess = 1.d0/740 
!vvap_guess = 1.d0/219
!T_guess = 296
!          CALL New_Rap3D(3,vliq,vvap,Tsat,&
!      &        res, Niter, flag, p,in_2,&
!      &        vliq_guess,vvap_guess, T_guess)

!print*, 'vliq   = ', 1/vliq
!print*, 'vvap   = ', 1/vvap
!print*, 'Tsat   = ', Tsat

INTEGER :: i
INTEGER, PARAMETER :: N = 100
REAL(8) :: press,Tsat,vv,vl,ev,el,ratio
REAL(8), DIMENSION(N) :: duL_dp, duV_dp, dvL_dp, dvV_dp, du_dp_x, dv_dp_x
REAL(8), DIMENSION(N) :: alpha,qual,cHEM
REAL(8), DIMENSION(N) :: rho  

press = 5.d6

! Using table determine el ev vl vv at pressure press 

CALL MAKE_GRID()
CALL satprop(3, press, Tsat, vv, vl, ev, el)

ratio = ((ev - el)/(vv - vl))   ! (J/kg)/(m3/kg)

DO i = 1, N
     alpha(i) = dfloat(i) / N 
     qual(i) = ( alpha(i)/vv ) / ( alpha(i)/vv + (1-alpha(i))/vl )
     CALL satderiv(3, press, duL_dp(i), duV_dp(i), dvL_dp(i), dvV_dp(i))
     rho(i) = alpha(i)/(qual(i)+1.d-10)/vv
     du_dp_x(i) = qual(i) * duV_dp(i) + (1.d0 - qual(i)) * duL_dp(i)
     dv_dp_x(i) = qual(i) * dvV_dp(i) + (1.d0 - qual(i)) * dvL_dp(i)
     cHEM(i) = SQRT((press + ratio)/(rho(i)**2*(du_dp_x(i) - ratio * dv_dp_x(i))))
     PRINT*, i, alpha(i), qual(i), cHEM(i)
ENDDO

PRINT*, 'rhol   = ', 1/vl
PRINT*, 'rhov   = ', 1/vv
PRINT*, 'Tsat   = ', Tsat

! Write sound velocity

!OPEN (UNIT = 22,             FILE = 'sound.txt', &
!&   FORM = 'formatted',     ACTION = 'write',   &
!&   STATUS = 'replace')
!
!DO i = 1, N
!     WRITE(22,*) alpha(i), qual(i), cHEM(i)
!ENDDO
 
!CLOSE(UNIT = 22)

END PROGRAM sound
