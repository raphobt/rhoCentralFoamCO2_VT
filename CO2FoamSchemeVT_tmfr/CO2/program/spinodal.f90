PROGRAM spinodal
     
     USE def_constants
     USE non_linear_solvers
     USE properties
     IMPLICIT NONE

! Without table

!      INTEGER::Niter,flag
!      REAL(8)::vliq,vvap,Tsat,res,p,in_2,vliq_guess,vvap_guess,T_guess

!press = 6.d6
!vliq_guess = 1.d0/740 
!vvap_guess = 1.d0/219
!T_guess = 296
!          CALL New_Rap3D(3,vliq,vvap,Tsat,&
!      &        res, Niter, flag, press, in_2,&
!      &        vliq_guess,vvap_guess, T_guess)

!print*, 'vliq   = ', 1/vliq
!print*, 'vvap   = ', 1/vvap
!print*, 'Tsat   = ', Tsat

INTEGER :: i, Niter, exitflag
REAL(8) :: in_2, resnorm, press, T_guess_sat, v_guess_sat
REAL(8) :: u_min, utest_max, delta
REAL(8), DIMENSION(NNN_sat_LL) :: T_liq_meta, v_liq_meta, y_mesh_sat_LL, T_sat_sat, v_Lsat_LL


!--------------------------------------------------------------------------
!
! Construct vector y_mesh_LL(physical domain), x_mesh_LL(tansformed
! domain)
!
!-------------------------------------------------------------------------
        utest_max = 1.0_pr*e_cr ! interal energy smaller than critical internal energy
        u_min = e_tri_L       
!
! array y_mesh_sat_LL
        delta             = (utest_max - u_min)/(NNN_sat_LL-1)
        y_mesh_sat_LL     = u_min + (/(i*delta, i=0,NNN_sat_LL-1)/)

! META: compute spinodal liquid between  -290.004 - -190.3 (kj/kg)



!e = -427.18 - -290.004  (kj/kg) p = 0
      T_guess_sat = 263.0
      v_guess_sat = 1.1e-3_pr
 DO i = 1, 172
    CALL New_Rap2D(4, T_liq_meta(i), v_liq_meta(i), &
              & resnorm, Niter, exitflag, y_mesh_sat_LL(i),in_2,&
              & T_guess_sat,v_guess_sat)
!
   IF (resnorm > 1e-5) THEN
       PRINT*,'LL_sat_spin-Convergence failed calculating T_sat_sat at spinodal(',i,')'
       PRINT*,'resnorm     = ',resnorm ,'N_iter = ',Niter
       STOP
   ENDIF
      T_guess_sat = T_liq_meta(i)
      v_guess_sat = v_liq_meta(i)+0.1e-4_pr
      CALL pressure(T_liq_meta(i),v_liq_meta(i), press)
   PRINT*, i, T_liq_meta(i), v_liq_meta(i), y_mesh_sat_LL(i), press

    T_liq_meta(i) = T_sat_sat(i)
    v_liq_meta(i) = v_Lsat_LL(i)
 ENDDO


!e = -290.004 - -201.476 (kj/kg) dp/dv_T = 0
      T_guess_sat = 275.0
      v_guess_sat = v_tri_L + 1e-4_pr
 DO i = 173, 284
!
        CALL New_Rap2D(3, T_liq_meta(i), v_liq_meta(i), &
              & resnorm, Niter, exitflag, y_mesh_sat_LL(i),in_2,&
              & T_guess_sat,v_guess_sat)
!
   IF (resnorm > 1e-5) THEN
       print*,'LL_sat_spin-Convergence failed calculating T_sat_sat at spinodal(',i,')'
       print*,'resnorm     = ',resnorm ,'N_iter = ',Niter
       STOP
   ENDIF
   IF (i<200) THEN
      T_guess_sat = T_liq_meta(i)
      v_guess_sat = v_liq_meta(i)+1e-4_pr
!   ELSEIF (i>=284) THEN
!      T_guess_sat = T_liq_meta(i)+0.1_pr
!      v_guess_sat = v_liq_meta(i)
   ELSE 
      T_guess_sat = T_liq_meta(i)
      v_guess_sat = v_liq_meta(i)-1e-4_pr
  
   ENDIF
      CALL pressure(T_liq_meta(i),v_liq_meta(i), press)
   PRINT*, i, T_liq_meta(i), v_liq_meta(i), y_mesh_sat_LL(i), press
 ENDDO


!e = -201.476 - -190.3 (kj/kg) saturation line
   DO i = 285, NNN_sat_LL
      T_liq_meta(i) = T_sat_sat(i)
      v_liq_meta(i) = v_Lsat_LL(i)
      PRINT*, i, T_liq_meta(i), v_liq_meta(i), y_mesh_sat_LL(i), press
 ENDDO
!   v_liq_meta(NNN_sat_LL)  = 1_pr/rho_cr

!print*, 'start p_max spline'
! P_max

PRINT*, ' '
PRINT*, 'To put in write.f90 to plot in (e,v) or other (python) along with sat curve!'
PRINT*, ' '

!OPEN (UNIT = 22,             FILE = 'spinodal.txt', &
!&   FORM = 'formatted',     ACTION = 'write',   &
!&   STATUS = 'replace')

!DO i = 1, N
!     WRITE(22,*) 
!ENDDO
 
!CLOSE(UNIT = 22)

END PROGRAM spinodal
