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
      INTEGER :: Niter, exitflag      
      REAL(pr):: T_in, v_in, p_in, e_in, h_in, cv, cp, s_in, c, v_l, v_v, v_l_Span, v_v_Span, &
     &           resnorm, guess_1, guess_2, T_span, v_span, e_span,helmho1
      REAL(pr):: p2, e2, cv2, cp2, s2, c2,helmho2,Helmholtz,eguess_in,out3,out2,press
      REAL(pr):: T1,v1,s_entro1,e1,p1,u1,press1,h1
      REAL(pr):: T2,v2,s_entro2,u2,press2,h2
!
! TEST 1 pressure, e, cv, cp, c, s, helmholtz
!     v = 0.00124937  ! corresponding to 1 MPa (Table 35, page 1568)
!      v = 1_pr/ 801.62_pr
     T_in = 294.15d0
     p_in = 1.883d6
!     eguess = -262219.0
     eguess_in = 1.0/30 !524.06 See SW article table to estimate e_guess <-> rho
!!!      CALL property (T,v,p,e,cv,cp,s,c)
!!!      helmho1 = Helmholtz(v,T)
 !    T1= 42+273.15
 !    v1= 1/582.22
 !    p1= 10e6 !9.0
 !    u1= 0 

 !    T2= 10.0+273.15
 !    v2= 1/161.8
 !    p2= 2.1e6 !6.1
 !    u2= 0
!
!!!        print*, 'test pressure -- T,v', T,v
!     CALL pressure(T,v,p2)
!      CALL entropy(T,v,s2)
!      CALL heat_cap_v(T,v,cv2)
!      CALL heat_cap_p(T,v,cp2)
!      CALL inter_energy(T,v,e2)
!      CALL sound_speed(T,v,c2)
!      CALL helmho(T,v,helmho2)

!       CALL MAKE_GRID()
!       CALL eos_1d(3, e, e2, resnorm, Niter,&
!     &             exitflag, p, eguess, v, p2) 

       CALL New_Rap1D(3, v_in, out2, resnorm, Niter,&
     &                             exitflag, p_in, eguess_in, T_in, out3)     
       CALL inter_energy(T_in,v_in,e2)
       CALL pressure(T_in,v_in,press)
     
      print*, 
      print*,'----------INLET OPERATING CONDITIONS-----------'
      print*, 
      print*, 'v_in   = ', v_in
      print*, 'rho_in = ', 1.0/v_in
      print*, '----------------------------------------------'
      print*, 'p_in         = ', press*1e-6, 'MPa (iterative)'
      print*, 'T_in         = ', T_in, 'K'
      print*, '----------------------------------------------'

       CALL entropy(T_in,v_in,s_in)
       CALL inter_energy(T_in,v_in,e_in)
       CALL pressure(T_in,v_in,p_in)
       h_in = e_in+p_in*v_in
     
      print*, 'e_in         = ', e2, 'J/kg'
      print*, '(e+e_ref)_in = ', e2+5.d5, 'J/kg'
      print*, 's_in         = ', s_in, 'J/kg'
      print*, 'h_in         = ', h_in, 'J/kg'
      print*, '----------------------------------------------'
      print*, 
 
! COLD OUTLET !
 
     T_in = 290.05d0
     p_in = 1.641d6
     eguess_in = 1.0/30 !524.06 See SW article table to estimate e_guess <-> rho

       CALL New_Rap1D(3, v_in, out2, resnorm, Niter,&
     &                             exitflag, p_in, eguess_in, T_in, out3)     
       CALL inter_energy(T_in,v_in,e2)
       CALL pressure(T_in,v_in,press)
     
      print*, 
      print*,'----------COLD OUTLET OPERATING CONDITIONS-----------'
      print*, 
      print*, 'v_cout   = ', v_in
      print*, 'rho_cout = ', 1.0/v_in
      print*, '----------------------------------------------'
      print*, 'p_cout         = ', press*1e-6, 'MPa (iterative)'
      print*, 'T_cout         = ', T_in, 'K'
      print*, '----------------------------------------------'

       CALL entropy(T_in,v_in,s_in)
       CALL inter_energy(T_in,v_in,e_in)
       CALL pressure(T_in,v_in,p_in)
       h_in = e_in+p_in*v_in
     
      print*, 'e_cout         = ', e2, 'J/kg'
      print*, '(e+e_ref)_cout = ', e2+5.d5, 'J/kg'
      print*, 's_cout         = ', s_in, 'J/kg'
      print*, 'h_cout         = ', h_in, 'J/kg'
      print*, '----------------------------------------------'
      print*, 
 
! HOT OUTLET !
 
     T_in = 293.35d0
     p_in = 1.629d6
     eguess_in = 1.0/30 !524.06 See SW article table to estimate e_guess <-> rho

       CALL New_Rap1D(3, v_in, out2, resnorm, Niter,&
     &                             exitflag, p_in, eguess_in, T_in, out3)     
       CALL inter_energy(T_in,v_in,e2)
       CALL pressure(T_in,v_in,press)
     
      print*, 
      print*,'----------HOT OUTLET OPERATING CONDITIONS-----------'
      print*, 
      print*, 'v_hout   = ', v_in
      print*, 'rho_hout = ', 1.0/v_in
      print*, '----------------------------------------------'
      print*, 'p_hout         = ', press*1e-6, 'MPa (iterative)'
      print*, 'T_hout         = ', T_in, 'K'
      print*, '----------------------------------------------'

       CALL entropy(T_in,v_in,s_in)
       CALL inter_energy(T_in,v_in,e_in)
       CALL pressure(T_in,v_in,p_in)
       h_in = e_in+p_in*v_in
     
      print*, 'e_hout         = ', e2, 'J/kg'
      print*, '(e+e_ref)_hout = ', e2+5.d5, 'J/kg'
      print*, 's_hout         = ', s_in, 'J/kg'
      print*, 'h_hout         = ', h_in, 'J/kg'
      print*, '----------------------------------------------'
      print*, 


 !      CALL entropy(T1,v1,s_entro1)
 !      CALL inter_energy(T1,v1,e1)
 !      CALL pressure(T1,v1,press1)
 !      h1 = e1+p1*v1

 !      CALL entropy(T2,v2,s_entro2)
 !      CALL inter_energy(T2,v2,e2)
 !      CALL pressure(T2,v2,press2)
 !      h2 = e2+p2*v2

 !     print*,'------------------Primary inlet-------------------------------------'
 !     print*, 's1= ', s_entro1,' h1= ', h1, "p1= ", press1, " ",p1
 !     print*,'------------------Secondary inlet----------------------------------------'
 !     print*, 's2= ', s_entro2,' h2= ', h2, "p2= ", press2, " ",p2
 !     print*,'----------------------------------------------------------'
 !     print*, 'exergy primary inlet= ', (h1-h2)+0.5*(u1*u1-u2*u2)-T2*(s_entro1-s_entro2)
 !     print*,'----------------------------------------------------------'

!      print*, ' cv ', cv2*1e-3_pr, ' error ', cv2*1e-3_pr - 0.68217_pr 
!      print*,'----------------------------------------------------------'
!      print*, ' cp  ' , cp2*1e-3_pr,' error ', cp2*1e-3_pr - 0.92089_pr
!      print*,'----------------------------------------------------------'
!      print*, ' s  ' , s2*1e-3_pr,' error ', s2*1e-3_pr + 0.44964
!      print*,'----------------------------------------------------------'
!      print*, ' c  ' , c2,' error ', c2 - 262.43_pr
!      print*,'----------------------------------------------------------'
!!      print*, ' helmholtz  ' , helmho2,' error ', helmho1*1e3_pr-helmho2
!!      print*,'----------------------------------------------------------'
!             
! TEST 3 - MARCO   saturation curve
!
!     T = 300_pr                                  ! K
! from table 34, page 1561. 
! At 300 K the saturated specific volumes are:
!     v_l_Span = 1_pr / 679.24_pr 
!     v_v_Span = 1_pr / 268.58_pr
!     
!     guess_1 =  v_l_Span*0.8_pr 
!     guess_2 =  v_v_Span*1.2_pr 
!
!
!      CALL New_Rap2D(1, v_l, v_v, &
!     &           resnorm, Niter, exitflag, T, 0d0, guess_1, guess_2)
!     
!      CALL pressure(T,v_l,p)
!
!      
!      print*, '----------------TEST 3-----------------------'
!      print*, 'resnorm, Niter, exitflag', resnorm, Niter, exitflag
!      print*,'-------------------------------------------------------------------------------'
!      print*, 'GUESS:        rho_l = ', 1_pr/guess_1,  ' rho_v= ',1_pr/guess_2
!      print*, 'From article: rho_l = ', 1_pr/v_l_Span, ' rho_v= ',1_pr/v_v_Span
!      print*, 'Newton:       rho_l = ', 1_pr/v_l,      ' rho_v= ',1_pr/v_v
!      print*,'--------------------------------------------------------------------------------'
!      print*, 'pressure in MPa, test 3', p*1e-6_pr
!             
! TEST 4 - MARCO    calculation of p and e corresponding to assigned T and V
!
! from table 34, page 1567. 
! At 300 K, 1 MPa the internal energy and the specific volumes are:
!     e_Span  = -61.765e3_pr 
!     v_Span  = 1_pr / 18.579_pr
!     T_Span       = 300_pr
!     
!     guess_1 = T*0.8_pr 
!     guess_2 = v_Span*0.8_pr 
!
!
!      CALL New_Rap2D(2, T, v, &
!     &           resnorm, Niter, exitflag, 1e6_pr, e_Span, guess_1, guess_2)
!     
!      CALL pressure(T,v,p)
!
!      
!      print*, '----------------TEST 4-----------------------'
!      print*, 'resnorm, Niter, exitflag', resnorm, Niter, exitflag
!     print*,'-------------------------------------------------------------------------------'
!      print*, 'GUESS:        T = ', guess_1,  ' v= ', guess_2
!      print*, 'From article: T = ', T_Span,   ' v= ', v_Span
!      print*, 'Newton:       T = ', T,        ' v= ', v
!      print*,'--------------------------------------------------------------------------------'
!      print*, 'pressure in MPa, test 4', p*1e-6_pr
!
!
END PROGRAM test_pressure
