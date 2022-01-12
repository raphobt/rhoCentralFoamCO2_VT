PROGRAM sound
     
     USE non_linear_solvers
     USE properties
     USE Grid
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
INTEGER, PARAMETER :: N = 100
REAL(8) :: resnorm, out2, out3 
REAL(8) :: press, Tsat, vv, vl, ev, el, ratio
REAL(8) :: Tml, vml, vml_guess, y
!REAL(8), DIMENSION(N) :: dvV_dp_T, dvml_dp_T,&
!                         dp_dv_Tv, dp_dv_Tml
REAL(8), DIMENSION(N) :: ds_dr_v, ds_dr_ml,&
                         ds_dT_v, ds_dT_ml,&
                         dr_dp_v, dr_dp_ml,&
                         dT_dp_v, dT_dp_ml,&
                         dr_dp_T_x, dT_dp_r_x, ds_dr_T_x, ds_dT_r_x
REAL(8), DIMENSION(N) :: duL_dp, duV_dp, dvL_dp, dvV_dp, du_dp_x, dv_dp_x ! HEM
REAL(8), DIMENSION(N) :: alpha,qual,qualEq,qualDEM,cHEM,cHRM,cDEM
REAL(8), DIMENSION(N) :: rho

press = 5.d6

!!!!!!!!!!!!!!!!!!!!!-HEM-!!!!!!!!!!!!!!!!!!!!!
CALL MAKE_GRID()
CALL satprop(3, press, Tsat, vv, vl, ev, el)
ratio = (ev - el)/(vv - vl)  ! (J/kg)/(m3/kg)

!!!!!!!!!!!!!!!!!!!!!-HRM-!!!!!!!!!!!!!!!!!!!!!
Tml = Tsat + 2.d1
vml_guess = 1.d0/120.d0 ! Respect to Tml
CALL New_Rap1D(3, vml, out2, resnorm, Niter, exitflag, press, vml_guess, Tml, out3) ! We get vml
!CALL dpdv_T(Tsat, vv, dp_dv_Tv)
!CALL dpdv_T(Tml, vml, dp_dv_Tml)
!dvV_dp_T = 1.d0/dp_dv_Tv
!dvml_dp_T = 1.d0/dp_dv_Tml

!!!!!!!!!!!!!!!!!!!!-DEM-!!!!!!!!!!!!!!!!!!!!!!
y       = 5.d-1
qualEq  = (vv - vml) / (vv - vl)  ! Stable phase quality in DEM
qualDEM = qualEq * y           ! Saturated vapor quality in DEM



DO i = 1, N
     alpha(i)     = dfloat(i) / N 
     qual(i)      = ( alpha(i)/vv ) / ( alpha(i)/vv + (1-alpha(i))/vl )

     ! HEM
     CALL satderiv(3, press, duL_dp(i), duV_dp(i), dvL_dp(i), dvV_dp(i))

     rho(i)       = alpha(i)/(qual(i)+1.d-10)/vv
     du_dp_x(i)   = qual(i) * duV_dp(i) + (1.d0 - qual(i)) * duL_dp(i)
     dv_dp_x(i)   = qual(i) * dvV_dp(i) + (1.d0 - qual(i)) * dvL_dp(i)

     cHEM(i)      = SQRT((press + ratio)/(rho(i)*rho(i)*(du_dp_x(i) - ratio * dv_dp_x(i))))

     ! HRM
     CALL dsdr_T(Tsat, vv, ds_dr_v(i))
     CALL dsdr_T(Tml, vml, ds_dr_ml(i))
     CALL dsdT_r(Tsat, vv, ds_dT_v(i))
     CALL dsdT_r(Tml, vml, ds_dT_ml(i))
     
     CALL drdp_T(Tsat, vv, dr_dp_v(i)) 
     CALL drdp_T(Tml, vml, dr_dp_ml(i)) 
     CALL dTdp_r(Tsat, vv, dT_dp_v(i))
     CALL dTdp_r(Tml, vml, dT_dp_ml(i))

     dr_dp_T_x(i) = qual(i) * dr_dp_v(i) + (1.d0 - qual(i)) * dr_dp_ml(i)
     dT_dp_r_x(i) = qual(i) * dT_dp_v(i) + (1.d0 - qual(i)) * dT_dp_ml(i)
     ds_dr_T_x(i) = qual(i) * ds_dr_v(i) + (1.d0 - qual(i)) * ds_dr_ml(i)
     ds_dT_r_x(i) = qual(i) * ds_dT_v(i) + (1.d0 - qual(i)) * ds_dT_ml(i)
     
     cHRM(i) = SQRT( 1.d0/dr_dp_T_x(i) - 1.d0/dT_dp_r_x(i) * ds_dr_T_x(i) / ds_dT_r_x(i) ) 
     
     ! DEM

     PRINT*, i, alpha(i), qual(i), cHEM(i), cHRM(i)
ENDDO

PRINT*, 'rhol   = ', 1.d0/vl
PRINT*, 'rhoml  = ', 1.d0/vml
PRINT*, 'rhov   = ', 1.d0/vv
PRINT*, 'Tsat   = ', Tsat

! Write sound velocity

OPEN (UNIT = 22,             FILE = 'sound.txt', &
&   FORM = 'formatted',     ACTION = 'write',   &
&   STATUS = 'replace')

DO i = 1, N
     WRITE(22,*) alpha(i), qual(i), cHEM(i), cHRM(i) !, cDEM(i)
ENDDO
 
CLOSE(UNIT = 22)

END PROGRAM sound
