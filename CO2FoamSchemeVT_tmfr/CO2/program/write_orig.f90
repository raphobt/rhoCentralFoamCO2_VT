PROGRAM write_solut
!      
      USE def_constants
      USE def_variables
      USE non_linear_solvers
      USE properties
      USE grid_functions
      USE Grid
      USE Interp_table
      IMPLICIT NONE


      INTEGER :: i,j,flag
      REAL(pr) :: cp
      REAL(pr) :: p_out, T_out, c_out, x_out, vp_out, xxx
      
      REAL(pr), DIMENSION(NNN_sat_LH) ::  PsatLH, TsatLH
      REAL(pr), DIMENSION(NNN_sat_LL) ::  PsatLL, TsatLL
      REAL(pr), DIMENSION(NNN_sat_R) ::  PsatR, TsatR
      

        CALL MAKE_GRID()
!---------------------------------------------------------------------------
!
!                     domain
!
!---------------------------------------------------------------------------
! Save two boundarys
!OPEN (UNIT = 22,             FILE = 'boundary_LH.txt', &
!&   FORM = 'formatted',     ACTION = 'write',   &
!&   STATUS = 'replace')
!
!DO i = 1, NNN_LH
!     WRITE(22,*), vvv_LH(i,1),vvv_LH(i,MMM_LH), y_mesh_LH(i)
!ENDDO

! CLOSE(UNIT = 22)

!
!
!
!!!!!!!!!!!!!!!! USE LOOK-UP TABLE TO OUTPUT: (Psat Tsat) = f (esat, vsat) !!!!!!!!!!!!!!!!!

p_out = 1.d0
T_out = 1.d0
c_out = 1.d0
x_out = 1.d0
vp_out  = 1.d0
xxx = 1.d0
flag  = 1

DO i = 1, NNN_sat_LL-1 ! avoid critical point
   CALL CO2BLLT_EQUI(p_out, T_out, c_out, x_out, vp_out, xxx, y_mesh_sat_LL(i), v_Lsat_LL(i), flag) 
   PsatLL(i) = p_out 
   TsatLL(i) = T_out
ENDDO

DO i = 2, NNN_sat_LH ! avoid critical point
   CALL CO2BLLT_EQUI(p_out, T_out, c_out, x_out, vp_out, xxx, y_mesh_sat_LH(i), v_Lsat_LH(i), flag) 
   PsatLH(i) = p_out 
   TsatLH(i) = T_out
ENDDO

DO i = 1, NNN_sat_R 
   CALL CO2BLLT_EQUI(p_out, T_out, c_out, x_out, vp_out, xxx, y_mesh_sat_R(i), v_Vsat(i), flag) 
   PsatR(i) = p_out 
   TsatR(i) = T_out
ENDDO

! Save spline
 
OPEN (UNIT = 22,             FILE = 'spline_LL.txt', &
&   FORM = 'formatted',     ACTION = 'write',   &
&   STATUS = 'replace')
!
DO i = 1, NNN_sat_LL-1 ! avoid critical point
      WRITE(22,*) v_Lsat_LL(i), v_Lpmax_LL(i),  y_mesh_sat_LL(i), v_liq_meta(i), PsatLL(i), TsatLL(i) 
ENDDO
!
 CLOSE(UNIT = 22)

OPEN (UNIT = 22,             FILE = 'spline_LH.txt', &
&   FORM = 'formatted',     ACTION = 'write',   &
&   STATUS = 'replace')

DO i = 2, NNN_sat_LH ! avoid critical point
      WRITE(22,*) v_Lsat_LH(i), v_Lpmax_LH(i),  y_mesh_sat_LH(i),  PsatLH(i), &
& TsatLH(i)
ENDDO

 CLOSE(UNIT = 22)

OPEN (UNIT = 22,             FILE = 'spline_HT.txt', &
&   FORM = 'formatted',     ACTION = 'write',   &
&   STATUS = 'replace')
!
DO i = 1, NNN_sat_HT
      WRITE(22,*) v_right_HT(i), v_left_HT(i),  y_mesh_sat_HT(i)
ENDDO

 CLOSE(UNIT = 22)

OPEN (UNIT = 22,             FILE = 'spline_R.txt', &
&   FORM = 'formatted',     ACTION = 'write',   &
&   STATUS = 'replace')
!
DO i = 1, NNN_sat_R
     WRITE(22,*) v_Vsat(i), v_Vpmin(i),  y_mesh_sat_R(i),  PsatR(i),  TsatR(i)
ENDDO

 CLOSE(UNIT = 22)
!STOP

! Save e, v, p, t, c
    OPEN (UNIT = 42,             FILE   = 'LH_p_T_c.txt', &
     &   FORM    = 'formatted',    ACTION = 'write',   &
     &   STATUS  = 'replace')
      DO i = 1, NNN_LH
         DO j = 1, MMM_LH
            CALL heat_cap_p(TTT_LH(i,j),vvv_LH(i,j),cp)
            WRITE(42,*) y_mesh_LH(i), vvv_LH(i,j), ppp_LH(i,j), &
     &                                          TTT_LH(i,j), ccc_LH(i,j),cp
         ENDDO
      ENDDO

      CLOSE(UNIT = 42)

    OPEN (UNIT = 42,             FILE   = 'LL_p_T_c.txt', &
     &   FORM    = 'formatted',    ACTION = 'write',   &
     &   STATUS  = 'replace')
      DO i = 1, NNN_LL
         DO j = 1, MMM_LL
            CALL heat_cap_p(TTT_LL(i,j),vvv_LL(i,j),cp)
            WRITE(42,*) y_mesh_LL(i), vvv_LL(i,j), ppp_LL(i,j), &
     &                                          TTT_LL(i,j), ccc_LL(i,j),cp
         ENDDO
      ENDDO

      CLOSE(UNIT = 42)
!
   OPEN (UNIT = 42,             FILE   = 'HT_p_T_c.txt', &
     &   FORM    = 'formatted',    ACTION = 'write',   &
     &   STATUS  = 'replace')
      DO i = 1, NNN_HT
         DO j = 1, MMM_HT
            CALL heat_cap_p(TTT_HT(i,j),vvv_HT(i,j),cp)
            WRITE(42,*) y_mesh_HT(i), vvv_HT(i,j), ppp_HT(i,j), &
     &                                          TTT_HT(i,j), ccc_HT(i,j),cp
         ENDDO
      ENDDO

      CLOSE(UNIT = 42)
!
   OPEN (UNIT = 42,             FILE   = 'R_p_T_c.txt', &
     &   FORM    = 'formatted',    ACTION = 'write',   &
     &   STATUS  = 'replace')
      DO i = 1, NNN_R
         DO j = 1, MMM_R
            CALL heat_cap_p(TTT_R(i,j),vvv_R(i,j),cp)
            WRITE(42,*) y_mesh_R(i), vvv_R(i,j), ppp_R(i,j), &
     &                                          TTT_R(i,j), ccc_R(i,j),cp
         ENDDO
      ENDDO

      CLOSE(UNIT = 42)





END PROGRAM write_solut
