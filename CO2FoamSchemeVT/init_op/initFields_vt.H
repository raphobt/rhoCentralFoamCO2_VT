#include <cmath> 
Info<< "Creating CO2 table ...\n" << endl;
//    //    Info<< "\n---------pp=" << c << endl;
    __grid_MOD_make_grid();
//    grid_construction_left_low_();
//    grid_construction_left_high_();
//    grid_construction_high_temperature_();
//    __interp_table_MOD_co2bllt_equi(&pp,&TT,&cc,&x_out,&a_out,&dummy,
//                                    &ee,&vv,&flag);
//    Info<< "\n pp,TT,cc" << pp << " " << TT << " " << cc << endl;
//    Info<< "\n xx,a,flag" << x_out<< " " << a_out << " " << flag << endl;
//*************************************************************************
//Info<< "\nInitialize conservative fields\n" << endl;
//
double ee=0.0; //-99772;
double vv=0.0; //1.0/128.4;
double pp=0.0;
double TT=0.0;
double cc=0.0;
double x_out;
double a_out;
double cp_out;
double cv_out;
double vp_out;
double dedr_out;
double lambda_out;
double dummy;
int    flag;
double vis;
const double e_ref=5.0e5; //avoir energie interne positive
scalar initial_time;
scalar M_max = 0.0;
initial_time = runTime.time().value();
//Info <<" time= " <<runTime.timeName()<<endl;
//Info <<" time= " <<initial_time<<endl;
//
if (initial_time < 0.00000001) 
{
//Info << "\n Smooth initial field \n" << endl;
Info << "\n Fluent initial field \n" << endl;
   const volVectorField&     center = mesh.C();
   scalar xcell;

/*   scalar e_in   = e_inlet + e_ref;
     scalar e_outh = e_outleth + e_ref;
     scalar e_outc = e_outletc + e_ref; */
 
/*     forAll(e,celli) 
   {
       e[celli] = e[celli] + e_ref;
       Info << e[celli] << endl;
   }
*/  
   
/*           xcell = center[celli].x();

       if ( (xcell > xColdIn) and (xcell < xHotIn) )
       {
              rho[celli]   = rho_inlet;
              e[celli]     = e_in;
             
              ee = e[celli] - e_ref;
              vv = 1.0/rho[celli];
              __interp_table_MOD_co2bllt_equi(&pp,&TT,&cc,&x_out,&a_out,&dummy,
                                           &ee,&vv,&flag);
              p[celli]  = pp;
              T[celli]  = TT;
              rhoU[celli] = rho[celli]*U[celli];
              rhoE[celli] = rho[celli]*(e[celli] + 0.5*magSqr(U[celli]));
       }
       else if ( xcell < xColdIn )
       {
              rho[celli] = rho_outletc - xcell/xColdIn * (rho_outletc - rho_inlet);
              e[celli]   = e_outc - xcell/xColdIn * (e_outc - e_in);
             
              ee = e[celli] - e_ref;
              vv = 1.0/rho[celli];
              __interp_table_MOD_co2bllt_equi(&pp,&TT,&cc,&x_out,&a_out,&dummy,
                                           &ee,&vv,&flag);
              p[celli]  = pp;
              T[celli]  = TT;
              rhoU[celli] = rho[celli]*U[celli];
              rhoE[celli] = rho[celli]*(e[celli] + 0.5*magSqr(U[celli]));

       }
       else if ( xcell > xHotIn )
       {
              rho[celli] = rho_inlet + (xcell-xHotIn)/(xHotOut-xHotIn) * (rho_outleth - rho_inlet);
              e[celli]   = e_in + (xcell-xHotIn)/(xHotOut-xHotIn) * (e_outh - e_in);
            
              ee = e[celli] - e_ref;
              vv = 1.0/rho[celli];
              __interp_table_MOD_co2bllt_equi(&pp,&TT,&cc,&x_out,&a_out,&dummy,
                                           &ee,&vv,&flag);
              p[celli]  = pp;
              T[celli]  = TT;
              rhoU[celli] = rho[celli]*U[celli];
              rhoE[celli] = rho[celli]*(e[celli] + 0.5*magSqr(U[celli]));

       }
   }
*/

}

//
//forAll(mesh.boundary(), patchID)
//        {
//         forAll (mesh.boundary()[patchID],facei)
//            {
//               TT = T.boundaryField()[patchID][facei];
//               vv = 1.0/rho.boundaryField()[patchID][facei];
//
//               __properties_MOD_inter_energy(&TT,&vv, &ee);
//
//               e.boundaryFieldRef()[patchID][facei] = ee + e_ref;
//            }
//        }
//
//rhoU = rho * U;
//rhoE = rho*(e + 0.5*magSqr(U));
//
//Info<< "\n e.boundary\n" << e.boundaryFieldRef() << endl;
//Info<< "\n rho\n" << rho << endl;
//Info<< "\n U\n"   << U   << endl;
//Info<< "\n e\n"   << e   << endl; 
