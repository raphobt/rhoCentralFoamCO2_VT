#include <cmath>
Info<< "\nUpdating inlet conditions" << endl;
rho.correctBoundaryConditions();
U.correctBoundaryConditions();
e.correctBoundaryConditions();
p.correctBoundaryConditions();
T.correctBoundaryConditions();
//
rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();
rhoE.boundaryFieldRef() == rho.boundaryField()*(
                           e.boundaryField() + 0.5*magSqr(U.boundaryField())
                           );
//
//
const volVectorField&     C = mesh.C();
scalar M_ii = 0.0;
scalar res2,res3,res4;
scalar dp_dv_u, dp_du_v;
//
label     adj1_x=0;//,adj2_x=0;
scalar dt = runTime.deltaTValue();
scalar p_averc=0, p_averh=0, numc=0, numh=0;//, rho_averh=0.0, rho_averc=0.0;
//scalar T_averc=0, T_averh=0; 
//
//
//
/*
// Patch identification
label inletPatchID = mesh.boundaryMesh().findPatchID("inlet");
label coldOutletPatchID = mesh.boundaryMesh().findPatchID("coldOutlet");
label hotOutletPatchID = mesh.boundaryMesh().findPatchID("hotOutlet");

// Initialization
scalar mfri=0.0, mfrc=0.0, mfrh=0.0, mfrh_tar=0.0, muc=0.0, areah=0.0, mfrc_tar=0.0, areac=0.0;
mfri     = gSum(rho.boundaryField()[inletPatchID]*mag(U.boundaryField()[inletPatchID])*mesh.magSf().boundaryField()[inletPatchID]);
mfrc     = gSum(rho.boundaryField()[coldOutletPatchID]*mag(U.boundaryField()[coldOutletPatchID])*mesh.magSf().boundaryField()[coldOutletPatchID]);
mfrh     = gSum(rho.boundaryField()[hotOutletPatchID]*mag(U.boundaryField()[hotOutletPatchID])*mesh.magSf().boundaryField()[hotOutletPatchID]);
mfrh_tar = (1-muc_tar)*mfri; //(1-0.64)*0.0542/6; // (1-muc_tar)*mfri; // muc_tar defined in OR_op/VT2.H
mfrc_tar = muc_tar*mfri;
muc      = mfrc/(mfri+1e-8);
areah    = gSum(mesh.magSf().boundaryField()[hotOutletPatchID]);
areac    = gSum(mesh.magSf().boundaryField()[coldOutletPatchID]);
*/
forAll (mesh.boundaryMesh(),ipatch) 
{
       word BCTypePhysical = mesh.boundaryMesh().physicalTypes()[ipatch];
       word BCType         = mesh.boundaryMesh().types()[ipatch];
       word BCName         = mesh.boundaryMesh().names()[ipatch];
       const UList<label> &bfaceCells = mesh.boundaryMesh()[ipatch].faceCells();
       
       //Info<< " " << BCTypePhysical <<" "<< BCType <<" "<< BCName << endl;

       if ( BCName == "hotOutlet" ) 
       {

        if ( BCTypePhysical == "NSCBC" )
        {
          // find max mach (defined in initFields.H)
          
          forAll ( bfaceCells , iface )
          {
                label  bcCell = bfaceCells[iface];
                M_ii = mag(U[bcCell])/c[bcCell];
                if (M_max_h < M_ii){ M_max_h = M_ii; }
          }
          
          forAll ( bfaceCells , iface ) 
          {
                 label  bcCell = bfaceCells[iface];
                 labelList adjI  = mesh.cellCells()[bcCell];
                 
                 // find adjacent cell
                 
              // scalar x_bc  = C[bcCell].x();
                 scalar xmin_adj=100.0;
                 
                 forAll (adjI,adjII)
                 {
                        scalar x_adj = C[adjI[adjII]].x();
                        if (x_adj<xmin_adj)
                        {
                           adj1_x   = adjI[adjII];
                           xmin_adj = x_adj;
                        }
                 }
                 //Info << adj1_x << endl;
                 
                 //compute gradient at boundary cells
                 
              // scalar deltax   =  C[bcCell].x()-C[adj1_x].x();
                 scalar drhodx   = (rho[bcCell] - rho[adj1_x]) / (C[bcCell].x()-C[adj1_x].x()); //gradrho[bcCell] & normal;
                 scalar dpdx     = (p[bcCell]   - p[adj1_x])   / (C[bcCell].x()-C[adj1_x].x()); //gradp[bcCell]   & normal;
                 scalar dudx     = (U[bcCell].x() - U[adj1_x].x()) / (C[bcCell].x()-C[adj1_x].x());
                 scalar dvdx     = (U[bcCell].y() - U[adj1_x].y()) / (C[bcCell].x()-C[adj1_x].x());
                 scalar dwdx     = (U[bcCell].z() - U[adj1_x].z()) / (C[bcCell].x()-C[adj1_x].x());
                 //Info << C[bcCell].x()-C[adj1_x].x() << endl;
                 //Info << rho[bcCell] - rho[adj1_x] << endl;
                 scalar p_i      = p[bcCell];
                 scalar rho_i    = rho[bcCell];
                 scalar v_i      = 1.0/rho_i;
                 scalar c_i      = c[bcCell];
                 scalar e_i      = e[bcCell] - e_ref;
                 scalar T_i      = T[bcCell];
                 scalar ux       = U[bcCell].x();
                 scalar uy       = U[bcCell].y();
                 scalar uz       = U[bcCell].z();
                 scalar rhoE_i   = rhoE[bcCell];
                 scalar M_i      = mag(U[bcCell])/c[bcCell];
                 p_averh         = p_i + p_averh;
               //  T_averh         = T_i + T_averh;
               //  rho_averh       = rho_i + rho_averh;
                 numh            = 1.0 + numh;
                 
              /* scalar cbc      = c.boundaryFieldRef()[ipatch][iface];
                 scalar mbc      = mag(U.boundaryFieldRef()[ipatch][iface])/cbc;
                 scalar p_i      = p[adj1_x];
                 scalar rho_i    = rho[adj1_x];
                 scalar c_i      = c[adj1_x];
                 scalar ux       = U[adj1_x].x();
                 scalar uy       = U[adj1_x].y();
                 scalar uz       = U[adj1_x].x();
                 scalar M_i      = mag(U[adj1_x])/c[adj1_x];
              */
                 if (M_max_h<1.0)
                 // if (M_i<1.0)
                 // if (mbc<1.0)
                 {  
                  // scalar K_i      = relax*0.25*c_i*(1.0-M_max_h*M_max_h)/lx;
                     scalar K_i      = relaxHot*0.25*c_i*(1.0-M_i*M_i)/lx;
                 
                     // characteristic waves
                 
                     scalar L1       = ux * (drhodx - dpdx/(c_i*c_i));
                     scalar L2       = ux * dvdx; 
                     scalar L3       = ux * dwdx;
                     scalar L4       = (ux+c_i) * ( dudx+dpdx/(rho_i*c_i) );
                     scalar L5       = K_i * (p_i - p_outleth);
                 
                     // characteristic variables
                 
                     scalar d1       =  L1 + rho_i*(L4+L5)/(2.0*c_i);
                     scalar d2       =  0.5*(L4-L5);
                     scalar d3       =  L2;
                     scalar d4       =  L3;
                     scalar d5       =  rho_i*c_i*0.5*(L4+L5);
                 
                     // compute derivatives
                 
                     //Info<<bcCell<<" ei= "<<e_i<<" vi= "<<v_i<<"T= "<<T_i<<"p= "<<p_i<<endl;
                 
                     __derivees_MOD_co2der(&dp_dv_u, &dp_du_v , &e_i, &v_i, &T_i, &p_i,&res2,&res3,&res4);
                 
                     //Info<<bcCell<<" ei= "<<e_i<<" vi= "<<v_i<<" dpdu= "<<dp_du_v<<"dpdv= "<<dp_dv_u<<endl;

                     scalar dvdr    =  -1.0 / (rho_i*rho_i);
                     scalar dedp_r  =   1.0 / dp_du_v;
                     scalar dedr_p  =  -dp_dv_u*dvdr / dp_du_v; 
                 
                     // conservative variables

                     scalar Qc1     = rho_i - dt*d1;
                     scalar Qc2     = rho_i*ux - dt*(ux*d1+rho_i*d2);                 
                     scalar Qc3     = rho_i*uy - dt*(uy*d1+rho_i*d3);                 
                     scalar Qc4     = rho_i*uz - dt*(uz*d1+rho_i*d4);                 
                     scalar Qc5     = rhoE_i - dt*(rhoE_i/rho_i*d1 + 
                                              rho_i*ux*d2 + rho_i*uy*d3 + rho_i*uz*d4 + 
                                              rho_i*(dedp_r*d5 + dedr_p*d1));

                     rho.boundaryFieldRef()[ipatch][iface]       = Qc1;
                     rhoU.boundaryFieldRef()[ipatch][iface].x()  = Qc2;
                     rhoU.boundaryFieldRef()[ipatch][iface].y()  = Qc3;
                     rhoU.boundaryFieldRef()[ipatch][iface].z()  = Qc4;
                     rhoE.boundaryFieldRef()[ipatch][iface]      = Qc5;

                     U.boundaryFieldRef()[ipatch][iface].x() = Qc2 / Qc1; 
                     U.boundaryFieldRef()[ipatch][iface].y() = Qc3 / Qc1; 
                     U.boundaryFieldRef()[ipatch][iface].z() = Qc4 / Qc1; 
                     e.boundaryFieldRef()[ipatch][iface]     = Qc5 / Qc1 - 0.5*magSqr(U.boundaryFieldRef()[ipatch][iface]); 
                     p.boundaryFieldRef()[ipatch][iface]     = p_i - dt * d5;
                 }
                 else
                 {
                     rho.boundaryFieldRef()[ipatch][iface]   = rho[bcCell];
                     U.boundaryFieldRef()[ipatch][iface]     = U[bcCell];
                     e.boundaryFieldRef()[ipatch][iface]     = e[bcCell];
                     p.boundaryFieldRef()[ipatch][iface]     = p[bcCell];
                     T.boundaryFieldRef()[ipatch][iface]     = T[bcCell];

                     rhoU.boundaryFieldRef()[ipatch][iface]  = rho[bcCell] * U[bcCell];
                     rhoE.boundaryFieldRef()[ipatch][iface]  = rho[bcCell] * (e[bcCell]+0.5*magSqr(U[bcCell]));               
                  /* rho[bcCell]   =  rho[adj1_x];
                     U[bcCell]     =  U[adj1_x];
                     e[bcCell]     =  e[adj1_x];
                     Info << "supersonic outlet "<<"Mach= " << M_i << endl; 
                  */
                 } //endif sub/supersonic
           
          } //end hot outlet NSCBC cells loop
          
   } //end hot outlet NSCBC type treatment

   
    /* rhoU.boundaryFieldRef()[ipatch] = rho.boundaryField()[ipatch]*U.boundaryField()[ipatch];
       rhoE.boundaryFieldRef()[ipatch] = rho.boundaryField()[ipatch]*
                                         ( e.boundaryField()[ipatch]+ 
                                           0.5*magSqr(U.boundaryField()[ipatch]) );
    */   
        else if (BCTypePhysical == "zeroExtrapolation")
        {
           Info << " outlet zero extrapolation not complete "<< endl;
           
            forAll ( bfaceCells , iface ) 
            {
                 // Info<< bfaceCells << " iface= "<< iface << endl; 
                             
                 label  bcCell = bfaceCells[iface];
                 labelList adjI  = mesh.cellCells()[bcCell];
                 
                 // find adjacent cell
                 
                 scalar xmin_adj=100.0;
                 forAll (adjI,adjII)
                 {
                        scalar x_adj = C[adjI[adjII]].x();
                        if (x_adj<xmin_adj)
                        {
                           adj1_x   = adjI[adjII];
                           xmin_adj = x_adj;
                        }
                 }
                 rho[bcCell] = rho[adj1_x];
                 U[bcCell]   = U[adj1_x];
                 p[bcCell]   = p[adj1_x];
                 e[bcCell]   = e[adj1_x];
                 T[bcCell]   = T[adj1_x];

             } //end hot outlet cells loop
          } //end hot outlet zeroExtrapolation type treatment
       
/*
if (fabs(1-fabs(muc/(1e-8+muc_tar))) > 0.05)
{
   Info<< "\n---> Adjusting hot pressure to match target cold mass fraction" << endl;
      
      if (p_outleth >= p_averh/(1e-8+numh))
      {
          Info<< "\n---> Hot pressure is lower than target" << endl;
          p_outleth = p_averh/(1e-8+numh) + 0.5 * rho_averh/(1e-8+numh) * (mfrh*mfrh - mfrh_tar*mfrh_tar) / (1e-8+rho_averh*areah*rho_averh*areah/(1e-8+numh)/(1e-8+numh));
      }
      else 
      {
          Info<< "\n---> Hot pressure is higher than target" << endl;
          p_outleth = p_averh/(1e-8+numh) - 0.5 * rho_averh/(1e-8+numh) * (mfrh*mfrh - mfrh_tar*mfrh_tar) / (1e-8+rho_averh*areah*rho_averh*areah/(1e-8+numh)/(1e-8+numh));
      }
}
else
{
   Info<< "\nTarget muc is matched within 5 %" << endl;
}
*/

} //end hot outlet BC update
       
       if ( BCName == "coldOutlet" ) 
       {

        if ( BCTypePhysical == "NSCBC" )
        {
          // find max mach (defined in initFields.H)
          
          forAll ( bfaceCells , iface )
          {
                label  bcCell = bfaceCells[iface];
                M_ii = mag(U[bcCell])/c[bcCell];
                if (M_max_c < M_ii){ M_max_c = M_ii; }
          }
          
          forAll ( bfaceCells , iface ) 
          {
                 label  bcCell = bfaceCells[iface];
                 labelList adjI  = mesh.cellCells()[bcCell];
                 
                 // find adjacent cell
                 
              // scalar x_bc  = C[bcCell].x();
                 scalar xmin_adj=-100.0;
                 
                 forAll (adjI,adjII)
                 {
                        scalar x_adj = C[adjI[adjII]].x();
                        if (x_adj>xmin_adj)
                        {
                           adj1_x   = adjI[adjII];
                           xmin_adj = x_adj;
                        }
                 }
                 //Info << C[adj1_x].x() << endl;
                 
                 //compute gradient at boundary cells
                 
              // scalar deltax   =  C[bcCell].x()-C[adj1_x].x();
                 scalar drhodx   = (rho[bcCell] - rho[adj1_x]) / (C[bcCell].x()-C[adj1_x].x()); //gradrho[bcCell] & normal;
                 scalar dpdx     = (p[bcCell]   - p[adj1_x])   / (C[bcCell].x()-C[adj1_x].x()); //gradp[bcCell]   & normal;
                 scalar dudx     = (U[bcCell].x() - U[adj1_x].x()) / (C[bcCell].x()-C[adj1_x].x());
                 scalar dvdx     = (U[bcCell].y() - U[adj1_x].y()) / (C[bcCell].x()-C[adj1_x].x());
                 scalar dwdx     = (U[bcCell].z() - U[adj1_x].z()) / (C[bcCell].x()-C[adj1_x].x());
                 //Info << C[bcCell].x()-C[adj1_x].x() << endl;
                 //Info << rho[bcCell] - rho[adj1_x] << endl;
                 scalar p_i      = p[bcCell];
                 scalar rho_i    = rho[bcCell];
                 scalar v_i      = 1.0/rho_i;
                 scalar c_i      = c[bcCell];
                 scalar e_i      = e[bcCell] - e_ref;
                 scalar T_i      = T[bcCell];
                 scalar ux       = U[bcCell].x();
                 scalar uy       = U[bcCell].y();
                 scalar uz       = U[bcCell].z();
                 scalar rhoE_i   = rhoE[bcCell];
                 scalar M_i      = mag(U[bcCell])/c[bcCell];
                 p_averc         = p_i + p_averc;
               //  T_averc         = T_i + T_averc;
               //  rho_averc       = rho_i + rho_averc;
                 numc            = 1.0 + numc;
                 
              /* scalar cbc      = c.boundaryFieldRef()[ipatch][iface];
                 scalar mbc      = mag(U.boundaryFieldRef()[ipatch][iface])/cbc;
                 scalar p_i      = p[adj1_x];
                 scalar rho_i    = rho[adj1_x];
                 scalar c_i      = c[adj1_x];
                 scalar ux       = U[adj1_x].x();
                 scalar uy       = U[adj1_x].y();
                 scalar uz       = U[adj1_x].x();
                 scalar M_i      = mag(U[adj1_x])/c[adj1_x];
              */
                 if (M_max_c<1.0)
                 // if (M_i<1.0)
                 // if (mbc<1.0)
                 {  
                  // scalar K_i      = relax*0.25*c_i*(1.0-M_max_c*M_max_c)/lx;
                     scalar K_i      = relaxCold*0.25*c_i*(1.0-M_i*M_i)/lx;
                 
                     // characteristic waves
                 
                     scalar L1       = ux * (drhodx - dpdx/(c_i*c_i));
                     scalar L2       = ux * dvdx; 
                     scalar L3       = ux * dwdx;
                     scalar L4       = (ux+c_i) * ( dudx+dpdx/(rho_i*c_i) );
                     scalar L5       = K_i * (p_i - p_outletc);
                 
                     // characteristic variables
                 
                     scalar d1       =  L1 + rho_i*(L4+L5)/(2.0*c_i);
                     scalar d2       =  0.5*(L4-L5);
                     scalar d3       =  L2;
                     scalar d4       =  L3;
                     scalar d5       =  rho_i*c_i*0.5*(L4+L5);
                 
                     // compute derivatives
                 
                     //Info<<bcCell<<" ei= "<<e_i<<" vi= "<<v_i<<"T= "<<T_i<<"p= "<<p_i<<endl;
                 
                     __derivees_MOD_co2der(&dp_dv_u, &dp_du_v , &e_i, &v_i, &T_i, &p_i,&res2,&res3,&res4);
                 
                     //Info<<bcCell<<" ei= "<<e_i<<" vi= "<<v_i<<" dpdu= "<<dp_du_v<<"dpdv= "<<dp_dv_u<<endl;

                     scalar dvdr    =  -1.0 / (rho_i*rho_i);
                     scalar dedp_r  =   1.0 / dp_du_v;
                     scalar dedr_p  =  -dp_dv_u*dvdr / dp_du_v; 
                 
                     // conservative variables

                     scalar Qc1     = rho_i - dt*d1;
                     scalar Qc2     = rho_i*ux - dt*(ux*d1+rho_i*d2);                 
                     scalar Qc3     = rho_i*uy - dt*(uy*d1+rho_i*d3);                 
                     scalar Qc4     = rho_i*uz - dt*(uz*d1+rho_i*d4);                 
                     scalar Qc5     = rhoE_i - dt*(rhoE_i/rho_i*d1 + 
                                              rho_i*ux*d2 + rho_i*uy*d3 + rho_i*uz*d4 + 
                                              rho_i*(dedp_r*d5 + dedr_p*d1));

                     rho.boundaryFieldRef()[ipatch][iface]       = Qc1;
                     rhoU.boundaryFieldRef()[ipatch][iface].x()  = Qc2;
                     rhoU.boundaryFieldRef()[ipatch][iface].y()  = Qc3;
                     rhoU.boundaryFieldRef()[ipatch][iface].z()  = Qc4;
                     rhoE.boundaryFieldRef()[ipatch][iface]      = Qc5;

                     U.boundaryFieldRef()[ipatch][iface].x() = Qc2 / Qc1; 
                     U.boundaryFieldRef()[ipatch][iface].y() = Qc3 / Qc1; 
                     U.boundaryFieldRef()[ipatch][iface].z() = Qc4 / Qc1; 
                     e.boundaryFieldRef()[ipatch][iface]     = Qc5 / Qc1 - 0.5*magSqr(U.boundaryFieldRef()[ipatch][iface]); 
                     p.boundaryFieldRef()[ipatch][iface]     = p_i - dt * d5;
                 }
                 else
                 {
                     rho.boundaryFieldRef()[ipatch][iface]   = rho[bcCell];
                     U.boundaryFieldRef()[ipatch][iface]     = U[bcCell];
                     e.boundaryFieldRef()[ipatch][iface]     = e[bcCell];
                     p.boundaryFieldRef()[ipatch][iface]     = p[bcCell];
                     T.boundaryFieldRef()[ipatch][iface]     = T[bcCell];

                     rhoU.boundaryFieldRef()[ipatch][iface]  = rho[bcCell] * U[bcCell];
                     rhoE.boundaryFieldRef()[ipatch][iface]  = rho[bcCell] * (e[bcCell]+0.5*magSqr(U[bcCell]));               
                  /* rho[bcCell]   =  rho[adj1_x];
                     U[bcCell]     =  U[adj1_x];
                     e[bcCell]     =  e[adj1_x];
                     Info << "supersonic outlet "<<"Mach= " << M_i << endl; 
                  */
                 } //endif sub/supersonic
                 
       } //end cold outlet cells loop
                    
   } //end cold outlet NSCBC type treatment
            
    /* rhoU.boundaryFieldRef()[ipatch] = rho.boundaryField()[ipatch]*U.boundaryField()[ipatch];
       rhoE.boundaryFieldRef()[ipatch] = rho.boundaryField()[ipatch]*
                                         ( e.boundaryField()[ipatch]+ 
                                           0.5*magSqr(U.boundaryField()[ipatch]) );
    */   
        else if (BCTypePhysical == "zeroExtrapolation")
        {
           Info << " outlet zero extrapolation not complete "<< endl;
           
            forAll ( bfaceCells , iface ) 
            {
                 // Info<< bfaceCells << " iface= "<< iface << endl; 
                             
                 label  bcCell = bfaceCells[iface];
                 labelList adjI  = mesh.cellCells()[bcCell];
                 
                 // find adjacent cell
                 
                 scalar xmin_adj=-100.0;
                 forAll (adjI,adjII)
                 {
                        scalar x_adj = C[adjI[adjII]].x();
                        if (x_adj<xmin_adj)
                        {
                           adj1_x   = adjI[adjII];
                           xmin_adj = x_adj;
                        }
                 }
                 rho[bcCell] = rho[adj1_x];
                 U[bcCell]   = U[adj1_x];
                 p[bcCell]   = p[adj1_x];
                 e[bcCell]   = e[adj1_x];
                 T[bcCell]   = T[adj1_x];

             } //end cold outlet cells loop
          } //end cold outlet zeroExtrapolation type treatment

/*
if (fabs(1-fabs(muc/(1e-8+muc_tar))) > 0.05)
{
   Info<< "\n---> Adjusting cold pressure to match target cold mass fraction" << endl;
      
      if (p_outletc >= p_averc/(1e-8+numc))
      {
          Info<< "\n---> Cold pressure is lower than target" << endl;
          p_outletc = p_averc/(1e-8+numc) + 0.5 * rho_averc/(1e-8+numc) * (mfrc*mfrc - mfrc_tar*mfrc_tar) / (1e-8+rho_averc*areac*rho_averc*areac/(1e-8+numc)/(1e-8+numc));
      }
      else 
      {
          Info<< "\n---> Cold pressure is higher than target" << endl;
          p_outletc = p_averc/(1e-8+numc) - 0.5 * rho_averc/(1e-8+numc) * (mfrc*mfrc - mfrc_tar*mfrc_tar) / (1e-8+rho_averc*areac*rho_averc*areac/(1e-8+numc)/(1e-8+numc));
      }
}
else
{
   Info<< "\nTarget muc is matched within 5 %" << endl;
}
*/
       } //end cold outlet BC update
       
       if ( BCName == "inlet" )
       {
          if ( BCTypePhysical == "NSCBC" ) 
	  {
              forAll ( bfaceCells , iface ) 
	      {
                 //Info<< bfaceCells << " iface= "<< iface << endl;             
                 
                 label bcCell = bfaceCells[iface];

                 rho.boundaryFieldRef()[ipatch][iface]   = rho_inlet;
                 e.boundaryFieldRef()[ipatch][iface]     = e_inlet + e_ref;
                 U.boundaryFieldRef()[ipatch][iface]     = U[bcCell-1];
                 T.boundaryFieldRef()[ipatch][iface]     = T_inlet;
                 p.boundaryFieldRef()[ipatch][iface]     = p_inlet;

                 rhoU.boundaryFieldRef()[ipatch][iface]  = rho_inlet * U[bcCell-1];
                 rhoE.boundaryFieldRef()[ipatch][iface]  = rho_inlet * ( e_inlet+e_ref+0.5*magSqr(U[bcCell-1]) );

              // Info << "rhoIN " << rho.boundaryFieldRef()[ipatch][iface] << endl;
              // Info << "eIN " << e.boundaryFieldRef()[ipatch][iface] << endl;
              // Info << "UIN " << U.boundaryFieldRef()[ipatch][iface] << endl;
              // Info << "TIN " << T.boundaryFieldRef()[ipatch][iface] << endl;
              // Info << "pIN " << p.boundaryFieldRef()[ipatch][iface] << endl;

              } //end inlet cells
              
           }      
       }

}
/*
// Target mass flow rate match

Info << "\nmfrh_tar 0.0032<-> "<< mfrh_tar << endl;
Info << "rho_averh 100<-> "<< rho_averh/(1e-8+numh) << endl;
Info << "rho_averc 100<-> "<< rho_averc/(1e-8+numc) << endl;
Info << "areah "<< areah << endl;
Info << "areac "<< areac << endl;

// Output boundary variables
Info << " " << endl;
Info << "cold outlet pressure "<< numc << " " << p_averc/(1e-8+numc) << " Pa" << " " << M_max_c << endl;
Info << "hot outlet pressure "<< numh << " " << p_averh/(1e-8+numh) << " Pa" << " " << M_max_h << endl;
Info << " " << endl;
Info << "cold outlet temperature "<< numc << " " << T_averc/(1e-8+numc) << " K" << endl;
Info << "hot outlet temperature "<< numh << " " << T_averh/(1e-8+numh) << " K" << endl;
Info << " " << endl;
Info << "inlet mass flow rate "<< mfri*6 << " kg/s (360deg)" << endl;
Info << "cold outlet mass flow rate "<< mfrc*6 << " kg/s (360deg)" << endl;
Info << "hot outlet mass flow rate "<< mfrh*6 << " kg/s (360deg)" << endl;
Info << " " << endl;
Info << "out/in mass fraction "<< (mfrc+mfrh)/(mfri+1e-8) << endl;
Info << "cold mass fraction "<< mfrc/(mfri+1e-8) << endl;
Info << " " << endl;
*/
