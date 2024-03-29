#include <cmath>
//************************** Output T,p,mass flux at inlet cout and hout *****************************//

// Patch identification
label inlet1PatchID = mesh.boundaryMesh().findPatchID("inlet1");
label inlet2PatchID = mesh.boundaryMesh().findPatchID("inlet2");
label inlet3PatchID = mesh.boundaryMesh().findPatchID("inlet3");
label inlet4PatchID = mesh.boundaryMesh().findPatchID("inlet4");
label inlet5PatchID = mesh.boundaryMesh().findPatchID("inlet5");
label inlet6PatchID = mesh.boundaryMesh().findPatchID("inlet6");
label coldOutletPatchID = mesh.boundaryMesh().findPatchID("coldOutlet");
label hotOutletPatchID = mesh.boundaryMesh().findPatchID("hotOutlet");

// Initialization
scalar mfri1    = fabs(gSum(rho.boundaryField()[inlet1PatchID]*mag(U.boundaryField()[inlet1PatchID])*mesh.magSf().boundaryField()[inlet1PatchID]));
scalar mfri2    = fabs(gSum(rho.boundaryField()[inlet2PatchID]*mag(U.boundaryField()[inlet2PatchID])*mesh.magSf().boundaryField()[inlet2PatchID]));
scalar mfri3    = fabs(gSum(rho.boundaryField()[inlet3PatchID]*mag(U.boundaryField()[inlet3PatchID])*mesh.magSf().boundaryField()[inlet3PatchID]));
scalar mfri4    = fabs(gSum(rho.boundaryField()[inlet4PatchID]*mag(U.boundaryField()[inlet4PatchID])*mesh.magSf().boundaryField()[inlet4PatchID]));
scalar mfri5    = fabs(gSum(rho.boundaryField()[inlet5PatchID]*mag(U.boundaryField()[inlet5PatchID])*mesh.magSf().boundaryField()[inlet5PatchID]));
scalar mfri6    = fabs(gSum(rho.boundaryField()[inlet6PatchID]*mag(U.boundaryField()[inlet6PatchID])*mesh.magSf().boundaryField()[inlet6PatchID]));
scalar mfrc     = gSum(rho.boundaryField()[coldOutletPatchID]*mag(U.boundaryField()[coldOutletPatchID])*mesh.magSf().boundaryField()[coldOutletPatchID]);
scalar mfrh     = gSum(rho.boundaryField()[hotOutletPatchID]*mag(U.boundaryField()[hotOutletPatchID])*mesh.magSf().boundaryField()[hotOutletPatchID]);
scalar muc      = mfrc/(mfri1+mfri2+mfri3+mfri4+mfri5+mfri6+1e-8);
scalar outsurin = (mfrc+mfrh)/(mfri1+mfri2+mfri3+mfri4+mfri5+mfri6+1e-8);

scalar p_in = 0.0, p_h = 0.0, p_c = 0.0;
scalar T_in = 0.0, T_h = 0.0, T_c = 0.0;
scalar nin = 0, nh = 0, nc = 0;

const UList<label> &bfaceCellsIn   = mesh.boundaryMesh()[inlet1PatchID].faceCells(); // pin & Tin are only assessed via inlet1!!
const UList<label> &bfaceCellsHot  = mesh.boundaryMesh()[hotOutletPatchID].faceCells();
const UList<label> &bfaceCellsCold = mesh.boundaryMesh()[coldOutletPatchID].faceCells();

forAll ( bfaceCellsIn, ifaceIn ) 
{
     label bcCellIn = bfaceCellsIn[ifaceIn];
     p_in = p_in + p[bcCellIn];
     T_in = T_in + T[bcCellIn];
     nin = nin + 1;
}
forAll ( bfaceCellsCold, ifaceCold ) 
{
     label bcCellCold = bfaceCellsCold[ifaceCold];
     p_c = p_c + p[bcCellCold];
     T_c = T_c + T[bcCellCold];
     nc = nc + 1;
}
forAll ( bfaceCellsHot, ifaceHot ) 
{
     label bcCellHot = bfaceCellsHot[ifaceHot];
     p_h = p_h + p[bcCellHot];
     T_h = T_h + T[bcCellHot];
     nh = nh + 1;
}

// Output boundary variables
Info << " " << endl;
Info << "inlet pressure "<< nin << " " << p_in/(1e-8+nin) << " Pa" << endl;
Info << " " << endl;
Info << "cold outlet pressure "<< nc << " " << p_c/(1e-8+nc) << " Pa" << endl;
Info << "hot outlet pressure "<< nh << " " << p_h/(1e-8+nh) << " Pa" << endl;
Info << " " << endl;
Info << "inlet temperature "<< nin << " " << T_in/(1e-8+nin) << " K" << endl;
Info << "cold outlet temperature "<< nc << " " << T_c/(1e-8+nc) << " K" << endl;
Info << "hot outlet temperature "<< nh << " " << T_h/(1e-8+nh) << " K" << endl;
Info << " " << endl;
Info << "inlet mass flow rate "<< mfri1+mfri2+mfri3+mfri4+mfri5+mfri6 << " kg/s (360deg)" << endl;
Info << "cold outlet mass flow rate "<< mfrc << " kg/s (360deg)" << endl;
Info << "hot outlet mass flow rate "<< mfrh << " kg/s (360deg)" << endl;
Info << " " << endl;
Info << "out/in mass fraction "<< outsurin << endl;
Info << "cold mass fraction "<< muc << endl;
Info << " " << endl;
