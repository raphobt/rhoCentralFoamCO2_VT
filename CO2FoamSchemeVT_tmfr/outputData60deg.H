//************************** Output T,p,mass flux at inlet cout and hout *****************************//

// Patch identification
label inletPatchID = mesh.boundaryMesh().findPatchID("inlet");
label coldOutletPatchID = mesh.boundaryMesh().findPatchID("coldOutlet");
label hotOutletPatchID = mesh.boundaryMesh().findPatchID("hotOutlet");

// Initialization
scalar mfri     = gSum(rho.boundaryField()[inletPatchID]*mag(U.boundaryField()[inletPatchID])*mesh.magSf().boundaryField()[inletPatchID]);
scalar mfrc     = gSum(rho.boundaryField()[coldOutletPatchID]*mag(U.boundaryField()[coldOutletPatchID])*mesh.magSf().boundaryField()[coldOutletPatchID]);
scalar mfrh     = gSum(rho.boundaryField()[hotOutletPatchID]*mag(U.boundaryField()[hotOutletPatchID])*mesh.magSf().boundaryField()[hotOutletPatchID]);
scalar muc      = mfrc/(mfri+1e-8);
scalar outsurin = (mfrc+mfrh)/(mfri+1e-8);

scalar p_in = 0.0, p_h = 0.0, p_c = 0.0;
scalar T_in = 0.0, T_h = 0.0, T_c = 0.0;
scalar nin = 0, nh = 0, nc = 0;

const UList<label> &bfaceCellsIn   = mesh.boundaryMesh()[inletPatchID].faceCells();
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
Info << "inlet mass flow rate "<< mfri*6 << " kg/s (360deg)" << endl;
Info << "cold outlet mass flow rate "<< mfrc*6 << " kg/s (360deg)" << endl;
Info << "hot outlet mass flow rate "<< mfrh*6 << " kg/s (360deg)" << endl;
Info << " " << endl;
Info << "out/in mass fraction "<< outsurin << endl;
Info << "cold mass fraction "<< muc << endl;
Info << " " << endl;
