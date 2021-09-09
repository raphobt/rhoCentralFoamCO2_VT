label inletPatchID = mesh.boundaryMesh().findPatchID("inlet");
label coldOutletPatchID = mesh.boundaryMesh().findPatchID("coldOutlet");
label hotOutletPatchID = mesh.boundaryMesh().findPatchID("hotOutlet");
scalar mfri=0.0, mfrc=0.0, mfrh=0.0;

//forAll(mesh.boundary()[inletPatchID], faceI)
//{
//      mfri = gSum(rho.boundaryField()[inletPatchID]*mag(U.boundaryField()[inletPatchID][faceI])*mesh.magSf().boundaryField()[inletPatchID][faceI]);
      mfri = gSum(rho.boundaryField()[inletPatchID]*mag(U.boundaryField()[inletPatchID])*mesh.magSf().boundaryField()[inletPatchID]);
//}
//forAll(mesh.boundary()[coldOutletPatchID], faceI)
//{
//      mfrc = gSum(rho.boundaryField()[coldOutletPatchID]*mag(U.boundaryField()[coldOutletPatchID][faceI])*mesh.magSf().boundaryField()[coldOutletPatchID][faceI]);
      mfrc = gSum(rho.boundaryField()[coldOutletPatchID]*mag(U.boundaryField()[coldOutletPatchID])*mesh.magSf().boundaryField()[coldOutletPatchID]);
//}
//forAll(mesh.boundary()[hotOutletPatchID], faceI)
//{
//      mfrh = gSum(rho.boundaryField()[hotOutletPatchID]*mag(U.boundaryField()[hotOutletPatchID][faceI])*mesh.magSf().boundaryField()[hotOutletPatchID][faceI]);
      mfrh = gSum(rho.boundaryField()[hotOutletPatchID]*mag(U.boundaryField()[hotOutletPatchID])*mesh.magSf().boundaryField()[hotOutletPatchID]);
//}
      
Info << " " << endl;
Info << "cold outlet pressure "<< numc << p_averc/(1e-8+numc) << " " << M_max << endl;
Info << "hot outlet pressure "<< numh << p_averh/(1e-8+numh) << " " << M_max << endl;
Info << " " << endl;
Info << "inlet mass flow rate "<< mfri*6 << " kg/s" << endl;
Info << "cold outlet mass flow rate "<< mfrc*6 << " kg/s" << endl;
Info << "hot outlet mass flow rate "<< mfrh*6 << " kg/s" << endl;
Info << " " << endl;
Info << "outflow mass fraction "<< (mfrc+mfrh)/(mfri+1e-10) << endl;
Info << "cold mass fraction "<< mfrc/(mfri+1e-10) << endl;
Info << " " << endl;
