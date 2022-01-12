label inletPatchID = mesh.boundaryMesh().findPatchID("inlet");
label coldOutletPatchID = mesh.boundaryMesh().findPatchID("coldOutlet");
label hotOutletPatchID = mesh.boundaryMesh().findPatchID("hotOutlet");
scalar mfri=0.0, mfrc=0.0, mfrh=0.0;
//scalar pc=0.0, ph=0.0, Tc=0.0, Th=0.0;
//scalar nc=0, nh=0;

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
//forAll(mesh.boundary()[coldOutletPatchID], faceI)
//{
//      pc = gSum(p.boundaryField()[coldOutletPatchID]);
//      nc = 1.0 + nc;
//}

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
Info << "out/in mass fraction "<< (mfrc+mfrh)/(mfri+1e-10) << endl;
Info << "cold mass fraction "<< mfrc/(mfri+1e-10) << endl;
Info << " " << endl;
