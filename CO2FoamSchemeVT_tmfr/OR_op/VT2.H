/*
// Inlet
   scalar e_inlet = -86347.8; //-198513.82; //-60273.5;
   scalar rho_inlet = 185.8; //582.22; //7.06;

// Hot outlet
//   scalar e_outleth = -71306.8; //-196343.42; //-60239.6;
//   scalar rho_outleth = 32.34; //6.97;

// Cold outlet
//   scalar e_outletc = -73894.6; //-196343.42; //-60239.6;
//   scalar rho_outletc = 33.12; //6.97;

// Vortex tube specific lengths
   scalar xColdIn  = 0.01647;
   scalar xHotIn   = xColdIn+0.0009;
   scalar xHotOut  = xHotIn+0.0816;
//   scalar rBaseNozzle = 0.00414;
*/
// BC conditions
   scalar p_inlet = 8e6; //10e6;
   scalar T_inlet = 336;

   scalar p_outleth  = 5.294e6;
   scalar T_outleth  = 304.7;  // For info
   
   scalar p_outletc  = 3.906e6;
   scalar T_outletc  = 301.7;  // For info
/*
// Characteristic length
   scalar lx = xHotOut; // Length of the vortex tube
   
// Relaxation factor
   scalar relaxHot = 10.0;
   scalar relaxCold = 10.0;
//
*/
Info<< "\n***************OR_vortextube Operating Condition******************\n" << endl;
Info<< "\nInlet: p = "<<p_inlet<<" T = "<<T_inlet<<"K\n"<< endl;
Info<< "\nHot outlet: p = "<<p_outleth<<" (T = "<<T_outleth<<"K)\n"<< endl;
Info<< "\nCold outlet: p = "<<p_outletc<<" (T = "<<T_outletc<<"K)\n"<< endl;
Info<< "\n****************************************************************\n" << endl;
