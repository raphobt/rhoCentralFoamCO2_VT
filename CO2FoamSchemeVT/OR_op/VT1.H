// Inlet
   scalar e_inlet = -60273.5;
   scalar rho_inlet = 7.05997;

// Hot outlet
   scalar e_outleth = -58027.3;
   scalar rho_outleth = 4.33378;

// Cold outlet
   scalar e_outletc = -65961.6;
   scalar rho_outletc = 2.49154;

// Vortex tube specific lengths
   scalar xColdIn  = 0.01647;
   scalar xHotIn   = 0.01747;
   scalar xHotOut  = 0.09907;

// Initial velocity
 /*scalar velx = 1e-6;
   scalar vely = -199.0;
   scalar velz = -199.0;*/

// BC conditions
   scalar p_inlet = 0.387e6;
   scalar T_inlet = 22.9+273.15;

   scalar p_outleth  = 0.241e6;
   scalar T_outleth  = 24.8+273.15;
   
   scalar p_outletc  = 0.133e6;
   scalar T_outletc  = 11.6+273.15;

// Characteristic lengths
   scalar lx = xHotOut;
   
// Relaxation factor
   scalar relax = 50;
//
//   
Info<< "\n***************OR_vortextube Operating Condition******************\n" << endl;
Info<< "\nInlet: p= "<<p_inlet<<" T= "<<T_inlet<<"K\n" << endl;
Info<< "\nHot outlet: p= "<<p_outleth<<" T= "<<T_outleth<<"K\n" << endl;
Info<< "\nCold outlet: p= "<<p_outletc<<" T= "<<T_outletc<<" relax= " << relax << endl;
Info<< "\n****************************************************************\n" << endl;
