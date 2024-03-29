/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::CO2InternalEnergy

Description
    Thermodynamics mapping class to expose the sensible internal energy function
    as the standard internal energy function e(T).

\*---------------------------------------------------------------------------*/

#ifndef CO2InternalEnergy_H
#define CO2InternalEnergy_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class CO2InternalEnergy Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class CO2InternalEnergy
{

public:

    // Constructors

        //- Construct
        CO2InternalEnergy()
        {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "CO2InternalEnergy";
        }

        // Fundamental properties

            static word name()
            {
                return "e";
            }

            //- Sensible Internal energy [J/kmol]
//            scalar he
//            (
//                const Thermo& thermo,
//                const scalar p,
//                const scalar T) const
//            {
//                return 1.0; //thermo.es(p, T);
//            }

            //- Heat capacity at constant volume [J/(kmol K)]
//            scalar cpv
//            (
//                const Thermo& thermo,
//                const scalar p,
//                const scalar T
//            ) const
//            {
//                return thermo.cv(p, T);
//            }

            //- cp/cv []
//            scalar cpBycpv
//            (
//                const Thermo& thermo,
//                const scalar p,
//                const scalar T
//            ) const
//            {
//                return thermo.gamma(p, T);
//            }


            //- Sensible internal energy [J/kg]
//            scalar HE
//            (
//                const Thermo& thermo,
//                const scalar p,
//                const scalar T
//            ) const
//            {
//                return 1.0; //thermo.Es(p, T);
//           }

            //- Temperature from sensible internal energy
            //  given an initial temperature T0
//            scalar THE
//            (
//                const Thermo& thermo,
//                const scalar e,
//                const scalar p,
//                const scalar T0
//            ) const
//            {
//                return thermo.TEs(e, p, T0);
//            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
