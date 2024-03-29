/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::CO2Thermo

Description
    Basic thermodynamic properties based on density

SourceFiles
    CO2Thermo.C

\*---------------------------------------------------------------------------*/

#ifndef CO2Thermo_H
#define CO2Thermo_H

#include "fluidThermo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class CO2Thermo Declaration
\*---------------------------------------------------------------------------*/

class CO2Thermo
:
    public fluidThermo
{

protected:

    // Protected data

        //- Density field [kg/m^3]
        //  Named 'CO2Thermo' to avoid (potential) conflict with solver density
        volScalarField rho_;

        //- Compressibility [s^2/m^2]
//        volScalarField psi_;

        //- Dynamic viscosity [kg/m/s]
        volScalarField mu_;


    // Protected Member Functions

        //- Construct as copy (not implemented)
        CO2Thermo(const CO2Thermo&);


public:

    //- Runtime type information
    TypeName("CO2Thermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        CO2Thermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Constructors

        //- Construct from mesh and phase name
        CO2Thermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh, dictionary and phase name
        CO2Thermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );


    //- Selector
    static autoPtr<CO2Thermo> New
    (
        const fvMesh&,
        const word& phaseName=word::null
    );


    //- Destructor
    virtual ~CO2Thermo();


    // Member functions

        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho();

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const;


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual tmp<volScalarField> mu() const;
//            virtual tmp<volScalarField> mu(const scalar rho, const scalar T);

            //- Dynamic viscosity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> mu(const label patchi) const;
//            virtual tmp<scalarField> mu(const label patchi,const scalar rho, const scalar T);
//
            virtual volScalarField&  mu();
//
            virtual tmp<volScalarField> nu() const;
            virtual tmp<scalarField> nu(const label patchi) const;
//            virtual volScalarField&  nu();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
