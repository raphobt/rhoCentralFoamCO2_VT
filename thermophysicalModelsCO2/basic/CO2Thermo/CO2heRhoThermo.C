/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "CO2heRhoThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::CO2heRhoThermo<BasicPsiThermo, MixtureType>::calculate()
{
//    const scalarField& hCells = this->he();
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
//    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& rhoCells = this->rho();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();
//Info<<"rho "<<rhoCells<<endl;
    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

//        TCells[celli] = mixture_.THE
//        (
//            hCells[celli],
//            pCells[celli],
//            TCells[celli]
//        );

//        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
//        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);

        muCells[celli] = 0.0;// mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = 0.0;//mixture_.alphah(pCells[celli], TCells[celli]);
    }

//    volScalarField::Boundary& pBf =
//        this->p_.boundaryFieldRef();
//
    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();
//
//    volScalarField::Boundary& psiBf =
//        this->psi_.boundaryFieldRef();
//
//    volScalarField::Boundary& rhoBf =
//        this->rho_.boundaryFieldRef();
//
//    volScalarField::Boundary& heBf =
//        this->he().boundaryFieldRef();
//
    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();
//
    volScalarField::Boundary& alphaBf =
        this->alpha_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
//        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
//        fvPatchScalarField& ppsi = psiBf[patchi];
//        fvPatchScalarField& prho = rhoBf[patchi];
//        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

//                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

//                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
//                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = 0.0;//mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = 0.0;//mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

//                pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);

//                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
//                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = 0.0;//mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = 0.0;//mixture_.alphah(pp[facei], pT[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::CO2heRhoThermo<BasicPsiThermo, MixtureType>::CO2heRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo_CO2<BasicPsiThermo, MixtureType>(mesh, phaseName)
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::CO2heRhoThermo<BasicPsiThermo, MixtureType>::~CO2heRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::CO2heRhoThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}


// ************************************************************************* //
