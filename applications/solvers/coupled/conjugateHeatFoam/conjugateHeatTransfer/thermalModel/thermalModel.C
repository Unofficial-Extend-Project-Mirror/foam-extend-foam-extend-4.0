/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 Hrvoje Jasak
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    thermalModel

Description
    Thermal  material properties for solids.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "thermalModel.H"
#include "volFields.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermalModel, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalModel::thermalModel(const volScalarField& T)
:
    IOdictionary
    (
        IOobject
        (
            "thermalProperties",
            T.time().constant(),
            T.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    T_(T),
    lawPtr_(thermalLaw::New("law", T_, subDict("thermal")))
{
    {
        PtrList<entry> entries(subDict("thermal").lookup("gaps"));
        gapPtr_.setSize(entries.size());

        forAll (gapPtr_, gapI)
        {
            gapPtr_.set
            (
                gapI,
                thermalGap::New
                (
                    entries[gapI].keyword(),
                    T,
                    entries[gapI].dict()
                )
            );
        }
    }

    {
        PtrList<entry> entries(subDict("thermal").lookup("sources"));
        sourcePtr_.setSize(entries.size());

        forAll (sourcePtr_, sourceI)
        {
            sourcePtr_.set
            (
                sourceI,
                thermalSource::New
                (
                    entries[sourceI].keyword(),
                    T,
                    entries[sourceI].dict()
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermalModel::modifyResistance
(
    surfaceScalarField& kf
) const
{
    forAll(gapPtr_, gapI)
    {
        gapPtr_[gapI].modifyResistance(kf);
    }
}

tmp<fvScalarMatrix> thermalModel::laplacian(volScalarField& T)
{
    word kScheme ("laplacian(k,T)");

    surfaceScalarField kf = fvc::interpolate(lawPtr_->k());

    modifyResistance(kf);

    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix( fvm::laplacian(kf, T, kScheme) )
    );
}

tmp<volScalarField> thermalModel::S() const
{
    tmp<volScalarField> tsource
    (
        new volScalarField
        (
            IOobject
            (
                "heatSource",
                T_.time().timeName(),
                T_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimensionedScalar
            (
                "zero",
                dimEnergy/dimTime/dimVolume,
                scalar(0.0)
            )
        )
    );

    forAll(sourcePtr_, sourceI)
    {
        sourcePtr_[sourceI].addSource(tsource());
    }

    return tsource;
}

bool thermalModel::read()
{
    if (regIOobject::read())
    {
        lawPtr_ = thermalLaw::New("law", T_, subDict("thermal"));

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
