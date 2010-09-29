/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    scalarTransportPOD

Description

\*---------------------------------------------------------------------------*/

#include "scalarTransportPOD.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(scalarTransportPOD, 0);

    addToRunTimeSelectionTable
    (
        PODODE,
        scalarTransportPOD,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::scalarTransportPOD::calcOrthoBase() const
{
    if (orthoBasePtr_)
    {
        FatalErrorIn
        (
            "scalarTransportPOD::calcOrthoBase()"
        )   << "Orthogonal base already calculated"
            << abort(FatalError);
    }

    // Create ortho-normal base
    scalar accuracy = readScalar(dict().lookup("accuracy"));

    // Get times list
    Time& runTime = const_cast<Time&>(mesh().time());

    // Remember time index to restore it after the scan
    label origTimeIndex = runTime.timeIndex();

    instantList Times = runTime.times();

    // Create a list of snapshots
    PtrList<volScalarField> fields(Times.size());

    label nSnapshots = 0;

    forAll (Times, i)
    {
        if (Times[i].value() < SMALL || Times[i] == runTime.constant())
        {
            Info << "Skipping time " << Times[i] << endl;

            continue;
        }

        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        IOobject phiHeader
        (
            phiName_,
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ
        );

        if (phiHeader.headerOk())
        {
            Info<< "    Reading " << phiName_ << endl;

            fields.set(nSnapshots, new volScalarField(phiHeader, mesh()));

            // Rename the field
            fields[nSnapshots].rename(phiName_ + name(i));
            nSnapshots++;
        }
        else
        {
            Info<< "    No " << phiName_ << endl;
        }
    }

    // Reset time index to initial state
    runTime.setTime(Times[origTimeIndex], origTimeIndex);

    // Resize snapshots
    if (nSnapshots < 2)
    {
        FatalErrorIn
        (
            "scalarTransportPOD::calcOrthoBase()"
        )   << "Insufficient number of snapshots: " << nSnapshots
            << abort(FatalError);
    }

    Info << "Number of snapshots: " << nSnapshots << endl;

    fields.setSize(nSnapshots);

    // Create ortho-normal base for transported variable
    orthoBasePtr_ = new scalarPODOrthoNormalBase(fields, accuracy);
}


void Foam::scalarTransportPOD::calcDerivativeCoeffs() const
{
    if (derivativeMatrixPtr_)
    {
        FatalErrorIn
        (
            "void scalarTransportPOD::calcDerivativeCoeffs() const"
        )   << "Derivative matrix already calculated"
            << abort(FatalError);
    }

    // Calculate coefficients for differential equation
    // Get times list
    Time& runTime = const_cast<Time&>(this->mesh().time());

    // Remember time index to restore it
    label origTimeIndex = runTime.timeIndex();

    // Read diffusivity

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading diffusivity D\n" << endl;

    dimensionedScalar DT
    (
        transportProperties.lookup("DT")
    );

    // Find velocity field

    word Uname(this->dict().lookup("velocity"));

    instantList Times = runTime.times();

    volVectorField* Uptr = NULL;

    forAll (Times, i)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        IOobject Uheader
        (
            Uname,
            runTime.timeName(),
            this->mesh(),
            IOobject::MUST_READ
        );

        if (Uheader.headerOk())
        {
            Info<< "    Reading " << Uname << endl;

            Uptr = new volVectorField(Uheader, this->mesh());
            break;
        }
        else
        {
            Info<< "    No " << Uname << endl;
        }
    }

    // Reset time index to initial state
    runTime.setTime(Times[origTimeIndex], origTimeIndex);

    if (!Uptr)
    {
        FatalErrorIn
        (
            "void scalarTransportPOD::calcDerivativeCoeffs() const"
        )   << "Cannot find velocity: " << Uname
            << abort(FatalError);
    }

    volVectorField& U = *Uptr;

    // Create derivative matrix

    const scalarPODOrthoNormalBase& b = orthoBase();

    derivativeMatrixPtr_ = new scalarSquareMatrix(b.baseSize(), 0.0);
    scalarSquareMatrix& derivative = *derivativeMatrixPtr_;

    for (label i = 0; i < b.baseSize(); i++)
    {
        const volScalarField& snapI = b.orthoField(i);

        volVectorField gradSnapI = fvc::grad(snapI);

        for (label j = 0; j < b.baseSize(); j++)
        {
            const volScalarField& snapJ = b.orthoField(j);

            volVectorField gradSnapJ = fvc::grad(snapJ);

            derivative[i][j] =
                DT.value()*POD::projection(fvc::div(gradSnapJ), snapI)
              - POD::projection((U & gradSnapJ), snapI);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::scalarTransportPOD::scalarTransportPOD
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    PODODE(mesh, dict),
    phiName_(dict.lookup("field")),
    coeffs_(),
    derivativeMatrixPtr_(NULL),
    orthoBasePtr_(NULL),
    fieldPtr_(NULL)
{
    // Grab coefficients from the first snapshot of the ortho-normal base
    coeffs_.setSize(orthoBase().baseSize());

    const scalarRectangularMatrix& orthoBaseCoeffs =
        orthoBase().interpolationCoeffs();

    forAll (coeffs_, i)
    {
        coeffs_[i] = orthoBaseCoeffs[0][i];
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::scalarTransportPOD::~scalarTransportPOD()
{
    deleteDemandDrivenData(derivativeMatrixPtr_);

    clearBase();
    clearFields();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::scalarTransportPOD::nEqns() const
{
    return coeffs().size();
}


Foam::scalarField& Foam::scalarTransportPOD::coeffs()
{
    return coeffs_;
}


const Foam::scalarField& Foam::scalarTransportPOD::coeffs() const
{
    return coeffs_;
}


void Foam::scalarTransportPOD::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    if (!derivativeMatrixPtr_)
    {
        calcDerivativeCoeffs();
    }

   const scalarSquareMatrix& derivative = *derivativeMatrixPtr_;

    forAll (dydx, i)
    {
        dydx[i] = 0;

        forAll (y, j)
        {
            dydx[i] += derivative[i][j]*y[j];
        }
    }
}


void Foam::scalarTransportPOD::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    derivatives(x, y, dfdx);
    dfdy = 0;
}


const Foam::scalarPODOrthoNormalBase&
Foam::scalarTransportPOD::orthoBase() const
{
    if (!orthoBasePtr_)
    {
        calcOrthoBase();
    }

    return *orthoBasePtr_;
}


void Foam::scalarTransportPOD::clearBase() const
{
    deleteDemandDrivenData(orthoBasePtr_);
}


const Foam::volScalarField& Foam::scalarTransportPOD::field() const
{
    if (!fieldPtr_)
    {
        updateFields();
    }

    return *fieldPtr_;
}


void Foam::scalarTransportPOD::updateFields() const
{
    if (!fieldPtr_)
    {
        // Allocate field
        fieldPtr_ =
            new volScalarField
            (
                IOobject
                (
                    phiName_ + "POD",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar
                (
                    "zero",
                    orthoBase().orthoField(0).dimensions(),
                    0
                )
            );
    }

    volScalarField& phi = *fieldPtr_;

    const scalarPODOrthoNormalBase& b = orthoBase();

    phi = dimensionedScalar("zero", b.orthoField(0).dimensions(), 0);

    forAll (coeffs_, i)
    {
        phi += coeffs_[i]*b.orthoField(i);
    }
}


void Foam::scalarTransportPOD::clearFields() const
{
    deleteDemandDrivenData(fieldPtr_);
}


void Foam::scalarTransportPOD::write() const
{
    // Recalculate field and force a write
    updateFields();
    field().write();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
