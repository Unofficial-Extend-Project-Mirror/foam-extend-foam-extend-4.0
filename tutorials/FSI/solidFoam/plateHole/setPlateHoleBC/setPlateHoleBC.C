/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "setPlateHoleBC.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "surfaceFields.H"

#include "ValuePointPatchField.H"
#include "tractionDisplacementFvPatchVectorField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setPlateHoleBC, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        setPlateHoleBC,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Foam::symmTensor Foam::setPlateHoleBC::plateHoleStress(const vector& C) const
{
    tensor sigma = tensor::zero;

    scalar T = 10000;
    scalar a = 0.5;

    scalar r = ::sqrt(sqr(C.x()) + sqr(C.y()));
    scalar theta = atan2(C.y(), C.x());

    sigma.xx() =
        T
       *(
           1.0
         - (sqr(a)/sqr(r))*(3*cos(2*theta)/2 + cos(4*theta))
         + (3*pow(a,4)/(2*pow(r,4)))*cos(4*theta)
        );

    sigma.yy() =
        T
       *(
         - (sqr(a)/sqr(r))*(cos(2*theta)/2 - cos(4*theta))
         - (3*pow(a,4)/(2*pow(r,4)))*cos(4*theta)
        );

    sigma.xy() =
        T
       *(
         - (sqr(a)/sqr(r))*(sin(2*theta)/2 + sin(4*theta))
         + (3*pow(a,4)/(2*pow(r,4)))*sin(4*theta)
        );

    sigma.yx() = sigma.xy();

//     coordinateSystem cs("polarCS", C, vector(0, 0, 1), C/mag(C));

//     sigma.xx() =
//         T*(1 - sqr(a)/sqr(r))/2
//       + T*(1 + 3*pow(a,4)/pow(r,4) - 4*sqr(a)/sqr(r))*cos(2*theta)/2;

//     sigma.xy() =
//       - T*(1 - 3*pow(a,4)/pow(r,4) + 2*sqr(a)/sqr(r))*sin(2*theta)/2;

//     sigma.yx() = sigma.xy();

//     sigma.yy() =
//         T*(1 + sqr(a)/sqr(r))/2
//       - T*(1 + 3*pow(a,4)/pow(r,4))*cos(2*theta)/2;

//     // Transformation to global coordinate system
//     sigma = ((cs.R() & sigma) & cs.R().T());

    symmTensor S = symmTensor::zero;

    S.xx() = sigma.xx();
    S.xy() = sigma.xy();
    S.yy() = sigma.yy();

    return S;
}


Foam::vector Foam::setPlateHoleBC::plateHoleDisplacement(const vector& C) const
{
    vector displacement = vector::zero;

    scalar T = 10000;
    scalar a = 0.5;

    scalar E = 2e11;
    scalar nu = 0.3;

    scalar mu = E/(2*(1.0 + nu));
    scalar kappa = (3.0-nu)/(1.0+nu);

    scalar r = ::sqrt(sqr(C.x()) + sqr(C.y()));
    scalar theta = atan2(C.y(), C.x());

    displacement.x() =
        (a*T/(8*mu))
       *(
           (r/a)*(kappa+1)*cos(theta)
         + (2*a/r)*((1+kappa)*cos(theta) + cos(3*theta))
         - (2*pow(a,3)/pow(r,3))*cos(3*theta)
        );

    displacement.y() =
        (a*T/(8*mu))
       *(
           (r/a)*(kappa-3)*sin(theta)
         + (2*a/r)*((1-kappa)*sin(theta) + sin(3*theta))
         - (2*pow(a,3)/pow(r,3))*sin(3*theta)
        );

    return displacement;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setPlateHoleBC::setPlateHoleBC
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion)
{
//     if (Pstream::parRun())
//     {
//         FatalErrorIn("setPlateHoleBC::setPlateHoleBC(...)")
//             << "setPlateHoleBC objec function "
//                 << "is not implemented for parallel run"
//                 << abort(FatalError);
//     }

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::setPlateHoleBC::start()
{
    setBC();

    return false;
}


bool Foam::setPlateHoleBC::execute()
{
    calcError();

    return false;
}


void Foam::setPlateHoleBC::setBC()
{
    Info << "Update boundary conditions" << endl;

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    volVectorField& D =
        const_cast<volVectorField&>(mesh.lookupObject<volVectorField>("D"));

    forAll(D.boundaryField(), patchI)
    {
        if
        (
            D.boundaryField()[patchI].type()
         == tractionDisplacementFvPatchVectorField::typeName
        )
        {
            tractionDisplacementFvPatchVectorField& patchD =
                refCast<tractionDisplacementFvPatchVectorField>
                (
                    D.boundaryField()[patchI]
                );

            vectorField nf = patchD.patch().nf();

            const vectorField& Cf = mesh.Cf().boundaryField()[patchI];

            forAll(patchD.traction(), faceI)
            {
                vector curC(Cf[faceI].x(), Cf[faceI].y(), 0);
                vector curN = nf[faceI];

                if (mesh.boundary()[patchI].name() == "hole")
                {
                    curC /= mag(curC);
                    curC *= 0.5;

                    curN = -curC/mag(curC);
                }

                patchD.traction()[faceI] = (nf[faceI]&plateHoleStress(curC));
            }
        }
    }
}


void Foam::setPlateHoleBC::calcError() const
{
    Info << "Calculating errors" << endl;

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    // Cell displacement
    {
        const volVectorField& D =
            mesh.lookupObject<volVectorField>("D");

        const vectorField& DI = D.internalField();

        const vectorField& C = mesh.C().internalField();

        scalarField DError(D.size(), 0);

        forAll(DError, cellI)
        {
            vector curR = vector(C[cellI].x(), C[cellI].y(), 0);
            vector curDa = plateHoleDisplacement(curR);

            DError[cellI] = mag(DI[cellI] - curDa);
        }

        Info << "DError, max : " << gMax(DError) << endl;
    }


    // Point displacement
    {
        const pointVectorField& pointD =
            mesh.lookupObject<pointVectorField>("pointD");

        const vectorField& pointDI = pointD.internalField();

        const vectorField& C = mesh.points();

        scalarField pointDError(pointDI.size(), 0);

        forAll(pointDError, pointI)
        {
            vector curR = vector(C[pointI].x(), C[pointI].y(), 0);
            vector curDa = plateHoleDisplacement(curR);

            pointDError[pointI] = mag(pointDI[pointI] - curDa);
        }

        Info << "pointDError, max : " << gMax(pointDError) << endl;
    }


    if (time_.outputTime())
    {
        // Stress error field for plate with hole case
        volScalarField sigmaXXErr
        (
            IOobject
            (
                "sigmaXXErr",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0)
        );

        volScalarField sigmaXYErr
        (
            IOobject
            (
                "sigmaXYErr",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0)
        );

        volScalarField sigmaYYErr
        (
            IOobject
            (
                "sigmaYYErr",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0)
        );

        scalarField& sigmaXXErrI = sigmaXXErr.internalField();
        scalarField& sigmaXYErrI = sigmaXYErr.internalField();
        scalarField& sigmaYYErrI = sigmaYYErr.internalField();

        const volSymmTensorField& sigma =
            mesh.lookupObject<volSymmTensorField>("sigma");

        const symmTensorField& sigmaI = sigma.internalField();
        const vectorField& C = mesh.C().internalField();

//         scalar maxSigmaXX = max(mag(sigma.component(symmTensor::XX))).value();
//         scalar maxSigmaXY = max(mag(sigma.component(symmTensor::XY))).value();
//         scalar maxSigmaYY = max(mag(sigma.component(symmTensor::YY))).value();

        forAll(sigmaI, cellI)
        {
            vector curR = vector(C[cellI].x(), C[cellI].y(), 0);
            symmTensor curSigmaA = plateHoleStress(curR);

            sigmaXXErrI[cellI] =
                mag(sigmaI[cellI].xx() - curSigmaA.xx());
//                 mag(sigmaI[cellI].xx() - curSigmaA.xx())/maxSigmaXX;

            sigmaXYErrI[cellI] =
                mag(sigmaI[cellI].xy() - curSigmaA.xy());
//                 mag(sigmaI[cellI].xy() - curSigmaA.xy())/maxSigmaXY;

            sigmaYYErrI[cellI] =
                mag(sigmaI[cellI].yy() - curSigmaA.yy());
//                 mag(sigmaI[cellI].yy() - curSigmaA.yy())/maxSigmaYY;
        }

        Info << "sigmaXXErr, max : " << gMax(sigmaXXErr) << endl;
        Info << "sigmaXYErr, max : " << gMax(sigmaXYErr) << endl;
        Info << "sigmaYYErr, max : " << gMax(sigmaYYErr) << endl;

        sigmaXXErr.write();
        sigmaXYErr.write();
        sigmaYYErr.write();
    }
}


bool Foam::setPlateHoleBC::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
