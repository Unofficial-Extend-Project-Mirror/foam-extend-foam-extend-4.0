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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "bubbleHistory.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "boundBox.H"
#include "freeSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bubbleHistory, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        bubbleHistory,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bubbleHistory::bubbleHistory
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    V0_(SMALL),
    historyFilePtr_(NULL)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    if (Pstream::master())
    {
        if (!time_.processorCase())
        {
            mkDir
            ( 
                time_.path()
               /"history"
               /time_.timeName()
            );
        
            historyFilePtr_ = 
                new OFstream
                (
                    time_.path()
                   /"history"
                   /time_.timeName()
                   /"history.dat"
                );
        }
        else
        {
            mkDir
            ( 
                time_.path()/".."/"history"
               /time_.timeName()
            );
        
            historyFilePtr_ = 
                new OFstream
                (
                    time_.path()/".."
                   /"history"
                   /time_.timeName()
                   /"history.dat"
                );
        }

        (*historyFilePtr_) 
            << "Time" << tab 
                << "Cx" << tab
                << "Cy" << tab
                << "Cz" << tab
                << "Ux" << tab
                << "Uy" << tab
                << "Uz" << tab
                << "ax" << tab
                << "ay" << tab
                << "az" << tab
                << "Fx" << tab
                << "Fy" << tab
                << "Fz" << tab
                << "D" << tab
                << "Lx" << tab
                << "Ly" << tab
                << "Lz" << tab
                << "CD" << tab
                << "CVM" << tab
                << "RE" << tab
                << "WE" << tab
                << "V/V0" << tab
                << "A" << tab
                << "B" << tab
                << "C" << endl;
    }

    const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName_);
    
    const volScalarField& fluidIndicator = 
        mesh.lookupObject<volScalarField>("fluidIndicator");

    V0_ = gSum((1 - fluidIndicator.internalField())*mesh.V().field());

    const freeSurface& fs = 
        mesh.lookupObject<freeSurface>("freeSurfaceProperties");

    if (!fs.twoFluids())
    {
        V0_ = 
          - gSum
            (
                mesh.Cf().boundaryField()[fs.aPatchID()]
              & mesh.Sf().boundaryField()[fs.aPatchID()]
            )/3;

        if (mesh.nGeometricD() != 3)
        {
            FatalErrorIn("bubbleHistory::")
                << "One-fluid bubble centroid calc "
                    << "is not implemented for 2d mesh"
                    << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::bubbleHistory::start()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const volScalarField& fluidIndicator = 
        mesh.lookupObject<volScalarField>("fluidIndicator");

    const freeSurface& fs = 
        mesh.lookupObject<freeSurface>("freeSurfaceProperties");

    scalar V = gSum((1 - fluidIndicator.internalField())*mesh.V().field());

    if (!fs.twoFluids())
    {
        V =
          - gSum
            (
                mesh.Cf().boundaryField()[fs.aPatchID()]
                & mesh.Sf().boundaryField()[fs.aPatchID()]
            )/3;

        if (mesh.nGeometricD() != 3)
        {
            FatalErrorIn("bubbleHistory::")
                << "One-fluid bubble centroid calc "
                    << "is not implemented for 2d mesh"
                    << abort(FatalError);
        }
    }

    const IOdictionary& mrf = 
        mesh.lookupObject<IOdictionary>("movingReferenceFrame");

    dimensionedVector C(mrf.lookup("XF"));
    dimensionedVector U(mrf.lookup("UF"));
    dimensionedVector a(mrf.lookup("aF"));

    vector F = fs.totalViscousForce() + fs.totalPressureForce();

    vector dragDir;
    scalar Uref;

    if(mag(U.value()) > SMALL)
    {
        dragDir = -U.value()/mag(U.value());
        Uref = mag(U.value());
    }
    else
    {
        dragDir = fs.g().value()/(mag(fs.g().value()) + SMALL);
        Uref = SMALL;
    }

    scalar dragForce = (dragDir&F);

    vector liftForce =  transform(I - dragDir*dragDir, F);

    scalar Deq = 2*pow(3*V/(4*M_PI), 1.0/3.0);

    scalar Ug = -(U.value()&(fs.g().value()/(mag(fs.g().value()) + SMALL)));

    scalar CD = 
        (4.0/3.0)*(fs.rhoFluidA().value() - fs.rhoFluidB().value())
       *mag(fs.g().value())*Deq
       /(fs.rhoFluidA().value()*mag(U.value())*Ug + SMALL);

    scalar ag = -(a.value()&(fs.g().value()/(mag(fs.g().value()) + SMALL)));

    scalar CVM = 
        (fs.rhoFluidA().value() - fs.rhoFluidB().value())*mag(fs.g().value())
       /(fs.rhoFluidA().value()*ag + SMALL)
      - (fs.rhoFluidB().value()/(fs.rhoFluidA().value() + SMALL));

    scalar RE = Ug*Deq*fs.rhoFluidA().value()/(fs.muFluidA().value() + SMALL);

    scalar WE = fs.rhoFluidA().value()*sqr(Ug)*Deq
       /(fs.cleanInterfaceSurfTension().value() + SMALL);

    boundBox box(mesh.C().boundaryField()[fs.aPatchID()]);

    if (Pstream::master())
    {
        historyFilePtr_->precision(12);

        (*historyFilePtr_) << time_.value() << tab 
            << C.value().x() << tab
            << C.value().y() << tab
            << C.value().z() << tab
            << U.value().x() << tab
            << U.value().y() << tab
            << U.value().z() << tab
            << a.value().x() << tab
            << a.value().y() << tab
            << a.value().z() << tab
            << F.x() << tab
            << F.y() << tab
            << F.z() << tab
            << dragForce << tab
            << liftForce.x() << tab
            << liftForce.y() << tab
            << liftForce.z() << tab
            << CD << tab
            << CVM << tab
            << RE << tab
            << WE << tab
            << mag(1.0 - V/V0_) << tab
            << (box.max().x()-box.min().x())/2 << tab
            << (box.max().y()-box.min().y())/2 << tab
            << (box.max().z()-box.min().z())/2 << endl;

        return true;
    }

    return false;
}


bool Foam::bubbleHistory::execute()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const volScalarField& fluidIndicator = 
        mesh.lookupObject<volScalarField>("fluidIndicator");

    const freeSurface& fs = 
        mesh.lookupObject<freeSurface>("freeSurfaceProperties");

    scalar V = gSum((1 - fluidIndicator.internalField())*mesh.V().field());

    if (!fs.twoFluids())
    {
        V =
          - gSum
            (
                mesh.Cf().boundaryField()[fs.aPatchID()]
                & mesh.Sf().boundaryField()[fs.aPatchID()]
            )/3;

        if (mesh.nGeometricD() != 3)
        {
            FatalErrorIn("bubbleHistory")
                << "One-fluid bubble centroid calc "
                    << "is not implemented for 2d mesh"
                    << abort(FatalError);
        }
    }

    const IOdictionary& mrf = 
        mesh.lookupObject<IOdictionary>("movingReferenceFrame");

    dimensionedVector C(mrf.lookup("XF"));
    dimensionedVector U(mrf.lookup("UF"));
    dimensionedVector a(mrf.lookup("aF"));


    vector F = fs.totalViscousForce() + fs.totalPressureForce();

    vector dragDir;
    scalar Uref;

    if(mag(U.value()) > SMALL)
    {
        dragDir = -U.value()/mag(U.value());
        Uref = mag(U.value());
    }
    else
    {
        dragDir = fs.g().value()/(mag(fs.g().value()) + SMALL);
        Uref = SMALL;
    }

    scalar dragForce = (dragDir&F);

    vector liftForce =  transform(I - dragDir*dragDir, F);

    scalar Deq = pow(3*V/(4*M_PI), 1.0/3.0);

    scalar Ug = -(U.value()&(fs.g().value()/(mag(fs.g().value()) + SMALL)));

    scalar CD = 
        (4.0/3.0)*(fs.rhoFluidA().value() - fs.rhoFluidB().value())
       *mag(fs.g().value())*Deq
       /(fs.rhoFluidA().value()*mag(U.value())*Ug + SMALL);

    scalar ag = -(a.value()&(fs.g().value()/(mag(fs.g().value()) + SMALL)));

    scalar CVM = 
        (fs.rhoFluidA().value() - fs.rhoFluidB().value())*mag(fs.g().value())
       /(fs.rhoFluidA().value()*ag + SMALL)
      - (fs.rhoFluidB().value()/(fs.rhoFluidA().value() + SMALL));

    scalar RE = Ug*Deq*fs.rhoFluidA().value()/(fs.muFluidA().value() + SMALL);

    scalar WE = fs.rhoFluidA().value()*sqr(Ug)*Deq
       /(fs.cleanInterfaceSurfTension().value() + SMALL);

    boundBox box(mesh.C().boundaryField()[fs.aPatchID()]);

    if (Pstream::master())
    {
        historyFilePtr_->precision(12);

        (*historyFilePtr_) << time_.value() << tab 
            << C.value().x() << tab
            << C.value().y() << tab
            << C.value().z() << tab
            << U.value().x() << tab
            << U.value().y() << tab
            << U.value().z() << tab
            << a.value().x() << tab
            << a.value().y() << tab
            << a.value().z() << tab
            << F.x() << tab
            << F.y() << tab
            << F.z() << tab
            << dragForce << tab
            << liftForce.x() << tab
            << liftForce.y() << tab
            << liftForce.z() << tab
            << CD << tab
            << CVM << tab
            << RE << tab
            << WE << tab
            << mag(1.0 - V/V0_) << tab
            << (box.max().x()-box.min().x())/2 << tab
            << (box.max().y()-box.min().y())/2 << tab
            << (box.max().z()-box.min().z())/2 << endl;

        return true;
    }

    return false;
}


bool Foam::bubbleHistory::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
