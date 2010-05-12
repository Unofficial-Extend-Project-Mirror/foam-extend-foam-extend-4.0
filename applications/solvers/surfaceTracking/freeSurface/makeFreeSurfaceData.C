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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "freeSurface.H"
#include "primitivePatchInterpolation.H"
#include "faceTetPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void freeSurface::makeInterpolators()
{
    if (debug)
    {
        Info<< "freeSurface::makeInterpolators(): "
            << "making patch to patch interpolator"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if 
    (
        interpolatorABPtr_ ||  
        interpolatorBAPtr_
    )
    {
        FatalErrorIn("freeSurface::makeInterpolators()")
            << "patch to patch interpolators already exists"
            << abort(FatalError);
    }


    if(aPatchID() == -1)
    {
        FatalErrorIn("freeSurface::makeInterpolators()")
            << "Free surface patch A not defined."
            << abort(FatalError);
    }


    if(bPatchID() == -1)
    {
        FatalErrorIn("freeSurface::makeInterpolators()")
            << "Free surface patch B not defined."
            << abort(FatalError);
    }

    interpolatorBAPtr_ = new patchToPatchInterpolation
    (
        mesh().boundaryMesh()[bPatchID()],
        mesh().boundaryMesh()[aPatchID()],
        intersection::VISIBLE
        // intersection::HALF_RAY
    );


    interpolatorABPtr_ = new patchToPatchInterpolation
    (
        mesh().boundaryMesh()[aPatchID()],
        mesh().boundaryMesh()[bPatchID()],
        intersection::VISIBLE
        // intersection::HALF_RAY
    );
}


void freeSurface::makeControlPoints()
{
    if (debug)
    {
        Info<< "freeSurface::makeControlPoints(): "
            << "making control points"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (controlPointsPtr_)
    {
        FatalErrorIn("freeSurface::makeControlPoints()")
            << "patch to patch interpolators already exists"
            << abort(FatalError);
    }


    IOobject controlPointsHeader
    (
        "controlPoints",
        time().timeName(),
        mesh(),
        IOobject::MUST_READ
    );

    if (controlPointsHeader.headerOk())
    {
        controlPointsPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "controlPoints",
                    time().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
    }
    else
    {
        controlPointsPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "controlPoints",
                    time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                aMesh().areaCentres().internalField()
            );
    }
}


void freeSurface::makeMotionPointsMask()
{
    if (debug)
    {
        Info<< "freeSurface::makeMotionPointsMask() : "
            << "making motion points mask"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (motionPointsMaskPtr_)
    {
        FatalErrorIn("freeSurface::motionPointsMask()")
            << "motion points mask already exists"
            << abort(FatalError);
    }


    if(aPatchID() == -1)
    {
        FatalErrorIn("freeSurface::makeMotionPointsMask()")
            << "Free surface patch A not defined."
            << abort(FatalError);
    }


    motionPointsMaskPtr_ = new scalarField
    (
        mesh().boundaryMesh()[aPatchID()].nPoints(),
        1.0
    );
}


void freeSurface::makeDirections()
{
    if (debug)
    {
        Info<< "freeSurface::makeDirections() : "
            << "making displacement directions for points and "
            << "control points"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if 
    (
        pointsDisplacementDirPtr_ ||  
        facesDisplacementDirPtr_
    )
    {
        FatalErrorIn("freeSurface::makeDirections()")
            << "points and control points displacement directions "
            << "already exists"
            << abort(FatalError);
    }


    if(aPatchID() == -1)
    {
        FatalErrorIn("freeSurface::makeMotionPointsMask()")
            << "Free surface patch A not defined."
            << abort(FatalError);
    }


    pointsDisplacementDirPtr_ = 
        new vectorField
        (
            mesh().boundaryMesh()[aPatchID()].nPoints(),
            -(g_/mag(g_)).value()
        );


    facesDisplacementDirPtr_ = 
        new vectorField
        (
            mesh().boundaryMesh()[aPatchID()].size(),
            -(g_/mag(g_)).value()
        );


    updateDisplacementDirections();
}


void freeSurface::makeTotalDisplacement()
{
    if (debug)
    {
        Info<< "freeSurface::makeTotalDisplacement() : "
            << "making zero total points displacement"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (totalDisplacementPtr_)
    {
        FatalErrorIn("freeSurface::makeTotalDisplacement()")
            << "total points displacement already exists"
            << abort(FatalError);
    }


    totalDisplacementPtr_ =
        new vectorIOField
        (
            IOobject
            (
                "totalDisplacement",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            vectorField
            (
                mesh().boundaryMesh()[aPatchID()].nPoints(), 
                vector::zero
            )
        );
}
 

void freeSurface::readTotalDisplacement()
{
    if (debug)
    {
        Info<< "freeSurface::readTotalDisplacement() : "
            << "reading total points displacement if present"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (totalDisplacementPtr_)
    {
        FatalErrorIn("freeSurface::makeTotalDisplacement()")
            << "total points displacement already exists"
            << abort(FatalError);
    }


    if
    (
        IOobject
        (
            "totalDisplacement",
            time().timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ).headerOk()
    )
    {
        totalDisplacementPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "totalDisplacement",
                    time().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )     
            );
    }       
}


void freeSurface::makeFaMesh()
{
    if (debug)
    {
        Info<< "freeSurface::makeFaMesh() : "
            << "making finite area mesh"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (aMeshPtr_)
    {
        FatalErrorIn("freeSurface::makeFaMesh()")
            << "finite area mesh already exists"
            << abort(FatalError);
    }


    aMeshPtr_ = new faMesh(mesh());
}


void freeSurface::makeMeshMotionSolver()
{
    if (debug)
    {
        Info<< "freeSurface::makeMeshMotionSolver() : "
            << "making mesh motion solver"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (mSolverPtr_)
    {
        FatalErrorIn("freeSurface::makeMeshMotionSolver()")
            << "mesh motion solver already exists"
            << abort(FatalError);
    }

    mSolverPtr_ = dynamic_cast<tetDecompositionMotionSolver*>
    (
        motionSolver::New(mesh()).ptr()
    );
}

void freeSurface::makePatchPointInterpolators()
{
    if (interpolatorAPtr_ || interpolatorBPtr_)
    {
        FatalErrorIn("freeSurface::makePatchPointInterpolators()")
            << "mesh motion solver already exists"
            << abort(FatalError);
    }

    // Make interpolators
    if(aPatchID() != -1)
    {
        interpolatorAPtr_ =
            new tetPolyPatchInterpolation
            (
                refCast<const faceTetPolyPatch>
                (
                    meshMotionSolver().motionU().boundaryField()
                    [aPatchID()].patch()
                )
            );
    }

    if(bPatchID() != -1)
    {
        interpolatorBPtr_ =
            new tetPolyPatchInterpolation
            (
                refCast<const faceTetPolyPatch>
                (
                    meshMotionSolver().motionU().boundaryField()
                    [bPatchID()].patch()
                )
            );
    }
}



void freeSurface::makeUs()
{
    if (debug)
    {
        Info<< "freeSurface::makeUs() : "
            << "making free-surface velocity field"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (UsPtr_)
    {
        FatalErrorIn("freeSurface::makeUs()")
            << "free-surface velocity field already exists"
            << abort(FatalError);
    }


    UsPtr_ = new areaVectorField
    (
        IOobject
        (
            "Us",
            time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh(),
        dimensioned<vector>("Us", dimVelocity, vector::zero),
        zeroGradientFaPatchVectorField::typeName
    );
}


void freeSurface::makePhis()
{
    if (debug)
    {
        Info<< "freeSurface::makePhis() : "
            << "making free-surface fluid flux"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (phisPtr_)
    {
        FatalErrorIn("freeSurface::makePhis()")
            << "free-surface fluid flux already exists"
            << abort(FatalError);
    }


    phisPtr_ = new edgeScalarField
    (
        IOobject
        (
            "phis",
            time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearEdgeInterpolate(Us()) & aMesh().Le()
    );
}


void freeSurface::makeSurfactConc()
{
    if (debug)
    {
        Info<< "freeSurface::makeSurfactConc() : "
            << "making free-surface surfactant concentration field"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (surfactConcPtr_)
    {
        FatalErrorIn("freeSurface::makeUs()")
            << "free-surface surfactant concentratio field already exists"
            << abort(FatalError);
    }


    surfactConcPtr_ = new areaScalarField
    (
        IOobject
        (
            "Cs",
            time().timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh()
    );
}


void freeSurface::makeSurfaceTension()
{
    if (debug)
    {
        Info<< "freeSurface::makeSurfaceTension() : "
            << "making surface tension field"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (surfaceTensionPtr_)
    {
        FatalErrorIn("freeSurface::makeSurfaceTension()")
            << "surface tension field already exists"
            << abort(FatalError);
    }


    surfaceTensionPtr_ = new areaScalarField
    (
        IOobject
        (
            "surfaceTension",
            time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cleanInterfaceSurfTension()
      + surfactant().surfactR()*
        surfactant().surfactT()*
        surfactant().surfactSaturatedConc()*
        log(1.0 - surfactantConcentration()/
        surfactant().surfactSaturatedConc())
    );
}


void freeSurface::makeSurfactant()
{
    if (debug)
    {
        Info<< "freeSurface::makeSurfactant() : "
            << "making surfactant properties"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (surfactantPtr_)
    {
        FatalErrorIn("freeSurface::makeSurfactant()")
            << "surfactant properties already exists"
            << abort(FatalError);
    }


    const dictionary& surfactProp = 
        this->subDict("surfactantProperties");

    surfactantPtr_ = new surfactantProperties(surfactProp);
}

void freeSurface::makeFluidIndicator()
{
    if (debug)
    {
        Info<< "freeSurface::makeFluidIndicator() : "
            << "making fluid indicator"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (fluidIndicatorPtr_)
    {
        FatalErrorIn("freeSurface::makeFluidIndicator()")
            << "fluid indicator already exists"
            << abort(FatalError);
    }

    fluidIndicatorPtr_ = new volScalarField
    (
        IOobject
        (
            "fluidIndicator",
            time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("1", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField& fluidIndicator = *fluidIndicatorPtr_;

    if (twoFluids())
    {
        // find start cell
        label pointOnShadowPatch =
            mesh().boundaryMesh()[bPatchID()][0][0];

        label startCell = mesh().pointCells()[pointOnShadowPatch][0];


        // get cell-cells addressing
        const labelListList& cellCells = mesh().cellCells();

        SLList<label> slList(startCell);

        while (slList.size())
        {
            label curCell = slList.removeHead();

            if (fluidIndicator[curCell] == 1)
            {
                fluidIndicator[curCell] = 0.0;

                for (int i = 0;i < cellCells[curCell].size();i++)
                {
                    slList.append(cellCells[curCell][i]);
                }
            }
        }
    }

    fluidIndicator.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const tetPolyPatchInterpolation& freeSurface::interpolatorA()
{
    if (!interpolatorAPtr_)
    {
        makePatchPointInterpolators();
    }
    
    return *interpolatorAPtr_;
}


const tetPolyPatchInterpolation& freeSurface::interpolatorB()
{
    if (!interpolatorBPtr_)
    {
        makePatchPointInterpolators();
    }
    
    return *interpolatorBPtr_;
}


const patchToPatchInterpolation& freeSurface::interpolatorAB()
{
    if (!interpolatorABPtr_)
    {
        makeInterpolators();
    }
    
    return *interpolatorABPtr_;
}


const patchToPatchInterpolation& freeSurface::interpolatorBA()
{
    if (!interpolatorBAPtr_)
    {
        makeInterpolators();
    }
    
    return *interpolatorBAPtr_;
}


vectorField& freeSurface::controlPoints()
{
    if (!controlPointsPtr_)
    {
        makeControlPoints();
    }

    return *controlPointsPtr_;
}


scalarField& freeSurface::motionPointsMask()
{
    if (!motionPointsMaskPtr_)
    {
        makeMotionPointsMask();
    }

    return *motionPointsMaskPtr_;
}


vectorField& freeSurface::pointsDisplacementDir()
{
    if (!pointsDisplacementDirPtr_)
    {
        makeDirections();
    }

    return *pointsDisplacementDirPtr_;
}


vectorField& freeSurface::facesDisplacementDir()
{
    if (!facesDisplacementDirPtr_)
    {
        makeDirections();
    }

    return *facesDisplacementDirPtr_;
}


vectorField& freeSurface::totalDisplacement()
{
    if (!totalDisplacementPtr_)
    {
        makeTotalDisplacement();
    }

    return *totalDisplacementPtr_;
}


faMesh& freeSurface::aMesh()
{
    if (!aMeshPtr_)
    {
        makeFaMesh();
    }
    
    return *aMeshPtr_;
}


tetDecompositionMotionSolver& freeSurface::meshMotionSolver()
{
    if (!mSolverPtr_)
    {
        makeMeshMotionSolver();
    }

    return *mSolverPtr_;
}


areaVectorField& freeSurface::Us()
{
    if (!UsPtr_)
    {
        makeUs();
    }
    
    return *UsPtr_;
}


edgeScalarField& freeSurface::Phis()
{
    if (!phisPtr_)
    {
        makePhis();
    }
    
    return *phisPtr_;
}


areaScalarField& freeSurface::surfactantConcentration()
{
    if (!surfactConcPtr_)
    {
        makeSurfactConc();
    }
    
    return *surfactConcPtr_;
}


areaScalarField& freeSurface::surfaceTension()
{
    if (!surfaceTensionPtr_)
    {
        makeSurfaceTension();
    }
    
    return *surfaceTensionPtr_;
}


const surfactantProperties& freeSurface::surfactant()
{
    if (!surfactantPtr_)
    {
        makeSurfactant();
    }
    
    return *surfactantPtr_;
}


const volScalarField& freeSurface::fluidIndicator()
{    
    if (!fluidIndicatorPtr_)
    {
        makeFluidIndicator();
    }

    return *fluidIndicatorPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
