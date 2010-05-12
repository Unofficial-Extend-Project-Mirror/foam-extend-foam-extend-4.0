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

Description

\*---------------------------------------------------------------------------*/

#include "freeSurface.H"
#include "volFields.H"
#include "primitivePatchInterpolation.H"
#include "fixedGradientFvPatchFields.H"
#include "emptyFaPatch.H"
#include "demandDrivenData.H"
#include "EulerDdtScheme.H"
#include "CrankNicholsonDdtScheme.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(freeSurface, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void freeSurface::clearOut()
{
    deleteDemandDrivenData(interpolatorAPtr_);
    deleteDemandDrivenData(interpolatorBPtr_);
    deleteDemandDrivenData(interpolatorABPtr_);
    deleteDemandDrivenData(interpolatorBAPtr_);
    deleteDemandDrivenData(controlPointsPtr_);
    deleteDemandDrivenData(motionPointsMaskPtr_);
    deleteDemandDrivenData(pointsDisplacementDirPtr_);
    deleteDemandDrivenData(facesDisplacementDirPtr_);
    deleteDemandDrivenData(totalDisplacementPtr_);
    deleteDemandDrivenData(aMeshPtr_);
    deleteDemandDrivenData(mSolverPtr_);
    deleteDemandDrivenData(UsPtr_);
    deleteDemandDrivenData(phisPtr_);
    deleteDemandDrivenData(surfactConcPtr_);
    deleteDemandDrivenData(surfaceTensionPtr_);
    deleteDemandDrivenData(surfactantPtr_);
    deleteDemandDrivenData(fluidIndicatorPtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

freeSurface::freeSurface
(
    fvMesh& m,
    const volScalarField& rho,
    volVectorField& Ub, 
    volScalarField& Pb, 
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "freeSurfaceProperties",
            Ub.mesh().time().constant(),
            Ub.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(m),
    rho_(rho),
    U_(Ub),
    p_(Pb),
    phi_(phi),
    curTimeIndex_(Ub.mesh().time().timeIndex()),
    twoFluids_
    (
        this->lookup("twoFluids")
    ),
    normalMotionDir_
    (
        this->lookup("normalMotionDir")
    ),
    freeSurfaceSmoothing_
    (
        this->lookup("freeSurfaceSmoothing")
    ),
    cleanInterface_
    (
        this->lookup("cleanInterface")
    ),
    aPatchID_(-1),
    bPatchID_(-1),
    muFluidA_
    (
        this->lookup("muFluidA")
    ),
    muFluidB_
    (
        this->lookup("muFluidB")
    ),
    rhoFluidA_
    (
        this->lookup("rhoFluidA")
    ),
    rhoFluidB_
    (
        this->lookup("rhoFluidB")
    ),
    g_(this->lookup("g")),
    cleanInterfaceSurfTension_
    (
        this->lookup("surfaceTension")
    ),
    fixedFreeSurfacePatches_
    (
        this->lookup("fixedFreeSurfacePatches")
    ),
    interpolatorAPtr_(NULL),
    interpolatorBPtr_(NULL),
    interpolatorABPtr_(NULL),
    interpolatorBAPtr_(NULL),
    controlPointsPtr_(NULL),
    motionPointsMaskPtr_(NULL),
    pointsDisplacementDirPtr_(NULL),
    facesDisplacementDirPtr_(NULL),
    totalDisplacementPtr_(NULL),
    aMeshPtr_(NULL),
    mSolverPtr_(NULL),
    UsPtr_(NULL),
    phisPtr_(NULL),
    surfactConcPtr_(NULL),
    surfaceTensionPtr_(NULL),
    surfactantPtr_(NULL),
    fluidIndicatorPtr_(NULL)
{
    // Make finite area mesh

    if (!aMeshPtr_)
    {
        makeFaMesh();
    }


    // Make mesh motion solver

    if (!mSolverPtr_)
    {
        makeMeshMotionSolver();
    }


    // Create free-surface velocity field

    if (!UsPtr_)
    {
        makeUs();
    }


    // Create surfactant concentration and fluid flux field

    if(!cleanInterface())
    {
        if (!surfactConcPtr_)
        {
            makeSurfactConc();
        }

        if (!phisPtr_)
        {
            makePhis();
        }
    }


    if(Pstream::master())
    {    
        // Detect the free surface patch

        forAll (mesh().boundary(), patchI)
        {
            if(mesh().boundary()[patchI].name() == "freeSurface")
            {
                aPatchID_ = patchI;
                
                Info<< "Found free surface patch. ID: " << aPatchID_
                    << endl;
            }
        }

        if(aPatchID() == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Free surface patch not defined.  Please make sure that "
                << " the free surface patches is named as freeSurface"
                << abort(FatalError);
        }


        // Detect the free surface shadow patch

        if (twoFluids())
        {
            forAll (mesh().boundary(), patchI)
            {
                if(mesh().boundary()[patchI].name() == "freeSurfaceShadow")
                {
                    bPatchID_ = patchI;
                    
                    Info<< "Found free surface shadow patch. ID: "
                        << bPatchID_ << endl;
                }
            }

            if(bPatchID() == -1)
            {
                FatalErrorIn("freeSurface::freeSurface(...)")
                    << "Free surface shadow patch not defined. "
                    << "Please make sure that the free surface shadow patch "
                    << "is named as freeSurfaceShadow."
                    << abort(FatalError);
            }
        }


        // Define free-surface points displacement directions

        if (!pointsDisplacementDirPtr_)
        {
            makeDirections();
        }


        // Define motion points displacement mask

        if (!motionPointsMaskPtr_)
        {
            makeMotionPointsMask();
        }

        
        // Mark free surface boundary points 
        // which do not belonge to empty patch

        forAll(aMesh().boundary(), patchI)
        {
            if 
            (
                aMesh().boundary()[patchI].type()
             != emptyFaPatch::typeName
            )
            {
                labelList patchPoints =
                    aMesh().boundary()[patchI].pointLabels();

                forAll(patchPoints, pointI)
                {
                    motionPointsMask()[patchPoints[pointI]] = 0.0;
                }                
            }
        }


        // Define control points position

        if (!controlPointsPtr_)
        {
            makeControlPoints();
        }


        // Read free-surface points total displacement if present

        if (!totalDisplacementPtr_)
        {
            readTotalDisplacement();
        }


        // Define variable surface tension calculation

        if (!cleanInterface())
        {
            // Read and create surfactant properties

            if (!surfactantPtr_)
            {
                makeSurfactant();
            }
            

            if (!surfaceTensionPtr_)
            {
                makeSurfaceTension();
            }
        }
    }
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

freeSurface::~freeSurface()
{
    clearOut();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void freeSurface::updateDisplacementDirections()

{
    if(Pstream::master())
    {    
        if(normalMotionDir())
        {
            pointsDisplacementDir() = aMesh().pointAreaNormals();
            
            facesDisplacementDir() =
                aMesh().faceAreaNormals().internalField();


            // Correction of control points postion
            const vectorField& Cf = aMesh().areaCentres().internalField();

            controlPoints() = 
                facesDisplacementDir()*
                (facesDisplacementDir()&(controlPoints() - Cf))
              + Cf;
        }    
    }
}


bool freeSurface::movePoints()
{
    if(Pstream::master())
    {
        scalarField sweptVolCorr = 
            phi().boundaryField()[aPatchID()]
          - fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()];

        word ddtScheme
        (
            mesh().ddtScheme("ddt(" + rho().name() + ',' + U().name()+')')
        );

        if 
        (
            ddtScheme
         == fv::CrankNicholsonDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= (1.0/2.0)*time().deltaT().value();
        }
        else if
        (
            ddtScheme
         == fv::EulerDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= time().deltaT().value();
        }
        else if
        (
            ddtScheme
         == fv::backwardDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= (2.0/3.0)*time().deltaT().value();
        }   
        else
        {
            FatalErrorIn("freeSurface::movePoints()")
                << "Unsupported temporal differencing scheme : "
                << ddtScheme
                << abort(FatalError);
        }


        const scalarField& Sf = aMesh().S();
        const vectorField& Nf = aMesh().faceAreaNormals().internalField();

        scalarField deltaH =
            sweptVolCorr/(Sf*(Nf & facesDisplacementDir()));
        

        pointField displacement = pointDisplacement(deltaH);


        // Move only free-surface points

        pointField newMeshPoints = mesh().points();

        const labelList& meshPointsA = 
            mesh().boundaryMesh()[aPatchID()].meshPoints();

        forAll (displacement, pointI)
        {
            newMeshPoints[meshPointsA[pointI]] += displacement[pointI];
        }

        if(twoFluids_)
        {
            const labelList& meshPointsB = 
                mesh().boundaryMesh()[bPatchID_].meshPoints();

            pointField displacementB =
                interpolatorAB().pointInterpolate
                (
                    displacement
                );
        
            forAll (displacementB, pointI)
            {
                newMeshPoints[meshPointsB[pointI]] += displacementB[pointI]; 
            }
        }

        mesh().movePoints(newMeshPoints);
       
        aMesh().movePoints(mesh().points());


        // Update total displacement field

        if(totalDisplacementPtr_ && (curTimeIndex_ < time().timeIndex()))
        {
            FatalErrorIn("freeSurface::movePoints()")
                << "Total displacement of free surface points "
                << "from previous time step is not absorbed by the mesh."
                << abort(FatalError);
        }
        else if (curTimeIndex_ < time().timeIndex())
        {
            totalDisplacement() = displacement;

            curTimeIndex_ = time().timeIndex();
        }
        else
        {
            totalDisplacement() += displacement;
        }
    }

    return true;
}


bool freeSurface::moveMeshPointsForOldFreeSurfDisplacement()
{
    bool moveMesh = false;

    if(Pstream::master() && totalDisplacementPtr_)
    {
        moveMesh = true;
        Pstream::scatter(moveMesh);

        pointField newPoints = mesh().points();

        const labelList& meshPointsA = 
            mesh().boundaryMesh()[aPatchID()].meshPoints();

        forAll (totalDisplacement(), pointI)
        {
            newPoints[meshPointsA[pointI]] -= totalDisplacement()[pointI]; 
        }

        meshMotionSolver().motionU().boundaryField()[aPatchID()] ==
            interpolatorA().pointToPointInterpolate
            (
                totalDisplacement()/time().deltaT().value()
            );


        if(twoFluids_)
        {
            const labelList& meshPointsB = 
                mesh().boundaryMesh()[bPatchID()].meshPoints();

            pointField totDisplacementB =
                interpolatorAB().pointInterpolate
                (
                    totalDisplacement()
                );
        
            forAll (totDisplacementB, pointI)
            {
                newPoints[meshPointsB[pointI]] -= totDisplacementB[pointI]; 
            }
            
            meshMotionSolver().motionU().boundaryField()[bPatchID()] == 
                interpolatorB().pointToPointInterpolate
                (
                    totDisplacementB/time().deltaT().value()
                );
        }


        mesh().movePoints(newPoints);


        deleteDemandDrivenData(totalDisplacementPtr_);
    }        

    if(moveMesh)
    {
        //-- Update finite volume mesh

        mesh().movePoints(meshMotionSolver().newPoints());
    }


    if(Pstream::master() && moveMesh)
    {
        aMesh().movePoints(mesh().points());

        meshMotionSolver().motionU().boundaryField()[aPatchID()] == 
            vector::zero;

        if(twoFluids_)
        {
            meshMotionSolver().motionU().boundaryField()[bPatchID()] ==
                vector::zero;
        }        
    }

    return true;
}


bool freeSurface::moveMeshPoints()
{
    if(Pstream::master())
    {
        scalarField sweptVolCorr = 
            phi().boundaryField()[aPatchID()]
          - fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()];

        word ddtScheme
        (
            mesh().ddtScheme("ddt(" + rho().name() + ',' + U().name()+')')
        );

        if 
        (
            ddtScheme
         == fv::CrankNicholsonDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= (1.0/2.0)*time().deltaT().value();
        }
        else if
        (
            ddtScheme
         == fv::EulerDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= time().deltaT().value();
        }
        else if
        (
            ddtScheme
         == fv::backwardDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= (2.0/3.0)*time().deltaT().value();
        }   
        else
        {
            FatalErrorIn("freeSurface::movePoints()")
                << "Unsupported temporal differencing scheme : "
                << ddtScheme
                << abort(FatalError);
        }


        const scalarField& Sf = aMesh().S();
        const vectorField& Nf = aMesh().faceAreaNormals().internalField();

        scalarField deltaH =
            sweptVolCorr/(Sf*(Nf & facesDisplacementDir()));


        pointField displacement = pointDisplacement(deltaH);


        //-- Set mesh motion boundary conditions

        meshMotionSolver().motionU().boundaryField()[aPatchID()] ==
            interpolatorA().pointToPointInterpolate
            (
                displacement/time().deltaT().value()
            );


        if (twoFluids())
        {
            meshMotionSolver().motionU().boundaryField()[bPatchID()] ==
                interpolatorB().pointToPointInterpolate
                (
                    interpolatorAB().pointInterpolate
                    (
                        displacement/time().deltaT().value()
                    )
                );
        }
    }


    //-- Update finite volume mesh

    mesh().movePoints(meshMotionSolver().newPoints());


    //-- Update finite area mesh

    if(Pstream::master())
    {
        aMesh().movePoints(mesh().points());


        meshMotionSolver().motionU().boundaryField()[aPatchID()] == 
            vector::zero;

        if(twoFluids_)
        {
            meshMotionSolver().motionU().boundaryField()[bPatchID()] ==
                vector::zero;
        }        
    }


    return true;
}


bool freeSurface::moveMeshPoints(const scalarField& delta)
{
    if(Pstream::master())
    {
        pointField displacement = delta*pointsDisplacementDir();


        //-- Set mesh motion boundary conditions

        meshMotionSolver().motionU().boundaryField()[aPatchID()] ==
            interpolatorA().pointToPointInterpolate
            (
                displacement/time().deltaT().value()
            );


        if (twoFluids())
        {
            meshMotionSolver().motionU().boundaryField()[bPatchID()] ==
                interpolatorB().pointToPointInterpolate
                (
                    interpolatorAB().pointInterpolate
                    (
                        displacement/time().deltaT().value()
                    )
                );
        }
    }


    //-- Update finite volume mesh

    mesh().movePoints(meshMotionSolver().newPoints());


    //-- Update finite area mesh

    if(Pstream::master())
    {
        aMesh().movePoints(mesh().points());


        meshMotionSolver().motionU().boundaryField()[aPatchID()] == 
            vector::zero;

        if(twoFluids_)
        {
            meshMotionSolver().motionU().boundaryField()[bPatchID()] ==
                vector::zero;
        }        
    }


    return true;
}


bool freeSurface::movePoints(scalar displ)
{
    if(Pstream::master())
    {
        controlPoints() += displ*facesDisplacementDir();

        pointField displacement = displ*pointsDisplacementDir();


        // Move only free-surface points

        pointField newMeshPoints = mesh().points();

        const labelList& meshPointsA = 
            mesh().boundaryMesh()[aPatchID()].meshPoints();

        forAll (displacement, pointI)
        {
            newMeshPoints[meshPointsA[pointI]] += displacement[pointI]; 
        }

        if(twoFluids())
        {
            const labelList& meshPointsB = 
                mesh().boundaryMesh()[bPatchID_].meshPoints();

            pointField displacementB =
                interpolatorAB().pointInterpolate
                (
                    displacement
                );
        
            forAll (displacementB, pointI)
            {
                newMeshPoints[meshPointsB[pointI]] += displacementB[pointI]; 
            }
        }

        mesh().movePoints(newMeshPoints);
       
        aMesh().movePoints(mesh().points());


        // Update total displacement field

        if(totalDisplacementPtr_ && (curTimeIndex_ < time().timeIndex()))
        {
            FatalErrorIn("freeSurface::movePoints()")
                << "Total displacement of free surface points "
                << "from previous time step is not absorbed by the mesh."
                << abort(FatalError);
        }
        else if (curTimeIndex_ < time().timeIndex())
        {
            totalDisplacement() = displacement;
            
            curTimeIndex_ = time().timeIndex();
        }
        else
        {
            totalDisplacement() += displacement;
        }
    }

    return true;
}


void freeSurface::correctBoundaryConditions()
{
    correctVelocity();
    correctVelocityGradient();
    correctSurfactantConcentration();
    correctPressure();
}


void freeSurface::correctVelocity()
{
    if(Pstream::master())
    {
        // Correct velocity boundary condition at the free-surface

         vectorField nA = mesh().boundary()[aPatchID()].nf();

        if(twoFluids())
        {
            vectorField nB = mesh().boundary()[bPatchID()].nf();

            // Calculate normal component of free surface velocity

            vectorField UnPB = interpolatorBA().faceInterpolate
            (
                U().boundaryField()[bPatchID()].patchInternalField()
            );
            UnPB = nA*(nA & UnPB);

            vectorField UnPA = 
                nA*(nA & U().boundaryField()[aPatchID()].patchInternalField());

            scalarField DnB = 
                interpolatorBA().faceInterpolate
                (
                    mesh().boundary()[bPatchID()].deltaCoeffs()
                );

            scalarField DnA = mesh().boundary()[bPatchID()].deltaCoeffs();

            Us().internalField() = (DnA*UnPA + DnB*UnPB)/(DnA + DnB);
            Us().correctBoundaryConditions();



            // Calculate the tangential component of free surface velocity

            vectorField UtPA = 
                U().boundaryField()[aPatchID()].patchInternalField();
            UtPA -= nA*(nA & UtPA);

            vectorField UtPB = interpolatorBA().faceInterpolate
            (
                U().boundaryField()[bPatchID()].patchInternalField()
            );
            UtPB -= nA*(nA & UtPB);

            DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

            DnB = interpolatorBA().faceInterpolate
            (
                mesh().boundary()[bPatchID()].deltaCoeffs()
            );

            vectorField sGradVn =
                fac::grad(Us()&aMesh().faceAreaNormals())
                ().internalField();

            vectorField UtFs = 
                muFluidA().value()*DnA*UtPA
              + muFluidB().value()*DnB*UtPB
              - (muFluidA().value() - muFluidB().value())*sGradVn;

            if(!cleanInterface())
            {
                UtFs += fac::grad(surfaceTension())().internalField();
            }

            UtFs /= muFluidA().value()*DnA + muFluidB().value()*DnB + SMALL;
                

            Us().internalField() += UtFs;
            Us().correctBoundaryConditions();


            U().boundaryField()[bPatchID()] == 
                interpolatorAB().faceInterpolate(UtFs)
              + nB*fvc::meshPhi(rho(),U())().boundaryField()[bPatchID()]/
                mesh().boundary()[bPatchID()].magSf();
        }
        else
        {
            // Calculate normal component of free surface velocity

            Us().internalField() = 
                nA*phi().boundaryField()[aPatchID()]/
                mesh().boundary()[aPatchID()].magSf();
            Us().correctBoundaryConditions();


            // Calculate the tangential component of free surface velocity

            vectorField UtPA = 
                U().boundaryField()[aPatchID()].patchInternalField();
            UtPA -= nA*(nA & UtPA);

            scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

            vectorField sGradVn =
                fac::grad(Us()&aMesh().faceAreaNormals())
                ().internalField();

            vectorField UtFs = 
                muFluidA().value()*DnA*UtPA
              - muFluidA().value()*sGradVn;

            if(!cleanInterface())
            {
                UtFs += fac::grad(surfaceTension())().internalField();
            }

            UtFs /= (muFluidA().value()*DnA + SMALL);
                

            Us().internalField() += UtFs;
            Us().correctBoundaryConditions();


//             Us().internalField() = U().boundaryField()[aPatchID()];
//             Us().correctBoundaryConditions();
        }
    }
}


void freeSurface::correctVelocityGradient()
{
    if(Pstream::master())
    {
        // Correct normal velocity gradient at the free-surface

        vectorField nA = mesh().boundary()[aPatchID()].nf();

        if(twoFluids())
        {
            vectorField VtFs = Us().internalField();
            VtFs -= nA*(nA & VtFs);

            vectorField VtPA =
                U().boundaryField()[aPatchID()].patchInternalField();
            VtPA -= nA*(nA & VtPA);
            

            vectorField nGradU = 
                (VtFs - VtPA)*mesh().boundary()[aPatchID()].deltaCoeffs();


            vectorField VnFs = Us().internalField();
            VnFs = nA*(nA & VnFs);

            vectorField VnPA =
                U().boundaryField()[aPatchID()].patchInternalField();
            VnPA = nA*(nA & VnPA);


            nGradU += 
                (VnFs - VnPA)*mesh().boundary()[aPatchID()].deltaCoeffs();


//             scalarField nnGradU = - fac::div(Us())().internalField();
            
//             nGradU += nnGradU*nA;


            
//             vectorField nGradU =
//               - muFluidB().value()/(muFluidA().value() + SMALL)*
//                 (
//                     (I - nA*nA)&
//                     interpolatorBA().faceInterpolate
//                     (
//                         U().boundaryField()[bPatchID()].snGrad()
//                     )
//                 );


//             if(true)
//             {
//                 nGradU -= (muFluidA().value() - muFluidB().value())/
//                     (muFluidA().value() + SMALL)*
//                     fac::grad(Us()&aMesh().faceAreaNormals())
//                     ().internalField();
//             }


//             if(normalShearStress())
//             {
//                 scalarField nnGradU = - fac::div(Us())().internalField();

//                 nGradU += 
//                     nnGradU*aMesh().faceAreaNormals().internalField();
//             }


//             if(!cleanInterface())
//             {
//                 nGradU += fac::grad(*surfaceTensionPtr_)().internalField()/
//                     muFluidA().value();
//             }


            if
            (
                U().boundaryField()[aPatchID()].type()
             == fixedGradientFvPatchField<vector>::typeName
            )
            {
                fixedGradientFvPatchField<vector>& aU =
                    refCast<fixedGradientFvPatchField<vector> >
                    (
                        U().boundaryField()[aPatchID()]
                    );

                aU.gradient() = nGradU;
            }
            else
            {
                Info << "Bounary condition on " << U().name() 
                    <<  " for freeSurface patch is " 
                    << U().boundaryField()[aPatchID()].type() 
                    << ", instead " 
                    << fixedGradientFvPatchField<vector>::typeName 
                    << '.' << endl;
            }
        }
        else
        {
            vectorField VtFs = Us().internalField();
            VtFs -= nA*(nA & VtFs);

            vectorField VtPA =
                U().boundaryField()[aPatchID()].patchInternalField();
            VtPA -= nA*(nA & VtPA);
            

            vectorField nGradU = 
                (VtFs - VtPA)*mesh().boundary()[aPatchID()].deltaCoeffs();


            vectorField VnFs = Us().internalField();
            VnFs = nA*(nA & VnFs);

            vectorField VnPA =
                U().boundaryField()[aPatchID()].patchInternalField();
            VnPA = nA*(nA & VnPA);


            nGradU += 
                (VnFs - VnPA)*mesh().boundary()[aPatchID()].deltaCoeffs();


//             scalarField nnGradU = - fac::div(Us())().internalField();
            
//             nGradU += nnGradU*nA;



//             vectorField nGradU
//             (
//                 mesh().boundaryMesh()[aPatchID()].size(),
//                 vector::zero
//             );

//             if(true)
//             {
//                 nGradU -=
//                     fac::grad(Us()&aMesh().faceAreaNormals())
//                     ().internalField();
//             }

//             if(true)
//             {
//                 scalarField nnGradU = - fac::div(Us())().internalField();

//                 nGradU += 
//                     nnGradU*aMesh().faceAreaNormals().internalField();
//             }


//             if(!cleanInterface())
//             {
//                 nGradU += fac::grad(*surfaceTensionPtr_)().internalField()/
//                     muFluidA().value();
//             }


            if
            (
                U().boundaryField()[aPatchID()].type()
             == fixedGradientFvPatchField<vector>::typeName
            )
            {
                fixedGradientFvPatchField<vector>& aU =
                    refCast<fixedGradientFvPatchField<vector> >
                    (
                        U().boundaryField()[aPatchID()]
                    );

                aU.gradient() = nGradU;
            }
            else
            {
                Info << "Bounary condition on " << U().name() 
                    <<  " for freeSurface patch is " 
                    << U().boundaryField()[aPatchID()].type() 
                    << ", instead " 
                    << fixedGradientFvPatchField<vector>::typeName 
                    << '.'<< endl;
            }
        }    
    }
}


void freeSurface::correctPressure()
{
    if(Pstream::master())
    {
        // Correct pressure boundary condition at the free-surface
        
        vectorField nA = mesh().boundary()[aPatchID()].nf();

        if(twoFluids())
        {
            scalarField pA =
                interpolatorBA().faceInterpolate
                (
                    p().boundaryField()[bPatchID()]
                );


            const scalarField& K = aMesh().faceCurvatures().internalField();
            
            Info << "Free surface curvature: min = " << min(K)
                << ", max = " << max(K)
                << ", average = " << average(K) << endl << flush;


            if(cleanInterface())
            {
                pA -= cleanInterfaceSurfTension().value()*(K - average(K));
            }
            else
            {
                scalarField surfTensionK =
                    surfaceTension().internalField()*K;
                
                pA -= surfTensionK - average(surfTensionK);
            }


            if(true)
            {
//                 scalarField nnGradU = - fac::div(Us())().internalField();

                scalarField nnGradU = 
                    nA&U().boundaryField()[aPatchID()].snGrad();

                pA += 2.0*(muFluidA().value() - muFluidB().value())
                    *nnGradU;
            }
            

            pA -= (rhoFluidA().value() - rhoFluidB().value())*
                (
                    g_.value()
                    & (
                        (mesh().C().boundaryField()[aPatchID()])
                      - average(mesh().C().boundaryField()[aPatchID()])
                    )
                );


            p().boundaryField()[aPatchID()] == pA;
        }
        else
        {
            vector R0 = average(mesh().C().boundaryField()[aPatchID()]);

            scalarField pA =
              - rhoFluidA().value()*
                (
                    g_.value()
                    & (
                        mesh().C().boundaryField()[aPatchID()]
                      - R0
                    )
                );


            const scalarField& K = aMesh().faceCurvatures().internalField();
        
            Info << "Free surface curvature: min = " << min(K)
                << ", max = " << max(K) << ", average = " << average(K) 
                << endl;


            if(cleanInterface())
            {
                pA -= cleanInterfaceSurfTension().value()*(K - average(K));
            }
            else
            {
                scalarField surfTensionK =
                    surfaceTension().internalField()*K;
                
                pA -= surfTensionK - average(surfTensionK);
            }


            if(true)
            {
//                 scalarField nnGradU = - fac::div(Us())().internalField();

                scalarField nnGradU = 
                    nA&U().boundaryField()[aPatchID()].snGrad();

                pA += 2.0*muFluidA().value()*nnGradU;
            }


            p().boundaryField()[aPatchID()] == pA;
        }
    }

    // Set modified pressure at patches with fixed apsolute
    // pressure

    vector R0 = average(mesh().C().boundaryField()[aPatchID()]);

    for (int patchI=0; patchI < p().boundaryField().size(); patchI++)
    {
        if 
        (
            p().boundaryField()[patchI].type()
         == fixedValueFvPatchScalarField::typeName
        )
        {
            if (patchI != aPatchID())
            {
                p().boundaryField()[patchI] == 
                  - rho().boundaryField()[patchI]*
                    (g_.value()&(mesh().C().boundaryField()[patchI] - R0));
            }
        }
    }
}


void freeSurface::correctSurfaceFlux()
{
    if(Pstream::master())
    {
        Phis() = fac::interpolate(Us()) & aMesh().Le();
    }
}


void freeSurface::correctSurfactantConcentration()
{
    if(!cleanInterface())
    {
        Info << "Correct surfactant concentration" << endl << flush;
        
        correctSurfaceFlux();

        faScalarMatrix CsEqn
        (
            fam::ddt(surfactantConcentration())
          + fam::div(Phis(), surfactantConcentration())
          - fam::laplacian
            (
                surfactant().surfactDiffusion(), 
                surfactantConcentration()
            )
        );

        CsEqn.solve();

        if(Pstream::master())
        {
            Info << "Correct surface tension" << endl;

            surfaceTension() =
                cleanInterfaceSurfTension()
              + surfactant().surfactR()*
                surfactant().surfactT()*
                surfactant().surfactSaturatedConc()*
                log(1.0 - surfactantConcentration()/
                surfactant().surfactSaturatedConc());

            if(neg(min(surfaceTension().internalField())))
            {
                FatalErrorIn
                (
                    "void freeSurface::correctSurfactantConcentration()"
                ) 
                << "Surface tension has negative value" 
                << abort(FatalError);
            }
        }
    }
}


vector freeSurface::totalPressureForce()
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& P = p().boundaryField()[aPatchID()];

    vectorField pressureForces = S*P*n;

    return sum(pressureForces);
}


vector freeSurface::totalViscousForce()
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    vectorField snGradU =
        U_.boundaryField()[aPatchID()].snGrad();

    vectorField viscousForces = ((I-n*n)&snGradU);

    if(true)
    {
        viscousForces +=
            (fac::grad(Us())&aMesh().faceAreaNormals())().internalField();

        viscousForces -=
            (2.0*aMesh().faceAreaNormals()*fac::div(Us()))
            ().internalField();
    }

    viscousForces *= -muFluidA().value()*S;

    return sum(viscousForces);
}


tmp<scalarField> freeSurface::undulationIndicator()
{
    tmp<scalarField> tUndulation
    (
        new scalarField
        (
            aMesh().nFaces(),
            0.0
        )
    );

    scalarField& undulation = tUndulation();

    primitivePatchInterpolation patchInterpolator
        (
            mesh().boundaryMesh()[aPatchID()]
        );

    undulation = 
        asin(mag(
            aMesh().faceAreaNormals().internalField()^
            patchInterpolator.pointToFaceInterpolate
            (
                aMesh().pointAreaNormals()
            )
        ))*180.0/M_PI;

    return tUndulation;
}


void freeSurface::smooth()
{
    if
    (
        Pstream::master()
     && !meshMotionSolver().twoDMotion()
     && freeSurfaceSmoothing()
    )
    {
        Info << "Smooting free surface ...";

        scalarField deltaH = facesDisplacementDir()&
        (
            aMesh().areaCentres().internalField()
          - controlPoints()
        );

        pointField displacement = pointDisplacement(deltaH);


        const faceList& faces = aMesh().faces();
        const pointField& points = aMesh().points();


        pointField newPoints = points + displacement;

        scalarField sweptVol(faces.size(), 0.0);

        forAll(faces, faceI)
        {
            sweptVol[faceI] = -faces[faceI].sweptVol(points, newPoints);
        }

        vectorField faceArea(faces.size(), vector::zero);

        forAll (faceArea, faceI)
        {
            faceArea[faceI] = faces[faceI].normal(newPoints);
        }

        forAll(deltaH, faceI)
        {
            deltaH[faceI] = sweptVol[faceI]/
                (faceArea[faceI] & facesDisplacementDir()[faceI]);
        }


        displacement = pointDisplacement(deltaH);


        // Move only free-surface points

        pointField newMeshPoints = mesh().points();

        const labelList& meshPointsA = 
            mesh().boundaryMesh()[aPatchID()].meshPoints();

        forAll (displacement, pointI)
        {
            newMeshPoints[meshPointsA[pointI]] += displacement[pointI]; 
        }

        if(twoFluids_)
        {
            const labelList& meshPointsB = 
                mesh().boundaryMesh()[bPatchID_].meshPoints();
            
            pointField displacementB =
                interpolatorAB().pointInterpolate
                (
                    displacement
                );
        
            forAll (displacementB, pointI)
            {
                newMeshPoints[meshPointsB[pointI]] += displacementB[pointI]; 
            }
        }

        mesh().movePoints(newMeshPoints);
       
        aMesh().movePoints(mesh().points());


        // Update total displacement field

        if(totalDisplacementPtr_ && (curTimeIndex_ < time().timeIndex()))
        {
            FatalErrorIn("freeSurface::movePoints()")
                << "Total displacement of free surface points "
                << "from previous time step is not absorbed by the mesh."
                << abort(FatalError);
        }
        else if (curTimeIndex_ < time().timeIndex())
        {
            totalDisplacement() = displacement;
        
            curTimeIndex_ = time().timeIndex();
        }
        else
        {
            totalDisplacement() += displacement;
        }

        Info << "done" << endl << flush;
    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
