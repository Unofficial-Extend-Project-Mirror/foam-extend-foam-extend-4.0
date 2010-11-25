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
#include "transformField.H"

#include "emptyFaPatch.H"
#include "wedgeFaPatch.H"
#include "wallFvPatch.H"

#include "EulerDdtScheme.H"
#include "CrankNicholsonDdtScheme.H"
#include "backwardDdtScheme.H"

#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "fixedValueTetPolyPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "twoDPointCorrector.H"

#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "zeroGradientCorrectedFvPatchFields.H"
#include "fixedGradientCorrectedFvPatchFields.H"
#include "fixedValueCorrectedFvPatchFields.H"

#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(freeSurface, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void freeSurface::clearOut()
{
    deleteDemandDrivenData(interpolatorABPtr_);
    deleteDemandDrivenData(interpolatorBAPtr_);
    deleteDemandDrivenData(controlPointsPtr_);
    deleteDemandDrivenData(motionPointsMaskPtr_);
    deleteDemandDrivenData(pointsDisplacementDirPtr_);
    deleteDemandDrivenData(facesDisplacementDirPtr_);
    deleteDemandDrivenData(totalDisplacementPtr_);
    deleteDemandDrivenData(aMeshPtr_);
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
    dynamicFvMesh& m,
    const volScalarField& rho,
    volVectorField& Ub, 
    volScalarField& Pb, 
    const surfaceScalarField& sfPhi
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
    phi_(sfPhi),
    curTimeIndex_(Ub.mesh().time().timeIndex()),
    twoFluids_
    (
        this->lookup("twoFluids")
    ),
    normalMotionDir_
    (
        this->lookup("normalMotionDir")
    ),
    motionDir_(0, 0, 0),
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
    pointNormalsCorrectionPatches_
    (
        this->lookup("pointNormalsCorrectionPatches")
    ),
    nFreeSurfCorr_
    (
        readInt(this->lookup("nFreeSurfaceCorrectors"))
    ),
    smoothing_(false),
    interpolatorABPtr_(NULL),
    interpolatorBAPtr_(NULL),
    controlPointsPtr_(NULL),
    motionPointsMaskPtr_(NULL),
    pointsDisplacementDirPtr_(NULL),
    facesDisplacementDirPtr_(NULL),
    totalDisplacementPtr_(NULL),
    aMeshPtr_(NULL),
    UsPtr_(NULL),
    phisPtr_(NULL),
    surfactConcPtr_(NULL),
    surfaceTensionPtr_(NULL),
    surfactantPtr_(NULL),
    fluidIndicatorPtr_(NULL)
{
    //Read motion direction
    if (!normalMotionDir_)
    {
        motionDir_ = vector(this->lookup("motionDir"));
        motionDir_ /= mag(motionDir_) + SMALL;
    }

    // Set point normal correction patches
    boolList& correction = aMesh().correctPatchPointNormals();

    forAll(pointNormalsCorrectionPatches_, patchI)
    {
        word patchName = pointNormalsCorrectionPatches_[patchI];

        label patchID = aMesh().boundary().findPatchID(patchName);

        if(patchID == -1)
        {
            FatalErrorIn
            (
                "freeSurface::freeSurface(...)"
            )   << "Patch name for point normals correction does not exist"
                << abort(FatalError);
        }

        correction[patchID] = true;
    }

    // Clear geometry
    aMesh().movePoints(aMesh().points());


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


    // Mark free surface boundary points 
    // which belonge to processor patches
    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == processorFaPatch::typeName
        )
        {
            const labelList& patchPoints =
                aMesh().boundary()[patchI].pointLabels();

            forAll(patchPoints, pointI)
            {
                motionPointsMask()[patchPoints[pointI]] = -1;
            }                
        }
    }


    // Mark fixed free surface boundary points 
    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID = 
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& patchPoints =
            aMesh().boundary()[fixedPatchID].pointLabels();

        forAll(patchPoints, pointI)
        {
            motionPointsMask()[patchPoints[pointI]] = 0;
        }
    }


    // Mark free-surface boundary point 
    // at the axis of 2-D axisymmetic cases
    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == wedgeFaPatch::typeName
        )
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

            if(wedgePatch.axisPoint() > -1)
            {
                motionPointsMask()[wedgePatch.axisPoint()] = 0;

                Info << "Axis point: " 
                    << wedgePatch.axisPoint()
                    << "vector: " 
                    << aMesh().points()[wedgePatch.axisPoint()] << endl;
            }
        }
    }


    // Read free-surface points total displacement if present
    readTotalDisplacement();


    // Read control points positions if present
    controlPoints();


    // Check if smoothing switch is set
    if (this->found("smoothing"))
    {
        smoothing_ = Switch(this->lookup("smoothing"));
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
    if(normalMotionDir())
    {
        // Update point displacement correction
        pointsDisplacementDir() = aMesh().pointAreaNormals();

        // Correcte point displacement direction 
        // at the "centerline" symmetryPlane which represents the axis
        // of an axisymmetric case
        forAll(aMesh().boundary(), patchI)
        {
            if(aMesh().boundary()[patchI].type() == wedgeFaPatch::typeName)
            {
                const wedgeFaPatch& wedgePatch =
                    refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);
		
                vector axis = wedgePatch.axis();

                label centerLinePatchID = 
                    aMesh().boundary().findPatchID("centerline");
	
                if(centerLinePatchID != -1)
                {
                    const labelList& pointLabels = 
                        aMesh().boundary()[centerLinePatchID].pointLabels();
                    
                    forAll(pointLabels, pointI)
                    {
                        vector dir = 
                            pointsDisplacementDir()[pointLabels[pointI]];
                    
                        dir = (dir&axis)*axis;
                        dir /= mag(dir);
                
                        pointsDisplacementDir()[pointLabels[pointI]] = dir;
                    }
                }
                else
                {
                    Info << "Warning: centerline polyPatch does not exist. " 
                        << "Free surface points displacement directions "
                        << "will not be corrected at the axis (centerline)" 
                        << endl; 
                }
            
                break;   
            }
        }

        // Update face displacement direction
        facesDisplacementDir() =
            aMesh().faceAreaNormals().internalField();

        // Correction of control points postion
        const vectorField& Cf = aMesh().areaCentres().internalField();

        controlPoints() = 
            facesDisplacementDir()
           *(facesDisplacementDir()&(controlPoints() - Cf))
          + Cf;
    }
}


bool freeSurface::predictPoints()
{
    // Smooth interface

    if (smoothing_)
    {
        controlPoints() = aMesh().areaCentres().internalField();
        movePoints(scalarField(controlPoints().size(), 0));
        movePoints(-fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()]);
    }

    for
    (
        int freeSurfCorr=0;
        freeSurfCorr<nFreeSurfCorr_;
        freeSurfCorr++
    )
    {
        movePoints(phi_.boundaryField()[aPatchID()]);
    }

    return true;
}


bool freeSurface::correctPoints()
{
    for
    (
        int freeSurfCorr=0;
        freeSurfCorr<nFreeSurfCorr_;
        freeSurfCorr++
    )
    {
        movePoints(phi_.boundaryField()[aPatchID()]);
    }

    return true;
}


bool freeSurface::movePoints(const scalarField& interfacePhi)
{
    pointField newMeshPoints = mesh().points();

    scalarField sweptVolCorr = 
        interfacePhi
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
        sweptVolCorr *= (1.0/2.0)*DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::EulerDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::backwardDdtScheme<vector>::typeName
    )
    {
        if (DB().timeIndex() == 1)
        {
            sweptVolCorr *= DB().deltaT().value();
        }
        else
        {
            sweptVolCorr *= (2.0/3.0)*DB().deltaT().value();
        }
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

    // Update total displacement field

    if(totalDisplacementPtr_ && (curTimeIndex_ < DB().timeIndex()))
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Total displacement of free surface points "
                << "from previous time step is not absorbed by the mesh."
                << abort(FatalError);
    }
    else if (curTimeIndex_ < DB().timeIndex())
    {
        totalDisplacement() = displacement;

        curTimeIndex_ = DB().timeIndex();
    }
    else
    {
        totalDisplacement() += displacement;
    }

    twoDPointCorrector twoDPointCorr(mesh());

    twoDPointCorr.correctPoints(newMeshPoints);

    mesh().movePoints(newMeshPoints);

    aMesh().movePoints(mesh().points());


    // Move correctedFvPatchField fvSubMeshes

    forAll(U().boundaryField(), patchI)
    {
        if
        (
            (
                U().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<vector>::typeName
            )
        )
        {
            correctedFvPatchField<vector>& pU =
                refCast<correctedFvPatchField<vector> >
                (
                    U().boundaryField()[patchI]
                );

            pU.movePatchSubMesh();
        }
    }

    forAll(p().boundaryField(), patchI)
    {
        if
        (
            (
                p().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<scalar>::typeName
            )
        )
        {
            correctedFvPatchField<scalar>& pP =
                refCast<correctedFvPatchField<scalar> >
                (
                    p().boundaryField()[patchI]
                );

            pP.movePatchSubMesh();
        }
    }

    return true;
}


bool freeSurface::moveMeshPointsForOldFreeSurfDisplacement()
{
    if(totalDisplacementPtr_)
    {
        pointField newPoints = mesh().points();

        const labelList& meshPointsA = 
            mesh().boundaryMesh()[aPatchID()].meshPoints();

        forAll (totalDisplacement(), pointI)
        {
            newPoints[meshPointsA[pointI]] -= totalDisplacement()[pointI]; 
        }


        // Check mesh motion solver type 
        bool feMotionSolver = 
            mesh().objectRegistry::foundObject<tetPointVectorField>
            (
                "motionU"
            );
        bool fvMotionSolver =
            mesh().objectRegistry::foundObject<pointVectorField>
            (
                "pointMotionU"
            );

        if (feMotionSolver)
        {
            tetPointVectorField& motionU =
                const_cast<tetPointVectorField&>
                (
                    mesh().objectRegistry::
                    lookupObject<tetPointVectorField>
                    (
                        "motionU"
                    )
                );

            fixedValueTetPolyPatchVectorField& motionUaPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[aPatchID()]
                );

            tetPolyPatchInterpolation tppiAPatch
            (
                refCast<const faceTetPolyPatch>
                (
                    motionUaPatch.patch()
                )
            );

            motionUaPatch ==
                tppiAPatch.pointToPointInterpolate
                (
                    totalDisplacement()/DB().deltaT().value()
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
                    newPoints[meshPointsB[pointI]] -= 
                        totDisplacementB[pointI]; 
                }
            
                fixedValueTetPolyPatchVectorField& motionUbPatch =
                    refCast<fixedValueTetPolyPatchVectorField>
                    (
                        motionU.boundaryField()[bPatchID()]
                    );

                tetPolyPatchInterpolation tppiBPatch
                (
                    refCast<const faceTetPolyPatch>(motionUbPatch.patch())
                );

                motionUbPatch == 
                    tppiBPatch.pointToPointInterpolate
                    (
                        totDisplacementB/DB().deltaT().value()
                    );
            }
        }
        else if (fvMotionSolver)
        {
            pointVectorField& motionU =
                const_cast<pointVectorField&>
                (
                    mesh().objectRegistry::
                    lookupObject<pointVectorField>
                    (
                        "pointMotionU"
                    )
                );

            fixedValuePointPatchVectorField& motionUaPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    motionU.boundaryField()[aPatchID()]
                );

            motionUaPatch ==
                totalDisplacement()/DB().deltaT().value();

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
                    newPoints[meshPointsB[pointI]] -= 
                        totDisplacementB[pointI]; 
                }
            
                fixedValuePointPatchVectorField& motionUbPatch =
                    refCast<fixedValuePointPatchVectorField>
                    (
                        motionU.boundaryField()[bPatchID()]
                    );

                motionUbPatch == 
                    totDisplacementB/DB().deltaT().value();
            }
        }

        twoDPointCorrector twoDPointCorr(mesh());

        twoDPointCorr.correctPoints(newPoints);

        mesh().movePoints(newPoints);

        deleteDemandDrivenData(totalDisplacementPtr_);

        mesh().update();

        aMesh().movePoints(mesh().points());

        // Move correctedFvPatchField fvSubMeshes

        forAll(U().boundaryField(), patchI)
        {
            if
            (
                (
                    U().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<vector>::typeName
                )
                ||
                (
                    U().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<vector>::typeName
                )
                ||
                (
                    U().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<vector>::typeName
                )
            )
            {
                correctedFvPatchField<vector>& aU =
                    refCast<correctedFvPatchField<vector> >
                    (
                        U().boundaryField()[patchI]
                    );

                aU.movePatchSubMesh();
            }
        }

        forAll(p().boundaryField(), patchI)
        {
            if
            (
                (
                    p().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<scalar>::typeName
                )
                ||
                (
                    p().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<scalar>::typeName
                )
                ||
                (
                    p().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<scalar>::typeName
                )
            )
            {
                correctedFvPatchField<scalar>& aP =
                    refCast<correctedFvPatchField<scalar> >
                    (
                        p().boundaryField()[patchI]
                    );

                aP.movePatchSubMesh();
            }
        }
    }

    return true;
}


bool freeSurface::moveMeshPoints()
{
        scalarField sweptVolCorr = 
            phi_.boundaryField()[aPatchID()]
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
            sweptVolCorr *= (1.0/2.0)*DB().deltaT().value();
        }
        else if
        (
            ddtScheme
         == fv::EulerDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= DB().deltaT().value();
        }
        else if
        (
            ddtScheme
         == fv::backwardDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= (2.0/3.0)*DB().deltaT().value();
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

        tetPointVectorField& motionU =
            const_cast<tetPointVectorField&>
            (
                mesh().objectRegistry::
                lookupObject<tetPointVectorField>
                (
                    "motionU"
                )
            );

        fixedValueTetPolyPatchVectorField& motionUaPatch =
            refCast<fixedValueTetPolyPatchVectorField>
            (
                motionU.boundaryField()[aPatchID()]
            );

        tetPolyPatchInterpolation tppiAPatch
        (
            refCast<const faceTetPolyPatch>
            (
                motionUaPatch.patch()
            )
        );

        motionUaPatch ==
            tppiAPatch.pointToPointInterpolate
            (
                displacement/DB().deltaT().value()
            );

        if (twoFluids())
        {
            fixedValueTetPolyPatchVectorField& motionUbPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[bPatchID()]
                );

            tetPolyPatchInterpolation tppiBPatch
            (
                refCast<const faceTetPolyPatch>(motionUbPatch.patch())
            );

            motionUbPatch == 
                tppiBPatch.pointToPointInterpolate
                (
                    interpolatorAB().pointInterpolate
                    (
                        displacement/DB().deltaT().value()
                    )
                );
        }

        mesh().update();

        aMesh().movePoints(mesh().points());


        // Move correctedFvPatchField fvSubMeshes

        forAll(U().boundaryField(), patchI)
        {
            if
            (
                (
                    U().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<vector>::typeName
                )
                ||
                (
                    U().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<vector>::typeName
                )
                ||
                (
                    U().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<vector>::typeName
                )
            )
            {
                correctedFvPatchField<vector>& aU =
                    refCast<correctedFvPatchField<vector> >
                    (
                        U().boundaryField()[patchI]
                    );

                aU.movePatchSubMesh();
            }
        }

        forAll(p().boundaryField(), patchI)
        {
            if
            (
                (
                    p().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<scalar>::typeName
                )
                ||
                (
                    p().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<scalar>::typeName
                )
                ||
                (
                    p().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<scalar>::typeName
                )
            )
            {
                correctedFvPatchField<scalar>& aP =
                    refCast<correctedFvPatchField<scalar> >
                    (
                        p().boundaryField()[patchI]
                    );

                aP.movePatchSubMesh();
            }
        }

    return true;
}


void freeSurface::updateBoundaryConditions()
{
    updateVelocity();
    updateSurfactantConcentration();
    updatePressure();
}


void freeSurface::updateVelocity()
{
    if(twoFluids())
    {
        vectorField nA = mesh().boundary()[aPatchID()].nf();

        vectorField nB = mesh().boundary()[bPatchID()].nf();

        scalarField DnB = interpolatorBA().faceInterpolate
        (
            mesh().boundary()[bPatchID()].deltaCoeffs()
        );

        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();


        vectorField UtPA = 
            U().boundaryField()[aPatchID()].patchInternalField();

        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            UtPA += aU.corrVecGrad();
        }

        UtPA -= nA*(nA & UtPA);


        vectorField UtPB = interpolatorBA().faceInterpolate
        (
            U().boundaryField()[bPatchID()].patchInternalField()
        );

        if
        (
            U().boundaryField()[bPatchID()].type()
         == fixedValueCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedValueCorrectedFvPatchField<vector>& bU =
                refCast<fixedValueCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[bPatchID()]
                );

            UtPB += interpolatorBA().faceInterpolate(bU.corrVecGrad());
        }

        UtPB -= nA*(nA & UtPB);

        vectorField UtFs = 
            muFluidA().value()*DnA*UtPA 
          + muFluidB().value()*DnB*UtPB;

        vectorField UnFs = 
            nA*phi_.boundaryField()[aPatchID()]
           /mesh().boundary()[aPatchID()].magSf();

        Us().internalField() += UnFs - nA*(nA&Us().internalField());
        correctUsBoundaryConditions();

        UtFs -= (muFluidA().value() - muFluidB().value())*
            (fac::grad(Us())&aMesh().faceAreaNormals())().internalField();


        vectorField tangentialSurfaceTensionForce(nA.size(), vector::zero);

        if(!cleanInterface())
        {
            tangentialSurfaceTensionForce = 
                surfaceTensionGrad()().internalField();
        }
        else
        {
            vectorField surfaceTensionForce =
                cleanInterfaceSurfTension().value()
               *fac::edgeIntegrate
                (
                    aMesh().Le()*aMesh().edgeLengthCorrection()
                )().internalField();

            tangentialSurfaceTensionForce =
                surfaceTensionForce 
              - cleanInterfaceSurfTension().value()
               *aMesh().faceCurvatures().internalField()*nA;
        }

        UtFs += tangentialSurfaceTensionForce;

        UtFs /= muFluidA().value()*DnA + muFluidB().value()*DnB + VSMALL;

        Us().internalField() = UnFs + UtFs;
        correctUsBoundaryConditions();

        // Store old-time velocity field U()
        U().oldTime();

        U().boundaryField()[bPatchID()] == 
            interpolatorAB().faceInterpolate(UtFs)
          + nB*fvc::meshPhi(rho(),U())().boundaryField()[bPatchID()]/
            mesh().boundary()[bPatchID()].magSf();

        if
        (
            p().boundaryField()[bPatchID()].type()
         == fixedGradientFvPatchField<scalar>::typeName
        )
        {
            fixedGradientFvPatchField<scalar>& pB =
                refCast<fixedGradientFvPatchField<scalar> >
                (
                    p().boundaryField()[bPatchID()]
                );

            pB.gradient() = 
               - rhoFluidB().value()
                *(
                     nB&fvc::ddt(U())().boundaryField()[bPatchID()]
                 );
        }


        // Update fixedGradient boundary condition on patch A

        vectorField nGradU =
            muFluidB().value()*(UtPB - UtFs)*DnA
          + tangentialSurfaceTensionForce
          - muFluidA().value()*nA*fac::div(Us())().internalField()
          + (muFluidB().value() - muFluidA().value())
           *(fac::grad(Us())().internalField()&nA);

        nGradU /= muFluidA().value() + VSMALL;
        

        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else if
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
            FatalErrorIn("freeSurface::updateVelocity()")
                << "Bounary condition on " << U().name() 
                    <<  " for freeSurface patch is " 
                    << U().boundaryField()[aPatchID()].type() 
                    << ", instead " 
                    << fixedGradientCorrectedFvPatchField<vector>::typeName 
                    << " or "
                    << fixedGradientFvPatchField<vector>::typeName 
                    << abort(FatalError);
        }
    }
    else
    {
        vectorField nA = aMesh().faceAreaNormals().internalField();

        vectorField UnFs =
            nA*phi_.boundaryField()[aPatchID()]
           /mesh().boundary()[aPatchID()].magSf();

        // Correct normal component of free-surface velocity
        Us().internalField() += UnFs - nA*(nA&Us().internalField());
        correctUsBoundaryConditions();

        vectorField tangentialSurfaceTensionForce(nA.size(), vector::zero);

        if(!cleanInterface())
        {
             tangentialSurfaceTensionForce = 
                 surfaceTensionGrad()().internalField();
        }
        else
        {
            vectorField surfaceTensionForce =
                cleanInterfaceSurfTension().value()
               *fac::edgeIntegrate
                (
                    aMesh().Le()*aMesh().edgeLengthCorrection()
                )().internalField();

            tangentialSurfaceTensionForce =
                surfaceTensionForce 
              - cleanInterfaceSurfTension().value()
               *aMesh().faceCurvatures().internalField()*nA;

            if (muFluidA().value() < SMALL)
            {
                tangentialSurfaceTensionForce = vector::zero;
            }
        }

        vectorField tnGradU =
            tangentialSurfaceTensionForce/(muFluidA().value() + VSMALL)
          - (fac::grad(Us())&aMesh().faceAreaNormals())().internalField();

        vectorField UtPA =
            U().boundaryField()[aPatchID()].patchInternalField();
        UtPA -= nA*(nA & UtPA);

        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

        vectorField UtFs = UtPA + tnGradU/DnA;
        
        Us().internalField() = UtFs + UnFs;
        correctUsBoundaryConditions();

        vectorField nGradU =
            tangentialSurfaceTensionForce/(muFluidA().value() + VSMALL)
          - nA*fac::div(Us())().internalField()
          - (fac::grad(Us())().internalField()&nA);

        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else if
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
            FatalErrorIn("freeSurface::updateVelocity()")
                << "Bounary condition on " << U().name() 
                    <<  " for freeSurface patch is " 
                    << U().boundaryField()[aPatchID()].type() 
                    << ", instead " 
                    << fixedGradientCorrectedFvPatchField<vector>::typeName 
                    << " or "
                    << fixedGradientFvPatchField<vector>::typeName 
                    << abort(FatalError);
        }
    }
}


void freeSurface::updatePressure()
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

        Info << "Free surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K)
            << ", average = " << gAverage(K) << endl << flush;

        if(cleanInterface())
        {   
//             pA -= cleanInterfaceSurfTension().value()*(K - gAverage(K));
            pA -= cleanInterfaceSurfTension().value()*K;
        }
        else
        {
            scalarField surfTensionK =
                surfaceTension().internalField()*K;
                
            pA -= surfTensionK - gAverage(surfTensionK);
        }

        pA -= 2.0*(muFluidA().value() - muFluidB().value())
            *fac::div(Us())().internalField();

//         vector R0 = gAverage(mesh().C().boundaryField()[aPatchID()]);
        vector R0 = vector::zero;

        pA -= (rhoFluidA().value() - rhoFluidB().value())*
            (
                g_.value()
              & (
                    mesh().C().boundaryField()[aPatchID()]
                  - R0
                )
            );

        p().boundaryField()[aPatchID()] == pA;
    }
    else
    {
//         vector R0 = gAverage(mesh().C().boundaryField()[aPatchID()]);
        vector R0 = vector::zero;

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
        
        Info << "Free surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K) << ", average = " << gAverage(K) 
            << endl;
        
        if(cleanInterface())
        {
//             pA -= cleanInterfaceSurfTension().value()*(K - gAverage(K));
            pA -= cleanInterfaceSurfTension().value()*K;
        }
        else
        {
            scalarField surfTensionK =
                surfaceTension().internalField()*K;
            
            pA -= surfTensionK - gAverage(surfTensionK);
        }

        pA -= 2.0*muFluidA().value()*fac::div(Us())().internalField();

        p().boundaryField()[aPatchID()] == pA;
    }


    // Set modified pressure at patches with fixed apsolute
    // pressure

//     vector R0 = gAverage(mesh().C().boundaryField()[aPatchID()]);
    vector R0 = vector::zero;

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
                  - rho().boundaryField()[patchI]
                   *(g_.value()&(mesh().C().boundaryField()[patchI] - R0));
            }
        }
    }
}


void freeSurface::updateSurfaceFlux()
{
    Phis() = fac::interpolate(Us()) & aMesh().Le();
}


void freeSurface::updateSurfactantConcentration()
{
    if(!cleanInterface())
    {
        Info << "Correct surfactant concentration" << endl << flush;
        
        updateSurfaceFlux();        

        // Crate and solve the surfactanta transport equation
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


        if(surfactant().soluble())
        {
            const scalarField& C =
                mesh().boundary()[aPatchID()]
               .lookupPatchField<volScalarField, scalar>("C");

            areaScalarField Cb
            (
                IOobject
                (
                    "Cb",
                    DB().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                aMesh(),
                dimensioned<scalar>("Cb", dimMoles/dimVolume, 0),
                zeroGradientFaPatchScalarField::typeName
            );

            Cb.internalField() = C;
            Cb.correctBoundaryConditions();

            CsEqn += 
                fam::Sp
                (
                    surfactant().surfactAdsorptionCoeff()*Cb
                  + surfactant().surfactAdsorptionCoeff()
                   *surfactant().surfactDesorptionCoeff(),
                    surfactantConcentration()
                )
              - surfactant().surfactAdsorptionCoeff()
               *Cb*surfactant().surfactSaturatedConc();
        }

        CsEqn.solve();

        Info << "Correct surface tension" << endl;

        surfaceTension() =
            cleanInterfaceSurfTension()
          + surfactant().surfactR()
           *surfactant().surfactT()
           *surfactant().surfactSaturatedConc()
           *log(1.0 - surfactantConcentration()
           /surfactant().surfactSaturatedConc());
        
        if(neg(min(surfaceTension().internalField())))
        {
            FatalErrorIn
            (
                "void freeSurface::correctSurfactantConcentration()"
            ) 
                << "Surface tension is negative" 
                    << abort(FatalError);
        }
    }
}


void freeSurface::correctUsBoundaryConditions()
{
    forAll(Us().boundaryField(), patchI)
    {
        if
        (
            Us().boundaryField()[patchI].type()
         == calculatedFaPatchVectorField::typeName
        )
        {
            vectorField& pUs = Us().boundaryField()[patchI];

            pUs = Us().boundaryField()[patchI].patchInternalField();

            label ngbPolyPatchID = 
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            if(ngbPolyPatchID != -1)
            {
                if
                (
                    (
                        U().boundaryField()[ngbPolyPatchID].type() 
                     == slipFvPatchVectorField::typeName
                    )
                 ||
                    (
                        U().boundaryField()[ngbPolyPatchID].type() 
                     == symmetryFvPatchVectorField::typeName
                    )
                )
                {
                    vectorField N = 
                        aMesh().boundary()[patchI].ngbPolyPatchFaceNormals();

                    pUs -= N*(N&pUs);
                }
            }
        }
    }

    Us().correctBoundaryConditions();
}


vector freeSurface::totalPressureForce() const
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& P = p().boundaryField()[aPatchID()];

    vectorField pressureForces = S*P*n;

    return gSum(pressureForces);
}


vector freeSurface::totalViscousForce() const
{
    const scalarField& S = aMesh().S();
    const vectorField& n = aMesh().faceAreaNormals().internalField();

    vectorField nGradU =
        U().boundaryField()[aPatchID()].snGrad();

    vectorField viscousForces = 
      - muFluidA().value()*S
       *(
            nGradU
          + (fac::grad(Us())().internalField()&n)
          - (n*fac::div(Us())().internalField())
        );

    return gSum(viscousForces);
}


vector freeSurface::totalSurfaceTensionForce() const
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& K = aMesh().faceCurvatures().internalField();

    vectorField surfTensionForces(n.size(), vector::zero);

    if(cleanInterface())
    {
        surfTensionForces =
            S*cleanInterfaceSurfTension().value()
           *fac::edgeIntegrate
            (
                aMesh().Le()*aMesh().edgeLengthCorrection()
            )().internalField();
    }
    else
    {
        surfTensionForces *= surfaceTension().internalField()*K;
    }

    return gSum(surfTensionForces);
}


void freeSurface::initializeControlPointsPosition()
{
    scalarField deltaH = scalarField(aMesh().nFaces(), 0.0);

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
}


scalar freeSurface::maxCourantNumber()
{
    scalar CoNum = 0;

    if(cleanInterface())
    {
        const scalarField& dE =aMesh().lPN();

        CoNum = gMax
        (
            DB().deltaT().value()/
            sqrt
            (
                rhoFluidA().value()*dE*dE*dE/
                2.0/M_PI/(cleanInterfaceSurfTension().value() + SMALL)
            )
        );
    }
    else
    {
        scalarField sigmaE = 
            linearEdgeInterpolate(surfaceTension())().internalField()
          + SMALL;

        const scalarField& dE =aMesh().lPN();
        
        CoNum = gMax
        (
            DB().deltaT().value()/
            sqrt
            (
                rhoFluidA().value()*dE*dE*dE/
                2.0/M_PI/sigmaE
            )
        );
    }

    return CoNum;
}


void freeSurface::updateProperties()
{
    muFluidA_ = dimensionedScalar(this->lookup("muFluidA"));

    muFluidB_ = dimensionedScalar(this->lookup("muFluidB"));

    rhoFluidA_ = dimensionedScalar(this->lookup("rhoFluidA"));

    rhoFluidB_ = dimensionedScalar(this->lookup("rhoFluidB"));

    g_ = dimensionedVector(this->lookup("g"));

    cleanInterfaceSurfTension_ = 
        dimensionedScalar(this->lookup("surfaceTension"));
}


void freeSurface::writeVTK() const
{
    aMesh().patch().writeVTK
    (
        DB().timePath()/"freeSurface",
        aMesh().patch(),
        aMesh().patch().points()
    );
}


void freeSurface::writeVTKControlPoints()
{
    // Write patch and points into VTK
    fileName name(DB().timePath()/"freeSurfaceControlPoints");
    OFstream mps(name + ".vtk");

    mps << "# vtk DataFile Version 2.0" << nl
        << name << ".vtk" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << controlPoints().size() << " float" << nl;

    forAll(controlPoints(), pointI)
    {
        mps << controlPoints()[pointI].x() << ' '
            << controlPoints()[pointI].y() << ' '
            << controlPoints()[pointI].z() << nl;
    }

    // Write vertices
    mps << "VERTICES " << controlPoints().size() << ' ' 
        << controlPoints().size()*2 << nl;

    forAll(controlPoints(), pointI)
    {
        mps << 1 << ' ' << pointI << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
