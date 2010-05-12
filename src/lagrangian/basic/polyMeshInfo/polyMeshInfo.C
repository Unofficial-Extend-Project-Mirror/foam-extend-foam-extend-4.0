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

\*---------------------------------------------------------------------------*/

#include "polyMeshInfo.H"
#include "wedgePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "emptyPolyPatch.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyMeshInfo, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyMeshInfo::setCentrePoint()
{
    boundBox bounds = mesh_.bounds();

    centrePoint_ = 0.5*(bounds.max() + bounds.min());
}


void Foam::polyMeshInfo::setEmptyComponent(const vector& dir)
{
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        if (dir[cmpt] > 1.0e-6)
        {
            emptyComponent_ = cmpt;
        }
    }
}


void Foam::polyMeshInfo::queryWedge()
{
    label patchId[4];

    bool symmPlaneExists = false;

    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchi];
        if (isA<wedgePolyPatch>(patch))
        {
            if (patch.size())
            {
                nWedge_++;
                if (nWedge_ > 4)
                {
                    break;
                }
                patchId[nWedge_-1] = patchi;
           }
        }
        else if (isA<symmetryPolyPatch>(patch))
        {
            symmPlaneExists = true;
        }
    }

    if (nWedge_ != 0 && nWedge_ != 2 && nWedge_ != 4)
    {
        FatalErrorIn("void polyMeshInfo::queryWedge() const")
            << "Number of wedge patches " << nWedge_ << " is incorrect, "
               "should be 0, 2 or 4"
            << exit(FatalError);
    }
    else
    {
        if (symmPlaneExists)
        {
            wedgeAngle_ = mathematicalConstant::pi;
        }
        else
        {
            if (nWedge_)
            {
               if (nWedge_ > 2)
                {
                    WarningIn("void polyMeshInfo::queryWedge() const")
                        << "Only configured for cases with 2 wedges. "
                        << "Lagrangian calculations may result in errors"
                        << endl;
                }
               // Get the vector normal to wedge patch1
                const wedgePolyPatch& patch1 = refCast<const wedgePolyPatch>
                    (mesh_.boundaryMesh()[patchId[0]]);
                vector n1 = patch1.patchNormal();

                // Get the vector normal to wedge patch2
                const wedgePolyPatch& patch2 = refCast<const wedgePolyPatch>
                    (mesh_.boundaryMesh()[patchId[1]]);
                vector n2 = patch2.patchNormal();

                // Calculate the angle swept between the patches
                const scalar arcCos = (n1 & n2)/mag(n1);
                wedgeAngle_ = mathematicalConstant::pi - acos(arcCos);

                // Get the centre normal
                centreNormal_ = patch1.centreNormal();

                // Get the wedge axis
                wedgeAxis_  = patch1.axis();
            }
        }
    }
}


void Foam::polyMeshInfo::queryDirections()
{
    vector dirVec = vector::zero;

    forAll(mesh_.boundaryMesh(), patchi)
    {
        if (isA<emptyPolyPatch>(mesh_.boundaryMesh()[patchi]))
        {
            if (mesh_.boundaryMesh()[patchi].size())
            {
                nEmpty_++;
                dirVec +=
                    sum(cmptMag(mesh_.boundaryMesh()[patchi].faceAreas()));
            }
        }
    }


    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        directions_[cmpt] = 1;
    }

    if (nEmpty_)
    {
        reduce(dirVec, sumOp<vector>());

        dirVec /= mag(dirVec);

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            if (dirVec[cmpt] > 1.0e-6)
            {
                directions_[cmpt] = -1;
                emptyComponent_ = cmpt;
            }
            else
            {
                directions_[cmpt] = 1;
            }
        }

        // Set the patch normal
        centreNormal_ = dirVec;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyMeshInfo::polyMeshInfo
(
    const polyMesh& mesh
)
:
    mesh_(mesh),
    nGeometricD_(0),
    nSolutionD_(0),
    centreNormal_(vector::zero),
    emptyComponent_(0),
    centrePoint_(vector::zero),
    nEmpty_(0),
    nWedge_(0),
    wedgeAxis_(vector::zero),
    wedgeAngle_(0)
{
    setCentrePoint();

    queryWedge();

    queryDirections();

    nSolutionD_ = cmptSum(directions_ + Vector<label>::one)/2;

    nGeometricD_ = nSolutionD_ - nWedge_/2;

    if (nGeometricD_ == 2)
    {
        setEmptyComponent(centreNormal_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyMeshInfo::~polyMeshInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::polyMeshInfo::nGeometricD() const
{
    return nGeometricD_;
}


Foam::label Foam::polyMeshInfo::nSolutionD() const
{
    return nSolutionD_;
}


bool Foam::polyMeshInfo::caseIs2dWedge() const
{
    return (nWedge_ > 0);
}


bool Foam::polyMeshInfo::caseIs2dSlab() const
{
    return (nEmpty_ > 0);
}


bool Foam::polyMeshInfo::caseIs2d() const
{
    return (nGeometricD_ == 2);
}


bool Foam::polyMeshInfo::caseIs3d() const
{
    return (nGeometricD_ == 3);
}


Foam::vector Foam::polyMeshInfo::wedgeAxis() const
{
    if (nWedge_)
    {
        return wedgeAxis_;
    }
    else
    {
        WarningIn("Foam::polyMeshInfo::wedgeAxis()")
            << "wedgeAxis() requested, but case is not of wedge type. "
            << "Returning zero vector" << endl;
        return vector::zero;
    }
}


Foam::scalar Foam::polyMeshInfo::wedgeAngle() const
{
    if (nWedge_)
    {
        return wedgeAngle_;
    }
    else
    {
        WarningIn("Foam::polyMeshInfo::wedgeAngle()")
            << "wedgeAngle() requested, but case is not of wedge type. "
            << "Returning zero" << endl;
        return 0.0;
    }
}


Foam::vector Foam::polyMeshInfo::centreNormal() const
{
    if (nGeometricD_ == 2)
    {
        return centreNormal_;
    }
    else
    {
        WarningIn("Foam::polyMeshInfo::centreNormal()")
            << "centreNormal() requested, but case is not 2-D. "
            << "Returning zero vector" << endl;
        return vector::zero;
    }
}


Foam::label Foam::polyMeshInfo::emptyComponent() const
{
    return emptyComponent_;
}


Foam::vector Foam::polyMeshInfo::centrePoint() const
{
    return centrePoint_;
}


// ************************************************************************* //
