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

Description

\*---------------------------------------------------------------------------*/

#include "wedgeFaPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "faBoundaryMesh.H"
#include "wedgePolyPatch.H"
#include "polyMesh.H"
#include "faMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(wedgeFaPatch, 0);
addToRunTimeSelectionTable(faPatch, wedgeFaPatch, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void wedgeFaPatch::findAxisPoint() const
{
    // Find axis point

    labelList ptLabels = pointLabels();
    
    labelListList ptEdges = pointEdges();

    const vectorField& points = boundaryMesh().mesh().points();

    const scalarField& magL = magEdgeLengths();

    forAll(ptEdges, pointI)
    {
        if( ptEdges[pointI].size() == 1 )
        {
            scalar r = mag((I-axis()*axis())&points[ptLabels[pointI]]);

            if( r < magL[ptEdges[pointI][0]] )
            {
                axisPoint_ = ptLabels[pointI];
                break;
            }
        }
    }

    axisPointChecked_ = true;
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

//- Construct from polyPatch
wedgeFaPatch::wedgeFaPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const faBoundaryMesh& bm
)
:
    faPatch(name, dict, index, bm),
    wedgePolyPatchPtr_(NULL),
    axisPoint_(-1),
    axisPointChecked_(false)
{
    if(ngbPolyPatchIndex() == -1)
    {
        FatalErrorIn
        (
            "wedgeFaPatch::wedgeFaPatch(const word&, const dictionary&, const label, const faBoundaryMesh&)"
        )   << "Neighbour polyPatch index is not specified for faPatch "
            << this->name() << exit(FatalError);
    }

    if
    (
        bm.mesh()().boundaryMesh()[ngbPolyPatchIndex()].type()
     == wedgePolyPatch::typeName
    )
    {
        const wedgePolyPatch& wedge =
            refCast<const wedgePolyPatch>
            (
                bm.mesh()().boundaryMesh()[ngbPolyPatchIndex()]
            );

        wedgePolyPatchPtr_ = &wedge;
    }
    else
    {
        FatalErrorIn
        (
            "wedgeFaPatch::wedgeFaPatch(const word&, const dictionary&, const label, const faBoundaryMesh&)"
        )   << "Neighbour polyPatch is not of type "
            << wedgePolyPatch::typeName
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
