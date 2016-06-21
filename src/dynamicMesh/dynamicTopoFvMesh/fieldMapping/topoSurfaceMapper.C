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

Class
    topoSurfaceMapper

Description
    Implementation of the topoSurfaceMapper class

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "topoMapper.H"
#include "mapPolyMesh.H"
#include "topoSurfaceMapper.H"
#include "demandDrivenData.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Clear out local storage
void topoSurfaceMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);
    deleteDemandDrivenData(insertedFaceLabelsPtr_);
}


//- Calculate the insertedFaceLabels list
void topoSurfaceMapper::calcInsertedFaceLabels() const
{
    if (insertedFaceLabelsPtr_)
    {
        FatalErrorIn
        (
            "void topoSurfaceMapper::calcInsertedFaceLabels() const"
        )   << " Inserted labels has already been calculated."
            << abort(FatalError);
    }

    // Allocate for inserted face labels
    label nInsertedFaces = 0;

    insertedFaceLabelsPtr_ = new labelList(size(), -1);
    labelList& insertedFaces = *insertedFaceLabelsPtr_;

    label nIntFaces = size();

    // Loop through the facesFromFaces map, and ensure that
    // inserted internal faces are not mapped from any parent faces.
    const List<objectMap>& fff = mpm_.facesFromFacesMap();

    forAll(fff, objectI)
    {
        const objectMap& fffI = fff[objectI];

        // Only pick internal faces
        if (fffI.index() < nIntFaces)
        {
            insertedFaces[nInsertedFaces++] = fffI.index();
        }
    }

    // Shorten inserted faces to actual size
    insertedFaces.setSize(nInsertedFaces);
}


//- Calculate addressing for mapping
void topoSurfaceMapper::calcAddressing() const
{
    if
    (
        directAddrPtr_
     || interpolationAddrPtr_
    )
    {
        FatalErrorIn("void topoSurfaceMapper::calcAddressing() const")
            << "Addressing already calculated."
            << abort(FatalError);
    }

    if (direct())
    {
        // Direct addressing, no weights - slice to size
        directAddrPtr_ =
        (
            new labelList
            (
                labelList::subList(mpm_.faceMap(), size())
            )
        );

        labelList& addr = *directAddrPtr_;

        // Loop through the addressing list,
        // and set dummy addressing for newly created faces.
        forAll(addr, faceI)
        {
            if (addr[faceI] < 0)
            {
                addr[faceI] = 0;
            }
        }
    }
    else
    {
        FatalErrorIn("void topoSurfaceMapper::calcAddressing() const")
            << " Request for interpolative addressing"
            << " on internal surface fields. The surfaceMapper"
            << " is not designed to do this. Interpolate volFields"
            << " to surfaceFields instead."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
topoSurfaceMapper::topoSurfaceMapper
(
    const mapPolyMesh& mpm,
    const topoMapper& mapper
)
:
    mesh_(mpm.mesh()),
    mpm_(mpm),
    tMapper_(mapper),
    direct_(false),
    sizeBeforeMapping_(mpm.nOldInternalFaces()),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL),
    insertedFaceLabelsPtr_(NULL)
{
    // Fetch offset sizes from topoMapper
    const labelList& sizes = tMapper_.faceSizes();

    // Add offset sizes
    if (sizes.size())
    {
        forAll(sizes, pI)
        {
            sizeBeforeMapping_ += sizes[pI];
        }
    }

    // Calculate the insertedFaces list
    calcInsertedFaceLabels();

    // Set to direct mapping.
    direct_ = true;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

topoSurfaceMapper::~topoSurfaceMapper()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return size
label topoSurfaceMapper::size() const
{
    return mesh_.nInternalFaces();
}


//- Return size before mapping
label topoSurfaceMapper::sizeBeforeMapping() const
{
    return sizeBeforeMapping_;
}


//- Is the mapping direct
bool topoSurfaceMapper::direct() const
{
    return direct_;
}


//- Return direct addressing
const unallocLabelList& topoSurfaceMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const unallocLabelList& "
            "topoSurfaceMapper::directAddressing() const"
        )   << "Requested direct addressing for an interpolative mapper."
            << abort(FatalError);
    }

    if (!directAddrPtr_)
    {
        calcAddressing();
    }

    return *directAddrPtr_;
}


//- Return interpolation addressing
const labelListList& topoSurfaceMapper::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& "
            "topoSurfaceMapper::addressing() const"
        )   << "Requested interpolative addressing for a direct mapper."
            << abort(FatalError);
    }

    if (!interpolationAddrPtr_)
    {
        calcAddressing();
    }

    return *interpolationAddrPtr_;
}


//- Return weights
const scalarListList& topoSurfaceMapper::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList& "
            "topoSurfaceMapper::weights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!weightsPtr_)
    {
        calcAddressing();
    }

    return *weightsPtr_;
}


//- Are there any inserted faces
bool topoSurfaceMapper::insertedObjects() const
{
    return insertedObjectLabels().size();
}


//- Return list of inserted faces
const labelList& topoSurfaceMapper::insertedObjectLabels() const
{
    if (!insertedFaceLabelsPtr_)
    {
        calcInsertedFaceLabels();
    }

    return *insertedFaceLabelsPtr_;
}


//- Return flux flip map
const labelHashSet& topoSurfaceMapper::flipFaceFlux() const
{
    return mpm_.flipFaceFlux();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void topoSurfaceMapper::operator=(const topoSurfaceMapper& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("void topoSurfaceMapper::operator=")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
