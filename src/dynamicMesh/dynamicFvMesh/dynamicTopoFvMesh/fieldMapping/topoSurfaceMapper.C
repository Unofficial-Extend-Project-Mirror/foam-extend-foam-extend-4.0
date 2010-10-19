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
    topoSurfaceMapper

Description
    Implementation of the topoSurfaceMapper class

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*----------------------------------------------------------------------------*/

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
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL),
    insertedFaceLabelsPtr_(NULL)
{
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
    return mpm_.nOldInternalFaces();
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


//- Map the internal field
template <class Type>
void topoSurfaceMapper::mapInternalField
(
    const word& fieldName,
    Field<Type>& iF
) const
{
    if (iF.size() != sizeBeforeMapping())
    {
        FatalErrorIn
        (
            "\n\n"
            "void topoSurfaceMapper::mapInternalField<Type>\n"
            "(\n"
            "    Field<Type>& iF\n"
            ") const\n"
        )  << "Incompatible size before mapping." << nl
           << " Field: " << fieldName << nl
           << " Field size: " << iF.size() << nl
           << " map size: " << sizeBeforeMapping() << nl
           << abort(FatalError);
    }

    // Map the internal field
    iF.autoMap(*this);

    // Flip the flux
    const labelList flipFaces = flipFaceFlux().toc();

    forAll (flipFaces, i)
    {
        if (flipFaces[i] < iF.size())
        {
            iF[flipFaces[i]] *= -1.0;
        }
        else
        {
            FatalErrorIn
            (
                "\n\n"
                "void topoSurfaceMapper::mapInternalField<Type>\n"
                "(\n"
                "    Field<Type>& iF\n"
                ") const\n"
            )  << "Cannot flip boundary face fluxes." << nl
               << " Field: " << fieldName << nl
               << " Field size: " << iF.size() << nl
               << " Face flip index: " << flipFaces[i] << nl
               << abort(FatalError);
        }
    }
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
