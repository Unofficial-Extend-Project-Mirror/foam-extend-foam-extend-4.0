/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Description
    FV surface mapper.

\*---------------------------------------------------------------------------*/

#include "faAreaMapper.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faAreaMapper::calcAddressing() const
{
    if
    (
        newFaceLabelsPtr_
     || newFaceLabelsMapPtr_
     || directAddrPtr_
     || interpolationAddrPtr_
     || weightsPtr_
     || insertedObjectLabelsPtr_
    )
    {
        FatalErrorIn("void faAreaMapper::calcAddressing() const)")
            << "Addressing already calculated"
            << abort(FatalError);
    }

    // Mapping

    const label oldNInternal = mpm_.nOldInternalFaces();

    // Calculate new face labels

    // Copy old face labels
    const labelList& oldFaces = mesh_.faceLabels();

    // Prepare a list of new face labels and (preliminary) addressing
    // Note: dimensioned to number of boundary faces of polyMesh
    // HJ, 10/Aug/2011
    newFaceLabelsPtr_ = new labelList
    (
        mesh_().nFaces() - mesh_().nInternalFaces(),
        -1
    );
    labelList& newFaceLabels = *newFaceLabelsPtr_;

    newFaceLabelsMapPtr_ = new labelList
    (
        mesh_().nFaces() - mesh_().nInternalFaces(),
        -1
    );
    labelList& newFaceLabelsMap = *newFaceLabelsMapPtr_;
    label nNewFaces = 0;

    Info<< "Old face list size: " << oldFaces.size()
        << " estimated new size " << newFaceLabels.size() << endl;

    // Get reverse face map
    const labelList& reverseFaceMap = mpm_.reverseFaceMap();

    // Pick up live old faces
    forAll (oldFaces, faceI)
    {
        if (reverseFaceMap[oldFaces[faceI]] > -1)
        {
            // Face is live, add it and record addressing
            newFaceLabels[nNewFaces] = reverseFaceMap[oldFaces[faceI]];
            newFaceLabelsMap[nNewFaces] = faceI;

            nNewFaces++;
        }
    }

    // Assemble the maps
    if (direct())
    {
        Info<< "Direct"<< endl;
        // Direct mapping: no further faces to add.  Resize list
        newFaceLabels.setSize(nNewFaces);

        directAddrPtr_ = new labelList(newFaceLabels.size());
        labelList& addr = *directAddrPtr_;

        // Adjust for creation of a boundary face from an internal face
        forAll (addr, faceI)
        {
            if (newFaceLabelsMap[faceI] < oldNInternal)
            {
                addr[faceI] = 0;
            }
            else
            {
                addr[faceI] = newFaceLabelsMap[faceI];
            }
        }
    }
    else
    {
        // There are further faces to add.  Prepare interpolation addressing
        // and weights to full size
        interpolationAddrPtr_ = new labelListList(newFaceLabels.size());
        labelListList& addr = *interpolationAddrPtr_;

        weightsPtr_ = new scalarListList(newFaceLabels.size());
        scalarListList& w = *weightsPtr_;

        // Insert single addressing and weights
        for (label addrI = 0; addrI < nNewFaces; addrI++)
        {
            addr[addrI] = labelList(1, newFaceLabelsMap[addrI]);
            w[addrI] = scalarList(1, scalar(1));
        }

        // Pick up faces from points, edges and faces where the origin
        // Only map from faces which were previously in the faMesh, using
        // fast lookup

        // Set of faces previously in the mesh
        labelHashSet oldFaceLookup(oldFaces);

        // Go through faces-from lists and add the ones where all
        // old face labels belonged to the faMesh

        const List<objectMap>& ffp = mpm_.facesFromPointsMap();

        forAll (ffp, ffpI)
        {
            // Get addressing
            const labelList& mo = ffp[ffpI].masterObjects();

            // Check if master objects are in faMesh
            labelList validMo(mo.size());
            label nValidMo = 0;

            forAll (mo, moI)
            {
                if (oldFaceLookup.found(mo[moI]))
                {
                    validMo[nValidMo] = oldFaceLookup[mo[moI]];
                    nValidMo++;
                }
            }

            if (nValidMo > 0)
            {
                // Some objects found: add face and interpolation to list
                newFaceLabels[nNewFaces] = ffp[ffpI].index();

                // No old face available
                newFaceLabelsMap[nNewFaces] = -1;

                // Map from masters, uniform weights
                addr[nNewFaces] = validMo;
                w[nNewFaces] = scalarList(validMo.size(), 1.0/validMo.size());

                nNewFaces++;
            }
        }

        const List<objectMap>& ffe = mpm_.facesFromEdgesMap();

        forAll (ffe, ffeI)
        {
            // Get addressing
            const labelList& mo = ffe[ffeI].masterObjects();

            // Check if master objects are in faMesh
            labelList validMo(mo.size());
            label nValidMo = 0;

            forAll (mo, moI)
            {
                if (oldFaceLookup.found(mo[moI]))
                {
                    validMo[nValidMo] = oldFaceLookup[mo[moI]];
                    nValidMo++;
                }
            }

            if (nValidMo > 0)
            {
                // Some objects found: add face and interpolation to list
                newFaceLabels[nNewFaces] = ffe[ffeI].index();

                // No old face available
                newFaceLabelsMap[nNewFaces] = -1;

                // Map from masters, uniform weights
                addr[nNewFaces] = validMo;
                w[nNewFaces] = scalarList(validMo.size(), 1.0/validMo.size());

                nNewFaces++;
            }
        }

        const List<objectMap>& fff = mpm_.facesFromFacesMap();

        forAll (fff, fffI)
        {
            // Get addressing
            const labelList& mo = fff[fffI].masterObjects();

            // Check if master objects are in faMesh
            labelList validMo(mo.size());
            label nValidMo = 0;

            forAll (mo, moI)
            {
                if (oldFaceLookup.found(mo[moI]))
                {
                    validMo[nValidMo] = oldFaceLookup[mo[moI]];
                    nValidMo++;
                }
            }

            if (nValidMo > 0)
            {
                // Some objects found: add face and interpolation to list
                newFaceLabels[nNewFaces] = fff[fffI].index();

                // No old face available
                newFaceLabelsMap[nNewFaces] = -1;

                // Map from masters, uniform weights
                addr[nNewFaces] = validMo;
                w[nNewFaces] = scalarList(validMo.size(), 1.0/validMo.size());

                nNewFaces++;
            }
        }

        // All faces collected.  Reset sizes of lists
        newFaceLabels.setSize(nNewFaces);
        newFaceLabelsMap.setSize(nNewFaces);
        addr.setSize(nNewFaces);
        w.setSize(nNewFaces);
        Info<< "addr: " << addr << nl
            << "w: " << w << endl;
    }

    // Inserted objects cannot appear in the new faMesh as they have no master
    // HJ, 10/Aug/2011
    insertedObjectLabelsPtr_ = new labelList(0);
}


void Foam::faAreaMapper::clearOut()
{
    deleteDemandDrivenData(newFaceLabelsPtr_);
    deleteDemandDrivenData(newFaceLabelsMapPtr_);

    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);

    deleteDemandDrivenData(insertedObjectLabelsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::faAreaMapper::faAreaMapper
(
    const faMesh& mesh,
    const mapPolyMesh& mpm
)
:
    mesh_(mesh),
    mpm_(mpm),
    insertedFaces_(false),
    direct_(false),
    sizeBeforeMapping_(mesh.nFaces()),
    newFaceLabelsPtr_(NULL),
    newFaceLabelsMapPtr_(NULL),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL),
    insertedObjectLabelsPtr_(NULL)
{
    // Check for possibility of direct mapping
    if
    (
        mpm_.facesFromPointsMap().empty()
     && mpm_.facesFromEdgesMap().empty()
     && mpm_.facesFromFacesMap().empty()
    )
    {
        direct_ = true;
    }
    else
    {
        direct_ = false;
    }

    // Inserted objects not suported: no master
    // HJ, 10/Aug/2011
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faAreaMapper::~faAreaMapper()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::faAreaMapper::newFaceLabels() const
{
    if (!newFaceLabelsPtr_)
    {
        calcAddressing();
    }

    return *newFaceLabelsPtr_;
}


const Foam::labelList& Foam::faAreaMapper::newFaceLabelsMap() const
{
    if (!newFaceLabelsMapPtr_)
    {
        calcAddressing();
    }

    return *newFaceLabelsMapPtr_;
}


const Foam::unallocLabelList& Foam::faAreaMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const unallocLabelList& faAreaMapper::"
            "directAddressing() const"
        )   << "Requested direct addressing for an interpolative mapper."
            << abort(FatalError);
    }

    if (!directAddrPtr_)
    {
        calcAddressing();
    }

    return *directAddrPtr_;
}


const Foam::labelListList& Foam::faAreaMapper::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& faAreaMapper::addressing() const"
        )   << "Requested interpolative addressing for a direct mapper."
            << abort(FatalError);
    }

    if (!interpolationAddrPtr_)
    {
        calcAddressing();
    }

    return *interpolationAddrPtr_;
}


const Foam::scalarListList& Foam::faAreaMapper::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList& faAreaMapper::weights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!weightsPtr_)
    {
        calcAddressing();
    }

    return *weightsPtr_;
}


const Foam::labelList& Foam::faAreaMapper::insertedObjectLabels() const
{
    if (!insertedObjectLabelsPtr_)
    {
        calcAddressing();
    }

    return *insertedObjectLabelsPtr_;
}


// ************************************************************************* //
