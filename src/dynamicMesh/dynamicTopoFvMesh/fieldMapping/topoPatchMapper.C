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
    topoPatchMapper

Description
    Implementation of the topoPatchMapper class

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "IOmanip.H"
#include "meshOps.H"
#include "topoMapper.H"
#include "mapPolyMesh.H"
#include "topoPatchMapper.H"
#include "demandDrivenData.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Clear out local storage
void topoPatchMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);
    deleteDemandDrivenData(insertedFaceLabelsPtr_);
    deleteDemandDrivenData(insertedFaceIndexMapPtr_);
    deleteDemandDrivenData(insertedFaceAddressingPtr_);
    deleteDemandDrivenData(areasPtr_);
    deleteDemandDrivenData(centresPtr_);
}


//- Calculate the insertedFaceLabels list
void topoPatchMapper::calcInsertedFaceAddressing() const
{
    if (insertedFaceLabelsPtr_ || insertedFaceAddressingPtr_)
    {
        FatalErrorIn
        (
            "void topoPatchMapper::calcInsertedFaceAddressing() const"
        )   << " Inserted labels has already been calculated."
            << abort(FatalError);
    }

    // Information from the old patch
    const label oldPatchSize = sizeBeforeMapping();
    const label oldPatchStart = mpm_.oldPatchStarts()[patch_.index()];

    if (oldPatchSize == 0)
    {
        // Nothing to map from. Return empty.
        insertedFaceLabelsPtr_ = new labelList(0);
        insertedFaceIndexMapPtr_ = new labelList(0);
        insertedFaceAddressingPtr_ = new labelListList(0);
        return;
    }

    // Allocate for inserted face labels and addressing
    label nInsertedFaces = 0;

    insertedFaceLabelsPtr_ = new labelList(patchSize(), -1);
    labelList& insertedFaces = *insertedFaceLabelsPtr_;

    insertedFaceIndexMapPtr_ = new labelList(patchSize(), -1);
    labelList& insertedFacesMap = *insertedFaceIndexMapPtr_;

    insertedFaceAddressingPtr_ = new labelListList(patchSize(), labelList(0));
    labelListList& insertedAddressing = *insertedFaceAddressingPtr_;

    // Fetch current patch start
    label pStart = patch_.patch().start();

    // Fetch the current boundary
    const polyBoundaryMesh& boundary = mpm_.mesh().boundaryMesh();

    // Loop through the facesFromFaces map, and ensure that
    // inserted faces are only mapped from faces on the same patch.
    const List<objectMap>& fff = mpm_.facesFromFacesMap();

    forAll(fff, objectI)
    {
        const objectMap& fffI = fff[objectI];

        // Only pick boundary faces in this patch
        if (boundary.whichPatch(fffI.index()) == patch_.index())
        {
            if (fffI.masterObjects().empty())
            {
                // Write out for post-processing
                meshOps::writeVTK
                (
                    mpm_.mesh(),
                    "patchFaceInsError_"
                  + Foam::name(fffI.index()),
                    labelList(1, fffI.index()),
                    2,
                    mpm_.mesh().points(),
                    List<edge>(0),
                    mpm_.mesh().faces()
                );

                FatalErrorIn
                (
                    "void topoPatchMapper::"
                    "calcInsertedFaceAddressing() const"
                )   << " Mapping for inserted boundary face is incorrect."
                    << " Found an empty masterObjects list."
                    << nl << " Face: " << fffI.index()
                    << nl << " Patch: " << patch_.name()
                    << abort(FatalError);
            }
            else
            {
                // Make an entry for the inserted label,
                // and renumber addressing to patch.
                insertedFaces[nInsertedFaces] = fffI.index() - pStart;

                // Add a mapping entry for facesFromFaces indices
                insertedFacesMap[nInsertedFaces] = objectI;

                // Make an entry for addressing
                labelList& addr = insertedAddressing[nInsertedFaces];

                // Check for illegal addressing
                addr = fffI.masterObjects();

                forAll(addr, faceI)
                {
                    if (addr[faceI] < 0 || addr[faceI] >= oldPatchSize)
                    {
                        // Write out for post-processing
                        meshOps::writeVTK
                        (
                            mpm_.mesh(),
                            "patchFacePatchError_"
                          + Foam::name(fffI.index()),
                            labelList(1, fffI.index()),
                            2,
                            mpm_.mesh().points(),
                            List<edge>(0),
                            mpm_.mesh().faces()
                        );

                        FatalErrorIn
                        (
                            "void topoPatchMapper::"
                            "calcInsertedFaceAddressing() const"
                        )
                            << "Addressing into another patch is not allowed."
                            << nl << " Patch face index: " << faceI
                            << nl << " Patch: " << patch_.name()
                            << nl << " fffI.index: " << fffI.index()
                            << nl << " pStart: " << pStart
                            << nl << " addr: " << addr
                            << nl << " addr[faceI]: " << addr[faceI]
                            << nl << " oldPatchStart: " << oldPatchStart
                            << nl << " oldPatchSize: " << oldPatchSize
                            << abort(FatalError);
                    }
                }

                nInsertedFaces++;
            }
        }
    }

    // Shorten inserted faces to actual size
    insertedFaces.setSize(nInsertedFaces);
    insertedFacesMap.setSize(nInsertedFaces);
    insertedAddressing.setSize(nInsertedFaces);
}


//- Calculate addressing for mapping
void topoPatchMapper::calcAddressing() const
{
    if (directAddrPtr_ || interpolationAddrPtr_)
    {
        FatalErrorIn("void topoPatchMapper::calcAddressing() const")
            << "Addressing already calculated."
            << abort(FatalError);
    }

    // Information from the old patch
    const label oldPatchSize = sizeBeforeMapping();
    const label oldPatchStart = mpm_.oldPatchStarts()[patch_.index()];
    const label oldPatchEnd = oldPatchStart + oldPatchSize;

    // Assemble the maps: slice to patch
    if (direct())
    {
        if (oldPatchSize == 0)
        {
            // Nothing to map from. Return empty.
            directAddrPtr_ = new labelList(0);
            return;
        }

        // Direct mapping - slice to size
        directAddrPtr_ =
        (
            new labelList(patch_.patch().patchSlice(mpm_.faceMap()))
        );

        labelList& addr = *directAddrPtr_;

        // Shift to local patch indices.
        // Also, check mapping for hits into other patches / internal faces.
        forAll(addr, faceI)
        {
            if
            (
                addr[faceI] >= oldPatchStart
             && addr[faceI] < oldPatchEnd
            )
            {
                addr[faceI] -= oldPatchStart;
            }
            else
            {
                // Relax addressing requirement for
                // processor patch faces. These require
                // cell-to-face interpolation anyway.
                if (isA<processorPolyPatch>(patch_.patch()))
                {
                    // Artificially map from face[0] of this patch.
                    addr[faceI] = 0;
                    continue;
                }

                // Write out for post-processing
                meshOps::writeVTK
                (
                    mpm_.mesh(),
                    "patchFacePatchError_"
                  + Foam::name(faceI),
                    labelList(1, faceI),
                    2,
                    mpm_.mesh().points(),
                    List<edge>(0),
                    mpm_.mesh().faces()
                );

                FatalErrorIn
                (
                    "void topoPatchMapper::calcAddressing() const"
                )
                    << "Addressing into another patch is not allowed."
                    << nl << " Patch: " << patch_.name()
                    << nl << " Patch index: " << patch_.index()
                    << nl << " Patch face index: " << faceI
                    << nl << " addr[faceI]: " << addr[faceI]
                    << nl << " oldPatchStart: " << oldPatchStart
                    << nl << " oldPatchSize: " << oldPatchSize
                    << nl << " oldPatchEnd: " << oldPatchEnd
                    << nl << " nInserted: " << insertedObjectLabels().size()
                    << abort(FatalError);
            }
        }
    }
    else
    {
        if (oldPatchSize == 0)
        {
            // Nothing to map from. Return empty.
            interpolationAddrPtr_ = new labelListList(0);
            return;
        }

        // Interpolative addressing
        interpolationAddrPtr_ = new labelListList(patchSize(), labelList(0));
        labelListList& addr = *interpolationAddrPtr_;

        // Fetch the list of inserted faces / addressing
        const labelList& insertedFaces = insertedObjectLabels();
        const labelListList& insertedAddressing = insertedFaceAddressing();

        // Make entries
        forAll(insertedFaces, faceI)
        {
            addr[insertedFaces[faceI]] = insertedAddressing[faceI];
        }

        // Do mapped faces. Note that this can already be set by insertedFaces
        // so check if addressing size still zero.
        const labelList::subList fm = patch_.patch().patchSlice(mpm_.faceMap());

        forAll(fm, faceI)
        {
            if (fm[faceI] > -1 && addr[faceI].size() == 0)
            {
                // Mapped from a single face
                label oldFace = fm[faceI];

                if
                (
                    oldFace >= oldPatchStart
                 && oldFace < oldPatchEnd
                )
                {
                    oldFace -= oldPatchStart;
                }
                else
                {
                    FatalErrorIn
                    (
                        "void topoPatchMapper::calcAddressing() const"
                    )
                        << "Addressing into another patch is not allowed."
                        << nl << " Patch: " << patch_.name()
                        << nl << " Patch index: " << patch_.index()
                        << nl << " Patch face index: " << faceI
                        << nl << " faceMap[faceI]: " << oldFace
                        << nl << " oldPatchStart: " << oldPatchStart
                        << nl << " oldPatchSize: " << oldPatchSize
                        << nl << " oldPatchEnd: " << oldPatchEnd
                        << abort(FatalError);
                }

                addr[faceI] = labelList(1, oldFace);
            }
        }

        // Check if we missed anything
        forAll(addr, faceI)
        {
            if (addr[faceI].empty())
            {
                // Relax addressing requirement for
                // processor patch faces. These require
                // cell-to-face interpolation anyway.
                if (isA<processorPolyPatch>(patch_.patch()))
                {
                    // Artificially map from face[0] of this patch.
                    addr[faceI] = labelList(1, 0);

                    continue;
                }

                label oldFace = (faceI < fm.size() ? fm[faceI] : -1);

                FatalErrorIn
                (
                    "void topoPatchMapper::calcAddressing() const"
                )
                    << "Addressing is missing." << nl
                    << " Patch face index: " << faceI << nl
                    << " nInsertedFaces: " << insertedFaces.size() << nl
                    << " faceMap[faceI]: " << oldFace << nl
                    << " patch faceMap size: " << fm.size() << nl
                    << " Patch: " << patch_.name() << nl
                    << " polyPatch: " << nl << patch_.patch() << nl
                    << " faceMap size: " << mpm_.faceMap().size() << nl
                    << abort(FatalError);
            }
        }
    }
}


//- Calculate inverse-distance weights for interpolative mapping
void topoPatchMapper::calcInverseDistanceWeights() const
{
    if (weightsPtr_)
    {
        FatalErrorIn
        (
            "void topoPatchMapper::calcInverseDistanceWeights() const"
        )
            << "Weights already calculated."
            << abort(FatalError);
    }

    // Fetch interpolative addressing
    const labelListList& addr = addressing();

    if (addr.empty())
    {
        // Nothing to map from. Return empty.
        weightsPtr_ = new scalarListList(0);
        return;
    }

    // Allocate memory
    weightsPtr_ = new scalarListList(patchSize());
    scalarListList& w = *weightsPtr_;

    // Obtain face-centre information from old/new meshes
    const vectorField& oldCentres = tMapper_.patchCentres(patch_.index());
    const vectorField& newCentres = patch_.patch().faceCentres();

    forAll(addr, faceI)
    {
        const labelList& mo = addr[faceI];

        // Do mapped faces
        if (mo.size() == 1)
        {
            w[faceI] = scalarList(1, 1.0);
        }
        else
        {
            // Map from masters, inverse-distance weights
            scalar totalWeight = 0.0;
            w[faceI] = scalarList(mo.size(), 0.0);

            forAll (mo, oldFaceI)
            {
                w[faceI][oldFaceI] =
                (
                    1.0/stabilise
                    (
                        magSqr
                        (
                            newCentres[faceI]
                          - oldCentres[mo[oldFaceI]]
                        ),
                        VSMALL
                    )
                );

                totalWeight += w[faceI][oldFaceI];
            }

            // Normalize weights
            scalar normFactor = (1.0/totalWeight);

            forAll (mo, oldFaceI)
            {
                w[faceI][oldFaceI] *= normFactor;
            }
        }
    }
}


//- Calculate intersection weights for conservative mapping
void topoPatchMapper::calcIntersectionWeightsAndCentres() const
{
    if (areasPtr_ || centresPtr_)
    {
        FatalErrorIn
        (
            "void topoPatchMapper::"
            "calcIntersectionWeightsAndCentres() const"
        )
            << "Weights already calculated."
            << abort(FatalError);
    }

    // Fetch interpolative addressing
    const labelListList& addr = addressing();

    if (addr.empty())
    {
        // Nothing to map from. Return empty.
        areasPtr_ = new List<scalarList>(0);
        centresPtr_ = new List<vectorList>(0);
        return;
    }

    // Allocate memory
    areasPtr_ = new List<scalarList>(patchSize(), scalarList(0));
    List<scalarList>& a = *areasPtr_;

    centresPtr_ = new List<vectorList>(patchSize(), vectorList(0));
    List<vectorList>& x = *centresPtr_;

    // Obtain stored face-centres
    const vectorField& faceCentres = tMapper_.patchCentres(patch_.index());

    // Fetch maps
    const List<vectorField>& mapFaceCentres = tMapper_.faceCentres();
    const List<scalarField>& mapFaceWeights = tMapper_.faceWeights();

    // Fill in maps first
    const labelList& insertedFaces = insertedObjectLabels();
    const labelList& insertedFacesMap = insertedObjectMap();

    forAll(insertedFaces, faceI)
    {
        a[insertedFaces[faceI]] = mapFaceWeights[insertedFacesMap[faceI]];
        x[insertedFaces[faceI]] = mapFaceCentres[insertedFacesMap[faceI]];
    }

    // Now do mapped faces
    forAll(addr, faceI)
    {
        const labelList& mo = addr[faceI];

        // Check if this is indeed a mapped face
        if (mo.size() == 1 && x[faceI].empty() && a[faceI].empty())
        {
            x[faceI] = vectorList(1, faceCentres[mo[0]]);
            a[faceI] = scalarList(1, 1.0);
        }
    }

    // Final check to ensure everything went okay
    forAll(a, faceI)
    {
        if (a[faceI].empty())
        {
            FatalErrorIn
            (
                "void topoPatchMapper::"
                "calcIntersectionWeightsAndCentres() const"
            )
                << "Weights / centres is missing."
                << nl << " Patch face index: " << faceI
                << abort(FatalError);
        }
    }
}


const List<scalarList>& topoPatchMapper::intersectionWeights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const List<scalarList>& "
            "topoPatchMapper::intersectionWeights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!areasPtr_)
    {
        calcIntersectionWeightsAndCentres();
    }

    return *areasPtr_;
}


const List<vectorList>& topoPatchMapper::intersectionCentres() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const List<vectorList>& "
            "topoPatchMapper::intersectionCentres() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!centresPtr_)
    {
        calcIntersectionWeightsAndCentres();
    }

    return *centresPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
topoPatchMapper::topoPatchMapper
(
    const fvPatch& patch,
    const mapPolyMesh& mpm,
    const topoMapper& mapper
)
:
    patch_(patch),
    mpm_(mpm),
    tMapper_(mapper),
    direct_(false),
    sizeBeforeMapping_(0),
    conservative_(false),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL),
    insertedFaceLabelsPtr_(NULL),
    insertedFaceIndexMapPtr_(NULL),
    insertedFaceAddressingPtr_(NULL),
    areasPtr_(NULL),
    centresPtr_(NULL)
{
    // Compute sizeBeforeMapping.
    // - This needs to be done before insertedObjects
    //   is computed to determine direct mapping
    if (isA<emptyPolyPatch>(patch_.patch()))
    {
        sizeBeforeMapping_ = 0;
    }
    else
    {
        label patchIndex = patch_.index();
        label totalSize = mpm_.oldPatchSizes()[patchIndex];

        // Fetch offset sizes from topoMapper
        const labelListList& sizes = tMapper_.patchSizes();

        // Add offset sizes
        if (sizes.size())
        {
            // Fetch number of physical patches
            label nPhysical = sizes[0].size();

            if (patchIndex < nPhysical)
            {
                forAll(sizes, pI)
                {
                    totalSize += sizes[pI][patchIndex];
                }
            }
        }

        sizeBeforeMapping_ = totalSize;
    }

    // Check for the possibility of direct mapping
    if (insertedObjects())
    {
        direct_ = false;
    }
    else
    {
        direct_ = true;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

topoPatchMapper::~topoPatchMapper()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return the polyPatch size
label topoPatchMapper::patchSize() const
{
    return patch_.patch().size();
}


//- Return size
label topoPatchMapper::size() const
{
    return patch_.size();
}


//- Return size before mapping
label topoPatchMapper::sizeBeforeMapping() const
{
    return sizeBeforeMapping_;
}


//- Is the mapping direct
bool topoPatchMapper::direct() const
{
    return direct_;
}


//- Return direct addressing
const unallocLabelList& topoPatchMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const unallocLabelList& "
            "topoPatchMapper::directAddressing() const"
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
const labelListList& topoPatchMapper::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& "
            "topoPatchMapper::addressing() const"
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
const scalarListList& topoPatchMapper::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList& "
            "topoPatchMapper::weights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (conservative_)
    {
        if (!areasPtr_ && !centresPtr_)
        {
            calcIntersectionWeightsAndCentres();
        }

        return *areasPtr_;
    }

    if (!weightsPtr_)
    {
        calcInverseDistanceWeights();
    }

    return *weightsPtr_;
}


//- Are there any inserted faces
bool topoPatchMapper::insertedObjects() const
{
    return insertedObjectLabels().size();
}


//- Return list of inserted faces
const labelList& topoPatchMapper::insertedObjectLabels() const
{
    if (!insertedFaceLabelsPtr_)
    {
        calcInsertedFaceAddressing();
    }

    return *insertedFaceLabelsPtr_;
}


//- Return addressing map for inserted faces
const labelList& topoPatchMapper::insertedObjectMap() const
{
    if (!insertedFaceIndexMapPtr_)
    {
        calcInsertedFaceAddressing();
    }

    return *insertedFaceIndexMapPtr_;
}


//- Return addressing for inserted faces
const labelListList& topoPatchMapper::insertedFaceAddressing() const
{
    if (!insertedFaceAddressingPtr_)
    {
        calcInsertedFaceAddressing();
    }

    return *insertedFaceAddressingPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void topoPatchMapper::operator=(const topoPatchMapper& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("void topoPatchMapper::operator=")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
