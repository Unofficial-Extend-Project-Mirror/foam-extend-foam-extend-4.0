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
    topoPatchMapper

Description
    Implementation of the topoPatchMapper class

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*----------------------------------------------------------------------------*/

#include "IOmanip.H"
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
    const label oldPatchSize = mpm_.oldPatchSizes()[patch_.index()];
    const label oldPatchStart = mpm_.oldPatchStarts()[patch_.index()];
    const label oldPatchEnd = oldPatchStart + oldPatchSize;

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

                // Renumber addressing to patch.
                // Also, check mapping for hits into
                // other patches / internal faces.
                addr = fffI.masterObjects();

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
                        FatalErrorIn
                        (
                            "void topoPatchMapper::"
                            "calcInsertedFaceAddressing() const"
                        )
                            << "Addressing into another patch is not allowed."
                            << nl << " Patch face index: " << faceI
                            << nl << " addr[faceI]: " << addr[faceI]
                            << nl << " oldPatchStart: " << oldPatchStart
                            << nl << " oldPatchSize: " << oldPatchSize
                            << nl << " oldPatchEnd: " << oldPatchEnd
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
    const label oldPatchSize = mpm_.oldPatchSizes()[patch_.index()];
    const label oldPatchStart = mpm_.oldPatchStarts()[patch_.index()];
    const label oldPatchEnd = oldPatchStart + oldPatchSize;

    // Assemble the maps: slice to patch
    if (direct())
    {
        // Direct mapping - slice to size
        directAddrPtr_ = new labelList(patch_.patch().patchSlice(mpm_.faceMap()));

        labelList& addr = *directAddrPtr_;

        // Shift to local patch indices.
        // Also, check mapping for hits into other patches / internal faces.
        forAll (addr, faceI)
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
                FatalErrorIn
                (
                    "void topoPatchMapper::calcAddressing() const"
                )
                    << "Addressing into another patch is not allowed."
                    << nl << " Patch face index: " << faceI
                    << nl << " addr[faceI]: " << addr[faceI]
                    << nl << " oldPatchStart: " << oldPatchStart
                    << nl << " oldPatchSize: " << oldPatchSize
                    << nl << " oldPatchEnd: " << oldPatchEnd
                    << abort(FatalError);
            }
        }
    }
    else
    {
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
        const labelList& fm = patch_.patch().patchSlice(mpm_.faceMap());

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
                FatalErrorIn
                (
                    "void topoPatchMapper::calcAddressing() const"
                )
                    << "Addressing is missing." << nl
                    << " Patch face index: " << faceI << nl
                    << " nInsertedFaces: " << insertedFaces.size() << nl
                    << " faceMap: " << fm[faceI] << nl
                    << " Patch: " << patch_.name() << nl
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

    // Allocate memory
    areasPtr_ = new List<scalarField>(patchSize(), scalarField(0));
    List<scalarField>& a = *areasPtr_;

    centresPtr_ = new List<vectorField>(patchSize(), vectorField(0));
    List<vectorField>& x = *centresPtr_;

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
            x[faceI] = vectorField(1, faceCentres[mo[0]]);
            a[faceI] = scalarField(1, 1.0);
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


const List<scalarField>& topoPatchMapper::intersectionWeights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const List<scalarField>& "
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


const List<vectorField>& topoPatchMapper::intersectionCentres() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const List<vectorField>& "
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
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL),
    insertedFaceLabelsPtr_(NULL),
    insertedFaceIndexMapPtr_(NULL),
    insertedFaceAddressingPtr_(NULL),
    areasPtr_(NULL),
    centresPtr_(NULL)
{
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
    if (patch_.type() == "empty")
    {
        return 0;
    }

    return mpm_.oldPatchSizes()[patch_.index()];
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


//- Map the patch field
template <class Type>
void topoPatchMapper::mapPatchField
(
    const word& fieldName,
    Field<Type>& pF
) const
{
    // To invoke inverse-distance weighting, use this:
    // pF.autoMap(*this);

    // Check for possibility of direct mapping
    if (direct())
    {
        pF.autoMap(*this);
    }
    else
    {
        if (pF.size() != sizeBeforeMapping())
        {
            FatalErrorIn
            (
                "\n\n"
                "void topoCellMapper::mapPatchField<Type>\n"
                "(\n"
                "    const word& fieldName,\n"
                "    Field<Type>& iF\n"
                ") const\n"
            )  << "Incompatible size before mapping." << nl
               << " Field: " << fieldName << nl
               << " Field size: " << pF.size() << nl
               << " map size: " << sizeBeforeMapping() << nl
               << abort(FatalError);
        }

        // Fetch addressing
        const labelListList& pAddressing = addressing();
        const List<scalarField>& wF = intersectionWeights();

        // Compute the integral of the source field
        Type intSource = sum(pF * tMapper_.patchAreas(patch_.index()));

        // Copy the original field
        Field<Type> fieldCpy(pF);

        // Resize to current dimensions
        pF.setSize(size());

        // Map the patch field
        forAll(pF, faceI)
        {
            const labelList& addr = pAddressing[faceI];

            pF[faceI] = pTraits<Type>::zero;

            // Accumulate area-weighted interpolate
            forAll(addr, faceJ)
            {
                pF[faceI] +=
                (
                    wF[faceI][faceJ] * fieldCpy[addr[faceJ]]
                );
            }
        }

        // Compute the integral of the target field
        const polyPatch& ppI = mpm_.mesh().boundaryMesh()[patch_.index()];

        Type intTarget = sum(pF * mag(ppI.faceAreas()));

        if (polyMesh::debug)
        {
            int oldP = Info().precision();

            // Compare the global integral
            Info << " Field : " << fieldName
                 << " Patch : " << ppI.name()
                 << " integral errors : "
                 << setprecision(15)
                 << " source : " << mag(intSource)
                 << " target : " << mag(intTarget)
                 << " norm : "
                 << (mag(intTarget - intSource) / (mag(intSource) + VSMALL))
                 << setprecision(oldP)
                 << endl;
        }
    }
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
