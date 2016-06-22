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
    topoCellMapper

Description
    Implementation of the topoCellMapper class

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "topoMapper.H"
#include "mapPolyMesh.H"
#include "topoCellMapper.H"
#include "demandDrivenData.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Clear out local storage
void topoCellMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);
    deleteDemandDrivenData(insertedCellLabelsPtr_);
    deleteDemandDrivenData(volumesPtr_);
    deleteDemandDrivenData(centresPtr_);
}


//- Calculate addressing for interpolative mapping
void topoCellMapper::calcAddressing() const
{
    if (directAddrPtr_ || interpolationAddrPtr_)
    {
        FatalErrorIn("void topoCellMapper::calcAddressing() const")
            << "Addressing already calculated."
            << abort(FatalError);
    }

    // Allocate for inserted cell labels
    label nInsertedCells = 0;

    insertedCellLabelsPtr_ = new labelList(mesh_.nCells(), -1);
    labelList& insertedCells = *insertedCellLabelsPtr_;

    if (direct())
    {
        // Direct addressing, no weights
        directAddrPtr_ = new labelList(mpm_.cellMap());
    }
    else
    {
        // Interpolative addressing
        interpolationAddrPtr_ = new labelListList(mesh_.nCells());
        labelListList& addr = *interpolationAddrPtr_;

        const List<objectMap>& cfc = mpm_.cellsFromCellsMap();

        forAll (cfc, cfcI)
        {
            // Get addressing
            const labelList& mo = cfc[cfcI].masterObjects();

            label cellI = cfc[cfcI].index();

            if (addr[cellI].size() > 0)
            {
                FatalErrorIn("void topoCellMapper::calcAddressing() const")
                    << "Master cell " << cellI
                    << " mapped from cells " << mo
                    << " is already destination for mapping."
                    << abort(FatalError);
            }

            // Set master objects
            addr[cellI] = mo;
        }

        // Do mapped cells. Note that this can already be set by cellsFromCells
        // so check if addressing size still zero.
        const labelList& cm = mpm_.cellMap();

        forAll (cm, cellI)
        {
            // Mapped from a single cell
            if (cm[cellI] > -1 && addr[cellI].empty())
            {
                addr[cellI] = labelList(1, cm[cellI]);
            }

            // Check for inserted cells without any addressing
            if (cm[cellI] < 0 && addr[cellI].empty())
            {
                insertedCells[nInsertedCells++] = cellI;
            }
        }
    }

    // Shorten inserted cells to actual size
    insertedCells.setSize(nInsertedCells);

    if (nInsertedCells)
    {
        FatalErrorIn("void topoCellMapper::calcAddressing() const")
            << " Found " << nInsertedCells << " which are"
            << " not mapped from any parent cells." << nl
            << " List: " << nl
            << insertedCells
            << abort(FatalError);
    }
}


//- Calculate inverse-distance weights for interpolative mapping
void topoCellMapper::calcInverseDistanceWeights() const
{
    if (weightsPtr_)
    {
        FatalErrorIn
        (
            "void topoCellMapper::calcInverseDistanceWeights() const"
        )
            << "Weights already calculated."
            << abort(FatalError);
    }

    // Fetch interpolative addressing
    const labelListList& addr = addressing();

    // Allocate memory
    weightsPtr_ = new scalarListList(size());
    scalarListList& w = *weightsPtr_;

    // Obtain cell-centre information from old/new meshes
    const vectorField& oldCentres = tMapper_.internalCentres();
    const vectorField& newCentres = mesh_.cellCentres();

    forAll(addr, cellI)
    {
        const labelList& mo = addr[cellI];

        // Do mapped cells
        if (mo.size() == 1)
        {
            w[cellI] = scalarList(1, 1.0);
        }
        else
        {
            // Map from masters, inverse-distance weights
            scalar totalWeight = 0.0;
            w[cellI] = scalarList(mo.size(), 0.0);

            forAll (mo, oldCellI)
            {
                w[cellI][oldCellI] =
                (
                    1.0/stabilise
                    (
                        magSqr
                        (
                            newCentres[cellI]
                          - oldCentres[mo[oldCellI]]
                        ),
                        VSMALL
                    )
                );

                totalWeight += w[cellI][oldCellI];
            }

            // Normalize weights
            scalar normFactor = (1.0/totalWeight);

            forAll (mo, oldCellI)
            {
                w[cellI][oldCellI] *= normFactor;
            }
        }
    }
}


//- Calculate intersection weights for conservative mapping
void topoCellMapper::calcIntersectionWeightsAndCentres() const
{
    if (volumesPtr_ || centresPtr_)
    {
        FatalErrorIn
        (
            "void topoCellMapper::"
            "calcIntersectionWeightsAndCentres() const"
        )
            << "Weights already calculated."
            << abort(FatalError);
    }

    // Fetch interpolative addressing
    const labelListList& addr = addressing();

    // Allocate memory
    volumesPtr_ = new List<scalarField>(size(), scalarField(0));
    List<scalarField>& v = *volumesPtr_;

    centresPtr_ = new List<vectorField>(size(), vectorField(0));
    List<vectorField>& x = *centresPtr_;

    // Obtain stored cell-centres
    const vectorField& cellCentres = tMapper_.internalCentres();

    // Fetch maps
    const List<objectMap>& cfc = mpm_.cellsFromCellsMap();
    const List<vectorField>& mapCellCentres = tMapper_.cellCentres();
    const List<scalarField>& mapCellWeights = tMapper_.cellWeights();

    // Fill in maps first
    forAll(cfc, indexI)
    {
        x[cfc[indexI].index()] = mapCellCentres[indexI];
        v[cfc[indexI].index()] = mapCellWeights[indexI];
    }

    // Now do mapped cells
    forAll(addr, cellI)
    {
        const labelList& mo = addr[cellI];

        // Check if this is indeed a mapped cell
        if (mo.size() == 1 && x[cellI].empty() && v[cellI].empty())
        {
            x[cellI] = vectorField(1, cellCentres[mo[0]]);
            v[cellI] = scalarField(1, 1.0);
        }
    }
}


const List<scalarField>& topoCellMapper::intersectionWeights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const List<scalarField>& "
            "topoCellMapper::intersectionWeights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!volumesPtr_)
    {
        calcIntersectionWeightsAndCentres();
    }

    return *volumesPtr_;
}


const List<vectorField>& topoCellMapper::intersectionCentres() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const List<vectorField>& "
            "topoCellMapper::intersectionCentres() const"
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
topoCellMapper::topoCellMapper
(
    const mapPolyMesh& mpm,
    const topoMapper& mapper
)
:
    mesh_(mpm.mesh()),
    mpm_(mpm),
    tMapper_(mapper),
    direct_(false),
    sizeBeforeMapping_(mpm.nOldCells()),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL),
    insertedCellLabelsPtr_(NULL),
    volumesPtr_(NULL),
    centresPtr_(NULL)
{
    // Fetch offset sizes from topoMapper
    const labelList& sizes = tMapper_.cellSizes();

    // Add offset sizes
    if (sizes.size())
    {
        forAll(sizes, pI)
        {
            sizeBeforeMapping_ += sizes[pI];
        }
    }

    // Check for possibility of direct mapping
    if
    (
        (min(mpm_.cellMap()) > -1)
     && mpm_.cellsFromPointsMap().empty()
     && mpm_.cellsFromEdgesMap().empty()
     && mpm_.cellsFromFacesMap().empty()
     && mpm_.cellsFromCellsMap().empty()
    )
    {
        direct_ = true;
    }
    else
    {
        direct_ = false;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

topoCellMapper::~topoCellMapper()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return size
label topoCellMapper::size() const
{
    return mesh_.nCells();
}


//- Return size before mapping
label topoCellMapper::sizeBeforeMapping() const
{
    return sizeBeforeMapping_;
}


//- Is the mapping direct
bool topoCellMapper::direct() const
{
    return direct_;
}


//- Return direct addressing
const unallocLabelList& topoCellMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const unallocLabelList& "
            "topoCellMapper::directAddressing() const"
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
const labelListList& topoCellMapper::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& "
            "topoCellMapper::addressing() const"
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
const scalarListList& topoCellMapper::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList& "
            "topoCellMapper::weights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!weightsPtr_)
    {
        calcInverseDistanceWeights();
    }

    return *weightsPtr_;
}


//- Are there any inserted cells
bool topoCellMapper::insertedObjects() const
{
    return insertedObjectLabels().size();
}


//- Return list of inserted cells
const labelList& topoCellMapper::insertedObjectLabels() const
{
    if (!insertedCellLabelsPtr_)
    {
        calcAddressing();
    }

    return *insertedCellLabelsPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void topoCellMapper::operator=(const topoCellMapper& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("void topoCellMapper::operator=")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
