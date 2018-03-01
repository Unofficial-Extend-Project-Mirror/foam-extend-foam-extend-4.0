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

#include "oversetMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::oversetMesh::interpolate
(
    Field<Type>& accF,
    const Field<Type>& cellF,
    const word& fieldName
) const
{
    // Check sizes
    if
    (
        accF.size() != this->acceptorCells().size()
     || cellF.size() != this->mesh().nCells()
    )
    {
        FatalErrorIn
        (
            "void oversetMesh::donorToAcceptor\n"
            "(\n"
            "    Field<Type>& accF,\n"
            "    const Field<Type>& ,cellF\n"
            ") const"
        )   << "Size of fields does not correspond to interpolation" << nl
            << "Source field size = " << cellF.size()
            << " mesh size = " << this->mesh().nCells()
            << " target field size = " << accF.size()
            << " acceptor list size = " << this->acceptorCells().size()
            << abort(FatalError);
    }

    // Create a copy of the field to interpolate and distribute its donor data
    // across processors.
    Field<Type> donorField = cellF;

    // Note: mapDistribute::distribute changes the size of the field it operates
    // on to correspond to the size of received donor data
    map().distribute(donorField);

    // Get remote donor to local acceptor map. Allows easy indexing of donor
    // that came from a (possibly) remote processor via its processor number and
    // cell number.
    const List<labelField>& remoteDAAddressing =
        remoteDonorToLocalAcceptorAddr();

    // Get overset interpolation for this field
    const oversetInterpolation& interpolation =
        this->interpolationScheme(fieldName);

    // Get interpolation weights for all donors for a given local acceptor in a
    // given region
    const oversetInterpolation::ListScalarFieldField& weights =
        interpolation.weights();

    // Note: accF field indexed by number of acceptors (region-wise
    // incremental, see oversetMesh::calcCellClassification() for details)
    label aI = 0;

    // Loop through all regions to account for all acceptors
    forAll (regions_, regionI)
    {
        // Get acceptors for this region
        const donorAcceptorList& curAcceptors = regions_[regionI].acceptors();

        // Get weights for this region
        const oversetInterpolation::ScalarFieldField& regionWeights =
            weights[regionI];

        // Loop through all acceptors of this region
        forAll (curAcceptors, regionAccI)
        {
            // Get necessary acceptor information
            const donorAcceptor& curDA = curAcceptors[regionAccI];

            // Get necessary donor information
            const label donorCellI = curDA.donorCell();

            // Get weights for this acceptor
            const scalarField& w = regionWeights[regionAccI];

            // Get remote donor to local acceptor addressing for this donor
            // processor. Note: it is assumed that all donors for an acceptor
            // come from the same processor (though it can be remote processor)
            const labelField& donorProcAddr =
                remoteDAAddressing[curDA.donorProcNo()];

            // Combine donor contributions for this acceptor
            // Note that accF is addressed by aI instead of acceptor cell
            // indices. See oversetMesh::calcCellClassification for details
            accF[aI] = pTraits<Type>::zero;

            // Master donor contribution (first entry of weights corresponds to
            // the weight for master donor)
            accF[aI] += w[0]*donorField[donorProcAddr[donorCellI]];

            // Neighbouring donor contributions
            const donorAcceptor::DynamicLabelList& nbrDonors =
                curDA.extendedDonorCells();

            forAll (nbrDonors, nbrDonorI)
            {
                // Get extended neighbour cell label
                const label nbrDonorCellI = nbrDonors[nbrDonorI];

                // Note indexing weights with + 1 offset since the first entry
                // is for the master donor and is already accounted for
                accF[aI] += w[nbrDonorI + 1]
                   *donorField[donorProcAddr[nbrDonorCellI]];
            }

            // Increment acceptor index
            ++aI;
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::oversetMesh::interpolate
(
    const Field<Type>& cellF,
    const word& fieldName
) const
{
    tmp<Field<Type> > tresult(new Field<Type>(this->acceptorCells().size()));
    Field<Type>& result = tresult();

    interpolate(result, cellF, fieldName);

    return tresult;
}


// ************************************************************************* //
