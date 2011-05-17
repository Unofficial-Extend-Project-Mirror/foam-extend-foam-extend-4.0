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

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Conservatively map the internal field
template <class Type, class gradType>
void topoCellMapper::mapInternalField
(
    const word& fieldName,
    const Field<gradType>& gF,
    Field<Type>& iF
) const
{
    if (iF.size() != sizeBeforeMapping() || gF.size() != sizeBeforeMapping())
    {
        FatalErrorIn
        (
            "\n\n"
            "void topoCellMapper::mapInternalField<Type>\n"
            "(\n"
            "    const word& fieldName,\n"
            "    const Field<gradType>& gF,\n"
            "    Field<Type>& iF\n"
            ") const\n"
        )  << "Incompatible size before mapping." << nl
           << " Field: " << fieldName << nl
           << " Field size: " << iF.size() << nl
           << " Gradient Field size: " << gF.size() << nl
           << " map size: " << sizeBeforeMapping() << nl
           << abort(FatalError);
    }

    // If we have direct addressing, map and bail out
    if (direct())
    {
        iF.autoMap(*this);
        return;
    }

    // Fetch addressing
    const labelListList& cAddressing = addressing();
    const List<scalarField>& wC = intersectionWeights();
    const List<vectorField>& xC = intersectionCentres();

    // Fetch geometry
    const vectorField& centres = tMapper_.internalCentres();

    // Copy the original field
    Field<Type> fieldCpy(iF);

    // Resize to current dimensions
    iF.setSize(size());

    // Map the internal field
    forAll(iF, cellI)
    {
        const labelList& addr = cAddressing[cellI];

        iF[cellI] = pTraits<Type>::zero;

        forAll(addr, cellJ)
        {
            // Accumulate volume-weighted Taylor-series interpolate
            iF[cellI] +=
            (
                wC[cellI][cellJ] *
                (
                    fieldCpy[addr[cellJ]]
                  + (
                        gF[addr[cellJ]] &
                        (
                            xC[cellI][cellJ]
                          - centres[addr[cellJ]]
                        )
                    )
                )
            );
        }
    }
}


} // End namespace Foam

// ************************************************************************* //
