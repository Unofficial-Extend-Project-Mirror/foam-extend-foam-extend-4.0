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

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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


} // End namespace Foam

// ************************************************************************* //
