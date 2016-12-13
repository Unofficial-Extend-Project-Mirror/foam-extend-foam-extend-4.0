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

#include "GeometricField.H"
#include "volMesh.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::solutionControl::storePrevIter() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    HashTable<const GeoField*>
        flds(mesh_.objectRegistry::lookupClass<GeoField>());

    forAllIter(typename HashTable<const GeoField*>, flds, iter)
    {
        const GeoField& fld = *iter();

        const word& fName = fld.name();

        size_t prevIterField = fName.find("PrevIter");

        if
        (
            (prevIterField == word::npos)
         && mesh_.solutionDict().relaxField(fName)
        )
        {
            if (debug)
            {
                Info<< algorithmName_ << ": storing previous iter for "
                    << fName << endl;
            }

            fld.storePrevIter();
        }
    }
}


// ************************************************************************* //
