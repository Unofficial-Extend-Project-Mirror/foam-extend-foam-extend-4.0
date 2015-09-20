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

\*---------------------------------------------------------------------------*/

#include "readFields.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class GeoField>
void readFields
(
    const vtkMesh& vMesh,
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields,
    PtrList<GeoField>& fields
)
{
    // Search list of objects for volScalarFields
    IOobjectList fieldObjects(objects.lookupClass(GeoField::typeName));

    // Construct the vol scalar fields
    label nFields = fields.size();
    fields.setSize(nFields + fieldObjects.size());

    for
    (
        IOobjectList::iterator iter = fieldObjects.begin();
        iter != fieldObjects.end();
        ++iter
    )
    {
        if (selectedFields.empty() || selectedFields.found(iter()->name()))
        {
            fields.set
            (
                nFields,
                vMesh.interpolate
                (
                    GeoField
                    (
                        *iter(),
                        mesh
                    )
                )
            );
            nFields++;
        }
    }

    fields.setSize(nFields);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
