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

#include "lagrangianFieldDecomposer.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void lagrangianFieldDecomposer::readFields
(
    const label cloudI,
    const IOobjectList& lagrangianObjects,
    PtrList<PtrList<IOField<Type> > >& lagrangianFields
)
{
    // Search list of objects for lagrangian fields
    IOobjectList lagrangianTypeObjects
    (
        lagrangianObjects.lookupClass(IOField<Type>::typeName)
    );

    lagrangianFields.set
    (
        cloudI,
        new PtrList<IOField<Type> >
        (
            lagrangianTypeObjects.size()
        )
    );

    label lagrangianFieldi=0;
    forAllIter(IOobjectList, lagrangianTypeObjects, iter)
    {
        lagrangianFields[cloudI].set
        (
            lagrangianFieldi++,
            new IOField<Type>(*iter())
        );
    }
}


template<class Type>
tmp<IOField<Type> > lagrangianFieldDecomposer::decomposeField
(
    const word& cloudName,
    const IOField<Type>& field
) const
{
    // Create and map the internal field values
    Field<Type> procField(field, particleIndices_);

    // Create the field for the processor
    return tmp<IOField<Type> >
    (
        new IOField<Type>
        (
            IOobject
            (
                field.name(),
                procMesh_.time().timeName(),
                cloud::prefix/cloudName,
                procMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procField
        )
    );
}


template<class GeoField>
void lagrangianFieldDecomposer::decomposeFields
(
    const word& cloudName,
    const PtrList<GeoField>& fields
) const
{
    if (particleIndices_.size())
    {
        forAll (fields, fieldI)
        {
            decomposeField(cloudName, fields[fieldI])().write();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
