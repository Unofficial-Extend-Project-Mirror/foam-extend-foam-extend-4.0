/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
#include "surfaceMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::fieldValues::fieldValueDelta::applyOperation
(
    const Type& value1,
    const Type& value2
) const
{
    Type result = pTraits<Type>::zero;

    switch (operation_)
    {
        case opAdd:
        {
            result = value1 + value2;
            break;
        }
        case opSubtract:
        {
            result = value1 - value2;
            break;
        }
        case opMin:
        {
            result = min(value1, value2);
            break;
        }
        case opMax:
        {
            result = max(value1, value2);
            break;
        }
        case opAverage:
        {
            result = 0.5*(value1 + value2);
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Type Foam::fieldValues::fieldValueDelta::applyOperation"
                "("
                    "const Type&, "
                    "const Type&"
                ") const"
            )
                << "Unable to process operation "
                << operationTypeNames_[operation_]
                << abort(FatalError);
        }
    }

    return result;
}


template<class Type>
void Foam::fieldValues::fieldValueDelta::processFields(bool& found)
{
    typedef GeometricField<Type, fvPatchField, volMesh> vf;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sf;

    const wordList& fields1 = source1Ptr_->fields();

    const dictionary& results1 = source1Ptr_->resultDict();
    const dictionary& results2 = source2Ptr_->resultDict();

    Type r1(pTraits<Type>::zero);
    Type r2(pTraits<Type>::zero);

    forAll(fields1, i)
    {
        const word& fieldName = fields1[i];

        if
        (
            (obr_.foundObject<vf>(fieldName) || obr_.foundObject<sf>(fieldName))
         && results2.found(fieldName)
        )
        {
            results1.lookup(fieldName) >> r1;
            results2.lookup(fieldName) >> r2;

            Type result = applyOperation(r1, r2);

            if (log_) Info<< "    " << operationTypeNames_[operation_]
                << "(" << fieldName << ") = " << result
                << endl;

            if (Pstream::master())
            {
                file()<< tab << result;
            }

            found = true;
        }
    }
}


// ************************************************************************* //
