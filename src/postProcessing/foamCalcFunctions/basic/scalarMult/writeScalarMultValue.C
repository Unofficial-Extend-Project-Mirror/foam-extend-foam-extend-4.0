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

template<class Type>
void Foam::calcTypes::scalarMult::writeScalarMultValue
(
    const IOobject& baseHeader,
    const string& valueStr,
    const fvMesh& mesh,
    bool& processed
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (baseHeader.headerClassName() == fieldType::typeName)
    {
        if (resultName_ == "")
        {
            resultName_ = baseHeader.name() + "_scalarMult";
        }

        scalar value;
        IStringStream(valueStr)() >> value;

        Info<< "    Reading " << baseHeader.name() << endl;
        fieldType baseField(baseHeader, mesh);

        fieldType newField
        (
            IOobject
            (
                resultName_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ
            ),
            baseField
        );

        Info<< "    Calculating " << resultName_ << endl;

        newField == value*baseField;
        newField.write();

        processed = true;
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
