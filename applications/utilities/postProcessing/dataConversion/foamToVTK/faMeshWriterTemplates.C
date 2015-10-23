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

#include "faMeshWriter.H"
#include "writeFuns.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::faMeshWriter::write
(
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    os_ << fld.name() << ' ' << pTraits<Type>::nComponents << ' '
        << aMesh_.nFaces() << " float" << std::endl;

    DynamicList<floatScalar> fField
    (
        pTraits<Type>::nComponents*aMesh_.nFaces()
    );

    writeFuns::insert(fld.internalField(), fField);
    writeFuns::write(os_, binary_, fField);
}


template<class Type>
void Foam::faMeshWriter::write
(
    const PtrList<GeometricField<Type, faPatchField, areaMesh> >& flds
)
{
    forAll(flds, fieldI)
    {
        write(flds[fieldI]);
    }
}


template<class Type>
void Foam::faMeshWriter::write
(
    const PrimitivePatchInterpolation<indirectPrimitivePatch>& pInter,
    const PtrList<GeometricField<Type, faPatchField, areaMesh> >& flds
)
{
    forAll(flds, fieldI)
    {
        write(pInter, flds[fieldI]);
    }
}


template<class Type>
void Foam::faMeshWriter::write
(
    const PrimitivePatchInterpolation<indirectPrimitivePatch>& pInter,
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    os_ << fld.name() << ' ' << pTraits<Type>::nComponents << ' '
        << aMesh_.nPoints() << " float" << std::endl;

    DynamicList<floatScalar> fField
    (
        pTraits<Type>::nComponents*aMesh_.nPoints()
    );

    writeFuns::insert
    (
        pInter.faceToPointInterpolate(fld.internalField())(),
        fField
    );
    writeFuns::write(os_, binary_, fField);
}


// ************************************************************************* //
