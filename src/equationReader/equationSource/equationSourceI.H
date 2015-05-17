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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
inline const UPtrList<const Type>& equationSource<Type>::singles() const
{
    return singles_;
}


template<class Type>
inline const wordList& equationSource<Type>::singleNames() const
{
    return singleNames_;
}


template<class Type>
inline const PtrList<dimensionSet>&
    equationSource<Type>::singleDimensions() const
{
    return singleDimensions_;
}


template<class Type>
inline const PtrList<UPtrList<const Field<Type> > >&
    equationSource<Type>::fields() const
{
    return fields_;
}


template<class Type>
inline const PtrList<dimensionSet>&
    equationSource<Type>::fieldDimensions() const
{
    return fieldDimensions_;
}


template<class Type>
inline const wordList& equationSource<Type>::fieldNames() const
{
    return fieldNames_;
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
