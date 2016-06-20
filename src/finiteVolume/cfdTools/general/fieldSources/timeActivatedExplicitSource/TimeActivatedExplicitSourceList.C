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

#include "TimeActivatedExplicitSourceList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::TimeActivatedExplicitSourceList<Type>::TimeActivatedExplicitSourceList
(
    const word& name,
    const fvMesh& mesh,
    const dimensionSet& dimensions,
    const wordList& fieldNames
)
:
    IOPtrList<TimeActivatedExplicitSource<Type> >
    (
        IOobject
        (
            name + "SourceProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        typename TimeActivatedExplicitSource<Type>::iNew(mesh, fieldNames)
    ),
    name_(name),
    mesh_(mesh),
    dimensions_(dimensions),
    fieldNames_(fieldNames)
{}


template<class Type>
Foam::TimeActivatedExplicitSourceList<Type>::TimeActivatedExplicitSourceList
(
    const word& name,
    const fvMesh& mesh,
    const dimensionSet& dimensions,
    const word& fieldName
)
:
    IOPtrList<TimeActivatedExplicitSource<Type> >
    (
        IOobject
        (
            name + "SourceProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        typename TimeActivatedExplicitSource<Type>::iNew
        (
            mesh,
            IStringStream('(' + fieldName + ')')()
        )
    ),
    name_(name),
    mesh_(mesh),
    dimensions_(dimensions),
    fieldNames_(1, fieldName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh> >
Foam::TimeActivatedExplicitSourceList<Type>::Su(const label fieldI)
{
    tmp<DimensionedField<Type, volMesh> > tSu
    (
       new DimensionedField<Type, volMesh>
       (
           IOobject
           (
               name_ + "Source_" + fieldNames_[fieldI],
               mesh_.time().timeName(),
               mesh_,
               IOobject::NO_READ,
               IOobject::NO_WRITE
           ),
           mesh_,
           dimensioned<Type>("zero", dimensions_, pTraits<Type>::zero)
       )
    );

    DimensionedField<Type, volMesh>& Su = tSu();

    forAll(*this, i)
    {
        this->operator[](i).addToField(Su, fieldI);
    }

    return tSu;
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh> >
Foam::TimeActivatedExplicitSourceList<Type>::SuTot()
{
    tmp<DimensionedField<Type, volMesh> > tSuTot
    (
       new DimensionedField<Type, volMesh>
       (
           IOobject
           (
               name_ + "TotalSource",
               mesh_.time().timeName(),
               mesh_,
               IOobject::NO_READ,
               IOobject::NO_WRITE
           ),
           mesh_,
           dimensioned<Type>("zero", dimensions_, pTraits<Type>::zero)
       )
    );

    DimensionedField<Type, volMesh>& SuTot = tSuTot();

    forAll(fieldNames_, fieldI)
    {
        forAll(*this, sourceI)
        {
            this->operator[](sourceI).addToField(SuTot, fieldI);
        }
    }

    return tSuTot;
}


template<class Type>
bool Foam::TimeActivatedExplicitSourceList<Type>::readData(Istream& is)
{
    this->clear();

    IOPtrList<TimeActivatedExplicitSource<Type> > newSources
    (
        IOobject
        (
            name_ + "TimeActivatedExplicitSource",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        typename TimeActivatedExplicitSource<Type>::iNew(mesh_, fieldNames_)
    );

    this->transfer(newSources);

    return is.good();
}


template<class Type>
bool Foam::TimeActivatedExplicitSourceList<Type>::writeData(Ostream& os) const
{
    // Write size of list
    os << nl << this->size();

    // Write beginning of contents
    os << nl << token::BEGIN_LIST;

    // Write list contents
    forAll(*this, i)
    {
        os << nl;
        this->operator[](i).writeData(os);
    }

    // Write end of contents
    os << token::END_LIST << token::END_STATEMENT << nl;

    // Check state of IOstream
    return os.good();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const TimeActivatedExplicitSourceList<Type>& sources
)
{
    sources.writeData(os);
    return os;
}


// ************************************************************************* //
