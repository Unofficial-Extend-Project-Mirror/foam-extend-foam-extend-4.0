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

Description
    Abstract base class for finite volume calculus div schemes.

\*---------------------------------------------------------------------------*/

#include "fv.H"
#include "HashTable.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<divScheme<Type> > divScheme<Type>::New
(
    const fvMesh& mesh,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        Info<< "divScheme<Type>::New(const fvMesh&, Istream&) : "
               "constructing divScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorIn
        (
            "divScheme<Type>::New(const fvMesh&, Istream&)",
            schemeData
        )   << "Div scheme not specified" << endl << endl
            << "Valid div schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(schemeName);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "divScheme<Type>::New(const fvMesh&, Istream&)",
            schemeData
        )   << "Unknown div scheme "
            << schemeName << nl << nl
            << "Valid div schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
divScheme<Type>::~divScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
>
divScheme<Type>::fvmUDiv
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    FatalErrorIn
    (
        "tmp<BlockLduSystem> divScheme<Type>::fvmDiv\n"
        "(\n"
        "    const GeometricField<Type, fvPatchField, volMesh>&"
        ")\n"
    )   << "Implicit div operator currently defined only for Gauss linear. "
        << abort(FatalError);

    typedef typename innerProduct<vector, Type>::type DivType;

    tmp<BlockLduSystem<vector, DivType> > tbs
    (
        new BlockLduSystem<vector, DivType>(vf.mesh())
    );

    return tbs;
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
>
divScheme<Type>::fvmUDiv
(
    const surfaceScalarField& flux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    FatalErrorIn
    (
        "tmp<BlockLduSystem> divScheme<Type>::fvmDiv\n"
        "(\n"
        "    surfaceScalarField&"
        "    const GeometricField<Type, fvPatchField, volMesh>&"
        ")\n"
    )   << "Implicit div operator currently defined only for Gauss linear. "
        << abort(FatalError);

    typedef typename innerProduct<vector, Type>::type DivType;

    tmp<BlockLduSystem<vector, DivType> > tbs
    (
        new BlockLduSystem<vector, DivType>(vf.mesh())
    );

    return tbs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
