/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Abstract base class for finite area calculus laplacian schemes.

\*---------------------------------------------------------------------------*/

#include "fa.H"
#include "HashTable.H"
#include "linearEdgeInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<laplacianScheme<Type> > laplacianScheme<Type>::New
(
    const faMesh& mesh,
    Istream& schemeData
)
{
    if (fa::debug)
    {
        Info<< "laplacianScheme<Type>::New(const faMesh&, Istream&) : "
               "constructing laplacianScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorIn
        (
            "laplacianScheme<Type>::New(const faMesh&, Istream&)",
            schemeData
        )   << "Laplacian scheme not specified" << endl << endl
            << "Valid laplacian schemes are :" << endl
            << IstreamConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    word schemeName(schemeData);

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(schemeName);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "laplacianScheme<Type>::New(const faMesh&, Istream&)",
            schemeData
        )   << "unknown laplacian scheme " << schemeName << endl << endl
            << "Valid laplacian schemes are :" << endl
            << IstreamConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
laplacianScheme<Type>::~laplacianScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<faMatrix<Type> >
laplacianScheme<Type>::famLaplacian
(
    const areaScalarField& gamma,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return famLaplacian(tinterpGammaScheme_().interpolate(gamma)(), vf);
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
laplacianScheme<Type>::facLaplacian
(
    const areaScalarField& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return facLaplacian(tinterpGammaScheme_().interpolate(gamma)(), vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
