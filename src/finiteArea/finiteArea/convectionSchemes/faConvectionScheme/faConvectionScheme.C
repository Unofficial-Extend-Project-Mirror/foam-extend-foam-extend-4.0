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
    Abstract base class for finite area calculus convection schemes.

\*---------------------------------------------------------------------------*/

#include "faConvectionScheme.H"
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
tmp<convectionScheme<Type> > convectionScheme<Type>::New
(
    const faMesh& mesh,
    const edgeScalarField& faceFlux,
    Istream& schemeData
)
{
    if (fa::debug)
    {
        Info<< "convectionScheme<Type>::New"
               "(const faMesh&, const edgeScalarField&, Istream&) : "
               "constructing convectionScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorIn
        (
            "convectionScheme<Type>::New"
            "(const faMesh&, const edgeScalarField&, Istream&)",
            schemeData
        )   << "Convection scheme not specified" << endl << endl
            << "Valid convection schemes are :" << endl
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
            "convectionScheme<Type>::New"
            "(const faMesh&, const edgeScalarField&, Istream&)",
            schemeData
        )   << "unknown convection scheme " << schemeName << endl << endl
            << "Valid convection schemes are :" << endl
            << IstreamConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, faceFlux, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
convectionScheme<Type>::~convectionScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
