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
    Abstract base class for finite volume calculus adjConvection schemes.

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
adjConvectionScheme<Type>::adjConvectionScheme(const adjConvectionScheme& cs)
:
    refCount(),
    mesh_(cs.mesh_)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<adjConvectionScheme<Type> > adjConvectionScheme<Type>::New
(
    const fvMesh& mesh,
    const volVectorField& Up,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        Info<< "adjConvectionScheme<Type>::New"
               "(const fvMesh&, const volVectorField&, Istream&) : "
               "constructing adjConvectionScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorIn
        (
            "adjConvectionScheme<Type>::New"
            "(const fvMesh&, const volVectorField&, Istream&)",
            schemeData
        )   << "AdjConvection scheme not specified" << endl << endl
            << "Valid adjConvection schemes are :" << endl
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
            "adjConvectionScheme<Type>::New"
            "(const fvMesh&, const volVectorField&, Istream&)",
            schemeData
        )   << "unknown adjConvection scheme " << schemeName << endl << endl
            << "Valid adjConvection schemes are :" << endl
            << IstreamConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, Up, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
adjConvectionScheme<Type>::~adjConvectionScheme()
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void adjConvectionScheme<Type>::operator=(const adjConvectionScheme<Type>& cs)
{
    if (this == &cs)
    {
        FatalErrorIn
        (
            "adjConvectionScheme<Type>::operator="
            "(const adjConvectionScheme<Type>&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
