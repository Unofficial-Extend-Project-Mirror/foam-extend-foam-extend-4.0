/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    Abstract base class for finite volume calculus gradient schemes.

\*---------------------------------------------------------------------------*/

#include "fv.H"
#include "HashTable.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<gradScheme<Type> > gradScheme<Type>::New
(
    const fvMesh& mesh,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        Info<< "gradScheme<Type>::New(Istream& schemeData) : "
               "constructing gradScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorIn
        (
            "gradScheme<Type>::New(Istream& schemeData)",
            schemeData
        )   << "Grad scheme not specified" << endl << endl
            << "Valid grad schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    word schemeName(schemeData);

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(schemeName);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "gradScheme<Type>::New(Istream& schemeData)",
            schemeData
        )   << "unknown grad scheme " << schemeName << endl << endl
            << "Valid grad schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
gradScheme<Type>::~gradScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<blockVectorMatrix> gradScheme<Type>::fvmGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    FatalIOErrorIn
    (
        "tmp<blockVectorMatrix> fvmGrad\n",
        "(\n"
        "    GeometricField<Type, fvPatchField, volMesh>&"
        ")\n"
    )   << "Implicit gradient operator currently defined only for Gauss grad."
        << abort(FatalIOError);

    tmp<blockVectorMatrix> tbm
    (
        new blockVectorMatrix
        (
           vf.mesh()
        )
    );

    return tbm;
}


template<class Type>
void gradScheme<Type>::correctBoundaryConditions
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >& gGrad
)
{
    forAll (vsf.boundaryField(), patchi)
    {
        if (!vsf.boundaryField()[patchi].coupled())
        {
            vectorField n = vsf.mesh().boundary()[patchi].nf();

            gGrad.boundaryField()[patchi] += n*
            (
                vsf.boundaryField()[patchi].snGrad()
              - (n & gGrad.boundaryField()[patchi])
            );
        }
    }

    // Note: coupled boundaries provide patchNeighbourField, which is only
    // updated on correct boundary conditions.  Therefore, evaluateCoupled()
    // should be called here. HJ, Apr/2013
    gGrad.boundaryField().evaluateCoupled();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
