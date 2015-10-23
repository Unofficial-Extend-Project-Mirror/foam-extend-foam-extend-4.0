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

Description
    Abstract base class for finite volume calculus gradient schemes.

\*---------------------------------------------------------------------------*/

#include "fv.H"
#include "HashTable.H"
#include "primitiveFields.H"
#include "objectRegistry.H"

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
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::gradScheme<Type>::grad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    if (!this->mesh().changing() && this->mesh().schemesDict().cache(name))
    {
        if (!mesh().objectRegistry::template foundObject<GradFieldType>(name))
        {
            if (fvSchemes::debug)
            {
                Info << "Cache: Calculating and caching " << name
                    << " originating from " << vsf.name()
                    << " event No. " << vsf.eventNo()
                    << endl;
            }
            tmp<GradFieldType> tgGrad = calcGrad(vsf, name);
            regIOobject::store(tgGrad.ptr());
        }

        if (fvSchemes::debug)
        {
            Info << "Cache: Retrieving " << name
                << " originating from " << vsf.name()
                << " event No. " << vsf.eventNo()
                << endl;
        }
        GradFieldType& gGrad = const_cast<GradFieldType&>
        (
            mesh().objectRegistry::template lookupObject<GradFieldType>(name)
        );

        if (gGrad.upToDate(vsf.name()))
        {
            return gGrad;
        }
        else
        {
            if (fvSchemes::debug)
            {
                Info << "Cache: Deleting " << name
                    << " originating from " << vsf.name()
                    << " event No. " << vsf.eventNo()
                    << endl;
            }
            gGrad.release();
            delete &gGrad;

            if (fvSchemes::debug)
            {
                Info << "Cache: Recalculating " << name
                    << " originating from " << vsf.name()
                    << " event No. " << vsf.eventNo()
                    << endl;
            }
            tmp<GradFieldType> tgGrad = calcGrad(vsf, name);

            if (fvSchemes::debug)
            {
                Info << "Cache: Storing " << name
                    << " originating from " << vsf.name()
                    << " event No. " << vsf.eventNo()
                    << endl;
            }
            regIOobject::store(tgGrad.ptr());
            GradFieldType& gGrad = const_cast<GradFieldType&>
            (
                mesh().objectRegistry::template lookupObject<GradFieldType>
                (
                    name
                )
            );

            return gGrad;
        }
    }
    else
    {
        if (mesh().objectRegistry::template foundObject<GradFieldType>(name))
        {
            GradFieldType& gGrad = const_cast<GradFieldType&>
            (
                mesh().objectRegistry::template lookupObject<GradFieldType>
                (
                    name
                )
            );

            if (gGrad.ownedByRegistry())
            {
                if (fvSchemes::debug)
                {
                    Info << "Cache: Deleting " << name
                        << " originating from " << vsf.name()
                        << " event No. " << vsf.eventNo()
                        << endl;
                }
                gGrad.release();
                delete &gGrad;
            }
        }

        if (fvSchemes::debug)
        {
            Info << "Cache: Calculating " << name
                << " originating from " << vsf.name()
                << " event No. " << vsf.eventNo()
                << endl;
        }
        return calcGrad(vsf, name);
    }
}

template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::gradScheme<Type>::grad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf
) const
{
    return grad(vsf, "grad(" + vsf.name() + ')');
}


template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::gradScheme<Type>::grad
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvsf
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    tmp<GradFieldType> tgrad = grad(tvsf());
    tvsf.clear();
    return tgrad;
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename outerProduct<vector, Type>::type>
>
gradScheme<Type>::fvmGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    FatalErrorIn
    (
        "tmp<BlockLduSystem> gradScheme<Type>::fvmGrad\n"
        "(\n"
        "    GeometricField<Type, fvPatchField, volMesh>&"
        ")\n"
    )   << "Implicit gradient operator currently defined only for Gauss linear "
        << "and leastSquares (cell and face limiters are optional)."
        << abort(FatalError);

    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<BlockLduSystem<vector, GradType> > tbs
    (
        new BlockLduSystem<vector, GradType>(vf.mesh())
    );

    return tbs;
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
