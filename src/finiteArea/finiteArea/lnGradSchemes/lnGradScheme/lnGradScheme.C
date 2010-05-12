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
    Abstract base class for lnGrad schemes.

\*---------------------------------------------------------------------------*/

#include "fa.H"
#include "lnGradScheme.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<lnGradScheme<Type> > lnGradScheme<Type>::New
(
    const faMesh& mesh,
    Istream& schemeData
)
{
    if (fa::debug)
    {
        Info<< "lnGradScheme<Type>::New(const faMesh&, Istream&)"
               " : constructing lnGradScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorIn
        (
            "lnGradScheme<Type>::New(const faMesh&, Istream&)",
            schemeData
        )   << "Discretisation scheme not specified"
            << endl << endl
            << "Valid schemes are :" << endl
            << MeshConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    word schemeName(schemeData);

    typename MeshConstructorTable::iterator constructorIter =
        MeshConstructorTablePtr_->find(schemeName);

    if (constructorIter == MeshConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "lnGradScheme<Type>::New(const faMesh&, Istream&)",
            schemeData
        )   << "Unknown discretisation scheme " << schemeName
            << endl << endl
            << "Valid schemes are :" << endl
            << MeshConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return constructorIter()(mesh, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
lnGradScheme<Type>::~lnGradScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return the face-lnGrad of the given cell field
//  with the given weigting factors
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
lnGradScheme<Type>::lnGrad
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const tmp<edgeScalarField>& tdeltaCoeffs,
    const word& lnGradName
)
{
    const faMesh& mesh = vf.mesh();

    // construct GeometricField<Type, faePatchField, edgeMesh>
    tmp<GeometricField<Type, faePatchField, edgeMesh> > tssf
    (
        new GeometricField<Type, faePatchField, edgeMesh>
        (
            IOobject
            (
                lnGradName + vf.name() + ')',
                vf.instance(),
                vf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            vf.dimensions()*tdeltaCoeffs().dimensions()
        )
    );
    GeometricField<Type, faePatchField, edgeMesh>& ssf = tssf();

    // set reference to difference factors array
    const scalarField& deltaCoeffs = tdeltaCoeffs().internalField();

    // owner/neighbour addressing
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    forAll(owner, faceI)
    {
        ssf[faceI] =
            deltaCoeffs[faceI]*(vf[neighbour[faceI]] - vf[owner[faceI]]);
    }

    forAll(vf.boundaryField(), patchI)
    {
        ssf.boundaryField()[patchI] = vf.boundaryField()[patchI].snGrad();
    }

    return tssf;
}


//- Return the face-lnGrad of the given cell field
//  with explicit correction
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
lnGradScheme<Type>::lnGrad
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    tmp<GeometricField<Type, faePatchField, edgeMesh> > tsf
        = lnGrad(vf, deltaCoeffs(vf));

    if (corrected())
    {
        tsf() += correction(vf);
    }

    return tsf;
}


//- Return the face-lnGrad of the given cell field
//  with explicit correction
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
lnGradScheme<Type>::lnGrad
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
) const
{
    tmp<GeometricField<Type, faePatchField, edgeMesh> > tinterpVf
        = lnGrad(tvf());
    tvf.clear();
    return tinterpVf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
