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

#include "facEdgeIntegrate.H"
#include "faMesh.H"
#include "zeroGradientFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
edgeIntegrate
(
    const GeometricField<Type, faePatchField, edgeMesh>& ssf
)
{
    const faMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, faPatchField, areaMesh> > tvf
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "edgeIntegrate("+ssf.name()+')',
                ssf.instance(),
                ssf.db()
            ),
            mesh,
            dimensioned<Type>
            (
                "0",
                ssf.dimensions()/dimArea,
                pTraits<Type>::zero
            ),
            zeroGradientFaPatchField<Type>::typeName
        )
    );
    GeometricField<Type, faPatchField, areaMesh>& vf = tvf();


    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    forAll(owner, faceI)
    {
        vf[owner[faceI]] += ssf[faceI];
        vf[neighbour[faceI]] -= ssf[faceI];
    }

    forAll(mesh.boundary(), patchI)
    {
        const unallocLabelList& pEdgeFaces =
            mesh.boundary()[patchI].edgeFaces();

        const faePatchField<Type>& pssf = ssf.boundaryField()[patchI];

        forAll(mesh.boundary()[patchI], faceI)
        {
            vf[pEdgeFaces[faceI]] += pssf[faceI];
        }
    }

    vf.internalField() /= mesh.S();
    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
edgeIntegrate
(
    const tmp<GeometricField<Type, faePatchField, edgeMesh> >& tssf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > tvf
    (
        fac::edgeIntegrate(tssf())
    );
    tssf.clear();
    return tvf;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
edgeSum
(
    const GeometricField<Type, faePatchField, edgeMesh>& ssf
)
{
    const faMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, faPatchField, areaMesh> > tvf
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "edgeSum("+ssf.name()+')',
                ssf.instance(),
                ssf.db()
            ),
            mesh,
            dimensioned<Type>("0", ssf.dimensions(), pTraits<Type>::zero),
            zeroGradientFaPatchField<Type>::typeName
        )
    );
    GeometricField<Type, faPatchField, areaMesh>& vf = tvf();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    forAll(owner, faceI)
    {
        vf[owner[faceI]] += ssf[faceI];
        vf[neighbour[faceI]] += ssf[faceI];
    }

    forAll(mesh.boundary(), patchI)
    {
        const unallocLabelList& pEdgeFaces =
            mesh.boundary()[patchI].edgeFaces();

        const faePatchField<Type>& pssf = ssf.boundaryField()[patchI];

        forAll(mesh.boundary()[patchI], faceI)
        {
            vf[pEdgeFaces[faceI]] += pssf[faceI];
        }
    }

    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> > edgeSum
(
    const tmp<GeometricField<Type, faePatchField, edgeMesh> >& tssf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > tvf =
        edgeSum(tssf());
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
