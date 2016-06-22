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

Class
    pointVolInterpolation

Description

\*---------------------------------------------------------------------------*/

#include "pointVolInterpolation.H"
#include "volFields.H"
#include "pointFields.H"
#include "primitiveMesh.H"
#include "emptyFvPatch.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::pointVolInterpolation::interpolate
(
    const GeometricField<Type, pointPatchField, pointMesh>& pf,
    GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (debug)
    {
        Info<< "pointVolInterpolation::interpolate("
            << "const GeometricField<Type, pointPatchField, pointMesh>&, "
            << "GeometricField<Type, fvPatchField, volMesh>&) : "
            << "interpolating field from points to cells"
            << endl;
    }

    const FieldField<Field, scalar>& weights = volWeights();
    const labelListList& cellPoints = vf.mesh().cellPoints();

    // Multiply pointField by weighting factor matrix to create volField
    forAll(cellPoints, cellI)
    {
        vf[cellI] = pTraits<Type>::zero;

        const labelList& curCellPoints = cellPoints[cellI];

        forAll(curCellPoints, cellPointI)
        {
            vf[cellI] +=
                weights[cellI][cellPointI]*pf[curCellPoints[cellPointI]];
        }
    }


    // Interpolate patch values: over-ride the internal values for the points
    // on the patch with the interpolated point values from the faces
    const fvBoundaryMesh& bm = vMesh().boundary();

    const PtrList<primitivePatchInterpolation>& pi = patchInterpolators();
    forAll (bm, patchI)
    {
        // If the patch is empty, skip it
        if
        (
            bm[patchI].type() != emptyFvPatch::typeName
        )
        {
            vf.boundaryField()[patchI] =
                pi[patchI].pointToFaceInterpolate
                (
                    pf.boundaryField()[patchI].patchInternalField()
                );
        }
    }

    vf.correctBoundaryConditions();

    if (debug)
    {
        Info<< "pointVolInterpolation::interpolate("
            << "const GeometricField<Type, pointPatchField, pointMesh>&, "
            << "GeometricField<Type, fvPatchField, volMesh>&) : "
            << "finished interpolating field from points to cells"
            << endl;
    }
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::pointVolInterpolation::interpolate
(
    const GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    // Construct tmp<pointField>
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "pointVolInterpolate(" + pf.name() + ')',
                pf.instance(),
                pf.db()
            ),
            vMesh(),
            pf.dimensions()
        )
    );

    // Perform interpolation
    interpolate(pf, tvf());

    return tvf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::pointVolInterpolation::interpolate
(
    const tmp<GeometricField<Type, pointPatchField, pointMesh> >& tpf
) const
{
    // Construct tmp<volField>
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf =
        interpolate(tpf());
    tpf.clear();
    return tvf;
}


// ************************************************************************* //
