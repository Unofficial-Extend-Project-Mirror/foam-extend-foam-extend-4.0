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

\*---------------------------------------------------------------------------*/

#include "edgeCorrectedVolPointInterpolation.H"
#include "emptyFvPatch.H"

// Interpolate from volField to pointField
// using inverse distance weighting with boundary correction
// given the field and its gradient
template<class Type>
void Foam::edgeCorrectedVolPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >& vGradf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Info<< "edgeCorrectedVolPointInterpolation::interpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    // Check if the gradient field is dimensionally correct to be the gradient
    // of vf
    if (vf.dimensions()/dimLength != vGradf.dimensions())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "void Foam::edgeCorrectedVolPointInterpolation::interpolate\n"
            "(\n"
            "    const GeometricField<Type, fvPatchField, volMesh>& vf,\n"
            "    const GeometricField\n"
            "    <\n"
            "        typename outerProduct<vector, Type>::type,\n"
            "        fvPatchField,\n"
            "        volMesh\n"
            "    >& vGradf,\n"
            "    GeometricField<Type, pointPatchField, pointMesh>& pf\n"
            ") const\n"
        )   << "For the interpolation to work correctly, vGradf field "
            << "should represent the gradient of vf.  Dimensionally, "
            << "this is not the case"
            << abort(FatalError);
    }

    // Do "normal" interpolation
    volPointInterpolation::interpolate(vf, pf);

    // Do the correction

    GeometricField<Type, pointPatchField, pointMesh>  pfCorr
    (
        IOobject
        (
            "edgeCorrectedVolPointInterpolate(" + vf.name() + ")Corr",
            vf.instance(),
            vf.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensioned<Type>("zero", pf.dimensions(), pTraits<Type>::zero),
        pf.boundaryField().types()
    );

    const labelList& ptc = boundaryPoints();

    const vectorListList& extraVecs = extrapolationVectors();

    const FieldField<Field, scalar>& w = pointBoundaryWeights();

    const labelListList& PointFaces = vMesh().pointFaces();

    const fvBoundaryMesh& bm = vMesh().boundary();

    forAll (ptc, pointI)
    {
        const label curPoint = ptc[pointI];

        const labelList& curFaces = PointFaces[curPoint];

        label fI = 0;

        // Go through all the faces
        forAll (curFaces, faceI)
        {
            if (!vMesh().isInternalFace(curFaces[faceI]))
            {
                // This is a boundary face.  If not in the empty patch
                // or coupled calculate the extrapolation vector
                label patchID =
                    vMesh().boundaryMesh().whichPatch(curFaces[faceI]);

                if
                (
                    !isA<emptyFvPatch>(vMesh().boundary()[patchID])
                 && !vMesh().boundary()[patchID].coupled()
                )
                {
                    label faceInPatchID =
                        bm[patchID].patch().whichFace(curFaces[faceI]);

                    pfCorr[curPoint] +=
                        w[pointI][fI]*
                        (
                            extraVecs[pointI][fI]
                          & vGradf.boundaryField()[patchID][faceInPatchID]
                        );

                    fI++;
                }
            }
        }
    }

    // Update coupled boundaries
    forAll (pfCorr.boundaryField(), patchI)
    {
        if (pfCorr.boundaryField()[patchI].coupled())
        {
            pfCorr.boundaryField()[patchI].initAddField();
        }
    }

    forAll (pfCorr.boundaryField(), patchI)
    {
        if (pfCorr.boundaryField()[patchI].coupled())
        {
            pfCorr.boundaryField()[patchI].addField(pfCorr.internalField());
        }
    }

    pfCorr.correctBoundaryConditions();
    pf += pfCorr;
    // Replace extrapolated values in pf
//     forAll (ptc, pointI)
//     {
//         pf[ptc[pointI]] = pfCorr[ptc[pointI]];
//     }

    pf.correctBoundaryConditions();

    if (debug)
    {
        Info<< "volPointInterpolation::interpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "finished interpolating field from cells to points"
            << endl;
    }
}


// Interpolate volField returning pointField
// using inverse distance weighting with boundary correction
// given the field and its gradient
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh> >
Foam::edgeCorrectedVolPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >& vGradf
) const
{
    // Construct tmp<pointField>
    tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf
    (
        new GeometricField<Type, pointPatchField, pointMesh>
        (
            IOobject
            (
                "edgeCorrectedVolPointInterpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            pMesh(),
            vf.dimensions()
        )
    );

    // Perform interpolation
    interpolate(vf, vGradf, tpf());

    return tpf;
}


// Interpolate tmp<volField> returning pointField
// using inverse distance weighting with boundary correction
// given the field and its gradient
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh> >
Foam::edgeCorrectedVolPointInterpolation::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >& vGradf
) const
{
    // Construct tmp<pointField>
    tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf =
        interpolate(tvf(), vGradf);
    tvf.clear();

    return tpf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh> >
Foam::edgeCorrectedVolPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const tmp<GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    > >& tvGradf
) const
{
    // Construct tmp<pointField>
    tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf =
        interpolate(vf, tvGradf());
    tvGradf.clear();

    return tpf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh> >
Foam::edgeCorrectedVolPointInterpolation::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const tmp<GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    > >& tvGradf
) const
{
    // Construct tmp<pointField>
    tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf =
        interpolate(tvf(), tvGradf());
    tvf.clear();
    tvGradf.clear();

    return tpf;
}


// ************************************************************************* //
