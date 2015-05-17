// Do "normal" interpolation
//volPointInterpolation::interpolate(vf, pf);

//- interpolation is done just before this file
pointVectorField& pf = pointDU;


// Do the correction
//GeometricField<Type, pointPatchField, pointMesh>  pfCorr
/*
pointVectorField pfCorr
(
    IOobject
    (
        //  "edgeCorrectedVolPointInterpolate(" + vf.name() + ")Corr",
        "edgeCorrectedVolPointInterpolate(" + DU.name() + ")Corr",
        //vf.instance(),
        DU,
        pMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    pMesh,
    //dimensioned<Type>("zero", pf.dimensions(), pTraits<Type>::zero),
    dimensionedVector("zero", pf.dimensions(), vector::zero),
    pf.boundaryField().types()
);
*/

pointVectorField pfCorr
(
    IOobject
    (
        "pointDUcorr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    pMesh,
    dimensionedVector("vector", dimLength, vector::zero),
    "calculated"
);

//const labelList& ptc = boundaryPoints();
#include "findBoundaryPoints.H"

//const FieldField<Field, vector>& extraVecs = extrapolationVectors(); ///*********
#include "calculateExtrapolationVectors.H"

//const FieldField<Field, scalar>& w = pointBoundaryWeights(); ///*********
#include "calculatePointBoundaryWeights.H"

//Info << "w: " << w << endl;

const labelListList& PointFaces = mesh.pointFaces();

//const fvBoundaryMesh& bm = mesh.boundary(); // declared in findBoundaryPoints.H

forAll (ptc, pointI)
{
    const label curPoint = ptc[pointI];

    const labelList& curFaces = PointFaces[curPoint];

    label fI = 0;

    // Go through all the faces
    forAll (curFaces, faceI)
    {
        if (!mesh.isInternalFace(curFaces[faceI]))
        {
            // This is a boundary face.  If not in the empty patch
            // or coupled calculate the extrapolation vector
            label patchID =
                mesh.boundaryMesh().whichPatch(curFaces[faceI]);

            if
            (
                !isA<emptyFvPatch>(mesh.boundary()[patchID])
                && !mesh.boundary()[patchID].coupled()
            )
            {
                label faceInPatchID =
                    bm[patchID].patch().whichFace(curFaces[faceI]);

                pfCorr[curPoint] +=
                    w[pointI][fI]*
                    (
                        extraVecs[pointI][fI]
                      & gradDU.boundaryField()[patchID][faceInPatchID]
                    );

                fI++;
            }
        }
    }
}

// Update coupled boundaries
/*
forAll (pfCorr.boundaryField(), patchI)
{
    if (pfCorr.boundaryField()[patchI].coupled())
    {
        pfCorr.boundaryField()[patchI].initAddField();
    }
}
*/

/*
forAll (pfCorr.boundaryField(), patchI)
{
    if (pfCorr.boundaryField()[patchI].coupled())
    {
        pfCorr.boundaryField()[patchI].addField(pfCorr.internalField());
    }
}
*/


//Info << "pfCorr: " << pfCorr << endl;
pfCorr.correctBoundaryConditions();

//pfCorr.write();

//Info << "pf: " << pf << endl;
pf += pfCorr;

//- this was already commented
// Replace extrapolated values in pf
//     forAll (ptc, pointI)
//     {
//         pf[ptc[pointI]] = pfCorr[ptc[pointI]];
//     }

pf.correctBoundaryConditions();

//pf.write();
