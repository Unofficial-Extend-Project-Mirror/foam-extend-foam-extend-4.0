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

\*---------------------------------------------------------------------------*/

#include "solidInterfaceTL.H"
#include "volFields.H"
#include "pointFields.H"
#include "globalPointPatch.H"
#include "emptyPolyPatch.H"
#include "volSurfaceMapping.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
void solidInterfaceTL::volToPointInterpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Info<< "solidInterfaceTL::volToPointInterpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    Field<Type>& pfI = pf.internalField();

    forAll(subMeshes(), meshI)
    {
        GeometricField<Type, fvPatchField, volMesh> subVf =
            subMeshes()[meshI].interpolate(vf);

        GeometricField<Type, pointPatchField, pointMesh> subPf =
            subMeshes()[meshI].interpolate(pf);

        subMeshVolToPoint()[meshI].interpolate(subVf, subPf);

        const Field<Type>& subPfI = subPf.internalField();

        // Map point field from sub-mesh to global mesh
        const labelList& pointMap = subMeshes()[meshI].pointMap();
        forAll(pointMap, pointI)
        {
            pfI[pointMap[pointI]] = subPfI[pointI];
        }
    }

    pf.correctBoundaryConditions();
}


template<class Type>
void solidInterfaceTL::volToPointInterpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const Field<Type>& interfaceVf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Info<< "solidInterfaceTL::volToPointInterpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    Field<Type>& pfI = pf.internalField();
    pfI = pTraits<Type>::zero;

    const labelList& noMat = pointNumOfMaterials();

    const fvMesh& mesh_ = D_.mesh();

    const labelList& spLabels =
        mesh_.globalData().sharedPointLabels();

    const labelList& spAddressing =
        mesh_.globalData().sharedPointAddr();

    List<List<Map<Type> > > glData(Pstream::nProcs());
    forAll(glData, procI)
    {
        glData[procI] =
            List<Map<Type> >(mesh_.globalData().nGlobalPoints(), Map<Type>());
    }

    forAll(subMeshes(), meshI)
    {
        // Sub mesh

        GeometricField<Type, fvPatchField, volMesh> subVf =
            subMeshes()[meshI].interpolate(vf);

        // Correct volume field at the interface
        const labelList& patchMap = subMeshes()[meshI].patchMap();
        label interfacePatchIndex = -1;
        forAll(patchMap, patchI)
        {
            if (patchMap[patchI] == -1)
            {
                interfacePatchIndex = patchI;
                break;
            }
        }
        Field<Type>& interfacePatchVf =
            subVf.boundaryField()[interfacePatchIndex];

        const labelList& fm = subMeshes()[meshI].faceMap();

        label interfacePatchStart =
            subMeshes()[meshI].subMesh().boundaryMesh()
            [
                interfacePatchIndex
            ].start();

        forAll(interfacePatchVf, faceI)
        {
            label curInterFace =
                findIndex(faces(), fm[interfacePatchStart + faceI]);

            interfacePatchVf[faceI] = interfaceVf[curInterFace];
        }

        GeometricField<Type, pointPatchField, pointMesh> subPf =
            subMeshes()[meshI].interpolate(pf);

        subMeshVolToPoint()[meshI].interpolate(subVf, subPf);

        const Field<Type>& subPfI = subPf.internalField();

        // Map point field from sub-mesh to global mesh
        const labelList& pointMap = subMeshes()[meshI].pointMap();
        forAll(pointMap, pointI)
        {
            label curMeshPoint = pointMap[pointI];

            bool sharedPoint(findIndex(spLabels, curMeshPoint) != -1);

            if (sharedPoint)
            {
                label k = findIndex(spLabels, curMeshPoint);
                label curSpIndex = spAddressing[k];
                glData[Pstream::myProcNo()][curSpIndex].insert
                    (
                        meshI,
                        subPfI[pointI]
                    );
            }
            else
            {
                pfI[curMeshPoint] += subPfI[pointI]/noMat[curMeshPoint];
            }
        }
    }

    Pstream::gatherList(glData);
    Pstream::scatterList(glData);

    // Gloabal points
    if (mesh_.globalData().nGlobalPoints())
    {
        for (label k=0; k<mesh_.globalData().nGlobalPoints(); k++)
        {
            label curSpIndex = findIndex(spAddressing, k);

            if (curSpIndex != -1)
            {
                List<label> matN(subMeshes().size(), 0);
                List<Type> matAvg(subMeshes().size(), pTraits<Type>::zero);

                forAll(glData, procI)
                {
                    const Map<Type>& curProcGlData = glData[procI][k];
                    for (label i=0; i<subMeshes().size(); i++)
                    {
                        if (curProcGlData.found(i))
                        {
                            matAvg[i] += curProcGlData[i];
                            matN[i]++;
                        }
                    }
                }

                label nMat = 0;
                Type avg = pTraits<Type>::zero;
                forAll(matAvg, matI)
                {
                    if (matN[matI])
                    {
                        matAvg[matI] /= matN[matI];
                        avg += matAvg[matI];
                        nMat++;
                    }
                }
                avg /= nMat;

                label curMeshPoint = spLabels[curSpIndex];
                pfI[curMeshPoint] = avg;
            }
        }
    }

    pf.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
