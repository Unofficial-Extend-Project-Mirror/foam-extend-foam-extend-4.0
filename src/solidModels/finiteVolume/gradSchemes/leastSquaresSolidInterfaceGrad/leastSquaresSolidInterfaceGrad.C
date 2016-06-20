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

#include "leastSquaresSolidInterfaceGrad.H"
#include "leastSquaresSolidInterfaceVectors.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "GeometricField.H"
#include "zeroGradientFvPatchField.H"
#include "solidInterface.H"
#include "IOReferencer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
leastSquaresSolidInterfaceGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vsf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tlsGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                name,
                vsf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "zero",
                vsf.dimensions()/dimLength,
                pTraits<GradType>::zero
            ),
            zeroGradientFvPatchField<GradType>::typeName
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& lsGrad = tlsGrad();

    // Get reference to least square vectors
    const leastSquaresSolidInterfaceVectors& lsv =
        leastSquaresSolidInterfaceVectors::New(mesh);

    const surfaceVectorField& ownLs = lsv.pVectors();
    const surfaceVectorField& neiLs = lsv.nVectors();

    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nei = mesh.neighbour();

    // interface fields
    const solidInterface& solInt =
      mesh.objectRegistry::lookupObject<IOReferencer<solidInterface> >
        ("solidInterface")();
    const labelList& interfaceFacesMap = solInt.indicatorFieldMap();
    const labelList& interfaceProcPatchMap = solInt.processorPatchMap();
    const labelListList& interfaceProcPatchFacesMap =
        solInt.processorPatchFacesMap();
    const vectorField& interfaceDispVec = solInt.interfaceDisplacement();
    const List<vectorField>& procInterDispVec =
        solInt.processorInterfaceDisplacement();

    // interface disp must be vector so
    // this is a quick hack
    Field<Type> interfaceDisp(interfaceDispVec.size(), pTraits<Type>::zero);
    if (pTraits<Type>::typeName != vector::typeName)
      {
          FatalError
              << "leastSquaresSolidInterface::grad() may only be used "
              "with a volVectorField" << exit(FatalError);
      }
    for (direction cmpt = 0; cmpt < pTraits<vector>::nComponents; cmpt++)
      {
          interfaceDisp.replace
              (
                  cmpt,
                  interfaceDispVec.component(cmpt)
                  );
      }
    List<Field<Type> > procInterDisp
        (procInterDispVec.size(), Field<Type>(0, pTraits<Type>::zero));
    forAll(procInterDisp, patchi)
      {
          procInterDisp[patchi].setSize
              (procInterDispVec[patchi].size(), pTraits<Type>::zero);

          for (direction cmpt = 0; cmpt < pTraits<vector>::nComponents; cmpt++)
          {
              procInterDisp[patchi].replace
                  (
                      cmpt,
                      procInterDispVec[patchi].component(cmpt)
                      );
          }
      }

    forAll(own, facei)
    {
        register label ownFaceI = own[facei];
        register label neiFaceI = nei[facei];

        if (interfaceFacesMap[facei] > -SMALL)
        {
            label interfaceFacei = interfaceFacesMap[facei];
            if (interfaceFacei == -1)
            {
                FatalError
                    << "leastSquaresSolidInterface::grad()"
                    << "face " << facei << " is not on the interface"
                    << exit(FatalError);
            }

            // for interface faces, we use the face centre value
            // instead of the neighbour cell centre
            Type ownDeltaVsf = interfaceDisp[interfaceFacei] - vsf[ownFaceI];
            Type neiDeltaVsf = vsf[neiFaceI] - interfaceDisp[interfaceFacei];

            lsGrad[ownFaceI] += ownLs[facei]*ownDeltaVsf;
            lsGrad[neiFaceI] -= neiLs[facei]*neiDeltaVsf;
        }
        else
        {
            // standard method
            Type deltaVsf = vsf[neiFaceI] - vsf[ownFaceI];

            lsGrad[ownFaceI] += ownLs[facei]*deltaVsf;
            lsGrad[neiFaceI] -= neiLs[facei]*deltaVsf;
        }
    }

    // Boundary faces
    forAll(vsf.boundaryField(), patchi)
    {
        const fvsPatchVectorField& patchOwnLs = ownLs.boundaryField()[patchi];

        const unallocLabelList& faceCells =
            lsGrad.boundaryField()[patchi].patch().faceCells();

        if (vsf.boundaryField()[patchi].coupled())
        {
            Field<Type> neiVsf =
                vsf.boundaryField()[patchi].patchNeighbourField();

            forAll(neiVsf, patchFaceI)
            {
                // philipc: special treatment of solid inter faces
                if (interfaceProcPatchFacesMap[patchi][patchFaceI] > -SMALL)
                {
                    label curPatch = interfaceProcPatchMap[patchi];
                    label curFace =
                        interfaceProcPatchFacesMap[patchi][patchFaceI];
                    if (curPatch == -1 || curFace == -1)
                    {
                        FatalError
                            << "leastSquaresSolidInterface::grad()"
                            << "proc face " << patchFaceI
                            << " is not on the interface"
                            << exit(FatalError);
                    }
                    // face is on solid interface
                    lsGrad[faceCells[patchFaceI]] +=
                        patchOwnLs[patchFaceI]
                        *(procInterDisp[curPatch][curFace]
                          - vsf[faceCells[patchFaceI]]);
                }
                else
                {
                    // standard method
                    lsGrad[faceCells[patchFaceI]] +=
                        patchOwnLs[patchFaceI]
                        *(neiVsf[patchFaceI] - vsf[faceCells[patchFaceI]]);
                }
            }
        }
        else
        {
            const fvPatchField<Type>& patchVsf = vsf.boundaryField()[patchi];

            forAll(patchVsf, patchFaceI)
            {
                lsGrad[faceCells[patchFaceI]] +=
                     patchOwnLs[patchFaceI]
                    *(patchVsf[patchFaceI] - vsf[faceCells[patchFaceI]]);
            }
        }
    }


    lsGrad.correctBoundaryConditions();
    gaussGrad<Type>::correctBoundaryConditions(vsf, lsGrad);

    return tlsGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
