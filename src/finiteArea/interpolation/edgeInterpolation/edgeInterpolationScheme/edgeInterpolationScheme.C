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
    Abstract base class for edge interpolation schemes.

\*---------------------------------------------------------------------------*/

#include "edgeInterpolationScheme.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "faPatchFields.H"
#include "coupledFaPatchField.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<edgeInterpolationScheme<Type> > edgeInterpolationScheme<Type>::New
(
    const faMesh& mesh,
    Istream& schemeData
)
{
    if (edgeInterpolation::debug)
    {
        Info<< "edgeInterpolationScheme<Type>::New(const faMesh&, Istream&)"
               " : constructing edgeInterpolationScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorIn
        (
            "edgeInterpolationScheme<Type>::New(const faMesh&, Istream&)",
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
            "edgeInterpolationScheme<Type>::New(const faMesh&, Istream&)",
            schemeData
        )   << "Unknown discretisation scheme " << schemeName
            << endl << endl
            << "Valid schemes are :" << endl
            << MeshConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return constructorIter()(mesh, schemeData);
}


// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<edgeInterpolationScheme<Type> > edgeInterpolationScheme<Type>::New
(
    const faMesh& mesh,
    const edgeScalarField& faceFlux,
    Istream& schemeData
)
{
    if (edgeInterpolation::debug)
    {
        Info<< "edgeInterpolationScheme<Type>::New"
               "(const faMesh&, const edgeScalarField&, Istream&) : "
               "constructing edgeInterpolationScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorIn
        (
            "edgeInterpolationScheme<Type>::New"
            "(const faMesh&, const edgeScalarField&, Istream&)",
            schemeData
        )   << "Discretisation scheme not specified"
            << endl << endl
            << "Valid schemes are :" << endl
            << MeshConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    word schemeName(schemeData);

    typename MeshFluxConstructorTable::iterator constructorIter =
        MeshFluxConstructorTablePtr_->find(schemeName);

    if (constructorIter == MeshFluxConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "edgeInterpolationScheme<Type>::New"
            "(const faMesh&, const edgeScalarField&, Istream&)",
            schemeData
        )   << "Unknown discretisation scheme " << schemeName
            << endl << endl
            << "Valid schemes are :" << endl
            << MeshFluxConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return constructorIter()(mesh, faceFlux, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
edgeInterpolationScheme<Type>::~edgeInterpolationScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return the face-interpolate of the given cell field
//  with the given owner and neighbour weighting factors
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
edgeInterpolationScheme<Type>::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const tmp<edgeScalarField>& tlambdas,
    const tmp<edgeScalarField>& tys
)
{
    if (edgeInterpolation::debug)
    {
        Info<< "edgeInterpolationScheme<Type>::uncorrectedInterpolate"
               "(const GeometricField<Type, faPatchField, areaMesh>&, "
               "const tmp<edgeScalarField>&, "
               "const tmp<edgeScalarField>&) : "
               "interpolating areaTypeField from cells to faces "
               "without explicit correction"
            << endl;
    }

    const edgeScalarField& lambdas = tlambdas();
    const edgeScalarField& ys = tys();

    const Field<Type>& vfi = vf.internalField();
    const scalarField& lambda = lambdas.internalField();
    const scalarField& y = ys.internalField();

    const faMesh& mesh = vf.mesh();
    const unallocLabelList& P = mesh.owner();
    const unallocLabelList& N = mesh.neighbour();

    tmp<GeometricField<Type, faePatchField, edgeMesh> > tsf
    (
        new GeometricField<Type, faePatchField, edgeMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    GeometricField<Type, faePatchField, edgeMesh>& sf = tsf();

    Field<Type>& sfi = sf.internalField();

    for (register label fi=0; fi<P.size(); fi++)
    {
        // ZT, 22/Apr/2003
        const tensorField& curT = mesh.edgeTransformTensors()[fi];

        const tensor& Te = curT[0];
        const tensor& TP = curT[1];
        const tensor& TN = curT[2];

        sfi[fi] = 
            transform
            (
                Te.T(),
                lambda[fi]*transform(TP, vfi[P[fi]])
              + y[fi]*transform(TN, vfi[N[fi]])
            );
    }


    // Interpolate across coupled patches using given lambdas and ys

    forAll (lambdas.boundaryField(), pi)
    {
        const faePatchScalarField& pLambda = lambdas.boundaryField()[pi];
        const faePatchScalarField& pY = ys.boundaryField()[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            label size = vf.boundaryField()[pi].patch().size();
            label start = vf.boundaryField()[pi].patch().start();

            Field<Type> pOwnVf = vf.boundaryField()[pi].patchInternalField();
            Field<Type> pNgbVf = vf.boundaryField()[pi].patchNeighbourField();

            Field<Type>& pSf = sf.boundaryField()[pi];

            for (label i=0; i<size; i++)
            {
                const tensorField& curT = 
                    mesh.edgeTransformTensors()[start + i];

                const tensor& Te = curT[0];
                const tensor& TP = curT[1];
                const tensor& TN = curT[2];

                pSf[i] =
                    transform
                    (
                        Te.T(),
                        pLambda[i]*transform(TP, pOwnVf[i])
                      + pY[i]*transform(TN, pNgbVf[i])
                    );
            }

//             sf.boundaryField()[pi] =
//                 pLambda*vf.boundaryField()[pi].patchInternalField()
//               + pY*vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sf.boundaryField()[pi] = vf.boundaryField()[pi];
        }
    }

    tlambdas.clear();
    tys.clear();

    return tsf;
}


//- Return the face-interpolate of the given cell field
//  with the given weigting factors
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
edgeInterpolationScheme<Type>::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const tmp<edgeScalarField>& tlambdas
)
{
    if (edgeInterpolation::debug)
    {
        Info<< "edgeInterpolationScheme<Type>::interpolate"
               "(const GeometricField<Type, faPatchField, areaMesh>&, "
               "const tmp<edgeScalarField>&) : "
               "interpolating areaTypeField from cells to faces "
               "without explicit correction"
            << endl;
    }

    const edgeScalarField& lambdas = tlambdas();

    const Field<Type>& vfi = vf.internalField();
    const scalarField& lambda = lambdas.internalField();

    const faMesh& mesh = vf.mesh();
    const unallocLabelList& P = mesh.owner();
    const unallocLabelList& N = mesh.neighbour();

    tmp<GeometricField<Type, faePatchField, edgeMesh> > tsf
    (
        new GeometricField<Type, faePatchField, edgeMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    GeometricField<Type, faePatchField, edgeMesh>& sf = tsf();

    Field<Type>& sfi = sf.internalField();

    for (register label eI = 0; eI < P.size(); eI++)
    {
        // ZT, 22/Apr/2003
        const tensorField& curT = mesh.edgeTransformTensors()[eI];

        const tensor& Te = curT[0];
        const tensor& TP = curT[1];
        const tensor& TN = curT[2];

        sfi[eI] = 
            transform
            (
                Te.T(),
                lambda[eI]*transform(TP, vfi[P[eI]]) 
              + (1 - lambda[eI])*transform(TN, vfi[N[eI]])
            );
    }


    // Interpolate across coupled patches using given lambdas

    forAll (lambdas.boundaryField(), pi)
    {
        const faePatchScalarField& pLambda = lambdas.boundaryField()[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            label size = vf.boundaryField()[pi].patch().size();
            label start = vf.boundaryField()[pi].patch().start();

            Field<Type> pOwnVf = vf.boundaryField()[pi].patchInternalField();
            Field<Type> pNgbVf = vf.boundaryField()[pi].patchNeighbourField();

            Field<Type>& pSf = sf.boundaryField()[pi];

            for (label i=0; i<size; i++)
            {
                const tensorField& curT =
                    mesh.edgeTransformTensors()[start + i];

                const tensor& Te = curT[0];
                const tensor& TP = curT[1];
                const tensor& TN = curT[2];

                pSf[i] =
                    transform
                    (
                        Te.T(),
                        pLambda[i]*transform(TP, pOwnVf[i])
                      + (1 - pLambda[i])*transform(TN, pNgbVf[i])
                    );
            }

//             tsf().boundaryField()[pi] =
//                 pLambda*vf.boundaryField()[pi].patchInternalField()
//              + (1 - pLambda)*vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sf.boundaryField()[pi] = vf.boundaryField()[pi];
        }
    }

    tlambdas.clear();

    return tsf;
}


//- Return the euclidian edge-interpolate of the given area field
//  with the given weigting factors
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
edgeInterpolationScheme<Type>::euclidianInterpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const tmp<edgeScalarField>& tlambdas
)
{
    if (edgeInterpolation::debug)
    {
        Info<< "edgeInterpolationScheme<Type>::euclidianInterpolate"
               "(const GeometricField<Type, faPatchField, areaMesh>&, "
               "const tmp<edgeScalarField>&) : "
               "interpolating areaTypeField from cells to faces "
               "without explicit correction"
            << endl;
    }

    const edgeScalarField& lambdas = tlambdas();

    const Field<Type>& vfi = vf.internalField();
    const scalarField& lambda = lambdas.internalField();

    const faMesh& mesh = vf.mesh();
    const unallocLabelList& P = mesh.owner();
    const unallocLabelList& N = mesh.neighbour();

    tmp<GeometricField<Type, faePatchField, edgeMesh> > tsf
    (
        new GeometricField<Type, faePatchField, edgeMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    GeometricField<Type, faePatchField, edgeMesh>& sf = tsf();

    Field<Type>& sfi = sf.internalField();

    for (register label eI = 0; eI < P.size(); eI++)
    {
        sfi[eI] = lambda[eI]*vfi[P[eI]] + (1 - lambda[eI])*vfi[N[eI]];
    }


    // Interpolate across coupled patches using given lambdas

    forAll (lambdas.boundaryField(), pi)
    {
        const faePatchScalarField& pLambda = lambdas.boundaryField()[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            tsf().boundaryField()[pi] =
                pLambda*vf.boundaryField()[pi].patchInternalField()
             + (1.0 - pLambda)*vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sf.boundaryField()[pi] = vf.boundaryField()[pi];
        }
    }

    tlambdas.clear();

    return tsf;
}


//- Return the face-interpolate of the given cell field
//  with explicit correction
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
edgeInterpolationScheme<Type>::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    if (edgeInterpolation::debug)
    {
        Info<< "edgeInterpolationScheme<Type>::interpolate"
               "(const GeometricField<Type, faPatchField, areaMesh>&) : "
            << "interpolating areaTypeField from cells to faces"
            << endl;
    }

    tmp<GeometricField<Type, faePatchField, edgeMesh> > tsf
        = interpolate(vf, weights(vf));

    if (corrected())
    {
        tsf() += correction(vf);
    }

    return tsf;
}

//- Return the euclidian edge-interpolate of the given area field
//  without explicit correction
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
edgeInterpolationScheme<Type>::euclidianInterpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    if (edgeInterpolation::debug)
    {
        Info<< "edgeInterpolationScheme<Type>::interpolate"
               "(const GeometricField<Type, faPatchField, areaMesh>&) : "
            << "interpolating areaTypeField from cells to faces"
            << endl;
    }

    tmp<GeometricField<Type, faePatchField, edgeMesh> > tsf
        = euclidianInterpolate(vf, weights(vf));

    return tsf;
}

//- Return the face-interpolate of the given cell field
//  with explicit correction
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
edgeInterpolationScheme<Type>::interpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
) const
{
    tmp<GeometricField<Type, faePatchField, edgeMesh> > tinterpVf
        = interpolate(tvf());
    tvf.clear();
    return tinterpVf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
