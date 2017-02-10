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

#include "steadyInertialDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<volScalarField> steadyInertialDdtScheme<Type>::CorDeltaT() const
{
    // Collect face delta t and pick the smallest for the cell
    surfaceScalarField cofrDeltaT = CofrDeltaT();

    tmp<volScalarField> tcorDeltaT
    (
        new volScalarField
        (
            IOobject
            (
                "CorDeltaT",
                cofrDeltaT.instance(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("CorDeltaT", cofrDeltaT.dimensions(), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& corDeltaT = tcorDeltaT();

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    forAll(owner, faceI)
    {
        corDeltaT[owner[faceI]] =
            max(corDeltaT[owner[faceI]], cofrDeltaT[faceI]);

        corDeltaT[neighbour[faceI]] =
            max(corDeltaT[neighbour[faceI]], cofrDeltaT[faceI]);
    }

    forAll(corDeltaT.boundaryField(), patchi)
    {
        const fvsPatchScalarField& pcofrDeltaT =
            cofrDeltaT.boundaryField()[patchi];

        const fvPatch& p = pcofrDeltaT.patch();
        const unallocLabelList& faceCells = p.patch().faceCells();

        forAll(pcofrDeltaT, patchFacei)
        {
            corDeltaT[faceCells[patchFacei]] = max
            (
                corDeltaT[faceCells[patchFacei]],
                pcofrDeltaT[patchFacei]
            );
        }
    }

    corDeltaT.correctBoundaryConditions();

    return tcorDeltaT;
}


template<class Type>
tmp<surfaceScalarField> steadyInertialDdtScheme<Type>::CofrDeltaT() const
{
    const objectRegistry& registry = this->mesh();

    const surfaceScalarField& phi =
        registry.lookupObject<surfaceScalarField>(phiName_);

    if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
    {
        // Calculate face velocity
        surfaceScalarField magFaceU = mag(phi)/mesh().magSf();

        // Calculate face Co number from local face velocity stabilised to
        // avoid U = 0 case.  Min velocity is assumed to be 1/1000 of the max
        // velocity.  HJ, 8/Mar/2010
        dimensionedScalar faceUlimit =
            Foam::max
            (
                0.001*max(magFaceU),  // Field max will do global reduction
                dimensionedScalar("one", dimVelocity, 1.0)
            );

        return max(magFaceU, faceUlimit)*
            mesh().surfaceInterpolation::deltaCoeffs()/maxCo_;
    }
    else if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
    {
        const volScalarField& rho =
            registry.lookupObject<volScalarField>(rhoName_);

        surfaceScalarField magFaceU  =
            mag(phi)/(fvc::interpolate(rho)*mesh().magSf());

        // Calculate face Co number from local face velocity stabilised to
        // avoid U = 0 case.  Min velocity is assumed to be 1/1000 of the max
        // velocity.  HJ, 8/Mar/2010
        dimensionedScalar faceUlimit = 0.001*max(magFaceU);

        return max(magFaceU, faceUlimit)*
            mesh().surfaceInterpolation::deltaCoeffs()/maxCo_;
    }
    else
    {
        FatalErrorIn("steaddyInertialDdtScheme<Type>::CofrDeltaT() const")
            << "Incorrect dimensions of phi: " << phi.dimensions()
            << abort(FatalError);

        return tmp<surfaceScalarField>(NULL);
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
steadyInertialDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    volScalarField rDeltaT = CorDeltaT();

    IOobject ddtIOobject
    (
        "ddt("+dt.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh> > tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<Type>::zero
                )
            )
        );

        tdtdt().internalField() =
            rDeltaT.internalField()*dt.value()*(1.0 - mesh().V0()/mesh().V());

        return tdtdt;
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<Type>::zero
                ),
                calculatedFvPatchField<Type>::typeName
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
steadyInertialDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    volScalarField rDeltaT = CorDeltaT();

    IOobject ddtIOobject
    (
        "ddt("+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                rDeltaT.internalField()*
                (
                    vf.internalField()
                  - vf.prevIter().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.boundaryField()*
                (
                    vf.boundaryField() - vf.prevIter().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*(vf - vf.prevIter())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
steadyInertialDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    volScalarField rDeltaT = CorDeltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.internalField()*rho.value()*
                (
                    vf.internalField()
                  - vf.prevIter().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.boundaryField()*rho.value()*
                (
                    vf.boundaryField() - vf.prevIter().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*rho*(vf - vf.prevIter())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
steadyInertialDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    volScalarField rDeltaT = CorDeltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.internalField()*
                (
                    rho.internalField()*vf.internalField()
                  - rho.internalField()
                   *vf.prevIter().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.boundaryField()*
                (
                    rho.boundaryField()*vf.boundaryField()
                  - rho.boundaryField()
                   *vf.prevIter().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*rho*(vf - vf.prevIter())
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type> >
steadyInertialDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm();

    scalarField rDeltaT = CorDeltaT()().internalField();

    fvm.diag() = rDeltaT*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*vf.prevIter().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT*vf.prevIter().internalField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
steadyInertialDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    scalarField rDeltaT = CorDeltaT()().internalField();

    fvm.diag() = rDeltaT*rho.value()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.prevIter().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.prevIter().internalField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
steadyInertialDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    scalarField rDeltaT = CorDeltaT()().internalField();

    fvm.diag() = rDeltaT*rho.internalField()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.internalField()
            *vf.prevIter().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.internalField()
            *vf.prevIter().internalField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<typename steadyInertialDdtScheme<Type>::fluxFieldType>
steadyInertialDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    IOobject ddtIOobject
    (
        "ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                mesh(),
                dimensioned<typename flux<Type>::type>
                (
                    "0",
                    rA.dimensions()*phi.dimensions()/dimTime,
                    pTraits<typename flux<Type>::type>::zero
                )
            )
        );
    }
    else
    {
        volScalarField rDeltaT = CorDeltaT();

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                this->fvcDdtPhiCoeff(U, phi)*
                (
                    fvc::interpolate(rDeltaT*rA)*phi
                  - (fvc::interpolate(rDeltaT*rA*U.prevIter()) & mesh().Sf())
                )
            )
        );
    }
}


template<class Type>
tmp<typename steadyInertialDdtScheme<Type>::fluxFieldType>
steadyInertialDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    IOobject ddtIOobject
    (
        "ddtPhiCorr("
      + rA.name() + ',' + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                mesh(),
                dimensioned<typename flux<Type>::type>
                (
                    "0",
                    rA.dimensions()*rho.dimensions()*phi.dimensions()/dimTime,
                    pTraits<typename flux<Type>::type>::zero
                )
            )
        );
    }
    else
    {
        volScalarField rDeltaT = CorDeltaT();

        if
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == dimVelocity*dimArea
        )
        {
            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    ddtIOobject,
                    this->fvcDdtPhiCoeff(U, phi)
                   *(
                        fvc::interpolate(rDeltaT*rA*rho)*phi
                      - (
                            fvc::interpolate(rDeltaT*rA*rho*U.prevIter()
                        )
                      & mesh().Sf())
                    )
                )
            );
        }
        else if
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == dimDensity*dimVelocity*dimArea
        )
        {
            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    ddtIOobject,
                    this->fvcDdtPhiCoeff
                    (
                        U.prevIter(),
                        phi/fvc::interpolate(rho)
                    )
                   *(
                        fvc::interpolate(rDeltaT*rA*rho)
                       *phi/fvc::interpolate(rho)
                      - (
                            fvc::interpolate
                            (
                                rDeltaT*rA*rho*U.prevIter()
                            ) & mesh().Sf()
                        )
                    )
                )
            );
        }
        else if
        (
            U.dimensions() == dimDensity*dimVelocity
         && phi.dimensions() == dimDensity*dimVelocity*dimArea
        )
        {
            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    ddtIOobject,
                    this->fvcDdtPhiCoeff(rho, U, phi)*
                    (
                        fvc::interpolate(rDeltaT*rA)*phi
                      - (
                            fvc::interpolate(rDeltaT*rA*U.prevIter()) & mesh().Sf()
                        )
                    )
                )
            );
        }
        else
        {
            FatalErrorIn
            (
                "steadyInertialDdtScheme<Type>::fvcDdtPhiCorr"
            )   << "dimensions of phi are not correct"
                << abort(FatalError);

            return fluxFieldType::null();
        }
    }
}


template<class Type>
tmp<surfaceScalarField> steadyInertialDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>&
)
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "meshPhi",
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("0", dimVolume/dimTime, 0.0)
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
