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
    
\*---------------------------------------------------------------------------*/

#include "EulerFaDdtScheme.H"
#include "facDiv.H"
#include "faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
EulerFaDdtScheme<Type>::facDdt
(
    const dimensioned<Type> dt
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+dt.name()+')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (mesh().moving())
    {
        tmp<GeometricField<Type, faPatchField, areaMesh> > tdtdt
        (
            new GeometricField<Type, faPatchField, areaMesh>
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
            rDeltaT.value()*dt.value()*(1.0 - mesh().S0()/mesh().S());

        return tdtdt;
    }
    else
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<Type>::zero
                ),
                calculatedFaPatchField<Type>::typeName
            )
        );
    }
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
EulerFaDdtScheme<Type>::facDdt0
(
    const dimensioned<Type> dt
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+dt.name()+')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    tmp<GeometricField<Type, faPatchField, areaMesh> > tdtdt0
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            ddtIOobject,
            mesh(),
            -rDeltaT*dt
        )
    );

    if (mesh().moving())
    {
        tdtdt0().internalField() =
            (-rDeltaT.value()*dt.value())*mesh().S0()/mesh().S();
    }

    return tdtdt0;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
EulerFaDdtScheme<Type>::facDdt
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+vf.name()+')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    vf.internalField()
                  - vf.oldTime().internalField()*mesh().S0()/mesh().S()
                ),
                rDeltaT.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                rDeltaT*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
EulerFaDdtScheme<Type>::facDdt0
(
    const GeometricField<Type, faePatchField, edgeMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt0("+vf.name()+')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, faePatchField, edgeMesh> >
        (
            new GeometricField<Type, faePatchField, edgeMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                (-rDeltaT.value())
                   *vf.oldTime().internalField()*mesh().S0()/mesh().S(),
                (-rDeltaT.value())
                   *vf.oldTime().boundaryField()
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, faePatchField, edgeMesh> >
        (
            new GeometricField<Type, faePatchField, edgeMesh>
            (
                ddtIOobject,
                (-rDeltaT)*vf.oldTime()
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
EulerFaDdtScheme<Type>::facDdt0
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt0("+vf.name()+')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                (-rDeltaT.value())*vf.oldTime().internalField(),
                (-rDeltaT.value())*vf.oldTime().boundaryField()
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                (-rDeltaT)*vf.oldTime()
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
EulerFaDdtScheme<Type>::facDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*rho.value()*
                (
                    vf.internalField()
                  - vf.oldTime().internalField()*mesh().S0()/mesh().S()
                ),
                rDeltaT.value()*rho.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                rDeltaT*rho*(vf - vf.oldTime())
            )
        );
    }
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
EulerFaDdtScheme<Type>::facDdt0
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt0("+rho.name()+','+vf.name()+')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                (-rDeltaT.value())*rho.value()*
                    vf.oldTime().internalField()*mesh().S0()/mesh().S(),
                (-rDeltaT.value())*rho.value()*
                    vf.oldTime().boundaryField()
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                (-rDeltaT)*rho*vf.oldTime()
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
EulerFaDdtScheme<Type>::facDdt
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    rho.internalField()*vf.internalField()
                  - rho.oldTime().internalField()
                   *vf.oldTime().internalField()*mesh().S0()/mesh().S()
                ),
                rDeltaT.value()*
                (
                    rho.boundaryField()*vf.boundaryField()
                  - rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                rDeltaT*(rho*vf - rho.oldTime()*vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
EulerFaDdtScheme<Type>::facDdt0
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt0("+rho.name()+','+vf.name()+')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                  - rho.oldTime().internalField()
                   *vf.oldTime().internalField()*mesh().S0()/mesh().S()
                ),
                rDeltaT.value()*
                (
                  - rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                (-rDeltaT)*rho.oldTime()*vf.oldTime()
            )
        );
    }
}

template<class Type>
tmp<faMatrix<Type> >
EulerFaDdtScheme<Type>::famDdt
(
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type> > tfam
    (
        new faMatrix<Type>
        (
            vf,
            vf.dimensions()*dimArea/dimTime
        )
    );

    faMatrix<Type>& fam = tfam();

    scalar rDeltaT = 1.0/mesh().time().deltaT().value();

    fam.diag() = rDeltaT*mesh().S();
    
    if (mesh().moving())
    {
        fam.source() = rDeltaT*vf.oldTime().internalField()*mesh().S0();
    }
    else
    {
        fam.source() = rDeltaT*vf.oldTime().internalField()*mesh().S();
    }

    return tfam;
}


template<class Type>
tmp<faMatrix<Type> >
EulerFaDdtScheme<Type>::famDdt
(
    const dimensionedScalar& rho,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type> > tfam
    (
        new faMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimArea/dimTime
        )
    );
    faMatrix<Type>& fam = tfam();

    scalar rDeltaT = 1.0/mesh().time().deltaT().value();

    fam.diag() = rDeltaT*rho.value()*mesh().S();
    
    if (mesh().moving())
    {
        fam.source() = rDeltaT
            *rho.value()*vf.oldTime().internalField()*mesh().S0();
    }
    else
    {
        fam.source() = rDeltaT
            *rho.value()*vf.oldTime().internalField()*mesh().S();
    }

    return tfam;
}


template<class Type>
tmp<faMatrix<Type> >
EulerFaDdtScheme<Type>::famDdt
(
    const areaScalarField& rho,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type> > tfam
    (
        new faMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimArea/dimTime
        )
    );
    faMatrix<Type>& fam = tfam();

    scalar rDeltaT = 1.0/mesh().time().deltaT().value();

    fam.diag() = rDeltaT*rho.internalField()*mesh().S();

    if (mesh().moving())
    {
        fam.source() = rDeltaT
            *rho.oldTime().internalField()
            *vf.oldTime().internalField()*mesh().S0();
    }
    else
    {
        fam.source() = rDeltaT
            *rho.oldTime().internalField()
            *vf.oldTime().internalField()*mesh().S();
    }

    return tfam;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
