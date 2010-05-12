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

#include "steadyStateFaDdtScheme.H"
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
steadyStateFaDdtScheme<Type>::facDdt
(
    const dimensioned<Type> dt
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh> >
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "ddt("+dt.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                dt.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
steadyStateFaDdtScheme<Type>::facDdt0
(
    const dimensioned<Type> dt
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh> >
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "ddt("+dt.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                dt.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
steadyStateFaDdtScheme<Type>::facDdt
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh> >
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "ddt("+vf.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                vf.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
steadyStateFaDdtScheme<Type>::facDdt0
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh> >
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "ddt0("+vf.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                vf.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}


template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
steadyStateFaDdtScheme<Type>::facDdt0
(
    const GeometricField<Type, faePatchField, edgeMesh>& vf
)
{
    return tmp<GeometricField<Type, faePatchField, edgeMesh> >
    (
        new GeometricField<Type, faePatchField, edgeMesh>
        (
            IOobject
            (
                "ddt0("+vf.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                vf.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
steadyStateFaDdtScheme<Type>::facDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh> >
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "ddt("+rho.name()+','+vf.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                rho.dimensions()*vf.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
steadyStateFaDdtScheme<Type>::facDdt0
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh> >
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "ddt0("+rho.name()+','+vf.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                rho.dimensions()*vf.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
steadyStateFaDdtScheme<Type>::facDdt
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh> >
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "ddt("+rho.name()+','+vf.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                rho.dimensions()*vf.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
steadyStateFaDdtScheme<Type>::facDdt0
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh> >
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "ddt0("+rho.name()+','+vf.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                rho.dimensions()*vf.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}

template<class Type>
tmp<faMatrix<Type> >
steadyStateFaDdtScheme<Type>::famDdt
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

    return tfam;
}


template<class Type>
tmp<faMatrix<Type> >
steadyStateFaDdtScheme<Type>::famDdt
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

    return tfam;
}


template<class Type>
tmp<faMatrix<Type> >
steadyStateFaDdtScheme<Type>::famDdt
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

    return tfam;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
