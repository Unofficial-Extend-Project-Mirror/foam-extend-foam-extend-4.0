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

#include "fvcAdjDiv.H"
#include "fvMesh.H"
#include "fvcSurfaceIntegrate.H"
#include "adjConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
adjDiv
(
    const volVectorField& Up,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::adjConvectionScheme<Type>::New
    (
        vf.mesh(),
        Up,
        vf.mesh().schemesDict().divScheme(name)
    )().fvcAdjDiv(Up, vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
adjDiv
(
    const tmp<volVectorField>& tUp,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > AdjDiv
    (
        fvc::adjDiv(tUp(), vf, name)
    );
    tUp.clear();
    return AdjDiv;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
adjDiv
(
    const volVectorField& Up,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > AdjDiv
    (
        fvc::adjDiv(Up, tvf(), name)
    );
    tvf.clear();
    return AdjDiv;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
adjDiv
(
    const tmp<volVectorField>& tUp,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > AdjDiv
    (
        fvc::adjDiv(tUp(), tvf(), name)
    );
    tUp.clear();
    tvf.clear();
    return AdjDiv;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
adjDiv
(
    const volVectorField& Up,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvc::adjDiv
    (
        Up, vf, "adjDiv(" + Up.name() + ',' + vf.name() + ')'
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
adjDiv
(
    const tmp<volVectorField>& tUp,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > AdjDiv
    (
        fvc::adjDiv(tUp(), vf)
    );
    tUp.clear();
    return AdjDiv;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
adjDiv
(
    const volVectorField& Up,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > AdjDiv
    (
        fvc::adjDiv(Up, tvf())
    );
    tvf.clear();
    return AdjDiv;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
adjDiv
(
    const tmp<volVectorField>& tUp,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > AdjDiv
    (
        fvc::adjDiv(tUp(), tvf())
    );
    tUp.clear();
    tvf.clear();
    return AdjDiv;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
