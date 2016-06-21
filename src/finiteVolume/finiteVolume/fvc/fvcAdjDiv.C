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
