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

InNamespace
    Foam::fvc

Description
    Surface Interpolation

SourceFiles
    surfaceInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceInterpolate_H
#define surfaceInterpolate_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    //- Return weighting factors for scheme given from Istream
    template<class Type>
    static tmp<surfaceInterpolationScheme<Type> > scheme
    (
        const surfaceScalarField& faceFlux,
        Istream& schemeData
    );

    //- Return weighting factors for scheme given by name in dictionary
    template<class Type>
    static tmp<surfaceInterpolationScheme<Type> > scheme
    (
        const surfaceScalarField& faceFlux,
        const word& name
    );


    //- Return weighting factors for scheme given from Istream
    template<class Type>
    static tmp<surfaceInterpolationScheme<Type> > scheme
    (
        const fvMesh& mesh,
        Istream& schemeData
    );

    //- Return weighting factors for scheme given by name in dictionary
    template<class Type>
    static tmp<surfaceInterpolationScheme<Type> > scheme
    (
        const fvMesh& mesh,
        const word& name
    );


    //- Interpolate field onto faces using scheme given by Istream
    template<class Type>
    static tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
    (
        const GeometricField<Type, fvPatchField, volMesh>& tvf,
        const surfaceScalarField& faceFlux,
        Istream& schemeData
    );

    //- Interpolate field onto faces using scheme given by name in fvSchemes
    template<class Type>
    static tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
    (
        const GeometricField<Type, fvPatchField, volMesh>& tvf,
        const surfaceScalarField& faceFlux,
        const word& name
    );

    //- Interpolate field onto faces using scheme given by name in fvSchemes
    template<class Type>
    static tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
    (
        const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
        const surfaceScalarField& faceFlux,
        const word& name
    );

    //- Interpolate field onto faces using scheme given by name in fvSchemes
    template<class Type>
    static tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
    (
        const GeometricField<Type, fvPatchField, volMesh>& tvf,
        const tmp<surfaceScalarField>& faceFlux,
        const word& name
    );

    //- Interpolate field onto faces using scheme given by name in fvSchemes
    template<class Type>
    static tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
    (
        const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
        const tmp<surfaceScalarField>& faceFlux,
        const word& name
    );


    //- Interpolate field onto faces using scheme given by Istream
    template<class Type>
    static tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
    (
        const GeometricField<Type, fvPatchField, volMesh>& tvf,
        Istream& schemeData
    );

    //- Interpolate field onto faces using scheme given by name in fvSchemes
    template<class Type>
    static tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
    (
        const GeometricField<Type, fvPatchField, volMesh>& tvf,
        const word& name
    );

    //- Interpolate field onto faces using scheme given by name in fvSchemes
    template<class Type>
    static tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
    (
        const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
        const word& name
    );


    //- Interpolate tmp field onto faces using central differencing
    template<class Type>
    static tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
    (
        const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
    );

    //- Interpolate field onto faces using central differencing
    template<class Type>
    static tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
    (
        const GeometricField<Type, fvPatchField, volMesh>& tvf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "surfaceInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
