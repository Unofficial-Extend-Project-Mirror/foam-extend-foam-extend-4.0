/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Namespace
    fac

Description

SourceFiles
    edgeInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef edgeInterpolate_H
#define edgeInterpolate_H

#include "tmp.H"
#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "edgeInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fa functions Declaration
\*---------------------------------------------------------------------------*/

namespace fac
{
    //- Return weighting factors for scheme given from Istream
    template<class Type>
    static tmp<edgeInterpolationScheme<Type> > scheme
    (
        const edgeScalarField& faceFlux,
        Istream& schemeData
    );

    //- Return weighting factors for scheme given by name in dictionary
    template<class Type>
    static tmp<edgeInterpolationScheme<Type> > scheme
    (
        const edgeScalarField& faceFlux,
        const word& name
    );


    //- Return weighting factors for scheme given from Istream
    template<class Type>
    static tmp<edgeInterpolationScheme<Type> > scheme
    (
        const faMesh& mesh,
        Istream& schemeData
    );

    //- Return weighting factors for scheme given by name in dictionary
    template<class Type>
    static tmp<edgeInterpolationScheme<Type> > scheme
    (
        const faMesh& mesh,
        const word& name
    );


    //- Interpolate field onto faces using scheme given by Istream
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, edgeMesh> > interpolate
    (
        const GeometricField<Type, faPatchField, areaMesh>& tvf,
        const edgeScalarField& faceFlux,
        Istream& schemeData
    );

    //- Interpolate field onto faces using scheme given by name in faSchemes
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, edgeMesh> > interpolate
    (
        const GeometricField<Type, faPatchField, areaMesh>& tvf,
        const edgeScalarField& faceFlux,
        const word& name
    );

    //- Interpolate field onto faces using scheme given by name in faSchemes
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, edgeMesh> > interpolate
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
        const edgeScalarField& faceFlux,
        const word& name
    );

    //- Interpolate field onto faces using scheme given by name in faSchemes
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, edgeMesh> > interpolate
    (
        const GeometricField<Type, faPatchField, areaMesh>& tvf,
        const tmp<edgeScalarField>& faceFlux,
        const word& name
    );

    //- Interpolate field onto faces using scheme given by name in faSchemes
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, edgeMesh> > interpolate
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
        const tmp<edgeScalarField>& faceFlux,
        const word& name
    );


    //- Interpolate field onto faces using scheme given by Istream
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, edgeMesh> > interpolate
    (
        const GeometricField<Type, faPatchField, areaMesh>& tvf,
        Istream& schemeData
    );

    //- Interpolate field onto faces using scheme given by name in faSchemes
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, edgeMesh> > interpolate
    (
        const GeometricField<Type, faPatchField, areaMesh>& tvf,
        const word& name
    );

    //- Interpolate field onto faces using scheme given by name in faSchemes
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, edgeMesh> > interpolate
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
        const word& name
    );


    //- Interpolate tmp field onto faces using central differencing
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, edgeMesh> > interpolate
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
    );

    //- Interpolate field onto faces using central differencing
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, edgeMesh> > interpolate
    (
        const GeometricField<Type, faPatchField, areaMesh>& tvf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "edgeInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
