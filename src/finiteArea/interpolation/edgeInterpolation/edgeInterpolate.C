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

Description

\*---------------------------------------------------------------------------*/

#include "edgeInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<edgeInterpolationScheme<Type> > scheme
(
    const edgeScalarField& faceFlux,
    Istream& streamData
)
{
    return edgeInterpolationScheme<Type>::New
    (
        faceFlux.mesh(),
        faceFlux,
        streamData
    );
}


// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<edgeInterpolationScheme<Type> > scheme
(
    const edgeScalarField& faceFlux,
    const word& name
)
{
    return edgeInterpolationScheme<Type>::New
    (
        faceFlux.mesh(),
        faceFlux,
        faceFlux.mesh().schemesDict().interpolationScheme(name)
    );
}


// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<edgeInterpolationScheme<Type> > scheme
(
    const faMesh& mesh,
    Istream& streamData
)
{
    return edgeInterpolationScheme<Type>::New
    (
        mesh,
        streamData
    );
}


// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<edgeInterpolationScheme<Type> > scheme
(
    const faMesh& mesh,
    const word& name
)
{
    return edgeInterpolationScheme<Type>::New
    (
        mesh,
        mesh.schemesDict().interpolationScheme(name)
    );
}


// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const edgeScalarField& faceFlux,
    Istream& schemeData
)
{
#   ifdef DEBUGInterpolations
    if (edgeInterpolation::debug)
    {
        Info<< "interpolate"
            << "(const GeometricField<Type, faPatchField, areaMesh>&, "
            << "const edgeScalarField&, Istream&) : "
            << "interpolating GeometricField<Type, faPatchField, areaMesh> "
            << endl;
    }
#   endif

    return scheme<Type>(faceFlux, schemeData)().interpolate(vf);
}


// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const edgeScalarField& faceFlux,
    const word& name
)
{
#   ifdef DEBUGInterpolations
    if (edgeInterpolation::debug)
    {
        Info<< "interpolate"
            << "(const GeometricField<Type, faPatchField, areaMesh>&, "
            << "const edgeScalarField&, const word&) : "
            << "interpolating GeometricField<Type, faPatchField, areaMesh> "
            << "using " << name
            << endl;
    }
#   endif

    return scheme<Type>(faceFlux, name)().interpolate(vf);
}

// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
interpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
    const edgeScalarField& faceFlux,
    const word& name
)
{
    tmp<GeometricField<Type, faePatchField, edgeMesh> > tsf =
        interpolate(tvf(), faceFlux, name);

    tvf.clear();

    return tsf;
}

// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const tmp<edgeScalarField>& tFaceFlux,
    const word& name
)
{
    tmp<GeometricField<Type, faePatchField, edgeMesh> > tsf =
        interpolate(vf, tFaceFlux(), name);

    tFaceFlux.clear();

    return tsf;
}

// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
interpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
    const tmp<edgeScalarField>& tFaceFlux,
    const word& name
)
{
    tmp<GeometricField<Type, faePatchField, edgeMesh> > tsf =
        interpolate(tvf(), tFaceFlux(), name);

    tvf.clear();
    tFaceFlux.clear();

    return tsf;
}


// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    Istream& schemeData
)
{
#   ifdef DEBUGInterpolations
    if (edgeInterpolation::debug)
    {
        Info<< "interpolate"
            << "(const GeometricField<Type, faPatchField, areaMesh>&, "
            << "Istream&) : "
            << "interpolating GeometricField<Type, faPatchField, areaMesh> "
            << endl;
    }
#   endif

    return scheme<Type>(vf.mesh(), schemeData)().interpolate(vf);
}

// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
#   ifdef DEBUGInterpolations
    if (edgeInterpolation::debug)
    {
        Info<< "interpolate"
            << "(const GeometricField<Type, faPatchField, areaMesh>&, "
            << "const word&) : "
            << "interpolating GeometricField<Type, faPatchField, areaMesh> "
            << "using " << name
            << endl;
    }
#   endif

    return scheme<Type>(vf.mesh(), name)().interpolate(vf);
}

// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
interpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faePatchField, edgeMesh> > tsf =
        interpolate(tvf(), name);

    tvf.clear();

    return tsf;
}


// Interpolate field onto faces using central differencing
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{

#   ifdef DEBUGInterpolations
    if (edgeInterpolation::debug)
    {
        Info<< "interpolate"
            << "(const GeometricField<Type, faPatchField, areaMesh>&) : "
            << "interpolating GeometricField<Type, faPatchField, areaMesh> "
            << "using run-time selected scheme"
            << endl;
    }
#   endif

    return interpolate(vf, "interpolate(" + vf.name() + ')');
}


// Interpolate field onto faces using central differencing
template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
interpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    tmp<GeometricField<Type, faePatchField, edgeMesh> > tsf =
        interpolate(tvf());
    tvf.clear();
    return tsf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
