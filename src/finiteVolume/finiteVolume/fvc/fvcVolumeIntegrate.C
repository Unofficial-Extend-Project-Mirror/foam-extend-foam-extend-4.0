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

#include "fvcVolumeIntegrate.H"
#include "fvMesh.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> >
volumeIntegrate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return vf.mesh().V()*vf.internalField();
}

template<class Type>
tmp<Field<Type> >
volumeIntegrate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    tmp<Field<Type> > tvivf = tvf().mesh().V()*tvf().internalField();
    tvf.clear();
    return tvivf;
}


template<class Type>
dimensioned<Type>
domainIntegrate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return dimensioned<Type>
    (
        "domainIntegrate(" + vf.name() + ')',
        dimVol*vf.dimensions(),
        gSum(fvc::volumeIntegrate(vf))
    );
}

template<class Type>
dimensioned<Type>
domainIntegrate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    dimensioned<Type> integral = domainIntegrate(tvf());
    tvf.clear();
    return integral;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
