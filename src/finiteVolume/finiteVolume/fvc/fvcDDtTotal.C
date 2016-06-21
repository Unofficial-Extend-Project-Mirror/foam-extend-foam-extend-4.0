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

Description


\*---------------------------------------------------------------------------*/

#include "fvcDDtTotal.H"
#include "fvcDiv.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
DDt
(
    const surfaceScalarField& phi,
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > ddtDivPhiPsi
        = fvc::ddt(psi) + fvc::div(phi, psi);

    if (phi.mesh().moving())
    {
        return ddtDivPhiPsi - fvc::div(phi + phi.mesh().phi())*psi;
    }
    else
    {
        return ddtDivPhiPsi - fvc::div(phi)*psi;
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
DDt
(
    const tmp<surfaceScalarField>& tphi,
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > DDtPsi
    (
        fvc::DDt(tphi(), psi)
    );
    tphi.clear();
    return DDtPsi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
