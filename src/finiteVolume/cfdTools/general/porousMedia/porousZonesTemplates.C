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

#include "porousZones.H"
#include "volFields.H"
#include "fvm.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::porousZones::modifyDdt(fvMatrix<Type>& m) const
{
    forAll(*this, i)
    {
        operator[](i).modifyDdt(m);
    }
}


// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::porousZones::ddt
(
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tres = fvm::ddt(vf);
    modifyDdt(tres());
    return tres;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::porousZones::ddt
(
    const geometricOneField&,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tres = fvm::ddt(vf);
    modifyDdt(tres());
    return tres;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::porousZones::ddt
(
    const dimensionedScalar& rho,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tres = fvm::ddt(rho,vf);
    modifyDdt(tres());
    return tres;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::porousZones::ddt
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tres = fvm::ddt(rho,vf);
    modifyDdt(tres());
    return tres;
}

// ************************************************************************* //
