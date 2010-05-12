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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Namespace
    fac

Description
    Calculate the lnGrad of the given areaField.

SourceFiles
    facLnGrad.C

\*---------------------------------------------------------------------------*/


#ifndef facLnGrad_H
#define facLnGrad_H

#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fac functions Declaration
\*---------------------------------------------------------------------------*/

namespace fac
{
    template<class Type>
    tmp<GeometricField<Type, faePatchField, edgeMesh> > lnGrad
    (
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word& name
    );

    template<class Type>
    tmp<GeometricField<Type, faePatchField, edgeMesh> > lnGrad
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&,
        const word& name
    );

    template<class Type>
    tmp<GeometricField<Type, faePatchField, edgeMesh> > lnGrad
    (
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faePatchField, edgeMesh> > lnGrad
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "facLnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
