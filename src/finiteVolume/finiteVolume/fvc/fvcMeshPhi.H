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
    Calculate the mesh motion flux and convert fluxes from absolute to relative
    and back.

SourceFiles
    fvcMeshPhi.C

\*---------------------------------------------------------------------------*/


#ifndef fvcMeshPhi_H
#define fvcMeshPhi_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    tmp<surfaceScalarField> meshPhi
    (
        const volVectorField& U
    );

    tmp<surfaceScalarField> meshPhi
    (
        const dimensionedScalar& rho,
        const volVectorField& U
    );

    tmp<surfaceScalarField> meshPhi
    (
        const volScalarField& rho,
        const volVectorField& U
    );


    //- Make the given flux relative
    void makeRelative
    (
        surfaceScalarField& phi,
        const volVectorField& U
    );

    //- Make the given flux relative
    void makeRelative
    (
        surfaceScalarField& phi,
        const dimensionedScalar& rho,
        const volVectorField& U
    );

    //- Make the given flux relative
    void makeRelative
    (
        surfaceScalarField& phi,
        const volScalarField& rho,
        const volVectorField& U
    );


    //- Make the given flux absolute
    void makeAbsolute
    (
        surfaceScalarField& phi,
        const volVectorField& U
    );

    //- Make the given flux absolute
    void makeAbsolute
    (
        surfaceScalarField& phi,
        const dimensionedScalar& rho,
        const volVectorField& U
    );

    //- Make the given flux absolute
    void makeAbsolute
    (
        surfaceScalarField& phi,
        const volScalarField& rho,
        const volVectorField& U
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
