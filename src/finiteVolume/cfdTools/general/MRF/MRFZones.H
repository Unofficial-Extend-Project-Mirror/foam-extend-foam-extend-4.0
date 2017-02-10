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

Class
    Foam::MRFZones

Description
    Container class for a set of MRFZones with the MRFZone member functions
    implemented to loop over the functions for each MRFZone.

SourceFiles
    MRFZones.C

\*---------------------------------------------------------------------------*/

#ifndef MRFZones_H
#define MRFZones_H

#include "MRFZone.H"
#include "IOPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class MRFZones Declaration
\*---------------------------------------------------------------------------*/

class MRFZones
:
    public IOPtrList<MRFZone>
{
    // Private data

        //- Mesh reference
        const fvMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        MRFZones(const MRFZones&);

        //- Disallow default bitwise assignment
        void operator=(const MRFZones&);


public:

    // Constructors

        //- Construct from fvMesh
        MRFZones(const fvMesh& mesh);


    // Member Functions

        //- Return raw correction flux
        tmp<surfaceScalarField> fluxCorrection() const;

        //- Return pseudo mesh flux
        tmp<surfaceScalarField> meshPhi() const;

        // Incompressible MRF

            //- Add the Coriolis force contribution to the momentum equation
            void addCoriolis(fvVectorMatrix& UEqn) const;

            //- Make the given absolute flux relative within the MRF region
            void relativeFlux(surfaceScalarField& phiAbs) const;

            //- Make the given relative flux absolute within the MRF region
            void absoluteFlux(surfaceScalarField& phiRel) const;


        // Compressible MRF

            //- Add the Coriolis force contribution to the momentum equation
            //  Use for equations with prefactor fields, eg. density.
            //  This is typically visible in the ddt term
            //  HJ, 12/Dec/2009
            void addCoriolis
            (
                const volScalarField& rho,
                fvVectorMatrix& UEqn
            ) const;

            //  Use for fluxes with prefactor fields, eg. density.
            //  For consistency, requires face-interpolated prefactor
            //  HJ. 12/Dec/2009

            //- Make the given relative mass-flux absolute within the
            //  MRF region
            void absoluteFlux
            (
                const surfaceScalarField& rho,
                surfaceScalarField& phi
            ) const;

            //- Make the given absolute flux relative within the MRF region
            void relativeFlux
            (
                const surfaceScalarField& rhof,
                surfaceScalarField& phiRel
            ) const;


        //- Make the given absolute velocity relative within the MRF region
        void relativeVelocity(volVectorField& Uabs) const;

        //- Make the given relative velocity absolute within the MRF region
        void absoluteVelocity(volVectorField& Urel) const;


        //- Correct the boundary velocity for the rotation of the MRF region
        void correctBoundaryVelocity(volVectorField& U) const;

        //- Compute source term for volScalarFields
        tmp<volScalarField> Su(const volScalarField& phi) const;

        //- Compute source term for volVectorFields
        tmp<volVectorField> Su(const volVectorField& phi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
