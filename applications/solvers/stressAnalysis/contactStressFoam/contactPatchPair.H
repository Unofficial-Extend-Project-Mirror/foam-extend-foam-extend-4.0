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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    contactPatchPair

Description
    
SourceFiles
    contactPatchPair.C
    contactPatchPairSlavePressure.C
    contactPatchPairUpdateContact.C

\*---------------------------------------------------------------------------*/

#ifndef contactPatchPair_H
#define contactPatchPair_H

#include "fvMesh.H"
#include "primitivePatchInterpolation.H"
#include "patchToPatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class contactPatchPair Declaration
\*---------------------------------------------------------------------------*/

class contactPatchPair
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Index of the directionMixed patch
        const label masterPatchIndex_;

        //- Index of the gradient patch
        const label slavePatchIndex_;

        //- Master patch interpolator
        primitivePatchInterpolation masterInterpolate_;

        //- Slave patch interpolator
        primitivePatchInterpolation slaveInterpolate_;

        //- Patch to patch interpolator
        patchToPatchInterpolation patchToPatchInterpolate_;

        //- Contact tolerance
        const scalar tol_;

        //- Touch fraction
        scalarField touchFraction_;

        //- Slave displacement
        vectorField slaveDisplacement_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        contactPatchPair(const contactPatchPair&);

        //- Disallow default bitwise assignment
        void operator=(const contactPatchPair&);


public:

    // Constructors

        //- Construct from components
        contactPatchPair
        (
            const fvMesh& m,
            const label master,
            const label slave,
            const scalar tolerance
        );


    // Member Functions

        //- Update contact data
        void updateContact
        (
            const volVectorField& disp
        );

        //- Return contact fraction for patch
        const scalarField& touchFraction() const
        {
            return touchFraction_;
        }

        //- Return slave patch displacement
        const vectorField& slaveDisplacement() const
        {
            return slaveDisplacement_;
        }

        //- Return slave patch pressure
        tmp<scalarField> slavePressure
        (
            const scalarField& traction
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
