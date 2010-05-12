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
    finiteRotation

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.
    Updated by Hrvoje Jasak

SourceFiles
    finiteRotation.C

\*---------------------------------------------------------------------------*/

#ifndef finiteRotation_H
#define finiteRotation_H

#include "HamiltonRodriguezRot.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                         Class finiteRotation Declaration
\*---------------------------------------------------------------------------*/

class finiteRotation
{
    // Private data

        //- Initial rotation
        const HamiltonRodriguezRot eInitial_;

        //- Rotational tensor
        tensor rotTensor_;

        //- Incremental rotational tensor
        tensor rotIncrementTensor_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        finiteRotation(const finiteRotation&);

        //- Disallow default bitwise assignment
        void operator=(const finiteRotation&);


        //- Calculate unit rotation vector from given rotation tensor
        static vector rotVector(const tensor& rotT);

        //- Calculate rotation angle from given rotation tensor
        static scalar rotAngle(const tensor& rotT);


public:

    // Constructors

        //- Construct from cooordinate rotation
        explicit finiteRotation(const HamiltonRodriguezRot& rot);

        //- Construct from rotation vector and angle
        finiteRotation
        (
            const vector& r,
            const scalar& angle
        );


    // Destructor

        ~finiteRotation();


    // Member Functions

        //- Update rotation
        void updateRotation(const HamiltonRodriguezRot& rot);

        //- Return initial quaternions
        const HamiltonRodriguezRot& eInitial() const;

        //- Return current quaternions
        HamiltonRodriguezRot eCurrent() const;

        //-Transformations to inertial coordinate system
        const tensor& rotTensor() const;

        //- Return rotation vector
        vector rotVector() const;

        //- Return rotation angle
        scalar rotAngle() const;

        // Transformations between new and previous rotation
        const tensor& rotIncrementTensor() const;

        //- Return rotation vector increment
        vector rotIncrementVector() const;

        //- Return rotation angle increment
        scalar rotIncrementAngle() const;

        // Return average rotational velocity vector
        vector omegaAverage(const scalar deltaT) const;

        // Return mean rotation for last time step
        tensor rotTensorAverage() const;

        // Return mean rotation vector for last time step
        vector rotVectorAverage() const;

        // Return average rotational velocity vector in absolute coord. system
        vector omegaAverageAbsolute(const scalar deltaT) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
