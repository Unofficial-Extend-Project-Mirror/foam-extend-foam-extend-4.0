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
#include "foamTime.H"

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
        HamiltonRodriguezRot eInitial_;

        //- Rotational tensor
        tensor rotTensor_;

        //- Incremental rotational tensor
        tensor rotIncrementTensor_;


    // Private Member Functions

        //- Calculate unit rotation vector from given rotation tensor
        static vector rotVector(const tensor& rotT);

        //- Calculate rotation angle from given rotation tensor
        static scalar rotAngle(const tensor& rotT);

        //- Calculate Euler angles (x-y-z (roll-pitch-yaw) convention) from
        //  given rotation tensor. Reference: Mike Day, Insomniac Games,
        //  Extracting Euler Angles from a Rotation Matrix.
        static vector eulerAngles(const tensor& rotT);


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

        //- Return Euler angles (see static function eulerAngles for convention)
        vector eulerAngles() const;

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
