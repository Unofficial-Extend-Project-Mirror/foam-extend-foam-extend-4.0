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
    mixerRotor

Description
    Single rotating region

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    mixerRotor.C

\*---------------------------------------------------------------------------*/

#ifndef mixerRotor_H
#define mixerRotor_H

#include "cylindricalCS.H"
#include "Switch.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class regionSplit;
class polyTopoChanger;

/*---------------------------------------------------------------------------*\
                        Class mixerRotor Declaration
\*---------------------------------------------------------------------------*/

class mixerRotor
{
    // Private data

        //- Name
        const word name_;

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Coordinate system
        cylindricalCS cs_;

        // - Rotational speed in rotations per minute (rpm)
        const scalar rpm_;

        //- Name of moving slider patch
        const word movingSliderName_;

        //- Name of static slider patch
        const word staticSliderName_;

        //- Rotating region marker point
        point rotatingRegionMarker_;

        //- Invert motion mask: marked region is stationary and the rest
        //  of the domain is moving
        const Switch invertMotionMask_;

        //- Use topological sliding.  Other option is GGI
        const Switch useTopoSliding_;

        //- Attach-detach action
        Switch attachDetach_;

        //- Markup field for points.  Moving points marked with 1
        mutable scalarField* movingPointsMaskPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        mixerRotor(const mixerRotor&);

        //- Disallow default bitwise assignment
        void operator=(const mixerRotor&);


        //- Return coordinate system
        const cylindricalCS& cs() const
        {
            return cs_;
        }

        //- Calculate moving mask
        void calcMovingMask() const;

        //- Return moving points mask
        const scalarField& movingPointsMask() const;

        //- Clear moving points mask
        void clearPointMask();


public:

    // Constructors

        //- Construct from  dictionary
        mixerRotor
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~mixerRotor();


    // Member Functions

        //- Return true if topological sliding is used
        bool useTopoSliding() const
        {
            return useTopoSliding_;
        }

        //- Return true if attach-detach action is used
        bool attachDetach() const
        {
            return attachDetach_;
        }

        //- Return accumulative point motion
        virtual tmp<vectorField> pointMotion() const;

        //- Add mixer rotor zones
        void addZones
        (
            DynamicList<pointZone*>& pz,
            DynamicList<faceZone*>& fz,
            DynamicList<cellZone*>& cz,
            const regionSplit& rs
        );

        //- Add mixer rotor  modifiers
        void addModifiers
        (
            polyTopoChanger& tc,
            label& nextI
        );

        //- Update topology
        void updateTopology();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
