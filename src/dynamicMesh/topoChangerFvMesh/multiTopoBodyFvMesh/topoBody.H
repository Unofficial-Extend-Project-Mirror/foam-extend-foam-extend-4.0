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
    topoBody

Description
    Single moving body with layer addition/removal interfaces around it

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    topoBody.C

\*---------------------------------------------------------------------------*/

#ifndef topoBody_H
#define topoBody_H

#include "polyMesh.H"
#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyTopoChanger;

/*---------------------------------------------------------------------------*\
                            Class topoBody Declaration
\*---------------------------------------------------------------------------*/

class topoBody
{
    // Private data

        //- Name
        const word name_;

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Name of moving cell zone
        const word movingCellsName_;

        //- Name of layer face zone
        const wordList layerFacesNames_;

        //- Min layering thickness
        const scalar minThickness_;

        //- Max layering thickness
        const scalar maxThickness_;

        //- Motion control function
        autoPtr<solidBodyMotionFunction> SBMFPtr_;

        //- Invert motion mask: marked region is stationary and the rest
        //  of the domain is moving
        const Switch invertMotionMask_;

        //- Markup field for points.  Moving points marked with 1
        mutable scalarField* movingPointsMaskPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        topoBody(const topoBody&);

        //- Disallow default bitwise assignment
        void operator=(const topoBody&);


        //- Calculate moving mask
        void calcMovingMask() const;

        //- Return moving points mask
        const scalarField& movingPointsMask() const;

        //- Clear moving points mask
        void clearPointMask();


public:

    // Constructors

        //- Construct from  dictionary
        topoBody
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~topoBody();


    // Member Functions

        //- Return accumulative point motion
        virtual tmp<vectorField> pointMotion() const;

        //- Add layering zones
        void addZones
        (
            DynamicList<pointZone*>& pz,
            DynamicList<faceZone*>& fz,
            DynamicList<cellZone*>& cz
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
