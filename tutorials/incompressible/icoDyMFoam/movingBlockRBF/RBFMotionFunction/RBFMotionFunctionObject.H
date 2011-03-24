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
    RBFMotionFunctionObject

Description
    FunctionObject specifying motion of control points using RBF

Author
    Frank Bos, TU Delft.  All rights reserved.

SourceFiles
    RBFMotionFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef RBFMotionFunctionObject_H
#define RBFMotionFunctionObject_H

#include "functionObject.H"
#include "dictionary.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{


/*---------------------------------------------------------------------------*\
                     Class RBFMotionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class RBFMotionFunctionObject
:
    public functionObject
{
    // Private data

        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

        scalar rotationAmplitude_;
        scalar rotationFrequency_;
        vector translationAmplitude_;
        vector translationFrequency_;
        vector initialRotationOrigin_;

        //- Stationary points
        pointField statPoints_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        RBFMotionFunctionObject(const RBFMotionFunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const RBFMotionFunctionObject&);


public:

    //- Runtime type information
    TypeName("RBFMotion");


    // Constructors

        //- Construct from components
        RBFMotionFunctionObject
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute();

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
