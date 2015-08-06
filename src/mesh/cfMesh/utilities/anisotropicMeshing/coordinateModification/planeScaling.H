/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Class
    planeScaling

Description
    Applies scaling to points on the positive side of the plane within
    the scaling distance.

SourceFiles
    planeScaling.C

\*---------------------------------------------------------------------------*/

#ifndef planeScaling_H
#define planeScaling_H

#include "coordinateModification.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class planeScaling Declaration
\*---------------------------------------------------------------------------*/

class planeScaling
:
    public coordinateModification
{
    // Private data
        //- origin of the plane
        point origin_;

        //- normal vector
        vector normal_;

        //- scaling distance
        scalar scalingDistance_;

        //- scaling factor
        scalar scalingFactor_;

public:

    //- Runtime type information
    TypeName("plane");


    // Constructors

        //- Null construct
        planeScaling();

        //- Construct from name, origin, normal, translation distance and scale
        planeScaling
        (
            const word& name,
            const point& original,
            const vector& normal,
            const scalar scalingDistance,
            const scalar scalingFactor
        );

        //- Construct from dictionary
        planeScaling(const word& name, const dictionary& dict);

        //- Construct and return a clone
        virtual autoPtr<coordinateModification> clone
        (
            const planeScaling& pt
        ) const
        {
            return autoPtr<coordinateModification>
            (
                new planeScaling
                (
                    pt.name(),
                    pt.origin_,
                    pt.normal_,
                    pt.scalingDistance_,
                    pt.scalingFactor_
                )
            );
        }

    // Member Functions

        //- return the origin of the plane
        virtual point origin() const;

        //- translate the object to the modified coordinates
        //- this is needed for backward transformation
        virtual void translateAndModifyObject(const vector&);

        //- calculate the displacement vector for plane translation
        virtual vector displacement(const point&) const;

        //- calculate the displacement vector for plane translation
        virtual vector backwardDisplacement(const point&) const;

        //- can this modification object be combined with other ones
        virtual bool combiningPossible() const;

        //- return that "bounding planes" of the scaling region for
        //- the given object
        virtual void boundingPlanes(PtrList<plane>&) const;

    //- Return as dictionary of entries
        dictionary dict(bool ignoreType = false) const;

      // Write

        //- Write
        void write(Ostream&) const;

        //- Write dictionary
        void writeDict(Ostream&, bool subDict = true) const;

    // Member Operators

        //- assign from dictionary
        void operator=(const dictionary&);

    // IOstream Operators

        Ostream& operator<<(Ostream&) const;
        friend Ostream& operator<<(Ostream&, const planeScaling&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
