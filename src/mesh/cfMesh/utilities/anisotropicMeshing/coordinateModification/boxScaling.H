/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    boxScaling

Description
    Checks if a box is contained inside the box object

SourceFiles
    boxScaling.C

\*---------------------------------------------------------------------------*/

#ifndef boxScaling_H
#define boxScaling_H

#include "coordinateModification.H"
#include "point.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class boxScaling Declaration
\*---------------------------------------------------------------------------*/

class boxScaling
:
    public coordinateModification
{
    // Private data
        //- centre of the box
        point centre_;

        //- length of box sides
        vector lengthVec_;

        //- scaling factors in all directions
        vector scaleVec_;

        //- min point
        point pMin_;

        //- max point
        point pMax_;

    // Private member functions
        //- calculate bounding box points
        void calculateBndBox();

public:

    //- Runtime type information
    TypeName("box");


    // Constructors

        //- Null construct
        boxScaling();

        //- Construct from name, cell size, centre and sizes and scaling factors
        boxScaling
        (
            const word& name,
            const point& centre,
            const scalar lengthX,
            const scalar lengthY,
            const scalar lengthZ,
            const scalar scaleX = 1.0,
            const scalar scaleY = 1.0,
            const scalar scaleZ = 1.0
        );

        //- Construct from dictionary
        boxScaling(const word& name, const dictionary& dict);

        //- Construct and return a clone
        virtual autoPtr<coordinateModification> clone
        (
            const boxScaling& bs
        ) const
        {
            return autoPtr<coordinateModification>
            (
                new boxScaling
                (
                    bs.name(),
                    bs.centre_,
                    bs.lengthVec_.x(),
                    bs.lengthVec_.y(),
                    bs.lengthVec_.z(),
                    bs.scaleVec_.x(),
                    bs.scaleVec_.y(),
                    bs.scaleVec_.z()
                )
            );
        }

    // Member Functions

        //- return the centre of the box
        virtual point origin() const;

        //- translate the object to the modified coordinates
        //- this is needed for backward transformation
        virtual void translateAndModifyObject(const vector&);

        //- calculate the displacement vector for box scaling
        virtual vector displacement(const point&) const;

        //- calculate the displacement vector for box scaling
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
        friend Ostream& operator<<(Ostream&, const boxScaling&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
