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
    boxRefinement

Description
    Checks if a box is contained inside the box object

SourceFiles
    boxRefinement.C

\*---------------------------------------------------------------------------*/

#ifndef boxRefinement_H
#define boxRefinement_H

#include "objectRefinement.H"
#include "point.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class boxRefinement Declaration
\*---------------------------------------------------------------------------*/

class boxRefinement
:
    public objectRefinement
{
    // Private data
        //- centre of the box
        point centre_;
    
        //- length of box sides
        scalar lengthX_;
        scalar lengthY_;
        scalar lengthZ_;
    
        //- angles to the global coordinate system
        //scalar angleX_;
        //scalar angleY_;
        //scalar angleZ_;

public:

    //- Runtime type information
    TypeName("box");


    // Constructors

        //- Null construct
        boxRefinement();

        //- Construct from name, cell size, centre and sizes
        boxRefinement
        (
            const word& name,
            const scalar cellSize,
            const point& centre,
            const scalar lengthX,
            const scalar lengthY,
            const scalar lengthZ
            //const scalar angleX = 0.,
            //const scalar angleY = 0.,
            //const scalar angleZ = 0.
        );

        //- Construct from dictionary
        boxRefinement(const word& name, const dictionary& dict);
        
        //- Construct and return a clone
        virtual autoPtr<objectRefinement> clone
        (
            const boxRefinement& sr
        ) const
        {
            return autoPtr<objectRefinement>
            (
                new boxRefinement
                (
                    sr.name(),
                    sr.cellSize(),
                    sr.centre_,
                    sr.lengthX_,
                    sr.lengthY_,
                    sr.lengthZ_
                    //sr.angleX_,
                    //sr.angleY_,
                    //sr.angleZ_
                )
            );
        }
    
    // Member Functions
    
        //- check if a boundBox intersects or is inside the object
        bool intersectsObject(const boundBox&) const;
        
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
