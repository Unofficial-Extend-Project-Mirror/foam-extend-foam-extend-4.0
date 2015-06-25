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
    hollowConeRefinement

Description
    Refine boxes with centres contained inside the specified hollow cone

SourceFiles
    hollowConeRefinement.C

\*---------------------------------------------------------------------------*/

#ifndef hollowConeRefinement_H
#define hollowConeRefinement_H

#include "objectRefinement.H"
#include "point.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class coneRefinement Declaration
\*---------------------------------------------------------------------------*/

class hollowConeRefinement
:
    public objectRefinement
{
    // Private data
        //- start point and the radii
        point p0_;
        scalar r0Outer_;
        scalar r0Inner_;
    
        //- end point and the radii
        point p1_;
        scalar r1Outer_;
        scalar r1Inner_;

public:

    //- Runtime type information
    TypeName("hollowCone");


    // Constructors

        //- Null construct
        hollowConeRefinement();

        //- Construct from name, cell size, centre and radii
        hollowConeRefinement
        (
            const word& name,
            const scalar cellSize,
            const direction additionalRefLevels,
            const point& p0,
            const scalar radius0Outer,
            const scalar radius0Inner,
            const point& p1,
            const scalar radius1Outer,
            const scalar radius1Inner
        );

        //- Construct from dictionary
        hollowConeRefinement(const word& name, const dictionary& dict);
        
        //- Construct and return a clone
        virtual autoPtr<objectRefinement> clone
        (
            const hollowConeRefinement& sr
        ) const
        {
            return autoPtr<objectRefinement>
            (
                new hollowConeRefinement
                (
                    sr.name(),
                    sr.cellSize(),
                    sr.additionalRefinementLevels(),
                    sr.p0_,
                    sr.r0Outer_,
                    sr.r0Inner_,
                    sr.p1_,
                    sr.r1Outer_,
                    sr.r1Inner_
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
