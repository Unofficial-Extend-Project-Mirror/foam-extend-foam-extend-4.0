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
    tetIntersection

Description
    Implementation of the tetrahedron intersection algorithm given in:

    D.H. Eberly, 3D Game Engine Design: A Practical Approach to Real-time
    Computer Graphics, Morgan Kaufmann, 2001.

    Geometric Tools, LLC
    Distributed under the Boost Software License, Version 1.0.
    http://www.boost.org/LICENSE_1_0.txt

Implemented by
    Sandeep Menon
    University of Massachusetts Amherst

SourceFiles
    tetIntersectionI.H

\*---------------------------------------------------------------------------*/

#ifndef tetIntersection_H
#define tetIntersection_H

#include "point.H"
#include "label.H"
#include "Tuple2.H"
#include "FixedList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class tetIntersection Declaration
\*---------------------------------------------------------------------------*/

class tetIntersection
{
    // Private data

        //- Const reference to clipping tetrahedron
        const FixedList<point, 4>& clipTet_;

        //- Hessian-normal plane definition
        typedef Tuple2<vector, scalar> hPlane;

        FixedList<hPlane, 4> clipPlanes_;

        //- Magnitude of clipping tetrahedron
        scalar clipTetMag_;

        //- Tetrahedra used as temporaries
        DynamicList<FixedList<point, 4> > inside_;

        //- All intersection tets
        DynamicList<FixedList<point, 4> > allTets_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetIntersection(const tetIntersection&);

        //- Disallow default bitwise assignment
        void operator=(const tetIntersection&);

        //- Compute clip-planes
        inline void computeClipPlanes();

        //- Split and decompose
        inline void splitAndDecompose
        (
            const label tetPlaneIndex,
            FixedList<point, 4>& tetra,
            DynamicList<FixedList<point, 4> >& decompTets
        ) const;

public:

    // Constructors

        //- Construct from components
        inline tetIntersection(const FixedList<point, 4>& clipTet);


    // Destructor

        inline ~tetIntersection();


    // Member Functions

        //- Return magnitude of clipping tetrahedron
        inline scalar clipTetMag() const;

        //- Evaluate for intersections against input tetrahedron
        inline bool evaluate(const FixedList<point, 4>& subjectTet);

        //- Return intersections
        inline const DynamicList<FixedList<point, 4> >& getIntersection() const;

        //- Evaluate and return volume / centroid
        inline void getVolumeAndCentre(scalar& volume, vector& centre) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tetIntersectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
