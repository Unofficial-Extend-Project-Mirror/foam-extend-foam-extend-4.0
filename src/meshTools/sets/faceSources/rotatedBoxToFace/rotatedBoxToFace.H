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
    Foam::rotatedBoxToFace

Description
    A topoSetSource to select faces based on face centres inside
    rotated/skewed box (parallelopiped?).

    Box defined as origin and i,j,k vectors.
    E.g. box rotated 45 degrees around z-axis with sizes sqrt(0.2^2+0.2^2)
    (and extra large, 200 in z direction):
    @verbatim
       origin   ( 0.4 0.4 -100);
       i        ( 0.2 0.2    0);
       j        (-0.2 0.2    0);
       k        ( 0.0 0.0  100);
    @endverbatim

SourceFiles
    rotatedBoxToFace.C

\*---------------------------------------------------------------------------*/

#ifndef rotatedBoxToFace_H
#define rotatedBoxToFace_H

#include "topoSetSource.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class rotatedBoxToFace Declaration
\*---------------------------------------------------------------------------*/

class rotatedBoxToFace
:
    public topoSetSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;


        //- skewed box
        const vector origin_;
        const vector i_;
        const vector j_;
        const vector k_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("rotatedBoxToFace");

    // Constructors

        //- Construct from components
        rotatedBoxToFace
        (
            const polyMesh& mesh,
            const vector& origin,
            const vector& i,
            const vector& j,
            const vector& k
        );

        //- Construct from dictionary
        rotatedBoxToFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        rotatedBoxToFace(const polyMesh& mesh, Istream&);


    // Destructor

        virtual ~rotatedBoxToFace();


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

