/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    Foam::faceMesh

Description
    Storage for surface mesh i.e. points and faces.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef faceMesh_H
#define faceMesh_H

#include "PrimitivePatchTemplate.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class faceMesh Declaration
\*---------------------------------------------------------------------------*/

class faceMesh
:
    public PrimitivePatch<face, Foam::List, pointField>
{
    // Private member functions

        PrimitivePatch<face, Foam::List, pointField> read(Istream& is)
        {
            pointField points(is);
            faceList faces(is);
            return PrimitivePatch<face, Foam::List, pointField>(faces, points);
        }


public:

    // Constructors

        //- Construct from components
        faceMesh(const faceList& faces, const pointField& points)
        :
            PrimitivePatch<face, Foam::List, pointField>(faces, points)
        {}

        //- Construct from Istream
        faceMesh(Istream& is)
        :
            PrimitivePatch<face, Foam::List, pointField>(read(is))
        {}


    // Member Functions

        void flip()
        {
            forAll(*this, i)
            {
                face& f = operator[](i);
                f = f.reverseFace();
            }
            clearOut();
        }

    // IOstream Operators

        friend Ostream& operator<<(Ostream& os, const faceMesh& fm)
        {
            return os
                << fm.points()
                << token::NL
                << static_cast<PrimitivePatch<face, Foam::List, pointField> >
                   (fm);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
