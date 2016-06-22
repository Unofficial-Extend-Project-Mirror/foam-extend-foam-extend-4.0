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
    Foam::ensightPartFaces

Description
    An implementation of ensightPart to hold volume mesh faces.

SourceFiles
    ensightPartFaces.C

\*---------------------------------------------------------------------------*/

#ifndef ensightPartFaces_H
#define ensightPartFaces_H

#include "ensightPart.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ensightPartFaces Declaration
\*---------------------------------------------------------------------------*/

class ensightPartFaces
:
    public ensightPart
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const ensightPartFaces&);

        //- track points used
        virtual localPoints calcLocalPoints() const;

        //- element connectivity
        virtual void writeConnectivity
        (
            ensightGeoFile& os,
            const string& key,
            const labelList& idList,
            const labelList& pointMap
        ) const;


protected:

    //- addressable ensight element types
    enum elemType
    {
        tria3Elements,
        quad4Elements,
        nsidedElements
    };

    // Static data members
    static List<word> elemTypes_;


public:

    //- Runtime type information
    TypeName("ensightFaces");

    // Constructors

        //- Construct empty part with number and description
        ensightPartFaces(label partNumber, const string& partDescription);

        //- Construct from polyMesh and polyPatch
        ensightPartFaces
        (
            label partNumber,
            const polyMesh&,
            const polyPatch&
        );

        //- Construct as copy
        ensightPartFaces(const ensightPartFaces&);

        //- Construct from Istream
        ensightPartFaces(Istream&);

        //- Construct on freestore from Istream
        static autoPtr<ensightPartFaces> New(Istream& is)
        {
            return autoPtr<ensightPartFaces>(new ensightPartFaces(is));
        }


    //- Destructor
    virtual ~ensightPartFaces();


    // Member Functions

        //- static listing of the element types
        virtual List<word> const& elementTypes() const
        {
            return elemTypes_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
