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
    Foam::ensightPartCells

Description
    An implementation of ensightPart to hold volume mesh cells.

SourceFiles
    ensightPartCells.C

\*---------------------------------------------------------------------------*/

#ifndef ensightPartCells_H
#define ensightPartCells_H

#include "ensightPart.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ensightPartCells Declaration
\*---------------------------------------------------------------------------*/

class ensightPartCells
:
    public ensightPart
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const ensightPartCells&);

        //- classify the cell types
        void classify(const labelList& idLabels = labelList::null());

        //- track points used
        virtual localPoints calcLocalPoints() const;

        //- track the points used
        // virtual void makeLocalPointMap();

        //- element connectivity
        virtual void writeConnectivity
        (
            ensightGeoFile& os,
            const string& key,
            const labelList& idList,
            const labelList& pointMap
        ) const;


protected:

    //- addressable Ensight element types
    enum elemType
    {
        tetra4Elements,
        pyramid5Elements,
        penta6Elements,
        hexa8Elements,
        nfacedElements
    };

    // Static data members
    static List<word> elemTypes_;


public:

    //- Runtime type information
    TypeName("ensightCells");

    // Constructors

        //- Construct empty part with number and description
        ensightPartCells(label partNumber, const string& partDescription);

        //- Construct from polyMesh without zones
        ensightPartCells(label partNumber, const polyMesh&);

        //- Construct from polyMesh and list of (non-zoned) cells
        ensightPartCells
        (
            label partNumber,
            const polyMesh&,
            const labelList&
        );

        //- Construct from polyMesh and cellZone
        ensightPartCells
        (
            label partNumber,
            const polyMesh&,
            const cellZone&
        );

        //- Construct as copy
        ensightPartCells(const ensightPartCells&);

        //- Construct from Istream
        ensightPartCells(Istream&);

        //- Construct on freestore from Istream
        static autoPtr<ensightPartCells> New(Istream& is)
        {
            return autoPtr<ensightPartCells>(new ensightPartCells(is));
        }


    //- Destructor
    virtual ~ensightPartCells();


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
