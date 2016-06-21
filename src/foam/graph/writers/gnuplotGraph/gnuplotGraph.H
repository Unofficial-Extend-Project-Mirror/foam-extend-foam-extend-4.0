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
    Foam::gnuplotGraph

Description
    Output in gnuplot (http://www.gnuplot.info) format

SourceFiles
    gnuplotGraph.C

\*---------------------------------------------------------------------------*/

#ifndef gnuplotGraph_H
#define gnuplotGraph_H

#include "graph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class gnuplotGraph Declaration
\*---------------------------------------------------------------------------*/

class gnuplotGraph
:
    public graph::writer
{

public:

    //- Runtime type information
    TypeName("gnuplot");

    //- FileName extension  for this graph format
    static const word ext_;


    // Constructors

        //- Construct null
        gnuplotGraph()
        {}


    // Destructor

        ~gnuplotGraph()
        {}


    // Member Functions

        // Access

            //- Return the appropriate fileName extension
            //  for this graph format
            const word& ext() const
            {
                return ext_;
            }


        // Write

            void write(const graph&, Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
