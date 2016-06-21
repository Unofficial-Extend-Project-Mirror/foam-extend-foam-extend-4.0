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
    Foam::cloud

Description
    A cloud is a collection of lagrangian particles

SourceFiles
    cloud.C

\*---------------------------------------------------------------------------*/

#ifndef cloud_H
#define cloud_H

#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                            Class cloud Declaration
\*---------------------------------------------------------------------------*/

class cloud
:
    public objectRegistry
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        cloud(const cloud&);

        //- Disallow default bitwise assignment
        void operator=(const cloud&);


public:

    //- Runtime type information
    TypeName("cloud");

        //- The prefix to local: %lagrangian
        static const word prefix;

        //- The default cloud name: %defaultCloud
        static word defaultName;


    // Constructors

        //- Construct for the given objectRegistry and named cloud instance
        cloud(const objectRegistry&, const word& cloudName = "");


    //- Destructor
        virtual ~cloud();


    // Member Functions

        // Edit

            //- Remap the cells of particles corresponding to the
            //  mesh topology change
            virtual void autoMap(const mapPolyMesh&) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
