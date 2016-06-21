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
    Foam::upwindCFCCellToFaceStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef upwindCFCCellToFaceStencilObject_H
#define upwindCFCCellToFaceStencilObject_H

#include "extendedUpwindCellToFaceStencil.H"
#include "CFCCellToFaceStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class upwindCFCCellToFaceStencilObject Declaration
\*---------------------------------------------------------------------------*/

class upwindCFCCellToFaceStencilObject
:
    public MeshObject<fvMesh, upwindCFCCellToFaceStencilObject>,
    public extendedUpwindCellToFaceStencil
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        upwindCFCCellToFaceStencilObject
        (
            const upwindCFCCellToFaceStencilObject&
        );

        //- Disallow default bitwise assignment
        void operator=(const upwindCFCCellToFaceStencilObject&);


public:

    TypeName("upwindCFCCellToFaceStencil");

    // Constructors

        //- Construct from uncompacted face stencil
        explicit upwindCFCCellToFaceStencilObject
        (
            const fvMesh& mesh,
            const bool pureUpwind,
            const scalar minOpposedness
        )
        :
            MeshObject<fvMesh, upwindCFCCellToFaceStencilObject>(mesh),
            extendedUpwindCellToFaceStencil
            (
                CFCCellToFaceStencil(mesh),
                pureUpwind,
                minOpposedness
            )
        {
            if (extendedCellToFaceStencil::debug)
            {
                Info<< "Generated off-centred stencil " << type()
                    << nl << endl;
                writeStencilStats(Info, ownStencil(), ownMap());
            }
        }


    // Destructor

        virtual ~upwindCFCCellToFaceStencilObject()
        {}


    // Member functions

        //- Delete the stencil when the mesh moves
        virtual bool movePoints() const
        {
            deleteObject();
            return true;
        }

        //- Delete the stencil when mesh updates
        virtual bool updateMesh(const mapPolyMesh&) const
        {
            deleteObject();
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
