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

#ifndef pureUpwindCFCCellToFaceStencilObject_H
#define pureUpwindCFCCellToFaceStencilObject_H

#include "extendedUpwindCellToFaceStencil.H"
#include "CFCCellToFaceStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class pureUpwindCFCCellToFaceStencilObject Declaration
\*---------------------------------------------------------------------------*/

class pureUpwindCFCCellToFaceStencilObject
:
    public MeshObject<fvMesh, pureUpwindCFCCellToFaceStencilObject>,
    public extendedUpwindCellToFaceStencil
{

public:

    TypeName("pureUpwindCFCCellToFaceStencil");

    // Constructors

        //- Construct from uncompacted face stencil
        explicit pureUpwindCFCCellToFaceStencilObject
        (
            const fvMesh& mesh
        )
        :
            MeshObject<fvMesh, pureUpwindCFCCellToFaceStencilObject>(mesh),
            extendedUpwindCellToFaceStencil(CFCCellToFaceStencil(mesh))
        {
            if (extendedCellToFaceStencil::debug)
            {
                Info<< "Generated pure upwind stencil " << type()
                    << nl << endl;
                writeStencilStats(Info, ownStencil(), ownMap());
            }
        }


    // Destructor

        virtual ~pureUpwindCFCCellToFaceStencilObject()
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
