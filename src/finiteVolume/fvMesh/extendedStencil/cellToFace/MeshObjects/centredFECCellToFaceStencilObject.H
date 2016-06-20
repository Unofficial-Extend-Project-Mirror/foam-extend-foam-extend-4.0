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
    Foam::centredFECCellToFaceStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef centredFECCellToFaceStencilObject_H
#define centredFECCellToFaceStencilObject_H

#include "extendedCentredCellToFaceStencil.H"
#include "FECCellToFaceStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class centredFECCellToFaceStencilObject Declaration
\*---------------------------------------------------------------------------*/

class centredFECCellToFaceStencilObject
:
    public MeshObject<fvMesh, centredFECCellToFaceStencilObject>,
    public extendedCentredCellToFaceStencil
{

public:

    TypeName("centredCFCCellToFaceStencil");

    // Constructors

        //- Construct from uncompacted face stencil
        explicit centredFECCellToFaceStencilObject
        (
            const fvMesh& mesh
        )
        :
            MeshObject<fvMesh, centredFECCellToFaceStencilObject>(mesh),
            extendedCentredCellToFaceStencil(FECCellToFaceStencil(mesh))
        {
            if (extendedCellToFaceStencil::debug)
            {
                Info<< "Generated centred stencil " << type()
                    << nl << endl;
                writeStencilStats(Info, stencil(), map());
            }
        }


    // Destructor

        virtual ~centredFECCellToFaceStencilObject()
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
