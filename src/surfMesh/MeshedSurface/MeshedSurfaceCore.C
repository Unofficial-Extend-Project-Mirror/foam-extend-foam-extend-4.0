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

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "foamTime.H"
#include "MeshedSurface.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{

    // specialization: from face -> triFace
    template<>
    void Foam::MeshedSurface<triFace>::transcribe(MeshedSurface<face>& surf)
    {
        // first triangulate
        surf.triangulate();
        this->storedPoints().transfer(surf.storedPoints());
        this->storedZones().transfer(surf.storedZones());

        // transcribe from face -> triFace
        List<face>&    origFaces = surf.storedFaces();
        List<triFace>  newFaces(origFaces.size());
        forAll(origFaces, faceI)
        {
            newFaces[faceI] = triFace
            (
                static_cast<const UList<label>&>(origFaces[faceI])
            );
        }
        surf.clear();

        this->storedFaces().transfer(newFaces);
    }


    // specialization: from face -> face
    template<>
    void Foam::MeshedSurface<face>::transcribe(MeshedSurface<face>& surf)
    {
        this->transfer(surf);
    }


}  // end of namespace Foam


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
