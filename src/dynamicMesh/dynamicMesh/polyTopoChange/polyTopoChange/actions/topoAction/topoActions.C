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

#include "topoAction.H"
#include "polyAddPoint.H"
#include "polyAddFace.H"
#include "polyAddCell.H"
#include "polyModifyPoint.H"
#include "polyModifyFace.H"
#include "polyModifyCell.H"
#include "polyRemovePoint.H"
#include "polyRemoveFace.H"
#include "polyRemoveCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(topoAction, 0);

    defineTypeNameAndDebug(polyAddPoint, 0);
    defineTypeNameAndDebug(polyModifyPoint, 0);
    defineTypeNameAndDebug(polyRemovePoint, 0);

    defineTypeNameAndDebug(polyAddFace, 0);
    defineTypeNameAndDebug(polyModifyFace, 0);
    defineTypeNameAndDebug(polyRemoveFace, 0);

    defineTypeNameAndDebug(polyAddCell, 0);
    defineTypeNameAndDebug(polyModifyCell, 0);
    defineTypeNameAndDebug(polyRemoveCell, 0);
}


// ************************************************************************* //
