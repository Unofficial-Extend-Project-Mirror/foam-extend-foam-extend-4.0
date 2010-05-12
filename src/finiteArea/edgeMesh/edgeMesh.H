/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    edgeMesh

Description
    Mesh data needed to do the Finite Area discretisation.

\*---------------------------------------------------------------------------*/

#ifndef edgeMesh_H
#define edgeMesh_H

#include "GeoMesh.H"
#include "faMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class edgeMesh Declaration
\*---------------------------------------------------------------------------*/

class edgeMesh
:
    public GeoMesh<faMesh>
{

public:

    explicit edgeMesh(const faMesh& mesh)
    :
        GeoMesh<faMesh>(mesh)
    {}

    label size() const
    {
        return size(mesh_);
    }

    static label size(const Mesh& mesh)
    {
        return mesh.nInternalEdges();
    }

    const edgeVectorField& C()
    {
        return mesh_.edgeCentres();
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
