/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    Foam::vtkMesh

Description
    Encapsulation of VTK mesh data. Holds mesh or meshsubset and
    polyhedral-cell decomposition on it.

SourceFiles
    vtkMesh.C

\*---------------------------------------------------------------------------*/

#ifndef vtkMesh_H
#define vtkMesh_H

#include "vtkTopo.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;

/*---------------------------------------------------------------------------*\
                           Class vtkMesh Declaration
\*---------------------------------------------------------------------------*/

class vtkMesh
:
    public fvMesh
{
    // Private data

        //- Mesh subset
        fvMeshSubset subsetMesh_;

        //- Current cellSet (or empty)
        const word setName_;

        //- Current decomposition of topology
        mutable autoPtr<vtkTopo> topoPtr_;



    // Private Member Functions

        //- Disallow default bitwise copy construct
        vtkMesh(const vtkMesh&);

        //- Disallow default bitwise assignment
        void operator=(const vtkMesh&);


public:

    // Constructors

        //- Construct from components
        vtkMesh(const IOobject& io, const word& setName = "");


    // Destructor

        virtual ~vtkMesh();


    // Member Functions

        // Access

            const fvMeshSubset& subsetMesh() const
            {
                return subsetMesh_;
            }

            //- Check if running subMesh
            bool useSubMesh() const
            {
                return setName_.size();
            }

            //- topology
            const vtkTopo& topo() const
            {
                if (topoPtr_.empty())
                {
                    topoPtr_.reset(new vtkTopo(mesh()));
                }
                return topoPtr_();
            }

            //- Access either mesh or submesh
            const fvMesh& mesh() const
            {
                if (useSubMesh())
                {
                    return subsetMesh_.subMesh();
                }
                else
                {
                    return *this;
                }
            }

            //- Access either pointMesh of base or pointMesh of subset
            const pointMesh& pMesh() const
            {
                if (useSubMesh())
                {
                    return pointMesh::New(subsetMesh_.subMesh());
                }
                else
                {
                    return pointMesh::New(*this);
                }
            }

            //- Number of field cells
            label nFieldCells() const
            {
                return topo().cellTypes().size();
            }

            //- Number of field points
            label nFieldPoints() const
            {
                return mesh().nPoints() + topo().addPointCellLabels().size();
            }


        // Edit

            //- Read mesh
            polyMesh::readUpdateState readUpdate();


            //- Map volume field (does in fact do very little interpolation;
            //  just copied from fvMeshSubset)
            template<class GeoField>
            tmp<GeoField> interpolate(const GeoField& fld) const
            {
                if (useSubMesh())
                {
                    tmp<GeoField> subFld = subsetMesh_.interpolate(fld);
                    subFld().rename(fld.name());
                    return subFld;
                }
                else
                {
                    return fld;
                }
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
