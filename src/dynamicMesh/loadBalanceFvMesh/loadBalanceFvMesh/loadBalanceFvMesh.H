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
    Foam::loadBalanceFvMesh

Description
    Load balancing fvMesh

SourceFiles
    loadBalanceFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef loadBalanceFvMesh_H
#define loadBalanceFvMesh_H

#include "dynamicFvMesh.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class loadBalanceFvMesh Declaration
\*---------------------------------------------------------------------------*/

class loadBalanceFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Dictionary of motion control parameters
        dictionary dynamicMeshCoeffs_;

        //- Imbalance trigger 0 < imbalance < 1
        scalar imbalanceTrigger_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        loadBalanceFvMesh(const loadBalanceFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const loadBalanceFvMesh&);

        //- Check balance.  Return true if re-balancing is required
        bool checkLoadBalance(const scalarField& weights) const;

        //- Send fields to other processor
        template<class GeoField, class Decomposer>
        void sendFields
        (
            const HashTable<const GeoField*>& geoFields,
            const Decomposer& decomposer,
            Ostream& toProc
        ) const;

        //- Insert local fields
        template<class GeoField, class Decomposer>
        void insertFields
        (
            const HashTable<const GeoField*>& geoFields,
            const Decomposer& decomposer,
            List<PtrList<GeoField> >& localFields
        ) const;

        //- Receive fields from other processor
        template<class GeoMesh, class GeoField>
        void receiveFields
        (
            const label procIndex,
            List<PtrList<GeoField> >& receivedFields,
            const GeoMesh& procMesh,
            Istream& fromProc
        ) const;

        //- Resize and rebuild fields
        template<class GeoField, class Reconstructor>
        void rebuildFields
        (
            const HashTable<const GeoField*>& geoFields,
            const Reconstructor& reconstructor,
            const List<PtrList<GeoField> >& receivedFields,
            const boolList& patchesToReplace
        ) const;


public:

    //- Runtime type information
    TypeName("loadBalanceFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit loadBalanceFvMesh(const IOobject& io);


    //- Destructor
    virtual ~loadBalanceFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "loadBalanceFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
