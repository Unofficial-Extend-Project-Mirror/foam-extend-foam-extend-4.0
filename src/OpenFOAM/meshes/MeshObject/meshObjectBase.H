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
    meshObjectBase

Description
    Abstract base-class for dynamic mesh objects used to automate
    their allocation to the mesh database and mesh-related updates

SourceFiles
    meshObjectBase.C

\*---------------------------------------------------------------------------*/

#ifndef meshObjectBase_H
#define meshObjectBase_H

#include "typeInfo.H"
#include "HashTable.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                         Class meshObjectBase Declaration
\*---------------------------------------------------------------------------*/

class meshObjectBase
{
public:

    //- Runtime type information
    TypeName("meshObject");

    // Static functions

        //- Update topology on all mesh objects
        template<class Mesh>
        static void allUpdateTopology(const Mesh& mesh, const mapPolyMesh& mpm)
        {
            HashTable<const meshObjectBase*> tbl =
                mesh.objectRegistry::lookupClass<meshObjectBase>();

            for
            (
                HashTable<const meshObjectBase*>::iterator iter =
                    tbl.begin();
                iter != tbl.end();
                ++iter
            )
            {
                const meshObjectBase& obj = *(iter());

                if (Mesh::debug)
                {
                    Info << "Updating object " << obj.type() << endl;
                }

                // obj.movePoints();
                obj.updateMesh(mpm);
            }

            if (Mesh::debug)
            {
                Info << "Done update topology" << endl;
            }
        }


        //- Move points on all mesh objects
        template<class Mesh>
        static void allMovePoints(const Mesh& mesh)
        {
            if (mesh.moving())
            {
                HashTable<const meshObjectBase*> tbl =
                    mesh.objectRegistry::lookupClass<meshObjectBase>();

                for
                (
                    HashTable<const meshObjectBase*>::iterator iter =
                        tbl.begin();
                    iter != tbl.end();
                    ++iter
                )
                {
                    const meshObjectBase& obj = *(iter());

                    if (Mesh::debug)
                    {
                        Info << "Moving object " << obj.type() << endl;
                    }

                    obj.movePoints();
                }

                if (Mesh::debug)
                {
                    Info << "Done moving" << endl;
                }
            }
        }

        //- Move points on all mesh objects
        template<class Mesh>
        static void allDelete(const Mesh& mesh)
        {
            HashTable<const meshObjectBase*> tbl =
                mesh.objectRegistry::lookupClass<meshObjectBase>();

            for
            (
                HashTable<const meshObjectBase*>::iterator iter =
                    tbl.begin();
                iter != tbl.end();
                ++iter
            )
            {
                const meshObjectBase& obj = *(iter());

                if (Mesh::debug)
                {
                    Info << "Deleting object " << obj.type() << endl;

                    obj.deleteObject();
                }
            }
        }


    // Destructor

        virtual ~meshObjectBase()
        {}


    // Member Functions

        //- Update after mesh motion
        virtual bool movePoints() const = 0;

        //- Update after topology change
        virtual bool updateMesh(const mapPolyMesh&) const = 0;

        //- Delete object from database
        virtual bool deleteObject() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
