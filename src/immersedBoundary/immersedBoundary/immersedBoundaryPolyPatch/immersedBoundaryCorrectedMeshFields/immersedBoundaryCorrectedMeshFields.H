/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::immersedBoundaryCorrectedMeshFields

Description
    MeshObjet derived class for global corrected data in immersed boundary poly
    patch.

SourceFiles
    immersedBoundaryCorrectedMeshFields.C

Author
    Inno Gatin, FMENA, Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryCorrectedMeshFields_H
#define immersedBoundaryCorrectedMeshFields_H

#include "MeshObject.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class immersedBoundaryCorrectedMeshFields Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryCorrectedMeshFields
:
    public MeshObject<polyMesh, immersedBoundaryCorrectedMeshFields>
{
    // Private data

        //- Corrected cell centres
        mutable vectorField* correctedCellCentresPtr_;

        //- Corrected face centres
        mutable vectorField* correctedFaceCentresPtr_;

        //- Corrected cell volumes
        mutable scalarField* correctedCellVolumesPtr_;

        //- Corrected face areas
        mutable vectorField* correctedFaceAreasPtr_;

        //- Current topo change index
        mutable label curTopoIndex_;


public:

    // Declare name of the class and its debug switch
    TypeName("immersedBoundaryCorrectedMeshFields");


    // Constructors

        //- Construct given an polyMesh
        explicit immersedBoundaryCorrectedMeshFields(const polyMesh&);

    // Destructor
        virtual ~immersedBoundaryCorrectedMeshFields();

    //- Clear data
        void clearOut(const label topoChangeIndex) const;

    // Access

        // Corrected mesh geometry

            //- Corrected cell centres
            vectorField& correctedCellCentres() const;

            //- Corrected face centres
            vectorField& correctedFaceCentres() const;

            //- Corrected cell volumes
            scalarField& correctedCellVolumes() const;

            //- Corrected face areas
            vectorField& correctedFaceAreas() const;


    // Virtual dummy functions: immersed boundary motion is handled elswere

    virtual bool movePoints() const
    {
        return true;
    }

    virtual bool updateMesh(const mapPolyMesh&) const
    {
        return true;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
