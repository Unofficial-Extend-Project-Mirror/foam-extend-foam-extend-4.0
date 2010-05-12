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
    coupleManager

Description
    Object to handle the coupling of region patches. It can be queried to
    return the neighbour information.

SourceFiles
    coupleManager.C

\*---------------------------------------------------------------------------*/

#ifndef coupleManager_H
#define coupleManager_H

#include "Ostream.H"
#include "dictionary.H"
#include "fvPatch.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coupleManager Declaration
\*---------------------------------------------------------------------------*/

class coupleManager
{
    // Private data

        //- Reference to the local fvPatch
        const fvPatch& patch_;

        //- Name of neighbour region
        word neighbourRegionName_;

        //- Name of patch on the neighbour region
        word neighbourPatchName_;

        //- Name of field on the neighbour region
        word neighbourFieldName_;

        //- Reference to the local region
        const fvMesh& localRegion_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const coupleManager&);


public:

    // Constructors

        //- Construct from fvPatch
        coupleManager(const fvPatch& patch);

        //- Construct from fvPatch and dictionary
        coupleManager(const fvPatch& patch, const dictionary& dict);

        //- Copy constructor
        coupleManager(const coupleManager& cm);


    // Destructor

        ~coupleManager();


    // Member Functions

        // Access

            //- Return a reference to the local patch
            inline const fvPatch& patch() const;

            //- Return the name of the neighbour region
            inline const word& neighbourRegionName() const;

            //- Return the name of the patch on the neighbour region
            inline const word& neighbourPatchName() const;

            //- Return the name of the field on the neighbour region
            inline const word& neighbourFieldName() const;

            //- Return a reference to the neighbour mesh
            inline const fvMesh& neighbourRegion() const;

            //- Return the neighbour patch ID
            inline label neighbourPatchID() const;

            //- Return a reference to the neighbour patch
            inline const fvPatch& neighbourPatch() const;

            //- Return a reference to the neighbour patch field
            template<class Type>
            inline const fvPatchField<Type>& neighbourPatchField() const;

            //- Am I owner (= first to evaluate) of this region interface?
            bool regionOwner() const;

            //- Check that the couple is valid
            void checkCouple() const;


        // Edit

            //- Return the name of the neighbour region
            word& neighbourRegionName();

            //- Return the name of the patch on the neighbour region
            word& neighbourPatchName();

            //- Return the name of the field on the neighbour region
            word& neighbourFieldName();


        // Write

           //- Write couple to obj file
           void coupleToObj() const;

           //- Write entries for patches
           void writeEntries(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "coupleManagerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
