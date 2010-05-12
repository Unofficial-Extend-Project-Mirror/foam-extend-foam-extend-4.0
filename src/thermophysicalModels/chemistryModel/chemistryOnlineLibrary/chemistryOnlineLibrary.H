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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description

    Implementation of the ISAT (In-situ adaptive tabulation), for chemistry 
    calculation. For forward information:
    
    S.B. Pope (1997) 
    ``Computationally efficient implementation of combustion chemistry using 
    in situ adaptive tabulation,'' 
    Combustion Theory and Modelling,  1, 41--63.

\*---------------------------------------------------------------------------*/

#ifndef chemistryOnlineLibrary_H
#define chemistryOnlineLibrary_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IOstream.H"
#include "dictionary.H"
#include "Switch.H"
#include "Matrix.H"
#include "scalarField.H"
#include "binaryTree.H" 

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class chemistryOnlineLibrary Declaration
\*---------------------------------------------------------------------------*/
class chemistryModel;

class chemistryOnlineLibrary
{
    // Private data
   
        //- Dictionary for data lookup
        const dictionary& onlineDict_;
        
        //- Reference to the chemistryModel
        const chemistryModel& chemistry_;

        //- Generate online lib (on/off)
        Switch online_;

        //- List of the stored 'points' organized in a binary tree
        binaryTree chemisTree_;
    
        //- Tolerance for species and temperature
        scalar tolerance_;

        //- List of scale factors for species, temperature and pressure
        scalarField scaleFactor_;
      
        //- List of scale factors for species, temperature and pressure
        scalarField solutionScaleFactor_;
      
        //- how to deal with temperature
        Switch logT_;
        
        //- Clean or not to Clean?
        Switch clean_;
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        chemistryOnlineLibrary(const chemistryOnlineLibrary&);

        //- Disallow default bitwise assignment
        void operator=(const chemistryOnlineLibrary&);

        


public:

    // Constructors

        //- Construct from dictionary
        chemistryOnlineLibrary
        (
            const chemistryModel& chemistry,
            const dictionary& chemistryProperties
        );


    // Selectors

    // Destructor

        ~chemistryOnlineLibrary();


    // Member Functions

        // Access

        inline binaryTree& chemisTree() 
        {
            return chemisTree_;
        }

        inline const Switch online() const
        {
            return online_;
        }
        inline const scalarField& scaleFactor() const
        {
            return scaleFactor_;
        }

        inline const Switch clean() const
        {
            return clean_;
        }
        
        // Database

        //- Find the closest stored composition in the binary tree
        chemPoint* findClosest(const scalarField& v0);
        
        //- Return the size of the binary tree
        inline label size()
        {
            return chemisTree_.size();
        }
        
        // Edit

        //- Add a new point to the binary tree
        void add
        (
            const scalarField& v0,    
            const scalarField& r, 
            const scalarField& Wi,
            const scalar rhoi,
            const scalar& Tr,
            const scalar& p,
            const scalar& deltaT
        );

        //- Calculate the new concentration on the basis of a stored data        
        const scalarField calcNewC
        (
            const chemPoint& c0
        );

        //- Clear the library
        void clear();

        //- Scan for debugging purposes
        void scan()
        {
            chemisTree_.scan();
        }
        
        // IOstream Operators

        friend Istream& operator>>(Istream&, chemistryOnlineLibrary&);
        friend Ostream& operator<<(Ostream&, const chemistryOnlineLibrary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "chemistryOnlineLibraryI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
