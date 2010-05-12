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
    Foam::leastSquaresFaVectors

Description
    Least-squares gradient scheme vectors for the Finite Area method

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    leastSquaresFaVectors.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresFaVectors_H
#define leastSquaresFaVectors_H

#include "MeshObject.H"
#include "faMesh.H"
#include "edgeFieldsFwd.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                     Class leastSquaresFaVectors Declaration
\*---------------------------------------------------------------------------*/

class leastSquaresFaVectors
:
    public MeshObject<faMesh, leastSquaresFaVectors>
{
    // Private data

        //- Least-squares gradient vectors
        mutable edgeVectorField* pVectorsPtr_;
        mutable edgeVectorField* nVectorsPtr_;


    // Private member functions

        //- Construct Least-squares gradient vectors
        void makeLeastSquaresVectors() const;


public:

    // Declare name of the class and its debug switch
    TypeName("leastSquaresFaVectors");


    // Constructors

        //- Construct given an faMesh
        explicit leastSquaresFaVectors(const faMesh&);


    // Destructor

        virtual ~leastSquaresFaVectors();


    // Member functions

        //- Return reference to owner least square vectors
        const edgeVectorField& pVectors() const;

        //- Return reference to neighbour least square vectors
        const edgeVectorField& nVectors() const;


        //- Delete the least square vectors when the mesh moves
        virtual bool movePoints() const;

        //- Delete the least square vectors when the mesh updates
        virtual bool updateMesh(const mapPolyMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
