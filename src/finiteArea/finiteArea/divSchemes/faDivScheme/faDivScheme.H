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
    divScheme

Description
    Abstract base class for div schemes.

SourceFiles
    faDivScheme.C
    faDivSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef faDivScheme_H
#define faDivScheme_H

#include "tmp.H"
#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "linearEdgeInterpolation.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class faMatrix;

class faMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                           Class divScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class divScheme
:
    public refCount
{

protected:

    // Protected data

        const faMesh& mesh_;
        tmp<edgeInterpolationScheme<Type> > tinterpScheme_;


    // Private Member Functions

        //- Disallow copy construct
        divScheme(const divScheme&);

        //- Disallow default bitwise assignment
        void operator=(const divScheme&);


public:

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            divScheme,
            Istream,
            (const faMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh
        divScheme(const faMesh& mesh)
        :
            mesh_(mesh),
            tinterpScheme_(new linearEdgeInterpolation<Type>(mesh))
        {}

        //- Construct from mesh and Istream
        divScheme(const faMesh& mesh, Istream& is)
        :
            mesh_(mesh),
            tinterpScheme_(edgeInterpolationScheme<Type>::New(mesh, is))
        {}


    // Selectors

        //- Return a pointer to a new divScheme created on freestore
        static tmp<divScheme<Type> > New
        (
            const faMesh& mesh,
            Istream& schemeData
        );


    // Destructor

        virtual ~divScheme();


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const
        {
            return mesh_;
        }

        virtual tmp
        <
            GeometricField
            <typename innerProduct<vector, Type>::type, faPatchField, areaMesh>
        > facDiv
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFaDivTypeScheme(SS, Type)                                          \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
divScheme<Type>::addIstreamConstructorToTable<SS<Type> >                       \
    add##SS##Type##IstreamConstructorToTable_;


#define makeFaDivScheme(SS)                                                    \
                                                                               \
makeFaDivTypeScheme(SS, vector)                                                \
makeFaDivTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "faDivScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
