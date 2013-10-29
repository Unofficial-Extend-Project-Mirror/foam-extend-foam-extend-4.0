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
    Foam::SymmTensor4thOrder

Description
    Templated 3D symmetric fourth order tensor derived from VectorSpace
    adding construction from 9 components,
    element access using xxxx(), xxyy() etc. member functions
    and double dot product with a second order tensor.

    It is assumed that the 4th order tensor takes the form (in 6x6 matrix form):
    / xxxx xxyy xxzz                \
    |      yyyy yyzz                |
    |           zzzz                |
    |                xyxy           |
    |                     yzyz      |
    \                          zxzx /
    where the tensor is symmetric and unspecified entries are zero.

SourceFiles
    SymmTensor4thOrderI.H

Author
    Philip Cardiff UCD
\*---------------------------------------------------------------------------*/

#ifndef SymmTensor4thOrder_H
#define SymmTensor4thOrder_H

#include "VectorSpace.H"
#include "SymmTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SymmTensor4thOrder Declaration
\*---------------------------------------------------------------------------*/

template <class Cmpt>
class SymmTensor4thOrder
:
    public VectorSpace<SymmTensor4thOrder<Cmpt>, Cmpt, 9>
{

public:

    //- Equivalent type of labels used for valid component indexing
    typedef SymmTensor4thOrder<label> labelType;


    // Member constants

        enum
        {
            rank = 2 // check: Rank of SymmTensor4thOrder is 4
        };


    // Static data members

        static const char* const typeName;
        static const char* componentNames[];

        static const SymmTensor4thOrder zero;
        static const SymmTensor4thOrder one;
        static const SymmTensor4thOrder max;
        static const SymmTensor4thOrder min;


    //- Component labeling enumeration
        enum components { XXXX, XXYY, XXZZ, YYYY, YYZZ, ZZZZ, XYXY, YZYZ, ZXZX };


    // Constructors

        //- Construct null
        inline SymmTensor4thOrder();

        //- Construct given VectorSpace
        inline SymmTensor4thOrder(const VectorSpace<SymmTensor4thOrder<Cmpt>, Cmpt, 9>&);

        //- Construct given SymmTensor (for compilation)
        inline SymmTensor4thOrder(const SymmTensor<Cmpt>&);

        //- Construct given the nine components
        inline SymmTensor4thOrder
        (
            const Cmpt txxxx, const Cmpt txxyy, const Cmpt txxzz,
	                      const Cmpt tyyyy, const Cmpt tyyzz,
	                                        const Cmpt tzzzz,
	                                                   	  const Cmpt txyxy,
	                                        	                            const Cmpt tyzyz,
	                                        	                                              const Cmpt tzxzx
        );

        //- Construct from Istream
        SymmTensor4thOrder(Istream&);


    // Member Functions

        // Access

            inline const Cmpt& xxxx() const;
            inline const Cmpt& xxyy() const;
            inline const Cmpt& xxzz() const;
            inline const Cmpt& yyyy() const;
            inline const Cmpt& yyzz() const;
            inline const Cmpt& zzzz() const;
            inline const Cmpt& xyxy() const;
            inline const Cmpt& yzyz() const;
            inline const Cmpt& zxzx() const;

            inline Cmpt& xxxx();
            inline Cmpt& xxyy();
            inline Cmpt& xxzz();
            inline Cmpt& yyyy();
            inline Cmpt& yyzz();
            inline Cmpt& zzzz();
            inline Cmpt& xyxy();
            inline Cmpt& yzyz();
            inline Cmpt& zxzx();


        //- Assign to a SymmTensor
        //- maybe I need this to compile
        inline void operator=(const SymmTensor<Cmpt>&);
};


// template<class Cmpt>
// class typeOfRank<Cmpt, 4>
// {
// public:

//     typedef SymmTensor4thOrder<Cmpt> type;
// };


// template<class Cmpt>
// class symmTypeOfRank<Cmpt, 4>
// {
// public:

//     typedef SymmTensor4thOrder<Cmpt> type;
// };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "SymmTensor4thOrderI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
