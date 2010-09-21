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
    tetMetrics

Description
    Header file of tetrahedral mesh-quality metrics.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    tetMetrics.C

\*---------------------------------------------------------------------------*/

#ifndef tetMetrics_H
#define tetMetrics_H

#include "tetMetric.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Metric Class Declarations
\*---------------------------------------------------------------------------*/

class Knupp
:
    public tetMetric
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        Knupp(const Knupp&);

        //- Disallow default bitwise assignment
        void operator=(const Knupp&);


public:

        //- Runtime type information
        TypeName("Knupp");

        static scalar metric
        (
            const point& p0,
            const point& p1,
            const point& p2,
            const point& p3
        );


    // Destructor

        virtual ~Knupp()
        {}
};


class Dihedral
:
    public tetMetric
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        Dihedral(const Dihedral&);

        //- Disallow default bitwise assignment
        void operator=(const Dihedral&);

        // Enumeration for tets
        static label tetEnum[6][4];


public:

        //- Runtime type information
        TypeName("Dihedral");

        static scalar metric
        (
            const point& p0,
            const point& p1,
            const point& p2,
            const point& p3
        );


    // Destructor

        virtual ~Dihedral()
        {}
};


class cubicMeanRatio
:
    public tetMetric
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        cubicMeanRatio(const cubicMeanRatio&);

        //- Disallow default bitwise assignment
        void operator=(const cubicMeanRatio&);

    public:

        //- Runtime type information
        TypeName("cubicMeanRatio");

        static scalar metric
        (
            const point& p0,
            const point& p1,
            const point& p2,
            const point& p3
        );


    // Destructor

        virtual ~cubicMeanRatio()
        {}
};


class Frobenius
:
    public tetMetric
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        Frobenius(const Frobenius&);

        //- Disallow default bitwise assignment
        void operator=(const Frobenius&);


public:

        //- Runtime type information
        TypeName("Frobenius");

        static scalar metric
        (
            const point& p0,
            const point& p1,
            const point& p2,
            const point& p3
        );


    // Destructor

        virtual ~Frobenius()
        {}
};


class PGH
:
    public tetMetric
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        PGH(const PGH&);

        //- Disallow default bitwise assignment
        void operator=(const PGH&);


public:

        //- Runtime type information
        TypeName("PGH");

        static scalar metric
        (
            const point& p0,
            const point& p1,
            const point& p2,
            const point& p3
        );


    // Destructor

        virtual ~PGH()
        {}
};


class CSG
:
    public tetMetric
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        CSG(const CSG&);

        //- Disallow default bitwise assignment
        void operator=(const CSG&);


public:

        //- Runtime type information
        TypeName("CSG");

        static scalar metric
        (
            const point& p0,
            const point& p1,
            const point& p2,
            const point& p3
        );


    // Destructor

        virtual ~CSG()
        {}
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
