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
    Foam::writer

Description
    Base class for graphics format writing. Entry points are
    - write(..). \n
      Write to an Ostream a table of points with corresponding values.
    - write(scalar/vector/sphericalTensor/symmTensor/tensor). \n
      Write single scalar/vector/sphericalTensor/symmTensor/tensor.
      Default is to write space separated components.

    Example:
    @verbatim
        // Construct writer of xmgr type
        autoPtr<writer<scalar> > scalarFormatter(writer<scalar>::New("xmgr"));

        // Output list of points and corresponding values
        scalarFormatter().write
        (
            coordSet
            (
                points,         // sample coordinates
                "someLine",     // name of coordSet
                "distance",     // write coordinates as distance to refPoint
                points[0]       // reference point
            ),
            "U.component(0)",   // name of values
            vals                // values
        );
    @endverbatim

SourceFiles
    writer.C

\*---------------------------------------------------------------------------*/

#ifndef writer_H
#define writer_H

#include "fileName.H"
#include "wordList.H"
#include "vector.H"
#include "sphericalTensor.H"
#include "diagTensor.H"
#include "symmTensor.H"
#include "tensor.H"
#include "symmTensor4thOrder.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class coordSet;

/*---------------------------------------------------------------------------*\
                           Class writer Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class writer
{

protected:

    //- Generates filename from coordSet and sampled fields
    fileName getBaseName(const coordSet&, const wordList&) const;

    void writeCoord(const coordSet&, const label sampleI, Ostream&) const;

    //- Writes single-column ascii write. Column 1 is coordSet coordinate,
    //  columns 2 is the value. Uses write() function
    //  to write coordinate in correct format.
    void writeTable(const coordSet&, const List<Type>&, Ostream&) const;

    //- Writes multi-column ascii write. Column 1 is coordSet coordinate,
    //  columns 2..n are the values. Uses write() function
    //  to write coordinate in correct format.
    void writeTable
    (
        const coordSet&,
        const List<const List<Type>*>&,
        Ostream& os
    ) const;


public:

    //- Runtime type information
    TypeName("writer");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            writer,
            word,
            (),
            ()
        );


    // Selectors

        //- Return a reference to the selected writer
        static autoPtr<writer> New(const word& writeFormat);


    // Constructors

        //- Construct null
        writer();


    //- Destructor
    virtual ~writer() = 0;


    // Member Functions

        //- Generate file name with correct extension
        virtual fileName getFileName
        (
            const coordSet&,
            const wordList&
        ) const = 0;

        //- General entry point for writing.
        //  The data is organized in a set of point with one or more values
        //  per point
        virtual void write
        (
            const coordSet&,
            const wordList&,
            const List<const Field<Type>*>&,
            Ostream&
        ) const = 0;

        //- General entry point for writing of multiple coordSets.
        //  Each coordSet (track) has same data variables.
        //  The data is per variable, per track, per point of track.
        //  If writeTracks adds connecting lines (wherever applicable)
        virtual void write
        (
            const bool writeTracks,
            const PtrList<coordSet>&,
            const wordList& valueSetNames,
            const List<List<Field<Type> > >&,
            Ostream&
        ) const = 0;

        //- Write scalar as ascii
        virtual Ostream& write(const scalar, Ostream&) const;

        template<class VSType>
        Ostream& writeVS(const VSType&, Ostream&) const;

        //- Write vector. Tab separated ascii
        virtual Ostream& write(const vector&, Ostream&) const;

        //- Write sphericalTensor. Tab separated ascii
        virtual Ostream& write(const sphericalTensor&, Ostream&) const;

        //- Write diagTensor. Tab separated ascii
        virtual Ostream& write(const diagTensor&, Ostream&) const;

        //- Write symmTensor. Tab separated ascii
        virtual Ostream& write(const symmTensor&, Ostream&) const;

        //- Write tensor. Tab separated ascii
        virtual Ostream& write(const tensor&, Ostream&) const;

        //- Write symmTensor4thOrder. Tab separated ascii
        virtual Ostream& write(const symmTensor4thOrder&, Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "writer.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
