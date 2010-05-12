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
    Foam::fieldAverage

Description
    Calculates the field averages given list of fieldAverageItems, e.g.

    @verbatim
    fieldAverage1
    {
        // Type of functionObject
        type fieldAverage;

        // Where to load it from (if not already in solver)
        functionObjectLibs ("libfieldAverage.so");

        // Fields to be probed. runTime modifiable!
        fields
        (
            U
            {
                mean            on;
                prime2Mean      on;
                base            time;
            }
            p
            {
                mean            on;
                prime2Mean      on;
                base            time;
            }
        );
    @endverbatim

    Member function calcAverages() calculates the averages.

    Member function fieldAverage::write() calls calcAverages(). Average
    field names are constructed by concatenating the base field with the
    averaging type, e.g.
    - base field, U
    - arithmetic mean field, UMean
    - prime-squared field, UPrime2Mean

    Information regarding the number of averaging steps, and total averaging
    time are written on a (base) per-field basis to the
    fieldAveragingProperties dictionary, located in \<time\>/uniform

SourceFiles
    fieldAverage.C
    fieldAverageTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fieldAverage_H
#define fieldAverage_H

#include "volFieldsFwd.H"
#include "pointFieldFwd.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class fieldAverageItem;
class OFstream;
template<class Type>
class List;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                         Class fieldAverage Declaration
\*---------------------------------------------------------------------------*/

class fieldAverage
{
protected:

    // Private data

        //- Name of this set of field averages.
        word name_;

        //- Database this class is registered to
        const objectRegistry& obr_;

        //- On/off switch
        bool active_;

        //- Time at last call, prevents repeated averaging
        label prevTimeIndex_;

        //- resetOnOutput flag
        Switch resetOnOutput_;

        //- List of field average items, describing what averages to be
        //  calculated and output
        List<fieldAverageItem> faItems_;

        // File and field name extensions

            //- Mean average
            static const word EXT_MEAN;

            //- Prime-squared average
            static const word EXT_PRIME2MEAN;


        // Lists of averages

            // Arithmetic mean fields
            wordList meanScalarFields_;
            wordList meanVectorFields_;
            wordList meanSphericalTensorFields_;
            wordList meanSymmTensorFields_;
            wordList meanTensorFields_;

            // Prime-squared fields - applicable to volVectorFields only
            wordList prime2MeanScalarFields_;
            wordList prime2MeanSymmTensorFields_;


        // Counters

            //- Iteration steps counter
            List<label> totalIter_;

            //- Total time counter
            List<scalar> totalTime_;


    // Private Member Functions

        // Initialisation routines

            //- Checkout fields (causes deletion) from the database
            void checkoutFields(const wordList&) const;

            //- Reset size of lists (clear existing values)
            void resetLists(const label nItems);

            //- Intitialise averaging. Check requested field averages are
            //  valid, and populate field lists
            void initialise();

            //- Add mean average field to PtrList
            template<class Type>
            void addMeanField(const label, wordList&) const;

            //- Add prime-squared average field to PtrList
            template<class Type1, class Type2>
            void addPrime2MeanField
            (
                const label,
                const wordList&,
                wordList&
            ) const;


        // Calculation functions

            //- Main calculation routine
            virtual void calcAverages();

            //- Calculate mean average fields
            template<class Type>
            void calculateMeanFields(const wordList&) const;

            //- Add mean-squared field value to prime-squared mean field
            template<class Type1, class Type2>
            void addMeanSqrToPrime2Mean
            (
                const wordList&,
                const wordList&
            ) const;

            //- Calculate prime-squared average fields
            template<class Type1, class Type2>
            void calculatePrime2MeanFields
            (
                const wordList&,
                const wordList&
            ) const;


        // I-O

            //- Write averages
            virtual void writeAverages() const;

            //- Write fields
            template<class Type>
            void writeFieldList(const wordList&) const;

            //- Write averaging properties - steps and time
            void writeAveragingProperties() const;

            //- Read averaging properties - steps and time
            void readAveragingProperties();


        // Functions to be over-ridden from IOoutputFilter class

            //- Update mesh
            virtual void updateMesh(const mapPolyMesh&);

            //- Move points
            virtual void movePoints(const Field<point>&);


        //- Disallow default bitwise copy construct
        fieldAverage(const fieldAverage&);

        //- Disallow default bitwise assignment
        void operator=(const fieldAverage&);


public:

    //- Runtime type information
    TypeName("fieldAverage");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        fieldAverage
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor

        virtual ~fieldAverage();


    // Member Functions

        //- Return name of the set of field averages
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the field average data
        virtual void read(const dictionary&);

        //- Calculate the field average data and write
        virtual void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fieldAverageTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
