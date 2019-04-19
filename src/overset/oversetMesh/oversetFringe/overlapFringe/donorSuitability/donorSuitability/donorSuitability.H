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
    Foam::donorSuitability::donorSuitability

Description
    An abstract base class for donor suitability algorithms. This class is
    used by oversetFringe algorithm to optimise overlap given donor suitability
    function.

    Donor suitability criteria is given in isDonorSuitable member function. The
    threshold value represents the admissible relative difference between donor
    and acceptor donor suitability value. For example, if using cellVolumes
    donor suitability, specifying threshold = 0.1 means that the volume of
    smaller cell (either donor or acceptor) needs to be within 10% of the
    other cell's volume.

SourceFiles
    donorSuitability.C
    newAcceptorSet.C

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

Contributor
    Josip Katavic, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef donorSuitability_H
#define donorSuitability_H

#include "dictionary.H"
#include "labelList.H"
#include "donorAcceptorList.H"
#include "scalarField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class oversetFringe;

namespace donorSuitability
{

/*---------------------------------------------------------------------------*\
                         Class donorSuitability Declaration
\*---------------------------------------------------------------------------*/

class donorSuitability
{
    // Private data

        //- Reference to oversetFringe algorithm
        const oversetFringe& oversetFringe_;

        //- Subdictionary for acceptor set specification
        const dictionary coeffDict_;

        //- Donor Suitability Function (DSF), a list of scalar fields for each
        //  processor. In case of serial run, it is a list with a single
        //  entry. Each field is defined in cell centres and is used
        //  for overlap minimisation.
        List<scalarField> dsf_;

        //- Threshold criterion for donor suitability
        scalar threshold_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        donorSuitability(const donorSuitability&);

        //- Disallow default bitwise assignment
        void operator=(const donorSuitability&);


protected:

    // Protected Member Functions

        // Access functions

            //- Return access to the donor suitability function
            List<scalarField>& donorSuitabilityFunction()
            {
                return dsf_;
            }

            //- Return access to the threshold criterion
            scalar& threshold()
            {
                return threshold_;
            }


        // Constructor helper for derived classes

            //- Perform gather scatter of the donor suitability function for
            //  parallel runs
            void combineDonorSuitabilityFunction(const scalarField& localDsf);


public:

    //- Runtime type information
    TypeName("donorSuitability");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            donorSuitability,
            dictionary,
            (
                const oversetFringe& oversetFringeAlgorithm,
                const dictionary& dict
            ),
            (oversetFringeAlgorithm, dict)
        );


    // Constructors

        //- Construct from components
        donorSuitability
        (
            const oversetFringe& oversetFringeAlgorithm,
            const dictionary& dict
        );


    // Selectors

        //- Return an autoPtr to the selected donorSuitability
        static autoPtr<donorSuitability> New
        (
            const oversetFringe& oversetFringeAlgorithm,
            const dictionary& dict
        );


    //- Destructor
    //  Note: the destructor is pure virtual to make this class abstract
    virtual ~donorSuitability() = 0;


    // Member Functions

        // Access

            //- Return a const reference to oversetFringe algorithm
            inline const oversetFringe& oversetFringeAlgorithm() const
            {
                return oversetFringe_;
            }

            //- Return coefficient dictionary
            inline const dictionary& coeffDict() const
            {
                return coeffDict_;
            }

            //- Return threshold
            inline const scalar& threshold() const
            {
                return threshold_;
            }


        // Overlap minimisation via Donor Suitability Function

            //- Donor/acceptor pair suitability fraction.
            //  Subtract smaller DSF from the greater. Divide difference with
            //  greater DSF to get relative difference. Return relative
            //  difference (a value which is always between 0 and 1. Higher
            //  value means better suitability).
            scalar suitabilityFraction
            (
                const donorAcceptor& daPair
            ) const
            {
                // Check whether the donor is valid for this pair
                if (!daPair.donorFound())
                {
                    FatalErrorIn
                    (
                        "scalar donorSuitability::suitabilityFraction"
                        "\n("
                        "\n    const donorAcceptor& daPair"
                        "\n) const"
                    )   << "Acceptor: " << daPair.acceptorCell()
                        << " on processor: " << daPair.acceptorProcNo()
                        << " did not find donor candidate."
                        << nl
                        << "Additional information: "
                        << nl
                        <<  (
                                daPair.withinBB()
                              ? "   Within bounding box"
                              : "   Not within bounding box"
                            )
                        << nl << tab
                        << "Donor index: " << daPair.donorCell()
                        << nl << tab
                        << "Donor processor index: " << daPair.donorProcNo()
                        << nl << tab
                        << "Number of extended donors: "
                        << daPair.extendedDonorCells().size()
                        << nl << nl
                        << "Please review your fringe assembly settings"
                        << " (or try using adaptiveOverlap fringe algorithm)."
                        << abort(FatalError);
                }

                // Get donor suitability function in acceptor/donor pair
                const scalar& dsfAcceptor =
                    dsf_[daPair.acceptorProcNo()][daPair.acceptorCell()];
                const scalar& dsfDonor =
                    dsf_[daPair.donorProcNo()][daPair.donorCell()];

                // Note: Donor Suitability Function should be non-negative, no
                // need to take the absolute value
                if (dsfAcceptor > dsfDonor)
                {
                    return (1 - (dsfAcceptor - dsfDonor)/(dsfAcceptor + SMALL));
                }
                else
                {
                    return (1 - (dsfDonor - dsfAcceptor)/(dsfDonor + SMALL));
                }
            }


            //- Is this donor/acceptor pair suitable?
            //  if sutabilityFraction > threshold
            bool isDonorSuitable
            (
                const donorAcceptor& daPair
            ) const
            {
                return (suitabilityFraction(daPair)) > threshold_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace donorSuitability
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
