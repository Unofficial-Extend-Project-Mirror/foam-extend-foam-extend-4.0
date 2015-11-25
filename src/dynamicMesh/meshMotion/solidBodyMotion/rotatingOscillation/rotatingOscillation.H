#ifndef rotatingOscillation_H
#define rotatingOscillation_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class rotatingOscillation Declaration
\*---------------------------------------------------------------------------*/

class rotatingOscillation
:
    public solidBodyMotionFunction
{
    // Private data

        //- Centre of gravity
        point origin_;

        //- Amplitude
        vector amplitude_;

        //- Motion period
        scalar period_;


    // Private Member Functions

        //- Disallow copy construct
        rotatingOscillation(const rotatingOscillation&);

        //- Disallow default bitwise assignment
        void operator=(const rotatingOscillation&);
	
	//- Calculate position
        vector calcPosition(const scalar t) const;

	//- Calculate tranformation
        septernion calcTransformation(const scalar t) const;


public:

    //- Runtime type information
    TypeName("rotatingOscillation");


    // Constructors

        //- Construct from components
        rotatingOscillation
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

 /*       //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new rotatingOscillation
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }
*/

    //- Destructor
    virtual ~rotatingOscillation();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;
	
	//- Return the solid-body motion velocity
        virtual septernion velocity() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
