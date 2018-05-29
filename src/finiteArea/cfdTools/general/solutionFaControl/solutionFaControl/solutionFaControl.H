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
    Foam::solutionFaControl

Description
    Base class for solution control classes

\*---------------------------------------------------------------------------*/

#ifndef solutionFaControl_H
#define solutionFaControl_H

#include "faMesh.H"
#include "edgeMesh.H"
#include "faePatchField.H"
#include "faMatrices.H"
#include "wordRe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solutionFaControl Declaration
\*---------------------------------------------------------------------------*/

class solutionFaControl
:
    public IOobject
{
public:

    struct fieldData
    {
        wordRe name;
        scalar absTol;
        scalar relTol;
        scalar initialResidual;
    };


protected:

    // Protected data

        //- Reference to the mesh database
        faMesh& mesh_;

        //- List of residual data per field
        List<fieldData> residualControl_;

        //- The dictionary name, e.g. SIMPLE, PIMPLE
        const word algorithmName_;


        // Solution controls

            //- Maximum number of non-orthogonal correctors
            label nNonOrthCorr_;

            //- Flag to indicate to solve for momentum
            bool momentumPredictor_;

            //- Flag to indicate to solve using transonic algorithm
            bool transonic_;

            //- Flag to indicate to relax pressure using the
            //  "consistent" approach of SIMPLEC
            bool consistent_;


        // Eareaution

            //- Current corrector loop index
            label corr_;

            //- Current non-orthogonal corrector loop index
            label corrNonOrtho_;


    // Protected Member Functions

        //- Read controls from faSolution dictionary
        virtual void read(const bool absTolOnly);

        //- Read controls from faSolution dictionary
        virtual void read();

        //- Return index of field in residualControl_ if present
        virtual label applyToField
        (
            const word& fieldName,
            const bool useRegEx = true
        ) const;

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied() = 0;

        //- Store previous iteration fields
        virtual void storePrevIterFields() const;

        //- Store previous iteration field for area<Type>Fields
        template<class Type>
        void storePrevIter() const;

        template<class Type>
        void maxTypeResidual
        (
            const word& fieldName,
            ITstream& data,
            scalar& firstRes,
            scalar& lastRes
        ) const;

        scalar maxResidual
        (
            const word& fieldName,
            ITstream& data,
            scalar& lastRes
        ) const;


private:

        //- Disallow default bitwise copy construct
        solutionFaControl(const solutionFaControl&);

        //- Disallow default bitwise assignment
        void operator=(const solutionFaControl&);


public:


    // Static Data Members

        //- Run-time type information
        TypeName("solutionFaControl");


    // Constructors

        //- Construct from mesh
        solutionFaControl(faMesh& mesh, const word& algorithmName);


    //- Destructor
    virtual ~solutionFaControl();


    // Member Functions

        // Access

            //- Return the solution dictionary
            inline const dictionary& dict() const;

            //- Current corrector loop index
            inline label corr() const;

            //- Current non-orthogonal corrector index
            inline label corrNonOrtho() const;


        // Solution control

            //- Maximum number of non-orthogonal correctors
            inline label nNonOrthCorr() const;

            //- Helper function to identify final non-orthogonal iteration
            inline bool finalNonOrthogonalIter() const;

            //- Flag to indicate to solve for momentum
            inline bool momentumPredictor() const;

            //- Flag to indicate to solve using transonic algorithm
            inline bool transonic() const;

            //- Flag to indicate to relax pressure using the
            //  "consistent" approach of SIMPLEC
            inline bool consistent() const;


        // Eareaution

            //- Main control loop
            virtual bool loop() = 0;

            //- Non-orthogonal corrector loop
            inline bool correctNonOrthogonal();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solutionFaControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "solutionFaControlTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
