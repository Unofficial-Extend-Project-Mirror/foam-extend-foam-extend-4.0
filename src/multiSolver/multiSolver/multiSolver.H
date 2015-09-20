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
    Foam::multiSolver

Description
    Manages multiple solvers within one master-solver.  Allows for multiple
    definitions of time.  Works by changing the case directory as required by
    each solver.

SourceFiles
    multiSolverI.H
    multiSolver.C
    multiSolverSetControls.C
    multiSolverSetInitialSolverDomain.C
    multiSolverSetNextSolverDomain.C
    multiSolverTimeFunctions.C

Author
    David L. F. Gaden
\*---------------------------------------------------------------------------*/

#ifndef multiSolver_H
#define multiSolver_H

#include "multiSolverVersion.H"
#include "multiTime.H"
#include "dummyControlDict.H"
#include "timeClusterList.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multiSolver Declaration
\*---------------------------------------------------------------------------*/

class multiSolver
{

    // Private data

        // A dummy controlDict, required to initialize multiDictRegistry
        dummyControlDict dcd_;

        // A mini-objectRegistry for the dictionaries - done this way to allow
        // run-time modification of dictionaries
        multiTime multiDictRegistry_;

        // Main control dictionary for multiSolver
        IOdictionary multiControlDict_;

        // All multiDicts contained in the case directory
        PtrList<IOdictionary> multiDicts_;

        // Reference to multiSolverControl subdictionary in multiControlDict
        dictionary& multiSolverControl_;

        // Reference to solverDomains subdictionary in multiControlDict
        dictionary& solverDomains_;

        // Current solverDomain dictionary
        dictionary currentSolverDomainDict_;

        // List of all the solver domain prefixes
        wordList prefixes_;

        // True when setInitialSolverDomain has been run
        bool initialized_;

        // Set by the solver to force a global end
        bool forcedEnd_;

        // True when superLoop++ just happened, but the previous solverDomain
        // has not yet been saved (to prevent the first solverDomain from
        // saving to the *next* superLoop.
        bool noSaveSinceSuperLoopIncrement_;

        // True if this is the lowest number on a local shared drive - parallel
        // only.
        bool manageLocalRoot_;


    // Private member functions

        // Set manageLocalRoot_ flags on all processors (parallel only)
        void setUpParallel();

        // Ensure all processors are synchronized (parallel only)
        void synchronizeParallel() const;

        // Load initial data and controls
        // *** Located in multiSolverSetInitialSolverDomain.C ***
        void setInitialSolverDomain(const word& solverDomainName);

        // Load data and controls of a subsequent solver domain
        // *** Located in multiSolverSetNextSolverDomain.C ***
        void setNextSolverDomain(const word& solverDomainName);

        // Sets controls from the multiSolverControl subdictionary in the
        // multiControlDict
        // *** Located in multiSolverSetControls.C ***
        void setMultiSolverControls();

        // Sets controls from the solverDomains subdictionary in the
        // multiControlDict
        // *** Located in multiSolverSetControls.C ***
        void setSolverDomainControls(const word& solverDomainName);

        // Use stopAt and finalStopAt settings to determine the endTime for
        // the current solverDomain.  Returns the controlDict's stopAt word,
        // which may differ from that in multiControlDict.
        word setLocalEndTime();

        // If the current solver domain should not start, returns true.
        bool checkGlobalEnd() const;

        // Looks for case/multiSolver/solverDomain/initial/0.  Fail if missing.
        void checkTimeDirectories() const;

        // Creates a solverDomain's dictionary from a multiDict format:
        // -Merges "default" if available
        // -Merge "sameAs" or "multiLoad" dicts, recursively
        // -Merges inputDict.subDict("solverDomainName") last
        // outputDict - overwrites any existing entries
        // inputDict - the "multiSolver" subdict
        void buildDictionary
        (
            dictionary& outputDict,
            const dictionary& inputDict,
            const word& solverDomainName
        );

        // Change all catalogued multiDicts to another solverDomain
        void swapDictionaries(const word& solverDomainName);

        // When setting up for a different solverDomain, the boundary
        // conditions are allowed to change.  This function copies all valid
        // case/intoSolverDomain/inital/0 files, and overwrites the
        // internalField with those found in the corresponding file in
        // dataSourcePath.  The result is placed in case/[time].
        void swapBoundaryConditions
        (
            const fileName& dataSourcePath,
            const word& intoSolverDomain
        );

        // Read multiDict files into the multiDictRegistry.  Loads:
        //  - case/system/multi*
        //  - case/system/[local]/multi*, where [local] is anything
        //  - case/constant/multi*
        //  - case/constant/[local]/multi*
        void readAllMultiDicts();

        // Scan a directory for multi*; read them into the multiDictRegistry
        void readMultiDictDirectory
        (
            const fileName& sourcePath,
            const word& local = word::null
        );

        // Rereads modified dictionaries and sets the controls
        void readIfModified();

        // Converts a fileName with the naming convention:
        //      solverDomain@superLoop@globalOffset@globalIndex@fileName
        // to a timeCluster
        timeCluster parseConditionedFile
        (
            const word& pcFile,
            const instant& inst
        ) const;


public:

    // multiSolverControl enumerations

        //- Read initial data control options
        enum initialStartFromControls
        {
            misFirstTime,
            misFirstTimeInStartDomain,
            misFirstTimeInStartDomainInStartSuperLoop,
            misStartTime,
            misStartTimeInStartDomain,
            misStartTimeInStartDomainInStartSuperLoop,
            misLatestTime,
            misLatestTimeInStartDomain,
            misLatestTimeInStartDomainInStartSuperLoop
        };

        //- Final stop at control options
        enum finalStopAtControls
        {
            mfsEndTime,
            mfsEndTimeInEndDomain,
            mfsEndTimeInEndDomainInEndSuperLoop,
            mfsEndSuperLoop,
            mfsWriteNow,
            mfsNoWriteNow,
            mfsNextWrite
        };

    // solverDomains enumerations

        //- Time value start from control options
        enum startFromControls
        {
            mtsFirstTime,
            mtsStartTime,
            mtsLatestTimeThisDomain,
            mtsLatestTimeAllDomains
        };

        //- Stop at control options
        enum stopAtControls
        {
            msaEndTime,
            msaNoWriteNow,
            msaWriteNow,
            msaNextWrite,
            msaIterations,
            msaSolverSignal,
            msaElapsedTime
        };

protected:

    // multiSolverControl data

        label superLoop_;

        word currentSolverDomain_;

        static const NamedEnum<initialStartFromControls, 9>
            initialStartFromControlsNames_;
        initialStartFromControls initialStartFrom_;

        scalar initialStartTime_;

        word startDomain_;

        label startSuperLoop_;

        static const NamedEnum<finalStopAtControls, 7>
            finalStopAtControlsNames_;
        finalStopAtControls finalStopAt_;

        word endDomain_;

        scalar finalEndTime_;

        label endSuperLoop_;

        bool multiDictsRunTimeModifiable_;

        scalar globalTimeOffset_;

        //- Each solverDomain / superLoop combination is assigned a unique
        //  index that increments chronologically
        label globalIndex_;


    // solverDomains data
    // This data is transient, changing between solver domains

        static const NamedEnum<startFromControls, 4>
            startFromControlsNames_;
        startFromControls startFrom_;

        scalar startTime_;

        static const NamedEnum<stopAtControls, 7> stopAtControlsNames_;
        stopAtControls stopAt_;

        scalar endTime_;

        label purgeWriteSuperLoops_;

        scalar deltaT_;

        label iterations_;

        scalar elapsedTime_;

        wordList storeFields_;

public:

    TypeName("multiSolver");

    //- The default multiSolver dictionary name
    static word multiControlDictName;

/* Not implemented yet
    //- Indicates whether the 'solverSignal' option for the finalStopAt setting
    //  in the multiControlDict is permitted.  Default false.  Set this to true
    //  in the global space of a solver that does support solverSignal.
    static bool supportsSolverSignal;
*/

    // Constructors

        //- Construct given the multiControlDict
        multiSolver
        (
            const dictionary& dict,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            bool showSplash = true
        );

        //- Construct reading the multiControlDict from file
        multiSolver
        (
            const word& multiControlDictName,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            bool showSplash = true
        );


    // Destructor
        virtual ~multiSolver();


    // Member functions

        // Access

            //- Version number
            word version() const;

            // Database

            inline const Time& multiDictRegistry() const;

            inline const IOdictionary& multiControlDict() const;

            inline const wordList& prefixes() const;

            // multiSolverControl data

            inline const label& superLoop() const;

            inline const word& currentSolverDomain() const;

            inline const initialStartFromControls& initialStartFrom() const;

            inline const word& startDomain() const;

            inline const scalar& initialStartTime() const;

            inline const finalStopAtControls& finalStopAt() const;

            inline const word& endDomain() const;

            inline const scalar& finalEndTime() const;

            inline const label& startSuperLoop() const;

            inline const label& endSuperLoop() const;

            inline const bool& multiDictsRunTimeModifiable() const;

            inline const scalar& globalTimeOffset() const;

            // Write permission
            inline scalar& globalTimeOffset();

            inline const label& globalIndex() const;

            // solverDomains data
            inline const startFromControls& startFrom() const;

            inline const stopAtControls& stopAt() const;

            inline const scalar& startTime() const;

            inline const scalar& endTime() const;
            inline scalar& endTime();

            inline const label& iterations() const;

            inline const scalar& elapsedTime() const;

            inline const wordList& storeFields() const;

            inline const label& purgeWriteSuperLoops() const;

            inline const scalar& deltaT() const;
            inline scalar& deltaT();


        // Solver (and pre/post-processor) interface functions

            // Pre-condition the directory for decomposePar or reconstructPar
            void preCondition(const word& processor = word::null);

            // Post-condition the directory after decomposePar or
            // reconstructPar
            void postCondition(const word& processor = word::null);

            // Return the source data that multiControlDict currently points to
            // Returns "default" solverDomain if initial conditions are to be
            // used and no solverDomain is set yet.
            timeCluster initialDataSource() const;

            // Switch to another solver domain
            void setSolverDomain(const word& solverDomainName);

            // Switch to another solver domain for post-processing only
            void setSolverDomainPostProcessing(const word& solverDomainName);

            // Stop the run at the next setSolverDomain
            inline void setFinished();

            // Archive the last set of data into case/multiSolver
            void finalize();

            // Increment the superLoop (prefix)
            multiSolver& operator++();

            // Increment the superLoop (postfix)
            multiSolver& operator++(int);


        // Check

            //- Return true if run should continue
            bool run() const;

            //- Return true if end of run
            bool end() const;


        // Time functions
            // The multiSolver time directory structure should have the form:
            //  case
            //  '-multiSolver
            //    |-prefix1 {solverDomain}
            //    | |-initial  {initial directory, superLoop -1}
            //    | |-0 {superLoop}
            //    | |-1 {superLoop}
            //    | | |-multiSolverTime {auto-generated file}
            //    | | '-0.01, etc.. {time directories}
            //    | '-2, etc..
            //    |-prefix2, etc..

            // Create a list of all superLoops in a directory, (usually in
            // case/prefix). Only looks for integer directory names, does not
            // check for valid time subdirectories.
            // *** Located in multiSolverTimeFunctions.C ***
            static labelList findSuperLoops(const fileName& path);

            // Find the closest global time to a given value in a
            // timeClusterList. Assumes timeClusters do not overlap global time
            // values (as they shouldn't).  If exact is true, this function
            // throws a FatalError when no exact match is found.
            // *** Located in multiSolverTimeFunctions.C ***
            static timeCluster findClosestGlobalTime
            (
                const scalar value,
                const timeClusterList& tcl,
                const bool& exact = false
            );

            // Find the closest local time to a given value in a
            // timeClusterList. Unlike global times, local times can overlap.
            // If overlaps exist, it uses only the latest superloop.  If exact
            // is true, this function throws a FatalError when no exact match
            // is found.
            // *** Located in multiSolverTimeFunctions.C ***
            static timeCluster findClosestLocalTime
            (
                const scalar value,
                const timeClusterList& tcl,
                const bool& exact = false
           );

            // Find the latest global time
            // *** Located in multiSolverTimeFunctions.C ***
            static timeCluster findLatestGlobalTime
            (
                const timeClusterList& tcl
            );

            // Find the latest global time
            // *** Located in multiSolverTimeFunctions.C ***
            static timeCluster findLatestLocalTime
            (
                const timeClusterList& tcl
            );

            // Find the timeCluster with the given globalIndex.  Returns an
            // empty timeCluster if not found.
            // *** Located in multiSolverTimeFunctions.C ***
            static timeCluster findGlobalIndex
            (
                const label& index,
                const timeClusterList& tcl
            );

            // Adds all earlier times in the same superLoop directory to the
            // time cluster
            void includePreviousTimes
            (
                timeCluster& tc
            ) const;

            // Find the path to a specific entry in a time cluster
            // *** Located in multiSolverTimeFunctions.C ***
            fileName findInstancePath
            (
                const timeCluster& tc,
                const label& index
            ) const;

            // Find the largest superLoop
            // *** Located in multiSolverTimeFunctions.C ***
            static label findMaxSuperLoopValue(const timeClusterList& tcl);

            // Find the timeClusterList index for the timeClusterList that has
            // the largest superLoop
            // *** Located in multiSolverTimeFunctions.C ***
            static labelList findMaxSuperLoopIndices
            (
                const timeClusterList& tcl
            );

            // Checks if any of the time ranges overlap one another in a
            // timeClusterList.  (If startTime = previous end, this is okay.)
            // True means they do not overlap.
            // *** Located in multiSolverTimeFunctions.C ***
            static bool nonOverlapping
            (
                const timeClusterList& tcl,
                const bool useGlobalTime = false
            );

            // Maps the time directories in a single superLoop directory
            // Include a processor name, and it uses the processorN directory
            // *** Located in multiSolverTimeFunctions.C ***
            timeCluster readSuperLoopTimes
            (
                const word& solverDomain,
                const label superLoop,
                const word& processor = word::null
            ) const;

            // Maps the time directories in a single solverDomain
            // Include a processor name, and it uses the processorN directory
            // *** Located in multiSolverTimeFunctions.C ***
            timeClusterList readSolverDomainTimes
            (
                const word& solverDomain,
                const word processor = word::null
            ) const;

            // Create a snapshot of all the multiSolver time directories
            // Give it a processor name, and it searches instead in the
            // processor directory
            // *** Located in multiSolverTimeFunctions.C ***
            timeClusterList readAllTimes
            (
                const word processor = word::null
            ) const;

            // Move or copy all output time directories that are catalogued by
            // a timeClusterList to case/[timeValue] format.  Returns false if
            // time values overlap when forced to useGlobalTime.
            // loadStoreFields will copy storeFields into every time folder
            // even though they do not change.
            // *** Located in multiSolverTimeFunctions.C ***
            bool loadTimeClusterList
            (
                const timeClusterList& tcl,
                const bool useGlobalTime = true,
                const bool loadStoreFields = true
            );

            // Move all the time directories from sourcePath to archivePath
            // *** Located in multiSolverTimeFunctions.C ***
            static void archiveTimeDirs
            (
                const fileName& sourcePath,
                const fileName& archivePath,
                const label& purgeWrite
            );

            // Delete all time directories in path, do not delete "constant"
            // *** Located in multiSolverTimeFunctions.C ***
            static void purgeTimeDirs(const fileName& path);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiSolverI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
