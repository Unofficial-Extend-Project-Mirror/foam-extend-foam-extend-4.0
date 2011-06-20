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

#include "Time.H"
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
        Time multiDictRegistry_;

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
        
        // True when the end condition has been met
        bool finished_;
        
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
        // which may differ from that in multiControlDict.  Sets finished_ to
        // true if the end condition will be met.
        word setLocalEndTime();

        // Looks for case/multiSolver/solverDomain/initial/0.  Fail if missing.
        void checkTimeDirectories() const;

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
        //      solverDomain@superLoop@globalOffset@fileName
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
            const word& constantName = "constant"
        );

        //- Construct reading the multiControlDict from file
        multiSolver
        (
            const word& multiControlDictName,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant"
        );


    // Destructor
        virtual ~multiSolver();


    // Member functions
    
        // Access

            // Database
            
            inline const Time& multiDictRegistry() const;

            inline const IOdictionary multiControlDict() const;

        
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
            
            // solverDomains data
            inline const startFromControls& startFrom() const;
            
            inline const stopAtControls& stopAt() const;
            
            inline const scalar& startTime() const;

            inline const scalar& endTime() const;
            
            inline const label& iterations() const;
            
            inline const scalar& elapsedTime() const;
            
            inline const wordList& storeFields() const;
            
            inline const label& purgeWriteSuperLoops() const;

            inline const scalar& deltaT() const;
            

        // Solver (and pre/post-processor) interface functions

            // Pre-condition the directory for decomposePar or reconstructPar
            void preCondition(const word& processor = word::null);

            // Post-condition the directory after decomposePar or
            // reconstructPar
            void postCondition(const word& processor = word::null);

            // Switch to another solver domain
            void setSolverDomain(const word& solverDomainName);

            // Switch to another solver domain for post-processing only
            void setSolverDomainPostProcessing(const word& solverDomainName);

            // Stop the run at the next setSolverDomain
            inline void setFinished();
            
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
            //    | | '-multiSolverTime {auto-generated dictionary}
            //    | '-2, etc..
            //    |-prefix2, etc..
            //
            // *** All time functions are located in ***
            // ***    multiSolverTimeFunctions.C     ***

            // Create a list of all superLoops in a directory, (usually in
            // case/prefix). Only looks for integer directory names, does not
            // check for valid time subdirectories.
            static labelList findSuperLoops(const fileName& path);

            // Find the closest global time to a given value in a
            // timeClusterList. Assumes timeClusters do not overlap global time
            // values (as they shouldn't).  If exact is true, this function
            // throws a FatalError when no exact match is found.
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
            static timeCluster findClosestLocalTime
            (
                const scalar value,
                const timeClusterList& tcl,
                const bool& exact = false
           );

            // Find the latest global time
            static timeCluster findLatestGlobalTime
            (
                const timeClusterList& tcl
            );

            // Find the latest global time
            static timeCluster findLatestLocalTime
            (
                const timeClusterList& tcl
            );

            // Find the path to a specific entry in a time cluster
            fileName findInstancePath
            (
                const timeCluster& tc,
                const label& index
            ) const;

            // Find the largest superLoop
            static label findMaxSuperLoopValue(const timeClusterList& tcl);

            // Find the timeClusterList index for the timeClusterList that has
            // the largest superLoop
            static labelList findMaxSuperLoopIndices(const timeClusterList& tcl);

            // Checks if any of the time ranges overlap one another in a
            // timeClusterList.  (If startTime = previous end, this is okay.)
            // True means they do not overlap.
            static bool nonOverlapping
            (
                const timeClusterList& tcl,
                const bool useGlobalTime = false
            );

            // Maps the time directories in a single superLoop directory
            // Include a processor name, and it uses the processorN directory
            timeCluster readSuperLoopTimes
            (
                const word& solverDomain,
                const label superLoop,
                const word& processor = word::null
            ) const;
            
            // Maps the time directories in a single solverDomain
            // Include a processor name, and it uses the processorN directory
            timeClusterList readSolverDomainTimes
            (
                const word& solverDomain,
                const word processor = word::null
            ) const;

            // Create a snapshot of all the multiSolver time directories
            // Give it a processor name, and it searches instead in the
            // processor directory
            timeClusterList readAllTimes
            (
                const word processor = word::null
            ) const;

            // Move or copy all output time directories that are catalogued by
            // a timeClusterList to case/[timeValue] format.  Returns false if
            // time values overlap when forced to useGlobalTime.
            // loadStoreFields will copy storeFields into every time folder
            // even though they do not change.
            bool loadTimeClusterList
            (
                const timeClusterList& tcl,
                const bool useGlobalTime = true,
                const bool loadStoreFields = true
            ) const;

            // Move all the time directories from sourcePath to archivePath
            static void archiveTimeDirs
            (
                const fileName& sourcePath,
                const fileName& archivePath,
                const label& purgeWrite
            );

            // Delete all time directories in path, do not delete "constant"
            static void purgeTimeDirs(const fileName& path);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiSolverI.H"

#endif
