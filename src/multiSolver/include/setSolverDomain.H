        } // previous solver domain goes out of scope
        multiRun.setSolverDomain(solverDomain);

// Clear defines that may interfere with other solver domains 
#undef createPhi_H
#undef createPhiV_H
#undef initContinuityErrs_H

        if (multiRun.run())
        { // next solver domain comes into scope
