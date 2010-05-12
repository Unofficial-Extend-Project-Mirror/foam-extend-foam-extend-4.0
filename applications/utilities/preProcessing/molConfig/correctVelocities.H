for (molN = totalMols; molN < totalMols + totalZoneMols; molN++)
{

    // Remove bulk momentum introduced by random numbers and add
    // desired bulk velocity

    // For systems with molecules of significantly differing masses, this may
    // need to be an iterative process or employ a better algorithm for
    // removing an appropriate share of the excess momentum from each molecule.

    initialVelocities(molN) += bulkVelocity - momentumSum/totalZoneMols/mass;
}

// momentumSum = vector::zero;
//
// for (molN = totalMols; molN < totalMols + totalZoneMols; molN++)
// {
//     momentumSum += mass*initialVelocities(molN);
// }
//
// Info << "Check momentum adjustment: " << momentumSum << endl;
