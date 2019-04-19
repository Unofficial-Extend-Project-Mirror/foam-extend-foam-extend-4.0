// Check convergence
if (maxResidual < convergenceCriterion)
{
    Info<< "reached convergence criterion: " << convergenceCriterion << endl;

    // Move to the next time-step
    break;
}

