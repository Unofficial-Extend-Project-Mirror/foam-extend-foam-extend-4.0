// use tetrahedral decomposition of the engine mesh

{
    tetPointVectorField& motionU = mSolver.motionU();

    vector pistonVel =
            piston().cs().axis()*engineTime_.pistonSpeed().value();

    // Set piston velocity
    if (piston().patchID().active())
    {

        if (debug)
        {
            Info << "Piston velocity: " << pistonVel;
        }

        componentMixedTetPolyPatchVectorField& pp =
            refCast<componentMixedTetPolyPatchVectorField>
            (
                motionU.boundaryField()[piston().patchID().index()]
            );

        pp.refValue() = pistonVel;
        
    }
    
    motionU.correctBoundaryConditions();
        
}
