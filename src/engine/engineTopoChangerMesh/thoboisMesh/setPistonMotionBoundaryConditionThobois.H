// use tetrahedral decomposition of the engine mesh

    vector pistonVel =
            piston().cs().axis()*engineTime_.pistonSpeed().value();

{
    tetPointVectorField& motionU = mSolver.motionU();
    
    Info << "setting the piston velocity" << endl;

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

        pp.refValue() = vector::zero;
        
    }
        
}
