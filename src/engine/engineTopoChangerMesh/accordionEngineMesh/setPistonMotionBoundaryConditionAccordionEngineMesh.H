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

        motionU.boundaryField()[piston().patchID().index()] == pistonVel;
//        motionU.boundaryField()[piston().patchID().index()] == vector::zero;

    }
        
}
