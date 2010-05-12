// use tetrahedral decomposition of the engine mesh

{
    
    Info << "setting valve motion b.c." << endl;

    tetPointVectorField& motionU = mSolver.motionU();

    // Set valve velocity
    forAll (valves_, valveI)
    {
        
        Info << "Valve n. " << valveI + 1 << " velocity = " <<  valves_[valveI].curVelocity() << endl;

        vector valveVel =
            valves_[valveI].curVelocity()*valves_[valveI].cs().axis();
        
        // If valve is present in geometry, set the motion
        if (valves_[valveI].bottomPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            motionU.boundaryField()[valves_[valveI].bottomPatchID().index()] ==
                valveVel;

//            if (debug)
            {
                Info<< "Valve " << valveI << " lift: "
                    << valves_[valveI].curLift()
                    << " velocity: " << valves_[valveI].curVelocity()
                    << endl;
            }
        }

        if (valves_[valveI].poppetPatchID().active())
        {
                            
            // Top of the valve does not move
            motionU.boundaryField()[valves_[valveI].poppetPatchID().index()] ==
                valveVel;
        }

        if (valves_[valveI].sidePatchID().active())
        {
            // Top of the valve does not move
            motionU.boundaryField()[valves_[valveI].sidePatchID().index()] ==
                valveVel;
        }
        
        if (valves_[valveI].detachInPortPatchID().active())
        {
            motionU.boundaryField()[valves_[valveI].detachInPortPatchID().index()] == vector::zero;
        }

        if (valves_[valveI].detachInCylinderPatchID().active())
        {
            motionU.boundaryField()[valves_[valveI].detachInCylinderPatchID().index()] == vector::zero;
        }

		
	}

//  Setting the boundary position

    {
        label cylinderHeadIndex = boundaryMesh().findPatchID(cylinderHeadName_);
        
        // Top of the valve does not move
        motionU.boundaryField()[cylinderHeadIndex] ==
            vector::zero;
    }
    
    Info << "valve motion boundary conditions set" << endl;
    
}
