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
        
        if
        (
	        moveDetach_ 
	        && 
            valves_[valveI].detachInPortPatchID().active()
            && 
            valves_[valveI].detachInCylinderPatchID().active()
        )
        {
        
            label detachInPortIndex = valves_[valveI].detachInPortPatchID().index();
            label detachInCylIndex = valves_[valveI].detachInCylinderPatchID().index();

            const coordinateSystem& cSysV = valves_[valveI].cs();
                        
            const pointField& detachInPortGlobal = 
                    motionU.boundaryField()[detachInPortIndex].patch().localPoints();
            
            const pointField& detachInCylGlobal = 
                    motionU.boundaryField()[detachInCylIndex].patch().localPoints();

            scalarField detachInPortLocal =  cSysV.localPosition(detachInPortGlobal)().component(vector::Z);               
            scalarField detachInCylLocal =  cSysV.localPosition(detachInCylGlobal)().component(vector::Z);               
                
            scalar maxDetachInPort = max(detachInPortLocal);
            scalar minDetachInPort = min(detachInPortLocal);

            scalar maxDetachInCyl = max(detachInCylLocal);
            scalar minDetachInCyl = min(detachInCylLocal);
                
            Info << "maxDetachInPort = " << maxDetachInPort << endl;
            Info << "minDetachInPort = " << minDetachInPort << endl;
            Info << "maxDetachInCyl = " << maxDetachInCyl << endl;
            Info << "minDetachInCyl = " << minDetachInCyl << endl;
				
				
            motionU.boundaryField()[detachInPortIndex] == valveVel*(maxDetachInPort -
            detachInPortLocal)/(maxDetachInPort-minDetachInPort);
            
            motionU.boundaryField()[detachInCylIndex] == valveVel*(maxDetachInCyl -
            detachInCylLocal)/(maxDetachInCyl-minDetachInCyl);
        }
        else
	    {
/*        	if (valves_[valveI].detachInPortPatchID().active())
        	{
            		motionU.boundaryField()[valves_[valveI].detachInPortPatchID().index()] == vector::zero;
        	}
        	Info << "detach in port" << endl;

        	if (valves_[valveI].detachInCylinderPatchID().active())
        	{	
            		motionU.boundaryField()[valves_[valveI].detachInCylinderPatchID().index()] == vector::zero;
       	 	}
        	Info << "detach in cyl" << endl;
*/	        if
	        (
	            !moveDetach_
  	            &&
                valves_[valveI].detachInPortPatchID().active()
                && 
                valves_[valveI].detachInCylinderPatchID().active()
	        )

	        {
        	    motionU.boundaryField()[valves_[valveI].detachInPortPatchID().index()] == vector::zero;
        	    motionU.boundaryField()[valves_[valveI].detachInCylinderPatchID().index()] == vector::zero;
	        }
		
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
