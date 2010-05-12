// use tetrahedral decomposition of the engine mesh

{
    
    Info << "setting valve motion b.c." << endl;

    tetPointVectorField& motionU = mSolver.motionU();

    
    // set piston motion boundary conditions
    
    // Set piston velocity
    if (piston().patchID().active())
    {
        vector pistonVel =
            piston().cs().axis()*engineTime_.pistonSpeed().value();

        if (debug)
        {
            Info << "Piston velocity: " << pistonVel;
        }

        motionU.boundaryField()[piston().patchID().index()] == pistonVel;

    }
    
    // Set valve velocity
    forAll (valves_, valveI)
    {
        
        vector valveVel =
            valves_[valveI].curVelocity()*valves_[valveI].cs().axis();
        
        Info << "Valve n. " << valveI + 1 << " velocity = " <<  valves_[valveI].curVelocity() << endl;
        
        // If valve is present in geometry, set the motion
        if (valves_[valveI].bottomPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            motionU.boundaryField()[valves_[valveI].bottomPatchID().index()] ==
                valveVel;

        }

        if (valves_[valveI].poppetPatchID().active())
        {
                            
            if(valves_[valveI].curLift() < valves_[valveI].deformationLift())
            {
                // Top of the valve does not move
                motionU.boundaryField()[valves_[valveI].poppetPatchID().index()] ==
                    valveVel;
            }
            else
            {
                // Top of the valve does not move
                motionU.boundaryField()[valves_[valveI].poppetPatchID().index()] ==
                    vector::zero;
            }
        
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

        if 
        (
            valves_[valveI].curtainInPortPatchID().active()
            &&
            valves_[valveI].poppetPatchID().active()
            &&
            valves_[valveI].curtainInCylinderPatchID().active()
            &&
            valves_[valveI].bottomPatchID().active()
        )
        {

//new

            label curtainInCylIndex = valves_[valveI].curtainInCylinderPatchID().index();
            label bottomIndex = valves_[valveI].bottomPatchID().index();
            
//            label poppetIndex = valves_[valveI].poppetPatchID().index();
//            label curtainInCylIndex = valves_[valveI].curtainInCylinderPatchID().index();
            
            const coordinateSystem& cSysV = valves_[valveI].cs();
                        
            const pointField& curtainInCylGlobal = 
                    motionU.boundaryField()[curtainInCylIndex].patch().localPoints();

            const pointField& bottomGlobal = 
                    motionU.boundaryField()[bottomIndex].patch().localPoints();
            
            scalarField curtainInCylLocal = cSysV.localPosition(curtainInCylGlobal)().component(vector::Z);               

            scalarField bottomLocal = cSysV.localPosition(bottomGlobal)().component(vector::Z);
                    
                
            scalar maxCurtain = max(curtainInCylLocal);

            scalar minCurtain = min(curtainInCylLocal);
//            scalar minCurtain = max(bottomLocal);
            scalar maxBottom = max(bottomLocal);

            
            Info << "minCurtain = " << minCurtain << endl;
            Info << "maxBottom = " << max(bottomLocal) << endl;
            
                				
//            motionU.boundaryField()[curtainInCylIndex] == (pos(curtainInCylLocal - minCurtain))*valveVel*(maxCurtain -
//            curtainInCylLocal)/(maxCurtain-minCurtain);
            

/*            
            forAll(motionU.boundaryField()[curtainInCylIndex], pointi)
            {
                if(curtainInCylLocal[pointi] > maxBottom)
                {
                    motionU.boundaryField()[curtainInCylIndex][pointi] 
                    == 
                    valveVel*
                    (maxCurtain - curtainInCylLocal[pointi])
                    /
                    (maxCurtain-maxBottom);
                }
                else
                {
                    motionU.boundaryField()[curtainInCylIndex][pointi] 
                    == 
                    valveVel;
                }
            }
*/


/*
            motionU.boundaryField()[curtainInCylIndex] == valveVel*(maxCurtain -
            curtainInCylLocal)/(maxCurtain-minCurtain);
*/
            motionU.boundaryField()[curtainInCylIndex] == pos(curtainInCylLocal - maxBottom) * valveVel*(maxCurtain -
            curtainInCylLocal)/(maxCurtain-minCurtain) + (1.0 - pos(curtainInCylLocal - maxBottom)) * valveVel;

//            Info << "motionU" << motionU.boundaryField()[curtainInCylIndex] << endl;
            
            Info << "valveVel" << valveVel << endl;
                        
//            Info << "motionU.boundaryField()[curtainInCylIndex]" << 
//            motionU.boundaryField()[curtainInCylIndex] << endl;

        }
		
	}

//  Setting the boundary position
    Info << "valve motion boundary conditions set" << endl;
    
}
