// use tetrahedral decomposition of the engine mesh

{
    
    Info << "setting valve motion b.c." << endl;

    tetPointVectorField& motionU = mSolver.motionU();


    // Set valve velocity
    forAll (valves_, valveI)
    {
        
        vector valveVel =
            valves_[valveI].curVelocity()*valves_[valveI].cs().axis();

        if(valves_[valveI].curLift() > valves_[valveI].deformationLift())
        {
            valveVel = vector::zero;
        }
        
        Info << "Valve n. " << valveI + 1 << " velocity = " <<  valves_[valveI].curVelocity() << endl;

        
        // If valve is present in geometry, set the motion
        if (valves_[valveI].bottomPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            motionU.boundaryField()[valves_[valveI].bottomPatchID().index()] ==
                vector::zero;

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
            {
                Info<< "Valve " << valveI << " lift: "
                    << valves_[valveI].curLift()
                    << " velocity: " << valves_[valveI].curVelocity()
                    << endl;
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

        // If valve is present in geometry, set the motion
        if 
        (
            valves_[valveI].curtainInPortPatchID().active()
            &&
            valves_[valveI].poppetPatchID().active()
            &&
            valves_[valveI].curtainInCylinderPatchID().active()
        )
        {


            label curtainInPortIndex = valves_[valveI].curtainInPortPatchID().index();
            label poppetIndex = valves_[valveI].poppetPatchID().index();
            label curtainInCylIndex = valves_[valveI].curtainInCylinderPatchID().index();
            
            const coordinateSystem& cSysV = valves_[valveI].cs();
                        
            const pointField& curtainInPortGlobal = 
                    motionU.boundaryField()[curtainInPortIndex].patch().localPoints();
            
            const pointField& poppetGlobal = 
                    motionU.boundaryField()[poppetIndex].patch().localPoints();

            const pointField& curtainInCylGlobal = 
                    motionU.boundaryField()[curtainInCylIndex].patch().localPoints();

            scalarField curtainInPortLocal =  cSysV.localPosition(curtainInPortGlobal)().component(vector::Z);               
            scalarField poppetLocal =  cSysV.localPosition(poppetGlobal)().component(vector::Z);               
            scalarField curtainInCylLocal =  cSysV.localPosition(curtainInCylGlobal)().component(vector::Z);               
                
            scalar maxCurtain = max(curtainInCylLocal);
            scalar minPoppet = min(poppetLocal);
                
            Info << "maxCurtain = " << maxCurtain << endl;
            Info << "minPoppet = " << minPoppet << endl;
				
            motionU.boundaryField()[curtainInPortIndex] == (pos(curtainInPortLocal - minPoppet))*valveVel*(maxCurtain -
            curtainInPortLocal)/(maxCurtain-minPoppet);
            
               

//            if (debug)
            {
                Info<< "Valve " << valveI << " lift: "
                    << valves_[valveI].curLift()
                    << " velocity: " << valves_[valveI].curVelocity()
                    << endl;
            }
        }
		
	}

//  Setting the boundary position
    Info << "valve motion boundary conditions set" << endl;
    
}
