/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    reactingCentralFoam

Description
    Solver for combustion with chemical reactions with the hybrid KT/PISO
    implicit scheme for advection

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulenceModel.H"
#include "psiCombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "coupledFvsPatchFields.H"
#include "cellQuality.H"
#include "customMULES.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"

    #include "readCourantType.H"
    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);
    Info<< "\nStarting time loop\n" << endl;
    #include "createSurfaceFields.H"
    #include "markBadQualityCells.H"

    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    pimpleControl pimple(mesh);
    
    

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "acousticCourantNo.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
	rho.oldTime();
	p.oldTime();
	U.oldTime();
	h.oldTime();
	psi.oldTime();
	forAll(Y, iCmpt)
	{
	    Y[iCmpt].oldTime();
	}

        {
            fvScalarMatrix rhoEqn
            (
                fvm::ddt(rho) + fvc::div(phi)
                ==
                fvOptions(rho)
            );
            
            fvOptions.constrain(rhoEqn);
            
            rhoEqn.solve();
            
            fvOptions.correct(rho);
        }
	
        while (pimple.loop())
        {
	    // Create limiter field for mass and energy fluxes
	    scalarField allFacesLambda(mesh.nFaces(), 1.0);
	    slicedSurfaceScalarField lambdaCoeffs
	    (
		IOobject
		(
		    "lambdaCoeffs",
		    mesh.time().timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE,
		    false
		),
		mesh,
		dimless,
		allFacesLambda,
		false   // Use slices for the couples
	    );
	    
	    #include "YEqn.H"
	    #include "UEqn.H"
	    #include "hEqn.H"
	    
	    // --- Pressure corrector loop
	    while (pimple.correct())
	    {
		#include "pEqn.H"
	    }
	    
	    #include "updateKappa.H"
	    
	    K = 0.5*magSqr(U);
	    EkChange = fvc::ddt(rho,K) + fvc::div(phiPos,K) + fvc::div(phiNeg,K);
	    dpdt = fvc::ddt(p);
	    
	    if (pimple.turbCorr())
	    {
		turbulence->correct();
	    }
	}

	if(runTime.write())
	{
	    psi.write();
	    rho.write();
	}

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
