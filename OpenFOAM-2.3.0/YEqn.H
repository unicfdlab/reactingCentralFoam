{
    /*
     *
     * Create fluxes and sources
     *
     */
    reaction->correct();
    dQ = reaction->dQ();
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);
    const label nComponents = Y.size();

    PtrList<fvScalarMatrix> SuSp(nComponents);
    
    forAll(Y, iCmpt)
    {
	volScalarField& Yi = Y[iCmpt];
	if (Yi.name() != inertSpecie)
	{
	    word YiPosScheme ("div(phiPos,Yi)");
	    word YiNegScheme ("div(phiNeg,Yi)");
	    
	    phiY.set
	    (
		iCmpt,
		fvc::flux
		(
		    phiPos,
		    Yi,
		    YiPosScheme
		)
		+
		fvc::flux
		(
		    phiNeg,
		    Yi,
		    YiNegScheme
		)
	    );
	    
	    SuSp.set
	    (
		iCmpt,
		new fvScalarMatrix
		(
		    reaction->R(Yi)
		    +
		    fvOptions(rho,Yi)
		)
	    );
	    
	    mDCf[iCmpt] = fvc::interpolate(turbulence->muEff());
	}
	else
	{
	    inertIndex = iCmpt;
	    
	    forAll(diffusiveFlux[iCmpt].internalField(), iFace)
	    {
		diffusiveFlux[iCmpt].internalField()[iFace] = 0.0;
	    }
	    forAll(diffusiveFlux[iCmpt].boundaryField(), iPatch)
	    {
		forAll(diffusiveFlux[iCmpt].boundaryField()[iPatch], iFace)
		{
		    diffusiveFlux[iCmpt].boundaryField()[iPatch][iFace] = 0.0;
		}
	    }
	}
    }


    /*
     *
     * Calculate the limiter
     *
     */
    forAll(Y, iCmpt)
    {
	volScalarField& Yi = Y[iCmpt];
	if ( Yi.name() != inertSpecie )
	{
	    surfaceScalarField& rhoPhiYCorr = phiY[iCmpt];
	
	    mulesWithDiffusionImplicitLimiter
	    (
		rho,
		Yi,
		phiPos,
		phiNeg,
		allFacesLambda,
		rhoPhiYCorr,
		diffusiveFlux[iCmpt],
		mDCf[iCmpt],
		SuSp[iCmpt]
	    );
	}
    }
    
    /*
     *
     * Solve for components
     *
     */
    upwind<scalar> UDsPos(mesh, phiPos);
    upwind<scalar> UDsNeg(mesh, phiNeg);
    
    forAll(Y, iCmpt)
    {
	volScalarField& Yi = Y[iCmpt];
	if (Yi.name() != inertSpecie)
	{
	    surfaceScalarField& rhoPhiYCorr = phiY[iCmpt];
	    fvScalarMatrix YConvection
	    (
		fvm::ddt(rho,Yi)
		+
		fv::gaussConvectionScheme<scalar>(mesh, phiPos, UDsPos).fvmDiv(phiPos, Yi)
		+
		fv::gaussConvectionScheme<scalar>(mesh, phiNeg, UDsNeg).fvmDiv(phiNeg, Yi)
		==
		SuSp[iCmpt]
	    );
	
	    fvScalarMatrix YDiffusion
	    (
		- fvm::laplacian(mDCf[iCmpt], Yi)
	    );

	    fvScalarMatrix YiEqn
	    (
		YConvection + YDiffusion +  fvc::div(lambdaCoeffs * rhoPhiYCorr)
	    );
	    
	    YiEqn.relax();
	    
	    fvOptions.constrain(YiEqn);
	    
	    YiEqn.solve(mesh.solver("Yi"));
	    
	    fvOptions.correct(Yi);
	    
	    Yi.max(0.0);
	    Yt += Yi;
	
	    phiY[iCmpt] = YConvection.flux() + lambdaCoeffs * rhoPhiYCorr;
	    
	    diffusiveFlux[iCmpt] = YDiffusion.flux();
	    diffusiveFlux[inertIndex] -= diffusiveFlux[iCmpt];
	    
	    Info << "max/min " << Yi.name() << ": " << max(Yi).value() << "/" << min(Yi).value() << endl;
	}
    }
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].correctBoundaryConditions();
    Y[inertIndex].max(0.0);
}
