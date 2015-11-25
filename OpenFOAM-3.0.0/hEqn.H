{
    /*
     *  Bounded flux of enthalpy or internal energy (u*gamma)
     */
    upwind<scalar> UDsPos(mesh, phiPos);
    upwind<scalar> UDsNeg(mesh, phiNeg);

    fvScalarMatrix hConvection
    (
	fvm::ddt(rho,h)
	+
	fv::gaussConvectionScheme<scalar>(mesh, phiPos, UDsPos).fvmDiv(phiPos, h)
	+
	fv::gaussConvectionScheme<scalar>(mesh, phiNeg, UDsNeg).fvmDiv(phiNeg, h)
    );
    
    /*
     *
     * Explicit flux correction, bounded with MULES
     *
     */
    surfaceScalarField hFluxCorr
    (
	lambdaCoeffs
	*
	(
	    fvc::flux
	    (
		phiPos,
		h,
		"div(phiPos," + h.name() + ")"
	    )
	    +
	    fvc::flux
	    (
		phiNeg,
		h,
		"div(phiNeg," + h.name() + ")"
	    )
	    -
	    hConvection.flux()
	)
    );

    /*
     *
     * Diffusive fluxes of energy due to mass diffusion
     *
     */
    PtrList<volScalarField> hi (Y.size());
    forAll (hi, iCmpt)
    {
	hi.set
	(
	    iCmpt,
	    new volScalarField
	    (
		IOobject
		(
		    "h" + Y[iCmpt].name(),
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh,
		dimEnergy/dimMass
	    )
	);
	
	scalarField& hiIF = hi[iCmpt].internalField();
	const scalarField& pIF  = p.internalField();
	const scalarField& TIF  = T.internalField();
	
	forAll(hiIF, iCell)
	{
	    hiIF[iCell] = thermo.composition().Hs(iCmpt, pIF[iCell], TIF[iCell]);
	}
	
	forAll(hi[iCmpt].boundaryField(), iPatch)
	{
	    fvPatchScalarField& hip = hi[iCmpt].boundaryField()[iPatch];
	    const fvPatchScalarField& pp = p.boundaryField()[iPatch];
	    const fvPatchScalarField& Tp = T.boundaryField()[iPatch];
	    forAll(hip, iFace)
	    {
		hip[iFace] = thermo.composition().Hs(iCmpt, pp[iFace], Tp[iFace]);
	    }
	}
    }
    
    surfaceScalarField dzetaPhi
    (
	"dzetaPhi",
	fvc::flux
	(
	    diffusiveFlux[0],
	    hi[0],
	    "div(rhoi*Uri,hi)"
	)
    );
    
    for (label iCmpt = 1; iCmpt < Y.size(); iCmpt++)
    {
	dzetaPhi +=
	fvc::flux
	(
	    diffusiveFlux[iCmpt],
	    hi[iCmpt],
	    "div(rhoi*Uri,hi)"
	);
    }

    volScalarField alphahEff ("alphahEff", turbulence->alphaEff());
    volScalarField Cp ("Cp", thermo.Cp());

    fvScalarMatrix hEqn
    (
	hConvection + fvc::div(hFluxCorr) + fvc::div(dzetaPhi)
	+
	EkChange
	==
	dpdt
	+
	fvm::laplacian(turbulence->alphaEff(), h)
	-
	fvc::laplacian(alphahEff*T,Cp)
	+
	fvc::div( ((-turbulence->devRhoReff()) & U) )
	+
	rho*(U&g)
	+
	reaction->Sh()
	+
	fvOptions(rho, h)
    );

    hEqn.relax();

    fvOptions.constrain(hEqn);

    hEqn.solve();

    fvOptions.correct(h);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
