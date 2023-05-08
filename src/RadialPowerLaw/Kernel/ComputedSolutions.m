(* ::Package:: *)

BeginPackage["ComputedSolutions`"]

Needs["DescriptionUtilities`"];

edgeSolution::usage = "Gives MK adn MM~ edge solutions given by Detournay data (zero buoyancy \[CapitalDelta]\[Gamma]=0). Input is the edge as string, the
fluid index n, the time t (list or scalar) and the location \[Rho].";
insideSolution::usage = "Gives the solution inside the rectangle in function of the trajectory parameter \[Phi] (zero buoyancy \[CapitalDelta]\[Gamma]=0). Input is the
fluid index n, the value of \[Phi], the time t (list or scalar) and the location \[Rho].";
solutionTime::usage = "Gives the interval of time where each solution is valid (zero buoyancy \[CapitalDelta]\[Gamma]=0). Sequence is edge MK, MMT or Inside and the fluid index n.";

Begin["`Private`"]

Get[$Path[[10]]<>"/Documents/Wolfram Mathematica/PyFrac_Mathematica_postprocessing/RadialPowerLaw/Kernel/MKEdgePowerLaw.txt"];
Get[$Path[[10]]<>"/Documents/Wolfram Mathematica/PyFrac_Mathematica_postprocessing/RadialPowerLaw/Kernel/InsidePowerLaw.txt"];
Get[$Path[[10]]<>"/Documents/Wolfram Mathematica/PyFrac_Mathematica_postprocessing/RadialPowerLaw/Kernel/MMtEdgePowerLaw.txt"];

edgeSolution[V_?StringQ,n_,t_,rho_]:=
	Module[
			{Scalar,Profile},
			
			If[
				Length[t]!=0,
				{
				Scalar = 
		     		{
		     			Time -> t, 
	   					FractureLength -> Transpose[{t,length[V][ToString[n]][t]}], 
	   					FractureVelocity -> Transpose[{t,fracvelocity[V][ToString[n]][t]}],
	   					FractureEfficiency -> Transpose[{t,efficiency[V][ToString[n]][t]}],
	   					InletOpening-> Transpose[{t,Omega1[V][ToString[n]][t]}],
   						InletPressure-> Transpose[{t,Pi1[V][ToString[n]][t]}],
   						InletFlux-> 1.(*Transpose[{t,Psi1[V][ToString[n]][t]}]*)
	   				};
	   			If[
	   				Length[rho]!=0,
	   				Profile =
   					{
   						Time -> #, 
     					Opening -> Transpose[{rho, Omegax[V][ToString[n]][#, rho]}], 
     					NetPressure -> Transpose[{rho, Pix[V][ToString[n]][#, rho]}],
     					FluidFlux -> Transpose[{rho, Psix[V][ToString[n]][#, rho]}]
     				} & /@ t,
     				Profile =
   					{
   						Time -> #, 
     					Opening -> {rho, Omegax[V][ToString[n]][#, rho]}, 
     					NetPressure -> {rho, Pix[V][ToString[n]][#, rho]},
     					FluidFlux -> {rho, Psix[V][ToString[n]][#, rho]}
     				} & /@ t
	   			]
				},
	   			{
	   			Scalar = 
		     		{
		     			Time -> t, 
	   					FractureLength -> {t,length[V][ToString[n]][t]},
	   					FractureVelocity -> {t,fracvelocity[V][ToString[n]][t]},
	   					FractureEfficiency -> {t,efficiency[V][ToString[n]][t]},
	   					InletOpening-> {t,Omega1[V][ToString[n]][t]},
   						InletPressure-> {t,Pi1[V][ToString[n]][t]},
   						InletFlux-> {t,Psi1[V][ToString[n]][t]}
	   				},
	   			If[
	   				Length[rho]!=0,
	   				Profile =
   					{
   						Time -> #, 
     					Opening -> Transpose[{rho, Omegax[V][ToString[n]][#, rho]}], 
     					NetPressure -> Transpose[{rho, Pix[V][ToString[n]][#, rho]}],
     					FluidFlux -> Transpose[{rho, Psix[V][ToString[n]][#, rho]}]
     				} & /@ Table[t,{1}],
     				Profile =
   					{
   						Time -> #, 
     					Opening -> {rho, Omegax[V][ToString[n]][#, rho]}, 
     					NetPressure -> {rho, Pix[V][ToString[n]][#, rho]},
     					FluidFlux -> {rho, Psix[V][ToString[n]][#, rho]}
     				} & /@ Table[t,{1}]
	   			]
	   			}
			]; 			
			
			{
 				( # -> (# /. Scalar) ) & /@ {Time, FractureLength, FractureVelocity, FractureEfficiency,InletOpening,InletPressure,InletFlux},
  				Profile
  			}
	
	];
	
insideSolution[n_,phi_,t_,rho_]:=
	Module[
		{Scalar,Profile},
		
		If[
			Length[t]!=0,
			{
			Scalar = 
	    	{
	     		Time -> t, 
   				FractureLength -> Transpose[{t,length["Inside"][ToString[n]][phi, t]}],
   				FractureVelocity -> Transpose[{t,fracvelocity["Inside"][ToString[n]][phi, t]}], 
   				FractureEfficiency -> Transpose[{t,efficiency["Inside"][ToString[n]][phi, t]}],
   				InletOpening-> Transpose[{t,Omega1["Inside"][ToString[n]][phi, t]}],
   				InletPressure-> Transpose[{t,Pi1["Inside"][ToString[n]][phi, t]}],
   				InletFlux-> Transpose[{t,Psi1["Inside"][ToString[n]][phi, t]}]
   			},
   			If[
   				Length[rho]!=0,
   				Profile = 
   					{
   						Time -> #, 
     					Opening -> Transpose[{rho, Omegax["Inside"][ToString[n]][phi, #, rho]}], 
     					NetPressure -> Transpose[{rho, Pix["Inside"][ToString[n]][phi, #, rho]}],
     					FluidFlux -> Transpose[{rho, Psix["Inside"][ToString[n]][phi, #, rho]}]
     				} & /@ t,
     			Profile = 
   					{
   						Time -> #, 
     					Opening -> {rho, Omegax["Inside"][ToString[n]][phi, #, rho]}, 
     					NetPressure -> {rho, Pix["Inside"][ToString[n]][phi, #, rho]},
     					FluidFlux -> {rho, Psix["Inside"][ToString[n]][phi, #, rho]}
     				} & /@ t
   			]
			},
			{
			Scalar = 
	    		{
	     			Time -> t, 
   					FractureLength -> {t,length["Inside"][ToString[n]][phi, t]},
   					FractureVelocity -> {t,fracvelocity["Inside"][ToString[n]][phi, t]}, 
   					FractureEfficiency -> {t,efficiency["Inside"][ToString[n]][phi, t]},
   					InletOpening-> {t,Omega1["Inside"][ToString[n]][phi, t]},
   					InletPressure-> {t,Pi1["Inside"][ToString[n]][phi, t]},
   					InletFlux-> {t,Psi1["Inside"][ToString[n]][phi, t]}
   				},
   			If[
   				Length[rho]!=0,
     			Profile = 
   				{
   					Time -> #, 
     				Opening -> Transpose[{rho, Omegax["Inside"][ToString[n]][phi, #, rho]}], 
     				NetPressure -> Transpose[{rho, Pix["Inside"][ToString[n]][phi, #, rho]}],
     				FluidFlux -> Transpose[{rho, Psix["Inside"][ToString[n]][phi, #, rho]}]
     			} & /@ Table[t,{1}],
     			Profile = 
   				{
   					Time -> #, 
     				Opening -> {rho, Omegax["Inside"][ToString[n]][phi, #, rho]}, 
     				NetPressure -> {rho, Pix["Inside"][ToString[n]][phi, #, rho]},
     				FluidFlux -> {rho, Psix["Inside"][ToString[n]][phi, #, rho]}
     			} & /@ Table[t,{1}]
   			]
			}
		];
		
		{
 			( # -> (# /. Scalar) ) & /@ {Time, FractureLength, FractureVelocity, FractureEfficiency,InletOpening,InletPressure,InletFlux},
  			Profile
  		}

	];

solutionTime[V_?StringQ,n_]=
	Switch[
			V,
			"MK",
				{InitialTime["MK"][ToString[n]], FinalTime["MK"][ToString[n]]},
			"MMt",
				{InitialTime["MMt"][ToString[n]],FinalTime["MMt"][ToString[n]]},
			"Inside",
				{InitialTime["Inside"][ToString[n]],FinalTime["Inside"][ToString[n]]}
		];

End[]

EndPackage[]






