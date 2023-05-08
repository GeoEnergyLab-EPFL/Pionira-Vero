(* ::Package:: *)

(* Mathematica Package *)

BeginPackage["DescriptionUtilities`"]
(* Exported symbols added here with SymbolName::usage *)  

(*
The goal of this package is to set up a number of structtures and Protected name for the descriptions:
i) of the solution of the problem : Structure of 'type' HFdescription
ii) of the model definition : geometry, kernels etc.
iii) the input parameters of the model

and a number of functions associated 

It also provides functions for  Handling Structure construct.

*)

(* Field values *)
Time::usage ="Symbol reserved for current time "
Opening::usage = "Symbol reserved for fracture opening"
FluidPressure::usage ="Symbol reserved for fluid pressure in fracture"
NetPressure::usage ="Symbol reserved for net pressure in fracture"
NormalStress::usage = " Symbol reserved for normal stress to the fracture"
FluidFlux::usage = "Symbol reserved for fluid flux"
PfGradient::usage="Symbol reseverd for fluid pressure gradient"
WellborePressure::usage = " Symbol reserved for fluid wellbore pressure"
FractureLength::usage="Symbol reserved for fracture length"
FractureVelocity::usage="Symbol reserved for fracture velocity" 
FractureEfficiency::usage= "Symbol reserved for fracture efficiency"
FractureVolume::usage = " Symbol reserved for fracture volume"

InletFlux::usage = " Symbol reserved for the flux entering the fracture"
InletPressure::usage = " Symbol reserved for the net pressure at the inlet of the fracture"
InletOpening::usage = " Symbol reserved for the opening at the inlet of the fracture"
OpeningIncrement::usage = " Symbol reserved for opening increment"

LeakOff::usage = " Symbol reserved for the cumulated volume leaked off the fracture at the current time "
InjectedVolume::usage = " Symbol reserverd for the cumulated injected volume inside the fracture at the current time"
LeakOffRate::usage = " Symbol reserved for leak off rate"
LeakOffVolume::usage = " Symbol lreserved for leak off volume in elts over a time step "
TriggerTime::usage = "symbol reserved for the leak-off trigger time "
TipAsymptote::usage = " symbol reserved for tip asymptote description "

CurrentGrid::usage = " Symbol reserved for list of mid-points grid cell coordinates on the fixed grid"
CurrentTimeStep::usage = " Symbol reserved for the current time step of the simulation (in the corresponding results list)"

gamma::usage = "Symbol reserved for dimensionless fracture length or radius "
    

(* Model definition *)
Geometry::usage = "Symbol reserved for fracture geometries ... Radial, RadialWellbore, PlaneStrain, PlaneStrainWellbore "
FlowFactor::usage = " Symbol reserved for the flow factor coefficient"
NormalStressFunction::usage = " Symbol reserved for the Normal Stress Function"
TipFunction::usage = "symbol reserved for the name of the tip volume function"

(* Parameters *)
NIndex::usage = " Symbol reserved for fluid consistency index"
KConsistency::usage = " Symbol reserved for fluid consistency coeffiecient k " 
YieldStress::usage = " Symbol reserved for fluid yield stress "

FluidProperties::usage = "Symbol reserved for fluid properties "

SigmaO::usage = " Symbol reserved for initial constant clamping stress normal to the fracture"


(* Dimensionless coefficients *)
\[ScriptCapitalK]::usage = "Dimensionless toughness coefficient."
\[ScriptCapitalM]::usage = "Dimensionless viscosity coefficient."
\[ScriptCapitalC]::usage = "Dimensionless leak-off coefficient."
\[ScriptCapitalV]::usage = "Dimensionless storage coefficient."
\[ScriptCapitalB]::usage = "Dimensionless buoyancy coefficient."


(* Buoyancy related parameters *)
\[CapitalDelta]\[Gamma]::usage = "Buoyancy of the gravity driven fracture."
Vh::usage = "Head volume for a buoyant fracture."
Vt::usage = "Tail volume for a buoyant fracture."


WellboreRadius::usage = " Symbol reserved for wellbore radis "

Ep::usage = " Symbole reserved for the plane strain Young's modulus Ep=E/(1-nu^2) "
Emod::usage = " Symbole reserved for the Young's modulus"
\[Nu]::usage = "Symbol reserved for the Poisson's ratio"
Kp::usage = "Symbol reserved for the equivalent toughness Kp = Sqrt[32/Pi] Kic "
KIc::usage = "Symbol reserved for the fracture toughness "
Mp::usage = "Symbol reserved for the fluid dimensionless viscosity"
\[Mu]p::usage = "Symbol reserved for the fluid dimensionless viscosity"
\[Mu]::usage = "Symbol reserved for the fluid viscosity"
Vstar::usage = "Symbol reserved for the critical crack velocity in sub-critical crack growth"
Nsubcritical::usage = " Symbol reserved for the sub-critical growth power law index"
n::usage = "Symbol reserved for the fluid index n."

Qo::usage =" Symbol reserved for the (total) injected flux "
Vo::usage =" Symbol reserved for the (total) injected fluid"
ts::usage =" Time of the shut-in of the pump."
InjectionRate::usage =" Symbol reserved for the (total) injected flux "
Qin::usage = "Symbol reserved for the flux entereing the fracture"
Cp::usage = "Symbol reserved for Total Perf Friction factor  Pa / (m^(2+d)/s)"
U::usage = " Symbol reserved for injection system compressibility " 
Cl::usage = "Symbol reserved for Carter Leak off coefficient m t^-1/2 "
S2::usage = "Symbol reserved for minimum stress" 
S1::usage = "Symbol reserved for max stress"

theta::usage = "Symbol reserved for the inclination angle of the defect from the wellbore in plane-strain"

lo::usage = " Symbol reserved for initial defect length "

AxialPosition::usage = " Symbol reserved for postion of the axial positions of different fracs in an array ";

PoissonRatio::usage = "Symbol reserved for Poisson's ratio"  (*   *)

WellboreRadius::usage = "Option for the kernel TensileKernelWellboreSym"

InitialHydraulicOpening::usage = " Symbol reseverd for initial hydraulic width, must be a list"

HFarray::usage = "Symbol reserved for the case of multiple fractures - used as switch in different part of the code"

(* Scales *)
Scales::usage = "Symbol reserved for the Scales definition"
DimensionlessNumbers::usage = "Symbol reserved for the dimensionless Number "
wstar::usage = " Symbol reserved for the opening scale "
pstar::usage = " Symbol reserved for the pressure scale "
Lstar::usage = " Symbol reserved for the length scale "
tstar::usage = " Symbol reserved for the characteristic time scale"
qstar::usage= " Symbol reserved for the flux scale "
Lstarv::usage = " Symbol reserved for the vertical length scale in buoyancy scalings "
Lstarh::usage = " Symbol reserved for the horizontal length scale in buoyancy scalings "

(* dimensionless numbers & symbol*)
Kbar::usage = " Symbol reserved for dimensionless toughness"
Mbar::usage = " Symbol reserved for dimensionless viscosity"
Cbar::usage = " Symbol reserved for dimensionless leak-off"

(* Functions*)
CreateHFDescription::usage = " CreateHFDescription[Name,ListRule] " 
DuplicateDescription::usage = "DuplicateDescription[Name_Symbol, NewName_Symbol, ListSymbol_List] duplicate a description Name@ListSymbol[[;;]]"
AppendValues::usage = "AppendValues[Name_Symbol, Args_List, Val_List] Append values to some of the args " 
IncreaseValues::usage = "IncreaseValues[Name_Symbol, Args_List,Incval_List] increment of a given value Name@Args by Incval, Args and Incval must be list of symbol, resp. values" 
ChangeValues::usage ="ChangeValues[Name_Symbol, Args_List,Val_List]  change the value of Name@Args with Val, Args and Val must be list of symbol, resp. values"

Begin["`Private`"] (* Begin Private Context *) 

Protect[#] & /@ {Opening, FluidPressure, NetPressure, NormalStress, FluidFlux, PfGradient, WellborePressure, FractureLength,
 InletFlux,LeakOff,InjectedVolume,TriggerTime,TipAsymptote,Geometry,LeakOffRate,FractureVolume,WellboreRadius,
 FractureEfficiency,LeakOffVolume} ;

Protect[#] & /@ {NIndex,KConsistency,YieldStress,SigmaO,Geometry,NormalStressFunction,FlowFactor,PoissonRatio,
Ep,Kp,Qin,Qo,S1,S2,Cl,lo,AxialPosition,Nsubcritical,Vstar,InitialHydraulicOpening};

Protect[#] & /@ {wstar,pstar,Lstar,qstar,Kbar,Mbar,Scales,DimensionlessNumbers,gamma}

CreateHFDescription[Name_,ListRule_]:= Module[{},
(* Name is a symbol *)
(* ListRule is a list of rule *)	
	(Name[#] = # /. ListRule) & /@ {Opening,FluidPressure,NormalStress,FluidFlux,FractureLength} ;
];

(* ------------------------- *)   
DuplicateDescription[Name_Symbol, NewName_Symbol] := (DownValues[NewName]=DownValues[Name] /. Name -> NewName );

 (* ------------------------- *)
AppendValues[Name_Symbol, Args_List, Val_List] :={(Name[Args[[#]]]=Flatten[Append[Name[Args[[#]]], Val[[#]]]]) & /@ Range[Length[Val]];} ;

(* ------------------------- *)
ChangeValues[Name_Symbol, Args_List,Val_List] := {(Name[Args[[#]]] = Val[[#]]) & /@ Range[Length[Val]];};

(* ------------------------- *)
IncreaseValues[Name_Symbol, Args_List,Incval_List] := {(Name[Args[[#]]] += Incval[[#]]) & /@ Range[Length[Incval]];};


End[] (* End Private Context *)

EndPackage[]
