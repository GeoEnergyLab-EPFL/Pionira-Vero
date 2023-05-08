(* ::Package:: *)

(* --------------------------------------------- 

Wolfram Language Package created by Andreas M\[ODoubleDot]ri
on Wednesday 18th November 2020
EPFL ENAC IIC GEL
GC B1 385 (Batiment GC)
Station 18, 1015 Lausanne, Switzerland
Phone: +41 (0)21 693 44 70
andreas.mori@epfl.ch
Geo-energy Laboratory EPFL 

--------------------------------------------- *)
BeginPackage["UtilityForScalings`"];

inputDataTransformation::usage =  "Transform the input data according to the input such that it can be used for simulations.";

Needs["DescriptionUtilities`"]

Begin["`Private`"];


inputDataTransformation[inpData_,prime_:True] := Module[{Kprime,Mprime,Eprime,Cprime,
Voint,tsint,Qoint,alpha,intn,dg},

(* Decide on Kprime *)
If[KeyExistsQ[inpData,KIc]&&prime,
Kprime = Sqrt[32/\[Pi]]KIc/.inpData//N,
If[KeyExistsQ[inpData,KIc],
Kprime = KIc/.inpData,
Kprime = Kp/.inpData]
];

(* Decide on n *)
If[KeyExistsQ[inpData,n],
alpha = (2^(n + 1) (2 n + 1)^n)/n^n;
intn = n/.inpData;,
alpha = 12;
intn = 1];

(* Decide on \[Mu]' *)
If[KeyExistsQ[inpData,\[Mu]],
Mprime = alpha \[Mu]/.inpData,
Mprime = \[Mu]p/.inpData];

(* Decide on E' *)
If[KeyExistsQ[inpData,Emod],
Eprime = Emod/(1-\[Nu]^2)/.inpData,
Eprime = Ep/.inpData];

(* Decide on C' *)
If[KeyExistsQ[inpData,Cl],
Cprime = 2 Cl/.inpData,
Cprime = Cp/.inpData];

(* Decide on Vo *)
If[KeyExistsQ[inpData,Vo],
Voint = Vo/.inpData,
Voint = Qo ts/.inpData];

(* Decide on Qo *)
Qoint = Qo/.inpData;

(* Decide on ts *)
tsint = ts/.inpData;

(* Decide on \[CapitalDelta]\[Gamma] *)
dg = \[CapitalDelta]\[Gamma]/.inpData;

{Ep-> Eprime, Vo -> Voint, Qo -> Qoint, ts -> tsint, Mp -> Mprime, Kp -> Kprime,
Cp-> Cprime,n -> intn, \[CapitalDelta]\[Gamma]-> dg}

];


End[]; (* End Private Context *)

EndPackage[];
