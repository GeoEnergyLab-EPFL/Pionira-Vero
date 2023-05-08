(* ::Package:: *)

(* Mathematica Package *)

BeginPackage["PostProcessRadialPL`"]
(* Exported symbols added here with SymbolName::usage *)  

Needs["Units`"];
 
Needs["DescriptionUtilities`"];

makePlotsHFResults::usage = "Function creating plots from the results structure of a HF"
backToDimensional::usage = "Back to dimensional results from dimensionless results structure of a HF"
makeTablesHFResults::usage = "Function creating Tables from the results structure of a HF"


Begin["`Private`"] (* Begin Private Context *) 

Protect[#]&/@ {Dimensionless,Animation}

Options[LegendsRulesHF]:={FontFamily-> "Times",  FontSize-> Larger}
LegendsRulesHF[Dimensionless_,opts: OptionsPattern[]]:=Module[{legends},
	(*     LEGENDS FOR HF VARIABLES , dimensional or dimensionless*)
	   If[Dimensionless,
   	(* Dimensionless *)
    legends = {
    	Time -> 
       Style[ "\[Tau]" , FilterRules[{opts},{FontFamily , FontSize }]], 
       FractureVelocity -> 
       Style[ "\[Gamma]" , FilterRules[{opts},{FontFamily , FontSize }]], 
      Opening -> 
       Style[ "\[CapitalOmega]" , FontFamily -> "Times" , 
        FilterRules[{opts},{FontFamily , FontSize }]], 
      WellborePressure -> 
       Style[ "\!\(\*SubscriptBox[\(\[CapitalPi]\), \(wb\)]\)" , 
       FilterRules[{opts},{FontFamily , FontSize }]], 
      EffectiveFlux -> 
       Style[ "\!\(\*SubscriptBox[\(\[CapitalPsi]\), \(eff\)]\)" , 
       FilterRules[{opts},{FontFamily , FontSize }]], 
      FluidPressure -> 
       Style[ "\!\(\*SubscriptBox[\(\[CapitalPi]\), \(f\)]\)"  , 
        FilterRules[{opts},{FontFamily , FontSize }]], 
        NetPressure -> Style[ "\[CapitalPi]"  , 
        FilterRules[{opts},{FontFamily , FontSize }]],        
      FractureLength -> 
       Style[ "\[Gamma]" ,FilterRules[{opts},{FontFamily , FontSize }]],
       FractureEfficiency -> 
       Style[ "\[Eta]" ,FilterRules[{opts},{FontFamily , FontSize }]]
       , 
       InletOpening -> Style[ "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(o\)]\)" , 
       FilterRules[{opts},{FontFamily , FontSize }]],
       InletPressure -> Style[ "\!\(\*SubscriptBox[\(\[CapitalPi]\), \(o\)]\)" , 
       FilterRules[{opts},{FontFamily , FontSize }]],
      X -> 
       Style[ "\[Xi]", FilterRules[{opts},{FontFamily , FontSize }]], 
      CohesiveZoneLength -> 
       Style[ "\!\(\*SubscriptBox[\(\[Gamma]\), \(\(coh\)\(.\)\)]\)", FilterRules[{opts},{FontFamily , FontSize }]], 
       CohesiveFraction -> Style[   "\!\(\*SubscriptBox[\(\[Xi]\), \(\(coh\)\(.\)\)]\)", FilterRules[{opts},{FontFamily , FontSize }]], 
      Energy -> 
       Style["E", FilterRules[{opts},{FontFamily , FontSize }]] 
       },       
       (*        Dimensional
        *)
    legends = {Time -> 
       Style[  "t" ,FilterRules[{opts},{FontFamily , FontSize }]],
       FractureVelocity -> 
       Style[ "l" , FilterRules[{opts},{FontFamily , FontSize }]],  
      Opening -> 
       Style[  "w" ,FilterRules[{opts},{FontFamily , FontSize }]], 
      WellborePressure -> 
       Style[  "\!\(\*SubscriptBox[\(p\), \(wb\)]\)" , 
        FilterRules[{opts},{FontFamily , FontSize }]], 
       EffectiveFlux -> 
       Style[ "\!\(\*SubscriptBox[\(Q\), \(eff\)]\)" , 
       FilterRules[{opts},{FontFamily , FontSize }]], 
     fluidpressure -> 
       Style[  "\!\(\*SubscriptBox[\(p\), \(f\)]\)", 
       	FilterRules[{opts},{FontFamily , FontSize }]] ,
       	netpressure ->  Style["p", FilterRules[{opts},{FontFamily , FontSize }]], 
      FractureLength -> 
       Style[ "l" ,  BaseStyle -> {FontFamily -> "Times", FontSize -> Larger} ,
        FilterRules[{opts},{FontFamily , FontSize }]],
       FractureEfficiency -> 
       Style[ "\[Eta]" ,FilterRules[{opts},{FontFamily , FontSize }]]
       ,
       InletOpening -> Style[ "\!\(\*SubscriptBox[\(w\), \(o\)]\)" , 
       FilterRules[{opts},{FontFamily , FontSize }]],
       InletPressure -> Style[ "\!\(\*SubscriptBox[\(p\), \(o\)]\)" , 
       	FilterRules[{opts},{FontFamily , FontSize }]],
      X -> Style[ "x",   FilterRules[{opts},{FontFamily , FontSize }]], 
      CohesiveZoneLength -> 
       Style[  "\!\(\*SubscriptBox[\(l\), \(\(coh\)\(.\)\)]\)", 
        FilterRules[{opts},{FontFamily , FontSize }]], 
      CohesiveFraction -> Style[   "\!\(\*SubscriptBox[\(\[Xi]\), \(\(coh\)\(.\)\)]\)", FilterRules[{opts},{FontFamily , FontSize }]],
      Energy -> 
       Style[ "E" ,FilterRules[{opts},{FontFamily , FontSize }]]}
    ]
	
];

Options[makePlotsHFResults] := {HistoryProfile -> True, 
   Dimensionless -> True, PlotStyle -> Automatic, Filling -> False, Animation-> True,FontFamily-> "Times",
   FontSize-> Larger,LogLog-> False,MovingMesh-> False,PlotMarkers-> None,PlotRange-> All};
makePlotsHFResults[results_, opts :   OptionsPattern[]] :=
  Module[{resultsField, resultsScal, legends, animationTime, 
    anOpening, anPressure, maxX, maxO, maxPf,ScalPlots,ScalFields,PressFields,i},
      
    
   legends=LegendsRulesHF[OptionValue[Dimensionless],FilterRules[{opts},{FontFamily , FontSize }]];
    
   If[OptionValue[HistoryProfile],
    {resultsField = results[[2]]; resultsScal = results[[1]];},
    resultsScal = results];

	(* get all the results rules *)

    ScalFields=Intersection[resultsScal[[;;,1]],
    	{FractureLength, FractureEfficiency, WellborePressure, EffectiveFlux,
   		InletOpening, InletPressure,CohesiveZoneLength,CohesiveFraction,EnergyInput,
   		CompressibilityEnergy, FractureEnergy, ViscousEnergy,FractureVelocity}];  
	 
	 
 
    (* SCALAR plots *)
    If[OptionValue[LogLog],
    ScalPlots = ( # -> ListLogLogPlot[# /. resultsScal, Joined -> True, 
     Frame -> True, 
     FrameLabel -> {Time /. legends, # /. legends}, 
       BaseStyle -> {FontFamily -> "Times", FontSize -> Larger},
       FilterRules[{opts}, Options[ListLogLogPlot]]      
     ] )  &  /@ Cases[ScalFields,Except[FractureEfficiency]];
     
     AppendTo[ScalPlots, FractureEfficiency-> ListLogLinearPlot[FractureEfficiency /. resultsScal, Joined -> True, 
     Frame -> True, 
     FrameLabel -> {Time /. legends, FractureEfficiency /. legends}, 
       BaseStyle -> {FontFamily -> "Times", FontSize -> Larger},
       FilterRules[{opts}, Options[ListLogLogPlot]]      
     ]
     ]
      	,
   ScalPlots = ( # -> ListPlot[# /. resultsScal, Joined -> True, 
     Frame -> True,  
     FrameLabel -> {Time /. legends, # /. legends}, 
       BaseStyle -> {FontFamily -> "Times", FontSize -> Larger}, 
     FilterRules[{opts}, Options[ListPlot]] 
     ] )  &  /@ ScalFields;
    ];
       
   	
   If[OptionValue[HistoryProfile],

	PressFields=Intersection[resultsField[[ 1, ;; , 1]] ,{NetPressure,FluidPressure}];

    (* CREATE ANIMATION *)
    animationTime = 
     Time /. resultsField[[#]] & /@ Range[Length[resultsField ]];
   	 
	If[OptionValue[MovingMesh], 
	maxX=1.05, 
    maxX = Max[Transpose[(FractureLength /. resultsScal)][[2]]] 1.05];
    
    maxO = Max[Transpose[Opening /. resultsField[[-1]]][[2]]] ;
    maxX = Max[Transpose[Opening /. resultsField[[-1]]][[1]]] 1.05 ;
  
   (* Print[maxX];*)
    If[maxO > 1, maxO = Ceiling[maxO]];
    If[Intersection[ScalFields,{WellborePressure}]!={},    	 
    	maxPf =Max[Transpose[(WellborePressure /. resultsScal)][[2]]] // Ceiling;,
    	maxPf =Max[Transpose[(InletPressure /. resultsScal)][[2]]] // Ceiling;
    ];
    
    If[OptionValue[Animation],
    {
    	anOpening =Opening ->   ListAnimate[ListPlot[Opening /. resultsField[[ #]], Joined -> True, 
        PlotMarkers -> None, PlotRange -> {{0., maxX}, {0., maxO}}, 
        FrameLabel -> {X /. legends, Opening /. legends}, 
        Frame -> True, FilterRules[{opts}, Options[ListPlot]]] & /@ 
      Range[Length[resultsField]] ,DefaultDuration-> 5
    (*animationTime[[-1]]*), AnimationRunning->True,Deployed-> True, BaseStyle -> {FontFamily -> "Times", FontSize -> Larger}] ;
    
     (* output either fluidpressure or netpressure, by default net pressure *)
    
    anPressure = Table[ PressFields[[i]]-> ListAnimate[ 
     ListPlot[PressFields[[i]] /. resultsField[[ #]], Joined -> True, 
        PlotMarkers -> None, 
        PlotRange -> {{0., maxX}, {-maxPf/2, maxPf}}, 
        FrameLabel -> {X /. legends, PressFields[[i]] /. legends}, 
        Frame -> True, FilterRules[{opts}, Options[ListPlot]]] & /@ 
      Range[Length[resultsField]],DefaultDuration-> 5
     (*animationTime[[-1]]*), AnimationRunning->True,Deployed-> True],{i,1,Length[PressFields]}] ;
     
    },
    {
    (* Series of LIST PLOT *)
     	anOpening =Opening ->
     ListPlot[(Opening /. resultsField[[#]]) & /@ 
      Range[Length[resultsField]], Joined -> True, 
        PlotMarkers -> None, PlotRange -> {{0., maxX}, {0., maxO}}, 
        FrameLabel -> {X /. legends, Opening /. legends},  
        Frame -> True, BaseStyle -> {FontFamily -> "Times", FontSize -> Larger} , FilterRules[{opts}, Options[ListPlot]] ](*,DefaultDuration->
    animationTime[[-1]], AnimationRunning->True,Deployed-> True]*);
    
     
    anPressure = Table[  PressFields[[i]]-> 
     ListPlot[( PressFields[[i]] /. resultsField[[ #]]) & /@ 
      Range[Length[resultsField]], Joined -> True, 
        PlotMarkers -> None, 
        PlotRange -> {{0., maxX}, {-maxPf/2, maxPf}}, 
        FrameLabel -> {X /. legends,  PressFields[[i]] /. legends}, 
        Frame -> True, BaseStyle -> {FontFamily -> "Times", FontSize -> Larger}, FilterRules[{opts}, Options[ListPlot]]],{i,1,Length[PressFields]}] ;

    }];
	{ScalPlots,anOpening,     anPressure} // Flatten 
    ,
    ScalPlots
    ]
   ];

(**************************----********************************)
Options[makeTablesHFResults]:= {HistoryProfile -> True, Background -> ColorData[18][7],Dimensionless -> True, FontSize-> Larger ,FontFamily-> "Times"}
makeTablesHFResults[results_,opts :OptionsPattern[]]:=Module[{legends,resultsField,resultsScal,ScalFields,PressFields,auxFrac,auxInlet,auxOpPro,auxPiPro,i,Nt},
	
	Off[Symbol::argx];
	Off[Intersection::heads];
	
   legends=LegendsRulesHF[OptionValue[Dimensionless],FilterRules[{opts},{FontFamily , FontSize }]];
 	   
	   (* Display results in Tables *)
 	If[OptionValue[HistoryProfile],     
    	{
    	resultsField = results[[2]]; 
    	resultsScal = results[[1]];},
    	resultsScal = results];
	(* get all the results rules *)
	 
    ScalFields=Intersection[resultsScal[[;;,1]],
   	{FractureLength,FractureVelocity,FractureEfficiency,WellborePressure,
   		InletOpening,InletPressure,CohesiveZoneLength,EnergyInput, CompressibilityEnergy, FractureEnergy, ViscousEnergy}];  
	 
	PressFields=Intersection[results[[2, 1, ;; , 1]] ,{NetPressure,FluidPressure}];

	(* BY CONVENTION WE MAKE 2 TABLES for Scalar Variables *)
	 
	(* one with Time,frac. length, frac efficiency *)
	
	auxFrac=Grid[
		Prepend[Transpose[{Time/.resultsScal,Transpose[FractureLength/.resultsScal][[2]],Transpose[FractureEfficiency/.resultsScal][[2]]}],{Time/.legends,FractureLength/.legends,FractureEfficiency/.legends}],
		Frame-> All,Background-> OptionValue[Background]
		];

	(* one with Time, inlet opg., inlet pressure *)
	
	auxInlet=Grid[
		Prepend[
			Transpose[{Time/.resultsScal,Flatten[Transpose[InletOpening/.resultsScal][[2]]],Flatten[Transpose[InletPressure/.resultsScal][[2]]]}] 
				,{Time/.legends,InletOpening/.legends,InletPressure/.legends}],
		Frame-> All,Background-> OptionValue[Background]
		];

   Nt=Length[resultsField];
  
	If[OptionValue[HistoryProfile],
		 
	auxOpPro=Table[ (ToString[(Time/.legends)] <> "=" <>
	    ToString[EngineeringForm[(Time/.resultsField[[i]]),3]] )->   
	Grid[Prepend[ (Opening/.resultsField[[i]]),{X/.legends,Opening/.legends}],
		Frame-> All,Background-> OptionValue[Background]
		] ,{i,1,Nt,2}] ;
		auxPiPro=
	  Table[ (ToString[(Time/.legends)] <> "=" <>
	    ToString[EngineeringForm[(Time/.resultsField[[i]]),3]] ) ->  
	Grid[Prepend[
		  (NetPressure/.resultsField[[i]])  
			,{X/.legends,NetPressure/.legends}],
		Frame-> All,Background-> OptionValue[Background]
		],{i,1,Nt,2}]  ;
		 
		{auxFrac,auxInlet,auxOpPro,auxPiPro}		
		,
		{auxFrac,auxInlet}		
	]
];



(**************************----********************************)
Options[backToDimensional] := {HistoryProfile -> True,MovingMesh-> False}
backToDimensional[results_, Scales_, Units_, OptionsPattern[]] := 
 Module[{aux,PiAux,OpAux,LAux,Lcoef,pcoef,wcoef,tcoef,resultsScal,ScalFields,i,qcoef,QeffAux,Velaux,Effaux},
  (*note this function can also be used  to express to switch the results back to another scaling if need be. *)
  (* Scales is a replacement rule of the form tstar-> , Lstar->  etc. defining the scales in the scaling in which  
  the 'dimensionless results' have been computed *)
  
  tcoef = tstar /. Scales /. Units; 
  wcoef = wstar /. Scales /. Units;
  Lcoef = Lstar /. Scales /. Units;
  pcoef = pstar /. Scales /. Units;
  qcoef = qstar /.Scales /. Units;
  
  resultsScal = results[[1]];
  ScalFields = resultsScal[[;;,1]];

  Intersection[{WellborePressure,InletPressure},ScalFields];

  PiAux =  ( # ->  ( Transpose[ {tcoef,pcoef} Transpose[# /.  resultsScal] ])  ) & /@ Intersection[{WellborePressure,InletPressure},ScalFields];
  OpAux  = ( # ->  ( Transpose[ {tcoef,wcoef} Transpose[# /.  resultsScal] ])  ) & /@ Intersection[{InletOpening},ScalFields];   
  LAux  = ( # ->  ( Transpose[ {tcoef,Lcoef} Transpose[# /.  resultsScal] ])  ) & /@  Intersection[{FractureLength,CohesiveZoneLength},ScalFields];   
  QeffAux = ( # ->  ( Transpose[ {tcoef,qcoef} Transpose[# /.  resultsScal] ])  ) & /@  Intersection[{EffectiveFlux},ScalFields];
  Velaux = (# -> (Transpose[{tcoef,Lcoef/tcoef} Transpose[# /.resultsScal] ])  )& /@ Intersection[{FractureVelocity},ScalFields];
  Effaux = (# -> (Transpose[{tcoef,1} Transpose[# /.resultsScal] ])  )& /@ Intersection[{FractureEfficiency},ScalFields];
  
  
	 aux= { (Time -> tcoef (Time /.resultsScal) ) ,PiAux, OpAux, LAux,QeffAux,Velaux,Effaux} //Flatten;
	  

  If[OptionValue[MovingMesh],Lcoef=1]; (* keep the scale w.r to fracture length for the absciss of the profile in that case *)
  If[OptionValue[HistoryProfile],
    {aux,
    Table[
     {Time -> tcoef ( Time /. results[[2, i]] ) ,
      Opening -> ({Lcoef,wcoef} # & /@  (Opening /. 
           results[[2, i]])),
      NetPressure -> ({Lcoef,pcoef} # & /@  (NetPressure /. 
           results[[2, i]]))
           (*,
      Fluidpressure ->  ({Lcoef,pcoef} # & /@  (Fluidpressure /. 
           results[[2, i]])),
       (* if exist *)    
      ClampingStress -> ({Lcoef,pcoef} # & /@  (ClampingStress /. 
           results[[2, i]])),
      CohesiveStress -> ({Lcoef,pcoef} # & /@  (CohesiveStress /. 
           results[[2, i]]))
           *)
           }, {i, 1, Length[results[[2]]]}]
    },
   aux]
  ]

End[] (* End Private Context *)

EndPackage[]
