(* ::Package:: *)

(* Mathematica Package for KGD SOLUTIONS IMPERMEABLE MEDIUM - zero leak off*)

BeginPackage["KGDSolutions`"]

Needs["ProtectedDefinitions`"];

(* Exported symbols added here with SymbolName::usage *)  

(*Needs["DescriptionUtilities`"];*)
KGDScaling::usage = "   viscosity or toughness scaling " 

ScalingTranslation::usage = "  " 
ZeroToughnessSolution::usage = "KGD self-similar Solution of a HF propagating under zero Toughness (M-vertex) M-scaling ";
ZeroToughnessSolution3Terms::usage = "KGD self-similar Solution of a HF propagating under zero Toughness (M-vertex) M-scaling - 3 firt terms only "
SmallToughnessSolution::usage = "KGD self-similar solution of a HF propagating with a small toughness (near M-vertex)  M-scaling, valid up to dimensionless toughness of 1";
ZeroViscositySolution::usage = "KGD self-similar solution of a HF propagating with a Large toughness (zero viscosity) (K-vertex) K-scaling";
SmallViscositySolution::usage = "KGD self-similar solution of a HF propagating with a small viscosity (near K-vertex)   K-scaling, 
									valid up to dimensionless viscosity of 0.01 (for opg+pressre), 0.1 for frac length";

SmallViscositySolution1::usage = " KGD self-similar solution of a HF propagating with a small viscosity (near K-vertex)   K-scaling,  linear expansion" 
  


Begin["`Private`"] (* Begin Private Context *) 


Protect[#] & /@ {Opening,  NetPressure, FractureLength} ;

(*Protect[#] & /@ {mup,Ep,Kp, Qo};*)

Protect[#] & /@ {Kbar,Mbar,gamma};

Protect[#] & /@ {wstar,pstar,Lstar,qstar};

(*----------------------------------------------------------------------------------*) 
(*         Scaling      *)
(*Needs to put Scaling as an options in as per the KGDWellboreScaling packages*)
KGDScaling["Viscosity"][Ep_,Qo_,mup_,Kp_,t_] =
{Scales-> {Lstar -> (Ep Qo^3 t^4/mup)^(1/6), wstar -> (mup/(Ep t))^ (1/3) (Ep Qo^3 t^4/mup)^(1/6)  ,pstar -> Ep (mup/(Ep t))^ (1/3) , qstar -> Qo },
    DimensionlessNumbers-> {Kbar-> Kp / (Ep^3 mup Qo )^(1/4) (*,Rbar-> 1,Ebar-> 1*), Mbar-> 1}  };


KGDScaling["Toughness"][Ep_,Qo_,mup_,Kp_,t_] =
{Scales-> {Lstar -> (Ep Qo  t/  ( Kp) )^(2/3), wstar -> (( Kp)^4/(Ep^4 Qo t))^(1/3) (Ep Qo  t/  ( Kp) )^(2/3)   ,pstar -> Ep (( Kp)^4/(Ep^4 Qo t))^(1/3) , qstar -> Qo },
    DimensionlessNumbers-> {Kbar-> 1 (*,Rbar-> 1,Ebar-> 1*), Mbar->  mup Ep^3 Qo /  (( Kp)^4) }  };


Options[ScalingTranslation] = {ToScaling -> "Toughness", DimensionlessNumber -> 0.1 };
ScalingTranslation[ opts: OptionsPattern[]] :=
    Module[ {M,K,rep},
        If[ OptionValue[ToScaling]== "Toughness",
            {M = OptionValue[DimensionlessNumber];
             K = M^-4;
             rep = {Kbar->K, gfact-> M^(-8/3), omfact -> M^8, pfact-> M^(16/3) , Lfact-> M^(8/3) };
                        
            },
            If[ OptionValue[ToScaling]=="Viscosity",
                {K = OptionValue[DimensionlessNumber];
                 M = K^(-1/4);
                 rep = {Mbar-> M,gfact->K^-(2/3) , omfact -> K^2,  pfact ->  K^4/3, Lfact-> K^(2/3)};
                
                },
                {Print["error in option ToScaling"];
                 Abort[]}
            ]
        ];
        rep
    ] 


(* --------------------------------------------------------------------------------*) 
(* BE CAREFUL ALL OPENINGS are Omegabar = Omega/gamma *)

(* Zero TOUGHNESS *)
Ombstst[xi_]=4Sqrt[1-xi^2]+2 xi^2 Log[Abs[(1-Sqrt[1-xi^2])/(1+Sqrt[1-xi^2])]];
Ombst[xi_,j_]=(1-xi^2)^(2/3) GegenbauerC[2j-2,2/3-1/2,xi];

Pistst[xi_]=2 - \[Pi] Abs[xi];
Pist1[xi_]=\[Alpha]/(2 \[Pi]) Beta[1/2,\[Alpha]]Hypergeometric2F1[1/2-\[Alpha], 1,1/2, xi^2]   /.{\[Alpha]->2/3};
Pist[xi_,j_]=((2\[Alpha]-1)(2 j-1))/(4\[Pi] (j-1+\[Alpha])) Beta[1/2-j,\[Alpha]+j](\[Alpha] xi^2 Hypergeometric2F1[5/2-j-\[Alpha], j,3/2, xi^2]  -1/2 Hypergeometric2F1[3/2-j-\[Alpha], j-1,1/2, xi^2])/.{\[Alpha]->2/3};
(* Series coefficients from *) 
AjsNewtonian={1.61750,0.39650,1.00297 10^-2,1.06179 10^-2,-2.97947 10^-3,-0.74268  10^-3,-1.34675  10^-3};
BjsNewtonian = 0.06858;

Ombar0[xi_]= BjsNewtonian Ombstst[xi] +Sum[AjsNewtonian[[j]] Ombst[xi,j],{j,1,7}];
Pim0[xi_]=BjsNewtonian Pistst[xi] + AjsNewtonian[[1]] Pist1[xi] +Sum[AjsNewtonian[[j]] Pist[xi,j],{j,2,7}] ;
(* Normal Stress ahead ! *)
Pist1ahead[xi_]=-2/(20Sqrt[\[Pi]]) Gamma[8/3]/Gamma[13/6] Hypergeometric2F1[2/3,7/6,13/6,1/ xi^2] / (Abs[xi]^(4/3) (xi^2-1)^(1/3));
Pistahead[xi_,j_]=Gamma[2j-5/3]/(Gamma[1/6]Gamma[2j-5/6]) ((6j-5)(xi^2-1)Hypergeometric2F1[j-1/3,j+1/6,2j-5/6,1/ xi^2]-4xi^2 Hypergeometric2F1[j-5/6,j-1/3,2j-5/6,1/ xi^2])/(Abs[2 xi]^(2j-2/3) (xi^2-1)^(1/3));
Piststahead[xi_]=2-2 Abs[xi]ArcCot[Sqrt[xi^2-1]];

NetStressAhead[xi_]= BjsNewtonian Piststahead[xi]+ AjsNewtonian[[1]] Pist1ahead[xi]+ Sum[AjsNewtonian[[j]] Pistahead[xi,j],{j,2,7}] ;


ZeroToughnessSolution[xi_] := {gamma -> 0.61524, 
    Omegabar -> If[Abs[xi]<1,Ombar0[xi],0],Opening ->  If[Abs[xi]<1,0.61524 Ombar0[xi] ,0],
    NetPressure -> If[Abs[xi]<1,Pim0[xi],NetStressAhead[xi]]
     };

(* in viscosity scaling  3 terms only in the series *)
ZeroToughnessSolution3Terms[xi_] :=
    {gamma -> 0.61524, 
    Omegabar -> 
    a01 (1 - xi^2)^(2/3) + a02 (1 - xi^2)^(5/3) + 
    a03 (2. Sqrt[1 - xi^2] + 
    1. xi^2 Log[    Abs[(1 - Sqrt[1 - xi^2])/(1 + Sqrt[1 - xi^2])]]), 
    Opening -> 0.61524(
    a01 (1 - xi^2)^(2/3) + a02 (1 - xi^2)^(5/3) + 
    a03 (2. Sqrt[1 - xi^2] + 
    1. xi^2 Log[    Abs[(1 - Sqrt[1 - xi^2])/(1 + Sqrt[1 - xi^2])]])),
    
    NetPressure -> 
    b01 Hypergeometric2F1[-1/6, 1., 1/2., xi^2] + 
    b02   Hypergeometric2F1[-7./6, 1, 1/2., xi^2]  + 
    b03 (2 - \[Pi] Abs[xi])} /. {a01 -> 3.^(1/2), a02 -> -0.15601,a03 -> 0.13264, b01 -> 0.475449, b02 -> -0.061178, b03 -> 0.066322};


(* SMALL TOUGHNESS *)
FirstOrderCorrectionSmallToughness[xi_] = {gamma -> - 0.17475, 
    Omegabar -> 
     a11 (1 - xi^2)^h + (a12 + a13 xi^2) (1 - xi^2)^h + 
      a14 (2 (1 - xi^2)^(1/2) + xi^2 Log[Abs[(1 - Sqrt[1 - xi^2])/(1 + Sqrt[1 - xi^2])]])
    , 
    NetPressure -> 
     b11  Hypergeometric2F1[c11, 1., 1/2., xi^2] + 
      b12   Hypergeometric2F1[c12, 1, 1/2., xi^2]   + 
      b13   Hypergeometric2F1[c13, 2, 3/2., xi^2] xi^2 + 
      b14  (2 - \[Pi] Abs[xi])} /. {h -> 0.13867, a11 -> 0.908354,a12 -> 0.025574, a13 -> -0.083814, a14 -> -0.09095,b11 -> 0.170654, 
    b12 -> 0.017132, b13 -> -0.039015, b14 -> -0.045476, c11 -> 0.36133, c12 -> -1.63867,c13-> -0.638673};

SmallToughnessSolution[xi_, K_] :=
    Module[ {epsilon, B  , b  , F0, F1,aux},
          B = 0.1076;b = 3.16796;
         epsilon = B K^b;         
         F0 = ZeroToughnessSolution[xi];
         F1 = FirstOrderCorrectionSmallToughness[xi];
        aux = ( # -> (# /. F0) + epsilon (# /. F1)  )& /@ {gamma, NetPressure, Omegabar } ;
        
        Join[aux,{Opening-> (Omegabar gamma) /.aux} ]    
    ]
   
      
(*   --------------------------------------------------------------------------------*) 
(* LARGE TOUGHNESS *)
(* in toughness scaling from Garagash EFM 2006 *)

ZeroViscositySolution[xi_] := {gamma -> 2./\[Pi]^(2/3),    NetPressure -> Pi^(1/3)/8.   + 0 xi, 
    Omegabar ->  Pi^(1/3)/2. Sqrt[1 - xi^2] , Opening ->2./\[Pi]^(2/3) Pi^(1/3)/2. Sqrt[1 - xi^2] }; 
 
 
FirstOrderCorrectionSmallViscosity[xi_] :=    {gamma ->  -32 (1. +    6 Log[2])/(9 Pi^(5/3)) , 
    NetPressure ->     8/(3. Pi^(
    2/3)) (1/24 + Log[4 Sqrt[1 - xi^2]] - 
    3/4 ( xi ArcCos[xi])/ ((1 - xi^2)^(1/2.))), 
    Omegabar ->  
    8/(3 Pi^(2/3.)) (2 Pi - 4 xi ArcSin[xi] - (5/6 - Log[2]) Sqrt[1 - xi^2] - 
    3/2 Log[(1 + Sqrt[1 - xi^2])^(1 + Sqrt[1 - xi^2])/(1 - Sqrt[1 - xi^2])^(
    1 - Sqrt[1 - xi^2])])   };

         
SmallViscositySolution[xi_, M_] :=
    Module[ {delta, Msi = 0.0333, F0, F1,aux},
         delta = M/Sqrt[1 + M/Msi];
         F0 = ZeroViscositySolution[xi];
         F1 = FirstOrderCorrectionSmallViscosity[xi];
        aux =  # -> ((# /. F0) + delta (# /. F1) ) & /@ {gamma, NetPressure, Omegabar  }   ;
           Join[aux,{Opening-> (Omegabar gamma) /.aux} ]
    ]         

SmallViscositySolution1[xi_, M_] :=
    Module[ {delta, F0, F1,aux},
         delta = M;
         F0 = ZeroViscositySolution[xi];
         F1 = FirstOrderCorrectionSmallViscosity[xi];
       aux = (  # -> (# /. F0) + delta (# /. F1) ) & /@ {gamma, NetPressure, Opening  };
              Join[aux,{Opening-> (Omegabar gamma) /.aux} ]    
    ]

(*	--------------------------------------------------------------------------------	*)  


End[] (* End Private Context *)

EndPackage[]
