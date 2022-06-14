(* Mathematica Package for KGD SOLUTIONS
 Newtonian Fracturing fluid
- zero leak off*)

BeginPackage["KGDSolutions`"]

Needs["ProtectedDefinitions`"];

KGDScaling::usage = "   viscosity or toughness scaling "

ScalingTranslation::usage = "  " 
StorageZeroToughnessSolution::usage = "KGD self-similar Solution of a HF propagating under zero Toughness (M-vertex) M-scaling ";
StorageSmallToughnessSolution::usage = "KGD self-similar solution of a HF propagating with a small toughness (near M-vertex)  M-scaling, valid up to dimensionless toughness of 1";
StorageZeroViscositySolution::usage = "KGD self-similar solution of a HF propagating with a Large toughness (zero viscosity) (K-vertex) K-scaling";
StorageSmallViscositySolution::usage = "KGD self-similar solution of a HF propagating with a small viscosity (near K-vertex)   K-scaling,
									valid up to dimensionless viscosity of 0.01 (for opg+pressre), 0.1 for frac length";
// StorageSmallViscositySolution1::usage = " KGD self-similar solution of a HF propagating with a small viscosity (near K-vertex)   K-scaling,  linear expansion"

Begin["`Private`"] (* Begin Private Context *)

(*----------------------------------------------------------------------------------*) 
(*         Scalings      *)

KGDScaling["ViscosityStorage"][Ep_,Qo_,mup_,Kp_,t_] =
{Scales-> {Lstar -> (Ep Qo^3 t^4/mup)^(1/6), wstar -> (mup/(Ep t))^ (1/3) (Ep Qo^3 t^4/mup)^(1/6)  ,pstar -> Ep (mup/(Ep t))^ (1/3) , qstar -> Qo },
    DimensionlessNumbers-> {Kbar-> Kp / (Ep^3 mup Qo )^(1/4) (*,Rbar-> 1,Ebar-> 1*), Mbar-> 1}  };

KGDScaling["ToughnessStorage"][Ep_,Qo_,mup_,Kp_,t_] =
{Scales-> {Lstar -> (Ep Qo  t/  Kp )^(2/3), wstar -> (( Kp)^4/(Ep^4 Qo t))^(1/3) (Ep Qo  t/  ( Kp) )^(2/3),
pstar -> Ep (( Kp)^4/(Ep^4 Qo t))^(1/3) , qstar -> Qo },
    DimensionlessNumbers-> {Kbar-> 1 , Mbar->  mup Ep^3 Qo /  (( Kp)^4) }  };


Options[ScalingTranslation] = {ToScaling -> "ToughnessStorage", DimensionlessNumber -> 0.1 };
ScalingTranslation[ opts: OptionsPattern[]] :=
    Module[ {M,K,rep},
        If[ OptionValue[ToScaling]== "ToughnessStorage",
            {M = OptionValue[DimensionlessNumber];
             K = M^-4;
             rep = {Kbar->K, gfact-> M^(-8/3), omfact -> M^8, pfact-> M^(16/3) , Lfact-> M^(8/3) };
            },
            If[ OptionValue[ToScaling]=="ViscosityStorage",
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

(* STORAGE DOMINATED SOLUTIONS *)
(* Omegabar= Omega / gamma *)
(* --------------------------------------------------------------------------------*) 
(* SMALL TOUGHNESS *) (* in viscosity scaling *)

(* Taken from Garagash & Detournay - JAM (2005) - appendix D*)
StorageZeroToughnessSolution[xi_] :=
    {gamma -> 0.61524, 
    Omegabar -> 
    a01 (1 - xi^2)^(2/3) + a02 (1 - xi^2)^(5/3) + 
    a03 (2. Sqrt[1 - xi^2] + 
    1. xi^2 Log[    Abs[(1 - Sqrt[1 - xi^2])/(1 + Sqrt[1 - xi^2])]]), 
    Opening -> 0.61524( a01 (1 - xi^2)^(2/3) + a02 (1 - xi^2)^(5/3) + a03 (2. Sqrt[1 - xi^2] + 1. xi^2 Log[    Abs[(1 - Sqrt[1 - xi^2])/(1 + Sqrt[1 - xi^2])]])),
    NetPressure -> 
    b01 Hypergeometric2F1[-1/6, 1., 1/2., xi^2] + 
    b02   Hypergeometric2F1[-7./6, 1, 1/2., xi^2]  + 
    b03 (2 - \[Pi] Abs[xi])} /. {a01 -> 3.^(1/2), a02 -> -0.15601,a03 -> 0.13264, b01 -> 0.475449, b02 -> -0.061178, b03 -> 0.066322};

StorageFirstOrderCorrectionSmallToughness[xi_] = {gamma -> - 0.17475,
    Omegabar -> 
     a11 (1 - xi^2)^h + (a12 + a13 xi^2) (1 - xi^2)^h + 
      a14 (2 (1 - xi^2)^(1/2) + xi^2 Log[Abs[(1 - Sqrt[1 - xi^2])/(1 + Sqrt[1 - xi^2])]])
    , NetPressure ->
     b11  Hypergeometric2F1[c11, 1., 1/2., xi^2] + 
      b12   Hypergeometric2F1[c12, 1, 1/2., xi^2]   + 
      b13   Hypergeometric2F1[c13, 2, 3/2., xi^2] xi^2 + 
      b14  (2 - \[Pi] Abs[xi])} /. {h -> 0.13867, a11 -> 0.908354,a12 -> 0.025574, a13 -> -0.083814, a14 -> -0.09095,b11 -> 0.170654, b12 -> 0.017132, b13 -> -0.039015, b14 -> -0.045476, c11 -> 0.36133, c12 -> -1.63867,c13-> -0.638673};

SmallToughnessSolution[xi_, K_] :=
    Module[ {epsilon, B  , b  , F0, F1,aux},
          B = 0.1076;
          b = 3.16796;
         epsilon = B K^b;         
         F0 = StorageZeroToughnessSolution[xi];
         F1 = StorageFirstOrderCorrectionSmallToughness[xi];
         aux = ( # -> (# /. F0) + epsilon (# /. F1)  )& /@ {gamma, NetPressure, Omegabar } ;
        Join[aux,{Opening-> (Omegabar gamma) /.aux} ]
    ]


(*   --------------------------------------------------------------------------------*) 
(* LARGE TOUGHNESS *) (* in toughness scaling from Garagash EFM 2006 *)

StorageZeroViscositySolution[xi_] := {gamma -> 2./\[Pi]^(2/3),    NetPressure -> Pi^(1/3)/8.   + 0 xi,
    Omegabar ->  Pi^(1/3)/2. Sqrt[1 - xi^2] , Opening ->2./\[Pi]^(2/3) Pi^(1/3)/2. Sqrt[1 - xi^2] }; 

StorageFirstOrderCorrectionSmallViscosity[xi_] :=    {gamma ->  -32 (1. + 6 Log[2])/(9 Pi^(5/3)) ,
    NetPressure ->     8/(3. Pi^(2/3)) (1/24 + Log[4 Sqrt[1 - xi^2]] -     3/4 ( xi ArcCos[xi])/ ((1 - xi^2)^(1/2.))),
    Omegabar ->  
    8/(3 Pi^(2/3.)) (2 Pi - 4 xi ArcSin[xi] - (5/6 - Log[2]) Sqrt[1 - xi^2] - 
    3/2 Log[(1 + Sqrt[1 - xi^2])^(1 + Sqrt[1 - xi^2])/(1 - Sqrt[1 - xi^2])^(1 - Sqrt[1 - xi^2])])   };

StorageSmallViscositySolution[xi_, M_] :=
    Module[ {delta, Msi = 0.0333, F0, F1,aux},
         delta = M/Sqrt[1 + M/Msi];
         F0 = StorageZeroViscositySolution[xi];
         F1 = StorageFirstOrderCorrectionSmallViscosity[xi];
        aux =  # -> ((# /. F0) + delta (# /. F1) ) & /@ {gamma, NetPressure, Omegabar  }   ;
           Join[aux,{Opening-> (Omegabar gamma) /.aux} ]
    ]         

StorageSmallViscositySolution1[xi_, M_] :=
    Module[ {delta, F0, F1,aux},
         delta = M;
         F0 = StorageZeroViscositySolution[xi];
         F1 = StorageFirstOrderCorrectionSmallViscosity[xi];
       aux = (  # -> (# /. F0) + delta (# /. F1) ) & /@ {gamma, NetPressure, Opening  };
              Join[aux,{Opening-> (Omegabar gamma) /.aux} ]    
    ]
(*	--------------------------------------------------------------------------------	*)
(* Lag / zero leak-off solution(s)...*)
(* early time Garagash 2006 solution *)

(* Lecampion - Detournay 2007 num. solution *)

(*	--------------------------------------------------------------------------------	*)
(* Leak-off solutions  *)

KGDScaling["ToughnessStorageLeakOffEdge"][Ep_,Qo_,mup_,Kp_,Cp_,t_] =
{Scales-> {Lstar ->  Kp^2 Qo^2 /(Ep^2 Cp^4), wstar -> (Cp^2/Qo) * Kp^2 Qo^2 /(Ep^2 Cp^4) ,
pstar ->  Ep * Cp^2/Qo  , qstar -> Qo ,tstar->Kp^4 Qo^2/ (Ep^4 Cp^6)  }, DimensionlessNumbers-> {Kbar-> 1 , Mbar->  mup Ep^3 Qo /  (Kp^4) }, Cbar-> 1  };

(* Large Toughness / Large Leak-off solution (Ktilde vertex ) - Bunger et al IJF (2006)*)
(* along the KKtilde edge - solution for large time / large leak-off *)
(* valid for \tau > 0.1 in the time-based edge scaling *)
auxKKtilde={ \[Alpha]-> -1/4, \[Beta]->1/2, gkt0->0.3183,gkt1-> -5.706*0.01,gkt2-> 1.341*0.01,gkt3-> -3.131*0.001,gkt4-> 6.368*0.0001}
solKKtildeEdge[tau_]= tau^(\[Beta]) * ( gkt0 + gkt1 * tau^(\[Alpha]) +gkt2 * tau^(\[Alpha]*2) + gkt3 * tau^(\[Alpha]*3) + gkt4 * tau^(\[Alpha]*4)    ) /.auxKKtilde
SmallViscosityLeakOffSolution[xi_,tau_]:={gamma-> solKKtildeEdge[tau],
Opening-> 2^(-1/2) (1-xi^2)^(1/2) (solKKtildeEdge[tau])^(1/2), NetPressure-> 2^(-5/2) (solKKtildeEdge[tau])^(-1/2)
}

(* viscosity / leak-off dominated   Adachi & Detournay (2008) *)
(* to be implemented *)





End[] (* End Private Context *)

EndPackage[]
