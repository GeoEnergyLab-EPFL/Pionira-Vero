


BeginPackage["ProtectedDefinitions`"]

omfact::usage = " factor of correspondance between 2 scalings for opening "
pfact::usage = " factor of correspondance between 2 scalings for net pressure "
Lfact::usage = " factor of correspondance between 2 scalings for length "
gfact::usage = " factor of correspondance between 2 scalings for dim groups "

Opening::usage = " Symbol for fracture opening "
NetPressure::usage ="Symbol reserved for net pressure in fracture"

Omegabar::usage = " Symbol reserved for the ratio between dimensionless opening and dimensionless length"
gamma::usage = "Symbol for dimensionless fracture length "

Scales::usage = " Symbol for Scales"
DimensionlessNumbers::usage = " Symbol for dimensionless numbers "

Qo::usage =" Symbol reserved for the (total) injected flux "
Ep::usage = " Symbole reserved for the plane strain Young's modulus Ep=E/(1-nu^2) "
Kp::usage = "Symbol reserved for the equivalent toughness Kp = Sqrt[32/Pi] Kic "
mup::usage = "Symbol reserved for the equivalent viscosity 12 muf "
Kic::usage = "Symbol reserved for the mode I fracture toughness "

(* Charateristic scales *)
wstar::usage = " Symbol reserved for the opening scale "
pstar::usage = " Symbol reserved for the pressure scale "
Lstar::usage = " Symbol reserved for the length scale "
qstar::usage= " Symbol reserved for the flux scale "
tstar::usage = "Symbol reserved for a time scale"

(* dimensionless numbers & symbol*)
Kbar::usage = " Symbol reserved for dimensionless toughness"
Mbar::usage = " Symbol reserved for dimensionless viscosity"
Cbar::usage = " Symbol reserved for dimensionless leak-off"
Sbar::usage = " Symbol reserved for dimensionless storage"
Obar::usage = " Symbol reserved for dimensionless lag"


Begin["`Private`"] (* Begin Private Context *)

Protect[#] & /@ {Opening,  NetPressure, FractureLength} ;

Protect[#] & /@ {mup,Ep,Kp, Qo};

Protect[#] & /@ {Kbar,Mbar,gamma};

Protect[#] & /@ {wstar,pstar,Lstar,qstar};




(*	--------------------------------------------------------------------------------	*)


End[] (* End Private Context *)

EndPackage[]
