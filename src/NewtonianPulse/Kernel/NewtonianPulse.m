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
BeginPackage["NewtonianPulse`"];

pulseVertexScalings::usage =  "Gives the vertex scales in the case of a pulse injection.
Input sequence is regime as string \[Element] [M, K, Mt, Kt, Mhh, Mht, MKh, MKt, Khh, Kht], the input data containing the information,
time (either a list or a scalar), and boolean deciding if Kp is used (default else set to false).";

pulseVertexSolutions::usage = "Gives the vertex solution of the \!\(\*SuperscriptBox[\(M\), \([V]\)]\)-vertex and \!\(\*SuperscriptBox[\(K\), \([V]\)]\)-vertex (zero buoyancy \[CapitalDelta]\[Gamma]).
Input sequence is regime as string \[Element] [M, K], the input data containing the information, value of \[Rho] dimensionless radius,
time (either a list or a scalar), boolean saying if the dimensional solution is wished (default else set to false), and 
boolean deciding if Kp is used (default else set to false).";

arrestRadius::usage = "Gives the solutions for the arrest radius (zero buoyancy \[CapitalDelta]\[Gamma]).
Input sequence is regime as string \[Element] [K, Kt, Mt] (first is toughness arrest radius, second and third viscosity),
the input data containing the information, boolean saying if the dimensional solution is wished (default else set to false), and 
boolean deciding if Kp is used (default else set to false)."

pulseTimeScales::usage = "Gives the transition timescales between the solutions. Input sequence is inputData (array of properties),
a boolean deciding if Kp is used (default else set to false), and a boolean deciding if the dimensional form (default) or dimensionless
form (set it to true). Output sequence is {phi,tmk,tmmt,tmtkt,tkkt,tmmh,tmhkh,tmhkh (head)}."

pulseTransitionMKScales::usage = "Gives the time-independent characteristic scales - call sequence is pulseTransitionMKScales[inpData,prime] with inpData an Array containg the material properties
and prime a boolean (true if Kp is used [default], false to use KIc)"

pulseTransitionMMhScales::usage =	"Gives the time-independent characteristic scales - call sequence is transitionMMhScales[inpData,prime] with inpData an Array containg the material properties
and prime a boolean (true if Kp is used [default], false to use KIc)"

pulseTransitionMhKhTailScales::usage =	"Gives the time-independent characteristic scales - call sequence is pulseTransitionMhKhTailScales[inpData,prime] with inpData an Array containg the material properties
and prime a boolean (true if Kp is used [default], false to use KIc)"

pulseTransitionMhKhHeadScales::usage =	"Gives the time-independent characteristic scales - call sequence is pulseTransitionMhKhTailScales[inpData,prime] with inpData an Array containg the material properties
and prime a boolean (true if Kp is used [default], false to use KIc)"

Needs["DescriptionUtilities`"]
Needs["UtilityForScalings`"]

Begin["`Private`"];


pulseVertexScalings[V_?StringQ,inpData_,t_,prime_:True] :=Module[{data},

data = inputDataTransformation[inpData,prime];

Switch[V,
	"M",
			{Lstar->(Ep^(1/9) t^(1/9) Vo^(1/3))/Mp^(1/9),
			wstar->(Mp^(2/9) Vo^(1/3))/(Ep^(2/9) t^(2/9)),
			pstar->(Ep^(2/3) Mp^(1/3))/t^(1/3),
			\[ScriptCapitalK] -> (Kp t^(5/18))/(Ep^(13/18) Mp^(5/18) Vo^(1/6)),
			\[ScriptCapitalC] -> (Cp Ep^(2/9) t^(13/18))/(Mp^(2/9) Vo^(1/3)),
			\[ScriptCapitalB] -> (t^(4/9) Vo^(1/3) \[CapitalDelta]\[Gamma])/(Ep^(5/9) Mp^(4/9))
			} /.data//N,
	"K",
			{Lstar->(Ep^(2/5) Vo^(2/5))/Kp^(2/5),
			wstar->(Kp^(4/5) Vo^(1/5))/Ep^(4/5),
			pstar->Kp^(6/5)/(Ep^(1/5) Vo^(1/5)),
			\[ScriptCapitalM] -> (Ep^(13/5) Mp Vo^(3/5))/(Kp^(18/5) t),
			\[ScriptCapitalC] -> (Cp Ep^(4/5) Sqrt[t])/(Kp^(4/5) Vo^(1/5)),
			\[ScriptCapitalB] -> (Ep^(3/5) Vo^(3/5) \[CapitalDelta]\[Gamma])/Kp^(8/5)
			} /.data//N,
	"Kt",
			{Lstar->Sqrt[Vo]/(Sqrt[Cp] t^(1/4)),
			wstar->(Kp Vo^(1/4))/(Cp^(1/4) Ep t^(1/8)),
			pstar->(Cp^(1/4) Kp t^(1/8))/Vo^(1/4),
			\[ScriptCapitalM] -> (Sqrt[Cp] Ep^3 Mp Sqrt[Vo])/(Kp^4 t^(3/4)),
			\[ScriptCapitalV] -> (Cp^(5/4) Ep t^(5/8))/(Kp Vo^(1/4))
			} /.data//N,
	"Mt",
			{Lstar->Sqrt[Vo]/(Sqrt[Cp] t^(1/4)),
			wstar->(Mp^(1/4) Vo^(3/8))/(Cp^(1/8) Ep^(1/4) t^(5/16)),
			pstar->(Cp^(3/8) Ep^(3/4) Mp^(1/4))/(t^(1/16) Vo^(1/8)),
			\[ScriptCapitalK] -> (Kp t^(3/16))/(Cp^(1/8) Ep^(3/4) Mp^(1/4) Vo^(1/8)),
			\[ScriptCapitalV] -> (Cp^(9/8) Ep^(1/4) t^(13/16))/(Mp^(1/4) Vo^(3/8))
			} /.data//N,
	"Mhh",
			{Lstarv->(Ep^(11/24) Vo^(1/8) Mp^(1/6))/(t^(1/6) \[CapitalDelta]\[Gamma]^(5/8)),
			Lstarh->(Ep^(11/24) Vo^(1/8) Mp^(1/6))/(t^(1/6) \[CapitalDelta]\[Gamma]^(5/8)),
			wstar->(Vo^(1/4) Mp^(1/3))/(Ep^(1/12) t^(1/3) \[CapitalDelta]\[Gamma]^(1/4)),
			pstar->(Ep^(11/24) Vo^(1/8) \[CapitalDelta]\[Gamma]^(3/8) Mp^(1/6))/t^(1/6),
			\[ScriptCapitalK] -> (Kp t^(1/4))/(Ep^(11/16) Vo^(3/16) \[CapitalDelta]\[Gamma]^(1/16) Mp^(1/4)),
			Vh -> (Ep^(5/6) Vo^(1/2) Mp^(2/3))/(t^(2/3) \[CapitalDelta]\[Gamma]^(3/2))
			} /.data//N,
	"Mht",
			{Lstarv->(t^(1/3) Vo^(1/2) \[CapitalDelta]\[Gamma]^(1/2))/(Ep^(1/6) Mp^(1/3)),
			Lstarh->(Ep^(1/4) Vo^(1/4))/(\[CapitalDelta]\[Gamma]^(1/4)),
			wstar->(Vo^(1/4) Mp^(1/3))/(Ep^(1/12) t^(1/3) \[CapitalDelta]\[Gamma]^(1/4)),
			pstar->(Ep^(2/3) Mp^(1/3))/(t^(1/3)),
			\[ScriptCapitalK] -> (Kp t^(1/3) \[CapitalDelta]\[Gamma]^(1/8))/(Ep^(19/24) Vo^(1/8) Mp^(1/3)),
			Vt -> Vo - (Ep^(5/6) Vo^(1/2) Mp^(2/3))/(t^(2/3) \[CapitalDelta]\[Gamma]^(3/2))
			} /.data//N,
	"MKh",
			{Lstarv->Kp^(2/3)/\[CapitalDelta]\[Gamma]^(2/3),
			Lstarh->Kp^(2/3)/\[CapitalDelta]\[Gamma]^(2/3),
			wstar->Kp^(4/3)/(Ep \[CapitalDelta]\[Gamma]^(1/3)),
			pstar->Kp^(2/3) \[CapitalDelta]\[Gamma]^(1/3),
			\[ScriptCapitalM] -> (Ep^(6/5) Vo^(2/3) \[CapitalDelta]\[Gamma]^(4/15) Mp^(2/5))/(Kp^(28/15) Qo^(4/15) t^(2/3)),
			Vh -> Kp^(8/3)/(Ep \[CapitalDelta]\[Gamma]^(5/3))
			} /.data//N,
	"MKt",
			{Lstarv-> (Kp^(4/5) t^(1/3) Vo^(2/3) \[CapitalDelta]\[Gamma]^(3/5))/(Ep^(4/5) Qo^(4/15) Mp^(3/5)),
			Lstarh->(Ep^(6/5) Qo^(2/5) Mp^(2/5))/(Kp^(6/5) \[CapitalDelta]\[Gamma]^(2/5)),
			wstar->(Kp^(2/5) Vo^(1/3) Mp^(1/5))/(t^(1/3) Ep^(2/5) Qo^(2/15) \[CapitalDelta]\[Gamma]^(1/5)),
			pstar->(Kp^(8/5) Vo^(1/3) \[CapitalDelta]\[Gamma]^(1/5))/(Ep^(3/5) Qo^(8/15) t^(1/3) Mp^(1/5)),
			Vt -> Vo - Kp^(8/3)/(Ep \[CapitalDelta]\[Gamma]^(5/3))
			} /.data//N,
	"Khh",
			{Lstarv->Kp^(2/3)/\[CapitalDelta]\[Gamma]^(2/3),
			Lstarh->Kp^(2/3)/\[CapitalDelta]\[Gamma]^(2/3),
			wstar->Kp^(4/3)/(Ep \[CapitalDelta]\[Gamma]^(1/3)),
			pstar->Kp^(2/3) \[CapitalDelta]\[Gamma]^(1/3),
			\[ScriptCapitalM] -> (Ep^3 Vo \[CapitalDelta]\[Gamma]^(2/3) Mp)/(Kp^(14/3) t),
			Vh -> Kp^(8/3)/(Ep \[CapitalDelta]\[Gamma]^(5/3))
			} /.data//N,
	"Kht",
			{Lstarv->(t^(1/3) Vo^(2/3) \[CapitalDelta]\[Gamma]^(7/9))/(Kp^(4/9) Mp^(1/3)),
			Lstarh->Kp^(2/3)/\[CapitalDelta]\[Gamma]^(2/3),
			wstar->(Vo^(1/3) Mp^(1/3))/(Kp^(2/9) t^(1/3) \[CapitalDelta]\[Gamma]^(1/9)),
			pstar->(Ep Vo^(1/3) Mp^(1/3) \[CapitalDelta]\[Gamma]^(5/9))/(Kp^(8/9) t^(1/3)),
			Vt -> Vo-Kp^(8/3)/(Ep \[CapitalDelta]\[Gamma]^(5/3))
			} /.data//N
]
];


pulseTimeScales[inpData_,prime_:True,dimless_:False] := 
Module[
   		{data,phi,tmk,tmmt,tmtkt,tkkt,tmb,tmhkhkbh, tmbkbhead, tmbkb},
   		
   		data = inputDataTransformation[inpData,prime];
   		
   		If[dimless,
   		{(Cp^2 Ep^(21/5) \[Mu]p Vo^(1/5))/Kp^(26/5),
   		(Ep^(13/5) Vo^(3/5) \[Mu]p)/Kp^(18/5),
   		(\[Mu]p^(4/13) Vo^(6/13))/(Cp^(18/13) Ep^(4/13)),
   		(Cp^(2/3) Ep^4 \[Mu]p^(4/3) Vo^(2/3))/Kp^(16/3),
   		(Kp^(8/5) Vo^(2/5))/(Cp^2 Ep^(8/5)),
   		(Ep^(5/4) \[Mu]p)/(Vo^(3/4) \[CapitalDelta]\[Gamma]^(9/4)),
   		((Ep^(11/4) Vo^(3/4) \[CapitalDelta]\[Gamma]^(1/4) \[Mu]p)/Kp^4),
   		(Ep^(19/8) Vo^(3/8) \[Mu]p)/(Kp^3 \[CapitalDelta]\[Gamma]^(3/8))
   		}//PowerExpand
   		
   ,
   		phi =
   		 If[
   			(Cp/.data)==0,
   			0,
   			((Cp^2 Ep^(21/5) Mp Vo^(1/5))/Kp^(26/5))/.data
   		];
   
      	(* Times along the edges of the rectangular space MKK~M~ *)
   		
   		tmk = (Ep^(13/5) Vo^(3/5) Mp)/Kp^(18/5)/.data;
   		tmmt = If[(Cp/.data)== 0, Infinity, (Mp^(4/13) Vo^(6/13))/(Cp^(18/13) Ep^(4/13))/.data];
   		tmtkt = If[(Cp/.data)==0,0,(Cp^(2/3) Ep^4 Mp^(4/3) Vo^(2/3))/Kp^(16/3)/.data];	
		   tkkt = If[(Cp/.data) == 0, Infinity, (Kp^(8/5) Vo^(2/5))/(Cp^2 Ep^(8/5))/.data];
		   tmb = If[(\[CapitalDelta]\[Gamma]/.data) == 0, Infinity, (Ep^(5/4) Mp)/(Vo^(3/4) \[CapitalDelta]\[Gamma]^(9/4))/.data];
		   tmbkbhead = If[(\[CapitalDelta]\[Gamma]/.data) == 0, Infinity, ((Ep^(11/4) Vo^(3/4) \[CapitalDelta]\[Gamma]^(1/4) Mp)/Kp^4)/.data];
		   tmbkb= If[(\[CapitalDelta]\[Gamma]/.data) == 0, Infinity, (Ep^(19/8) Vo^(3/8) Mp)/(Kp^3 \[CapitalDelta]\[Gamma]^(3/8))/.data];
			
   		{phi, tmk, tmmt, tmtkt, tkkt, tmb, tmbkbhead, tmbkb}
   		
   		]
   
];


pulseVertexSolutions[V_?StringQ,inpData_,\[Rho]_,t_,dimensional_:True,prime_:True] := Module[{data},

data = inputDataTransformation[inpData,prime];

Switch[V,
	"M",
			If[dimensional,
			{FractureLength-> 0.8359862276159556` Lstar,
			Opening-> InterpolatingFunction[{{-0.9999802608561371, 0.9999802608561371}}, {5, 7, 0, {499}, {4}, 0, 0, 0, 0, Automatic, {}, {}, False}, CompressedData["
1:eJwtl3lcjWkbxyMzVJgiSimlkQxZkjX5MSNCKkIJmddMifFqxhZCkiWyFjVR
lH0LNTGKXGVaMK2ipEV1TtvpnPO8GlGJ3j5z3X/U86lOz3I/9/X7fb/mq/0W
eXfX0NDo1vX1T9e3mOTYM8WdEgVqqNd5dh3vmo/YGPpFoo6otYlrP0t0QDdm
dMsniWxcDiqs2yWSJj9R9G+VKDs11fZqi0RFOZRX3yzRlx/k9EqSKK5isO5v
SolajMc9TW6QaKXN43d35F3nUU6MdKuSaEW47K9LZRK9za8OvlgsUUbot3KX
Qola8zc1X/+76xi7PjUhU6LupXsWryGJTDQ+Zj9/IFG3oQqr8gSJCuPmHjp7
Q6IHOqWfel2UaLl8zKlvz0o0Ymuyjzqs6zrRSSE+h7vu68h43aNBErnHGOuv
3tZ1f1rDbtVu6Hqu+2/+Z+Atkeu7qE+tnhKZWofVhbh2XTc9Iyfboet4qcfz
e1Mlemr2oX3RGIn2hfbddd6i67zt4RvOGkik+z5aPVtHoipTXZMLX9R0Ubpo
FP9OTUsXV/bylavp4Ep97aJiNe2vXTtT/VRNGycUvkhOUdOca9qyCbfU9GD+
pcjVMWr6JrBt4IxjalIftwnJ262mzsWz+vf1U9PJLKnqi5earE+nG8U5d33u
sV7Te3s1XTDR3NU+Sk0lnXs7EozVdPy2Y4SxtpoqP7/0n9SqIvfi3/K/qlOR
Q+q4nANFKpJ7SjEpaSoK++K0JS5eRa1mScF2Z1Q0Wb74U8gBFWmHDfvm8EYV
qQxqmqd7qcgxSKfs6lwVhfvVd2TYqmjficBjEUNUNL1fXNZgbRUdHlBe7v5e
SYamwVqOlUpK23NjvypbSSX5w0MdE5QU5nBy8bIzSlo/yt7ELFhJlx1mmsf8
oqRn1htjCtyUtIhW56fYKclj008tyy2U5J4W5JykrSRd/4NGGe+aaLdJe1To
6yaqmmPUV4eaqEWqefD95Saye96rwia0ibyi2lLLfm2ivwLLztsvbSJblyPZ
HnZNdL84doe1WRPtG/1z58MeTTSme0l0j0YF1b2+faZnroLali+ZmXlXQdvP
1pfan1KQ9sIJt3/zV1CLWVqv/3gqaFVpim1vewX9eafX7k1DFOS9c9q4yG4K
uvo+8qy/rJHu1WjrDMhspEErbPL8rjSSwwd325CDjbS77dejXr6N1HzS5Lt/
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

"], {Developer`PackedArrayForm, CompressedData["
1:eJwV1AOXEAwUBNBs27XZtm1udm22bdu2bdu2bdt95t13zp35BfMCgroEdg4Z
IkSIn8HhgisUoQlDWMIRnghEJBKRiUJUohGdGMQkFrGJQ1ziEZ8EJCQRiUlC
UpKRnBQEkJJUpCYNaUlHejKQkUxkJgtZyUZ2cpCTXOQmD3nJR34KUJBCFKYI
RSlGcUpQklKUpgxlKUd5KlCRSlSmClWpRnVqEEhNalGbOtSlHvVpQEMa0Zgm
NKUZzQmiBS1pRWva0JZ2tKcDHelEZ7rQlW50pwc96UVv+tCXfvRnAAMZxGCG
MJRhDGcEIxnFaMYwlnGMZwITmcRkpjCVaUxnBjOZxWzmMJd5zGcBC1nEYpaw
lGUsZwUrWcVq1rCWdaxnAxvZxGa2sJVtbGcHO9nFbvawl33s5wAHOcRhjnCU
YxznBCc5xWnOcJZznOcCF7nEZa5wlWtc5wY3ucVt7nCXe9znAQ95xGOe8JRn
POcFL3nFa97wlne85wMf+cRnvvCVb3znBz/5hV/5jd/5gz/5i7/5h3/5j+Dx
hyQUoQlDWMIRnghEJBKRiUJUohGdGMQkFrGJQ1ziEZ8EJCQRiUlCUpKRnBQE
kJJUpCYNaUlHejKQkUxkJgtZyUZ2cpCTXOQmD3nJR34KUJBCFKYIRSlGcUpQ
klKUpgxlKUd5KlCRSlSmClWpRnVqEEhNalGbOtSlHvVpQEMa0ZgmNKUZzQmi
BS1pRWva0JZ2tKcDHelE8APvQle60Z0e9KQXvelDX/rRnwEMZBCDGcJQhjGc
EYxkFKMZw1jGMZ4JTGQSk5nCVKYxnRnMZBazmcNc5jGfBSxkEYtZwlKWsZwV
rGQVq1nDWtaxng1sZBOb2cJWtrGdHexkF7vZw172sZ8DHOQQhznCUY5xnBOc
5BSnOcNZznGeC1zkEpe5wlWucZ0b3OQWt7nDXe5xnwc85BGPecJTnvGcF7zk
Fa95w1ve8Z4PfOQTn/nCV77xnR/8D2hl9us=
"], CompressedData["
1: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
+fEyA2552Mj/Ihjw3H2dTaKnr+lmx3PiGZTcu14wdJfBF+06f9Y+YHD+2Ce3
YxkMVr8W5eXkMnhvN7G2tJDByKn7zM+lDLQ3rzo4UD59fzK26ktrGGz8XB/9
pp6BanDg61AOgxXi4TU9TQzESu1O8VoZbHviH/ugY/r5Qny4Jl0MlrRnWUv1
MFiQecD2Yh+DTi3H7kI+g6ihhUnKIgbVNyuXOzAM2h/bXTCfvl4VXOGUCxh0
0O819Ho0vZ9L/89+niJ9vgJ9H/v9s+h+2Ptj79eA7v8cPY8mPR/7vG/p+TVo
Pc7R+sjTehXR+vFpPdn1Zde7gNaf3Q92f9j9ukr7J6Y7s5+qtL/sfvfT/vvS
PLDzwc7LK5ofdp7Y+cqgeWPnT47mkZ3PaprX8zS/B2me2fk+SvM+QvPPngf2
fLDn5QGdn0N0ntjzxZ63HXT+fmDPI51P9rzOpvPLnucIOt88Ou9hdP6LKQ/Y
fPClvNhL+RFIeVJJ+cKhvGHz5zzlEZtP4ZRXKpRfbpRnYZRvhpR3KpR/NZSH
bD6epLw0p/xMpjxl8/UO5S2bv52Ux0WUz82U12x+p1Oe11O+11HeG1H+s31g
Sv1wnPoilPqD7RO2X5ZQ38hQ/9RRH82jfvpAfaVE/RVLfbad+q2f+s6C+u8P
6kO2H9dSX16k/txOfcr26yrqW0/q37XUx47Uz/uorx9Rf8+hPmf7fRn1fSf1
vwbxwDfEB4uIFxYQP7A8oUp8sZR4g+WPEOKRLuKT3cQresQvs4lnQHxTT7zD
8k8J8ZA68RGHeGkN8dNe4qnvia92EW8lEn+xPCYgPltOvNZC/LaJeO4u8R3L
e7nEf/nEg0eID9OJF0OIH1mevEB82UO8qUP8OUg8Wk18up54VZv41ZR4VpL4
9jHx7ibiX2ni4dvEx7rEyxHEz5uIp3WIr4OIt/WJvzcQj7N8bkO83kD8vpJ4
fhHxfTXxvljtDP//RT5gSX7A+gLrD7XkEx/IL/71mvGNCPIPd/KRJPITe/KV
TvIX1meyyG8ekO9sJP+xIR/qJz8yIF9KIH9KJZ9qI7/KJt9yJ/+aTT62hfzs
LPlaC/nbS/I5DvmdHvneb+R/cuSDtuSHVuSL7eSPy8knpcgvC8k3K8k/D5OP
ishPJcdmfHUL+avG1IzPypLffku+e4n8V0g+fJT8+CfyZdnsGX/OJZ/eT35t
uG/Gt/cmz/i3q/WMj/8fpp3tFA==
"]}, {Automatic}][\[Rho]] wstar,
			NetPressure-> InterpolatingFunction[{{-0.9999802608561371, 0.9999802608561371}}, {5, 7, 0, {499}, {4}, 0, 0, 0, 0, Automatic, {}, {}, False}, CompressedData["
1:eJwtl3lcjWkbxyMzVJgiSimlkQxZkjX5MSNCKkIJmddMifFqxhZCkiWyFjVR
lH0LNTGKXGVaMK2ipEV1TtvpnPO8GlGJ3j5z3X/U86lOz3I/9/X7fb/mq/0W
eXfX0NDo1vX1T9e3mOTYM8WdEgVqqNd5dh3vmo/YGPpFoo6otYlrP0t0QDdm
dMsniWxcDiqs2yWSJj9R9G+VKDs11fZqi0RFOZRX3yzRlx/k9EqSKK5isO5v
SolajMc9TW6QaKXN43d35F3nUU6MdKuSaEW47K9LZRK9za8OvlgsUUbot3KX
Qola8zc1X/+76xi7PjUhU6LupXsWryGJTDQ+Zj9/IFG3oQqr8gSJCuPmHjp7
Q6IHOqWfel2UaLl8zKlvz0o0Ymuyjzqs6zrRSSE+h7vu68h43aNBErnHGOuv
3tZ1f1rDbtVu6Hqu+2/+Z+Atkeu7qE+tnhKZWofVhbh2XTc9Iyfboet4qcfz
e1Mlemr2oX3RGIn2hfbddd6i67zt4RvOGkik+z5aPVtHoipTXZMLX9R0Ubpo
FP9OTUsXV/bylavp4Ep97aJiNe2vXTtT/VRNGycUvkhOUdOca9qyCbfU9GD+
pcjVMWr6JrBt4IxjalIftwnJ262mzsWz+vf1U9PJLKnqi5earE+nG8U5d33u
sV7Te3s1XTDR3NU+Sk0lnXs7EozVdPy2Y4SxtpoqP7/0n9SqIvfi3/K/qlOR
Q+q4nANFKpJ7SjEpaSoK++K0JS5eRa1mScF2Z1Q0Wb74U8gBFWmHDfvm8EYV
qQxqmqd7qcgxSKfs6lwVhfvVd2TYqmjficBjEUNUNL1fXNZgbRUdHlBe7v5e
SYamwVqOlUpK23NjvypbSSX5w0MdE5QU5nBy8bIzSlo/yt7ELFhJlx1mmsf8
oqRn1htjCtyUtIhW56fYKclj008tyy2U5J4W5JykrSRd/4NGGe+aaLdJe1To
6yaqmmPUV4eaqEWqefD95Saye96rwia0ibyi2lLLfm2ivwLLztsvbSJblyPZ
HnZNdL84doe1WRPtG/1z58MeTTSme0l0j0YF1b2+faZnroLali+ZmXlXQdvP
1pfan1KQ9sIJt3/zV1CLWVqv/3gqaFVpim1vewX9eafX7k1DFOS9c9q4yG4K
uvo+8qy/rJHu1WjrDMhspEErbPL8rjSSwwd325CDjbS77dejXr6N1HzS5Lt/
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

"], {Developer`PackedArrayForm, CompressedData["
1:eJwV1AOXEAwUBNBs27XZtm1udm22bdu2bdu2bdt95t13zp35BfMCgroEdg4Z
IkSIn8HhgisUoQlDWMIRnghEJBKRiUJUohGdGMQkFrGJQ1ziEZ8EJCQRiUlC
UpKRnBQEkJJUpCYNaUlHejKQkUxkJgtZyUZ2cpCTXOQmD3nJR34KUJBCFKYI
RSlGcUpQklKUpgxlKUd5KlCRSlSmClWpRnVqEEhNalGbOtSlHvVpQEMa0Zgm
NKUZzQmiBS1pRWva0JZ2tKcDHelEZ7rQlW50pwc96UVv+tCXfvRnAAMZxGCG
MJRhDGcEIxnFaMYwlnGMZwITmcRkpjCVaUxnBjOZxWzmMJd5zGcBC1nEYpaw
lGUsZwUrWcVq1rCWdaxnAxvZxGa2sJVtbGcHO9nFbvawl33s5wAHOcRhjnCU
YxznBCc5xWnOcJZznOcCF7nEZa5wlWtc5wY3ucVt7nCXe9znAQ95xGOe8JRn
POcFL3nFa97wlne85wMf+cRnvvCVb3znBz/5hV/5jd/5gz/5i7/5h3/5j+Dx
hyQUoQlDWMIRnghEJBKRiUJUohGdGMQkFrGJQ1ziEZ8EJCQRiUlCUpKRnBQE
kJJUpCYNaUlHejKQkUxkJgtZyUZ2cpCTXOQmD3nJR34KUJBCFKYIRSlGcUpQ
klKUpgxlKUd5KlCRSlSmClWpRnVqEEhNalGbOtSlHvVpQEMa0ZgmNKUZzQmi
BS1pRWva0JZ2tKcDHelE8APvQle60Z0e9KQXvelDX/rRnwEMZBCDGcJQhjGc
EYxkFKMZw1jGMZ4JTGQSk5nCVKYxnRnMZBazmcNc5jGfBSxkEYtZwlKWsZwV
rGQVq1nDWtaxng1sZBOb2cJWtrGdHexkF7vZw172sZ8DHOQQhznCUY5xnBOc
5BSnOcNZznGeC1zkEpe5wlWucZ0b3OQWt7nDXe5xnwc85BGPecJTnvGcF7zk
Fa95w1ve8Z4PfOQTn/nCV77xnR/8D2hl9us=
"], CompressedData["
1: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
886gPxsGoU9/GkS/1opA3xYT/1Z8gD7+f9AFESs=
"]}, {Automatic}][\[Rho]] pstar}/.pulseVertexScalings[V,inpData,t,prime]/.data//N,
			{FractureLength-> 0.8359862276159556`,
			Opening-> InterpolatingFunction[{{-0.9999802608561371, 0.9999802608561371}}, {5, 7, 0, {499}, {4}, 0, 0, 0, 0, Automatic, {}, {}, False}, CompressedData["
1:eJwtl3lcjWkbxyMzVJgiSimlkQxZkjX5MSNCKkIJmddMifFqxhZCkiWyFjVR
lH0LNTGKXGVaMK2ipEV1TtvpnPO8GlGJ3j5z3X/U86lOz3I/9/X7fb/mq/0W
eXfX0NDo1vX1T9e3mOTYM8WdEgVqqNd5dh3vmo/YGPpFoo6otYlrP0t0QDdm
dMsniWxcDiqs2yWSJj9R9G+VKDs11fZqi0RFOZRX3yzRlx/k9EqSKK5isO5v
SolajMc9TW6QaKXN43d35F3nUU6MdKuSaEW47K9LZRK9za8OvlgsUUbot3KX
Qola8zc1X/+76xi7PjUhU6LupXsWryGJTDQ+Zj9/IFG3oQqr8gSJCuPmHjp7
Q6IHOqWfel2UaLl8zKlvz0o0Ymuyjzqs6zrRSSE+h7vu68h43aNBErnHGOuv
3tZ1f1rDbtVu6Hqu+2/+Z+Atkeu7qE+tnhKZWofVhbh2XTc9Iyfboet4qcfz
e1Mlemr2oX3RGIn2hfbddd6i67zt4RvOGkik+z5aPVtHoipTXZMLX9R0Ubpo
FP9OTUsXV/bylavp4Ep97aJiNe2vXTtT/VRNGycUvkhOUdOca9qyCbfU9GD+
pcjVMWr6JrBt4IxjalIftwnJ262mzsWz+vf1U9PJLKnqi5earE+nG8U5d33u
sV7Te3s1XTDR3NU+Sk0lnXs7EozVdPy2Y4SxtpoqP7/0n9SqIvfi3/K/qlOR
Q+q4nANFKpJ7SjEpaSoK++K0JS5eRa1mScF2Z1Q0Wb74U8gBFWmHDfvm8EYV
qQxqmqd7qcgxSKfs6lwVhfvVd2TYqmjficBjEUNUNL1fXNZgbRUdHlBe7v5e
SYamwVqOlUpK23NjvypbSSX5w0MdE5QU5nBy8bIzSlo/yt7ELFhJlx1mmsf8
oqRn1htjCtyUtIhW56fYKclj008tyy2U5J4W5JykrSRd/4NGGe+aaLdJe1To
6yaqmmPUV4eaqEWqefD95Saye96rwia0ibyi2lLLfm2ivwLLztsvbSJblyPZ
HnZNdL84doe1WRPtG/1z58MeTTSme0l0j0YF1b2+faZnroLali+ZmXlXQdvP
1pfan1KQ9sIJt3/zV1CLWVqv/3gqaFVpim1vewX9eafX7k1DFOS9c9q4yG4K
uvo+8qy/rJHu1WjrDMhspEErbPL8rjSSwwd325CDjbS77dejXr6N1HzS5Lt/
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

"], {Developer`PackedArrayForm, CompressedData["
1:eJwV1AOXEAwUBNBs27XZtm1udm22bdu2bdu2bdt95t13zp35BfMCgroEdg4Z
IkSIn8HhgisUoQlDWMIRnghEJBKRiUJUohGdGMQkFrGJQ1ziEZ8EJCQRiUlC
UpKRnBQEkJJUpCYNaUlHejKQkUxkJgtZyUZ2cpCTXOQmD3nJR34KUJBCFKYI
RSlGcUpQklKUpgxlKUd5KlCRSlSmClWpRnVqEEhNalGbOtSlHvVpQEMa0Zgm
NKUZzQmiBS1pRWva0JZ2tKcDHelEZ7rQlW50pwc96UVv+tCXfvRnAAMZxGCG
MJRhDGcEIxnFaMYwlnGMZwITmcRkpjCVaUxnBjOZxWzmMJd5zGcBC1nEYpaw
lGUsZwUrWcVq1rCWdaxnAxvZxGa2sJVtbGcHO9nFbvawl33s5wAHOcRhjnCU
YxznBCc5xWnOcJZznOcCF7nEZa5wlWtc5wY3ucVt7nCXe9znAQ95xGOe8JRn
POcFL3nFa97wlne85wMf+cRnvvCVb3znBz/5hV/5jd/5gz/5i7/5h3/5j+Dx
hyQUoQlDWMIRnghEJBKRiUJUohGdGMQkFrGJQ1ziEZ8EJCQRiUlCUpKRnBQE
kJJUpCYNaUlHejKQkUxkJgtZyUZ2cpCTXOQmD3nJR34KUJBCFKYIRSlGcUpQ
klKUpgxlKUd5KlCRSlSmClWpRnVqEEhNalGbOtSlHvVpQEMa0ZgmNKUZzQmi
BS1pRWva0JZ2tKcDHelE8APvQle60Z0e9KQXvelDX/rRnwEMZBCDGcJQhjGc
EYxkFKMZw1jGMZ4JTGQSk5nCVKYxnRnMZBazmcNc5jGfBSxkEYtZwlKWsZwV
rGQVq1nDWtaxng1sZBOb2cJWtrGdHexkF7vZw172sZ8DHOQQhznCUY5xnBOc
5BSnOcNZznGeC1zkEpe5wlWucZ0b3OQWt7nDXe5xnwc85BGPecJTnvGcF7zk
Fa95w1ve8Z4PfOQTn/nCV77xnR/8D2hl9us=
"], CompressedData["
1:eJw1l3k4lekfxqVCoTSi0iRLIxSmEjXRXQYVmphImyFGMQnTz0xpGUuRmVQo
GTFEEqWyG2uUnaZj55DlbC9FllTSmH6uy/f9x3txjnPe93m+z31/PqrOXt+7
zhITE3s7/eOhpQL3Suc2zLlnd+pC4hEcsAs+6LDCC1Iucb7OfqfAW2E32Jzr
h8ismuJQ54vwKVOXFpf+HUsiFYLnSVxFV4W45pZr4XDL7Zk8bXwDXyn0Ro9N
3MRXX8qvcmmIBu9fpUMZObGw8pzKX5cVj/FRqTXlDgloquyZWDWUCNsxaWFM
QBIiDbNtDVWSMdS9/PbJZ/dgwS96MOqVCptbFXNTVR4gJ9dHIqsjDZvHjyXk
3HyEvfd3DpXqpENxi5yM+BcZ8Ll15ETdxwwsLTXeYszPROnsSsziZGFM7NWo
ZUk2BsVUHsc/zkG42bDsD4m5sM4J7xmKyoOY06649qt/4yde1yqvkHyUb4zj
G18oQCZnIuW8fyGuf73bdPi3IkQMzJ9d6FcMXcfAHyMCSjBm9q19ZNAT7GvS
XbDZrBSGDduUY6pL8ejt82VuVmXIPpP4Wa2hDMEBG/Za7HuKhl0ftG26nsI0
aXxWivMzROVo6VQOPIMHR7t//c/liPI0k3f7WA4x7+/13wRUYEq4UKAlXYkm
vUbVhBuVCFYTHilQrsIpj2eLPqdUwVxgE/9yQzWujQttBkqqMVnDPRpoUYN6
i1t3vmmtgcK8YLOTzrXQOPwotPxNLTqU3p25drYOVsF+Y+5S9Tj6x4tLSZH1
0HVQ6RepPYd29FzFs+nPETh3MFFu6z/4Wrn6/du66WuAe17lgRf485f1tkLm
BdaZd0v99w0HcWIyU5qBHMwrNJ/cWcOBgfCOQ4BcA9Ky9h5PtW/A3Cvng/hx
Ddg5otAtLWqAVkGrzQ6dRmxITVP09GlE7VU9taTCRgyHbqorFG+CtqO6Zseu
JryMdQ1vCWuC/8896Ghrguq+fZnlys24FBGknuvajMc3h/SD0prhMvSx/7u3
zSgueL1idHML8mqG86/6t2Dp/GMDnypbkChMEP4g2wqDrRJzbu1tRf3G38uq
o1uxTHJHUVVPK/r8+Aeyv2qDfsSy/SeOtyHJ0illcUYbLGW1q4Let2G/fuu7
ri3tsD89x2nCvx3a7kq+4xXtqOu4HJc9vwMqn+w1bPd0oOCM8Yn86x1Yo1I3
NtXWgRcOgyMSX3KxTdibOurIhYSm99Nbd7hYt6g2ewXDhe1n1UvO2p3gLmhy
8z/RiSCH2Y626Z34FLXEaGKsE4qVar/YbOyC2tlrcv871QWlT9L2ZvldkLdy
d+RPdkHxgefgeqOXMJpaHbD7/EuoOSrvWVryEu0m1sGJ/73EEkmP8cat3TD0
+M6v4rduSBZaOLmUdCPzKCORP9U9/Txd67ONenBk8eI17md7sG6ZRH15fg/c
o7w92j/0QE2vL7xCqRfOGS1LzmzuxTalX+L/te9FcHHkh/2/9kKvvDEt4EYv
gj5GPPLP7MVlH56MFacXlikbjN4M9YIxSHM4LN2Hrac9ViZo9iHGN2/qoVkf
ltZql1xx7sP7wIhcfb8+zFUI3Zce04fny80MPuT1QdbpXalc8/T7bptpjQ73
IT3MPj1RmoeKk/nqSqt5uGynNeVswoP7Abt2bwceTPzdn313mocQLfPy4Qge
jE08FV0e8vBc+Y+EmCoeOA+vlsX38ZB5golx/cTDSTHvyMnFfLSOGJju0eUj
/8mJTSd38MF1HNWzdeJjnfj9bjFfPuJj3y3yCOcjS+BbcSeVj3urde/cKuPD
LETe2qGDDx872Tm9I3ysnfgrWVdKgKr3ceHGKwVQ1b2/WsFAAIPjVffSrAQI
XztLXtZFgGOWst6avtOvX8yum39NgLA3773vJgmg6GNqM7dAgHMVF3jKL6av
1glpk3wBXni8Vrw2IUC5z2Rlv4wQAdstVSdUhLDN/zqzRl+IM17DfOud039/
xKjFHBKiqNmzONpTiPB22YH9AULw3YsbOdeFuD449qt0shCrXqWZTOYJMac2
ZfhujRApRQrJEp1CmGdyo3UGhTjspG0hMyWEmfhB0zRZERap6BhPrhChTjTL
YL6uCA6Obkb/GIuQ9oT7hfVuEd7GxKpcPCyCuPqAne9xEc44cW+rnRFhtMH5
8e8hIrgKri1PvClC2Ktxm7NJIhx6s11SMlMEv9itq62eiFD795rgHfUiGJbd
lZhqF2HhZb1ZzkIRkv9budFvVIRMMfUJiykRnu7pv8KRYhD+SU4gvZhBv2sO
d95KBopyhgUFWgxiZJ9kK+szOHv4xu1NWxkMrF+psXAnA43FJ8YjbRioXIiS
rz/I4LOO18F7LgxCv3RKNfJgcNJRve9Xn+nPTSlwdjvHIMzk8ZjURQb56zav
+fEyA2552Mj/Ihjw3H2dTaKnr+lmx3PiGZTcu14wdJfBF+06f9Y+YHD+2Ce3
YxkMVr8W5eXkMnhvN7G2tJDByKn7zM+lDLQ3rzo4UD59fzK26ktrGGz8XB/9
pp6BanDg61AOgxXi4TU9TQzESu1O8VoZbHviH/ugY/r5Qny4Jl0MlrRnWUv1
MFiQecD2Yh+DTi3H7kI+g6ihhUnKIgbVNyuXOzAM2h/bXTCfvl4VXOGUCxh0
0O819Ho0vZ9L/89+niJ9vgJ9H/v9s+h+2Ptj79eA7v8cPY8mPR/7vG/p+TVo
Pc7R+sjTehXR+vFpPdn1Zde7gNaf3Q92f9j9ukr7J6Y7s5+qtL/sfvfT/vvS
PLDzwc7LK5ofdp7Y+cqgeWPnT47mkZ3PaprX8zS/B2me2fk+SvM+QvPPngf2
fLDn5QGdn0N0ntjzxZ63HXT+fmDPI51P9rzOpvPLnucIOt88Ou9hdP6LKQ/Y
fPClvNhL+RFIeVJJ+cKhvGHz5zzlEZtP4ZRXKpRfbpRnYZRvhpR3KpR/NZSH
bD6epLw0p/xMpjxl8/UO5S2bv52Ux0WUz82U12x+p1Oe11O+11HeG1H+s31g
Sv1wnPoilPqD7RO2X5ZQ38hQ/9RRH82jfvpAfaVE/RVLfbad+q2f+s6C+u8P
6kO2H9dSX16k/txOfcr26yrqW0/q37XUx47Uz/uorx9Rf8+hPmf7fRn1fSf1
vwbxwDfEB4uIFxYQP7A8oUp8sZR4g+WPEOKRLuKT3cQresQvs4lnQHxTT7zD
8k8J8ZA68RGHeGkN8dNe4qnvia92EW8lEn+xPCYgPltOvNZC/LaJeO4u8R3L
e7nEf/nEg0eID9OJF0OIH1mevEB82UO8qUP8OUg8Wk18up54VZv41ZR4VpL4
9jHx7ibiX2ni4dvEx7rEyxHEz5uIp3WIr4OIt/WJvzcQj7N8bkO83kD8vpJ4
fhHxfTXxvljtDP//RT5gSX7A+gLrD7XkEx/IL/71mvGNCPIPd/KRJPITe/KV
TvIX1meyyG8ekO9sJP+xIR/qJz8yIF9KIH9KJZ9qI7/KJt9yJ/+aTT62hfzs
LPlaC/nbS/I5DvmdHvneb+R/cuSDtuSHVuSL7eSPy8knpcgvC8k3K8k/D5OP
ishPJcdmfHUL+avG1IzPypLffku+e4n8V0g+fJT8+CfyZdnsGX/OJZ/eT35t
uG/Gt/cmz/i3q/WMj/8fpp3tFA==
"]}, {Automatic}][\[Rho]],
			NetPressure-> InterpolatingFunction[{{-0.9999802608561371, 0.9999802608561371}}, {5, 7, 0, {499}, {4}, 0, 0, 0, 0, Automatic, {}, {}, False}, CompressedData["
1:eJwtl3lcjWkbxyMzVJgiSimlkQxZkjX5MSNCKkIJmddMifFqxhZCkiWyFjVR
lH0LNTGKXGVaMK2ipEV1TtvpnPO8GlGJ3j5z3X/U86lOz3I/9/X7fb/mq/0W
eXfX0NDo1vX1T9e3mOTYM8WdEgVqqNd5dh3vmo/YGPpFoo6otYlrP0t0QDdm
dMsniWxcDiqs2yWSJj9R9G+VKDs11fZqi0RFOZRX3yzRlx/k9EqSKK5isO5v
SolajMc9TW6QaKXN43d35F3nUU6MdKuSaEW47K9LZRK9za8OvlgsUUbot3KX
Qola8zc1X/+76xi7PjUhU6LupXsWryGJTDQ+Zj9/IFG3oQqr8gSJCuPmHjp7
Q6IHOqWfel2UaLl8zKlvz0o0Ymuyjzqs6zrRSSE+h7vu68h43aNBErnHGOuv
3tZ1f1rDbtVu6Hqu+2/+Z+Atkeu7qE+tnhKZWofVhbh2XTc9Iyfboet4qcfz
e1Mlemr2oX3RGIn2hfbddd6i67zt4RvOGkik+z5aPVtHoipTXZMLX9R0Ubpo
FP9OTUsXV/bylavp4Ep97aJiNe2vXTtT/VRNGycUvkhOUdOca9qyCbfU9GD+
pcjVMWr6JrBt4IxjalIftwnJ262mzsWz+vf1U9PJLKnqi5earE+nG8U5d33u
sV7Te3s1XTDR3NU+Sk0lnXs7EozVdPy2Y4SxtpoqP7/0n9SqIvfi3/K/qlOR
Q+q4nANFKpJ7SjEpaSoK++K0JS5eRa1mScF2Z1Q0Wb74U8gBFWmHDfvm8EYV
qQxqmqd7qcgxSKfs6lwVhfvVd2TYqmjficBjEUNUNL1fXNZgbRUdHlBe7v5e
SYamwVqOlUpK23NjvypbSSX5w0MdE5QU5nBy8bIzSlo/yt7ELFhJlx1mmsf8
oqRn1htjCtyUtIhW56fYKclj008tyy2U5J4W5JykrSRd/4NGGe+aaLdJe1To
6yaqmmPUV4eaqEWqefD95Saye96rwia0ibyi2lLLfm2ivwLLztsvbSJblyPZ
HnZNdL84doe1WRPtG/1z58MeTTSme0l0j0YF1b2+faZnroLali+ZmXlXQdvP
1pfan1KQ9sIJt3/zV1CLWVqv/3gqaFVpim1vewX9eafX7k1DFOS9c9q4yG4K
uvo+8qy/rJHu1WjrDMhspEErbPL8rjSSwwd325CDjbS77dejXr6N1HzS5Lt/
HBvJvPPq8tkjGulouuZmD61GSszKyRzR2EBO4bMvpTxtoMZflv+kc62BcgPO
ug8+2EDW7XbFKu8GurNg4IQ9sxpo1eirz4uGNlB5Zh91g0YDfczV/kiV9XQ9
vtN+2aN6WntKZXwvqp6qwzcrirbWk/6DF5r33eqp0kT2ZPnYrmNPk+DM3vVU
8S4s9p+GOjp24Z1/XUYd5dgF+pyLraPZc42eG++sowu5k7utcq+jQc7163xt
6uhwebbP5D51pG+mtyS/vpaMVhyNsH5SS41boo64RddSZqtF1MyttXTc8JZp
q0stNSXO2OM/opY+OQeap3WvpXPG4T4vyuTkpByfdDtJTvtsL/stOSqnR93u
fcnyllP/9ynZmtPlFDAhclzvgXIK80u4WKWSkWHax/1BmTLSulUxQ4qWkbfe
cvvvNsvotEGT3Ga+jA6YhsZpDZWR/63kitutNbTSc0sPi/waGkV2QT6Xa+jm
sl5p2wNqyGC7m/XKhTUU9K3eRL3hNdSgcXdOREc1GU00LFAWVpPXohv2A65W
U/6o6Z3f7KymW1nm+ypcq8ltpOG03cOqqZ+LR2RzWxXJlaFf2+VVkWJ1Xrv7
hSp6dvRm7dytXb+/YDpFb14Vbdc6uea2SRVN8k7M/qHmLdlUbfz+dchbOnIq
1dJ4zFtafJ1W9HlVSVs1/XPv7KikXMWd1ZpmleQY372nVmYFvUfh5fS1FdRZ
MWvk6L4VdO/y/H0uieUUND7inOXScmqZeto3qa2MVCb+hW3RZRQTNOzZR5SR
Vs33M+7UvKEFFlZW5vvfkL15zUbn4W/I0n5M/8nPSokad3ZUrislnwm1Q2f1
LqWdi0YF+Ma/pmmmj9udnF/T3cs3I5vVJTTo+jVHz+MlFGGi9dX+MSXkEf84
Y1N+MfWduHTPcL9ikjm0jI3rW0zlEe4F1fGvKDEkx0Ph9Ioco2wzHza9pINl
X3TdD7+k8vWWk1OtXpJz45hJUlYRGVZa6qh/LqJ+dmcTUroXUcSTc1ZLY1/Q
pG1PNqTbv6DeffWDP78pJJmrsbfOtkIyCi/RV+oX0uyBqWGxCQWUvCnx7Qjn
AjrRc3C73qt8Gh5iWbqrTz6N3F71Y5RDHi3M3hLhsyuXdLtXBNYn5dDqIS8s
TJV/01v7wkO9Lf6m8J/1UpKWPachMwNTDU88o2rZnzFTs57Snh2a3oM7sul0
L1fTtHHZVGA66tXwNVm04ET66SXRmbTTsmD9/MIMCnC7tkbn6wwqrRl2KnLq
X+Rw/YZG84YntNdm0KPBF9Np/ob8V4NK0ki+PG5lk3Yake/VbMtdj8nAyrp5
uPwRqWZprKyc95Cmb+o70jkhmd4bdAQGGDygpPBBe9btuk8LQz1ch8iS6J/R
Mweec/yD/gis/6iIT6CA8Cpj7f53ae9dv0Wd0i3Ku/H02Pnk63TvdcC0z3uv
0NqxtTfNnC6Sl8dH634DYql4jpNV0sMoemx04cduC8Kp7vH9RVMeHaJfP5yM
8Xy0nXqu0B9veeXktM/8M6zp37+jmj+Pdv5/+PH5sIPPj2y+Hqr4+gjn+8E+
vj884vtFB98/PPl5kMrPh8/8vJjEz48PvB7Q5fXBA14vVPH6wYXXE4G8vpjJ
641iXn/48/tAEL8fzOP3hRx+fzjJ7xO7+P2ilt83LPj94zzvByh5f8CP9wsG
8P6BF+8nTOL9hfG83xDN+w8PeT/CifcnzHi/opH3L/rxfsYU3t/4nfc79Hn/
w5jnAa48H6jkecEhnh/M43lCEs8XKnneUMvzB12eR3jyfOJ3nlcY8/wikecZ
03m+sZvnHb48/0jnPIAV5wPAeQEXzg/ocJ4gmvMFTZw3aOf8wQHOIzzkfEK3
yn/zCh84vzCP8wz5nG/YzHkHF84/HOA8hC3nI6ZwXmIH5ydknKd4yvmKRs5b
yDh/ocd5DFfOZ8RzXqOA8xurOM9hyvkOBec9gjn/YcR9gHjuB4zmvoAX9we2
cZ/gMPcLfue+gS/3D3pzH8GI+wmnuK+wi/sL+txneMj9hmDuO7hw/+E89yE6
uB+h5L7ECe5PiD6F6FeIvkV/7l8c4T6GEfczLnFfYw73N3K5z3Gc+x2i7yH6
H4IHMJD5ADLmBaxjfsAN5gkIvoDgDXgxf0DwCEYxnyCPeQUK5hcsYJ7BH8w3
ELwDwT8QPATBRxC8BEPmJwiewjXmKwjeguAveDGPQfAZtJjXsI35De3Mc6hn
vsNY5j0cYP7DA+ZBTGQ+RAbzIn5kfsQ05kkIvoTgTQj+RD/mUSxjPsUy5lW4
Mb/iOfMsrjDf4r/MuxD8C8HDEHyMQczLCGV+huBpCL6G4G0I/obgcQg+xyTm
dbQxvyOceR61zPcQvA/B/xA+AOEHEL4A4Q/QZZ/AaPYLhLFvQPgHJPYRCD/B
n+wrmM3+AuEzEH4D4TsQ/gPhQxB+BOFLEP6EYPYpZLNfQfgWhH9B+BiEn+Ez
+xo62d/gwT4H4XcQvgfhf1jBPohk9kMIX4TwR5iyT0KT/RJt7JsQ/gnhoxB+
CuGrEP4K4bMQfgvhuxD+C+HDEH4M4csQ/gzh0xB+DeHbEP4N4eP4P7gplQg=

"], {Developer`PackedArrayForm, CompressedData["
1:eJwV1AOXEAwUBNBs27XZtm1udm22bdu2bdu2bdt95t13zp35BfMCgroEdg4Z
IkSIn8HhgisUoQlDWMIRnghEJBKRiUJUohGdGMQkFrGJQ1ziEZ8EJCQRiUlC
UpKRnBQEkJJUpCYNaUlHejKQkUxkJgtZyUZ2cpCTXOQmD3nJR34KUJBCFKYI
RSlGcUpQklKUpgxlKUd5KlCRSlSmClWpRnVqEEhNalGbOtSlHvVpQEMa0Zgm
NKUZzQmiBS1pRWva0JZ2tKcDHelEZ7rQlW50pwc96UVv+tCXfvRnAAMZxGCG
MJRhDGcEIxnFaMYwlnGMZwITmcRkpjCVaUxnBjOZxWzmMJd5zGcBC1nEYpaw
lGUsZwUrWcVq1rCWdaxnAxvZxGa2sJVtbGcHO9nFbvawl33s5wAHOcRhjnCU
YxznBCc5xWnOcJZznOcCF7nEZa5wlWtc5wY3ucVt7nCXe9znAQ95xGOe8JRn
POcFL3nFa97wlne85wMf+cRnvvCVb3znBz/5hV/5jd/5gz/5i7/5h3/5j+Dx
hyQUoQlDWMIRnghEJBKRiUJUohGdGMQkFrGJQ1ziEZ8EJCQRiUlCUpKRnBQE
kJJUpCYNaUlHejKQkUxkJgtZyUZ2cpCTXOQmD3nJR34KUJBCFKYIRSlGcUpQ
klKUpgxlKUd5KlCRSlSmClWpRnVqEEhNalGbOtSlHvVpQEMa0ZgmNKUZzQmi
BS1pRWva0JZ2tKcDHelE8APvQle60Z0e9KQXvelDX/rRnwEMZBCDGcJQhjGc
EYxkFKMZw1jGMZ4JTGQSk5nCVKYxnRnMZBazmcNc5jGfBSxkEYtZwlKWsZwV
rGQVq1nDWtaxng1sZBOb2cJWtrGdHexkF7vZw172sZ8DHOQQhznCUY5xnBOc
5BSnOcNZznGeC1zkEpe5wlWucZ0b3OQWt7nDXe5xnwc85BGPecJTnvGcF7zk
Fa95w1ve8Z4PfOQTn/nCV77xnR/8D2hl9us=
"], CompressedData["
1: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
886gPxsGoU9/GkS/1opA3xYT/1Z8gD7+f9AFESs=
"]}, {Automatic}][\[Rho]]}//N],
	"K",
			If[dimensional,
			{FractureLength-> (3/(\[Pi] Sqrt[2]))^(2/5) Lstar,
			Opening-> (3/(8 \[Pi]))^(1/5) (1-\[Rho]^2)^(1/2) wstar,
			NetPressure-> \[Pi]/8 (\[Pi]/12)^(1/5) pstar}/.pulseVertexScalings[V,inpData,t,prime]/.data//N,
			{FractureLength-> (3/(\[Pi] Sqrt[2]))^(2/5),
			Opening-> (3/(8 \[Pi]))^(1/5) (1-\[Rho]^2)^(1/2),
			NetPressure-> \[Pi]/8 (\[Pi]/12)^(1/5)}//N]
]
];


arrestRadius[V_?StringQ,inpData_,dimensional_:True,prime_:True] := Module[{data},

data = inputDataTransformation[inpData,prime];

Switch[V,
	"K",
			If[dimensional,
			(3/(\[Pi] Sqrt[2]))^(2/5) Lstar/.pulseVertexScalings[V,inpData,t,prime]/.data//N,
			(3/(\[Pi] Sqrt[2]))^(2/5)//N],
	"Kt",
			If[dimensional,
			(0.52177177375563246783358808245421881975`17.747909562467935 Ep^(1/13) Vo^(5/13))/(Cp^(2/13) Mp^(1/13))/.data//N,
			0.52177177375563246783358808245421881975`17.747909562467935],
	"Mt",
			If[dimensional,
			(0.52177177375563246783358808245421881975`17.747909562467935 Ep^(1/13) Vo^(5/13))/(Cp^(2/13) Mp^(1/13))/.data//N,
			0.52177177375563246783358808245421881975`17.747909562467935]
]
];


pulseTransitionMKScales[inpData_,prime_:True] := 
Module[{data},
   
	data = inputDataTransformation[inpData,prime];
         
   		{
	   	   Lstar-> (Ep^(2/5) Vo^(2/5))/Kp^(2/5),
			  wstar->(Kp^(4/5) Vo^(1/5))/Ep^(4/5),
			  pstar-> Kp^(6/5)/(Ep^(1/5) Vo^(1/5)),
	      	qstar -> (Kp^(18/5) Vo^(2/5))/(Ep^(13/5) Mp),
	      	tstar -> (Ep^(13/5) Vo^(3/5) Mp)/Kp^(18/5)
      	} /.data//N
];


pulseTransitionMMhScales[inpData_,prime_:True] := 
Module[{data},
   
	data = inputDataTransformation[inpData,prime];
         
   		{
	   		Lstar-> (Ep^(1/4) Vo^(1/4))/\[CapitalDelta]\[Gamma]^(1/4),
			  wstar->(Sqrt[Vo] Sqrt[\[CapitalDelta]\[Gamma]])/Sqrt[Ep],
			  pstar-> Ep^(1/4) Vo^(1/4) \[CapitalDelta]\[Gamma]^(3/4),
	      	qstar -> (Vo^(7/4) \[CapitalDelta]\[Gamma]^(9/4))/(Ep^(5/4) Mp),
	      	tstar -> (Ep^(5/4) Mp)/(Vo^(3/4) \[CapitalDelta]\[Gamma]^(9/4))
      	} /.data//N
];


pulseTransitionMhKhTailScales[inpData_,prime_:True] := 
Module[{data},
   
	data = inputDataTransformation[inpData,prime];
         
   		{
	   	   Lstar-> Ep^(1/4)Vo^(1/4)/\[CapitalDelta]\[Gamma]^(1/4),
			  wstar->Kp^(4/3)/(Ep \[CapitalDelta]\[Gamma]^(1/3)),
			  pstar-> Kp^(4/3)/(Ep^(1/4) Vo^(1/4) \[CapitalDelta]\[Gamma]^(1/12)),
	      	tstar -> ((Ep^(11/4) Vo^(3/4) \[CapitalDelta]\[Gamma]^(1/4) Mp)/Kp^4)
      	} /.data//N
];


pulseTransitionMhKhHeadScales[inpData_,prime_:True] := 
Module[{data},
   
	data = inputDataTransformation[inpData,prime];
         
   		{
	   	   Lstar-> Kp^(2/3)/\[CapitalDelta]\[Gamma]^(2/3),
			  wstar->Kp^(4/3)/(Ep \[CapitalDelta]\[Gamma]^(1/3)),
			  pstar-> Kp^(2/3) \[CapitalDelta]\[Gamma]^(1/3),
	      	tstar -> ((Ep^(11/4) Vo^(3/4) \[CapitalDelta]\[Gamma]^(1/4) Mp)/Kp^4)
      	} /.data//N
];


End[]; (* End Private Context *)

EndPackage[];
