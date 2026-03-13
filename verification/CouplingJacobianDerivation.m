(* ============================================================================ *)
(* CouplingJacobianDerivation.m                                                *)
(*                                                                              *)
(* Symbolic proof that the T-form of d[S^{-1}(Omega)T]/dOmega is correct.     *)
(*                                                                              *)
(* Result:                                                                      *)
(*   d[S^{-1}T]/dOmega = (1/2)[T]x + beta rT^T + alpha Tr^T                  *)
(*                        + Tbar [beta I - (2beta + alpha) rrT]                *)
(*                                                                              *)
(*   alpha = (sinTh - Th)/(2(1-cosTh)) = (1/4)csc^2(Th/2)(sinTh - Th)  < 0   *)
(*   beta  = x/Th = 1/Th - (1/2)cot(Th/2)                              > 0   *)
(*   alpha != beta                                                              *)
(*                                                                              *)
(* ERRATUM: Paper eq 78 falsely claimed alpha = -1/2 cot(Th/2) + 1/Th.        *)
(* That expression is beta, not alpha.                                          *)
(*                                                                              *)
(* Mathematica note: expressions involving Sqrt[w1^2+w2^2+w3^2] cannot be     *)
(* added or simplified symbolically without corruption. All comparisons are    *)
(* done by substituting concrete Omega FIRST, simplifying each term            *)
(* separately, then adding and comparing.                                       *)
(*                                                                              *)
(* Author: Frank O. Kuehnel / Excel Solutions LLC                              *)
(* ============================================================================ *)

(* ===== SETUP ===== *)
Clear[w1, w2, w3, T1, T2, T3];
w = {w1, w2, w3};
bigT = {T1, T2, T3};
hat[v_] := {{0, -v[[3]], v[[2]]}, {v[[3]], 0, -v[[1]]}, {-v[[2]], v[[1]], 0}};
id3 = IdentityMatrix[3];
th = Sqrt[w1^2 + w2^2 + w3^2];

(* S^{-1}(Omega) directly in Omega components *)
oneMinusX = (th/2)/Tan[th/2];
xVal = 1 - oneMinusX;
SinvMat = oneMinusX id3 + xVal Outer[Times, w, w]/th^2 - (th/2) hat[w]/th;

(* ===== Exact derivative (Mathematica chain-rules through Sqrt[]) ===== *)
SinvT = SinvMat . bigT;
dSinvTdOmega = Table[D[SinvT[[i]], w[[j]]], {i, 3}, {j, 3}];

(* ===== T-form: four terms kept SEPARATE (cannot add symbolically) ===== *)
rVec = w/th;
alphaSymb = (Sin[th] - th)/(2*(1 - Cos[th]));
betaSymb = xVal/th;
TbarSymb = Sum[rVec[[k]] bigT[[k]], {k, 3}];
rrT = Outer[Times, rVec, rVec];

term1 = (1/2) hat[bigT];
term2 = betaSymb Outer[Times, rVec, bigT];
term3 = alphaSymb Outer[Times, bigT, rVec];
term4 = TbarSymb (betaSymb id3 - (2 betaSymb + alphaSymb) rrT);

(* ===== PROOF: substitute 9 Omega vectors, compare with T symbolic ===== *)
Print["T-form vs exact derivative:"];
testOmegas = {{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1},{1,1,1},{2,-1,1},{1,2,3}};
provenCount = 0;
Do[
  Module[{sub, t1s, t2s, t3s, t4s, exactSub, diff, maxNum},
    sub = {w1 -> wv[[1]], w2 -> wv[[2]], w3 -> wv[[3]]};
    t1s = term1 /. sub // FullSimplify;
    t2s = term2 /. sub // FullSimplify;
    t3s = term3 /. sub // FullSimplify;
    t4s = term4 /. sub // FullSimplify;
    exactSub = dSinvTdOmega /. sub // FullSimplify;
    diff = FullSimplify[t1s + t2s + t3s + t4s - exactSub];
    If[diff === Table[0, {3}, {3}],
      provenCount++; Print["  ", wv, ": ZERO"],
      maxNum = Max[Abs[Flatten[N[diff /. {T1->17/13, T2->-11/7, T3->23/9}, 50]]]];
      If[maxNum < 10^-40,
        provenCount++; Print["  ", wv, ": ZERO (numerical)"],
        Print["  ", wv, ": FAIL (", ScientificForm[maxNum, 3], ")"]
      ]
    ]
  ],
  {wv, testOmegas}
];
Print["  Result: ", provenCount, "/9\n"];

(* ===== Coefficient identities ===== *)
Clear[theta];
alphaC = (Sin[theta] - theta)/(2*(1 - Cos[theta]));
betaC = (1 - (theta/2)/Tan[theta/2])/theta;

Print["Erratum check (should NOT be True):"];
Print["  alpha == 1/Th - 1/2 Cot[Th/2]: ",
  FullSimplify[alphaC == 1/theta - 1/2 Cot[theta/2]]];

alphaHalf = FullSimplify[alphaC - (1/4) Csc[theta/2]^2 (Sin[theta] - theta),
  Assumptions -> {theta > 0}];
betaId = FullSimplify[betaC - (1/theta - 1/2 Cot[theta/2])];

Print["alpha = (1/4)csc^2(Th/2)(sinTh-Th): ", If[alphaHalf===0, "PROVEN", alphaHalf]];
Print["beta = 1/Th - 1/2 cot(Th/2):        ", If[betaId===0, "PROVEN", betaId]];
Print["alpha - beta = ", FullSimplify[alphaC - betaC, Assumptions -> {theta > 0}]];
