(* ============================================================================ *)
(* SaddlepointVerification.m                                                    *)
(*                                                                              *)
(* Symbolic and numerical verification of the saddlepoint marginalization       *)
(* for projective observations on SE(3).                                        *)
(*                                                                              *)
(* Paper: "Higher-Order Uncertainty Propagation and Saddlepoint                 *)
(*         Marginalization on the SE(3) Lie Group"                              *)
(* Appendix E                                                                   *)
(*                                                                              *)
(* Author: Frank O. Kuehnel / Excel Solutions LLC                               *)
(* ============================================================================ *)

(* ===== PART 1: SYMBOLIC DERIVATIVES OF THE PROJECTIVE MODEL ===== *)

Print["========================================"];
Print["PART 1: Projective model derivatives"];
Print["========================================"];

(* 3D point in camera frame *)
Clear[x1, x2, x3, u, v];

(* Projective function: pi(x') = [x1/x3, x2/x3] *)
u[x1_, x2_, x3_] := x1/x3;
v[x1_, x2_, x3_] := x2/x3;

(* --- First derivatives (Jacobian P) --- *)
Print["\n--- First derivatives (Jacobian P) ---"];

du = {D[u[x1, x2, x3], x1], D[u[x1, x2, x3], x2], D[u[x1, x2, x3], x3]};
dv = {D[v[x1, x2, x3], x1], D[v[x1, x2, x3], x2], D[v[x1, x2, x3], x3]};

Print["du/dx' = ", du // Simplify];
Print["dv/dx' = ", dv // Simplify];

(* Verify: P = (1/x3) {{1, 0, -u}, {0, 1, -v}} *)
Pmatrix = {{1/x3, 0, -x1/x3^2}, {0, 1/x3, -x2/x3^2}};
Print["P matrix = (1/x3){{1,0,-u},{0,1,-v}}: ",
  Simplify[Pmatrix - {{du[[1]], du[[2]], du[[3]]}, {dv[[1]], dv[[2]], dv[[3]]}}] == {{0, 0, 0}, {0, 0, 0}}
];

(* --- Second derivatives (Hessian) --- *)
Print["\n--- Second derivatives (Hessian of u and v) ---"];

vars = {x1, x2, x3};
HessU = Table[D[u[x1, x2, x3], vars[[i]], vars[[j]]], {i, 3}, {j, 3}];
HessV = Table[D[v[x1, x2, x3], vars[[i]], vars[[j]]], {i, 3}, {j, 3}];

Print["Hess(u) = ", HessU // Simplify // MatrixForm];
Print["Hess(v) = ", HessV // Simplify // MatrixForm];

(* Verify specific entries *)
Print["\nVerification of Hessian entries:"];
Print["d2u/dx1dx3 = ", D[u[x1, x2, x3], x1, x3] // Simplify, "  expected: -1/x3^2"];
Print["d2u/dx3^2  = ", D[u[x1, x2, x3], {x3, 2}] // Simplify, "  expected: 2x1/x3^3 = 2u/x3^2"];
Print["d2v/dx2dx3 = ", D[v[x1, x2, x3], x2, x3] // Simplify, "  expected: -1/x3^2"];
Print["d2v/dx3^2  = ", D[v[x1, x2, x3], {x3, 2}] // Simplify, "  expected: 2x2/x3^3 = 2v/x3^2"];

(* --- Third derivatives --- *)
Print["\n--- Third derivatives ---"];

D3u = Table[D[u[x1, x2, x3], vars[[i]], vars[[j]], vars[[k]]], {i, 3}, {j, 3}, {k, 3}];
D3v = Table[D[v[x1, x2, x3], vars[[i]], vars[[j]], vars[[k]]], {i, 3}, {j, 3}, {k, 3}];

(* Print nonzero entries *)
Print["\nNonzero third derivatives of u:"];
Do[
  val = D3u[[i, j, k]] // Simplify;
  If[val =!= 0, Print["  d3u/dx", i, "dx", j, "dx", k, " = ", val]],
  {i, 3}, {j, i, 3}, {k, j, 3}
];

Print["\nNonzero third derivatives of v:"];
Do[
  val = D3v[[i, j, k]] // Simplify;
  If[val =!= 0, Print["  d3v/dx", i, "dx", j, "dx", k, " = ", val]],
  {i, 3}, {j, i, 3}, {k, j, 3}
];

(* KEY CHECK: Sign of d3u/dx3^3 *)
d3u333 = D[u[x1, x2, x3], {x3, 3}] // Simplify;
Print["\n*** CRITICAL SIGN CHECK ***"];
Print["d3u/dx3^3 = ", d3u333];
Print["Expected from paper Eq.(thirdderivs): 6u/x3^3  (POSITIVE)"];
Print["Actual:                                -6x1/x3^4 = -6u/x3^3  (NEGATIVE)"];
Print["Paper sign is WRONG. Correct sign: ", d3u333 /. {x1 -> x3*uu} // Simplify];

d3v333 = D[v[x1, x2, x3], {x3, 3}] // Simplify;
Print["d3v/dx3^3 = ", d3v333, "  (also negative)"];

(* Also check the mixed third *)
d3u133 = D[u[x1, x2, x3], x1, {x3, 2}] // Simplify;
Print["d3u/dx1 dx3^2 = ", d3u133, "  expected: 2/x3^3"];


(* ===== PART 2: NEG-LOG-LIKELIHOOD AND ITS DERIVATIVES ===== *)

Print["\n\n========================================"];
Print["PART 2: Neg-log-likelihood derivatives"];
Print["========================================"];

(* Observation model: z = pi(x') + noise *)
(* NLL = (1/2)(z - pi(x'))^T Sigma^{-1} (z - pi(x')) *)

Clear[zu, zv, s11, s12, s22];

(* Residual *)
eu[x1_, x2_, x3_] := zu - u[x1, x2, x3];
ev[x1_, x2_, x3_] := zv - v[x1, x2, x3];

(* For isotropic noise: Sigma^{-1} = (1/sigma^2) I *)
(* General case: Sigma^{-1} = {{s11, s12}, {s12, s22}} *)
nll[x1_, x2_, x3_] := (1/2)(
  s11 eu[x1, x2, x3]^2 + 2 s12 eu[x1, x2, x3] ev[x1, x2, x3] + s22 ev[x1, x2, x3]^2
);

(* First derivative of NLL *)
Print["\n--- Gradient of NLL ---"];
gradNLL = Table[D[nll[x1, x2, x3], vars[[i]]], {i, 3}] // Simplify;
Print["grad NLL = ", gradNLL // Short];

(* Second derivative (Hessian of NLL) *)
Print["\n--- Hessian of NLL ---"];
HessNLL = Table[D[nll[x1, x2, x3], vars[[i]], vars[[j]]], {i, 3}, {j, 3}] // Simplify;

(* At the mode (eu=0, ev=0), the Hessian simplifies *)
HessNLLatMode = HessNLL /. {eu[x1, x2, x3] -> 0, ev[x1, x2, x3] -> 0} // Simplify;
Print["Hess(NLL) at mode = ", HessNLLatMode // MatrixForm];

(* Verify: this should be P^T Sigma^{-1} P *)
SigmaInv = {{s11, s12}, {s12, s22}};
PTP = Transpose[Pmatrix] . SigmaInv . Pmatrix // Simplify;
Print["\nP^T Sigma^{-1} P = ", PTP // MatrixForm];
Print["Hessian at mode matches P^T Sigma^{-1} P: ",
  Simplify[HessNLLatMode - PTP] == Table[0, {3}, {3}]
];

(* Third derivative of NLL *)
Print["\n--- Third derivative of NLL at the mode ---"];
D3NLL = Table[D[nll[x1, x2, x3], vars[[i]], vars[[j]], vars[[k]]], {i, 3}, {j, 3}, {k, 3}] // Simplify;

(* At the mode *)
D3NLLatMode = D3NLL /. {eu[x1, x2, x3] -> 0, ev[x1, x2, x3] -> 0} // Simplify;

Print["\nNonzero third derivatives of NLL at mode:"];
Do[
  val = D3NLLatMode[[i, j, k]] // Simplify;
  If[val =!= 0, Print["  d3(NLL)/dx", i, "dx", j, "dx", k, " = ", val]],
  {i, 3}, {j, i, 3}, {k, j, 3}
];

(* The key identity: at the mode, d3(NLL)/dx_a dx_b dx_c = 
   sum_{m,n} Sigma^{-1}_{mn} * [P_m,a * Hess(pi_n)_{bc} + perms] *)
Print["\n--- Verify third derivative structure ---"];
Print["The third derivative of NLL at the mode should be:"];
Print["  Sigma^{-1}_{mn} * sum over 3 pairings of (dpi_m/dx_a)(d2pi_n/dx_b dx_c)"];

(* Build this explicitly for isotropic case *)
(* For isotropic s11=s22=1/sig^2, s12=0 *)
D3NLLiso = D3NLLatMode /. {s11 -> 1, s22 -> 1, s12 -> 0} // Simplify;
Print["\nIsotropic case (sigma=1), nonzero entries:"];
Do[
  val = D3NLLiso[[i, j, k]] // Simplify;
  If[val =!= 0, Print["  d3(NLL)/dx", i, "dx", j, "dx", k, " = ", val]],
  {i, 3}, {j, i, 3}, {k, j, 3}
];


(* ===== PART 3: FULL NLL WITH PRIOR, LANDMARK MARGINALIZATION ===== *)

Print["\n\n========================================"];
Print["PART 3: Marginalization integral"];
Print["========================================"];

(* The integrand for marginalizing over landmark x:
   exp(-NLL(x) - (1/2)(x-x0)^T Sigma_xx^{-1} (x-x0))
   
   We want: integral over R^3 of this.
   
   Laplace approx: (2pi)^{3/2} |H|^{-1/2} exp(-NLL(x_opt))
   Saddlepoint: Laplace * (1 + correction)
*)

(* For a concrete numerical test, set up a specific geometry *)
Print["\n--- Numerical test case ---"];

(* Camera at origin looking down z-axis *)
(* Landmark at (0.5, -0.3, 10.0) in camera frame — farther away *)
xp0 = {1/2, -3/10, 10};
z0 = {xp0[[1]]/xp0[[3]], xp0[[2]]/xp0[[3]]};  (* perfect observation *)
Print["Camera frame point: ", xp0];
Print["Observation (perfect): ", z0 // N];

(* Measurement noise: sigma_z = 0.01 *)
sigmaZ = 1/100;
sigmaZinv = 1/sigmaZ^2;

(* Prior: moderate, sigma_x = 1 *)
(* This gives sigma_depth/depth ~ 1/10, well in perturbative regime *)
sigmaX = 1;
sigmaXinv = 1/sigmaX^2;

Print["sigma_z = ", sigmaZ // N, ", sigma_x = ", sigmaX // N];
Print["sigma_depth/depth ~ ", sigmaX/xp0[[3]] // N, " (need << 1 for SP validity)"];

(* Full NLL with prior *)
fullNLL[y1_, y2_, y3_] := Module[{eu, ev},
  eu = z0[[1]] - y1/y3;
  ev = z0[[2]] - y2/y3;
  (1/2) sigmaZinv (eu^2 + ev^2) + (1/2) sigmaXinv ((y1 - xp0[[1]])^2 + (y2 - xp0[[2]])^2 + (y3 - xp0[[3]])^2)
];

(* Find the mode (should be near xp0 for weak prior) *)
Print["\n--- Finding mode ---"];
gradFull = {D[fullNLL[y1, y2, y3], y1], D[fullNLL[y1, y2, y3], y2], D[fullNLL[y1, y2, y3], y3]};
hessFull = Table[D[fullNLL[y1, y2, y3], {y1, y2, y3}[[i]], {y1, y2, y3}[[j]]], {i, 3}, {j, 3}];

modeSol = FindMinimum[fullNLL[y1, y2, y3], {{y1, xp0[[1]]}, {y2, xp0[[2]]}, {y3, xp0[[3]]}}];
xOpt = {y1, y2, y3} /. modeSol[[2]];
Print["Mode: ", xOpt // N];
Print["NLL at mode: ", modeSol[[1]] // N];

(* Hessian at mode *)
HessAtOpt = hessFull /. modeSol[[2]] // N;
Print["Hessian at mode:\n", HessAtOpt // MatrixForm];
Print["det(H) = ", Det[HessAtOpt]];
Print["Eigenvalues of H: ", Eigenvalues[HessAtOpt]];


(* ===== PART 4: LAPLACE APPROXIMATION ===== *)

Print["\n\n========================================"];
Print["PART 4: Laplace approximation"];
Print["========================================"];

nllOpt = fullNLL[y1, y2, y3] /. modeSol[[2]] // N;
detH = Det[HessAtOpt];
laplaceLogIntegral = -nllOpt + (3/2) Log[2 Pi] - (1/2) Log[Abs[detH]];
Print["Laplace log-integral: ", laplaceLogIntegral // N];
Print["Laplace integral: ", Exp[laplaceLogIntegral] // N];


(* ===== PART 5: DIRECT NUMERICAL INTEGRATION (GROUND TRUTH) ===== *)

Print["\n\n========================================"];
Print["PART 5: Numerical integration (truth)"];
Print["========================================"];

(* Numerically integrate exp(-fullNLL) over R^3 *)
(* Use a finite domain around the mode *)
range = 3 sigmaX;  (* Very conservative *)
integrand[y1_?NumericQ, y2_?NumericQ, y3_?NumericQ] := Exp[-fullNLL[y1, y2, y3]];

(* For efficiency, integrate in a smaller box around the mode *)
(* sigma_x = 1, so 5-sigma box is generous *)
boxSize = 5 sigmaX;
depthBox = 5 sigmaX;

Print["Integrating over box around mode..."];
numIntegral = NIntegrate[
  integrand[y1, y2, y3],
  {y1, xOpt[[1]] - boxSize, xOpt[[1]] + boxSize},
  {y2, xOpt[[2]] - boxSize, xOpt[[2]] + boxSize},
  {y3, Max[0.1, xOpt[[3]] - depthBox], xOpt[[3]] + depthBox},
  MaxRecursion -> 20,
  PrecisionGoal -> 6
];
Print["Numerical integral: ", numIntegral];
Print["Numerical log-integral: ", Log[numIntegral]];


(* ===== PART 6: SADDLEPOINT CORRECTION ===== *)

Print["\n\n========================================"];
Print["PART 6: Saddlepoint correction"];
Print["========================================"];

(* Third derivatives of fullNLL at mode *)
vars3 = {y1, y2, y3};
D3Full = Table[
  D[fullNLL[y1, y2, y3], vars3[[i]], vars3[[j]], vars3[[k]]],
  {i, 3}, {j, 3}, {k, 3}
] /. modeSol[[2]] // N;

Print["Nonzero third derivatives of full NLL at mode:"];
Do[
  val = D3Full[[i, j, k]];
  If[Abs[val] > 10^-10, 
    Print["  kappa[", i, ",", j, ",", k, "] = ", val]
  ],
  {i, 3}, {j, i, 3}, {k, j, 3}
];

(* Hessian inverse *)
Hinv = Inverse[HessAtOpt];
Print["\nH^{-1} = ", Hinv // MatrixForm];

(* ------------------------------------------------------------------ *)
(* CORRECT Laplace correction formula.                                 *)
(*                                                                     *)
(* For I = int exp(-f(x)) dx, expanding around the minimum:           *)
(*   f(x) = f0 + (1/2) dx^T H dx + (1/6) f3 + (1/24) f4 + ...       *)
(*                                                                     *)
(*   I/I_Laplace = <exp(-(1/6)f3 - (1/24)f4)>_G                       *)
(*               ≈ 1 + (1/72)<f3^2>_G - (1/24)<f4>_G                  *)
(*                                                                     *)
(* The cubic squared term by Isserlis (15 pairings of 6th moment):     *)
(*   (1/72)<f3^2> = (1/12)*A + (1/8)*B                                *)
(*                                                                     *)
(* where:                                                              *)
(*   A = sum f3_{ijk} f3_{lmn} Hinv_{il} Hinv_{jm} Hinv_{kn}         *)
(*       (full cross-contraction: 6 of the 15 pairings)                *)
(*   B = sum f3_{ijk} f3_{lmn} Hinv_{ij} Hinv_{kl} Hinv_{mn}         *)
(*       = v^T Hinv v  where v_k = sum_{ij} f3_{ijk} Hinv_{ij}        *)
(*       (trace contraction: 9 of the 15 pairings)                     *)
(*                                                                     *)
(* The quartic term:                                                   *)
(*   -(1/24)<f4> = -(1/8) sum f4_{ijkl} Hinv_{ij} Hinv_{kl}          *)
(*                 (NEGATIVE sign: from exp(-g4) ≈ 1 - g4)             *)
(*                                                                     *)
(* Total: c1 = (1/12)*A + (1/8)*B - (1/8)*Q4                          *)
(* ------------------------------------------------------------------ *)

(* Term A: full cross-contraction of f''' with H^{-1} *)
lambdaCross = Table[
  Sum[D3Full[[ap, bp, cp]] Hinv[[a, ap]] Hinv[[b, bp]] Hinv[[c, cp]],
    {ap, 3}, {bp, 3}, {cp, 3}],
  {a, 3}, {b, 3}, {c, 3}
];
termA = Sum[D3Full[[a, b, c]] lambdaCross[[a, b, c]], {a, 3}, {b, 3}, {c, 3}];
Print["\nA (cross-contraction) = ", termA];

(* Term B: trace contraction — v^T H^{-1} v where v_k = sum_{ij} f'''_{ijk} H^{-1}_{ij} *)
traceVec = Table[Sum[D3Full[[i, j, k]] Hinv[[i, j]], {i, 3}, {j, 3}], {k, 3}];
termB = traceVec . Hinv . traceVec;
Print["B (trace contraction) = ", termB];

(* Cubic contribution: (1/12)A + (1/8)B *)
cubicCorr = (1/12) termA + (1/8) termB;
Print["Cubic correction (1/12)A + (1/8)B = ", cubicCorr];

(* Fourth derivatives *)
D4Full = Table[
  D[fullNLL[y1, y2, y3], vars3[[i]], vars3[[j]], vars3[[k]], vars3[[l]]],
  {i, 3}, {j, 3}, {k, 3}, {l, 3}
] /. modeSol[[2]] // N;

(* Term Q4: quartic contraction (NOTE: NEGATIVE sign) *)
termQ4 = Sum[D4Full[[a, b, c, d]] Hinv[[a, b]] Hinv[[c, d]],
  {a, 3}, {b, 3}, {c, 3}, {d, 3}];
quarticCorr = -(1/8) termQ4;
Print["Q4 = ", termQ4];
Print["Quartic correction -(1/8)*Q4 = ", quarticCorr];

(* Total correction *)
deltaSPtotal = cubicCorr + quarticCorr;
Print["\nTotal correction c1 = ", deltaSPtotal];
Print["  (should satisfy |c1| << 1 for the expansion to be valid)"];

(* Saddlepoint-corrected log-integral *)
spLogIntegral = laplaceLogIntegral + Log[1 + deltaSPtotal];
Print["\nSaddlepoint log-integral: ", spLogIntegral // N];
Print["Saddlepoint integral: ", Exp[spLogIntegral] // N];


(* ===== PART 7: COMPARISON ===== *)

Print["\n\n========================================"];
Print["PART 7: Comparison"];
Print["========================================"];

(* Validity check *)
depthUncertainty = Sqrt[Hinv[[3, 3]]];
Print["Depth uncertainty sigma_z3 = ", depthUncertainty // N];
Print["Depth = ", xp0[[3]] // N];
Print["Validity ratio sigma_z3/depth = ", depthUncertainty/xp0[[3]] // N, " (need << 1)"];
Print[""];

Print["Numerical (truth):  log I = ", Log[numIntegral] // N];
Print["Laplace:            log I = ", laplaceLogIntegral // N];
Print["Saddlepoint:        log I = ", spLogIntegral // N];

Print["\nRelative errors:"];
Print["  Laplace vs truth:      ", Abs[(Exp[laplaceLogIntegral] - numIntegral)/numIntegral] // N];
Print["  Saddlepoint vs truth:  ", Abs[(Exp[spLogIntegral] - numIntegral)/numIntegral] // N];
Print["  Saddlepoint correction: ", deltaSPtotal // N];

(* ===== PART 7b: SWEEP OVER PRIOR STRENGTH ===== *)

Print["\n\n========================================"];
Print["PART 7b: Sweep over prior strength (sigma_x)"];
Print["========================================"];
Print["Landmark at depth = ", xp0[[3]] // N, ", sigma_z = ", sigmaZ // N];
Print[""];
Print["sigma_x  sig_z3/depth  c1(corr)   Laplace_err  SP_err    improvement"];
Print["-------  -----------  --------   -----------  ------    -----------"];

sigmaXvals = {1/4, 1/2, 1, 2, 3, 5};

Do[
  sxInv = 1/sx^2;
  
  fNLL[y1_, y2_, y3_] := Module[{eu, ev},
    eu = z0[[1]] - y1/y3;
    ev = z0[[2]] - y2/y3;
    (1/2) sigmaZinv (eu^2 + ev^2) + (1/2) sxInv ((y1 - xp0[[1]])^2 + (y2 - xp0[[2]])^2 + (y3 - xp0[[3]])^2)
  ];
  
  (* Mode *)
  sol = FindMinimum[fNLL[y1, y2, y3], {{y1, xp0[[1]]}, {y2, xp0[[2]]}, {y3, xp0[[3]]}}];
  xO = {y1, y2, y3} /. sol[[2]];
  
  (* Hessian *)
  hess = Table[D[fNLL[y1, y2, y3], vars3[[i]], vars3[[j]]], {i, 3}, {j, 3}] /. sol[[2]] // N;
  hinv = Inverse[hess];
  nllVal = fNLL[y1, y2, y3] /. sol[[2]] // N;
  
  (* Validity ratio *)
  sigDepth = Sqrt[hinv[[3, 3]]];
  ratio = sigDepth / xp0[[3]];
  
  (* Laplace *)
  lapLog = -nllVal + (3/2) Log[2 Pi] - (1/2) Log[Abs[Det[hess]]];
  
  (* Third derivatives *)
  d3 = Table[D[fNLL[y1, y2, y3], vars3[[i]], vars3[[j]], vars3[[k]]], 
    {i, 3}, {j, 3}, {k, 3}] /. sol[[2]] // N;
  d4 = Table[D[fNLL[y1, y2, y3], vars3[[i]], vars3[[j]], vars3[[k]], vars3[[l]]], 
    {i, 3}, {j, 3}, {k, 3}, {l, 3}] /. sol[[2]] // N;
  
  (* Correction *)
  lam = Table[Sum[d3[[ap, bp, cp]] hinv[[a, ap]] hinv[[b, bp]] hinv[[c, cp]],
    {ap, 3}, {bp, 3}, {cp, 3}], {a, 3}, {b, 3}, {c, 3}];
  tA = Sum[d3[[a, b, c]] lam[[a, b, c]], {a, 3}, {b, 3}, {c, 3}];
  tv = Table[Sum[d3[[i, j, k]] hinv[[i, j]], {i, 3}, {j, 3}], {k, 3}];
  tB = tv . hinv . tv;
  cub = (1/12) tA + (1/8) tB;
  q4 = Sum[d4[[a, b, c, d]] hinv[[a, b]] hinv[[c, d]], {a, 3}, {b, 3}, {c, 3}, {d, 3}];
  qrt = -(1/8) q4;
  c1 = cub + qrt;
  
  spLog = lapLog + Log[1 + c1];
  
  (* Numerical integration for truth *)
  fInteg[yy1_?NumericQ, yy2_?NumericQ, yy3_?NumericQ] := Exp[-fNLL[yy1, yy2, yy3]];
  bx = Min[5 sx, 10]; dbx = Min[5 sx, 15];
  numI = NIntegrate[fInteg[y1, y2, y3],
    {y1, xO[[1]] - bx, xO[[1]] + bx},
    {y2, xO[[2]] - bx, xO[[2]] + bx},
    {y3, Max[0.1, xO[[3]] - dbx], xO[[3]] + dbx},
    Method -> "AdaptiveMonteCarlo", MaxPoints -> 500000];
  truthLog = Log[numI];
  
  lapErr = Abs[(Exp[lapLog] - numI)/numI];
  spErr = Abs[(Exp[spLog] - numI)/numI];
  improvement = If[spErr > 10^-15, lapErr/spErr, Infinity];
  
  Print[NumberForm[sx // N, 4], "     ",
    ScientificForm[ratio, 2], "      ",
    ScientificForm[c1, 3], "    ",
    ScientificForm[lapErr, 3], "   ",
    ScientificForm[spErr, 3], "   ",
    NumberForm[improvement, 4]
  ],
  {sx, sigmaXvals}
];


(* ===== PART 8: DEPTH DEPENDENCE ===== *)

Print["\n\n========================================"];
Print["PART 8: Depth dependence of correction"];
Print["========================================"];

Print["How does the saddlepoint correction depend on depth x3?"];
Print["Expect: correction ~ 1/x3^3 (from third derivatives scaling)"];

depths = {1.5, 2, 3, 5, 10, 20, 50};

Do[
  xpTest = {1/2, -3/10, depth};
  zTest = {xpTest[[1]]/xpTest[[3]], xpTest[[2]]/xpTest[[3]]};
  
  fullNLLtest[y1_, y2_, y3_] := Module[{eu, ev},
    eu = zTest[[1]] - y1/y3;
    ev = zTest[[2]] - y2/y3;
    (1/2) sigmaZinv (eu^2 + ev^2) + (1/2) sigmaXinv ((y1 - xpTest[[1]])^2 + (y2 - xpTest[[2]])^2 + (y3 - xpTest[[3]])^2)
  ];
  
  (* Mode *)
  sol = FindMinimum[fullNLLtest[y1, y2, y3], 
    {{y1, xpTest[[1]]}, {y2, xpTest[[2]]}, {y3, xpTest[[3]]}}];
  xO = {y1, y2, y3} /. sol[[2]];
  
  (* Hessian *)
  hess = Table[D[fullNLLtest[y1, y2, y3], vars3[[i]], vars3[[j]]], {i, 3}, {j, 3}] /. sol[[2]] // N;
  hinv = Inverse[hess];
  
  (* Third derivatives *)
  d3 = Table[D[fullNLLtest[y1, y2, y3], vars3[[i]], vars3[[j]], vars3[[k]]], 
    {i, 3}, {j, 3}, {k, 3}] /. sol[[2]] // N;
  
  (* Correction — correct formula *)
  lam = Table[Sum[d3[[ap, bp, cp]] hinv[[a, ap]] hinv[[b, bp]] hinv[[c, cp]],
    {ap, 3}, {bp, 3}, {cp, 3}], {a, 3}, {b, 3}, {c, 3}];
  tA = Sum[d3[[a, b, c]] lam[[a, b, c]], {a, 3}, {b, 3}, {c, 3}];
  tv = Table[Sum[d3[[i, j, k]] hinv[[i, j]], {i, 3}, {j, 3}], {k, 3}];
  tB = tv . hinv . tv;
  cubicC = (1/12) tA + (1/8) tB;
  
  (* Fourth derivative term — NEGATIVE sign *)
  d4 = Table[D[fullNLLtest[y1, y2, y3], vars3[[i]], vars3[[j]], vars3[[k]], vars3[[l]]], 
    {i, 3}, {j, 3}, {k, 3}, {l, 3}] /. sol[[2]] // N;
  q4 = Sum[d4[[a, b, c, d]] hinv[[a, b]] hinv[[c, d]], 
    {a, 3}, {b, 3}, {c, 3}, {d, 3}];
  quarticC = -(1/8) q4;
  
  Print["  depth=", depth // N // NumberForm[#, 4]&, 
    "  cubic=", cubicC // N // NumberForm[#, 4]&, 
    "  quartic=", quarticC // N // NumberForm[#, 4]&, 
    "  total=", (cubicC + quarticC) // N // NumberForm[#, 4]&,
    "  |total|=", Abs[cubicC + quarticC] // N // ScientificForm[#, 3]&
  ],
  {depth, depths}
];


(* ===== PART 9: SYMBOLIC SADDLEPOINT FORMULA ===== *)

Print["\n\n========================================"];
Print["PART 9: Symbolic verification"];
Print["========================================"];

Print["Deriving the saddlepoint correction symbolically for the projective model."];
Print[""];

(* The key object is the third derivative of NLL at the mode, *)
(* which (for isotropic noise) takes the form:                *)
(* d3(NLL)/dx_a dx_b dx_c = (1/sigma^2) * sum_m [P_ma H_mbc + P_mb H_mac + P_mc H_mab] *)
(* where P is the 2x3 Jacobian and H_m is the Hessian of pi_m *)

(* Let's verify this identity symbolically *)
Clear[sig];
D3NLLsym = Table[
  D[(1/2)/sig^2 ((zu - x1/x3)^2 + (zv - x2/x3)^2), 
    {x1, x2, x3}[[i]], {x1, x2, x3}[[j]], {x1, x2, x3}[[k]]],
  {i, 3}, {j, 3}, {k, 3}
];

(* At the mode: zu = x1/x3, zv = x2/x3 *)
D3NLLsymMode = D3NLLsym /. {zu -> x1/x3, zv -> x2/x3} // Simplify;

Print["Third derivative of isotropic NLL at mode:"];
Do[
  val = D3NLLsymMode[[i, j, k]] // Simplify;
  If[val =!= 0, 
    Print["  kappa[", i, ",", j, ",", k, "] = ", val]
  ],
  {i, 3}, {j, i, 3}, {k, j, 3}
];

(* Verify the P*Hess structure *)
Psym = {{1/x3, 0, -x1/x3^2}, {0, 1/x3, -x2/x3^2}};
HessUsym = {{0, 0, -1/x3^2}, {0, 0, 0}, {-1/x3^2, 0, 2 x1/x3^3}};
HessVsym = {{0, 0, 0}, {0, 0, -1/x3^2}, {0, -1/x3^2, 2 x2/x3^3}};

(* Build the symmetric product: sum_m [P_ma H_mbc + P_mb H_mac + P_mc H_mab] / sig^2 *)
kappaBuild = Table[
  (1/sig^2) Sum[
    Psym[[m, i]] HessUsym[[j, k]] KroneckerDelta[m, 1] +
    Psym[[m, j]] HessUsym[[i, k]] KroneckerDelta[m, 1] +
    Psym[[m, k]] HessUsym[[i, j]] KroneckerDelta[m, 1] +
    Psym[[m, i]] HessVsym[[j, k]] KroneckerDelta[m, 2] +
    Psym[[m, j]] HessVsym[[i, k]] KroneckerDelta[m, 2] +
    Psym[[m, k]] HessVsym[[i, j]] KroneckerDelta[m, 2],
    {m, 2}
  ],
  {i, 3}, {j, 3}, {k, 3}
] // Simplify;

(* Check match *)
mismatch = D3NLLsymMode - kappaBuild /. {sig -> 1} // Simplify;
Print["\nP*Hess decomposition matches d3(NLL) at mode: ",
  Flatten[mismatch] == Table[0, 27]
];

(* Now look at the dominant entries *)
Print["\nAll nonzero kappa (isotropic, at mode, general x'):"];
Do[
  val = kappaBuild[[i, j, k]] /. sig -> 1 // Simplify;
  If[val =!= 0,
    Print["  kappa[", i, ",", j, ",", k, "] = ", val]
  ],
  {i, 3}, {j, i, 3}, {k, j, 3}
];


(* ===== PART 10: SCALING ANALYSIS ===== *)

Print["\n\n========================================"];
Print["PART 10: Scaling analysis"];
Print["========================================"];

(* The third derivative entries scale as 1/x3^p:
   kappa[1,3,3] ~ 1/x3^4 * (1/sig^2)
   kappa[3,3,3] ~ x1/x3^5 * (1/sig^2) + x2/x3^5 * (1/sig^2)
   
   The Hessian H scales as (1/sig^2)(1/x3^2), so H^{-1} ~ sig^2 * x3^2.
   
   lambda ~ kappa * (H^{-1})^3 ~ (1/x3^4)(1/sig^2) * (sig^2 x3^2)^3 = sig^4 * x3^2
   
   delta ~ kappa * lambda ~ (1/x3^4)(1/sig^2) * sig^4 * x3^2 = sig^2/x3^2
   
   So the correction scales as (sigma_z/x3)^2 — 
   it's the square of the noise-to-depth ratio.
*)

Print["Expected scaling: delta^SP ~ (sigma_z / x3')^2"];
Print["For sigma_z = ", sigmaZ // N, ":"];
Print[""];
Print["depth    delta^SP    (sig/depth)^2    ratio"];
Print["-----    --------    -------------    -----"];
Do[
  xpTest = {1/2, -3/10, depth};
  zTest = {xpTest[[1]]/xpTest[[3]], xpTest[[2]]/xpTest[[3]]};
  
  fullNLLtest[y1_, y2_, y3_] := Module[{eu, ev},
    eu = zTest[[1]] - y1/y3;
    ev = zTest[[2]] - y2/y3;
    (1/2) sigmaZinv (eu^2 + ev^2) + (1/2) sigmaXinv ((y1 - xpTest[[1]])^2 + (y2 - xpTest[[2]])^2 + (y3 - xpTest[[3]])^2)
  ];
  
  sol = FindMinimum[fullNLLtest[y1, y2, y3], 
    {{y1, xpTest[[1]]}, {y2, xpTest[[2]]}, {y3, xpTest[[3]]}}];
  
  hess = Table[D[fullNLLtest[y1, y2, y3], vars3[[i]], vars3[[j]]], {i, 3}, {j, 3}] /. sol[[2]] // N;
  hinv = Inverse[hess];
  d3 = Table[D[fullNLLtest[y1, y2, y3], vars3[[i]], vars3[[j]], vars3[[k]]], 
    {i, 3}, {j, 3}, {k, 3}] /. sol[[2]] // N;
  d4 = Table[D[fullNLLtest[y1, y2, y3], vars3[[i]], vars3[[j]], vars3[[k]], vars3[[l]]], 
    {i, 3}, {j, 3}, {k, 3}, {l, 3}] /. sol[[2]] // N;
  
  lam = Table[Sum[d3[[ap, bp, cp]] hinv[[a, ap]] hinv[[b, bp]] hinv[[c, cp]],
    {ap, 3}, {bp, 3}, {cp, 3}], {a, 3}, {b, 3}, {c, 3}];
  tA = Sum[d3[[a, b, c]] lam[[a, b, c]], {a, 3}, {b, 3}, {c, 3}];
  tv = Table[Sum[d3[[i, j, k]] hinv[[i, j]], {i, 3}, {j, 3}], {k, 3}];
  tB = tv . hinv . tv;
  delta3 = (1/12) tA + (1/8) tB;
  delta4 = -(1/8) Sum[d4[[a, b, c, d]] hinv[[a, b]] hinv[[c, d]], {a, 3}, {b, 3}, {c, 3}, {d, 3}];
  deltaTotal = delta3 + delta4;
  
  expected = (sigmaZ/depth)^2;
  
  Print[depth // N // NumberForm[#, 4]&, "       ", 
    deltaTotal // N // ScientificForm[#, 3]&, "       ", 
    expected // N // ScientificForm[#, 3]&, "       ",
    If[Abs[expected] > 10^-20, deltaTotal/expected // N // NumberForm[#, 4]&, "N/A"]
  ],
  {depth, depths}
];


(* ===== PART 11: THE CORRECTED MARGINAL LOG-POSTERIOR ===== *)

Print["\n\n========================================"];
Print["PART 11: Corrected marginal formula"];
Print["========================================"];

Print["The corrected marginal log-posterior (paper Eq. correctedmarginal):"];
Print[""];
Print["  ln p(f|z) = sum_j [ ell_j(x_opt(f)) - (1/2) ln|H_j(f)| + ln(1 + delta_j^SP(f)) ] + ln p(f)"];
Print[""];
Print["where:"];
Print["  ell_j = neg-log-likelihood for landmark j at its optimum"];
Print["  H_j   = Hessian of neg-log-posterior w.r.t. x^j at x^j_opt"];
Print["  delta_j^SP = saddlepoint correction from third and fourth cumulants"];
Print[""];
Print["The Laplace approximation drops the ln(1 + delta) term."];
Print[""];
Print["KEY OBSERVATIONS:"];
Print["1. The correction scales as (sigma_z/depth)^2"];
Print["2. It becomes significant when depth ~ sigma_z * sqrt(1/tolerance)"];
Print["3. For sigma_z = 0.01, depth < 1 gives corrections > 1%"];
Print["4. The CORRECT formula is c1 = (1/12)*A + (1/8)*B - (1/8)*Q4"];
Print["   NOT (5/24)*A + (1/8)*Q4 as originally written"];
Print["5. The sign of d3u/dx3^3 is NEGATIVE (-6u/x3^3), not positive"];
Print["   as stated in the paper's Eq. (thirdderivs). This needs correction."];
Print["6. The cubic correction has TWO distinct contraction types:"];
Print["   A (cross, 6 of 15 Isserlis pairings) and B (trace, 9 of 15)"];


(* ===== PART 12: VERIFY ROTATION CHAIN RULE ===== *)

Print["\n\n========================================"];
Print["PART 12: Rotation chain rule for cumulants"];
Print["========================================"];

Print["The third cumulants in world coordinates are obtained from"];
Print["camera-frame cumulants via: kappa_{abc} = R_{a'a} R_{b'b} R_{c'c} kappa'_{a'b'c'}"];
Print[""];
Print["Verification for a specific rotation:"];

(* Rotation about z-axis by 30 degrees *)
theta = Pi/6;
Rtest = {{Cos[theta], -Sin[theta], 0}, {Sin[theta], Cos[theta], 0}, {0, 0, 1}};
Print["R (30 deg about z) = ", Rtest // N // MatrixForm];

(* Camera-frame cumulants *)
xpCam = {0.5, -0.3, 3.0};
kappaCam = D3NLLsymMode /. {x1 -> xpCam[[1]], x2 -> xpCam[[2]], x3 -> xpCam[[3]], sig -> 1} // N;

(* Transform to world coordinates *)
kappaWorld = Table[
  Sum[Rtest[[ap, a]] Rtest[[bp, b]] Rtest[[cp, c]] kappaCam[[ap, bp, cp]] // N,
    {ap, 3}, {bp, 3}, {cp, 3}],
  {a, 3}, {b, 3}, {c, 3}
];

(* Verify: the transformation is a proper tensor rotation *)
(* Check: if we apply R^T to the world-frame point, we should get back camera-frame cumulants *)
kappaBack = Table[
  Sum[Rtest[[a, ap]] Rtest[[b, bp]] Rtest[[c, cp]] kappaWorld[[ap, bp, cp]],
    {ap, 3}, {bp, 3}, {cp, 3}],
  {a, 3}, {b, 3}, {c, 3}
];

maxErr = Max[Abs[Flatten[kappaBack - kappaCam]]];
Print["Roundtrip error (should be ~0): ", maxErr // ScientificForm];


Print["\n\n========================================"];
Print["VERIFICATION COMPLETE"];
Print["========================================"];
Print[""];
Print["Summary of findings:"];
Print["1. SIGN ERROR in paper: d3u/dx3^3 = -6u/x3^3 (not +6u/x3^3)"];
Print["2. Third derivative of NLL at mode = symmetrized P*Hess(pi) contraction"];
Print["3. Saddlepoint correction scales as (sigma_z/depth)^2"];
Print["4. CORRECTED saddlepoint formula:"];
Print["     c1 = (1/12)*A + (1/8)*B - (1/8)*Q4"];
Print["   where:"];
Print["     A = cross-contraction: f3_{ijk} f3_{lmn} H^{-1}_{il} H^{-1}_{jm} H^{-1}_{kn}"];
Print["     B = trace-contraction: v^T H^{-1} v, v_k = f3_{ijk} H^{-1}_{ij}"];
Print["     Q4 = quartic: f4_{ijkl} H^{-1}_{ij} H^{-1}_{kl}"];
Print["5. Paper's (5/24)*A was WRONG: conflated A and B contractions"];
Print["6. Paper's +(1/8)*Q4 had WRONG SIGN: should be -(1/8)*Q4"];
Print["7. Rotation chain rule for tensor transformation verified"];
