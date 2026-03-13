(* ============================================================================ *)
(* CouplingJacobianDerivation.m                                                *)
(*                                                                              *)
(* Symbolic derivation and verification of d[S^{-1}(Omega) T]/dOmega          *)
(* from first principles, settling the correct closed-form formula.            *)
(*                                                                              *)
(* Author: Frank O. Kuehnel / Excel Solutions LLC                              *)
(* ============================================================================ *)

Print["================================================================"];
Print["Symbolic derivation of d[S^{-1}(Omega) T] / dOmega"];
Print["================================================================"];

(* ===== SETUP ===== *)
Clear[theta, r1, r2, r3, T1, T2, T3];

r = {r1, r2, r3};
bigT = {T1, T2, T3};

hat[w_] := {{0, -w[[3]], w[[2]]}, {w[[3]], 0, -w[[1]]}, {-w[[2]], w[[1]], 0}};
id3 = IdentityMatrix[3];
H = hat[r];
H2 = H . H // Simplify;
rrT = Outer[Times, r, r];

(* Constraint: |r| = 1 *)
unitrule = r1^2 + r2^2 + r3^2 -> 1;

(* S^{-1}(Omega) = (1-x) I + x r r^T - (theta/2) H *)
(* where 1-x = (theta/2) cot(theta/2) *)
oneMinusXsymb = (theta/2) / Tan[theta/2];
xVal = 1 - oneMinusXsymb;

Sinv = oneMinusXsymb id3 + xVal rrT - (theta/2) H;

Print["\nS^{-1} = (1-x)I + x r r^T - (theta/2) H"];
Print["where 1-x = (theta/2) cot(theta/2)"];

(* ===== STEP 1: Compute S^{-1} T symbolically ===== *)
SinvT = Sinv . bigT // Simplify;
Print["\nS^{-1} T = ", SinvT];


(* ===== STEP 2: Differentiate w.r.t. Omega_j = theta * r_j ===== *)
(*
   We need d[S^{-1} T]/dOmega_j where Omega = theta * r.
   Using chain rule: d/dOmega_j = (dr_k/dOmega_j)(d/dr_k) + (dtheta/dOmega_j)(d/dtheta)
   
   Since Omega = theta * r, with |r|=1:
   - dtheta/dOmega_j = r_j  (theta = |Omega|, so dtheta/dOmega_j = Omega_j/theta = r_j)
   - dr_k/dOmega_j = (delta_{kj} - r_k r_j) / theta
   
   So: d/dOmega_j = r_j (d/dtheta) + (1/theta)*(delta_{kj} - r_k r_j)(d/dr_k)
*)

Print["\n================================================================"];
Print["STEP 2: Symbolic differentiation"];
Print["================================================================"];

(* Derivative of S^{-1} T w.r.t. theta *)
dSinvTdtheta = D[SinvT, theta];

(* Derivative of S^{-1} T w.r.t. r_k *)
dSinvTdr = Table[D[SinvT, r[[k]]], {k, 3}];

(* Assemble: [d(S^{-1}T)/dOmega]_{i,j} = r_j [d(S^{-1}T)/dtheta]_i 
                                           + (1/theta) sum_k (delta_{kj} - r_k r_j) [d(S^{-1}T)/dr_k]_i *)

dSinvTdOmega = Table[
  r[[j]] dSinvTdtheta[[i]] + (1/theta) Sum[(KroneckerDelta[k, j] - r[[k]] r[[j]]) dSinvTdr[[k, i]], {k, 3}],
  {i, 3}, {j, 3}
];

(* Simplify with |r|=1 *)
dSinvTdOmega = dSinvTdOmega /. unitrule // FullSimplify;

Print["Symbolic d[S^{-1}T]/dOmega computed (3x3 matrix)."];
Print["(Expression is large; will verify numerically below.)"];


(* ===== STEP 3: Numerical verification at test point ===== *)

Print["\n================================================================"];
Print["STEP 3: Numerical verification"];
Print["================================================================"];

thetaNum = Sqrt[0.5^2 + 0.3^2 + 0.7^2] // N;
rNum = {0.5, -0.3, 0.7}/thetaNum;
tExp = {1.0, -0.5, 0.3};
omegaNum = thetaNum rNum;

(* Compute T = S(Omega) t *)
rodriguesExp[omega_] := Module[{th = Norm[omega], rr, HH},
  If[th < 10^-12, Return[IdentityMatrix[3]]];
  rr = omega/th; HH = hat[rr];
  id3 + Sin[th] HH + (1 - Cos[th]) HH . HH
];

SmatrixFunc[omega_] := Module[{th = Norm[omega], rr, HH, RR},
  If[th < 10^-12, Return[IdentityMatrix[3]]];
  rr = omega/th; HH = hat[rr];
  RR = rodriguesExp[omega];
  IdentityMatrix[3] + HH . HH + (1/th)*(IdentityMatrix[3] - RR) . HH
];

SinvFunc[omega_] := Module[{th = Norm[omega], rr, HH, omx},
  If[th < 10^-12, Return[IdentityMatrix[3]]];
  rr = omega/th; HH = hat[rr];
  omx = (th/2)/Tan[th/2];
  omx id3 + (1 - omx) Outer[Times, rr, rr] - (th/2) HH
];

bigTnum = SmatrixFunc[omegaNum] . tExp // N;
Print["Test point:"];
Print["  omega = ", omegaNum];
Print["  T = S(omega).t = ", bigTnum];

(* Evaluate symbolic result *)
numrule = {theta -> thetaNum, r1 -> rNum[[1]], r2 -> rNum[[2]], r3 -> rNum[[3]],
           T1 -> bigTnum[[1]], T2 -> bigTnum[[2]], T3 -> bigTnum[[3]]};

dSinvTsymbolic = dSinvTdOmega /. numrule // N;
Print["\nd[S^{-1}T]/dOmega (symbolic):"];
Print[dSinvTsymbolic // MatrixForm];

(* FD ground truth *)
h = 10^-8;
dSinvTfd = Table[
  Module[{omP, omM},
    omP = omegaNum; omP[[j]] += h;
    omM = omegaNum; omM[[j]] -= h;
    (SinvFunc[omP] . bigTnum - SinvFunc[omM] . bigTnum)/(2 h)
  ],
  {j, 3}
] // Transpose;

Print["\nd[S^{-1}T]/dOmega (FD):"];
Print[dSinvTfd // MatrixForm];

errSymbolic = Max[Abs[Flatten[dSinvTsymbolic - dSinvTfd]]];
Print["\nSymbolic vs FD max error: ", ScientificForm[errSymbolic, 3]];


(* ===== STEP 4: Compare ALL candidate formulas numerically ===== *)

Print["\n================================================================"];
Print["STEP 4: Compare candidate closed-form formulas (numerical)"];
Print["================================================================"];

omxN = (thetaNum/2)/Tan[thetaNum/2] // N;
xN = 1 - omxN;
alphaN = (Sin[thetaNum] - thetaNum)/(2*(1 - Cos[thetaNum])) // N;
TbarN = rNum . bigTnum;
tbarN = rNum . tExp;
xOverTh = xN/thetaNum;
HNum = hat[rNum]; H2Num = HNum . HNum;
SNum = SmatrixFunc[omegaNum] // N;
SinvNegNum = SinvFunc[-omegaNum] // N;

(* --- Formula A: NEW (T-form, alpha-based) --- *)
formulaNew = 0.5 hat[bigTnum] + xOverTh Outer[Times, rNum, bigTnum] 
  + alphaN Outer[Times, bigTnum, rNum] 
  + TbarN*(xOverTh IdentityMatrix[3] - (2 xOverTh + alphaN) Outer[Times, rNum, rNum]);

errNew = Max[Abs[Flatten[formulaNew - dSinvTfd]]];
Print["Formula A — NEW (T-form, alpha): max err = ", ScientificForm[errNew, 3]];

(* --- Formula B: OLD structure with t (original paper Eq. 73 bottom) --- *)
formulaOldt = 0.5 hat[tExp] 
  + (xN/thetaNum)*(Outer[Times, rNum, tExp] - tbarN H2Num + Outer[Times, tExp, rNum])
  + (tbarN/thetaNum)*(SNum - SinvNegNum);

errOldt = Max[Abs[Flatten[formulaOldt - dSinvTfd]]];
Print["Formula B — OLD (t-form, S-S^{-1}): max err = ", ScientificForm[errOldt, 3]];

(* --- Formula C: OLD structure with T (naive t->T swap, WRONG) --- *)
formulaOldT = 0.5 hat[bigTnum] 
  + (xN/thetaNum)*(Outer[Times, rNum, bigTnum] - TbarN H2Num + Outer[Times, bigTnum, rNum])
  + (TbarN/thetaNum)*(SNum - SinvNegNum);

errOldT = Max[Abs[Flatten[formulaOldT - dSinvTfd]]];
Print["Formula C — OLD with T (naive swap): max err = ", ScientificForm[errOldT, 3]];


(* ===== STEP 5: Symbolic verification — both forms vs exact ===== *)

Print["\n================================================================"];
Print["STEP 5: Symbolic verification — both forms"];
Print["================================================================"];

alphaSymb = (Sin[theta] - theta)/(2*(1 - Cos[theta]));
xSymb = 1 - oneMinusXsymb;
xOverThSymb = xSymb/theta;
TbarSymb = Sum[r[[k]] bigT[[k]], {k, 3}];

(* Formula A (T-form) built symbolically *)
formulaNewSymb = (1/2) hat[bigT] + xOverThSymb Outer[Times, r, bigT]
  + alphaSymb Outer[Times, bigT, r]
  + TbarSymb*(xOverThSymb id3 - (2 xOverThSymb + alphaSymb) rrT);

(* Check: Formula A = exact symbolic derivative? *)
diffA = formulaNewSymb - dSinvTdOmega /. unitrule // FullSimplify;
Print["Formula A (T-form) - exact derivative:"];
diffAnum = diffA /. numrule // N;
Print["  Max numerical residual: ", Max[Abs[Flatten[diffAnum]]] // ScientificForm[#, 3]&];

isZeroA = Simplify[diffA, Assumptions -> {theta > 0, r1^2 + r2^2 + r3^2 == 1}];
Print["  Simplified: ", isZeroA // MatrixForm];

(* Formula B (t-form) — need to express t symbolically *)
(* t = S^{-1}(Omega) T, so we use t directly as a separate variable *)
Clear[t1, t2, t3];
tVec = {t1, t2, t3};
tbarSymb = Sum[r[[k]] tVec[[k]], {k, 3}];

(* S and S^{-1}(-Omega) symbolically *)
Rsymb = id3 + Sin[theta] H + (1 - Cos[theta]) H2;
Ssymb = id3 + H2 + (1/theta)*(id3 - Rsymb) . H // Simplify;
SinvNegSymb = oneMinusXsymb id3 + xSymb rrT + (theta/2) H;  (* S^{-1}(-Omega) = J_wr(Omega) *)

formulaOldSymb = (1/2) hat[tVec] 
  + (xSymb/theta)*(Outer[Times, r, tVec] - tbarSymb H2 + Outer[Times, tVec, r])
  + (tbarSymb/theta)*(Ssymb - SinvNegSymb);

(* To compare with exact, substitute t = S^{-1} T and T = S t *)
(* Actually: the exact derivative is d[S^{-1}T]/dOmega with T fixed.
   If the old formula with t is correct, then after substituting 
   T = S(Omega) t (where t is the exp coord), it should give the same result.
   
   The key insight: in the old formula, t is NOT held fixed — it is the 
   exponential coordinate which is the fundamental free variable.
   T = S(Omega) t, so when we vary Omega with t fixed, T varies too.
   
   But d[S^{-1}T]/dOmega with T FIXED is what we need for J_t.
   
   So the question is: does the old formula with t give d[S^{-1}(Omega) S(Omega) t]/dOmega 
   (which would be 0 since S^{-1}S = I), or does it give d[S^{-1}(Omega)]/dOmega . T?
   
   Let's just substitute T = S.t into Formula A and compare with Formula B. *)

(* Substitute T = S t into Formula A *)
Tsub = Ssymb . tVec // Simplify;
TbarSub = Simplify[Sum[r[[k]] Tsub[[k]], {k, 3}] /. unitrule];

formulaNewWithT = formulaNewSymb /. {T1 -> Tsub[[1]], T2 -> Tsub[[2]], T3 -> Tsub[[3]]};
formulaNewWithT = formulaNewWithT /. unitrule // Simplify;

(* Compare Formula A(T=St) with Formula B(t) *)
diffAB = formulaNewWithT - formulaOldSymb /. unitrule // FullSimplify;
Print["\n--- Key test: Formula A(T=St) - Formula B(t) ---"];
diffABnum = diffAB /. {theta -> thetaNum, r1 -> rNum[[1]], r2 -> rNum[[2]], r3 -> rNum[[3]],
                        t1 -> tExp[[1]], t2 -> tExp[[2]], t3 -> tExp[[3]]} // N;
Print["  Max numerical residual: ", Max[Abs[Flatten[diffABnum]]] // ScientificForm[#, 3]&];

diffABsimp = Simplify[diffAB, Assumptions -> {theta > 0, r1^2 + r2^2 + r3^2 == 1}];
Print["  Simplified: ", diffABsimp // MatrixForm];


(* ===== STEP 6: What went wrong in the original Rust implementation? ===== *)

Print["\n================================================================"];
Print["STEP 6: Diagnosing the original Rust bug"];
Print["================================================================"];

(* If both formulas are correct, the original Rust bug was NOT in the formula.
   Let's check: the original Rust code used t_exp (the exponential coordinate),
   but did it compute T̄ = r·t or T̄ = r·T? And which hat[] did it use? *)

Print["Original Rust j_coupling used:"];
Print["  Term 1: (1/2) hat(t_exp)  — this is [t]×, CORRECT for t-form"];
Print["  Term 2: (x/Θ)(r t^T - tbar H^2 + t r^T) — tbar = r.t, CORRECT for t-form"];
Print["  Term 3: (tbar/Θ)(S(Omega) - S^{-1}(-Omega)) — CORRECT for t-form"];
Print[""];
Print["If the t-form formula is symbolically correct, the bug must be elsewhere."];
Print["Checking numerically: the old formula with t at the test point..."];

(* The old formula with t, evaluated correctly *)
Print["\nOld formula (t-form) evaluated carefully:"];
Print[formulaOldt // MatrixForm];
Print["FD ground truth:"];
Print[dSinvTfd // MatrixForm];
Print["Error: ", ScientificForm[errOldt, 3]];

(* Element-by-element comparison *)
Print["\nElement-by-element:"];
Do[
  Print["  [", i, ",", j, "] old_t=", NumberForm[formulaOldt[[i, j]], 8],
    "  FD=", NumberForm[dSinvTfd[[i, j]], 8],
    "  err=", ScientificForm[Abs[formulaOldt[[i, j]] - dSinvTfd[[i, j]]], 2]],
  {i, 3}, {j, 3}
];


(* ===== STEP 7: Summary ===== *)

Print["\n================================================================"];
Print["STEP 7: Summary"];
Print["================================================================"];
Print[""];
Print["d[S^{-1}T]/dOmega has two equivalent forms (Paper Eq. 73):"];
Print[""];
Print["  T-form: (1/2)[T]× + (x/Θ)rT^T + αTr^T + T̄[(x/Θ)I - (2x/Θ+α)rr^T]"];
Print["  t-form: (1/2)[t]× + (x/Θ)(rt^T - t̄H² + tr^T) + (t̄/Θ)(S - S^{-1}(-Ω))"];
Print[""];
Print["  where T = S(Ω)t, T̄ = r·T, t̄ = r·t"];
Print["  α = (sinΘ - Θ)/(2(1-cosΘ))"];
Print[""];
Print["Numerical errors vs FD:"];
Print["  T-form (new):  ", ScientificForm[errNew, 3]];
Print["  t-form (old):  ", ScientificForm[errOldt, 3]];
Print["  T naively in old structure: ", ScientificForm[errOldT, 3], " (WRONG)"];
Print[""];
Print["Symbolic equivalence: Formula A(T=St) - Formula B(t) = ", 
  If[diffABsimp == Table[0, {3}, {3}], "ZERO (proven)", "check above"]];

