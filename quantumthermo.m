
(* ::Package:: *)


(* {{{ *)BeginPackage["Quantumthermo`",{"Carlos`", "Quantum`","Quantum2`"}];
Daviesgen::usage="Daviesgen[H_, \[Beta]_, \[Gamma]m_] gives the Davies generator"
MatrixCommutator::usage="Gives the conmutator"
timefordaviesdecomposition::usage="timefordaviesdecomposition[channel_?MatrixQ, Davies_?MatrixQ,tmin_?NumericQ, tmax_?NumericQ, step_?NumericQ, tol_?NumericQ]  Gives the time in an interval where the Channel can be decomposed"
DaviesMapPauliBasis[H_, \[Beta]_, t_, \[Gamma]_]::usage="DaviesMapPauliBasis[H_, \[Beta]_, t_, \[Gamma]_]  gives the Davies map at time t in Pauli Basis"
DaviesMapUnitBasis[H_, \[Beta]_, t_, \[Gamma]_]::usage="DaviesMapUnitBasis[H_, \[Beta]_, t_, \[Gamma]_]  gives the Davies map at time t in Unit Basis"



Begin["`Private`"];

(* {{{ Implementación privada *)

Daviesgen[H_, \[Beta]_, \[Gamma]m_][\[Rho]_] :=
  Module[{\[Sigma]p, \[Sigma]m, vals, vecs, \[Omega]},
   {vals, vecs} = H // Eigensystem;
   \[Omega] = vals[[2]] - vals[[1]];
   vecs = Normalize /@ Orthogonalize[vecs];
   \[Sigma]p = Proyector[vecs[[2]], vecs[[1]]];
   \[Sigma]m = Proyector[vecs[[1]], vecs[[2]]];
   -I MatrixCommutator[
      H, \[Rho]] + \[Gamma]m Exp[-\[Beta] \[Omega]] (\[Sigma]p . \
\[Rho] . \[Sigma]m -
       1/2 (\[Sigma]m . \[Sigma]p . \[Rho] + \[Rho] . \[Sigma]m . \
\[Sigma]p)) + \[Gamma]m (\[Sigma]m . \[Rho] . \[Sigma]p -
       1/2 (\[Sigma]p . \[Sigma]m . \[Rho] + \[Rho] . \[Sigma]p . \
\[Sigma]m))];

Daviesgen[H_, \[Beta]_, \[Gamma]m_, \[Gamma]deph_][\[Rho]_] :=
  Module[{\[Sigma]p, \[Sigma]m, vals, vecs, \[Omega]},
   {vals, vecs} = H // Eigensystem;
   vecs = Normalize /@ Orthogonalize[vecs];
   \[Omega] = vals[[2]] - vals[[1]];
   \[Sigma]p = Proyector[vecs[[2]], vecs[[1]]];
   \[Sigma]m = Proyector[vecs[[1]], vecs[[2]]];
   -I MatrixCommutator[
      H, \[Rho]] + \[Gamma]m Exp[-\[Beta] \[Omega]] (\[Sigma]p . \
\[Rho] . \[Sigma]m -
       1/2 (\[Sigma]m . \[Sigma]p . \[Rho] + \[Rho] . \[Sigma]m . \
\[Sigma]p)) + \[Gamma]m (\[Sigma]m . \[Rho] . \[Sigma]p -
       1/2 (\[Sigma]p . \[Sigma]m . \[Rho] + \[Rho] . \[Sigma]p . \
\[Sigma]m)) + \[Gamma]deph (Proyector[vecs[[1]]] . \[Rho] .
        Proyector[vecs[[1]]] +
       Proyector[vecs[[2]]] . \[Rho] . Proyector[vecs[[2]]] - \[Rho])];

MatrixCommutator[A_?MatrixQ, B_?MatrixQ] := A . B - B . A;

  timefordaviesdecomposition[channel_?MatrixQ, Davies_?MatrixQ,
      tmin_?NumericQ, tmax_?NumericQ, step_?NumericQ, tol_?NumericQ] :=
    Module[{choi, f, data, signChangeIntervals, tL, tR, tm,
        eps = 10^-10, maxAttempts = 10, attempt = 0, tMax = tmax},
      choi[t_] := Module[{c},
          c = (MatrixExp[-t Davies] . channel) // FromPauliToUnit //
                 Reshuffle // Chop;
          (c + ConjugateTranspose[c])/2];
      f[t_?NumericQ] := Min[Re@Eigenvalues[choi[t]]];
      If[f[tmin] <= eps, Return[tmin]];
      signChangeIntervals = {};
      While[attempt < maxAttempts,
        data = Table[{t, f[t]}, {t, tmin, tMax, step}];
        signChangeIntervals =

     Select[Partition[data, 2,
       1], #[[1, 2]] > eps && #[[2, 2]] < -eps &];
        If[Length[signChangeIntervals] > 0,
          {tL, tR} = signChangeIntervals[[1]][[All, 1]];
          Break[]];
        tMax *= 2; attempt++];
      If[Length[signChangeIntervals] == 0,
        nFailedDecompositions += 1;
        Return[Null]];
      While[tR - tL > tol, tm = (tL + tR)/2;
        If[f[tm] > eps, tL = tm, tR = tm]];
      (tL + tR)/2];

DaviesMapPauliBasis[H_, \[Beta]_, t_, \[Gamma]_] :=
  MatrixExp[t QuantumMapInPauliBasis[Lth[H, \[Beta], \[Gamma]]]] //
   Re;
DaviesMapUnitBasis[H_, \[Beta]_, t_, \[Gamma]_] :=
  MatrixExp[t QuantumMapInUnitBasis[Lth[H, \[Beta], \[Gamma]]]] //
   Re;

(* }}} *)

End[];

EndPackage[];
