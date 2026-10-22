(* ::Package:: *)

(* ::Package:: *)
(**)


(* {{{ *)BeginPackage["Quantumthermo`",{"Carlos`", "Quantum`","Quantum2`"}];
Daviesgen::usage="Daviesgen[H_, \[Beta]_, \[Gamma]m_] gives the Davies generator"
MatrixCommutator::usage="Gives the conmutator"
timefordaviesdecomposition::usage="timefordaviesdecomposition[channel_?MatrixQ, Davies_?MatrixQ,tmin_?NumericQ, tmax_?NumericQ, step_?NumericQ, tol_?NumericQ]  Gives the time in an interval where the Channel can be decomposed"
DaviesMapPauliBasis[H_, \[Beta]_, t_, \[Gamma]_]::usage="DaviesMapPauliBasis[H_, \[Beta]_, t_, \[Gamma]_]  gives the Davies map at time t in Pauli Basis"
DaviesMapUnitBasis[H_, \[Beta]_, t_, \[Gamma]_]::usage="DaviesMapUnitBasis[H_, \[Beta]_, t_, \[Gamma]_]  gives the Davies map at time t in Unit Basis"
PureDensityMatrix[dim_Integer]::usage="PureDensityMatrix[dim_Integer] generates a density matrix that represents purestates of dim d"
BuildTimeFinder[Davies_?MatrixQ]::usage="BuildTimeFinder[\[IndentingNewLine]Davies_?MatrixQ] given a Davies generator builds a function to call the time for a set of data"
Ottocycle[Hh_,Hc_,\[Beta]h_,\[Beta]c_,th_,tc_,\[Tau]h_,\[Tau]c_,\[Gamma]h_,\[Gamma]c_,\[Tau]h_,\[Tau]c_]::usage="Ottocycle[Hh_,Hc_,\[Beta]h_,\[Beta]c_,th_,tc_,\[Tau]h_,\[Tau]c_,\[Gamma]h_,\[Gamma]c_,\[Tau]h_,\[Tau]c_] gives the Otto cycle represented in a Channel"
Ottoperfectquench[Hh_,Hc_,\[Beta]h_,\[Beta]c_,th_,tc_,\[Gamma]h_,\[Gamma]c_]::usage="Ottoperfectquench[Hh_,Hc_,\[Beta]h_,\[Beta]c_,th_,tc_,\[Gamma]h_,\[Gamma]c_] gives the Otto cycle considering a perfect quench"


Begin["`Private`"];

(* {{{ Implementaci\[OAcute]n privada *)

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
  MatrixExp[t QuantumMapInPauliBasis[Daviesgen[H, \[Beta], \[Gamma]]]] //
   Re;
DaviesMapUnitBasis[H_, \[Beta]_, t_, \[Gamma]_] :=
  MatrixExp[t QuantumMapInUnitBasis[Daviesgen[H, \[Beta], \[Gamma]]]] //
   Re;
PureDensityMatrix[dim_Integer]:=Module[{psi},psi=RandomState[dim];
Outer[Times,psi,Conjugate[psi]]];
pauliMatrices={IdentityMatrix[2],PauliMatrix[1],PauliMatrix[2],PauliMatrix[3]};
matrixToPauli[rho_]:=Table[Tr[rho . pauliMatrices[[i]]],{i,4}];
pauliToMatrix[c_]:=1/2 Sum[c[[i]] pauliMatrices[[i]],{i,4}];
Channelidxy[pid_,pz_]:=pid DiagonalMatrix[{1,1,1,1}]+pz DiagonalMatrix[{1,-1,-1,1}]+(1-pz-pid)DiagonalMatrix[{1,0,0,-1(*1/2*)}];
Channelidyz[pid_,px_]:=pid DiagonalMatrix[{1,1,1,1}]+px DiagonalMatrix[{1,1,-1,-1}]+(1-px-pid)DiagonalMatrix[{1,-1(*1/2*),0,0}];
nFailedDecompositions=0;
nMalformedInput=0;

BuildTimeFinder[
Davies_?MatrixQ]:=Module[{muD,vecsD,Pmat,PmatInv,negRates,maxDecayRate,tOverflowSafe},{muD,vecsD}=Eigensystem[N[Davies]];
Pmat=Transpose[vecsD];
PmatInv=Inverse[Pmat];
negRates=Select[Re[muD],#<0&];
maxDecayRate=If[Length[negRates]>0,Max[Abs[negRates]],1.];
tOverflowSafe=0.9 (Log[$MaxMachineNumber]/maxDecayRate);
Function[{channel,tmin,tmax,step,tol},Module[{choi,f,eps=10^-10,tCur,fPrev,fCur,tPrev,tL,tR,fL,found=False,currentTmax=Min[tmax,tOverflowSafe],maxAttempts=8,attempt=0,fAtTmin},If[!(MatrixQ[channel,NumericQ]&&Dimensions[channel]==={4,4}),nMalformedInput+=1;tmin,choi[t_]:=Module[{c,ev},ev=Exp[-t muD];
c=Chop[Reshuffle[FromPauliToUnit[Pmat . (ev*PmatInv) . channel]]];
(c+ConjugateTranspose[c])/2];
f[t_?NumericQ]:=Module[{val},val=Quiet[Min[Re@Eigenvalues[choi[t]]]];
If[NumericQ[val]&&FreeQ[val,Overflow|Underflow|Indeterminate|ComplexInfinity|DirectedInfinity],val,$Failed]];
fAtTmin=f[tmin];
If[Chop[fAtTmin]<=eps,tmin,tPrev=tmin;fPrev=fAtTmin;
For[tCur=tmin+step,tCur<=currentTmax&&!found,tCur+=step,fCur=f[tCur];
If[fCur===$Failed,Break[]];
If[Chop[fCur]==0,tL=tR=tCur;found=True,If[fPrev>eps&&fCur<-eps,tL=tPrev;tR=tCur;found=True]];
tPrev=tCur;fPrev=fCur];
While[!found&&attempt<maxAttempts&&currentTmax<tOverflowSafe,currentTmax=Min[currentTmax*2,tOverflowSafe];
tPrev=tmin;fPrev=fAtTmin;
For[tCur=tmin+step,tCur<=currentTmax&&!found,tCur+=step,fCur=f[tCur];
If[fCur===$Failed,Break[]];
If[Chop[fCur]==0,tL=tR=tCur;found=True,If[fPrev>eps&&fCur<-eps,tL=tPrev;tR=tCur;found=True]];
tPrev=tCur;fPrev=fCur];
attempt++];
If[!found,nFailedDecompositions+=1;tmin,If[tL==tR,tL,fL=f[tL];
While[tR-tL>tol,Module[{tm=(tL+tR)/2,fm},fm=f[tm];
Which[fm===$Failed,tR=tm,Chop[fm]==0,tL=tR=tm,fL*fm<0,tR=tm,True,tL=tm;fL=fm]]];
(tL+tR)/2]]]]]]];
Ottocycle[Hh_,Hc_,\[Beta]h_,\[Beta]c_,th_,tc_,\[Tau]h_,\[Tau]c_,\[Gamma]h_,\[Gamma]c_,\[Tau]h_,\[Tau]c_]:=Module[{pauliOps,S,Ut1,Ut2,Ut1sup,Ut2sup,pauliunit1,pauliunitt2,lcpauli,lhpauli,paso1,paso2,paso3,paso4,paso5,paso6,paso7,otto,ideal,U\[Tau]1,U\[Tau]2,U\[Tau]1sup,U\[Tau]2sup,pauliunit\[Tau]1,pauliunit\[Tau]2},
pauliOps={IdentityMatrix[2],PauliMatrix[1],PauliMatrix[2], PauliMatrix[3]}/Sqrt[2];
S=Transpose[Table[Flatten[Transpose[pauliOps[[k]]]],{k,1,4}]];
Ut1=MatrixExp[-I Hc tc]//FullSimplify;
     Ut2=MatrixExp[-I Hh th//FullSimplify];
U\[Tau]1=MatrixExp[-I Hc \[Tau]c]//FullSimplify;   
     U\[Tau]2=MatrixExp[-I Hc \[Tau]h]//FullSimplify;   
U\[Tau]1sup=KroneckerProduct[Conjugate[U\[Tau]1],U\[Tau]1]//FullSimplify; 
     U\[Tau]2sup=KroneckerProduct[Conjugate[U\[Tau]2],U\[Tau]2]//FullSimplify; 
Ut1sup=KroneckerProduct[Conjugate[Ut1],Ut1]//FullSimplify;
     Ut2sup=KroneckerProduct[Conjugate[Ut2],Ut2]//FullSimplify;
pauliunit\[Tau]1= ConjugateTranspose[S] . U\[Tau]1sup . S // FullSimplify;
pauliunit\[Tau]2= ConjugateTranspose[S] . U\[Tau]2sup . S // FullSimplify;
pauliunit1= ConjugateTranspose[S] . Ut1sup . S // FullSimplify;
pauliunitt2= ConjugateTranspose[S] . Ut2sup . S // FullSimplify;
lcpauli=DaviesMapPauliBasis[Hc,\[Beta]c,tc,\[Gamma]c]//FullSimplify;
lhpauli=DaviesMapPauliBasis[Hh,\[Beta]h,th,\[Gamma]h]//FullSimplify;
paso1=pauliunit1;
paso2=paso1 . IdentityMatrix[4]//FullSimplify;
paso3=paso2 . pauliunitt2//FullSimplify;
paso4=paso3 . lhpauli//FullSimplify;
paso5=paso4 . pauliunit\[Tau]2//FullSimplify;
paso6=paso5 . IdentityMatrix[4]//FullSimplify;
paso7=paso6 . pauliunit\[Tau]1//FullSimplify;
otto=paso7 . lcpauli//FullSimplify]

Ottoperfectquench[Hh_,Hc_,\[Beta]h_,\[Beta]c_,th_,tc_,\[Gamma]h_,\[Gamma]c_]:=Module[{pauliOps,S,Ut1,Ut2,Ut1sup,Ut2sup,pauliunit1,pauliunitt2,lcpauli,lhpauli,paso1,paso2,paso3,paso4,ideal},
pauliOps={IdentityMatrix[2],PauliMatrix[1],PauliMatrix[2], PauliMatrix[3]}/Sqrt[2];
S=Transpose[Table[Flatten[Transpose[pauliOps[[k]]]],{k,1,4}]];
(*Ut1=MatrixExp[-I Hc tc]//FullSimplify;
     Ut2=MatrixExp[-I Hh th//FullSimplify];
Ut1sup=KroneckerProduct[Conjugate[Ut1],Ut1]//FullSimplify;
     Ut2sup=KroneckerProduct[Conjugate[Ut2],Ut2]//FullSimplify;*);
pauliunit1= ConjugateTranspose[S] . Ut1sup . S // FullSimplify;
pauliunitt2= ConjugateTranspose[S] . Ut2sup . S // FullSimplify;
lcpauli=DaviesMapPauliBasis[Hc,\[Beta]c,tc,\[Gamma]c]//FullSimplify;
lhpauli=DaviesMapPauliBasis[Hh,\[Beta]h,th,\[Gamma]h]//FullSimplify;
paso1=IdentityMatrix[4];
paso2=paso1 . lhpauli;
paso3=paso2 . IdentityMatrix[4]//FullSimplify;
paso4=paso3 . lcpauli//FullSimplify;
Chop[ideal=paso4//FullSimplify,10^-10]]
(* }}} *)
(*random1*)

End[];

EndPackage[];
