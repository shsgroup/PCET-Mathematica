(* ::Package:: *)

(* Mathematica Package *)

(*
MathPCET
*)

(*
Nonadiabatic PCET: rate expressions etc.
*)

BeginPackage["MathPCET`"]

Needs["NumericalCalculus`"]

MorsePotential::usage="MorsePotential[DE,Beta,x] ->
Returns Morse potential in kcal/mol at x (Bohr) with the minimum at x=0.
Parameters:
DE is the dissociation energy (kcal/mol),
Beta is the parameter (1/Bohr) related to the force constant at the minimum.";

HarmonicPotential::usage="HarmonicPotential[Mass,Frequency,x]
Returns Harmonic potential in kcal/mol at x (Bohr) with the minimum at x=0.
Parameters:
Mass is a mass of the particle (Daltons),
Frequency is the frequency (1/cm).";

MorsePotentialInverted::usage="MorsePotentialInverted[DE,Beta,x]
Returns inverted Morse potential in kcal/mol at x (Bohr) with the minimum at x=0.
Parameters:
DE is the dissociation energy (kcal/mol),
Beta is the parameter (1/Bohr) related to the force constant at the minimum.";

MorseBeta::usage="MorseBeta[Omega,DE]
Returns the Morse beta parameter (1/Bohr) corresponding to potential \
for the vibration with frequency Omega (1/cm).
Parameters:
DE is the dissociation energy in kcal/mol,
Omega is the frequency (1/cm).
Options:
M -> mass of the particle in Daltons.";

MorseEnergy::usage="MorseEnergy[n,Beta,DE]
Returns the energy (kcal/mol) of the Morse oscillator relative to the bottom of the potential well.
Parameters:
n = 0, 1, 2,... - quantum number,
Beta - Morse parameter (in 1/Bohr),
DE - dissociation energy (in kcal/mol).
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).";

MorseBound::usage="MorseBound[Beta,DE]
Returns number of bound states for the Morse oscillator.
Parameters:
Beta - beta parameter (1/Bohr),
DE - dissociation energy (kcal/mol),
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).";

MorseWavefunctionLeft::usage="MorseWavefunctionLeft[n,Beta,DE,x]
Returns wavefunction for the Morse oscillator in units of (Bohr)^(-1/2).
Parameters:
n = 0, 1, 2,... - quantum number,
Beta - Morse parameter (in 1/Bohr),
DE - dissociation energy (in kcal/mol),
x - coordinate (in Bohrs) relative to the minimum of the potential.
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).";

MorseWavefunctionRight::usage="MorseWavefunctionRight[n,Beta,DE,x]
Returns wavefunction for the mirrored Morse oscillator in units of (Bohr)^(-1/2).
Parameters:
n = 0, 1, 2,... - quantum number,
Beta - Morse parameter (in 1/Bohr),
DE - dissociation energy (in kcal/mol),
x - coordinate (in Bohrs) relative to the minimum of the potential.
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).";

HOWavefunction::usage="HOWavefunction[n_,Freq_,x_]
Returns wavefunction for the Harmonic oscillator in units of (Bohr)^(-1/2).
Parameters:
n - quantum number (0 for ground state),
Freq - frequency (1/cm),
x - coordinate (Bohrs),
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).";

HOEnergy::usage="HOEnergy[freq_,n_]
Returns the energy (kcal/mol) level of a quantum harmonic oscillator.
Parameters:
freq - frequency (1/cm),
n - quantum number (0 for ground state).";

QHarmonic::usage="Coming Soon!";

QMorseExact::usage="Coming Soon!";

QMorseApprox::usage="Coming Soon!";

MarcusReactantMorse::usage="MarcusReactant[n_,beta_,de_,lambda_,X_]
Marcus parabola for a given reactant vibronic state.
Returns free energy in kcal/mol.
Parameters:
n = 0, 1, 2,... - vibrational quantum number,
beta - Morse parameter (in 1/Bohr),
de - dissociation energy (in kcal/mol),
lambda - reorganization energy (kcal/mol),
X - energy gap reaction coordinate (in kcal/mol).
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).";

MarcusProductMorse::usage="MarcusProduct[n_,beta_,de_,lambda_,X_]
Marcus parabola for a given product vibronic state.
Returns free energy in kcal/mol.
Parameters:
n = 0, 1, 2,... - vibrational quantum number,
beta - Morse parameter (in 1/Bohr),
de - dissociation energy (in kcal/mol),
lambda - reorganization energy (kcal/mol),
X - energy gap reaction coordinate (in kcal/mol).
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).";

HarmonicOverlap::usage="HarmonicOverlap[n1_,n2_,f1_,f2_,d_]
Overlap between two harmonic oscillator wavefunctions
[J.-L. Chang, J. Mol. Spectrosc. 232 (2005) 102-104]
Parameters:
n1 - quantum number for the oscillator on the left (0 for ground state),
n2 - quantum number for the oscillator on the right (0 for ground state),
f1 - frequency for the oscillator on the left (1/cm),
f2 - frequency for the oscillator on the right (1/cm),
d - distance (Bohr) between the minima of the displaced harmonic potentials.
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).";

MorseOverlapSym::usage="MorseOverlap0[DE_,Beta_,nu1_,nu2_,d_]
Morse overlap for identical mirrored Morse potentials.
[modified expression from J. C. Lopez, A. L. Rivera, Yu. F. Smirnov, A. Frank, arXiv:physics/0109017v].
Parameters:
DE - dissociation energy (kcal/mol),
Beta - beta parameter (1/Bohr),
nu1 - quantum number for the oscillator on the left (0 for ground state),
nu2 - quantum number for the oscillator on the right (0 for ground state),
d - distance between the minima of the displaced and mirrored Morse potentials.
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).
Acc -> accuracy (number of digits after decimal point), Default is 50.";

MorseOverlap::usage="MorseOverlap[DE1_,Beta1_,DE2_,Beta2_,nu1_,nu2_,d_]
Morse overlap for general mirrored Morse potentials
[modified expression from J. C. Lopez, A. L. Rivera, Yu. F. Smirnov, A. Frank, arXiv:physics/0109017v]
Parameters:
DE1 - dissociation energy (kcal/mol),
DE2 - dissociation energy (kcal/mol),
Beta1 - beta parameter (1/Bohr),
Beta2 - beta parameter (1/Bohr),
nu1 - quantum number for the oscillator on the left (0 for ground state),
nu2 - quantum number for the oscillator on the right (0 for ground state),
d - distance between the minima of the displaced and mirrored Morse potentials.
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).
Acc -> accuracy (number of digits after decimal point), Default is 50.";

NMorseOverlap::usage="NMorseOverlap[DE1_,Beta1_,DE2_,Beta2_,nu1_,nu2_,d_]
Numerical Morse overlap for general mirrored Morse potentials.
Parameters:
DE1 - dissociation energy (kcal/mol),
DE2 - dissociation energy (kcal/mol),
Beta1 - beta parameter (1/Bohr),
Beta2 - beta parameter (1/Bohr),
nu1 - quantum number for the oscillator on the left (0 for ground state),
nu2 - quantum number for the oscillator on the right (0 for ground state),
d - distance between the minima of the displaced and mirrored Morse potentials.
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).";

NAlphaMorse::usage="NAlphaMorse[DE1_,Beta1_,DE2_,Beta2_,i_,j_,d_]
Distance dependence parameter: Morse overlap for general mirrored Morse potentials.
Numerical logarithmic derivative of the numerical overlap integrals.
Parameters:
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
i, j - quantum numbers for the left and right oscillators, respectively,
d - distance between the minima of the morse potentials (Bohr).
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).
Step -> step size in numerical differentiation, Default is 0.0001 Bohr.";

AlphaMorseSym::usage="AlphaMorse0[DE_,Beta_,nu1_,nu2_,d_]
Distance dependence parameter: Morse overlap for identical mirrored Morse potentials.
Analytical derivative of analytical overlap integrals (for identical mirrored potentials).
Parameters:
DE - dissociation energy for the oscillator (kcal/mol),
Beta - beta parameter for the oscillator on the left (1/Bohr),
nu1, nu2 - quantum numbers for the left and right oscillators, respectively,
d - distance between the minima of the morse potentials (Bohr).
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).
Acc -> accuracy (number of digits after decimal point), Default is 50.";

AlphaMorse::usage="AlphaMorse0[DE1_,Beta1_,DE2_,Beta2_,nu1_,nu2_,d_]
Distance dependence parameter: Morse overlap for identical mirrored Morse potentials.
Analytical derivative of analytical overlap integrals (for identical mirrored potentials).
Parameters:
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
nu1, nu2 - quantum numbers for the left and right oscillators, respectively,
d - distance between the minima of the morse potentials (Bohr).
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).
Acc -> accuracy (number of digits after decimal point), Default is 50.";

PClassical::usage="PClassical[T_, R_, MR_, FR_, x_]:
Classical distribution function for harmonic oscillator.
Parameters:
T - temperature (K),
R - equilibrium value of coordinate (Bohr),
MR - mass of the oscillator (Dalton),
FR - frequency of the oscillator (1/cm),
x - coordinate (Bohr).";

PQuantal::usage="PQuantal[T_, R_, MR_, FR_, x_]:
Quantal distribution function for harmonic oscillator.
Parameters:
T - temperature (K),
R - equilibrium value of coordinate (Bohr),
MR - mass of the oscillator (Dalton),
FR - frequency of the oscillator (1/cm),
x - coordinate (Bohr).";

HarmonicFranckCondonAveraged::usage="HarmonicFranckCondonAveraged[mu_,nu_,f1_,f2_,d_,T_,MR_,FR_]:
Averaged Franck-Condon factor for two harmonic wavefunctions.
Parameters:
mu, nu - quantum numbers for left and right wavefunctions, respectively,
f1, f2 - frequancies of left and right oscillators, respectively (1/cm),
d - equilibrium distance between the minima of left and right oscillators (Bohr),
T - temperature (K),
MR - mass of the R-oscillator (Dalton),
FR - frequency of the R-oscillator (1/cm).
Options:
M->mass of the particle in Daltons (Default - mass of the proton),
Distribution->
Classical or Quantal (for R-oscillator distribution function).";

MorseFranckCondonAveraged::usage="HarmonicFranckCondonAveraged[mu_,nu_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,MR_,FR_]:
Averaged Franck-Condon factor for two Morse wavefunctions.
Parameters:
mu, nu - quantum numbers for left and right wavefunctions, respectively,
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of left and right oscillators (Bohr),
T - temperature (K),
MR - mass of the R-oscillator (Dalton),
FR - frequency of the R-oscillator (1/cm).
Options:
M->mass of the particle in Daltons (Default - mass of the proton),
Overlap->
Numerical or Analytical (Default - Numerical),
Distribution->
Classical or Quantal (for R-oscillator distribution function).";

Sigma2RClassical::usage="Sigma2RClassical[MR_,Freq_,T_]
Mean square displacement of the classical oscillator.
Parameters:
MR - reduced mass of the oscillator (Daltons),
Freq - frequency of the oscillator (1/cm),
T - temperature (K).
Returns variance in Bohr^2.";

Sigma2RQuantum::usage="Sigma2RQuantum[MR_,Freq_,T_]
Mean square displacement of the quantum oscillator.
Parameters:
MR - reduced mass of the oscillator (Daltons),
Freq - frequency of the oscillator (1/cm),
T - temperature (K).
Returns variance in Bohr^2.";

CRClassical::usage="CRClassical[MR_,Freq_,T_,t_]
Time correlation function of the classical oscillator.
Parameters:
MR - reduced mass of the oscillator (Daltons),
Freq - frequency of the oscillator (1/cm),
T - temperature (K),
t - time (atomic units).
Returns correlation function in Bohr^2.";

CRQuantum::usage="CRQuantum[MR_,Freq_,T_,t_]
Time correlation function of the quantum oscillator.
Parameters:
MR - reduced mass of the oscillator (Daltons),
Freq - frequency of the oscillator (1/cm),
T - temperature (K),
t - time (atomic units).
Returns correlation function in Bohr^2.";

NLambdaAlphaMorse::usage="NLambdaAlphaMorse[MR_,DE1_,Beta1_,DE2_,Beta2_,i_,j_,d_]
Coupling reorganization energy with numerical alpha.
Expression for general mirrored Morse potentials.
Parameters:
MR - reduced mass of the donor-acceptor mode (Dalton),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
i, j - quantum numbers for the left and right oscillators, respectively,
d - distance between the minima of the morse potentials (Bohr).
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).
Step -> step size in finite differences method, Default is 0.0001 Bohr.
Returns the coupling reorganization energy in kcal/mol.";

LambdaAlphaMorseSym::usage="LambdaAlphaMorseSym[MR_,DE_,Beta_,nu1_,nu2_,d_]
Coupling reorganization energy with analytical alpha (for identical mirrored Morse potentials).
Expression using analytical derivative of analytical overlap integrals (for identical mirrored potentials).
Parameters:
DE - dissociation energy for the oscillator (kcal/mol),
Beta - beta parameter for the oscillator on the left (1/Bohr),
nu1, nu2 - quantum numbers for the left and right oscillators, respectively,
d - distance between the minima of the morse potentials (Bohr).
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).
Acc -> accuracy (number of digits after decimal point), Default is 50.";

LambdaAlphaMorse::usage="LambdaAlphaMorse[MR_,DE_,Beta_,nu1_,nu2_,d_]
Coupling reorganization energy with analytical alpha (for general mirrored Morse potentials).
Expression using analytical derivative of analytical overlap integrals.
Parameters:
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
nu1, nu2 - quantum numbers for the left and right oscillators, respectively,
d - distance between the minima of the morse potentials (Bohr).
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).
Acc -> accuracy (number of digits after decimal point), Default is 50.";

LambdaAlphaInput::usage="LambdaAlphaInput[MR_,alpha_]
Coupling reorganization energy for constant alpha given as input parameter.
Parameters:
MR - reduced mass of the donor-acceptor mode (Daltons),
alpha - alpha parameter (1/Bohr).
Returns coupling reorganization energy in kcal/mol.";

SolventDampingClassical::usage="SolventDampingClassical[lambda_,T_,t_]
Classical solvent damping term (Gaussian damping for high temperature and short time).
Parameters:
lambda - solvent reorganization energy (kcal/mol),
T - temperature (K),
t - time (atomic units).
Return solvent damping Gaussian (unitless).";

CoherentTermMorse::usage="CoherentTermMorse[DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,nu1_,nu2_,t_]
Coherent (oscillatory) term.
Parameters:
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
nu1, nu2 - quantum numbers for the left and right oscillators, respectively,
d - distance between the minima of the morse potentials (Bohr),
t - time (atomic units).
Options:
M -> mass of the particle in Daltons (Default - mass of the proton),
Alpha -> method for alpha calculation:
Numerical (default), Analytical, AnalyticalSym, or <value> for fixed input.";

CoherentTermFGH::usage="CoherentTermFGH[fghList_,DG_,lambda_,MR_,Freq_,nu1_,nu2_,t_]
Coherent (oscillatory) term.
Parameters:
fghList - nested list with the results of FGH calculation,
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
nu1, nu2 - quantum numbers for the left and right oscillators, respectively,
t - time (atomic units).
Options:
Alpha -> method for alpha calculation:
Numerical (default) <value> for fixed input.";

CoherentTermFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

RmodeTermMorse::usage="Term originating from the donor-acceptor mode vibrational mode.
RmodeTerm[MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,nu1_,nu2_,t_].
Parameters:
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - distance between the minima of the morse potentials (Bohr),
T - temperature (K),
nu1, nu2 - quantum numbers for the left and right oscillators, respectively,
t - time (atomic units).
Options:
M -> mass of the particle in Daltons (Default - mass of the proton),
Alpha -> method for alpha calculation:
Numerical (default), Analytical, AnalyticalSym, or <value> for fixed input.";

RmodeTermFGH::usage="Term originating from the donor-acceptor mode vibrational mode.
RmodeTermFGH[fghList_,MR_,Freq_,T_,nu1_,nu2_,t_].
Parameters:
fghList - nested list with the results of FGH calculation,
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
T - temperature (K),
nu1, nu2 - quantum numbers for the left and right oscillators, respectively,
t - time (atomic units).
Options:
Alpha -> method for alpha calculation:
Numerical (default) or <value> for fixed input.";

RmodeTermFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

RealFluxQuantumMorse::usage="Total nonadiabatic quantum flux correlation function.
RealFluxQuantum[DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,nu1_,nu2_,t_]
Parameters:
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of the morse potentials (Bohr),
T - temperature (K),
nu1, nu2 - quantum numbers for the left and right oscillators, respectively,
t - time (atomic units).
Options:
M -> mass of the particle in Daltons (Default - mass of the proton),
Alpha -> method for alpha calculation:
Numerical (default), Analytical, AnalyticalSym, or <value> for fixed input.";

RealFluxQuantumFGH::usage="Total nonadiabatic quantum flux correlation function.
RealFluxQuantumFGH[fghList_,DG_,lambda_,MR_,Freq_,T_,nu1_,nu2_,t_]
Parameters:
fghList - nested list with the results of FGH calculation,
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
T - temperature (K),
nu1, nu2 - quantum numbers for the left and right oscillators, respectively,
t - time (atomic units).
Options:
Alpha -> method for alpha calculation:
Numerical (default) or <value> for fixed input.";

PartialRateQfluxMorse::usage="Partial quantum nonadiabatic first order rate constant k(mu->nu) (1/sec).
PartialRateQflux[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mu_,nu_]
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of the morse potentials (Bohr),
T - temperature (K),
mu, nu - quantum numbers for the left and right oscillators, respectively.
Options:
M -> mass of the particle in Daltons (Default - mass of the proton),
Overlap -> method for overlap calculation:
Numerical (default), Analytical, or AnalyticalSym,
Alpha -> method for alpha calculation:
Numerical (default), Analytical, AnalyticalSym, or <value> for fixed input.";

PartialRateQfluxFGH::usage="Partial quantum nonadiabatic first order rate constant k(mu->nu) (1/sec).
PartialRateQfluxFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mu_,nu_]
Parameters:
fghList - nested list with the results of FGH calculation,
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
T - temperature (K),
mu, nu - quantum numbers for the left and right oscillators, respectively.
Options:
Alpha -> method for alpha calculation:
Numerical (default) or <value> for fixed input.";

PartialRateQfluxFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

TotalRateQfluxMorse::usage="Total quantum nonadiabatic first order rate constant (1/sec).
TotalRateQflux[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_]
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of the morse potentials (Bohr),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states.
Options:
M -> mass of the particle in Daltons (Default - mass of the proton),
Overlap -> method for overlap calculation:
Numerical (default), Analytical, or AnalyticalSym,
Alpha -> method for alpha calculation:
Numerical (default), Analytical, AnalyticalSym, or <value> for fixed input.";

TotalRateQfluxFGH::usage="Total quantum nonadiabatic first order rate constant (1/sec).
TotalRateQfluxFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_]
Parameters:
fghList - nested list with the results of FGH calculation,
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states.
Options:
Alpha -> method for alpha calculation:
Numerical (default) or <value> for fixed input.";

TotalRateQfluxFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

KIEQfluxMorse::usage="Kinetic isotope effect (KIE) for quantum nonadiabatic rate constant.
KIEQflux[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_]
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of the morse potentials (Bohr),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states.
Options:
M1 -> mass of the particle for the rate in the numerator (default: mass of the Protium),
M2 -> mass of the particle for the rate in the denominator (default: mass of the Deuterium),
Overlap -> method for alpha calculation:
Numerical (default), Analytical, or AnalyticalSym,
Alpha -> method for alpha calculation:
Numerical (default), Analytical, AnalyticalSym, or <value> for fixed input.";

KIEQfluxFGH::usage="Kinetic isotope effect (KIE) for quantum nonadiabatic rate constant.
KIEQfluxFGH[fghList1_,fghList2_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_]
Parameters:
fghList1 - nested list with the results of FGH calculation for isotope 1,
fghList2 - nested list with the results of FGH calculation for isotope 2,
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states.
Options:
Alpha1 -> method for alpha calculation for the isotope 1:
Numerical (default) <value> for fixed input,
Alpha2 -> method for alpha calculation for the isotope 2:
Numerical (default) <value> for fixed input.";

TableRateQfluxMorse::usage="Table with rate channel contributions for quantum nonadiabatic rate constant.
TableRateQflux[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_]
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor - acceptor mode (Dalton),
Freq - donor - acceptor mode frequency (1/cm),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of the morse potentials (Bohr),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states.
Options:
M -> mass of the particle in Daltons (Default - mass of the proton),
Overlap -> method for alpha calculation:
Numerical (default), Analytical, or AnalyticalSym,
Alpha -> method for alpha calculation:
Numerical (default), Analytical, AnalyticalSym, or <value> for fixed input.";

TableRateQfluxFGH::usage="Table with rate channel contributions for quantum nonadiabatic rate constant.
TableRateQfluxFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_]
Parameters:
fghList - nested list with the results of FGH calculation,
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor - acceptor mode (Dalton),
Freq - donor - acceptor mode frequency (1/cm),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states.
Options:
Alpha -> method for alpha calculation:
Numerical (default) or <value> for fixed input.";

TableRateQfluxFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

PartialRateRfixedMorse::usage="Partial nonadiabatic first order rate constant k(mu->nu) (1/sec) for a fixed proton donor-acceptor distance (Morse potentials).
PartialRateRfixed[V_,DG_,lambda_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mu_,nu_]
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of the morse potentials (Bohr),
T - temperature (K),
mu, nu - quantum numbers for the left and right oscillators, respectivly.
Options:
M -> mass of the particle in Daltons (Default - mass of the proton).
Overlap -> method for overlap calculation:
Numerical (default), Analytical, or AnalyticalSym.
Return partial rate constant in sec^-1.";

PartialRateRfixedFGH::usage="Partial nonadiabatic first order rate constant k(mu->nu) (1/sec) for a fixed proton donor-acceptor distance (general potentials).
PartialRateRfixed[fghList_,V_,DG_,lambda_,T_,mu_,nu_]
Parameters:
fghList - nested list with the results of FGH calculation,
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
T - temperature (K),
mu, nu - quantum numbers for the left and right oscillators, respectivly.
Options: None.
Return partial rate constant in sec^-1.";

PartialRateRfixedFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

TotalRateRfixedMorse::usage="Total nonadiabatic first order rate constant (1/sec) for a fixed proton donor-acceptor distance.
TotalRateRfixed[V_,DG_,lambda_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_].
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of the morse potentials (Bohr),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states. 
Options: 
M -> mass of the particle in Daltons (Default - mass of the proton). 
Overlap -> method for overlap calculation, 
Numerical (default), Analytical, or AnalyticalSym.";

TotalRateRfixedFGH::usage="Total nonadiabatic first order rate constant (1/sec) for a fixed proton donor-acceptor distance (general potentials).
TotalRateRfixed[fghList_,V_,DG_,lambda_,T_,mumax_,numax_].
Parameters:
fghList - nested list with the results of FGH calculation,
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states. 
Options: None.";

TotalRateRfixedFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

KIERfixedMorse::usage="Kinetic isotope effect for nonadiabatic rate constant for a fixed proton donor-acceptor distance (Morse potentials).
KIERfixed[V_,DG_,lambda_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_].
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of the morse potentials (Bohr),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states. 
Options: 
M1 -> mass of the particle for the rate in the numerator (default: mass of the Protium),
M2 -> mass of the particle for the rate in the denominator (default: mass of the Deuterium),
Overlap -> method for overlap calculation: 
Numerical (default), Analytical, or AnalyticalSym.";

KIERfixedFGH::usage="Kinetic isotope effect for nonadiabatic rate constant for a fixed proton donor-acceptor distance (general potentials).
KIERfixedFGH[fghList1_,fghList2_,V_,DG_,lambda_,T_,mumax_,numax_].
Parameters:
fghList1 - nested list with the results of FGH calculation for isotope 1,
fghList2 - nested list with the results of FGH calculation for isotope 2,
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states. 
Options: None.";

TableRateRfixedMorse::usage="Table with rate channel contributions (for a fixed proton donor-acceptor distance and Morse potentials).
TableRateRfixedMorse[V_,DG_,lambda_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_].
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
M - mass of the particle (Dalton),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of the morse potentials (Bohr),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states. 
Options: 
M -> mass of the particle in Daltons (Default - mass of the proton). 
Overlap -> method for overlap calculation: 
Numerical (default), Analytical, or AnalyticalSym.";

TableRateRfixedFGH::usage="Table with rate channel contributions (for a fixed proton donor-acceptor distance and general potentials).
TableRateRfixedFGH[fghList_,V_,DG_,lambda_,T_,mumax_,numax_].
Parameters:
fghList - nested list with the results of FGH calculation,
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states. 
Options: None.";

TableRateRfixedFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

PartialRateUKHarmonic::usage="PartialRateUKHarmonic[V_,DG_,lambda_,MR_,FR_,f1_,f2_,d_,T_,mu_,nu_]:
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
f1 - frequency of the left (proton donor) oscillator,
f2 - frequency of the right (proton acceptor) oscillator,
MR - reduced mass of the donor-acceptor mode (Dalton),
FR - donor-acceptor mode frequency (1/cm),
d - equilibrium distance between the minima of the proton potentials (Bohr),
T - temperature (K),
mu, nu - quantum numbers for the left and right oscillators, respectively.
Options:
M->mass of the particle in Daltons (Default - mass of the proton),
Distribution->
Classical or Quantal (for R-oscillator distribution function).
Return partial rate constant in sec^-1.";

TotalRateUKHarmonic::usage="TotalRateUKHarmonic[V_,DG_,lambda_,MR_,FR_,f1_,f2_,d_,T_,mumax_,numax_]:
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
f1 - frequency of the left (proton donor) oscillator,
f2 - frequency of the right (proton acceptor) oscillator,
MR - reduced mass of the donor-acceptor mode (Dalton),
FR - donor-acceptor mode frequency (1/cm),
d - equilibrium distance between the minima of the proton potentials (Bohr),
T - temperature (K),
mumax, numax - highest quantum numbers for the left and right oscillators states, respectively.
Options:
M->mass of the particle in Daltons (Default - mass of the proton),
Distribution->
Classical or Quantal (for R-oscillator distribution function).
Return total rate constant in sec^-1.";

KIEUKHarmonic::usage="KIEUKHarmonic[V_,DG_,lambda_,MR_,FR_,f1_,f2_,d_,T_,mumax_,numax_]:
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
f1 - frequency of the left (proton donor) protium oscillator,
f2 - frequency of the right (proton acceptor) protium oscillator,
MR - reduced mass of the donor-acceptor mode (Dalton),
FR - donor-acceptor mode frequency (1/cm),
d - equilibrium distance between the minima of the proton potentials (Bohr),
T - temperature (K),
mumax, numax - highest quantum numbers for the left and right oscillators states, respectively.
Options:
M1,M2->masses of the particle in Daltons (Default - protium and deuterium),
Distribution->
Classical or Quantal (for R-oscillator distribution function).";

TableRateUKHarmonic::usage="TableRateUKHarmonic[V_,DG_,lambda_,MR_,FR_,f1_,f2_,d_,T_,mumax_,numax_]:
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
f1 - frequency of the left (proton donor) oscillator,
f2 - frequency of the right (proton acceptor) oscillator,
MR - reduced mass of the donor-acceptor mode (Dalton),
FR - donor-acceptor mode frequency (1/cm),
d - equilibrium distance between the minima of the proton potentials (Bohr),
T - temperature (K),
mumax, numax - highest quantum numbers for the left and right oscillators states, respectively.
Options:
M->mass of the particle in Daltons (Default - protium),
Distribution->
Classical or Quantal (for R-oscillator distribution function).";

PartialRateUKMorse::usage="PartialRateUKMorse[V_,DG_,lambda_,MR_,FR_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mu_,nu_]:
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
DE1, DE2 - Dissociation energies for reactant and product Morse potentials (kcal/mol),
Beta1, Beta2 - Beta parameters for reactant and product Morse potenmtials (1/Bohr),
MR - reduced mass of the donor-acceptor mode (Dalton),
FR - donor-acceptor mode frequency (1/cm),
d - equilibrium distance between the minima of the proton potentials (Bohr),
T - temperature (K),
mu, nu - quantum numbers for the left and right oscillators, respectively.
Options:
M->mass of the particle in Daltons (Default - mass of the proton),
Overlap->Numerical/Analytical (for Morse overlaps),
Distribution->
Classical or Quantal (for R-oscillator distribution function).
Return partial rate constant in sec^-1.";

TotalRateUKMorse::usage="TotalRateUKMorse[V_,DG_,lambda_,MR_,FR_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_]:
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
DE1, DE2 - Dissociation energies for reactant and product Morse potentials (kcal/mol),
Beta1, Beta2 - Beta parameters for reactant and product Morse potenmtials (1/Bohr),
MR - reduced mass of the donor-acceptor mode (Dalton),
FR - donor-acceptor mode frequency (1/cm),
d - equilibrium distance between the minima of the proton potentials (Bohr),
T - temperature (K),
mumax, numax - highest quantum numbers for the left and right oscillators states, respectively.
Options:
M->mass of the particle in Daltons (Default - mass of the proton),
Overlap->Numerical/Analytical (for Morse overlaps),
Distribution->
Classical or Quantal (for R-oscillator distribution function).
Return total rate constant in sec^-1.";

KIEUKMorse::usage="KIEUKMorse[V_,DG_,lambda_,MR_,FR_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_]:
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
DE1, DE2 - Dissociation energies for reactant and product Morse potentials (kcal/mol),
Beta1, Beta2 - Beta parameters for reactant and product Morse potenmtials (1/Bohr),
MR - reduced mass of the donor-acceptor mode (Dalton),
FR - donor-acceptor mode frequency (1/cm),
d - equilibrium distance between the minima of the proton potentials (Bohr),
T - temperature (K),
mumax, numax - highest quantum numbers for the left and right oscillators states, respectively.
Options:
M1,M2->masses of the particle in Daltons (Default - protium and deuterium),
Overlap->Numerical/Analytical (for Morse overlaps),
Distribution->
Classical or Quantal (for R-oscillator distribution function).";

TableRateUKMorse::usage="TableRateUKHarmonic[V_,DG_,lambda_,MR_,FR_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_]:
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
DE1, DE2 - Dissociation energies for reactant and product Morse potentials (kcal/mol),
Beta1, Beta2 - Beta parameters for reactant and product Morse potenmtials (1/Bohr),
MR - reduced mass of the donor-acceptor mode (Dalton),
FR - donor-acceptor mode frequency (1/cm),
d - equilibrium distance between the minima of the proton potentials (Bohr),
T - temperature (K),
mumax, numax - highest quantum numbers for the left and right oscillators states, respectively.
Options:
M->mass of the particle in Daltons (Default - protium),
Overlap->Numerical/Analytical (for Morse overlaps),
Distribution->
Classical or Quantal (for R-oscillator distribution function).";

PartialRateRhighTMorse::usage="Partial nonadiabatic first order rate constant k(mu->nu) (1/sec) in high-T limit for proton donor-acceptor mode (Morse potentials).
PartialRateRhighT[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mu_,nu_].
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of the morse potentials (Bohr),
T - temperature (K),
mu, nu - quantum numbers for the left and right oscillators, respectively. 
Options: 
M -> mass of the particle in Daltons (Default - mass of the proton). 
Overlap -> method for overlap calculation, 
Numerical (default), Analytical, or AnalyticalSym, 
Alpha -> method for alpha calculation, 
Numerical (default), Analytical, AnalyticalSym, or <value> for fixed input.";

PartialRateRhighTFGH::usage="Partial nonadiabatic first order rate constant k(mu->nu) (1/sec) in high-T limit for proton donor-acceptor mode (general potentials).
PartialRateRhighTFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mu_,nu_].
Parameters:
fghList - nested list with the results of FGH calculation,
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
T - temperature (K),
mu, nu - quantum numbers for the left and right oscillators, respectively. 
Options: 
Alpha -> method for alpha calculation, 
Numerical (default), or <value> for fixed input.";

PartialRateRhighTFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

TotalRateRhighTMorse::usage="Total nonadiabatic first order rate constant (1/sec) in high-T limit for proton donor-acceptor mode (Morse potentials).
TotalRateRhighT[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_].
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of the morse potentials (Bohr),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states. 
Options: 
M -> mass of the particle in Daltons (Default - mass of the proton). 
Overlap -> method for overlap calculation, 
Numerical (default), Analytical, or AnalyticalSym, 
Alpha -> method for alpha calculation, 
Numerical (default), Analytical, AnalyticalSym, or <value> for fixed input.";

TotalRateRhighTFGH::usage="Total nonadiabatic first order rate constant (1/sec) in high-T limit for proton donor-acceptor mode (Morse potentials).
TotalRateRhighTFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_].
Parameters:
fghList - nested list with the results of FGH calculation,
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states. 
Options: 
Alpha -> method for alpha calculation, 
Numerical (default) or <value> for fixed input.";

TotalRateRhighTFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

KIERhighTMorse::usage="Kinetic isotope effect for nonadiabatic rate constant in high-T limit for proton donor-acceptor mode (Morse potentials).
KIERhighT[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_].
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of the morse potentials (Bohr),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states. 
Options: 
M1 -> mass of the particle for the rate in the numerator (default: mass of the Protium),
M2 -> mass of the particle for the rate in the denominator (default: mass of the Deuterium),
Overlap -> method for overlap calculation, 
Numerical (default), Analytical, or AnalyticalSym, 
Alpha -> method for alpha calculation, 
Numerical (default), Analytical, AnalyticalSym, or <value> for fixed input.";

KIERhighTFGH::usage="Kinetic isotope effect for nonadiabatic rate constant in high-T limit for proton donor-acceptor mode (general potentials).
KIERhighTFGH[fghList1_,fghList2_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_].
Parameters:
fghList1 - nested list with the results of FGH calculation for isotope 1,
fghList2 - nested list with the results of FGH calculation for isotope 2,
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states. 
Options: 
Alpha -> method for alpha calculation, 
Numerical (default) or <value> for fixed input.";

TableRateRhighTMorse::usage="Table with rate channel contributions (in high-T limit for proton donor-acceptor mode) (Morse potentials).
TableRateRhighT[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_].
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor - acceptor mode (Dalton),
Freq - donor - acceptor mode frequency (1/cm),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of the morse potentials (Bohr),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states. 
Options: 
M -> mass of the particle in Daltons (Default - mass of the proton). 
Overlap -> method for overlap calculation, 
Numerical (default), Analytical, or AnalyticalSym, 
Alpha -> method for alpha calculation, 
Numerical (default), Analytical, AnalyticalSym, or <value> for fixed input.";

TableRateRhighTFGH::usage="Table with rate channel contributions (in high-T limit for proton donor-acceptor mode) (general potentials).
TableRateRhighTFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_].
Parameters:
fghList - nested list with the results of FGH calculation,
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor - acceptor mode (Dalton),
Freq - donor - acceptor mode frequency (1/cm),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states. 
Options: 
Alpha -> method for alpha calculation, 
Numerical (default) or <value> for fixed input.";

TableRateRhighTFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

PartialRateRlowTMorse::usage="Partial nonadiabatic first order rate constant k(mu->nu) (1/sec) in low-T limit for proton donor-acceptor mode (Morse potentials).
PartialRateRlowTMorse[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mu_,nu_].
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of the morse potentials (Bohr),
T - temperature (K),
mu, nu - quantum numbers for the left and right oscillators, respectively. 
Options: 
M -> mass of the particle in Daltons (Default - mass of the proton). 
Overlap -> method for overlap calculation, 
Numerical (default), Analytical, or AnalyticalSym, 
Alpha -> method for alpha calculation, 
Numerical (default), Analytical, AnalyticalSym, or <value> for fixed input.";

PartialRateRlowTFGH::usage="Partial nonadiabatic first order rate constant k(mu->nu) (1/sec) in low-T limit for proton donor-acceptor mode (general potentials).
PartialRateRlowTFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mu_,nu_].
Parameters:
fghList - nested list with the results of FGH calculation,
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
T - temperature (K),
mu, nu - quantum numbers for the left and right oscillators, respectively. 
Options: 
Alpha -> method for alpha calculation, 
Numerical (default) or <value> for fixed input.";

PartialRateRlowTFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

TotalRateRlowTMorse::usage="Total nonadiabatic first order rate constant (1/sec) in low-T limit for proton donor-acceptor mode (Morse potentials).
TotalRateRlowT[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_].
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of the morse potentials (Bohr),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states.
Options: 
M -> mass of the particle in Daltons (Default - mass of the proton). 
Overlap -> method for overlap calculation, 
Numerical (default), Analytical, or AnalyticalSym, 
Alpha -> method for alpha calculation, 
Numerical (default), Analytical, AnalyticalSym, or <value> for fixed input.";

TotalRateRlowTFGH::usage="Total nonadiabatic first order rate constant (1/sec) in low-T limit for proton donor-acceptor mode (general potentials).
TotalRateRlowTFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_].
Parameters:
fghList - nested list with the results of FGH calculation,
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states.
Options: 
Alpha -> method for alpha calculation, 
Numerical (default) or <value> for fixed input.";

TotalRateRlowTFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

KIERlowTMorse::usage="Kinetic isotope effect for nonadiabatic rate constant in low-T limit for proton donor-acceptor mode (Morse potentials).
KIERlowT[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_].
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of the morse potentials (Bohr),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states. 
Options: 
M1 -> mass of the particle for the rate in the numerator (default: mass of the Protium),
M2 -> mass of the particle for the rate in the denominator (default: mass of the Deuterium),
Overlap -> method for overlap calculation, 
Numerical (default), Analytical, or AnalyticalSym, 
Alpha -> method for alpha calculation, 
Numerical (default), Analytical, AnalyticalSym, or <value> for fixed input.";

KIERlowTFGH::usage="Kinetic isotope effect for nonadiabatic rate constant in low-T limit for proton donor-acceptor mode (general potentials).
KIERlowTFGH[fghList1_,fghList2_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_].
Parameters:
fghList1 - nested list with the results of FGH calculation for isotope 1,
fghList2 - nested list with the results of FGH calculation for isotope 2,
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor-acceptor mode (Dalton),
Freq - donor-acceptor mode frequency (1/cm),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states. 
Options: 
Alpha -> method for alpha calculation, 
Numerical (default) or <value> for fixed input.";

TableRateRlowTMorse::usage="Table with rate channel contributions (in low-T limit for proton donor-acceptor mode) (Morse potentials).
TableRateRlowTMorse[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_].
Parameters:
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor - acceptor mode (Dalton),
Freq - donor - acceptor mode frequency (1/cm),
DE1 - dissociation energy for the oscillator on the left (kcal/mol),
Beta1 - beta parameter for the oscillator on the left (1/Bohr),
DE2 - dissociation energy for the oscillator on the right (kcal/mol),
Beta2 - beta parameter for the oscillator on the right (1/Bohr),
d - equilibrium distance between the minima of the morse potentials (Bohr),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states.
Options: 
M -> mass of the particle in Daltons (Default - mass of the proton). 
Overlap -> method for overlap calculation, 
Numerical (default), Analytical, or AnalyticalSym, 
Alpha -> method for alpha calculation, 
Numerical (default), Analytical, AnalyticalSym, or <value> for fixed input.";

TableRateRlowTFGH::usage="Table with rate channel contributions (in low-T limit for proton donor-acceptor mode) (general potentials).
TableRateRlowTFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_].
Parameters:
fghList - nested list with the results of FGH calculation,
V - electronic coupling (kcal/mol),
DG - reaction free energy (free energy bias) (kcal/mol),
lambda - solvent reorganization energy (kcal/mol),
MR - reduced mass of the donor - acceptor mode (Dalton),
Freq - donor - acceptor mode frequency (1/cm),
T - temperature (K),
mumax - highest quantum number for the reactant vibronic states,
numax - highest quantum number for the product vibronic states.
Options: 
Alpha -> method for alpha calculation, 
Numerical (default) or <value> for fixed input.";

TableRateRlowTFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

FGHWavefunctions::usage="Fourier Grid Hamiltonian calculations of the wavefunctions and overlap integrals.
FGHWavefunctions[file_,nGridPoints_,nStates_].
Parameters:
file - name of the file with the potential data (three columns: x(A) Ureactant(Hartree) Uproduct(Hartree));
nGridPoints - number of grid points (should be power of two);
nStates - number of vibrational states in each well to calculate.
Options:
M -> mass of the particle in Daltons (default is mass of a proton).
Output:
Returns a nested list of the following structure:
(
 (reactantEnergies_list, reactantPotential_matrix, reactantWavefunctions_matrix),
 (productEnergies_list, productPotential_matrix, productWavefunctions_matrix),
 overlap_matrix,
 alpha_matrix
)
where
xxxEnergies_list - is a list (vector) of length nStates;
xxxPotential_matrix - a 2D-list (table) of the splined potential [nGridPoints by 2];
xxxWavefunctions_matrix - a 2D-list (matrix)  [nStates by nGridPoints];
overlap_matrix - 2D-list (matrix) of overlap integrals [nStates by nStates];
alpha_matrix - 2D-list (matrix) of alpha-parameters [nStates by nStates].";

FGHWavefunctions::inputdata="The file `1` does not exist";
FGHWavefunctions::noexec="Executable not found";
BsplineCommand="/usr/local/bspline/bin/fgh_bspline.bin";

FGHVibronicCouplings::usage="Fourier Grid Hamiltonian calculations of the wavefunctions, overlap integrals, and vibronic couplings (electronically nonadiabatic, electronically adiabatic, and semiclassical approximations).
FGHVibronicCouplings[file_,nGridPoints_,nStates_].
Parameters:
file - name of the file with the potential data (four columns: x(A) Ureactant(kcal/mol) Uproduct(kcal/mol) ElectronicCoupling (kcal/mol));
nGridPoints - number of grid points (should be power of two);
nStates - number of vibrational states in each well to calculate.
Options:
M -> mass of the particle in Daltons (Default - mass of a proton).
Output:
Returns a nested list of the following structure:
(
 (reactantEnergiesList, reactantPotentialMatrix, reactantWavefunctionsMatrix),
 (productEnergiesList, productPotentialMatrix, productWavefunctionsMatrix),
 overlapMatrix,
 diabaticCouplingMatrix,
 productCouplingMatrix,
 semiclassicalCouplingMatrix,
 halfTunnelingSplittingMatrix,
 taueMatrix,
 taupMatrix,
 adiabaticityParameterMatrix,
 kappaMatrix
)
where
xEnergiesList - is a list (vector) of length nStates;
xPotentialMatrix - a 2D-list (table) of the splined potential [nGridPoints by 2];
xWavefunctionsMatrix - a 2D-list (matrix)  [nStates by nGridPoints];
overlapMatrix - 2D-list (matrix) of overlap integrals [nStates by nStates];
diabaticCouplingMatrix - 2D-list (matrix) of vibronic couplings \!\(\*SubscriptBox[
StyleBox[\"V\",\nFontSlant->\"Italic\"], \(\[Mu]\[Nu]\)]\)=<\!\(\*SubscriptBox[\(\[Chi]\), \(\[Mu]\)]\)|\!\(\*SuperscriptBox[
StyleBox[\"V\",\nFontSlant->\"Italic\"], \(el\)]\)(x)|\!\(\*SubscriptBox[\(\[Chi]\), \(\[Nu]\)]\)> (kcal/mol) [nStates by nStates].;
productCouplingMatrix - 2D-list (matrix) of vibronic couplings \!\(\*SubscriptBox[
StyleBox[\"V\",\nFontSlant->\"Italic\"], \(\[Mu]\[Nu]\)]\)=\!\(\*SuperscriptBox[
StyleBox[\"V\",\nFontSlant->\"Italic\"], \(el\)]\)<\!\(\*SubscriptBox[\(\[Chi]\), \(\[Mu]\)]\)|\!\(\*SubscriptBox[\(\[Chi]\), \(\[Nu]\)]\)> (kcal/mol) [nStates by nStates].;
semiclassicalCouplingMatrix - 2D-list (matrix) of semiclassical vibronic couplings (kcal/mol) (Georgievskii-Stuchebrukhov) [nStates by nStates].;
halfTunnelingSplittingMatrix - 2D-list (matrix) of half tunneling splittings (kcal/mol) (electronically adiabatic approximation) [nStates by nStates].;
taueMatrix - 2D-list (matrix) of electronic transition times (fs) [nStates by nStates].;
taupMatrix - 2D-list (matrix) of proton tunneling times (fs) [nStates by nStates].;
adiabaticityParameterMatrix - 2D-list (matrix) of adiabaticity parameters p [nStates by nStates].;
kappaMatrix - 2D-list (matrix) of prefactors \[Kappa] [nStates by nStates].";

FGHVibronicCouplings::inputdata="The file `1` does not exist";
FGHVibronicCouplings::noexec="Executable not found";
KappaCouplingCommand="/usr/local/kappa_coupling/bin/kappa_coupling.bin";

GridSolverND::usage="GridSolverND[vGrid, delta, n]
The input parameters:
vGrid - N-dimensional matrix of the potential values (atomic units) on the N-dimensional grid;
delta - grid spacing (in Bohrs), the same in each dimension;
n - number of lowest eigenvalues (and eigenvectors) desired;
Options:
M -> mass of the particle in Daltons (default is mass of a proton, MassH)
Output:
Returns a nested list of the following structure (energies are in atomic units):
(
   (energyLevel_1, energyLevel_2, ..., energyLevel_n),
   (
      (wavefunctionGrid_1, wavefunctionGrid_2, ..., wavefunctionGrid_n)
   )
)
";

GenGridSolverND::usage="GenGridSolverND[vGrid, delta, n]
The input parameters:
vGrid - N-dimensional matrix of the potential values (atomic units) on the N-dimensional grid;
delta - grid spacing (in Bohrs), the same in each dimension;
n - number of lowest eigenvalues (and eigenvectors) desired;
Options:
M -> mass of the particle in Daltons (default is mass of a proton)
DiffOrder -> accuracy order for the numerical second derivative on the grid (default is 4)
Output:
Returns a nested list of the following structure (energies are in atomic units):
{
   (energyLevel_1, energyLevel_2, ..., energyLevel_n),
   (
      (wavefunctionGrid_1, wavefunctionGrid_2, ..., wavefunctionGrid_n)
   )
)
";

FGHMax::usage="FGHMax[potx,nGrid,nRoots]
Input parameters:
potx - table with the 1D potential: {x,U(x)}, x in Angstroms, U(X) in Hartrees;
nGrid - number of grid points (arbitrary);
nRoots - number of lowest eigenstates inquired.
Options:
M -> particle mass in Daltons (amu).
Output:
Returns a nested list of the following structure (energies in kcal/mol, coordinates in Angstroms):
{
 {energy(1), energy(2),...,energy(nRoots)},
 {{x(1),pot(1)},{x(2),pot(2)},...,{x(nGrid),pot(nGrid)}},
 {
   {{x(1),wf[1](1)},{x(2),wf[1](2)},...,{x(nGrid),wf[1](nGrid)}},...,
   {{x(1),wf[nRoots](1)},{x(2),wf[nRoots](2)},...,{x(nGrid),wf[nRoots](nGrid)}}
 },
}
";

Fermi::usage="Fermi[eps,T]
Parameters:
eps - epsilon in Hartrees;
T - temperature in Kelvin;
Returns fermi distribution value at a prescribed energy and temperature";

AlphaHarmonic::usage="AlphaHarmonic[mu,nu,omega1,omega2,d]
Parameters:
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
d - distance between the minima of the displaced potentials in bohr;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Returns the alpha parameter in inverse Bohr for the harmonic oscillator approximation";

LambdaR::usage="LambdaR[MR,omega,dR]
Parameters:
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equilibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
Returns the R mode reorganization energy in kcal/mol";

NLambdaMorse::usage="NLambdaMorse[lambda,MR,omega,dR,DE1,Beta1,DE2,Beta2,mu,nu,d]
Parameters:
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in Bohrs;
Options:
M -> mass of the particle in Daltons (Default - mass of a proton);
Step -> step to take in numerical differentiaion (Default - 0.0001);
Returns the total reorganization energy in kcal/mol";

NLambdaMorsePlus::usage="NLambdaMorsePlus[lambda,MR,omega,dR,DE1,Beta1,DE2,Beta2,mu,nu,d]
Parameters:
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - the temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Step -> step to take in numerical differentiaion (Default - 0.0001);
Returns the adjusted reorganization energy for oxidation processes in kcal/mol";

NLambdaMorseMinus::usage="NLambdaMorseMinus[lambda,MR,omega,dR,DE1,Beta1,DE2,Beta2,mu,nu,d]
Parameters:
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - the temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Step -> step to take in numerical differentiaion (Default - 0.0001);
Returns the adjusted reorganization energy for reductive processes in kcal/mol";

LambdaAlphaHarmonic::usage="LambdaAlphaHarmonic[MR,mu,nu,omega1,omega2,d]
Parameters:
MR - reduced mass of the donor-acceptor mode in Daltons;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
d - distance between the minima of the displaced potentials in bohr;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Returns the coupling reorganization energy in kcal/mol for the harmonic oscillator approximation";

LambdaHarmonic::usage="LambdaHarmonic[lambda,MR,dR,omega1,omega2,mu,nu,d]
Paramters:
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Returns the total reorganization energy in kcal/mol for the harmonic oscillator approximation";

LambdaHarmonicPlus::usage="LambdaHarmonicPlus[lambda,MR,dR,omega,d,omega1,omega2,mu,nu,T]
Parameters:
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
d - distance between the minima of the displaced potentials in bohr;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
T - temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Returns the adjusted reorganization energy in kcal/mol for the oxidative process for the harmonic oscillator approximation";

LambdaHarmonicMinus::usage="LambdaHarmonicMinus[lambda,MR,dR,omega,d,omega1,omega2,mu,nu,T]
Parameters:
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
d - distance between the minima of the displaced potentials in bohr;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
T - temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Returns the adjusted reorganization energy in kcal/mol for the reductive process for the harmonic oscillator approximation";

BoltzmannMorse::usage="BoltzmannMorse[Beta1,DE1,mu,T]
Parameters:
Beta1 - beta parameter for oscillator in inverse Bohr;
DE1 - dissociation energy of the oscillator in kcal/mol;
mu - quantum number;
T - temperature in Kelvin;
Returns the Boltzmann weight per quantum number for a Morse approximation";

BoltzmannHarmonic::usage="BoltzmannHarmonic[omega1,mu,T]
Parameters:
omega1 - frequency of the harmonic oscillator in wavenumbers;
mu - quantum number;
T - temperature in Kelvin;
Returns the Boltzmann weight per quantum number for a harmonic approximation";

dGMorse::usage="dGMorse[eta,eps,Beta1,DE1,Beta2,DE2,T,mu,nu]
Parameters:
eta - overpotential in volts;
eps - energy in Hartrees (used in Fermi distribution);
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
T - temperature in Kelvin;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
Returns the reaction free energy for a Morse approximation in Hartrees";

dGHarmonic::usage="dGHarmonic[eta,eps,omega1,omega2,T,mu,nu]
Parameters:
eta - overpotential in volts;
eps - energy in Hartrees (used in Fermi distribution);
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
Returns the reaction free energy for a harmonic approximation in Hartrees";

AnodicRateConstantMorse::usage="AnodicRateConstantMorse[Vel,eta,lambda,MR,omega,dR,Beta1,DE1,Beta2,DE2,mumax,numax,d,T]
Parameters:
Vel - electronic coupling in kcal/mol;
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Step -> step to take in numerical differentiaion (Default - 0.0001);
Returns the anodic (oxidation) PCET rate constant (sec^-1) for Morse proton potentials";

CathodicRateConstantMorse::usage="AnodicRateConstantMorse[Vel,eta,lambda,MR,omega,dR,Beta1,DE1,Beta2,DE2,mumax,numax,d,T]
Parameters:
Vel - electronic coupling in kcal/mol;
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Step -> step to take in numerical differentiaion (Default - 0.0001);
Returns the cathodic (reduction) PCET rate constant (sec^-1) for Morse proton potentials";

AnodicRateConstantHarmonic::usage="AnodicRateConstantHarmonic[Vel,eta,lambda,MR,dR,omega,omega1,omega2,mumax,numax,d,T]
Parameters:
Vel - electronic coupling in kcal/mol;
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Returns the anodic (oxidation) PCET rate constant (sec^-1) for harmonic proton potentials";

CathodicRateConstantHarmonic::usage="CathodicRateConstantHarmonic[Vel,eta,lambda,MR,dR,omega,omega1,omega2,mumax,numax,d,T]
Parameters:
Vel - electronic coupling in kcal/mol;
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Returns the cathodic (reduction) PCET rate constant (sec^-1) for harmonic proton potentials";

AnodicRateConstantFGH::usage="AnodicRateConstantFGH[fghlist,Vel,eta,lambda,MR,dR,omega,mumax,numax,d,T]
Parameters:
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
Vel - electronic coupling in kcal/mol;
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Returns the anodic (oxidation) PCET rate constant (sec^-1) for arbitrary proton potentials";

AnodicRateConstantFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

CathodicRateConstantFGH::usage="CathodicRateConstantFGH[fghlist,Vel,eta,lambda,MR,dR,omega,mumax,numax,d,T]
Parameters:
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
Vel - electronic coupling in kcal/mol;
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Returns the cathodic (reduction) PCET rate constant (sec^-1) for arbitrary proton potentials";

CathodicRateConstantFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

AnodicTransferCoefficientMorse::usage="Coming soon!";

CathodicTransferCoefficientMorse::usage="Coming soon!";

AnodicTransferCoefficientHarmonic::usage="Coming soon!";

CathodicTransferCoefficientHarmonic::usage="Coming soon!";

AnodicTransferCoefficientFGH::usage="Coming soon!";

CathodicTransferCoefficientFGH::usage="Coming soon!";

AnodicCurrentDensityMorse::usage="AnodicCurrentDensityMorse[eta,lambda,MR,omega,dR,Beta1,DE1,Beta2,DE2,mumax,numax,d,T]
Parameters:
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Step -> step to take in numerical differentiaion (Default - 0.0001);
Returns the current density for oxidation in the Morse approximation";

CathodicCurrentDensityMorse::usage="CathodicCurrentDensityMorse[eta,lambda,MR,omega,dR,Beta1,DE1,Beta2,DE2,mumax,numax,d,T]
Parameters:
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Step -> step to take in numerical differentiaion (Default - 0.0001);
Returns the current density for reduction in the Morse approximation";

TotalCurrentDensityMorse::usage="TotalCurrentDensityMorse[eta,lambda,MR,omega,dR,Beta1,DE1,Beta2,DE2,mumax,numax,d,T,Cox,Cred]
Parameters:
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Cox - concentration of the oxidized complex (any consistent units);
Cred - concentration of the reduced complex (any consistent units);
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Step -> step to take in numerical differentiaion (Default - 0.0001);
Returns the total current density for the Morse approximation";

AnodicCurrentDensityHarmonic::usage="AnodicCurrentDensityHarmonic[eta,lambda,MR,dR,omega,omega1,omega2,mumax,numax,d,T]
Parameters:
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Returns the current density for oxidation in the harmonic approximation";

CathodicCurrentDensityHarmonic::usage="CathodicCurrentDensityHarmonic[eta,lambda,MR,dR,omega,omega1,omega2,mumax,numax,d,T]
Parameters:
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the donor-acceptor mode in wavenumbers;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Returns the current density for reduction in the harmonic approximation";

TotalCurrentDensityHarmonic::usage="TotalCurrentDensityHarmonic[eta,lambda,MR,dR,omega,omega1,omega2,mumax,numax,d,T,Cox,Cred]
Parameters:
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Cox - concentration of the oxidized complex (any consistent units);
Cred - concentration of the reduced complex (any consistent units);
M -> mass of the particle in Daltons (Default - mass of the proton);
Returns the total current density for the harmonic approximation";

AnodicCurrentDensityLowTMorse::usage="AnodicCurrentDensityLowTMorse[eta,lambda,MR,omega,dR,Beta1,DE1,Beta2,DE2,mumax,numax,d,T]
Parameters:
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Step -> step to take in numerical differentiaion (Default - 0.0001);
Returns the current density for oxidation in the Morse approximation for Low T";

CathodicCurrentDensityLowTMorse::usage="CathodicCurrentDensityLowTMorse[eta,lambda,MR,omega,dR,Beta1,DE1,Beta2,DE2,mumax,numax,d,T]
Parameters:
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Step -> step to take in numerical differentiaion (Default - 0.0001);
Returns the current density for reduction in the Morse approximation for Low T";

TotalCurrentDensityLowTMorse::usage="TotalCurrentDensityLowTMorse[eta,lambda,MR,omega,dR,Beta1,DE1,Beta2,DE2,mumax,numax,d,T,Cox,Cred]
Parameters:
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Cox - concentration of the oxidized complex (any consistent units);
Cred - concentration of the reduced complex (any consistent units);
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Step -> step to take in numerical differentiaion (Default - 0.0001);
Returns the total current density for the Morse approximation for Low T";

AnodicCurrentDensityLowTHarmonic::usage="AnodicCurrentDensityLowTHarmonic[eta,lambda,MR,dR,omega,omega1,omega2,mumax,numax,d,T]
Parameters:
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Returns the current density for oxidation in the harmonic approximation for Low T";

CathodicCurrentDensityLowTHarmonic::usage="CathodicCurrentDensityLowTHarmonic[eta,lambda,MR,dR,omega,omega1,omega2,mumax,numax,d,T]
Parameters:
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
M -> mass of the particle in Daltons (Default - mass of the proton);
Returns the current density for reduction in the harmonic approximation for Low T";

TotalCurrentDensityLowTHarmonic::usage="CurrentDensityTotalHarmonicT[eta,lambda,MR,dR,omega,omega1,omega2,mumax,numax,d,T,Cox,Cred]
Parameters:
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Cox - concentration of the oxidized complex (any consistent units);
Cred - concentration of the reduced complex (any consistent units);
M -> mass of the particle in Daltons (Default - mass of the proton);
Returns the total current density for the harmonic approximation for Low T";

AnodicCurrentDensityFGH::usage="AnodicCurrentDensityFGH[fghlist,eta,lambda,MR,dR,omega,mumax,numax,d,T]
Parameters:
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
Returns the current density for oxidation for general proton potentials";

AnodicCurrentDensityFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

CathodicCurrentDensityFGH::usage="CathodicCurrentDensityFGH[fghlist,eta,lambda,MR,dR,omega,mumax,numax,d,T]
Parameters:
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
Returns the current density for reduction for general proton potentials";

CathodicCurrentDensityFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

TotalCurrentDensityFGH::usage="TotalCurrentDensityFGH[fghlist,eta,lambda,MR,dR,omega,mumax,numax,d,T,Cox,Cred]
Parameters:
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Cox - concentration of the oxidized complex (any consistent units);
Cred - concentration of the reduced complex (any consistent units);
Returns the total current density for general proton potentials";

AnodicCurrentDensityLowTFGH::usage="AnodicCurrentDensityLowTFGH[fghlist,eta,lambda,MR,dR,omega,mumax,numax,d,T]
Parameters:
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
Returns the current density for oxidation for general proton potentials in LowT";

AnodicCurrentDensityLowTFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

CathodicCurrentDensityLowTFGH::usage="CathodicCurrentDensityLowTFGH[fghlist,eta,lambda,MR,dR,omega,mumax,numax,d,T]
Parameters:
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Options:
Returns the current density for reduction for general proton potentials in LowT";

CathodicCurrentDensityLowTFGH::numberofstates=
"Number of reactant or product states chosen exceeds number of states available in fghList.";

TotalCurrentDensityLowTFGH::usage="TotalCurrentDensityLowTFGH[fghlist,eta,lambda,MR,dR,omega,mumax,numax,d,T,Cox,Cred]
Parameters:
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Cox - concentration of the oxidized complex (any consistent units);
Cred - concentration of the reduced complex (any consistent units);
Returns the total current density for general proton potentials in LowT";


Pathway1::usage="Pathway1[Eta,e1,pKa1,e2,Lambda1,Lambda2,pKaAcid]
Returns Pathway1 free energy diagram. red lines denotes monometallic mechanism,
blue lines denotes bimetallic mechanism,black lines denote states that are applicable
to both mechanisms.
Parameters:
Eta is applied potiential (eV);
e1 is Co(II)/Co(I) reduction potiential(ev);
pKa1 is pKa of Co(III)H;
e2 is Co(III)/Co(II) reduction potiential(ev);
Lambda1 is Co(II) to Co(I) reorganization energy(eV);
Lambda2 is Co(III) to Co(II) reorganization energy(eV);
pKaAcid is Acid pKa;
Options:
Ref -> reference electrode (Default - HA/H2)";


Pathway2::usage="Pathway2[Eta,e1,pKa1,e3,pKa2,e4,Lambda1,Lambda3,pKaAcid]
Returns Pathway2 free energy diagram. red lines denotes monometallic mechanism,
blue lines denotes bimetallic mechanism,black lines denote states that are applicable
to both mechanisms.
Parameters:
Eta is applied potiential (eV);
e1 is Co(II)/Co(I) reduction potiential(ev);
pKa1 is pKa of Co(III)H;
e3 is Co(III)H/Co(II)H reduction potiential(ev);
pKa2 is pKa of Co(II)H;
e3 is Co(III)H/Co(II)H reduction potiential(ev);
Lambda1 is Co(II) to Co(I) reorganization energy(eV);
Lambda3 is Co(III)H to Co(II)H reorganization energy(eV);
pKaAcid is Acid pKa;
Options:
Ref -> reference electrode (Default - HA/H2)";


Pathway3::usage="Pathway3[Eta,e1,Lambda1,e4,Lambda4,pKa2,pKaAcid]
Returns Pathway3 free energy diagram. red lines denotes monometallic mechanism,
blue lines denotes bimetallic mechanism,black lines denote states that are applicable
to both mechanisms.
Parameters:
Eta is applied potiential (eV),
e1 is Co(II)/Co(I) reduction potiential(ev);
Lambda1 is Co(II) to Co(I) reorganization energy(eV);
e4 is Co(I)/Co(0) reduction potiential(ev);
Lambda4 is Co(I) to Co(0) reorganization energy(eV);
pKa2 is pKa of Co(II)H;
pKaAcid is Acid pKa,
Options:
Ref -> reference electrode (Default - HA/H2)";


(*
Constants
*)

r=0.0000861688;
t=298.15;

(*
Unit conversion factors:
*)

bohr2a=0.529177`10;
a2bohr=1/bohr2a;
au2kcal=627.5095`10;
kcal2au=1/au2kcal;
ev2cm=8065.54477`10;
cm2ev = 1/ev2cm ;
au2ev=27.2113`10;
ev2au = 1/au2ev;
cm2au=cm2ev*ev2au;     
au2cm=1/cm2au;
au2ps=2.4189`10*10^(-5);
ps2au=1/au2ps;
au2fs=au2ps*1000;
fs2au=1/au2fs;

(*
Physical constants:
*)

hbar=1;
hbarps=0.047685388`10/(au2kcal*\[Pi]);
Dalton=1822.8900409014022`15;
MassH=1.0072756064562605`15;
MassD=2.0135514936645316`15;
MassT=3.0160492`15;
MassE=1/Dalton;
mH=MassH*Dalton;
mD=MassD*Dalton;
mT=MassT*Dalton;
mE=1;
kb=3.16683`10*10^-6;

strcm=ToString[Superscript["cm","-1"],TraditionalForm];
strr0=ToString[Subscript[R,0],TraditionalForm];
strpmu=ToString[Subscript[P,"\[Mu]"],TraditionalForm];
stramunu=ToString[Subscript["\[Alpha]","\[Mu]\[Nu]"],TraditionalForm];
strdg00=ToString[Subscript[\[CapitalDelta]G,"00"],TraditionalForm];
strdg000=ToString[Subsuperscript[\[CapitalDelta]G,"00","0"],TraditionalForm];
strdgmunu=ToString[Subscript[\[CapitalDelta]G,"\[Mu]\[Nu]"],TraditionalForm];
strdg0munu=ToString[Subsuperscript[\[CapitalDelta]G,"\[Mu]\[Nu]","0"],TraditionalForm];
strdgdd=ToString[Superscript[\[CapitalDelta]G,"\[DoubleDagger]"],TraditionalForm];
strdgddmunu=ToString[Subsuperscript[\[CapitalDelta]G,"\[Mu]\[Nu]","\[DoubleDagger]"],TraditionalForm];
strs2=ToString[Superscript[S,2],TraditionalForm];
strs2munu=ToString[Subsuperscript[S,"\[Mu]\[Nu]",2],TraditionalForm];

Begin["`Private`"]

(*
Harmonic potentials
*)

(*
Potential (in kcal/mol):
M - mass (Daltons);
Freq - frequency (1/cm);
x - coordinate (Bohr) (potential minimum is at x=0).
*)

HarmonicPotential[Freq_,x_,OptionsPattern[{M->MassH}]]:=Module[{mu,f,pot},
	mu=OptionValue[M]*Dalton;
	f=Freq*cm2au;
	pot=mu*f^2*x^2/2;
	pot*au2kcal
];

(*
Morse potentials
*)

(*
Potential (in kcal/mol):
DE - dissociation energy (kcal/mol);
beta - beta parameter (1/Bohr);
x - coordinate (Bohr) (potential minimum is at x=0).
*)

MorsePotential[DE_,beta_,x_]:=Module[{DEE,pot},
	DEE=DE*kcal2au;
	pot=DEE*(1-Exp[-beta*x])^2;
	pot*au2kcal
];

MorsePotentialInverted[DE_,beta_,x_]:=Module[{DEE,pot},
	DEE=DE*kcal2au;
	pot=DEE*(1-Exp[beta*x])^2;
	pot*au2kcal
];

(*
Beta parameter (1/Bohr):
M - mass of the particle (Daltons),
Omega - frequency (1/cm),
DE - dissociation energy (kcal/mol).
*)

MorseBeta[Omega_,DE_,OptionsPattern[{M->MassH}]]:=Module[{mu,DEA,wA},
	mu=OptionValue[M]*Dalton;
	DEA=DE*kcal2au;
	wA=Omega*cm2au;
	wA*Sqrt[mu/(2*DEA)]
];

(*
Energies and wavefunctions for a Morse oscillator
from [J. P. Dahl, M. Springborg, J. Chem. Phys. 88, 4535 (1987)] 
(typo in Eq.(38): factor Sqrt[Alpha] is missing)
*)

(*
Energies of the discrete states
*)

(*
n = 0, 1, 2,... - quantum number;
Beta - Morse parameter (in 1/Bohr);
DE - dissociation energy (in kcal/mol);
M - mass of the particle (in Daltons);
MorseEnergy returns energy in kcal/mol relative to the bottom of the potential well.
*)

MorseEnergy[n_,Beta_,DE_,OptionsPattern[{M->MassH}]]:=Module[{DEA,Lambda,mu},
	DEA=DE*kcal2au;
	mu=OptionValue[M]*Dalton;
	Lambda=Sqrt[2*mu*DEA]/(hbar*Beta);
	au2kcal*hbar*((n+1/2)-1/(2*Lambda)*(n+1/2)^2)*Sqrt[(2*DEA*Beta^2)/mu]
];

(*
Number of bound states:
Beta - beta parameter (1/Bohr);
DE - dissociation energy (kcal/mol);
M - mass (Daltons).
*)

MorseBound[Beta_,DE_,OptionsPattern[{M->MassH}]]:=Module[{DEA,Lambda,mu},
	DEA=DE*kcal2au;
	mu=OptionValue[M]*Dalton;
	Lambda=Sqrt[2*mu*DEA]/(hbar*Beta);
	IntegerPart[Lambda-1/2]
];

(*
Wavefunctions
*)

(*
n = 0, 1, 2,... - quantum number;
Beta - Morse parameter (in 1/Bohr);
DE - dissociation energy (in kcal/mol);
Mu - mass of the particle (in Daltons);
x - coordinate (in Bohrs) relative to the minimum of the potential.
Returns wavefunction in units of (Bohr)^(-1/2).
*)

MorseWavefunctionLeft[n_,Beta_,DE_,x_,OptionsPattern[{M->MassH}]]:=Module[{Mu,DEA,Lambda,Xi,wnorm},
	Mu=OptionValue[M]*Dalton;
	DEA=DE*kcal2au;
	Lambda=Sqrt[2*Mu*DEA]/(hbar*Beta);
	Xi=2*Lambda*Exp[-Beta*x];
	wnorm=Sqrt[(Beta*(2*Lambda-2*n-1)*Gamma[n+1])/Gamma[2*Lambda-n]];
	wnorm*Exp[-(Xi/2)]*Xi^(Lambda-n-1/2)*LaguerreL[n,2*Lambda-2*n-1,Xi]
];

MorseWavefunctionRight[n_,Beta_,DE_,x_,OptionsPattern[{M->MassH}]]:=Module[{Mu,DEA,Lambda,Xi,wnorm},
	Mu=OptionValue[M]*Dalton;
	DEA=DE*kcal2au;
	Lambda=Sqrt[2*Mu*DEA]/(hbar*Beta);
	Xi=2*Lambda*Exp[Beta*x];
	wnorm=Sqrt[(Beta*(2*Lambda-2*n-1)*Gamma[n+1])/Gamma[2*Lambda-n]];
	wnorm*Exp[-(Xi/2)]*Xi^(Lambda-n-1/2)*LaguerreL[n,2*Lambda-2*n-1,Xi]
];

(*
Energies and wavefunctions for a Harmonic oscillator
*)

(*
freq - frequency (1/cm);
n - quantum number (0 for ground state).
Returns energy in kcal/mol;
*)

HOEnergy[freq_,n_]:=Module[{freqau,eau},
	freqau=freq*cm2au;
	eau=freqau*(n+1/2);
	eau*au2kcal
];

(*
n - quantum number (0 for ground state);
M - mass (Daltons);
Freq - frequency (1/cm);
x - coordinate (Bohrs);
(Potential minimum is at x=0).
Returns wavefunction in units of (Bohr)^(-1/2).
*)

HOWavefunction[n_,Freq_,x_,OptionsPattern[{M->MassH}]]:=Module[{mu,f,wnorm,w1,w2},
	mu=OptionValue[M]*Dalton;
	f=Freq*cm2au;
	wnorm=(mu*f/(\[Pi]*hbar))^(1/4)*2^(-n/2)*(n!)^(-1/2);
	w1=Exp[-mu*f*x^2/(2*hbar)];
	w2=HermiteH[n,x*Sqrt[mu*f/hbar]];
	wnorm*w1*w2
];

(*
Vibrational partition function for Harmonic oscillator
*)

QHarmonic[omega_,T_]:=Module[{},
	1/(Exp[(omega*cm2au)/(2*kb*T)]-Exp[-(omega*cm2au)/(2*kb*T)])
];

(*
Vibrational partition function for Morse oscillator
*)

QMorseExact[Beta_,DE_,T_,OptionsPattern[{M->MassH}]]:=Module[{mass,mumax,e0},
	mass=OptionValue[M];
	mumax=MorseBound[Beta,DE,M->mass];
	e0=MorseEnergy[0,Beta,DE,M->mass];
	Exp[-e0*kcal2au/(kb*T)]*(1 + Sum[Exp[-kcal2au*(MorseEnergy[mu,Beta,DE,M->mass]-e0)/(kb*T)],{mu,1,mumax}])
];

QMorseApprox[Beta_,DE_,T_,OptionsPattern[{M->MassH}]]:=Module[{mass,m,DEA,omega,theta,xe,a,b,nmax,e0,Q0,QN,sqa,cN},
	mass=OptionValue[M];
	m = mass*Dalton;
	DEA = DE*kcal2au;
	nmax = MorseBound[Beta,DE,M->mass];
	xe = Sqrt[Beta/(8*m*DEA)];
	omega = Sqrt[(2*DEA*Beta^2)/m];
	theta = omega/(kb*T);
	a = xe*theta;
	b = theta*(1-xe);
	e0=MorseEnergy[0,Beta,DE,M->mass]*kcal2au;
	sqa = Sqrt[a];
	cN = Exp[a*(nmax+1)^2-b*(nmax+1)];
	Q0 = 1/(1-Exp[-b]) + DawsonF[b/(2*sqa)]/sqa - 1/b;
	QN = cN*(1/(Exp[2*a*(nmax+1)-b]-1) + DawsonF[(2*a*(nmax+1)-b)/(2*sqa)]/sqa - 1/(2*a*(nmax+1)-b));
	Exp[-e0/(kb*T)]*(Q0 + QN)
];

(*
Marcus parabolas for PCET reaction with Morse potentials
*)

(*
n = 0, 1, 2,... - vibrational quantum number;
beta - Morse parameter (in 1/Bohr);
de - dissociation energy (in kcal/mol);
m - mass of the particle (in Daltons);
lambda - reorganization energy (kcal/mol);
X - energy gap reaction coordinate (in kcal/mol).
Returns free energy in kcal/mol.
*)

MarcusReactantMorse[n_,beta_,de_,lambda_,X_,OptionsPattern[{M->MassH}]]:=Module[{fe,en,e0},
	fe=1/(4*lambda)*(X+lambda)^2;
	e0=MorseEnergy[0,beta,de,M->OptionValue[M]];
	en=MorseEnergy[n,beta,de,M->OptionValue[M]];
	fe+en-e0
];

MarcusProductMorse[n_,beta_,de_,lambda_,X_,OptionsPattern[{M->MassH}]]:=Module[{fe,en,e0},
	fe=1/(4*lambda)*(X-lambda)^2;
	e0=MorseEnergy[0,beta,de,M->OptionValue[M]];
	en=MorseEnergy[n,beta,de,M->OptionValue[M]];
	fe+en-e0
];

(*
Overlap integrals between harmonic oscillator wavefunctions
*)

(*
Overlap between two harmonic oscillator wavefunctions
[J.-L. Chang, J. Mol. Spectrosc. 232 (2005) 102-104]:
n1 - quantum number for the oscillator on the left (0 for ground state);
n2 - quantum number for the oscillator on the right (0 for ground state);
f1 - frequency for the oscillator on the left (1/cm);
f2 - frequency for the oscillator on the right (1/cm);
m - mass of the particle (Daltons);
d - distance between the minima of the displaced harmonic potentials.
*)

Kaux[i_,j_,a_,b_]:=If[OddQ[i+j],0,(i+j-1)!!/(a+b)^((i+j)/2)];

HarmonicOverlap[n1_,n2_,f1_,f2_,d_,OptionsPattern[{M->MassH}]]:=Module[{o1,o2,mu,a1,a2,s,A,Prefactor,b1,b2,i,j},
	mu=OptionValue[M]*Dalton;
	o1=f1*cm2au;
	o2=f2*cm2au;
	a1=(mu*o1)/hbar;
	a2=(mu*o2)/hbar;
	s=(a1*a2*d^2)/(a1+a2);
	A=(2*Sqrt[a1*a2])/(a1+a2);
	Prefactor=Sqrt[(A*Exp[-s])/(2^(n1+n2)*n1!*n2!)];
	b1=-((d*a2*Sqrt[a1])/(a1+a2));
	b2=(d*a1*Sqrt[a2])/(a1+a2);
	Prefactor*Sum[Binomial[n1,i]*Binomial[n2,j]*HermiteH[n1-i,b1]*HermiteH[n2-j,b2]*(2*Sqrt[a1])^i*(2*Sqrt[a2])^j*Kaux[i,j,a1,a2],{j,0,n2},{i,0,n1}]
];

(*
Overlap integrals between morse oscillator wavefunctions
*)

(*
Case I: identical mirrored Morse potentials
[modified expression from J. C. Lopez, A. L. Rivera, Yu. F. Smirnov, A. Frank, arXiv:physics/0109017v]
*)

(*
AccN - accuracy (number of digits after decimal point), usually 20 is enough;
M - mass of the particle (Daltons);
DE - dissociation energy (kcal/mol);
Beta - beta parameter (1/Bohr);
nu1 - quantum number for the oscillator on the left (0 for ground state);
nu2 - quantum number for the oscillator on the right (0 for ground state);
d - distance between the minima of the displaced and mirrored Morse potentials.
*)

MorseOverlapSym[DE_,Beta_,nu1_,nu2_,d_,OptionsPattern[{Acc->100,M->MassH}]]:=
Module[{AccN,DEA,R,mu,beta,lambda,l1,l2,j1,j2,y1,y2,Fcsum,fac0a,fac0b,fac0,fac1a,fac1b,fac1,fac2,fac3,Fcaux,Lambda,res},
(*-----------------------------------------------------------------------*)
(* This program computes the Franck-Condon factor of two Morse functions *)
(* for the mirrored potentials for the case when Subscript[\[Beta], 1]=Subscript[\[Beta], 2]=\[Beta] *)
(* Modified version of Lopez-Rivera-Smirnov-Frank *)
(* Usage: *)
(* FC0H[R,nu1,nu2] *)
(*---------------------------------------------------------------------*)
	AccN=OptionValue[Acc];
	DEA=SetAccuracy[DE*kcal2au,AccN];
	beta=SetAccuracy[Beta,AccN];
	R=SetAccuracy[d,AccN];
	mu=SetAccuracy[OptionValue[M]*Dalton,AccN];
	lambda=SetAccuracy[Sqrt[2*mu*DEA]/(beta*hbar),AccN];
	j1=SetAccuracy[lambda-1/2,AccN];
	j2=SetAccuracy[lambda-1/2,AccN];
	y1=SetAccuracy[(2*j1+1)*Exp[-beta*R/2],AccN];
	y2=SetAccuracy[(2*j2+1)*Exp[-beta*R/2],AccN];
	Fcsum=0;
	fac0a=SetAccuracy[2*Sqrt[FullSimplify[nu1!*(j1-nu1)*nu2!*(j2-nu2)/(Gamma[2*j1-nu1+1]*Gamma[2*j2-nu2+1])]],AccN];
	fac0b=SetAccuracy[y1^(j1-nu1)*y2^(j2-nu2),AccN];
	fac0=fac0a*fac0b;
	Do[(*Do loops for summations*)
		fac1a=SetAccuracy[(2*j1+1)^(l1)*(2*j2+1)^(l2),AccN];
		fac1b=SetAccuracy[Exp[-beta*R*(l1)/2]*Exp[-beta*R*(l2)/2],AccN];
		fac1=SetAccuracy[fac1a*fac1b,AccN];
		fac2=SetAccuracy[(-1)^(l1+l2)/(l1!*l2!),AccN];
		fac3=SetAccuracy[FullSimplify[Binomial[2*j1-nu1,nu1-l1]*Binomial[2*j2-nu2,nu2-l2]],AccN];
		Fcaux=SetAccuracy[fac0*fac1*fac2*fac3,AccN];
		Lambda=SetAccuracy[j1-nu1+l1-(j2-nu2+l2),AccN];
		res=SetAccuracy[2(y1/y2)^(-(Lambda/2))BesselK[-Lambda,Sqrt[y1*y2]],AccN];
		Fcsum=Fcsum+SetAccuracy[Fcaux*res,AccN],{l1,0,nu1},{l2,0,nu2}
	];
	Fcsum
];

(*
Case II: general mirrored Morse potentials
[modified expression from J. C. Lopez, A. L. Rivera, Yu. F. Smirnov, A. Frank, arXiv:physics/0109017v]
*)

(*
AccN - accuracy (number of digits after decimal point), usually 20 is enough;
M - mass of the particle (Daltons);
DE1 - dissociation energy (kcal/mol);
DE2 - dissociation energy (kcal/mol);
Beta1 - beta parameter (1/Bohr);
Beta2 - beta parameter (1/Bohr);
nu1 - quantum number for the oscillator on the left (0 for ground state);
nu2 - quantum number for the oscillator on the right (0 for ground state);
d - distance between the minima of the displaced and mirrored Morse potentials.
*)

MorseOverlap[DE1_,Beta1_,DE2_,Beta2_,nu1_,nu2_,d_,OptionsPattern[{Acc->100,M->MassH}]]:=
Module[{Acc1,AccN,DE1A,DE2A,R,mu,lambda1,lambda2,l1,l2,j1,beta1,j2,beta2,y1,y2,a,Fcsum,fac0a,fac0b,fac0,fac1a,fac1b,fac1,fac2,fac3,Fcaux,Lambda,aux,peak,xatmax,fatmax,integr,res,finteg,finteg2,smallestx,x,X},
(*-----------------------------------------------------------------------*)
(* This program computes the Franck-Condon factor of two Morse functions *)
(* for the mirrored potentials *)
(* Modified version of Lopez-Rivera-Smirnov-Frank *)
(*---------------------------------------------------------------------*)
(* Precision and accuracy goals in numerical integration: *)
	AccN=OptionValue[Acc];
	Acc1=15;
	DE1A=SetAccuracy[DE1*kcal2au,AccN];
	DE2A=SetAccuracy[DE2*kcal2au,AccN];
	beta1=SetAccuracy[Beta1,AccN];
	beta2=SetAccuracy[Beta2,AccN];
	R=SetAccuracy[d,AccN];
	mu=SetAccuracy[OptionValue[M]*Dalton,AccN];
	lambda1=SetAccuracy[Sqrt[2*mu*DE1A]/(beta1*hbar),AccN];
	lambda2=SetAccuracy[Sqrt[2*mu*DE2A]/(beta2*hbar),AccN];
	j1=SetAccuracy[lambda1-1/2,AccN];
	j2=SetAccuracy[lambda2-1/2,AccN];
	a=SetAccuracy[beta2/beta1,AccN];
	y1=SetAccuracy[(2*j1+1)*Exp[-beta1*R/2],AccN];
	y2=SetAccuracy[(2*j2+1)*Exp[-beta2*R/2],AccN];
	Fcsum=0;
	fac0a=SetAccuracy[2*Sqrt[a*FullSimplify[nu1!*(j1-nu1)*nu2!*(j2-nu2)/(Gamma[2*j1-nu1+1]*Gamma[2*j2-nu2+1])]],AccN];
	fac0b=SetAccuracy[y1^(j1-nu1)*y2^(j2-nu2),AccN];
	fac0=fac0a*fac0b;
	Do[(*Do loops for summations*)
		fac1a=SetAccuracy[(2*j1+1)^(l1)*(2*j2+1)^(l2),AccN];
		fac1b=SetAccuracy[Exp[-beta1*R*(l1)/2]*Exp[-beta2*R*(l2)/2],AccN];
		fac1=SetAccuracy[fac1a*fac1b,AccN];
		fac2=SetAccuracy[(-1)^(l1+l2)/(l1!*l2!),AccN];
		fac3=SetAccuracy[FullSimplify[Binomial[2*j1-nu1,nu1-l1]*Binomial[2*j2-nu2,nu2-l2]],AccN];
		Fcaux=SetAccuracy[fac0*fac1*fac2*fac3,AccN];
		(* ===========================================================*)
		(* Integration Section *)
		(* Integral Ia.Calculation with scaling *)
		(* ===========================================================*)
		Lambda=SetAccuracy[j1-nu1+l1-a*(j2-nu2+l2),AccN];
		(* Integrand *)
		integr=SetAccuracy[x^(Lambda-1)*Exp[-(y1*x+y2*x^(-a))/2],AccN];
		finteg[xx_]:=xx^(Lambda-1)*Exp[-(y1*xx+y2*xx^(-a))/2];
		(* Find smallest x before underflow occurs *)
		Off[General::"unfl"];
		smallestx=0.1`100;
		While[finteg[smallestx]>$MinNumber,smallestx=smallestx/2];
		smallestx=smallestx*2;
		On[General::unfl];
		(*Looking for position of maximum*)
		aux=FindMaximum[integr,{x,smallestx+1,smallestx,\[Infinity]},WorkingPrecision->AccN];
		peak=aux[[2]][[1]][[2]];
		xatmax=SetAccuracy[peak,AccN];
		(*Integrand value at maximum*)
		fatmax=SetAccuracy[finteg[peak],AccN];
		finteg2=SetAccuracy[(1/fatmax)*(X*xatmax)^(Lambda-1)*Exp[-(y1*(X*xatmax)+y2*(X*xatmax)^(-a))/2],AccN];
		res=SetAccuracy[fatmax*xatmax*NIntegrate[finteg2,{X,0,Infinity},PrecisionGoal->Acc1,AccuracyGoal->Acc1,Method->"DoubleExponential"],AccN];
		(*Adding terms*)
		Fcsum=Fcsum+SetAccuracy[Fcaux*res,AccN],{l1,0,nu1},{l2,0,nu2}
	];
	(*Output*)
	Fcsum
];

(*
Case III: Numerical overlap integrals for general mirrored Morse potentials.
*)

(*
M - mass of the particle (Daltons);
DE1 - dissociation energy (kcal/mol);
DE2 - dissociation energy (kcal/mol);
Beta1 - beta parameter (1/Bohr);
Beta2 - beta parameter (1/Bohr);
nu1 - quantum number for the oscillator on the left (0 for ground state);
nu2 - quantum number for the oscillator on the right (0 for ground state);
d - distance between the minima of the displaced and mirrored Morse potentials.
*)

NMorseOverlap[DE1_,Beta1_,DE2_,Beta2_,nu1_,nu2_,d_,OptionsPattern[{M->MassH}]]:=Module[{Mass,mu,DE1A,DE2A,X,integrand},
	Mass=OptionValue[M];
	mu=Mass*Dalton;
	DE1A=DE1*kcal2au;
	DE2A=DE2*kcal2au;
	integrand[x_]:=MorseWavefunctionLeft[nu1,Beta1,DE1,x,M->Mass]*MorseWavefunctionRight[nu2,Beta2,DE2,x-d,M->Mass];
	NIntegrate[integrand[X],{X,-\[Infinity],\[Infinity]},Method->{"DoubleExponential","SymbolicProcessing"->0},Compiled->False]
];

(*
Distance dependence of overlap integrals
*)

(*
Numerical logarithmic derivative of the numerical overlap integrals:
M - mass of the particle (Dalton);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
i, j - quantum numbers for the left and right oscillators, respectively;
d - distance between the minima of the morse potentials (Bohr).
*)

NAlphaMorse[DE1_,Beta1_,DE2_,Beta2_,i_,j_,d_,OptionsPattern[{M->MassH,Step->0.0001}]]:=Module[{h,S0,S1,S2},
	h=OptionValue[Step];
	S0=NMorseOverlap[DE1,Beta1,DE2,Beta2,i,j,d,M->OptionValue[M]];
	S2=NMorseOverlap[DE1,Beta1,DE2,Beta2,i,j,d+h,M->OptionValue[M]];
	S1=NMorseOverlap[DE1,Beta1,DE2,Beta2,i,j,d-h,M->OptionValue[M]];
	-(1/S0)*(S2-S1)/(2*h)
];

(*
Analytical derivative of analytical overlap integrals (for identical mirrored potentials):
AccN - accuracy
M - mass of the particle (Dalton);
DE - dissociation energy for the oscillator (kcal/mol);
Beta - beta parameter for the oscillator (1/Bohr);
nu1, nu2 - quantum numbers for the left and right oscillators, respectively;
d - distance between the minima of the morse potentials (Bohr).
*)

AlphaMorseSym[DE_,Beta_,nu1_,nu2_,d_,OptionsPattern[{Acc->50,M->MassH}]]:=Module[{x,s},
	s[r_]=MorseOverlapSym[DE,Beta,nu1,nu2,r,Acc->OptionValue[Acc],M->OptionValue[M]];
	-(1/s[x])*D[s[x]]/.x->d
];

(*
Analytical derivative of analytical overlap integrals (for general mirrored potentials):
AccN - accuracy
M - mass of the particle (Dalton);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
nu1, nu2 - quantum numbers for the left and right oscillators, respectively;
d - distance between the minima of the morse potentials (Bohr).
*)

AlphaMorse[DE1_,Beta1_,DE2_,Beta2_,nu1_,nu2_,d_,OptionsPattern[{Acc->50,M->MassH}]]:=Module[{x,s},
	s[r_]=MorseOverlap[DE1,Beta1,DE2,Beta2,nu1,nu2,r,Acc->OptionValue[Acc],M->OptionValue[M]];
	-(1/s[x])*D[s[x]]/.x->d
];

(*
FGH - method for energies and wavefunctions of a particle in an arbitrary one - dimensional potential
*)

(*
The following function calculates energies and wavefunctions of the particle in two external one-dimensional potentials as well as the ovelap integral matrix and matrix of alpha-parameters.
*)

FGHWavefunctions[file_,nGridPoints_,nStates_,OptionsPattern[{M->MassH}]]:=
Module[{mass,massString,filename,nG,nS,error,command,reactantEnergies,productEnergies,reactantPotential,productPotential,reactantWavefunctions,productWavefunctions,overlapMatrix,alphaMatrix},
filename=StringSplit[file,"."][[1]];
nG=nGridPoints;
(* Check if nStates < nGridPoints, otherwise set nStates to default value of 10 *)
If[nStates<nGridPoints,nS=nStates,nS=10];
mass=OptionValue[M]*Dalton;
massString=ToString[FortranForm[mass]];
(* run external executable to generate all the data *)
command=BsplineCommand<>" "<>file<>" "<>ToString[nG]<>" "<>massString<>" "<>ToString[nS]<>" > /dev/null";
error=Run[command];
If[error!=0,Message[FGHWavefunctions::noexec];Abort[]];
(* Read energies from disk *)
reactantEnergies=ReadList[filename<>"_reactant_en.dat",Real];
productEnergies=ReadList[filename<>"_product_en.dat",Real];
(* Read splined potentials from disk *)
reactantPotential=Part[ReadList[filename<>"_bspline.dat",{Number,Number,Number}],All,{1,2}];
productPotential=Part[ReadList[filename<>"_bspline.dat",{Number,Number,Number}],All,{1,3}];
(* Read wavefunctions from disk *)
reactantWavefunctions={};
productWavefunctions={};
Do[reactantWavefunctions=Append[reactantWavefunctions,Part[ReadList[filename<>"_reactant_wf.dat",Table["Real",{k,1,1+nS}]],All,1+i]],{i,1,nS}];
Do[productWavefunctions=Append[productWavefunctions,Part[ReadList[filename<>"_product_wf.dat",Table["Real",{k,1,1+nS}]],All,1+i]],{i,1,nS}];
(* Read overlap matrix from disk *)
overlapMatrix=ReadList[filename<>"_overlaps.dat",Table["Real",{k,1,nS}]];
(* Read alpha matrix from disk *)
alphaMatrix=ReadList[filename<>"_alphas.dat",Table["Real",{k,1,nS}]];
(* clean the data files in the current directory *)
DeleteFile[FileNames[filename<>"_*.dat"]];
(* Now gather all the results in one nested list *)
{{reactantEnergies,reactantPotential,reactantWavefunctions},{productEnergies,productPotential,productWavefunctions},overlapMatrix,alphaMatrix}
];

(*
The following function calculates energies and wavefunctions of the particle in two diabatic one-dimensional potentials as well as the ovelap integrals and vibronic couplings.
*)

FGHVibronicCouplings[file_,nGridPoints_,nStates_,OptionsPattern[{M->MassH}]]:=
Module[{mass,massString,filename,nG,nS,error,command,reactantEnergies,productEnergies,reactantPotential,productPotential,reactantWavefunctions,productWavefunctions,overlapMatrix,diabCouplingMatrix,productCouplingMatrix,kappaCouplingMatrix,halfTunnelingSplittingMatrix,taueMatrix,taupMatrix,pMatrix,kappaMatrix},
filename=StringSplit[file,"."][[1]];
nG=nGridPoints;
(* Check if nStates < nGridPoints, otherwise set nStates to default value of 10 *)
If[nStates<nGridPoints,nS=nStates,nS=10];
mass=OptionValue[M]*Dalton;
massString=ToString[FortranForm[mass]];
(* run external executable to generate all the data *)
command=KappaCouplingCommand<>" "<>file<>" "<>ToString[nG]<>" "<>massString<>" "<>ToString[nS]<>" > /dev/null";
error=Run[command];
If[error!=0,Message[FGHVibronicCouplings::noexec];Abort[]];
(* Read energies from disk *)
reactantEnergies=ReadList[filename<>"_reactant_en.dat",Real];
productEnergies=ReadList[filename<>"_product_en.dat",Real];
(* Read splined potentials from disk *)
reactantPotential=ReadList[filename<>"_left.dat",Table["Real",{k,1,2+2*nS}]][[All,{1,2}]];
productPotential=Part[ReadList[filename<>"_right.dat",Table["Real",{k,1,2+2*nS}]],All,{1,2}];
(* Read wavefunctions from disk *)
reactantWavefunctions={};
productWavefunctions={};
Do[reactantWavefunctions=Append[reactantWavefunctions,Part[ReadList[filename<>"_left.dat",Table["Real",{k,1,2+2*nS}]],All,2+2*i]],{i,1,nS}];
Do[productWavefunctions=Append[productWavefunctions,Part[ReadList[filename<>"_right.dat",Table["Real",{k,1,2+2*nS}]],All,2+2*i]],{i,1,nS}];
(* Read other matrices from disk files *)
overlapMatrix=ReadList[filename<>"_overlaps.dat",Table["Real",{k,1,nS}]];
diabCouplingMatrix=ReadList[filename<>"_vibcouplings_diab.dat",Table["Real",{k,1,nS}]];
productCouplingMatrix=ReadList[filename<>"_vibcouplings_product.dat",Table["Real",{k,1,nS}]];
kappaCouplingMatrix=ReadList[filename<>"_vibcouplings_kappa.dat",Table["Real",{k,1,nS}]];
halfTunnelingSplittingMatrix=ReadList[filename<>"_tunneling_splittings.dat",Table["Real",{k,1,nS}]]/2;
taueMatrix=ReadList[filename<>"_tau_e.dat",Table["Real",{k,1,nS}]];
taupMatrix=ReadList[filename<>"_tau_p.dat",Table["Real",{k,1,nS}]];
pMatrix=ReadList[filename<>"_p_adiabaticity.dat",Table["Real",{k,1,nS}]];
kappaMatrix=ReadList[filename<>"_kappa.dat",Table["Real",{k,1,nS}]];
(* clean the data files in the current directory *)
DeleteFile[FileNames[filename<>"_*.dat"]];
(* Now gather all the results in one nested list *)
{{reactantEnergies,reactantPotential,reactantWavefunctions},{productEnergies,productPotential,productWavefunctions},overlapMatrix,diabCouplingMatrix,productCouplingMatrix,kappaCouplingMatrix,halfTunnelingSplittingMatrix,taueMatrix,taupMatrix,pMatrix,kappaMatrix}
];

(*
DVR for an arbitrary N-dimensional potential on the grid: energy levels and wavefunctions
*)

GridSolverND[vGrid_, \[Delta]_, n_, OptionsPattern[{M -> MassH}]] := 
  Module[{Vmin,vGrid0,tuples, dims, ndim, \[Mu], T, V, h, en, v},
     \[Mu] = OptionValue[M]*Dalton;
     dims = Dimensions[vGrid];
     Vmin=Min[vGrid];
     vGrid0=vGrid-Vmin;
     ndim = Length[dims];
     tuples = Tuples[Table[Range[1, dims[[idim]]], {idim, 1, ndim}]];
     T = -(1/(2 \[Mu] \[Delta]^2)) AdjacencyMatrix[GridGraph[dims]];
     V = DiagonalMatrix[SparseArray[(vGrid0[[##]] & @@@ tuples) + ndim/(\[Mu] \[Delta]^2)]];
     h = N[T + V];
     en = Reverse@Eigenvalues[h, -n, Method -> {"Arnoldi", "Criteria" -> "RealPart"}];
     v = Reverse@Eigenvectors[h, -n, Method -> {"Arnoldi", "Criteria" -> "RealPart"}];
     {en+Vmin, Table[ArrayReshape[v[[i]], dims], {i, 1, Length[v]}]}
];

(*
DVR for an arbitrary N-dimensional potential on the grid: energy levels and wavefunctions.
General approach with an arbitrary difference order for the numerical second derivative
*)

GenGridSolverND[vGrid_,\[Delta]_,nRoots_,OptionsPattern[{M->MassH,DiffOrder->4}]]:=
   Module[{Vmin,vGrid0,\[Mu],dims,ndim,Ranges,tuples,T,V,h,en,v,nOrder},
      nOrder=OptionValue[DiffOrder];
      \[Mu]=OptionValue[M]*Dalton;
      dims=Dimensions[vGrid];
      Vmin=Min[vGrid];
      vGrid0=vGrid-Vmin;
      ndim=Length[dims];
      Ranges=Table[Range[1,dims[[idim]]],{idim,1,ndim}];
      tuples=Tuples[Ranges];
      T=-(1/(2 \[Mu] \[Delta]^2)) Sum[NDSolve`FiniteDifferenceDerivative[i,Ranges,"DifferenceOrder"->nOrder]["DifferentiationMatrix"],{i,DiagonalMatrix[Table[2,{i,1,ndim}]]}];
      V=DiagonalMatrix[SparseArray[(vGrid0[[##]]&@@@tuples)]];
      h=N[T+V];
      en=Reverse@Re@Eigenvalues[h,-nRoots,Method->{"Arnoldi","Criteria"->"RealPart"}];
      v=Reverse@Re@Eigenvectors[h,-nRoots,Method->{"Arnoldi","Criteria"->"RealPart"}];
      {en+Vmin,Table[ArrayReshape[v[[i]],dims],{i,1,Length[v]}]}
];

(*
DVR direct, without Fourier TRansform
(original code by Max Secor, 2019)
*)
FGHMax[potx_,nGrid_,nRoots_,OptionsPattern[{M->MassH}]]:=Module[
{\[Mu],potSpline,potxSplined,pot,xGrid,wavefunctions,k,xLeft,xRight,dx,KE,PE,Ham,energies,vectors},

\[Mu]=OptionValue[M]*Dalton;
potSpline=Interpolation[potx, Method->"Spline",InterpolationOrder->3];

(* calculate potential on the grid *)
xLeft=Min[potx[[All,1]]];
xRight=Max[potx[[All,1]]];
dx=(xRight-xLeft)/(nGrid-1);
xGrid=Table[x,{x,xLeft,xRight,dx}];
pot=Table[potSpline[x],{x,xLeft,xRight,dx}];
potxSplined=Table[{x,potSpline[x]au2kcal},{x,xLeft,xRight,dx}];

(*Construct kinetic and potential energy matrices*)
k=Pi/(dx a2bohr);
KE=Table[
1/(2\[Mu]) If[
i==j,k^2/3,
((2*k^2)/\[Pi]^2)*(-1)^(i-j)/(i-j)^2
],
{i,0,nGrid-1},{j,0,nGrid-1}];
PE=DiagonalMatrix[pot];

(*Construct the Hamiltionian matrix and diagonalize*)
Ham=N[KE+PE];
energies = Reverse@Eigenvalues[Ham, -nRoots, Method -> {"Arnoldi", "Criteria" -> "RealPart"}];
vectors = Reverse@Eigenvectors[Ham, -nRoots, Method -> {"Arnoldi", "Criteria" -> "RealPart"}];
wavefunctions=Table[Table[{xGrid[[k]],vectors[[i,k]]},{k,1,nGrid}],{i,1,nRoots}];

(*Output*)
{energies*au2kcal,potxSplined,wavefunctions}
];

(*
Averaged Franck-Condon factors
*)

(*
Distribution functions for R-oscillator
*)

PQuantal[T_,R_,MR_,FR_,x_]:=Module[{M,Omega},
   M=MR*Dalton;
   Omega=FR*cm2au;
   Sqrt[(M*Omega/\[Pi])*Tanh[Omega/(2*kb*T)]]*Exp[-M*Omega*Tanh[Omega/(2*kb*T)]*(x-R)^2]
];

PClassical[T_,R_,MR_,FR_,x_]:=Module[{M,Omega},
   M=MR*Dalton;
   Omega=FR*cm2au;
   Sqrt[(M*Omega^2)/(2*\[Pi]*kb*T)]*Exp[-(1/(2*kb*T))*M*Omega^2*(x-R)^2]
];

(*
Franck-Condon Factors
*)

HarmonicFranckCondonAveraged[mu_,nu_,f1_,f2_,d_,T_,MR_,FR_,OptionsPattern[{M->MassH,Distribution->"Classical"}]]:=Module[{m,o1,o2,distr,integrand},
   m=OptionValue[M];
   distr=OptionValue[Distribution];
   o1=f1*cm2au;
   o2=f2*cm2au;
   integrand[x_]:=Switch[distr,
                        "Classical",HarmonicOverlap[mu,nu,f1,f2,x,M->m]^2*PClassical[T,d,MR,FR,x],
                        "Quantal",  HarmonicOverlap[mu,nu,f1,f2,x,M->m]^2*PQuantal[T,d,MR,FR,x],
                        _,          HarmonicOverlap[mu,nu,f1,f2,x,M->m]^2*PClassical[T,d,MR,FR,x]];
   NIntegrate[integrand[X],{X,0,\[Infinity]}]
];

MorseFranckCondonAveraged[mu_,nu_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,MR_,FR_,OptionsPattern[{M->MassH,Overlap->"Numerical",Distribution->"Classical"}]]:=
Module[{m,sw,PR,s,distr,integrand},
   m=OptionValue[M];
   distr=OptionValue[Distribution];
   sw=OptionValue[Overlap];
   PR[x_]:=Switch[distr,
                  "Classical",PClassical[T,d,MR,FR,x],
                  "Quantal",  PQuantal[T,d,MR,FR,x],
                  _,          PClassical[T,d,MR,FR,x]];
   s[x_]:=Switch[sw,
                 "Numerical",     NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,x,M->m],
                 "Analytical",    MorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,x,M->m],
                 "AnalyticalSym", MorseOverlapSym[DE1,Beta1,mu,nu,x,M->m],
                 _,               NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,x,M->m]];
   NIntegrate[Abs[s[x]]^2*PR[x],{x,0,\[Infinity]}]
];


(*
Time correlation functions of the R-oscillator (proton donor-acceptor mode)
*)

(*
Classical oscillator
*)

(*
Variance (mean square displacement): \[Sigma]^2=C (0).
MR - reduced mass of the oscillator (Daltons);
Freq - frequency of the oscillator (1/cm);
T - temperature (K).
Returns variance in Bohr^2.
*)

Sigma2RClassical[MR_,Freq_,T_]:=Module[{mu,f},
	mu=MR*Dalton;
	f=Freq*cm2au;
	kb*T/(mu*f^2)
];

(*
Time correlation function C (t):
MR - reduced mass of the oscillator (Daltons);
Freq - frequency of the oscillator (1/cm);
T - temperature (K);
t - time (atomic units).
Returns correlation function in Bohr^2.
*)

CRClassical[MR_,Freq_,T_,t_]:=Module[{mu,f},
	mu=MR*Dalton;
	f=Freq*cm2au;
	kb*T*Cos[f*t]/(mu*f^2)
];

(*
Quantum oscillator
*)

(*
Variance (mean square displacement) C (0):
MR - reduced mass of the oscillator (Daltons);
Freq - frequency of the oscillator (1/cm);
T - temperature (K).
Returns variance in Bohr^2.
*)

Sigma2RQuantum[MR_,Freq_,T_]:=Module[{mu,f},
	mu=MR*Dalton;
	f=Freq*cm2au;
	hbar*Coth[hbar*f/(2*kb*T)]/(2*mu*f)
];

(*
Time correlation function C (t):
MR - reduced mass of the oscillator (Daltons);
Freq - frequency of the oscillator (1/cm);
T - temperature (K);
t - time (atomic units).
Returns correlation function in Bohr^2.
*)

CRQuantum[MR_,Freq_,T_,t_]:=Module[{mu,f,c1,c2},
	mu=MR*Dalton;
	f=Freq*cm2au;
	c1=hbar^2/(2*mu*f);
	c2=Coth[hbar*f/(2*kb*T)]*Cos[f*t]+I*Sin[f*t];
	c1*c2
];

(*
Coupling reorganization energy (lambda_alpha)
*)

(*
Expression for general mirrored Morse potentials:
MR - reduced mass of the donor-acceptor mode (Dalton);
M - mass of the particle (Dalton);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
i, j - quantum numbers for the left and right oscillators, respectively;
d - distance between the minima of the morse potentials (Bohr).
Returns the coupling reorganization energy in kcal/mol.
*)

NLambdaAlphaMorse[MR_,DE1_,Beta1_,DE2_,Beta2_,i_,j_,d_,OptionsPattern[{M->MassH,Step->0.0001}]]:=Module[{mur,alpha},
	mur=MR*Dalton;
	alpha=NAlphaMorse[DE1,Beta1,DE2,Beta2,i,j,d,M->OptionValue[M],Step->OptionValue[Step]];
	au2kcal*hbar^2*alpha^2/(2*mur)
];

(*
Expression using analytical derivative of analytical overlap integrals (for identical mirrored Morse potentials):
AccN - accuracy
MR - reduced mass of the donor-acceptor mode (Dalton);
M - mass of the particle (Dalton);
DE - dissociation energy for the oscillator (kcal/mol);
Beta - beta parameter for the oscillator on the left (1/Bohr);
nu1, nu2 - quantum numbers for the left and right oscillators, respectively;
d - distance between the minima of the morse potentials (Bohr).
*)

LambdaAlphaMorseSym[MR_,DE_,Beta_,nu1_,nu2_,d_,OptionsPattern[{Acc->50,M->MassH}]]:=Module[{mur,alpha},
	mur=MR*Dalton;
	alpha=AlphaMorseSym[DE,Beta,nu1,nu2,d,Acc->OptionValue[Acc],M->OptionValue[M]];
	au2kcal*hbar^2*alpha^2/(2*mur)
];

(*
Expression using analytical derivative of analytical overlap integrals (for general mirrored Morse potentials):
*)

LambdaAlphaMorse[MR_,DE1_,Beta1_,DE2_,Beta2_,nu1_,nu2_,d_,OptionsPattern[{Acc->50,M->MassH}]]:=Module[{mur,alpha},
	mur=MR*Dalton;
	alpha=AlphaMorse[DE1,Beta1,DE2,Beta2,nu1,nu2,d,Acc->OptionValue[Acc],M->OptionValue[M]];
	au2kcal*hbar^2*alpha^2/(2*mur)
];

(*
Coupling reorganization energy for constant alpha given as input parameter.
MR - reduced mass of the donor-acceptor mode (Daltons),
alpha - alpha parameter (1/Bohr).
Returns coupling reorganization energy in kcal/mol.
*)

LambdaAlphaInput[MR_,alpha_]:=Module[{mur},
	mur=MR*Dalton;
	au2kcal*hbar^2*alpha^2/(2*mur)
];

(*
Nonadiabatic flux probability correlation function
*)

(*
Flux components
*)

(*
Classical solvent damping term (Gaussian damping for high temperature and short time):
lambda - solvent reorganization energy (kcal/mol);
T - temperature (K);
t - time (atomic units).
Return solvent damping Gaussian (unitless).
*)

SolventDampingClassical[lambda_,T_,t_]:=Module[{lau},
	lau=lambda*kcal2au;
	Exp[-lau*kb*T*t^2/hbar^2]
];

(*
Coherent (oscillatory) term:
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
M - mass of the particle (Dalton);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
i, j - quantum numbers for the left and right oscillators, respectively;
d - distance between the minima of the morse potentials (Bohr);
t - time (atomic units).
*)

CoherentTermMorse[DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,nu1_,nu2_,t_,OptionsPattern[{M->MassH,Alpha->"Numerical"}]]:=Module[{sw,dg,la,f,e2,e20,e1,e10,lalpha},
	sw=OptionValue[Alpha];
	dg=DG*kcal2au;
	la=lambda*kcal2au;
	f=Freq*cm2au;
	e10=MorseEnergy[0,Beta1,DE1,M->OptionValue[M]]*kcal2au;
	e1=MorseEnergy[nu1,Beta1,DE1,M->OptionValue[M]]*kcal2au;
	e20=MorseEnergy[0,Beta2,DE2,M->OptionValue[M]]*kcal2au;
	e2=MorseEnergy[nu2,Beta2,DE2,M->OptionValue[M]]*kcal2au;
	lalpha=Switch[Head[sw],
		Integer,kcal2au*LambdaAlphaInput[MR,sw],
		Real,kcal2au*LambdaAlphaInput[MR,sw],
		Rational,kcal2au*LambdaAlphaInput[MR,sw],
		String,kcal2au*Switch[sw,
			"Numerical",NLambdaAlphaMorse[MR,DE1,Beta1,DE2,Beta2,nu1,nu2,d,M->OptionValue[M]],
			"Analytical",LambdaAlphaMorse[MR,DE1,Beta1,DE2,Beta2,nu1,nu2,d,M->OptionValue[M]],
			"AnalyticalSym",LambdaAlphaMorseSym[MR,DE1,Beta1,nu1,nu2,d,M->OptionValue[M]],
			_,NLambdaAlphaMorse[MR,DE1,Beta1,DE2,Beta2,nu1,nu2,d,M->OptionValue[M]]],
		_,kcal2au*NLambdaAlphaMorse[MR,DE1,Beta1,DE2,Beta2,nu1,nu2,d,M->OptionValue[M]]
	];
	Cos[(dg+la+e2-e20-e1+e10)*t+lalpha/(f*hbar)*Sin[f*t]]
];

(*
Coherent (oscillatory) term:
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
i, j - quantum numbers for the left and right oscillators, respectively;
t - time (atomic units).
*)

CoherentTermFGH[fghList_,DG_,lambda_,MR_,Freq_,mu_,nu_,t_,OptionsPattern[{Alpha->"Numerical"}]]:=
Module[{sa,dg,la,f,e2,e20,e1,e10,lalpha,Ns},
        Ns = Length[fghList[[1,1]]];
        If[mu+1>Ns||nu+1>Ns,Message[CoherentTermFGH::numberofstates];Abort[]];
	sa=OptionValue[Alpha];
	dg=DG*kcal2au;
	la=lambda*kcal2au;
	f=Freq*cm2au;
	e10=fghList[[1,1,1]]*kcal2au;
	e1=fghList[[1,1,mu+1]]*kcal2au;
	e20=fghList[[2,1,1]]*kcal2au;
	e2=fghList[[2,1,nu+1]]*kcal2au;
	lalpha=Switch[Head[sa],
		Integer,kcal2au*LambdaAlphaInput[MR,sa],
		Real,kcal2au*LambdaAlphaInput[MR,sa],
		Rational,kcal2au*LambdaAlphaInput[MR,sa],
		String,Switch[sa,
			"Numerical",(fghList[[4]][[mu+1,nu+1]]/a2bohr)^2/(2*MR*Dalton),
			          _,(fghList[[4]][[mu+1,nu+1]]/a2bohr)^2/(2*MR*Dalton)],
		_,(fghList[[4]][[mu+1,nu+1]]/a2bohr)^2/(2*MR*Dalton)
	];
	Cos[(dg+la+e2-e20-e1+e10)*t+lalpha/(f*hbar)*Sin[f*t]]
];

(*
Term originating from the donor-acceptor vibrational mode:
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
M - mass of the particle (Dalton);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - distance between the minima of the morse potentials (Bohr);
T - temperature (K);
nu1, nu2 - quantum numbers for the left and right oscillators, respectively;
t - time (atomic units).

*)

RmodeTermMorse[MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,nu1_,nu2_,t_,OptionsPattern[{M->MassH,Alpha->"Numerical"}]]:=Module[{sw,f,lalpha},
	sw=OptionValue[Alpha];
	f=Freq*cm2au;
	lalpha=Switch[Head[sw],
		Integer,kcal2au*LambdaAlphaInput[MR,sw],
		Real,kcal2au*LambdaAlphaInput[MR,sw],
		Rational,kcal2au*LambdaAlphaInput[MR,sw],
		String,kcal2au*Switch[sw,
			"Numerical",NLambdaAlphaMorse[MR,DE1,Beta1,DE2,Beta2,nu1,nu2,d,M->OptionValue[M]],
			"Analytical",LambdaAlphaMorse[MR,DE1,Beta1,DE2,Beta2,nu1,nu2,d,M->OptionValue[M]],
			"AnalyticalSym",LambdaAlphaMorseSym[MR,DE1,Beta1,nu1,nu2,d,M->OptionValue[M]],
			_,NLambdaAlphaMorse[MR,DE1,Beta1,DE2,Beta2,nu1,nu2,d,M->OptionValue[M]]],
		_,kcal2au*NLambdaAlphaMorse[MR,DE1,Beta1,DE2,Beta2,nu1,nu2,d,M->OptionValue[M]]
	];
	Exp[lalpha/(hbar*f)*Coth[(hbar*f)/(2*kb*T)]*(Cos[f*t]+1)]
];

(*
Term originating from the donor-acceptor vibrational mode:
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
T - temperature (K);
nu1, nu2 - quantum numbers for the left and right oscillators, respectively;
t - time (atomic units).

*)

RmodeTermFGH[fghList_,MR_,Freq_,T_,mu_,nu_,t_,OptionsPattern[{Alpha->"Numerical"}]]:=
Module[{sa,f,lalpha,Ns},
        Ns = Length[fghList[[1,1]]];
        If[mu+1>Ns||nu+1>Ns,Message[RmodeTermFGH::numberofstates];Abort[]];
	sa=OptionValue[Alpha];
	f=Freq*cm2au;
	lalpha=Switch[Head[sa],
		Integer,kcal2au*LambdaAlphaInput[MR,sa],
		Real,kcal2au*LambdaAlphaInput[MR,sa],
		Rational,kcal2au*LambdaAlphaInput[MR,sa],
		String,Switch[sa,
			"Numerical",(fghList[[4]][[mu+1,nu+1]]/a2bohr)^2/(2*MR*Dalton),
			          _,(fghList[[4]][[mu+1,nu+1]]/a2bohr)^2/(2*MR*Dalton)],
		_,(fghList[[4]][[mu+1,nu+1]]/a2bohr)^2/(2*MR*Dalton)
	];
	Exp[lalpha/(hbar*f)*Coth[(hbar*f)/(2*kb*T)]*(Cos[f*t]+1)]
];

(*
An finally, the real part of the total nonadiabatic flux correlation function:
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
M - mass of the particle (Dalton);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
nu1, nu2 - quantum numbers for the left and right oscillators, respectively;
t - time (atomic units).
*)

RealFluxQuantumMorse[DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,nu1_,nu2_,t_,OptionsPattern[{M->MassH,Alpha->"Numerical"}]]:=Module[{sterm,cterm,rterm},
	sterm=SolventDampingClassical[lambda,T,t];
	cterm=CoherentTermMorse[DG,lambda,MR,Freq,DE1,Beta1,DE2,Beta2,d,nu1,nu2,t,M->OptionValue[M],Alpha->OptionValue[Alpha]];
	rterm=RmodeTermMorse[MR,Freq,DE1,Beta1,DE2,Beta2,d,T,nu1,nu2,t,M->OptionValue[M],Alpha->OptionValue[Alpha]];
	sterm*cterm*rterm
];

(*
The real part of the total nonadiabatic flux correlation function (general potentials):
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
T - temperature (K);
nu1, nu2 - quantum numbers for the left and right oscillators, respectively;
t - time (atomic units).
*)

RealFluxQuantumFGH[fghList_,DG_,lambda_,MR_,Freq_,T_,mu_,nu_,t_,OptionsPattern[{Alpha->"Numerical"}]]:=Module[{sterm,cterm,rterm},
	sterm=SolventDampingClassical[lambda,T,t];
	cterm=CoherentTermFGH[fghList,DG,lambda,MR,Freq,mu,nu,t,Alpha->OptionValue[Alpha]];
	rterm=RmodeTermFGH[fghList,MR,Freq,T,mu,nu,t,Alpha->OptionValue[Alpha]];
	sterm*cterm*rterm
];

(*
Exact nonadiabatic rate constant and kinetic isotope effect (KIE):
numerical integration of the quantum flux expression
*)

(*
Rate constant
*)

(*
Partial nonadiabatic first order rate constant k (mu->nu) (1/sec):
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
M - mass of the particle (Dalton);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mu, nu - quantum numbers for the left and right oscillators, respectively.
*)

PartialRateQfluxMorse[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mu_,nu_,OptionsPattern[{M->MassH,Overlap->"Numerical",Alpha->"Numerical"}]]:=Module[{sw,sa,time,rterm0,reflux,s,s2,prefactor,v2,fluxint},
	sw=OptionValue[Overlap];
	sa=OptionValue[Alpha];
	rterm0=RmodeTermMorse[MR,Freq,DE1,Beta1,DE2,Beta2,d,T,mu,nu,0,M->OptionValue[M],Alpha->sa];
	reflux[t_]:=SetAccuracy[RealFluxQuantumMorse[DG,lambda,MR,Freq,DE1,Beta1,DE2,Beta2,d,T,mu,nu,t,M->OptionValue[M],Alpha->sa]/rterm0,30];
	s=Switch[sw,
		"Numerical",NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
		"Analytical",MorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
		"AnalyticalSym",MorseOverlapSym[DE1,Beta1,mu,nu,d,M->OptionValue[M]],
		_,NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]];
	s2=s^2;
	prefactor=10^12/au2ps;
	v2=(V*kcal2au)^2;
	fluxint=NIntegrate[reflux[time],{time,0,\[Infinity]},Method->{"DoubleExponential","SymbolicProcessing"->0},WorkingPrecision->$MachinePrecision];
	2*prefactor*v2*s2*rterm0*fluxint
];

(*
Partial nonadiabatic first order rate constant k (mu->nu) (1/sec) (general potentials):
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
T - temperature (K);
mu, nu - quantum numbers for the left and right oscillators, respectively.
*)

PartialRateQfluxFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mu_,nu_,OptionsPattern[{Alpha->"Numerical"}]]:=
Module[{sa,time,rterm0,reflux,s,s2,prefactor,v2,fluxint,Ns},
        Ns = Length[fghList[[1,1]]];
        If[mu+1>Ns||nu+1>Ns,Message[PartialRateQfluxFGH::numberofstates];Abort[]];
	sa=OptionValue[Alpha];
	rterm0=RmodeTermFGH[fghList,MR,Freq,T,mu,nu,0,Alpha->sa];
	reflux[t_]:=SetAccuracy[RealFluxQuantumFGH[fghList,DG,lambda,MR,Freq,T,mu,nu,t,Alpha->sa]/rterm0,30];
	s=fghList[[3]][[mu+1,nu+1]];
	s2=s^2;
	prefactor=10^12/au2ps;
	v2=(V*kcal2au)^2;
	fluxint=NIntegrate[reflux[time],{time,0,\[Infinity]},Method->{"DoubleExponential","SymbolicProcessing"->0},WorkingPrecision->$MachinePrecision];
	2*prefactor*v2*s2*rterm0*fluxint
];

(*
Total nonadiabatic first order rate constant k (1/sec):
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
M - mass of the particle (Dalton);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TotalRateQfluxMorse[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_,OptionsPattern[{M->MassH,Overlap->"Numerical",Alpha->"Numerical"}]]:=
Module[{Emu,Z,Pmu,mu,nu},
	Array[Emu,mumax+1,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*MorseEnergy[mu,Beta1,DE1,M->OptionValue[M]],{mu,0,mumax}];
	Z = Exp[Emu[0]/(kb*T)]*QMorseExact[Beta1,DE1,T,M->OptionValue[M]];
	Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
	Sum[Pmu[mu]*Sum[PartialRateQfluxMorse[V,DG,lambda,MR,Freq,DE1,Beta1,DE2,Beta2,d,T,mu,nu,M->OptionValue[M],Overlap->OptionValue[Overlap],Alpha->OptionValue[Alpha]],{nu,0,numax}],{mu,0,mumax}]
];

(*
Total nonadiabatic first order rate constant k (1/sec):
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TotalRateQfluxFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_,OptionsPattern[{Alpha->"Numerical"}]]:=
Module[{Emu,Z,Pmu,mu,nu,Ns},
        Ns = Length[fghList[[1,1]]];
        If[mumax+1>Ns||numax+1>Ns,Message[TotalRateQfluxFGH::numberofstates];Abort[]];
	Array[Emu,Ns,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*fghList[[1,1,mu+1]],{mu,0,Ns-1}];
	Z=Sum[Exp[-Emu[mu]/(kb*T)],{mu,0,Ns-1}];
	Do[Pmu[mu]=Exp[-Emu[mu]/(kb*T)]/Z,{mu,0,mumax}];
	Sum[Pmu[mu]*Sum[PartialRateQfluxFGH[fghList,V,DG,lambda,MR,Freq,T,mu,nu,Alpha->OptionValue[Alpha]],{nu,0,numax}],{mu,0,mumax}]
];

(*
Kinetic isotope effect
*)

(*
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

KIEQfluxMorse[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_,OptionsPattern[{M1->MassH,M2->MassD,Overlap->"Numerical",Alpha->"Numerical"}]]:=Module[{rate1,rate2},
	rate1=TotalRateQfluxMorse[V,DG,lambda,MR,Freq,DE1,Beta1,DE2,Beta2,d,T,mumax,numax,M->OptionValue[M1],Overlap->OptionValue[Overlap],Alpha->OptionValue[Alpha]];
	rate2=TotalRateQfluxMorse[V,DG,lambda,MR,Freq,DE1,Beta1,DE2,Beta2,d,T,mumax,numax,M->OptionValue[M2],Overlap->OptionValue[Overlap],Alpha->OptionValue[Alpha]];
	rate1/rate2
];

(*
fghList1 - nested list with energies, wavefunctions, overlap integrals and alphas for isotope 1;
fghList2 - nested list with energies, wavefunctions, overlap integrals and alphas for isotope 2;
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

KIEQfluxFGH[fghList1_,fghList2_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_,OptionsPattern[{Alpha1->"Numerical",Alpha2->"Numerical"}]]:=Module[{rate1,rate2},
	rate1=TotalRateQfluxFGH[fghList1,V,DG,lambda,MR,Freq,T,mumax,numax,Alpha->OptionValue[Alpha1]];
	rate2=TotalRateQfluxFGH[fghList2,V,DG,lambda,MR,Freq,T,mumax,numax,Alpha->OptionValue[Alpha2]];
	rate1/rate2
];

(*
Table with rate channel contributions
*)

(*
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor - acceptor mode (Dalton);
Freq - donor - acceptor mode frequency (1/cm);
M - mass of the particle (Dalton);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TableRateQfluxMorse[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_,OptionsPattern[{M->MassH,Overlap->"Numerical",Alpha->"Numerical"}]]:=Module[{headings,dg0,la,totalrate,sw,sa,s,s2,alpha,dg,Ea,w,Emu,Enu,Pmu,Z,kmunu,info,mu,nu,Mass},
	Mass=OptionValue[M];
	sw=OptionValue[Overlap];
	sa=OptionValue[Alpha];
	headings={"\[Mu]","\[Nu]",strpmu,strdgmunu,strdgddmunu,strs2munu,stramunu,"exp(-\[Beta]"<>strdgddmunu<>")","%"};
	info="Isotope: "<>ToString[NumberForm[Mass,3]]<>"amu; "<>strdg00<>"="<>ToString[NumberForm[DG,4]]<>"kcal/mol; \[Lambda]="<>ToString[NumberForm[lambda,4]]<>"kcal/mol; M="<>ToString[NumberForm[MR,3]]<>"amu; \[CapitalOmega]="<>ToString[NumberForm[Freq,4]]<>strcm;
	dg0=DG*kcal2au;
	la=lambda*kcal2au;
	totalrate=TotalRateQfluxMorse[V,DG,lambda,MR,Freq,DE1,Beta1,DE2,Beta2,d,T,mumax,numax,M->OptionValue[M],Overlap->OptionValue[Overlap],Alpha->OptionValue[Alpha]];
	Array[s2,{mumax+1,numax+1},{0,0}];
	Array[alpha,{mumax+1,numax+1},{0,0}];
	Array[dg,{mumax+1,numax+1},{0,0}];
	Array[Ea,{mumax+1,numax+1},{0,0}];
	Array[w,{mumax+1,numax+1},{0,0}];
	Array[Emu,mumax+1,0];
	Array[Enu,numax+1,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*MorseEnergy[mu,Beta1,DE1,M->OptionValue[M]],{mu,0,mumax}];
	Do[Enu[nu]=kcal2au*MorseEnergy[nu,Beta2,DE2,M->OptionValue[M]],{nu,0,numax}];
	Z = Exp[Emu[0]/(kb*T)]*QMorseExact[Beta1,DE1,T,M->OptionValue[M]];
	Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
	Do[
		s=Switch[sw,
			"Numerical",NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
			"Analytical",MorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
			"AnalyticalSym",MorseOverlapSym[DE1,Beta1,mu,nu,d,M->OptionValue[M]],
			_,NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]];
		s2[mu,nu]=s^2;
	        alpha[mu,nu]=Switch[Head[sa],
		        Integer,sa/bohr2a,
		        Real,sa/bohr2a,
		        Rational,sa/bohr2a,
		        String,Switch[sa,
			        "Numerical",NAlphaMorse[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]/bohr2a,
			        "Analytical",AlphaMorse[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]/bohr2a,
			        "AnalyticalSym",AlphaMorseSym[DE1,Beta1,mu,nu,d,M->OptionValue[M]]/bohr2a,
			        _,NAlphaMorse[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]/bohr2a
			        ],
		        _,NAlphaMorse[DE1,Beta1,DE2,Beta2,nu1,nu2,d,M->OptionValue[M]]
	                ];
		dg[mu,nu]=dg0+Enu[nu]-Enu[0]-(Emu[mu]-Emu[0]);
		Ea[mu,nu]=(dg[mu,nu]+la)^2/(4*la);
		kmunu=PartialRateQfluxMorse[V,DG,lambda,MR,Freq,DE1,Beta1,DE2,Beta2,d,T,mu,nu,M->OptionValue[M],Overlap->OptionValue[Overlap],Alpha->OptionValue[Alpha]];
		w[mu,nu]=Pmu[mu]*kmunu*100/totalrate,{mu,0,mumax},{nu,0,numax}
	];
	Grid[Join[{{info,SpanFromLeft}},{headings},Flatten[Table[{mu,nu,ScientificForm[Pmu[mu],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[dg[mu,nu]*au2kcal,{6,3}],PaddedForm[Ea[mu,nu]*au2kcal,{6,3}],ScientificForm[s2[mu,nu],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[alpha[mu,nu],{6,3}],ScientificForm[Exp[-Ea[mu,nu]/(kb*T)],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[w[mu,nu],{7,3}]},{mu,0,mumax},{nu,0,numax}],1]],Frame->All,Spacings->{Automatic,1}]
];


(*
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor - acceptor mode (Dalton);
Freq - donor - acceptor mode frequency (1/cm);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TableRateQfluxFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_,OptionsPattern[{Alpha->"Numerical"}]]:=
Module[{headings,dg0,la,totalrate,sa,s,s2,alpha,dg,Ea,w,Emu,Enu,Pmu,Z,kmunu,info,mu,nu,Ns},
        Ns = Length[fghList[[1,1]]];
        If[mumax+1>Ns||numax+1>Ns,Message[TableRateQfluxFGH::numberofstates];Abort[]];
	sa=OptionValue[Alpha];
	headings={"\[Mu]","\[Nu]",strpmu,strdgmunu,strdgddmunu,strs2munu,stramunu,"exp(-\[Beta]"<>strdgddmunu<>")","%"};
	info=strdg00<>"="<>ToString[NumberForm[DG,4]]<>"kcal/mol; \[Lambda]="<>ToString[NumberForm[lambda,4]]<>"kcal/mol; M="<>ToString[NumberForm[MR,3]]<>"amu; \[CapitalOmega]="<>ToString[NumberForm[Freq,4]]<>strcm;
	dg0=DG*kcal2au;
	la=lambda*kcal2au;
	totalrate=TotalRateQfluxFGH[fghList,V,DG,lambda,MR,Freq,T,mumax,numax,Alpha->OptionValue[Alpha]];
	Array[s2,{mumax+1,numax+1},{0,0}];
	Array[alpha,{mumax+1,numax+1},{0,0}];
	Array[dg,{mumax+1,numax+1},{0,0}];
	Array[Ea,{mumax+1,numax+1},{0,0}];
	Array[w,{mumax+1,numax+1},{0,0}];
	Array[Emu,Ns,0];
	Array[Enu,Ns,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*fghList[[1,1,mu+1]],{mu,0,Ns-1}];
	Do[Enu[nu]=kcal2au*fghList[[2,1,nu+1]],{nu,0,Ns-1}];
	Z=Sum[Exp[-(Emu[mu]-Emu[0])/(kb*T)],{mu,1,Ns-1}];
	Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
	Do[
		s=fghList[[3]][[mu+1,nu+1]];
		s2[mu,nu]=s^2;
	        alpha[mu,nu]=Switch[Head[sa],
		        Integer,sa/bohr2a,
		        Real,sa/bohr2a,
		        Rational,sa/bohr2a,
		        String,Switch[sa,
			        "Numerical",fghList[[4]][[mu+1,nu+1]],
			        _,fghList[[4]][[mu+1,nu+1]]
			        ],
		        _,fghList[[4]][[mu+1,nu+1]]
	                ];
		dg[mu,nu]=dg0+Enu[nu]-Enu[0]-(Emu[mu]-mu[0]);
		Ea[mu,nu]=(dg[mu,nu]+la)^2/(4*la);
		kmunu=PartialRateQfluxFGH[fghList,V,DG,lambda,MR,Freq,T,mu,nu,Alpha->OptionValue[Alpha]];
		w[mu,nu]=Pmu[mu]*kmunu*100/totalrate,{mu,0,mumax},{nu,0,numax}
	];
	Grid[Join[{{info,SpanFromLeft}},{headings},Flatten[Table[{mu,nu,ScientificForm[Pmu[mu],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[dg[mu,nu]*au2kcal,{6,3}],PaddedForm[Ea[mu,nu]*au2kcal,{6,3}],ScientificForm[s2[mu,nu],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[alpha[mu,nu],{6,3}],ScientificForm[Exp[-Ea[mu,nu]/(kb*T)],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[w[mu,nu],{7,3}]},{mu,0,mumax},{nu,0,numax}],1]],Frame->All,Spacings->{Automatic,1}]
];

(*
Nonadiabatic rate constant and KIE for fixed proton donor-acceptor distance
*)

(*
Rate constant
*)

(*
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
M - mass of the particle (Dalton);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mu, nu - quantum numbers for the left and right oscillators, respectively.
Return partial rate constant in sec^-1.
*)

PartialRateRfixedMorse[V_,DG_,lambda_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mu_,nu_,OptionsPattern[{M->MassH,Overlap->"Numerical"}]]:=Module[{sw,s,s2,prefactor,v2,la,dg,sq,e10,e1,e20,e2,Ea},
	sw=OptionValue[Overlap];
	s=Switch[sw,
		"Numerical",NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
		"Analytical",MorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
		"AnalyticalSym",MorseOverlapSym[DE1,Beta1,mu,nu,d,M->OptionValue[M]],
		_,NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]];
	s2=s^2;
	prefactor=10^12/au2ps;
	v2=(V*kcal2au)^2;
	la=lambda*kcal2au;
	dg=DG*kcal2au;
	sq=Sqrt[\[Pi]/(kb*T*la)];
	e10=kcal2au*MorseEnergy[0,Beta1,DE1,M->OptionValue[M]];
	e1=kcal2au*MorseEnergy[mu,Beta1,DE1,M->OptionValue[M]];
	e20=kcal2au*MorseEnergy[0,Beta2,DE2,M->OptionValue[M]];
	e2=kcal2au*MorseEnergy[nu,Beta2,DE2,M->OptionValue[M]];
	Ea=(dg+la+e2-e20-e1+e10)^2/(4*la);
	prefactor*v2*s2*sq*Exp[-Ea/(kb*T)]
];


(*
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
T - temperature (K);
mu, nu - quantum numbers for the left and right oscillators, respectively.
Return partial rate constant in sec^-1.
*)

PartialRateRfixedFGH[fghList_,V_,DG_,lambda_,T_,mu_,nu_]:=
Module[{s,s2,prefactor,v2,la,dg,sq,e10,e1,e20,e2,Ea,Ns},
    Ns = Length[fghList[[1,1]]];
    If[mu+1>Ns||nu+1>Ns,Message[PartialRateRfixedFGH::numberofstates];Abort[]];
    s=fghList[[3]][[mu+1,nu+1]];
    s2=s^2;
    prefactor=10^12/au2ps;
    v2=(V*kcal2au)^2;
    la=lambda*kcal2au;
    dg=DG*kcal2au;
    sq=Sqrt[\[Pi]/(kb*T*la)];
    e10=kcal2au*fghList[[1,1,1]];
    e1=kcal2au*fghList[[1,1,mu+1]];
    e20=kcal2au*fghList[[2,1,1]];
    e2=kcal2au*fghList[[2,1,nu+1]];
    Ea=(dg+la+e2-e20-e1+e10)^2/(4*la);
    prefactor*v2*s2*sq*Exp[-Ea/(kb*T)]
];


(*
Total nonadiabatic first order rate constant k (1/sec):
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
M - mass of the particle (Dalton);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TotalRateRfixedMorse[V_,DG_,lambda_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_,OptionsPattern[{M->MassH,Overlap->"Numerical"}]]:=
Module[{Emu,Z,Pmu,mu,nu},
	Array[Emu,mumax+1,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*MorseEnergy[mu,Beta1,DE1,M->OptionValue[M]],{mu,0,mumax}];
	Z = Exp[Emu[0]/(kb*T)]*QMorseExact[Beta1,DE1,T,M->OptionValue[M]];
	Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
	Sum[Pmu[mu]*Sum[PartialRateRfixedMorse[V,DG,lambda,DE1,Beta1,DE2,Beta2,d,T,mu,nu,M->OptionValue[M],Overlap->OptionValue[Overlap]],{nu,0,numax}],{mu,0,mumax}]
];

(*
Total nonadiabatic first order rate constant k (1/sec):
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TotalRateRfixedFGH[fghList_,V_,DG_,lambda_,T_,mumax_,numax_]:=
Module[{Emu,Z,Pmu,mu,nu,Ns},
        Ns = Length[fghList[[1,1]]];
        If[mumax+1>Ns||numax+1>Ns,Message[TotalRateRfixedFGH::numberofstates];Abort[]];
	Array[Emu,Ns,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*fghList[[1,1,mu+1]],{mu,0,Ns-1}];
	Z=Sum[Exp[-(Emu[mu]-Emu[0])/(kb*T)],{mu,0,Ns-1}];
	Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
	Sum[Pmu[mu]*Sum[PartialRateRfixedFGH[fghList,V,DG,lambda,T,mu,nu],{nu,0,numax}],{mu,0,mumax}]
];

(*
Kinetic isotope effect
*)

(*
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

KIERfixedMorse[V_,DG_,lambda_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_,OptionsPattern[{M1->MassH,M2->MassD,Overlap->"Numerical"}]]:=Module[{rate1,rate2},
	rate1=TotalRateRfixedMorse[V,DG,lambda,DE1,Beta1,DE2,Beta2,d,T,mumax,numax,M->OptionValue[M1],Overlap->OptionValue[Overlap]];
	rate2=TotalRateRfixedMorse[V,DG,lambda,DE1,Beta1,DE2,Beta2,d,T,mumax,numax,M->OptionValue[M2],Overlap->OptionValue[Overlap]];
	rate1/rate2
];


(*
fghList1 - nested list with energies, wavefunctions, overlap integrals and alphas for isotope 1;
fghList2 - nested list with energies, wavefunctions, overlap integrals and alphas for isotope 2;
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

KIERfixedFGH[fghList1_,fghList2_,V_,DG_,lambda_,T_,mumax_,numax_]:=Module[{rate1,rate2},
	rate1=TotalRateRfixedFGH[fghList1,V,DG,lambda,T,mumax,numax];
	rate2=TotalRateRfixedFGH[fghList2,V,DG,lambda,T,mumax,numax];
	rate1/rate2
];

(*
Table with rate channel contributions
*)

(*
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
M - mass of the particle (Dalton);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TableRateRfixedMorse[V_,DG_,lambda_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_,OptionsPattern[{M->MassH,Overlap->"Numerical"}]]:=
Module[{headings,dg0,la,totalrate,sw,s,s2,dg,Ea,w,Emu,Enu,Pmu,Z,kmunu,info,Mass},
	Mass=OptionValue[M];
	sw=OptionValue[Overlap];
	dg0=DG*kcal2au;
	headings={"\[Mu]","\[Nu]",strpmu,strdgmunu,strdgddmunu,strs2munu,"exp(-\[Beta]"<>strdgddmunu<>")","%"};
	info="Isotope: "<>ToString[NumberForm[Mass,3]]<>"amu; "<>strdg00<>"="<>ToString[NumberForm[DG,4]]<>"kcal/mol; \[Lambda]="<>ToString[NumberForm[lambda,4]]<>"kcal/mol";
	la=lambda*kcal2au;
	totalrate=TotalRateRfixedMorse[V,DG,lambda,DE1,Beta1,DE2,Beta2,d,T,mumax,numax,M->OptionValue[M],Overlap->OptionValue[Overlap]];
	Array[s2,{mumax+1,numax+1},{0,0}];
	Array[dg,{mumax+1,numax+1},{0,0}];
	Array[Ea,{mumax+1,numax+1},{0,0}];
	Array[w,{mumax+1,numax+1},{0,0}];
	Array[Emu,mumax+1,0];
	Array[Enu,numax+1,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*MorseEnergy[mu,Beta1,DE1,M->OptionValue[M]],{mu,0,mumax}];
	Do[Enu[nu]=kcal2au*MorseEnergy[nu,Beta2,DE2,M->OptionValue[M]],{nu,0,numax}];
	Z = Exp[Emu[0]/(kb*T)]*QMorseExact[Beta1,DE1,T,M->OptionValue[M]];
	Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
	Do[
		s=Switch[sw,
			"Numerical",NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
			"Analytical",MorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
			"AnalyticalSym",MorseOverlapSym[DE1,Beta1,mu,nu,d,M->OptionValue[M]],
			_,NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]];
		s2[mu,nu]=s^2;
		dg[mu,nu]=dg0+Enu[nu]-Enu[0]-(Emu[mu]-Emu[0]);
		Ea[mu,nu]=(dg[mu,nu]+la)^2/(4*la);
		kmunu=PartialRateRfixedMorse[V,DG,lambda,DE1,Beta1,DE2,Beta2,d,T,mu,nu,M->OptionValue[M],Overlap->OptionValue[Overlap]];
		w[mu,nu]=Pmu[mu]*kmunu*100/totalrate,{mu,0,mumax},{nu,0,numax}
	];
	Grid[Join[{{info,SpanFromLeft}},{headings},Flatten[Table[{mu,nu,ScientificForm[Pmu[mu],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[dg[mu,nu]*au2kcal,{6,3}],PaddedForm[Ea[mu,nu]*au2kcal,{6,3}],ScientificForm[s2[mu,nu],6,NumberFormat->(Row[{#1,"e",#3}]&)],ScientificForm[Exp[-Ea[mu,nu]/(kb*T)],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[w[mu,nu],{7,3}]},{mu,0,mumax},{nu,0,numax}],1]],Frame->All,Spacings->{Automatic,1}]
];

(*
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TableRateRfixedFGH[fghList_,V_,DG_,lambda_,T_,mumax_,numax_]:=
Module[{headings,dg0,la,totalrate,s,s2,dg,Ea,w,Emu,Enu,Pmu,Z,kmunu,info,Ns},
        Ns = Length[fghList[[1,1]]];
        If[mumax+1>Ns||numax+1>Ns,Message[TableRateRfixedFGH::numberofstates];Abort[]];
	dg0=DG*kcal2au;
	headings={"\[Mu]","\[Nu]",strpmu,strdgmunu,strdgddmunu,strs2munu,"exp(-\[Beta]"<>strdgddmunu<>")","%"};
	info=strdg00<>"="<>ToString[NumberForm[DG,4]]<>"kcal/mol; \[Lambda]="<>ToString[NumberForm[lambda,4]]<>"kcal/mol";
	la=lambda*kcal2au;
	totalrate=TotalRateRfixedFGH[fghList,V,DG,lambda,T,mumax,numax];
	Array[s2,{mumax+1,numax+1},{0,0}];
	Array[dg,{mumax+1,numax+1},{0,0}];
	Array[Ea,{mumax+1,numax+1},{0,0}];
	Array[w,{mumax+1,numax+1},{0,0}];
	Array[Emu,Ns,0];
	Array[Enu,Ns,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*fghList[[1,1,mu+1]],{mu,0,Ns-1}];
	Do[Enu[nu]=kcal2au*fghList[[2,1,nu+1]],{nu,0,Ns-1}];
	Z=Sum[Exp[-(Emu[mu]-Emu[0])/(kb*T)],{mu,0,Ns-1}];
	Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
	Do[
	    s=fghList[[3]][[mu+1,nu+1]];
	    s2[mu,nu]=s^2;
	    dg[mu,nu]=dg0+Enu[nu]-Enu[0]-(Emu[mu]-Emu[0]);
	    Ea[mu,nu]=(dg[mu,nu]+la)^2/(4*la);
	    kmunu=PartialRateRfixedFGH[fghList,V,DG,lambda,T,mu,nu];
	    w[mu,nu]=Pmu[mu]*kmunu*100/totalrate,{mu,0,mumax},{nu,0,numax}
	];
	Grid[Join[{{info,SpanFromLeft}},{headings},Flatten[Table[{mu,nu,ScientificForm[Pmu[mu],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[dg[mu,nu]*au2kcal,{6,3}],PaddedForm[Ea[mu,nu]*au2kcal,{6,3}],ScientificForm[s2[mu,nu],6,NumberFormat->(Row[{#1,"e",#3}]&)],ScientificForm[Exp[-Ea[mu,nu]/(kb*T)],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[w[mu,nu],{7,3}]},{mu,0,mumax},{nu,0,numax}],1]],Frame->All,Spacings->{Automatic,1}]
];





(*
Nonadiabatic rate constant and KIE with averaged Franck-Condon factors (Ulstrup-Kuznetsov model)
*)

(*
Rate constant for harmonic proton potentials
*)

(*
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
M - mass of the particle (Dalton);
f1 - frequency of the left (proton donor) oscillator;
f2 - frequency of the right (proton acceptor) oscillator;
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mu, nu - quantum numbers for the left and right oscillators, respectively.
Return partial rate constant in sec^-1.
*)

PartialRateUKHarmonic[V_,DG_,lambda_,MR_,FR_,f1_,f2_,d_,T_,mu_,nu_,OptionsPattern[{M->MassH,Distribution->"Classical"}]]:=Module[{m,distr,prefactor,v2,la,fc,dg,sq,e10,e1,e20,e2,Ea},
   m=OptionValue[M];
   distr=OptionValue[Distribution];
   fc=HarmonicFranckCondonAveraged[mu,nu,f1,f2,d,T,MR,FR,M->m,Distribution->distr];
	prefactor=10^12/au2ps;
	v2=(V*kcal2au)^2;
	la=lambda*kcal2au;
	dg=DG*kcal2au;
	sq=Sqrt[\[Pi]/(kb*T*la)];
	e10=kcal2au*HOEnergy[0,f1];
	e1=kcal2au*HOEnergy[mu,f1];
	e20=kcal2au*HOEnergy[0,f2];
	e2=kcal2au*HOEnergy[nu,f2];
	Ea=(dg+la+e2-e20-e1+e10)^2/(4*la);
	prefactor*v2*fc*sq*Exp[-Ea/(kb*T)]
];

TotalRateUKHarmonic[V_,DG_,lambda_,MR_,FR_,f1_,f2_,d_,T_,mumax_,numax_,OptionsPattern[{M->MassH,Distribution->"Classical"}]]:=
Module[{m,distr,Emu,Z,Pmu,mu,nu},
   m=OptionValue[M];
   distr=OptionValue[Distribution];
	Array[Emu,mumax+1,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*HOEnergy[mu,f1],{mu,0,mumax}];
	Z = Exp[Emu[0]/(kb*T)]*QHarmonic[f1,T];
	Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
	Sum[Pmu[mu]*Sum[PartialRateUKHarmonic[V,DG,lambda,MR,FR,f1,f2,d,T,mu,nu,M->m,Distribution->distr],{nu,0,numax}],{mu,0,mumax}]
];

(*
Kinetic isotope effect
*)

(*
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
FR - donor-acceptor mode frequency (1/cm);
f1 - frequency of the protium oscillator on the left (1/cm);
f2 - frequency of the protium oscillator on the right (1/cm);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

KIEUKHarmonic[V_,DG_,lambda_,MR_,FR_,f1_,f2_,d_,T_,mumax_,numax_,OptionsPattern[{M1->MassH,M2->MassD,Distribution->"Classical"}]]:=Module[{fx1,fx2,fy1,fy2,mx,my,rate1,rate2},
   (* assume that the input frequencies are for protium oscillator *)
   mx=OptionValue[M1];
   my=OptionValue[M2];
   fx1=f1*Sqrt[MassH/mx];
   fx2=f2*Sqrt[MassH/mx];
   fy1=f1*Sqrt[MassH/my];
   fy2=f2*Sqrt[MassH/my];
	rate1=TotalRateUKHarmonic[V,DG,lambda,MR,FR,fx1,fx2,d,T,mumax,numax,M->mx,Distribution->OptionValue[Distribution]];
	rate2=TotalRateUKHarmonic[V,DG,lambda,MR,FR,fy1,fy2,d,T,mumax,numax,M->my,Distribution->OptionValue[Distribution]];
	rate1/rate2
];

(*
Table with rate channel contributions
*)

(*
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor - acceptor mode (Dalton);
FR - donor - acceptor mode frequency (1/cm);
M - mass of the particle (Dalton);
f1 - frequency for the oscillator on the left (1/cm);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TableRateUKHarmonic[V_,DG_,lambda_,MR_,FR_,f1_,f2_,d_,T_,mumax_,numax_,OptionsPattern[{M->MassH,Distribution->"Classical"}]]:=
Module[{distr,headings,dg0,la,totalrate,s,s2,alpha,dg,Ea,w,Emu,Enu,Pmu,Z,kmunu,info,Mass},
	distr=OptionValue[Distribution];
	Mass=OptionValue[M];
	headings={"\[Mu]","\[Nu]",strpmu,strdgmunu,strdgddmunu,"<"<>strs2munu<>">","exp(-\[Beta]"<>strdgddmunu<>")","%"};
	info="Isotope: "<>ToString[NumberForm[Mass,3]]<>"amu; "<>strdg00<>"="<>ToString[NumberForm[DG,4]]<>"kcal/mol; \[Lambda]="<>ToString[NumberForm[lambda,4]]<>"kcal/mol; M="<>ToString[NumberForm[MR,3]]<>"amu; \[CapitalOmega]="<>ToString[NumberForm[FR,4]]<>strcm;
	dg0=DG*kcal2au;
	la=lambda*kcal2au;
	totalrate=TotalRateUKHarmonic[V,DG,lambda,MR,FR,f1,f2,d,T,mumax,numax,M->OptionValue[M],Distribution->distr];
	Array[fc,{mumax+1,numax+1},{0,0}];
	Array[dg,{mumax+1,numax+1},{0,0}];
	Array[Ea,{mumax+1,numax+1},{0,0}];
	Array[w,{mumax+1,numax+1},{0,0}];
	Array[Emu,mumax+1,0];
	Array[Enu,numax+1,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*HOEnergy[mu,f1],{mu,0,mumax}];
	Do[Enu[nu]=kcal2au*HOEnergy[nu,f2],{nu,0,numax}];
	Z = Exp[Emu[0]/(kb*T)]*QHarmonic[f1,T];
	Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
	Do[
                fc[mu,nu]=HarmonicFranckCondonAveraged[mu,nu,f1,f2,d,T,MR,FR,M->Mass,Distribution->distr];
		dg[mu,nu]=dg0+Enu[nu]-Enu[0]-(Emu[mu]-Emu[0]);
		Ea[mu,nu]=(dg[mu,nu]+la)^2/(4*la);
		kmunu=PartialRateUKHarmonic[V,DG,lambda,MR,FR,f1,f2,d,T,mu,nu,M->Mass,Distribution->distr];
		w[mu,nu]=Pmu[mu]*kmunu*100/totalrate,{mu,0,mumax},{nu,0,numax}
	];
	Grid[Join[{{info,SpanFromLeft}},{headings},Flatten[Table[{mu,nu,ScientificForm[Pmu[mu],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[dg[mu,nu]*au2kcal,{6,3}],PaddedForm[Ea[mu,nu]*au2kcal,{6,3}],ScientificForm[fc[mu,nu],6,NumberFormat->(Row[{#1,"e",#3}]&)],ScientificForm[Exp[-Ea[mu,nu]/(kb*T)],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[w[mu,nu],{7,3}]},{mu,0,mumax},{nu,0,numax}],1]],Frame->All,Spacings->{Automatic,1}]
];

(*
Rate constant for Morse proton potentials
*)

(*
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
M - mass of the particle (Dalton);
f1 - frequency of the left (proton donor) oscillator;
f2 - frequency of the right (proton acceptor) oscillator;
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mu, nu - quantum numbers for the left and right oscillators, respectively.
Return partial rate constant in sec^-1.
*)

PartialRateUKMorse[V_,DG_,lambda_,MR_,FR_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mu_,nu_,OptionsPattern[{M->MassH,Overlap->"Numerical",Distribution->"Classical"}]]:=Module[{m,sw,distr,prefactor,v2,la,fc,dg,sq,e10,e1,e20,e2,Ea},
   m=OptionValue[M];
   sw=OptionValue[Overlap];
   distr=OptionValue[Distribution];
   fc=MorseFranckCondonAveraged[mu,nu,DE1,Beta1,DE2,Beta2,d,T,MR,FR,M->m,Overlap->sw,Distribution->distr];
   prefactor=10^12/au2ps;
   v2=(V*kcal2au)^2;
   la=lambda*kcal2au;
   dg=DG*kcal2au;
   sq=Sqrt[\[Pi]/(kb*T*la)];
   e10=kcal2au*MorseEnergy[0,Beta1,DE1,M->m];
   e1=kcal2au*MorseEnergy[mu,Beta1,DE1,M->m];
   e20=kcal2au*MorseEnergy[0,Beta2,DE2,M->m];
   e2=kcal2au*MorseEnergy[nu,Beta2,DE2,M->m];
   Ea=(dg+la+e2-e20-e1+e10)^2/(4*la);
   prefactor*v2*fc*sq*Exp[-Ea/(kb*T)]
];

TotalRateUKMorse[V_,DG_,lambda_,MR_,FR_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_,OptionsPattern[{M->MassH,Overlap->"Numerical",Distribution->"Classical"}]]:=
Module[{m,sw,distr,Emu,Z,Pmu,mu,nu},
   m=OptionValue[M];
   sw=OptionValue[Overlap];
   distr=OptionValue[Distribution];
   Array[Emu,mumax+1,0];
   Array[Pmu,mumax+1,0];
   Do[Emu[mu]=kcal2au*MorseEnergy[mu,Beta1,DE1,M->m],{mu,0,mumax}];
   Z = Exp[Emu[0]/(kb*T)]*QMorseExact[Beta1,DE1,T,M->m];
   Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
   Sum[Pmu[mu]*Sum[PartialRateUKMorse[V,DG,lambda,MR,FR,DE1,Beta1,DE2,Beta2,d,T,mu,nu,M->m,Overlap->sw,Distribution->distr],{nu,0,numax}],{mu,0,mumax}]
];

(*
Kinetic isotope effect
*)

(*
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
FR - donor-acceptor mode frequency (1/cm);
f1 - frequency of the protium oscillator on the left (1/cm);
f2 - frequency of the protium oscillator on the right (1/cm);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

KIEUKMorse[V_,DG_,lambda_,MR_,FR_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_,OptionsPattern[{M1->MassH,M2->MassD,Overlap->"Numerical",Distribution->"Classical"}]]:=Module[{sw,fx1,fx2,fy1,fy2,mx,my,rate1,rate2},
   (* assume that the input frequencies are for protium oscillator *)
   mx=OptionValue[M1];
   my=OptionValue[M2];
   sw=OptionValue[Overlap];
   rate1=TotalRateUKMorse[V,DG,lambda,MR,FR,DE1,Beta1,DE2,Beta2,d,T,mumax,numax,M->mx,Overlap->sw,Distribution->OptionValue[Distribution]];
   rate2=TotalRateUKMorse[V,DG,lambda,MR,FR,DE1,Beta1,DE2,Beta2,d,T,mumax,numax,M->my,Overlap->sw,Distribution->OptionValue[Distribution]];
   rate1/rate2
];

(*
Table with rate channel contributions
*)

(*
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor - acceptor mode (Dalton);
FR - donor - acceptor mode frequency (1/cm);
M - mass of the particle (Dalton);
f1 - frequency for the oscillator on the left (1/cm);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TableRateUKMorse[V_,DG_,lambda_,MR_,FR_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_,OptionsPattern[{M->MassH,Overlap->"Numerical",Distribution->"Classical"}]]:=
Module[{sw,distr,headings,dg0,la,totalrate,s,s2,alpha,dg,Ea,w,Emu,Enu,Pmu,Z,kmunu,fc,info,Mass},
   sw=OptionValue[Overlap];
   distr=OptionValue[Distribution];
   Mass=OptionValue[M];
   headings={"\[Mu]","\[Nu]",strpmu,strdgmunu,strdgddmunu,"<"<>strs2munu<>">","exp(-\[Beta]"<>strdgddmunu<>")","%"};
   info="Isotope: "<>ToString[NumberForm[Mass,3]]<>"amu; "<>strdg00<>"="<>ToString[NumberForm[DG,4]]<>"kcal/mol; \[Lambda]="<>ToString[NumberForm[lambda,4]]<>"kcal/mol; M="<>ToString[NumberForm[MR,3]]<>"amu; \[CapitalOmega]="<>T
   dg0=DG*kcal2au;
   la=lambda*kcal2au;
   totalrate=TotalRateUKHarmonic[V,DG,lambda,MR,FR,DE1,Beta1,DE2,Beta2,d,T,mumax,numax,M->OptionValue[M],Overlap->sw,Distribution->distr];
   Array[fc,{mumax+1,numax+1},{0,0}];
   Array[dg,{mumax+1,numax+1},{0,0}];
   Array[Ea,{mumax+1,numax+1},{0,0}];
   Array[w,{mumax+1,numax+1},{0,0}];
   Array[Emu,mumax+1,0];
   Array[Enu,numax+1,0];
   Array[Pmu,mumax+1,0];
   Do[Emu[mu]=kcal2au*MorseEnergy[mu,Beta1,DE1,M->Mass],{mu,0,mumax}];
   Do[Enu[nu]=kcal2au*MorseEnergy[nu,Beta2,DE2,M->Mass],{nu,0,numax}];
   Z = Exp[Emu[0]/(kb*T)]*QMorseExact[Beta1,DE1,T,M->Mass];
   Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
   Do[
      fc[mu,nu]=MorseFranckCondonAveraged[mu,nu,DE1,Beta1,DE2,Beta2,d,T,MR,FR,M->Mass,Overlap->sw,Distribution->distr];
      dg[mu,nu]=dg0+Enu[nu]-Enu[0]-(Emu[mu]-Emu[0]);
      Ea[mu,nu]=(dg[mu,nu]+la)^2/(4*la);
      kmunu=PartialRateUKMorse[V,DG,lambda,MR,FR,DE1,Beta1,DE2,Beta2,d,T,mu,nu,M->Mass,Overlap->sw,Distribution->distr];
      w[mu,nu]=Pmu[mu]*kmunu*100/totalrate,{mu,0,mumax},{nu,0,numax}
   ];
   Grid[Join[{{info,SpanFromLeft}},{headings},Flatten[Table[{mu,nu,ScientificForm[Pmu[mu],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[dg[mu,nu]*au2kcal,{6,3}],PaddedForm[Ea[mu,nu]*au2kcal,{6,3}],ScientificForm[fc[mu,nu],6,NumberFormat->(Row[{#1,"e",#3}]&)],ScientificForm[Exp[-Ea[mu,nu]/(kb*
T)],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[w[mu,nu],{7,3}]},{mu,0,mumax},{nu,0,numax}],1]],Frame->All,Spacings->{Automatic,1}]
];


(*
Nonadiabatic rate constant and KIE in high-temperature limit for proton donor-acceptor mode
*)

(*
Rate constant
*)

(*
Partial nonadiabatic first order rate constant k (mu->nu) (1/sec):
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
M - mass of the particle (Dalton);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mu, nu - quantum numbers for the left and right oscillators, respectively.
*)

PartialRateRhighTMorse[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mu_,nu_,OptionsPattern[{M->MassH,Overlap->"Numerical",Alpha->"Numerical"}]]:=Module[{sw,sa,s,s2,prefactor,v2,la,lalpha,latot,f,dg,rterm,sq,e10,e1,e20,e2,Ea},
	sw=OptionValue[Overlap];
	sa=OptionValue[Alpha];
	s=Switch[sw,
		"Numerical",NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
		"Analytical",MorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
		"AnalyticalSym",MorseOverlapSym[DE1,Beta1,mu,nu,d,M->OptionValue[M]],
		_,NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]];
	s2=s^2;
	prefactor=10^12/au2ps;
	v2=(V*kcal2au)^2;
	la=lambda*kcal2au;
	lalpha=Switch[Head[sa],
		Integer,kcal2au*LambdaAlphaInput[MR,sa],
		Real,kcal2au*LambdaAlphaInput[MR,sa],
		Rational,kcal2au*LambdaAlphaInput[MR,sa],
		String,kcal2au*Switch[sa,
			"Numerical",NLambdaAlphaMorse[MR,DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
			"Analytical",LambdaAlphaMorse[MR,DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
			"AnalyticalSym",LambdaAlphaMorseSym[MR,DE1,Beta1,mu,nu,d,M->OptionValue[M]],
			_,NLambdaAlphaMorse[MR,DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]],
		_,kcal2au*NLambdaAlphaMorse[MR,DE1,Beta1,DE2,Beta2,nu1,nu2,d,M->OptionValue[M]]
	];
	latot=la+lalpha;
	f=Freq*cm2au;
	dg=DG*kcal2au;
	rterm=Exp[4*kb*T*lalpha/f^2];
	sq=Sqrt[\[Pi]/(kb*T*latot)];
	e10=kcal2au*MorseEnergy[0,Beta1,DE1,M->OptionValue[M]];
	e1=kcal2au*MorseEnergy[mu,Beta1,DE1,M->OptionValue[M]];
	e20=kcal2au*MorseEnergy[0,Beta2,DE2,M->OptionValue[M]];
	e2=kcal2au*MorseEnergy[nu,Beta2,DE2,M->OptionValue[M]];
	Ea=(dg+latot+e2-e20-e1+e10)^2/(4*latot);
	prefactor*v2*s2*rterm*sq*Exp[-Ea/(kb*T)]
];

(*
Partial nonadiabatic first order rate constant k (mu->nu) (1/sec):
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
T - temperature (K);
mu, nu - quantum numbers for the left and right oscillators, respectively.
*)

PartialRateRhighTFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mu_,nu_,OptionsPattern[{Alpha->"Numerical"}]]:=
Module[{sa,s,s2,prefactor,v2,la,lalpha,latot,f,dg,rterm,sq,e10,e1,e20,e2,Ea,Ns},
        Ns = Length[fghList[[1,1]]];
        If[mu+1>Ns||nu+1>Ns,Message[PartialRateRhighTFGH::numberofstates];Abort[]];
	sa=OptionValue[Alpha];
	s=fghList[[3]][[mu+1,nu+1]];
	s2=s^2;
	prefactor=10^12/au2ps;
	v2=(V*kcal2au)^2;
	la=lambda*kcal2au;
	lalpha=Switch[Head[sa],
		Integer,kcal2au*LambdaAlphaInput[MR,sa],
		Real,kcal2au*LambdaAlphaInput[MR,sa],
		Rational,kcal2au*LambdaAlphaInput[MR,sa],
		String,Switch[sa,
			"Numerical",(fghList[[4]][[mu+1,nu+1]]/a2bohr)^2/(2*MR*Dalton),
			          _,(fghList[[4]][[mu+1,nu+1]]/a2bohr)^2/(2*MR*Dalton)],
		_,(fghList[[4]][[mu+1,nu+1]]/a2bohr)^2/(2*MR*Dalton)
	];
	latot=la+lalpha;
	f=Freq*cm2au;
	dg=DG*kcal2au;
	rterm=Exp[4*kb*T*lalpha/f^2];
	sq=Sqrt[\[Pi]/(kb*T*latot)];
	e10=kcal2au*fghList[[1,1,1]];
	e1=kcal2au*fghList[[1,1,mu+1]];
	e20=kcal2au*fghList[[2,1,1]];
	e2=kcal2au*fghList[[2,1,nu+1]];
	Ea=(dg+latot+e2-e20-e1+e10)^2/(4*latot);
	prefactor*v2*s2*rterm*sq*Exp[-Ea/(kb*T)]
];

(*
Total nonadiabatic first order rate constant k (1/sec):
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
M - mass of the particle (Dalton);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TotalRateRhighTMorse[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_,OptionsPattern[{M->MassH,Overlap->"Numerical",Alpha->"Numerical"}]]:=
Module[{Emu,Z,Pmu,mu,nu},
	Array[Emu,mumax+1,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*MorseEnergy[mu,Beta1,DE1,M->OptionValue[M]],{mu,0,mumax}];
	Z = Exp[Emu[0]/(kb*T)]*QMorseExact[Beta1,DE1,T,M->OptionValue[M]];
	Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
	Sum[Pmu[mu]*Sum[PartialRateRhighTMorse[V,DG,lambda,MR,Freq,DE1,Beta1,DE2,Beta2,d,T,mu,nu,M->OptionValue[M],Overlap->OptionValue[Overlap],Alpha->OptionValue[Alpha]],{nu,0,numax}],{mu,0,mumax}]
];

(*
Total nonadiabatic first order rate constant k (1/sec):
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TotalRateRhighTFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_,OptionsPattern[{Alpha->"Numerical"}]]:=
Module[{Emu,Z,Pmu,mu,nu,Ns},
        Ns = Length[fghList[[1,1]]];
        If[mumax+1>Ns||numax+1>Ns,Message[TotalRateRhighTFGH::numberofstates];Abort[]];
	Array[Emu,Ns,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*fghList[[1,1,mu+1]],{mu,0,Ns-1}];
	Z=Sum[Exp[-(Emu[mu]-Emu[0])/(kb*T)],{mu,0,Ns-1}];
	Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
	Sum[Pmu[mu]*Sum[PartialRateRhighTFGH[fghList,V,DG,lambda,MR,Freq,T,mu,nu,Alpha->OptionValue[Alpha]],{nu,0,numax}],{mu,0,mumax}]
];

(*
Kinetic isotope effect
*)

(*
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

KIERhighTMorse[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_,OptionsPattern[{M1->MassH,M2->MassD,Overlap->"Numerical",Alpha->"Numerical"}]]:=Module[{rate1,rate2},
	rate1=TotalRateRhighTMorse[V,DG,lambda,MR,Freq,DE1,Beta1,DE2,Beta2,d,T,mumax,numax,M->OptionValue[M1],Overlap->OptionValue[Overlap],Alpha->OptionValue[Alpha]];
	rate2=TotalRateRhighTMorse[V,DG,lambda,MR,Freq,DE1,Beta1,DE2,Beta2,d,T,mumax,numax,M->OptionValue[M2],Overlap->OptionValue[Overlap],Alpha->OptionValue[Alpha]];
	rate1/rate2
];

(*
fghList1 - nested list with energies, wavefunctions, overlap integrals and alphas for isotope 1;
fghList2 - nested list with energies, wavefunctions, overlap integrals and alphas for isotope 2;
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

KIERhighTFGH[fghList1_,fghList2_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_,OptionsPattern[{Alpha1->"Numerical",Alpha2->"Numerical"}]]:=Module[{rate1,rate2},
	rate1=TotalRateRhighTFGH[fghList1,V,DG,lambda,MR,Freq,T,mumax,numax,Alpha->OptionValue[Alpha1]];
	rate2=TotalRateRhighTFGH[fghList2,V,DG,lambda,MR,Freq,T,mumax,numax,Alpha->OptionValue[Alpha2]];
	rate1/rate2
];

(*
Table with rate channel contributions
*)

(*
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor - acceptor mode (Dalton);
Freq - donor - acceptor mode frequency (1/cm);
M - mass of the particle (Dalton);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TableRateRhighTMorse[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_,OptionsPattern[{M->MassH,Overlap->"Numerical",Alpha->"Numerical"}]]:=
Module[{sw,sa,headings,dg0,la,totalrate,s,s2,alpha,dg,Ea,w,Emu,Enu,Pmu,Z,kmunu,info,Mass},
	sw=OptionValue[Overlap];
	sa=OptionValue[Alpha];
	Mass=OptionValue[M];
	headings={"\[Mu]","\[Nu]",strpmu,strdgmunu,strdgddmunu,strs2munu,stramunu,"exp(-\[Beta]"<>strdgddmunu<>")","%"};
	info="Isotope: "<>ToString[NumberForm[Mass,3]]<>"amu; "<>strdg00<>"="<>ToString[NumberForm[DG,4]]<>"kcal/mol; \[Lambda]="<>ToString[NumberForm[lambda,4]]<>"kcal/mol; M="<>ToString[NumberForm[MR,3]]<>"amu; \[CapitalOmega]="<>ToString[NumberForm[Freq,4]]<>strcm;
	dg0=DG*kcal2au;
	la=lambda*kcal2au;
	totalrate=TotalRateRhighTMorse[V,DG,lambda,MR,Freq,DE1,Beta1,DE2,Beta2,d,T,mumax,numax,M->OptionValue[M],Overlap->OptionValue[M],Alpha->OptionValue[Alpha]];
	Array[s2,{mumax+1,numax+1},{0,0}];
	Array[dg,{mumax+1,numax+1},{0,0}];
	Array[Ea,{mumax+1,numax+1},{0,0}];
	Array[w,{mumax+1,numax+1},{0,0}];
	Array[Emu,mumax+1,0];
	Array[Enu,numax+1,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*MorseEnergy[mu,Beta1,DE1,M->OptionValue[M]],{mu,0,mumax}];
	Do[Enu[nu]=kcal2au*MorseEnergy[nu,Beta2,DE2,M->OptionValue[M]],{nu,0,numax}];
	Z = Exp[Emu[0]/(kb*T)]*QMorseExact[Beta1,DE1,T,M->OptionValue[M]];
	Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
	Do[
		s=Switch[sw,
			"Numerical",NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
			"Analytical",MorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
			"AnalyticalSym",MorseOverlapSym[DE1,Beta1,mu,nu,d,M->OptionValue[M]],
			_,NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]];
		s2[mu,nu]=s^2;
	        alpha[mu,nu]=Switch[Head[sa],
		        Integer,sa/bohr2a,
		        Real,sa/bohr2a,
		        Rational,sa/bohr2a,
		        String,Switch[sa,
			        "Numerical",NAlphaMorse[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]/bohr2a,
			        "Analytical",AlphaMorse[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]/bohr2a,
			        "AnalyticalSym",AlphaMorseSym[DE1,Beta1,mu,nu,d,M->OptionValue[M]]/bohr2a,
			        _,NAlphaMorse[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]/bohr2a
			        ],
		        _,NAlphaMorse[DE1,Beta1,DE2,Beta2,nu1,nu2,d,M->OptionValue[M]]
	                ];
		dg[mu,nu]=dg0+Enu[nu]-Enu[0]-(Emu[mu]-Emu[0]);
		Ea[mu,nu]=(dg[mu,nu]+la)^2/(4*la);
		kmunu=PartialRateRhighTMorse[V,DG,lambda,MR,Freq,DE1,Beta1,DE2,Beta2,d,T,mu,nu,M->OptionValue[M],Overlap->sw,Alpha->sa];
		w[mu,nu]=Pmu[mu]*kmunu*100/totalrate,{mu,0,mumax},{nu,0,numax}
	];
	Grid[Join[{{info,SpanFromLeft}},{headings},Flatten[Table[{mu,nu,ScientificForm[Pmu[mu],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[dg[mu,nu]*au2kcal,{6,3}],PaddedForm[Ea[mu,nu]*au2kcal,{6,3}],ScientificForm[s2[mu,nu],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[alpha[mu,nu],{6,3}],ScientificForm[Exp[-Ea[mu,nu]/(kb*T)],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[w[mu,nu],{7,3}]},{mu,0,mumax},{nu,0,numax}],1]],Frame->All,Spacings->{Automatic,1}]
];

(*
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor - acceptor mode (Dalton);
Freq - donor - acceptor mode frequency (1/cm);
M - mass of the particle (Dalton);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TableRateRhighTFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_,OptionsPattern[{Alpha->"Numerical"}]]:=
Module[{sa,headings,dg0,la,totalrate,s,s2,alpha,dg,Ea,w,Emu,Enu,Pmu,Z,kmunu,info,Ns},
        Ns = Length[fghList[[1,1]]];
        If[mumax+1>Ns||numax+1>Ns,Message[TableRateRhighTFGH::numberofstates];Abort[]];
	sa=OptionValue[Alpha];
	headings={"\[Mu]","\[Nu]",strpmu,strdgmunu,strdgddmunu,strs2munu,stramunu,"exp(-\[Beta]"<>strdgddmunu<>")","%"};
	info=strdg00<>"="<>ToString[NumberForm[DG,4]]<>"kcal/mol; \[Lambda]="<>ToString[NumberForm[lambda,4]]<>"kcal/mol; M="<>ToString[NumberForm[MR,3]]<>"amu; \[CapitalOmega]="<>ToString[NumberForm[Freq,4]]<>strcm;
	dg0=DG*kcal2au;
	la=lambda*kcal2au;
	totalrate=TotalRateRhighTFGH[fghList,V,DG,lambda,MR,Freq,T,mumax,numax,Alpha->OptionValue[Alpha]];
	Array[s2,{mumax+1,numax+1},{0,0}];
	Array[dg,{mumax+1,numax+1},{0,0}];
	Array[Ea,{mumax+1,numax+1},{0,0}];
	Array[w,{mumax+1,numax+1},{0,0}];
	Array[Emu,Ns,0];
	Array[Enu,Ns,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*fghList[[1,1,mu+1]],{mu,0,Ns-1}];
	Do[Enu[nu]=kcal2au*fghList[[2,1,nu+1]],{nu,0,Ns-1}];
	Z=Sum[Exp[-(Emu[mu]-Emu[0])/(kb*T)],{mu,0,Ns-1}];
	Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
	Do[
		s=fghList[[3]][[mu+1,nu+1]];
		s2[mu,nu]=s^2;
	        alpha[mu,nu]=Switch[Head[sa],
		        Integer,sa/bohr2a,
		        Real,sa/bohr2a,
		        Rational,sa/bohr2a,
		        String,Switch[sa,
			        "Numerical",fghList[[4]][[mu+1,nu+1]],
			        _,fghList[[4]][[mu+1,nu+1]]
			        ],
		        _,fghList[[4]][[mu+1,nu+1]]
	                ];
		dg[mu,nu]=dg0+Enu[nu]-Enu[0]-(Emu[mu]-Emu[0]);
		Ea[mu,nu]=(dg[mu,nu]+la)^2/(4*la);
		kmunu=PartialRateRhighTFGH[fghList,V,DG,lambda,MR,Freq,T,mu,nu,Alpha->sa];
		w[mu,nu]=Pmu[mu]*kmunu*100/totalrate,{mu,0,mumax},{nu,0,numax}
	];
	Grid[Join[{{info,SpanFromLeft}},{headings},Flatten[Table[{mu,nu,ScientificForm[Pmu[mu],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[dg[mu,nu]*au2kcal,{6,3}],PaddedForm[Ea[mu,nu]*au2kcal,{6,3}],ScientificForm[s2[mu,nu],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[alpha[mu,nu],{6,3}],ScientificForm[Exp[-Ea[mu,nu]/(kb*T)],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[w[mu,nu],{7,3}]},{mu,0,mumax},{nu,0,numax}],1]],Frame->All,Spacings->{Automatic,1}]
];

(*
Nonadiabatic rate constant and KIE in low-temperature limit for proton donor-acceptor mode
*)

(*
Rate constant
*)

(*
Partial nonadiabatic first order rate constant k (mu->nu) (1/sec):
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
M - mass of the particle (Dalton);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mu, nu - quantum numbers for the left and right oscillators, respectively.
*)

PartialRateRlowTMorse[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mu_,nu_,OptionsPattern[{M->MassH,Overlap->"Numerical",Alpha->"Numerical"}]]:=Module[{sw,sa,s,s2,prefactor,v2,la,lalpha,f,dg,rterm,sq,e10,e1,e20,e2,Ea},
	sw=OptionValue[Overlap];
	sa=OptionValue[Alpha];
	s=Switch[sw,
		"Numerical",NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
		"Analytical",MorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
		"AnalyticalSym",MorseOverlapSym[DE1,Beta1,mu,nu,d,M->OptionValue[M]],
		_,NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]];
	s2=s^2;
	prefactor=10^12/au2ps;
	v2=(V*kcal2au)^2;
	la=lambda*kcal2au;
	lalpha=Switch[Head[sa],
		Integer,kcal2au*LambdaAlphaInput[MR,sa],
		Real,kcal2au*LambdaAlphaInput[MR,sa],
		Rational,kcal2au*LambdaAlphaInput[MR,sa],
		String,kcal2au*Switch[sa,
			"Numerical",NLambdaAlphaMorse[MR,DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
			"Analytical",LambdaAlphaMorse[MR,DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
			"AnalyticalSym",LambdaAlphaMorseSym[MR,DE1,Beta1,mu,nu,d,M->OptionValue[M]],
			_,NLambdaAlphaMorse[MR,DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]],
		_,kcal2au*NLambdaAlphaMorse[MR,DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]
	];
	f=Freq*cm2au;
	dg=DG*kcal2au;
	rterm=Exp[lalpha*Coth[hbar*f/(2*kb*T)]/(hbar*f)];
	sq=Sqrt[\[Pi]/(kb*T*la)];
	e10=kcal2au*MorseEnergy[0,Beta1,DE1,M->OptionValue[M]];
	e1=kcal2au*MorseEnergy[mu,Beta1,DE1,M->OptionValue[M]];
	e20=kcal2au*MorseEnergy[0,Beta2,DE2,M->OptionValue[M]];
	e2=kcal2au*MorseEnergy[nu,Beta2,DE2,M->OptionValue[M]];
	Ea=(dg+la+e2-e20-e1+e10)^2/(4*la);
	prefactor*v2*s2*rterm*sq*Exp[-Ea/(kb*T)]
];

(*
Partial nonadiabatic first order rate constant k (mu->nu) (1/sec):
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
T - temperature (K);
mu, nu - quantum numbers for the left and right oscillators, respectively.
*)

PartialRateRlowTFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mu_,nu_,OptionsPattern[{Alpha->"Numerical"}]]:=
Module[{sa,s,s2,prefactor,v2,la,lalpha,f,dg,rterm,sq,e10,e1,e20,e2,Ea,Ns},
        Ns = Length[fghList[[1,1]]];
        If[mu+1>Ns||nu+1>Ns,Message[PartialRateRlowTFGH::numberofstates];Abort[]];
	sa=OptionValue[Alpha];
	s=fghList[[3]][[mu+1,nu+1]];
	s2=s^2;
	prefactor=10^12/au2ps;
	v2=(V*kcal2au)^2;
	la=lambda*kcal2au;
	lalpha=Switch[Head[sa],
		Integer,kcal2au*LambdaAlphaInput[MR,sa],
		Real,kcal2au*LambdaAlphaInput[MR,sa],
		Rational,kcal2au*LambdaAlphaInput[MR,sa],
		String,Switch[sa,
			"Numerical",(fghList[[4]][[mu+1,nu+1]]/a2bohr)^2/(2*MR*Dalton),
			          _,(fghList[[4]][[mu+1,nu+1]]/a2bohr)^2/(2*MR*Dalton)],
		_,(fghList[[4]][[mu+1,nu+1]]/a2bohr)^2/(2*MR*Dalton)
	];
	f=Freq*cm2au;
	dg=DG*kcal2au;
	rterm=Exp[lalpha*Coth[hbar*f/(2*kb*T)]/(hbar*f)];
	sq=Sqrt[\[Pi]/(kb*T*la)];
	e10=kcal2au*fghList[[1,1,1]];
	e1=kcal2au*fghList[[1,1,mu+1]];
	e20=kcal2au*fghList[[2,1,1]];
	e2=kcal2au*fghList[[2,1,nu+1]];
	Ea=(dg+la+e2-e20-e1+e10)^2/(4*la);
	prefactor*v2*s2*rterm*sq*Exp[-Ea/(kb*T)]
];

(*
Total nonadiabatic first order rate constant k (1/sec):
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
M - mass of the particle (Dalton);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TotalRateRlowTMorse[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_,OptionsPattern[{M->MassH,Overlap->"Numerical",Alpha->"Numerical"}]]:=
Module[{Emu,Z,Pmu,mu,nu},
	Array[Emu,mumax+1,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*MorseEnergy[mu,Beta1,DE1,M->OptionValue[M]],{mu,0,mumax}];
	Z = Exp[Emu[0]/(kb*T)]*QMorseExact[Beta1,DE1,T,M->OptionValue[M]];
	Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
	Sum[Pmu[mu]*Sum[PartialRateRlowTMorse[V,DG,lambda,MR,Freq,DE1,Beta1,DE2,Beta2,d,T,mu,nu,M->OptionValue[M],Overlap->OptionValue[Overlap],Alpha->OptionValue[Alpha]],{nu,0,numax}],{mu,0,mumax}]
];

(*
Total nonadiabatic first order rate constant k (1/sec):
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TotalRateRlowTFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_,OptionsPattern[{Alpha->"Numerical"}]]:=
Module[{Emu,Z,Pmu,mu,nu,Ns},
        Ns = Length[fghList[[1,1]]];
        If[mumax+1>Ns||numax+1>Ns,Message[TotalRateRlowTFGH::numberofstates];Abort[]];
	Array[Emu,Ns,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*fghList[[1,1,mu+1]],{mu,0,Ns-1}];
	Z=Sum[Exp[-(Emu[mu]-Emu[0])/(kb*T)],{mu,0,Ns-1}];
	Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
	Sum[Pmu[mu]*Sum[PartialRateRlowTFGH[fghList,V,DG,lambda,MR,Freq,T,mu,nu,Alpha->OptionValue[Alpha]],{nu,0,numax}],{mu,0,mumax}]
];

(*
kinetic isotope effect
*)

(*
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

KIERlowTMorse[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_,OptionsPattern[{M1->MassH,M2->MassD,Overlap->"Numerical",Alpha->"Numerical"}]]:=Module[{rate1,rate2},
	rate1=TotalRateRlowTMorse[V,DG,lambda,MR,Freq,DE1,Beta1,DE2,Beta2,d,T,mumax,numax,M->OptionValue[M1],Overlap->OptionValue[Overlap],Alpha->OptionValue[Alpha]];
	rate2=TotalRateRlowTMorse[V,DG,lambda,MR,Freq,DE1,Beta1,DE2,Beta2,d,T,mumax,numax,M->OptionValue[M2],Overlap->OptionValue[Overlap],Alpha->OptionValue[Alpha]];
	rate1/rate2
];

(*
fghList1 - nested list with energies, wavefunctions, overlap integrals and alphas for isotope 1;
fghList2 - nested list with energies, wavefunctions, overlap integrals and alphas for isotope 2;
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor-acceptor mode (Dalton);
Freq - donor-acceptor mode frequency (1/cm);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

KIERlowTFGH[fghList1_,fghList2_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_,OptionsPattern[{Alpha1->"Numerical",Alpha2->"Numerical"}]]:=Module[{rate1,rate2},
	rate1=TotalRateRlowTFGH[fghList1,V,DG,lambda,MR,Freq,T,mumax,numax,Alpha->OptionValue[Alpha1]];
	rate2=TotalRateRlowTFGH[fghList2,V,DG,lambda,MR,Freq,T,mumax,numax,Alpha->OptionValue[Alpha2]];
	rate1/rate2
];

(*
Table with rate channel contributions
*)

(*
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor - acceptor mode (Dalton);
Freq - donor - acceptor mode frequency (1/cm);
M - mass of the particle (Dalton);
DE1 - dissociation energy for the oscillator on the left (kcal/mol);
Beta1 - beta parameter for the oscillator on the left (1/Bohr);
DE2 - dissociation energy for the oscillator on the right (kcal/mol);
Beta2 - beta parameter for the oscillator on the right (1/Bohr);
d - equilibrium distance between the minima of the morse potentials (Bohr);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TableRateRlowTMorse[V_,DG_,lambda_,MR_,Freq_,DE1_,Beta1_,DE2_,Beta2_,d_,T_,mumax_,numax_,OptionsPattern[{M->MassH,Overlap->"Numerical",Alpha->"Numerical"}]]:=
Module[{sw,sa,headings,dg0,la,totalrate,s,s2,alpha,dg,Ea,w,Emu,Enu,Pmu,Z,kmunu,info,Mass},
	sw=OptionValue[Overlap];
	sa=OptionValue[Alpha];
	Mass=OptionValue[M];
	headings={"\[Mu]","\[Nu]",strpmu,strdgmunu,strdgddmunu,strs2munu,stramunu,"exp(-\[Beta]"<>strdgddmunu<>")","%"};
	info="Isotope: "<>ToString[NumberForm[Mass,3]]<>"amu; "<>strdg00<>"="<>ToString[NumberForm[DG,4]]<>"kcal/mol; \[Lambda]="<>ToString[NumberForm[lambda,4]]<>"kcal/mol; M="<>ToString[NumberForm[MR,3]]<>"amu; \[CapitalOmega]="<>ToString[NumberForm[Freq,4]]<>strcm;
	dg0=DG*kcal2au;
	la=lambda*kcal2au;
	totalrate=TotalRateRlowTMorse[V,DG,lambda,MR,Freq,DE1,Beta1,DE2,Beta2,d,T,mumax,numax,M->MassH,Overlap->OptionValue[Overlap],Alpha->OptionValue[Alpha]];
	Array[s2,{mumax+1,numax+1},{0,0}];
	Array[alpha,{mumax+1,numax+1},{0,0}];
	Array[dg,{mumax+1,numax+1},{0,0}];
	Array[Ea,{mumax+1,numax+1},{0,0}];
	Array[w,{mumax+1,numax+1},{0,0}];
	Array[Emu,mumax+1,0];
	Array[Enu,numax+1,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*MorseEnergy[mu,Beta1,DE1,M->OptionValue[M]],{mu,0,mumax}];
	Do[Enu[nu]=kcal2au*MorseEnergy[nu,Beta2,DE2,M->OptionValue[M]],{nu,0,numax}];
	Z = Exp[Emu[0]/(kb*T)]*QMorseExact[Beta1,DE1,T,M->OptionValue[M]];
	Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
	Do[
		s=Switch[sw,
			"Numerical",NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
			"Analytical",MorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],
			"AnalyticalSym",MorseOverlapSym[DE1,Beta1,mu,nu,d,M->OptionValue[M]],
			_,NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]
			];
		s2[mu,nu]=s^2;
	        alpha[mu,nu]=Switch[Head[sa],
		        Integer,sa/bohr2a,
		        Real,sa/bohr2a,
		        Rational,sa/bohr2a,
		        String,Switch[sa,
			        "Numerical",NAlphaMorse[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]/bohr2a,
			        "Analytical",AlphaMorse[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]/bohr2a,
			        "AnalyticalSym",AlphaMorseSym[DE1,Beta1,mu,nu,d,M->OptionValue[M]]/bohr2a,
			        _,NAlphaMorse[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]]/bohr2a
			        ],
		        _,NAlphaMorse[DE1,Beta1,DE2,Beta2,nu1,nu2,d,M->OptionValue[M]]
	                ];
		dg[mu,nu]=dg0+Enu[nu]-Enu[0]-(Emu[mu]-Emu[0]);
		Ea[mu,nu]=(dg[mu,nu]+la)^2/(4*la);
		kmunu=PartialRateRlowTMorse[V,DG,lambda,MR,Freq,DE1,Beta1,DE2,Beta2,d,T,mu,nu,M->MassH,Overlap->OptionValue[Overlap],Alpha->OptionValue[Alpha]];
		w[mu,nu]=Pmu[mu]*kmunu*100/totalrate,{mu,0,mumax},{nu,0,numax}
	];
	Grid[Join[{{info,SpanFromLeft}},{headings},Flatten[Table[{mu,nu,ScientificForm[Pmu[mu],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[dg[mu,nu]*au2kcal,{6,3}],PaddedForm[Ea[mu,nu]*au2kcal,{6,3}],ScientificForm[s2[mu,nu],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[alpha[mu,nu],{6,3}],ScientificForm[Exp[-Ea[mu,nu]/(kb*T)],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[w[mu,nu],{7,3}]},{mu,0,mumax},{nu,0,numax}],1]],Frame->All,Spacings->{Automatic,1}]
];



(*
fghList - nested list with energies, wavefunctions, overlap integrals and alphas;
V - electronic coupling (kcal/mol);
DG - reaction free energy (free energy bias) (kcal/mol);
lambda - solvent reorganization energy (kcal/mol);
MR - reduced mass of the donor - acceptor mode (Dalton);
Freq - donor - acceptor mode frequency (1/cm);
T - temperature (K);
mumax - highest quantum number for the reactant vibronic states;
numax - highest quantum number for the product vibronic states.
*)

TableRateRlowTFGH[fghList_,V_,DG_,lambda_,MR_,Freq_,T_,mumax_,numax_,OptionsPattern[{Alpha->"Numerical"}]]:=
Module[{sa,headings,dg0,la,totalrate,s,s2,alpha,dg,Ea,w,Emu,Enu,Pmu,Z,kmunu,info,Ns},
        Ns = Length[fghList[[1,1]]];
        If[mumax+1>Ns||numax+1>Ns,Message[TableRateRlowTFGH::numberofstates];Abort[]];
	sa=OptionValue[Alpha];
	headings={"\[Mu]","\[Nu]",strpmu,strdgmunu,strdgddmunu,strs2munu,stramunu,"exp(-\[Beta]"<>strdgddmunu<>")","%"};
	info=strdg00<>"="<>ToString[NumberForm[DG,4]]<>"kcal/mol; \[Lambda]="<>ToString[NumberForm[lambda,4]]<>"kcal/mol; M="<>ToString[NumberForm[MR,3]]<>"amu; \[CapitalOmega]="<>ToString[NumberForm[Freq,4]]<>strcm;
	dg0=DG*kcal2au;
	la=lambda*kcal2au;
	totalrate=TotalRateRlowTFGH[fghList,V,DG,lambda,MR,Freq,T,mumax,numax,Alpha->OptionValue[Alpha]];
	Array[s2,{mumax+1,numax+1},{0,0}];
	Array[alpha,{mumax+1,numax+1},{0,0}];
	Array[dg,{mumax+1,numax+1},{0,0}];
	Array[Ea,{mumax+1,numax+1},{0,0}];
	Array[w,{mumax+1,numax+1},{0,0}];
	Array[Emu,Ns,0];
	Array[Enu,Ns,0];
	Array[Pmu,mumax+1,0];
	Do[Emu[mu]=kcal2au*fghList[[1,1,mu+1]],{mu,0,Ns-1}];
	Do[Enu[nu]=kcal2au*fghList[[2,1,nu+1]],{nu,0,Ns-1}];
	Z=Sum[Exp[-(Emu[mu]-Emu[0])/(kb*T)],{mu,0,Ns-1}];
	Do[Pmu[mu]=Exp[-(Emu[mu]-Emu[0])/(kb*T)]/Z,{mu,0,mumax}];
	Do[
		s=fghList[[3]][[mu+1,nu+1]];
		s2[mu,nu]=s^2;
	        alpha[mu,nu]=Switch[Head[sa],
		        Integer,sa/bohr2a,
		        Real,sa/bohr2a,
		        Rational,sa/bohr2a,
		        String,Switch[sa,
			        "Numerical",fghList[[4]][[mu+1,nu+1]],
			        _,          fghList[[4]][[mu+1,nu+1]]
			        ],
		        _,fghList[[4]][[mu+1,nu+1]]
	                ];
		dg[mu,nu]=dg0+Enu[nu]-Enu[0]-(Emu[mu]-Emu[0]);
		Ea[mu,nu]=(dg[mu,nu]+la)^2/(4*la);
		kmunu=PartialRateRlowTFGH[fghList,V,DG,lambda,MR,Freq,T,mu,nu,Alpha->OptionValue[Alpha]];
		w[mu,nu]=Pmu[mu]*kmunu*100/totalrate,{mu,0,mumax},{nu,0,numax}
	];
	Grid[Join[{{info,SpanFromLeft}},{headings},Flatten[Table[{mu,nu,ScientificForm[Pmu[mu],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[dg[mu,nu]*au2kcal,{6,3}],PaddedForm[Ea[mu,nu]*au2kcal,{6,3}],ScientificForm[s2[mu,nu],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[alpha[mu,nu],{6,3}],ScientificForm[Exp[-Ea[mu,nu]/(kb*T)],6,NumberFormat->(Row[{#1,"e",#3}]&)],PaddedForm[w[mu,nu],{7,3}]},{mu,0,mumax},{nu,0,numax}],1]],Frame->All,Spacings->{Automatic,1}]
];

(*
Fermi Distribution
*)

(*
eps - epsilon in Hartrees;
T - temperature in Kelvin;
Returns fermi distribution value at a prescribed energy and temperature;
*)

Fermi[eps_,T_] := Module[{},
	1/(1+Exp[eps/(kb*T)])
];

(*
Harmonic Alpha
*)

(*
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
d - distance between the minima of the displaced potentials in bohr;
Returns the alpha parameter in inverse Bohr for the harmonic oscillator approximation;
*)

AlphaHarmonic[mu_,nu_,omega1_,omega2_,d_,OptionsPattern[{M->MassH}]]:=Module[{s,x},
	s=HarmonicOverlap[mu,nu,omega1,omega2,d,M->OptionValue[M]];
	-(1/s)*ND[HarmonicOverlap[mu,nu,omega1,omega2,x,M->OptionValue[M]],x,d]
];

(*
Reorganization Energies
*)

(*
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equilibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
Returns the R mode reorganization energy in kcal/mol; 
*)

LambdaR[MR_,omega_,dR_]:=Module[{},
	au2kcal*MR*Dalton*(omega*cm2au)^2*dR^2/2
];

(*
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
Returns the total reorganization energy in kcal/mol; 
*)

NLambdaMorse[lambda_,MR_,omega_,dR_,DE1_,Beta1_,DE2_,Beta2_,mu_,nu_,d_,OptionsPattern[{M->MassH,Step->0.0001}]]:= Module[{},
	lambda+LambdaR[MR,omega,dR]+NLambdaAlphaMorse[MR,DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M],Step->OptionValue[Step]]
];

(*
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - the temperature in Kelvin;
Returns the adjusted reorganization energy for oxidation processes in kcal/mol; 
*)

NLambdaMorsePlus[lambda_,MR_,omega_,dR_,DE1_,Beta1_,DE2_,Beta2_,mu_,nu_,d_,T_,OptionsPattern[{M->MassH,Step->0.0001}]]:= Module[{},
	NLambdaMorse[lambda,MR,omega,dR,DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M],Step->OptionValue[Step]]+2*NAlphaMorse[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M],Step->OptionValue[Step]]*dR*kb*T*au2kcal
];

(*
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - the temperature in Kelvin;
Returns the adjusted reorganization energy for reductive processes in kcal/mol; 
*)

NLambdaMorseMinus[lambda_,MR_,omega_,dR_,DE1_,Beta1_,DE2_,Beta2_,mu_,nu_,d_,T_,OptionsPattern[{M->MassH,Step->0.0001}]]:=Module[{},
	NLambdaMorse[lambda,MR,omega,dR,DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M],Step->OptionValue[Step]]-2*NAlphaMorse[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M],Step->OptionValue[Step]]*dR*kb*T*au2kcal
];

(*
MR - reduced mass of the donor-acceptor mode in Daltons;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
d - distance between the minima of the displaced potentials in bohr;
Returns the coupling reorganization energy in kcal/mol for the harmonic oscillator approximation;
*)

LambdaAlphaHarmonic[MR_,mu_,nu_,omega1_,omega2_,d_,OptionsPattern[M->MassH]]:=Module[{mur,alpha},
	mur=MR*Dalton;
	alpha=AlphaHarmonic[mu,nu,omega1,omega2,d,M->OptionValue[M]];
	au2kcal*hbar^2*alpha^2/(2*mur)
];

(*
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
Returns the total reorganization energy in kcal/mol for the harmonic oscillator approximation;
*)

LambdaHarmonic[lambda_,MR_,dR_,omega_,omega1_,omega2_,mu_,nu_,d_,OptionsPattern[M->MassH]]:=Module[{},
lambda+LambdaR[MR,omega,dR]+LambdaAlphaHarmonic[MR,mu,nu,omega1,omega2,d,M->OptionValue[M]]
];

(*
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
d - distance between the minima of the displaced potentials in bohr;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
T - temperature in Kelvin;
Returns the adjusted reorganization energy in kcal/mol for the oxidative process for the harmonic oscillator approximation;
*)

LambdaHarmonicPlus[lambda_,MR_,dR_,omega_,d_,omega1_,omega2_,mu_,nu_,T_,OptionsPattern[M-> MassH]]:=Module[{},
LambdaHarmonic[lambda,MR,dR,omega,omega1,omega2,mu,nu,d,M->OptionValue[M]]+2*AlphaHarmonic[mu,nu,omega1,omega2,d,M->OptionValue[M]]*dR*kb*T*au2kcal
];

(*
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
d - distance between the minima of the displaced potentials in bohr;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
T - temperature in Kelvin;
Returns the adjusted reorganization energy in kcal/mol for the reductive process for the harmonic oscillator approximation;
*)

LambdaHarmonicMinus[lambda_,MR_,dR_,omega_,d_,omega1_,omega2_,mu_,nu_,T_,OptionsPattern[M->MassH]]:=Module[{},
	LambdaHarmonic[lambda,MR,dR,omega,omega1,omega2,mu,nu,d,M-> OptionValue[M]]-2*AlphaHarmonic[mu,nu,omega1,omega2,d,M->OptionValue[M]]*dR*kb*T*au2kcal
];

(*
Boltzmann Weights
*)

(*
Beta1 - beta parameter for oscillator in inverse Bohr;
DE1 - dissociation energy of the oscillator in kcal/mol;
mu - quantum number;
T - temperature in Kelvin;
Returns the Boltzmann weight per quantum number for a Morse approximation;
*)

BoltzmannMorse[Beta1_,DE1_,mu_,T_,OptionsPattern[{M->MassH}]]:=
Module[{mass,Z},
	mass=OptionValue[M];
	Z = QMorseExact[Beta1,DE1,T,M->mass];
	Exp[-kcal2au*MorseEnergy[mu,Beta1,DE1,M->mass]/(kb*T)]/Z
];

(*
omega1 - frequency of the harmonic oscillator in wavenumbers;
mu - quantum number;
T - temperature in Kelvin;
Returns the Boltzmann weight per quantum number for a harmonic approximation;
*)

BoltzmannHarmonic[omega1_,mu_,T_]:=Module[{},
	Exp[-(kcal2au*HOEnergy[omega1,mu])/(kb*T)]*(Exp[(omega1*cm2au*hbar)/(2*kb*T)]-Exp[-(omega1*cm2au*hbar)/(2*kb*T)])
];

(*
Gibbs Free Energy
*)

(*
eta - overpotential in volts;
eps - energy in Hartrees (used in Fermi distribution);
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
T - temperature in Kelvin;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
Returns the reaction free energy for a Morse approximation in Hartrees;
*)

dGMorse[eta_,eps_,Beta1_,DE1_,Beta2_,DE2_,T_,mu_,nu_,OptionsPattern[{M->MassH}]]:=
Module[{mass,QI,QII},
	mass=OptionValue[M];
	QI = QMorseExact[Beta1,DE1,T,M->mass];
	QII = QMorseExact[Beta2,DE2,T,M->mass];
	kcal2au*(MorseEnergy[nu,Beta2,DE2,M->mass]-MorseEnergy[mu,Beta1,DE1,M->mass])+kb*T*(Log[QII/QI]) + eps - eta*ev2au
];

(*
eta - overpotential in volts;
eps - energy in Hartrees (used in Fermi distribution);
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mu - quantum number of oscillator on the left;
nu - quantum number of oscillator on the right;
Returns the reaction free energy for a harmonic approximation in Hartrees;
*)

dGHarmonic[eta_,eps_,omega1_,omega2_,T_,mu_,nu_]:=Module[{QI,QII},
	QI = 1/(Exp[(omega1*cm2au)/(2*kb*T)]-Exp[-(omega1*cm2au)/(2*kb*T)]);
	QII = 1/(Exp[(omega2*cm2au)/(2*kb*T)]-Exp[-(omega2*cm2au)/(2*kb*T)]);
	omega2*(nu+1/2)*cm2au - omega1*(mu+1/2)*cm2au + kb*T*Log[QII/QI] + eps - eta*ev2au
];

(*
Heterogeneous rate constants
*)

(*
Anodic (oxidation) rate constants at a fixed distance (sec^-1) - High temperature approximation for the donor-acceptor mode
*)

AnodicRateConstantMorse[Vel_,eta_,lambda_,MR_,dR_,omega_,Beta1_,DE1_,Beta2_,DE2_,mumax_,numax_,d_,T_,OptionsPattern[{M->MassH,Step->0.0001}]]:=
Module[{eps,s,alpha,lamalph,x,h,S0,S2,S1,lam,SS},
	Array[alpha,{mumax+1,numax+1},0];
	Array[lamalph,{mumax+1,numax+1},0];
	h=OptionValue[Step];
	Array[S0,{mumax+1,numax+1},0];
	Array[S2,{mumax+1,numax+1},0];
	Array[S1,{mumax+1,numax+1},0];
	Do[S0[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[S2[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d+h,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[S1[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d-h,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=-(1/S0[mu,nu])*(S2[mu,nu]-S1[mu,nu])/(2*h),{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	lam=kcal2au*(lambda+LambdaR[MR,omega,dR]);
	SS=Sum[
	    BoltzmannMorse[Beta1,DE1,mu,T,M->OptionValue[M]]*S0[mu,nu]^2*
		Exp[4*lamalph[mu,nu]*kb*T/(hbar^2*(omega*cm2au)^2)]*
		Sqrt[Pi/(kb*T*(lam+lamalph[mu,nu]))]*
		NIntegrate[
			(1-Fermi[eps,T])*
			Exp[-(dGMorse[eta,eps,Beta1,DE1,Beta2,DE2,T,mu,nu,M->OptionValue[M]] + lam + lamalph[mu,nu] + 2*alpha[mu,nu]*dR*kb*T)^2/
				(4*kb*T*(lam + lamalph[mu,nu]))],
	    {eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}];
	10^12*(Vel*kcal2au)^2*SS/au2ps
];

AnodicRateConstantHarmonic[Vel_,eta_,lambda_,MR_,dR_,omega_,omega1_,omega2_,mumax_,numax_,d_,T_,OptionsPattern[{M->MassH}]]:=
Module[{eps,s,alpha,lamalph,x,h,ov,lam,SS},
	Array[alpha,{mumax+1,numax+1},0];
	Array[lamalph,{mumax+1,numax+1},0];
	Array[ov,{mumax+1,numax+1},0];
	Do[ov[mu,nu]=HarmonicOverlap[mu,nu,omega1,omega2,d,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=-(1/ov[mu,nu])*ND[HarmonicOverlap[mu,nu,omega1,omega2,x,M->OptionValue[M]],x,d],{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	lam=kcal2au*(lambda+LambdaR[MR,omega,dR]);
	SS=Sum[
		BoltzmannHarmonic[omega1,mu,T]*ov[mu,nu]^2*
		Exp[4*lamalph[mu,nu]*kb*T/(hbar^2*(omega*cm2au)^2)]*
		Sqrt[Pi/(kb*T*(lam+lamalph[mu,nu]))]*
		NIntegrate[
			(1-Fermi[eps,T])*
			Exp[-(dGHarmonic[eta,eps,omega1,omega2,T,mu,nu] + lam + lamalph[mu,nu] + 2*alpha[mu,nu]*dR*kb*T)^2/
			(4*kb*T*(lam+lamalph[mu,nu]))],
		{eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}];
	10^12*(Vel*kcal2au)^2*SS/au2ps
];

AnodicRateConstantFGH[fghlist_,Vel_,eta_,lambda_,MR_,dR_,omega_,mumax_,numax_,T_]:=
Module[{eps,alpha,lamalph,ov,lam,energymu,energynu,Z1,Z2,SS,Ns},
        Ns = Length[fghlist[[1,1]]];
        If[mumax+1>Ns||numax+1>Ns,Message[AnodicRateConstantFGH::numberofstates];Abort[]];
	Array[alpha,{mumax+1,numax+1},0];
	Array[lamalph,{mumax+1,numax+1},0];
	Array[ov,{mumax+1,numax+1},0];
	Array[energymu,Ns,0];
	Array[energynu,Ns,0];
	Do[ov[mu,nu]=fghlist[[3]][[mu+1,nu+1]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=fghlist[[4]][[mu+1,nu+1]],{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	Do[energymu[mu]=kcal2au*fghlist[[1]][[1,mu+1]],{mu,0,Ns-1}];
	Do[energynu[nu]=kcal2au*fghlist[[2]][[1,nu+1]],{nu,0,Ns-1}];
	Z1 = Sum[Exp[-energymu[mu]/(kb*T)],{mu,0,Ns-1}];
	Z2 = Sum[Exp[-energynu[nu]/(kb*T)],{nu,0,Ns-1}];
	lam = kcal2au*(lambda + LambdaR[MR,omega,dR]);
	SS=Sum[
		(Exp[-energymu[mu]/(kb*T)]/Z1)*
		ov[mu,nu]^2*
		Exp[4*lamalph[mu,nu]*kb*T/(hbar^2*(omega*cm2au)^2)]*
		Sqrt[Pi/(kb*T*(lam+lamalph[mu,nu]))]*
		NIntegrate[
			(1-Fermi[eps,T])*
			Exp[-(energynu[nu]-energymu[mu] + kb*T*Log[Z2/Z1] + eps - eta*ev2au + lam + lamalph[mu,nu] + 2*alpha[mu,nu]*dR*kb*T)^2/
			(4*kb*T*(lam+lamalph[mu,nu]))],
		{eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}];
	10^12*(Vel*kcal2au)^2*SS/au2ps
];

(*
Cathodic (reduction) rate constants at a fixed distance (sec^-1) - High temperature approximation for the donor-acceptor mode
*)

CathodicRateConstantMorse[Vel_,eta_,lambda_,MR_,dR_,omega_,Beta1_,DE1_,Beta2_,DE2_,mumax_,numax_,d_,T_,OptionsPattern[{M->MassH,Step->0.0001}]]:=
Module[{eps,s,alpha,lamalph,x,h,S0,S2,S1,lam,SS},
	Array[alpha,{mumax+1,numax+1},{0,0}];
	Array[lamalph,{mumax+1,numax+1},{0,0}];
	h=OptionValue[Step];
	Array[S0,{mumax+1,numax+1},{0,0}];
	Array[S2,{mumax+1,numax+1},{0,0}];
	Array[S1,{mumax+1,numax+1},{0,0}];
	Do[S0[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[S2[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d+h,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[S1[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d-h,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=-(1/S0[mu,nu])*(S2[mu,nu]-S1[mu,nu])/(2*h),{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	lam=kcal2au*(lambda+LambdaR[MR,omega,dR]);
	SS=Sum[
	    BoltzmannMorse[Beta2,DE2,nu,T,M->OptionValue[M]]*S0[mu,nu]^2*
		Exp[-2*alpha[mu,nu]*dR]*
		Exp[4*lamalph[mu,nu]*kb*T/(hbar^2*(omega*cm2au)^2)]*
		Sqrt[Pi/(kb*T*(lam+lamalph[mu,nu]))]*
		NIntegrate[
			Fermi[eps,T]*
			Exp[-(-dGMorse[eta,eps,Beta1,DE1,Beta2,DE2,T,mu,nu,M->OptionValue[M]] + lam + lamalph[mu,nu] - 2*alpha[mu,nu]*dR*kb*T)^2/
			(4*kb*T*(lam + lamalph[mu,nu]))],
	    {eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}];
	10^12*(Vel*kcal2au)^2*SS/au2ps
];

CathodicRateConstantHarmonic[Vel_,eta_,lambda_,MR_,dR_,omega_,omega1_,omega2_,mumax_,numax_,d_,T_,OptionsPattern[{M->MassH}]]:=
Module[{eps,s,alpha,lamalph,x,h,ov,lam,SS},
	Array[alpha,{mumax+1,numax+1},0];
	Array[lamalph,{mumax+1,numax+1},0];
	Array[ov,{mumax+1,numax+1},0];
	Do[ov[mu,nu]=HarmonicOverlap[mu,nu,omega1,omega2,d,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=-(1/ov[mu,nu])*ND[HarmonicOverlap[mu,nu,omega1,omega2,x,M->OptionValue[M]],x,d],{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	lam=kcal2au*(lambda+LambdaR[MR,omega,dR]);
	SS=Sum[
		BoltzmannHarmonic[omega1,nu,T]*ov[mu,nu]^2*
		Exp[-2*alpha[mu,nu]*dR]*
		Exp[4*lamalph[mu,nu]*kb*T/(hbar^2*(omega*cm2au)^2)]*
		Sqrt[Pi/(kb*T*(lam+lamalph[mu,nu]))]*
		NIntegrate[
			(Fermi[eps,T])*
			Exp[-(-dGHarmonic[eta,eps,omega1,omega2,T,mu,nu]+lam+lamalph[mu,nu]-2*alpha[mu,nu]*dR*kb*T)^2/
			(4*kb*T*(lam+lamalph[mu,nu]))],
		{eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}];
	10^12*(Vel*kcal2au)^2*SS/au2ps
];

CathodicRateConstantFGH[fghlist_,Vel_,eta_,lambda_,MR_,dR_,omega_,mumax_,numax_,T_]:=
Module[{eps,alpha,lamalph,ov,lam,energymu,energynu,Z1,Z2,SS,Ns},
        Ns = Length[fghlist[[1,1]]];
        If[mumax+1>Ns||numax+1>Ns,Message[CathodicRateConstantFGH::numberofstates];Abort[]];
	Array[alpha,{mumax+1,numax+1},0];
	Array[lamalph,{mumax+1,numax+1},0];
	Array[ov,{mumax+1,numax+1},0];
	Array[energymu,Ns,0];
	Array[energynu,Ns,0];
	Do[ov[mu,nu]=fghlist[[3]][[mu+1,nu+1]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=fghlist[[4]][[mu+1,nu+1]],{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	Do[energymu[mu]=kcal2au*fghlist[[1]][[1,mu+1]],{mu,0,Ns-1}];
	Do[energynu[nu]=kcal2au*fghlist[[2]][[1,nu+1]],{nu,0,Ns-1}];
	Z1 = Sum[Exp[-energymu[mu]/(kb*T)],{mu,0,Ns-1}];
	Z2 = Sum[Exp[-energynu[nu]/(kb*T)],{nu,0,Ns-1}];
	lam=kcal2au*(lambda+LambdaR[MR,omega,dR]);
	SS=Sum[
		(Exp[-energynu[nu]/(kb*T)]/Z2)*
		ov[mu,nu]^2*
		Exp[-2*alpha[mu,nu]*dR]*
		Exp[4*lamalph[mu,nu]*kb*T/(hbar^2*(omega*cm2au)^2)]*
		Sqrt[Pi/(kb*T*(lam+lamalph[mu,nu]))]*
		NIntegrate[
			Fermi[eps,T]*
			Exp[-(-energynu[nu] + energymu[mu] - kb*T*Log[Z2/Z1] - eps + eta*ev2au + lam + lamalph[mu,nu] - 2*alpha[mu,nu]*dR*kb*T)^2/
			(4*kb*T*(lam+lamalph[mu,nu]))],
		{eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}];
	10^12*(Vel*kcal2au)^2*SS/au2ps
];

(*
Transfer coefficients - High temperature approximation for the donor-acceptor mode
*)

AnodicTransferCoefficientMorse[Vel_,eta_,lambda_,MR_,dR_,omega_,Beta1_,DE1_,Beta2_,DE2_,mumax_,numax_,d_,T_,OptionsPattern[{Dstep->0.001,M->MassH,Step->0.0001}]]:=
Module[{Mass,St,jleft,jright,\[Delta]},
	Mass=OptionValue[M];
	St=OptionValue[Step];
	\[Delta]=OptionValue[Dstep];
	\[Eta]=eta;
	jleft=Log[AnodicRateConstantMorse[Vel,\[Eta]-\[Delta],lambda,MR,dR,omega,Beta1,DE1,Beta2,DE2,mumax,numax,d,T,M->Mass,Step->St]];
	jright=Log[AnodicRateConstantMorse[Vel,\[Eta]+\[Delta],lambda,MR,dR,omega,Beta1,DE1,Beta2,DE2,mumax,numax,d,T,M->Mass,Step->St]];
	-kb*T*(jright-jleft)/(2*\[Delta]*ev2au)
];

CathodicTransferCoefficientMorse[Vel_,eta_,lambda_,MR_,dR_,omega_,Beta1_,DE1_,Beta2_,DE2_,mumax_,numax_,d_,T_,OptionsPattern[{Dstep->0.001,M->MassH,Step->0.0001}]]:=
Module[{Mass,St,jleft,jright,\[Delta]},
	Mass=OptionValue[M];
	St=OptionValue[Step];
	\[Delta]=OptionValue[Dstep];
	\[Eta]=eta;
	jleft=Log[CathodicRateConstantMorse[Vel,\[Eta]-\[Delta],lambda,MR,dR,omega,Beta1,DE1,Beta2,DE2,mumax,numax,d,T,M->Mass,Step->St]];
	jright=Log[CathodicRateConstantMorse[Vel,\[Eta]+\[Delta],lambda,MR,dR,omega,Beta1,DE1,Beta2,DE2,mumax,numax,d,T,M->Mass,Step->St]];
	-kb*T*(jright-jleft)/(2*\[Delta]*ev2au)
];

AnodicTransferCoefficientHarmonic[Vel_,eta_,lambda_,MR_,dR_,omega_,omega1_,omega2_,mumax_,numax_,d_,T_,OptionsPattern[{M->MassH,Dstep->0.001}]]:=
Module[{Mass,jleft,jright,\[Delta]},
	Mass=OptionValue[M];
	\[Delta]=OptionValue[Dstep];
	\[Eta]=eta;
	jleft=Log[AnodicRateConstantHarmonic[Vel,\[Eta]-\[Delta],lambda,MR,dR,omega,omega1,omega2,mumax,numax,d,T,M->Mass]];
	jright=Log[AnodicRateConstantHarmonic[Vel,\[Eta]+\[Delta],lambda,MR,dR,omega,omega1,omega2,mumax,numax,d,T,M->Mass]];
	-kb*T*(jright-jleft)/(2*\[Delta]*ev2au)
];

CathodicTransferCoefficientHarmonic[Vel_,eta_,lambda_,MR_,dR_,omega_,omega1_,omega2_,mumax_,numax_,d_,T_,OptionsPattern[{M->MassH,Dstep->0.001}]]:=
Module[{Mass,jleft,jright,\[Delta]},
	Mass=OptionValue[M];
	\[Delta]=OptionValue[Dstep];
	\[Eta]=eta;
	jleft=Log[CathodicRateConstantHarmonic[Vel,\[Eta]-\[Delta],lambda,MR,dR,omega,omega1,omega2,mumax,numax,d,T,M->Mass]];
	jright=Log[CathodicRateConstantHarmonic[Vel,\[Eta]+\[Delta],lambda,MR,dR,omega,omega1,omega2,mumax,numax,d,T,M->Mass]];
	-kb*T*(jright-jleft)/(2*\[Delta]*ev2au)
];

AnodicTransferCoefficientFGH[fghlist_,Vel_,eta_,lambda_,MR_,dR_,omega_,mumax_,numax_,T_,OptionsPattern[{Dstep->0.001}]]:=
Module[{jleft,jright,\[Delta]},
	\[Delta]=OptionValue[Dstep];
	\[Eta]=eta;
	jleft=Log[AnodicRateConstantFGH[fghlist,Vel,\[Eta]-\[Delta],lambda,MR,dR,omega,mumax,numax,T]];
	jright=Log[AnodicRateConstantFGH[fghlist,Vel,\[Eta]+\[Delta],lambda,MR,dR,omega,mumax,numax,T]];
	-kb*T*(jright-jleft)/(2*\[Delta]*ev2au)
];

CathodicTransferCoefficientFGH[fghlist_,Vel_,eta_,lambda_,MR_,dR_,omega_,mumax_,numax_,T_,OptionsPattern[{Dstep->0.001}]]:=
	Module[{jleft,jright,\[Delta]},
	\[Delta]=OptionValue[Dstep];
	\[Eta]=eta;
	jleft=Log[CathodicRateConstantFGH[fghlist,Vel,\[Eta]-\[Delta],lambda,MR,dR,omega,mumax,numax,T]];
	jright=Log[CathodicRateConstantFGH[fghlist,Vel,\[Eta]+\[Delta],lambda,MR,dR,omega,mumax,numax,T]];
	-kb*T*(jright-jleft)/(2*\[Delta]*ev2au)
];

(*
Current Densities - High Temperature Approximation
*)

(*
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Returns the current density for the Morse approximation - Oxidation or Reduction;
*)

AnodicCurrentDensityMorse[eta_,lambda_,MR_,dR_,omega_,Beta1_,DE1_,Beta2_,DE2_,mumax_,numax_,d_,T_,OptionsPattern[{M->MassH,Step->0.0001}]]:=
Module[{s,alpha,lamalph,x,h,S0,S2,S1,lam},
	Array[alpha,{mumax+1,numax+1},0];
	Array[lamalph,{mumax+1,numax+1},0];
	h=OptionValue[Step];
	Array[S0,{mumax+1,numax+1},0];
	Array[S2,{mumax+1,numax+1},0];
	Array[S1,{mumax+1,numax+1},0];
	Do[S0[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[S2[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d+h,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[S1[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d-h,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=-(1/S0[mu,nu])*(S2[mu,nu]-S1[mu,nu])/(2*h),{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	lam=kcal2au*(lambda+LambdaR[MR,omega,dR]);
	Sum[
	    BoltzmannMorse[Beta1,DE1,mu,T,M->OptionValue[M]]*S0[mu,nu]^2*
		Exp[4*lamalph[mu,nu]*kb*T/(hbar^2*(omega*cm2au)^2)]*
		Sqrt[Pi/(kb*T*(lam+lamalph[mu,nu]))]*
		NIntegrate[
			(1-Fermi[eps,T])*
			Exp[-(dGMorse[eta,eps,Beta1,DE1,Beta2,DE2,T,mu,nu,M->OptionValue[M]] + lam + lamalph[mu,nu] + 2*alpha[mu,nu]*dR*kb*T)^2/
				(4*kb*T*(lam + lamalph[mu,nu]))],
	    {eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}]
];

CathodicCurrentDensityMorse[eta_,lambda_,MR_,dR_,omega_,Beta1_,DE1_,Beta2_,DE2_,mumax_,numax_,d_,T_,OptionsPattern[{M->MassH,Step->0.0001}]]:=
Module[{eps,s,alpha,lamalph,x,h,S0,S2,S1,lam},
	Array[alpha,{mumax+1,numax+1},{0,0}];
	Array[lamalph,{mumax+1,numax+1},{0,0}];
	h=OptionValue[Step];
	Array[S0,{mumax+1,numax+1},{0,0}];
	Array[S2,{mumax+1,numax+1},{0,0}];
	Array[S1,{mumax+1,numax+1},{0,0}];
	Do[S0[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[S2[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d+h,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[S1[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d-h,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=-(1/S0[mu,nu])*(S2[mu,nu]-S1[mu,nu])/(2*h),{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	lam=kcal2au*(lambda+LambdaR[MR,omega,dR]);
	Sum[
	    BoltzmannMorse[Beta2,DE2,nu,T,M->OptionValue[M]]*S0[mu,nu]^2*
		Exp[-2*alpha[mu,nu]*dR]*
		Exp[4*lamalph[mu,nu]*kb*T/(hbar^2*(omega*cm2au)^2)]*
		Sqrt[Pi/(kb*T*(lam+lamalph[mu,nu]))]*
		NIntegrate[
			Fermi[eps,T]*
			Exp[-(-dGMorse[eta,eps,Beta1,DE1,Beta2,DE2,T,mu,nu,M->OptionValue[M]] + lam + lamalph[mu,nu] - 2*alpha[mu,nu]*dR*kb*T)^2/
			(4*kb*T*(lam + lamalph[mu,nu]))],
	    {eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}]
];

(*
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Cox - concentration of the oxidized complex (any consistent units);
Cred - concentration of the reduced complex (any consistent units);
Returns the total current density for the Morse approximation;
*)

TotalCurrentDensityMorse[eta_,lambda_,MR_,dR_,omega_,Beta1_,DE1_,Beta2_,DE2_,mumax_,numax_,d_,T_,Cox_,Cred_,OptionsPattern[{M->MassH,Step->0.0001}]]:=
Module[{Conc},
	Conc=Cred*Exp[(eta*ev2au)/(kb*T)]/Cox-1;
	Cox*CathodicCurrentDensityMorse[eta,lambda,MR,dR,omega,Beta1,DE1,Beta2,DE2,mumax,numax,d,T,M->OptionValue[M],Step->OptionValue[Step]]*Conc
];

(*
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Returns the current density for the harmonic approximation - Reduction or Oxidation;
*)

AnodicCurrentDensityHarmonic[eta_,lambda_,MR_,dR_,omega_,omega1_,omega2_,mumax_,numax_,d_,T_,OptionsPattern[{M->MassH}]]:=
Module[{eps,s,alpha,lamalph,x,h,ov,lam},
	Array[alpha,{mumax+1,numax+1},0];
	Array[lamalph,{mumax+1,numax+1},0];
	Array[ov,{mumax+1,numax+1},0];
	Do[ov[mu,nu]=HarmonicOverlap[mu,nu,omega1,omega2,d,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=-(1/ov[mu,nu])*ND[HarmonicOverlap[mu,nu,omega1,omega2,x,M->OptionValue[M]],x,d],{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	lam=kcal2au*(lambda+LambdaR[MR,omega,dR]);
	Sum[
		BoltzmannHarmonic[omega1,mu,T]*ov[mu,nu]^2*
		Exp[4*lamalph[mu,nu]*kb*T/(hbar^2*(omega*cm2au)^2)]*
		Sqrt[Pi/(kb*T*(lam+lamalph[mu,nu]))]*
		NIntegrate[
			(1-Fermi[eps,T])*
			Exp[-(dGHarmonic[eta,eps,omega1,omega2,T,mu,nu] + lam + lamalph[mu,nu] + 2*alpha[mu,nu]*dR*kb*T)^2/
			(4*kb*T*(lam+lamalph[mu,nu]))],
		{eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}]
];

CathodicCurrentDensityHarmonic[eta_,lambda_,MR_,dR_,omega_,omega1_,omega2_,mumax_,numax_,d_,T_,OptionsPattern[{M->MassH,Step->0.0001}]]:=
Module[{eps,s,alpha,lamalph,x,h,ov,lam},
	Array[alpha,{mumax+1,numax+1},0];
	Array[lamalph,{mumax+1,numax+1},0];
	Array[ov,{mumax+1,numax+1},0];
	Do[ov[mu,nu]=HarmonicOverlap[mu,nu,omega1,omega2,d,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=-(1/ov[mu,nu])*ND[HarmonicOverlap[mu,nu,omega1,omega2,x,M->OptionValue[M]],x,d],{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	lam=kcal2au*(lambda+LambdaR[MR,omega,dR]);
	Sum[
		BoltzmannHarmonic[omega1,nu,T]*ov[mu,nu]^2*
		Exp[-2*alpha[mu,nu]*dR]*
		Exp[4*lamalph[mu,nu]*kb*T/(hbar^2*(omega*cm2au)^2)]*
		Sqrt[Pi/(kb*T*(lam+lamalph[mu,nu]))]*
		NIntegrate[
			Fermi[eps,T]*
			Exp[-(-dGHarmonic[eta,eps,omega1,omega2,T,mu,nu]+lam+lamalph[mu,nu]-2*alpha[mu,nu]*dR*kb*T)^2/
			(4*kb*T*(lam+lamalph[mu,nu]))],
		{eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}]
];

(*
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Cox - concentration of the oxidized complex (any consistent units);
Cred - concentration of the reduced complex (any consistent units);
Returns the total current density for the harmonic approximation;
*)

TotalCurrentDensityHarmonic[eta_,lambda_,MR_,dR_,omega_,omega1_,omega2_,mumax_,numax_,d_,T_,Cox_,Cred_,OptionsPattern[{M->MassH}]]:=Module[{Conc},
	Conc=Cred*Exp[(eta*ev2au)/(kb*T)]/Cox-1;
	Cox*CathodicCurrentDensityHarmonic[eta,lambda,MR,dR,omega,omega1,omega2,mumax,numax,d,T]*Conc
];

(*
Current Densities - Low Temperature
*)

(*
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Returns the current density for the Morse approximation - Oxidation or Reduction;
*)

AnodicCurrentDensityLowTMorse[eta_,lambda_,MR_,dR_,omega_,Beta1_,DE1_,Beta2_,DE2_,mumax_,numax_,d_,T_,OptionsPattern[{M->MassH,Step->0.0001}]]:=
Module[{eps,s,alpha,lamalph,x,h,S0,S2,S1,lam},
	Array[alpha,{mumax+1,numax+1},0];
	Array[lamalph,{mumax+1,numax+1},0];
	h=OptionValue[Step];
	Array[S0,{mumax+1,numax+1},0];
	Array[S2,{mumax+1,numax+1},0];
	Array[S1,{mumax+1,numax+1},0];
	Do[S0[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[S2[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d+h,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[S1[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d-h,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=-(1/S0[mu,nu])*(S2[mu,nu]-S1[mu,nu])/(2*h),{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	lam =kcal2au*lambda;
	Sum[
		BoltzmannMorse[Beta1,DE1,mu,T,M->OptionValue[M]]*
		S0[mu,nu]^2*
		Exp[-alpha[mu,nu]*dR]*
		Exp[(lamalph[mu,nu]-kcal2au*LambdaR[MR,omega,dR])/(hbar*omega*cm2au)]*
		Sqrt[Pi/(kb*T*lam)]*
		NIntegrate[
			(1-Fermi[eps,T])*
			Exp[-(dGMorse[eta,eps,Beta1,DE1,Beta2,DE2,T,mu,nu,M->OptionValue[M]] + lam)^2/(4*kb*T*lam)],
		{eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}]
];

CathodicCurrentDensityLowTMorse[eta_,lambda_,MR_,dR_,omega_,Beta1_,DE1_,Beta2_,DE2_,mumax_,numax_,d_,T_,OptionsPattern[{M->MassH,Step->0.0001}]]:=
Module[{eps,s,alpha,lamalph,x,h,S0,S2,S1,lam},
	Array[alpha,{mumax+1,numax+1},{0,0}];
	Array[lamalph,{mumax+1,numax+1},{0,0}];
	h=OptionValue[Step];
	Array[S0,{mumax+1,numax+1},{0,0}];
	Array[S2,{mumax+1,numax+1},{0,0}];
	Array[S1,{mumax+1,numax+1},{0,0}];
	Do[S0[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[S2[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d+h,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[S1[mu,nu]=NMorseOverlap[DE1,Beta1,DE2,Beta2,mu,nu,d-h,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=-(1/S0[mu,nu])*(S2[mu,nu]-S1[mu,nu])/(2*h),{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	lam =kcal2au*lambda;
	Sum[
		BoltzmannMorse[Beta2,DE2,nu,T,M->OptionValue[M]]*
		S0[mu,nu]^2*
		Exp[-alpha[mu,nu]*dR]*
		Exp[(lamalph[mu,nu]-kcal2au*LambdaR[MR,omega,dR])/(hbar*omega*cm2au)]*
		Sqrt[Pi/(kb*T*lam)]*
		NIntegrate[
			Fermi[eps,T]*
			Exp[-(-dGMorse[eta,eps,Beta1,DE1,Beta2,DE2,T,mu,nu,M->OptionValue[M]]+lam)^2/(4*kb*T*lam)],
		{eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}]
];

(*
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
omega - frequency of the coupling between the oscillators in wavenumbers;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
Beta1 - beta parameter for oscillator on the left in inverse Bohr;
DE1 - dissociation energy of the oscillator on the left in kcal/mol;
Beta2 - beta parameter for oscillator on the right in inverse Bohr;
DE2 - dissociation energy of the oscillator on the right in kcal/mol;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Cox - concentration of the oxidized complex (any consistent units);
Cred - concentration of the reduced complex (any consistent units);
Returns the total current density for the Morse approximation;
*)

TotalCurrentDensityLowTMorse[eta_,lambda_,MR_,dR_,omega_,Beta1_,DE1_,Beta2_,DE2_,mumax_,numax_,d_,T_,Cox_,Cred_,OptionsPattern[{M->MassH,Step->0.0001}]]:=
Module[{Conc},
	Conc=Cred*Exp[(eta*ev2au)/(kb*T)]/Cox-1;
	Cox*CathodicCurrentDensityLowTMorse[eta,lambda,MR,dR,omega,Beta1,DE1,Beta2,DE2,mumax,numax,d,T,M->OptionValue[M],Step->OptionValue[Step]]*Conc
];

(*
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Returns the current density for the harmonic approximation - Reduction or Oxidation;
*)

AnodicCurrentDensityLowTHarmonic[eta_,lambda_,MR_,dR_,omega_,omega1_,omega2_,mumax_,numax_,d_,T_,OptionsPattern[{M->MassH}]]:=
Module[{eps,s,alpha,lamalph,x,h,ov,lam},
	Array[alpha,{mumax+1,numax+1},0];
	Array[lamalph,{mumax+1,numax+1},0];
	Array[ov,{mumax+1,numax+1},0];
	Do[ov[mu,nu]=HarmonicOverlap[mu,nu,omega1,omega2,d,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=-(1/ov[mu,nu])*ND[HarmonicOverlap[mu,nu,omega1,omega2,x,M->OptionValue[M]],x,d],{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	lam =kcal2au*lambda;
	Sum[
		BoltzmannHarmonic[omega1,mu,T]*
		ov[mu,nu]^2*
		Exp[-alpha[mu,nu]*dR]*
		Exp[(lamalph[mu,nu]-kcal2au*LambdaR[MR,omega,dR])/(hbar*omega*cm2au)]*
		Sqrt[Pi/(kb*T*lam)]*
		NIntegrate[
			(1-Fermi[eps,T])*
			Exp[-(dGHarmonic[eta,eps,omega1,omega2,T,mu,nu]+lam)^2/(4*kb*T*lam)],
		{eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}]
];

CathodicCurrentDensityLowTHarmonic[eta_,lambda_,MR_,dR_,omega_,omega1_,omega2_,mumax_,numax_,d_,T_,OptionsPattern[{M->MassH,Step->0.0001}]]:=
Module[{eps,s,alpha,lamalph,x,h,ov,lam},
	Array[alpha,{mumax+1,numax+1},0];
	Array[lamalph,{mumax+1,numax+1},0];
	Array[ov,{mumax+1,numax+1},0];
	Do[ov[mu,nu]=HarmonicOverlap[mu,nu,omega1,omega2,d,M->OptionValue[M]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=-(1/ov[mu,nu])*ND[HarmonicOverlap[mu,nu,omega1,omega2,x,M->OptionValue[M]],x,d],{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	lam=kcal2au*lambda;
	Sum[
		BoltzmannHarmonic[omega1,nu,T]*
		ov[mu,nu]^2*
		Exp[-alpha[mu,nu]*dR]*
		Exp[(lamalph[mu,nu]-kcal2au*LambdaR[MR,omega,dR])/(hbar*omega*cm2au)]*
		Sqrt[Pi/(kb*T*lam)]*
		NIntegrate[
			Fermi[eps,T]*
			Exp[-(-dGHarmonic[eta,eps,omega1,omega2,T,mu,nu]+lam)^2/(4*kb*T*lam)],
		{eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}]
];

(*
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
omega1 - frequency of oscillator on the left in wavenumbers;
omega2 - frequency of oscillator on the right in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
d - distance between the minima of the displaced potentials in bohr;
T - temperature in Kelvin;
Cox - concentration of the oxidized complex (any consistent units);
Cred - concentration of the reduced complex (any consistent units);
Returns the total current density for the harmonic approximation;
*)

TotalCurrentDensityLowTHarmonic[eta_,lambda_,MR_,dR_,omega_,omega1_,omega2_,mumax_,numax_,d_,T_,Cox_,Cred_,OptionsPattern[{M->MassH}]]:=
Module[{Conc},
	Conc=Cred*Exp[(eta*ev2au)/(kb*T)]/Cox-1;
	Cox*CathodicCurrentDensityLowTHarmonic[eta,lambda,MR,dR,omega,omega1,omega2,mumax,numax,d,T]*Conc
];

(*
Current Densities - Alpha and overlaps from FGH calculations for arbitrary potentials
*)

(*
fghlist - nested list with energies, wavefunctions, overlap integrals and alphas;
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
T - temperature in Kelvin;
Returns the current density for general proton potentials - Reduction or Oxidation;
*)

(* we are here *)

AnodicCurrentDensityFGH[fghlist_,eta_,lambda_,MR_,dR_,omega_,mumax_,numax_,T_]:=
Module[{eps,alpha,lamalph,ov,lam,energymu,energynu,Z1,Z2,Ns},
        Ns = Length[fghlist[[1,1]]];
        If[mumax+1>Ns||numax+1>Ns,Message[AnodicCurrentDensityFGH::numberofstates];Abort[]];
	Array[alpha,{mumax+1,numax+1},0];
	Array[lamalph,{mumax+1,numax+1},0];
	Array[ov,{mumax+1,numax+1},0];
	Array[energymu,Ns,0];
	Array[energynu,Ns,0];
	Do[ov[mu,nu]=fghlist[[3]][[mu+1,nu+1]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=fghlist[[4]][[mu+1,nu+1]],{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	Do[energymu[mu]=kcal2au*fghlist[[1]][[1,mu+1]],{mu,0,Ns-1}];
	Do[energynu[nu]=kcal2au*fghlist[[2]][[1,nu+1]],{nu,0,Ns-1}];
	Z1 = Sum[Exp[-energymu[mu]/(kb*T)],{mu,0,Ns-1}];
	Z2 = Sum[Exp[-energynu[nu]/(kb*T)],{nu,0,Ns-1}];
	lam = kcal2au*(lambda + LambdaR[MR,omega,dR]);
	Sum[
		(Exp[-energymu[mu]/(kb*T)]/Z1)*
		ov[mu,nu]^2*
		Exp[4*lamalph[mu,nu]*kb*T/(hbar^2*(omega*cm2au)^2)]*
		Sqrt[Pi/(kb*T*(lam+lamalph[mu,nu]))]*
		NIntegrate[
			(1-Fermi[eps,T])*
			Exp[-(energynu[nu]-energymu[mu] + kb*T*Log[Z2/Z1] + eps - eta*ev2au + lam + lamalph[mu,nu] + 2*alpha[mu,nu]*dR*kb*T)^2/
			(4*kb*T*(lam+lamalph[mu,nu]))],
		{eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}]
];

CathodicCurrentDensityFGH[fghlist_,eta_,lambda_,MR_,dR_,omega_,mumax_,numax_,T_]:=
Module[{eps,alpha,lamalph,ov,lam,energymu,energynu,Z1,Z2,Ns},
        Ns = Length[fghlist[[1,1]]];
        If[mumax+1>Ns||numax+1>Ns,Message[CathodicCurrentDensityFGH::numberofstates];Abort[]];
	Array[alpha,{mumax+1,numax+1},0];
	Array[lamalph,{mumax+1,numax+1},0];
	Array[ov,{mumax+1,numax+1},0];
	Array[energymu,Ns,0];
	Array[energynu,Ns,0];
	Do[ov[mu,nu]=fghlist[[3]][[mu+1,nu+1]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=fghlist[[4]][[mu+1,nu+1]],{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	Do[energymu[mu]=kcal2au*fghlist[[1]][[1,mu+1]],{mu,0,Ns-1}];
	Do[energynu[nu]=kcal2au*fghlist[[2]][[1,nu+1]],{nu,0,Ns-1}];
	Z1=Sum[Exp[-energymu[mu]/(kb*T)],{mu,0,Ns-1}];
	Z2=Sum[Exp[-energynu[nu]/(kb*T)],{nu,0,Ns-1}];
	lam=kcal2au*(lambda+LambdaR[MR,omega,dR]);
	Sum[
		(Exp[-(energynu[nu])/(kb*T)]/Z2)*
		ov[mu,nu]^2*
		Exp[-2*alpha[mu,nu]*dR]*
		Exp[4*lamalph[mu,nu]*kb*T/(hbar^2*(omega*cm2au)^2)]*
		Sqrt[Pi/(kb*T*(lam+lamalph[mu,nu]))]*
		NIntegrate[
			Fermi[eps,T]*
			Exp[-(-energynu[nu] + energymu[mu] - kb*T*Log[Z2/Z1] - eps + eta*ev2au + lam + lamalph[mu,nu] - 2*alpha[mu,nu]*dR*kb*T)^2/
			(4*kb*T*(lam+lamalph[mu,nu]))],
		{eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}]
];

(*
fghlist - nested list with energies, wavefunctions, overlap integrals and alphas;
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
T - temperature in Kelvin;
Cox - concentration of the oxidized complex (any consistent units);
Cred - concentration of the reduced complex (any consistent units);
Returns the total current density for general proton potentials;
*)

TotalCurrentDensityFGH[fghlist_,eta_,lambda_,MR_,dR_,omega_,mumax_,numax_,T_,Cox_,Cred_]:=
Module[{Conc},
	Conc=Cred*Exp[(eta*ev2au)/(kb*T)]/Cox-1;
	Cox*CathodicCurrentDensityFGH[fghlist,eta,lambda,MR,dR,omega,mumax,numax,T]*Conc
];

(*
fghlist - nested list with energies, wavefunctions, overlap integrals and alphas;
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
T - temperature in Kelvin;
Returns the current density for general proton potentials in the low temperature approximation - Reduction or Oxidation;
*)

AnodicCurrentDensityLowTFGH[fghlist_,eta_,lambda_,MR_,dR_,omega_,mumax_,numax_,T_]:=
Module[{eps,alpha,lamalph,ov,lam,energymu,energynu,Z1,Z2,Ns},
        Ns = Length[fghlist[[1,1]]];
        If[mumax+1>Ns||numax+1>Ns,Message[AnodicCurrentDensityLowTFGH::numberofstates];Abort[]];
	Array[alpha,{mumax+1,numax+1},0];
	Array[lamalph,{mumax+1,numax+1},0];
	Array[ov,{mumax+1,numax+1},0];
	Array[energymu,Ns,0];
	Array[energynu,Ns,0];
	Do[ov[mu,nu]=fghlist[[3]][[mu+1,nu+1]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=fghlist[[4]][[mu+1,nu+1]],{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	Do[energymu[mu]=kcal2au*fghlist[[1]][[1,mu+1]],{mu,0,Ns-1}];
	Do[energynu[nu]=kcal2au*fghlist[[2]][[1,nu+1]],{nu,0,Ns-1}];
	Z1 = Sum[Exp[-energymu[mu]/(kb*T)],{mu,0,Ns-1}];
	Z2 = Sum[Exp[-energynu[nu]/(kb*T)],{nu,0,Ns-1}];
	lam=kcal2au*lambda;
	Sum[
		(Exp[-energymu[mu]/(kb*T)]/Z1)*
		ov[mu,nu]^2*
		Exp[-alpha[mu,nu]*dR]*
		Exp[(lamalph[mu,nu]-kcal2au*LambdaR[MR,omega,dR])/(hbar*omega*cm2au)]*
		Sqrt[Pi/(kb*T*lam)]*
		NIntegrate[
			(1-Fermi[eps,T])*Exp[-(energynu[nu] - energymu[mu] + kb*T*Log[Z2/Z1] + eps - eta*ev2au+lam)^2/
			(4*kb*T*lam)],
		{eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}]
];

CathodicCurrentDensityLowTFGH[fghlist_,eta_,lambda_,MR_,dR_,omega_,mumax_,numax_,T_]:=
Module[{eps,alpha,lamalph,ov,lam,energymu,energynu,Z1,Z2,Ns},
        Ns = Length[fghlist[[1,1]]];
        If[mumax+1>Ns||numax+1>Ns,Message[CathodicCurrentDensityLowTFGH::numberofstates];Abort[]];
	Array[alpha,{mumax+1,numax+1},0];
	Array[lamalph,{mumax+1,numax+1},0];
	Array[ov,{mumax+1,numax+1},0];
	Array[energymu,Ns,0];
	Array[energynu,Ns,0];
	Do[ov[mu,nu]=fghlist[[3]][[mu+1,nu+1]],{mu,0,mumax},{nu,0,numax}];
	Do[alpha[mu,nu]=fghlist[[4]][[mu+1,nu+1]],{mu,0,mumax},{nu,0,numax}];
	Do[lamalph[mu,nu]=hbar^2*alpha[mu,nu]^2/(2*MR*Dalton),{mu,0,mumax},{nu,0,numax}];
	Do[energymu[mu]=kcal2au*fghlist[[1]][[1,mu+1]],{mu,0,Ns-1}];
	Do[energynu[nu]=kcal2au*fghlist[[2]][[1,nu+1]],{nu,0,Ns-1}];
	Z1 = Sum[Exp[-energymu[mu]/(kb*T)],{mu,0,Ns-1}];
	Z2 = Sum[Exp[-energynu[nu]/(kb*T)],{nu,0,Ns-1}];
	lam=kcal2au*lambda;
	Sum[
		(Exp[-energynu[nu]/(kb*T)]/Z2)*
		ov[mu,nu]^2*
		Exp[-alpha[mu,nu]*dR]*
		Exp[(lamalph[mu,nu]-kcal2au*LambdaR[MR,omega,dR])/(hbar*omega*cm2au)]*
		Sqrt[Pi/(kb*T*lam)]*
		NIntegrate[
			Fermi[eps,T]*Exp[-(-energynu[nu] + energymu[mu] - kb*T*Log[Z2/Z1] - eps + eta*ev2au + lam)^2/
			(4*kb*T*lam)],
		{eps,-Infinity,Infinity}],
	{nu,0,numax},{mu,0,mumax}]
];

(*
fghlist - nested list with energies, wavefunctions, overlap integrals and alphas;
eta - overpotential in volts;
lambda - solvent reorganization energy in kcal/mol;
MR - reduced mass of the donor-acceptor mode in Daltons;
dR - distance between the equlibrium value of the R mode at product state and the reactant state (Rnu-Rmu) in Bohr;
omega - frequency of the coupling between the oscillators in wavenumbers;
mumax - maximum quantum number of oscillator on the left;
numax - maximum quantum number of oscillator on the right;
T - temperature in Kelvin;
Cox - concentration of the oxidized complex (any consistent units);
Cred - concentration of the reduced complex (any consistent units);
Returns the total current density for general proton potentials;
*)

TotalCurrentDensityLowTFGH[fghlist_,eta_,lambda_,MR_,dR_,omega_,mumax_,numax_,T_,Cox_,Cred_]:=
Module[{Conc},
	Conc=Cred*Exp[(eta*ev2au)/(kb*T)]/Cox-1;
	Cox*CathodicCurrentDensityLowTFGH[fghlist,eta,lambda,MR,dR,omega,mumax,numax,T]*Conc
];



Pathway1[Eta_,e1_,pKa1_,e2_,Lambda1_,Lambda2_,pKaAcid_,OptionsPattern[{Ref->"HA/H2"}]]:=Module[{a,b,c,d,f,g,ref,e0},
ref=OptionValue[Ref];
e0=Switch[ref,
"HA/H2",Log[10.]*r*t*pKaAcid,
"NHE",0.0,
"SCE",0.2,
_,Log[10.]*r*t*pKaAcid
];
a=-e1+Eta;
b=Log[10.]*r*t*(pKaAcid-pKa1);
c=Log[10.]*r*t*(pKaAcid+pKa1)+(e2+e1-2e0);
d=-e2+Eta;
f=(a+Lambda1)^2/(4Lambda1);
g=(d+Lambda2)^2/(4Lambda2);
Graphics[{
{FontSize->16,FontFamily-> "Arial",Text["Pathway 1",{0.75,1.85}]},
{FontSize->12,FontFamily-> "Arial", Text["2Co(II)",{0.5,.3}]},
{FontSize->12,FontFamily-> "Arial", Text["+ 2HA",{0.5,.2}]},
{FontSize->12,FontFamily-> "Arial", Text[Superscript["+ 2e","-"],{0.5,.1}]},
{FontSize->12,FontFamily-> "Arial",Text["Co(II)",{3,a+0.4}]},
{FontSize->12,FontFamily-> "Arial",Text["+ Co(I\!\(\*SuperscriptBox[\()\), \(-\)]\)",{3,a+0.3}]},
{FontSize->12,FontFamily-> "Arial",Text["+ 2HA",{3,a+0.2}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ e","-"],{3,a+0.1}]},
{FontSize->12,FontFamily-> "Arial",Text["Co(II)",{4.875,a+b-0.1}]},
{FontSize->12,FontFamily-> "Arial",Text["+ Co(III)H",{4.875,a+b-0.2}]},
{FontSize->12,FontFamily-> "Arial",Text["+ HA +" Superscript["A","-"],{4.875,a+b-0.3}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ e","-"],{4.875,a+b-0.4}]},
{FontSize->12,FontFamily-> "Arial",Text["2Co(I\!\(\*SuperscriptBox[\()\), \(-\)]\)",{5.5,2a-0.1}]},
{FontSize->12,FontFamily-> "Arial",Text["+ 2HA",{5.5,2a-0.2}]},
{FontSize->12,FontFamily-> "Arial",Text["Co(II)",{6.75,a+b+c-0.1}]},
{FontSize->12,FontFamily-> "Arial",Text["+ Co(III\!\(\*SuperscriptBox[\()\), \(+\)]\)",{6.75,a+b+c-0.2}]},
{FontSize->12,FontFamily-> "Arial",Text[Subscript["+ H",2]Superscript["+ 2A","-"],{6.75,a+b+c-0.3}]},
{FontSize->12,FontFamily-> "Arial",Text["Co(I\!\(\*SuperscriptBox[\()\), \(-\)]\)",{7,2a+b+0.3}]},
{FontSize->12,FontFamily-> "Arial",Text["+ Co(III)H",{7,2a+b+0.2}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ HA + A","-"],{7,2a+b+0.1}]},
{FontSize->12,FontFamily-> "Arial",Text["2Co(III)H",{8.5,2a+2b-0.1}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ 2A","-"],{8.5,2a+2b-0.2}]},
{FontSize->12,FontFamily-> "Arial",Text["2Co(II)",{10,a+b+c+d+0.2}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ 2A","-"]Subscript["+ H",2],{10,a+b+c+d+0.1}]},
{FontSize->12,FontFamily->"Arial",Text["Reference"ref,{0.75,-0.3}]},
{Thick,Line[{{0,0},{1,0}}]},
{Dashed,Line[{{1,0},{1.5,f}}]},
{Dotted,Line[{{1.5,f},{2,f}}]},
{Dashed,Line[{{2,f},{2.5,a}}]},
{Thick,Line[{{2.5,a},{3.5,a}}]},
{Dashed,Red,Line[{{3.5,a},{4.375,a+b}}]},
{Dashed,Blue,Line[{{3.5,a},{4,a+f}}]},
{Thick,Red,Line[{{4.375,a+b},{5.375,a+b}}]},
{Dotted,Blue,Line[{{4,a+f},{4.5,a+f}}]},
{Dashed,Red,Line[{{5.375,a+b},{6.25,a+b+c}}]},
{Dashed,Blue,Line[{{4.5,a+f},{5,2a}}]},
{Thick,Red,Line[{{6.25,a+b+c},{7.25,a+b+c}}]},
{Thick,Blue,Line[{{5,2a},{6,2a}}]},
{Dashed,Red,Line[{{7.25,a+b+c},{8.125,a+b+c+g}}]},
{Dashed,Blue,Line[{{6,2a},{6.5,2a+b}}]},
{Dotted,Red,Line[{{8.125,a+b+c+g},{8.625,a+b+c+g}}]},
{Thick,Blue,Line[{{6.5,2a+b},{7.5,2a+b}}]},
{Dashed,Red,Line[{{8.625,a+b+c+g},{9.5,a+b+c+d}}]},
{Dashed,Blue,Line[{{7.5,2a+b},{8,2a+2b}}]},
{Thick,Blue,Line[{{8,2a+2b},{9,2a+2b}}]},
{Dashed,Blue,Line[{{9,2a+2b},{9.5,a+b+c+d}}]},
{Thick,Line[{{9.5,a+b+c+d},{10.5,a+b+c+d}}]}
},Axes->False,Frame->True,FrameLabel->{"Reaction Coordinate","Relative Free Energy (eV)"},
LabelStyle->(FontFamily->"Arial"),FrameTicks->{False,True,False,True},AspectRatio->0.5,ImageSize->{800},
BaseStyle->{FontSize->16},PlotRange-> {{-0.25,10.75},{-0.5,2}}]
];






Pathway2[Eta_,e1_,pKa1_,e3_,pKa2_,e4_,Lambda1_,Lambda3_,pKaAcid_,OptionsPattern[{Ref->"HA/H2"}]]:=Module[{a,b,h,i,j,f,k,l,ref,e0},
ref=OptionValue[Ref];
e0=Switch[ref,
"HA/H2",Log[10.]*r*t*pKaAcid,
"NHE",0.0,
"SCE",0.2,
_,Log[10.]*r*t*pKaAcid
];
a=-e1+Eta;
b=Log[10.]*r*t*(pKaAcid-pKa1);
h=-e3+Eta;
i=Log[10.]*r*t*(pKaAcid+pKa2)+(e1+e4-2e0);
j=2Log[10.]*r*t*pKa2+2(e4-e0);
f=(a+Lambda1)^2/(4Lambda1);
k=(h+Lambda3)^2/(4Lambda3);
l=(-a+Lambda1)^2/(4Lambda1);
      Graphics[{
{FontSize->16,FontFamily-> "Arial",Text["Pathway 2",{-3.5,2.8}]},
{FontSize->12,FontFamily-> "Arial", Text["2Co(II)",{-4.5,0.45}]},
{FontSize->12,FontFamily-> "Arial", Text["+ 2HA",{-4.5,0.3}]},
{FontSize->12,FontFamily-> "Arial", Text[Superscript["+ 4e","-"],{-4.5,0.15}]},
{FontSize->12,FontFamily-> "Arial",Text["Co(II)",{-2,a+0.52}]},
{FontSize->12,FontFamily-> "Arial",Text["+ Co(I\!\(\*SuperscriptBox[\()\), \(-\)]\)",{-2,a+0.39}]},
{FontSize->12,FontFamily-> "Arial",Text["+ 2HA",{-2,a+0.26}]},
{FontSize->12,FontFamily-> "Arial", Text[Superscript["+ 3e","-"],{-2,a+0.13}]},
{FontSize->12,FontFamily-> "Arial",Text["2Co(I\!\(\*SuperscriptBox[\()\), \(-\)]\)",{0.5,2a-0.13}]},
{FontSize->12,FontFamily-> "Arial",Text["+ 2HA",{0.5,2a-0.26}]},
{FontSize->12,FontFamily-> "Arial", Text[Superscript["+ 2e","-"],{0.5,2a-0.39}]},
{FontSize->12,FontFamily-> "Arial",Text["Co(III)H",{2,2a+b+0.65}]},
{FontSize->12,FontFamily-> "Arial",Text["+ Co(I\!\(\*SuperscriptBox[\()\), \(-\)]\)",{2,2a+b+0.52}]},
{FontSize->12,FontFamily-> "Arial",Text["+ HA",{2,2a+b+0.39}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ A","-"],{2,2a+b+0.26}]},
{FontSize->12,FontFamily-> "Arial", Text[Superscript["+ 2e","-"],{2,2a+b+0.13}]},
{FontSize->12,FontFamily-> "Arial",Text["Co(II)",{2.375,a+b-0.13}]},
{FontSize->12,FontFamily-> "Arial",Text["+ Co(III)H",{2.37,a+b-0.26}]},
{FontSize->12,FontFamily-> "Arial",Text["+ HA +" Superscript["A","-"],{2.375,a+b-0.39}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ 3e","-"],{2.375,a+b-0.52}]},
{FontSize->12,FontFamily-> "Arial",Text["2Co(III)H",{3.45,2a+2b+0.39}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ 2A","-"],{3.45,2a+2b+0.26}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ 2e","-"],{3.45,2a+2b+0.13}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["Co(II)H","-"],{5.8,2a+2b+h+0.52}]},
{FontSize->12,FontFamily-> "Arial",Text["+ Co(III)H",{5.8,2a+2b+h+0.39}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ 2A","-"],{5.8,2a+2b+h+0.26}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ e","-"],{5.8,2a+2b+h+0.13}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["2Co(II)H","-"],{8.5,2a+2b+2h+0.26}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ 2A","-"],{8.5,2a+2b+2h+0.13}]},
{FontSize->12,FontFamily-> "Arial",Text["2Co(I\!\(\*SuperscriptBox[\()\), \(-\)]\)",{10,2a+2b+2h+j-0.13}]},
{FontSize->12,FontFamily->"Arial",Text[Subscript["+ H",2]Superscript["+ 2A","-"],{10,2a+2b+2h+j-0.26}]},
{FontSize->12,FontFamily-> "Arial",Text["Co(II)",{10.625,a+b+h-0.13}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ Co(II)H","-"],{10.625,a+b+h-0.26}]},
{FontSize->12,FontFamily-> "Arial",Text["+ HA +" Superscript["A","-"],{10.625,a+b+h-0.39}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ 2e","-"],{10.625,a+b+h-0.52}]},
{FontSize->12,FontFamily-> "Arial",Text["Co(I\!\(\*SuperscriptBox[\()\), \(-\)]\)",{12.5,a+2b+2h+j+0.65}]},
{FontSize->12,FontFamily-> "Arial",Text["+ Co(II)",{12.5,a+2b+2h+j+0.52}]},
{FontSize->12,FontFamily->"Arial",Text[Subscript["+ H",2],{12.5,a+2b+2h+j+0.39}]},
{FontSize->12,FontFamily->"Arial",Text[Superscript["+ 2A","-"],{12.5,a+2b+2h+j+0.26}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ e","-"],{12.5,a+2b+2h+j+0.13}]},
{FontSize->12,FontFamily-> "Arial",Text["2Co(II)",{15,a+b+h+i+0.45}]},
{FontSize->12,FontFamily->"Arial",Text[Subscript["+ H",2]Superscript["+ 2A","-"],{15,a+b+h+i+0.30}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ 2e","-"],{15,a+b+h+i+0.15}]},
{FontSize->12,FontFamily->"Arial",Text["Reference"ref,{-3.5,-0.3}]},
{Thick,Line[{{-5,0},{-4,0}}]},
{Dashed,Line[{{-4,0},{-3.5,f}}]},
{Dotted,Line[{{-3.5,f},{-3,f}}]},
{Dashed,Line[{{-3,f},{-2.5,a}}]},
{Thick,Line[{{-2.5,a},{-1.5,a}}]},
{Dashed,Red,Line[{{-1.5,a},{1.875,a+b}}]},
{Dashed,Blue,Line[{{-1.5,a},{-1,a+f}}]},
{Dotted,Blue,Line[{{-1,a+f},{-0.5,a+f}}]},
{Dashed,Blue,Line[{{-0.5,a+f},{0,2a}}]},
{Thick,Blue,Line[{{0,2a},{1,2a}}]},
{Dashed,Blue,Line[{{1,2a},{1.5,2a+b}}]},
{Thick,Blue,Line[{{1.5,2a+b},{2.5,2a+b}}]},
{Thick,Red,Line[{{1.875,a+b},{2.875,a+b}}]},
{Dashed,Blue,Line[{{2.5,2a+b},{3,2a+2b}}]},
{Dashed,Red,Line[{{2.875,a+b},{6.25,a+b+k}}]},
{Thick,Blue,Line[{{3,2a+2b},{4,2a+2b}}]},
{Dashed,Blue,Line[{{4,2a+2b},{4.5,2a+2b+k}}]},
{Dotted,Blue,Line[{{4.5,2a+2b+k},{5,2a+2b+k}}]},
{Dashed,Blue,Line[{{5,2a+2b+k},{5.5,2a+2b+h}}]},
{Thick,Blue,Line[{{5.5,2a+2b+h},{6.5,2a+2b+h}}]},
{Dotted,Red,Line[{{6.25,a+b+k},{6.75,a+b+k}}]},
{Dashed,Blue,Line[{{6.5,2a+2b+h},{7,2a+2b+2k}}]},
{Dashed,Red,Line[{{6.75,a+b+k},{10.125,a+b+h}}]},
{Dotted,Blue,Line[{{7,2a+2b+2k},{7.5,2a+2b+2k}}]},
{Dashed,Blue,Line[{{7.5,2a+2b+2k},{8,2a+2b+2h}}]},
{Thick,Blue,Line[{{8,2a+2b+2h},{9,2a+2b+2h}}]},
{Dashed,Blue,Line[{{9,2a+2b+2h},{9.5,2a+2b+2h+j}}]},
{Thick,Blue,Line[{{9.5,2a+2b+2h+j},{10.5,2a+2b+2h+j}}]},
{Thick,Red,Line[{{10.125,a+b+h},{11.125,a+b+h}}]},
{Dashed,Blue,Line[{{10.5,2a+2b+2h+j},{11,2a+2b+2h+j+l}}]},
{Dotted,Blue,Line[{{11,2a+2b+2h+j+l},{11.5,2a+2b+2h+j+l}}]},
{Dashed,Red,Line[{{11.125,a+b+h},{14.5,a+b+h+i}}]},
{Dashed,Blue,Line[{{11.5,2a+2b+2h+j+l},{12,a+2b+2h+j}}]},
{Thick,Blue,Line[{{12,a+2b+2h+j},{13,a+2b+2h+j}}]},
{Dashed,Blue,Line[{{13,a+2b+2h+j},{13.5,a+2b+2h+j+l}}]},
{Dotted,Blue,Line[{{13.5,a+2b+2h+j+l},{14,a+2b+2h+j+l}}]},
{Dashed,Blue,Line[{{14,a+2b+2h+j+l},{14.5,a+b+h+i}}]},
{Thick,Line[{{14.5,a+b+h+i},{15.5,a+b+h+i}}]}
},Axes->False,Frame->True,FrameLabel->{"Reaction Coordinate","Relative Free Energy (eV)"},
FrameTicks->{False,True,False,True},AspectRatio->0.5,ImageSize->{800},BaseStyle->{FontSize->16},
LabelStyle->(FontFamily->"Arial"),PlotRange->{{-5.25,16},{-1,4}}]
];


Pathway3[Eta_,e1_,Lambda1_,e4_,Lambda4_,pKa2_,pKaAcid_,OptionsPattern[{Ref->"HA/H2"}]]:=Module[{a,f,n,m,p,j,i,l,ref,e0},
ref=OptionValue[Ref];
e0=Switch[ref,
"HA/H2",Log[10.]*r*t*pKaAcid,
"NHE",0.0,
"SCE",0.2,
_,Log[10.]*r*t*pKaAcid
];
a=-e1+Eta;
f=(a+Lambda1)^2/(4Lambda1);
m=(n+Lambda4)^2/(4Lambda4);
n=-e4+Eta;
p=Log[10.]*r*t*(pKaAcid-pKa2);
j=2Log[10.]*r*t*pKa2+2(e4-e0);
i=Log[10.]*r*t*(pKaAcid+pKa2)+(e1+e4-2e0);
l=(-a+Lambda1)^2/(4Lambda1);
Graphics[{
{FontSize->16,FontFamily-> "Arial",Text["Pathway 3",{-3.5,4.2}]},
{FontSize->12,FontFamily-> "Arial",Text["2Co(II)",{-4.5,0.6}]},
{FontSize->12,FontFamily-> "Arial",Text["+ 2HA",{-4.5,0.4}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ 4e","-"],{-4.5,0.2}]},
{FontSize->12,FontFamily-> "Arial",Text["Co(II)",{-2,a+0.8}]},
{FontSize->12,FontFamily-> "Arial",Text["+ Co(I\!\(\*SuperscriptBox[\()\), \(-\)]\)",{-2,a+0.6}]},
{FontSize->12,FontFamily-> "Arial",Text["+ 2HA",{-2,a+0.4}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ 3e","-"],{-2,a+0.2}]},
{FontSize->12,FontFamily-> "Arial",Text["2Co(I\!\(\*SuperscriptBox[\()\), \(-\)]\)",{0.4,2a+0.6}]},
{FontSize->12,FontFamily-> "Arial",Text["+ 2HA",{0.4,2a+0.4}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ 2e","-"],{0.4,2a+0.2}]},
{FontSize->12,FontFamily-> "Arial",Text["Co(0\!\(\*SuperscriptBox[\()\), \(\(2\)\(-\)\)]\)",{3,2a+n+0.8}]},
{FontSize->12,FontFamily-> "Arial",Text["+ Co(I\!\(\*SuperscriptBox[\()\), \(-\)]\)",{3,2a+n+0.6}]},
{FontSize->12,FontFamily-> "Arial",Text["+ 2HA",{3,2a+n+0.4}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ e","-"],{3,2a+n+0.2}]},
{FontSize->12,FontFamily-> "Arial",Text["2Co(0\!\(\*SuperscriptBox[\()\), \(\(2\)\(-\)\)]\)",{5.25,2a+2n-0.2}]},
{FontSize->12,FontFamily-> "Arial",Text["+ 2HA",{5.25,2a+2n-0.4}]},
{FontSize->12,FontFamily-> "Arial",Text["Co(II)",{6.25,a+n-0.2}]},
{FontSize->12,FontFamily-> "Arial",Text["+ Co(0\!\(\*SuperscriptBox[\()\), \(\(2\)\(-\)\)]\)",{6.25,a+n-0.4}]},
{FontSize->12,FontFamily-> "Arial",Text["+ 2HA",{6.25,a+n-0.6}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ 2e","-"],{6.25,a+n-0.8}]},
{FontSize->12,FontFamily-> "Arial",Text["Co(0\!\(\*SuperscriptBox[\()\), \(\(2\)\(-\)\)]\)",{7.5,2a+2n+p+0.6}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ Co(II)H","-"],{7.5,2a+2n+p+0.4}]},
{FontSize->12,FontFamily-> "Arial",Text["+ HA +" Superscript["A","-"],{7.5,2a+2n+p+0.2}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["2Co(II)H","-"],{8.5,2a+2n+2p-0.2}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ 2A","-"],{8.5,2a+2n+2p-0.4}]},
{FontSize->12,FontFamily-> "Arial",Text["2Co(I\!\(\*SuperscriptBox[\()\), \(-\)]\)",{10,2a+2n+2p+j+0.4}]},
{FontSize->12,FontFamily->"Arial",Text[Subscript["+ H",2]Superscript["+ 2A","-"],{10,2a+2n+2p+j+0.2}]},
{FontSize->12,FontFamily-> "Arial",Text["Co(II)",{10.625,a+n+p-0.2}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ Co(II)H","-"],{10.625,a+n+p-0.4}]},
{FontSize->12,FontFamily-> "Arial",Text["+ HA +" Superscript["A","-"],{10.625,a+n+p-0.6}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ 2e","-"],{10.625,a+n+p-0.8}]},
{FontSize->12,FontFamily-> "Arial",Text["Co(I\!\(\*SuperscriptBox[\()\), \(-\)]\)",{12.75,a+2n+2p+j+l+0.8}]},
{FontSize->12,FontFamily-> "Arial",Text["+ Co(II)",{12.75,a+2n+2p+j+l+0.6}]},
{FontSize->12,FontFamily->"Arial",Text[Subscript["+ H",2]Superscript["+ 2A","-"],{12.75,a+2n+2p+j+l+0.4}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ e","-"],{12.75,a+2n+2p+j+l+0.2}]},
{FontSize->12,FontFamily-> "Arial",Text["2Co(II)",{15,a+n+p+i+0.8}]},
{FontSize->12,FontFamily->"Arial",Text[Subscript["+ H",2]Superscript["+ 2A","-"],{15,a+n+p+i+0.6}]},
{FontSize->12,FontFamily-> "Arial",Text[Superscript["+ 2e","-"],{15,a+n+p+i+0.4}]},
{FontSize->12,FontFamily->"Arial",Text["Reference"ref,{0.75,-0.3}]},
{Thick,Line[{{-5,0},{-4,0}}]},
{Dashed,Line[{{-4,0},{-3.5,f}}]},
{Dotted,Line[{{-3.5,f},{-3,f}}]},
{Dashed,Line[{{-3,f},{-2.5,a}}]},
{Thick,Line[{{-2.5,a},{-1.5,a}}]},
{Dashed,Red,Line[{{-1.5,a},{1.875,a+m}}]},
{Dashed,Blue,Line[{{-1.5,a},{-1,a+f}}]},
{Dotted,Blue,Line[{{-1,a+f},{-0.5,a+f}}]},
{Dashed,Blue,Line[{{-0.5,a+f},{0,2a}}]},
{Thick,Blue,Line[{{0,2a},{1,2a}}]},
{Dashed,Blue,Line[{{1,2a},{1.5,2a+m}}]},
{Dotted,Blue,Line[{{1.5,2a+m},{2,2a+m}}]},
{Dotted,Red,Line[{{1.875,a+m},{2.375,a+m}}]},
{Dashed,Blue,Line[{{2,2a+m},{2.5,2a+n}}]},
{Dashed,Red,Line[{{2.375,a+m},{5.75,a+n}}]},
{Thick,Blue,Line[{{2.5,2a+n},{3.5,2a+n}}]},
{Dashed,Blue,Line[{{3.5,2a+n},{4,2a+n+m}}]},
{Dotted,Blue,Line[{{4,2a+n+m},{4.5,2a+n+m}}]},
{Dashed,Blue,Line[{{4.5,2a+n+m},{5,2a+2n}}]},
{Thick,Blue,Line[{{5,2a+2n},{6,2a+2n}}]},
{Thick,Red,Line[{{5.75,a+n},{6.75,a+n}}]},
{Dashed,Blue,Line[{{6,2a+2n},{6.5,2a+2n+p}}]},
{Thick,Blue,Line[{{6.5,2a+2n+p},{7.5,2a+2n+p}}]},
{Dashed,Red,Line[{{6.75,a+n},{10.125,a+n+p}}]},
{Dashed,Blue,Line[{{7.5,2a+2n+p},{8,2a+2n+2p}}]},
{Thick,Blue,Line[{{8,2a+2n+2p},{9,2a+2n+2p}}]},
{Dashed,Blue,Line[{{9,2a+2n+2p},{9.5,2a+2n+2p+j}}]},
{Thick,Blue,Line[{{9.5,2a+2n+2p+j},{10.5,2a+2n+2p+j}}]},
{Dashed,Blue,Line[{{10.5,2a+2n+2p+j},{11,2a+2n+2p+j+l}}]},
{Thick,Red,Line[{{10.125,a+n+p},{11.125,a+n+p}}]},
{Dotted,Blue,Line[{{11,2a+2n+2p+j+l},{11.5,2a+2n+2p+j+l}}]},
{Dashed,Red,Line[{{11.125,a+n+p},{14.5,a+n+p+i}}]},
{Dashed,Blue,Line[{{11.5,2a+2n+2p+j+l},{12,a+2n+2p+j+l}}]},
{Thick,Blue,Line[{{12,a+2n+2p+j+l},{13,a+2n+2p+j+l}}]},
{Dashed,Blue,Line[{{13,a+2n+2p+j+l},{13.5,a+2n+2p+j+2l}}]},
{Dotted,Blue,Line[{{13.5,a+2n+2p+j+2l},{14,a+2n+2p+j+2l}}]},
{Dashed,Blue,Line[{{14,a+2n+2p+j+2l},{14.5,a+n+p+i}}]},
{Thick,Line[{{14.5,a+n+p+i},{15.5,a+n+p+i}}]}
},Axes->False,Frame->True,FrameLabel->{"Reaction Coordinate","Relative Free Energy (eV)"},
FrameTicks->{False,True,False,True},AspectRatio->0.5,ImageSize->{800},BaseStyle->{FontSize->16},
LabelStyle->(FontFamily->"Arial"),PlotRange->{{-5.25,16},{-0.5,4.5}}]
];

(*
Closing statements
*)

End[]

EndPackage[]
