(* ::Package:: *)

BeginPackage["PourbaixDiagram`"]
Needs["ComputationalGeometry`"]

listedge::usage="create list of known edges to be used in Pourbaix input";
GraphPourbaix::usage=
"PourbaixDiagram[specieslist,valueslist,ne,T,pHmin,pHmax,E0min,E0max]
Returns Pourbaix (phase) diagram that relates potential vs. pH, listing the most stable species for a given set of parameters:

--Input parameters:
* specieslist - list of species given in a matrix/nested list form. If some species is not found , use  '*'  instead.
Example: { {M(III), M(II), M(I)}, {M(III)H, M(II)H, * }, {*, M(II)H2, M(I)H2},..}  
 left to right : decreasing oxidation state
 up to down: increasing protonation
* valueslist - list of values connecting a pair of species in a nested list form
Format: { {<E0 from M(II)to M(I)>, M(II), M(I)}, {<pKa of MH to M>, M, MH} }
Example: { {-1.2, M(III), M(II)}, {10, M(III),M(III)H}, {-0.8,  M(II)H2, M(I)H2 },..}
Here, -1.2 is the reduction potential corresponding to M(III)/M(III) and 10 is pKa corresponding to the deprotonation of M(III)H to form M(III).
* ne - number of electrons transfered in reduction steps;
* T - the absolute temperature in K;
* pHmin,pHmax - pH range of the diagram;
* E0min,E0max - potential range of the diagram;
-Options:
* RefIndex - Index of the desired reference species in the flattened specieslist. Reference species is the one with zero relative free energy, 
w.r.t which all other free energies are calculated. For instance, if the reference spcies you want is M(II)H2 in the spcieslist given by
{ {M(III), M(II), M(I)}, {M(III)H, M(II)H, x }, {x, M(II)H2, M(I)H2},..} , then use RefIndex->8 (see next line to know why) in the input function. 
       1       2       3       4       5      6      7       8
Default value of RefIndex is 1. 
* Ref - Name of the reference electrode (default: Fc+/Fc);
* Solvent - Name of the solvent;
* WriteFiles (True/False) - if True, write out files with the coordinates of the triple points, points on the edges, and vertices of stability regions 
 (NOTE: Do not use subscripts and superscripts in the labels if you want to see clean items in the output files);
* ColorScheme - Name of the color scheme to use for regions of the Pourbaix diagram (default : \"Rainbow\").
  The following schemes can be used:  \n "<>ToString[ColorData["Gradients"]];
showgraph::usage="Show the graph which is used to generate the Pourbaix Diagram";
printcomments::usage="Print comments on the given data (in case of insufficient data";
printdata::usage="Print the coordinates for triple points, edge points and the stability regions";

Begin["`Private`"];


listedge[alledgelist_]:=
Module[
{valueslist,x},
valueslist={};
x={};
For[i=1,i<=Length[alledgelist],i++,
If[MemberQ[alledgelist[[i]],Null],Null,AppendTo[valueslist,alledgelist[[i]]]];
];

Return[valueslist]
];


GraphPourbaix[specieslist_,species_,verticals_,horizontals_,ne_,T_,pHmin_,pHmax_,E0min_,E0max_,
OptionsPattern[{RefIndex->1,Ref-> "Fc+/Fc",Solvent->  "MeCN",Tol->  0.0,ColorScheme->  "Rainbow"}]]:= 
Module[{speciesall,ncolumnlist,valueslist,edgelist,nspeciesrow,eid,Nspecies,Nvalues,vcoord,ps,es,Rln10,F,lS,con,speciesnull,
path,dGi,speciesLabels,nS,id,valId,dG,G,Stable,uniqueTriples,triplePoints,sol,solpH,solE0,ln10,
edgePoints,uniquePairs,regions,point3,tLabel,regionPoints,uniqueRegions,cornerPoints,notEmpty,ref,solvent,tol,mm,
writefiles,colorscheme,vertexFile,vertexChannel,regionFile,regionChannel,a,b,c,i,j,k,iT,iP,ri},

(*Constants*)
F=23.0605488670097;
ln10=Log[10.0];
Rln10=0.00198720414566658*ln10;
(* options *)
ref=OptionValue[Ref];
solvent=OptionValue[Solvent];
tol=OptionValue[Tol];
writefiles=OptionValue[WriteFiles];
colorscheme=OptionValue[ColorScheme];
ri=OptionValue[RefIndex];
(* Defining variables and lists *)
valueslist=Join[verticals,horizontals];
nspeciesrow=Length[specieslist]; 
ncolumnlist=Table[Length[specieslist[[i]]],{i,1,nspeciesrow}]; 
lS=Length[species];
speciesnull={};
For[i=1,i<=lS,i++,
If[species[[i]]==Null,AppendTo[speciesnull,species[[i]]]
];
];
speciesall={};
For[i=1,i<=lS,i++,
If[MemberQ[speciesnull,species[[i]]], Null, AppendTo[speciesall,species[[i]]]
];
];

Nspecies=Length[speciesall];
Nvalues=Length[valueslist]; 
edgelist=Table[valueslist[[i,2]] \[UndirectedEdge] valueslist[[i,3]],{i,1,Nvalues}];

eid={};
pid={};
For[i=1,i<=Nvalues,i++,
If[MemberQ[horizontals,valueslist[[i]]],AppendTo[eid,1];AppendTo[pid,0],
AppendTo[eid,0];AppendTo[pid,1]]
];


For[i=1,i<=nspeciesrow,i++,
For[j=1,j<=Length[specieslist[[i]]],j++,
For[k=1,k<=lS,k++,
If[MemberQ[specieslist[[i]],species[[k]]],ps[species[[k]]]=i];
If[MemberQ[specieslist[[All,j]],species[[k]]],es[species[[k]]]=j];
]]];




G=Graph[speciesall,edgelist,VertexLabels->Placed["Name",Center]];
dGi[i_,pH_,E0_]:= If[eid[[i]]== 1,If[es[valueslist[[i,2]]]<es[valueslist[[i,3]]],-ne*F*(valueslist[[i,1]]-E0),ne*F*(valueslist[[i,1]]-E0)],
If[ps[valueslist[[i,2]]]<ps[valueslist[[i,3]]],Rln10*T*(pH-valueslist[[i,1]]),-Rln10*T*(pH-valueslist[[i,1]])]];

speciesLabels={};
For[i=1,i<=  Length[speciesall],i++,
If[Length[FindShortestPath[G,speciesall[[ri]],speciesall[[i]]]]!= 0 ,AppendTo[speciesLabels,speciesall[[i]]] 
];
];

nS=Length[speciesLabels];
path[i_]:=FindShortestPath[G,speciesLabels[[ri]],speciesLabels[[i]]];
For[i=1,i<=nS,i++,
For[j=1,j<Length[path[i]] ,j++,
For[k=1,k<=Nvalues,k++,
If[path[i][[j]]\[UndirectedEdge] path[i][[j+1]]===edgelist[[k]] \[Or] path[i][[j+1]]\[UndirectedEdge] path[i][[j]]=== edgelist[[k]] ,valId[i,j]=k]
];
];
];
dG[pH_,E0_]:=Table[Sum[dGi[valId[i,j],pH,E0] ,{j,1,Length[path[i]]-1}],{i,1,nS}];

(* This logical function defines a stability region for species i *)
Stable[pH_,E0_,i_]:=Ordering[dG[pH,E0],1][[1]]== i;

(* Find all triple points on the diagram *)
triplePoints={};
uniqueTriples=Subsets[Range[nS],{3}];
 For[iT=1,iT<Length[uniqueTriples]+1,iT++,
{i,j,k}=uniqueTriples[[iT]];
sol=Solve[dG[x,y][[i]]== dG[x,y][[j]]== dG[x,y][[k]]&&pHmin<x<pHmax&&E0min<y<E0max,{x,y}];
If[sol== {},Continue[]];
{solpH,solE0}=sol[[1,All,2]];
If[Stable[solpH,solE0,i]\[Or]Stable[solpH,solE0,j]\[Or]Stable[solpH,solE0,k],
AppendTo[triplePoints,{{solpH,solE0},speciesLabels[[{i,j,k}]]}]
];
];
(*
*)
(* Find all edge points on the diagram *)
edgePoints={};
uniquePairs=Subsets[Range[nS],{2}];

(* Bottom edge of the frame *)
For[iP=1,iP<Length[uniquePairs]+1,iP++,
{i,j}=uniquePairs[[iP]];
sol=Solve[dG[x,E0min][[i]]== dG[x,E0min][[j]]&&pHmin<x<pHmax,x];
If[sol== {},Continue[]];
solpH=sol[[1,1,2]];
If[Stable[solpH,E0min,i]\[Or]Stable[solpH,E0min,j],
AppendTo[edgePoints,{{solpH,E0min},speciesLabels[[{i,j}]]}]
];
];

(* Top edge of the frame *)
For[iP=1,iP<Length[uniquePairs]+1,iP++,
{i,j}=uniquePairs[[iP]];
sol=Solve[dG[x,E0max][[i]]== dG[x,E0max][[j]]&&pHmin<x<pHmax,x];
If[sol== {},Continue[]];
solpH=sol[[1,1,2]];
If[Stable[solpH,E0max,i]\[Or]Stable[solpH,E0max,j],
AppendTo[edgePoints,{{solpH,E0max},speciesLabels[[{i,j}]]}]
];
];

(* Left edge of the frame *)
For[iP=1,iP<Length[uniquePairs]+1,iP++,
{i,j}=uniquePairs[[iP]];
sol=Solve[dG[pHmin,y][[i]]== dG[pHmin,y][[j]]&&E0min<y<E0max,y];
If[sol== {},Continue[]];
solE0=sol[[1,1,2]];
If[Stable[pHmin,solE0,i]\[Or]Stable[pHmin,solE0,j],
AppendTo[edgePoints,{{pHmin,solE0},speciesLabels[[{i,j}]]}]
];
];

(* Right edge of the frame *)
For[iP=1,iP<Length[uniquePairs]+1,iP++,
{i,j}=uniquePairs[[iP]];
sol=Solve[dG[pHmax,y][[i]]== dG[pHmax,y][[j]]&&E0min<y<E0max,y];
If[sol== {},Continue[]];
solE0=sol[[1,1,2]];
If[Stable[pHmax,solE0,i]\[Or]Stable[pHmax,solE0,j],
AppendTo[edgePoints,{{pHmax,solE0},speciesLabels[[{i,j}]]}]
];
];

(* Build stability regions as collection of vertices selected from triple points and edgepoints *)

regions={};

(* Loop over triple points *)
Do[
point3=triplePoints[[iTPoint]];
(* Loop over labels of each triple point *)
Do[
tLabel=point3[[2,iTLabel]];
regionPoints={};
(* loop over edge points *)
Do[
If[tLabel== edgePoints[[iP2]][[2,1]]\[Or]tLabel== edgePoints[[iP2]][[2,2]],
AppendTo[regionPoints,edgePoints[[iP2]][[1]]]
]
,{iP2,1,Length[edgePoints]}];
(* loop over triple points *)
Do[
If[tLabel== triplePoints[[iP3]][[2,1]]\[Or]tLabel== triplePoints[[iP3]][[2,2]]\[Or]tLabel== triplePoints[[iP3]][[2,3]],
AppendTo[regionPoints,triplePoints[[iP3]][[1]]]
]
,{iP3,1,Length[triplePoints]}];
AppendTo[regions,{regionPoints,tLabel}]
,{iTLabel,1,3}]
,{iTPoint,1,Length[triplePoints]}];

(* Get rid of duplicate regions *)
uniqueRegions=Union[regions,SameTest-> (#1[[2]]==#2[[2]]&)]; 

(* Add corner points *)
cornerPoints={{pHmin,E0min},{pHmin,E0max},{pHmax,E0max},{pHmax,E0min}};
(* loop over four corner points *)
Do[
(* loop over unique regions *)
Do[
regionPoints=uniqueRegions[[iRegion,1]];
If[MemberQ[regionPoints[[All,1]],cornerPoints[[iCorner,1]]]&&MemberQ[regionPoints[[All,2]],cornerPoints[[iCorner,2]]],
uniqueRegions=Insert[uniqueRegions,cornerPoints[[iCorner]],{iRegion,1,-1}]
]
,{iRegion,1,Length[uniqueRegions]}]
,{iCorner,1,4}];

(*Print["\nUNIQUE REGIONS"];*)
(* Scan[Print,Sort[uniqueRegions]]; *)
(*Print["\nAREAS OF UNIQUE REGIONS:"];
Scan[Print,Sort[Table[{ConvexHullArea[uniqueRegions[[i,1]][[ConvexHull[uniqueRegions[[i,1]]]]]],uniqueRegions[[i,2]]},{i,1,Length[uniqueRegions]}]]];
*)
(* Mark the empty regions *)
notEmpty=Table[MemberQ[uniqueRegions[[All,2]],speciesLabels[[i]]],{i,1,nS}];
(*Print["Species not having stability regions in the given range of pH and \!\(\*SuperscriptBox[\(E\), \(0\)]\)"];
For[i=1,i<= nS,i++, 
If[notEmpty[[i]],Null,Print[speciesLabels[[i]]];]];*)
Show[
Graphics[Table[{EdgeForm[Directive[Thick,White]],Opacity[0.5],ColorData[colorscheme][iRegion/Length[uniqueRegions]],
                Polygon[uniqueRegions[[iRegion,1]][[ConvexHull[uniqueRegions[[iRegion,1]]]]]]},
                {iRegion,1,Length[uniqueRegions]}]],
Graphics[Table[Text[uniqueRegions[[i,2]],{Total[uniqueRegions[[i,1]][[All,1]]],Total[uniqueRegions[[i,1]][[All,2]]]}/
Length[uniqueRegions[[i,1]]],{0,0},BaseStyle->  {FontFamily->  "Helvetica",FontSize->  14,Plain}],{i,1,Length[uniqueRegions]}]],
ImageSize-> Large,AspectRatio-> 0.75,Frame->True,
FrameLabel->{{"Potential [V] vs. "<>ref<>" in "<>solvent,Null},{"pH",Text[Style["Pourbaix Diagram",Bold,FontFamily->"Helvetica"]]}},
(*FrameLabel->{{"Potential [V] ",Null},{"pH",Text[Style["Pourbaix Diagram",Bold,FontFamily->"Helvetica"]]}},*)
BaseStyle->{FontSize->16,FontFamily->"Helvetica"}]
];


showgraph[specieslist_,species_,verticals_,horizontals_]:=
Module[{speciesall,ncolumnlist,valueslist,edgelist,nspeciesrow,eid,Nspecies,Nvalues,vcoord,ps,es,Rln10,F,lS,con,speciesnull,pid,id},
(* Defining variables and lists *)
valueslist=Join[verticals,horizontals];
nspeciesrow=Length[specieslist]; 
ncolumnlist=Table[Length[specieslist[[i]]],{i,1,nspeciesrow}]; 
lS=Length[species];
speciesnull={};
For[i=1,i<=lS,i++,
If[species[[i]]==Null,AppendTo[speciesnull,species[[i]]]
];
];
speciesall={};
For[i=1,i<=lS,i++,
If[MemberQ[speciesnull,species[[i]]], Null, AppendTo[speciesall,species[[i]]]
];
];

Nspecies=Length[speciesall];
Nvalues=Length[valueslist]; 
edgelist=Table[valueslist[[i,2]] \[UndirectedEdge] valueslist[[i,3]],{i,1,Nvalues}];

eid={};
pid={};
For[i=1,i<=Nvalues,i++,
If[MemberQ[horizontals,valueslist[[i]]],AppendTo[eid,1];AppendTo[pid,0],
AppendTo[eid,0];AppendTo[pid,1]]
];

id={};
(* Generating the  list of edgelabels*)
For[i=1,i<=Nvalues,i++,
If[eid[[i]]== 1,AppendTo[id,"E0 = "],AppendTo[id,"pKa = "]];
];
(* Generating the Vertex Coordinates*) 
vcoord={};
For[i=1,i<= nspeciesrow,i++,
For[j=1,j<= ncolumnlist[[i]],j++,
If[specieslist[[i,j]]===0,Null,AppendTo[vcoord,{(j-1),-(i-1)}]];
];
];
(* Make a graph from given input *)
Graph[speciesall,edgelist,
EdgeLabels-> Table[edgelist[[i]]->  Placed[id[[i]] <>ToString[valueslist[[i,1]]],0.5],{i,1,Nvalues}],
EdgeLabelStyle-> Directive[Italic,12,FontFamily-> "Helvetica"], EdgeStyle-> Thick,
VertexCoordinates-> vcoord,VertexLabels-> Placed["Name",Center],VertexLabelStyle-> Directive[Bold,12,FontFamily-> "Helvetica"],
VertexSize-> 0.30,VertexShapeFunction-> None,ImageSize-> Large]
];


printcomments[specieslist_,species_,verticals_,horizontals_]:=
Module[{speciesall,ncolumnlist,valueslist,edgelist,nspeciesrow,eid,Nspecies,Nvalues,vcoord,ps,es,Rln10,F,lS,con,speciesnull,pid,id,G,comments},
(* Defining variables and lists *)
valueslist=Join[verticals,horizontals];
nspeciesrow=Length[specieslist]; 
ncolumnlist=Table[Length[specieslist[[i]]],{i,1,nspeciesrow}]; 
lS=Length[species];
speciesnull={};
For[i=1,i<=lS,i++,
If[species[[i]]==Null,AppendTo[speciesnull,species[[i]]]
];
];
speciesall={};
For[i=1,i<=lS,i++,
If[MemberQ[speciesnull,species[[i]]], Null, AppendTo[speciesall,species[[i]]]
];
];

Nspecies=Length[speciesall];
Nvalues=Length[valueslist]; 
edgelist=Table[valueslist[[i,2]] \[UndirectedEdge] valueslist[[i,3]],{i,1,Nvalues}];

eid={};
pid={};
For[i=1,i<=Nvalues,i++,
If[MemberQ[horizontals,valueslist[[i]]],AppendTo[eid,1];AppendTo[pid,0],
AppendTo[eid,0];AppendTo[pid,1]]
];

id={};
(* Generating the  list of edgelabels*)
For[i=1,i<=Nvalues,i++,
If[eid[[i]]== 1,AppendTo[id,"E0 = "],AppendTo[id,"pKa = "]];
];
(* Generating the Vertex Coordinates*) 
vcoord={};
For[i=1,i<= nspeciesrow,i++,
For[j=1,j<= ncolumnlist[[i]],j++,
If[specieslist[[i,j]]===0,Null,AppendTo[vcoord,{(j-1),-(i-1)}]];
];
];
(* Make a graph from given input *)
G=Graph[speciesall,edgelist];
comments={};
For[i=1,i<=Length[speciesall],i++,
If[Length[FindShortestPath[G,speciesall[[1]],speciesall[[i]]]]===0,AppendTo[comments,speciesall[[i]]]
];];
If[Length[comments]>0,AppendTo[comments,"Isolated Species that need more data to be included in the diagram:"]];
MatrixForm[Reverse[comments]]
];

printdata[specieslist_,species_,verticals_,horizontals_,ne_,T_,pHmin_,pHmax_,E0min_,E0max_,
OptionsPattern[{RefIndex->1,Ref-> "Fc+/Fc",Solvent->  "MeCN",Tol->  0.0,ColorScheme->  "Rainbow"}]]:= 
Module[{speciesall,ncolumnlist,valueslist,edgelist,nspeciesrow,eid,Nspecies,Nvalues,vcoord,ps,es,Rln10,F,lS,con,speciesnull,
path,dGi,speciesLabels,nS,id,valId,dG,G,Stable,uniqueTriples,triplePoints,sol,solpH,solE0,ln10,
edgePoints,uniquePairs,regions,point3,tLabel,regionPoints,uniqueRegions,cornerPoints,notEmpty,ref,solvent,tol,mm,
writefiles,colorscheme,vertexFile,vertexChannel,regionFile,regionChannel,a,b,c,i,j,k,iT,iP,ri,data},

(*Constants*)
F=23.0605488670097;
ln10=Log[10.0];
Rln10=0.00198720414566658*ln10;
(* options *)
ref=OptionValue[Ref];
solvent=OptionValue[Solvent];
tol=OptionValue[Tol];
writefiles=OptionValue[WriteFiles];
colorscheme=OptionValue[ColorScheme];
ri=OptionValue[RefIndex];
(* Defining variables and lists *)
valueslist=Join[verticals,horizontals];
nspeciesrow=Length[specieslist]; 
ncolumnlist=Table[Length[specieslist[[i]]],{i,1,nspeciesrow}]; 
lS=Length[species];
speciesnull={};
For[i=1,i<=lS,i++,
If[species[[i]]==Null,AppendTo[speciesnull,species[[i]]]
];
];
speciesall={};
For[i=1,i<=lS,i++,
If[MemberQ[speciesnull,species[[i]]], Null, AppendTo[speciesall,species[[i]]]
];
];

Nspecies=Length[speciesall];
Nvalues=Length[valueslist]; 
edgelist=Table[valueslist[[i,2]] \[UndirectedEdge] valueslist[[i,3]],{i,1,Nvalues}];

eid={};
pid={};
For[i=1,i<=Nvalues,i++,
If[MemberQ[horizontals,valueslist[[i]]],AppendTo[eid,1];AppendTo[pid,0],
AppendTo[eid,0];AppendTo[pid,1]]
];


For[i=1,i<=nspeciesrow,i++,
For[j=1,j<=Length[specieslist[[i]]],j++,
For[k=1,k<=lS,k++,
If[MemberQ[specieslist[[i]],species[[k]]],ps[species[[k]]]=i];
If[MemberQ[specieslist[[All,j]],species[[k]]],es[species[[k]]]=j];
]]];




G=Graph[speciesall,edgelist,VertexLabels->Placed["Name",Center]];
dGi[i_,pH_,E0_]:= If[eid[[i]]== 1,If[es[valueslist[[i,2]]]<es[valueslist[[i,3]]],-ne*F*(valueslist[[i,1]]-E0),ne*F*(valueslist[[i,1]]-E0)],
If[ps[valueslist[[i,2]]]<ps[valueslist[[i,3]]],Rln10*T*(pH-valueslist[[i,1]]),-Rln10*T*(pH-valueslist[[i,1]])]];

speciesLabels={};
For[i=1,i<=  Length[speciesall],i++,
If[Length[FindShortestPath[G,speciesall[[ri]],speciesall[[i]]]]!= 0 ,AppendTo[speciesLabels,speciesall[[i]]] 
];
];

nS=Length[speciesLabels];
path[i_]:=FindShortestPath[G,speciesLabels[[ri]],speciesLabels[[i]]];
For[i=1,i<=nS,i++,
For[j=1,j<Length[path[i]] ,j++,
For[k=1,k<=Nvalues,k++,
If[path[i][[j]]\[UndirectedEdge] path[i][[j+1]]===edgelist[[k]] \[Or] path[i][[j+1]]\[UndirectedEdge] path[i][[j]]=== edgelist[[k]] ,valId[i,j]=k]
];
];
];
dG[pH_,E0_]:=Table[Sum[dGi[valId[i,j],pH,E0] ,{j,1,Length[path[i]]-1}],{i,1,nS}];

(* This logical function defines a stability region for species i *)
Stable[pH_,E0_,i_]:=Ordering[dG[pH,E0],1][[1]]== i;
(* Find all triple points on the diagram *)
triplePoints={};
uniqueTriples=Subsets[Range[nS],{3}];
 For[iT=1,iT<Length[uniqueTriples]+1,iT++,
{i,j,k}=uniqueTriples[[iT]];
sol=Solve[dG[x,y][[i]]== dG[x,y][[j]]== dG[x,y][[k]]&&pHmin<x<pHmax&&E0min<y<E0max,{x,y}];
If[sol== {},Continue[]];
{solpH,solE0}=sol[[1,All,2]];
If[Stable[solpH,solE0,i]\[Or]Stable[solpH,solE0,j]\[Or]Stable[solpH,solE0,k],
AppendTo[triplePoints,{{solpH,solE0},speciesLabels[[{i,j,k}]]}];
];
];
(*
*)
(* Find all edge points on the diagram *)
edgePoints={};
uniquePairs=Subsets[Range[nS],{2}];

(* Bottom edge of the frame *)
For[iP=1,iP<Length[uniquePairs]+1,iP++,
{i,j}=uniquePairs[[iP]];
sol=Solve[dG[x,E0min][[i]]== dG[x,E0min][[j]]&&pHmin<x<pHmax,x];
If[sol== {},Continue[]];
solpH=sol[[1,1,2]];
If[Stable[solpH,E0min,i]\[Or]Stable[solpH,E0min,j],
AppendTo[edgePoints,{{solpH,E0min},speciesLabels[[{i,j}]]}]
];
];

(* Top edge of the frame *)
For[iP=1,iP<Length[uniquePairs]+1,iP++,
{i,j}=uniquePairs[[iP]];
sol=Solve[dG[x,E0max][[i]]== dG[x,E0max][[j]]&&pHmin<x<pHmax,x];
If[sol== {},Continue[]];
solpH=sol[[1,1,2]];
If[Stable[solpH,E0max,i]\[Or]Stable[solpH,E0max,j],
AppendTo[edgePoints,{{solpH,E0max},speciesLabels[[{i,j}]]}]
];
];

(* Left edge of the frame *)
For[iP=1,iP<Length[uniquePairs]+1,iP++,
{i,j}=uniquePairs[[iP]];
sol=Solve[dG[pHmin,y][[i]]== dG[pHmin,y][[j]]&&E0min<y<E0max,y];
If[sol== {},Continue[]];
solE0=sol[[1,1,2]];
If[Stable[pHmin,solE0,i]\[Or]Stable[pHmin,solE0,j],
AppendTo[edgePoints,{{pHmin,solE0},speciesLabels[[{i,j}]]}]
];
];

(* Right edge of the frame *)
For[iP=1,iP<Length[uniquePairs]+1,iP++,
{i,j}=uniquePairs[[iP]];
sol=Solve[dG[pHmax,y][[i]]== dG[pHmax,y][[j]]&&E0min<y<E0max,y];
If[sol== {},Continue[]];
solE0=sol[[1,1,2]];
If[Stable[pHmax,solE0,i]\[Or]Stable[pHmax,solE0,j],
AppendTo[edgePoints,{{pHmax,solE0},speciesLabels[[{i,j}]]}]
];
];

(* Build stability regions as collection of vertices selected from triple points and edgepoints *)

regions={};

(* Loop over triple points *)
Do[
point3=triplePoints[[iTPoint]];
(* Loop over labels of each triple point *)
Do[
tLabel=point3[[2,iTLabel]];
regionPoints={};
(* loop over edge points *)
Do[
If[tLabel== edgePoints[[iP2]][[2,1]]\[Or]tLabel== edgePoints[[iP2]][[2,2]],
AppendTo[regionPoints,edgePoints[[iP2]][[1]]]
]
,{iP2,1,Length[edgePoints]}];
(* loop over triple points *)
Do[
If[tLabel== triplePoints[[iP3]][[2,1]]\[Or]tLabel== triplePoints[[iP3]][[2,2]]\[Or]tLabel== triplePoints[[iP3]][[2,3]],
AppendTo[regionPoints,triplePoints[[iP3]][[1]]]
]
,{iP3,1,Length[triplePoints]}];
AppendTo[regions,{regionPoints,tLabel}]
,{iTLabel,1,3}]
,{iTPoint,1,Length[triplePoints]}];

(* Get rid of duplicate regions *)
uniqueRegions=Union[regions,SameTest-> (#1[[2]]==#2[[2]]&)]; 

(* Add corner points *)
cornerPoints={{pHmin,E0min},{pHmin,E0max},{pHmax,E0max},{pHmax,E0min}};
(* loop over four corner points *)
Do[
(* loop over unique regions *)
Do[
regionPoints=uniqueRegions[[iRegion,1]];
If[MemberQ[regionPoints[[All,1]],cornerPoints[[iCorner,1]]]&&MemberQ[regionPoints[[All,2]],cornerPoints[[iCorner,2]]],
uniqueRegions=Insert[uniqueRegions,cornerPoints[[iCorner]],{iRegion,1,-1}]
]
,{iRegion,1,Length[uniqueRegions]}]
,{iCorner,1,4}];

(*Print["\nUNIQUE REGIONS"];*)
(* Scan[Print,Sort[uniqueRegions]]; *)
(*Print["\nAREAS OF UNIQUE REGIONS:"];
Scan[Print,Sort[Table[{ConvexHullArea[uniqueRegions[[i,1]][[ConvexHull[uniqueRegions[[i,1]]]]]],uniqueRegions[[i,2]]},{i,1,Length[uniqueRegions]}]]];
*)
(* Mark the empty regions *)
notEmpty=Table[MemberQ[uniqueRegions[[All,2]],speciesLabels[[i]]],{i,1,nS}];
(*Print["Species not having stability regions in the given range of pH and \!\(\*SuperscriptBox[\(E\), \(0\)]\)"];
For[i=1,i<= nS,i++, 
If[notEmpty[[i]],Null,Print[speciesLabels[[i]]];]];*)

data={"=============== Coordinate Data ===============","---- Triple Points"};
For[i=1,i<=Length[triplePoints],i++,
AppendTo[data,triplePoints[[i]]];
];
AppendTo[data,"---- Edge Points"];
For[i=1,i<=Length[edgePoints],i++,
AppendTo[data,edgePoints[[i]]];
];
AppendTo[data,"---- Region Points"];
For[i=1,i<=Length[regionPoints],i++,
AppendTo[data,uniqueRegions[[i]]];
];
TableForm[data,TableSpacing->{1,2}]
];


End[]

EndPackage[]


(* ::Input:: *)
(**)
