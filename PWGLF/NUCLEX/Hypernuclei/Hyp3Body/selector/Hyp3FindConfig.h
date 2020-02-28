#ifndef HYP3FINDCONFIG_H
#define HYP3FINDCONFIG_H
//number of bin in ct
const int kNBinCt = 10;
//number of bin in pt
const int kNBinPt = 10;
//number of variables for the cuts
const int kNcuts = 3;
//number of variations for each cut
const int kNvariations = 5;
//all the sets of cut [variables of the cut][cut set][particle[deuton,proton,pion]]
const float kCuts[kNcuts][kNvariations][3] = {{{3,3,3},{3.5,3.5,3.5},{4.5,4.5,4.5},{5.5,5.5,5.5},{6.5,6.5,6.5}},
{{70,70,70},{70,70,70},{75,75,75},{125,125,125},{150,150,150}},
{{0,0,0},{1,1,1},{2,2,2},{3,3,3},{4,4,4}}};
//chi2 max for chi_2progs,chi3_prong,chi2_topology
const float kMaxKFchi2[3] = {5,5,5};
#endif
