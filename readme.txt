
Installation:

```
mamba env create -n coffea-env python=3.10
mamba activate coffea-env
pip install coffea --upgrade --pre

```

Some naming conventions:

dType :
- E.g. data, signal, GJets, ...

dataset : 
- Subsets of a dType including the era/year
- E.g. data : (EGamma_2018A, EGamma_2018B, ...)
       GJets: (GJets_HT-100To200_2018, ...)
       signal: (BkkToGRadionToGGG_M1-140_R0-0p7_2018, ...)

The samples file contains datasets with their associated formats (MiniAOD, NanoAOD) and their access methods (das, vast, hadoop).