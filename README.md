# MRIRelaxometryNLM
This repository contains the Matlab / Octave code of a nonlocal means filter for MRI relaxometry.
Its use is described in 
    L. Starke, K. Tabelow, T. Niendorf & A. Pohlmann,
    "Denoising for Improved Parametric MRI of the Kidney: Protocol for Nonlocal Means Filtering",
    Springer Protocols book "Preclinical MRI of the Kidney", (2021)

Contains:

- nlmFilter2D.m             : nlm filter for 2D
- stackNlmFilter.m          : nlm filter for a stack of relaxometry images (dim 1&2: spatial, dim 3: relaxometry)
- nlmExample.m              : example running nlmFilter2D.m and stackNlmFilter.m on synthetic data

For implementation details and references, see the above publication, especially notes 1 and 6.

Everything is licensed under GNU GPLv3 

Ludger Starke; Max Delbrück Center for Molecular Medicine in the Helmholtz Association, Berlin; 21-01-25

Contact me at Ludger.Starke@mdc-berlin.de for feedback or support

