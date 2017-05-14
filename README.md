# Stevens_and_Latimer_2015_GCB
Full workflow for: Stevens, J. T., and A. M. Latimer. 2015. Snowpack, fire, and forest disturbance: interactions affect montane invasions by non-native shrubs. Global Change Biology 21:2379-2393
R Code README 11/21/14

Workflow:

1) Run entire body of code in BroomStatModels.R; first for spp=“cs” and second for spp=“sj”. This produces parameter estimates and associated standard deviations for every treatment combination, which can then be used to run demographic models in CanopyProjections.R and FireProjections.R

If running Canopy analysis (e.g. to make MS Fig 3):

2) Run CanopyProjections.R Lines 1-212 (#1 - #4a) for spp=“cs” (to where the code says “STOP”)
3) Run CanopyProjections.R Lines 1-212 for spp=“sj”
4) Run CanopyProjections.R Lines 214-end to make figure

If running Fire analysis (e.g. to make MS Fig 4):

2) Run CanopyProjections.R Lines 1-98 (#1 - #2b) for spp=“cs”
3) Run FireProjections.R Lines 1-95 (#1 - #4) for spp=“cs”
4) Run CanopyProjections.R Lines 1-98 (#1 - #2b) for spp=“sj”
5) Run FireProjections.R Lines 1-95 (#1 - #4) for spp=“sj”
6) Run FireProjections.R Lines 99-end
