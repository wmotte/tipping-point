# tipping-point
Determine tipping point in structural-functional connectivities.

## Study Background
Understanding how structural and functional networks are related will help to unravel the complex organization of brain networks. These networks are complex because the topological behavior is difficult to explain by the constituent parts, such as local connectivity or regional synchronization. A simple linear correlation analysis between structural and functional connectivity may not completely cover the complex nature of the structure-function relationship and ignores potential non-linear components. Therefore, we mapped the brain-wide structure-function relationship considering potential non-linear behavior using an additive model.

## Methods 
We determined structural connectivities with diffusion-weighted MRI and functional connectivities with resting-state fMRI, both the human and rat brain. To determine the structure-function relationship without an a priori determined linear model, we fitted generalized additive models in R [mgcv](https://www.rdocumentation.org/packages/mgcv), using five knots and the restricted maximum likelihood method with functional connectivity dependent on structural connectivity. Based on these generalized additive models, we proposed a novel measure of structure-function correlations, which we named the structure-function tipping point. To determine this tipping point, we calculated the first derivative of the generalized additive model and its corresponding 95% confidence interval based on 10,000 iterations. This first derivative represents the change in functional connectivity for each increasing step of structural connectivity. Structure-function tipping points are structural connectivity values where the increasing step of structural connectivity is associated with a change in functional connectivity that is significantly higher or lower than zero (the 95% confidence interval of the first derivative does not include zero).

## Script + example
As input the script needs structural and functional connectivity matrices. The provided script enables the user to plot individual datapoints of structural and functional connectivity, to plot the generalized additive model and its first derivative including the 95% confidence interval. The script also enables to determine the tipping point, the structural connectivity value from which structural and functional connectivity become significantly associated (the 95% confidence interval does not include 0). 

A non-linear fit between structural and functional connectivity in the human brain:  
![alt text](https://github.com/wmotte/tipping-point/raw/master/images/fit.png "Example of fit")

The individual datapoints are shown as blue circles; the blue line represents the generalized additive model fit and the orange line represents the first derivative including shades that represent its 95% confidence interval. The dotted vertical line represents the structure-function tipping point. 
