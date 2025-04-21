# RCS

This code was developed to perform regional curve standardization on geoduck (<em>Panopea generosa</em>)

## Features

*   Build a custom regional curve based on the rwl and po (pith offset) files
*   RCS detrending

-----

## Requirements


R, the language, is availalble from CRAN. R Studio is an IDE that will work using R and make your workflow much easier.

**R**

https://cran.r-project.org

 **R Studio**

 https://www.rstudio.com

-------------


## Installation

Create a new project in R Studio and start with a fresh workspace:

![Workspace](https://www.dropbox.com/s/wg2w9ag7hwi7knd/1_fresh.png?raw=1)

Install _remotes_ in the console window:
 
    install.packages("remotes")

Use _remotes_ to install the RCS package from github:

    remotes::install_github("DaveEdge1/RCS")

Load the _RCS_ package:

    library(RCS)

And that's it! You have successfully installed RCS and are ready to start working.

## Core functions

### **robustRC()**

**Purpose:** 
Created a regional curve by 1) estimating the average ontogenetic growth and 2) fitting a curve to the average values. Curve fitting is by dplR methods or a custom age-varying spline method, ```tvSpline()```. 

**Parameters:**
> **rwlFile**
> 
> ring width file as created by read.rwl in dplR
>
> **poFile**
> 
> two-column file created by dplR
>
> **truncRC**
>
> limit curve calculation to regions of a designated minimum sample depth
>
> **tvSpline**
>
> TRUE/FALSE, use an age-varying spline (recommended fit method, dMethod must be set to NULL)
>
> **tvRange**
>
> min and max spline lengths in nyrs (if using tvSpline), eg. c(3,80)
>
> **tvStiff**
>
> set the age at which the spline length will max out (if using tvSpline)

**Returns:**
> 
> **RC**
> 
> regional curve by age
>
> **rcInfo**
>
> additional details from RC creation
>
> **AgeAligned**
>
> the growth increment series aligned by ontogenetic age
>
> **plot**
>
> visualization of RC creation

### **newRCS()**

**Purpose:** 
Detrend the growth increment series using the custom RC created from **robustRC**

**Parameters:**
> **rwlFile**
> 
> ring width file as created by read.rwl in dplR
> 
> **poFile**
> 
> two-column file created by dplR
>
> **truncRC**
>
> limit curve calculation to regions of a designated minimum sample depth
>
> **ratios**
>
> T/F, detrend by division
>
> **rcIn**
>
> custom regional curve, the ```RC``` output from ```robustRC()```
>
> ageMin
>
> minimum ontogentic age of samples for inclusion in chronology
>
> ageMax
>
> maximum ontogentic age of samples for inclusion in chronology

If you use this code, please cite the original publication using the citation in the right-hand margin.
  number={9},
  pages={e2021PA004291},
  year={2021},
  publisher={Wiley Online Library}
}
