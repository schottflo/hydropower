# Project Hydropower

## Content
In Project Hydropower I co-developed a time series regression model for the discharge of an East African river to study profitability of hydropower projects and identify optimization 
potentials of the turbine set-up. In this repo I shared one general approach to model the discharge of any river (plus the corresponding electricity production).

## Requirements
The user needs to convert the discharge values of interest into a .rds file and place them in the working directory.

The user needs to acquire weather variables and save them as a .rds file. Those weather variables are meant to be used in the time series regression, thus they should be
associated with the river discharge to produce sensible results. The corresponding files need to start with "Extracted_" and need to be placed in the working directory.

## Scripts
Discharge_Analysis: The file gives an overview how to conduct the described analysis for any river. It is NOT runnable. It contains gaps or placeholders to make future analysis easier for the user.

Functions: The file which contains all functions used in Discharge_Analysis.

## Notes
Some functions are tailored to the specific setup of the problem given to my team and me. Be aware that the used weather variables were only available from 01.01.1979 to 31.12.2018 and
that the electricity production is subject to a specific power transformation function.

## Co-authors
Sven Kohler,
Alexander Arzt,
Martin Mächler
