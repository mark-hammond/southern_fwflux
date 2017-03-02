RUNOFF GENERATOR

This MATLAB function takes information about ice shelves and icebergs and produces a runoff file "runoff.bin".

Details of the data files and function are below. Example iceberg, ice shelf and river runoff data files are in this directory.

If you don't want the new river runoff (linked below, and used in runoff_NVND.bin) then just remove line 49:
"totalfluxmyr = totalfluxmyr + RiverRunoffGen(mask,DXC,DYC,XC,YC,cmap);"
or substitute your own river runoff (such as the upper portion of the old SOSE runoff)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DATA FORMAT

Ice shelf data format: Name | Basal Melt Rate in Gt/yr | Facing Angle | Opening Angle
Names must be accepted by the BEDMAP2 Toolbox for MATLAB
Facing angles are defined relative to the north-pointing vector
See "Facing Angle Demo.pdf" for a visual explanation.

River runoff data format:   No m2s_ratio lonm   latm   area(km2) Vol(km3/yr)  nyr  yrb yre  elev(m) CT CN River_Name OCN Station_Name
(from http://www.cgd.ucar.edu/cas/catalog/surface/dai-runoff/)

Iceberg tracks from http://www.scp.byu.edu/data/iceberg/database1.html should be placed in an "Iceberg Tracks" folder and have their extensions changed to ".dat" (blame MATLAB).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


FUNCTION DEFINITION

function [ totalfluxmyr totalvol a b c ] = RunoffFileGen_river( mask,DXC,DYC,XC,YC,smallfrac,basalvol,calvingvol,cmap)
%RunoffFileGen Produces a binary Antarctic freshwater runoff flux file
%   1) XC/YC: Grid cell centers
%   2) DXC/DYC: X/Y cell sizes
%   3) mask: Mask defining land (0) and ocean (1) on above grid
%   4) smallfrac: Fraction of "small" (normally distributed from coast
%   rather than following tracking database distribution) icebergs
%   5) basalvol/calvingvol: Total yearly volume of freshwater from basal
%   melting and iceberg calving, in Gt/yr
%   6) cmap: colourmap
 
% You'll need the "Bedmap2 ToolBox for Matlab" installed for ice shelf
% locations and shapes.
 
% This function should be run in a directory containing a folder 
% "Iceberg Tracks" with location files from 
% http://www.scp.byu.edu/data/iceberg/database1.html; I used QSCAT.
 
% The iceberg track files need to have their extensions changed to .dat 
% to work nicely with MATLAB's import function.
 
% It also need the riverrunoff.txt file (data from
% http://www.cgd.ucar.edu/cas/catalog/surface/dai-runoff/)

