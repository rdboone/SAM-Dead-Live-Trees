Read Me for Tree Mortality data prep
Larissa Yocom
Sept. 2017

MortalitySites.csv = File with site information

Original files that were modified to create MortalitySites.csv were from Maxime Caillaret and Katie Ireland: Treedata.csv and KatieAspenTreeData.csv.
Sites were selected where the 1st mortality factor was identified as drought.
Sites were numbered in order of decreasing longitude (roughly, eastern Europe to western North America).
116 sites are represented in the file.

Columns:
SiteNum		numbers 1-116, ordered by longitude in decreasing order
site		site name as given by original researcher
species		tree species
longitude	longitude
latitude	latitude
reference	text which is a combination of tree species, lead author name and date (style from Maxime Caillaret)



SpeciesTable.csv = File with species numbers (1-23) and names


RingWidths.csv = File with ring width information

Original files that were modified to create RingWidths.csv were from Maxime Caillaret and Katie Ireland: RWdata.csv and KatieRW.csv.
Cores were selected from only the sites in MortalitySites.csv (drought identified as #1 mortality factor).
Data was processed to give site number, tree number, calendar year, tree age, ring width, and year number for each record; 314,500 records total
Tree age was calculated differently for trees that came with an estimated age vs. trees with no estimated age. For trees with no estimated age,
the first tree ring was assigned age 1 and other tree rings followed in sequence up through the last measured tree ring. For trees with estimated ages,
the difference between the ring number of the last tree ring and the estimated age was added to the ring number of the first (1) through last tree ring. 
For example, tree x has tree rings 1-120. Estimated age is 127. The difference is 7. Tree ring 1 is assigned age 8, tree ring 120 is assigned age 127, etc. 

Columns:
SiteNum		site number, 1-116, matching numbers in MortalitySites.csv
TreeNum		tree number, 1-4262, in order of site so that first tree numbers are in site 1 and last tree numbers are in site 116
Year		calendar year, starting in 1905 because high-resolution gridded CRU data are available from 1901 (for SAM model to always have a few prior years)
Age		age for each tree ring (see above for explanation of how this was calculated)
Width		measured ring width
Species		tree species number (1-23)
YearNum		years numbered from 1-110 (1905-2014)