
# dragmine.py

Dragmine - Python 3 program which mines repetitive sequences (spidroins) from protein genome datasets. 

## SYNOPSIS

networkCreate.py mandatory arguments: [-input proteinFile.fa ]  
  optional arguments: [-out results] [-nmer int] [-scramble y or n] [-chop y or n] [-order e n or f]

## DESCRIPTION

Designed to mine for spider silk proteins (spidroins), Dragmine has been designed to rank protein sequences in order of their repetitiveness. This is achieved by first generating an n-mer (k-mer) profile of each sequence. The n-mer profile is then used to calculate an evenness score (based on the Shannon's entropy equation) for each protein, where highly repetitive protein sequences receive the lowest scores. These results are output to a csv table, which can be directly interpreted or used for further statistical analysis.

## VERSIONS

dragmine.py Last Updated 28 August 2019

## AUTHORS

Charles Bannister

	
