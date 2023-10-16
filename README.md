# Relative_volumes_tripartite
Supplementar material to our work. This repository contains codes that can be used to reproduce the results found in the paper
"Estimating the volumes of correlation sets in causal networks".

There we are interested in using non-convex optimization techniques by means of the Gurobi package as well as convex ones
to estimate the volume occupied by typical correlation sets in tripartite causal networks. In this last case, we make use 
of the Inflation package for Python with its underlying Navascués-Pironio-Acín hierarchy and other methods that can be 
found in the literature, in particular the employment of Finner and Entropic inequalities. Below one will find a short 
guide on the libraries at disposition in this repository.

## gurobi_codes
This one includes functions that can be used to solve optimization tasks with quadratic constraints for the Bilocal and 
Evans scenarios. Given a behaviour "p" in the form of a list with the appropriate number of entries depending on the 
scenario, functions are going to return feasibility status and distances when applicable upon request. Below is a list 
of the functions one will find there:

b_compatibility(p,args): used to test locality and bilocality of a given behaviour p for the Bilocal scenario;

e_compatiblity(p,args): used to test compatibility of a given behaviour p with the Evans scenario, including 
interventions if requested;

e_NSI_test(p): test whether p for the Evans scenario is compatible with NSI conditions.

## inflation_codes
There is a single function here to test the compatibility of a given behaviour "p" with the (inflated) scenario of interest, 
among the possibilities: "bilocal", "evans" and "triangle", by means of a direct application of the Inflation package: 
inflation_test(p,args)

## inequalities_test
To be added. It contains functions to assess whether a given triangle behaviour "p" satisfies the inequality constraints of 
interest.

finner_test(p): compatibility with the Finner inequality;

entropic_test(p): compatibility with the entropic inequalities at disposal in the literature for this scenario.
