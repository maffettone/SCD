# Side Chain Decorator
In this case the goal in designing a side chain decoration algorithm was not to access the correct orientations of the side chains, but to create a sensible starting structure for further refinement. As clashing between chains and atomic overlap produces divergent potentials in MD, avoiding this common pitfall was a priority.
The side chain decoration program (SCD) was built using the modern Fortran standard and the same objects and methods utilized in PIPA. A modular design allows SCD to act as a standalone program that can be implemented in the workflow of other methods, as well as integrate directly into the PIPA output. 

The major steps are:
1. input of data and generating coordinates
2. appending side chains from the AMBER94 library
3. generating k- dimensional oriented polytopes (k-DOPs) which contain all rotamers
4. generating an adjacency matrix of feasible interactions
5. running Monte Carlo minimization using hard spheres potentials
6. outputting the results


## Used in combination with [PIPA](https://github.com/maffettone/PIPA)
Bayesian approach for information poor structure determination as highlighted by protein tertiary structure determination from an NMR prior and total scattering based likelihood. The algorithm is a Bayesian reverse Monte Carlo. 

The program, PIPA, has been developed in a proof-of-concept study to demonstrate its potential with the protein melittin. Using this protein with simulated data, I parameterize the program for use and discuss the utility of particular types of prior information. When considering an ensemble of structure refinements, the algorithm effectively discriminates sensible structure predictions from other ensemble members by using the $\chi^2$ parameter. In considering the models in isolation, the approach suffers from a lack of precision, but could potentially be used as a starting configuration in molecular dynamics. This opens the opportunity for inclusion of other information sources, while still demonstrating the pair distribution function as having sufficient information content to determine the overall tertiary structure of a protein. Combined with sparse chemical shift data, this work demonstrates a facile solution phase probe into protein structural solution, bringing the community one step closer to effective understanding of challenging proteins, especially membrane proteins. 

Acta Cryst. (2017), A70, C654
