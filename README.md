Welcome to Genome Size Evolution project!
==============

Scientific motivation:
-------------------

This is an implementation of a model presented in the paper:

Bentkowski P, Van Oosterhout C, Mock T (2015) *A model of genome size
evolution for prokaryotes in stable and fluctuating environments*.
Genome Biol Evol evv148. http://dx.doi.org/10.1093/gbe/evv148 
(it's Open Access)

The "Master" branch includes the newest stable development. For the
ORIGINAL code from the publication check this branch:

   **As_in_Bentkowski_et_al_2015_GBE**

 Software used:
-------------------

For writing this project I have used:

* NetBeans for coding, https://netbeans.org/
* Doxygen and Doxywizard for generating documentation:
      http://www.stack.nl/~dimitri/doxygen/

Compilation:
----------

To compile one needs GNU Scientific Library, http://www.gnu.org/software/gsl/ ,
installed. It's in any GNU/Linux repos. To read the documentation you need to
download the thing first.

In the directory '00_CurrentCode' one can find nbproject directory. If you running
NetBeans as well, you can load a NetBeans project from there.

Documentation:
-------------

Directory 'Doc' contains documentation in HTML format. Find file 'index.htm' and
run it in a browser to view the docs. Do not alter files in Doc manually, as they
are generated automatically by Doxygen. File 'Doxyfile' contains the set-up for
Doxygen (you mignt need to tweek some dir paths at your local workstation -
Doxywizard is a GUI front-end to Doxygen one might find useful).

Fairly accurate description what the model does can be found in the documentation
(directory 'Doc') or in Chapter 2 " Construction of the model and testing the
parameter space" of Piotr Bentkowski's PhD thesis:
https://ueaeprints.uea.ac.uk/50553/

Auxiliary code:
-----------

Directory '01_PyScrips' contains Python scripts used for plotting and analysis
of the output files. They are rather poorly documented - sorry :-(

* SigleRuns - scripts for processing single simulation.
* MultiRuns - scripts for processing and comparing a larger number of simulations,
often varying in parameter set-up between each other.
* GenotypeAnalysys - scripts for counting mutations and evolution related stuff.

Acknowledgements:
--------------

All simulations were carried out on the High Performance Computing Cluster supported
by the Research and Specialist Computing Support service at the University of East Anglia.
  http://rscs.uea.ac.uk/high-performance-computing