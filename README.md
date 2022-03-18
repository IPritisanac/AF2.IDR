# AF2.IDR
Code and data for "Systematic identification of conditionally folded intrinsically disordered regions by AlphaFold2"

## (1) System Requirements
The collection of Python scripts were written in Python 3.8.8 and do not require external packages beyond those that are included in The Python Standard Library (https://docs.python.org/3/library/). In a couple of instances, the external Python package BioPython (https://biopython.org/) is used, as well as the external software programs DSSP and MAFFT, as indicated in the header of the file containing the import statements. The scripts were run on macOS and Windows platforms and do not require specific hardware features. 

## (2) Installation guide
Installation requires a standard Python 3 installation, e.g.: https://docs.anaconda.com/anaconda/install/index.html. Once Python 3 has been installed, each script can be run as outlined in the README files. For the Python scripts that require external libraries/software, installation of the following libraries can be achieved as follows. DSSP is most easily installed via Anaconda on macOS and Linux systems (https://anaconda.org/salilab/dssp) and MAFFT can be installed on macOS, Linux, and Windows platforms (https://mafft.cbrc.jp/alignment/software/). For the one script that requires DSSP (read_AF_DSSP.py), we created an alternative script that runs without this feature (read_AF.py). The intsallation time is less than 5 minutes.

## (3) Demo
Please see each of the sub-directories and the README files for explanations on how to run the demo scripts and recapitulate the results shown in the preprint. The longest run time is for read_AF.py in FIG1, which takes ~20 minutes on a Windows laptop with an Intel Core i5-8250U CPU at 1.60 GHz and 8 GB of RAM. On a Linux system that has an i9 core with 64 GB of RAM, the run time is reduced to approximately 2 minutes.

## (4) Instructions on how to use
Please see the README files for explanations on how to use the software and reproduce the plots in the preprint. 
