# Global_OCpetro_Oxidation

This repository contains the code and data for the manuscript "Global CO2 release by rock organic carbon oxidation negates silicate weathering" by Jesse R. Zondervan, Robert G. Hilton, Fiona Clubb, Mathieu Dellinger, Tobias Roylands, Mateja Ogrič.

This repository contains the code only and a link to the python environment and the data used to generate the final global OCpetro oxidation model.

The code was developed in a python environment detailed in anaconda-project.yml, with every recursive dependency down to the individual build in anaconda-project-locked.yml. The code file is "Glob_newmethod_parr_globalresidual.py" in this repository.

This code should be reproducible indefinitely, without depending on online package repositories, as long as nothing depends on system libraries. Both the commands and the environment have been captured into a fully locked anaconda project with all conda packages unpacked and included using anaconda-project --pack-envs. Note that only packages for running the code on Linux can be unpacked in this way; building on other platforms (Windows, Mac) will still require access to repositories.

A permanent DOI was created through Zenodo and can be found here: 

## Data
Data needed to run this code is uploaded and available on Zenodo. To run the code with available data in a functional environment, it is advised to download the files, code and environment from Zenodo. Anaconda-project can be used to unpack the zip file and to run the model. Instructions can be found by searching for anaconda-project online, or directly via https://anaconda-project.readthedocs.io/en/latest/user-guide/tasks/create-project-archive.html?highlight=unarchive#extracting-the-archive-file [last accessed 26/01/2023]

##Notes
This code was run on an HPC environment with a job submitter called SLURM. As such, the code will run according to a slurm job array with numbers from 1-10000 (10,000 monte carlo simulations). sbatch --array=1-10000:1 job_script_file_name.sh

###example:
sbatch --array=1-10000:1 run_Glob_OCpetro_model.sh

An example of a job script file has been appended. Please not that the details of this job script depend on your machine or HPC.
