# Global_OCpetro_Oxidation

This repository contains the code and data for the manuscript "Global CO2 release by rock organic carbon oxidation negates silicate weathering" by Jesse R. Zondervan, Robert G. Hilton, Fiona Clubb, Mathieu Dellinger, Tobias Roylands, Mateja Ogrič.

This repository contains the code only and a link to the python environment and data used to generate the final global OCpetro oxidation model. In addition, for transparency it contains code used to generate the submodels of denudation and OCpetro stocks which are used as inputs for the OCpetro oxidation model. The outputs of these submodels are available on Zenodo as data inputs for the OCpetro oxidation model, which means there is no need to rerun this code when attempting to run the final model.

The code was developed in a python environment detailed in anaconda-project.yml, with every recursive dependency down to the individual build in anaconda-project-locked.yml. The code file is "Glob_newmethod_parr_globalresidual.py" in this repository.

This code should be reproducible indefinitely, without depending on online package repositories. Both the commands and the environment have been captured into a fully locked anaconda project with all conda packages unpacked and included using anaconda-project --pack-envs. Note that only packages for running the code on Linux can be unpacked in this way; building on other platforms (Windows, Mac) will still require access to repositories.

A permanent DOI for this Github entry was created through Zenodo.

[![DOI](https://zenodo.org/badge/592910202.svg)](https://zenodo.org/badge/latestdoi/592910202)

## Data
Data needed to run this code is uploaded and available on Zenodo. To run the code with available data in a functional environment, it is advised to download the files, code and environment from Zenodo. Anaconda-project can be used to unpack the zip file and to run the model. Instructions can be found by searching for anaconda-project online, or directly via https://anaconda-project.readthedocs.io/en/latest/user-guide/tasks/create-project-archive.html?highlight=unarchive#extracting-the-archive-file [last accessed 26/01/2023]

## Notes
This code was run on an HPC environment with a job submitter called SLURM. As such, the code will run according to a slurm job array with numbers from 1-100 (10,000 monte carlo simulations). The command to run the Monte Carlo simulation as 10,000 seperate jobs is done like this: "sbatch --array=1-10000:1 job_script_file_name.sh". Note that the version of code uploaded here is set to run 100 simulations ("sbatch --array=1-100:1 job_script_file_name.sh"). When running 10,000 simulations, please change line 46 to "quantile = float(os.getenv('SLURM_ARRAY_TASK_ID'))/10000."

Whilst it is possible to run this code on a single machine such as a personal computer, the user is warned that it takes 24 core hours per simulation to run. For example, a typical 4-core laptop would need 6 hours to run one simulation. Now calculate how many 10,000 would take...

To run one simulation, line 46 ("quantile = float(os.getenv('SLURM_ARRAY_TASK_ID'))/100.") can be replaced with ("quantile = float(number between 0 and 1)"). Outputs will be saved for each simulation, which can rack up a lot of space, unless you specifically put in lines to delete these from the disk, or, in the case of the example job script for HPC usage, exclude the files when moving data from the node that ran the job.

### Example:

sbatch --array=1-100:1 run_Glob_OCpetro_model.sh    #note that this runs 100 simulations.

An example of a job script file has been appended. Please note that the details of this job script depend on your machine or HPC system. Please consult your HPC support or platform's (Linux, Mac, Windows) command prompt instructions.
