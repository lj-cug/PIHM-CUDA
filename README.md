# PIHM_CUDA
PIHM_CUDA is developed basing on the Serial and OpenMP versions of PIHM developed in  Multi-Modular Penn State Integrated Hydrologic Model (MM-PIHM), which is a physically based watershed model with multiple optional modules.
MM-PIHM is the **sweetest** PIHM, ever!

PIHM is a spatially-distributed, physically based hydrologic model.
PIHM-FBR adds fractured bedrock hydrology to PIHM to simulate deep groundwater processes.
Flux-PIHM adds a land surface model (adapted from the Noah land surface model) to PIHM for
the simulation of land surface processes.
Flux-PIHM-BGC couples Flux-PIHM with a terrestrial ecosystem model (adapted from Biome-BGC) that enables the simulation of carbon and nitrogen cycles.
The source code for the reactive transport (RT) module will be provided in future releases.

MM-PIHM is open source software licensed under the MIT License.
All bug reports and feature requests should be submitted using the [Issues](https://github.com/PSUmodeling/MM-PIHM/issues) page.

# Notes on CUDA version of PIHM (LI JIAN)

It should be noted that it's on heavy development stage. Now I use Visual Studio 2013 and Nvidia CUDA 9.5 to compile the host and device codes.
The CUDA version of SUNDIALS CVODE should be compiled first, then linked with the compiled host and device codes for PIHM_CUDA.

# Introduction to the directory in the repository

benchmark: The example files used for test the serial, openMP and CUDA C versions of MM-PIHM.

cvode: I used version 2.6.2, August 2015, the other versions should be tested. The library can be compiled using CMAKE, and have been configured in the sub-directory "build_vs";

doc: The related documentaries for PIHM and development recording.

src: The host and device codes for PIHM_CUDA.

util: The Shell Script for checking versions of PIHM.

pihm_cuda_lib: The Visual Studio 2013 Project to compile the device code of PIHM_CUDA, the shared library file will be generated.

pihm_cvode3.2.1_cuda: The Visual Studio 2013 Project to compile the cuda version PIHM.

pihm_cvode3.2.1_openmp: The Visual Studio 2013 Project to compile the openmp version PIHM.

pihm_cvode3.2.1_serial: The Visual Studio 2013 Project to compile the serial version PIHM.

# TODO
I need to edit a Makefile to compile all versions of CVDOE and PIHM in Linux OS.


# MM-PIHM Usage

# Installing CVODE

MM-PIHM uses the SUNDIALS CVODE v2.9.0 implicit solvers.
The CVODE Version 2.9.0 source code is provided with the MM-PIHM package for users' convenience.
SUNDIALS (:copyright:2012--2016) is copyrighted software produced at the Lawrence Livermore National Laboratory.
A SUNDIALS copyright note can be found in the `cvode` directory.

If you already have CVODE v2.9.0 installed, you can edit the Makefile and point `CVODE_PATH` to your CVODE directory.
Otherwise, you need to install CVODE before compiling MM-PIHM, by doing

```shell
$ make cvode
```

in your MM-PIHM directory.

Currently CMake (version 2.8.1 or higher) is the only supported method of CVODE installation.
If CMake is not available on your system, the CMake Version 3.7.2 binary for Linux (or Mac OS, depending on your OS) will be downloaded from [http://www.cmake.org](http://www.cmake.org) automatically when you choose to `make cvode`.

### Installing MM-PIHM

Once CVODE is installed, you can compile MM-PIHM models from the MM-PIHM directory by doing

```shell
$ make [model]
```

The `[model]` should be replaced by the name of model that you want to compile, which could be `pihm`, `pihm-fbr`, `flux-pihm`, or `flux-pihm-bgc`.

The command

```shell
$ make clean
```

will clean the executables and object files.

Note: If you want to switch from one MM-PIHM model to another one, you must `make clean` first.

A help message will appear if you run `make`.

#### Installation options

By default, MM-PIHM is paralleled using OpenMP, which significantly improves the computational efficiency of MM-PIHM models, especially Flux-PIHM and Flux-PIHM-BGC.
CVODE, however, is not implemented using OpenMP by default.
According to CVODE document, CVODE state variables (i.e., unknowns) "should be of length at least 100, 000 before the overhead associated with creating and using the threads is made up by the parallelism in the vector calculations".
In other words, you should using OpenMP for CVODE if your model domain has about 30, 000 or more model grids.
If you do want to test using OpenMP for CVODE, you can compile MM-PIHM models using

```shell
$ make CVODE_OMP=on [model]
```

Note that in order to use OpenMP for CVODE, you also need to turn on the OPENMP_ENABLE option when using CMake to install CVODE.

You can also turn off OpenMP for MM-PIHM (NOT RECOMMENDED):

```shell
$ make OMP=off [model]
```

By default, PIHM compilation is optimized using the `-O2` gcc option.
If you wish to debug using gdb, you may want to use

```shell
$ make DEBUG=on [model]
```

which will compile using `-O0` gcc option.

### Running MM-PIHM

#### Setting up OpenMP environment

To optimize PIHM efficiency, you need to set the number of threads in OpenMP.
For example, in command line

```shell
$ export OMP_NUM_THREADS=20
```

The command above will enable MM-PIHM model simulations using twenty (20) OpenMP threads.

If you use a PBS script, you must require the right number of ppn (processor cores per node) before setting the number of threads.
The ppn should be the same as the number of threads you want to use.
For example, your PBS script may look like

```shell
#PBS -l nodes=1:ppn=8
#PBS -l walltime=1:00:00
#PBS -j oe
#PBS -l pmem=1gb

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=8
./pihm example
```

#### Running MM-PIHM models

Now you can run MM-PIHM models using:

```shell
$ ./[model] [-c] [-d] [-t] [-V] [-v] [-o dir_name] [project]
```

where `[model]` is the installed executable, `[project]` is the name of the project, and `[-cdotVv]` are optional parameters.

The optional `-c` parameter will turn on the elevation correction mode.
Surface elevation of all model grids will be checked, and changed if needed before simulation, to avoid surface sinks.

The optional `-d` parameter will turn on the debug mode.
In debug mode, helpful information is displayed on screen and a CVODE log file will be produced.

The optional `-t` parameter will turn on Tecplot output.

The `-V` parameter will display model version.
Note that model will quit after displaying the version information.
No simulation will be performed when using the `-V` parameter.

The optional `-v` parameter will turn on the verbose mode.

The optional `-o` parameter will specify the name of directory to store model output.
All model output variables will be stored in the `output/dir_name` directory when `-o` option is used.
If `-o` parameter is not used, model output will be stored in a directory named after the project and the system time when the simulation is executed.

Example input files are provided with each release.
For a description of input files, please refer to the [*User's Guide*](https://github.com/PSUmodeling/MM-PIHM/releases/download/v0.11.0-alpha/guide.pdf "Guide")
 that can be downloaded from the release page.

### Penn State Users

The Penn State Lion-X clusters support both batch job submissions and interactive jobs.
The clusters for batch job submissions, Lion-XX clusters, usually support eight processors per node for OMP jobs.
The cluster for interactive jobs, hammer, supports twenty processors per node for OMP jobs, but is limited to short jobs only.