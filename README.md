# GSFLOW-GRASS

***Python toolkit using GRASS GIS to generate inputs for the USGS hydrologic model "GSFLOW".***

GSFLOW is the US Geological Survey's combined groundwater (MODFLOW) and surface water (PRMS) model. GSFLOW-GRASS builds a stream network and its "hydrologic response units" -- here, watershed sub-basins -- that combine to form a surface-water hydrologic network (PRMS). It then links these to a regular finite difference grid to compute groundwater flow using MODFLOW.

These instructions are meant to take a user familiar with computers but new to (or a beginner with) GSFLOW and GRASS GIS through the basics of how to get GSFLOW-GRASS to work. *Please leave a message if you have trouble working with GSFLOW-GRASS; your comments could assist both you and the more general improvement of this documentation.*

When you use GSFLOW-GRASS, please contact us; a reference (Ng et al.) will be noted here as soon as the paper is made available.

This manual is written in the style of a quick(-ish) start guide that allows users to get GSFLOW up and running for their target watershed using mostly default options in our toolkit and without modifying GRASS-GIS or Python scripts. Additional customization is possible by more advanced users by editing our GRASS GIS and input file building scripts to extend the set of parameters incorporated **\todo{Andy: make section on adding GRASS parameters; perhaps just a small part at the end with appropriate links}**.

## Required Software

* **GSFLOW 1.2**
* **GRASS GIS 7.3+** and extensions (described below); **7.4** is the stable version at the time of writing
* **Python 2.7.X**
* **GSFLOW-GRASS Toolkit** (this software)

### Installing GSFLOW

***Download and install GSFLOW 1.2.0 or 1.2.2***

Obtain the source code from:
https://github.com/UMN-Hydro/GSFLOW-1.2.0
compile and install it. For windows, you can also download the executable file already compiled from the USGS website: https://water.usgs.gov/ogw/gsflow/#downloads.

***Hoping to have a better integration with GSFLOW v1.2.2, so not writing much more in the way of instructions! 
May need tp remove link to USGS website if 1.2.2 is not available and we don't want to point people towards 1.2.1***

### Installing GRASS GIS

***Download and install GRASS GIS 7.3+***

Two options:
* Cross-platform binaries:
https://grass.osgeo.org/download/software/
* Instructions to build from source:
  * https://grasswiki.osgeo.org/wiki/Compile_and_Install
  * https://grasswiki.osgeo.org/wiki/Compile_and_Install_Ubuntu

If you choose to compile GRASS GIS from source, we have used these configuration flags many times on Ubuntu (`configure_ubuntu.sh`):

```configure
CFLAGS="-O2 -Wall" LDFLAGS="-s" ./configure \
--enable-largefile=yes \
--with-nls \
--with-cxx \
--with-readline \
--with-pthread \
--with-proj-share=/usr/share/proj \
--with-geos=/usr/bin/geos-config \
--with-wxwidgets \
--with-cairo \
--with-opengl-libs=/usr/include/GL \
--with-freetype=yes --with-freetype-includes="/usr/include/freetype2/" \
--with-postgres=yes --with-postgres-includes="/usr/include/postgresql" \
--with-sqlite=yes \
--with-mysql=yes --with-mysql-includes="/usr/include/mysql" \
--with-odbc=no \
--with-netcdf=/usr/bin/nc-config
```

### Installing Python

GSFLOW-GRASS has been tested on **Python 2.7**, and should work (with a few possible changes) on future versions of Python 2.X. It has not been tested on Python 3.X.

In order to run properly, GSFLOW-GRASS requires the following Python modules (in addition to those that come with GRASS GIS v7.3 or greater, above):
* `numpy`
* `matplotlib`
* `pandas`
* `osgeo` (also listed under `gdal`)

*For users who are new to Python, follow these directions to install the Python interpreters onto your computer.*

#### Linux

Use your package manager to download and install the required Python packages. For Debian/Ubuntu, it will be something like:

```bash
# Basic packages
sudo apt-get install \
python python-numpy python-scipy \
python-setuptools python-matplotlib \
python-gdal

# pip (recommended for automatic installs via setuptools)
sudo apt-get install python-pip

# iPython console -- very useful (optional)
sudo apt-get install ipython

# Sypder IDE (I don't personally use it but many others like it: optional)
sudo apt-get install spyder
```

#### Windows and Mac

We recommend using **Anaconda**, which comes with most of the Python modules that you might need for the execusion of the GSFLOW-GRASS codes (including numpy, matplotlib and pandas but not including osgeo). Additional modules may be installed with either "conda" (the Anaconda package manager) or "pip" (the Python package manager), for example:

```bash
# Anaconda
conda install python-numpy
# Pip
pip install numpy
```
For Windows, you need to use conda through the Anaconda Prompt (After installing Anaconda, go to the windows start menu, search programs and files for "Anaconda Prompt"), it will not work through the regular command prompt. Installing osgeo/gdal may also downgrade certain packages that you can re-upgrade afterwords. For example (in the Anaconda prompt):

```bash
# In the Anaconda Prompt
conda install gdal
conda upgrade numpy
```
Anaconda also comes with the Spyder development environment which is useful if you want to individually run or edit the python codes.

### Installing FFMPEG

(Optional)

In order to make movies, you need an active copy of ffmpeg on your computer. For Linux, this is simple:
```sh
sudo apt-get install ffmpeg
```
For Windows or Mac, download via https://www.ffmpeg.org/download.html. You can select the correct version for your operating system under “More downloading options” and “Get the packages” and follow the links therein. Ffmpeg is installed by adding the executable file to your system's or user path variable. Windows installation instructions can be found here https://www.wikihow.com/Install-FFmpeg-on-Windows.

## Directory Structure

**GSFLOW-GRASS** has four main directories:

* **domain_builder** holds GRASS GIS commands and associated code to turn a DEM into an input domain of Hydrologic Response Units (HRUs) and stream segments.
* **input_file_builder** holds code that turns the outputs from **domain_builder** and user-specified parameters in the **settings** file (see Step 1, below) into the set of input files that are required by GSFLOW
*  **Run** holds the files to run GSFLOW using these inputs
* **visualization** holds scripts to build plots and movies based on the outputs from both GSFLOW and the GRASS GIS domain builder

In addition, the "figures" directory holds images used for this README.

<!---

This seems like a lot of work to maintain when they could just look at the actual directory structure.

\todo{This isn't actual directory structure now.  Change after finalizing.  Will need to make sure everything still runs, i.e., can find dependent files etc}

* `GSFLOW-GRASS -> domain_builder` GRASS GIS Python scripts to build the watershed and MODFLOW geometries
* `GSFLOW-GRASS -> input_file_builder` Python scripts to create input files for GSFLOW, including the control file, paramers file for PRMS, and input file for MODFLOW (NWT)
  * `GSFLOW_print_controlfile_current.py`
  * `GSFLOW_print_PRMSparamfile_current.py`
  * `print_MODFLOW_inputs_res_NWT_current.py`
  * `-> MODFLOW_scripts` \todo{consider doing away with this sub-directory?}
    * `MODFLOW_NWT_lib_current.py` \todo{modify to read in down-gradient cell of discharge point}
* `Create_hydcond_array.py`
* `GSFLOW-GRASS -> Run ->`
  * `settings_template.ini`
  * `settings_test.py`
  * `run_Python_GSFLOW_current.sh` **\todo{re-name as go-GSFLOW.sh}**
  * `run_Python_GSFLOW_current_Win.bat` **\todo{re-name as go-GSFLOW.bat}**
* `GSFLOW-GRASS -> visualization`
  * `Plot_MODFLOW_inputs.py`
  * `plot_gsflow_csv.py`
  * `gsflow_csv_table.py`
  * `Plot_MODFLOW_3D_head3.py`
  * `Plot_MODFLOW_3D_uzf.py`
\todo{Andy: put GRASS-GIS scripts for plotting HRU-distributed variables here?  Or above with domain-builder scripts?}
* `GSFLOW-GRASS -> examples`
  * `Shullcas`
  * `Santa Rosa` **\todo{space here?}**
-->

## Pre-processing
***Generating inputs to GSFLOW using Python and GRASS***

### Step 1: Customize the Settings File

The **Settings** file holds user-defined information that defines how GSFLOW is set up and will run.

Use `settings_template.ini` in the 'Run' folder as a template for creating your own **Settings** File, which can have any name. **Boldface** options are required. This files includes:

#### "settings" section

| Option             | Description
| ------------------ | ---------
| **proj_name**      | Label for the project; no spaces
| **gsflow_exe**     | Full pathname for GSFLOW executable
| **gsflow_path_simdir** | Full pathname for location where GSFLOW simulation<br>directory should go.
| **fl_print_climate_hru** | **1** to print spatially uniform climate data over all HRU's<br>using climate data from file specified in *climate_data_file*.<br>**0** if user already has pre-existing HRU-distributed climate files.
| climate_data_file  | Only for *fl_print_climate_hru*=1:<br>Name of file containing climate data for single weather station site,<br>to be uniformly distributed over all HRU's using<br>`GSFLOW_print_data_climatehru_files1_metric.py`.<br>If *fl_print_climate_hru*=0, this entry can be omitted;<br>if it is included anyway, it will be ignored.<br>
| climate_hru_dir    | **Only for *fl_print_climate_hru*=0**:<br>Name of directory with pre-existing climate_hru data files <br>containing HRU-distributed climate inputs:<br>**tmin.day**, **tmax.day**, **precip.day**, and **empty.day**.<br>See GSFLOW manual or example files in example cases<br>(e.g., in Shullcas -> inputs -> PRMS_GSFLOW) for format of climate_hru data files.<br>**If *fl_print_climate_hru*=1**, this entry can be omitted;<br>if it is included anyway, it will be ignored.<br>
| **sw_1spinup_2restart** | **1** for spin-up run starting from steady-state MODFLOW period<br>**2** for restart run starting from states saved in the below files
| restart_PRMSfil | optional: for restart runs (sw_1spinup_2restart=2)<br>full pathname of file that is saved under ``save_vars_to_file''<br>in the GSFLOW control file during a previous run.<br>This entry won't be used (but should still be entered) if *sw_1spinup_2restart*=1<br>for startup runs.
| restart_MODfil | optional: for restart runs (sw_1spinup_2restart=2)<br>Full pathname of file that is saved under ``IWRT'' in the MODFLOW name file<br>during a previous run. This entry won't be used (but should still be entered)<br>if sw_1spinup_2restart= 1 for startup runs

#### "custom_params" section

| Option             | Description
| ------------------ | ---------
| **fl_create_hydcond**  | **1** to implement Python script to create spatially distributed hydraulic conductivity.<br>**0** to use values or pre-existing file entered in *hydcond*<br>**\todo{Crystal: implement this in Create_hydcond.py}**
| **hydcond**            | For uniform hydraulic conductivity within each layer in the saturated domain:<br>enter value(s) (in [m/d]), using comma-separated list for multiple layer domains,<br>starting with top layer.  For spatially distributed values: Enter file name containing<br>array of values (in [m/d]); if *fl_create_hydcond*=1, contents of file will be created<br>using the Python script called in the Go-GSFLOW File (see below description of<br>Go-GSFLOW File).<br>User may need to adjust hydraulic conductivity values to reach numerically<br>convergent results and to match observations.
| finf | Optional: Only for spin-up runs; this entry is ignored (but should still be entered)<br>for restart runs. For uniform infiltration to the unsaturated zone over the watershed:<br>enter a single value (in [m/d]).  For spatially distributed values: enter file name<br>containing array of values (in [m/d]).[br] User may need to adjust this value to reach<br>numerically convergent results and for reasonable start of transient results.

#### "domain" section

| Option             | Description
| ------------------ | ---------
| **start_date**     | Start date of simulation, format: YYYY-MM-DD
| **end_date**       | End date of simulation, format: YYYY-MM-DD
| **NLAY**           | Integer number of vertical layers
| **DZ**             | Layer thicknesses (in [m]).  For multiple layers (NLAY>1), enter comma-separated list,<br>starting from top layer.

#### GRASS section

| Option             | Description
| ------------------ | ---------
| **DEM_file_path_to_import**   | DEM for import; if blank, assumes that DEM has already been imported<br>and that the associated initial calculations (flow accumulation, offmap<br> flow converted to NULL cells) have been completed.
| **threshold_drainage_area_meters2** | Threshold drainage area in square meters at which flow is considered to create<br>a channel
| **MODFLOW_grid_resolution_meters** | Target cell side length in meters for MODFLOW grid; side lengths will not be<br>exactly this long, as the nearest value to create an integer number<br>of cells in the domain will be chosen.
| **outlet_point_x** | Pour point approximate x (Easting) position; the nearest stream<br>segment to this point is chosen as the true pour point.
| **outlet_point_y** | Pour point approximate y (Northing) position; the nearest stream<br>segment to this point is chosen as the true pour point.
| **icalc** | Method selector for hydraulic geometry computation<br>**0:** Constant<br>**1:** Manning's equation with the wide channel assumption.<br>**2:** Manning's Equation.<br>**3:** Power-law relationship between width, depth, velocity,<br/>and discharge, per Leopold and Maddock (1953)
| gisdb | Optional: Directory that holds grass GIS locations.<br>Typically `~/grassdata`<br>Not currently used.<br>(Will be used to run this while starting GRASS in the background)
| version | Optional: GRASS GIS version number without any "." characters.<br>We used **73**<br>Option is not currently used.<br>(Will be used to run this while starting GRASS in the background)


### Step 2. Running GRASS GIS to build the model domain

#### Launch GRASS GIS and create your location

1. Launch GRASS GIS
2. Create a folder to hold your GRASS locations. This is typically called **"grassdata"** and placed in your home directory.
3. Click on "New", and follow the prompts.
  * We recommend naming the **"Project Location"** the same as **"proj_name"**.
  * No "Location Title" is needed.
  * We suggest that you **"Read projection and datum terms from a georeferenced data file"** to set the coordinate system. Each GRASS GIS location has only one coordinate system. Your DEM **must be in a projected coordinate system**: we do not test our codes using geographic (lat/lon) coordinate systems.

![GRASS GIS start-up screen](figures/GRASS_startup_screen.png)

##### Choosing pour points (if needed)

If you need to choose your pour point manually, we recommend that either you (a) import the DEM and find it now using GRASS, or (b) use another program like QGIS to find this location. You should still keep the path of the DEM in **DEM_file_path_to_import** until you have finished the first run of the code. Keeping this in during later runs will not cause problems; it will just take extra time as the DEM is re-imported and re-corrected.

#### Start GRASS GIS session in your location

With your location selected, click **"Start GRASS session"**.

#### Install required GRASS GIS extensions

***First time running GRASS GIS only, or after GRASS extension updates***

Either using the terminal or clicking on the "Console" tab in the GRASS GIS Layer Manager, run install all of the necessary extensions for GSFLOW-GRASS.

In the terminal, this may be done as follows:

On Linux, cd `GSFLOW-GRASS/domain-builder` and type:
```sh
sh install_extensions.sh
```

On Windows, likewise cd to `GSFLOW-GRASS\domain-builder` and type:
```bat
install_extensions
```

Perhaps more easily, on either operating system, just copy and paste the contents of `install_extensions.sh` into the terminal or GRASS GIS console. We have also pasted these here for convenience.
```sh
# From us
g.extension v.gsflow.export
g.extension v.gsflow.gravres
g.extension v.gsflow.grid
g.extension v.gsflow.hruparams
g.extension v.gsflow.reaches
g.extension v.gsflow.segments
g.extension r.gsflow.hydrodem
g.extension v.stream.inbasin
g.extension v.stream.network
g.extension r.cell.area

# From others
g.extension r.stream.basins
g.extension r.hydrodem
```

#### Create the domain inputs for input-file-builder scripts

Either using the terminal (Linux) or clicking on the "Console" tab in the GRASS GIS Layer Manager (Linux or Windows), run `workflow_GRASS.py`. For example, if GSFLOW-GRASS is in your "models" folder:

**TO DO: MUST INTEGRATE THIS WITH CRYSTAL'S PYTHON SETTINGS.INI READER, SO DON'T HAVE TO BE IN THE SAME LOCATION (ADDPATH)**

```sh
python ~/models/GSFLOW-GRASS/domain_builder/workflow_GRASS.sh
```

**TO DO: how to manage flow to the ocean with null cells? see v.gsflow.grid**

Time will pass and a lot of text will go past on the screen. If it ends with "Done.", regardless of warning/error messages about adding fields to shapefiles. If it does not end with "Done.", please contact us!

Once this has finished check our **"gsflow_simdir"** for a **"GIS"** subfolder that contains the outputs of your work here. The files there will be automatically read in during Step 4.

Pat yourself on the back! The GRASS portion is complete.

### Step 3: Customize the Go-GSFLOW File to set input-file-builder options

The Go-GSFLOW File (`go-GSFLOW.sh` on Linux and `go-GSFLOW.bat` on Windows) in the 'Run' folder, is for pre-processing and running GSFLOW.

**\todo{Crystal: change file name; currently `run_Python_GSFLOW_current.sh.`}**

At the top of the file, the user should customize:

* `toolkit_dir` full pathname for location of GSFLOWGRASS_toolkit
* `settings_file` Full pathname for Settings File (customized in Step 1).

The rest of the file will execute scripts to: (1) set up certain model inputs (climate forcing data and hydraulic conductivity) according to the *Settings* File, (2) create GSFLOW input files (GSFLOW control file, PRMS parameter file, and MODFLOW input files), and (3) run GSFLOW.  In the default implementation, the user does not need to change the bottom part of the file with Python scripts.  However, certain lines may be commented out or changed, as described in the following:

* `Create_hydcond_array.py`: If *fl_create_hydcond*=1 in Settings File, this script creates spatially distributed hydraulic conductivity values; see top of this script to select from options.  This line may be changed to a different script name if the user writes their own script for creating spatially distributed hydraulic conductivity. If *fl_create_hydcond*=0, this line may be left in; nothing will be done in the script.
* `GSFLOW_print_data_climatehru_files1_metric.py` If *fl_print_climate_hru*=1 in Settings File, this script creates climate_hru files with spatially uniform conditions, based on data from *climate_data_file* in the Settings File.  This line may be changed to a different script name if the user writes their own script for creating spatially distributed climate inputs. If *fl_print_climate*=0, this line may be left in; nothing will be done in the script.
* `GSFLOW_print_controlfile_current.py`: This script creates GSFLOW control file.  This line generally should be included, but it may be commented out if the user has already run the script previously and will be using the same file in its same location.
* `GSFLOW_print_PRMSparamfile_current.py`: This script creates the PRMS parameter file.  This line generally should be included, but it may be commented out if the user has already run the script previously and will be using the same file in its same location.
* `print_MODFLOW_inputs_res_NWT_current.py`: This script creates all the MODFLOW input files.  This line generally should be included, but it may be commented out if the user has already run the script previously and will be using the same files in their same location.
* `run_GSFLOW.py`: This script executes the GSFLOW model.  This line generally should be included, but it may be commented out if the user only wishes to create the input files without running the model.

Note that the above Python scripts can also be run independently using Python, outside of the Go-GSFLOW File; just be sure to include the *Settings File* name as an argument.

### Step 4. Optional steps

For the default implementation, the user can proceed to Step 5.  However, extra steps are needed if the user has specified any of the following:

* **Settings File, fl_create_hydcond=1**: Set options at the top of the `Create_hydcond_array.py` input-file-builder script for different spatial distribution configurations.  Other steps may be needed if the user replaces this script with their own to create spatially distributed hydraulic conductivity fields.
* **Settings File, fl_print_climate=1**: Create the file specified in climate_data_file, which should have climate data time series from one weather station for daily minimum temperature, maximum temperature, and precipitation.  See example problems for the file format (e.g., in `Shullcas -> UserData ->` ), which has the following format:
  * line 1: comment
  * line 2: `tmax` for daily maximum temperature, `1` for number of weather stations
  * line 3: `tmin` for daily minimum temperature, `1` for number of weather stations
  * line 4: `precip` for daily precipitation, `1` for number of weather stations
  * line 5: `####################################` to indicate start of data
  * line 6: `YYYY Month Day 0 0 0 (tmax value) (tmin value) (precip value)`, etc. for all dates in daily time series
All temperature data in this file are assumed to be in [&deg;C], and precipitation data in [mm/d] (these are eventually converted to [&deg;C] and [in/d] for the PRMS model component).  This file can be expanded to include relative humidity (in [%]) (used if Penman-Monteith option is selected for the potential ET module) and solar radiation (in [MJ/m<sup>2</sup>]) if available.  Other steps may be needed if the user replaces this script with their own to create spatially distributed hydraulic conductivity fields.
* **Settings File, fl_print_climate_hru=0**: Make sure climate_hru_dir is directory with pre-existing climate_hru data files containing HRU-distributed climate inputs: **tmin.day**, **tmax.day**, **precip.day**, and **empty.day**. See GSFLOW manual or example files in example cases (e.g., in Shullcas -> inputs -> PRMS_GSFLOW) for format of climate_hru data files.

*Changing other input parameters*: See "Advanced Customization" below.

### Step 5. Running GSFLOW
The pre-processing and GSFLOW model execution can be carried out by entering the Go-GSFLOW at the command line:
* Linux prompt: ./go-GSFLOW.sh (user needs to press Enter at the end, to return to command prompt)
* Windows command prompt: .\go-GSFLOW.bat (or by double clicking on the .bat file in windows explorer)

Output files will be located in the gsflow_simdir (specified in Step 1 in the *Settings File*) -> *outputs*.

After running GSFLOW, the user should check the following:
* Model ran to completion: GSFLOW will print "Normal termination of simulation".  For Windows, this will print to the screen.  For Linux, this will appear at the end of *out.txt* created in the directory where go-GSFLOW.sh is executed (all output to screen is redirected to out.txt).
* Model converged: Check gsflow_simdir -> *control* -> *gsflow.log*.  This will summarize convergence performance.
* Model ran with good water balance: Check gsflow_simdir -> *outputs* -> *PRMS_GSFLOW* -> *gsflow.out*.  This will report instances when there are large water balance discrepancies (indicated with "WARNING...")
* MODFLOW component ran smoothly: Check gsflow_simdir -> *outputs* -> *MODFLOW_NWT* -> *test.lst*.  This will include convergence and run details specific to MODFLOW, which are only summarized in *gsflow.log*.


### Step 6. Visualization
Our toolkit includes Python scripts in `GSFLOW-toolkit -> visualization` for graphically depicting major GSFLOW inputs and outputs.  Each of these scripts can be run in Python using the following syntax at a Python console: *TODO: Andy - checking my wording, do you "run at a Python console"?*

```bash
run (visualization_script).py (Settings File full path)
```

By default, the visualization scripts will plot output files in directories specified by *Settings File*; this facilitates visualizing results directly after running GSFLOW with this toolkit.  However, the user can over-ride these file locations in the section `*** CHANGE FILE NAMES AS NEEDED` in order to plot arbitrary output files (not based on *Settings File*).  However, our visualization scripts assume certain output file formats and are not guaranteed to work with GSFLOW output files generated outside of this toolkit.

For other main plotting options, the user should use the section at the top of each visualization script indicated with `*** SET THE FOLLOWING:...`.  Other more detailed customizations (e.g., figure formatting etc) can be manually made in the rest of the script.

The visualization scripts in the toolkit include the following:

* Visualization scripts for plotting GSFLOW inputs (inputs that are created by our toolkit):
  * `plotBasin.py`: Creates figure showing HRU and stream discretization.  User specifies HRU id's to be highlighted.
  * `Plot_MODFLOW_inputs.py`: Creates various figures showing spatially discretized MODFLOW inputs: active grid cells, top elevation of each layer, bottom elevation of each layer, and hydraulic conductivity for each layer.
* Visualization scripts for plotting GSFLOW outputs:
  * `plotHRUvars.py`: Creates movie of HRU-discretized output variables.
  * `plotSegmentDischarge.py`: Creates movie of stream segment-discretized discharge values.
  * `Plot_MODFLOW_3D_head3.py`: Creates movie of spatially discretized hydraulic head.  Also options to show water table depth and change-in-head over print time steps.
  * `Plot_MODFLOW_3D_uzf.py`: Creates movie of spatially discretized unsaturated-zone outputs, such as recharge from the unsaturated to the saturated zone.
  * `plot_gsflow_csv.py`: Creates a time series figure with basin-total variables.

## Advanced Customization

**Changing other parameters in GSFLOW input files:**

Our toolkit is set up to easily change hydraulic conductivity, climate, and infiltration inputs through the *Settings* File.  To change other model input parameters (described in the GSFLOW manual), the user must locate those entries in the Python input-file-builder scripts and edit the values there.  These scripts are in Toolkit_GSFLOW ->input_file_builder and include: 
* `GSFLOW_print_controlfile_current.py`: Builds GSFLOW control file, which controls model options.  See commented Section headings to make changes.
* `GSFLOW_print_PRMSparamfile_current.py`: Builds PRMS parameter file, which contains all (non-stream) surface properties in "Section 2: Parameters."  While any of these may be customized, those of particular interest are commented with "# *** CHANGE FOR SPECIFIC SITE"
* `MODFLOW_NWT_lib_current.py`: Library of functions to build the various MODFLOW input files (used in `print_MODFLOW_inputs_res_NWT_current.py.`  See individual functions to change input parameters for the different corresponding MODFLOW packages.
* Subdirectories of `domain_builder` hold the GRASS GIS commands. If the user needs to expand the set of variables exported under these commands, one should:
  1. Open the desired `.py` file and follow the pattern shown to edit its parser inputs/outputs (formatted comments near the top) and the Python code inside that first parses these comments and then runs GRASS GIS commands
  2. Ensure that you have a C compiler installed (gcc)
  3. Run `make MODULE_TOPDIR=<PATH_TO_YOUR_GRASS_GIS_INSTALL>` to compile the updated GRASS GIS commands. See and/or `makeall.sh`, in which this path is provided by the user, for a simple way to recompile many commands.

**If you make modifications, please contact us.** We would be happy to include your contributions in the main code base, and any changes to the GRASS GIS commands to the main GRASS GIS extensions (add-ons) repository: this will help future users!
