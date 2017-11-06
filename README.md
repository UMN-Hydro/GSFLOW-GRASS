# GSFLOW-GRASS

***Python toolkit using GRASS GIS to generate inputs for the USGS hydrologic model "GSFLOW".***

GSFLOW is the US Geological Survey's combined groundater (MODFLOW) and surface water (PRMS) model. GSFLOW-GRASS builds a stream network and its "hydrologic response units" -- here, watershed sub-basins -- that combine to form a surface-water hydrologic network (PRMS). It then links these to a regular finite difference grid to compute groundwater flow using MODFLOW.

These instructions are meant to take an user familiar with computers but new to (or a beginner with) GSFLOW and GRASS GIS through the basics of how to get GSFLOW-GRASS to work. *Please leave a message if you have trouble working with GSFLOW-GRASS; your comments could assist both you and the more general improvement of this documentation.*

When you use GSFLOW-GRASS, please contact us; a reference (Ng et al.) will be noted here as soon as the paper is made available.

This manual is written in the style of a quick(-ish) start guide that allows users to get GSFLOW up and running for their target watershed using mostly default options in our toolkit and without modifying GRASS-GIS or Python scripts. Additional customization is possible by more advanced users by editing our GRASS GIS and input file building scripts to extend the set of parameters incorporated **\todo{Andy: make section on adding GRASS parameters; perhaps just a small part at the end with appropriate links}**.

## Required Software

* **GSFLOW 1.2**
* **GRASS GIS 7.3+** and extensions (described below)
* **Python 2.7.X**
* **GSFLOW-GRASS Toolkit** (this software)

## Directory Structure:
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

## Pre-Processing

### Step 1: Customize the Settings File

Use `settings_template.ini` as a template for creating your own Settings File, which can have any name. **Boldface** options are required. This files includes:

**\todo{MAKE EVERYTHING REQUIRED BOLD}**

#### "settings" section

| Option             | Description
| ------------------ | ---------
| **proj_name**      | Label for the project; no spaces
| **gsflow_exe**     | Full pathname for GSFLOW executable
| **gsflow_path_simdir** | Full pathname for location where GSFLOW simulation<br>directory should go.
| **fl_print_climate_hru** | **1** to print spatially uniform climate data over all HRU's<br>using climate data from file specified in *climate_data_file*.<br>**0** if user already has pre-existing HRU-distributed climate files.
| climate_data_file  | Only for *fl_print_climate_hru*=1:<br>Name of file containing climate data for single weather station site,<br>to be uniformly distributed over all HRU's using<br>`GSFLOW_print_data_climatehru_files1_metric.py`.<br>If *fl_print_climate_hru*=0, this entry can be omitted;<br>if it is included anyway, it will be ignored.<br>**\todo{create Settings file parser so that it doesn't use this in unless fl_print_climate_hru=1}**
| climate_hru_dir    | **Only for *fl_print_climate_hru*=0**:<br>Name of directory with pre-existing climate_hru data files <br>containing HRU-distributed climate inputs:<br>**tmin.day**, **tmax.day**, **precip.day**, and **empty.day**.<br>See GSFLOW manual or example files in example cases<br>(e.g., in Shullcas -> inputs -> PRMS_GSFLOW) for format of climate_hru data files.<br>**If *fl_print_climate_hru*=1**, this entry can be omitted;<br>if it is included anyway, it will be ignored.<br>**\todo{create Settings file parser so that it doesn't use this in unless fl_print_climate_hru=1.  Add something to Settings parser file to copy over climate_hru files to input directory.**
| **sw_1spinup_2restart** | **1** for spin-up run starting from steady-state MODFLOW period<br>**2** for restart run starting from states saved in the below files
| restart_PRMSfil | optional: for restart runs (sw_1spinup_2restart=2)<br>full pathname of file that is saved under ``save_vars_to_file''<br>in the GSFLOW control file during a previous run.<br>This entry won't be used (but should still be entered) if *sw_1spinup_2restart*=1<br>for startup runs.
| restart_MODfil | optional: for restart runs (sw_1spinup_2restart=2)<br>Full pathname of file that is saved under ``IWRT'' in the MODFLOW name file<br>during a previous run. This entry won't be used (but should still be entered)<br>if sw_1spinup_2restart= 1 for startup runs

#### "custom_params" section

| Option             | Description
| ------------------ | ---------
| fl_create_hydcond  | **1** to implement Python script to create spatially distributed hydraulic conductivity.<br>Otherwise, use values or pre-existing file entered in *hydcond*<br>**\todo{implement this in Create_hydcond.py}**
| hydcond            | For uniform hydraulic conductivity within each layer in the saturated domain<br>enter value(s) (in [m/d]), using comma-separated list for multiple layer domains,<br>starting with top layer.  For spatially distributed values: Enter file name containing<br>array of values (in [m/d]); if *fl_create_hydcond*=1, contents of file will be created<br>using the Python script called in the Go-GSFLOW File (see below description of<br>Go-GSFLOW File).<br>User may need to adjust hydraulic conductivity values to reach numerically<br>convergent results and to match observations.
| finf | Optional: Only for spin-up runs; this entry is ignored (but should still be entered)<br>for restart runs. For uniform infiltration to the unsaturated zone over the watershed:<br>enter a single value (in [m/d]).  For spatially distributed values, enter file name<br>containing array of values (in [m/d]).  May need to adjust this value to reach<br>numerically convergent results and for reasonable start of transient results.

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
| **threshold_drainage_area_meters2** | Threshold drainage area at which flow is considered to create<br>a channel
| **MODFLOW_grid_resolution_meters** | Target cell side length for MODFLOW grid; side lengths will not be<br>exactly this long, as the nearest value to create an integer number<br>of cells in the domain will be chosen.
| **outlet_point_x** | Pour point approximate x (Easting) position; the nearest stream<br>segment to this point is chosen as the true pour point.
| **outlet_point_y** | Pour point approximate y (Northing) position; the nearest stream<br>segment to this point is chosen as the true pour point.
| **icalc** | Method selector for hydraulic geometry computation<br>**0:** Constant<br>**1:** Manning's equation with the wide channel assumption.<br>**2:** Manning's Equation.<br>**3:** Power-law relationship between width, depth, velocity,<br/>and discharge, per Leopold and Maddock (1953)
| gisdb | Optional: Directory that holds grass GIS locations.<br>Typically `~/grassdata`<br>Not currently used.<br>(Will be used to run this while starting GRASS in the background)
| version | Optional: GRASS GIS version number without any "." characters.<br>We used **73**<br>Option is not currently used.<br>(Will be used to run this while starting GRASS in the background)


\todo{how to run GRASS module? From Pre-Processing + Run File? Or separately?}

### Step 2: Customize the Go-GSFLOW File
\todo{increment list number if need separate step for GRASS module}

The Go-GSFLOW File (`go-GSFLOW.sh` on Linux and `go-GSFLOW.bat` on Windows) is executed for pre-processing and running GSFLOW.

\todo{change file name; currently `run_Python_GSFLOW_current.sh.'  Also, update Windows batch file}

At the top of the file, the user should customize:

* `preproc_dir` full pathname for location of GSFLOWGRASS_toolkit] \todo{scripts might end up in directory w/in preproc_dir; make changes accordingly.}
* `settings_file` Full pathname for Settings File (customized in Step 1).

The rest of the file will execute pre-processor scripts to set up certain inputs (climate forcing data and hydraulic conductivity) according to the Settings File, and to create GSFLOW input files (control file, PRMS parameter file, and MODFLOW input files); it then runs GSFLOW.  In the default implementation, the user does not need to change the bottom part of the file with Python scripts.  However, certain lines may be commented out or changed, as described in the following:

* `python2.7 \${preproc_dir}settings_test.py}` Imports all the values set in Settings File.  This line must always be included.
* `python2.7 \${preproc_dir}Create_hydcond_array.py`: If *fl_create_hydcond*=1 in Settings File, this script creates spatially distributed hydraulic conductivity values; see top of script to select from options.  This line may be changed to a different script name if the user writes their own script for creating spatially distributed hydraulic conductivity. If *fl_create_hydcond*=0, this line may be left in; nothing will be done in the script.
* `python2.7 \${preproc_dir}GSFLOW_print_data_climatehru_files1_metric.py` If *fl_print_climate_hru*=1 in Settings File, this script creates climate_hru files with spatially uniform conditions, based on data from *climate_data_file* in Settings File.  This line may be changed to a different script name if the user writes their own script for creating spatially distributed climate inputs. If *fl_print_climate*=0, this line may be left in; nothing will be done in the script.
* `python2.7 \${preproc_dir}GSFLOW_print_controlfile_current.py`: This script creates GSFLOW control file.  This line generally should be included, but it may be commented out if the user has already run the script previously and will be using the same file in its same location.
* `python2.7 \${preproc_dir}GSFLOW_print_PRMSparamfile_current.py`: This script creates the PRMS parameter file.  This line generally should be included, but it may be commented out if the user has already run the script previously and will be using the same file in its same location.
* `python2.7* \${preproc_dir}print_MODFLOW_inputs_res_NWT_current.py`: This script creates all the MODFLOW input files.  This line generally should be included, but it may be commented out if the user has already run the script previously and will be using the same files in their same location.
* `python2.7* \${preproc_dir}run_GSFLOW.py`: This script executes the GSFLOW model.  This line generally should be included, but it may be commented out if the user only wishes to create the input files without running the model.

### Step 3. Optional steps

For the default implementation, the user can proceed to Step 3.  However, extra steps are needed if the user has specified any of the following:

* **Settings File, fl_create_hydcond=1**: Set options at the top of `\${preproc_dir}Create_hydcond_array.py` script for different spatial distribution configurations.  Other steps may be needed if the user replaces this script with their own to create spatially distributed hydraulic conductivity fields.
* **Settings File, fl_print_climate=1**: Create file specified climate_data_file, which should have climate data time series from one weather station for daily minimum temperature, maximum temperature, and precipitation.  See example problems for the file format (e.g., `in Shullcas -> UserData ->` ), which should be the following:
  * line 1: comment
  * line 2: `tmax` for daily maximum temperature, **1** for number of weather stations
  * line 3: `tmin` for daily minimum temperature, **1** for number of weather stations
  * line 4: `precip` for daily precipitation, **1** for number of weather stations
  * line 5: `####################################` to indicate start of data
  * line 6: `YYYY Month Day 0 0 0 (value for tmax) (value for tmin) (value for precip)`, etc. for all dates in daily time series

All temperature data in this file are assumed to be in [&deg;C], and precipitation data in [mm/d] (these are eventually converted to [&deg;C] and [in/d] for the PRMS model component).  This file can be expanded to include relative humidity (in [%]) (used if Penman-Monteith option is selected for the potential ET module) and solar radiation (in [MJ/m<sup>2</sup>]) if available.  Other steps may be needed if the user replaces this script with their own to create spatially distributed hydraulic conductivity fields.

Settings File, *fl_print_climate_hru*=0: **RESUME HERE**

### Step 4. Running GRASS GIS and generating output

#### Download and install GRASS GIS 7.3+

Two options:
* Cross-platform binaries:
https://grass.osgeo.org/download/software/
* Instructions to build from source:
  * https://grasswiki.osgeo.org/wiki/Compile_and_Install
  * https://grasswiki.osgeo.org/wiki/Compile_and_Install_Ubuntu

If you choose to compile GRASS GIS from source, A. Wickert has used these configuration flags many times on Ubuntu (`configure_ubuntu.sh`):

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

## Running GSFLOW
Before running GSFLOW, the user should:




## Visualization

## Advanced Customization

**\todo{Andy, this is where to put ways to customize GRASS GIS code; mostly links should be OK}**
**\todo{Also link this to Crystal's code}**

## More Details **Delete??**

### Pre-processing

#### GRASS domain-builder

#### Settings File

#### GSFLOW Control File builder

#### PRMS Parameter File builder

#### MODFLW Input Files builder

### Visualization Tools

#### Visualization Tool for Model Inputs

#### Visualization Tool for Model Outputs

###
