# gwcatpy
Python library to access and update data from LIGO-Virgo gravitational wave detections, built from [Gravitational Wave Open Science Centre](https://www.gw-openscience.org/) and [GraceDB](https://gracedb.ligo.org/latest/). For usage see gwcat.cardiffgravity.org

## Data output
For output see https://data.cardiffgravity.org/gwcat-data/
For javascript interaction with the data see https://gwcat.cardiffgravity.org/

##Usage (python)

### Importing the libarary
Make sure the directory containing the gwcatpy module is in your path. For example, if that is the current directory then the following should work
~~~~
import sys
sys.path.insert(0,'../') 
import gwcatpy
~~~~

### Loading the data
Download the data from https://data.cardiffgravity.org/gwcat-data/data/gwosc_gracedb.json
~~~~
GC = gwcatpy.GWCat(filein=<path-to-json-data-file>)
~~~~
This will make a data directory (if it doesn't already exist), and a `status.json` file (if that doesn't exist).

### Accessing the data

 * `GC = gwcatpy.GWCat(filein=<path-to-json-data-file>)`
 * `GC.data`: dictionary object containing event data, indexed by event name
 * `GC.links`: dictionary object containing links, indexed by event name
 * `GC.datadict`: dictionary object containing parameters, indexed by parameter name

 * `GC.events`: pandas Dataframe object containing event data
 * `GC.eventrefs`: pandas Dataframe object containing event links
 * `GC.events`: pandas Dataframe object containing parameter data

##### To export to json:
~~~~
GC.exportJson(fileout,[dir='', verbose=False])
~~~~
 * fileout = [string] filename of data file (containing data, parameters and links)
 * dir = [string, optional] directory to export data to (default=current directory)
 * verbose = [boolean, optional] set to print messages

##### To export to csv:
~~~~
GC.exportCSV(datafileout,[dictfileout=None, linksfileout=None, dir='', verbose=False, clearcols=True])
~~~~
 * datafileout = [string] filename of data file
 * dir = [string, optional] directory to export data to (default=current directory)
 * dictfileout = [string, optional] filename to export parameters to(if not present, no parameters output)
 * linksfileout = [string, optional] filename of links to (if not present, no links output)
 * verbose = [boolean, optional] set to print messages
 * clearcols = [boolean, optional] set to clear empty rows if present

##### To export to Excel:
~~~~
GC.exportExcel(fileout, [dir='', verbose=False])
~~~~
 * fileout = [string] filename of data file (containing data, parameters and links in separate sheets)
 * dir = [string, optional] directory to export data to (default=current directory)
 * verbose = [boolean, optional] set to print messages

### Produced data files (GWOSC + GraceDB):
 * [all data (JSON)](https://git.ligo.org/chris.north/gwcat-data-dev/-/blob/master/data/gwosc_gracedb.json) (_includes parameter definitions and referece links_)
 * [all data (JSON, minified))](https://git.ligo.org/chris.north/gwcat-data-dev/-/blob/master/data/gwosc_gracedb.json.min) (_includes parameter definitions and reference links_)
 * [data values (CSV)](https://git.ligo.org/chris.north/gwcat-data-dev/-/blob/master/data/gwosc_gracedb.csv)
 * [Data parameters (CSV)](https://git.ligo.org/chris.north/gwcat-data-dev/-/blob/master/data/parameters.csv)
 * [References/links (CSV)](https://git.ligo.org/chris.north/gwcat-data-dev/-/blob/master/data/gwosc_gracedb_links.csv)

### Sky localisation images:
 * Location: http://data.cardiffgravity.org/gwcat-data/data/png/\<filename\>
 * Filename formats for \<event\> (e.g. GW150914, S190510g etc.):
    * Mollweide: data/png/\<event\>_moll.png
        * (e.g. [data/png/S191105e_moll.png](http://data.cardiffgravity.org/gwcat-data/data/png/S191105e_moll.png))
    * Mollweide (rotated to centre on peak location): data/png/\<event\>_moll_rot.png
        * (e.g. [data/png/S191105e_moll_rot.png](http://data.cardiffgravity.org/gwcat-data/data/png/S191105e_moll_rot.png))
    * Cartesian fullsky: data/png/\<event\>_cart.png
        * (e.g. [data/png/S191105e_cart.png](http://data.cardiffgravity.org/gwcat-data/data/png/S191105e_cart.png)
    * Cartesian fullsky (rotated to centre on peak location): data/png/\<event\>_cart_rot.png
        * (e.g. [data/png/S191105e_cart_rot.png](http://data.cardiffgravity.org/gwcat-data/data/png/S191105e_cart_rot.png)
    * Cartesian zoomed onto peak location: data/png/\<event\>_cartzoom.png
        * (e.g. [data/png/S191105e_cartzoom.png](http://data.cardiffgravity.org/gwcat-data/data/png/S191105e_cartzoom.png)
 
### High-res sky localisation images:
 * Location: http://data.cardiffgravity.org/gwcat-data/data/gravoscope/\<event\>_8192.png
   * e.g. [data/gravoscope/S191105e_9192.png](http://data.cardiffgravity.org/gwcat-data/data/gravoscope/S191105e_9192.png)
 * 8192x4096 px
 * grayscale image

### Dependencies
* python>=3.6 (it almost certainly works in earlier versions of python3, and it may well work in python 2, possibly with minor changes)
* Contains the submodule gwcatpy, which has the following dependencies:
 * LIGO-specific dependencies
   * pesummary
   * ligo.gracedb
 * Other community-developed dependencies:
   * astropy
   * healpy
   * matplotlib
   * pandas
   * h5py
   * astropy_healpix
   * numpy
   * scipy
   * requests
   * json
   * bs4
   
### Processing

##### Data access
* gwcatpy.gwosc.getGWTC
    * Download JSON catalogue from from GWOSC
    * Access individual JSON files for each event
* gracedb.getSuperevents
    * Download superevents list
        * Access individual event records
        * Check version/date and redownload if necessary

##### Catalogue processing methods
See testcat.py for example script. Status is saved in data/status.json (note that data/ is not included in this repo)

* importGWTC:
    * Convert parameters to gwcat format
        * mostly parameter names from parameters > gwtc2_pe_<event>
        * FAR and SNR from parameters > gwtc2_gstlal_allsky_<event>
    * Add links for open-data, sky-map etc.
* importGraceDB:
    * Convert parameters to gwcat format
    * Add links for open-data, sky-map etc.
    * Download skymap fits file
* matchGraceDB:
    * identify GraceDB entries that already have confirmed detection entry
* removeCandidates:
    * remove GraceDB entries that have confirmed detection
* addRefs:
    * add references to sky-maps manually for GWTC1 (not included in GWOSC machine-readable data)
* updateH5:
    * dowload HDF (.hdf/.h5) files for GWOSC data
        * may need extracting from tarball, along with skymap
    * get parameters from PublicationSamples approximant and update database values
* setPrecision:
    * reduce precision of stored data
* updateMaps:
    * download maps that aren't already downloaded
    * Calculate 90% sky area
* plotMapPngs:
    * make range of PNG maps for skymaps
* makeGravoscopeTiles:
    * make tilesets for Gravoscope
* makeWaveforms:
    * make illustrative waveforms
* exportJson:
    * export to JSON (data/gwosc_gracedb.json)
* exportCSV:
    * export to CSV files
        * Data values in data/gwosc_gracedb.csv
        * Parameters in data/parameters.csv
        * Links in data/gwosc_gracedb_links.csv

The data can also be exported to jsonp format, which seems to more nicely with CORS.

### Data format

The data array contains an object for each event. Each event has a set of parameters (M1, M2, Mchirp, UTC, etc.), within which is an object containing the relevant values.

##### Data values

gwcat.data = [
	{
		"Event1 Parameter1": {
			"Event1 Parameter1 Value1": string/number/array
		},
		"Event1 Parameter2": {
			"Event1 Parameter2 Value1": string/array/value
			"Event1 Parameter2 Value2": string/array/value
		},...
	},
	{
		"Event2 Parameter1":{
			"Event2 Parameter1 Value1": string/number/array
		},...
	},...
]


The parameters are those recorded for each event. Note that not all are present for each event.

 * name: Event name (e.g. GW150914)
 * UTC: UTC time of detection (YYYY-MM-DDThh:mm:ss)
 * GPS: GPS time of detection
 * M1: Primary mass (Msun)
 * M2: Secondary mass (Msun)
 * Mchirp: Chirp mass (Msun)
 * Mtotal: Total mass (Msun)
 * Mfinal: Final mass (Msun)
 * Mratio: Mass ratio
 * chi: Effective inspiral spin
 * af: Final spin
 * DL: Luminosity distance (MPc)
 * z: Redshift
 * lpeak: Peak luminosity (1056 erg s-1)
 * Erad: Radiated energy (Msun c2)
 * FAR: False alarm rate (yr-1)
 * deltaOmega: Sky localization area (deg2)
 * rho: Signal-to-noise ratio

The values can be any of:

 * best: exact of best-fit value of parameter. Can be string (e.g. name, UTC), or number (e.g. masses, spins, GPS).
 * lower: a (numerical) lower limit on the parameter.
 * upper: a (numerical) upper limit on the parameter.
 * lim: a two-element array (of numberse) containing the range of plausible values (where applicable), in order [min, max]. Used where a best-fit value and corresponding error isn't appropriate.
 * err: a two-element array (of numberse) containing the errors on the "best" value (where applicable), in order [upper, lower]. Only accompanies a "best" value.

#### Links

For each event, links are stored as a list, with each entry being an object with various parameters

gwcat.links = {
     "Event1":[
         {
             "Parameter1": value
             "Parameter2": value
         },
         {
             "Parameter1": value
             "Parameter2": value
         },...
     ],
     "Event2":[
         {
             "Parameter1": value
             "Parameter2": value
         },
         {
             "Parameter1": value
             "Parameter2": value
         },...
     ],...
}

Parameters are:
 * url: url to link
 * text: description of link
 * type: type of link (designed to be machine-readable)
 * filetype: type of file (for distinduishing between h5 and tar files)