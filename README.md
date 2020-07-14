# gwcatpy
Data from LIGO-Virgo gravitational wave detections, built from [Gravitational Wave Open Science Centre](https://www.gw-openscience.org/) and [GraceDB](https://gracedb.ligo.org/latest/). For usage see gwcat.cardiffgravity.org

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
