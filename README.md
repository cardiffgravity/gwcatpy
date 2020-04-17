# gwcat-data
Data from LIGO-Virgo gravitational wave detections, built from [Gravitational Wave Open Science Centre](https://www.gw-openscience.org/) and [GraceDB](https://gracedb.ligo.org/latest/). For usage see gwcat.cardiffgravity.org

### Data files (GWOSC + GraceDB):
 * [all data (JSON)](data/gwosc_gracedb.json) (_includes parameter definitions and referece links_)
 * [all data (JSON, minified))](data/gwosc_gracedb.json.min) (_includes parameter definitions and reference links_)
 * [data values (CSV)](data/gwosc_gracedb.csv)
 * [Data parameters (CSV)](data/parameters.csv)
 * [References/links (CSV)](data/gwosc_gracedb_links.csv)

### Sky localisation images:
 * Location: http://data.cardiffgravity.org/gwcat-data/data/png/\<filename\>
 * Filename formats for \<event\> (e.g. GW150914, S190510g etc.):
    * Mollweide: data/png/\<event\>_moll.png
        * (e.g. [data/png/S191105e_moll.png](data/png/S191105e_moll.png))
    * Mollweide (rotated to centre on peak location): data/png/\<event\>_moll_rot.png
        * (e.g. [data/png/S191105e_moll_rot.png](data/png/S191105e_moll_rot.png))
    * Cartesian fullsky: data/png/\<event\>_cart.png
        * (e.g. [data/png/S191105e_cart.png](data/png/S191105e_cart.png)
    * Cartesian fullsky (rotated to centre on peak location): data/png/\<event\>_cart_rot.png
        * (e.g. [data/png/S191105e_cart_rot.png](data/png/S191105e_cart_rot.png)
    * Cartesian zoomed onto peak location: data/png/\<event\>_cartzoom.png
        * (e.g. [data/png/S191105e_cartzoom.png](data/png/S191105e_cartzoom.png)
 
 
### High-res sky localisation images:
 * Location: http://data.cardiffgravity.org/gwcat-data/data/gravoscope/\<event\>_8192.png
   * e.g. [data/gravoscope/S191105e_9192.png](data/gravoscope/S191105e_9192.png)
 * 8192x4096 px
 * grayscale image
