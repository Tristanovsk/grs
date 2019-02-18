# GRS algorithm package
## GRS (Glint Removal for Sentinel-2-like sensors)

...

see [Harmel et al., 2018](https://www.sciencedirect.com/science/article/pii/S0034425717304856)

![flowchart](images/flowchart_sunglint_S2.png)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

Download and install the [SNAP software](http://step.esa.int/main/download/). 
Configure the [SNAP-python interface](https://senbox.atlassian.net/wiki/spaces/SNAP/pages/50855941/Configure+Python+to+use+the+SNAP-Python+snappy+interface) 
and link the obtained `snappy` folder to your python site-packages as `esasnappy`. For example:

```bash
ln -s /FULL_PATH/.snap/snap-python/snappy /PATH_TO_LIB_PYTHON/lib/python3.6/site-packages/esasnappy
```


Compilers such gcc and gfortran are needed to install the package.

Compile all C and fortran files into shared libraries:

```
make
```

### Installing

To install the package:
```
python setup.py install
```

or 

```
python setup.py install --user
```

If the installation is successful, you should have:
```
$ grs
Usage:
   grs <input_file> [--sensor <sensor>] [-o <ofile>] [--odir <odir>] [--shape <shp>] [--wkt <wktfile>]   [--longlat <longmax,longmin,latmax,latmin> ]    [--altitude=alt] [--aerosol=DB] [--aeronet=<afile>]    [--aot550=aot] [--angstrom=ang] [--output param]   [--resolution=res] [--levname <lev>] [--no_clobber]
  grs -h | --help
  grs -v | --version
  
```

## Running the tests

Examples of output images

![image_output](images/Fig_valid_qualit_sea_scale.png)

## Deployment

Add additional notes about how to deploy this on a live system

## Contributing

Please contact [authors](tristan.harmel@ntymail.com) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Tristan Harmel** - *Initial work* - [contact](tristan.harmel@ntymail.com)

See also the list of [contributors](...) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* The [Step forum](http://forum.step.esa.int) and Marco Peters are acknowledged for their useful help to process Sentinel-2 data
with the snappy API.
* The authors are very grateful to Olivier Hagolle
for providing open source codes to perform gaseous absorption correction and massive Sentinel-2 data download.
