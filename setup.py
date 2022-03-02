# import ez_setup
# ez_setup.use_setuptools()

from setuptools import setup, find_packages

# from grs import __version__, __package__

__package__ = 'grs'
__version__ = '1.3.4'

setup(
    name=__package__,
    version=__version__,
    packages=find_packages(exclude=['build']),
    package_data={'': ['*.so', '*h', '*angles*']},
    #     # If any package contains *.txt files, include them:
    #     '': ['*.txt'],
    #     'lut': ['data/lut/*.nc'],
    #     'aux': ['data/aux/*']
    # },
    include_package_data=True,

    url='https://gitlab.irstea.fr/telquel-obs2co/satellite/grs',
    license='MIT',
    author='T. Harmel',
    author_email='tristan.harmel@ntymail.com',
    description='scientific code for atmosphere and sunglint correction of Sentinel-2 and Landsat series',

    # Dependent packages (distributions)
    install_requires=['cdsapi', 'pandas', 'geopandas', 'scipy', 'numpy', 'netCDF4',
                      'matplotlib', 'docopt',  'python-dateutil', 'pyproj',
                      'Shapely', 'xarray', 'dask', 'dask[array]'],

    # install_requires=['cdsapi','pandas==1.1.3','geopandas==0.9.0','scipy==1.5.3','numpy==1.19.1','netCDF4==1.5.5',
    #                  'matplotlib==3.3.3','docopt==0.6.2','GDAL==3.0.2','python-dateutil==2.8.1','pyproj==3.0.0.post1', 'Shapely==1.7.1', 'xarray'],

    entry_points={
        'console_scripts': [
            'grs = grs.run:main'
        ]}
)
