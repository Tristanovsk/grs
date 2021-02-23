# import ez_setup
# ez_setup.use_setuptools()

from setuptools import setup, find_packages
#from grs import __version__, __package__

__package__ = 'grs'
__version__ = '1.2.4'

setup(
    name=__package__,
    version=__version__,
    packages=find_packages(exclude=['build']),
    package_data={'':['*.so','*h','*angles*']},
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
    #install_requires=['pandas','geopandas','scipy','numpy','netCDF4',
    #                  'matplotlib','docopt','GDAL','python-dateutil','pyproj'],

   
    install_requires=['pandas==1.1.3','geopandas','scipy==1.5.3','numpy','netCDF4',
                      'matplotlib==3.3.3','docopt','GDAL','python-dateutil','pyproj'],


    entry_points={
          'console_scripts': [
              'grs = grs.run:main'
          ]}
)
