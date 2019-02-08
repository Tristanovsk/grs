# import ez_setup
# ez_setup.use_setuptools()

from setuptools import setup, find_packages
from grs import __version__, __package__

setup(
    name=__package__,
    version=__version__,
    packages=find_packages(),
    package_data={'':['*.so']},
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
    install_requires=['pandas','geopandas','scipy','numpy','netCDF4',
                      'matplotlib','docopt','GDAL','python-dateutil','pyproj'],

    entry_points={
          'console_scripts': [
              'grs = grs.run:main'
          ]}
)
