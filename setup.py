from setuptools import setup, find_packages

# read the contents of README file
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

__package__ = 'grs'
__version__ = '2.1.6'

setup(
    name=__package__,
    version=__version__,
    # other arguments omitted
    long_description=long_description,
    long_description_content_type='text/markdown',

    packages=find_packages(),
    package_data={
        # If any package contains *.txt files, include them:
        'grs': ['*.yml'],
        '': ['data/lut/gases/*.nc', 'data/aux/*.txt'],

    },
    include_package_data=True,

    url='',
    license='MIT',
    author='T. Harmel',
    author_email='tristan.harmel@gmail.com',
    description='scientific code for atmosphere and sunglint correction of Sentinel-2 and Landsat series',

    # Dependent packages (distributions)
    # install_requires=['cdsapi', 'pandas', 'geopandas', 'scipy', 'numba', 'numpy', 'netCDF4',
    #                  'matplotlib', 'docopt', 'python-dateutil', 'pyproj', 'xarray'],

    entry_points={
        'console_scripts': [
            'grs = grs.run:main'
        ]}
)
