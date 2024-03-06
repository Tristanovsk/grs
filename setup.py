from setuptools import setup, find_packages

__package__ = 'grs'
__version__ = '2.1.4'

setup(
    name=__package__,
    version=__version__,
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
    #install_requires=['cdsapi', 'pandas', 'geopandas', 'scipy', 'numba', 'numpy', 'netCDF4',
    #                  'matplotlib', 'docopt', 'python-dateutil', 'pyproj', 'xarray'],

    entry_points={
        'console_scripts': [
            'grs = grs.run:main'
        ]}
)
