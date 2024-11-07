"""
This function download the cams data used in grs from CDS API

The main program takes two argument : start  and end year
"""

import argparse
import os
from datetime import date, timedelta

import cdsapi
from pathlib import Path


def main(dic, i):
    data_type = 'cams-global-reanalysis-eac4'  # dic['mode']
    today = date.today() - timedelta(days=i)
    print(today)
    yesterday=today - timedelta(days=1)
    print(yesterday)
    d1 = today.strftime("%Y-%m-%d")
    
    # specify the period to catch data
    # try:
    odir = Path(dic['outdir'], yesterday.strftime("%Y"), yesterday.strftime("%m"), yesterday.strftime("%d"))
   
    print(str(odir))
    if not os.path.exists(odir):
        os.makedirs(odir)
    date_str = str(yesterday)+"/"+str(yesterday)
    print(date_str)

    c = cdsapi.Client()
    if dic['mode'] == 'reanalysis':
            data_type = 'cams-global-reanalysis-eac4'
            datafile =  odir + str(yesterday.strftime("%Y-%m-%d") + "-" + data_type + '.nc')
            print('processing ' + datafile + '...')
            c.retrieve(
                data_type,
                    {
                    'nocache': '456',
                    'format': 'netcdf',
                    'date': date_str,  # '2003-09-01/2003-09-30',
                    'time': [
                        '00:00', '03:00', '06:00',
                        '09:00', '12:00', '15:00',
                        '18:00', '21:00',
                    ],
                    'variable': [
                        '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_temperature',
                        'black_carbon_aerosol_optical_depth_550nm', 'dust_aerosol_optical_depth_550nm',
                        'mean_sea_level_pressure',
                        'organic_matter_aerosol_optical_depth_550nm', 'sea_salt_aerosol_optical_depth_550nm',
                        'sulphate_aerosol_optical_depth_550nm',
                        'surface_pressure', 'total_aerosol_optical_depth_1240nm',
                        'total_aerosol_optical_depth_469nm',
                        'total_aerosol_optical_depth_550nm', 'total_aerosol_optical_depth_670nm',
                        'total_aerosol_optical_depth_865nm',
                        'total_column_carbon_monoxide', 'total_column_methane', 'total_column_nitrogen_dioxide',
                        'total_column_ozone', 'total_column_water_vapour',
                    ],
               },
                datafile)
    else:
            data_type = 'cams-global-atmospheric-composition-forecasts'
            datafile = Path(odir, f"{yesterday.strftime('%Y-%m-%d')}-{data_type}.nc").resolve()
            if datafile.exists():
                print(f'!! {datafile} already exists !!')
            print(f'processing {datafile} ...')
            c.retrieve(
                data_type,
                {
                    'nocache': '456',
                    'date': date_str,
                    'type': 'forecast',
                    'format': 'netcdf',
                    'variable': [
                        '10m_u_component_of_wind', '10m_v_component_of_wind',
                        '2m_temperature',
                        'mean_sea_level_pressure', 'surface_pressure',
                        'ammonium_aerosol_optical_depth_550nm', 'black_carbon_aerosol_optical_depth_550nm',
                        'dust_aerosol_optical_depth_550nm',
                        'nitrate_aerosol_optical_depth_550nm', 'organic_matter_aerosol_optical_depth_550nm',
                        'sea_salt_aerosol_optical_depth_550nm',
                        'secondary_organic_aerosol_optical_depth_550nm', 'sulphate_aerosol_optical_depth_550nm',
                        'total_aerosol_optical_depth_1240nm',
                        'total_aerosol_optical_depth_469nm',
                        'total_aerosol_optical_depth_550nm',
                        'total_aerosol_optical_depth_670nm',
                        'total_aerosol_optical_depth_865nm',
                        'total_column_carbon_monoxide',
                        'total_column_methane',
                        'total_column_nitrogen_dioxide',
                        'total_column_ozone', 'total_column_water_vapour'],

                    'time': ['00:00', '12:00'],
                    'leadtime_hour': ['0', '3', '6', '9'],

                },
                str(datafile))


# except:
#    print('Error: not appropriate cams settings for download. Refers to ecmwf.')
#    sys.exit()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Download Cams datasets from CDS')
    parser.add_argument('mode',
                        help='choose `reanalysis` or `forecast` dataset')
    parser.add_argument('outdir',
                        help='output directory where the products will be downloaded',
                        default="/datalake/watcal/ECMWF/CAMS/")
    args = parser.parse_args()
    #for i in range(1, 100):
    main(vars(args), 2)

