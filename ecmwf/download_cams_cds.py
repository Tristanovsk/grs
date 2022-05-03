"""
This function download the cams data used in grs from CDS API

The main program takes two argument : start  and end year
"""

import os, sys

import argparse
from datetime import date
import calendar
import cdsapi


def main(dic):
    data_type = 'cams-global-reanalysis-eac4'  # dic['mode']
    # month = dic['month']
    year_start = int(dic['year_start'])
    year_end = int(dic['year_end'])
    # area = "90/-180/-90/180"

    # specify the period to catch data
    # try:
    for year in range(year_start, year_end + 1):
        odir = "/datalake/watcal/ECMWF/CAMS/" + str(year) + "/"
        if not os.path.exists(odir):
            os.makedirs(odir)
        for month in range(1, 13):
            numberOfDays = calendar.monthrange(year, month)[1]
            date = '-'.join((str(year), str(month).zfill(2), "01/" + str(year),
                             str(month).zfill(2), str(numberOfDays)))
            print(date)

            c = cdsapi.Client()
            if dic['mode'] == 'reanalysis':
                data_type = 'cams-global-reanalysis-eac4'
                datafile = odir + str(year) + '-' + str(month).zfill(
                    2) + '_month_' + data_type + '.nc'
                if os.path.exists(datafile):
                    continue
                print('processing ' + date + datafile + '...')

                c.retrieve(
                    data_type,
                    {
                        'format': 'netcdf',
                        'date': date,  # '2003-09-01/2003-09-30',
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
                datafile = odir + str(year) + '-' + str(month).zfill(
                    2) + '_month_' + data_type + '.nc'

                if os.path.exists(datafile):
                    print('!!' + datafile + 'already exists !!')
                    continue
                print('processing ' +date+ datafile + '...')

                c.retrieve(
                    data_type,
                    {   'nocache': '123',
                        'date': date,
                        'type': 'forecast',
                        'format': 'netcdf',
                        'variable': [
                            '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_temperature',
                            'mean_sea_level_pressure', 'surface_pressure',
                            'single_scattering_albedo_1020nm',
                            'single_scattering_albedo_1240nm',
                            'single_scattering_albedo_1640nm',
                            'single_scattering_albedo_2130nm',
                            'single_scattering_albedo_355nm',
                            'single_scattering_albedo_380nm',
                            'single_scattering_albedo_400nm',
                            'single_scattering_albedo_440nm',
                            'single_scattering_albedo_500nm',
                            'single_scattering_albedo_550nm',
                            'single_scattering_albedo_645nm',
                            'single_scattering_albedo_670nm',
                            'single_scattering_albedo_800nm',
                            'single_scattering_albedo_865nm',
                            'total_aerosol_optical_depth_1020nm',
                            'total_aerosol_optical_depth_1064nm',
                            'total_aerosol_optical_depth_1240nm',
                            'total_aerosol_optical_depth_1640nm',
                            'total_aerosol_optical_depth_2130nm',
                            'total_aerosol_optical_depth_355nm',
                            'total_aerosol_optical_depth_380nm',
                            'total_aerosol_optical_depth_400nm',
                            'total_aerosol_optical_depth_440nm',
                            'total_aerosol_optical_depth_469nm',
                            'total_aerosol_optical_depth_500nm',
                            'total_aerosol_optical_depth_550nm',
                            'total_aerosol_optical_depth_645nm',
                            'total_aerosol_optical_depth_670nm',
                            'total_aerosol_optical_depth_800nm',
                            'total_aerosol_optical_depth_865nm',
                            'total_column_carbon_monoxide', 'total_column_formaldehyde',
                            'total_column_hydroxyl_radical', 'total_column_methane', 'total_column_nitrogen_dioxide',
                            'total_column_ozone', 'total_column_propane', 'total_column_water_vapour',
                        ],
                        'leadtime_hour': [
                            '0', '12', '18',
                            '21', '3', '6',
                            '9',
                        ],
                        'time': '00:00',
                    },
                    datafile)


# except:
#    print('Error: not appropriate cams settings for download. Refers to ecmwf.')
#    sys.exit()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Download Cams datasets from CDS')
    parser.add_argument('year_start',
                        help='the start year of period to be downloaded ')
    parser.add_argument('year_end',
                        help='the end year of period to be downloaded')
    parser.add_argument('mode',
                        help='choose `reanalysis` or `forecast` dataset')
    args = parser.parse_args()

    main(vars(args))
