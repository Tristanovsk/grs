import sys
from ecmwfapi import ECMWFDataServer
import argparse
from datetime import date


"""
This function download the cams data used in grs 

The main program takes one argument : the mode cams_reanalysis or cams_forecast needed

"""


def main(dic):

        data_type=dic['mode']
        month=dic['month']
        area="90/-180/-90/180"
 
        if month!='all':
            today = date.today()
            d1 = today.strftime("%Y%m%d")

        server = ECMWFDataServer()

        #dataset will be download by default from 2017 to 2020
        step = '0'
        if data_type == 'cams_forecast':
            class_ = 'mc'
            dataset = 'cams_nrealtime'
            time = '00:00:00'
            step = '0/6/12/18'
            type = 'fc'


        # data will be download by default from 2000 to 2017
        elif data_type == 'cams-reanalysis':
            class_ = 'mc'
            dataset = 'cams_reanalysis'
            date='20150301/TO/20170101'
            time = '00:00:00/06:00:00/12:00:00/18:00:00'
            type = 'an'

        elif data_type == 'interim':
            class_ = 'ei'
            dataset = 'interim'
            type = 'an'
        else:
            print('Error: not appropriate dataset for ecmwf/cams download')
            sys.exit()

        import calendar


        #specify the period to catch data
        #try:
        for year in [2017, 2018]:
            for month in range(1, 12):
                numberOfDays = calendar.monthrange(year, month)[1]
                date=str(year)+str(month).zfill(2)+"01/TO/"+str(year)+str(month).zfill(2)+str(numberOfDays)
                print(date)
                server.retrieve({ \
                    'class': class_, \
                    'dataset': dataset, \
                    'date': date, \
                    'grid': "0.125/0.125", \
                    'levtype': 'sfc', \
                    'param': "125.210/137.128/151.128/165.128/166.128/167.128/206.128/207.210/213.210/214.210/215.210/216.210", \
                    'step': "0", \
                    'stream': 'oper', \
                    'time': time, \
                    'type': type, \
                    'format': 'netcdf',
                    'area': "90/-180/-90/180",
                    'target': "/work/ALT/swot/aval/OBS2CO/git/grs/grsdata/ECMWF/CAMS/"+str(year)+"-"+str(month).zfill(2)+"_month_"+data_type+".nc"})
        #except:
        #    print('Error: not appropriate cams settings for download. Refers to ecmwf.')
        #    sys.exit()


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Downlad Cams datasets.')
    parser.add_argument('mode',
                help='the cams_forecast or cams_reanalysis mode')
    parser.add_argument('month',
                help='the month from which data will be downloaded. Set all for downloading data from 2000 to today. Set a date in %Y%m%d format if you want to download data from this date until now.')

    args = parser.parse_args()

    main(vars(args))

