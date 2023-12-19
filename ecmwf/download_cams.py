"""
This function download the cams data used in grs from CDS API

The main program takes two argument : start  and end year
"""

import argparse
from calendar import monthrange
from datetime import datetime, timedelta
from typing import Optional

from ecmwfapi import ECMWFService
from cdsapi import Client
from pathlib import Path
from shutil import rmtree
import xarray as xr

ECMWF_START_DAY = 2
ECMWF_END_DAY = 1
COPERNICUS_START_DAY = 2
COPERNICUS_END_DAY = 1

ECMWF_REQUEST = {
    'class': "od",
    'date': None,
    'expver': "1",
    'levelist': "1/2/3/5/7/10/20/30/50/70/100/150/200/250/300/400/500/600/700/850/925/1000",
    'levtype': "pl",
    'param': "157.128",
    'step': "0",
    'stream': "oper",
    'time': "00/12",
    # 'time': ['00:00', '03:00', '06:00', '09:00', '12:00', '18:00', '21:00'],
    #    'leadtime_hour': ['0', '3', '6', '9'],
    'type': "fc",
    'grid': "0.4/0.4",
    'format': "netcdf",
}
COPERNICUS_FC_REQUEST = {
    'nocache': '456',
    'format': 'netcdf',
    'date': None,
    'time': ['00:00', '12:00'],
    'leadtime_hour': ['0', '3', '6', '9'],
    # 'leadtime_hour': ['0', '12', '18', '21', '3', '6', '9', ],
    # 'time': '00:00',
    'type': 'forecast',
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
}

COPERNICUS_EAC4_REQUEST = {
    'nocache': '456',
    'format': 'netcdf',
    'date': None,
    'time': ['00:00', '03:00', '06:00', '09:00', '12:00', '15:00', '18:00', '21:00'],
    'variable': [
        '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_temperature',
        'mean_sea_level_pressure', 'surface_pressure',
        'total_aerosol_optical_depth_469nm',
        'total_aerosol_optical_depth_550nm',
        'total_aerosol_optical_depth_670nm',
        'total_aerosol_optical_depth_865nm',
        'total_aerosol_optical_depth_1240nm',
        'black_carbon_aerosol_optical_depth_550nm',
        'dust_aerosol_optical_depth_550nm',
        'organic_matter_aerosol_optical_depth_550nm',
        'sea_salt_aerosol_optical_depth_550nm',
        'sulphate_aerosol_optical_depth_550nm',
        'total_column_carbon_monoxide', 'total_column_methane', 'total_column_nitrogen_dioxide',
        'total_column_ozone', 'total_column_water_vapour',
    ],
}


def valid_dir(outdir) -> Path:
    return Path(outdir).resolve(strict=True)


def valid_month(m) -> Optional[int]:
    try:
        if m is not None:
            int_m = int(m)
            if not 1 <= int_m <= 12:
                raise argparse.ArgumentTypeError(f"not a valid month: {m}. Please enter a value between 1 and 12")
            return int(m)
        else:
            return None
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a valid month: {m}. Please enter a value between 1 and 12")


def valid_year(y) -> Optional[int]:
    try:
        if y is not None:
            int_y = int(y)
            if not 2015 <= int_y <= datetime.now().year:
                raise argparse.ArgumentTypeError(
                    f"not a valid year: {y}. Please enter a value between 2015 and {datetime.now().year}")
            return int(y)
        else:
            return None
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a valid month: {y}. Please enter a value between 1 and 12")


def valid_date(s) -> Optional[datetime]:
    try:
        if s is not None:
            d = datetime.strptime(s, "%Y-%m-%d")
            if not datetime(year=2015, month=1, day=1) <= d <= datetime.now():
                raise argparse.ArgumentTypeError(
                    f"not a valid date: {s}. "
                    f"Please enter a date between 2015 and {datetime.now().strftime('%Y-%m-%d')}")
            return d
        else:
            return None
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a valid date: {s}")


def sdate(d: datetime) -> str:
    return d.strftime('%Y-%m-%d')


def get_request_date(args, source) -> tuple[datetime, datetime]:
    start_date = None
    end_date = None
    start_default = COPERNICUS_START_DAY
    end_default = COPERNICUS_END_DAY

    if args.month is not None:
        if args.year is None:
            y = datetime.now().year
        else:
            y = args.year
        start_date = datetime(year=y, month=args.month, day=1)
        end_date = datetime(year=y, month=args.month, day=monthrange(y, args.month)[1])

    if source == "copernicus":
        if not start_date and args.cstart is not None:
            start_date = args.cstart
        if not end_date and args.cend is not None:
            end_date = args.cend
    elif source == "ecmwf":
        start_default = ECMWF_START_DAY
        end_default = ECMWF_END_DAY
        if not start_date and args.estart is not None:
            start_date = args.estart
        if not end_date and args.eend is not None:
            end_date = args.eend

    if start_date is None:
        start_date = datetime.today() - timedelta(days=start_default)

    if end_date is None:
        end_date = datetime.today() - timedelta(days=end_default)

    if start_date > end_date:
        raise argparse.ArgumentTypeError(f"wrong date range "
                                         f"start : {sdate(start_date)} end : {sdate(end_date)}")
    return start_date, end_date


def retrieve_ecmwf(args) -> Path:
    estart, eend = get_request_date(args, source="ecmwf")
    ECMWF_REQUEST["date"] = f"{sdate(estart)}/to/{sdate(eend)}"
    target_name = f"{sdate(estart)}_{sdate(eend)}_relative-humidity-forecast.nc"

    if args.source == "both" and args.combine:
        target_p = Path(args.outdir, "tmp", target_name).resolve()
        if target_p.exists():
            target_p.unlink()
    else:
        target_p = Path(args.outdir, str(estart.year), str(estart.month).zfill(2), str(estart.day).zfill(2),
                        target_name).resolve()
    target_p.parent.mkdir(exist_ok=True, parents=True)
    if target_p.exists():
        if args.overwrite:
            print(f"overwriting existing file {target_p}")
        else:
            raise FileExistsError(str(target_p))

    server = ECMWFService("mars")
    server.execute(ECMWF_REQUEST, str(target_p))

    return target_p


def retrieve_copernicus(args) -> Path:
    if args.mode == "reanalisys":
        data_type = 'cams-global-reanalysis-eac4'
        cop_req = COPERNICUS_EAC4_REQUEST
    else:
        data_type = 'cams-global-atmospheric-composition-forecasts'
        cop_req = COPERNICUS_FC_REQUEST

    cstart, cend = get_request_date(args, source="copernicus")
    cop_req["date"] = f"{sdate(cstart)}/{sdate(cend)}"
    target_name = f"{sdate(cstart)}_{sdate(cend)}-{data_type}.nc"

    if args.source == "both" and args.combine:
        target_p = Path(args.outdir, "tmp", target_name).resolve()
        if target_p.exists():
            target_p.unlink()
    else:
        target_p = Path(args.outdir, str(cstart.year), str(cstart.month).zfill(2), str(cstart.day).zfill(2),
                        target_name).resolve()

    target_p.parent.mkdir(exist_ok=True, parents=True)
    if target_p.exists():
        if args.overwrite:
            print(f"overwriting existing file {target_p}")
        else:
            raise FileExistsError(str(target_p))
    Client().retrieve(data_type, cop_req, target_p)

    return target_p


def combine_cop_ecmwf(args, ecmwf_cams_p, copernicus_cams_p):
    ecmwf_ds = xr.open_dataset(ecmwf_cams_p)
    copernicus_ds = xr.open_dataset(copernicus_cams_p)
    combined_ds = xr.combine_by_coords([ecmwf_ds, copernicus_ds], combine_attrs="drop_conflicts")
    combined_ds.attrs["history"] = f"{ecmwf_ds.attrs['history']} | {copernicus_ds.attrs['history']}"

    estart, eend = get_request_date(args, source="ecmwf")
    cstart, cend = get_request_date(args, source="copernicus")
    start_date = min(estart, cstart)
    end_date = max(eend, cend)
    if args.mode == "reanalisys":
        data_type = 'cams-global-reanalysis-eac4'
    else:
        data_type = 'cams-global-atmospheric-composition-forecasts'

    target_name = f"{sdate(start_date)}_{sdate(end_date)}_{data_type}_relative_humidity.nc"
    target_p = Path(args.outdir, str(start_date.year), str(start_date.month).zfill(2), str(start_date.day).zfill(2),
                    target_name).resolve()
    target_p.parent.mkdir(exist_ok=True, parents=True)
    # combined_ds = combined_ds.interpolate_na(dim="time", fill_value="extrapolate")
    combined_ds.to_netcdf(target_p)


def main(args):
    copernicus_cams_p = None
    if args.month is not None or args.source in ["copernicus", "both"]:
        copernicus_cams_p = retrieve_copernicus(args)

    ecmwf_cams_p = None
    if args.source in ["ecmwf", "both"]:
        ecmwf_cams_p = retrieve_ecmwf(args)

    if args.source == "both" and args.source == "combine":
        combine_cop_ecmwf(args, ecmwf_cams_p, copernicus_cams_p)
        rmtree(Path(args.outdir, "tmp"), ignore_errors=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download Cams datasets from CDS")
    parser.add_argument("--mode", choices=["forecast", " "],
                        default="forecast",
                        help="for copernicus source, choose between `reanalysis` and `forecast` dataset")
    parser.add_argument("--source", choices=["ecmwf", "copernicus", "both"],
                        default="both",
                        help="choose `reanalysis` or `forecast` dataset")
    parser.add_argument("--outdir", "-o",
                        help="output directory where the products will be downloaded",
                        default="/work/datalake/watcal/ECMWF/CAMS/")
    parser.add_argument("--cstart",
                        help="copernicus starting date for dataset (default is today -7), format:YYYY-MM-DD",
                        default=None,
                        type=valid_date)
    parser.add_argument("--cend",
                        help="copernicus end date for dataset",
                        default=None,
                        type=valid_date)
    parser.add_argument("--estart",
                        help="ecmwf starting date for dataset (default is today -2), format:YYYY-MM-DD",
                        default=None,
                        type=valid_date)
    parser.add_argument("--eend",
                        help="ecmwf end date for dataset",
                        default=None,
                        type=valid_date)
    parser.add_argument("--month",
                        help="get data from ADS for the entire month. supersede any date provided",
                        type=valid_month)
    parser.add_argument("--year",
                        help="to use with --month to define the year. default to the current year",
                        type=valid_year)
    parser.add_argument("--overwrite", action="store_true",
                        help="any existing data will be overwritten")
    parser.add_argument("--combine", action="store_true",
                        help="used with --source=both will try to combine the data into one netcdf file")

    main(parser.parse_args())
