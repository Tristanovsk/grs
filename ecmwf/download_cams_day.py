"""
This function download the cams data used in grs from CDS API

The main program takes two argument : start  and end year
"""

import argparse
from datetime import datetime, timedelta
from typing import Optional

from ecmwfapi import ECMWFService
from cdsapi import Client
from pathlib import Path
from shutil import rmtree
import xarray as xr

ECMWF_START_DAY = 2
ECMWF_END_DAY = 1
COPERNICUS_START_DAY = 7
COPERNICUS_END_DAY = 6
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
    'type': "fc",
    'grid': "0.4/0.4",
    'format': "netcdf",
}
COPERNICUS_FC_REQUEST = {
    'nocache': '456',
    'format': 'netcdf',
    'date': None,
    'leadtime_hour': ['0', '12', '18', '21', '3', '6', '9', ],
    'time': '00:00',
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


def valid_date(s) -> Optional[datetime]:
    try:
        if s is not None:
            return datetime.strptime(s, "%Y-%m-%d")
        else:
            return None
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a valid date: {s}")


def sdate(d: datetime) -> str:
    return d.strftime('%Y-%m-%d')


def get_request_date(start_date, end_date, start_default, end_default) -> tuple[datetime, datetime]:
    if start_date is None:
        estart = datetime.today() - timedelta(days=start_default)
    else:
        estart = start_date

    if end_date is None:
        eend = datetime.today() - timedelta(days=end_default)
    else:
        eend = end_date

    if estart > eend:
        raise argparse.ArgumentTypeError(f"wrong date range for ECMWF "
                                         f"start : {sdate(estart)} end : {sdate(eend)}")
    return estart, eend


def retrieve_ecmwf(args) -> Path:
    estart, eend = get_request_date(args.estart, args.eend, ECMWF_START_DAY, ECMWF_END_DAY)
    ECMWF_REQUEST["date"] = f"{sdate(estart)}/to/{sdate(eend)}"
    target_name = f"{sdate(estart)}_{sdate(eend)}-relative-humidity-forecast.nc"

    if args.source == "both":
        target_p = Path(args.outdir, "tmp", target_name).resolve()
        if target_p.exists():
            target_p.unlink()
    else:
        target_p = Path(args.outdir, str(estart.year), str(estart.month), str(estart.day), target_name).resolve()
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

    cstart, cend = get_request_date(args.cstart, args.cend, COPERNICUS_START_DAY, COPERNICUS_END_DAY)
    cop_req["date"] = f"{sdate(cstart)}/{sdate(cend)}"
    target_name = f"{sdate(cstart)}_{sdate(cend)}-{data_type}.nc"

    if args.source != "both":
        target_p = Path(args.outdir, "tmp", target_name).resolve()
        if target_p.exists():
            target_p.unlink()
    else:
        target_p = Path(args.outdir, str(cstart.year), str(cstart.month), str(cstart.day), target_name).resolve()

    target_p.parent.mkdir(exist_ok=True)
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

    estart, eend = get_request_date(args.estart, args.eend, ECMWF_START_DAY, ECMWF_END_DAY)
    cstart, cend = get_request_date(args.cstart, args.cend, COPERNICUS_START_DAY, COPERNICUS_END_DAY)
    start_date = min(estart, cstart)
    end_date = max(eend, cend)
    if args.mode == "reanalisys":
        data_type = 'cams-global-reanalysis-eac4'
    else:
        data_type = 'cams-global-atmospheric-composition-forecasts'

    target_name = f"{sdate(start_date)}_{sdate(end_date)}_{data_type}_relative_humidity.nc"
    target_p = Path(args.outdir,
                    str(start_date.year), str(start_date.month), str(start_date.day), target_name).resolve()
    combined_ds.to_netcdf(target_p)


def main(args):
    ecmwf_cams_p = None
    if args.source in ["ecmwf", "both"]:
        ecmwf_cams_p = retrieve_ecmwf(args)

    copernicus_cams_p = None
    if args.source in ["copernicus", "both"]:
        copernicus_cams_p = retrieve_copernicus(args)

    if args.source == "both":
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
                        default="/datalake/watcal/ECMWF/CAMS/")
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
    parser.add_argument("--overwrite", action="store_true",
                        help="any existing data will be overwritten")

    main(parser.parse_args())
