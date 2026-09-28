from typing import Dict
from pprint import pformat
from pathlib import Path

import argparse
import xarray as xr
import numpy as np


def compare_attrs(attrs: Dict, attrs_ref: Dict):
    differences = {}
    # compare the attributes keys
    for attr in (attrs_ref.keys() - attrs.keys()):
        differences[attr] = f"attribute {attr} from Ref not in output"
    for attr in (attrs.keys() - attrs_ref.keys()):
        differences[attr] = f"attribute {attr} from output not in Ref"

    # compare the attributes values
    for attr in (attrs.keys() - differences.keys()):
        if attr == "processing_time":
            continue
        if isinstance(attrs[attr], np.ndarray):
            if not np.allclose(attrs[attr], attrs_ref[attr], equal_nan=True):
                differences[attr] = f"{attr} : ({attrs[attr]} | {attrs_ref[attr]})"
        elif attrs[attr] != attrs_ref[attr]:
            if Path(str(attrs[attr])).exists():
                if not Path(attrs[attr]).name == Path(attrs_ref[attr]).name:
                    differences[attr] = f"{attr} : ({attrs[attr]} | {attrs_ref[attr]})"
            else:
                differences[attr] = f"{attr} : ({attrs[attr]} | {attrs_ref[attr]})"

    return differences


def compare_nc(nc: xr.Dataset, nc_ref: xr.Dataset):
    # compare the list of coordinates
    if list(nc_ref.coords) != list(nc.coords):
        raise ValueError("Differences in the coordinates")

    # compare the list of data variables
    dv_differences = {}
    # compare the attributes keys
    for k in (nc_ref.keys() - nc.keys()):
        dv_differences[k] = f"data variable {k} from Ref not in output"
    for k in (nc.keys() - nc_ref.keys()):
        dv_differences[k] = f"data variable {k} from output not in Ref"

    # compare the global attributes
    attr_differences = compare_attrs(nc.attrs, nc_ref.attrs)

    # compare the values of the data variables
    for data_var in (nc.keys() - dv_differences.keys()):
        # compare the data variables attributes
        compare_attrs(nc[data_var].attrs, nc_ref[data_var].attrs)

        # compare the data variables values
        if not np.allclose(nc[data_var], nc_ref[data_var], equal_nan=True):
            diff_max = np.nanmax(np.abs(nc_ref[data_var].values - nc[data_var].values))
            dv_differences[data_var] = f"data variable {data_var} differences with ref up to {diff_max}"

    if attr_differences.keys() or dv_differences.keys():
        raise ValueError(pformat(attr_differences | dv_differences))


def main():

    parser = argparse.ArgumentParser(description="compare grs outputs")
    parser.add_argument("--nc-ref", help="reference cetcdf")
    parser.add_argument("--nc", help="grs output to compare")

    args = parser.parse_args()

    nc = xr.open_dataset(args.nc)
    nc_ref = xr.open_dataset(args.nc_ref)

    compare_nc(nc, nc_ref)
    print(f"{args.nc} file is equivalent to it's reference")


if __name__ == "__main__":
    main()
