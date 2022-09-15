
import os, glob
import numpy as np
import xarray as xr
import xml.etree.ElementTree as ET


from collections import namedtuple


opj =os.path.join
imageSAFE = "test/data/S2B_MSIL1C_20180927T103019_N0206_R108_T31TGK_20180927T143835.SAFE"
xml_file=glob.glob(opj(imageSAFE,'GRANULE','*','MTD_TL.xml'))[0]
with open(xml_file) as xml_:
    tree = ET.parse(xml_)
    root = tree.getroot()


# Internal parsing function for angular grids
def parse_angular_grid_node(node):
    values = []
    for c in node.find('Values_List'):
        values.append(np.array([float(t) for t in c.text.split()]))
        values_array = np.stack(values)
    return values_array

# Parse sun angles
Angles = namedtuple('Angles', 'zenith azimuth')

sun_angles = Angles(
    parse_angular_grid_node(
        root.find('.//Angles_Grids_List/Sun_Angles_Grids/Zenith')),
    parse_angular_grid_node(
        root.find(
            './/Angles_Grids_List/Sun_Angles_Grids/Azimuth')))

# Parse incidence angles
self.incidence_angles = {}
for b in root.find(
        './/Angles_Grids_List/Viewing_Incidence_Angles_Grids_List'
):
    if b.attrib['band_id'] != 'B1':
        band_key = self.Band(b.attrib['band_id'])
        band_dict = {}
        for d in b.findall('Viewing_Incidence_Angles_Grids'):
            det_key = self.Detector(int(d.attrib['detector_id']))
            zen = parse_angular_grid_node(d.find('Zenith'))
            az = parse_angular_grid_node(d.find('Azimuth'))
            band_dict[det_key] = Angles(zen, az)
        self.incidence_angles[band_key] = band_dict



