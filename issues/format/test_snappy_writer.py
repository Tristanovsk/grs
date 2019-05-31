import numpy as np
from esasnappy import Product, ProductIO, ProductData

formats = ([".dim", "BEAM-DIMAP"],
           [".h5", "HDF5"],
           [".img", "GDAL-HFA-WRITER"],
           [".jp2", "JPEG2000"],
           [".bin", "Generic Binary BSQ"],
           [".tif", "GeoTIFF+XML"],
           [".grd", "GDAL-GSBG-WRITER"],
           [".nc", "NetCDF-CF"],
           [".nc", "NetCDF-BEAM"],
           [".nc", "NetCDF4-BEAM"],
           [".tif", "GeoTIFF"],
           ["", "JP2"],
           [".tif", "GDAL-GTiff-WRITER"],
           [".tif", "GeoTIFF-BigTIFF"],
           [".csv", "CSV"],
           [".bmp", "GDAL-BMP-WRITER"],
           [".rst", "GDAL-RST-WRITER"],
           [".nc", "NetCDF4-CF"],
           [".hdr", "ENVI"])

for format in formats:
    print()
    print('---------------')
    print(format)
    file = 'issues/test' + format[0]

    product = Product('test', 'test', 100, 100)
    band = product.addBand('SZA', ProductData.TYPE_FLOAT32)
    band.setModified(True)
    band.setNoDataValue(np.nan)
    band.setNoDataValueUsed(True)
    product.getBand('SZA').setDescription('Solar zenith angle in deg.')

    writer = ProductIO.getProductWriter(format[1])

    product.setProductWriter(writer)
    try:
        #product.writeHeader(file)
        ProductIO.writeProduct(product, file, format[1])
        p = ProductIO.readProduct(file)
        try:

            print(p.getBand('SZA').getDescription())
            print('SUCCESS!!!!!!!')
        except:
            print('error for format ', format[1])
    except:
        print('no header allowed for', format[1])
    print()
