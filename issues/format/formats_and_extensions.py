from esasnappy import jpy

ProductIOPlugInManager = jpy.get_type('org.esa.snap.core.dataio.ProductIOPlugInManager')
ProductReaderPlugIn = jpy.get_type('org.esa.snap.core.dataio.ProductReaderPlugIn')
ProductWriterPlugIn = jpy.get_type('org.esa.snap.core.dataio.ProductWriterPlugIn')

read_plugins = ProductIOPlugInManager.getInstance().getAllReaderPlugIns()
write_plugins = ProductIOPlugInManager.getInstance().getAllWriterPlugIns()

print('Writer formats')
while write_plugins.hasNext():
    plugin = write_plugins.next()
    print('["'+ plugin.getDefaultFileExtensions()[0]+'","'+ plugin.getFormatNames()[0]+'"],')

print('Reader formats')
while read_plugins.hasNext():
    plugin = read_plugins.next()
    print('  ', plugin.getDefaultFileExtensions()[0], plugin.getFormatNames()[0])