
import sys

def grs_call(self,p):
        args,fjunk = p
        for arg in args:
            file_tbp, outfile, wkt, altitude, aerosol, aeronet_file, resolution, \
            aot550, angstrom, unzip, untar, startrow, angleonly = arg
            print('yop',file_tbp)
            #return
            try:
                from grs import grs_process
                grs_process.Process().execute(file_tbp, outfile, wkt, altitude=altitude, aerosol=aerosol,
                                              dem=None, aeronet_file=aeronet_file, resolution=resolution,
                                              aot550=aot550, angstrom=angstrom, unzip=unzip, untar=untar,
                                              startrow=startrow, angleonly=angleonly)
            except:
                print('-------------------------------')
                print('error for file  ', file_tbp, ' skip')
                print('-------------------------------')
                with open(fjunk, "a") as myfile:
                    myfile.write(file_tbp + ' error during grs \n')
                continue
        # here sys.exit instead of "return" to terminate and close snappy and free memory
        #sys.exit()
        return

if __name__ == '__main__':
    p = sys.argv[1]
    grs_call(p)