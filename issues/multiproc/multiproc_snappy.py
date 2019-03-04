
import sys
from multiprocessing import Pool



def call_snappy(i):


    print(i,' call of snappy')
    import esasnappy
    sys.exit()

if __name__ == '__main__':
    with Pool(processes=3) as pool:
        pool.map(call_snappy, range(10))

