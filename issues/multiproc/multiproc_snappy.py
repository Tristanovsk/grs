
import os, sys
from multiprocessing import Pool

def exit():
    os._exit(0)

def call_snappy(i):


    print(i,' call of snappy')
    sys.exit()
    import esasnappy
    exit()
    #sys.exit()

if __name__ == '__main__':
    with Pool(processes=3) as pool:
        processes = pool.map(call_snappy, range(30),1)
        print('processes ',processes)
        pool.close()
        # pool.join()
    print('end')