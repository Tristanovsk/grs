import multiprocessing as mp

def test(img):
    from grs import grs_process

    print(img)

imgs =['1','2','3','4','5','6','7','8']
p = mp.Pool(3)
p.map(test,imgs)