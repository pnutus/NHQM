"""import multiprocessing
import ctypes
import numpy as np

def f(i,j):
    return i*j
size = 10
shared_array_base = multiprocessing.Array(ctypes.c_double, size*size)
shared_array = np.ctypeslib.as_array(shared_array_base.get_obj())
shared_array = shared_array.reshape(size, size)

# No copy was made
assert shared_array.base.base is shared_array_base.get_obj()

# Parallel processing
def my_func(tup, f, def_param=shared_array):
    i,j = tup
    shared_array[i,j] = f(i,j)

if __name__ == '__main__':
    pool = multiprocessing.Pool(processes=4)
    pool.map(my_func, range(10), f)

    print shared_array
    
"""
    
    
import sharedmem as shm
import numpy as np
import multiprocessing as mp

def worker(q,arr):
    done = False
    while not done:
        cmd = q.get()
        if cmd == 'done':
            done = True
        elif cmd == 'data':
            ##Fake data. In real life, get data from hardware.
            rnd=np.random.randint(100)
            print('rnd={0}'.format(rnd))
            arr[:]=rnd
        q.task_done()

if __name__=='__main__':
    N=10
    arr=shm.zeros(N,dtype=np.uint8)
    q=mp.JoinableQueue()    
    proc = mp.Process(target=worker, args=[q,arr])
    proc.daemon=True
    proc.start()

    for i in range(3):
        q.put('data')
        # Wait for the computation to finish
        q.join()   
        print arr.shape
        print(arr)
    q.put('done')
    proc.join()    