from __future__ import division
import sys, time
from datetime import timedelta

# TODO
# timing should weigh first iterations less than subsequent

def update_progress(done_ratio, 
                    time_passed = None, time_left = None):
    """
    Takes a number in the range [0,1] and outputs a progressbar.
    """
    bar = '#' * int(done_ratio * 20)
    string = "\r[{0:20s}] {1:.1%}". \
                format(bar, done_ratio)
    if time_passed is not None:
        time_str = str(timedelta(seconds=time_passed))
        string += " Time elapsed: {}".format(time_str)
    if time_left is not None:
        time_str = str(timedelta(seconds=time_left))
        string += " Time left: {}".format(time_str)
    sys.stdout.write(string)
    sys.stdout.flush()
    if done_ratio >= 1:
        print # newline when done

def progress(iterator, timing = True):
    """
    Makes a progressbar to be used with a for-loop, 
    where each iteration takes the same amount of time.
    """
    length = len(iterator)
    t = time.time()
    time_passed = 0
    for (i, elem) in enumerate(iterator):
        if timing:
            time_passed += time.time() - t
            time_left = time_passed / (i+1) * (length - (i+1))
            t = time.time()
        update_progress((i+1)/length, 
                        time_passed = int(time_passed),
                        time_left = int(time_left))
        yield elem
        
