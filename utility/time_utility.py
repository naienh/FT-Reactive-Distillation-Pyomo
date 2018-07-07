#python3
import atexit
from time import time, strftime, localtime
from datetime import timedelta, datetime

def create_filename_time():
    return datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss')

start = time()
midway = time()

def secondsToStr(elapsed=None):
    if elapsed is None:
        return strftime("%Y-%m-%d %H:%M:%S", localtime())
    else:
        delta = timedelta(seconds=round(elapsed))
        return str(delta)

def log(s, elapsed=None):
    line = "="*108
    print(line)
    print(secondsToStr(), '-', s)
    if elapsed:
        print("Elapsed time:", elapsed)
    print(line)

def log_now():
    global midway
    current = time()
    elapsed = current-midway
    midway = time()
    log("Elapsed time from last log", secondsToStr(elapsed))

def log_end():
    end = time()
    elapsed = end-start
    log("End Program", secondsToStr(elapsed))

atexit.register(log_end)
log("Start Program")
