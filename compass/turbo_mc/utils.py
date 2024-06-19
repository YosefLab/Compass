import time

TIMEIT = False


def timeit(prefix, always_time: bool = False):
    def decorate(func):
        def call(*args, **kwargs):
            if TIMEIT or always_time:
                print(f"{prefix}::{func.__name__} ...")
                start_time = time.time()
            result = func(*args, **kwargs)
            if TIMEIT or always_time:
                tot_time = time.time() - start_time
                print(f"{prefix}::{func.__name__}, time = %.2f" % tot_time)
            return result
        return call
    return decorate
