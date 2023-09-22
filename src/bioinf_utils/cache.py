import gzip
import pickle
import hashlib

import numpy


def array_to_hasheable_tuple(array, subsample_size=10000):
    array_size = array.size
    if array_size <= subsample_size:
        b = array.flat[:]
    else:
        rng = numpy.random.RandomState(42)
        idx = rng.randint(low=0, high=array.size, size=subsample_size)
        b = array.flat[idx]
    b.flags.writeable = False
    return (array.shape, b.data.tobytes())


def hash_from_tuple(tuple_):
    str_tuple = tuple(map(str, tuple_))
    return hashlib.md5(" ".join(str_tuple).encode()).hexdigest()


class MissingCachedResult(RuntimeError):
    pass


def save_cache(value_, cache_path, use_gzip=False):
    if use_gzip:
        fhand = gzip.open(cache_path, "wb")
    else:
        fhand = open(cache_path, "wb")
    pickle.dump(value_, fhand)


def load_cache(cache_path):
    if not cache_path.exists():
        raise MissingCachedResult()
    try:
        return pickle.load(gzip.open(cache_path, "rb"))
    except gzip.BadGzipFile:
        return pickle.load(open(cache_path, "rb"))


def get_result(
    funct,
    cache_path,
    args=None,
    kwargs=None,
    use_gzip=False,
    update_cache=False,
):
    if not update_cache and cache_path.exists():
        return load_cache(cache_path)

    if args is None:
        args = tuple()
    if kwargs is None:
        kwargs = {}

    result = funct(*args, **kwargs)
    save_cache(result, cache_path=cache_path, use_gzip=use_gzip)
    return result
