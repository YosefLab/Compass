import os
import json

from ..globals import RESOURCE_DIR
PREPROCESS_CACHE_DIR = os.path.join(RESOURCE_DIR, 'COMPASS')
if not os.path.isdir(PREPROCESS_CACHE_DIR):
    os.mkdir(PREPROCESS_CACHE_DIR)

_cache = {}


def load(model):
    global _cache

    if model.name not in _cache:

        cache_file = os.path.join(PREPROCESS_CACHE_DIR,
                                  model.name + ".preprocess")

        if os.path.exists(cache_file):
            with open(cache_file) as fin:
                out = json.load(fin)

            _cache[model.name] = out

        else:
            _cache[model.name] = {}

    return _cache[model.name]


def save(model):
    global _cache

    cache_data = _cache[model.name]

    cache_file = os.path.join(PREPROCESS_CACHE_DIR,
                              model.name + ".preprocess")

    with open(cache_file, 'w') as fout:
        json.dump(cache_data, fout, indent=1)


def clear(model):
    global _cache

    _cache[model.name] = {}

    save(model)
