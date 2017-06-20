import os
import json

from ..globals import RESOURCE_DIR
PREPROCESS_CACHE_DIR = os.path.join(RESOURCE_DIR, 'COMPASS')
if not os.path.isdir(PREPROCESS_CACHE_DIR):
    os.mkdir(PREPROCESS_CACHE_DIR)

# Keys are tuple of (model.name, model.media)
_cache = {}


def load(model):
    global _cache

    if (model.name, model.media) not in _cache:

        cache_file = os.path.join(PREPROCESS_CACHE_DIR,
                                  model.name, model.media, "preprocess.json")

        if os.path.exists(cache_file):
            with open(cache_file) as fin:
                out = json.load(fin)

            _cache[(model.name, model.media)] = out

        else:
            _cache[(model.name, model.media)] = {}

    return _cache[(model.name, model.media)]


def save(model):
    global _cache

    cache_data = _cache[(model.name, model.media)]

    cache_dir = os.path.join(PREPROCESS_CACHE_DIR, model.name, model.media)
    os.makedirs(cache_dir)
    cache_file = os.path.join(cache_dir, 'preprocess.json')

    with open(cache_file, 'w') as fout:
        json.dump(cache_data, fout, indent=1)


def clear(model):
    global _cache

    _cache[(model.name, model.media)] = {}

    save(model)
