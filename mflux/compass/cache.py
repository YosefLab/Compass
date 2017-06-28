import os
import json

from ..globals import RESOURCE_DIR
PREPROCESS_CACHE_DIR = os.path.join(RESOURCE_DIR, 'COMPASS')
if not os.path.isdir(PREPROCESS_CACHE_DIR):
    os.mkdir(PREPROCESS_CACHE_DIR)

# Keys are tuple of (model, media)
_cache = {}


def load(model, media=None):
    global _cache

    if media is None:
        media = model.media

    if not isinstance(model, str):
        model = model.name

    if (model, media) not in _cache:

        cache_file = os.path.join(PREPROCESS_CACHE_DIR,
                                  model, media, "preprocess.json")

        if os.path.exists(cache_file):
            with open(cache_file) as fin:
                out = json.load(fin)

            _cache[(model, media)] = out

        else:
            _cache[(model, media)] = {}

    return _cache[(model, media)]


def save(model, media=None):
    global _cache

    if media is None:
        media = model.media

    if not isinstance(model, str):
        model = model.name

    cache_data = _cache[(model, media)]

    cache_dir = os.path.join(PREPROCESS_CACHE_DIR, model, media)

    if not os.path.isdir(cache_dir):
        os.makedirs(cache_dir)

    cache_file = os.path.join(cache_dir, 'preprocess.json')

    with open(cache_file, 'w') as fout:
        json.dump(cache_data, fout, indent=1)


def clear(model, media=None):
    global _cache

    if media is None:
        media = model.media

    if not isinstance(model, str):
        model = model.name

    _cache[(model, media)] = {}

    save(model, media)
