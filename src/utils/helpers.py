from argparse import Namespace

def dict_to_namespace(d):
    if isinstance(d, dict):
        return Namespace(**{k: dict_to_namespace(v) for k, v in d.items()})
    return d