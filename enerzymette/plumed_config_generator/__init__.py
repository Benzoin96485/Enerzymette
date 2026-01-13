def get_plumed_patch(key: str):
    if key == "sammt":
        from . import sammt
        return sammt.__file__
