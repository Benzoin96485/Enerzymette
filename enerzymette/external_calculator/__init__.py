def get_calculator_patch(key: str):
    if key == "uma":
        from . import uma
        return uma.__file__