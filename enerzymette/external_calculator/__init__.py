def get_calculator_patch(key: str):
    if key == "uma":
        from . import uma
        return uma.__file__
    elif key == "uma_turbo":
        from . import uma_turbo
        return uma_turbo.__file__