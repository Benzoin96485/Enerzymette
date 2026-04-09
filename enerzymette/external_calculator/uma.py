import os
from fairchem.core import pretrained_mlip, FAIRChemCalculator

def get_uma_calculator(
    model="uma-s-1p2",
    task="omol",
    device="cuda",
    cache_dir=os.getenv("FAIRCHEM_CACHE_DIR"),
    mode="turbo"
):
    _predict_unit = pretrained_mlip.get_predict_unit(
        model_name=model,
        device=device,
        cache_dir=cache_dir,
        inference_settings=mode
    )
    return FAIRChemCalculator(
        predict_unit=_predict_unit,
        task_name=task
    )
