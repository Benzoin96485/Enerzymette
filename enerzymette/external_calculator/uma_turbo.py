import os
from fairchem.core import pretrained_mlip, FAIRChemCalculator

_predict_unit = pretrained_mlip.get_predict_unit(
    model_name="uma-s-1p1", 
    device="cuda", 
    cache_dir=os.getenv("FAIRCHEM_CACHE_DIR"),
    inference_settings="turbo"
)
uma_calculator = FAIRChemCalculator(
    predict_unit=_predict_unit,
    task_name="omol",
)
