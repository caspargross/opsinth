# opsinth/__init__.py

import yaml
from pathlib import Path

VERSION = "0.1"

def load_variants():
    """Load variant definitions from YAML file."""
    yaml_path = Path(__file__).parent.parent / 'data' / 'variants.yaml'
    with open(yaml_path) as f:
        variants = yaml.safe_load(f)
    return variants['color_variants'], variants['splicing_variants']

# Load variants from YAML file
OPSIN_COLOR_VARIANTS, OPSIN_SPLICING_VARIANTS = load_variants()
