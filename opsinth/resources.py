"""
Default files and references
"""

OPN1MW_REF = "opsinth/resources/OPN1MW.cDNA.GRCh38.fa"
OPN1LW_REF = "opsinth/resources/OPN1LW.cDNA.GRCh38.fa"

DEFAULT_ANCHORS = "opsinth/resources/anchors.WGS.GRCh38.fa"
DEFAULT_REF = "opsinth/resources/GRCh38.chrX.roi.fa.gz"
DEFAULT_BED = "opsinth/resources/OPN1_region.GRCh38.bed"
DEFAULT_OUT = "opsinth"

"""
Module for handling resource loading and management.
"""

import yaml
import logging
from pathlib import Path

def get_data_dir() -> Path:
    """Get the path to the data directory."""
    return Path(__file__).parent.parent / 'data'

def load_variants():
    """
    Load variant definitions from YAML file.
    
    Returns:
        Tuple containing:
        - color_variants: List of color-determining variants
        - splicing_variants: List of splicing variants
    """
    yaml_path = get_data_dir() / 'variants.yaml'
    
    try:
        with open(yaml_path) as f:
            variants = yaml.safe_load(f)
            
        logging.debug(f"Loaded variants from {yaml_path}")
        return variants['color_variants'], variants['splicing_variants']
        
    except FileNotFoundError:
        logging.error(f"Variants file not found: {yaml_path}")
        raise
    except yaml.YAMLError as e:
        logging.error(f"Error parsing variants YAML: {e}")
        raise
    except KeyError as e:
        logging.error(f"Missing required variant type in YAML: {e}")
        raise

# Load variants when module is imported
OPSIN_EXON5_COLOR_VARIANTS, OPSIN_EXON3_SPLICING_VARIANTS = load_variants() 