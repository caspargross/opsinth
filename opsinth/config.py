"""
Module for handling resource loading and management.
"""

import os
import yaml
import logging
import pkg_resources

def get_resource_path(filename):
    """Get absolute path to a resource file."""
    return pkg_resources.resource_filename('opsinth', os.path.join('resources', filename))

def load_variants(variant_file):
    """Load variant definitions from YAML file."""
    yaml_path = get_resource_path(variant_file)
    
    try:
        with open(yaml_path) as f:
            variants = yaml.safe_load(f)
            
        logging.debug(f"Loaded variants from {yaml_path}")
        return variants
        
    except FileNotFoundError:
        logging.error(f"Variants file not found: {yaml_path}")
        raise
    except yaml.YAMLError as e:
        logging.error(f"Error parsing variants YAML: {e}")
        raise
    except KeyError as e:
        logging.error(f"Missing required variant type in YAML: {e}")
        raise

# Default paths configuration
DEFAULTS = {
    'out': 'opsinth_out',
    'bed': get_resource_path('OPN1_region.GRCh38.bed'),
    'ref': get_resource_path('GRCh38.chrX.roi.fa.gz'),
    'anchors': get_resource_path('anchors.WGS.GRCh38.fa'),
    'racon': 'racon',
    'n_polish_rounds': 6
}

# Reference sequences
REFERENCES = {
    'opn1mw_cdna': get_resource_path('OPN1MW.cDNA.GRCh38.fa'),
    'opn1lw_cdna': get_resource_path('OPN1LW.cDNA.GRCh38.fa'),
    'opn1mw': get_resource_path('OPN1MW.GRCh38.fa'),
    'opn1lw': get_resource_path('OPN1LW.GRCh38.fa')
}

VARIANTS = {
    'exon5_color_variants':  load_variants('haplotype_defining_variants.yaml').get('exon5_color_variants', []),
    'exon3_splicing_variants':  load_variants('haplotype_defining_variants.yaml').get('exon3_splicing_variants', [])
}

# Make everything available at module level
__all__ = [
	'DEFAULTS',
	'REFERENCES',
	'VARIANTS'
] 