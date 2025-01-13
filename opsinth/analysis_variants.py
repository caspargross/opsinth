import logging
from opsinth.utils import query_to_ref_pos
from opsinth.config import *


def genotype_known_variants(alignment_dict: dict, ref_seq: str, variants: list) -> list:
    """
    Find variants in an aligned sequence.
    
    Args:
        alignment_dict: Dictionary containing alignment information
        ref_seq: Reference sequence to compare against
        variants: List of variant dictionaries
        
    Returns:
        List of dictionaries with variant information
    """
    results = []
    
    logger = logging.getLogger('opsinth.analysis_variants')
    logger.debug(f"Reference sequence length: {len(ref_seq)}")
    logger.debug(f"Alignment region: {alignment_dict['ref_start']}-{alignment_dict['ref_end']}")
    logger.debug(f"CIGAR string: {alignment_dict['cigar']}")
    
    for var in variants:
        try:
            # Adjust cDNA position based on query start
            adjusted_cdna_pos = var['cdna_pos'] + alignment_dict['query_start']
            
            # Get reference position (0-based)
            ref_pos = query_to_ref_pos(
                adjusted_cdna_pos - 1,  # Convert to 0-based coordinate
                alignment_dict['cigar'],
                query_start=alignment_dict['query_start'],
                ref_start=alignment_dict['ref_start']
            )
            
            # Add 1 to get correct base
            ref_pos = ref_pos + 1
            
            logger.debug(f"Processing variant at cDNA pos {var['cdna_pos']} (adjusted to {adjusted_cdna_pos}):")
            logger.debug(f"Converted to reference pos: {ref_pos}")
            
            if ref_pos >= len(ref_seq):
                logger.warning(f"Position {ref_pos} is beyond reference sequence length {len(ref_seq)}")
                continue
                
            observed_base = ref_seq[ref_pos]
            is_alt = observed_base == var['alt_base']
            is_ref = observed_base == var['ref_base']
            
            # Determine the type of variant
            type = var['ref_type'] if is_ref else (var['alt_type'] if is_alt else "Unknown")
            
            result = {
                'cdna_pos': var['cdna_pos'],
                'ref_pos': ref_pos,
                'expected_ref': var['ref_base'],
                'expected_alt': var['alt_base'],
                'observed': observed_base,
                'aa_pos': var['aa_pos'],
                'aa_change': var['aa_change'],
                'is_alt': is_alt,
                'is_ref': is_ref,
                'ref_type': var['ref_type'],
                'alt_type': var['alt_type'],
                'type': type
            }
            results.append(result)
            
        except ValueError as e:
            logger.error(f"Could not map position {var['cdna_pos']}: {str(e)}")
            
    return results

def format_variant_classification(variants: list) -> str:
    """
    Format variant results as a readable table.
    
    Args:
        variants: List of variant dictionaries from genotype_known_variants
        
    Returns:
        Formatted string containing the variant table
    """
    if not variants:
        return "No variants found"
        
    lines = []
    lines.append("\nVariant Analysis:")
    lines.append("Pos cDNA | Pos Ref | AA Change | Expected(Ref/Alt) | Observed | Type")
    lines.append("-" * 70)
    
    for var in variants:
        variant_type = var['ref_type'] if var['is_ref'] else (var['alt_type'] if var['is_alt'] else "Unknown")
        lines.append(
            f"cDNA {var['cdna_pos']} |  {var['ref_pos']: >6} | {var['aa_change']}       | "
            f"{var['expected_ref']}/{var['expected_alt']}               | "
            f"{var['observed']}        | {variant_type}"
        )
    
    return "\n".join(lines)

