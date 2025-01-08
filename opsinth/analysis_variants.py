def find_variants(alignment_dict: dict, ref_seq: str, variants: list) -> list:
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
    
    logging.debug(f"Reference sequence length: {len(ref_seq)}")
    logging.debug(f"Alignment region: {alignment_dict['ref_start']}-{alignment_dict['ref_end']}")
    logging.debug(f"CIGAR string: {alignment_dict['cigar']}")
    
    for var in variants:
        try:
            # Get reference position (0-based)
            ref_pos = query_to_ref_pos(
                var['cdna_pos'] - 1,  # Convert to 0-based coordinate
                alignment_dict['cigar'],
                query_start=alignment_dict['query_start'],
                ref_start=alignment_dict['ref_start']
            )
            
            # Add 1 to get correct base
            ref_pos = ref_pos + 1
            
            logging.debug(f"Processing variant at cDNA pos {var['cdna_pos']}:")
            logging.debug(f"Converted to reference pos: {ref_pos}")
            
            if ref_pos >= len(ref_seq):
                logging.warning(f"Position {ref_pos} is beyond reference sequence length {len(ref_seq)}")
                continue
                
            observed_base = ref_seq[ref_pos]
            
            result = {
                'cdna_pos': var['cdna_pos'],
                'ref_pos': ref_pos,
                'expected_ref': var['ref_base'],
                'expected_alt': var['alt_base'],
                'observed': observed_base,
                'aa_pos': var['aa_pos'],
                'aa_change': var['aa_change'],
                'is_alt': observed_base == var['alt_base'],
                'is_ref': observed_base == var['ref_base'],
                'ref_type': var['ref_type'],
                'alt_type': var['alt_type']
            }
            results.append(result)
            
        except ValueError as e:
            logging.error(f"Could not map position {var['cdna_pos']}: {str(e)}")
            
    return results