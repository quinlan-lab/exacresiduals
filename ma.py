def match_alleles(alt,csq_alt,ref,impact,is_multi):
    """
    >>> match_alleles('GT','GT','G','frameshift_variant&feature_elongation',True)
    True
    >>> match_alleles('GT','T','G','regulatory_region_variant',True)
    False
    >>> match_alleles('GT','T','G','frameshift_variant&feature_elongation',False)
    True
    >>> match_alleles('T','-','TCAGGG','stop_gained&frameshift_variant', False) 
    True
    >>> match_alleles('AAC','-','A','frameshift_variant&feature_truncation', False) 
    True
    >>> match_alleles('AG','AG','A','frameshift_variant&feature_elongation', False)
    False
    >>> match_alleles('AG','-','A','frameshift_variant&missense_variant&feature_truncation', True)
    False
    >>> match_alleles('C','C','CGCGAGGACTCTGCCTCCCCA','frameshift_variant&feature_truncation',True)
    True
    >>> match_alleles('TCAGCTCC','CAGCTCC','T','stop_gained&frameshift_variant&feature_elongation',False)
    True
    >>> match_alleles('TCAGCTCC','CAGCTCC','T','stop_gained&frameshift_variant&feature_elongation',True)
    False
    >>> match_alleles('T','T','A','missense_variant',False)
    True
    >>> match_alleles('T','C','A','missense_variant',False)
    False
    >>> match_alleles('T','C','A','missense_variant',True)
    False
    """
    if "frameshift" in impact:
        if is_multi == True:
            if alt == csq_alt:
                return True
            else:
                return False
        else:
            if alt == csq_alt:
                if "truncation" in impact or "stop_gained" in impact or "elongation" in impact:
                    return False
            else:
                if "truncation" in impact or ("stop_gained" in impact and "elongation" not in impact):
                    if csq_alt == "-":
                        return True
                    else:
                        return False
                elif "elongation" in impact:
                    if alt == ref+csq_alt:
                        return True
                    else:
                        return False
                else:
                    return False
    else:
        if alt == csq_alt:
            return True
        else:
            return False

if __name__ == "__main__":
    import doctest
    print(doctest.testmod())
