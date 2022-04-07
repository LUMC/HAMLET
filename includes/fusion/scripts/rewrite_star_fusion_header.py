#!/usr/bin/env python3

import sys

# Excludeds the fields est_J and est_S for compatibility with crimson, see
# https://github.com/bow/crimson/pull/18
output_fields =  ['#FusionName', 'JunctionReadCount', 'SpanningFragCount',
        'SpliceType', 'LeftGene', 'LeftBreakpoint',
        'RightGene', 'RightBreakpoint', 'JunctionReads', 'SpanningFrags',
        'LargeAnchorSupport', 'FFPM', 'LeftBreakDinuc', 'LeftBreakEntropy',
        'RightBreakDinuc', 'RightBreakEntropy', 'annots']

if __name__ == '__main__':
    with open(sys.argv[1]) as fin:
        # Input header
        header = next(fin).strip().split('\t')

        # Write the output header
        print(*output_fields, sep='\t')
        
        for line in fin:
            d = {field: value for field, value in zip(header, line.strip().split('\t'))}
            print(*(d.get(field,'') for field in output_fields), sep='\t')

