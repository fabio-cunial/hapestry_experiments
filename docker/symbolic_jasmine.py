"""
Makes all DEL and INV calls symbolic to speed up Jasmine.
"""

import sys
import pysam
import truvari


if __name__ == '__main__':
    vcf = pysam.VariantFile(sys.argv[1])
    n_header = vcf.header.copy()
    out = pysam.VariantFile("/dev/stdout", 'w', header=n_header)
    for entry in vcf:
        if 'SVTYPE' in entry.info:
            if entry.info['SVTYPE'] == 'DEL':
                entry.ref = entry.ref[0]
                entry.alts = ['<DEL>']
            elif entry.info['SVTYPE'] == 'INV':
                entry.ref = entry.ref[0]
                entry.alts = ['<INV>']
        
        # Outputting
        entry.translate(n_header)
        try:
            out.write(entry)
        except Exception:
            sys.stderr.write(f"{entry}\n{type(entry)}\n")
