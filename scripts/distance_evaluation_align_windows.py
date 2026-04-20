#!/usr/bin/env python3
"""
Compute global alignments between sequence pairs in genomic windows using pywfa.

The script processes paired reference and query directories, each containing
subdirectories representing genomic windows. For each window pair, it computes
global alignment scores for all sequence pairings and selects the best scoring
assignment (bipartite matching of sequences).

Usage:
    python align_windows_pywfa.py <reference_dir> <query_dir> <output_csv>

This script was mostly written by Copilot with auto model.
"""

import argparse
import csv
import os
import sys
import re
from pathlib import Path
from typing import List, Tuple, Dict
import logging

try:
    from pywfa import WavefrontAligner
except ImportError:
    print("Error: pywfa is not installed. Install it with: pip install pywfa")
    sys.exit(1)


class ReferenceSequenceCountError(ValueError):
    """Raised when the reference FASTA has fewer than two sequences."""


class QuerySequenceCountError(ValueError):
    """Raised when the query FASTA has fewer than two sequences."""


def setup_logging(level=logging.INFO):
    """Configure logging."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)


def parse_bed_file(bed_path: str) -> Dict[str, List[Tuple[int, int]]]:
    """
    Parse a BED file and return confident regions organized by chromosome.
    
    BED format (0-based, half-open intervals [start, end)):
    - Column 1: chromosome
    - Column 2: start position
    - Column 3: end position
    - Columns 4+: optional fields (ignored)
    
    Args:
        bed_path: Path to BED file
        
    Returns:
        Dict mapping chromosome names to sorted lists of (start, end) tuples
    """
    confident_regions = {}
    
    try:
        with open(bed_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                fields = line.split('\t')
                if len(fields) < 3:
                    logging.warning(f"Skipping invalid BED line {line_num}: {line}")
                    continue
                
                try:
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    
                    if chrom not in confident_regions:
                        confident_regions[chrom] = []
                    confident_regions[chrom].append((start, end))
                except ValueError as e:
                    logging.warning(f"Skipping BED line {line_num} with invalid coordinates: {line}")
                    continue
        
        # Sort regions for each chromosome
        for chrom in confident_regions:
            confident_regions[chrom].sort()
        
        logging.info(f"Loaded confident regions from {bed_path}")
        for chrom, regions in confident_regions.items():
            logging.debug(f"  {chrom}: {len(regions)} confident regions")
        
        return confident_regions
    
    except IOError as e:
        logging.error(f"Error reading BED file {bed_path}: {e}")
        raise


def parse_window_name(window_name: str) -> Tuple[str, int, int]:
    """
    Parse a window name in the format 'chr1_212185028-212185584'.
    
    Args:
        window_name: Window name string
        
    Returns:
        Tuple of (chromosome, start, end)
    """
    # Match pattern like "chr1_212185028-212185584" or "chrX_1000-2000"
    match = re.match(r'^(.+?)_(\d+)-(\d+)$', window_name)
    if not match:
        raise ValueError(f"Invalid window name format: {window_name}")
    
    chrom, start_str, end_str = match.groups()
    return chrom, int(start_str), int(end_str)


def is_window_in_confident_region(window_name: str, confident_regions: Dict[str, List[Tuple[int, int]]]) -> bool:
    """
    Check if a window is completely contained in at least one confident region.
    
    Args:
        window_name: Window name (e.g., 'chr1_212185028-212185584')
        confident_regions: Dict of confident regions by chromosome
        
    Returns:
        True if window is contained in at least one confident region, False otherwise
    """
    try:
        chrom, window_start, window_end = parse_window_name(window_name)
    except ValueError as e:
        logging.warning(f"Could not parse window name {window_name}: {e}")
        return False
    
    if chrom not in confident_regions:
        return False
    
    # Check if window [window_start, window_end) is contained in any confident region
    for region_start, region_end in confident_regions[chrom]:
        if region_start <= window_start and window_end <= region_end:
            return True
    
    return False


def read_fasta(fasta_path: str) -> List[Tuple[str, str]]:
    """
    Parse a FASTA file and return a list of (header, sequence) tuples.
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        List of (header, sequence) tuples
    """
    sequences = []
    current_header = None
    current_seq = []
    
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                if not line:
                    continue
                if line.startswith('>'):
                    if current_header is not None:
                        sequences.append((current_header, ''.join(current_seq)))
                    current_header = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Don't forget the last sequence
            if current_header is not None:
                sequences.append((current_header, ''.join(current_seq)))
    
    except IOError as e:
        logging.error(f"Error reading FASTA file {fasta_path}: {e}")
        raise
    
    return sequences


def select_longest_by_haplotype(records: List[Tuple[str, str]], fasta_path: str) -> List[Tuple[str, str]]:
    """
    Select one sequence for haplotype 1 and one for haplotype 2, keeping the
    longest sequence for each haplotype.

    Expected header format: "chrZ_S_X Y", where X is 1 or 2.
    Haplotype is inferred from the first whitespace-delimited token, using the
    suffix after the final underscore.
    """
    best_by_haplotype = {}

    for header, sequence in records:
        token = header.split()[0]
        haplotype = token.rsplit('_', 1)[-1] if '_' in token else None
        if haplotype not in {"1", "2"}:
            continue

        current = best_by_haplotype.get(haplotype)
        if current is None or len(sequence) > len(current[1]):
            best_by_haplotype[haplotype] = (header, sequence)

    if "1" not in best_by_haplotype or "2" not in best_by_haplotype:
        raise ValueError(
            f"Expected at least one sequence with haplotype 1 and one with haplotype 2 in {fasta_path}"
        )

    return [best_by_haplotype["1"], best_by_haplotype["2"]]


def cigar_to_edit_distance(cigar_string: str) -> int:
    """
    Compute edit distance from a CIGAR string.
    Edit distance counts mismatches (X), insertions (I), and deletions (D).
    """
    edit_distance = 0
    for m in re.finditer(r'(\d+)([MIDNSHP=X])', cigar_string):
        count = int(m.group(1))
        op = m.group(2)
        if op in ('X', 'I', 'D'):
            edit_distance += count
    return edit_distance


def cigar_to_alignment_length(cigar_string: str) -> int:
    """
    Compute alignment length from a CIGAR string as the total number of
    operations, including matches.
    """
    alignment_length = 0
    for m in re.finditer(r'(\d+)([MIDNSHP=X])', cigar_string):
        alignment_length += int(m.group(1))
    return alignment_length


def compute_normalized_edit_distance(seq1: str, seq2: str, aligner: WavefrontAligner) -> float:
    """
    Compute edit distance between two sequences via WFA CIGAR string,
    normalized by the total alignment length in the CIGAR string.

    Args:
        seq1: Reference sequence
        seq2: Query sequence
        aligner: WavefrontAligner instance configured for full alignment

    Returns:
        Edit distance (mismatches + insertions + deletions) divided by the
        total number of CIGAR operations, including matches
    """
    try:
        result = aligner(seq1, seq2)
        edit_distance = cigar_to_edit_distance(result.cigarstring)
        alignment_length = cigar_to_alignment_length(result.cigarstring)
        if alignment_length == 0:
            raise ValueError(f"Alignment length is zero for CIGAR: {result.cigarstring}")
        if edit_distance > alignment_length:
            print(f"Edit distance {edit_distance} exceeds alignment length {alignment_length} for CIGAR: {result.cigarstring}")
            exit(1)
        return edit_distance / alignment_length
    except Exception as e:
        logging.error(f"Error computing alignment: {e}")
        exit(1)


def find_best_pairing(normalized_distances: Dict[Tuple[int, int], float]) -> Tuple[float, Tuple[int, int], Tuple[int, int]]:
    """
    Find the best pairing of 2 reference sequences with 2 query sequences.
    
    For a 2x2 bipartite matching, there are only 2 possible pairings:
    - (ref[0], query[0]) and (ref[1], query[1])
    - (ref[0], query[1]) and (ref[1], query[0])
    
    Args:
        normalized_distances: Dict mapping (ref_idx, query_idx) to normalized edit distance
        
    Returns:
        Tuple of (normalized_distance, ref_assignment, query_assignment)
        where ref_assignment[i] = j means ref sequence i pairs with query sequence j
    """
    # Pairing option 1: ref[0] -> query[0], ref[1] -> query[1]
    distance1 = normalized_distances.get((0, 0), 0) + normalized_distances.get((1, 1), 0)
    
    # Pairing option 2: ref[0] -> query[1], ref[1] -> query[0]
    distance2 = normalized_distances.get((0, 1), 0) + normalized_distances.get((1, 0), 0)
    
    if distance1 < distance2:
        return distance1, (0, 0), (1, 1)  # ref[0]->query[0], ref[1]->query[1]
    else:
        return distance2, (0, 1), (1, 0)  # ref[0]->query[1], ref[1]->query[0]


def process_window(window_name: str, ref_window_dir: str, query_window_dir: str) -> Tuple[str, int, float]:
    """
    Process a single genomic window pair and compute the min pairing score.
    
    Args:
        window_name: Name of the genomic window
        ref_window_dir: Path to reference window directory
        query_window_dir: Path to query window directory
        
    Returns:
        Tuple of (window_name, window_length, min_normalized_distance)
    """
    ref_fasta = os.path.join(ref_window_dir, 'sequences.fasta')
    query_fasta = os.path.join(query_window_dir, 'sequences.fasta')
    
    if not os.path.exists(ref_fasta):
        raise FileNotFoundError(f"Reference FASTA not found: {ref_fasta}")
    if not os.path.exists(query_fasta):
        raise FileNotFoundError(f"Query FASTA not found: {query_fasta}")
    
    # Read sequences
    ref_sequences = read_fasta(ref_fasta)
    query_sequences = read_fasta(query_fasta)
    
    if len(ref_sequences) > 2:
        ref_sequences = select_longest_by_haplotype(ref_sequences, ref_fasta)
        logging.info(
            "  Reference FASTA has >2 sequences; selected a longest hap 1 and a longest hap2 sequence"
        )
    elif len(ref_sequences) < 2:
        raise ReferenceSequenceCountError(
            f"Expected exactly 2 sequences in {ref_fasta}, found {len(ref_sequences)}"
        )
    if len(query_sequences) > 2:
        query_sequences = select_longest_by_haplotype(query_sequences, query_fasta)
        logging.info(
            "  Query FASTA has >2 sequences; selected a longest hap 1 and a longest hap2 sequence"
        )    
    elif len(query_sequences) < 2:
        raise QuerySequenceCountError(
            f"Expected exactly 2 sequences in {query_fasta}, found {len(query_sequences)}"
        )
    
    # Extract sequence strings
    ref_seqs = [seq for _, seq in ref_sequences]
    query_seqs = [seq for _, seq in query_sequences]
    
    logging.info(f"Processing window {window_name}")
    logging.debug(f"  Reference sequences: {len(ref_seqs[0])} bp, {len(ref_seqs[1])} bp")
    logging.debug(f"  Query sequences: {len(query_seqs[0])} bp, {len(query_seqs[1])} bp")
    
    # Compute all pairwise edit distances
    _, window_start, window_end = parse_window_name(window_name)
    window_length = window_end - window_start

    aligner = WavefrontAligner(span="end-to-end")
    normalized_distances = {}
    for i, ref_seq in enumerate(ref_seqs):
        for j, query_seq in enumerate(query_seqs):
            distance = compute_normalized_edit_distance(ref_seq, query_seq, aligner)
            normalized_distances[(i, j)] = distance
            logging.debug(f"  Normalized edit distance ref[{i}] vs query[{j}]: {distance:.6f}")

    # Find best pairing
    best_normalized_distance, _, _ = find_best_pairing(normalized_distances)
    logging.info(f"  Best pairing normalized edit distance: {best_normalized_distance:.6f}")

    return window_name, window_length, best_normalized_distance


def main():
    parser = argparse.ArgumentParser(
        description='Compute global alignments between sequence pairs in genomic windows using WFA',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python align_windows_pywfa.py ./reference_dir ./query_dir results.csv
  python align_windows_pywfa.py ./reference_dir ./query_dir results.csv --verbose
        """
    )
    
    parser.add_argument(
        'reference_dir',
        help='Path to reference directory containing genomic window subdirectories'
    )
    parser.add_argument(
        'query_dir',
        help='Path to query directory containing genomic window subdirectories'
    )
    parser.add_argument(
        'output_csv',
        help='Path to output CSV file'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    parser.add_argument(
        '--confident-bed',
        type=str,
        default=None,
        help='Path to BED file with confident regions. Only windows completely contained in these regions will be processed.'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = setup_logging(log_level)
    
    # Validate input directories
    ref_dir = Path(args.reference_dir)
    query_dir = Path(args.query_dir)
    
    if not ref_dir.is_dir():
        logger.error(f"Reference directory does not exist: {args.reference_dir}")
        sys.exit(1)
    
    if not query_dir.is_dir():
        logger.error(f"Query directory does not exist: {args.query_dir}")
        sys.exit(1)
    
    logger.info(f"Reference directory: {ref_dir}")
    logger.info(f"Query directory: {query_dir}")
    logger.info(f"Output CSV: {args.output_csv}")
    
    # Load confident regions if BED file is provided
    confident_regions = {}
    if args.confident_bed:
        if not os.path.exists(args.confident_bed):
            logger.error(f"Confident BED file does not exist: {args.confident_bed}")
            sys.exit(1)
        confident_regions = parse_bed_file(args.confident_bed)
        logger.info(f"Loaded confident regions from {args.confident_bed}")
    else:
        logger.info("No confident BED file specified; processing all windows")
    
    # Get list of window directories from reference
    ref_windows = sorted([d.name for d in ref_dir.iterdir() if d.is_dir()])
    logger.info(f"Found {len(ref_windows)} windows in reference directory")
    
    # Filter windows by confident regions if BED file is provided
    windows_rejected = 0
    if confident_regions:
        ref_windows_filtered = []
        for window_name in ref_windows:
            if is_window_in_confident_region(window_name, confident_regions):
                ref_windows_filtered.append(window_name)
            else:
                windows_rejected += 1
                logger.debug(f"Window {window_name} not in confident region; skipping")
        ref_windows = ref_windows_filtered
        logger.info(f"After filtering by confident regions: {len(ref_windows)} windows to process ({windows_rejected} rejected)")
    
    # Process each window
    results = []
    windows_processed = 0
    windows_failed = 0
    windows_failed_not_found = 0
    windows_failed_ref_lt2 = 0
    windows_failed_query_lt2 = 0
    
    for window_name in ref_windows:
        ref_window_path = ref_dir / window_name
        query_window_path = query_dir / window_name
        
        if not query_window_path.exists():
            logger.warning(f"Corresponding query window not found: {window_name}")
            windows_failed += 1
            windows_failed_not_found += 1
            continue
        
        try:
            window_name_result, window_length, score = process_window(
                window_name,
                str(ref_window_path),
                str(query_window_path)
            )
            results.append((window_name_result, window_length, score))
            windows_processed += 1
        except ReferenceSequenceCountError as e:
            logger.warning(f"Skipping window {window_name}: {e}")
            windows_failed += 1
            windows_failed_ref_lt2 += 1
            continue
        except QuerySequenceCountError as e:
            logger.warning(f"Skipping window {window_name}: {e}")
            windows_failed += 1
            windows_failed_query_lt2 += 1
            continue
        except Exception as e:
            logger.error(f"Error processing window {window_name}: {e}")
            windows_failed += 1
            continue
    
    # Write results to CSV
    try:
        output_path = Path(args.output_csv)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, lineterminator='\n')
            writer.writerow(['window', 'window_length', 'normalized_edit_distance'])
            writer.writerows(results)
        
        logger.info(f"Results written to {args.output_csv}")
    except IOError as e:
        logger.error(f"Error writing output CSV: {e}")
        sys.exit(1)
    
    # Summary
    logger.info(f"\n=== Summary ===")
    if confident_regions:
        logger.info(f"Windows rejected (not in confident region): {windows_rejected}")
    logger.info(f"Windows processed: {windows_processed}")
    logger.info(f"Windows failed (total): {windows_failed}")
    logger.info(f"Windows failed (ref window not found in query): {windows_failed_not_found}")
    logger.info(f"Windows failed (ref has <2 sequences): {windows_failed_ref_lt2}")
    logger.info(f"Windows failed (query has <2 sequences): {windows_failed_query_lt2}")
    logger.info(f"Total results: {len(results)}")
    
    if windows_failed > 0:
        sys.exit(1)


if __name__ == '__main__':
    main()
