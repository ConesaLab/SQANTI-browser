"""
Validation track building: STAR SJ, CAGE peaks, PolyA peaks, reference GTF to bigBed.
"""

from __future__ import annotations

import os
import logging
import shutil
import subprocess
from typing import Any, Optional

from src.constants import DEFAULT_VALIDATION_COLORS

logger = logging.getLogger(__name__)


class ValidationTrackBuilder:
    """Builds validation tracks (STAR SJ, CAGE, PolyA, reference) from optional input files."""

    def __init__(self, converter: Any) -> None:
        self.converter = converter

    def _get_validation_track_color(self, track_name: str) -> tuple[int, int, int]:
        """Get RGB tuple for a validation track (CAGE_peaks, polyA_peaks, star_sj, reference)."""
        custom = getattr(self.converter, 'custom_validation_colors', {}) or {}
        return custom.get(track_name, DEFAULT_VALIDATION_COLORS.get(track_name, (128, 128, 128)))

    def _pack_rgb(self, r: int, g: int, b: int) -> int:
        """Pack RGB to BED itemRgb integer."""
        return r * 65536 + g * 256 + b

    def create_star_sj_bigbed(self, star_sj_file: str) -> Optional[Any]:
        """Convert STAR splice junctions to bigBed"""
        logger.info("Converting STAR splice junctions to bigBed...")
        try:
            if self.converter.chrom_sizes_file:
                chrom_sizes = self.converter.chrom_sizes_file
            elif self.converter.two_bit_file:
                chrom_sizes = self.converter.extract_chrom_sizes_from_twobit()
            else:
                chrom_sizes = self.converter.extract_chrom_sizes()

            if not chrom_sizes:
                return None

            r, g, b = self._get_validation_track_color("star_sj")
            item_rgb = self._pack_rgb(r, g, b)

            bed_file = os.path.join(self.converter.temp_dir, "star_junctions.bed")
            with open(star_sj_file, 'r') as infile, open(bed_file, 'w') as outfile:
                for i, line in enumerate(infile):
                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        continue
                    chrom = parts[0]
                    start = int(parts[1]) - 1
                    end = int(parts[2])
                    strand_val = parts[3]
                    strand = '+' if strand_val == '1' else ('-' if strand_val == '2' else '.')
                    unique = int(parts[6])
                    multi = int(parts[7])
                    score = min(unique + multi, 1000)
                    name = f"JUNC{i+1}"
                    outfile.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{start}\t{end}\t{item_rgb}\n")

            sorted_bed = os.path.join(self.converter.temp_dir, "star_junctions.sorted.bed")
            env = os.environ.copy()
            env["LC_COLLATE"] = "C"
            subprocess.run(['sort', '-k1,1', '-k2,2n', bed_file, '-o', sorted_bed], check=True, env=env)

            bb_file = self.converter.output_dir / self.converter.genome / f"{self.converter.genome}_star_sj.bb"
            subprocess.run(['bedToBigBed', '-type=bed9', sorted_bed, chrom_sizes, str(bb_file)], check=True)
            logger.info(f"STAR bigBed created: {bb_file}")
            return bb_file
        except Exception as e:
            logger.error(f"Error converting STAR SJ: {e}")
            return None

    def create_cage_bigbed(self, cage_file: str) -> Optional[Any]:
        """Convert CAGE peaks BED to bigBed"""
        logger.info("Converting CAGE peaks to bigBed...")
        try:
            if self.converter.chrom_sizes_file:
                chrom_sizes = self.converter.chrom_sizes_file
            elif self.converter.two_bit_file:
                chrom_sizes = self.converter.extract_chrom_sizes_from_twobit()
            else:
                chrom_sizes = self.converter.extract_chrom_sizes()

            if not chrom_sizes:
                return None

            bed_file = os.path.join(self.converter.temp_dir, "cage_peaks.bed")
            valid_chroms = {}
            with open(chrom_sizes, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        valid_chroms[parts[0]] = int(parts[1])

            env = os.environ.copy()
            env["LC_COLLATE"] = "C"
            r, g, b = self._get_validation_track_color("CAGE_peaks")
            item_rgb = self._pack_rgb(r, g, b)
            filtered_bed = os.path.join(self.converter.temp_dir, "cage_peaks_filtered.bed")
            filtered_count = 0
            total_count = 0
            skipped_bounds = 0

            with open(cage_file, 'r') as infile, open(filtered_bed, 'w') as outfile:
                for line in infile:
                    total_count += 1
                    parts = line.rstrip('\n').split('\t')
                    if parts and len(parts) >= 3 and parts[0] in valid_chroms:
                        try:
                            chrom = parts[0]
                            start = int(parts[1])
                            end_pos = int(parts[2])
                            name = parts[3] if len(parts) > 3 else "."
                            score = parts[4] if len(parts) > 4 else "0"
                            strand = parts[5] if len(parts) > 5 else "."
                            if strand not in ('+', '-', '.'):
                                strand = "."
                            try:
                                score = int(float(score))
                                score = max(0, min(score, 1000))
                            except (ValueError, TypeError):
                                score = 0
                            if end_pos <= valid_chroms[chrom]:
                                outfile.write(f"{chrom}\t{start}\t{end_pos}\t{name}\t{score}\t{strand}\t{start}\t{end_pos}\t{item_rgb}\n")
                                filtered_count += 1
                            else:
                                skipped_bounds += 1
                        except (ValueError, IndexError):
                            continue

            if filtered_count == 0:
                logger.warning(f"No CAGE peaks kept after filtering. Total: {total_count}. Skipped out of bounds: {skipped_bounds}")
                return None
            if filtered_count < total_count:
                logger.warning(f"Filtered CAGE peaks to match genome: {filtered_count}/{total_count} peaks kept. ({skipped_bounds} out of bounds, {total_count - filtered_count - skipped_bounds} on unknown chromosomes)")

            sorted_bed = os.path.join(self.converter.temp_dir, "cage_peaks.sorted.bed")
            subprocess.run(['sort', '-k1,1', '-k2,2n', filtered_bed, '-o', sorted_bed], check=True, env=env)
            bb_file = self.converter.output_dir / self.converter.genome / f"{self.converter.genome}_cage_peaks.bb"
            subprocess.run(['bedToBigBed', '-type=bed9', '-tab', sorted_bed, chrom_sizes, str(bb_file)], check=True)
            logger.info(f"CAGE bigBed created: {bb_file}")
            return bb_file
        except Exception as e:
            logger.error(f"Error converting CAGE peaks: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return None

    def create_polya_bigbed(self, polya_file: str) -> Optional[Any]:
        """Convert PolyA peaks BED to bigBed"""
        logger.info("Converting PolyA peaks to bigBed...")
        try:
            if self.converter.chrom_sizes_file:
                chrom_sizes = self.converter.chrom_sizes_file
            elif self.converter.two_bit_file:
                chrom_sizes = self.converter.extract_chrom_sizes_from_twobit()
            else:
                chrom_sizes = self.converter.extract_chrom_sizes()
            if not chrom_sizes:
                return None

            valid_chroms = {}
            with open(chrom_sizes, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        valid_chroms[parts[0]] = int(parts[1])
            env = os.environ.copy()
            env["LC_COLLATE"] = "C"
            r, g, b = self._get_validation_track_color("polyA_peaks")
            item_rgb = self._pack_rgb(r, g, b)
            filtered_bed = os.path.join(self.converter.temp_dir, "polya_peaks_filtered.bed")
            filtered_count = 0
            total_count = 0
            skipped_bounds = 0

            with open(polya_file, 'r') as infile, open(filtered_bed, 'w') as outfile:
                for line in infile:
                    total_count += 1
                    parts = line.rstrip('\n').split('\t')
                    if parts and len(parts) >= 3 and parts[0] in valid_chroms:
                        try:
                            chrom = parts[0]
                            start = int(parts[1])
                            end_pos = int(parts[2])
                            name = parts[3] if len(parts) > 3 else "."
                            score_val = parts[4] if len(parts) > 4 else "0"
                            strand = parts[5] if len(parts) > 5 else "."
                            if strand not in ('+', '-', '.'):
                                strand = "."
                            try:
                                score_float = float(score_val)
                                score = min(int(score_float * 1000), 1000)
                                score = max(0, score)
                            except (ValueError, TypeError):
                                score = 0
                            if end_pos <= valid_chroms[chrom]:
                                outfile.write(f"{chrom}\t{start}\t{end_pos}\t{name}\t{score}\t{strand}\t{start}\t{end_pos}\t{item_rgb}\n")
                                filtered_count += 1
                            else:
                                skipped_bounds += 1
                        except (ValueError, IndexError):
                            continue

            if filtered_count == 0:
                logger.warning(f"No PolyA peaks kept after filtering. Total: {total_count}. Skipped out of bounds: {skipped_bounds}")
                return None
            if filtered_count < total_count:
                logger.warning(f"Filtered PolyA peaks to match genome: {filtered_count}/{total_count} peaks kept. ({skipped_bounds} out of bounds, {total_count - filtered_count - skipped_bounds} on unknown chromosomes)")

            sorted_bed = os.path.join(self.converter.temp_dir, "polya_peaks.sorted.bed")
            subprocess.run(['sort', '-k1,1', '-k2,2n', filtered_bed, '-o', sorted_bed], check=True, env=env)
            bb_file = self.converter.output_dir / self.converter.genome / f"{self.converter.genome}_polya_peaks.bb"
            subprocess.run(['bedToBigBed', '-type=bed9', '-tab', sorted_bed, chrom_sizes, str(bb_file)], check=True)
            logger.info(f"PolyA bigBed created: {bb_file}")
            return bb_file
        except Exception as e:
            logger.error(f"Error converting PolyA peaks: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return None

    def create_reference_bigbed(self, ref_gtf_file: str) -> Optional[Any]:
        """Convert reference GTF to bigBed for direct comparison"""
        logger.info("Converting reference GTF to bigBed...")
        try:
            if self.converter.chrom_sizes_file:
                chrom_sizes = self.converter.chrom_sizes_file
            elif self.converter.two_bit_file:
                chrom_sizes = self.converter.extract_chrom_sizes_from_twobit()
            else:
                chrom_sizes = self.converter.extract_chrom_sizes()
            if not chrom_sizes:
                return None

            ref_genepred = os.path.join(self.converter.temp_dir, "reference.genepred")
            logger.info("Converting reference GTF to GenePred...")
            subprocess.run(['gtfToGenePred', ref_gtf_file, ref_genepred], check=True)
            ref_bed = os.path.join(self.converter.temp_dir, "reference.bed")
            logger.info("Converting reference GenePred to BED...")
            subprocess.run(['genePredToBed', ref_genepred, ref_bed], check=True)

            valid_chroms = {}
            with open(chrom_sizes, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        valid_chroms[parts[0]] = int(parts[1])
            filtered_bed = os.path.join(self.converter.temp_dir, "reference_filtered.bed")
            filtered_count = 0
            total_count = 0
            skipped_bounds = 0
            skipped_chroms = 0

            with open(ref_bed, 'r') as infile, open(filtered_bed, 'w') as outfile:
                for line in infile:
                    total_count += 1
                    parts = line.split('\t')
                    if parts and len(parts) >= 3:
                        try:
                            chrom = parts[0]
                            start_pos = int(parts[1])
                            end_pos = int(parts[2])
                            if chrom not in valid_chroms:
                                skipped_chroms += 1
                                continue
                            if end_pos <= valid_chroms[chrom] and start_pos >= 0:
                                outfile.write(line)
                                filtered_count += 1
                            else:
                                skipped_bounds += 1
                        except (ValueError, IndexError):
                            continue

            if filtered_count == 0:
                logger.warning(f"No reference transcripts kept after filtering. Total: {total_count}. Skipped: {skipped_bounds} out of bounds, {skipped_chroms} on unknown chromosomes")
                return None
            if filtered_count < total_count:
                logger.warning(f"Filtered reference transcripts to match genome: {filtered_count}/{total_count} transcripts kept. ({skipped_bounds} out of bounds, {skipped_chroms} on unknown chromosomes)")

            env = os.environ.copy()
            env["LC_COLLATE"] = "C"
            ref_sorted_bed = os.path.join(self.converter.temp_dir, "reference.sorted.bed")
            subprocess.run(['sort', '-k1,1', '-k2,2n', filtered_bed, '-o', ref_sorted_bed], check=True, env=env)
            bb_file = self.converter.output_dir / self.converter.genome / f"{self.converter.genome}_reference.bb"
            subprocess.run(['bedToBigBed', '-type=bed12', '-tab', ref_sorted_bed, chrom_sizes, str(bb_file)], check=True)
            logger.info(f"Reference bigBed created: {bb_file}")
            return bb_file
        except Exception as e:
            logger.error(f"Error converting reference GTF: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return None
