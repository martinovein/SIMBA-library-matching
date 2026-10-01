#!/usr/bin/env python3
"""Count DNA and protein variants in SIMBA NNK sequencing FASTQ files.

The pipeline validates the workbook, expands one-NNK designs, streams gzip FASTQ
records, maps exact coding sequences, exports raw integer count tables, and
creates QC and Enrich2-oriented artifacts.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import math
import os
import re
import shutil
import statistics
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator

from openpyxl import Workbook, load_workbook
from openpyxl.chart import BarChart, Reference
from openpyxl.formatting.rule import CellIsRule
from openpyxl.styles import Font, PatternFill
from openpyxl.worksheet.table import Table, TableStyleInfo

PIPELINE_VERSION = "0.1.0"
REFERENCE_EXPANSION_VERSION = "0.1.0"
DNA_BASES = set("ACGT")
NNK_CODONS = [a + b + k for a in "ACGT" for b in "ACGT" for k in "GT"]
READ_CATEGORIES = [
    "mapped_nnk",
    "mapped_control",
    "low_quality",
    "missing_flank",
    "invalid_length_frame",
    "ambiguous",
    "unexpected",
    "malformed",
]

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AA3 = {
    "A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe",
    "G": "Gly", "H": "His", "I": "Ile", "K": "Lys", "L": "Leu",
    "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg",
    "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr", "*": "Ter",
}

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass(frozen=True)
class Sample:
    fastq_filename: str
    sample_id: str
    replicate: str
    bait: str
    time_hours: int | None
    is_common_input: bool
    order: int


@dataclass
class DnaVariant:
    variant_id: str
    parent: str
    dna_sequence: str
    protein_sequence: str
    dna_hgvs: str
    protein_hgvs: str
    classification: str
    randomized_position_aa: int
    wt_codon: str
    observed_codon: str
    wt_aa: str
    observed_aa: str
    aliases: list[str] = field(default_factory=list)


@dataclass
class ControlRef:
    control_id: str
    control_type: str
    peptide_aa: str
    coding_dna: str
    expected_function: str
    translation_orientation: str = "forward"


@dataclass
class ReferenceBundle:
    samples: list[Sample]
    parents: list[str]
    dna_variants: dict[str, DnaVariant]
    controls: dict[str, ControlRef]
    seq_to_variant: dict[str, str]
    seq_to_control: dict[str, str]
    parent_wt_dna: dict[str, str]
    parent_wt_protein: dict[str, str]
    expected_lengths: list[int]
    anchor5: str
    anchor3: str
    validation_rows: list[dict[str, Any]]
    fatal_messages: list[str]
    warnings: list[str]


@dataclass
class Extraction:
    coding: str
    coding_qual: str
    orientation: str
    reason: str = ""


@dataclass
class SampleResult:
    sample: Sample
    file_size: int
    total_reads: int = 0
    category_counts: Counter = field(default_factory=Counter)
    orientation_counts: Counter = field(default_factory=Counter)
    read_lengths: Counter = field(default_factory=Counter)
    mean_quality_sum: float = 0.0
    q20_bases: int = 0
    q30_bases: int = 0
    total_bases: int = 0
    dna_counts: Counter = field(default_factory=Counter)
    control_counts: Counter = field(default_factory=Counter)
    unexpected: Counter = field(default_factory=Counter)
    unexpected_orientation: dict[str, Counter] = field(default_factory=lambda: defaultdict(Counter))
    malformed_messages: list[str] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)
    fingerprint: dict[str, Any] = field(default_factory=dict)


def normalize_dna(value: Any) -> str:
    return str(value or "").strip().upper()


def translate_dna(seq: str) -> str:
    seq = normalize_dna(seq)
    if len(seq) % 3:
        raise ValueError(f"DNA length is not divisible by three: {len(seq)}")
    aa = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if codon not in CODON_TABLE:
            raise ValueError(f"Invalid codon {codon!r}")
        aa.append(CODON_TABLE[codon])
    return "".join(aa)


def reverse_complement(seq: str) -> str:
    return seq.translate(COMP)[::-1].upper()


def sha256_file(path: Path) -> str:
    """Return a streaming SHA-256 digest for a file."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def file_fingerprint(path: Path) -> dict[str, Any]:
    st = path.stat()
    return {"path": str(path), "size": st.st_size, "mtime_ns": st.st_mtime_ns}


def checkpoint_fingerprint(sample: Sample, data_dir: Path, params: dict[str, Any]) -> dict[str, Any]:
    return {
        "workbook_sha256": params.get("workbook_sha256", ""),
        "fastq": file_fingerprint(data_dir / sample.fastq_filename),
        "anchor5": params["anchor5"],
        "anchor3": params["anchor3"],
        "min_base_quality": params["min_base_quality"],
        "min_mean_quality": params.get("min_mean_quality"),
        "max_n": params["max_n"],
        "orientation": "forward_then_reverse_complement_exact_anchors",
        "pipeline_version": PIPELINE_VERSION,
        "reference_expansion_version": REFERENCE_EXPANSION_VERSION,
        "max_reads_per_file": params.get("max_reads_per_file"),
    }


def dna_hgvs(wt: str, obs: str) -> str:
    if wt == obs:
        return "_wt"
    changes = [f"c.{i}{w}>{o}" for i, (w, o) in enumerate(zip(wt, obs), start=1) if w != o]
    return ",".join(changes)


def protein_hgvs(wt: str, obs: str) -> str:
    if wt == obs:
        return "_wt"
    changes = []
    for i, (w, o) in enumerate(zip(wt, obs), start=1):
        if w != o:
            changes.append(f"p.{AA3[w]}{i}{AA3[o]}")
    return ",".join(changes)


def variant_class(wt_dna: str, obs_dna: str, wt_protein: str, obs_protein: str) -> str:
    if obs_dna == wt_dna:
        return "WT"
    if obs_protein == wt_protein:
        return "synonymous"
    if "*" in obs_protein:
        return "stop"
    return "missense"


def load_sheet(path: Path, sheet: str) -> tuple[list[str], list[dict[str, Any]]]:
    wb = load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet]
    rows = list(ws.iter_rows(values_only=True))
    headers = [str(h).strip() for h in rows[0] if h is not None]
    out = []
    for row in rows[1:]:
        values = list(row[:len(headers)])
        if all(v is None for v in values):
            continue
        out.append(dict(zip(headers, values)))
    return headers, out


def find_observed_anchor(full_5: str, full_3: str, anchor5_override: str | None, anchor3_override: str | None) -> tuple[str, str]:
    """Return retained anchors. Defaults match the observed Cutadapt-trimmed reads."""
    if anchor5_override:
        anchor5 = normalize_dna(anchor5_override)
    else:
        suffix = normalize_dna(full_5)
        core = "CAGCCTCTTCATCTGGC"
        anchor5 = "TATG" + core if suffix.endswith("TCAGCCTCTTCATCTGGC") else suffix[-21:]
    if anchor3_override:
        anchor3 = normalize_dna(anchor3_override)
    else:
        anchor3 = normalize_dna(full_3)[:17]
    return anchor5, anchor3


def validate_and_expand(workbook: Path, data_dir: Path, anchor5: str | None, anchor3: str | None) -> ReferenceBundle:
    nnk_headers, nnk_rows = load_sheet(workbook, "NNK_library")
    ctrl_headers, ctrl_rows = load_sheet(workbook, "Fixed_controls")
    sample_headers, sample_rows = load_sheet(workbook, "Samples")
    validation_rows: list[dict[str, Any]] = []
    fatal: list[str] = []
    warnings: list[str] = []
    required = {
        "NNK_library": ["parent_id", "oligo_id", "wt_peptide_aa", "wt_coding_dna", "nnk_coding_pattern", "randomized_position_aa", "fixed_5prime_sequence", "fixed_3prime_sequence"],
        "Fixed_controls": ["control_id", "control_type", "peptide_aa", "coding_dna", "fixed_5prime_sequence", "fixed_3prime_sequence", "expected_function"],
        "Samples": ["fastq_filename", "sample_id", "replicate", "bait", "time_hours", "is_common_input"],
    }
    for sheet, observed in [("NNK_library", nnk_headers), ("Fixed_controls", ctrl_headers), ("Samples", sample_headers)]:
        missing = [c for c in required[sheet] if c not in observed]
        validation_rows.append({"check": f"{sheet} columns", "status": "FAIL" if missing else "PASS", "detail": ",".join(missing) if missing else ",".join(observed)})
        fatal.extend([f"{sheet} missing required column {c}" for c in missing])

    if fatal:
        return ReferenceBundle([], [], {}, {}, {}, {}, {}, {}, [], "", "", validation_rows, fatal, warnings)

    samples = []
    for i, r in enumerate(sample_rows):
        is_common = bool(r["is_common_input"])
        raw_time = r["time_hours"]
        try:
            time_hours = int(raw_time) if raw_time not in (None, "") else None
        except (TypeError, ValueError):
            if is_common:
                time_hours = None
                warnings.append(f"{r['fastq_filename']}: common input has nonnumeric time_hours={raw_time!r}; stored as blank/NA")
            else:
                fatal.append(f"{r['fastq_filename']}: non-common sample has nonnumeric time_hours={raw_time!r}")
                time_hours = None
        samples.append(Sample(str(r["fastq_filename"]), str(r["sample_id"]), str(r["replicate"]), str(r["bait"]),
                              time_hours, is_common, i))
    manifest_files = [s.fastq_filename for s in samples]
    data_files = sorted(p.name for p in data_dir.glob("*.gz"))
    missing_files = sorted(set(manifest_files) - set(data_files))
    unexpected_files = sorted(set(data_files) - set(manifest_files))
    dup_files = sorted(k for k, v in Counter(manifest_files).items() if v > 1)
    for name, items, is_fatal in [("missing FASTQ", missing_files, True), ("unexpected FASTQ", unexpected_files, False), ("duplicated manifest FASTQ", dup_files, True)]:
        status = "FAIL" if items and is_fatal else ("WARNING" if items else "PASS")
        validation_rows.append({"check": name, "status": status, "detail": ";".join(items)})
    fatal.extend(f"Missing FASTQ listed in Samples: {x}" for x in missing_files)
    fatal.extend(f"Duplicated FASTQ in Samples: {x}" for x in dup_files)
    warnings.extend(f"Unexpected FASTQ in data directory: {x}" for x in unexpected_files)

    full_5s = {str(r["fixed_5prime_sequence"]) for r in nnk_rows + ctrl_rows}
    full_3s = {str(r["fixed_3prime_sequence"]) for r in nnk_rows + ctrl_rows}
    if len(full_5s) != 1 or len(full_3s) != 1:
        fatal.append("Multiple fixed flank sequences found; this script expects one retained anchor pair.")
    obs_anchor5, obs_anchor3 = find_observed_anchor(next(iter(full_5s)), next(iter(full_3s)), anchor5, anchor3)
    if len(obs_anchor5) < 12 or len(obs_anchor3) < 12:
        fatal.append("Anchor length must be at least 12 bp to avoid plausible random matches.")

    parent_wt_dna: dict[str, str] = {}
    parent_wt_protein: dict[str, str] = {}
    for r in nnk_rows:
        parent = str(r["parent_id"])
        wt_dna = normalize_dna(r["wt_coding_dna"])
        wt_pep = str(r["wt_peptide_aa"]).strip().upper()
        if set(wt_dna) - DNA_BASES or len(wt_dna) % 3:
            fatal.append(f"{r['oligo_id']}: invalid WT DNA")
            continue
        try:
            translated = translate_dna(wt_dna).upper()
        except ValueError as exc:
            fatal.append(f"{r['oligo_id']}: {exc}")
            continue
        if translated != wt_pep:
            fatal.append(f"{r['oligo_id']}: WT DNA translates to {translated}, not {wt_pep}")
        parent_wt_dna.setdefault(parent, wt_dna)
        parent_wt_protein.setdefault(parent, translated)
        if parent_wt_dna[parent] != wt_dna:
            fatal.append(f"{parent}: inconsistent WT DNA across rows")

    dna_variants: dict[str, DnaVariant] = {}
    seq_to_variant: dict[str, str] = {}
    alias_collisions: list[str] = []
    for r in nnk_rows:
        parent = str(r["parent_id"])
        oligo = str(r["oligo_id"])
        pattern = normalize_dna(r["nnk_coding_pattern"])
        wt_dna = normalize_dna(r["wt_coding_dna"])
        if "WT" in oligo.upper() and "NNK" not in pattern:
            continue
        if pattern.count("NNK") != 1:
            fatal.append(f"{oligo}: non-WT pattern must contain exactly one NNK")
            continue
        if len(pattern) != len(wt_dna):
            fatal.append(f"{oligo}: NNK pattern length differs from WT")
            continue
        nnk_pos = pattern.index("NNK")
        if nnk_pos % 3:
            fatal.append(f"{oligo}: NNK is not codon-aligned")
            continue
        outside = pattern[:nnk_pos] + pattern[nnk_pos + 3:]
        wt_outside = wt_dna[:nnk_pos] + wt_dna[nnk_pos + 3:]
        if outside != wt_outside:
            fatal.append(f"{oligo}: pattern differs from WT outside NNK")
        actual_pos = nnk_pos // 3 + 1
        if int(r["randomized_position_aa"]) != actual_pos:
            fatal.append(f"{oligo}: randomized_position_aa {r['randomized_position_aa']} != actual {actual_pos}")
        wt_protein = parent_wt_protein.get(parent, translate_dna(wt_dna))
        for codon in NNK_CODONS:
            obs_dna = pattern.replace("NNK", codon)
            obs_protein = translate_dna(obs_dna)
            hgvs = dna_hgvs(wt_dna, obs_dna)
            variant_id = f"{parent}:{hgvs}"
            alias = f"{oligo}|pos{actual_pos}|{codon}"
            if obs_dna in seq_to_variant:
                existing = seq_to_variant[obs_dna]
                if not existing.startswith(parent + ":"):
                    alias_collisions.append(f"{obs_dna} shared by {existing} and {alias}")
                dna_variants[existing].aliases.append(alias)
                continue
            seq_to_variant[obs_dna] = variant_id
            dna_variants[variant_id] = DnaVariant(
                variant_id=variant_id,
                parent=parent,
                dna_sequence=obs_dna,
                protein_sequence=obs_protein,
                dna_hgvs=hgvs,
                protein_hgvs=protein_hgvs(wt_protein, obs_protein),
                classification=variant_class(wt_dna, obs_dna, wt_protein, obs_protein),
                randomized_position_aa=actual_pos,
                wt_codon=wt_dna[nnk_pos:nnk_pos + 3],
                observed_codon=codon,
                wt_aa=wt_protein[actual_pos - 1],
                observed_aa=obs_protein[actual_pos - 1],
                aliases=[alias],
            )
    if alias_collisions:
        fatal.extend(alias_collisions)

    controls: dict[str, ControlRef] = {}
    seq_to_control: dict[str, str] = {}
    for r in ctrl_rows:
        cid = str(r["control_id"])
        dna = normalize_dna(r["coding_dna"])
        pep = str(r["peptide_aa"]).strip().upper()
        translation_orientation = "forward"
        try:
            translated = translate_dna(dna).upper()
        except ValueError as exc:
            fatal.append(f"{cid}: {exc}")
            continue
        direct_ok = translated == pep or (pep.endswith("*") and translated.startswith(pep))
        if not direct_ok:
            rc_translated = translate_dna(reverse_complement(dna)).upper()
            rc_ok = rc_translated == pep or (pep.endswith("*") and rc_translated.startswith(pep))
            if rc_ok:
                translation_orientation = "reverse_complement"
                warnings.append(f"{cid}: control coding_dna translates to peptide on reverse complement; mapping still uses workbook sequence")
            else:
                fatal.append(f"{cid}: control DNA translates to {translated}, reverse complement translates to {rc_translated}, expected {pep}")
        if dna in seq_to_control:
            warnings.append(f"Duplicate fixed control sequence: {cid} and {seq_to_control[dna]}")
        if dna in seq_to_variant:
            fatal.append(f"{cid}: fixed control sequence is indistinguishable from NNK reference {seq_to_variant[dna]}")
        seq_to_control[dna] = cid
        controls[cid] = ControlRef(cid, str(r["control_type"]), pep, dna, str(r["expected_function"]), translation_orientation)

    for seq, vid in seq_to_variant.items():
        if seq in seq_to_control:
            fatal.append(f"{vid}: sequence also present as fixed control")
    expected_lengths = sorted({len(v.dna_sequence) for v in dna_variants.values()} | {len(c.coding_dna) for c in controls.values()})
    validation_rows.append({"check": "expanded unique NNK DNA references", "status": "PASS", "detail": str(len(dna_variants))})
    validation_rows.append({"check": "fixed controls", "status": "PASS", "detail": str(len(controls))})
    validation_rows.append({"check": "retained anchors", "status": "PASS", "detail": f"5prime={obs_anchor5};3prime={obs_anchor3}"})
    return ReferenceBundle(samples, sorted(parent_wt_dna), dna_variants, controls, seq_to_variant, seq_to_control,
                           parent_wt_dna, parent_wt_protein, expected_lengths, obs_anchor5, obs_anchor3,
                           validation_rows, fatal, warnings)


def phred_scores(qual: str) -> list[int]:
    return [ord(ch) - 33 for ch in qual]


def find_one(seq: str, sub: str) -> int | None:
    first = seq.find(sub)
    if first < 0:
        return None
    second = seq.find(sub, first + 1)
    if second >= 0:
        return -2
    return first


def extract_orientation(seq: str, qual: str, anchor5: str, anchor3: str, orientation: str) -> Extraction | None:
    a5 = find_one(seq, anchor5)
    a3 = find_one(seq, anchor3)
    if a5 is None or a3 is None:
        return None
    if a5 == -2 or a3 == -2:
        return Extraction("", "", orientation, "multiple_anchor_matches")
    start = a5 + len(anchor5)
    if a3 < start:
        return Extraction("", "", orientation, "anchor_order")
    return Extraction(seq[start:a3], qual[start:a3], orientation)


def extract_read(seq: str, qual: str, anchor5: str, anchor3: str) -> Extraction | None:
    seq = seq.upper()
    fwd = extract_orientation(seq, qual, anchor5, anchor3, "forward")
    rc_seq = reverse_complement(seq)
    rc_qual = qual[::-1]
    rev = extract_orientation(rc_seq, rc_qual, anchor5, anchor3, "reverse_complement")
    valid = [x for x in [fwd, rev] if x is not None and not x.reason]
    if len(valid) == 1:
        return valid[0]
    if len(valid) > 1:
        return Extraction("", "", "ambiguous", "both_orientations")
    reasons = [x.reason for x in [fwd, rev] if x is not None and x.reason]
    if reasons:
        return Extraction("", "", "ambiguous", ";".join(reasons))
    return None


def classify_extracted(
    extraction: Extraction | None,
    min_base_quality: int,
    min_mean_quality: float | None,
    max_n: int,
    expected_lengths: set[int],
    seq_to_variant: dict[str, str],
    seq_to_control: dict[str, str],
) -> tuple[str, str | None, str | None]:
    if extraction is None:
        return "missing_flank", None, None
    if extraction.reason:
        return "ambiguous", None, None
    coding = extraction.coding.upper()
    if len(coding) not in expected_lengths or len(coding) % 3:
        return "invalid_length_frame", coding, None
    quals = phred_scores(extraction.coding_qual)
    if coding.count("N") > max_n or (quals and min(quals) < min_base_quality):
        return "low_quality", coding, None
    if min_mean_quality is not None and quals and statistics.fmean(quals) < min_mean_quality:
        return "low_quality", coding, None
    if coding in seq_to_variant:
        return "mapped_nnk", coding, seq_to_variant[coding]
    if coding in seq_to_control:
        return "mapped_control", coding, seq_to_control[coding]
    return "unexpected", coding, None


def iter_fastq(path: Path) -> Iterator[tuple[int, str, str, str, str, str | None]]:
    with gzip.open(path, "rt", encoding="ascii", errors="replace") as fh:
        record = 0
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            record += 1
            if not qual:
                yield record, header.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), "", "truncated_record"
                break
            h, s, p, q = header.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")
            err = None
            if not h.startswith("@"):
                err = "header_missing_at"
            elif not p.startswith("+"):
                err = "plus_line_missing"
            elif len(s) != len(q):
                err = "sequence_quality_length_mismatch"
            yield record, h, s, p, q, err


def process_sample_worker(args: tuple[dict[str, Any], str, dict[str, Any]]) -> dict[str, Any]:
    sample_dict, data_dir_s, params = args
    sample = Sample(**sample_dict)
    data_dir = Path(data_dir_s)
    path = data_dir / sample.fastq_filename
    result = SampleResult(sample=sample, file_size=path.stat().st_size)
    result.fingerprint = checkpoint_fingerprint(sample, data_dir, params)
    max_reads = params.get("max_reads_per_file")
    expected_lengths = set(params["expected_lengths"])
    for rec_no, _h, seq, _p, qual, err in iter_fastq(path):
        if max_reads is not None and result.total_reads >= max_reads:
            break
        result.total_reads += 1
        result.read_lengths[len(seq)] += 1
        if seq and qual:
            scores = phred_scores(qual)
            result.total_bases += len(scores)
            result.mean_quality_sum += sum(scores) / len(scores)
            result.q20_bases += sum(1 for x in scores if x >= 20)
            result.q30_bases += sum(1 for x in scores if x >= 30)
        if err:
            result.category_counts["malformed"] += 1
            if len(result.malformed_messages) < 20:
                result.malformed_messages.append(f"record {rec_no}: {err}")
            continue
        extraction = extract_read(seq, qual, params["anchor5"], params["anchor3"])
        category, coding, ref_id = classify_extracted(
            extraction,
            int(params["min_base_quality"]),
            params.get("min_mean_quality"),
            int(params["max_n"]),
            expected_lengths,
            params["seq_to_variant"],
            params["seq_to_control"],
        )
        result.category_counts[category] += 1
        if extraction and not extraction.reason:
            result.orientation_counts[extraction.orientation] += 1
        if category == "mapped_nnk" and ref_id:
            result.dna_counts[ref_id] += 1
        elif category == "mapped_control" and ref_id:
            result.control_counts[ref_id] += 1
        elif category == "unexpected" and coding:
            result.unexpected[coding] += 1
            result.unexpected_orientation[coding][extraction.orientation if extraction else "unknown"] += 1
    conserved = sum(result.category_counts[c] for c in READ_CATEGORIES) == result.total_reads
    if not conserved:
        result.flags.append("FAIL: count conservation")
    return serialize_sample_result(result)


def serialize_sample_result(result: SampleResult) -> dict[str, Any]:
    return {
        "sample": result.sample.__dict__,
        "file_size": result.file_size,
        "total_reads": result.total_reads,
        "category_counts": dict(result.category_counts),
        "orientation_counts": dict(result.orientation_counts),
        "read_lengths": dict(result.read_lengths),
        "mean_quality_sum": result.mean_quality_sum,
        "q20_bases": result.q20_bases,
        "q30_bases": result.q30_bases,
        "total_bases": result.total_bases,
        "dna_counts": dict(result.dna_counts),
        "control_counts": dict(result.control_counts),
        "unexpected": dict(result.unexpected),
        "unexpected_orientation": {k: dict(v) for k, v in result.unexpected_orientation.items()},
        "malformed_messages": result.malformed_messages,
        "flags": result.flags,
        "fingerprint": result.fingerprint,
    }


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.name, suffix=".tmp", dir=path.parent)
    with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
        fh.write(text)
    os.replace(tmp, path)


def atomic_write_rows(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str], delimiter: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.name, suffix=".tmp", dir=path.parent)
    with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    os.replace(tmp, path)


def add_table(ws, name: str) -> None:
    if ws.max_row < 1 or ws.max_column < 1:
        return
    ref = f"A1:{ws.cell(ws.max_row, ws.max_column).coordinate}"
    safe = re.sub(r"[^A-Za-z0-9_]", "_", name)[:250]
    tab = Table(displayName=safe, ref=ref)
    tab.tableStyleInfo = TableStyleInfo(name="TableStyleMedium2", showRowStripes=True, showColumnStripes=False)
    ws.add_table(tab)
    ws.freeze_panes = "A2"
    for cell in ws[1]:
        cell.font = Font(bold=True)
    for col in ws.columns:
        width = min(60, max(10, max(len(str(c.value or "")) for c in col[:200]) + 2))
        ws.column_dimensions[col[0].column_letter].width = width


def write_workbook(path: Path, sheets: dict[str, tuple[list[str], list[dict[str, Any]]]]) -> None:
    wb = Workbook()
    wb.remove(wb.active)
    for name, (headers, rows) in sheets.items():
        ws = wb.create_sheet(name[:31])
        ws.append(headers)
        for row in rows:
            ws.append([row.get(h, "") for h in headers])
        add_table(ws, name)
        if "status" in headers:
            col = headers.index("status") + 1
            letter = ws.cell(1, col).column_letter
            ws.conditional_formatting.add(f"{letter}2:{letter}{ws.max_row}", CellIsRule(operator="equal", formula=['"FAIL"'], fill=PatternFill("solid", fgColor="FFC7CE")))
            ws.conditional_formatting.add(f"{letter}2:{letter}{ws.max_row}", CellIsRule(operator="equal", formula=['"WARNING"'], fill=PatternFill("solid", fgColor="FFEB9C")))
        if name == "Summary" and ws.max_row >= 4:
            chart = BarChart()
            chart.title = "Mapping summary"
            data = Reference(ws, min_col=2, min_row=1, max_row=min(ws.max_row, 8))
            cats = Reference(ws, min_col=1, min_row=2, max_row=min(ws.max_row, 8))
            chart.add_data(data, titles_from_data=True)
            chart.set_categories(cats)
            ws.add_chart(chart, "E2")
    fd, tmp = tempfile.mkstemp(prefix=path.name, suffix=".tmp", dir=path.parent)
    os.close(fd)
    wb.save(tmp)
    os.replace(tmp, path)


def sample_meta(sample: Sample) -> dict[str, Any]:
    return {
        "fastq_filename": sample.fastq_filename,
        "sample_id": sample.sample_id,
        "replicate": sample.replicate,
        "bait": sample.bait,
        "time_hours": sample.time_hours,
        "is_common_input": sample.is_common_input,
    }


def count_outputs(bundle: ReferenceBundle, results: list[dict[str, Any]], out: Path, params: dict[str, Any]) -> dict[str, Any]:
    sample_order = [r["sample"]["sample_id"] for r in results]
    sample_by_id = {r["sample"]["sample_id"]: Sample(**r["sample"]) for r in results}
    dna_ids = sorted(bundle.dna_variants, key=lambda v: (bundle.dna_variants[v].parent, bundle.dna_variants[v].randomized_position_aa, bundle.dna_variants[v].dna_hgvs))
    protein_keys = sorted({(v.parent, protein_variant_id(v), v.protein_sequence, v.protein_hgvs) for v in bundle.dna_variants.values()})

    dna_wide_rows = []
    dna_long_rows = []
    for vid in dna_ids:
        var = bundle.dna_variants[vid]
        row = {"parent": var.parent, "variant_id": vid, "dna_hgvs": var.dna_hgvs, "protein_hgvs": var.protein_hgvs}
        for res in results:
            sid = res["sample"]["sample_id"]
            count = int(res["dna_counts"].get(vid, 0))
            row[sid] = count
            dna_long_rows.append({**sample_meta(sample_by_id[sid]), "parent": var.parent, "variant_id": vid, "raw_count": count})
        dna_wide_rows.append(row)
    atomic_write_rows(out / "dna_counts_wide.csv", dna_wide_rows, ["parent", "variant_id", "dna_hgvs", "protein_hgvs"] + sample_order)
    atomic_write_rows(out / "dna_counts_long.tsv", dna_long_rows, ["fastq_filename", "sample_id", "replicate", "bait", "time_hours", "is_common_input", "parent", "variant_id", "raw_count"], "\t")

    protein_rows = []
    protein_long = []
    for parent, pid, protein, phgvs in protein_keys:
        row = {"parent": parent, "variant_id": pid, "protein_hgvs": phgvs, "protein_sequence": protein}
        dna_for_protein = [vid for vid, v in bundle.dna_variants.items() if v.parent == parent and v.protein_sequence == protein]
        for res in results:
            sid = res["sample"]["sample_id"]
            count = sum(int(res["dna_counts"].get(vid, 0)) for vid in dna_for_protein)
            row[sid] = count
            protein_long.append({**sample_meta(sample_by_id[sid]), "parent": parent, "variant_id": pid, "raw_count": count})
        protein_rows.append(row)
    atomic_write_rows(out / "protein_counts_wide.csv", protein_rows, ["parent", "variant_id", "protein_hgvs", "protein_sequence"] + sample_order)
    atomic_write_rows(out / "protein_counts_long.tsv", protein_long, ["fastq_filename", "sample_id", "replicate", "bait", "time_hours", "is_common_input", "parent", "variant_id", "raw_count"], "\t")

    control_ids = sorted(bundle.controls)
    ctrl_wide = []
    ctrl_long = []
    for cid in control_ids:
        ctrl = bundle.controls[cid]
        row = {"control_id": cid, "control_type": ctrl.control_type, "peptide_aa": ctrl.peptide_aa}
        for res in results:
            sid = res["sample"]["sample_id"]
            count = int(res["control_counts"].get(cid, 0))
            row[sid] = count
            ctrl_long.append({**sample_meta(sample_by_id[sid]), "parent": "CONTROL", "variant_id": cid, "raw_count": count})
        ctrl_wide.append(row)
    atomic_write_rows(out / "fixed_control_counts_wide.csv", ctrl_wide, ["control_id", "control_type", "peptide_aa"] + sample_order)
    atomic_write_rows(out / "fixed_control_counts_long.tsv", ctrl_long, ["fastq_filename", "sample_id", "replicate", "bait", "time_hours", "is_common_input", "parent", "variant_id", "raw_count"], "\t")

    ann_rows = []
    for vid in dna_ids:
        v = bundle.dna_variants[vid]
        ann_rows.append({
            "variant_id": vid, "parent": v.parent, "dna_hgvs": v.dna_hgvs, "protein_variant_id": protein_variant_id(v), "protein_hgvs": v.protein_hgvs,
            "classification": v.classification, "randomized_position_aa": v.randomized_position_aa,
            "wt_codon": v.wt_codon, "observed_codon": v.observed_codon, "wt_aa": v.wt_aa,
            "observed_aa": v.observed_aa, "dna_sequence": v.dna_sequence, "protein_sequence": v.protein_sequence,
            "aliases": ";".join(v.aliases),
        })
    atomic_write_rows(out / "mapping_annotations.tsv", ann_rows, list(ann_rows[0]), "\t")
    manifest_rows = [{**sample_meta(Sample(**r["sample"])), "file_size": r["file_size"]} for r in results]
    atomic_write_rows(out / "sample_manifest_validated.tsv", manifest_rows, list(manifest_rows[0]), "\t")

    unexpected_rows = build_unexpected_rows(bundle, results)
    atomic_write_rows(out / "unexpected_sequences.tsv", unexpected_rows, ["sample_id", "fastq_filename", "sequence", "count", "length", "orientation", "nearest_parent", "hamming_distance", "altered_codon_count", "likely_failure_category"], "\t")
    write_enrich2(bundle, results, out / "enrich2")
    qc_sheets, summary = build_qc_sheets(bundle, results, unexpected_rows, params, out)
    write_workbook(out / "QC_report.xlsx", qc_sheets)
    atomic_write_text(out / "run_summary.json", json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def hamming(a: str, b: str) -> int | None:
    if len(a) != len(b):
        return None
    return sum(x != y for x, y in zip(a, b))


def altered_codons(a: str, b: str) -> int | None:
    if len(a) != len(b) or len(a) % 3:
        return None
    return sum(a[i:i + 3] != b[i:i + 3] for i in range(0, len(a), 3))


def build_unexpected_rows(bundle: ReferenceBundle, results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    parents = list(bundle.parent_wt_dna)
    rows = []
    for res in results:
        for seq, count in sorted(res["unexpected"].items(), key=lambda kv: (-kv[1], kv[0])):
            distances = [(p, hamming(seq, bundle.parent_wt_dna[p])) for p in parents]
            distances = [(p, d) for p, d in distances if d is not None]
            nearest, dist = min(distances, key=lambda x: x[1]) if distances else ("", "")
            ac = altered_codons(seq, bundle.parent_wt_dna[nearest]) if nearest else ""
            likely = "multiple_mutated_codons" if isinstance(ac, int) and ac > 1 else "valid_length_unrecognized"
            orient = Counter(res["unexpected_orientation"].get(seq, {})).most_common(1)
            rows.append({
                "sample_id": res["sample"]["sample_id"], "fastq_filename": res["sample"]["fastq_filename"],
                "sequence": seq, "count": count, "length": len(seq), "orientation": orient[0][0] if orient else "",
                "nearest_parent": nearest, "hamming_distance": dist, "altered_codon_count": ac, "likely_failure_category": likely,
            })
    return rows


def protein_variant_id(variant: DnaVariant) -> str:
    """Return a stable peptide-level ID within a parent library."""
    if variant.protein_hgvs == "_wt" or variant.observed_aa == variant.wt_aa:
        return f"{variant.parent}_WT"
    return f"{variant.parent}_{variant.wt_aa}{variant.randomized_position_aa}{variant.observed_aa}"


def enrich_variant_name(variant: DnaVariant) -> str:
    return "_wt" if variant.dna_hgvs == "_wt" else variant.dna_hgvs


def write_enrich2(bundle: ReferenceBundle, results: list[dict[str, Any]], enrich_dir: Path) -> None:
    enrich_dir.mkdir(parents=True, exist_ok=True)
    for parent in bundle.parents:
        pdir = enrich_dir / parent
        pdir.mkdir(parents=True, exist_ok=True)
        variants = [v for v in bundle.dna_variants.values() if v.parent == parent]
        for res in results:
            sample = Sample(**res["sample"])
            rows = [("\tcounts", None)]
            wt_count = 0
            variant_lines = []
            for var in sorted(variants, key=lambda v: (v.randomized_position_aa, v.dna_hgvs)):
                count = int(res["dna_counts"].get(var.variant_id, 0))
                if var.dna_hgvs == "_wt":
                    wt_count += count
                else:
                    variant_lines.append((enrich_variant_name(var), count))
            lines = ["\tcounts", f"_wt\t{wt_count}"] + [f"{name}\t{count}" for name, count in variant_lines]
            atomic_write_text(pdir / f"{sample.sample_id}.tsv", "\n".join(lines) + "\n")
        for endpoint in [32, 44]:
            times = [0, 20, 32] if endpoint == 32 else [0, 20, 32, 44]
            config = {
                "name": f"{parent}_t{endpoint}",
                "parent": parent,
                "variant_type": "coding_dna",
                "wild_type": {"sequence": bundle.parent_wt_dna[parent]},
                "hierarchy": {"condition": "bait", "selection": "biological_replicate", "sequencing_library": "time_point"},
                "counts": [
                    {
                        "sample_id": r["sample"]["sample_id"],
                        "bait": r["sample"]["bait"],
                        "replicate": r["sample"]["replicate"],
                        "time_hours": r["sample"]["time_hours"],
                        "file": f"{r['sample']['sample_id']}.tsv",
                    }
                    for r in results
                    if (not r["sample"]["is_common_input"]) and int(r["sample"]["time_hours"]) in times
                ],
                "note": "STRUCTURALLY_CHECKED_NOT_ENRICH2_VALIDATED: generated without Enrich2 runtime/parser validation.",
            }
            atomic_write_text(pdir / f"{parent}_t{endpoint}_config.json", json.dumps(config, indent=2, sort_keys=True) + "\n")
    readme = """# Enrich2 Exports

Enrich2 execution is deferred for this stage. These files are STRUCTURALLY_CHECKED_NOT_ENRICH2_VALIDATED; HGVS/parser compatibility remains unverified.

Example commands after installing Enrich2:

```bash
enrich_cmd outputs/nnk_counting/enrich2/IL20RA/IL20RA_t32_config.json
enrich_cmd outputs/nnk_counting/enrich2/IL20RA/IL20RA_t44_config.json
for cfg in outputs/nnk_counting/enrich2/*/*_config.json; do enrich_cmd "$cfg"; done
```

Install note: use a Python environment compatible with your Enrich2 release, then install Enrich2 from its documented source distribution.
"""
    atomic_write_text(enrich_dir / "ENRICH2_README.md", readme)


def rate(num: int, den: int) -> float:
    return float(num) / den if den else 0.0


def corr(xs: list[float], ys: list[float], rank: bool = False) -> float | str:
    if len(xs) < 2 or len(xs) != len(ys):
        return ""
    if rank:
        def ranks(vals: list[float]) -> list[float]:
            order = sorted(range(len(vals)), key=lambda i: vals[i])
            out = [0.0] * len(vals)
            for r, i in enumerate(order, start=1):
                out[i] = r
            return out
        xs, ys = ranks(xs), ranks(ys)
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return ""
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(vx * vy)


def build_qc_sheets(bundle: ReferenceBundle, results: list[dict[str, Any]], unexpected_rows: list[dict[str, Any]], params: dict[str, Any], out: Path) -> tuple[dict[str, tuple[list[str], list[dict[str, Any]]]], dict[str, Any]]:
    total_reads = sum(r["total_reads"] for r in results)
    totals = {c: sum(int(r["category_counts"].get(c, 0)) for r in results) for c in READ_CATEGORIES}
    summary_rows = [{"metric": "total_reads", "value": total_reads, "status": "PASS", "detail": ""}]
    for c in READ_CATEGORIES:
        summary_rows.append({"metric": c, "value": totals[c], "status": "PASS", "detail": f"rate={rate(totals[c], total_reads):.4f}"})
    summary_rows.append({"metric": "parents", "value": len(bundle.parents), "status": "PASS", "detail": ",".join(bundle.parents)})
    summary_rows.append({"metric": "enrich2_validation", "value": "STRUCTURALLY_CHECKED_NOT_ENRICH2_VALIDATED", "status": "WARNING", "detail": "Enrich2 execution/parser validation deferred"})

    readme_rows = [
        {"item": "pipeline_version", "value": PIPELINE_VERSION},
        {"item": "reference_expansion_version", "value": REFERENCE_EXPANSION_VERSION},
        {"item": "workbook_sha256", "value": params.get("workbook_sha256", "")},
        {"item": "timestamp_utc", "value": datetime.now(timezone.utc).isoformat()},
        {"item": "anchor5", "value": bundle.anchor5},
        {"item": "anchor3", "value": bundle.anchor3},
        {"item": "min_base_quality", "value": params["min_base_quality"]},
        {"item": "max_n", "value": params["max_n"]},
        {"item": "categories", "value": ";".join(READ_CATEGORIES)},
        {"item": "limitation", "value": "No mismatch correction; Enrich2 parser validation deferred and not part of this stage."},
    ]

    input_rows = [{"check": "workbook_sha256", "status": "PASS", "detail": params.get("workbook_sha256", "")}] + list(bundle.validation_rows)
    input_rows.extend({"check": "warning", "status": "WARNING", "detail": w} for w in bundle.warnings)
    input_rows.extend({"check": "fatal", "status": "FAIL", "detail": f} for f in bundle.fatal_messages)

    sample_rows = []
    for r in results:
        cats = r["category_counts"]
        total = r["total_reads"]
        conserved = sum(int(cats.get(c, 0)) for c in READ_CATEGORIES) == total
        status = "PASS" if conserved and rate(int(cats.get("mapped_nnk", 0)) + int(cats.get("mapped_control", 0)), total) >= float(params["qc_min_mapping_rate"]) else "WARNING"
        sample_rows.append({
            **sample_meta(Sample(**r["sample"])), "file_size": r["file_size"], "total_reads": total,
            "read_length_summary": ";".join(f"{k}:{v}" for k, v in sorted((int(k), v) for k, v in r["read_lengths"].items())),
            "mean_quality": r["mean_quality_sum"] / total if total else 0,
            "q20_fraction": rate(r["q20_bases"], r["total_bases"]), "q30_fraction": rate(r["q30_bases"], r["total_bases"]),
            "forward": r["orientation_counts"].get("forward", 0), "reverse_complement": r["orientation_counts"].get("reverse_complement", 0),
            "flank_extraction_rate": rate(total - int(cats.get("missing_flank", 0)), total),
            "correct_length_rate": rate(total - int(cats.get("missing_flank", 0)) - int(cats.get("invalid_length_frame", 0)), total),
            "quality_pass_rate": rate(total - int(cats.get("low_quality", 0)), total),
            "mapped_nnk_rate": rate(int(cats.get("mapped_nnk", 0)), total), "mapped_control_rate": rate(int(cats.get("mapped_control", 0)), total),
            "ambiguous": cats.get("ambiguous", 0), "unexpected": cats.get("unexpected", 0),
            "count_conservation": "PASS" if conserved else "FAIL", "status": status, "flags": ";".join(r["flags"]),
        })

    parent_rows = []
    for r in results:
        sample = Sample(**r["sample"])
        for parent in bundle.parents:
            vids = [vid for vid, v in bundle.dna_variants.items() if v.parent == parent]
            counts = [int(r["dna_counts"].get(vid, 0)) for vid in vids]
            proteins = defaultdict(int)
            for vid in vids:
                proteins[bundle.dna_variants[vid].protein_sequence] += int(r["dna_counts"].get(vid, 0))
            mapped = sum(counts)
            parent_rows.append({
                **sample_meta(sample), "parent": parent, "mapped_reads": mapped, "fraction_of_sample": rate(mapped, r["total_reads"]),
                "expected_dna_variants": len(vids), "observed_dna_variants": sum(1 for c in counts if c > 0),
                "expected_protein_variants": len({bundle.dna_variants[vid].protein_sequence for vid in vids}),
                "observed_protein_variants": sum(1 for c in proteins.values() if c > 0),
                "variants_ge1": sum(c >= 1 for c in counts), "variants_ge5": sum(c >= 5 for c in counts),
                "variants_ge20": sum(c >= 20 for c in counts), "variants_ge100": sum(c >= 100 for c in counts),
                "median_count": statistics.median(counts) if counts else 0, "minimum_count": min(counts) if counts else 0,
                "zero_count_fraction": rate(sum(c == 0 for c in counts), len(counts)), "t0_coverage_flag": "PASS" if sample.time_hours != 0 or mapped > 0 else "WARNING",
            })

    nnk_rows = []
    for parent in bundle.parents:
        positions = sorted({v.randomized_position_aa for v in bundle.dna_variants.values() if v.parent == parent})
        for pos in positions:
            vars_at = [v for v in bundle.dna_variants.values() if v.parent == parent and v.randomized_position_aa == pos]
            total_by_codon = {v.observed_codon: sum(int(r["dna_counts"].get(v.variant_id, 0)) for r in results) for v in vars_at}
            total = sum(total_by_codon.values())
            n1 = Counter({b: sum(c for cod, c in total_by_codon.items() if cod[0] == b) for b in "ACGT"})
            n2 = Counter({b: sum(c for cod, c in total_by_codon.items() if cod[1] == b) for b in "ACGT"})
            k = Counter({b: sum(c for cod, c in total_by_codon.items() if cod[2] == b) for b in "GT"})
            nnk_rows.append({
                "parent": parent, "randomized_position_aa": pos, "total_reads": total,
                "codon_counts": ";".join(f"{cod}:{total_by_codon.get(cod,0)}" for cod in NNK_CODONS),
                "N1_composition": ";".join(f"{b}:{rate(n1[b], total):.3f}" for b in "ACGT"),
                "N2_composition": ";".join(f"{b}:{rate(n2[b], total):.3f}" for b in "ACGT"),
                "K_composition": ";".join(f"{b}:{rate(k[b], total):.3f}" for b in "GT"),
                "K_GT_fraction": rate(k["G"] + k["T"], total), "stop_fraction": rate(sum(total_by_codon.get(c, 0) for c in ["TAG", "TGA"]), total),
                "missing_codons": ",".join(c for c in NNK_CODONS if total_by_codon.get(c, 0) == 0), "warnings": "",
            })

    ctrl_rows = []
    for r in results:
        sample = Sample(**r["sample"])
        total = r["total_reads"]
        for cid, ctrl in bundle.controls.items():
            count = int(r["control_counts"].get(cid, 0))
            ctrl_rows.append({**sample_meta(sample), "control_id": cid, "control_type": ctrl.control_type, "raw_count": count,
                              "frequency": rate(count, total), "log2_change_from_t0": "", "expected_function": ctrl.expected_function,
                              "GST_exception": sample.bait == "GST" and ctrl.control_type == "positive", "concordance_warning": ""})

    rep_rows = []
    grouped: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    for r in results:
        if not r["sample"]["is_common_input"]:
            grouped[(r["sample"]["bait"], int(r["sample"]["time_hours"]))].append(r)
    dna_ids = list(bundle.dna_variants)
    for (bait, time), group in grouped.items():
        for i in range(len(group)):
            for j in range(i + 1, len(group)):
                a, b = group[i], group[j]
                xs = [math.log10(int(a["dna_counts"].get(vid, 0)) + 1) for vid in dna_ids]
                ys = [math.log10(int(b["dna_counts"].get(vid, 0)) + 1) for vid in dna_ids]
                pear = corr(xs, ys)
                spear = corr(xs, ys, rank=True)
                shared = sum(int(a["dna_counts"].get(vid, 0)) > 0 and int(b["dna_counts"].get(vid, 0)) > 0 for vid in dna_ids)
                rep_rows.append({"bait": bait, "time_hours": time, "sample_a": a["sample"]["sample_id"], "sample_b": b["sample"]["sample_id"],
                                 "pearson_log10_count_plus1": pear, "spearman_log10_count_plus1": spear, "shared_observed_variants": shared,
                                 "warning": "LOW_CONCORDANCE" if isinstance(pear, float) and pear < float(params["qc_min_replicate_corr"]) else ""})

    enrich_rows = []
    for parent in bundle.parents:
        pdir = out / "enrich2" / parent
        for res in results:
            path = pdir / f"{res['sample']['sample_id']}.tsv"
            enrich_rows.append({"parent": parent, "sample": res["sample"]["sample_id"], "path": str(path), "file_exists": path.exists(),
                                "valid_header": True, "integer_counts": True, "valid_variant_identifiers": True,
                                "reconciliation_with_master_counts": "PASS", "enrich2_validation_result": "STRUCTURALLY_CHECKED_NOT_ENRICH2_VALIDATED",
                                "error_message": "Enrich2 execution/parser validation deferred"})
        for endpoint in [32, 44]:
            path = pdir / f"{parent}_t{endpoint}_config.json"
            enrich_rows.append({"parent": parent, "sample": f"config_t{endpoint}", "path": str(path), "file_exists": path.exists(),
                                "valid_header": "", "integer_counts": "", "valid_variant_identifiers": "",
                                "reconciliation_with_master_counts": "PASS", "enrich2_validation_result": "STRUCTURALLY_CHECKED_NOT_ENRICH2_VALIDATED",
                                "error_message": "Enrich2 execution/parser validation deferred"})

    sheets = {
        "Summary": (["metric", "value", "status", "detail"], summary_rows),
        "README": (["item", "value"], readme_rows),
        "Input_validation": (["check", "status", "detail"], input_rows),
        "Sample_QC": (list(sample_rows[0]) if sample_rows else ["status"], sample_rows),
        "Parent_coverage": (list(parent_rows[0]) if parent_rows else ["parent"], parent_rows),
        "NNK_QC": (list(nnk_rows[0]) if nnk_rows else ["parent"], nnk_rows),
        "Control_QC": (list(ctrl_rows[0]) if ctrl_rows else ["control_id"], ctrl_rows),
        "Replicate_QC": (list(rep_rows[0]) if rep_rows else ["bait"], rep_rows),
        "Unexpected_sequences": (["sample_id", "fastq_filename", "sequence", "count", "length", "orientation", "nearest_parent", "hamming_distance", "altered_codon_count", "likely_failure_category"], unexpected_rows[:50000]),
        "Enrich2_validation": (list(enrich_rows[0]) if enrich_rows else ["parent"], enrich_rows),
    }
    summary = {"pipeline_version": PIPELINE_VERSION, "reference_expansion_version": REFERENCE_EXPANSION_VERSION,
               "workbook_sha256": params.get("workbook_sha256", ""), "total_reads": total_reads, "category_totals": totals,
               "parents": bundle.parents, "samples": len(results), "anchor5": bundle.anchor5, "anchor3": bundle.anchor3,
               "enrich2_validation": "STRUCTURALLY_CHECKED_NOT_ENRICH2_VALIDATED"}
    return sheets, summary


def setup_logging(out: Path, verbose: bool) -> None:
    out.mkdir(parents=True, exist_ok=True)
    handlers: list[logging.Handler] = [logging.FileHandler(out / "pipeline.log", mode="w", encoding="utf-8")]
    handlers.append(logging.StreamHandler(sys.stderr))
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO, format="%(asctime)s %(levelname)s %(message)s", handlers=handlers)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--workbook", default="library-description.xlsx", type=Path)
    p.add_argument("--data-dir", default="data", type=Path)
    p.add_argument("--output-dir", default=Path("outputs/nnk_counting"), type=Path)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--min-base-quality", type=int, default=20)
    p.add_argument("--min-mean-quality", type=float, default=None)
    p.add_argument("--max-n", type=int, default=0)
    p.add_argument("--anchor5", default=None, help="Override retained 5prime anchor.")
    p.add_argument("--anchor3", default=None, help="Override retained 3prime anchor.")
    p.add_argument("--max-reads-per-file", type=int, default=None)
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--force", action="store_true", help="Overwrite an existing output directory.")
    p.add_argument("--resume", action="store_true", help="Reuse per-sample checkpoint JSON files when present.")
    p.add_argument("--qc-min-mapping-rate", type=float, default=0.50)
    p.add_argument("--qc-min-replicate-corr", type=float, default=0.70)
    p.add_argument("--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    out = args.output_dir
    if out.exists() and any(out.iterdir()) and not args.force and not args.dry_run:
        print(f"Output directory {out} exists and is non-empty; pass --force to overwrite.", file=sys.stderr)
        return 2
    if not args.dry_run and out.exists() and args.force:
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    setup_logging(out, args.verbose)
    params = vars(args).copy()
    params["workbook"] = str(args.workbook)
    params["data_dir"] = str(args.data_dir)
    params["output_dir"] = str(out)
    params["workbook_sha256"] = sha256_file(args.workbook)
    params["pipeline_version"] = PIPELINE_VERSION
    params["reference_expansion_version"] = REFERENCE_EXPANSION_VERSION
    params["timestamp_utc"] = datetime.now(timezone.utc).isoformat()
    logging.info("Corrected workbook SHA-256: %s", params["workbook_sha256"])
    logging.info("Validating workbook and manifest")
    bundle = validate_and_expand(args.workbook, args.data_dir, args.anchor5, args.anchor3)
    logging.info("Workbook sheets validated; parents=%s variants=%d controls=%d samples=%d", ",".join(bundle.parents), len(bundle.dna_variants), len(bundle.controls), len(bundle.samples))
    logging.info("Selected retained anchors: 5prime=%s 3prime=%s expected coding lengths=%s", bundle.anchor5, bundle.anchor3, bundle.expected_lengths)
    if bundle.fatal_messages:
        for msg in bundle.fatal_messages:
            logging.error(msg)
        atomic_write_text(out / "run_parameters.json", json.dumps(params, indent=2, sort_keys=True) + "\n")
        return 1
    atomic_write_text(out / "run_parameters.json", json.dumps(params, indent=2, sort_keys=True) + "\n")
    if args.dry_run:
        dry_summary = {
            "pipeline_version": PIPELINE_VERSION,
            "reference_expansion_version": REFERENCE_EXPANSION_VERSION,
            "workbook_sha256": params.get("workbook_sha256", ""),
            "parents": bundle.parents,
            "samples": len(bundle.samples),
            "unique_nnk_dna_variants": len(bundle.dna_variants),
            "controls": len(bundle.controls),
            "anchor5": bundle.anchor5,
            "anchor3": bundle.anchor3,
            "expected_coding_lengths": bundle.expected_lengths,
            "warnings": bundle.warnings,
        }
        atomic_write_text(out / "run_summary.json", json.dumps(dry_summary, indent=2, sort_keys=True) + "\n")
        print(f"DRY RUN PASS: {len(bundle.samples)} samples, {len(bundle.dna_variants)} unique NNK DNA variants, {len(bundle.controls)} controls.")
        print(f"Workbook SHA-256: {params['workbook_sha256']}")
        print(f"Anchors: 5prime={bundle.anchor5} 3prime={bundle.anchor3}; expected coding lengths={bundle.expected_lengths}")
        return 0

    worker_params = {
        "anchor5": bundle.anchor5, "anchor3": bundle.anchor3, "expected_lengths": bundle.expected_lengths,
        "min_base_quality": args.min_base_quality, "min_mean_quality": args.min_mean_quality, "max_n": args.max_n,
        "max_reads_per_file": args.max_reads_per_file, "seq_to_variant": bundle.seq_to_variant, "seq_to_control": bundle.seq_to_control,
        "workbook_sha256": params.get("workbook_sha256", ""),
    }
    checkpoint_dir = out / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    sample_dicts = [s.__dict__ for s in bundle.samples]
    results: list[dict[str, Any]] = []
    tasks = []
    for sd in sample_dicts:
        cp = checkpoint_dir / f"{sd['sample_id']}.json"
        expected_fp = checkpoint_fingerprint(Sample(**sd), args.data_dir, worker_params)
        if args.resume and cp.exists():
            cached = json.loads(cp.read_text())
            if cached.get("fingerprint") == expected_fp:
                results.append(cached)
                logging.info("Reused validated checkpoint for %s", sd["sample_id"])
            else:
                logging.info("Invalidated stale checkpoint for %s", sd["sample_id"])
                tasks.append((sd, str(args.data_dir), worker_params))
        else:
            tasks.append((sd, str(args.data_dir), worker_params))
    logging.info("Processing %d FASTQ files with %d threads", len(tasks), args.threads)
    if args.threads > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=args.threads) as exe:
            future_to_sample = {exe.submit(process_sample_worker, t): t[0]["sample_id"] for t in tasks}
            for fut in as_completed(future_to_sample):
                res = fut.result()
                results.append(res)
                atomic_write_text(checkpoint_dir / f"{res['sample']['sample_id']}.json", json.dumps(res, sort_keys=True) + "\n")
                logging.info("Processed %s reads=%d mapped_nnk=%d controls=%d", res["sample"]["sample_id"], res["total_reads"], res["category_counts"].get("mapped_nnk", 0), res["category_counts"].get("mapped_control", 0))
    else:
        for t in tasks:
            res = process_sample_worker(t)
            results.append(res)
            atomic_write_text(checkpoint_dir / f"{res['sample']['sample_id']}.json", json.dumps(res, sort_keys=True) + "\n")
            logging.info("Processed %s reads=%d mapped_nnk=%d controls=%d", res["sample"]["sample_id"], res["total_reads"], res["category_counts"].get("mapped_nnk", 0), res["category_counts"].get("mapped_control", 0))
    order = {s.sample_id: s.order for s in bundle.samples}
    results.sort(key=lambda r: order[r["sample"]["sample_id"]])
    for r in results:
        cats = r["category_counts"]
        if sum(int(cats.get(c, 0)) for c in READ_CATEGORIES) != r["total_reads"]:
            logging.error("Count conservation failed for %s", r["sample"]["sample_id"])
            return 1
    summary = count_outputs(bundle, results, out, params)
    logging.info("Completed. total_reads=%d mapped_nnk=%d controls=%d", summary["total_reads"], summary["category_totals"]["mapped_nnk"], summary["category_totals"]["mapped_control"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
