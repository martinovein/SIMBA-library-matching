import sys
from pathlib import Path
from collections import Counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import count_nnk_fastq as base
import count_nnk_fastq_tolerant as tol

A5='TATGCAGCCTCTTCATCTGGC'
A3='GGTGGAGGATCCGGAGG'
INS='ATGAAA'
CTRL='GGTCCA'
QUAL='I'


def classify(seq, max5=0, max3=0, insert=INS, q=None, refs=None, ctrls=None):
    if q is None: q=QUAL*len(seq)
    refs = refs or {INS:'v1'}
    ctrls = ctrls or {CTRL:'ctrl'}
    ext=tol.extract_read_tolerant(seq, q, A5, A3, max5, max3, {6})
    return ext, tol.classify_tolerant(ext, 20, None, 0, refs, ctrls)


def test_exact_anchors_behave_like_baseline():
    seq=A5+INS+A3; q=QUAL*len(seq)
    be=base.extract_read(seq,q,A5,A3)
    te=tol.extract_read_tolerant(seq,q,A5,A3,0,0,{6})
    assert be.coding == te.coding == INS
    assert classify(seq)[1] == ('mapped_nnk', INS, 'v1')


def test_one_5prime_anchor_substitution_can_be_rescued():
    seq='C'+A5[1:]+INS+A3
    ext, cls=classify(seq, max5=1)
    assert cls[0]=='mapped_nnk'
    assert ext.anchor5_mismatches==1 and ext.anchor5_positions==(1,)


def test_one_3prime_anchor_substitution_can_be_rescued():
    seq=A5+INS+A3[:-1]+'A'
    ext, cls=classify(seq, max3=1)
    assert cls[0]=='mapped_nnk'
    assert ext.anchor3_mismatches==1 and ext.anchor3_positions==(len(A3),)


def test_one_substitution_in_each_anchor_can_be_rescued():
    seq='C'+A5[1:]+INS+A3[:-1]+'A'
    ext, cls=classify(seq, max5=1, max3=1)
    assert cls[0]=='mapped_nnk'
    assert (ext.anchor5_mismatches, ext.anchor3_mismatches)==(1,1)


def test_mismatch_above_threshold_rejected():
    seq='CC'+A5[2:]+INS+A3
    ext, cls=classify(seq, max5=1)
    assert cls[0]=='missing_flank'


def test_anchor_insertion_not_rescued():
    seq=A5[:4]+'A'+A5[4:]+INS+A3
    ext, cls=classify(seq, max5=1)
    assert cls[0] in {'missing_flank','ambiguous'}
    assert cls[0] != 'mapped_nnk'


def test_anchor_deletion_not_rescued():
    seq=A5[:4]+A5[5:]+INS+A3
    ext, cls=classify(seq, max5=1)
    assert cls[0] in {'missing_flank','ambiguous'}
    assert cls[0] != 'mapped_nnk'


def test_53bp_insert_not_converted_to_54bp():
    seq=A5+('A'*53)+A3
    ext=tol.extract_read_tolerant(seq, QUAL*len(seq), A5,A3,1,1,{54})
    assert ext.reason == 'invalid_length_frame'
    assert tol.classify_tolerant(ext,20,None,0,{}, {})[0] == 'invalid_length_frame'


def test_55bp_insert_not_converted_to_54bp():
    seq=A5+('A'*55)+A3
    ext=tol.extract_read_tolerant(seq, QUAL*len(seq), A5,A3,1,1,{54})
    assert ext.reason == 'invalid_length_frame'
    assert tol.classify_tolerant(ext,20,None,0,{}, {})[0] == 'invalid_length_frame'


def test_expected_insert_assigned_after_tolerant_extraction():
    seq='C'+A5[1:]+INS+A3
    assert classify(seq, max5=1)[1] == ('mapped_nnk', INS, 'v1')


def test_unrecognized_insert_not_nearest_assigned():
    seq='C'+A5[1:]+'TTTCCC'+A3
    ext, cls=classify(seq, max5=1)
    assert cls[0]=='unexpected' and cls[2] is None


def test_multiple_equally_valid_anchor_placements_ambiguous():
    # Two possible 5prime anchors generate the same expected-length insert.
    seq=A5 + 'A' + A5 + INS + A3
    ext=tol.extract_read_tolerant(seq, QUAL*len(seq), A5, A3, 1, 0, {6, 6+len(A5)+1})
    assert ext.reason == 'multiple_anchor_candidates'


def test_reverse_quality_coordinates_correct():
    fwd=A5+INS+A3
    rc=base.reverse_complement(fwd)
    qual='5'*len(A5)+QUAL*len(INS)+'!'*len(A3)
    ext=tol.extract_read_tolerant(rc, qual[::-1], A5,A3,0,0,{6})
    assert ext.orientation=='reverse_complement'
    assert ext.coding_qual == QUAL*len(INS)


def test_terminal_count_conservation_and_controls_separate(tmp_path):
    sample=base.Sample('x.fastqsanger.gz','x','r','b',0,False,0)
    data=tmp_path
    import gzip
    reads=[A5+INS+A3, A5+CTRL+A3, 'C'+A5[1:]+INS+A3, A5+'TTTCCC'+A3]
    with gzip.open(data/'x.fastqsanger.gz','wt') as fh:
        for i,r in enumerate(reads): fh.write(f'@r{i}\n{r}\n+\n{QUAL*len(r)}\n')
    params={'anchor5':A5,'anchor3':A3,'expected_lengths':[6],'min_base_quality':20,'min_mean_quality':None,'max_n':0,'max_reads_per_file':None,'seq_to_variant':{INS:'v1'},'seq_to_control':{CTRL:'ctrl'},'workbook_sha256':'x','max_5prime_anchor_mismatches':1,'max_3prime_anchor_mismatches':0}
    res=tol.process_sample_worker((sample.__dict__, str(data), params))
    assert sum(res['category_counts'].get(k,0) for k in base.READ_CATEGORIES)==res['total_reads']==4
    assert res['dna_counts']=={'v1':2}
    assert res['control_counts']=={'ctrl':1}


def test_valid_reference_cannot_be_valid_length_unrecognized():
    seq='C'+A5[1:]+INS+A3
    ext, cls=classify(seq, max5=1)
    assert cls[0]=='mapped_nnk'
