import os
import pytest
from seqtools import DNASequence, RNASequence, AminoAcidSequence, filter_fastq


@pytest.fixture
def input_fastq(tmp_path):
    filepath = tmp_path / "input.fastq"
    with open(filepath, "w") as f:
        f.write("@good_read\nATGCATGC\n+\nIIIIIIII\n")
        f.write("@bad_read\nATGCATGC\n+\n!!!!!!!!\n")
    yield filepath

# при таком тесте у меня не удаляются output файлы после теста,
# но как сделать fixture на output не меняя код функции filter_fastq я не успел придумать
def test_filter_fastq(input_fastq):
    filter_fastq(str(input_fastq), "out.fastq", quality_threshold=20)
    with open(os.path.join("filtered", "out.fastq"), "r") as f:
        lines = f.readlines()
    assert lines == ["@good_read\n", "ATGCATGC\n", "+\n", "IIIIIIII\n"]


class TestDNASequence:
    def test_length(self):
        seq = DNASequence("ATGCTTCGA")
        assert len(seq) == 9

    def test_complement(self):
        seq = DNASequence("ATGC")
        assert str(seq.complement()) == "TACG"

    def test_reverse(self):
        seq = DNASequence("ATGAAAC")
        assert str(seq.reverse()) == "CAAAGTA"

    def test_invalid_dna_raises_error(self):
        with pytest.raises(ValueError):
            DNASequence("AUGC")


class TestAminoAcidSequence:
    def test_trypsin_sites(self):
        seq = AminoAcidSequence("AKPpMKrGc")
        assert seq.trypsin_sites() == [5, 6]


class TestRNASequence:
    def test_invalid_rna_raises_error(self):
        with pytest.raises(ValueError):
            RNASequence("ATGC")

    def test_reverse_complement_rna(self):
        seq = RNASequence("GCAU")
        assert str(seq.reverse_complement()) == "AUGC"