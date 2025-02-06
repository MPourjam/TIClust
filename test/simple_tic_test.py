import unittest
import sys
import os
import pytest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
import simple_tic as stic
import tic_helper as tich
import logging
import tempfile
import pathlib as pl


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class TestTaxonomy:

    def test_taxonomy_initialization(self):
        tax_str = "tax=kingdom;phylum;class;order;family;genus;species;"
        taxonomy = stic.Taxonomy(tax_str)
        # log tax_str
        logging.info(taxonomy.tax_str)
        assert taxonomy.tax_str == "kingdom;phylum;class;order;family;genus;species"
        assert taxonomy.kingdom == "kingdom"
        assert taxonomy.phylum == "phylum"
        assert taxonomy.class_ == "class"
        assert taxonomy.order == "order"
        assert taxonomy.family == "family"
        assert taxonomy.genus == "genus"
        assert taxonomy.species == "species"

    def test_get_tax_upto(self):
        tax_str = "tax=kingdom;phylum;class;order;family;genus;species;"
        taxonomy = stic.Taxonomy(tax_str)
        assert taxonomy.get_tax_upto('family').tax_str == "kingdom;phylum;class;order;family"

    def test_clean_taxonomy(self):
        tax_str = "tax=kingdom;phylum;unclass;unorder;family;genus;species;"
        taxonomy = stic.Taxonomy(tax_str)
        assert taxonomy.tax_list == ["kingdom", "phylum"]
        assert str(taxonomy) == "tax=kingdom;phylum;NA-Class;NA-Order;family;genus;species"
        assert taxonomy.get_tax_upto('family', ret_type = 'str') == "kingdom;phylum;NA-Class;NA-Order;NA-Family"
        # missed family level so we cut from invalid or missed level
        tax_str = "tax=kingdom;phylum;class;order;;genus;species;"
        taxonomy = stic.Taxonomy(tax_str)
        assert taxonomy.tax_list == ["kingdom", "phylum", "class", "order"]
        assert str(taxonomy) == "tax=kingdom;phylum;class;order;NA-Family;genus;species"

    def test_full_tax(self):
        tax_str = "tax=kingdom;phylum;class;order;family;genus;species;"
        taxonomy = stic.Taxonomy(tax_str)
        assert taxonomy.full_tax == "kingdom;phylum;class;order;family;genus;species"

    def test_last_known_level(self):
        tax_str = "tax=kingdom;phylum;class;order;family;genus;species;"
        taxonomy = stic.Taxonomy(tax_str)
        assert taxonomy.last_known_level == "species"
        tax_str = "tax=kingdom;phylum;class;order;family;"
        taxonomy = stic.Taxonomy(tax_str)
        assert taxonomy.last_known_level == "family"
        tax_str = "tax=kingdom;phylum;;order;family;"
        taxonomy = stic.Taxonomy(tax_str)
        assert taxonomy.last_known_level == "phylum"

    def test_last_known_tax(self):
        tax_str = "tax=kingdomA;phylumA;classA;orderA;familyA;genusA;speciesA;"
        taxonomy = stic.Taxonomy(tax_str)
        assert taxonomy.last_known_tax == "speciesA"
        tax_str = "tax=kingdomA;phylumA;classA;orderA;familyA;;"
        taxonomy = stic.Taxonomy(tax_str)
        assert taxonomy.last_known_tax == "familyA"
        tax_str = "tax=kingdomA;phylumA;classA;;familyA;"
        taxonomy = stic.Taxonomy(tax_str)
        assert taxonomy.last_known_tax == "classA"

    def test_is_known_upto(self):
        tax_str = "tax=kingdom;phylum;class;order;family;genus;species;"
        taxonomy = stic.Taxonomy(tax_str)
        assert taxonomy.is_known_upto("family")
        assert taxonomy.is_known_upto("species")
        tax_str = "tax=kingdom;phylum;;order;;genus;;"
        taxonomy = stic.Taxonomy(tax_str)
        assert not taxonomy.is_known_upto("family")
        assert taxonomy.is_known_upto("phylum")

    def test_set_level(self):
        tax_str = "tax=kingdom;phylum;class;order;family;genus;species;"
        taxonomy = stic.Taxonomy(tax_str)
        taxonomy.set_level("family", "new_family")
        assert taxonomy.family == "new_family"
        assert taxonomy.full_tax == "kingdom;phylum;class;order;new_family;genus;species"
        assert taxonomy.tax_str == "kingdom;phylum;class;order;new_family;genus;species"
        tax_str = "tax=kingdom;phylum;class;order;;;;"
        taxonomy = stic.Taxonomy(tax_str)
        taxonomy.set_level("family", "new_family")
        assert taxonomy.family == "new_family"
        assert taxonomy.tax_str == "kingdom;phylum;class;order;new_family"
        taxonomy.set_level("species", "new_species")
        assert taxonomy.species == "new_species"
        assert taxonomy.full_tax == "kingdom;phylum;class;order;new_family;NA-Genus;new_species"
        assert taxonomy.tax_str == "kingdom;phylum;class;order;new_family"
        with pytest.raises(ValueError):
            taxonomy.set_level("ZOTU", "new_ZOTU")
    
    def test_get_clean_taxonomy(self):
        tax_str = "tax=aBacteria;NA_phylum__aBacteria;NA_class__aBacteria;NA_order__aBacteria;"
        taxonomy = stic.Taxonomy(tax_str)
        assert taxonomy.get_clean_taxonomy(tax_str) == tax_str[:-1]
        full_tax = taxonomy.get_clean_taxonomy(tax_str, force_full_path=True)
        assert full_tax == tax_str + "NA-Family;NA-Genus;NA-Species"


class TestSeqID:

    def test_seq_id_initialization(self):
        header = ">seq1 some description"
        seq_id = stic.SeqID(header)
        assert seq_id.head_id == "seq1"
        header = ">seq1 ; some; description"
        seq_id = stic.SeqID(header)
        assert seq_id.head_id == "seq1"


class TestSeqHeader:

    def test_seq_header_initialization(self):
        header = ">seq1 tax=kingdom;phylum;class;order;family;genus;species;"
        seq_header = stic.SeqHeader(header)
        assert str(seq_header.seq_id) == "seq1"
        assert isinstance(seq_header.taxonomy, stic.Taxonomy)
        assert seq_header.taxonomy.tax_str == "kingdom;phylum;class;order;family;genus;species"
        assert str(seq_header.taxonomy) == "tax=kingdom;phylum;class;order;family;genus;species"


class TestSequence:

    def test_sequence_initialization(self):
        header = ">seq1 tax=kingdom;phylum;class;order;family;genus;species;"
        sequence = "ATCGATCGATCG"
        seq = stic.Sequence(header, sequence)
        assert str(seq.header.seq_id) == "seq1"
        assert seq.sequence == "ATCGATCGATCG"

    def test_is_sequence_correct(self):
        header = ">seq1 tax=kingdom;phylum;class;order;family;genus;species;"
        sequence = "ATCGATCGATCG"
        seq = stic.Sequence(header, sequence)
        assert seq.is_sequence_correct()
    
    def test_get_hash(self):
        header = ">seq1 tax=kingdom;phylum;class;order;family;genus;species;"
        sequence = "ATCGATCGATCG"
        seq = stic.Sequence(header, sequence)
        assert isinstance(seq.__hash__(), int)
        assert seq.__hash__() == hash(seq)
    
    def test_get_tax_upto(self):
        header = ">seq1 tax=kingdom;phylum;class;order;family;genus;species;"
        sequence = "ATCGATCGATCG"
        seq = stic.Sequence(header, sequence)
        assert seq.header.taxonomy.get_tax_upto('family', ret_type = "str") == "kingdom;phylum;class;order;family"

    def test_equal(self):
        header = ">seq1 tax=kingdom;phylum;class;order;family;genus;species;"
        sequence = "ATCGATCGATCG"
        seq1 = stic.Sequence(header, sequence)
        seq2 = stic.Sequence(header, sequence)
        assert seq1 == seq2
        header = ">seq2 tax=kingdom;phylum;class;order;family;genus;species;"
        seq3 = stic.Sequence(header, sequence)
        assert seq1 != seq3


class TestSequenceCluster:

    def setup_method(self):
        self.header1 = ">seq1 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence1 = "ATCGATCGATCG"
        self.seq1 = stic.Sequence(self.header1, self.sequence1)

        self.header2 = ">seq2 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence2 = "GCTAGCTAGCTA"
        self.seq2 = stic.Sequence(self.header2, self.sequence2)

        self.header3 = ">seq3 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence3 = "CGTACGTACGTA"
        self.seq3 = stic.Sequence(self.header3, self.sequence3)

        self.header4 = ">seq4 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence4 = "TACGTACGTACG"
        self.seq4 = stic.Sequence(self.header4, self.sequence4)

        self.header5 = ">seq4 tax=kingdom;phylum;class;order;;genus;species;"
        self.sequence5 = "TACGTACGTACG"
        self.seq5 = stic.Sequence(self.header5, self.sequence5)

        self.sequences = [self.seq1, self.seq2, self.seq3, self.seq4]
        self.non_homogenous_sequences = [self.seq1, self.seq2, self.seq3, self.seq5]

    def test_sequence_cluster_initialization(self):
        cluster = stic.SequenceCluster(self.sequences, centroid=self.seq1)
        assert cluster.centroid == self.seq1
        assert len(cluster.sequences) == 4
        assert cluster.homogenous_tax_group
        non_homogenous_sequences = stic.SequenceCluster(
            self.non_homogenous_sequences,
            centroid=self.seq1,
            force_homogeneity=False
        )
        assert not non_homogenous_sequences.homogenous_tax_group


    def test_add_sequence(self):
        cluster = stic.SequenceCluster(self.sequences, centroid=self.seq1)
        new_header = ">seq5 tax=kingdom;phylum;class;order;family;genus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        cluster.add_sequence(new_seq)
        assert len(cluster.sequences) == 5
        assert cluster.homogenous_tax_group

    def test_set_level(self):
        cluster = stic.SequenceCluster(self.sequences, centroid=self.seq1)
        cluster.set_level("family", "new_family")
        for seq in cluster.sequences:
            assert seq.header.taxonomy.family == "new_family"

    def test_closest_common_ancestor(self):
        cluster = stic.SequenceCluster(self.sequences, centroid=self.seq1)
        common_ancestor = cluster.closest_common_ancestor
        assert common_ancestor.tax_str == "kingdom;phylum;class;order;family;genus;species"
        non_homogenous_sequences = stic.SequenceCluster(
            self.non_homogenous_sequences,
            centroid=self.seq1,
            force_homogeneity=False
        )
        assert non_homogenous_sequences.closest_common_ancestor.tax_str == "kingdom;phylum;class;order"

    def test_is_homogenous_group(self):
        assert stic.SequenceCluster.is_homogenous_group(self.sequences, upto_level="species")

    def test_cluster_size(self):
        cluster = stic.SequenceCluster(self.sequences, centroid=self.seq1)
        assert cluster.cluster_size == 4

    def test_write_to_fasta(self, tmp_path):
        cluster = stic.SequenceCluster(self.sequences, centroid=self.seq1)
        output_file_path = tmp_path.joinpath("test_output.fasta")
        cluster.write_to_fasta(output_file_path)
        with open(output_file_path, 'r') as f:
            lines = f.readlines()
        assert len(lines) == 8  # 5 sequences, each with a header and sequence line
        output_file_path.unlink()  # Clean up the test file


class TestKingdomCluster:

    def setup_method(self):
        self.header1 = ">seq1 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence1 = "ATCGATCGATCG"
        self.seq1 = stic.Sequence(self.header1, self.sequence1)

        self.header2 = ">seq2 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence2 = "GCTAGCTAGCTA"
        self.seq2 = stic.Sequence(self.header2, self.sequence2)

        self.header3 = ">seq3 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence3 = "CGTACGTACGTA"
        self.seq3 = stic.Sequence(self.header3, self.sequence3)

        self.header4 = ">seq4 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence4 = "TACGTACGTACG"
        self.seq4 = stic.Sequence(self.header4, self.sequence4)

        self.sequences = [self.seq1, self.seq2, self.seq3, self.seq4]

    def test_kingdom_cluster_initialization(self):
        cluster = stic.KingdomCluster(self.sequences, centroid=self.seq1)
        assert cluster.centroid == self.seq1
        assert len(cluster.sequences) == 4
        assert cluster.homogenous_tax_group

    def test_add_sequence(self):
        cluster = stic.KingdomCluster(self.sequences, centroid=self.seq1)
        new_header = ">seq5 tax=kingdom;phylum;class;order;family;genus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        cluster.add_sequence(new_seq)
        assert len(cluster.sequences) == 5
        assert cluster.homogenous_tax_group
        new_header = ">seq5 tax=Anotherkingdom;phylum;class;order;family;genus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        with pytest.raises(ValueError):
            cluster.add_sequence(new_seq)

    def test_set_level(self):
        cluster = stic.KingdomCluster(self.sequences, centroid=self.seq1)
        cluster.set_level("family", "new_family")
        for seq in cluster.sequences:
            assert seq.header.taxonomy.family == "new_family"

    def test_closest_common_ancestor(self):
        cluster = stic.KingdomCluster(self.sequences, centroid=self.seq1)
        common_ancestor = cluster.closest_common_ancestor
        assert common_ancestor.tax_str == "kingdom;phylum;class;order;family;genus;species"

    def test_is_homogenous_group(self):
        assert stic.KingdomCluster.is_homogenous_group(self.sequences, upto_level="kingdom")

    def test_write_to_fasta(self, tmp_path):
        cluster = stic.KingdomCluster(self.sequences, centroid=self.seq1)
        output_file_path = tmp_path.joinpath("test_output.fasta")
        cluster.write_to_fasta(output_file_path)
        with open(output_file_path, 'r') as f:
            lines = f.readlines()
        assert len(lines) == 8  # 5 sequences, each with a header and sequence line
        output_file_path.unlink()  # Clean up the test file


class TestPhylumCluster:

    def setup_method(self):
        self.header1 = ">seq1 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence1 = "ATCGATCGATCG"
        self.seq1 = stic.Sequence(self.header1, self.sequence1)

        self.header2 = ">seq2 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence2 = "GCTAGCTAGCTA"
        self.seq2 = stic.Sequence(self.header2, self.sequence2)

        self.header3 = ">seq3 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence3 = "CGTACGTACGTA"
        self.seq3 = stic.Sequence(self.header3, self.sequence3)

        self.header4 = ">seq4 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence4 = "TACGTACGTACG"
        self.seq4 = stic.Sequence(self.header4, self.sequence4)

        self.sequences = [self.seq1, self.seq2, self.seq3, self.seq4]

    def test_phylum_cluster_initialization(self):
        cluster = stic.PhylumCluster(self.sequences, centroid=self.seq1)
        assert cluster.centroid == self.seq1
        assert len(cluster.sequences) == 4
        assert cluster.homogenous_tax_group

    def test_add_sequence(self):
        cluster = stic.PhylumCluster(self.sequences, centroid=self.seq1)
        new_header = ">seq5 tax=kingdom;phylum;class;order;family;genus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        cluster.add_sequence(new_seq)
        assert len(cluster.sequences) == 5
        assert cluster.homogenous_tax_group
        new_header = ">seq5 tax=kingdom;Anotherphylum;class;order;family;genus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        with pytest.raises(ValueError):
            cluster.add_sequence(new_seq)
        new_header = ">seq5 tax=Anotherkingdom;phylum;class;order;family;genus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        with pytest.raises(ValueError):
            cluster.add_sequence(new_seq)

    def test_set_level(self):
        cluster = stic.PhylumCluster(self.sequences, centroid=self.seq1)
        cluster.set_level("family", "new_family")
        for seq in cluster.sequences:
            assert seq.header.taxonomy.family == "new_family"

    def test_closest_common_ancestor(self):
        cluster = stic.PhylumCluster(self.sequences, centroid=self.seq1)
        common_ancestor = cluster.closest_common_ancestor
        assert common_ancestor.tax_str == "kingdom;phylum;class;order;family;genus;species"

    def test_is_homogenous_group(self):
        assert stic.PhylumCluster.is_homogenous_group(self.sequences, upto_level="phylum")

    def test_write_to_fasta(self, tmp_path):
        cluster = stic.PhylumCluster(self.sequences, centroid=self.seq1)
        output_file_path = tmp_path.joinpath("test_output.fasta")
        cluster.write_to_fasta(output_file_path)
        with open(output_file_path, 'r') as f:
            lines = f.readlines()
        assert len(lines) == 8  # 4 sequences, each with a header and sequence line
        output_file_path.unlink()


class TestClassCluster:
    
    def setup_method(self):
        self.header1 = ">seq1 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence1 = "ATCGATCGATCG"
        self.seq1 = stic.Sequence(self.header1, self.sequence1)

        self.header2 = ">seq2 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence2 = "GCTAGCTAGCTA"
        self.seq2 = stic.Sequence(self.header2, self.sequence2)

        self.header3 = ">seq3 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence3 = "CGTACGTACGTA"
        self.seq3 = stic.Sequence(self.header3, self.sequence3)

        self.header4 = ">seq4 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence4 = "TACGTACGTACG"
        self.seq4 = stic.Sequence(self.header4, self.sequence4)

        self.sequences = [self.seq1, self.seq2, self.seq3, self.seq4]

    def test_class_cluster_initialization(self):
        cluster = stic.ClassCluster(self.sequences, centroid=self.seq1)
        assert cluster.centroid == self.seq1
        assert len(cluster.sequences) == 4
        assert cluster.homogenous_tax_group

    def test_add_sequence(self):
        cluster = stic.ClassCluster(self.sequences, centroid=self.seq1)
        new_header = ">seq5 tax=kingdom;phylum;class;order;family;genus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        cluster.add_sequence(new_seq)
        assert len(cluster.sequences) == 5
        assert cluster.homogenous_tax_group
        new_header = ">seq5 tax=kingdom;phylum;Anotherclass;order;family;genus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        with pytest.raises(ValueError):
            cluster.add_sequence(new_seq)
        new_header = ">seq5 tax=Anotherkingdom;phylum;class;order;family;genus;species;"

    def test_set_level(self):
        cluster = stic.ClassCluster(self.sequences, centroid=self.seq1)
        cluster.set_level("family", "new_family")
        for seq in cluster.sequences:
            assert seq.header.taxonomy.family == "new_family"

    def test_closest_common_ancestor(self):
        cluster = stic.ClassCluster(self.sequences, centroid=self.seq1)
        common_ancestor = cluster.closest_common_ancestor
        assert common_ancestor.tax_str == "kingdom;phylum;class;order;family;genus;species"


class TestOrderCluster:
    
    def setup_method(self):
        self.header1 = ">seq1 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence1 = "ATCGATCGATCG"
        self.seq1 = stic.Sequence(self.header1, self.sequence1)

        self.header2 = ">seq2 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence2 = "GCTAGCTAGCTA"
        self.seq2 = stic.Sequence(self.header2, self.sequence2)

        self.header3 = ">seq3 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence3 = "CGTACGTACGTA"
        self.seq3 = stic.Sequence(self.header3, self.sequence3)

        self.header4 = ">seq4 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence4 = "TACGTACGTACG"
        self.seq4 = stic.Sequence(self.header4, self.sequence4)

        self.sequences = [self.seq1, self.seq2, self.seq3, self.seq4]

    def test_order_cluster_initialization(self):
        cluster = stic.OrderCluster(self.sequences, centroid=self.seq1)
        assert cluster.centroid == self.seq1
        assert len(cluster.sequences) == 4
        assert cluster.homogenous_tax_group

    def test_add_sequence(self):
        cluster = stic.OrderCluster(self.sequences, centroid=self.seq1)
        new_header = ">seq5 tax=kingdom;phylum;class;order;family;genus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        cluster.add_sequence(new_seq)
        assert len(cluster.sequences) == 5
        assert cluster.homogenous_tax_group
        new_header = ">seq5 tax=kingdom;phylum;class;Anotherorder;family;genus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        with pytest.raises(ValueError):
            cluster.add_sequence(new_seq)
        new_header = ">seq5 tax=Anotherkingdom;phylum;class;order;family;genus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        with pytest.raises(ValueError):
            cluster.add_sequence(new_seq)

    def test_set_level(self):
        cluster = stic.OrderCluster(self.sequences, centroid=self.seq1)
        cluster.set_level("family", "new_family")
        for seq in cluster.sequences:
            assert seq.header.taxonomy.family == "new_family"

    def test_closest_common_ancestor(self):
        cluster = stic.OrderCluster(self.sequences, centroid=self.seq1)
        common_ancestor = cluster.closest_common_ancestor
        assert common_ancestor.tax_str == "kingdom;phylum;class;order;family;genus;species"


class TestFamilyCluster:
    
    def setup_method(self):
        self.header1 = ">seq1 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence1 = "ATCGATCGATCG"
        self.seq1 = stic.Sequence(self.header1, self.sequence1)

        self.header2 = ">seq2 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence2 = "GCTAGCTAGCTA"
        self.seq2 = stic.Sequence(self.header2, self.sequence2)

        self.header3 = ">seq3 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence3 = "CGTACGTACGTA"
        self.seq3 = stic.Sequence(self.header3, self.sequence3)

        self.header4 = ">seq4 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence4 = "TACGTACGTACG"
        self.seq4 = stic.Sequence(self.header4, self.sequence4)

        self.sequences = [self.seq1, self.seq2, self.seq3, self.seq4]

    def test_family_cluster_initialization(self):
        cluster = stic.FamilyCluster(self.sequences, centroid=self.seq1)
        assert cluster.centroid == self.seq1
        assert len(cluster.sequences) == 4
        assert cluster.homogenous_tax_group

    def test_add_sequence(self):
        cluster = stic.FamilyCluster(self.sequences, centroid=self.seq1)
        new_header = ">seq5 tax=kingdom;phylum;class;order;family;genus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        cluster.add_sequence(new_seq)
        assert len(cluster.sequences) == 5
        assert cluster.homogenous_tax_group
        new_header = ">seq5 tax=kingdom;phylum;class;order;Anotherfamily;genus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        with pytest.raises(ValueError):
            cluster.add_sequence(new_seq)
        new_header = ">seq5 tax=Anotherkingdom;phylum;class;order;family;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        with pytest.raises(ValueError):
            cluster.add_sequence(new_seq)

    def test_set_level(self):
        cluster = stic.FamilyCluster(self.sequences, centroid=self.seq1)
        cluster.set_level("family", "new_family")
        for seq in cluster.sequences:
            assert seq.header.taxonomy.family == "new_family"

    def test_closest_common_ancestor(self):
        cluster = stic.FamilyCluster(self.sequences, centroid=self.seq1)
        common_ancestor = cluster.closest_common_ancestor
        assert common_ancestor.tax_str == "kingdom;phylum;class;order;family;genus;species"


class TestGenusCluster:
    
    def setup_method(self):
        self.header1 = ">seq1 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence1 = "ATCGATCGATCG"
        self.seq1 = stic.Sequence(self.header1, self.sequence1)

        self.header2 = ">seq2 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence2 = "GCTAGCTAGCTA"
        self.seq2 = stic.Sequence(self.header2, self.sequence2)

        self.header3 = ">seq3 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence3 = "CGTACGTACGTA"
        self.seq3 = stic.Sequence(self.header3, self.sequence3)

        self.header4 = ">seq4 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence4 = "TACGTACGTACG"
        self.seq4 = stic.Sequence(self.header4, self.sequence4)

        self.sequences = [self.seq1, self.seq2, self.seq3, self.seq4]

    def test_genus_cluster_initialization(self):
        cluster = stic.GenusCluster(self.sequences, centroid=self.seq1)
        assert cluster.centroid == self.seq1
        assert len(cluster.sequences) == 4
        assert cluster.homogenous_tax_group

    def test_add_sequence(self):
        cluster = stic.GenusCluster(self.sequences, centroid=self.seq1)
        new_header = ">seq5 tax=kingdom;phylum;class;order;family;genus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        cluster.add_sequence(new_seq)
        assert len(cluster.sequences) == 5
        assert cluster.homogenous_tax_group
        new_header = ">seq5 tax=kingdom;phylum;class;order;family;Anothergenus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        with pytest.raises(ValueError):
            cluster.add_sequence(new_seq)
        new_header = ">seq5 tax=Anotherkingdom;phylum;class;order;family;genus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        with pytest.raises(ValueError):
            cluster.add_sequence(new_seq)

    def test_set_level(self):
        cluster = stic.GenusCluster(self.sequences, centroid=self.seq1)
        cluster.set_level("family", "new_family")
        for seq in cluster.sequences:
            assert seq.header.taxonomy.family == "new_family"

    def test_closest_common_ancestor(self):
        cluster = stic.GenusCluster(self.sequences, centroid=self.seq1)
        common_ancestor = cluster.closest_common_ancestor
        assert common_ancestor.tax_str == "kingdom;phylum;class;order;family;genus;species"


class TestSpeciesCluster:
    
    def setup_method(self):
        self.header1 = ">seq1 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence1 = "ATCGATCGATCG"
        self.seq1 = stic.Sequence(self.header1, self.sequence1)

        self.header2 = ">seq2 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence2 = "GCTAGCTAGCTA"
        self.seq2 = stic.Sequence(self.header2, self.sequence2)

        self.header3 = ">seq3 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence3 = "CGTACGTACGTA"
        self.seq3 = stic.Sequence(self.header3, self.sequence3)

        self.header4 = ">seq4 tax=kingdom;phylum;class;order;family;genus;species;"
        self.sequence4 = "TACGTACGTACG"
        self.seq4 = stic.Sequence(self.header4, self.sequence4)

        self.sequences = [self.seq1, self.seq2, self.seq3, self.seq4]

    def test_species_cluster_initialization(self):
        cluster = stic.SpeciesCluster(self.sequences, centroid=self.seq1)
        assert cluster.centroid == self.seq1
        assert len(cluster.sequences) == 4
        assert cluster.homogenous_tax_group

    def test_add_sequence(self):
        cluster = stic.SpeciesCluster(self.sequences, centroid=self.seq1)
        new_header = ">seq5 tax=kingdom;phylum;class;order;family;genus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        cluster.add_sequence(new_seq)
        assert len(cluster.sequences) == 5
        assert cluster.homogenous_tax_group
        new_header = ">seq5 tax=kingdom;phylum;class;order;family;genus;Anotherspecies;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        with pytest.raises(ValueError):
            cluster.add_sequence(new_seq)
        new_header = ">seq5 tax=Anotherkingdom;phylum;class;order;family;genus;species;"
        new_sequence = "ATCGATCGATCG"
        new_seq = stic.Sequence(new_header, new_sequence)
        with pytest.raises(ValueError):
            cluster.add_sequence(new_seq)

    def test_set_level(self):
        cluster = stic.SpeciesCluster(self.sequences, centroid=self.seq1)
        cluster.set_level("family", "new_family")
        for seq in cluster.sequences:
            assert seq.header.taxonomy.family == "new_family"

    def test_closest_common_ancestor(self):
        cluster = stic.SpeciesCluster(self.sequences, centroid=self.seq1)
        common_ancestor = cluster.closest_common_ancestor
        assert common_ancestor.tax_str == "kingdom;phylum;class;order;family;genus;species"


class TestFastaFile:

    def setup_method(self):
        self.fasta_content = """>seq1 tax=kingdom;phylum;class;order;family;genus;species;
        ATCGATCGATCG
        >seq2 tax=kingdom;phylum;class;order;family;genus;species;
        GCTAGCTAGCTA
        >seq3 tax=kingdom;phylum;class;order;family;genus;species;
        CGTACGTACGTA
        >seq4 tax=kingdom;phylum;class;order;family;genus;species;
        TACGTACGTACG
        """
        self.fasta_file_path = pl.Path("test_fasta_file.fasta")
        with open(self.fasta_file_path, 'w') as f:
            f.write(self.fasta_content)
        self.fasta_file = stic.FastaFile(self.fasta_file_path)

    def teardown_method(self):
        self.fasta_file_path.unlink()

    def test_get_hash_table(self):
        hash_table = self.fasta_file.get_hash_table()
        assert isinstance(hash_table, dict)
        assert len(hash_table) == 4
        for key, value in hash_table.items():
            assert isinstance(key, int)
            assert isinstance(value, stic.SeqHeader)

    def test_get_seq_headers(self):
        headers = self.fasta_file.get_seq_headers()
        assert isinstance(headers, list)
        assert len(headers) == 4
        for header in headers:
            assert isinstance(header, stic.SeqHeader)

    def test_get_sequences(self):
        sequences = self.fasta_file.get_sequences()
        assert isinstance(sequences, list)
        assert len(sequences) == 4
        for seq in sequences:
            assert isinstance(seq, stic.Sequence)

    def test_filter_seq_by_hash(self):
        sequences = self.fasta_file.get_sequences()
        hash_list = [hash(seq) for seq in sequences[:2]]
        filtered_sequences = self.fasta_file.filter_seq_by_hash(hash_list)
        assert isinstance(filtered_sequences, list)
        assert len(filtered_sequences) == 2
        for seq in filtered_sequences:
            assert isinstance(seq, stic.Sequence)
            assert hash(seq) in hash_list

    def test_write_to_fasta_file(self, tmp_path):
        sequences = self.fasta_file.get_sequences()
        output_file_path = tmp_path.joinpath("output.fasta")
        self.fasta_file.write_to_fasta_file(output_file_path, sequences)
        with open(output_file_path, 'r') as f:
            lines = f.readlines()
        assert len(lines) == 8  # 4 sequences, each with a header and sequence line
        output_file_path.unlink()

    def test_get_seq_by_seq_id(self):
        seq_ids_list = ["seq1", "seq3"]
        sequences = self.fasta_file.get_seq_by_seq_id(seq_ids_list)
        assert isinstance(sequences, list)
        assert len(sequences) == 2
        for seq in sequences:
            assert isinstance(seq, stic.Sequence)
            assert str(seq.header.seq_id) in seq_ids_list


class TestTaxedFastaFile:

    def setup_method(self):
        self.fasta_content = """>seq1 tax=kingdom;phylum;class;order;family;genus;species;
        ATCGATCGATCG
        >seq2 tax=kingdom;phylum;class;order;family;genus;species;
        GCTAGCTAGCTA
        >seq3 tax=kingdom;phylum;class;order;family;genus;species;
        CGTACGTACGTA
        >seq4 tax=kingdom;phylum;class;order;family;genus;species;
        TACGTACGTACG
        """
        self.fasta_file_path = pl.Path("test_taxed_fasta_file.fasta")
        with open(self.fasta_file_path, 'w') as f:
            f.write(self.fasta_content)
        self.taxed_fasta_file = stic.TaxedFastaFile(self.fasta_file_path)

    def teardown_method(self):
        self.fasta_file_path.unlink()

    def test_get_tax_obj_list(self):
        tax_obj_list = self.taxed_fasta_file.tax_obj_list
        assert isinstance(tax_obj_list, list)
        assert len(tax_obj_list) == 4
        for tax in tax_obj_list:
            assert isinstance(tax, stic.Taxonomy)

    def test_filter_seq_by_tax(self):
        tax_str = "tax=kingdom;phylum;class;order;family;genus;species;"
        taxonomy = stic.Taxonomy(tax_str)
        sequences = self.taxed_fasta_file.filter_seq_by_tax(taxonomy)
        assert isinstance(sequences, list)
        assert all([seq.header.taxonomy == taxonomy for seq in sequences])
        assert all([isinstance(seq, stic.Sequence) for seq in sequences])
        assert len(sequences) == 4
        for seq in sequences:
            assert isinstance(seq, stic.Sequence)
            assert seq.header.taxonomy == taxonomy

    def test_get_tax_seq_map(self):
        tax_seq_map = self.taxed_fasta_file.get_tax_seq_map()
        assert isinstance(tax_seq_map, dict)
        assert len(tax_seq_map) == 1  # All taxa are the same
        for key, value in tax_seq_map.items():
            assert isinstance(key, stic.Taxonomy)
            assert isinstance(value, list)
            for seq in value:
                assert isinstance(seq, stic.Sequence)

    def test_filter_tax_set_at_last_known_level(self):
        taxonomies = self.taxed_fasta_file.filter_tax_set_at_last_known_level('species')
        assert isinstance(taxonomies, list)
        assert len(taxonomies) == 1  # All taxa are at the species level
        for tax in taxonomies:
            assert isinstance(tax, stic.Taxonomy)
            assert tax.last_known_level == 'species'


class TestTICUClust:

    def setup_method(self):
        self.usearch_bin = pl.Path('../bin/usearch11.0.667_i86linux64').absolute()
        self.uclust_work_dir = pl.Path('./Uclust-WD').absolute()
        self.uclust_work_dir.mkdir(exist_ok=True)
        self.ticuclust = stic.TICUClust(self.usearch_bin, self.uclust_work_dir)
        self.header1 = ">seq1 tax=kingdomA;phylumA;classA;orderA"
        self.sequence1 = "ATCGATCGATCG"
        self.seq1 = stic.Sequence(self.header1, self.sequence1)

        self.header2 = ">seq2 tax=kingdomA;phylumA;classA;orderA"
        self.sequence2 = "GCTAGCTAGCTA"
        self.seq2 = stic.Sequence(self.header2, self.sequence2)

        self.header3 = ">seq3 tax=kingdomA;phylumA;classA;orderA"
        self.sequence3 = "CGTACGTACGTA"
        self.seq3 = stic.Sequence(self.header3, self.sequence3)

        self.header4 = ">seq4 tax=kingdomA;phylumA;classA;orderA"
        self.sequence4 = "TACGTACGTACGAACTG"
        self.seq4 = stic.Sequence(self.header4, self.sequence4)

        self.sequences = [self.seq4, self.seq1, self.seq2, self.seq3]

    def test_run_uclust(self):
        with tempfile.TemporaryDirectory(dir=self.uclust_work_dir) as temp_cluster_dir:
            run_dir = pl.Path(temp_cluster_dir)
            input_fasta_path = run_dir.joinpath("input.fasta").absolute()
            sequence_cluster = stic.OrderCluster(self.sequences)
            sequence_cluster.write_to_fasta(input_fasta_path)
            file_pair_tuple = self.ticuclust.run_uclust(self.sequences, 0.987)
            in_fasta_file_path, centroid_seq_dict = file_pair_tuple
            assert in_fasta_file_path.exists()
            assert isinstance(centroid_seq_dict, dict)
            assert len(centroid_seq_dict) == 4
            # first element in the centroid_seq_dict should be self.seq4 as it
            # is the longest sequence
            assert list(centroid_seq_dict.keys())[0] == str(self.seq4.header)

    def test_parse_uc_file(self):
        uc_content = """\
            S\t0\t17\t*\t.\t*\t*\t*\tseq4 tax=kingdomA;phylumA;classA;orderA;NA-Family;NA-Genus;NA-Species\t*\n\
            S\t1\t12\t*\t.\t*\t*\t*\tseq3 tax=kingdomA;phylumA;classA;orderA;NA-Family;NA-Genus;NA-Species\t*\n\
            S\t2\t12\t*\t.\t*\t*\t*\tseq2 tax=kingdomA;phylumA;classA;orderA;NA-Family;NA-Genus;NA-Species\*\n\
            S\t3\t12\t*\t.\t*\t*\t*\tseq1 tax=kingdomA;phylumA;classA;orderA;NA-Family;NA-Genus;NA-Species\t*\n\
            C\t0\t1\t*\t*\t*\t*\t*\tseq4 tax=kingdomA;phylumA;classA;orderA;NA-Family;NA-Genus;NA-Species\t*\n\
            C\t1\t1\t*\t*\t*\t*\t*\tseq3 tax=kingdomA;phylumA;classA;orderA;NA-Family;NA-Genus;NA-Species\t*\n\
            C\t2\t1\t*\t*\t*\t*\t*\tseq2 tax=kingdomA;phylumA;classA;orderA;NA-Family;NA-Genus;NA-Species\t*\n\
            C\t3\t1\t*\t*\t*\t*\t*\tseq1 tax=kingdomA;phylumA;classA;orderA;NA-Family;NA-Genus;NA-Species\t*\n\
        """
        with tempfile.NamedTemporaryFile(delete=False) as uc_file:
            uc_file.write(uc_content.encode('utf-8'))
            uc_file_path = uc_file.name
        uc_dict = self.ticuclust.parse_uc_file(uc_file_path)
        assert isinstance(uc_dict, dict)
        assert len(uc_dict) == 4 # 4 sequences and each is one cluster
        assert list(uc_dict.values())[0][0] == ">seq4 tax=kingdomA;phylumA;classA;orderA;NA-Family;NA-Genus;NA-Species"
        uc_dict = self.ticuclust.parse_uc_file(uc_file_path, cut_tax=True)
        assert isinstance(uc_dict, dict)
        assert len(uc_dict) == 4
        assert list(uc_dict.values())[0][0] == ">seq4"
        pl.Path(uc_file_path).unlink()

    # TODO: check the rest of the methods
    def test_get_sequences_clusters(self):
        clusters = self.ticuclust.get_sequences_clusters(self.sequences, 0.987)
        assert isinstance(clusters, list)
        assert len(clusters) == 4 # 4 sequences, each in its own cluster
        assert isinstance(clusters[0], stic.SequenceCluster)

    def test_sort_seqs(self):
        with tempfile.TemporaryDirectory(dir=self.uclust_work_dir) as temp_cluster_dir:
            run_dir = pl.Path(temp_cluster_dir)
            input_fasta_path = run_dir.joinpath("input.fasta").absolute()
            sequence_cluster = stic.SequenceCluster(self.sequences, force_homogeneity=False)
            sequence_cluster.write_to_fasta(input_fasta_path)
            sorted_fasta_file = self.ticuclust.sort_seqs(str(input_fasta_path), by="length")
            assert sorted_fasta_file.exists()
            assert sorted_fasta_file.stat().st_size > 0
            # first sequence in the sorted file should be self.seq4 as it is the longest sequence
            with open(sorted_fasta_file, 'r') as f:
                lines = f.readlines()
                assert lines[0].strip() == str(self.seq4.header)
            sorted_fasta_file.unlink()


class TestTICAnalysis:

    def setup_method(self):
        self.fasta_content = """>seq1 tax=Bacteria;
        ATCGATCGATCG
        >seq2 tax=Bacteria;phylum;
        GCTAGCTAGCTA
        >seq3 tax=Bacteria;phylum;class;
        CGTACGTACGTA
        >seq4 tax=Bacteria;phylum;class;order;
        TACGTGGVACGTACG
        >seq5 tax=Bacteria;phylum;class;order;family;
        TACGTAAACGTACG
        >seq6 tax=Bacteria;phylum;class;order;family;genus;
        TACGTACCCGTACG
        >seq7 tax=Bacteria;phylum;class;order;family;genus;species;
        TACGTACGTTTTACG
        """
        self.fasta_file_path = pl.Path("test_tic_analysis.fasta")
        with open(self.fasta_file_path, 'w') as f:
            f.write(self.fasta_content)
        self.tic_analysis = stic.TICAnalysis(self.fasta_file_path)

    def teardown_method(self):
        self.fasta_file_path.unlink()

    def test_filter_tax_set_at_last_known_level(self):
        taxonomies = self.tic_analysis.filter_tax_set_at_last_known_level('species')
        assert isinstance(taxonomies, list)
        assert len(taxonomies) == 1  # All taxa are at the species level
        for tax in taxonomies:
            assert isinstance(tax, stic.Taxonomy)
            assert tax.last_known_level == 'species'

    def test_fill_upto_order(self):
        output_fasta_path = self.tic_analysis.fill_upto_order()
        assert output_fasta_path.exists()
        print(f"output_fasta_path: {output_fasta_path}")
        with open(output_fasta_path, 'r') as f:
            lines = f.readlines()
        # only 3 sequences have incomplete taxonomy upto the order level
        assert len(lines) == 6
        output_fasta_path.unlink()

    # TODO test the rest of methods
    def test_grow_taxonomy(self):
        tax_str = "tax=kingdom;phylum;class;order;family;genus;species;"
        taxonomy = stic.Taxonomy(tax_str)
        clusters = self.tic_analysis.grow_taxonomy(taxonomy, cluster_threshold=0.987)
        assert isinstance(clusters, list)
        for cluster in clusters:
            assert isinstance(cluster, stic.SequenceCluster)

    def test_complete_family_level(self):
        all_known_order_fasta_path = self.tic_analysis.fill_up_to_order()
        output_fasta_path = self.tic_analysis.complete_family_level(all_known_order_fasta_path)
        assert output_fasta_path.exists()
        with open(output_fasta_path, 'r') as f:
            lines = f.readlines()
        assert len(lines) == 8  # 4 sequences, each with a header and sequence line
        output_fasta_path.unlink()

    def test_complete_genus_level(self):
        all_known_order_fasta_path = self.tic_analysis.fill_up_to_order()
        all_known_family_fasta_path = self.tic_analysis.complete_family_level(all_known_order_fasta_path)
        output_fasta_path = self.tic_analysis.complete_genus_level(all_known_family_fasta_path)
        assert output_fasta_path.exists()
        with open(output_fasta_path, 'r') as f:
            lines = f.readlines()
        assert len(lines) == 8  # 4 sequences, each with a header and sequence line
        output_fasta_path.unlink()

    def test_complete_species_level(self):
        all_known_order_fasta_path = self.tic_analysis.fill_up_to_order()
        all_known_family_fasta_path = self.tic_analysis.complete_family_level(all_known_order_fasta_path)
        all_known_genus_fasta_path = self.tic_analysis.complete_genus_level(all_known_family_fasta_path)
        output_fasta_path = self.tic_analysis.complete_species_level(all_known_genus_fasta_path)
        assert output_fasta_path.exists()
        with open(output_fasta_path, 'r') as f:
            lines = f.readlines()
        assert len(lines) == 8  # 4 sequences, each with a header and sequence line
        output_fasta_path.unlink()
    
    def test_run(self):
        output_fasta_path = self.tic_analysis.run()
        assert output_fasta_path.exists()
        with open(output_fasta_path, 'r') as f:
            lines = f.readlines()
            print(lines)
        assert len(lines) == 14  # 3 sequence
        output_fasta_path.unlink()


if __name__ == '__main__':
    pytest.main()