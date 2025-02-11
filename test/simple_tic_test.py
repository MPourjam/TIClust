import unittest
import sys
import os
import pytest
import re
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from tic import simple_tic as stic
import logging
import tempfile
from collections import namedtuple as nt
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
    
    def test_taxonomy_init_unclassifed(self):
        tax_str = "tax=Unclassified"
        taxonomy = stic.Taxonomy(tax_str)
        assert taxonomy.tax_str == ""
        assert taxonomy.kingdom == "NA-Kingdom"
        assert taxonomy.phylum == "NA-Phylum"
        assert taxonomy.class_ == "NA-Class"
        assert taxonomy.order == "NA-Order"
        assert taxonomy.family == "NA-Family"
        assert taxonomy.genus == "NA-Genus"
        assert taxonomy.species == "NA-Species"
        assert taxonomy.full_tax == "NA-Kingdom;NA-Phylum;NA-Class;NA-Order;NA-Family;NA-Genus;NA-Species"
        assert taxonomy.last_known_level == ""
        assert taxonomy.orig_tax == "Unclassified"
        assert str(taxonomy) == "tax=NA-Kingdom;NA-Phylum;NA-Class;NA-Order;NA-Family;NA-Genus;NA-Species"
        assert taxonomy.get_level("kingdom") == "NA-Kingdom"

    def test_taxonomy_all_invalid(self):
        tax_str = "tax=;;;;;;"
        taxonomy = stic.Taxonomy(tax_str)
        assert taxonomy.tax_str == ""
        assert taxonomy.kingdom == "NA-Kingdom"
        assert taxonomy.phylum == "NA-Phylum"
        assert taxonomy.class_ == "NA-Class"
        assert taxonomy.order == "NA-Order"
        assert taxonomy.family == "NA-Family"
        assert taxonomy.genus == "NA-Genus"
        assert taxonomy.species == "NA-Species"
        assert taxonomy.full_tax == "NA-Kingdom;NA-Phylum;NA-Class;NA-Order;NA-Family;NA-Genus;NA-Species"
        assert taxonomy.last_known_level == ""
        assert taxonomy.orig_tax == ";;;;;;"
        assert str(taxonomy) == "tax=NA-Kingdom;NA-Phylum;NA-Class;NA-Order;NA-Family;NA-Genus;NA-Species"
        assert taxonomy.tax_list == ['']
        assert taxonomy.full_tax_list == [
            'NA-Kingdom',
            'NA-Phylum',
            'NA-Class',
            'NA-Order',
            'NA-Family',
            'NA-Genus',
            'NA-Species'
        ]
        assert taxonomy == stic.Taxonomy("tax=Unclassified;;;;;")
        assert bool(taxonomy) == False
    
    def test_taxonomy_partial_invalid(self):
        tax_str = "tax=Bacteria;Proteobacteria;Alphaproteobacteria;SAR11 clade;Clade III"
        taxonomy = stic.Taxonomy(tax_str)
        assert taxonomy.tax_str == "Bacteria;Proteobacteria;Alphaproteobacteria;SAR11 clade;Clade III"
        assert taxonomy.kingdom == "Bacteria"
        assert taxonomy.phylum == "Proteobacteria"
        assert taxonomy.class_ == "Alphaproteobacteria"
        assert taxonomy.order == "SAR11 clade"
        assert taxonomy.family == "Clade III"
        assert taxonomy.genus == "NA-Genus"
        assert taxonomy.species == "NA-Species"
        assert taxonomy.full_tax == "Bacteria;Proteobacteria;Alphaproteobacteria;SAR11 clade;Clade III;NA-Genus;NA-Species"
        assert taxonomy.last_known_level == "family"
        assert taxonomy.orig_tax == "Bacteria;Proteobacteria;Alphaproteobacteria;SAR11 clade;Clade III"
        assert str(taxonomy) == "tax=Bacteria;Proteobacteria;Alphaproteobacteria;SAR11 clade;Clade III;NA-Genus;NA-Species"
        assert taxonomy.is_known_upto("family")


    def test_taxonomy_eukaryotic(self):
        tax_str = "tax=Eukaryota;Amorphea;Obazoa;Opisthokonta;"\
            "Nucletmycea;Fungi;Dikarya;Basidiomycota;Pucciniomycotina;"\
            "Pucciniomycetes;Pucciniales;Pucciniaceae;Gymnosporangium"
        taxonomy = stic.Taxonomy(tax_str)
        assert taxonomy.tax_str == "Eukaryota;Amorphea;Obazoa;Opisthokonta;Nucletmycea;Fungi;Dikarya"
        assert taxonomy.kingdom == "Eukaryota"
        assert taxonomy.phylum == "Amorphea"
        assert taxonomy.class_ == "Obazoa"
        assert taxonomy.order == "Opisthokonta"
        assert taxonomy.family == "Nucletmycea"
        assert taxonomy.genus == "Fungi"
        assert taxonomy.species == "Dikarya"
        assert taxonomy.full_tax == "Eukaryota;Amorphea;Obazoa;Opisthokonta;Nucletmycea;Fungi;Dikarya"
        assert taxonomy.last_known_level == "species"
        assert taxonomy.get_level("kingdom") == "Eukaryota"

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
        assert seq_id.id_str == "seq1"
        header = ">seq1 ; some; description"
        seq_id = stic.SeqID(header)
        assert seq_id.id_str == "seq1"


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
        self.ticuclust = stic.TICUClust(self.uclust_work_dir)
        self.test_fasta_file_path = pl.Path("All_known_order_500_sequences.fasta")
        self.test_fasta_file = stic.TaxedFastaFile(self.test_fasta_file_path)
        self.sequences = self.test_fasta_file.get_sequences()

    @pytest.fixture
    def fix_run_uclust(self):
        # get a path a temporary directory
        tmp_dir = tempfile.TemporaryDirectory(dir=self.uclust_work_dir)
        run_dir = pl.Path(tmp_dir.name)
        input_fasta_path = run_dir.joinpath("input.fasta").absolute()
        sequence_cluster = stic.OrderCluster(self.sequences)
        sequence_cluster.write_to_fasta(input_fasta_path)
        uclust_obj = stic.TICUClust(self.uclust_work_dir)
        yield uclust_obj
        shutil.rmtree(tmp_dir)

    def test_run_uclust(self, fix_run_uclust):
        uclust_obj = fix_run_uclust
        file_pair_tuple = uclust_obj.run_uclust(self.sequences, 0.987)
        in_fasta_file_path, centroid_seq_dict = file_pair_tuple
        assert in_fasta_file_path.exists()
        assert isinstance(centroid_seq_dict, dict)
        assert len(centroid_seq_dict) == 375
        # first element in the centroid_seq_dict should be self.seq4 as it
        # is the longest sequence
        # count of clusters seqs should equal to the total number of sequences
        clusters_sequs = [seq for cent, clust in centroid_seq_dict.items() for seq in clust]
        assert len(clusters_sequs) == len(self.sequences)

    def test_parse_uc_file(self):
        uc_content = """S\t0\t17\t*\t.\t*\t*\t*\tseq4 tax=kingdomA;phylumA;classA;orderA;NA-Family;NA-Genus;NA-Species\t*
S\t1\t12\t*\t.\t*\t*\t*\tseq3 tax=kingdomA;phylumA;classA;orderA;NA-Family;NA-Genus;NA-Species\t*
S\t2\t12\t*\t.\t*\t*\t*\tseq2 tax=kingdomA;phylumA;classA;orderA;NA-Family;NA-Genus;NA-Species\t*
S\t3\t12\t*\t.\t*\t*\t*\tseq1 tax=kingdomA;phylumA;classA;orderA;NA-Family;NA-Genus;NA-Species\t*
C\t0\t1\t*\t*\t*\t*\t*\tseq4 tax=kingdomA;phylumA;classA;orderA;NA-Family;NA-Genus;NA-Species\t*
C\t1\t1\t*\t*\t*\t*\t*\tseq3 tax=kingdomA;phylumA;classA;orderA;NA-Family;NA-Genus;NA-Species\t*
C\t2\t1\t*\t*\t*\t*\t*\tseq2 tax=kingdomA;phylumA;classA;orderA;NA-Family;NA-Genus;NA-Species\t*
C\t3\t1\t*\t*\t*\t*\t*\tseq1 tax=kingdomA;phylumA;classA;orderA;NA-Family;NA-Genus;NA-Species\t*
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


seq_count_stats = nt(
    "seq_count_stats",
    [
        "total_seq_count",
        "bac_seq_count",
        "non_bac_seq_count",
        "new_fam_seq_count",
        "new_gen_seq_count",
        "new_spe_seq_count"
    ]
)

cluster_params = nt(
    "cluster_params",
    [
        "fam_thr",
        "gen_thr",
        "spe_thr"
    ]
)

fafi_stats = nt(
    "fafi_stats",
    [
        "fasta_file_path",
        "cluster_params",
        "seq_count_stats"
    ]
)

INPUT_FASTAs = [
    fafi_stats(
        fasta_file_path="Minor_pseudo_sequences.fasta",
        cluster_params=cluster_params(
            fam_thr=0.9,
            gen_thr=0.95,
            spe_thr=0.987
        ),
        seq_count_stats=seq_count_stats(
            total_seq_count=7,
            bac_seq_count=7,
            non_bac_seq_count=0,
            new_fam_seq_count=5,
            new_gen_seq_count=6,
            new_spe_seq_count=6
        )
    ),
    fafi_stats(
        fasta_file_path="All_known_order_500_sequences.fasta",
        cluster_params=cluster_params(
            fam_thr=0.9,
            gen_thr=0.95,
            spe_thr=0.987
        ),
        seq_count_stats=seq_count_stats(
            total_seq_count=500,
            bac_seq_count=500,
            non_bac_seq_count=0,
            new_fam_seq_count=33,  # simple tic logs: Completed genus level for 33 families
            new_gen_seq_count=123,  # TODO simple tic logs: Completed species level for 130 genera. But counting results in 123
            new_spe_seq_count=408
        )
    ),
    fafi_stats(
        fasta_file_path="Mixed_kingdoms_1K_V3-V4_sequences.fasta",
        cluster_params=cluster_params(
            fam_thr=0.9,
            gen_thr=0.95,
            spe_thr=0.987
        ),
        seq_count_stats=seq_count_stats(
            total_seq_count=1000,
            bac_seq_count=1000,
            non_bac_seq_count=0,
            new_fam_seq_count=60,  # TODO simple tic logs: Completed genus level for 123 families. But counting results in 60
            new_gen_seq_count=343,  # TODO simple tic logs: Completed species level for 340 genera. But counting results in 343
            new_spe_seq_count=564
        )
    )
]


@pytest.mark.parametrize(
    "in_fasta",
    INPUT_FASTAs
)
class TestTICAnalysis:

    @pytest.fixture(autouse=True)
    def setup_method(self, in_fasta):
        self.fasta_file_path = pl.Path(in_fasta.fasta_file_path).resolve()
        self.tic_analysis = stic.TICAnalysis(self.fasta_file_path)
        yield
        self.tic_analysis.cleanup(full=True)

    @pytest.fixture
    def fix_run(self):
        self.tic_analysis.run()
        yield self.tic_analysis
        self.tic_analysis.cleanup(full=True)

    @pytest.fixture
    def fix_parent_child_pairs(self):
        tax_levels = stic.Taxonomy.level_tax_map
        parent_child_pairs = []
        for i in range(len(tax_levels) - 2):
            for j in range(len(tax_levels) - 1, i, -1):
                parent_child_pairs.append((tax_levels[i], tax_levels[j]))
        return parent_child_pairs

    def test_preserve_counts(self, fix_run):
        tic_analysis = fix_run
        output_fasta_path = tic_analysis.tic_output_fasta_path
        with open(output_fasta_path, 'r') as f:
            lines = f.readlines()
        # sequence ids of input and output fasta files should be the same
        output_fasta = stic.TaxedFastaFile(output_fasta_path)
        bac_input_seqs = [] + tic_analysis.filter_bac_seq_last_kown_at("kingdom", flatten=True)
        bac_input_seqs += tic_analysis.filter_bac_seq_last_kown_at("phylum", flatten=True)
        bac_input_seqs += tic_analysis.filter_bac_seq_last_kown_at("class", flatten=True)
        bac_input_seqs += tic_analysis.filter_bac_seq_last_kown_at("order", flatten=True)
        bac_input_seqs += tic_analysis.filter_bac_seq_last_kown_at("family", flatten=True)
        bac_input_seqs += tic_analysis.filter_bac_seq_last_kown_at("genus", flatten=True)
        bac_input_seqs += tic_analysis.filter_bac_seq_last_kown_at("species", flatten=True)
        bac_input_seq_ids = [seq.header.seq_id for seq in bac_input_seqs]
        output_seq_ids = output_fasta.get_seq_ids()
        for seq_id in bac_input_seq_ids:
            assert seq_id in output_seq_ids
        assert len(bac_input_seq_ids) == len(output_seq_ids)
        # all sequences in the output file should have full taxonomy
        for taxa in output_fasta.tax_obj_list:
            assert taxa.last_known_level == "species"
        all_input_seqs = self.tic_analysis.fasta_file.get_seq_ids()
        non_bac_fasta = stic.TaxedFastaFile(self.tic_analysis.non_bact_fasta_path)
        non_bac_seq_ids = non_bac_fasta.get_seq_ids()
        tic_seq_processed_count = len(output_seq_ids) + len(non_bac_seq_ids)
        init_seq_count = len(all_input_seqs)
        assert tic_seq_processed_count == init_seq_count

    def test_append_non_bacteria_seqs(self, fix_run):
        tic_analysis = fix_run
        non_bact_fasta_path = tic_analysis.non_bact_fasta_path
        bact_counts = 0
        non_bact_counts = 0
        with open(tic_analysis.fasta_file.fasta_file_path, 'r') as f:
            line = f.readline().strip()
            while line:
                if line.startswith(">"):
                    if re.search(r"tax=Bacteria;", line, re.IGNORECASE):
                        bact_counts += 1
                    else:
                        non_bact_counts += 1
                line = f.readline().strip()
        # assert non_bact_fasta_path.exists()
        output_non_bact_counts = 0
        with open(non_bact_fasta_path, 'r') as f:
            line = f.readline().strip()
            while line:
                if line.startswith(">"):
                    output_non_bact_counts += 1
                line = f.readline().strip()
        assert non_bact_counts == output_non_bact_counts

    def test_taxonomy_completeness(self, fix_run):
        tic_analysis = fix_run
        output_fasta_path = tic_analysis.tic_output_fasta_path
        output_fasta = stic.TaxedFastaFile(output_fasta_path)
        tax_obj_list = output_fasta.tax_obj_list
        for tax in tax_obj_list:
            assert tax.last_known_level == "species"
        assert all([tax.last_known_level == "species" for tax in tax_obj_list])

    def test_taxonomy_consistency(self, fix_run, fix_parent_child_pairs):
        tic_analysis = fix_run
        output_fasta_path = tic_analysis.tic_output_fasta_path
        output_fasta = stic.TaxedFastaFile(output_fasta_path)
        tax_obj_list = output_fasta.tax_obj_list
        #
        phyl_king_d = {}
        for tax in tax_obj_list:
            for parent, child in fix_parent_child_pairs:
                phyl_king_d.get(tax.get_level(child), []).append(tax.get_level(parent))

        for child, parents in phyl_king_d.items():
            if len(parents) > 1:
                logging.warning(f"Taxonomy {child} has multiple parents: {parents}")
        assert all([len(parents) == 1 for parents in phyl_king_d.values()])

    # TODO merge following two test classes into this class. Use stats INPUT_FASTAs to check the results cluster counts.


@pytest.mark.minor_test
class TestTICAnalysisSmall:

    def setup_method(self):
        self.fasta_file_path = pl.Path(__file__).parent.joinpath("Minor_pseudo_sequences.fasta").resolve()
        self.tic_analysis = stic.TICAnalysis(self.fasta_file_path)

    def teardown_method(self):
        self.tic_analysis.cleanup(full=True)
    
    @pytest.mark.parametrize("level", ["phylum", "class", "order", "family", "genus", "species"])
    def test_filter_tax_set_at_last_known_level(self, level):
        taxonomies = self.tic_analysis.filter_tax_set_at_last_known_level(level)
        assert isinstance(taxonomies, list)
        assert len(taxonomies) == 1  # All taxa are at the species level
        for tax in taxonomies:
            assert isinstance(tax, stic.Taxonomy)
            assert tax.last_known_level == level

    @pytest.fixture
    def fix_fill_upto_order(self):
        all_known_order_fasta = self.tic_analysis.fill_upto_order()
        yield all_known_order_fasta
        all_known_order_fasta.unlink()

    @pytest.fixture
    def fix_complete_family_level(self, fix_fill_upto_order):
        all_known_family_fasta = self.tic_analysis.complete_family_level(fix_fill_upto_order)
        yield all_known_family_fasta
        all_known_family_fasta.unlink()

    @pytest.fixture
    def fix_complete_genus_level(self, fix_complete_family_level):
        all_known_genus_fasta = self.tic_analysis.complete_genus_level(fix_complete_family_level)
        yield all_known_genus_fasta
        all_known_genus_fasta.unlink()

    @pytest.fixture
    def fix_complete_species_level(self, fix_complete_genus_level):
        all_known_species_fasta = self.tic_analysis.complete_species_level(fix_complete_genus_level)
        yield all_known_species_fasta
        all_known_species_fasta.unlink()

    @pytest.fixture
    def fix_run(self):
        self.tic_analysis.run()
        yield self.tic_analysis
        self.tic_analysis.cleanup(full=True)

    def test_fill_upto_order(self, fix_fill_upto_order):
        with open(fix_fill_upto_order, 'r') as f:
            lines = f.readlines()
        # only 3 sequences have incomplete taxonomy upto the order level + 1 sequence with known order level
        assert len(lines) == 8

    def test_complete_family_level(self, fix_complete_family_level):
        with open(fix_complete_family_level, 'r') as f:
            lines = f.readlines()
        # 4 sequences with known order level but unknown family level + 1 sequence with known family level
        assert len(lines) == 10

    def test_complete_genus_level(self, fix_complete_genus_level):
        with open(fix_complete_genus_level, 'r') as f:
            lines = f.readlines()
        # 5 sequences with known family level but unknown genus level + 1 sequence with known genus level
        assert len(lines) == 12

    def test_complete_species_level(self, fix_complete_species_level):
        with open(fix_complete_species_level, 'r') as f:
            lines = f.readlines()
        # 6 sequences with known genus level but unknown species level + 1 sequence with known species level
        assert len(lines) == 14

    def test_run(self, fix_run):
        tic_analysis = fix_run
        output_fasta_path = tic_analysis.tic_output_fasta_path
        with open(output_fasta_path, 'r') as f:
            lines = f.readlines()
        # 7 sequences with known species level
        assert len(lines) == 14
        # sequence ids of input and output fasta files should be the same
        output_fasta = stic.TaxedFastaFile(output_fasta_path)
        input_seq_ids = tic_analysis.fasta_file.get_seq_ids()
        output_seq_ids = output_fasta.get_seq_ids()
        for seq_id in input_seq_ids:
            assert seq_id in output_seq_ids
        assert len(input_seq_ids) == len(output_seq_ids)
        # all sequences in the output file should have full taxonomy
        for taxa in output_fasta.tax_obj_list:
            assert taxa.last_known_level == "species"


@pytest.mark.medium_test
class TestTICAnalysisAllKnownOrder:
    REPEAT_PARAMETER = 1

    def setup_method(self):
        self.fasta_file_path = pl.Path(__file__).parent.joinpath("All_known_order_500_sequences.fasta").resolve()
        self.tic_analysis = stic.TICAnalysis(self.fasta_file_path)

    def teardown_method(self):
        self.tic_analysis.cleanup(full=True)

    @pytest.mark.parametrize("level", ["phylum", "class", "order", "family", "genus", "species"])
    def test_filter_tax_set_at_last_known_level(self, level):
        taxonomies = self.tic_analysis.filter_tax_set_at_last_known_level(level)
        assert isinstance(taxonomies, list)
        if level == "order":
            print(taxonomies)
            assert len(taxonomies) == 1  # All taxa are at the species level
        for tax in taxonomies:
            assert isinstance(tax, stic.Taxonomy)
            assert tax.last_known_level == level

    @pytest.fixture
    def fix_fill_upto_order(self):
        all_known_order_fasta = self.tic_analysis.fill_upto_order()
        yield all_known_order_fasta
        all_known_order_fasta.unlink()

    @pytest.fixture
    def fix_complete_family_level(self, fix_fill_upto_order):
        all_known_family_fasta = self.tic_analysis.complete_family_level(fix_fill_upto_order)
        yield all_known_family_fasta
        all_known_family_fasta.unlink()

    @pytest.fixture
    def fix_complete_genus_level(self, fix_complete_family_level):
        all_known_genus_fasta = self.tic_analysis.complete_genus_level(fix_complete_family_level)
        yield all_known_genus_fasta
        all_known_genus_fasta.unlink()

    @pytest.fixture
    def fix_complete_species_level(self, fix_complete_genus_level):
        all_known_species_fasta = self.tic_analysis.complete_species_level(fix_complete_genus_level)
        yield all_known_species_fasta
        # all_known_species_fasta.unlink()

    @pytest.fixture
    def fix_run(self):
        self.tic_analysis.run()
        yield self.tic_analysis
        self.tic_analysis.cleanup(full=True)

    def test_fill_upto_order(self, fix_fill_upto_order):
        with open(fix_fill_upto_order, 'r') as f:
            lines = f.readlines()
        assert len(lines) == 1000
        # sequence ids of input and output fasta files should be the same
        # Because all sequences input file have known taxonomy upto the order level
        # then all sequences in the output file of this should have the same taxonomy
        output_fasta = stic.TaxedFastaFile(fix_fill_upto_order)
        input_seq_ids = self.tic_analysis.fasta_file.get_seq_ids()
        output_seq_ids = output_fasta.get_seq_ids()
        for seq_id in input_seq_ids:
            assert seq_id in output_seq_ids
        assert len(input_seq_ids) == len(output_seq_ids)

    @pytest.mark.parametrize("i", range(REPEAT_PARAMETER))
    def test_complete_family_level(self, i, fix_complete_family_level):
        logging.info(f"Attempt {i + 1} of complete_family_level")
        with open(fix_complete_family_level, 'r') as f:
            lines = f.readlines()
        assert len(lines) == 1000
        # There should be 33 distinct FOTU in the output file
        output_fasta = stic.TaxedFastaFile(fix_complete_family_level)
        tax_obj_list = output_fasta.tax_obj_list
        family_set = set()
        for tax in tax_obj_list:
            family_set.add(tax.family)
        assert len(family_set) == 33

    @pytest.mark.parametrize("i", range(REPEAT_PARAMETER))
    def test_complete_genus_level(self, i, fix_complete_genus_level):
        logging.info(f"Attempt {i + 1} of complete_genus_level")
        with open(fix_complete_genus_level, 'r') as f:
            lines = f.readlines()
        assert len(lines) == 1000
        # There should be 130 distinct GOTU in the output file
        output_fasta = stic.TaxedFastaFile(fix_complete_genus_level)
        tax_obj_list = output_fasta.tax_obj_list
        genus_set = set()
        for tax in tax_obj_list:
            genus_set.add(tax.genus)
        assert len(genus_set) == 130

    @pytest.mark.parametrize("i", range(REPEAT_PARAMETER))
    def test_complete_species_level(self, i, fix_complete_species_level):
        logging.info(f"Attempt {i + 1} of complete_species_level")
        with open(fix_complete_species_level, 'r') as f:
            lines = f.readlines()
        assert len(lines) == 1000
        # There should be 1000 distinct SOTU in the output file
        output_fasta = stic.TaxedFastaFile(fix_complete_species_level)
        tax_obj_list = output_fasta.tax_obj_list
        species_set = set()
        for tax in tax_obj_list:
            species_set.add(tax.species)
        assert len(species_set) == 412


@pytest.mark.mixed_kingdoms_test
class TestTICAnalysisMixedKingdoms:

    def setup_method(self):
        self.fasta_file_path = pl.Path(__file__).parent.joinpath("Mixed_kingdoms_1K_V3-V4_sequences.fasta").resolve()
        self.tic_analysis = stic.TICAnalysis(self.fasta_file_path)

    def teardown_method(self):
        self.tic_analysis.cleanup(full=True)

    # define a fixture passable to other tests in this class
    @pytest.fixture
    def fix_fill_upto_order(self):
        all_known_order_fasta = self.tic_analysis.fill_upto_order()
        yield all_known_order_fasta
        all_known_order_fasta.unlink()

    @pytest.fixture
    def fix_complete_family_level(self, fix_fill_upto_order):
        all_knwon_family_file_path = self.tic_analysis.complete_family_level(fix_fill_upto_order)
        yield all_knwon_family_file_path
        all_knwon_family_file_path.unlink()

    @pytest.fixture
    def fix_complete_genus_level(self, fix_complete_family_level):
        all_known_genus_fasta_path = self.tic_analysis.complete_genus_level(fix_complete_family_level)
        yield all_known_genus_fasta_path
        all_known_genus_fasta_path.unlink()

    @pytest.fixture
    def fix_complete_species_level(self, fix_complete_genus_level):
        all_known_species_fasta_path = self.tic_analysis.complete_species_level(fix_complete_genus_level)
        yield all_known_species_fasta_path
        all_known_species_fasta_path.unlink()

    @pytest.fixture
    def fix_run(self):
        output_fasta_path = self.tic_analysis.run()
        yield self.tic_analysis
        self.tic_analysis.cleanup(full=True)

    def test_fill_upto_order(self, fix_fill_upto_order):
        with open(fix_fill_upto_order, 'r') as f:
            lines = f.readlines()
        assert len(lines) == 212

    def test_complete_family_level(self, fix_complete_family_level):
        with open(fix_complete_family_level, 'r') as f:
            lines = f.readlines()
        assert len(lines) == 522

    def test_complete_genus_level(self, fix_complete_genus_level):
        with open(fix_complete_genus_level, 'r') as f:
            lines = f.readlines()
        assert len(lines) == 1314

    def test_complete_species_level(self, fix_complete_species_level):
        with open(fix_complete_species_level, 'r') as f:
            lines = f.readlines()
        assert len(lines) == 1314


if __name__ == '__main__':
    pytest.main()
