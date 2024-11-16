import unittest
from simple_tic import Taxonomy, SeqID, SeqHeader, Sequence, ZOTUFASTA

class TestTaxonomy(unittest.TestCase):

    def test_taxonomy_initialization(self):
        tax_str = "tax=kingdom;phylum;class;order;family;genus;species;"
        taxonomy = Taxonomy(tax_str)
        self.assertEqual(taxonomy.tax_str, "kingdom;phylum;class;order;family;genus;species;")
        self.assertEqual(taxonomy.kingdom, "kingdom")
        self.assertEqual(taxonomy.phylum, "phylum")
        self.assertEqual(taxonomy.class_, "class")
        self.assertEqual(taxonomy.order, "order")
        self.assertEqual(taxonomy.family, "family")
        self.assertEqual(taxonomy.genus, "genus")
        self.assertEqual(taxonomy.species, "species")

    def test_get_tax_upto(self):
        tax_str = "tax=kingdom;phylum;class;order;family;genus;species;"
        taxonomy = Taxonomy(tax_str)
        self.assertEqual(taxonomy.get_tax_upto('family'), "kingdom___phylum___class___order___family")

    def test_clean_taxonomy(self):
        tax_str = "tax=kingdom;phylum;unclass;unorder;family;genus;species;"
        taxonomy = Taxonomy(tax_str)
        self.assertEqual(taxonomy.tax_list, ["kingdom", "phylum"])
        self.assertEqual(str(taxonomy), "kingdom___phylum")
        self.assertEqual(taxonomy.get_tax_upto('family'), "kingdom___phylum")
        tax_str = "tax=kingdom;phylum;class;order;;genus;species;"
        taxonomy = Taxonomy(tax_str)
        self.assertEqual(taxonomy.tax_list, ["kingdom", "phylum", "class", "order"])
        self.assertEqual(str(taxonomy), "kingdom___phylum___class___order")


class TestSeqID(unittest.TestCase):

    def test_seq_id_initialization(self):
        header = ">seq1 some description"
        seq_id = SeqID(header)
        self.assertEqual(seq_id.head_id, ">seq1")


class TestSeqHeader(unittest.TestCase):

    def test_seq_header_initialization(self):
        header = ">seq1 tax=kingdom;phylum;class;order;family;genus;species;"
        seq_header = SeqHeader(header)
        self.assertEqual(str(seq_header.seq_id), ">seq1")
        self.assertEqual(str(seq_header.taxonomy), "kingdom___phylum___class___order___family___genus___species")


class TestSequence(unittest.TestCase):

    def test_sequence_initialization(self):
        header = ">seq1 tax=kingdom;phylum;class;order;family;genus;species;"
        sequence = "ATCGATCGATCG"
        seq = Sequence(header, sequence)
        self.assertEqual(str(seq.header.seq_id), ">seq1")
        self.assertEqual(seq.sequence, "ATCGATCGATCG")

    def test_is_sequence_correct(self):
        header = ">seq1 tax=kingdom;phylum;class;order;family;genus;species;"
        sequence = "ATCGATCGATCG"
        seq = Sequence(header, sequence)
        self.assertTrue(seq.is_sequence_correct())


class TestZOTUFASTA(unittest.TestCase):

    def test_get_taxonomies(self):
        fasta_file = "path/to/fasta_file.fasta"
        zotu_fasta = ZOTUFASTA(fasta_file)
        taxonomies = zotu_fasta.get_taxonomies()
        self.assertIsInstance(taxonomies, list)

    def test_get_hash_table(self):
        fasta_file = "path/to/fasta_file.fasta"
        zotu_fasta = ZOTUFASTA(fasta_file)
        hash_table = zotu_fasta.get_hash_table()
        self.assertIsInstance(hash_table, dict)



if __name__ == '__main__':
    unittest.main()