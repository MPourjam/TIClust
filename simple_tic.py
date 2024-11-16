import re
from collections import defaultdict
from typing import List, Dict, Tuple
from concurrent.futures import ThreadPoolExecutor
import threading
import tempfile
import pathlib as pl
import subprocess


class Taxonomy:

    output_delimiter = '___'
    delimiter = ';'  # Default delimiter for taxonomy
    tax_regex = re.compile(r"tax=(?P<tax>([^;]+;)*[^;]+;?)$", re.IGNORECASE)
    level_tax_map = [
        'kingdom',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species',
        'ZOTU'
    ]
    invalid_tax_rgx = [
        re.compile(r"uncultured", re.IGNORECASE),
        re.compile(r"unclassified", re.IGNORECASE),
        re.compile(r"metagenome", re.IGNORECASE),
        re.compile(r"unidentified", re.IGNORECASE),
        re.compile(r"unknown", re.IGNORECASE),
        re.compile(r"unidentified", re.IGNORECASE),
        re.compile(r"unking", re.IGNORECASE),
        re.compile(r"unphyl", re.IGNORECASE),
        re.compile(r"unorder", re.IGNORECASE),
        re.compile(r"unclass", re.IGNORECASE),
        re.compile(r"unfamily", re.IGNORECASE),
        re.compile(r"ungenus", re.IGNORECASE),
        re.compile(r"unspec", re.IGNORECASE),
        # # we don't filter incertae sedis as it could have sublevels
        # re.compile(r"incertae", re.IGNORECASE),
    ]

    def __init__(self, tax_str: str):
        tax_reg_match = self.tax_regex.search(tax_str)
        if tax_reg_match:
            self.tax_str = tax_reg_match.group('tax')
        else:
            self.tax_str = ''
        self.tax_list = self.get_tax_list(self.delimiter)

    def get_tax_list(self, delimiter: str = None) -> List[str]:
        delimiter = delimiter if delimiter else self.delimiter
        clean_tax = self.get_clean_taxonomy(self.tax_str, delimiter)
        return clean_tax.split(delimiter)

    @classmethod
    def get_clean_taxonomy(cls, tax_str: str, delimiter: str = None) -> str:
        # NOTE needs testing
        delimiter = delimiter if delimiter else cls.delimiter
        clean_tax = []
        for curr_tax in tax_str.split(delimiter):
            for invalid_tax in cls.invalid_tax_rgx:
                invalid_tax_match = invalid_tax.search(curr_tax)
                if not invalid_tax_match and curr_tax:
                    clean_tax.append(curr_tax)
                    break

        return delimiter.join(clean_tax)

    @property
    def kingdom(self):
        return self._get_level('kingdom')

    @property
    def phylum(self):
        return self._get_level('phylum')

    @property
    def class_(self):
        return self._get_level('class')

    @property
    def order(self):
        return self._get_level('order')

    @property
    def family(self):
        return self._get_level('family')

    @property
    def genus(self):
        return self._get_level('genus')

    @property
    def species(self):
        return self._get_level('species')

    def _get_level(self, level: str) -> str:
        tax_level_ind = self.level_tax_map.index(level)
        if len(self.tax_list) > tax_level_ind:
            return self.tax_list[tax_level_ind]
        return f'NA-{str(level).capitalize()}'

    def get_tax_upto(self, level: str) -> str:
        tax_level_ind = self.level_tax_map.index(level) + 1
        sliced_tax = []
        for lev in range(tax_level_ind):
            cur_level = self.level_tax_map[lev]
            if self._get_level(cur_level):
                sliced_tax.append(self._get_level(cur_level))
            else:
                raise ValueError(f"Taxonomy level {cur_level} is not available for {self.tax_str}")
        return "___".join(sliced_tax)

    def __repr__(self):
        return self.get_tax_upto('ZOTU').replace(self.delimiter, self.output_delimiter)

    def __str__(self):
        return self.__repr__()
    
    def __hash__(self):
        return hash(self.tax_str)


class SeqID:

    seq_id_regex = re.compile(r"^(?P<seq_header>>\S+)\s*.*$", re.IGNORECASE)

    def __init__(self, header: str):
        seq_id_match = self.seq_id_regex.match(header)
        if not seq_id_match:
            raise ValueError(f"Incorrect header format for {seq_id_match}")
        self.head_id = str(seq_id_match.group('seq_header'))

    def __repr__(self):
        return self.head_id

    def __str__(self):
        return self.__repr__()


class SeqHeader:

    def __init__(self, header: str):
        self.seq_id = SeqID(header)
        self.taxonomy = Taxonomy(header)
        self.header = str(self.seq_id) + ' ' + str(self.taxonomy)

    def __repr__(self):
        return self.header

    def __str__(self):
        return self.__repr__()

    def __hash__(self):
        return hash(self.header)


class Sequence:

    seq_header_regex = re.compile(r"^>")
    # NOTE in some SILVA entries there are 'K' 'M', and 'N' characters in the sequence
    seq_content_regex = re.compile(r"^[A-Z]+$", re.IGNORECASE)

    def __init__(self, header: str, sequence: str):
        self.header = SeqHeader(header)
        self.sequence = sequence
        if not self.is_sequence_correct():
            raise ValueError(f"Incorrect sequence format for \n{self.header}\n{self.sequence}\n")

    def is_header_correct(self, header_regex: re.Pattern = None) -> bool:
        header_regex_to_check = header_regex if header_regex else self.seq_header_regex
        return bool(header_regex_to_check.match(str(self.header)))

    def is_sequence_correct(self) -> bool:
        return bool(self.seq_content_regex.match(self.sequence))
    
    def get_tax_seq_map(self) -> Dict[Taxonomy, 'Sequence']:
        return {self.header.taxonomy: self}

    def __hash__(self):
        return hash((str(self.header), self.sequence.lower()))

    def __repr__(self):
        return f"{self.header.seq_id} {self.sequence[:7]}...{self.sequence[-7:]}"

    def __str__(self):
        return self.__repr__()


class ZOTUFASTA:

    def __init__(self, fasta_file: pl.Path):
        self.fasta_file = fasta_file.absolute()
        pass

    def get_taxonomies(self) -> List[Tuple[SeqID, Taxonomy]]:
        taxonomies = []
        with open(self.fasta_file, 'r', encoding='utf-8') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    taxonomies.append(Taxonomy(line))
        return taxonomies

    def get_tax_seq_map(self) -> Dict[Taxonomy, Sequence]:
        tax_seq_map = defaultdict(list)
        with open(self.fasta_file, 'r', encoding='utf-8') as fasta:
            curr_line = fasta.readline().strip()
            seq_str = ''
            seq_header = ''
            while curr_line:
                if curr_line.startswith('>'):
                    if seq_str:
                        seq_obj = Sequence(seq_header, seq_str)
                        tax_seq_map.get(seq_obj.header.taxonomy).append(seq_obj)
                        seq_str = ''
                    seq_header = curr_line
                else:
                    seq_str += curr_line
                curr_line = fasta.readline().strip()
            if seq_str:
                seq_obj = Sequence(seq_header, seq_str)
                tax_seq_map[seq_obj.header.taxonomy] = seq_obj
        return tax_seq_map

    def get_hash_table(self) -> Dict[int, str]:
        # NOTE suprisingly, this does not help reduce memory usage as compared to get_sequences()
        hash_table = {}
        with open(self.fasta_file, 'r', encoding='utf-8') as fasta:
            curr_line = fasta.readline().strip()
            seq_str = ''
            seq_header = ''
            header_no = 0
            while curr_line:
                if curr_line.startswith('>'):
                    if seq_str:
                        seq_obj = Sequence(seq_header, seq_str)
                        seq_obj_hash = hash(seq_obj)
                        hash_table[seq_obj_hash] = seq_obj.header
                        seq_str = ''
                    seq_header = curr_line
                    header_no += 1
                else:
                    seq_str += curr_line
                curr_line = fasta.readline().strip()
            if seq_str:
                seq_obj = Sequence(seq_header, seq_str)
                seq_obj_hash = hash(seq_obj)
                hash_table[seq_obj_hash] = seq_obj.header
        return hash_table

    def get_sequences(self) -> List[Sequence]:
        sequences = []
        with open(self.fasta_file, 'r', encoding='utf-8') as fasta:
            curr_line = fasta.readline().strip()
            seq_str = ''
            seq_header = ''
            while curr_line:
                if curr_line.startswith('>'):
                    if seq_str:
                        seq_obj = Sequence(seq_header, seq_str)
                        sequences.append(seq_obj)
                        seq_str = ''
                    seq_header = curr_line
                else:
                    seq_str += curr_line
                curr_line = fasta.readline().strip()
            if seq_str:
                seq_obj = Sequence(seq_header, seq_str)
                sequences.append(seq_obj)
        return sequences

    def filter_sequences(self, hash_list: List[int]) -> List[Sequence]:
        """
        It returns a list of sequence objects whose hash values are in the hash_list.
        """
        sequences = []
        with open(self.fasta_file, 'r', encoding='utf-8') as fasta:
            curr_line = fasta.readline().strip()
            seq_str = ''
            seq_header = ''
            while curr_line:
                if curr_line.startswith('>'):
                    if seq_str:
                        seq_obj = Sequence(seq_header, seq_str)
                        seq_obj_hash = hash(seq_obj)
                        if seq_obj_hash in hash_list:
                            sequences.append(seq_obj)
                        seq_str = ''
                    seq_header = curr_line
                else:
                    seq_str += curr_line
                curr_line = fasta.readline().strip()
            if seq_str:
                seq_obj = Sequence(seq_header, seq_str)
                seq_obj_hash = hash(seq_obj)
                if seq_obj_hash in hash_list:
                    sequences.append(seq_obj)
        return sequences


class Species:

    def __init__(self):
        pass


class Genus:

    def __init__(self):
        pass


class Family:

    def __init__(self):
        pass


class Order:

    def __init__(self):
        pass


class Class:

    def __init__(self):
        pass


class Phylum:

    def __init__(self):
        pass


class Kingdom:

    def __init__(self):
        pass


class UClust:
    """
    Objects of this class will take zOTUs and cluster them using UClust.
    """
    uclust_bin = 'uclust'

    def __init__(
            self,
            zotus: List[Sequence],
            sim_threshold: float,
            uclust_bin: str = 'uclust',
            threads: int = 1,
            uclust_work_dir: pl.Path = pl.Path(tempfile.mkdtemp())):
        self.zotus = zotus
        self.sim_threshold = sim_threshold
        self.uclust_bin = uclust_bin
        self.threads = threads
        self.uclust_work_dir = pl.Path(pl.PurePath(uclust_work_dir)).absolute()
        pass

    def gather_zotus(self) -> pl.Path:
        pass

    def cluster(self) -> Dict[str, List[str]]:
        gathered_zotus = self.gather_zotus()
        centroids_file = self.uclust_work_dir / 'centroids.fasta'
        uc_file = self.uclust_work_dir / 'clusters.uc'
        subprocess.run(
            [
                *str(self.uclust_bin).strip().split(),
                str(gathered_zotus),
                '-id',
                self.sim_threshold,
                '-strand',
                'both',
                '-top_hits_only',
                '-centroids',
                str(centroids_file),
                '-uc',
                str(uc_file)
            ]
        )
        pass

    def parse_uc(self, uc_file: pl.Path) -> Dict[str, List[str]]:
        pass

    def write_uc(self, uc_file: pl.Path):
        pass

    def write_fasta(self, fasta_file: pl.Path):
        pass

    def write_clusters(self, clusters: Dict[str, List[str]]):
        pass


class FamilyCluster:
    """
    Objects of this class will take zOTUs with taxonomy known to Order level and cluster
        them to Family level.

    """
    prefix = 'FOTU'

    def __init__():
        pass


class GenusCluster:
    """
    Objects of this class will take zOTUs with taxonomy known to Family level and cluster
        them to Genus level.
    """
    prefix = 'GOTU'

    def __init__():
        pass


class SpeciesCluster:
    """
    Objects of this class will take zOTUs with taxonomy known to Genus level and cluster
        them to Species level.
    """
    prefix = 'SOTU'

    def __init__():
        pass


class FeatureTable:
    """
    Parent class for OTUTable and zOTUTable.
    """
    def __init__():
        pass


class OTUTable(FeatureTable):
    """
    Objects of this class will take zOTUs and cluster them to OTUs.
    """
    def __init__():
        pass


class zOTUTable(FeatureTable):
    """
    Objects of this class will take zOTUs and cluster them to OTUs.
    """

    def __init__():
        pass
