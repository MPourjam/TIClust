# pylint: disable=missing-docstring
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring
# pylint: disable=too-few-public-methods
# pylint: disable=no-method-argument
import re
from collections import defaultdict
from typing import List, Dict, Tuple, Set
from concurrent.futures import ThreadPoolExecutor
import threading
import tempfile
import pathlib as pl
import subprocess


class Taxonomy:

    output_delimiter = '___'
    delimiter = ';'  # Default delimiter for taxonomy
    tax_regex = re.compile(r"tax=(?P<tax>([^;]+;)*([^;]+)?;?)$", re.IGNORECASE)
    complete_taxonomy_length = 8
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

    def __init__(self, tax_str: str, delimiter: str = None):
        self.delimiter = delimiter if delimiter else self.delimiter
        tax_reg_match = self.tax_regex.search(tax_str)
        self.tax_str = ''
        if tax_reg_match:
            tax_str = tax_reg_match.group('tax')
            self.tax_str = self.get_clean_taxonomy(tax_str, self.delimiter)
        truncated_tax = self.delimiter.join(self._get_tax_list())
        self.tax_str = truncated_tax

    @property
    def tax_list(self) -> List[str]:
        return self._get_tax_list(self.delimiter)

    def _get_tax_list(self, delimiter: str = None) -> List[str]:
        delimiter = delimiter if delimiter else self.delimiter
        # Any taxonomy with more than complete_taxonomy_length levels is truncated
        return self.tax_str.split(delimiter)[:self.complete_taxonomy_length]

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

    def _set_level(self, level: str, value: str):
        tax_level_ind = self.level_tax_map.index(level)
        cur_tax_list = self.get_tax_upto('ZOTU').split(self.delimiter)
        if len(self.tax_list) > tax_level_ind:
            cur_tax_list[tax_level_ind] = value
        else:
            raise ValueError(f"Taxonomy level {level} is not available for {self.tax_str}")
        self.tax_str = self.delimiter.join(cur_tax_list)

    def get_tax_upto(self, level: str, only_knowns: bool = False) -> str:
        tax_level_ind = self.level_tax_map.index(level) + 1
        sliced_tax = []
        unknown_tax_reg = re.compile(
            r"NA-(Kingdom|Phylum|Class|Order|Family|Genus|Species)+",
            re.IGNORECASE
        )
        for lev in range(tax_level_ind):
            cur_level = self.level_tax_map[lev]
            cur_level = self._get_level(cur_level)
            if not only_knowns:
                sliced_tax.append(cur_level)
            elif only_knowns:
                unknown_tax_match = unknown_tax_reg.match(cur_level)
                if not unknown_tax_match:
                    sliced_tax.append(cur_level)
            elif not cur_level:
                raise ValueError(f"Taxonomy level {cur_level} is not available for {self.tax_str}")
        return self.delimiter.join(sliced_tax)

    def _is_tax_complete(self) -> bool:
        # When 'ToBeDone' is in the taxonomy, it is not complete and set for completion by TIC
        to_be_done = any(True for tax in self.tax_list if 'ToBeDone' in tax)
        complete_tax_levels = len(self.tax_list) == self.complete_taxonomy_length
        return complete_tax_levels and not to_be_done

    @property
    def last_known_tax_level(self) -> str:
        return self._get_last_known_tax_level()

    def _get_last_known_tax_level(self) -> str:
        length_known_tax = len(self.get_tax_upto('ZOTU', True).split(self.delimiter))
        return self.level_tax_map[length_known_tax - 1] if length_known_tax > 0 else ''

    def __repr__(self):
        return self.get_tax_upto('ZOTU').replace(self.delimiter, self.output_delimiter)

    def __str__(self):
        return self.__repr__().replace(self.output_delimiter, self.delimiter)

    def __hash__(self):
        return hash(self.tax_str)

    def __bool__(self):
        return bool(self.tax_str)

    def __eq__(self, other):
        # NOTE taxonomy should not be case-sensitive
        tax_list_lower = [tax.lower() for tax in self.tax_list]
        other_tax_list_lower = [tax.lower() for tax in other.tax_list]
        return tax_list_lower == other_tax_list_lower


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
        return f"{self.header}\n{self.sequence}"


class FastaFile:

    def __init__(self, fasta_file_path: pl.Path):
        self.fasta_file_path = fasta_file_path.absolute()

    def get_hash_table(self) -> Dict[int, str]:
        # NOTE suprisingly, this does not help reduce memory usage as compared to get_sequences()
        hash_table = {}
        with open(self.fasta_file_path, 'r', encoding='utf-8') as fasta:
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

    def get_seq_headers(self) -> List[SeqHeader]:
        headers = []
        with open(self.fasta_file_path, 'r', encoding='utf-8') as fasta:
            curr_line = fasta.readline().strip()
            while curr_line:
                if curr_line.startswith('>'):
                    headers.append(SeqHeader(curr_line))
                curr_line = fasta.readline().strip()
        return headers

    def get_sequences(self) -> List[Sequence]:
        sequences = []
        with open(self.fasta_file_path, 'r', encoding='utf-8') as fasta:
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

    def filter_seq_by_hash(self, hash_list: List[int]) -> List[Sequence]:
        """
        It returns a list of sequence objects whose hash values are in the hash_list.
        """
        sequences = []
        with open(self.fasta_file_path, 'r', encoding='utf-8') as fasta:
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

    def write_fasta_file(self, output_file_path: pl.Path = None, sequences: List[Sequence] = None):
        if not output_file_path:
            output_file_path = self.fasta_file_path
        with open(output_file_path, 'w', encoding='utf-8') as fasta:
            for seq in sequences:
                fasta.write(str(seq) + '\n')


class TreeNode:

    def __init__(self, node_name: str, parent_name: 'TreeNode' = None):
        self.node_name = node_name
        self.parent_name = parent_name
        self.children = set()

    def add_child(self, child: 'TreeNode'):
        self.children.add(child)

    def get_parents(self) -> Tuple['TreeNode']:
        parents = []
        parent = self.parent_name
        while parent:
            parents.append(parent)
            parent = parent.parent_name
        return tuple(parents)

    def __repr__(self):
        return self.node_name

    def __str__(self):
        return self.__repr__()

    def __hash__(self):
        return hash(self.node_name)


class TaxononomyTree:

    def __init__(self, taxonomies: List[Taxonomy]):
        self.tree_lineages = list(set(taxonomies))
        self.root = self.plant_tree()

    def plant_tree(self) -> TreeNode:
        root = TreeNode('root')
        for tax in self.tree_lineages:
            tax_levels = tax.tax_list
            parent = root
            for level in tax_levels:
                child = TreeNode(level, parent)
                parent.add_child(child)
                parent = child
        return root

    def _get_nodes(self) -> Set[TreeNode]:
        nodes = set()
        def get_nodes(node: TreeNode):
            nodes.add(node)
            for child in node.children:
                get_nodes(child)
        get_nodes(self.root)
        return nodes

    def add_lineage(self, taxonomy: Taxonomy):
        if taxonomy not in self.tree_lineages:
            self.tree_lineages.append(taxonomy)
            self.root = self.plant_tree()

    def node_exist(self, node_obj: TreeNode) -> bool:
        nodes = self._get_nodes()
        return any(True for node in nodes if node.node_name == node_obj.node_name)


class TaxedFastaFile(FastaFile):

    def __init__(self, fasta_file_path: pl.Path):
        super().__init__(fasta_file_path)
        self.taxonomies_objs_list = self.__get_seq_taxonomies()

    @property
    def tax_obj_list(self) -> List[Taxonomy]:
        return [header.taxonomy for header in self.get_seq_headers()]

    def __get_seq_taxonomies(self) -> List[Taxonomy]:
        headers = self.get_seq_headers()
        taxonomies = [header.taxonomy for header in headers]
        if any(True for tax in taxonomies if not tax):
            print("Sodme sequences do not have any taxonomy information.")
        return taxonomies

    def filter_seq_by_tax(self, taxonomy: Taxonomy) -> List[Sequence]:
        sequences = []
        with open(self.fasta_file_path, 'r', encoding='utf-8') as fasta:
            curr_line = fasta.readline().strip()
            seq_str = ''
            seq_header = ''
            while curr_line:
                if curr_line.startswith('>'):
                    if seq_str:
                        seq_obj = Sequence(seq_header, seq_str)
                        if seq_obj.header.taxonomy == taxonomy:
                            sequences.append(seq_obj)
                        seq_str = ''
                    seq_header = curr_line
                else:
                    seq_str += curr_line
                curr_line = fasta.readline().strip()
            if seq_str:
                seq_obj = Sequence(seq_header, seq_str)
                if seq_obj.header.taxonomy == taxonomy:
                    sequences.append(seq_obj)
        return sequences

    def get_tax_seq_map(self) -> Dict[Taxonomy, Sequence]:
        tax_seq_map = defaultdict(list)
        with open(self.fasta_file_path, 'r', encoding='utf-8') as fasta:
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

    def group_by_taxonomic_level(self, level: str) -> Dict[str, List[Sequence]]:
        groups = defaultdict(list)
        for seq in self.get_sequences():
            key = seq.header.taxonomy.get_tax_upto(level)
            if key:
                groups[key].append(seq)
        return groups


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
            zotus: List[Sequence],  # this could be a an argument of a method
            sim_threshold: float,
            uclust_bin: str = 'uclust',
            threads: int = 1,
            uclust_work_dir: pl.Path = pl.Path(tempfile.mkdtemp())):
        self.zotus = zotus
        self.sim_threshold = sim_threshold
        self.uclust_bin = uclust_bin
        self.threads = threads
        self.uclust_work_dir = pl.Path(pl.PurePath(uclust_work_dir)).absolute()
        # TODO

    def gather_zotus(self) -> pl.Path:
        # TODO
        pass

    def cluster(self) -> Dict[str, List[str]]:
        gathered_zotus = self.gather_zotus()
        centroids_file = self.uclust_work_dir / 'centroids.fasta'
        uc_file = self.uclust_work_dir / 'clusters.uc'
        # TODO
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
        # TODO

    def parse_uc(self, uc_file: pl.Path) -> Dict[str, List[str]]:
        pass

    def write_uc(self, uc_file: pl.Path):
        pass

    def write_fasta(self, fasta_file: pl.Path):
        pass

    def write_clusters(self, clusters: Dict[str, List[str]]):
        pass


class TIC:

    uclust_bin = pl.Path('./uclust.bin')
    def __init__(self, taxed_fasta_file_path: pl.Path):
        self.fasta_file = TaxedFastaFile(taxed_fasta_file_path)
        self.tax_tree = TaxononomyTree(self.fasta_file.tax_obj_list)

    def run(
        self,
        species_threshold: float = 0.97,
        genus_threshold: float = 0.95,
        family_threshold: float = 0.90,
        threads: int = 1):
        pass

    def group_by_taxonomic_level(self, level: str):
        return self.fasta_file.group_by_taxonomic_level(level)

    def group_by_taxonomic_level(self, level: str):
        pass

    def find_centroid(self, sequences: List[Sequence]):
        pass

    def cluster_sequences(self, similarity_threshold: float):
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


class ZOTUTable(FeatureTable):
    """
    Objects of this class will take zOTUs and cluster them to OTUs.
    """

    def __init__():
        pass
