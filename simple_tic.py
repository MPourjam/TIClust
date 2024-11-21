# pylint: disable=missing-docstring
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring
# pylint: disable=too-few-public-methods
# pylint: disable=no-method-argument
import re
from collections import defaultdict
from typing import List, Dict, Tuple, Set, Optional
from concurrent.futures import ThreadPoolExecutor
from collections import OrderedDict
import threading
import tempfile
import pathlib as pl
import subprocess
import shutil


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
    unknown_tax_reg = re.compile(
        r"NA-(Kingdom|Phylum|Class|Order|Family|Genus|Species|ZOTU)+",
        re.IGNORECASE
    )

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
        """
        It returns the taxonomy as a list of strings based on self.tax_str.
        NOTE self.tax_str is the cleaned original taxonomy string without implicit NA-levels.
        """
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
        return self.get_level('kingdom')

    @property
    def phylum(self):
        return self.get_level('phylum')

    @property
    def class_(self):
        return self.get_level('class')

    @property
    def order(self):
        return self.get_level('order')

    @property
    def family(self):
        return self.get_level('family')

    @property
    def genus(self):
        return self.get_level('genus')

    @property
    def species(self):
        return self.get_level('species')

    def get_level(self, level: str) -> str:
        tax_level_ind = self.level_tax_map.index(level)
        if len(self.tax_list) > tax_level_ind:
            return self.tax_list[tax_level_ind]
        return f'NA-{str(level).capitalize()}'

    def set_level(self, level: str, value: str):
        try:
            tax_level_ind = self.level_tax_map.index(level)
        except IndexError as ind_exc:
            raise IndexError(
                f'{ind_exc}: Taxonomy level {level} is not available for {self.tax_str}'
            ) from ind_exc
        except Exception as e:
            raise e
        cur_tax_list = self.get_tax_upto('ZOTU').tax_list
        if len(self.tax_list) > tax_level_ind:
            cur_tax_list[tax_level_ind] = value
        else:
            raise ValueError(f"Taxonomy level {level} is not available for {self.tax_str}")
        self.tax_str = self.delimiter.join(cur_tax_list)

    def is_known_upto(self, level: str) -> bool:
        """
        It gets the level and checks if the taxonomy is known up to the given level.
        """
        tax_level_ind = self.level_tax_map.index(level)
        return len(self.tax_list) > tax_level_ind

    def get_tax_upto(
        self,
        level: str,
        only_knowns: bool = False,
        top_down: bool = True) -> 'Taxonomy':
        """
        It returns taxonomy up to the given level.
            e.g: Kingdom;Phylum;Class;Order;Family;Genus;Species;ZOTU

        Args:
            level (str): The taxonomy level up to which the taxonomy is needed.
            only_knowns (bool): If True, only known taxonomies will be returned.
            top_down (bool): If True, the taxonomy will be returned from Kingdom to ZOTU.
        """
        tax_lv_ind = self.level_tax_map.index(level)
        sliced_tax = []
        tax_range  = range(tax_lv_ind + 1)
        if not top_down:
            tax_range = range(-1, tax_lv_ind - self.complete_taxonomy_length - 1, -1)
        for lev in tax_range:
            cur_level = self.level_tax_map[lev]
            cur_level = self.get_level(cur_level)
            if not only_knowns:
                sliced_tax.append(cur_level)
            elif only_knowns:
                unknown_tax_match = self.unknown_tax_reg.match(cur_level)
                if not unknown_tax_match:
                    sliced_tax.append(cur_level)
            elif not cur_level:
                raise ValueError(f"Taxonomy level {cur_level} is not available for {self.tax_str}")
        tax_ = "tax="
        if top_down:
            tax_ += self.delimiter.join(sliced_tax)
        else:
            tax_ += self.delimiter.join(sliced_tax[::-1])
        return Taxonomy(tax_)

    def _is_tax_complete(self) -> bool:
        # When 'ToBeDone' is in the taxonomy, it is not complete and set for completion by TIC
        to_be_done = any(True for tax in self.tax_list if 'NA-' in tax)
        complete_tax_levels = len(self.tax_list) == self.complete_taxonomy_length
        return complete_tax_levels and not to_be_done

    @property
    def last_known_level(self) -> str:
        known_tax = self.get_tax_upto('ZOTU', only_knowns=True)
        print(known_tax)
        known_tax_levels = known_tax.tax_list
        return self.level_tax_map[len(known_tax_levels) - 1]

    def __repr__(self):
        return f'{self.tax_str}'

    def __str__(self):
        return f"{self.tax_str}"

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


class SequenceGroup:
    """
    Collection of sequences with the same taxonomy.
    """
    __centroid: Optional[Sequence] = None
    __sequences: List[Sequence] = []
    homogenous_tax_group: bool = True

    def __init__(
        self,
        sequences: List[Sequence],
        centroid: Sequence = None,
        homogenous_tax_group: bool = True):
        self.homogenous_tax_group = homogenous_tax_group
        will_be_homogenous = self.__class__.is_homogenous_group(sequences)
        if self.homogenous_tax_group and not will_be_homogenous:
            raise ValueError("All sequences should have the same taxonomy.")
        self.__sequences = sequences
        self.__set_centroid(centroid)

    @property
    def sequences(self) -> List[Sequence]:
        return self.__sequences

    @sequences.setter
    def sequences(self, sequences: List[Sequence]):
        will_be_homogenous = self.__class__.is_homogenous_group(sequences)
        if self.homogenous_tax_group and not will_be_homogenous:
            raise ValueError("All sequences should have the same taxonomy.")
        self.__sequences = sequences
        self.__centroid = None

    @property
    def centroid(self) -> Sequence:
        return self.__centroid

    @centroid.setter
    def centroid(self, centroid: Sequence):
        self.__set_centroid(centroid)

    def __set_centroid(self, centroid: Sequence):
        if centroid not in self.sequences:
            raise ValueError("The centroid should be one of the sequences in the group.")
        self.__centroid = centroid

    @property
    def tax(self) -> Taxonomy:
        if self.__class__.is_homogenous_group(self.sequences):
            return self.sequences[0].header.taxonomy.get_tax_upto('ZOTU', only_knowns=True)
        print(f"The sequences do not have the same taxonomy in {self}")
        return None

    def add_sequence(self, sequence: Sequence):
        """
        It checks if the given sequence has the same taxonomy as the other members of the group.
        """
        will_be_homogenous = self.__class__.is_homogenous_group([sequence] + self.sequences)
        if self.homogenous_tax_group and not will_be_homogenous:
            raise ValueError(
                "The sequence does not have the same taxonomy as the other members of the group."
            )
        self.__sequences.append(sequence)
        self.__centroid = None


    @staticmethod
    def is_homogenous_group(sequences: List[Sequence], upto_level: str = "ZOTU") -> bool:
        taxs_list = [str(seq.header.taxonomy.get_tax_upto(upto_level)) for seq in sequences]
        taxs_set = set(taxs_list)
        return len(taxs_set) == 1

    def set_level(self, level: str, value: str):
        for seq in self.sequences:
            seq.header.taxonomy.set_level(level, value)

    def get_tax_upto(self, level: str, only_knowns: bool = False, top_down: bool = True) -> str:
        return self.sequences[0].header.taxonomy.get_tax_upto(level, only_knowns, top_down)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.tax}, ({len(self.sequences)})"

    def __str__(self):
        return self.__repr__()

    def __hash__(self):
        return hash([seq.__hash__ for seq in self.sequences])

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def __len__(self):
        return len(self.sequences)

    def __iter__(self):
        return iter(self.sequences)

    def __bool__(self):
        return bool(self.sequences)

    def write_to_fasta(self, output_file_path: pl.Path):
        with open(output_file_path, 'w', encoding='utf-8') as fasta:
            for seq in self.sequences:
                fasta.write(str(seq) + '\n')


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
        """
        It parses pairs of lines in the fasta file and returns a list of sequence objects.
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

    def __eq__(self, other):
        return self.node_name == other.node_name


class TaxononomyTree:

    delimiter = ';'

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

    def find_node(self, node_name: str) -> TreeNode:
        pass

    def traverse_tree(self, node: TreeNode = None) -> List[List[TreeNode]]:
        node = node if node else self.root
        result = [[node]]
        if not node.children:
            return result
        for child in node.children:
            for path in self.traverse_tree(child):
                result.append([node] + path)
        return result

    def get_lineages(self) -> List[Taxonomy]:
        lineages_list = self.traverse_tree()
        lineages_tax = []
        for lineage in lineages_list:
            tax_str = self.delimiter.join([node.node_name for node in lineage])
            tax = Taxonomy(tax_str)
            lineages_tax.append(tax)
        return lineages_tax

    def add_lineage(self, taxonomy: Taxonomy):
        if taxonomy not in self.tree_lineages:
            self.tree_lineages.append(taxonomy)

        for level in taxonomy.tax_list:
            if not level:
                continue
            if not any(True for node in self._get_nodes() if node.node_name == level):
                parent = self.root
                for node_name in taxonomy.tax_list:
                    child = TreeNode(node_name, parent)
                    parent.add_child(child)
                    parent = child

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

    def tax_set_at_known_level(self, level: str) -> List[Taxonomy]:
        """
        It returns a list of sequences that their knwon taxonomy is known up to the given level.

        Args:
            level (str): The taxonomy level up to which the taxonomy is needed.
            only_knowns (bool): If True, only known taxonomies will be returned.
        """
        taxs = self.tax_obj_list
        tax_list = []
        for tax in taxs:
            if tax.last_known_level == level:
                tax_list.append(tax)
        tax_list = list(set(tax_list))
        return tax_list


class UClust:
    """
    Objects of this class will take zOTUs and cluster them using UClust.
    """
    uclust_bin = pl.Path('./uclust').absolute()

    def __init__(
            self,
            sequences: List[Sequence],  # this could be a an argument of a method
            uclust_bin: pl.Path = pl.Path('./uclust').absolute(),
            uclust_work_dir: pl.Path = pl.Path('./uclust_work_dir').absolute()
            ):
        self.sequences = sequences
        self.uclust_bin = uclust_bin
        if not self.uclust_bin.exists() or not self.uclust_bin.is_file():
            raise FileNotFoundError(f"UClust binary not found at {self.uclust_bin}")
        if uclust_work_dir.is_file():
            raise FileNotFoundError(f"UClust work directory is a file at {uclust_work_dir}")
        # create a temporary directory for UClust in the uclust_work_dir
        with tempfile.TemporaryDirectory(dir=uclust_work_dir) as temp_cluster_dir:
            # write sequences to a fasta file in the temporary directory
            self.uclust_work_dir = pl.Path(temp_cluster_dir.name)
            self.sequences.write_fasta(self.uclust_work_dir / 'sequences.fasta')


    def cluster(
        self,
        sim_threshold: float,
        ) -> List[SequenceGroup]:
        """
        It clusters the sequences using UClust and returns a dictionary of clusters.
        """
        # get sure that the self.uclust_work_dir is empty
        shutil.rmtree(self.uclust_work_dir)
        centroids_file = self.uclust_work_dir / 'centroids.fasta'
        uc_file = self.uclust_work_dir / 'clusters.uc'
        # TODO
        subprocess.run(
            [
                str(self.uclust_bin),
                self.uclust_work_dir / 'sequences.fasta',
                '-id',
                sim_threshold,
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


class TICAnalysis:

    uclust_bin = pl.Path('./uclust.bin')
    threads = 1
    default_thresholds = OrderedDict({
        'kingdom': None,
        'phylum': None,
        'class': None,
        'order': None,
        'family': 0.90,
        'genus': 0.95,
        'species': 0.97,
    })

    def __init__(self, taxed_fasta_file_path: pl.Path):
        self.fasta_file = TaxedFastaFile(taxed_fasta_file_path)
        self.target_lieages = list(set(list(self.fasta_file.tax_obj_list)))
        self.tax_tree = TaxononomyTree(self.fasta_file.tax_obj_list)
        self.cluster_thresholds = self.default_thresholds

    def run(
        self,
        threads: int = 1,
        cluster_thresholds: Dict[str, float] = None):
        self.cluster_thresholds.update(cluster_thresholds)
        self.complete_tax_at_level('species', threads)

    def complete_tax_at_level(self, level: str, threads: int):
        """
        It takes all sequences of taxonomies known at the given level but unknown at the next level.
        Then tries to cluster them at the next level in a parallel manner.
        """
        tax_lineage_at_level = self.fasta_file.tax_set_at_known_level(level)
        # TODO

    def complete_known_genus(self, threads: int):
        known_genus_taxs = self.fasta_file.tax_set_at_known_level('genus')
        # start a thread pool
        with ThreadPoolExecutor(max_workers=threads) as executor:
            # create a list of tasks
            tasks = []
            for tax in known_genus_taxs:
                sequences = self.fasta_file.filter_seq_by_tax(tax)
                tasks.append(executor.submit(self.complete_known_genus, tax))
                # TODO

    def complete_known_family(self, threads: int):
        pass
        # TODO

    def complete_known_order(self, threads: int):
        pass
        # TODO
