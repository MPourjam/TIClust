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
    level_tax_map = {
        'kingdom': 0,
        'phylum': 1,
        'class': 2,
        'order': 3,
        'family': 4,
        'genus': 5,
        'species': 6,
        'ZOTU': 7
    }

    def __init__(self, tax_str: str, delimiter: str = ';'):
        self.tax_str = tax_str
        self.tax_list = self.parse_taxonomy(delimiter)
        self.delimiter = delimiter
        pass

    def parse_taxonomy(self, delimiter: str) -> List[str]:
        return self.tax_str.split(delimiter)

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

    def _get_level(self, level: str):
        tax_level_ind = self.level_tax_map[level]
        return self.tax_list[tax_level_ind] if len(self.tax_list) > tax_level_ind else f'NA_{str(level).capitalize()}'

    def get_tax_upto(self, level: str):
        tax_level_ind = self.level_tax_map[level] + 1
        sliced_tax = []
        for lev in range(tax_level_ind):
            sliced_tax.append(self._get_level(lev))
        return "___".join(sliced_tax)

    def __repr__(self):
        return self.tax_str.replace(self.delimiter, self.output_delimiter)

    def __str__(self):
        return self.__repr__()


class Sequence:

    seq_header_regex = re.compile(r"^>$")
    seq_content_regex = re.compile(r"^[ACUGTN]+$")

    def __init__(self, header: str, sequence: str):
        self.header = header
        self.sequence = sequence
        if not self.is_header_correct(header):
            raise ValueError(f"Incorrect header format for {self.header}")
        if not self.is_sequence_correct(sequence):
            raise ValueError(f"Incorrect sequence format for {self.header}")

    def is_header_correct(self, header_regex: re.Pattern = None) -> bool:
        header_regex_to_check = header_regex if header_regex else self.seq_header_regex
        return bool(header_regex_to_check.match(self.header))

    def is_sequence_correct(self) -> bool:
        return bool(self.seq_content_regex.match(self.sequence))

    def __hash__(self):
        return hash((self.header, self.sequence.lower()))


class ZOTU(Sequence):
    zotus_header_regex = re.compile(r"^>(?P<zotu_id>Zotu[0-9]+).*tax=(?P<tax>([^;]+;){7,8});?", re.IGNORECASE)

    def __init__(self, header: str, sequence: str):
        super().__init__(header, sequence)
        if not self.is_header_correct(header):
            raise ValueError(f"Incorrect header format for {self.header}")
        self.taxonomy = Taxonomy(self._get_taxonomy())
        self.zotu_id = self._get_zotu_id()

    def _get_taxonomy(self) -> str:
        parsed_tax_match = self.zotus_header_regex.match(self.header).group('tax')
        return str(parsed_tax_match)

    def is_header_correct(self, header) -> bool:
        return bool(self.zotus_header_regex.match(header))

    def __hash__(self) -> int:
        return super().__hash__()


class ZOTUFASTA:

    def __init__(self, fasta_file: pl.Path):
        self.fasta_file = fasta_file
        pass


class Species:

    def __init__():
        pass


class Genus:

    def __init__():
        pass


class Family:

    def __init__():
        pass


class Order:

    def __init__():
        pass


class Class:

    def __init__():
        pass


class Phylum:

    def __init__():
        pass


class Kingdom:

    def __init__():
        pass


class UClust:
    """
    Objects of this class will take zOTUs and cluster them using UClust.
    """
    uclust_bin = 'uclust'

    def __init__(
            self,
            zotus: List[ZOTU],
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
        # TODO

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