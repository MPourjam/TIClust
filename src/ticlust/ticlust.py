# pylint: disable=missing-docstring,missing-class-docstring,missing-function-docstring
# pylint: disable=too-few-public-methods,no-method-argument,too-many-lines
# Description: This is a simple version of the TIC module.
import re
import tempfile
import shutil
import logging
import pathlib as pl
import hashlib
from os import cpu_count
from collections import defaultdict, OrderedDict
from typing import List, Dict, Set, Tuple, Optional, Union
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
from .ticlust_helper.utils import system_sub, onelinefasta
from .ticlust_helper.config import tic_configs


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class Taxonomy:

    output_delimiter = '___'
    delimiter = ';'  # Default delimiter for taxonomy
    # TODO tax_tag could be defined at initialization level for more flexibility.
    # Add regex template and update the tax_regex to include the user-given tax_tag.
    tax_tag = 'tax='
    tax_regex = re.compile(r"(?P<tax_tag>tax=)(?P<tax>([^;]*;)*[^;]*;?)$", re.IGNORECASE)
    complete_taxonomy_length = 8
    level_tax_map = [
        'kingdom',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species'
        # 'ZOTU'
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
        # NOTE incertae sedis could have sublevels
        re.compile(r"incertae", re.IGNORECASE),
    ]
    unknown_tax_reg = re.compile(
        r"NA-(Kingdom|Phylum|Class|Order|Family|Genus|Species|ZOTU)+",
        re.IGNORECASE
    )

    hypo_full_tax_list = [
        "NA-Kingdom",
        "NA-Phylum",
        "NA-Class",
        "NA-Order",
        "NA-Family",
        "NA-Genus",
        "NA-Species"
    ]

    def __init__(self, tax_str: str, delimiter: str = None):
        """
        self.__tax_str: is string representation of taxnomoy without implicit NA-levels and
        any invalid taxonomies from kingdom to farthest continuous known level.
        self.__full_tax: is the full taxonomy string with all levels and implicit NA-levels.

        NOTE: if a taxonomy has several invalid or missed levels between known levels, then
        these taxonomies are considered as invalid. We assume if at any level a taxonomy is 
        known, its parents should also be known. Thus a taxonomy like 
        tax=kingdom;phylum;unclass;unorder;;genus;species;
        should be considered as invalid and only longest continuous known taxonomy should be 
        considered which is tax=kingdom;phylum;
        """
        self.delimiter = delimiter if delimiter else self.delimiter
        tax_reg_match = self.tax_regex.search(tax_str)
        # NOTE self.__tax_str is the cleaned original taxonomy string without implicit NA-levels.
        self.__tax_str = ''
        self.__full_tax = ''
        self.__orig_tax = ''
        if tax_reg_match:
            tax_str = tax_reg_match.group('tax')
            self.tax_tag = tax_reg_match.group('tax_tag') if tax_reg_match.group('tax_tag') else ''
            self.__orig_tax = tax_str
            self.full_tax = tax_str

    @property
    def full_tax(self) -> str:
        return self.delimiter.join(self.full_tax_list)

    @full_tax.setter
    def full_tax(self, full_tax: str, delimiter: str = None):
        # new full-level tax should not have any invalid taxonomies
        self.delimiter = delimiter if delimiter else self.delimiter
        full_tax_ = self.get_clean_taxonomy(full_tax, self.delimiter, force_full_path=True)
        tax_str_ = self.get_clean_taxonomy(full_tax, self.delimiter, force_full_path=False)
        self.__full_tax = full_tax_
        self.__tax_str = tax_str_


    @property
    def tax_str(self) -> str:
        return self.delimiter.join(self.tax_list)

    @tax_str.setter
    def tax_str(self, tax_str: str):
        self.__tax_str = self.get_clean_taxonomy(tax_str, self.delimiter, force_full_path=False)
        self.__full_tax = self.get_clean_taxonomy(tax_str, self.delimiter, force_full_path=True)

    @property
    def orig_tax(self) -> str:
        return self.__orig_tax

    @property
    def tax_list(self) -> List[str]:
        """
        It returns the taxonomy as a list of strings based on self.__tax_str.
        NOTE self.__tax_str is the cleaned original taxonomy string without implicit NA-levels.
        NOTE if there is not valid taxonomy, self.tax_list[0] is ''
        """
        return self._get_tax_list(self.delimiter)

    @property
    def full_tax_list(self) -> List[str]:
        """
        It returns the full taxonomy as a list of strings based on self.__full_tax.
        NOTE self.__full_tax is the full taxonomy string with all levels and implicit NA-levels.
        """
        return self._get_tax_list(self.delimiter, full_level=True)

    def _get_tax_list(self, delimiter: str = None, full_level: bool = False) -> List[str]:
        delimiter = delimiter if delimiter else self.delimiter
        # Any taxonomy with more than complete_taxonomy_length levels is truncated
        if full_level:
            return self.__full_tax.split(delimiter)[:self.complete_taxonomy_length]

        return self.__tax_str.split(delimiter)[:self.complete_taxonomy_length]

    @classmethod
    def get_clean_taxonomy(
            cls,
            tax_str: str,
            delimiter: str = None,
            force_full_path: bool = False) -> str:
        delimiter = delimiter if delimiter else cls.delimiter
        clean_full_path = cls.hypo_full_tax_list.copy()
        tax_list_ = tax_str.split(delimiter)
        for ind, _ in enumerate(clean_full_path):
            valid_tax = True
            curr_tax = tax_list_[ind] if ind < len(tax_list_) else ''
            for invalid_tax in cls.invalid_tax_rgx:
                invalid_tax_match = invalid_tax.search(curr_tax)
                if invalid_tax_match:
                    valid_tax = False
                    break
            if valid_tax and curr_tax:
                clean_full_path[ind] = curr_tax

        clean_longest_continuous_path = []
        for hypo_tax in clean_full_path:
            if cls.unknown_tax_reg.match(hypo_tax):
                break
            clean_longest_continuous_path.append(hypo_tax)
        if not force_full_path:
            return delimiter.join(clean_longest_continuous_path)

        return delimiter.join(clean_full_path)

    @property
    def kingdom(self):
        return self.full_tax_list[0]

    @property
    def phylum(self):
        return self.full_tax_list[1]

    @property
    def class_(self):
        return self.full_tax_list[2]

    @property
    def order(self):
        return self.full_tax_list[3]

    @property
    def family(self):
        return self.full_tax_list[4]

    @property
    def genus(self):
        return self.full_tax_list[5]

    @property
    def species(self):
        return self.full_tax_list[6]

    def get_level(self, level: str) -> str:
        tax_level_ind = self.level_tax_map.index(level)
        # if taxonomy is empty, self.tax_list[0] is ''
        if len(self.tax_list) > tax_level_ind and self.tax_list[0]:
            return self.tax_list[tax_level_ind]
        return f"NA-{str(level).capitalize()}"

    def set_level(self, level: str, value: str):
        # NOTE This is not consistent with the rest of the code
        # It is setting self.__tax_str with NA-levels which is not
        # supposed throughout the code.
        try:
            tax_level_ind = self.level_tax_map.index(level)
        except IndexError as ind_exc:
            raise IndexError(
                f'{ind_exc}: Taxonomy level {level} is not available for {self.__tax_str}'
            ) from ind_exc
        except Exception as e:
            raise e
        cur_tax_list = self.full_tax_list
        cur_tax_list[tax_level_ind] = value
        self.full_tax = self.delimiter.join(cur_tax_list)

    def is_known_upto(self, level: str) -> bool:
        """
        It gets the level and checks if the taxonomy is known up to the given level.
        """
        tax_level_ind = self.level_tax_map.index(level)
        return len(self.tax_list) > tax_level_ind and self.tax_list[0]

    def get_tax_upto(
            self,
            level: str,
            only_knowns: bool = False,
            top_down: bool = True,
            ret_type: str = 'Taxonomy') -> Union['Taxonomy', str]:
        """
        It returns taxonomy up to the given level.
            e.g: Kingdom;Phylum;Class;Order;Family;Genus;Species;ZOTU

        :param level: The taxonomy level up to which the taxonomy is needed.
        :param only_knowns: If True, only known taxonomies will be returned.
        :param top_down: If True, the taxonomy will be returned from Kingdom to ZOTU.
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
                else:
                    break
            elif not cur_level:
                raise ValueError(
                    f"Taxonomy level {cur_level} is not available for {self.__tax_str}"
                )
        tax_ = ""

        if top_down:
            tax_ += self.delimiter.join(sliced_tax)
        else:
            tax_ += self.delimiter.join(sliced_tax[::-1])
        if ret_type == 'Taxonomy':
            # NOTE tax_ is a string without tax_tag. Thus we add tax_tag to it.
            return Taxonomy("tax=" + tax_)
        if ret_type == 'str':
            return tax_
        raise ValueError(f"Invalid return type {ret_type}")

    def _is_tax_complete(self) -> bool:
        # When 'ToBeDone' is in the taxonomy, it is not complete and set for completion by TIC
        to_be_done = any(True for tax in self.full_tax_list if 'NA-' in tax)
        complete_tax_levels = len(self.tax_list) == self.complete_taxonomy_length
        return complete_tax_levels and not to_be_done

    @property
    def last_known_level(self) -> str:
        known_tax = self.get_tax_upto('species', only_knowns=True)
        # NOTE if the known taxonomy is continuous, the last known level is the last level.
        # Otherwise it is the level before the first unknown level.
        known_tax_levels = known_tax.tax_list
        if len(known_tax_levels) > 0 and known_tax_levels[0]:
            return self.level_tax_map[len(known_tax_levels) - 1]
        return ''

    @property
    def last_known_tax(self) -> str:
        known_tax = self.get_tax_upto('species', only_knowns=True)
        return known_tax.get_level(self.last_known_level)

    def __repr__(self):
        return f'<{self.__class__.__name__}>{self.tax_str}'

    def __str__(self):
        if self.kingdom.lower() == 'bacteria'.lower():
            tax_str_ = f'{self.tax_tag}{self.tax_str}'
        elif self.tax_str == '':
            tax_str_ = ''
        else:
            tax_str_ = f'{self.tax_tag}{self.orig_tax}'
        return tax_str_

    def __hash__(self):
        return hash(self.tax_str)

    def __bool__(self):
        return bool(self.tax_str)

    def __eq__(self, other):
        # NOTE taxonomy should not be case-sensitive
        tax_str_eq = self.tax_str.lower() == other.tax_str.lower()
        return tax_str_eq


class SeqID:

    seq_id_regex = re.compile(r"^>(?P<seq_id>[^\s;]+)", re.IGNORECASE)

    def __init__(self, header: str):
        seq_id_match = self.seq_id_regex.match(header)
        if not seq_id_match:
            raise ValueError(f"Incorrect header format for {seq_id_match}")
        self.id_str = str(seq_id_match.group('seq_id'))

    def __repr__(self):
        return f">{self.id_str}"

    def __str__(self):
        return f"{self.id_str}"

    def __eq__(self, other):
        return str(self.id_str) == str(other)


class SeqHeader:

    header_regex = re.compile(r"^>.*$", re.IGNORECASE)

    def __init__(self, header: str):
        if not self.header_regex.match(header):
            raise ValueError(f"Incorrect header format for {header}")
        self.seq_id = SeqID(header)
        self.taxonomy = Taxonomy(header)
        self.orig_header = header

    @property
    def header(self) -> str:
        seq_header = '>' + str(self.seq_id)
        seq_header += ' ' + str(self.taxonomy)
        return seq_header

    def __repr__(self):
        return self.header

    def __str__(self):
        return self.__repr__()

    def __hash__(self):
        return hash(self.header)

    def __eq__(self, other):
        return self.orig_header == other.orig_header


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
        return hash((str(self.header.orig_header), self.sequence.lower()))

    def __repr__(self):
        return f">{self.header.seq_id} {self.sequence[:7]}...{self.sequence[-7:]}"

    def __str__(self):
        return f"{self.header}\n{self.sequence}"

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()


class SequenceCluster:
    """
    Collection of sequences with the same taxonomy.
    """
    __centroid: Optional[Sequence] = None
    __sequences: List[Sequence] = []
    homogenous_tax_group: bool = True
    homogeneity_level: str = 'species'
    tax_levels: List[str] = [
        'kingdom',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species'
        # 'ZOTU'
    ]

    def __init__(
            self,
            sequences: List[Sequence],
            centroid: Sequence = None,
            homogeneity_level: str = 'species',
            force_homogeneity: bool = True):
        self.force_homogeneity = force_homogeneity
        self.homogeneity_level = homogeneity_level
        will_be_homogenous = self.__class__.is_homogenous_group(sequences, self.homogeneity_level)
        if self.force_homogeneity and not will_be_homogenous:
            raise ValueError("Forcing homogeneity but the sequences are not homogenous.")
        self.homogenous_tax_group = will_be_homogenous
        self.__sequences = sequences
        if not centroid and will_be_homogenous:
            self.__centroid = sequences[0]
        else:
            self.__centroid = centroid

    @property
    def sequences(self) -> List[Sequence]:
        return self.__sequences

    @sequences.setter
    def sequences(self, sequences: List[Sequence], centroid: Sequence = None):
        will_be_homogenous = self.__class__.is_homogenous_group(
            sequences, self.homogeneity_level
        )
        if self.force_homogeneity and not will_be_homogenous:
            raise ValueError("All sequences should have the same taxonomy.")
        self.homogenous_tax_group = will_be_homogenous
        self.__sequences = sequences
        self.centroid = centroid

    @property
    def centroid(self) -> Sequence:
        return self.__centroid

    @centroid.setter
    def centroid(self, centroid: Sequence):
        self.__set_centroid(centroid)

    def __set_centroid(self, centroid: Sequence):
        if centroid not in self.sequences:
            # if centroid is None then we do not end up here as None is element of all lists
            raise ValueError("The centroid should be one of the sequences in the group.")
        self.__centroid = centroid

    @property
    def taxas(self) -> List[Taxonomy]:
        return [seq.header.taxonomy for seq in self.sequences]

    @property
    def closest_common_ancestor(self) -> Taxonomy:
        """
        It returns the common ancestor of the sequences' taxonomies in the group.
        """
        return self.get_closest_common_ancestor()

    def get_closest_common_ancestor(self) -> Taxonomy:
        """
        It returns the closest common ancestor of the sequences' taxonomies in the group.
        # NOTE what if the sequences have different taxonomies even from kingdom level?
        # TIC does not work with non-bacterial sequences.
        """
        common_ancestor = []
        for tax_level in self.tax_levels:
            taxs_list = [
                seq.header.taxonomy.get_tax_upto(tax_level, only_knowns=True, ret_type='str')
                # if tax_level='Bacteria' and is unknown, we get ''
                for seq in self.sequences if seq.header.taxonomy.is_known_upto('kingdom')
            ]
            if len(set(taxs_list)) == 1:
                common_ancestor = taxs_list[0]
            else:
                break
        return Taxonomy('tax=' + common_ancestor)

    @property
    def tax(self) -> Taxonomy:
        is_homogenous_group = self.__class__.is_homogenous_group(
            self.sequences,
            upto_level=self.homogeneity_level
        )
        tax_ = Taxonomy('tax=')
        if len(self.sequences) == 0:
            return tax_
        if is_homogenous_group:
            tax_ = self.sequences[0].header.taxonomy.get_tax_upto('species', only_knowns=True)
        elif not is_homogenous_group:
            tax_ = self.closest_common_ancestor
        return tax_

    def add_sequence(self, sequence: Sequence):
        """
        It checks if the given sequence has the same taxonomy as the other members of the group.
        """
        will_be_homogenous = self.__class__.is_homogenous_group(
            [sequence] + self.sequences,
            upto_level=self.homogeneity_level
        )
        if self.force_homogeneity and not will_be_homogenous:
            raise ValueError(
                "The sequence does not have the same taxonomy as the other members of the group."
            )
        self.homogenous_tax_group = will_be_homogenous
        self.__sequences.append(sequence)
        self.__centroid = None

    @property
    def cluster_size(self) -> int:
        return len(self.sequences)

    def __str__(self):
        return f"<{self.__class__.__name__}>{self.tax}|size={len(self.sequences)}"

    def __repr__(self):
        return self.__str__()

    @staticmethod
    def is_homogenous_group(sequences: List[Sequence], upto_level: str = "species") -> bool:
        """
        It checks if all sequences in the group have the same taxonomy up to the given level.

        :param sequences: The list of sequences to be checked.
        :param upto_level: The taxonomy level up to which the taxonomy should be checked.
            default: species
        """
        taxs_list = [
            str(seq.header.taxonomy.get_tax_upto(upto_level))
            for seq in sequences if seq.header.taxonomy
        ]
        taxs_set = set(taxs_list)
        return len(taxs_set) == 1

    def set_level(self, level: str, value: str):
        # NOTE taxonomy.set_level is not logically consistent with the rest of the code
        if len(self.sequences) == 0:
            print("No sequences in the group.")
        for seq in self.sequences:
            seq.header.taxonomy.set_level(level, value)

    def get_tax_upto(self, level: str, only_knowns: bool = False, top_down: bool = True) -> str:
        if len(self.sequences) == 0:
            return ''
        return self.sequences[0].header.taxonomy.get_tax_upto(level, only_knowns, top_down)

    @property
    def min_len(self) -> int:
        return min(len(seq.sequence) for seq in self.sequences)

    @property
    def max_len(self) -> int:
        return max(len(seq.sequence) for seq in self.sequences)

    def __hash__(self):
        return hash([seq.__hash__ for seq in self.sequences])

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def __len__(self):
        return len(self.sequences)

    def __iter__(self):
        # iter over sequences. First sequence is the centroid.
        yield from self.sequences

    def __bool__(self):
        return bool(self.sequences)

    def write_to_fasta(self, output_file_path: pl.Path, mode: str = 'w'):
        with open(output_file_path, mode, encoding='utf-8') as fasta:
            if self.centroid:
                fasta.write(f"{str(self.centroid)}\n")
            for seq in self.sequences:
                if seq == self.centroid:
                    continue
                fasta.write(str(seq) + '\n')



class SpeciesCluster(SequenceCluster):

    def __init__(self, sequences: List[Sequence], centroid: Sequence = None):
        if not self.is_homogenous_group(sequences, upto_level='species'):
            raise ValueError(
                "All sequences should have the same taxonomy up to species level."
            )
        super().__init__(
            sequences,
            centroid,
            force_homogeneity=True,
            homogeneity_level='species'
        )


class GenusCluster(SequenceCluster):

    def __init__(self, sequences: List[Sequence], centroid: Sequence = None):
        if not self.is_homogenous_group(sequences, upto_level='genus'):
            raise ValueError(
                "All sequences should have the same taxonomy up to genus level."
            )
        super().__init__(
            sequences,
            centroid,
            force_homogeneity=True,
            homogeneity_level='genus'
        )


class FamilyCluster(SequenceCluster):

    def __init__(self, sequences: List[Sequence], centroid: Sequence = None):
        if not self.is_homogenous_group(sequences, upto_level='family'):
            raise ValueError(
                "All sequences should have the same taxonomy up to genus level."
            )
        super().__init__(
            sequences,
            centroid,
            force_homogeneity=True,
            homogeneity_level='family'
        )


class OrderCluster(SequenceCluster):

    def __init__(self, sequences: List[Sequence], centroid: Sequence = None):
        if not self.is_homogenous_group(sequences, upto_level='order'):
            raise ValueError(
                "All sequences should have the same taxonomy up to genus level."
            )
        super().__init__(
            sequences,
            centroid,
            force_homogeneity=True,
            homogeneity_level='order'
        )


class ClassCluster(SequenceCluster):

    def __init__(self, sequences: List[Sequence], centroid: Sequence = None):
        if not self.is_homogenous_group(sequences, upto_level='class'):
            raise ValueError(
                "All sequences should have the same taxonomy up to genus level."
            )
        super().__init__(
            sequences,
            centroid,
            force_homogeneity=True,
            homogeneity_level='class'
        )


class PhylumCluster(SequenceCluster):

    def __init__(self, sequences: List[Sequence], centroid: Sequence = None):
        if not self.is_homogenous_group(sequences, upto_level='phylum'):
            raise ValueError(
                "All sequences should have the same taxonomy up to genus level."
            )
        super().__init__(
            sequences,
            centroid,
            force_homogeneity=True,
            homogeneity_level='phylum'
        )


class KingdomCluster(SequenceCluster):

    def __init__(self, sequences: List[Sequence], centroid: Sequence = None):
        if not self.is_homogenous_group(sequences, upto_level='kingdom'):
            raise ValueError(
                "All sequences should have the same taxonomy up to genus level."
            )
        super().__init__(
            sequences,
            centroid,
            force_homogeneity=True,
            homogeneity_level='kingdom'
        )


class FastaFile:

    def __init__(self, fasta_file_path: pl.Path):
        self.fasta_file_path = fasta_file_path.absolute()

    def get_hash_table(self) -> Dict[int, str]:
        # NOTE suprisingly, this does not help reduce memory footprint
        # as compared to get_sequences()
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

    def write_to_fasta_file(
        self,
        output_file_path: pl.Path = None,
        sequences: List[Sequence] = None,
        mode: str = 'w'):
        if not output_file_path:
            output_file_path = self.fasta_file_path
        with open(output_file_path, mode, encoding='utf-8') as fasta:
            for seq in sequences:
                fasta.write(str(seq) + '\n')

    def get_seq_by_seq_id(self, seq_ids_list: list = []) -> List[Sequence]:
        output_list = []
        if not seq_ids_list:
            return output_list
        with open(self.fasta_file_path, 'r', encoding='utf-8') as fasta_fio:
            curr_line = fasta_fio.readline().strip()
            seq_str = ''
            seq_header = ''
            while curr_line:
                if curr_line.startswith('>'):
                    if seq_str:
                        seq_obj = Sequence(seq_header, seq_str)
                        if str(seq_obj.header.seq_id) in seq_ids_list:
                            output_list.append(seq_obj)
                        seq_str = ''
                    seq_header = curr_line
                else:
                    seq_str += curr_line
                curr_line = fasta_fio.readline().strip()
            if seq_str:
                seq_obj = Sequence(seq_header, seq_str)
                if str(seq_obj.header.seq_id) in seq_ids_list:
                    output_list.append(seq_obj)

        return output_list

    def get_seq_ids(self) -> List[str]:
        seq_ids = []
        with open(self.fasta_file_path, 'r', encoding='utf-8') as fasta:
            curr_line = fasta.readline().strip()
            while curr_line:
                if curr_line.startswith('>'):
                    seq_id = SeqID(curr_line)
                    seq_ids.append(str(seq_id))
                curr_line = fasta.readline().strip()
        return seq_ids


class TaxedFastaFile(FastaFile):

    @property
    def tax_obj_list(self) -> List[Taxonomy]:
        return self.__get_seq_taxonomies(uniq=False)

    @property
    def tax_obj_set(self) -> Set[Taxonomy]:
        return self.__get_seq_taxonomies(uniq=True)

    def __get_seq_taxonomies(self, uniq=True) -> List[Taxonomy]:
        headers = self.get_seq_headers()
        # NOTE we might need to omit the check for tax as sequences without
        # taxnomoies are planned to be written to output file untouched
        taxonomies = []
        for header in headers:
            tax = header.taxonomy
            taxonomies.append(tax)
        if uniq:
            return set(taxonomies)
        return taxonomies

    def filter_seq_by_tax(self, taxonomy: Taxonomy) -> List[Sequence]:
        sequences = []
        logging.debug("Filtering sequences by %s", taxonomy)
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
                        seq_list = tax_seq_map[seq_obj.header.taxonomy]
                        seq_list.append(seq_obj)
                        tax_seq_map[seq_obj.header.taxonomy] = seq_list
                        seq_str = ''
                    seq_header = curr_line
                else:
                    seq_str += curr_line
                curr_line = fasta.readline().strip()
            if seq_str:
                seq_obj = Sequence(seq_header, seq_str)
                seq_list = tax_seq_map[seq_obj.header.taxonomy]
                seq_list.append(seq_obj)
                tax_seq_map[seq_obj.header.taxonomy] = seq_list
        return tax_seq_map

    def filter_tax_set_at_last_known_level(self, level: str) -> List[Taxonomy]:
        """
        It returns a list of sequences that their knwon taxonomy is known up to the given level.

        :param level: The taxonomy level up to which the taxonomy is needed.
        :param only_knowns: If True, only known taxonomies will be returned.
        """
        taxs = self.tax_obj_set
        tax_list = []
        for tax in taxs:
            if tax.last_known_level == level:
                tax_list.append(tax)
        sorted(tax_list, key=lambda x: x.tax_str)
        return tax_list

    def write_seqs_last_known_at(self, level: str, dest_fasta_path: pl.Path, mode: str = 'w'):
        """
        It writes sequences whose last known level is the given level to the destination fasta file.

        :param level: The taxonomy level up to which the taxonomy is needed.
        :param only_knowns: If True, only known taxonomies will be returned.
        :param dest_fasta_path: The path to the destination fasta file.
        """
        taxs = self.filter_tax_set_at_last_known_level(level)
        with open(dest_fasta_path, mode, encoding='utf-8') as fasta:
            for tax in taxs:
                seqs = self.filter_seq_by_tax(tax)
                for seq in seqs:
                    fasta.write(str(seq) + '\n')

    def get_map_to_parent(
            self,
            parent_level: str,
            child_level: str) -> List[Tuple[str, str]]:
        """
        It returns a list of tuples of parent and child taxonomies.
        It works upto child_level = 'species'.

        :param parent_level: The parent taxonomy level.
        :param child_level: The child taxonomy level.

        :return: List[Tuple[str, str]]
        """
        taxs = self.tax_obj_set
        map_list = []
        for tax in taxs:
            if tax.is_known_upto(child_level):
                parent_tax = tax.get_level(parent_level)
                child_tax = tax.get_level(child_level)
                map_list.append((parent_tax, child_tax))
        sorted(map_list, key=lambda x: x[0])

        return list(set(map_list))

    def filter_seq_by_level_value(self, level: str, value: str) -> List[Sequence]:
        """
        It returns a list of sequences that have the given value at the given level.

        :param level: The taxonomy level to be checked.
        :param value: The value to be checked.

        :return: List[Sequence]
        """
        taxs = self.tax_obj_set
        seqs = []
        for tax in taxs:
            if tax.get_level(level) == value:
                logging.debug("Checking %s", tax)
                seqs += self.filter_seq_by_tax(tax)
        return seqs


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


class ZOTUTable:

    def __init__(self, zotu_table_file: pl.Path, sep: str = '\t'):
        self.zotu_table_file = pl.Path(zotu_table_file).resolve()
        self.sep = sep
        self.table_df = pd.read_csv(self.zotu_table_file, sep=sep, index_col=0)

    @property
    def tax_col_ind(self) -> int:
        cols_list = list(self.table_df.columns)
        for i, col in enumerate(cols_list):
            if 'taxonomy' in col.lower():
                return i
        return -1

    @property
    def zotus_ids(self) -> List[str]:
        return self.table_df.index.tolist()

    @property
    def sample_ids(self) -> List[str]:
        columns = list(self.table_df.columns)
        # tax column could be at any position
        if self.tax_col_ind != -1:
            columns.pop(self.tax_col_ind)
        return columns

    def get_zotu_counts(self, otu_id: str) -> List[int]:
        samples_ids = self.sample_ids
        counts = self.table_df.loc[otu_id, samples_ids].tolist()
        return counts

    def get_zotu_tax(self, otu_id: str) -> str:
        if self.tax_col_ind == -1:
            return ''
        otu_id_row_ind = self.table_df.index.get_loc(otu_id)
        return self.table_df.iloc[otu_id_row_ind, self.tax_col_ind]

    def collapse_zotu_groups(self, zotu_group: List[str]) -> Tuple[str, List[int], str]:
        """
        It combines the counts of the given zotus and returns the combined counts.
        """
        init_counts = [0 for _ in range(len(self.sample_ids))]
        group_tax = ''
        centroid = ''
        centroid_counts = init_counts
        group_counts = init_counts
        for zotu in zotu_group:
            otu_counts = self.get_zotu_counts(zotu)
            is_centroid = sum(otu_counts) >= sum(centroid_counts)
            group_counts = [sum(x) for x in zip(group_counts, otu_counts)]
            if is_centroid:
                group_tax = self.get_zotu_tax(zotu)
                centroid = zotu
                centroid_counts = otu_counts
        return centroid, group_counts, group_tax

    def sort_zotus_by_abundance(self, zotu_ids: list = None, reverse: bool = True) -> List[str]:
        if not zotu_ids:
            zotu_ids = self.zotus_ids
        zotu_counts = [sum(self.get_zotu_counts(zotu)) for zotu in zotu_ids]
        zotu_count_zip = list(zip(zotu_ids, zotu_counts))
        sorted(
            zotu_count_zip,
            key=lambda x: x[1],
            reverse=reverse
        )
        return [zotu for zotu, _ in zotu_count_zip]

    def total_count(self, zotu_ids: list = None) -> int:
        # We drop the taxonomy column. It could be at any position
        sampel_ids = self.sample_ids
        if zotu_ids is None:
            zotu_ids = self.zotus_ids
        else:
            assert all(zotu in self.zotus_ids for zotu in zotu_ids)
        total_count = self.table_df.loc[zotu_ids, sampel_ids].sum().sum()
        return int(total_count)


class TICUClust:
    """
    Objects of this class will take zOTUs and cluster them using UClust.
    """
    vsearch_bin = pl.Path(tic_configs["VSEARCH_BIN_PATH"]).resolve()
    def __init__(
            self,
            uclust_work_pd: pl.Path,
            vsearch_bin: pl.Path = None
            ):
        if vsearch_bin is None:
            self.vsearch_bin = pl.Path(tic_configs["VSEARCH_BIN_PATH"]).resolve()
        else:
            self.vsearch_bin = pl.Path(pl.PurePath(vsearch_bin)).absolute().resolve()

        if not self.vsearch_bin.exists() or not self.vsearch_bin.is_file():
            raise FileNotFoundError(f"Vsearch binary not found at {self.vsearch_bin}")
        if uclust_work_pd.is_file():
            raise FileNotFoundError(f"Vsearch work directory is a file at {uclust_work_pd}")
        uclust_work_pd.mkdir(parents=True, exist_ok=True)
        self.uclust_work_dir = uclust_work_pd

    def run_uclust(
            self,
            sequences: List[Sequence],
            cluster_threshold: float = 0.987,
            cut_tax: bool = False) -> Tuple[pl.Path, dict]:
        """
        It runs uclust_smallmem on the input fasta file and
        returns the uc file path and the centroids fasta file path

        :param input_fasta: str
            The path to the input fasta file
        :return: str
        """
         # create a temporary directory for UClust in the uclust_work_dir
        with tempfile.TemporaryDirectory(dir=self.uclust_work_dir) as temp_cluster_dir:
            run_dir = pl.Path(temp_cluster_dir).absolute().resolve()
        run_dir.mkdir(parents=True, exist_ok=True)
        sequence_cluster = SequenceCluster(sequences, force_homogeneity=False)
        # sequences with length below 32 will be discarded by vsearch
        if sequence_cluster.min_len < 32:
            logging.warning("Vsearch will discard sequences shorter than 32 nt.")
        input_fasta_path = run_dir.joinpath("input.fasta").absolute()
        sequence_cluster.write_to_fasta(input_fasta_path)
        sorted_seq_file = self.sort_seqs(str(input_fasta_path), by="length")
        centroids_file = input_fasta_path.parent.joinpath("cluster_centroids.fasta")
        uc_file = input_fasta_path.parent.joinpath("otu_clusters.uc")
        cmd_to_call = [
            str(self.vsearch_bin),
            "--cluster_smallmem",
            str(sorted_seq_file),
            "--centroids",
            str(centroids_file),
            "--uc",
            str(uc_file),
            "--id",
            str(cluster_threshold),
            "--strand",
            "both",
        ]
        system_sub(cmd_to_call, force_log=False, quiet=True)
        onelinefasta(centroids_file)
        centroid_cluster_dict = self.parse_uc_file(str(uc_file), cut_tax=cut_tax)
        return input_fasta_path, centroid_cluster_dict

    @staticmethod
    def parse_uc_file(
            uc_file: str,
            cut_tax: bool = False) -> dict:
        """
        It parsed uc files from -uclust_smallmem. It returns a dictionary with
        the centroid as the key and the members as the values. The centroid is
        always the first element in the list of members

        :param uc_file: str
            The path to the uc file
        :return: dict
        """
        uc_file = pl.Path(pl.PurePath(uc_file)).absolute()
        tax_reg = re.compile(r"\s(?P<tax_tag>tax=)?(?P<tax>([^;]+;)*([^;]+)?;?)", re.IGNORECASE)

        uc_dict = {}
        with open(uc_file, 'r', encoding='utf-8') as uc_h:
            line = uc_h.readline().strip()
            while line:
                if line.startswith("H"):
                    line_tokens = line.split("\t")
                    query = line_tokens[8]
                    target_centroid = line_tokens[9]
                    sim_ = round(float(line_tokens[3]), 5)
                    # omitting the taxonomic information
                    if cut_tax:
                        query = tax_reg.sub("", query).strip()
                        target_centroid = tax_reg.sub("", target_centroid).strip()
                    # prepending '>' to the centroid and query.
                    # As uclust removes the '>' from the query
                    target_centroid = '>' + target_centroid
                    query = '>' + query
                    matched_clusters = uc_dict.get(query, [])
                    matched_clusters.append((target_centroid, sim_))
                    uc_dict[query] = matched_clusters
                elif line.startswith("S"):
                    # case for singletons
                    line_tokens = line.split("\t")
                    centroid = line_tokens[8]
                    sim_ = round(float(100), 5)
                    # omitting the taxonomic information
                    if cut_tax:
                        centroid = tax_reg.sub("", centroid).strip()
                    # prepending '>' to the centroid
                    centroid = '>' + centroid
                    query = centroid
                    matched_clusters = uc_dict.get(query, [])
                    matched_clusters.append((centroid, sim_))
                    uc_dict[query] = matched_clusters
                line = uc_h.readline().strip()
        # some sequences might match more than one more centers
        # we will keep the one with the highest similarity
        cent_members_dict = {}
        for member, cent_matches in uc_dict.items():
            cent_matches = sorted(cent_matches, key=lambda x: x[1], reverse=True)
            cent_id = cent_matches[0][0]
            mems_list = cent_members_dict.get(cent_id, [])
            mems_list.append(member)
            cent_members_dict[cent_id] = mems_list

        return cent_members_dict

    def get_sequences_clusters(
            self,
            sequences: List[Sequence],
            cluster_threshold: float) -> List[SequenceCluster]:
        clusters = []
        input_fasta_path, cent_clust_dict = self.run_uclust(
            sequences,
            cluster_threshold,
            cut_tax=True
        )
        input_fasta = TaxedFastaFile(input_fasta_path)
        for _, seqs_id in cent_clust_dict.items():
            # First element in the seqs_id is the centroid
            seqs_to_get = [str(SeqHeader(seq_id).seq_id) for seq_id in seqs_id]
            seq_objs = input_fasta.get_seq_by_seq_id(seqs_to_get)
            if not seq_objs:
                raise ValueError(f"No sequence found for {seqs_to_get}")
            # pop the centroid
            centroid = seq_objs[0]
            last_knwon_tax_level = centroid.header.taxonomy.last_known_level
            cluster = SequenceCluster(seq_objs, centroid, last_knwon_tax_level)
            clusters.append(cluster)
        # deleting the intermediate files
        shutil.rmtree(input_fasta_path.parent)
        return clusters

    def sort_seqs(self, to_sort_fasta: str, by: str = 'size') -> str:
        to_sort_path = pl.Path(to_sort_fasta).resolve()
        dir_path = to_sort_path.parent
        sorted_fasta_file = dir_path / f"sorted_{by}.fasta"
        if by == 'size':
            cmd_to_call_list = [
                str(self.vsearch_bin),
                "--sortbysize",
                str(to_sort_path),
                "--output",  # vsearch equivalent of usearch --fastaout
                str(sorted_fasta_file),
            ]
            system_sub(cmd_to_call_list, force_log=False, quiet=True)
        elif by == 'length':
            cmd_to_call_list = [
                str(self.vsearch_bin),
                "--sortbylength",
                str(to_sort_path),
                "--output",  # vsearch equivalent of usearch --fastaout
                str(sorted_fasta_file),
            ]
            system_sub(cmd_to_call_list, force_log=False, quiet=True)
        else:
            raise ValueError("Sorting should be by either 'size' or 'length'")
        onelinefasta(sorted_fasta_file)
        return sorted_fasta_file


class TICAnalysis:

    threads = 1
    default_thresholds = OrderedDict({
        'kingdom': None,
        'phylum': None,
        'class': None,
        'order': None,
        'family': 0.90,
        'genus': 0.95,
        'species': 0.987,
    })
    default_output_fasta_name = "TIC-FullTaxonomy-Sequences.fasta"
    default_non_bact_fasta_name = "Non-Bacteria-Sequences.fasta"
    default_fotu_gotu_file_name = "Map-FOTU-GOTU.tab"
    default_gotu_sotu_file_name = "Map-GOTU-SOTU.tab"
    default_sotu_zotu_file_name = "Map-SOTU-ZOTU.tab"
    default_sotu_table_file_name = "SOTUs-Table.tab"
    default_sotu_fasta_file_name = "SOTUs-Seq.fasta"
    default_work_dir_name = "TIC-WD"
    default_uclust_work_dir_name = "Uclust-WD"

    def __init__(
        self,
        taxed_fasta_file_path: pl.Path,
        zotu_table_file: pl.Path = None):
        taxed_fasta_file_path = pl.Path(taxed_fasta_file_path).resolve()
        self.fasta_file = TaxedFastaFile(taxed_fasta_file_path)
        self.zotu_table_file = pl.Path(zotu_table_file).resolve() if zotu_table_file else None
        self.zotu_table: ZOTUTable = ZOTUTable(self.zotu_table_file) if zotu_table_file else None
        if self.zotu_table and not self.table_fasta_ids_match:
            raise ValueError("Sequence IDs zOTU table and the fasta file do not match.")
        self.tic_wd = self.fasta_file.fasta_file_path.parent / self.default_work_dir_name
        if self.tic_wd.exists():
            logging.warning(
                "TIC work directory already exists. Running TIC overwrites it."
            )
            shutil.rmtree(self.tic_wd)
        self.uclust_wd = self.tic_wd / self.default_uclust_work_dir_name

        # creating the directories
        self.tic_wd.mkdir(parents=True, exist_ok=True)
        self.uclust_wd.mkdir(parents=True, exist_ok=True)
        # creating the file paths
        self.non_bact_fasta_path = self.tic_wd / self.default_non_bact_fasta_name
        self.fotu_gotu_file_path = self.tic_wd / self.default_fotu_gotu_file_name
        self.gotu_sotu_file_path = self.tic_wd / self.default_gotu_sotu_file_name
        self.sotu_zotu_file_path = self.tic_wd / self.default_sotu_zotu_file_name
        self.tic_output_fasta_path = self.tic_wd / self.default_output_fasta_name
        self.sotu_table_path = self.tic_wd / self.default_sotu_table_file_name
        self.sotu_fasta_path = self.tic_wd / self.default_sotu_fasta_file_name
        # creating the lists
        self.fotu_gotu_list: List[Tuple[str, str]] = []
        self.gotu_sotu_list: List[Tuple[str, str]] = []
        self.sotu_zotu_list: List[Tuple[str, str]] = []
        # setting the default thresholds
        self.cluster_thresholds = self.default_thresholds.copy()

    def filter_tax_set_at_last_known_level(self, level: str) -> List[Taxonomy]:
        this_level_known_tax = self.fasta_file.filter_tax_set_at_last_known_level(level)
        this_level_known_tax = [
            tax for tax in this_level_known_tax
            if tax.kingdom == 'Bacteria'
        ]
        this_level_known_tax = list(set(this_level_known_tax))
        return this_level_known_tax

    def run(
            self,
            threads: int = int(cpu_count() * 0.75),
            cluster_thresholds_d: Dict[str, float] = None
            ) -> None:
        self.threads = threads
        if cluster_thresholds_d:
            self.cluster_thresholds.update(cluster_thresholds_d)
        logging.debug("Running TIC with the following thresholds: %s", self.cluster_thresholds)
        all_known_order_fasta_path = self.fill_upto_order()
        all_known_family_fasta_path = self.complete_family_level(all_known_order_fasta_path)
        all_known_genus_fasta_path = self.complete_genus_level(all_known_family_fasta_path)
        all_known_species_fasta_path = self.complete_species_level(all_known_genus_fasta_path)
        # self.tic_output_fasta = TaxedFastaFile(all_known_species_fasta_path)
        # completing maps
        self.write_fotu_gotu_map_to_file(all_known_genus_fasta_path)
        self.write_gotu_sotu_map_to_file(all_known_species_fasta_path)
        self.write_sotu_zotu_map_to_file(all_known_species_fasta_path)
        # writing the sotu table to file. Also SOTU to ZOTU map
        if self.zotu_table:
            # We rely on ZOTU table to infer the centroids
            # and the members of the SOTUs
            self.deflate_zotu_table(all_known_species_fasta_path)
            # writing the sotu fasta to file
            self.extract_sotu_fasta_file(all_known_species_fasta_path)
        # writing non-bacteria sequences to the output fasta
        self.append_non_bacteria_seqs(self.non_bact_fasta_path)
        # deleting the intermediate files
        all_known_order_fasta_path.unlink()
        all_known_family_fasta_path.unlink()
        all_known_genus_fasta_path.unlink()
        # log the output fasta file path
        logging.info("Output fasta wrote to %s", all_known_species_fasta_path)

        shutil.rmtree(self.uclust_wd, ignore_errors=True)
        return all_known_species_fasta_path

    def grow_taxonomy(
            self,
            tax_to_complete: Taxonomy,
            input_fasta: pl.Path,
            cluster_threshold: float = 0.987
            ) -> List[SequenceCluster]:
        input_fasta = TaxedFastaFile(pl.Path(input_fasta).absolute().resolve())
        uclust_obj = TICUClust(uclust_work_pd=self.uclust_wd)
        homo_level = tax_to_complete.last_known_level
        tax_seqs = input_fasta.filter_seq_by_tax(tax_to_complete)
        sequence_cluster = SequenceCluster(
            tax_seqs,
            homogeneity_level=homo_level
        )

        sub_clusters_list = uclust_obj.get_sequences_clusters(
            list(sequence_cluster),
            cluster_threshold
        )
        return sub_clusters_list

    def __get_fasta_file(self, file_stem_name: str, make_unique: bool = True) -> TaxedFastaFile:
        # Add some random parts to file_stem_name so that it does not overwrite other file
        if make_unique:
            file_stem_name = next(tempfile._get_candidate_names()) + "_" + file_stem_name
        target_fasta_file_path = self.tic_wd / f"{str(file_stem_name)}.fasta"
        # deleting the file if it exists
        if target_fasta_file_path.exists():
            target_fasta_file_path.unlink()
        target_fasta_file_path.touch()
        target_fasta_file_obj = TaxedFastaFile(target_fasta_file_path)
        return target_fasta_file_obj

    @property
    def table_fasta_ids_match(self):
        if self.zotu_table:
            # zotus in zotu_table and sequence IDs in fasta file should match
            zotu_ids = set(self.zotu_table.zotus_ids)
            fasta_ids = set(self.fasta_file.get_seq_ids())
            if zotu_ids == fasta_ids:
                return True
        return False

    def filter_bac_seq_last_kown_at(self, level: str, flatten: bool = False) -> List[Sequence]:
        last_known_ = self.filter_tax_set_at_last_known_level(level)
        last_known_bac = [
            tax for tax in last_known_
            if tax.kingdom == 'Bacteria'
        ]
        last_known_bac = list(set(last_known_bac))
        seqs = []
        for tax in last_known_bac:
            seqs.append(self.fasta_file.filter_seq_by_tax(tax))
        seqs = [seq_list for seq_list in seqs if seq_list]
        if flatten:
            seqs = [seq for seq_list in seqs for seq in seq_list]
        return seqs

    def fill_upto_order(self) -> pl.Path:
        """
        It completes the taxonomy upto order level for the sequences known upto kingdom level.
        NOTE: returned fasta file will have all the bacterial sequences filled upto order level
        """
        last_known_kingdoms = self.filter_bac_seq_last_kown_at('kingdom')
        clusters_to_fill = []
        for known_king_seq in last_known_kingdoms:
            seq_cluster = KingdomCluster(
                known_king_seq
            )
            kingdom_ = seq_cluster.centroid.header.taxonomy.kingdom
            seq_cluster.set_level('phylum', 'NA_phylum__' + kingdom_)
            seq_cluster.set_level('class', 'NA_class__' + kingdom_)
            seq_cluster.set_level('order', 'NA_order__' + kingdom_)
            clusters_to_fill.append(seq_cluster)

        last_known_pyl_bacteria = self.filter_bac_seq_last_kown_at('phylum')
        for known_phyl_seq in last_known_pyl_bacteria:
            seq_cluster = PhylumCluster(
                known_phyl_seq
            )
            phylum_ = seq_cluster.centroid.header.taxonomy.phylum
            seq_cluster.set_level('class', 'NA_class__' + phylum_)
            seq_cluster.set_level('order', 'NA_order__' + phylum_)
            clusters_to_fill.append(seq_cluster)

        last_known_class_bacteria = self.filter_bac_seq_last_kown_at('class')
        for known_class_seq in last_known_class_bacteria:
            seq_cluster = ClassCluster(
                known_class_seq
            )
            class_ = seq_cluster.centroid.header.taxonomy.class_
            seq_cluster.set_level('order', 'NA_order__' + class_)
            clusters_to_fill.append(seq_cluster)

        all_knwon_order_fasta = self.__get_fasta_file("All-Known-Order")
        seq_list_to_write = [
            seq_ for seq_cluster in clusters_to_fill for seq_ in seq_cluster
        ]
        all_knwon_order_fasta.write_to_fasta_file(sequences=seq_list_to_write, mode='a')
        # we append all originally known sequences upto order level to the output fasta
        last_known_orders_from_input_fasta = self.filter_bac_seq_last_kown_at('order', flatten=True)
        all_knwon_order_fasta.write_to_fasta_file(
            sequences=last_known_orders_from_input_fasta,
            mode='a'
        )

        return all_knwon_order_fasta.fasta_file_path

    def complete_family_level(self, all_known_order_fasta_path: pl.Path) -> pl.Path:
        """
        It completes the taxonomy upto family level for the sequences known upto order level.
        NOTE: All the sequences in the input fasta file should be known upto order level and 
        should be from 'Bacteria' kingdom
        """
        all_known_order_fasta = TaxedFastaFile(all_known_order_fasta_path)
        # Due to filter_bac_seq_last_kown_at, we are sure that the last known level is 'order'
        # and all the sequences are from 'Bacteria' kingdom
        last_known_order_taxs = all_known_order_fasta.filter_tax_set_at_last_known_level('order')
        # Below is for getting sure that we are working with only bacteria
        last_known_ord_bacteria = [
            tax for tax in last_known_order_taxs
            if tax.kingdom == 'Bacteria'
        ]
        args_list = [
            (tax, all_known_order_fasta.fasta_file_path, self.cluster_thresholds['family'])
            for tax in last_known_ord_bacteria
        ]
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            result_clusters = list(
                executor.map(
                    lambda args_tup: self.grow_taxonomy(*args_tup),
                    args_list
                )
            )
            # make sure that each thread is done
            for _ in result_clusters:
                pass
        all_known_family_fasta = self.__get_fasta_file("All-Known-Family")
        all_known_family_seq_to_write = []
        for ind, clusters_list in enumerate(result_clusters):
            ind += 1
            for inner_ind, cluster in enumerate(clusters_list):
                inner_ind += 1
                cluster.set_level('family', 'FOTU' + self.concat_nums(ind, inner_ind))
                all_known_family_seq_to_write += list(cluster)

        all_known_family_fasta.write_to_fasta_file(
            sequences=all_known_family_seq_to_write,
            mode='a'
        )

        logging.info("Completed family level for %s orders.", str(len(result_clusters)))
        # we append all originally known sequences upto family level to the output fasta
        last_known_families_from_input_fasta = self.filter_bac_seq_last_kown_at(
            'family',
            flatten=True
        )
        all_known_family_fasta.write_to_fasta_file(
            sequences=last_known_families_from_input_fasta,
            mode='a'
        )

        return all_known_family_fasta.fasta_file_path

    def complete_genus_level(self, all_known_family_fasta_path: pl.Path) -> pl.Path:
        """
        It completes the taxonomy upto genus level for the sequences known upto family level.
        NOTE: All the sequences in the input fasta file should be known upto family level and
        should be from 'Bacteria' kingdom
        """
        all_known_family_fasta = TaxedFastaFile(all_known_family_fasta_path)
        # Due to filter_bac_seq_last_kown_at, we are sure that the last known level is 'family'
        # and all the sequences are from 'Bacteria' kingdom
        last_known_family_taxs = all_known_family_fasta.filter_tax_set_at_last_known_level('family')
        # Below is for getting sure that we are working with only bacteria
        last_known_fam_bacteria = [
            tax for tax in last_known_family_taxs
            if tax.kingdom == 'Bacteria'
        ]
        args_list = [
            (tax, all_known_family_fasta.fasta_file_path, self.cluster_thresholds['genus'])
            for tax in last_known_fam_bacteria
        ]
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            result_clusters = list(
                executor.map(
                    lambda args_tup: self.grow_taxonomy(*args_tup),
                    args_list
                )
            )
            # make sure that each thread is done
            for _ in result_clusters:
                pass
        all_known_genus_fasta = self.__get_fasta_file("All-Known-Genus")
        all_known_genus_to_write = []
        for ind, clusters_list in enumerate(result_clusters):
            ind += 1
            for inner_ind, cluster in enumerate(clusters_list):
                inner_ind += 1
                gotu_id = 'GOTU' + self.concat_nums(ind, inner_ind)
                cluster.set_level('genus', gotu_id)
                all_known_genus_to_write += list(cluster)
        all_known_genus_fasta.write_to_fasta_file(sequences=all_known_genus_to_write, mode='a')

        logging.info("Completed genus level for %s families.", str(len(result_clusters)))
        # we append all originally known sequences upto genus level to the output fasta
        last_known_genus_from_input_fasta = self.filter_bac_seq_last_kown_at('genus', flatten=True)
        all_known_genus_fasta.write_to_fasta_file(
            sequences=last_known_genus_from_input_fasta,
            mode='a'
        )
        self.fotu_gotu_list = self.get_fotu_gotu_map(all_known_family_fasta.fasta_file_path)

        return all_known_genus_fasta.fasta_file_path

    def complete_species_level(self, all_known_genus_fasta_path: pl.Path) -> pl.Path:
        """
        It completes the taxonomy upto species level for the sequences known upto genus level.
        NOTE: All the sequences in the input fasta file should be known upto genus level and
        should be from 'Bacteria' kingdom
        """
        all_known_genus_fasta = TaxedFastaFile(all_known_genus_fasta_path)
        # Due to filter_bac_seq_last_kown_at, we are sure that the last known level is 'genus'
        # and all the sequences are from 'Bacteria' kingdom
        last_known_genus_taxs = all_known_genus_fasta.filter_tax_set_at_last_known_level('genus')
        # Below is for getting sure that we are working with only bacteria
        last_knwon_gen_bacteria = [
            tax for tax in last_known_genus_taxs
            if tax.kingdom == 'Bacteria'
        ]
        args_list = [
            (tax, all_known_genus_fasta.fasta_file_path, self.cluster_thresholds['species'])
            for tax in last_knwon_gen_bacteria
        ]
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            result_clusters = list(
                executor.map(
                    lambda args_tup: self.grow_taxonomy(*args_tup),
                    args_list
                )
            )
            # make sure that each thread is done
            for _ in result_clusters:
                pass
        all_known_species_fasta = TaxedFastaFile(self.tic_output_fasta_path)
        all_known_speceis_to_write = []
        for ind, clusters_list in enumerate(result_clusters):
            ind += 1
            for inner_ind, cluster in enumerate(clusters_list):
                inner_ind += 1
                sotu_id = 'SOTU' + self.concat_nums(ind, inner_ind)
                cluster.set_level('species', sotu_id)
                all_known_speceis_to_write += list(cluster)

        all_known_species_fasta.write_to_fasta_file(sequences=all_known_speceis_to_write, mode='a')
        logging.info("Completed species level for %s genera.", str(len(result_clusters)))

        # finally we append all the sequences known upto species level to the output fasta
        last_known_species_from_input_fasta = self.filter_bac_seq_last_kown_at(
            'species',
            flatten=True
        )
        all_known_species_fasta.write_to_fasta_file(
            sequences=last_known_species_from_input_fasta,
            mode='a'
        )
        # getting the maps
        self.gotu_sotu_list = self.get_gotu_sotu_map(all_known_species_fasta.fasta_file_path)
        self.sotu_zotu_list = self.get_sotu_zotu_map(all_known_species_fasta.fasta_file_path)

        return all_known_species_fasta.fasta_file_path

    def write_fotu_gotu_map_to_file(self, fasta_file_path: pl.Path) -> None:
        map_list = self.get_fotu_gotu_map(fasta_file_path)
        with open(self.fotu_gotu_file_path, 'w', encoding='utf-8') as map_file_io:
            map_file_io.write("FOTU\tGOTU\n")
            for fotu, gotu in map_list:
                map_file_io.write(f"{fotu}\t{gotu}\n")

    def write_gotu_sotu_map_to_file(self, fasta_file_path: pl.Path) -> None:
        map_list = self.get_gotu_sotu_map(fasta_file_path)
        with open(self.gotu_sotu_file_path, 'w', encoding='utf-8') as map_file_io:
            map_file_io.write("GOTU\tSOTU\n")
            for gotu, sotu in map_list:
                map_file_io.write(f"{gotu}\t{sotu}\n")

    def write_sotu_zotu_map_to_file(self, fasta_file_path: pl.Path) -> None:
        map_list = self.get_sotu_zotu_map(fasta_file_path)
        map_dict = self.__species_zotus_pair_to_dict(map_list)
        if self.zotu_table:
            # if there is a zotu_table given to TIC then the zotus
            # will be sorted by abundance.
            map_dict = {
                sotu: self.zotu_table.sort_zotus_by_abundance(zotus)
                for sotu, zotus in map_dict.items()
            }
        with open(self.sotu_zotu_file_path, 'w', encoding='utf-8') as map_file_io:
            map_file_io.write("SOTU\tZOTU\n")
            for sotu, zotu_li in map_dict.items():
                for zotu in zotu_li:
                    map_file_io.write(f"{sotu}\t{zotu}\n")

    @staticmethod
    def get_fotu_gotu_map(taxed_fasta_path: pl.Path) -> List[Tuple[str, str]]:
        taxed_fasta_obj = TaxedFastaFile(taxed_fasta_path)
        return taxed_fasta_obj.get_map_to_parent('family', 'genus')

    @staticmethod
    def get_gotu_sotu_map(taxed_fasta_path: pl.Path) -> List[Tuple[str, str]]:
        taxed_fasta_obj = TaxedFastaFile(taxed_fasta_path)
        return taxed_fasta_obj.get_map_to_parent('genus', 'species')

    @staticmethod
    def get_sotu_zotu_map(taxed_fasta_path: pl.Path) -> List[Tuple[str, str]]:
        """
        Given a fasta file path, it returns a list of tuples of SOTU and ZOTU.
        NOTE if value at species level is not uniuqe, in the whole taxonomic
        tree then we end of with multiple ZOTUs for a single SOTU. For example,
        if there are two species with the same name but different genus, then
        we will have two ZOTUs for a single species.
        """
        map_set = set()
        with open(taxed_fasta_path, 'r', encoding='utf-8') as fasta:
            curr_line = fasta.readline().strip()
            while curr_line:
                if curr_line.startswith('>'):
                    seq_header = SeqHeader(curr_line)
                    if seq_header.taxonomy.is_known_upto('species'):
                        map_set.add((str(seq_header.taxonomy.species), str(seq_header.seq_id)))
                    else:
                        logging.warning(
                            "Sequence %s is not known upto species level."\
                            " It won't be included in SOTU/ZOTU map.",
                            str(seq_header.seq_id)
                        )
                curr_line = fasta.readline().strip()
        map_list = list(map_set)
        sorted(map_list, key=lambda x: x[0])
        return map_list

    def append_non_bacteria_seqs(self, fasta_file_path: pl.Path) -> None:
        non_bacteria_seqs = []
        non_bacteria_seqs = self.fasta_file.filter_seq_by_level_value('kingdom', 'Archaea')
        non_bacteria_seqs += self.fasta_file.filter_seq_by_level_value('kingdom', 'Eukaryota')
        non_bacteria_seqs += self.fasta_file.filter_seq_by_level_value('kingdom', 'Viruses')
        non_bacteria_seqs += self.fasta_file.filter_seq_by_level_value('kingdom', 'Viroids')
        non_bacteria_seqs += self.fasta_file.filter_seq_by_level_value('kingdom', 'NA-Kingdom')
        logging.info(
            "Appending %s non-bacteria sequences to the output fasta.", str(len(non_bacteria_seqs))
        )
        with open(fasta_file_path, 'a', encoding='utf-8') as fasta:
            for seq in non_bacteria_seqs:
                fasta.write(str(seq) + '\n')

    @staticmethod
    def generate_id(input_str: str, prefix: str, length: int = 4) -> str:
        """
        It generates a unique id for the given string.
        NOTE There is a risk of hash collision.
        """
        hash_object = hashlib.sha256(str(input_str).lower().encode())
        hash_hex = hash_object.hexdigest()
        # Take the first 'length' digits of the hash
        short_hash_str = hash_hex[:length]
        # Convert the short hash to an integer
        short_hash_int = int(short_hash_str, 16)

        return str(prefix) + str(short_hash_int)

    def __species_zotus_pair_to_dict(self, sotu_zotu_map_pairs: List[Tuple[str, str]]) -> dict:
        sotu_to_zotu_map = {}
        for sotu, zotu in sotu_zotu_map_pairs:
            zotus = sotu_to_zotu_map.get(sotu, [])
            zotus.append(zotu)
            sotu = sotu if "SOTU" in sotu else self.generate_id(sotu, "SOTU0")
            sotu_to_zotu_map[sotu] = zotus
        return sotu_to_zotu_map

    def get_sotu_table_lines(
        self,
        sotu_zotu_map: List[Tuple[str, str]],
        known_species_fasta: pl.Path) -> Dict[str, Tuple[str, List[int], str]]:
        """
        It returns the lines of the SOTU table. Each line is a tuple of 
         SOTU, (centroid ZOTU, counts, taxonomy)
        """
        if self.zotu_table is None:
            logging.warning("ZOTU table is not provided. Skipping the SOTU table generation.")
            return {}

        known_species_fasta = pl.Path(known_species_fasta).resolve()
        known_species_fasta = TaxedFastaFile(known_species_fasta)
        sotu_to_zotu_map = self.__species_zotus_pair_to_dict(sotu_zotu_map)
        sotu_to_line_map = {}
        for sotu, zotus in sotu_to_zotu_map.items():
            collapsed_zotu_line = self.zotu_table.collapse_zotu_groups(zotus)
            sotu_to_line_map[sotu] = collapsed_zotu_line
        new_sotu_to_line_map = {}
        for sotu, (cent_z, counts, _) in sotu_to_line_map.items():
            new_sotu_to_line_map[cent_z] = (
                sotu,
                counts,
                known_species_fasta.get_seq_by_seq_id([cent_z])[0].header.taxonomy.tax_str
            )
        sotu_to_line_map = new_sotu_to_line_map
        return sotu_to_line_map

    def deflate_zotu_table(
        self,
        known_species_fasta: pl.Path) -> None:
        """
        It deflates the zOTU table to the SOTU table. SOTUs might be actual species or SOTUs
        """
        known_species_fasta = pl.Path(known_species_fasta).resolve()
        if self.zotu_table is None or not known_species_fasta.is_file():
            logging.warning(
                "ZOTU table or known species fasta file is not provided. Skipping the deflation."
            )
            return None

        species_zotus_map: List[Tuple[str, str]] = self.get_sotu_zotu_map(known_species_fasta)
        sotu_zotu_map = species_zotus_map
        sotu_to_line_map = self.get_sotu_table_lines(sotu_zotu_map, known_species_fasta)
        sotu_df = pd.DataFrame(
            columns= ["#SOTU"] + self.zotu_table.sample_ids + ["Taxonomy"]
        )
        for _, (sotu, counts, tax) in sotu_to_line_map.items():
            # write the line to the dataframe
            sotu_df.loc[sotu] = [sotu] + counts + [tax]
        if not sotu_df.empty:
            sotu_df.to_csv(self.sotu_table_path, sep='\t', index=False)
        else:
            logging.warning("No SOTUs found. The SOTU table won't be generated.")

    def extract_sotu_fasta_file(
        self,
        known_species_fasta: pl.Path) -> None:
        """
        It extracts the sequences of the SOTUs from the source fasta file and
         writes them to a new file.
        """
        seq_list = []
        known_species_fasta = pl.Path(known_species_fasta).resolve()
        known_species_fasta_obj = TaxedFastaFile(known_species_fasta)
        species_zotus_map: List[Tuple[str, str]] = self.get_sotu_zotu_map(known_species_fasta)
        sotu_zotu_map = species_zotus_map
        sotu_to_line_map = self.get_sotu_table_lines(sotu_zotu_map, known_species_fasta)
        for cent_z, (sotu, _, _) in sotu_to_line_map.items():
            curr_seq = known_species_fasta_obj.get_seq_by_seq_id([cent_z])[0]
            curr_seq.header.seq_id.id_str = sotu
            seq_list.append(curr_seq)

        self.fasta_file.write_to_fasta_file(
            sequences=seq_list,
            mode='w',
            output_file_path=self.sotu_fasta_path
            )

    @staticmethod
    def concat_nums(num_1: int, num_2: int) -> str:
        """Encodes two numbers into a single unique number."""
        assert num_1 > 0 and num_2 > 0, "Both numbers should be non-negative!"
        num_1_str = str(num_1)
        num_2_str = str(num_2)
        return f"{num_1_str}0{num_2_str}"

    def __del__(self):
        shutil.rmtree(self.uclust_wd, ignore_errors=True)

    def cleanup(self, full: bool = False):
        if full:
            shutil.rmtree(self.tic_wd, ignore_errors=True)
        shutil.rmtree(self.uclust_wd, ignore_errors=True)

    def __enter__(self):
        return self

    # when the object is deleted, the temporary directory should be deleted
    def __exit__(self, exc_type, exc_value, traceback) -> None:
        full_cleanup = False
        if exc_type is not None:
            logging.error("An error occurred while running TIC. Error: %s", exc_value)
            full_cleanup = True
        self.cleanup(full=full_cleanup)
