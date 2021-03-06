from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
import gzip
import logging
from pathlib import Path
import re

logger = logging.getLogger(__name__)


def compression_aware_open(pth):
    pth = Path(pth)
    if pth.suffix == '.gz':
        return gzip.open(pth, 'rt')
    else:
        return open(pth)


class FileMapping:
    def __init__(self, file_map_pth, catalog_pth, sample_list_pth, pth_replace_patterns = None):
        self.file_map_pth = file_map_pth
        self.catalog_pth = catalog_pth
        self.sample_list_pth = sample_list_pth
        self.pth_replace_patterns = pth_replace_patterns

        logger.info(f'... reading sample list')
        self.samples = self.read_samples()

        logger.info(f'... read catalog')
        self.uuid_to_catalog = self.read_catalog()

        # A reverse mapping from file path to its details from the catalog
        self.pth_to_catalog = {}

    def read_samples(self):
        """Read in the list of samples."""
        with open(self.sample_list_pth) as f:
            samples = [x for x in f.read().splitlines() if x]

        unique_samples = set(samples)
        # Find sample duplications
        if len(samples) != len(unique_samples):
            sample_counter = Counter(samples)
            dup_samples = sorted(
                s for s, c in sample_counter.items()
                if c > 1
            )
            logger.error(f"Found duplicated samples: {dup_samples}")
            raise ValueError("Duplicated samples in the given sample list")
        return unique_samples

    def read_catalog(self):
        uuid_to_catalog = {}
        with compression_aware_open(self.catalog_pth) as f:
            reader = csv.DictReader(f, dialect='excel-tab')
            for row in reader:
                sample = row['aliquot']
                if sample in self.samples:
                    uuid_to_catalog[row['UUID']] = row
        return uuid_to_catalog

    def iter_file_map(self):
        """Read file map and only return the matched samples."""
        with open(self.file_map_pth) as f:
            reader = csv.DictReader(f, dialect='excel-tab')
            for row in reader:
                yield row


class FASTQFileMapping(FileMapping):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        logger.info(f'... finding FASTQs')
        # FASTQ map structcure
        # {
        #     sample1_id: {R1: fq_pth, R2: fq_pth}
        #     sample2_id: { ... },
        #     ...
        # }
        self.fastq_map = self.build_fastq_map()

    def build_fastq_map(self):
        fastq_map = defaultdict(lambda: {'R1': None, 'R2': None})
        for row in self.iter_file_map():
            # Samples not selected will not be in the mapping
            if row['UUID'] not in self.uuid_to_catalog:
                continue
            correct_file_type = (
                row['experimental_strategy'] == 'RNA-Seq'
                and row['data_format'] == 'FASTQ'
            )
            if not correct_file_type:
                continue
            catalog = self.uuid_to_catalog[row['UUID']]
            sample = catalog['aliquot']

            fq_pth = Path(row['data_path'])
            read_strand = re.search(r'_(R[12])_\d+\.fastq\.gz$', fq_pth.name).group(1)
            fastq_map[sample][read_strand] = fq_pth
            self.pth_to_catalog[fq_pth] = catalog

        # Validate the mapping
        for sample, fqs in fastq_map.items():
            if not fqs:
                raise ValueError(f'{sample} has no FASTQs')
            if any(fq_pth is None for fq_pth in fqs.values()):
                raise ValueError(f'{sample} has unpaired FASTQs')

        samples_without_fq = set(self.samples) - set(fastq_map)
        if samples_without_fq:
            ValueError(
                f"These samples have no FASTQs: {' '.join(samples_without_fq)}"
            )
        return fastq_map


class BAMFileMapping(FileMapping):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        logger.info(f'... finding BAMs')
        # BAM map structcure
        # {
        #     sample1_id: {'genomic': bam_pth, 'transcriptome': bam_pth},
        #     sample2_id: { ... },
        #     ...
        # }
        self.bam_map = self.build_bam_map()

    def build_bam_map(self):
        bam_map = defaultdict(dict)
        for row in self.iter_file_map():
            # Samples not selected will not be in the mapping
            if row['UUID'] not in self.uuid_to_catalog:
                continue

            catalog = self.uuid_to_catalog[row['UUID']]
            sample = catalog['aliquot']
            correct_file_type = (
                row['experimental_strategy'] == 'RNA-Seq'
                and row['data_format'] == 'BAM'
                and row['reference'] == 'hg38'
                and catalog['result_type'] in ['genomic', 'transcriptome']
            )
            if not correct_file_type:
                continue

            pth = Path(row['data_path'])
            bam_map[sample][catalog['result_type']] = pth
            self.pth_to_catalog[pth] = catalog

        # Validate the mapping
        samples_without_bam = set(self.samples) - set(bam_map)
        if samples_without_bam:
            ValueError(
                f"These samples have no BAMs: {' '.join(samples_without_bam)}"
            )
        return bam_map