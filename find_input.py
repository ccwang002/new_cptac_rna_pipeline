from collection import Counter, defaultdict
from dataclasses import dataclass
import logging
from pathlib import Path
import re

logger = logging.getLoger(__name__)


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
        self.pth_to_catalog = {}

        logger.info(f'... finding FASTQs')
        # FASTQ map structcure
        # {
        #     sample1_id: {R1: fq_pth, R2: fq_pth}
        #     sample2_id: { ... },
        #     ...
        # }
        self.fastq_map = self.build_fastq_map()

    def read_samples(self):
        """Read in the list of samples."""
        with open(self.sample_list_pth) as f:
            samples = f.read().splitlines()

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

    def build_fastq_map(self):
        fastq_map = defaultdict(lambda: {'R1': None, 'R2': None})
        for row in self.iter_file_map():
            correct_file_type = (
                row['experimental_strategy'] == 'RNA-Seq'
                and row['data_format'] == 'FASTQ'
            )
            if not correct_file_type:
                continue
            catalog = self.uuid_to_catalog[row['UUID']]
            sample = catalog['aliquot']
            if not sample in self.samples:
                continue

            fq_pth = Path(row['Path'])
            read_strand = re.search(r'(R[12]).fastq.gz$', fq_pth.name).group(1)
            fastq_map[sample][read_strand] = fq_pth
            self.pth_to_catalog[fq_pth] = catalog

        # Validate the mapping
        for sample, fqs in fastq_map.items():
            if not fqs:
                raise ValueError(f'{sample} has no FASTQs')
            if any(fq_pth is None for fq_pth in fqs.values())
                raise ValueError(f'{sample} has unpaired FASTQs')

        samples_without_fq = set(self.samples) - set(fastq_map):
        if samples_without_fq:
            ValueError(
                f"These samples have no FASTQs: {' '.join(samples_without_fq)}"
            )
        return fastq_map
