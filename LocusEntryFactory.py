from pdb import set_trace as ST
from pprint import pprint as PP

from LocusEntry import LocusEntry
from IlluminaBeadArrayFiles import RefStrand

class LocusEntryFactory(object):
    """
    Class to create locus entries from BPM records
    """

    def __init__(self, vcf_record_factory, chrom_order, unsquash_duplicates, split_multiallelics, logger):
        """
        Create new locus entry factory

        Args:
            vcf_record_factory (VcfRecordFactory): Creates vcf._Record objects
            chrom_sort_function (func(string, int)): Function used to sort chromosomes
            unsquash_duplicates (bool): True to generate separate entries for duplicates
            logger (logging.Logger): Logger to report warnings/errors
        """
        self._vcf_record_factory = vcf_record_factory
        self._chrom_order = chrom_order
        self._unsquash_duplicates = unsquash_duplicates
        self._split_multiallelics = split_multiallelics
        self._logger = logger

    def create_locus_entries(self, bpm_reader):
        """
        Generate locus entries from BPM records without
        any sample information.

        Args:
            bpm_reader (BPMReader): Provides BPM records
            loci_to_filter (set(string)): Set of loci names to filter from the manifest

        Returns:
            list(LocusEntry): List of locus entries corresponding to BPM records
        """
        result = []
        for record_group in self._group_bpm_records(bpm_reader.get_bpm_records()):
            result.append(self._generate_locus_entry(record_group))
        return sorted(result, key=lambda entry: (self._chrom_order[str(entry.vcf_record.CHROM)], entry.vcf_record.POS, entry.vcf_record.REF))

    def _group_bpm_records(self, bpm_records):
        """
        Group BPM records into groups where all BPM records in a single
        group will be represented in the same VCF record

        The settings of the --unsquash-duplicates and --split-multiallelics
        options affect how BPM records are grouped.

        If there is only one BPM record for a given position, this record will
        be yielded as a singleton list.

        If neither --unsquash-duplicates nor --split-multiallelics is set, all
        the BPM records for a given position are yielded as a list.

        If either --unsquash-duplicates or --split-multiallelics are specified,
        and there multiple BPM records for a given position, these records will
        be grouped as described below, and these groups will be yielded one at
        a time.

        CASE 1: the position is biallelic (in other words, the same ALT allele
                appears in all the BPM records for this position).

                In this case, the --split-multiallelics flag is not applicable.

                Therefore, if the --unsquash-duplicates flag is True, the
                individual BPM records will be yielded as singleton lists, one
                at a time.  Otherewise, they will be yielded together as a
                single list.

        CASE 2: the position is multiallelic (in other words, more than one ALT
                allele appear among the BPM records for this position).

                If --split-multiallelics is not set, then all the records for
                this position will be yielded together as a single group.
                (Note, in particular, that in this case, the
                --unsquash-duplicates flag will have no effect, irrespective of
                the duplicates that may exist within the group.)

                If --split-multiallelics is set, and --unsquash-duplicates is
                also set, then the individual BPM records will be yielded as
                singleton lists, one at a time.

                If --split-multiallelics is set, but --unsquash-duplicates is
                not set, then the BPM records for a given position will be
                grouped according to their ALT allele (thus creating groups of
                duplicates), and these groups will be yielded one at a time.

        Args:
            bpm_records (BPMReader): Provides BPM records

        Yields:
            list(BPMRecord): Next group of BPM records
        """
        position2record = {}
        for record in bpm_records:
            position = (record.chromosome, record.pos, (record.get_indel_source_sequences(RefStrand.Plus)[1], record.is_deletion) if record.is_indel() else None)
            position2record.setdefault(position, []).append(record)

        for _, value in position2record.items():

            if (len(value) > 1 and
                (self._unsquash_duplicates or self._split_multiallelics)):

                refs = set()
                duplicates_grouped_by_alt = dict()

                for bpm_record in value:
                    ref, alt = bpm_record.plus_strand_alleles
                    refs.add(ref)
                    (duplicates_grouped_by_alt.setdefault(alt, [])
                     .append(bpm_record))

                assert len(refs) == 1

                nalts = len(duplicates_grouped_by_alt.value())

                if nalts == 1:

                    # POSITION IS BIALLELIC

                    # nalts == 1 means that there is only one alt allele for
                    # this position; in other words, this is a biallelic
                    # position; hence --split-multiallelics has no effect in
                    # this if-branch.

                    if self._unsquash_duplicates:
                        for bpm_record in value:
                            yield [bpm_record]
                    else:
                        yield value

                else:

                    # POSITION IS MULTIALLELIC

                    assert nalts > 1
                    # nalts > 1 means that there are multiple alt alleles for
                    # this position; in other words, this is a multiallelic
                    # position, and therefore --split-multiallelics has an
                    # effect in this if-branch.

                    if self._split_multiallelics:
                        if self._unsquash_duplicates:
                            for bpm_record in value:
                                yield [bpm_record]
                        else:
                            for group in duplicates_grouped_by_alt.values():
                                yield group
                    else:
                        assert self._unsquash_duplicates

                        # Since self._split_multiallelics is False here, we
                        # cannot split the current multiallelic position, and
                        # therefore the fact that self._unsquash_duplicates is
                        # True has no effect.

                        yield value

            else:
                yield value

    def _generate_locus_entry(self, bpm_record_group):
        """
        Generate a single VCF record fromm a group of BPM records

        Args:
            bpm_record_group (list(BPMRecord)): Group of BPM records for single site

        Returns:
            LocusEntry: LocusEntry for the site
        """
        return LocusEntry(bpm_record_group, self._vcf_record_factory.create_vcf_record(bpm_record_group))
