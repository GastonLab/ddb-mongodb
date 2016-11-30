#!/usr/bin/env python

import re
import sys
import utils
import cyvcf2
import getpass
import argparse
import vcf_parsing

from cyvcf2 import VCF
from datetime import datetime
from pymongo import MongoClient
from collections import defaultdict


def get_population_freqs(variant):
    freqs = {'esp_ea': variant.INFO.get('aaf_esp_ea') or -1,
             'esp_aa': variant.INFO.get('aaf_esp_aa') or -1,
             'esp_all': variant.INFO.get('aaf_esp_all') or -1,
             '1kg_amr': variant.INFO.get('aaf_1kg_amr') or -1,
             '1kg_eas': variant.INFO.get('aaf_1kg_eas') or -1,
             '1kg_sas': variant.INFO.get('aaf_1kg_sas') or -1,
             '1kg_afr': variant.INFO.get('aaf_1kg_afr') or -1,
             '1kg_eur': variant.INFO.get('aaf_1kg_eur') or -1,
             '1kg_all': variant.INFO.get('aaf_1kg_all') or -1,
             'exac_all': variant.INFO.get('aaf_exac_all') or -1,
             'adj_exac_all': variant.INFO.get('aaf_adj_exac_all') or -1,
             'adj_exac_afr': variant.INFO.get('aaf_adj_exac_afr') or -1,
             'adj_exac_amr': variant.INFO.get('aaf_adj_exac_amr') or -1,
             'adj_exac_eas': variant.INFO.get('aaf_adj_exac_eas') or -1,
             'adj_exac_fin': variant.INFO.get('aaf_adj_exac_fin') or -1,
             'adj_exac_nfe': variant.INFO.get('aaf_adj_exac_nfe') or -1,
             'adj_exac_oth': variant.INFO.get('aaf_adj_exac_oth') or -1,
             'adj_exac_sas': variant.INFO.get('aaf_adj_exac_sas') or -1}

    return freqs


def get_amplicon_data(variant):
    data = {'amplicon': variant.INFO.get('amplicon_target') or "None",
            'panel_amplicon': variant.INFO.get('panel_target') or "None",
            'intersect': variant.INFO.get('amplicon_intersect') or "None"}

    return data


if __name__ == "__main__":
    sample = ""
    library = ""

    client = MongoClient()
    db = client['ddb']

    parse_functions = {'mutect': vcf_parsing.parse_mutect_vcf_record,
                       'freebayes': vcf_parsing.parse_freebayes_vcf_record,
                       'vardict': vcf_parsing.parse_vardict_vcf_record,
                       'scalpel': vcf_parsing.parse_scalpel_vcf_record,
                       'platypus': vcf_parsing.parse_platypus_vcf_record,
                       'pindel': vcf_parsing.parse_pindel_vcf_record}

    caller_records = defaultdict(lambda: dict())

    sys.stdout.write("Parsing Caller VCF Files\n")
    vcf_parsing.parse_vcf("{}.mutect.normalized.vcf".format(sample), "mutect", caller_records)
    vcf_parsing.parse_vcf("{}.vardict.normalized.vcf".format(sample), "vardict", caller_records)
    vcf_parsing.parse_vcf("{}.freebayes.normalized.vcf".format(sample), "freebayes", caller_records)
    vcf_parsing.parse_vcf("{}.scalpel.normalized.vcf".format(sample), "scalpel", caller_records)
    vcf_parsing.parse_vcf("{}.platypus.normalized.vcf".format(sample), "platypus", caller_records)
    vcf_parsing.parse_vcf("{}.pindel.normalized.vcf".format(sample), "pindel", caller_records)

    annotated_vcf = "{}.vcfanno.snpEff.GRCh37.75.vcf".format(sample)

    sys.stdout.write("Parsing VCFAnno VCF\n")
    vcf = VCF(annotated_vcf)

    sys.stdout.write("Parsing VCFAnno VCF with CyVCF2\n")
    reader = cyvcf2.VCFReader(annotated_vcf)
    desc = reader["ANN"]["Description"]
    annotation_keys = [x.strip("\"'") for x in re.split("\s*\|\s*", desc.split(":", 1)[1].strip('" '))]

    # Filter out variants with minor allele frequencies above the threshold but
    # retain any that are above the threshold but in COSMIC or in ClinVar and not listed as benign.
    sys.stdout.write("Processing individual variants\n")

    for var in vcf:
        # Parsing VCF and creating data structures for Cassandra model
        callers = var.INFO.get('CALLERS').split(',')
        effects = utils.get_effects(var, annotation_keys)
        top_impact = utils.get_top_impact(effects)
        population_freqs = get_population_freqs(var)
        amplicon_data = get_amplicon_data(var)

        var_key_string = "GRCh37.75_{}_{}_{}_{}".format(var.CHROM, var.start, var.REF, var.ALT[0])
        library_var_key_string = "GRCh37.75_{}_{}_{}_{}_{}".format(library, var.CHROM, var.start, var.REF, var.ALT[0])

        caller_variant_data_dicts = defaultdict(dict)
        max_som_aaf = -1.00
        max_depth = -1
        min_depth = 100000000

        for caller in callers:
            caller_variant_data_dicts[caller] = parse_functions[caller](caller_records[caller][key])
            if float(caller_variant_data_dicts[caller]['AAF']) > max_som_aaf:
                max_som_aaf = float(caller_variant_data_dicts[caller]['AAF'])
            if int(caller_variant_data_dicts[caller]['DP']) < min_depth:
                min_depth = int(caller_variant_data_dicts[caller]['DP'])
            if int(caller_variant_data_dicts[caller]['DP']) > max_depth:
                max_depth = int(caller_variant_data_dicts[caller]['DP'])

        if min_depth == 100000000:
            min_depth = -1
