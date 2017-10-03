"""
Copyright Government of Canada 2017

Written by: Camy Tran, National Microbiology Laboratory,
            Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

from collections import defaultdict


def parse_genes_from_codon_variants_csv(csv):
    genes = defaultdict(lambda: defaultdict(lambda:
                        defaultdict(lambda: defaultdict(int))))

    with open(csv, "r") as f:
        for line in f:
            if line[0] != "#":
                (
                    gene, gene_start_end, nt_start,
                    nt_end, ref_codon, mutant_codon,
                    ref_aa, mutant_aa, coverage,
                    mutant_freq, mutant_type,
                    ns_count, s_count
                ) = line.rstrip().split(",")
                gene_start, gene_end = gene_start_end.split('-')

                gene_start = int(gene_start)
                gene_end = int(gene_end)
                nt_start = int(nt_start)
                nt_end = int(nt_end)
                mutant_freq = float(mutant_freq)
                ns_count = float(ns_count)
                s_count = float(s_count)

                genes[gene]['start'] = gene_start
                genes[gene]['end'] = gene_end

                if ns_count > 0:
                    genes[gene][nt_start-gene_start][
                        'NS'][ns_count] += mutant_freq/100.0
                if s_count > 0:
                    genes[gene][nt_start-gene_start][
                        'S'][s_count] += mutant_freq/100.0

    f.close()

    return genes
