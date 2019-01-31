"""
Copyright Government of Canada 2015

Written by: Eric Enns, National Microbiology Laboratory,
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

INSERTION = 1
DELETION = 2
REF_SKIP = 3
SOFT_CLIP = 4
HARD_CLIP = 5


def sam_alignment_to_padded_alignment(alignment, reference):
    ref_dna = reference.sub_seq(alignment.reference_start,
                                alignment.reference_end-1)
    query_dna = alignment.query_alignment_sequence

    pad_ref, pad_match, pad_query = '', '', ''

    for (operation, length) in alignment.cigartuples:
        if operation is SOFT_CLIP:
            """nothing needs to be done for soft clips"""
        elif operation is INSERTION:
            pad_ref += '-' * length
            pad_query += query_dna[0:length]
            query_dna = query_dna[length:]
            pad_match += ' ' * length
        elif operation is DELETION or operation is REF_SKIP:
            pad_ref += ref_dna[0:length]
            ref_dna = ref_dna[length:]
            pad_query += '-' * length
            pad_match += ' ' * length
        elif operation is HARD_CLIP:
            """nothing needs to be done for hard clips"""
        else:
            pad_ref += ref_dna[0:length]
            ref_dna = ref_dna[length:]
            pad_query += query_dna[0:length]
            query_dna = query_dna[length:]
            pad_match += '|' * length

    return (pad_ref, pad_match, pad_query)


def pairwise_alignment_to_differences(pad_ref, pad_query, ref_start):
    differences = dict()

    index = -1
    for i, c in enumerate(pad_ref):
        if c == '-':
            if ref_start + index not in differences.keys():
                differences[ref_start + index] = '.'
            differences[ref_start + index] += pad_query[i]
        else:
            index += 1
            if c is not pad_query[i]:
                differences[ref_start + index] = pad_query[i]

    return differences
