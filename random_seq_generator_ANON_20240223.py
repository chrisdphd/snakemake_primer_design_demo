#!/usr/local/bin/python3
# generates [a number] of randomers of [a fixed length] to run thru a filter (GC%, Tm, homopolymer check), and then output the desired number of filter-passing sequences, with properties

import random
import re
from Bio.SeqUtils import MeltingTemp as melttemp
from argparse import ArgumentParser



def generate_randomer(N, alphabet='ACGT'):
    return ''.join([random.choice(alphabet) for i in range(N)])

def dna_regex(dna):
    regexp = re.compile(r'AAAAA|TTTTTT|CCCC|GGGG')
    if regexp.search(dna):
      return 'matched'

def get_base_frequencies(dna):
    return {base: dna.count(base) / int(len(dna))
            for base in 'ATGC'}

def format_frequencies(frequencies):
    return ', '.join(['%s: %.3f' % (base, frequencies[base]) for base in frequencies])

def get_gcfreq(dna):
    return {base: dna.count(base) / int(len(dna))
            for base in 'GC'}

def format_gcpct(gc_freqs):
    cum = (sum([gc_freqs[base] for base in gc_freqs]))
    return '%.1f' % (cum * 100)

def calcTm(dnaSeq):
    tm = melttemp.Tm_NN(dnaSeq, Na = 50, Mg = 0.0, dNTPs = 0.0, saltcorr=5, dnac1=250, dnac2=0) ### for client 20180129
    return  '%.1f' % (tm)



def workflow(num_candidates, num_required, oligo_len, trim_5, prefix_tail5, suffix_tail3):

    candidates = 0
    num_passed = 0
    for cand_num in range(1, num_candidates+1, 1):
        if num_passed < num_required:
            candidates += 1
            dna = generate_randomer(oligo_len)
            subdna = dna[trim_5:] ### removes bases 0,1,2 or more from start (left side, or 5 prime end). End at 999 = i.e. effectively infinite, for this purpose

            homopol = dna_regex(dna)
            gc_freqs = get_gcfreq(dna)
            gc_pct = float(format_gcpct(gc_freqs))
            tm = float(calcTm(subdna))

            print ("Candidate:",dna)
            ### applying basic in-line tests/filter to the randomer:
            #if (36 <= gc_pct <= 55) and (42 <= tm <= 48) and homopol is not "matched": ### client barcodes?
            #if (36 <= gc_pct <= 55) and (66 <= tm <= 72) and homopol != "matched": ### client backbones
            if (20 <= gc_pct <= 80) and (50 <= tm <= 69) and homopol != "matched":

                frequencies = get_base_frequencies(dna)
                base_freqs = format_frequencies(frequencies)

                if prefix_tail5:
                    tail5 = prefix_tail5
                    dna = (tail5 + dna)
                if suffix_tail3:
                    tail3 = suffix_tail3
                    dna = (dna + tail3)
                final_length = len(dna)

                num_passed += 1

                print ("Passed:", str(num_passed) + '/' + str(cand_num), dna, 'length:', final_length, 'GC%:', gc_pct, 'Tm:', tm, 'randomer_base_freqs:', base_freqs)





def parseArgs():
    argparser = ArgumentParser(description='Creates RANDOMers and does some properties and filtering')
    # the arguments:
    argparser.add_argument('-nc', '--num_candidates', help='Number of CANDIDATE randomers to create', type=int, required=True)
    argparser.add_argument('-nr', '--num_required', help='(up to) Number of PASSING-FILTERS randomers to create', type=int, required=False)
    argparser.add_argument('-l', '--length', help='Length of randomers to create', type=int, required=True)
    argparser.add_argument('-t', '--trim', help='Num bases to trim from 5 prime end before getting Tm', type=int, default=0, required=False)
    argparser.add_argument('-p', '--prefix_tail5', help='Bases to add to 5 prime end of final randomer', type=str, default=False, required=False)
    argparser.add_argument('-s', '--suffix_tail3', help='Bases to add to 3 prime end of final randomer', type=str, default=False, required=False)

    return argparser.parse_args()


def main():
    args = parseArgs()

    num_candidates = args.num_candidates

    if args.num_required:
        num_required = args.num_required
    else:
        num_required = args.num_candidates    ### it will simply stop when candidates exhauseted

    oligo_len = args.length

    trim_5 = args.trim  ### i.e.  ignore 5' overhang for Tm measurement, because of client requirement. default is 0
    prefix_tail5 = args.prefix_tail5
    suffix_tail3 = args.suffix_tail3

    #if args.self:    # not implemented here
    #    self_dimer = "true"

    workflow(num_candidates, num_required, oligo_len, trim_5, prefix_tail5, suffix_tail3)


if __name__ == "__main__":
    main()
