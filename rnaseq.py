'''A class for normalizing RNAseq libraries'''

__author__ = 'Alex H Wagner'

import pandas
from biotk.range import Range
import sys


class RNASeq:

    def __init__(self, fpkm_df=None, trx_gtf=None, range_dict=None):
        self.gtf_model = None
        if trx_gtf:
            self.gtf_model = self.gtf_reader(trx_gtf)
        elif range_dict:
            self.gtf_model = range_dict
        if fpkm_df is not None:
            self.df = fpkm_df
            header = list(self.df)
            if header != ['id', 'FPKM']:
                raise ValueError('Expected dataframe header in form ["id", "FPKM"]. Got {0}'
                                 .format(header))
            #calculate pseudo-counts, TPMs

    def gtf_reader(gtf):
        d = dict()
        for i, line in enumerate(gtf):
            sys.stdout.write("\rEvaluating line {0} of gtf".format(i))
            sys.stdout.flush()
            a = line.split('\t')
            if a[2] != 'exon':
                continue
            b = {x: y for x, y in [z.split() for z in a[8].split(';')[:-1]]}
            gene = b['gene_id']
            trx = b['transcript_id']
            for id in [gene, trx]:
                try:
                    if d[id]['chr'] != a[0]:
                        raise 'Trans-scaffold id {0}'.format(id)
                    d[id]['start'].append(int(a[3]))
                    d[id]['end'].append(int(a[4]))
                except KeyError:
                    d[id] = {'chr': a[0],
                             'start': [int(a[3])],
                             'end': [int(a[4])]}
        sys.stdout.write("\n")
        for id in d:
            l = list(zip(d[id]['start'], d[id]['end']))
            d[id] = Range(l, d[id]['chr'])
        return d

if __name__ == '__main__':
    # import pickle
    # with open('data/transcripts.gtf') as f:
    #     d = RNASeq.gtf_reader(f)
    # with open('data/gtf_dict.pickle', 'wb') as f:
    #     pickle.dump(d, f)
    pass
