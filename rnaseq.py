"""A class for normalizing RNAseq libraries"""

__author__ = 'Alex H Wagner'

import pandas
from biotk.range import Range
import sys


class RNASeq:

    def __init__(self, fpkm_df=None, fpkm_tracking_file=None, trx_gtf_file=None, range_dict=None, read_length=None,
                 aligned_fragment_count=None, counts_df=None):
        self.gtf_model = None
        if range_dict:
            self.gtf_model = range_dict
        elif trx_gtf_file:
            self.gtf_model = self.gtf_reader(trx_gtf_file)
        else:
            raise AttributeError('No gtf dict or gtf file provided.')
        # if counts_df
        # elif count_tracking_file
        # convert following statement to elif
        if fpkm_df is not None:
            self.df = fpkm_df
            header = list(self.df)
            if header != ['id', 'FPKM']:
                raise ValueError('Expected dataframe header in form ["id", "FPKM"]. Got {0}'
                                 .format(', '.join(header)))
        elif fpkm_tracking_file:
            df = pandas.read_csv(fpkm_tracking_file, sep="\t")
            df = df[['tracking_id', 'FPKM']]
            df.columns = ['id', 'FPKM']
            self.df = df
        else:
            raise AttributeError('No fpkm dataframe or tracking file provided.')
        if read_length:
            self.rl = read_length
        else:
            raise AttributeError('No average read length provided.')
        if aligned_fragment_count is None:
            aligned_fragment_count = 1000000

        # calculate pseudo-counts
        # formulae come from http://lynchlab.uchicago.edu/publications/Wagner,%20Kin,%20and%20Lynch%20(2012).pdf
        self.df['length'] = [self.gtf_model[x].length for x in self.df['id']]
        self.df['pcounts'] = self.df['FPKM'] * self.df['length'] * aligned_fragment_count / 1000000000

    @staticmethod
    def gtf_reader(gtf):
        d = dict()
        for i, line in enumerate(gtf):
            if i % 1000 == 0:
                sys.stdout.write("\rEvaluating line {0} of gtf".format(i))
                sys.stdout.flush()
            a = line.split('\t')
            if a[2] != 'exon':
                continue
            b = {x: y.strip('"') for x, y in [z.split() for z in a[8].split(';')[:-1]]}
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
        sys.stdout.write("\rEvaluating line {0} of gtf\n".format(i))
        sys.stdout.flush()
        for id in d:
            l = list(zip(d[id]['start'], d[id]['end']))
            d[id] = Range(l, d[id]['chr'])
        return d

    def normalized_expr(self, index):
        """Return a dataframe of TPM normalized to only the indexed genes/transcripts"""
        pass


if __name__ == '__main__':
    # import pickle
    # with open('data/transcripts.gtf') as f:
    #     d = RNASeq.gtf_reader(f)
    # with open('data/gtf_dict.pickle', 'wb') as f:
    #     pickle.dump(d, f)
    pass
