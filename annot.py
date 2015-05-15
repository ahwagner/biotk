__author__ = 'Alex H Wagner'

import sqlite3
import csv
import pathlib
import warnings

class Gene:

    @staticmethod
    def import_database():
        # Hard-coded for now. Maybe address dynamic tables and updates through API in the future.
        c = Gene._conn.cursor()
        try:
            c.execute("DROP TABLE Gene;")
        except sqlite3.OperationalError:
            pass
        c.executescript("CREATE TABLE Gene ('ensembl_id','entrez_id','hgnc_symbol');"
                        "CREATE INDEX GeneEnsembl ON Gene(ensembl_id);"
                        "CREATE INDEX GeneHgnc ON Gene(hgnc_symbol);"
                        "CREATE INDEX GeneEntrez ON Gene(entrez_id);")
        with open(str(Gene._path / 'data'/ 'biomart_gene.tsv')) as f:
            dr = csv.DictReader(f, delimiter="\t")
            to_db = [(i['Ensembl Gene ID'], i['EntrezGene ID'], i['HGNC symbol']) for i in dr]
        c.executemany("INSERT INTO Gene ('ensembl_id','entrez_id','hgnc_symbol') VALUES (?, ?, ?);", to_db)
        c.close()
        Gene._conn.commit()

    annotation_types = ['ensembl_id', 'entrez_id', 'hgnc_symbol']
    _type_error_string = "Unknown id type. Valid options are:\n{0}".format("\n".join(annotation_types))
    _path = pathlib.Path(__file__).parent
    _conn = sqlite3.connect(str(_path / 'data' / 'annot.db'))

    def __init__(self, identifier, id_type):
        if id_type.lower() in Gene.annotation_types:
            self.annotations = {id_type: str(identifier)}
        else:
            raise TypeError(Gene._type_error_string)
        self._c = Gene._conn.cursor()
        self._lookup()

    def is_defined(self, item):
        return item.lower() in self.annotations

    def __getattr__(self, item):
        item = item.lower()
        if item in Gene.annotation_types:
            if not self.is_defined(item):
                return ''
            return self.annotations[item]
        else:
            raise TypeError(Gene._type_error_string)

    def _lookup(self):
        id_type = list(self.annotations)[0]
        id_value = self.annotations[id_type]
        sql = "SELECT * FROM Gene WHERE {0} = '{1}'".format(id_type, id_value)
        self._c.execute(sql)
        r = self._c.fetchall()
        if len(r) == 0:
            pass    # Do nothing if not found in table.
        elif len(r) > 1:
            for row in r:
                for i, v in enumerate(row):
                    try:
                        self.annotations[Gene.annotation_types[i]].append(v)
                    except (AttributeError, KeyError):
                        self.annotations[Gene.annotation_types[i]] = [v]
        else:
            self.annotations = dict(zip(Gene.annotation_types, r[0]))


if __name__ == "__main__":
    # g = Gene(1, 'entrez_id')
    # print(g.ensembl_id)
    Gene.import_database()