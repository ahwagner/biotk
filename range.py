__author__ = 'Alex H Wagner'

import operator as o

class Range:
    """A class encompassing a genomic range on a scaffold. Starts and stops are 1-based."""

    def __init__(self, sub_ranges=(), scaffold=None):
        sub_ranges = list(sub_ranges)
        if len(sub_ranges) == 2 and isinstance(sub_ranges[0], int) and isinstance(sub_ranges[1], int):
            sub_ranges = [sub_ranges]
        for i, r in enumerate(sub_ranges):
            if len(r) != 2 or not isinstance(r[0], int) or not isinstance(r[1], int):
                raise ValueError('Ranges should be a list of (start <int>, end <int>) iterables.\n'
                                 'Instead got {0}.'.format(sub_ranges))
            sub_ranges[i] = list(r)
        self.scaffold = scaffold
        self.sub_ranges = sub_ranges
        self.length = 0
        self.merge()

    def sort(self):
        self.sub_ranges = sorted(self.sub_ranges, key=o.itemgetter(0, 1))

    def __add__(self, other):
        if isinstance(other, Range):
            r = Range(self.sub_ranges + other.sub_ranges) # This isn't as efficient as doing a sorted merge. Fix later.
            if self.scaffold is None:
                r.scaffold = other.scaffold
            elif other.scaffold is None:
                r.scaffold = self.scaffold
            elif self.scaffold == other.scaffold:
                r.scaffold = self.scaffold
            else:
                raise ValueError('Scaffolds do not match.')
        elif len(other) == 2 and isinstance(other[0], int) and isinstance(other[1], int):
            o = Range(other, self.scaffold)
            r = self + o
        else:
            raise TypeError('No method for adding Range and elements of {0}.'.format(type(other)))
        return r

    def merge(self):
        if not self.sub_ranges:
            return
        self.sort()
        new = []
        last = list(self.sub_ranges.pop(0))
        while self.sub_ranges:
            nxt = self.sub_ranges.pop(0)
            if last[1] >= nxt[0]:
                last[1] = max(last[1], nxt[1])
            else:
                new.append(tuple(last))
                last = nxt
        new.append(tuple(last))
        self.sub_ranges = new
        self.length = 0
        for t in new:
            self.length += (t[1] - t[0] + 1)

    def __iter__(self):
        for r in self.sub_ranges:
            yield((self.scaffold, r[0], r[1]))


if __name__ == '__main__':
    r = Range((1, 100), 'chr1')
    r2 = r + (80, 160)
    for r in r2:
        print(r)
