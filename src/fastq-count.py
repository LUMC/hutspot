"""
fastq-count.py

Pythonic equivalent to https://github.com/sndrtj/fastq-count
"""
import click
import gzip
import json

from collections import namedtuple


class SimpleFastqRecord(namedtuple("SimpleFastqRecord",
                                   ["read_id", "sequence", "qualities"])):

    def __str__(self):
        return "{0}{1}+\n{2}".format(self.read_id.decode(),
                                     self.sequence.decode(),
                                     self.qualities.decode())

    def __bytes__(self):
        return self.read_id + self.sequence + b"+\n" + self.qualities

    @property
    def id(self):
        try:
            return self._id
        except AttributeError:
            self._id = self.read_id.strip().split()[0][1:]
        return self._id

    @property
    def seq(self):
        return self.sequence.strip().decode()


class SimpleFastqParser(object):
    """
    An iterator returning SimpleFastqRecord objects

    :arg handle: Any iterator returning lines of bytestrings,
    preferable in open file handle in rb mode.
    :return SimpleFastqRecord objects
    """
    def __init__(self, handle):
        self.__handle = handle
        self.__bucket = [None, None, None]

    def __iter__(self):
        return self

    def __next__(self):
        i = 0
        while i < 3:
            line = next(self.__handle)
            if line == b"+\n":
                continue
            self.__bucket[i] = line
            i += 1
        read = SimpleFastqRecord(self.__bucket[0], self.__bucket[1],
                                 self.__bucket[2])
        self.__bucket = [None, None, None]
        return read

    def next(self):  # python 2 compatibility
        return self.__next__()

    def close(self):
        self.__handle.close()


def count(handle_r1, handle_r2):
    """
    Count reads and bases for handles in rb mode
    :param handle:
    :return:
    """
    p = SimpleFastqParser(handle_r1)
    reads = 0
    bases = 0
    for r in p:
        reads += 1
        bases += len(r.seq)
    p2 = SimpleFastqParser(handle_r2)
    for r2 in p2:
        reads += 1
        bases += len(r2.seq)
    return {"reads": reads, "bases": bases}


@click.command()
@click.argument("r1",
                type=click.Path(dir_okay=False, readable=True, exists=True),
                required=True)
@click.argument("r2",
                type=click.Path(dir_okay=False, readable=True, exists=True),
                required=True)
def main(r1, r2):
    if r1.endswith(".gz"):
        r1_handle = gzip.open(r1, mode="rb")
    else:
        r1_handle = open(r1, mode="rb")

    if r2.endswith(".gz"):
        r2_handle = gzip.open(r2, mode="rb")
    else:
        r2_handle = open(r2, mode="rb")

    d = count(r1_handle, r2_handle)
    print(json.dumps(d))


if __name__ == "__main__":
    main()