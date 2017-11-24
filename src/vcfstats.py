import click
import cyvcf2
import numpy


class Sample(object):
    def __init__(self, name):
        self.name = name
        self.transversions = 0
        self.transitions = 0
        self.hom_ref = 0
        self.het = 0
        self.hom_alt = 0

    @property
    def ti_tv(self):
        if self.transversions > 0:
            return self.transitions/self.transversions
        return numpy.nan



@click.command()
@click.option("-i",
              "--input",
              type=click.Path(exists=True, dir_okay=False, readable=True),
              required=True,
              help="Input VCF file")
@click.option("-o",
              "--output",
              type=click.Path(dir_okay=False, writable=True),
              required=True,
              help="Output json")
def main(input, output):
    pass


if __name__ == "__main__":
    main()