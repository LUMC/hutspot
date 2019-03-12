#   hutspot - a DNAseq variant calling pipeline
#   Copyright (C) 2017-2019, Sander Bollen, Leiden University Medical Center
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Affero General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Affero General Public License for more details.
#
#   You should have received a copy of the GNU Affero General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
collect_stats.py
~~~~~~~~~~~~~~~~

:copyright: (c) 2017-2019 Sander Bollen
:copyright: (c) 2017-2019 Leiden University Medical Center
:license: AGPL-3.0
"""

import click
import json


def parse_json(path):
    with open(path) as handle:
        return json.load(handle)


@click.command()
@click.option("--vcfstats",
              type=click.Path(exists=True, dir_okay=False, readable=True),
              required=True,
              help="Path to vcfstats json")
@click.argument("collectstats",
              type=click.Path(exists=True, dir_okay=False, readable=True),
              nargs=-1,
              required=True)
def main(vcfstats, collectstats):
    v = parse_json(vcfstats)
    cs = [parse_json(x) for x in collectstats]
    d = {
        "sample_stats": cs,
        "multisample_vcfstats": v
    }
    print(json.dumps(d))


if __name__ == "__main__":
    main()
