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