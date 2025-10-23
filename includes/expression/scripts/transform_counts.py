import argparse


def main(countfile: str, strand: str, sample: str) -> None:
    columns = {"unstranded": 1, "forward": 2, "reverse": 3}
    # Select the column for the specified strand
    column = columns[strand]
    with open(countfile) as fin:
        # Skip the STAR summary headers
        for _ in range(4):
            next(fin)

        for line in fin:
            spline = line.strip("\n").split("\t")
            gene = spline[0]
            counts = spline[column]
            print(gene, counts, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--counts", required=True, help="STAR counts file")
    parser.add_argument(
        "--strand",
        choices=["unstranded", "forward", "reverse"],
        required=True,
        help="Strand of interest",
    )
    parser.add_argument("--sample", required=True, help="Sample name")

    args = parser.parse_args()

    main(args.counts, args.strand, args.sample)
