import argparse


def main(countfile, strand, sample):
    # Print the header
    print("X", sample, sample, sep="\t")
    columns = {"unstranded": 1, "forward": 2, "reverse": 3}
    # Select the column for the specified strand
    column = columns[strand]
    with open(countfile) as fin:
        # Skip the STAR summary headers
        for _ in range(4):
            next(fin)

        for i, line in enumerate(fin, 1):
            spline = line.strip("\n").split("\t")
            gene = spline[0]
            counts = spline[column]
            print(i, gene, counts, counts, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--count", required=True, help="STAR counts file")
    parser.add_argument(
        "--strand",
        choices=["unstranded", "forward", "reverse"],
        required=True,
        help="Strand of interest",
    )
    parser.add_argument("--sample", required=True, help="Sample name")

    args = parser.parse_args()

    main(args.count, args.strand, args.sample)
