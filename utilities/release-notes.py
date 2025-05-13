#!/usr/bin/env python3

import argparse
import re


def main(changelog_file, version):
    changes = changelog_by_release(changelog_file)

    for change in changes.get(version, list()):
        print(change)


def is_header(line):
    pattern = r"^\*+$"
    return re.match(pattern, line)


def changelog_by_release(changelog_file):
    results = dict()

    # Previous release, for when we reach the end of the changes
    prev_release = None
    changes = list()

    with open(changelog_file) as fin:
        for line in fin:
            if is_header(line):
                # Get the next release
                release = next(fin).strip("\n")
                # Skip next header line
                next(fin)

                if prev_release is not None:
                    results[prev_release] = changes

                prev_release = release
                changes = list()
            else:
                if change := line.strip("\n"):
                    changes.append(change)

        results[prev_release] = changes

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("changelog")
    parser.add_argument("version")

    args = parser.parse_args()
    main(args.changelog, args.version)
