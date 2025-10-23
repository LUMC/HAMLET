#!/usr/bin/env python3

import sys


def get_changelog_version() -> str:
    changelog = "CHANGELOG.rst"
    with open(changelog) as fin:
        for line in fin:
            # The next line holds the version number
            if line.startswith("******"):
                changelog_version = next(fin).strip("\n")
                return changelog_version
    return ""


def get_hamlet_version() -> str:
    hamlet = "common.smk"
    with open(hamlet) as fin:
        for line in fin:
            if line.startswith("PIPELINE_VERSION"):
                hamlet_version = line.split('"')[1]
                return hamlet_version
    return ""


def get_doc_version() -> str:
    config = "docs/source/conf.py"
    with open(config) as fin:
        for line in fin:
            if line.startswith("release = "):
                doc_version = line.split('"')[1]
                return doc_version
    return ""


if __name__ == "__main__":

    version: str = sys.argv[1]

    changelog_version = get_changelog_version()
    hamlet_version = get_hamlet_version()
    docs_version = get_doc_version()
    print(f"Expected version:  {version}")
    print(f"Changelog version: {changelog_version}")
    print(f"HAMLET version:    {hamlet_version}")
    print(f"Docs version:      {docs_version}")

    if (
        changelog_version == version
        and hamlet_version == version
        and docs_version == version
    ):
        exit(0)
    else:
        exit(1)
