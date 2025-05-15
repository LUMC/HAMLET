#!/usr/bin/env python3

import functools
import sys
from release_notes import changelog_by_release
import pytest
from collections import OrderedDict
import copy

@pytest.fixture()
def changelog():
    return OrderedDict(
        {"v1.0-dev": ["changes"],
         "v0.9": ["older", "changes"]
         }
    )

def test_set_new_version(changelog):
    new = set_version(changelog, "v2.0")
    assert next(iter(new)) == "v2.0"

def test_unset_dev(changelog):
    new = set_version(changelog, "v1.0")
    assert next(iter(new)) == "v1.0"
    assert len(new) == 2

def test_set_existing_version(changelog):
    new = set_version(changelog, "v1.0-dev")
    latest_version, changes = next(iter(new.items()))

    assert latest_version == "v1.0-dev"
    assert changes == ["changes"]

def set_version(changelog: OrderedDict, version: str) -> OrderedDict:
    """Set version of the ordered dict

    Either by adding a new entry, or removing -dev from the latest entry
    """
    new = copy.deepcopy(changelog)

    latest_version = next(iter(changelog))

    # If the version is already set, we do nothing
    if latest_version == version:
        pass
    elif latest_version == f"{version}-dev":
        dev_version, changes = new.popitem(last=False)
        new[version] = changes
    else:
        new[version] = list()
    new.move_to_end(version, last=False)
    return new

def set_version_changelog(new_version) -> None:
    changelog = "CHANGELOG.rst"
    releases = changelog_by_release(changelog)
    new = set_version(releases, new_version)

    # get the start of the file
    file_start = list()
    with open(changelog) as fin:
        for line in fin:
            if line.startswith("***"):
                break
            else:
                file_start.append(line.strip("\n"))

    with open(changelog, 'wt') as fout:
        write = functools.partial(print, end='\n', file=fout)
        # Write the start of the file
        for line in file_start:
            write(line)

        # Write the updated changelog
        for version, changes in new.items():
            write("*"*len(version))
            write(version)
            write("*"*len(version))
            write()
            for change in changes:
                write(change)
            write()


def set_version_hamlet(version):
    hamlet = "common.smk"
    with open(hamlet) as fin:
        lines = fin.readlines()
    with open(hamlet, 'wt') as fout:
        write = functools.partial(print, end='\n', file=fout)
        for line in lines:
            if line.startswith("PIPELINE_VERSION"):
                write(f'PIPELINE_VERSION = "{version}"')
            else:
                write(line, end='')


def set_version_docs(version) -> str:
    config = "docs/source/conf.py"
    with open(config) as fin:
        lines = fin.readlines()
    with open(config, 'wt') as fout:
        write = functools.partial(print, end='\n', file=fout)
        for line in lines:
            if line.startswith("release = "):
                write(f'release = "{version}"')
            else:
                write(line, end='')


if __name__ == "__main__":

    version: str = sys.argv[1]

    set_version_changelog(version)
    set_version_hamlet(version)
    set_version_docs(version)
