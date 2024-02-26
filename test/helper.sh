#!/usr/bin/env bash

set -euxo pipefail

ref=test/data/reference/hamlet-ref.fa
star=test/data/reference/hamlet-star/Genome

function setup {
  unxz -k ${ref}.xz
  unxz -k ${star}.xz
}

function cleanup {
  rm -f ${ref}
  rm -f ${star}
}

trap cleanup EXIT

setup
