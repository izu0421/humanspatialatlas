"""The query payload: the sample table in a form the browser can filter, with coverage.

The dashboard ships aggregates. A question like "how many Xenium breast cancer
samples are there" cannot be answered from an aggregate, so the per-sample table
travels too, dictionary-encoded: 13,652 human samples over 12 fields come to
about 550 kB of JSON and 80 kB over the wire.

COVERAGE TRAVELS WITH IT, per field, because this catalogue's fields are filled
in very unevenly and a count computed over a sparse field looks exactly like a
count computed over a complete one:

    assay            100%        sex               31%
    tissue            97%        age group         21%
    disease state     90%

So "how many female samples" has 4,257 records that could answer it and 9,395
that are silent, and an answer that does not say so is wrong by omission. This
project has shipped that bug four times in a week in other forms: a blank field
read as a measured negative.
"""
from __future__ import annotations

import pandas as pd

FIELDS = [
    ("assay", "assay"), ("tissue", "tissue"), ("disease_state", "disease state"),
    ("disease_level1", "disease"), ("disease_level2", "disease detail"),
    ("sex", "sex"), ("age_group", "age group"), ("preservation", "preservation"),
    ("source", "source"), ("accession", "study"), ("donor_id", "donor"),
    ("metadata_source", "metadata source"),
]
BLANK = {"", "unknown", "none", "nan", "na", "n/a", "not applicable", "not collected"}


def encode(series: pd.Series) -> dict:
    """Dictionary-encode one column, with every blank spelling mapped to index 0.

    Index 0 is always the empty string, so the browser can test "is this field
    known for this sample" with one integer comparison rather than a string
    match against a list of spellings that will grow.
    """
    s = series.fillna("").astype(str).str.strip()
    s = s.where(~s.str.lower().isin(BLANK), "")
    vals = [""] + sorted(v for v in set(s) if v)
    idx = {v: i for i, v in enumerate(vals)}
    return {"v": vals, "c": [idx[x] for x in s]}


def build(sm: pd.DataFrame) -> dict:
    """The query payload for the human samples in `sm`."""
    hum = sm[sm.species.isin(["homo_sapiens", "human_in_mouse_xenograft"])].copy()
    cols, cover = {}, {}
    for key, label in FIELDS:
        if key not in hum.columns:
            continue
        enc = encode(hum[key])
        cols[key] = {**enc, "label": label}
        cover[key] = int(sum(1 for c in enc["c"] if c))
    return {"n": int(len(hum)), "cols": cols, "coverage": cover,
            "n_nonhuman": int(len(sm) - len(hum))}
