"""
HBB Sanger Sequence Analyzer — minimal demo
=============================================
A simplified, working example of the kind of app at hbb-sanger.streamlit.app.

What it does:
1. Accepts a FASTA (.fa/.fasta) or raw pasted sequence
   (real tools also accept .ab1 chromatogram files via Bio.SeqIO's "abi" parser)
2. Aligns it to a reference HBB coding sequence
3. Calls point differences (simple substitution-level "mutations")
4. Flags any that match a small known-variant lookup table (e.g. HbS)
5. Displays results in a table + simple alignment view

Run locally with:  streamlit run app.py
Deploy by pushing this + requirements.txt to GitHub, then
connecting the repo at https://share.streamlit.io
"""

import streamlit as st
from Bio import Align

# ---------------------------------------------------------------------------
# 1. Reference data
# ---------------------------------------------------------------------------
# Real apps load a curated reference (e.g. NM_000518.5 HBB CDS) from a file.
# Here's a short illustrative stretch of the HBB coding sequence
# (exon 1, start of the beta-globin chain) for demonstration purposes.
HBB_REFERENCE = (
    "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAA"
    "GTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTT"
    "GGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTC"
)

# A tiny illustrative lookup of known pathogenic single-base changes.
# Real tools use a curated clinical variant database (e.g. ClinVar, HbVar).
KNOWN_VARIANTS = {
    20: {
        "wild_type": "A",
        "mutant": "T",
        "name": "HBB c.20A>T (p.Glu7Val)",
        "note": "Causes HbS — sickle cell disease (classic example)",
    },
}

st.set_page_config(page_title="HBB Sanger Analyzer (Demo)", layout="wide")
st.title("🧬 HBB Sanger Sequence Analyzer — Demo")
st.caption(
    "A minimal illustration of how a tool like hbb-sanger.streamlit.app "
    "is structured. Not for clinical use."
)

# ---------------------------------------------------------------------------
# 2. Input
# ---------------------------------------------------------------------------
st.subheader("1. Provide a sequence")

input_mode = st.radio("Input method", ["Paste sequence", "Upload FASTA file"], horizontal=True)

query_seq = None

if input_mode == "Paste sequence":
    pasted = st.text_area(
        "Paste a DNA sequence (A/C/G/T only)",
        placeholder="ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAA...",
        height=120,
    )
    if pasted.strip():
        query_seq = "".join(pasted.split()).upper()

else:
    uploaded = st.file_uploader("Upload a .fasta / .fa file", type=["fasta", "fa", "txt"])
    if uploaded:
        text = uploaded.read().decode("utf-8")
        lines = [l.strip() for l in text.splitlines() if l.strip() and not l.startswith(">")]
        query_seq = "".join(lines).upper()

# Demo button so reviewers can see it work without their own file
if st.button("Use example sequence (with HbS mutation)"):
    mutated = list(HBB_REFERENCE)
    mutated[20 - 1] = "T"  # introduce the example A>T at position 20
    query_seq = "".join(mutated)

# ---------------------------------------------------------------------------
# 3. Analysis
# ---------------------------------------------------------------------------
if query_seq:
    st.subheader("2. Alignment to HBB reference")

    if not set(query_seq) <= set("ACGTN"):
        st.error("Sequence contains characters other than A/C/G/T/N — please check input.")
        st.stop()

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(HBB_REFERENCE, query_seq)[0]
    ref_aln, qry_aln = str(alignment[0]), str(alignment[1])

    with st.expander("Show raw pairwise alignment"):
        st.code(str(alignment))

    # ------------------------------------------------------------------
    # 4. Variant calling — walk the alignment and report mismatches
    # ------------------------------------------------------------------
    st.subheader("3. Detected differences")

    diffs = []
    ref_pos = 0
    for ref_base, qry_base in zip(ref_aln, qry_aln):
        if ref_base != "-":
            ref_pos += 1
        if ref_base != qry_base and ref_base != "-" and qry_base != "-":
            diffs.append((ref_pos, ref_base, qry_base))

    if not diffs:
        st.success("No differences detected vs. reference.")
    else:
        rows = []
        for pos, ref_b, qry_b in diffs:
            known = KNOWN_VARIANTS.get(pos)
            rows.append({
                "Position": pos,
                "Reference": ref_b,
                "Observed": qry_b,
                "Known variant?": known["name"] if known else "—",
                "Clinical note": known["note"] if known else "Not in demo lookup table",
            })
        st.dataframe(rows, use_container_width=True)

        flagged = [r for r in rows if r["Known variant?"] != "—"]
        if flagged:
            st.warning(f"⚠️ {len(flagged)} known variant(s) detected — see table above.")

else:
    st.info("Paste a sequence, upload a FASTA file, or click the example button to begin.")

st.divider()
st.caption(
    "This is a teaching example: tiny reference fragment, toy variant table, "
    "simple pairwise alignment. A production tool would use a full reference "
    "genome/transcript, a real chromatogram (.ab1) parser, quality-score "
    "trimming, and a maintained clinical variant database."
)
