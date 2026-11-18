#!/usr/bin/env python3
"""
hr_clone_validate.py
=====================
In-silico validation of Vazyme C116 (ClonExpress II / Ultra) homologous-
recombination cloning.

Given N PCR fragments (N >= 2), each amplified by one primer pair from its own
SnapGene .dna template, the tool will:

  1. Perform in-silico PCR for each fragment -> N double-stranded amplicons.
     The 5' portion of a primer that does NOT match its template is reported
     as the added homology (HR) arm / overhang.
  2. Search for a circular assembly order (and per-fragment orientation) in
     which every consecutive junction (incl. last -> first) shares terminal
     homology (suffix of one amplicon == prefix of the next).
  3. Decide whether recombination can circularize all fragments
     (ClonExpress needs ~15-20 bp of terminal homology at each junction).
  4. Report each HR-arm sequence + length, the assembly order/orientations,
     the assembled circular product, and a recommended annealing temperature
     (nearest-neighbour Tm of each whole primer).

N fragments == N templates (templates may be the same file; in the built-in
demo all three fragments come from one .dna).

Dependencies: snapgene_reader, biopython.

It also does SELF-HR MUTAGENESIS (around-the-horn): one template, one divergent
primer pair whose single linear product self-circularises because its two ends
share an HR arm. Reports the arm, the circular product size, and classifies the
edit (point substitution / insertion / deletion) vs the template.

Usage:
  # built-in 3-fragment demo (tiles one plasmid into 3 amplicons):
  python hr_clone_validate.py

  # explicit fragments (one template each; LABEL@file.dna, or --template default):
  python hr_clone_validate.py --template pcDNA3-BTK-nLuc.dna \
      --fragment VECTOR  VecF=...  vecR=... \
      --fragment INSERT1 aF=...    aR=... \
      --fragment INSERT2@other.dna bF=... bR=...

  # self-HR mutagenesis (one template, one primer pair):
  python hr_clone_validate.py --template pcDNA3-BTK-nLuc.dna \
      --mutagenesis MutF=... MutR=...
  python hr_clone_validate.py --demo-mut          # built-in deletion demo
"""

import argparse
import os
import sys

from snapgene_reader import snapgene_file_to_dict
from Bio.SeqUtils import MeltingTemp as mt
try:                                   # biopython >= 1.80
    from Bio.SeqUtils import gc_fraction
    def _gc(s):
        return gc_fraction(s)
except ImportError:                    # biopython <= 1.79
    from Bio.SeqUtils import GC
    def _gc(s):
        return GC(s) / 100.0

# --------------------------------------------------------------------------- #
# basic sequence helpers
# --------------------------------------------------------------------------- #
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def load_template(path: str):
    """Return (sequence_upper, is_circular) from a SnapGene .dna file."""
    d = snapgene_file_to_dict(path)
    seq = d["seq"].upper()
    topo = (d.get("dna") or {}).get("topology", "linear")
    return seq, (topo == "circular")


# --------------------------------------------------------------------------- #
# in-silico PCR
# --------------------------------------------------------------------------- #
class AnnealSite:
    """A primer's 3'-anchored binding site, expressed on the template top strand."""

    def __init__(self, name, strand, top_start, top_end, anneal_seq, tail):
        self.name = name
        self.strand = strand            # '+' forward, '-' reverse
        self.top_start = top_start      # inclusive, on template top strand
        self.top_end = top_end          # exclusive
        self.anneal_seq = anneal_seq    # template-binding part of the primer (5'->3')
        self.tail = tail                # 5' overhang of the primer that does NOT match


def _find_3prime_anchored(top_primer, search, anchored, seed_len, tmpl_len):
    """
    Locate matches of `top_primer` against `search` (a possibly-doubled template),
    anchoring the primer's 3' end.

      anchored == 'right' : matching part is a SUFFIX of top_primer (forward case)
      anchored == 'left'  : matching part is a PREFIX of top_primer (reverse case)

    Returns list of (top_start, top_end, anneal_len) with starts in [0, tmpl_len).
    """
    hits = []
    eff_seed = min(seed_len, len(top_primer))

    if anchored == "right":
        seed = top_primer[-eff_seed:]
        j = search.find(seed)
        while j != -1:
            if j < tmpl_len:
                E = j + eff_seed                      # 3' end of primer on top strand
                k = 0
                while (j - 1 - k >= 0
                       and len(top_primer) - eff_seed - 1 - k >= 0
                       and top_primer[len(top_primer) - eff_seed - 1 - k] == search[j - 1 - k]):
                    k += 1
                anneal_len = eff_seed + k
                hits.append((E - anneal_len, E, anneal_len))
            j = search.find(seed, j + 1)
    else:  # anchored == 'left'
        seed = top_primer[:eff_seed]
        j = search.find(seed)
        while j != -1:
            if j < tmpl_len:
                S = j                                 # 3' end of primer on top strand (leftmost)
                k = 0
                while (S + eff_seed + k < len(search)
                       and eff_seed + k < len(top_primer)
                       and top_primer[eff_seed + k] == search[S + eff_seed + k]):
                    k += 1
                anneal_len = eff_seed + k
                hits.append((S, S + anneal_len, anneal_len))
            j = search.find(seed, j + 1)
    return hits


def map_primer(name, primer, template, is_circular, seed_len=15):
    """
    Find the best binding site of `primer` on `template` (either strand).
    Returns a list of AnnealSite (usually 1). Longest anneal wins.
    """
    primer = primer.upper()
    L = len(template)
    search = template + template if is_circular else template
    sites = []

    # forward: primer reads on top strand, suffix matches, 3' end on the right
    for s, e, alen in _find_3prime_anchored(primer, search, "right", seed_len, L):
        tail = primer[: len(primer) - alen]
        sites.append(AnnealSite(name, "+", s % L if is_circular else s,
                                e, primer[len(primer) - alen:], tail))

    # reverse: revcomp(primer) reads on top strand, prefix matches, 3' end on the left
    q = revcomp(primer)
    for s, e, alen in _find_3prime_anchored(q, search, "left", seed_len, L):
        tail_rc = q[alen:]                  # revcomp of the primer's 5' overhang
        tail = revcomp(tail_rc)
        sites.append(AnnealSite(name, "-", s, e, revcomp(q[:alen]), tail))

    sites.sort(key=lambda x: (x.top_end - x.top_start), reverse=True)
    return sites


def pcr(name, fwd, rev, template, is_circular, seed_len=15):
    """
    Run in-silico PCR. Returns dict describing the amplicon, or raises on failure.
    """
    f_sites = [s for s in map_primer(name + ":F", fwd, template, is_circular, seed_len)
               if s.strand == "+"]
    r_sites = [s for s in map_primer(name + ":R", rev, template, is_circular, seed_len)
               if s.strand == "-"]
    if not f_sites:
        raise ValueError(f"[{name}] forward primer does not anneal to template")
    if not r_sites:
        raise ValueError(f"[{name}] reverse primer does not anneal to template")

    L = len(template)
    best = None
    for f in f_sites:
        for r in r_sites:
            # express on a linear (doubled if circular) coordinate; forward 3' end
            # (f.top_end) must lie upstream of reverse 3' end (r.top_start)
            f_end = f.top_end
            r_start = r.top_start
            if is_circular and r_start < f.top_start:
                r_start += L
                r_end = r.top_end + L
            else:
                r_end = r.top_end
            if r_start < f_end:                 # reverse upstream of forward -> not a product
                if is_circular:
                    r_start += L
                    r_end += L
                else:
                    continue
            span = r_end - f.top_start
            if span <= 0 or span > (2 * L if is_circular else L):
                continue
            if best is None or span < best[0]:
                best = (span, f, r, r_end)

    if best is None:
        raise ValueError(f"[{name}] no productive forward/reverse orientation found")

    span, f, r, r_end = best
    dbl = template + template if is_circular else template
    middle = dbl[f.top_start:r_end]             # annealF + interior + annealR (template-derived)
    fwd_tail = f.tail
    rev_tail_rc = revcomp(rev.upper())[r.top_end - r.top_start:]  # right-side overhang
    amplicon = fwd_tail + middle + rev_tail_rc

    return {
        "name": name,
        "fwd_site": f,
        "rev_site": r,
        "amplicon": amplicon,
        "length": len(amplicon),
        "fwd_arm": fwd_tail,        # 5' HR arm contributed by forward primer
        "rev_arm": revcomp(rev_tail_rc),  # 5' HR arm contributed by reverse primer (as primer reads)
    }


# --------------------------------------------------------------------------- #
# recombination / HR-arm detection
# --------------------------------------------------------------------------- #
def longest_suffix_prefix(a, b, max_len=80):
    """Largest L (<=max_len) with a[-L:] == b[:L]."""
    lim = min(len(a), len(b), max_len)
    for L in range(lim, 0, -1):
        if a[-L:] == b[:L]:
            return L
    return 0


# --- Vazyme ClonExpress design window for homology arms -------------------
# The kit manual specifies 15-20 bp of homology per junction. ARM_MAX is a
# design limit, not a thermodynamic one: a longer arm still recombines, but it
# is outside the kit's validated range, so the report flags it loudly instead
# of passing it silently. Fix an over-spec arm by shortening the primer 5' tail
# or moving the junction -- never by just leaving it long.
ARM_MIN = 15
ARM_MAX = 20


def assemble_circular(frags, min_arm=ARM_MIN):
    """
    Find a circular HR assembly of N fragments (N >= 2).

    Each fragment may go in either orientation; the fragment order is searched.
    A valid circle requires every consecutive junction (incl. last->first) to
    share >= min_arm bp (suffix of one == prefix of the next).

    `frags` : list of dicts, each with 'name' and 'amplicon'.
    Returns a result dict with the best (highest total-overlap) assembly found.
    The first fragment is the anchor: kept in its given orientation and placed
    first, which removes rotational/reflection duplicates.
    """
    from itertools import permutations, product

    n = len(frags)
    if n < 2:
        raise ValueError("need at least two fragments to assemble a circle")
    seqs0 = [f["amplicon"] for f in frags]
    names = [f["name"] for f in frags]
    others = list(range(1, n))

    best = None
    for perm in permutations(others):
        order = [0] + list(perm)
        # anchor (index 0) orientation fixed to 'as-given'; flip the rest freely
        for bits in product((0, 1), repeat=n - 1):
            oris = [0] + list(bits)
            seqs = [revcomp(seqs0[fi]) if o else seqs0[fi]
                    for fi, o in zip(order, oris)]
            arms = []
            ok = True
            for k in range(n):
                ov = longest_suffix_prefix(seqs[k], seqs[(k + 1) % n])
                arms.append(ov)
                if ov < min_arm:
                    ok = False
                    break
            if not ok:
                continue
            score = sum(arms)
            if best is None or score > best["score"]:
                # build the circular sequence
                merged = seqs[0]
                for k in range(1, n):
                    merged += seqs[k][arms[k - 1]:]
                wrap = arms[-1]                       # junction last -> first
                circle = merged[:len(merged) - wrap]
                junctions = []
                for k in range(n):
                    a, b = order[k], order[(k + 1) % n]
                    junctions.append({
                        "from": names[a], "from_ori": oris[k],
                        "to": names[b], "to_ori": oris[(k + 1) % n],
                        "len": arms[k], "seq": seqs[k][-arms[k]:],
                    })
                best = {
                    "ok": True,
                    "score": score,
                    "order": [names[i] for i in order],
                    "orientations": dict(zip((names[i] for i in order),
                                            ("revcomp" if o else "as-given"
                                             for o in oris))),
                    "junctions": junctions,
                    "circle": circle,
                    "circle_len": len(circle),
                }

    if best is None:
        # report the best partial set of junctions for diagnostics (given order)
        seqs = seqs0
        diag = [{
            "from": names[k], "to": names[(k + 1) % n],
            "len": longest_suffix_prefix(seqs[k], seqs[(k + 1) % n]),
        } for k in range(n)]
        return {"ok": False, "order": names, "junctions_diag": diag}
    return best


# --------------------------------------------------------------------------- #
# self-HR mutagenesis (around-the-horn): one template, one divergent primer pair
# --------------------------------------------------------------------------- #
def longest_self_overlap(s, max_len=120):
    """Largest L (1 <= L < len(s), L <= max_len) with s[-L:] == s[:L]."""
    lim = min(len(s) - 1, max_len)
    for L in range(lim, 0, -1):
        if s[-L:] == s[:L]:
            return L
    return 0


def characterize_change(seq, circle, origin_t, origin_p):
    """
    Compare the circular mutagenesis product to the circular template by rotating
    both to a shared, guaranteed-unchanged origin (the forward primer's 3' anneal
    start), then trimming the common prefix/suffix to isolate the edited region.

    Returns (kind, tmpl_region, prod_region, start_pos, end_pos) where positions
    are 0-based template coordinates of the edited window.
    """
    L = len(seq)
    s = seq[origin_t:] + seq[:origin_t]
    c = circle[origin_p:] + circle[:origin_p]

    pre = 0
    while pre < min(len(s), len(c)) and s[pre] == c[pre]:
        pre += 1
    suf = 0
    while suf < min(len(s), len(c)) - pre and s[-1 - suf] == c[-1 - suf]:
        suf += 1

    dt = s[pre:len(s) - suf]      # template bases that were removed/replaced
    dp = c[pre:len(c) - suf]      # product bases that took their place
    start = (origin_t + pre) % L
    end = (origin_t + len(s) - suf) % L

    if not dt and not dp:
        kind = "none (silent / no change vs template)"
    elif dt and not dp:
        kind = f"deletion of {len(dt)} bp"
    elif dp and not dt:
        kind = f"insertion of {len(dp)} bp"
    elif len(dt) == len(dp):
        kind = f"substitution of {len(dt)} bp"
    else:
        kind = f"replacement {len(dt)} bp -> {len(dp)} bp"
    return kind, dt, dp, start, end


def mutagenesis(template, fwd, rev, min_arm=15, seed_len=15):
    """
    One template, one primer pair, self-HR. Linearise the plasmid by PCR (the
    around-the-horn product), then test whether the single amplicon's two ends
    share >= min_arm bp so it self-circularises. Returns a result dict.
    """
    seq, circ = template
    res = pcr("MUT", fwd[1], rev[1], seq, circ, seed_len)
    A = res["amplicon"]
    arm = longest_self_overlap(A)
    out = {
        "pcr": res,
        "amplicon": A,
        "amplicon_len": len(A),
        "arm_len": arm,
        "arm_seq": A[-arm:] if arm else "",
        "ok": arm >= min_arm,
    }
    if out["ok"]:
        circle = A[:-arm]
        out["circle"] = circle
        out["circle_len"] = len(circle)
        origin_t = res["fwd_site"].top_start
        origin_p = len(res["fwd_site"].tail)
        out["change"] = characterize_change(seq, circle, origin_t, origin_p)
    return out


# --------------------------------------------------------------------------- #
# thermodynamics
# --------------------------------------------------------------------------- #
def tm_nn(seq):
    """Nearest-neighbour Tm (SantaLucia 1998) at default 50 mM Na+, 25 nM primer."""
    try:
        return mt.Tm_NN(seq, nn_table=mt.DNA_NN3, dnac1=250, dnac2=0, Na=50)
    except Exception:
        return float("nan")


def primer_report(name, primer_seq, site):
    """Report Tm for the whole primer oligo (per user request: Tm is for the primer)."""
    a = site.anneal_seq
    return {
        "name": name,
        "primer_len": len(primer_seq),
        "match_len": len(a),
        "tail_seq": site.tail,
        "primer_tm": tm_nn(primer_seq),      # Tm of the full primer
        "primer_gc": 100 * _gc(primer_seq) if primer_seq else 0.0,
    }


def recommend_ta(tm_f, tm_r):
    """Common rule: Ta = min(primer Tm) - 5 C, clamped to a sane PCR window."""
    base = min(tm_f, tm_r) - 5.0
    return max(50.0, min(72.0, base))


# --------------------------------------------------------------------------- #
# reporting
# --------------------------------------------------------------------------- #
def run(specs, min_arm=ARM_MIN, seed_len=15, max_arm=ARM_MAX):
    """
    specs : list of fragment dicts, one per fragment / per template:
            {'label': str,
             'tpl':   (seq, is_circular),
             'tpl_name': str,
             'fwd':   (name, seq),
             'rev':   (name, seq)}
    N fragments == N templates (templates may point at the same file).
    """
    frags = []
    print("=" * 74)
    print(f"  IN-SILICO PCR   ({len(specs)} fragments / {len(specs)} templates)")
    print("=" * 74)
    for sp in specs:
        seq, circ = sp["tpl"]
        res = pcr(sp["label"], sp["fwd"][1], sp["rev"][1], seq, circ, seed_len)
        res.update(label=sp["label"], pair=(sp["fwd"], sp["rev"]),
                   tpl_name=sp["tpl_name"], tpl_len=len(seq), tpl_circ=circ)
        frags.append(res)

        print(f"\n[{sp['label']}]  template {sp['tpl_name']}  "
              f"primers {sp['fwd'][0]} / {sp['rev'][0]}")
        print(f"  amplicon length : {res['length']} bp")
        print(f"  fwd 3'-binds top strand {res['fwd_site'].top_start}..{res['fwd_site'].top_end}"
              f"  |  rev 3'-binds top strand "
              f"{res['rev_site'].top_start}..{res['rev_site'].top_end} (- strand)")
        for which, pr in (("F", res["fwd_site"]), ("R", res["rev_site"])):
            pname, pseq = sp["fwd"] if which == "F" else sp["rev"]
            rp = primer_report(pname, pseq, pr)
            tail = rp["tail_seq"] if rp["tail_seq"] else "(none, fully templated)"
            print(f"    {rp['name']:<10} primer={rp['primer_len']:>2}nt  "
                  f"template-match={rp['match_len']:>2}nt  "
                  f"Tm(primer)={rp['primer_tm']:5.1f}C  GC={rp['primer_gc']:4.1f}%  "
                  f"5'tail={tail}")

    print("\n" + "=" * 74)
    print("  ANNEALING TEMPERATURE  (whole-primer Tm, NN / SantaLucia'98)")
    print("=" * 74)
    for res in frags:
        (fn, fseq), (rn, rseq) = res["pair"]
        tf, tr = tm_nn(fseq), tm_nn(rseq)
        ta = recommend_ta(tf, tr)
        print(f"  {res['label']:<8} {fn} Tm={tf:5.1f}C | {rn} Tm={tr:5.1f}C "
              f"-> recommended Ta ~ {ta:4.1f}C")
    print("  (Ta = lower whole-primer Tm of the pair - 5 C, clamped to 50-72 C.)")

    print("\n" + "=" * 74)
    print(f"  RECOMBINATION  (Vazyme ClonExpress; arms must be "
          f"{min_arm}-{max_arm} bp per junction)")
    print("=" * 74)
    res = assemble_circular(
        [{"name": f["label"], "amplicon": f["amplicon"]} for f in frags], min_arm)

    if not res["ok"]:
        print("  Recombination: NO   -- could not form a circle with all junctions "
              f">= {min_arm} bp")
        print("  Best overlaps in the given order:")
        for j in res["junctions_diag"]:
            flag = "ok" if j["len"] >= min_arm else "TOO SHORT"
            print(f"    {j['from']} -> {j['to']}: {j['len']} bp  [{flag}]")
        print("=" * 74)
        return False

    print("  Recombination: YES  -- fragments will recombine into a circular plasmid")
    print(f"  Assembly order : {' -> '.join(res['order'])} -> (circular)")
    print("  Orientations   : "
          + ", ".join(f"{k}={v}" for k, v in res["orientations"].items()))
    print("\n  Junctions (HR arms):")
    over = []
    for i, j in enumerate(res["junctions"], 1):
        fo = "" if j["from_ori"] == 0 else "(rc)"
        to = "" if j["to_ori"] == 0 else "(rc)"
        if j["len"] > max_arm:
            over.append(i)
            flag = f"   <<< OVER SPEC (> {max_arm} bp)"
        else:
            flag = ""
        print(f"    J{i}: {j['from']}{fo} 3' -> {j['to']}{to} 5'   "
              f"{j['len']} bp   {j['seq']}{flag}")

    if over:
        js = ", ".join(f"J{i}" for i in over)
        print(f"\n  WARNING: {js} exceed the {min_arm}-{max_arm} bp ClonExpress "
              "design window.")
        print("  Shorten the 5' tail of the primer that writes that arm, or move "
              "the junction to a")
        print("  better-balanced stretch of sequence. Do not leave an arm long "
              "just to raise its Tm.")

    print(f"\n  Assembled circular product: {res['circle_len']} bp")
    # reference: if every fragment came from one circular template, compare to it
    circ_lens = {f["tpl_name"]: f["tpl_len"] for f in frags if f["tpl_circ"]}
    if len(circ_lens) == 1:
        name, ref = next(iter(circ_lens.items()))
        print(f"  (template {name} is {ref} bp circular; "
              f"delta = {res['circle_len'] - ref:+d} bp)")
    print("=" * 74)
    return True


def run_mutagenesis(template, tpl_name, pair, min_arm=ARM_MIN, seed_len=15,
                    max_arm=ARM_MAX):
    """Report a one-template, one-primer-pair self-HR mutagenesis."""
    fwd, rev = pair
    seq, circ = template
    out = mutagenesis(template, fwd, rev, min_arm, seed_len)
    res = out["pcr"]

    print("=" * 74)
    print("  SELF-HR MUTAGENESIS  (one template, one primer pair, around-the-horn)")
    print("=" * 74)
    print(f"\n  template : {tpl_name}  ({len(seq)} bp{', circular' if circ else ''})")
    print(f"  primers  : {fwd[0]} / {rev[0]}")
    print(f"  linear PCR product : {out['amplicon_len']} bp")
    print(f"  fwd 3'-binds top strand {res['fwd_site'].top_start}..{res['fwd_site'].top_end}"
          f"  |  rev 3'-binds top strand "
          f"{res['rev_site'].top_start}..{res['rev_site'].top_end} (- strand)")
    for which, pr in (("F", res["fwd_site"]), ("R", res["rev_site"])):
        pname, pseq = pair[0 if which == "F" else 1]
        rp = primer_report(pname, pseq, pr)
        tail = rp["tail_seq"] if rp["tail_seq"] else "(none)"
        print(f"    {rp['name']:<10} primer={rp['primer_len']:>2}nt  "
              f"template-match={rp['match_len']:>2}nt  "
              f"Tm(primer)={rp['primer_tm']:5.1f}C  GC={rp['primer_gc']:4.1f}%  5'tail={tail}")
    tf, tr = tm_nn(fwd[1]), tm_nn(rev[1])
    print(f"\n  recommended Ta ~ {recommend_ta(tf, tr):4.1f}C "
          f"(lower whole-primer Tm - 5 C)")

    print("\n" + "-" * 74)
    print(f"  Self-circularisation (end overlap must be "
          f"{min_arm}-{max_arm} bp):")
    if not out["ok"]:
        print(f"  RESULT: NO  -- the product's ends overlap only {out['arm_len']} bp; "
              "primers' 5' arms do not match. No self-HR.")
        print("=" * 74)
        return False
    print(f"  RESULT: YES -- product self-circularises via a "
          f"{out['arm_len']} bp HR arm")
    print(f"    HR arm = {out['arm_seq']}")
    if out["arm_len"] > max_arm:
        print(f"  WARNING: the arm is {out['arm_len']} bp, over the "
              f"{min_arm}-{max_arm} bp ClonExpress design window.")
        print("  Trim the overlapping 5' tails of the two divergent primers.")
    print(f"  Circular product: {out['circle_len']} bp "
          f"(template {len(seq)} bp; delta = {out['circle_len'] - len(seq):+d} bp)")

    kind, dt, dp, start, end = out["change"]
    print(f"\n  Edit vs template: {kind}")
    if dt or dp:
        print(f"    template coordinates ~ {start}..{end}")
        if dt:
            print(f"    removed/replaced (template): {dt if len(dt) <= 60 else dt[:57] + '...'}")
        if dp:
            print(f"    introduced (product)      : {dp if len(dp) <= 60 else dp[:57] + '...'}")
    print("=" * 74)
    return True


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _parse_named(item):
    if "=" not in item:
        raise argparse.ArgumentTypeError(f"expected NAME=SEQ, got {item!r}")
    n, s = item.split("=", 1)
    return (n, s.strip())


def _parse_fragment(tokens, default_tpl):
    """
    tokens = [LABEL[@template.dna], FWD=SEQ, REV=SEQ]
    Returns a fragment spec (template not yet loaded).
    """
    if len(tokens) != 3:
        raise argparse.ArgumentTypeError(
            "each --fragment needs: LABEL[@template.dna] FWD=SEQ REV=SEQ")
    head = tokens[0]
    if "@" in head:
        label, tpl = head.split("@", 1)
    else:
        label, tpl = head, default_tpl
    if not tpl:
        raise argparse.ArgumentTypeError(
            f"fragment {label}: no template (use LABEL@file.dna or --template)")
    return {"label": label, "tpl_path": tpl,
            "fwd": _parse_named(tokens[1]), "rev": _parse_named(tokens[2])}


# Demo: one .dna, tiled into 3 fragments amplified by 3 real in-silico primer
# pairs that share ~30 bp overlaps -> should reassemble the original plasmid.
# Resolved relative to this script so the demos work from any directory.
DEMO_TEMPLATE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "examples", "pcDNA3-BTK-nLuc.dna")


def build_tiling_demo(template_path, n=3, overlap=30, bind=25):
    seq, circ = load_template(template_path)
    L = len(seq)
    dbl = seq + seq
    bps = [round(i * L / n) for i in range(n)]
    specs = []
    for i in range(n):
        p = bps[i]
        q = bps[(i + 1) % n]
        if q <= p:
            q += L
        end = q + overlap
        fwd = (f"F{i + 1}", dbl[p:p + bind])
        rev = (f"R{i + 1}", revcomp(dbl[end - bind:end]))
        specs.append({"label": f"FRAG{i + 1}", "tpl": (seq, circ),
                      "tpl_name": template_path, "fwd": fwd, "rev": rev})
    return specs


def build_mut_demo(template_path, d1=2000, d2=2030, arm=18, bind=22):
    """
    Build a self-HR deletion-mutagenesis demo: divergent primers that delete the
    template region [d1, d2) and recreate the junction via a shared HR arm.
      F = [arm upstream of deletion] + [anneal downstream of deletion]
      R = revcomp([anneal upstream of deletion])
    The single PCR product self-circularises into the template minus (d2-d1) bp.
    """
    seq, circ = load_template(template_path)
    fwd = ("MutF", seq[d1 - arm:d1] + seq[d2:d2 + bind])
    rev = ("MutR", revcomp(seq[d1 - bind:d1]))
    return (seq, circ), template_path, (fwd, rev), (d1, d2)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Validate Vazyme ClonExpress HR cloning in silico "
                    "(2+ fragments; one template per fragment).")
    ap.add_argument("--template", help="default SnapGene .dna template for fragments "
                                       "that do not specify their own via LABEL@file.dna")
    ap.add_argument("--fragment", action="append", nargs="+", metavar="TOK",
                    help="a fragment: LABEL[@template.dna] FWD=SEQ REV=SEQ "
                         "(repeat for each fragment)")
    ap.add_argument("--mutagenesis", nargs=2, metavar="NAME=SEQ",
                    help="self-HR mutagenesis: one --template, one primer pair "
                         "FWD=SEQ REV=SEQ (divergent / around-the-horn)")
    ap.add_argument("--demo-mut", action="store_true",
                    help="run the built-in self-HR mutagenesis demo (one .dna)")
    ap.add_argument("--min-arm", type=int, default=ARM_MIN,
                    help=f"minimum HR-arm length per junction "
                         f"(default {ARM_MIN}, per the ClonExpress manual)")
    ap.add_argument("--max-arm", type=int, default=ARM_MAX,
                    help=f"maximum HR-arm length per junction before the report "
                         f"flags it as over spec (default {ARM_MAX}, per the "
                         f"ClonExpress manual)")
    ap.add_argument("--seed", type=int, default=15,
                    help="3'-anchor seed length for in-silico PCR (default 15)")
    args = ap.parse_args(argv)

    # --- mutagenesis mode: one template, one primer pair, self-HR ---
    if args.demo_mut:
        print("(--demo-mut -> built-in self-HR deletion mutagenesis on one .dna: "
              f"{DEMO_TEMPLATE}, deleting bp 2000..2030)\n")
        tpl, name, pair, _ = build_mut_demo(DEMO_TEMPLATE)
        ok = run_mutagenesis(tpl, name, pair, min_arm=args.min_arm,
                             seed_len=args.seed, max_arm=args.max_arm)
        return 0 if ok else 1
    if args.mutagenesis:
        if not args.template:
            ap.error("--mutagenesis needs --template (one template)")
        tpl = load_template(args.template)
        pair = (_parse_named(args.mutagenesis[0]), _parse_named(args.mutagenesis[1]))
        ok = run_mutagenesis(tpl, args.template, pair,
                             min_arm=args.min_arm, seed_len=args.seed,
                             max_arm=args.max_arm)
        return 0 if ok else 1

    # --- assembly mode: N fragments / N templates ---
    if not args.fragment:
        print("(no --fragment given -> running built-in 3-fragment demo on one .dna: "
              f"{DEMO_TEMPLATE})\n")
        specs = build_tiling_demo(DEMO_TEMPLATE, n=3)
    else:
        raw = [_parse_fragment(toks, args.template) for toks in args.fragment]
        cache = {}
        specs = []
        for r in raw:
            if r["tpl_path"] not in cache:
                cache[r["tpl_path"]] = load_template(r["tpl_path"])
            specs.append({"label": r["label"], "tpl": cache[r["tpl_path"]],
                          "tpl_name": r["tpl_path"], "fwd": r["fwd"], "rev": r["rev"]})

    ok = run(specs, min_arm=args.min_arm, seed_len=args.seed,
             max_arm=args.max_arm)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
