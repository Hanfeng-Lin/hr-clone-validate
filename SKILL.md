---
name: hr-clone-validate
description: >
  Validate Vazyme ClonExpress / C116-style homologous-recombination (HR) cloning
  in silico from SnapGene .dna files. Use when the user wants to check whether
  PCR fragments will recombine into a plasmid, find the HR (homology) arms at each
  junction, get the assembled product size, run an around-the-horn self-HR
  mutagenesis check (point mutation / insertion / deletion), or compute primer Tm
  and annealing temperature for HR cloning.
---

# hr-clone-validate skill

A command-line tool (`hr_clone_validate.py`) that simulates HR cloning from
SnapGene `.dna` templates and primer sequences.

## When to use

Invoke this skill when the user asks to:

- check whether a set of PCR fragments will **recombine / assemble** into a
  circular plasmid (Gibson / ClonExpress / In-Fusion style);
- find the **HR arm / overlap** shared between fragments at each junction;
- predict the **assembled plasmid size** and the assembly **order/orientation**;
- validate a **self-HR mutagenesis** (around-the-horn) design — one template, one
  divergent primer pair — and classify the edit (substitution / insertion /
  deletion);
- get **whole-primer Tm** and a **recommended annealing temperature** for HR
  cloning primers.

## Prerequisites

- Python 3.8+, `snapgene_reader`, `biopython` (`pip install -r requirements.txt`).
- The templates as SnapGene `.dna` files (paths can be relative, e.g. `../foo.dna`).
- The primer sequences (5′→3′). Case-insensitive; lowercase tails are fine.

## Core model

- **N fragments = N templates.** Each fragment is one PCR reaction with its own
  template and primer pair; templates may repeat (point at the same file).
- A primer's 3′ region anneals to the template; its 5′ non-templated part is the
  **HR arm**. Two fragments join where `suffix(A) == prefix(B)` ≥ `--min-arm` bp.
- For **mutagenesis**, the single linear product self-circularises through the
  overlap between its two ends.

## How to run

The script ships next to this `SKILL.md`. Run it by its full path inside the
skill directory (e.g. `~/.claude/skills/hr-clone-validate/hr_clone_validate.py`
or `~/.codex/skills/hr-clone-validate/hr_clone_validate.py`) while staying in the
user's working directory, so relative `.dna` paths resolve. Below,
`hr_clone_validate.py` stands for that full path.

Assembly (each `--fragment` = `LABEL[@template.dna] FWD=SEQ REV=SEQ`):

```bash
python hr_clone_validate.py \
    --fragment VECTOR@vector.dna  VecF=<seq> vecR=<seq> \
    --fragment INSERT@insert.dna  insF=<seq> insR=<seq>
```

Use `--template default.dna` for fragments written without an `@template`.

Self-HR mutagenesis (one template, one pair):

```bash
python hr_clone_validate.py --template plasmid.dna --mutagenesis MutF=<seq> MutR=<seq>
```

Demos (no user inputs needed):

```bash
python hr_clone_validate.py            # 3-fragment assembly demo (one .dna)
python hr_clone_validate.py --demo-mut # deletion mutagenesis demo (one .dna)
```

Options: `--min-arm N` (default 15, minimum overlap per junction),
`--seed N` (default 15, 3′-anchor seed for in-silico PCR). Exit code `0` =
valid circular product found, `1` = not.

## Interpreting the output

1. **IN-SILICO PCR** — per fragment: amplicon length, where each primer binds,
   primer length / template-match length / whole-primer Tm / GC, and the 5′ tail
   (= HR arm; `(none, fully templated)` means the whole primer matches template).
2. **ANNEALING TEMPERATURE** — per fragment, each primer's Tm and the recommended
   Ta (lower primer Tm − 5 °C, clamped 50–72 °C).
3. **RECOMBINATION** — YES/NO, the assembly order, per-fragment orientation
   (`as-given` / `revcomp`), every junction's HR-arm length + sequence, and the
   assembled circular size (compared to the template when all fragments share one
   circular template). On failure it prints the best overlap per junction so you
   can see which arm is too short.

## Guidance for the agent

- **Always confirm what the insert actually contains.** Primers may amplify only
  part of a template (e.g. a tag, not the whole gene). After a successful
  assembly, check the amplicon's template coordinates against the template's
  feature annotations and tell the user what the product really encodes — don't
  assume the construct name equals the insert.
- If a primer fails to anneal, `pcr()` raises; report which primer and on which
  template, and suggest checking the template choice or the 3′ binding region.
- If recombination is NO, point to the junction whose overlap is below
  `--min-arm` and suggest lengthening that primer's 5′ HR arm to ~15–20 bp.
- Vazyme ClonExpress recommends ~15–20 bp homology per junction; `--min-arm 15`
  is the default floor, raise it for a stricter check.
- This is a sequence-level prediction only — no PCR efficiency, primer-dimer, or
  secondary-structure modelling.
