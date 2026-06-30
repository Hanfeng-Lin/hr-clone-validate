# hr-clone-validate

In-silico validation of **Vazyme ClonExpress / C116-style homologous-recombination
(HR) cloning** directly from SnapGene `.dna` files.

Given your templates and primers, it performs in-silico PCR, finds the homology
(HR) arms that drive recombination, decides whether the fragments assemble into a
circular plasmid, and reports recommended annealing temperatures — all from the
command line, no SnapGene cloning wizard required.

It covers three workflows:

1. **Multi-fragment assembly** — N fragments (each with its own template + primer
   pair) joined into one circle. Handles 2-fragment inserts, 3-fragment (and
   higher) assemblies, and automatically figures out fragment **order** and
   **orientation**.
2. **Self-HR mutagenesis** (around-the-horn) — one template, one divergent primer
   pair; the single linear PCR product self-circularises because its two ends
   share an HR arm. Classifies the edit as a point **substitution**,
   **insertion**, or **deletion**.
3. **In-silico PCR + Tm** — for every primer: the template-binding region, the
   5′ non-templated tail (= the HR arm), whole-primer melting temperature
   (nearest-neighbour, SantaLucia 1998) and a recommended PCR annealing
   temperature.

---

## Requirements

- Python 3.8+
- [`snapgene_reader`](https://pypi.org/project/snapgene-reader/)
- [`biopython`](https://biopython.org/) (1.79+ works; uses `Tm_NN`)

```bash
pip install -r requirements.txt
```

---

## Quick start

```bash
# built-in 3-fragment demo (tiles one plasmid into 3 amplicons and reassembles it)
python hr_clone_validate.py

# built-in self-HR mutagenesis demo (deletes 30 bp from one plasmid)
python hr_clone_validate.py --demo-mut
```

---

## Usage

### Multi-fragment assembly

Each fragment is one PCR reaction with its own template. **N fragments = N
templates** (templates may point at the same file).

```bash
python hr_clone_validate.py \
    --fragment VECTOR  VecF=<seq>  vecR=<seq> \
    --fragment INSERT1 aF=<seq>    aR=<seq> \
    --template pcDNA3-BTK-nLuc.dna          # default template for fragments without @
```

Give a fragment its **own** template with `LABEL@path.dna`:

```bash
python hr_clone_validate.py \
    --fragment VECTOR@pcDNA3-BTK-HiBiT.dna   VecF=<seq>     vecR=<seq> \
    --fragment INSERT@../vav1/pRP-nLuc-vav1.dna BTKnlucF=<seq> nlucR=<seq>
```

- Each `--fragment` takes exactly three tokens: `LABEL[@template.dna]`,
  `FWDNAME=SEQ`, `REVNAME=SEQ`.
- Mix freely: fragments without `@` use the `--template` default.
- The assembler tries every order and per-fragment orientation; it reports the
  one that closes the circle with every junction ≥ `--min-arm` bp.

### Self-HR mutagenesis (one template, one primer pair)

```bash
python hr_clone_validate.py --template pcDNA3-BTK-nLuc.dna \
    --mutagenesis MutF=<seq> MutR=<seq>
```

The two primers diverge (point away from each other / around the plasmid). The
single linear product self-circularises through the HR arm formed by the primers'
overlapping 5′ ends. The tool reports the arm, the circular product size, and the
edit relative to the template.

### Options

| flag | default | meaning |
|------|---------|---------|
| `--template PATH` | – | default `.dna` template for fragments without their own `@file` |
| `--fragment TOK…` | – | one fragment: `LABEL[@template.dna] FWD=SEQ REV=SEQ` (repeatable) |
| `--mutagenesis F=SEQ R=SEQ` | – | self-HR mutagenesis mode (needs `--template`) |
| `--demo-mut` | – | run the built-in deletion-mutagenesis demo |
| `--min-arm N` | 15 | minimum HR-arm length required at each junction |
| `--seed N` | 15 | 3′-anchor seed length used by the in-silico PCR |

Exit code is `0` when a valid circular product is found, `1` otherwise (useful in
scripts).

---

## Example output

```
==========================================================================
  RECOMBINATION  (Vazyme ClonExpress; need >= 15 bp per junction)
==========================================================================
  Recombination: YES  -- fragments will recombine into a circular plasmid
  Assembly order : VECTOR -> INSERT -> (circular)
  Orientations   : VECTOR=as-given, INSERT=as-given

  Junctions (HR arms):
    J1: VECTOR 3' -> INSERT 5'   18 bp   GGATCGTCAAGTGGCGGT
    J2: INSERT 3' -> VECTOR 5'   18 bp   TGCGAACGCATTCTGGCG

  Assembled circular product: 7888 bp
```

Mutagenesis:

```
  RESULT: YES -- product self-circularises via a 18 bp HR arm
  Circular product: 7858 bp (template 7888 bp; delta = -30 bp)
  Edit vs template: deletion of 30 bp
    removed/replaced (template): CTGCAGGACTCATATCCAGGCTCAAATATC
```

---

## How it works

1. **In-silico PCR.** Each primer's 3′ end is anchored to the template (forward
   on the top strand, reverse on the bottom strand) and extended to find the
   longest contiguous match. Circular templates are searched on a doubled
   sequence so primers spanning the origin still bind. The amplicon = forward
   primer 5′ tail + templated interior + reverse-complemented reverse primer.
   The part of a primer that does **not** match the template is the added HR arm.
2. **Junction detection.** For assembly, the tool searches all fragment orders
   and orientations for a circle where every junction shares
   `suffix(fragment_i) == prefix(fragment_{i+1})` of at least `--min-arm` bp. For
   mutagenesis it looks for the largest self-overlap between the single product's
   two ends.
3. **Edit classification (mutagenesis).** The circular product and template are
   rotated to a shared, guaranteed-unchanged origin (the forward primer's 3′
   anneal start), then the common prefix/suffix is trimmed to isolate the edited
   window and classify it (substitution / insertion / deletion / replacement).
4. **Thermodynamics.** Whole-primer Tm via Biopython `Tm_NN` (SantaLucia 1998).
   Recommended annealing temperature = lower primer Tm of the pair − 5 °C,
   clamped to 50–72 °C.

---

## Notes & caveats

- Tm is computed for the **whole primer oligo**. When a primer is fully templated
  (e.g. a self-reconstitution test that reuses one plasmid for every fragment),
  the reported 5′ tail is `(none)` because nothing is non-templated — the real HR
  arms still appear in the junction report.
- The built-in demos auto-design fully-templated primers; for real cloning your
  inserts' foreign HR arms will appear as each primer's `5'tail`.
- This is a sequence-level check. It does not model PCR efficiency, secondary
  structure, primer dimers, or recombination efficiency — confirm long or
  GC-skewed primers experimentally.

---

## License

MIT (see `LICENSE`), unless your lab specifies otherwise.
