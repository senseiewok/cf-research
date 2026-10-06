# CFTR structures: what four published models contain and how we read them

What the four structure files of human CFTR that the lab's 3D model uses actually say about themselves, what UniProt adds, how they line up, and what we deliberately do not claim. The model and its tests are in [`tools/cftr-model/`](../tools/cftr-model/README.md). Every file and number below was read or computed on 2026-10-05 from the files in `sources/downloads/` (fetched with [`tools/sources/fetch_sources.py`](../tools/sources/README.md); the catalog entries are `rcsb-pdb-6msm`, `rcsb-pdb-5uak`, `rcsb-pdb-8eiq`, `rcsb-pdb-8ej1` and `uniprot-p13569`, and the manifest holds their SHA-256).

Not medical advice. A structure is a still picture of one arrangement of a protein under one set of experimental conditions. Nothing here describes how any medicine acts, and nothing here should change anyone's care.

## The four structures

| Entry | What the file calls it | Method | Resolution | Deposited | Paper (journal, year, PubMed) |
| --- | --- | --- | --- | --- | --- |
| 6MSM | Phosphorylated, ATP-bound human CFTR | cryo-EM | 3.2 Å | 2018-10-16 | PNAS 2018, 30459277 |
| 5UAK | Dephosphorylated, ATP-free human CFTR | cryo-EM | 3.87 Å | 2016-12-19 | Cell 2017, 28340353 |
| 8EIQ | Phosphorylated human delta F508 CFTR with the three components of Trikafta and ATP/Mg | cryo-EM | 3 Å | 2022-09-15 | Science 2022, 36264792 |
| 8EJ1 | Dephosphorylated human delta F508 CFTR | cryo-EM | 6.9 Å | 2022-09-16 | Science 2022, 36264792 |

`[verified]` Source: the header fields `_struct.title`, `_exptl.method`, `_em_3d_reconstruction.resolution`, `_pdbx_database_status.recvd_initial_deposition_date` and `_citation.*` of each file, read by command 2026-10-05. The wording in the second column paraphrases each file's own title.

8EJ1 has the lowest resolution of the four, which limits how much detail its model carries. `[verified]` (resolution field); the consequence for interpretation is our reading, `[hypothesis]`.

## What each file says about its own sequence

Each file carries a list of where its protein differs from the UniProt reference sequence. Read directly (`_struct_ref_seq_dif`):

- 6MSM, 8EIQ and 8EJ1 each record an engineered substitution at position 1371: glutamine in the file where UniProt has glutamate. 5UAK does not. `[verified]`
- 8EIQ and 8EJ1 lack phenylalanine 508. The 8EIQ file calls this a deletion; the 8EJ1 file calls it a variant and its entity description states a deletion of F508 together with the E1371Q substitution. `[verified]`
- 6MSM and 5UAK each list 9 expression-tag residues beyond the end of the protein; 8EIQ and 8EJ1 list none. These are not drawn. `[verified]`
- The files do not say why residue 1371 was changed. We do not say either. `[verified]` that they are silent.
- Our reading: because 8EIQ also carries E1371Q, the Trikafta components and ATP with magnesium, a difference between 6MSM and 8EIQ cannot be attributed to the deletion alone. `[hypothesis]`

## Numbering

The residue numbers in the model are the entries' author numbers. They equal UniProt's numbers: every one of the 4,609 residues the four models place was compared with the amino acid UniProt has at that number, and the only disagreements are the 3 declared substitutions at 1371 above. `[verified]` by `tools/cftr-model/build_traces.py`, which stops if one disagrees without its file declaring it.

For the two deletion entries the file's own reference table counts the deleted residue out: it lists the changed residue as 1370 in the file's sequence and 1371 in UniProt. The coordinates keep the UniProt number. `[verified]`

## What each model places

Each model places only part of the 1,480-residue protein. 6MSM places 1,181 (299 not placed), 5UAK 1,139 (341), 8EIQ 1,162 (318) and 8EJ1 1,127 (353). None of them places any residue of the stretch UniProt marks as a disordered regulatory region (residues 654 to 831) in the chain the model draws (chain A, model 1). Three of the files also hold a short separate fragment of unknown identity whose place in the sequence is not assigned: 17 residues in 6MSM (chain B), 19 in 5UAK (chain R; the file names it 'R domain') and 16 in 8EIQ (chain B); 8EJ1 has none. The model does not draw these fragments. `[verified]` (computed from the residue numbers that have an alpha-carbon atom in chain A, model 1).

In 8EIQ and 8EJ1 the chain runs from residue 507 straight to 509. The gap at 508 is the only one-residue gap in 8EIQ; its other gaps between placed stretches are 16 to 209 residues. `[verified]`

## How the models line up

We fit each model onto 6MSM using the alpha-carbon atoms in the two membrane-spanning parts (UniProt residues 81 to 365 and 859 to 1155) that 6MSM and the other model both place, at the same residue number and with the same amino acid, by a least-squares rigid rotation and shift (Horn's method). The remaining root-mean-square distance:

| Model | Residues used | Distance left |
| --- | --- | --- |
| 8EIQ | 566 | 0.84 Å |
| 5UAK | 557 | 5.01 Å |
| 8EJ1 | 547 | 4.87 Å |

`[verified]` Computed for the model; the method is proved on a known rotation in `tools/cftr-model/tests/test_data.py`. These numbers describe differences in shape between four separate experiments. They do not say why the shapes differ, they do not say that any state is open or closed, and they do not show what any compound does. `[hypothesis]` about how to read them, and the reason the page says so.

## What UniProt adds (entry P13569)

Entry version 286, annotation last updated 2026-09-02, sequence length 1,480. `[verified]`

- Four marked parts: two membrane-spanning (residues 81 to 365 and 859 to 1155) and two nucleotide-binding, which UniProt names as ABC transporter domains (423 to 646 and 1210 to 1443). Its keywords include ATP-binding and nucleotide-binding, and its ATP-binding site annotations fall inside those two parts. `[verified]`
- 12 transmembrane segments and 13 topological domains (cytoplasmic or extracellular), and 118 helix, strand and turn ranges (76 helix, 31 strand, 11 turn), each with the PDB entries UniProt cites. These are UniProt's annotations; the lab does not compute them from the four models. `[verified]`
- 209 natural-variant entries: 196 change one residue (at 165 different positions) and 13 span several residues. The entry for F508del carries a long description with several functional statements and many citations. `[verified]` as UniProt's wording and counts. The lab shows the position and UniProt's own words with attribution and the caveat below, and makes no statement of its own about any variant.
- UniProt's licence text applies CC BY 4.0 to its databases and adds that medical and genetic information is for research and education, not a substitute for medical advice. `[verified]` (help text, last modified 2024-12-18, read 2026-10-05).
- Not read here: CFTR2, the CF-specific authority on variant consequences, and any registry count of how common F508del is. The model therefore says nothing about frequency, severity or consequence. `[verified]` that they were not read.

## Licences and conduct

- The wwPDB policy, as stated on the RCSB policy page read directly on 2026-10-05, makes the archive's data files available under CC0 1.0 and encourages crediting the original authors. `[verified]`
- `files.rcsb.org` has no robots.txt (a 404) and its download service is documented for scripted downloads; the repo's fetcher took the four files one at a time with its project user agent. `[verified]`
- `data.rcsb.org/rest` is disallowed for all crawlers by that host's robots.txt (read 2026-10-05), so nothing here uses that API. A research sub-agent used it during the session before this was noticed; the slip is in the [ledger](../ledger/2026-10-04-first-session.md). `[verified]`

## What we do not claim

These were asked for or suggested during the work and are not stated by any source we read directly, so the model and its page leave them out. `[unverified]` unless noted.

- Where the pore and the gate are, how chloride passes, and whether any of the four is an open or a closed channel.
- How much F508del protein reaches the cell surface, and any percentage of people with the variant.
- How any medicine acts, or what the compounds in 8EIQ do.
- A membrane plane through the protein. The sides come from UniProt's topology annotation, and the page says so.

Sources that describe these (a patient-education page and a review summarised by a fetch tool) were read only through summaries, so they stay `[unverified]` until someone reads them on the page.

## How this was made and checked

The process record, including the blind reviews and what they found, is in [`proposals/2026-10-05-cftr-model-review.md`](../proposals/2026-10-05-cftr-model-review.md).
