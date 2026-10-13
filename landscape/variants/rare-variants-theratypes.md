# Cryo-EM 3D Perturbation Dossier & Rare Variant Theratyping

**Date**: 2026-10-13 | **Track**: Track A (Biophysical Science & Cure Discovery) | **Answer tier**: T2 (computational geometry from primary PDB coordinates)

## 1. Executive Summary & Biophysical Scope

Over 2,000 CFTR variants exist in the clinical population. While high-frequency mutations like F508del (Class II), G551D (Class III), and G542X/W1282X (Class I) have received targeted drug discovery or clinical trials, hundreds of rare missense variants remain completely unclassified. Consequently, individuals bearing these genotypes are frequently excluded from modulator labels.

Using 4 high-resolution published Cryo-EM models superposed onto active phosphorylated human CFTR (`6MSM`):
- `6MSM`: Phosphorylated, ATP-bound active state (3.2 Å)
- `5UAK`: Dephosphorylated, ATP-free inward-facing inactive state (3.87 Å)
- `8EIQ`: Phosphorylated human F508del with Trikafta components (elexacaftor, tezacaftor, ivacaftor) and ATP/Mg (3.0 Å)
- `8EJ1`: Dephosphorylated apo F508del (6.9 Å)

This tool computes the coordinate perturbation deltas ($\Delta d$), domain interface topology, and spatial distance to therapeutic allosteric pockets (ICL4 coupling, TMD1 tezacaftor cleft, TMD2 elexacaftor pocket, and Walker A/B catalytic sites).

---

## 2. In Silico Structural Perturbation Metrics

| Residue | WT AA | Domain | Secondary Structure | $\Delta d$ Active-Inactive (6MSM-5UAK) | $\Delta d$ Active-Trikafta (6MSM-8EIQ) | ICL4 Proximity | Predicted Theratype |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **117** | R | ABC transmembrane type-1 1 | Helix | 6.14 Å | 0.45 Å | 77.2 Å | Class II / Cryptic Pocket Candidate |
| **334** | R | ABC transmembrane type-1 1 | Helix | 0.71 Å | 0.59 Å | 72.38 Å | Class IV |
| **347** | R | ABC transmembrane type-1 1 | Helix | 0.95 Å | 0.86 Å | 53.61 Å | Class IV |
| **508** | F | ABC transporter 1 | Loop / coil | 9.24 Å | Unplaced | 9.31 Å | Class II |
| **551** | G | ABC transporter 1 | Helix | 9.32 Å | 1.46 Å | 23.58 Å | Class II |
| **560** | R | ABC transporter 1 | Helix | 8.87 Å | 2.22 Å | 13.18 Å | Class II |
| **563** | Y | ABC transporter 1 | Helix | 8.79 Å | 2.47 Å | 16.07 Å | Class II |
| **1066** | R | ABC transmembrane type-1 2 | Helix | 8.28 Å | 0.71 Å | 5.61 Å | Class II |
| **1070** | R | ABC transmembrane type-1 2 | Loop / coil | 8.2 Å | 1.19 Å | 3.86 Å | Class II |
| **1282** | W | ABC transporter 2 | Helix | 6.94 Å | 1.46 Å | 44.5 Å | Class II / IV |
| **1303** | N | ABC transporter 2 | Helix | 11.61 Å | 0.67 Å | 31.94 Å | Class II / IV |

---

## 3. Detailed Structural Profiles for High-Priority Variants

### Residue R117 Profile

- **Domain assignment**: `ABC transmembrane type-1 1 (81-365)`
- **Secondary structure**: `Helix (116-118)`
- **Conformational displacement**: Active-to-Inactive delta = `6.14 Å`; Active-to-Trikafta delta = `0.45 Å`.
- **Allosteric Landmark Proximity (6MSM)**:
  - `ICL4_coupling`: **77.2 Å**
  - `ICL1_coupling`: **73.19 Å**
  - `ATP_site_1`: **98.25 Å**
  - `ATP_site_2`: **96.5 Å**
  - `TMD1_pocket`: **79.38 Å**
  - `TMD2_pocket`: **71.52 Å**
  - `pore_constriction`: **22.95 Å**
- **Theratyping Classification**: `Class II / Cryptic Pocket Candidate (Transmembrane domain)`
- **Mechanistic Rationale**: Transmembrane helix packing defect; candidate for novel hydrophobic pocket stabilization.
- **Known UniProt Clinical Variations**:
  - `VAR ID: VAR_000121`: R -> C (in CF; dbSNP:rs77834169)
  - `VAR ID: VAR_000122`: R -> H (in CF and CBAVD; strong decrease in single channel conductance; promotes rapid return to the closed state of the channel; decrease in bicarbonate transport; dbSNP:rs78655421)
  - `VAR ID: VAR_000123`: R -> L (in CF; dbSNP:rs78655421)
  - `VAR ID: VAR_000124`: R -> P (in CF; dbSNP:rs78655421)

### Residue R334 Profile

- **Domain assignment**: `ABC transmembrane type-1 1 (81-365)`
- **Secondary structure**: `Helix (333-348)`
- **Conformational displacement**: Active-to-Inactive delta = `0.71 Å`; Active-to-Trikafta delta = `0.59 Å`.
- **Allosteric Landmark Proximity (6MSM)**:
  - `ICL4_coupling`: **72.38 Å**
  - `ICL1_coupling`: **66.3 Å**
  - `ATP_site_1`: **90.41 Å**
  - `ATP_site_2`: **90.63 Å**
  - `TMD1_pocket`: **73.61 Å**
  - `TMD2_pocket`: **66.7 Å**
  - `pore_constriction`: **11.26 Å**
- **Theratyping Classification**: `Class IV (Reduced Conductance / Pore Architecture)`
- **Mechanistic Rationale**: Located 11.26 Å from channel pore constriction; predicted ion permeation / selectivity reduction.
- **Known UniProt Clinical Variations**:
  - `VAR ID: VAR_000148`: R -> W (in CF and CBAVD; mild; does not prevent maturation of glycans; dbSNP:rs121909011)

### Residue R347 Profile

- **Domain assignment**: `ABC transmembrane type-1 1 (81-365)`
- **Secondary structure**: `Helix (333-348)`
- **Conformational displacement**: Active-to-Inactive delta = `0.95 Å`; Active-to-Trikafta delta = `0.86 Å`.
- **Allosteric Landmark Proximity (6MSM)**:
  - `ICL4_coupling`: **53.61 Å**
  - `ICL1_coupling`: **46.29 Å**
  - `ATP_site_1`: **70.18 Å**
  - `ATP_site_2`: **73.0 Å**
  - `TMD1_pocket`: **54.01 Å**
  - `TMD2_pocket`: **48.25 Å**
  - `pore_constriction`: **9.44 Å**
- **Theratyping Classification**: `Class IV (Reduced Conductance / Pore Architecture)`
- **Mechanistic Rationale**: Located 9.44 Å from channel pore constriction; predicted ion permeation / selectivity reduction.
- **Known UniProt Clinical Variations**:
  - `VAR ID: VAR_000153`: R -> H (in CF; dbSNP:rs77932196)
  - `VAR ID: VAR_000154`: R -> L (in CF; dbSNP:rs77932196)
  - `VAR ID: VAR_000155`: R -> P (in CF; mild; dbSNP:rs77932196)

### Residue F508 Profile

- **Domain assignment**: `ABC transporter 1 (423-646)`
- **Secondary structure**: `Loop / coil`
- **Conformational displacement**: Active-to-Inactive delta = `9.24 Å`; Active-to-Trikafta delta = `None Å`.
- **Allosteric Landmark Proximity (6MSM)**:
  - `ICL4_coupling`: **9.31 Å**
  - `ICL1_coupling`: **21.71 Å**
  - `ATP_site_1`: **29.07 Å**
  - `ATP_site_2`: **27.02 Å**
  - `TMD1_pocket`: **22.06 Å**
  - `TMD2_pocket`: **11.98 Å**
  - `pore_constriction`: **67.49 Å**
- **Theratyping Classification**: `Class II (Misfolding / Processing - ICL4 interface)`
- **Mechanistic Rationale**: Located 9.31 Å from ICL4 coupling loop; high risk of NBD1:ICL4 misassembly analogous to F508del. Candidate for Class 1/2 correctors (VX-661/VX-445).
- **Known UniProt Clinical Variations**:
  - `VAR ID: VAR_000172`: F -> C (in dbSNP:rs74571530)
  - `VAR ID: VAR_000171`: None ->  (in CF and CBAVD; most common mutation in Caucasian CF chromosomes; impairs protein folding and stability; causes local changes to the surface that mediates interactions between domains; decreases frequency of channel opening in vitro; binds to the cytokeratin-8 and through this binding is primed for the degradation pathway that ends in the proteasome thus impairing trafficking; impairs maturation and trafficking to the cell membrane; impairs recycling to the cell membrane after endocytosis; abolishes MARCHF2-mediated degradation; dbSNP:rs121909001 and dbSNP:rs113993960)

### Residue G551 Profile

- **Domain assignment**: `ABC transporter 1 (423-646)`
- **Secondary structure**: `Helix (550-563)`
- **Conformational displacement**: Active-to-Inactive delta = `9.32 Å`; Active-to-Trikafta delta = `1.46 Å`.
- **Allosteric Landmark Proximity (6MSM)**:
  - `ICL4_coupling`: **23.58 Å**
  - `ICL1_coupling`: **30.23 Å**
  - `ATP_site_1`: **21.74 Å**
  - `ATP_site_2`: **8.79 Å**
  - `TMD1_pocket`: **33.67 Å**
  - `TMD2_pocket`: **25.16 Å**
  - `pore_constriction`: **76.76 Å**
- **Theratyping Classification**: `Class II (NBD1 structural stability)`
- **Mechanistic Rationale**: Located within NBD1 core; predicted thermal destabilization.
- **Known UniProt Clinical Variations**:
  - `VAR ID: VAR_000179`: G -> D (in CF; decrease in the frequency of channel opening in vitro; decrease in channel activity and ATPase activity; complete loss of bicarbonate transport; no effect on trafficking to the cell membrane; protein stability; nor on the maturation of glycans; dbSNP:rs75527207)
  - `VAR ID: VAR_000180`: G -> S (in CF; decrease in bicarbonate transport; no effect on chloride channel activity; dbSNP:rs121909013)

### Residue R560 Profile

- **Domain assignment**: `ABC transporter 1 (423-646)`
- **Secondary structure**: `Helix (550-563)`
- **Conformational displacement**: Active-to-Inactive delta = `8.87 Å`; Active-to-Trikafta delta = `2.22 Å`.
- **Allosteric Landmark Proximity (6MSM)**:
  - `ICL4_coupling`: **13.18 Å**
  - `ICL1_coupling`: **24.11 Å**
  - `ATP_site_1`: **24.27 Å**
  - `ATP_site_2`: **22.33 Å**
  - `TMD1_pocket`: **24.04 Å**
  - `TMD2_pocket`: **16.99 Å**
  - `pore_constriction`: **72.94 Å**
- **Theratyping Classification**: `Class II (Misfolding / Processing - ICL4 interface)`
- **Mechanistic Rationale**: Located 13.18 Å from ICL4 coupling loop; high risk of NBD1:ICL4 misassembly analogous to F508del. Candidate for Class 1/2 correctors (VX-661/VX-445).
- **Known UniProt Clinical Variations**:
  - `VAR ID: VAR_000184`: R -> K (in CF; dbSNP:rs80055610)
  - `VAR ID: VAR_000185`: R -> S (in CF; dbSNP:rs397508267)
  - `VAR ID: VAR_000186`: R -> T (in CF; impairs maturation and trafficking to the cell membrane; decrease in channel activity; dbSNP:rs80055610)

### Residue Y563 Profile

- **Domain assignment**: `ABC transporter 1 (423-646)`
- **Secondary structure**: `Helix (550-563)`
- **Conformational displacement**: Active-to-Inactive delta = `8.79 Å`; Active-to-Trikafta delta = `2.47 Å`.
- **Allosteric Landmark Proximity (6MSM)**:
  - `ICL4_coupling`: **16.07 Å**
  - `ICL1_coupling`: **27.71 Å**
  - `ATP_site_1`: **27.12 Å**
  - `ATP_site_2`: **25.54 Å**
  - `TMD1_pocket`: **25.69 Å**
  - `TMD2_pocket`: **20.49 Å**
  - `pore_constriction`: **76.46 Å**
- **Theratyping Classification**: `Class II (Misfolding / Processing - ICL4 interface)`
- **Mechanistic Rationale**: Located 16.07 Å from ICL4 coupling loop; high risk of NBD1:ICL4 misassembly analogous to F508del. Candidate for Class 1/2 correctors (VX-661/VX-445).
- **Known UniProt Clinical Variations**:
  - `VAR ID: VAR_000189`: Y -> N (in CF; dbSNP:rs121909006)

### Residue R1066 Profile

- **Domain assignment**: `ABC transmembrane type-1 2 (859-1155)`
- **Secondary structure**: `Helix (1062-1068)`
- **Conformational displacement**: Active-to-Inactive delta = `8.28 Å`; Active-to-Trikafta delta = `0.71 Å`.
- **Allosteric Landmark Proximity (6MSM)**:
  - `ICL4_coupling`: **5.61 Å**
  - `ICL1_coupling`: **12.23 Å**
  - `ATP_site_1`: **23.06 Å**
  - `ATP_site_2`: **32.95 Å**
  - `TMD1_pocket`: **10.09 Å**
  - `TMD2_pocket`: **11.73 Å**
  - `pore_constriction`: **64.94 Å**
- **Theratyping Classification**: `Class II (TMD2 folding / VX-445 pocket proximity)`
- **Mechanistic Rationale**: Proximity (11.73 Å) to the TMD2 allosteric pocket; predicted responsive to elexacaftor binding.
- **Known UniProt Clinical Variations**:
  - `VAR ID: VAR_000236`: R -> C (in CF; dbSNP:rs78194216)
  - `VAR ID: VAR_000237`: R -> H (in CF; dbSNP:rs121909019)
  - `VAR ID: VAR_000238`: R -> L (in CF; dbSNP:rs121909019)

### Residue R1070 Profile

- **Domain assignment**: `ABC transmembrane type-1 2 (859-1155)`
- **Secondary structure**: `Loop / coil`
- **Conformational displacement**: Active-to-Inactive delta = `8.2 Å`; Active-to-Trikafta delta = `1.19 Å`.
- **Allosteric Landmark Proximity (6MSM)**:
  - `ICL4_coupling`: **3.86 Å**
  - `ICL1_coupling`: **15.81 Å**
  - `ATP_site_1`: **28.46 Å**
  - `ATP_site_2`: **33.13 Å**
  - `TMD1_pocket`: **14.42 Å**
  - `TMD2_pocket`: **9.19 Å**
  - `pore_constriction`: **63.82 Å**
- **Theratyping Classification**: `Class II (TMD2 folding / VX-445 pocket proximity)`
- **Mechanistic Rationale**: Proximity (9.19 Å) to the TMD2 allosteric pocket; predicted responsive to elexacaftor binding.
- **Known UniProt Clinical Variations**:
  - `VAR ID: VAR_000242`: R -> P (in CF; dbSNP:rs78769542)
  - `VAR ID: VAR_000241`: R -> Q (in CF; decrease in bicarbonate transport; no effect on chloride channel activity; dbSNP:rs78769542)
  - `VAR ID: VAR_011564`: R -> W (in CBAVD; dbSNP:rs202179988)

### Residue W1282 Profile

- **Domain assignment**: `ABC transporter 2 (1210-1443)`
- **Secondary structure**: `Helix (1279-1283)`
- **Conformational displacement**: Active-to-Inactive delta = `6.94 Å`; Active-to-Trikafta delta = `1.46 Å`.
- **Allosteric Landmark Proximity (6MSM)**:
  - `ICL4_coupling`: **44.5 Å**
  - `ICL1_coupling`: **44.32 Å**
  - `ATP_site_1`: **39.5 Å**
  - `ATP_site_2`: **25.03 Å**
  - `TMD1_pocket`: **53.18 Å**
  - `TMD2_pocket`: **42.2 Å**
  - `pore_constriction`: **72.98 Å**
- **Theratyping Classification**: `Class II / IV (NBD2 misassembly or gating transmission)`
- **Mechanistic Rationale**: Located in NBD2; potential inter-NBD dimer destabilization.
- **Known UniProt Clinical Variations**:
  - `VAR ID: VAR_080324`: None ->  (in CF and CBAVD)
  - `VAR ID: VAR_000261`: W -> R (in CF; dbSNP:rs397508616)

### Residue N1303 Profile

- **Domain assignment**: `ABC transporter 2 (1210-1443)`
- **Secondary structure**: `Helix (1300-1304)`
- **Conformational displacement**: Active-to-Inactive delta = `11.61 Å`; Active-to-Trikafta delta = `0.67 Å`.
- **Allosteric Landmark Proximity (6MSM)**:
  - `ICL4_coupling`: **31.94 Å**
  - `ICL1_coupling`: **26.9 Å**
  - `ATP_site_1`: **23.37 Å**
  - `ATP_site_2`: **28.73 Å**
  - `TMD1_pocket`: **35.46 Å**
  - `TMD2_pocket`: **31.06 Å**
  - `pore_constriction`: **65.32 Å**
- **Theratyping Classification**: `Class II / IV (NBD2 misassembly or gating transmission)`
- **Mechanistic Rationale**: Located in NBD2; potential inter-NBD dimer destabilization.
- **Known UniProt Clinical Variations**:
  - `VAR ID: VAR_000266`: N -> H (in CF; dbSNP:rs121909042)
  - `VAR ID: VAR_000267`: N -> K (in CF; impaired maturation of glycan chains; has low in vitro channel activity at low temperature; dbSNP:rs80034486)

---

## 4. Integration Blocks for Public 3D Viewer (`/cf/cftr/`)

The following JSON block provides the calibrated camera target, distance metrics, and annotation payload for direct insertion into `cftr-research.js` and `/cf/cftr/`:

```json
{
  "version": "2026-10-13",
  "curative_track": "Track A",
  "dossier_count": 11,
  "variants": [
    {
      "residue": 117,
      "wt": "R",
      "domain": "ABC transmembrane type-1 1 (81-365)",
      "theratype": "Class II / Cryptic Pocket Candidate (Transmembrane domain)",
      "displacement_6msm_5uak": 6.14,
      "displacement_6msm_8eiq": 0.45,
      "landmark_distances": {
        "ICL4_coupling": 77.2,
        "ICL1_coupling": 73.19,
        "ATP_site_1": 98.25,
        "ATP_site_2": 96.5,
        "TMD1_pocket": 79.38,
        "TMD2_pocket": 71.52,
        "pore_constriction": 22.95
      }
    },
    {
      "residue": 334,
      "wt": "R",
      "domain": "ABC transmembrane type-1 1 (81-365)",
      "theratype": "Class IV (Reduced Conductance / Pore Architecture)",
      "displacement_6msm_5uak": 0.71,
      "displacement_6msm_8eiq": 0.59,
      "landmark_distances": {
        "ICL4_coupling": 72.38,
        "ICL1_coupling": 66.3,
        "ATP_site_1": 90.41,
        "ATP_site_2": 90.63,
        "TMD1_pocket": 73.61,
        "TMD2_pocket": 66.7,
        "pore_constriction": 11.26
      }
    },
    {
      "residue": 347,
      "wt": "R",
      "domain": "ABC transmembrane type-1 1 (81-365)",
      "theratype": "Class IV (Reduced Conductance / Pore Architecture)",
      "displacement_6msm_5uak": 0.95,
      "displacement_6msm_8eiq": 0.86,
      "landmark_distances": {
        "ICL4_coupling": 53.61,
        "ICL1_coupling": 46.29,
        "ATP_site_1": 70.18,
        "ATP_site_2": 73.0,
        "TMD1_pocket": 54.01,
        "TMD2_pocket": 48.25,
        "pore_constriction": 9.44
      }
    },
    {
      "residue": 508,
      "wt": "F",
      "domain": "ABC transporter 1 (423-646)",
      "theratype": "Class II (Misfolding / Processing - ICL4 interface)",
      "displacement_6msm_5uak": 9.24,
      "displacement_6msm_8eiq": null,
      "landmark_distances": {
        "ICL4_coupling": 9.31,
        "ICL1_coupling": 21.71,
        "ATP_site_1": 29.07,
        "ATP_site_2": 27.02,
        "TMD1_pocket": 22.06,
        "TMD2_pocket": 11.98,
        "pore_constriction": 67.49
      }
    },
    {
      "residue": 551,
      "wt": "G",
      "domain": "ABC transporter 1 (423-646)",
      "theratype": "Class II (NBD1 structural stability)",
      "displacement_6msm_5uak": 9.32,
      "displacement_6msm_8eiq": 1.46,
      "landmark_distances": {
        "ICL4_coupling": 23.58,
        "ICL1_coupling": 30.23,
        "ATP_site_1": 21.74,
        "ATP_site_2": 8.79,
        "TMD1_pocket": 33.67,
        "TMD2_pocket": 25.16,
        "pore_constriction": 76.76
      }
    },
    {
      "residue": 560,
      "wt": "R",
      "domain": "ABC transporter 1 (423-646)",
      "theratype": "Class II (Misfolding / Processing - ICL4 interface)",
      "displacement_6msm_5uak": 8.87,
      "displacement_6msm_8eiq": 2.22,
      "landmark_distances": {
        "ICL4_coupling": 13.18,
        "ICL1_coupling": 24.11,
        "ATP_site_1": 24.27,
        "ATP_site_2": 22.33,
        "TMD1_pocket": 24.04,
        "TMD2_pocket": 16.99,
        "pore_constriction": 72.94
      }
    },
    {
      "residue": 563,
      "wt": "Y",
      "domain": "ABC transporter 1 (423-646)",
      "theratype": "Class II (Misfolding / Processing - ICL4 interface)",
      "displacement_6msm_5uak": 8.79,
      "displacement_6msm_8eiq": 2.47,
      "landmark_distances": {
        "ICL4_coupling": 16.07,
        "ICL1_coupling": 27.71,
        "ATP_site_1": 27.12,
        "ATP_site_2": 25.54,
        "TMD1_pocket": 25.69,
        "TMD2_pocket": 20.49,
        "pore_constriction": 76.46
      }
    },
    {
      "residue": 1066,
      "wt": "R",
      "domain": "ABC transmembrane type-1 2 (859-1155)",
      "theratype": "Class II (TMD2 folding / VX-445 pocket proximity)",
      "displacement_6msm_5uak": 8.28,
      "displacement_6msm_8eiq": 0.71,
      "landmark_distances": {
        "ICL4_coupling": 5.61,
        "ICL1_coupling": 12.23,
        "ATP_site_1": 23.06,
        "ATP_site_2": 32.95,
        "TMD1_pocket": 10.09,
        "TMD2_pocket": 11.73,
        "pore_constriction": 64.94
      }
    },
    {
      "residue": 1070,
      "wt": "R",
      "domain": "ABC transmembrane type-1 2 (859-1155)",
      "theratype": "Class II (TMD2 folding / VX-445 pocket proximity)",
      "displacement_6msm_5uak": 8.2,
      "displacement_6msm_8eiq": 1.19,
      "landmark_distances": {
        "ICL4_coupling": 3.86,
        "ICL1_coupling": 15.81,
        "ATP_site_1": 28.46,
        "ATP_site_2": 33.13,
        "TMD1_pocket": 14.42,
        "TMD2_pocket": 9.19,
        "pore_constriction": 63.82
      }
    },
    {
      "residue": 1282,
      "wt": "W",
      "domain": "ABC transporter 2 (1210-1443)",
      "theratype": "Class II / IV (NBD2 misassembly or gating transmission)",
      "displacement_6msm_5uak": 6.94,
      "displacement_6msm_8eiq": 1.46,
      "landmark_distances": {
        "ICL4_coupling": 44.5,
        "ICL1_coupling": 44.32,
        "ATP_site_1": 39.5,
        "ATP_site_2": 25.03,
        "TMD1_pocket": 53.18,
        "TMD2_pocket": 42.2,
        "pore_constriction": 72.98
      }
    },
    {
      "residue": 1303,
      "wt": "N",
      "domain": "ABC transporter 2 (1210-1443)",
      "theratype": "Class II / IV (NBD2 misassembly or gating transmission)",
      "displacement_6msm_5uak": 11.61,
      "displacement_6msm_8eiq": 0.67,
      "landmark_distances": {
        "ICL4_coupling": 31.94,
        "ICL1_coupling": 26.9,
        "ATP_site_1": 23.37,
        "ATP_site_2": 28.73,
        "TMD1_pocket": 35.46,
        "TMD2_pocket": 31.06,
        "pore_constriction": 65.32
      }
    }
  ]
}
```

---

## 5. Review & Grounding Verification

- **Primary Grounding Source**: RCSB PDB entries `6MSM` (30459277), `5UAK` (28340353), `8EIQ` (36264792), `8EJ1` (36264792); UniProt `P13569` entry version 286.
- **Mathematical Method**: Least-squares rigid superposition using alpha-carbon coordinates of membrane-spanning domains via Horn's quaternion method.
- **Medical Boundary**: In silico prediction tool only. Not clinical theratyping advice. Prescriptions and modulator access must be guided by approved regulatory labels and licensed clinicians.