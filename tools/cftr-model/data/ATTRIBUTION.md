# Attribution for CFTR Structure and Sequence Data

This directory holds the five canonical structure and sequence data files required for offline verification and execution of `tools/cftr-model` (resolving task **T-0132**).

All five files are verified byte-for-byte against the cryptographic SHA-256 digests recorded in [`sources/manifest.json`](../../../sources/manifest.json), and every structure entry below records the exact `_struct.title`, deposition date, resolution, and primary citation fields (`_citation.title`, `_citation.journal_abbrev`, `_citation.year`, `_citation.journal_volume`, `_citation.page_first`, `_citation.page_last`, `_citation.pdbx_database_id_DOI`, `_citation.pdbx_database_id_PubMed`) stored inside the mmCIF file.

---

## 1. RCSB Protein Data Bank (PDB) Structures (CC0 1.0)

Under the [RCSB PDB Data Usage Policy](https://www.rcsb.org/pages/usage-policy), all structure coordinate files in the PDB archive are dedicated to the public domain under the Creative Commons CC0 1.0 Universal Public Domain Dedication.

### 6MSM (`6MSM.cif.gz`)
- **Title**: Phosphorylated, ATP-bound human cystic fibrosis transmembrane conductance regulator (CFTR)
- **Method**: ELECTRON MICROSCOPY, 3.2 Å resolution
- **Deposited**: 2018-10-16
- **Authors**: Zhang Z, Liu F, Chen J
- **Primary Citation**: Zhang Z, Liu F, Chen J. *Molecular structure of the ATP-bound, phosphorylated human CFTR.* **Proc. Natl. Acad. Sci. U.S.A.** 2018; 115:12757-12762. DOI: [10.1073/pnas.1815287115](https://doi.org/10.1073/pnas.1815287115). PubMed: [30459277](https://pubmed.ncbi.nlm.nih.gov/30459277/).
- **License**: Creative Commons CC0 1.0 Universal Public Domain Dedication
- **SHA-256**: `82933988abc3ed24a343819d410b39dc0d86913eaba4426dd945e3273e34fa59`

### 5UAK (`5UAK.cif.gz`)
- **Title**: Dephosphorylated, ATP-free human cystic fibrosis transmembrane conductance regulator (CFTR)
- **Method**: ELECTRON MICROSCOPY, 3.87 Å resolution
- **Deposited**: 2016-12-19
- **Authors**: Liu F, Zhang Z, Csanády L, Gadsby DC, Chen J
- **Primary Citation**: Liu F, Zhang Z, Csanády L, Gadsby DC, Chen J. *Molecular Structure of the Human CFTR Ion Channel.* **Cell** 2017; 169:85-95.e8. DOI: [10.1016/j.cell.2017.02.024](https://doi.org/10.1016/j.cell.2017.02.024). PubMed: [28340353](https://pubmed.ncbi.nlm.nih.gov/28340353/).
- **License**: Creative Commons CC0 1.0 Universal Public Domain Dedication
- **SHA-256**: `8bf00f2612dd961d3b1ccd2f2193781a534ad44f32666a6de8da948afb8d84cf`

### 8EIQ (`8EIQ.cif.gz`)
- **Title**: The complex of phosphorylated human delta F508 cystic fibrosis transmembrane conductance regulator (CFTR) with Trikafta [elexacaftor (VX-445), tezacaftor (VX-661), ivacaftor (VX-770)] and ATP/Mg
- **Method**: ELECTRON MICROSCOPY, 3 Å resolution
- **Deposited**: 2022-09-15
- **Authors**: Fiedorczuk K, Chen J
- **Primary Citation**: Fiedorczuk K, Chen J. *Molecular structures reveal synergistic rescue of Delta 508 CFTR by Trikafta modulators.* **Science** 2022; 378:284-290. DOI: [10.1126/science.ade2216](https://doi.org/10.1126/science.ade2216). PubMed: [36264792](https://pubmed.ncbi.nlm.nih.gov/36264792/).
- **License**: Creative Commons CC0 1.0 Universal Public Domain Dedication
- **SHA-256**: `d13cbfc8ba6575ab6228f351578c6f9cc5b9dafe13e5913d9b27b3cacebc18d6`

### 8EJ1 (`8EJ1.cif.gz`)
- **Title**: Dephosphorylated human delta F508 cystic fibrosis transmembrane conductance regulator (CFTR)
- **Method**: ELECTRON MICROSCOPY, 6.9 Å resolution
- **Deposited**: 2022-09-16
- **Authors**: Fiedorczuk K, Chen J
- **Primary Citation**: Fiedorczuk K, Chen J. *Molecular structures reveal synergistic rescue of Delta 508 CFTR by Trikafta modulators.* **Science** 2022; 378:284-290. DOI: [10.1126/science.ade2216](https://doi.org/10.1126/science.ade2216). PubMed: [36264792](https://pubmed.ncbi.nlm.nih.gov/36264792/).
- **License**: Creative Commons CC0 1.0 Universal Public Domain Dedication
- **SHA-256**: `9582203b74c0add00a7917bed765cbc8a943fcc3a88ffbac8ae8dde7de2bc159`

---

## 2. UniProtKB Record (CC BY 4.0)

### Human CFTR (`p13569.json`)
- **Entry Accession**: `P13569`
- **Entry Name**: `CFTR_HUMAN`
- **Protein**: Cystic fibrosis transmembrane conductance regulator (ATP-binding cassette sub-family C member 7)
- **Gene**: *CFTR* (Chromosome 7q31.2)
- **Source**: UniProt Consortium (European Bioinformatics Institute, Swiss Institute of Bioinformatics, Protein Information Resource)
- **URL**: `https://rest.uniprot.org/uniprotkb/P13569.json`
- **License**: Creative Commons Attribution 4.0 International License ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/))
- **Attribution Statement**: CFTR sequence, topological domain boundaries, transmembrane helical assignments, and natural variant disease annotations are provided by the UniProt Consortium under CC BY 4.0.
- **SHA-256**: `c8feabfabe8690ef4a06af7de5b5f8dc3c8da554246c5ab3d403edf265d4a893`
