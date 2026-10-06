# DNA and variant MCP servers, APIs and datasets for CFTR work: a survey

Status: proposed, written 2026-10-05. Everything here is **unverified**. A research agent read the pages with a fetch tool that passes each page through a small summarising model, so quotations are relayed, not byte-checked, and some pages did not render or refused the request. Nothing was installed, cloned or run. Nothing here changes `sources/catalog.yaml`: a source enters the catalog only after a person reads its terms and robots.txt and the lab's own tools do the reading.

## Conflicts with what the catalog already records

| Tool | What it does | What the catalog says |
| --- | --- | --- |
| BioMCP, OpenCRAVAT (ClinGen lookups) | Reach ClinGen services on the user's behalf | `clingen`: robots.txt refuses all non-search-engine agents, so access is `forbidden`. The Allele Registry is a separate host (`reg.clinicalgenome.org`) whose robots.txt was not read |
| BioMCP, GenomeMCP, a third-party gnomAD MCP | Call the gnomAD GraphQL API | `gnomad`: robots.txt disallows API paths and no API terms were found, so access is `manual` |

So an MCP wrapper is not a way around a catalog entry. Admitting one means its calls must stay inside what the catalog allows for each source, and that is the board's T-0042 test.

## Candidates, by use

| Source | Kind | What it gives for CFTR | Terms as relayed | Read |
| --- | --- | --- | --- | --- |
| Ensembl REST (VEP, variant recoder, lookup) | API | Consequence of an HGVS variant on transcript and protein | Data "no restrictions"; a published limit of 55,000 requests an hour per IP; predictor plugins such as SpliceAI are non-commercial | partly (limit and disclaimer pages) |
| NCBI Variation Services, Datasets, dbSNP | API | HGVS or rsID to SPDI; reference records | NCBI places no restrictions; E-utilities at most 3 requests a second | partly (documentation did not render) |
| UniProt REST | API | CFTR (P13569) domains, topology and variant features | CC BY 4.0, from a mirror's licence text, not uniprot.org | partly |
| Open Targets Platform | API, dataset | Target-disease evidence, drugs and GWAS links; not a variant classifier | Data CC0 1.0, code Apache-2.0 | read |
| ClinVar bulk files | dataset | Offline classification set, no API call per query | Attribution requested; no explicit open licence on the page | read |
| UCSC Genome Browser API | API | Coordinates and tracks; not interpretation | Its two pages give different rate limits, so use the stricter: one hit per 15 seconds, 5,000 a day | read |
| NHGRI-EBI GWAS Catalog | API | Modifier-gene associations; not variant classification | EMBL-EBI terms; summary statistics CC0 | read |
| ClinPGx (PharmGKB) | API | Modulator pharmacogenomic annotations | CC BY-SA 4.0 (a share-alike constraint); official page returned nothing | partly |
| Monarch Initiative v3 | API | Gene-disease and phenotype links | CC BY per a snippet | partly |
| BioMCP | MCP | One server over dozens of public APIs including ClinVar, gnomAD, ClinGen, UniProt | MIT; upstream terms govern results; curl-pipe-shell install, so use PyPI or a pinned binary; high release churn | read |
| Open Targets MCP, BioContextAI knowledge base MCP | MCP | Wrappers over APIs above; BioContextAI adds Ensembl and InterPro | Apache-2.0; both say experimental or fair-use | read |
| UniProt MCP (individual maintainer) | MCP | CFTR protein entry | MIT; two stars | partly |
| Anthropic life-sciences marketplace | MCP | Literature and trials only, no variant data; several paid | Per vendor | partly |

Not recommended, with the reason given: OMIM (key tied to an account, research-only terms, the agreement page refused the read), CFTR1 at SickKids (certificate mismatch, terms unknown; do not bypass), the shared LOVD CFTR database (already `forbidden`), a hosted OpenCRAVAT endpoint (sends variant queries to a third party), a gnomAD MCP with no licence and few commits, GenomeMCP's optional Supabase persistence (sends data to a third party), and BioContextAI's Google Scholar and KEGG tools. GA4GH Beacon and refget were named but not researched. PanelApp's terms are a PDF that could not be read (the fetch tool saved it locally; a person can open it).

## Suggested order

1. Ensembl REST for the consequence of an HGVS variant: official source, no key, a stated limit the lab would cap far lower. Read its terms and robots.txt through the lab's own reader first.
2. ClinVar bulk `variant_summary` as an offline option beside the existing E-utilities provider.
3. UniProt REST for protein domains, after one human read of uniprot.org's licence page.
4. NCBI Variation Services for normalisation, after its documentation is read.
5. Wrappers last, and only BioMCP in a pinned release as a discovery aid under T-0042, with its ClinGen and gnomAD calls blocked until the catalog entries are revisited.

Board follow-ups (proposed): T-0070 (read the terms and robots.txt of the four sources in the first four steps with the lab's reader, and propose catalog entries), T-0071 (add the conflicts above to T-0042's acceptance test).

## Unknown

Whether BioMCP honours robots.txt or the lab's network rules; Ensembl's robots.txt; NCBI Variation and Datasets limits; the ClinGen Allele Registry's terms; uniprot.org's own licence text; the official terms of ClinPGx, Monarch and PanelApp; any published API terms for gnomAD; the current status of CFTR1; the OpenCRAVAT licence; and how the third-party MCP servers pin dependencies or run code at install.
