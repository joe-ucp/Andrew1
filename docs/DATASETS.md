# Acquire the original inputs

Place the six files below together in a read-only input directory outside the repository. Keep their original names. The pipeline reads archives directly or extracts selected contents into its own new working directory; do not hand-edit timestamps, units, field names, source rows or CSV encodings.

| Dataset | Exact file | Version / identity | Access |
|---|---|---|---|
| DCLP3 | `DCLP3 Public Dataset - Release 3 - 2022-08-04.zip` | Release 3, 4 August 2022 | [Jaeb 573](https://public.jaeb.org/dataset/573); registration and data-use terms |
| REPLACE-BG | `Replace-BG Dataset.zip` | Provider archive identified by the SHA-256 below; filename has no release number | [Jaeb 546](https://public.jaeb.org/dataset/546); registration and data-use terms |
| CLVer | `CLVerPublicDataset.zip` | Provider archive identified by SHA-256; includes protocol v5.2, 1 March 2022 | [Jaeb 591](https://public.jaeb.org/dataset/591); registration and data-use terms |
| TN10 metadata | `TN10_metadata.csv` | Zenodo record 17063348, version 1, 5 September 2025 | [Zenodo](https://zenodo.org/records/17063348); CC BY-SA 4.0 |
| TN10 expression | `TN10_raw_counts.csv` | Same fixed record/version | [Zenodo](https://zenodo.org/records/17063348); CC BY-SA 4.0 |
| AbATE signature | `JCI176403_Table6.pdf` | Supplementary Table 6 for JCI e176403 | [JCI Table 6](https://www.jci.org/articles/view/176403/sd/6); CC BY 4.0. Save under this filename. |

The provider, file bytes and checksum jointly specify the inputs. The CLVer protocol date is a document version, not an invented release date for the archive. If a provider supplies different bytes, retain them separately and investigate the difference before attempting this frozen reproduction; do not update the expected hashes to force acceptance.

| File | Bytes | SHA-256 |
|---|---:|---|
| DCLP3 ZIP | 195,304,936 | `6dbe8580fe38612622f0ebf5da58fed04b4ed556d1e06bf3a67813791222f5dc` |
| REPLACE-BG ZIP | 174,942,355 | `6d4bb4e47e71fff72cd3eb94c5a5e60c9f32e066606bb3883821debf528dc004` |
| CLVer ZIP | 86,726,379 | `8421a9a9ff7dd1485ce0e4bdb984b5ea037300a0e7216f9f445334e56d32d385` |
| TN10 metadata | 44,999 | `dabd7beff35f82aef231f9c9e0fd5b91974cdbc9f175edd8831c475824517586` |
| TN10 counts | 8,508,942 | `74085ba79ef79c44d5e25679e36f3615562c0a644700f2332426d71b942e1ee8` |
| JCI Table 6 | 210,945 | `93b00f6ccb5688eea28f936cfc8e72117e65f70b612c656531cea11c5b51f90c` |

Validate them before computation:

```text
python -I -B src/release.py inputs --inputs-root /path/to/t1d-inputs
```

The [machine-readable manifest](../provenance/source_data.json) contains the same identities. No registration credentials, authorization tokens or participant-level records are included in the release. Download permission is governed by the providers; an input checksum is an identity check, not a grant of permission.

## From originals to required intermediates

The [full pipeline](REPRODUCING.md#reproduce-from-original-datasets) regenerates every dependency required by `src/reconstruct_paper.py`:

1. `card001_data.py` and `card001_replacebg.py` read the original study tables, reconstruct source clocks and participant windows, and write canonical grids with missingness preserved.
2. `card001_run.py` trains the committed 63-fit development grid and selects operating points using development/validation data only; subsequent replay writes event, alarm and participant ledgers for both test populations.
3. `card003_data.py` reconstructs assays, clinical variables and forward CGM intervals from CLVer. The frozen DCLP3 policy is transferred by `card003_policy.py`; clinical replay and estimation produce the interval/event/alarm and fitted-model dependencies.
4. `card002_run.py` extracts the original gene signature from the PDF and rebuilds the 73-record analysis, primary contrast and 10,000 bootstrap draws from original TN10 CSVs. Revision code reconstructs curves, continuous interactions and follow-up bounds.
5. Warning diagnostics build the availability audit; paper diagnostics combine the regenerated ledgers and model summaries. Required paths and generating commands are recorded in [external_dependencies.json](../provenance/external_dependencies.json).

Original archives remain unchanged. Newly produced participant-level intermediates are for the authorized reviewer's local analysis and are excluded from the repository.
