# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][],
and this project adheres to [Semantic Versioning][].

[keep a changelog]: https://keepachangelog.com/en/1.0.0/
[semantic versioning]: https://semver.org/spec/v2.0.0.html

## [Unreleased]

### Added

- Requests-cache dependency for robust dataset downloads
- Documentation improvements and tutorial updates

### Fixed

- Parameter name in enrichment analysis function
- Typo in requests-cache dependency
- MS CSF dataset download reliability

## [0.1.2] -- 2025-12-22

### Added

- Lupus tutorial with data pre-processing
- Support for Cytokine-induced Immune Programs (CIP) signatures as input
- Support for custom gene-program signatures
- CIP gene signature data
- Well-bias filter in `load_human_cytokine_dict()`
- Manuscript plots and README abstract

### Fixed

- Indexing bug in `_get_senders()`
- Docstrings for `_get_genesets()`

## [0.1.1] -- 2025-12-17

### Fixed

- Bug in `get_senders()` cytokine communication function
- Updated tutorials

## [0.1.0] -- 2025-12-16

### Added

- Cytokine communication analysis (`get_one_senders_and_receivers`, `get_all_senders_and_receivers`)
- Circos plot for cell-cell communication (`plot_communication`)
- Data loading functions for MS CSF, Lupus, and cytokine info datasets
- Citation for the cytokine dictionary article

### Changed

- Restructured enrichment test and dataset modules

## [0.0.1] -- 2025-12-07

### Added

- Initial release
- Core enrichment analysis (`run_one_enrichment_test`, `run_all_enrichment_test`)
- Robustness testing (`check_robustness`, `get_robust_significant_results`)
- Heatmap plotting (`plot_significant_results`)
- Human Cytokine Dictionary data loading
- Multiple Sclerosis tutorial notebooks
