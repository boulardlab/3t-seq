# Changelog

## [1.2.0](https://github.com/boulardlab/3t-seq/compare/v1.1.0...v1.2.0) (2024-05-16)


### Features

* add rulegraph ([804bb1c](https://github.com/boulardlab/3t-seq/commit/804bb1cad06341c13b962b995688fca2bae0c854))
* add workflow documentation.md ([3fb9713](https://github.com/boulardlab/3t-seq/commit/3fb9713f8f97f951bec56a80e2d76c33a9d0d125))
* allow to download only RMSK annotation only for selected chromosomes ([92abc0a](https://github.com/boulardlab/3t-seq/commit/92abc0ae41810ab2ed094693fb1f94b3a8fd3db7))
* bump Docker environment version ([f899fd3](https://github.com/boulardlab/3t-seq/commit/f899fd3c2a656a2b4f1dc8e610bd27c1c52b6cb1))
* init implementation of unit tests - currently not all work ([d9eb2e3](https://github.com/boulardlab/3t-seq/commit/d9eb2e39649c6ac2cb264818425e506ae0c289d3))
* port to Snakemake 8 ([9a0a8af](https://github.com/boulardlab/3t-seq/commit/9a0a8af9a8d7c89787c29377c80247cc007b2336))
* refactor common functions. Now pipeline requires sample sheet to work ([d723b12](https://github.com/boulardlab/3t-seq/commit/d723b1231856fb59539b493e824b90485d1960fa))
* rename and cleanup integration tests folder to comply with WorkflowHub folder tree ([2e908ca](https://github.com/boulardlab/3t-seq/commit/2e908cae29debd028b5d8d57fbb6890c9c2f1ddf))
* updated Dockerfile ([85f4b02](https://github.com/boulardlab/3t-seq/commit/85f4b026cba90380d9f1f2addca5d9bd1077c8ce))


### Bug Fixes

* add resources/sample-sheet.csv and update integration tests config files accordingly ([d61b396](https://github.com/boulardlab/3t-seq/commit/d61b396d739ca46a20f457419e932e414b2e8d09))
* build_trna_matrix, sort columns according to the sample sheet order ([e97ee2e](https://github.com/boulardlab/3t-seq/commit/e97ee2ee77dbc7cbe7f44097bdb73e9786311a39))
* clean .gitattributes ([8637d9c](https://github.com/boulardlab/3t-seq/commit/8637d9c37abf5cc4d957b97a931c1f8d5c94d9d8))
* clean deseq2 rules ([bc33426](https://github.com/boulardlab/3t-seq/commit/bc33426915bb0897734b7a13898f05f1dfcbeb4c))
* cleanup Snakefile ([8c0f6df](https://github.com/boulardlab/3t-seq/commit/8c0f6df31cd2a65ed0da44bad6dc5cc3ed532303))
* coverage_tRNA now makes tmp folder in the correct folder ([9b22e04](https://github.com/boulardlab/3t-seq/commit/9b22e0446aef978d5d51b8883c82aba16c1d77c1))
* fix fastqc_markdup input function and command ([919d088](https://github.com/boulardlab/3t-seq/commit/919d088816faec4dbcff775b44d58145d1e5ec48))
* fix logging in download RMSK script ([dcf7f29](https://github.com/boulardlab/3t-seq/commit/dcf7f29c3718fa1ffd07faa9d9f3a0ad320c098b))
* make filename regex more general ([01eb56f](https://github.com/boulardlab/3t-seq/commit/01eb56fb3252257f832c9d2bbd06790ca55e6769))
* port all rules to new functions ([8488e66](https://github.com/boulardlab/3t-seq/commit/8488e668b73d5b1363d1ca7c690b8792a685f89f))
* protected file path as strings, remove comments ([157d98e](https://github.com/boulardlab/3t-seq/commit/157d98edd6d0654669df64b38ddf5c6f221e50fd))
* remove commented line ([57761e1](https://github.com/boulardlab/3t-seq/commit/57761e1508a83e5a9411ad534bb11e9fdc49fc5f))
* remove shadow directive. As of Snakemake 8.10.0 this is not working in conjunction with Apptainer and Slurm ([602b112](https://github.com/boulardlab/3t-seq/commit/602b112fa0ce439b2abc3566a2a89416f0ce0a94))
* remove unneeded files ([2afaea4](https://github.com/boulardlab/3t-seq/commit/2afaea4be65799788035e72e8d9f7a982b259a4d))
* SalmonTE, polish shell code ([13e16eb](https://github.com/boulardlab/3t-seq/commit/13e16eb54266940b1385a4ad72f7dbb337eff767))
* sorting of tRNA count matrix columns ([ac1e1e4](https://github.com/boulardlab/3t-seq/commit/ac1e1e4bcaedb9a789e599103418403364f4f4b3))
* SSL error on some platforms ([198cbd0](https://github.com/boulardlab/3t-seq/commit/198cbd013d63ab3873eac837fa1add6322f14351))
* STAR, properly handle log files ([ba0c2c7](https://github.com/boulardlab/3t-seq/commit/ba0c2c748333fd033c466140fc1d9339d56ddcb6))
* STARTE, properly hande log files ([fde3c36](https://github.com/boulardlab/3t-seq/commit/fde3c3636b3083f975235f9806f5c069b79e0632))
* trim rules now have stats files as logs, removed unparsable log files from output ([de2bed0](https://github.com/boulardlab/3t-seq/commit/de2bed08518ce4f3a475372b1c52e0eee41c771b))
* update GH Actions workflow to match new test folder tree, add unit tests ([0e2a836](https://github.com/boulardlab/3t-seq/commit/0e2a8365a0d4248420808d47d7dbdf927fa485e1))
* update Snakefile to align with new config ([efa5b57](https://github.com/boulardlab/3t-seq/commit/efa5b57acccd1feffcbc9be3aec15e7aea9ed55b))
* write lfc-shrinked table to the correct file ([e5e3c60](https://github.com/boulardlab/3t-seq/commit/e5e3c60913630f69c7583e268452d1fb254a5e7c))
* YTE template now is a param instead of an input file ([2aab702](https://github.com/boulardlab/3t-seq/commit/2aab7023c353a10ca7fbf8abdc3f52f2a6712095))

## [1.1.0](https://github.com/boulardlab/3t-seq/compare/v1.0.1...v1.1.0) (2024-03-12)


### Features

* containerize pipeline's output ([a3b4086](https://github.com/boulardlab/3t-seq/commit/a3b4086e50eb7e03f81a0f25ef83a00734c60346))
* enable caching for all references and STAR indexes ([bc99687](https://github.com/boulardlab/3t-seq/commit/bc996872682a72aa9228f01b13fd3971eb054a4e))
* test Snakemake report generation ([c2c0f46](https://github.com/boulardlab/3t-seq/commit/c2c0f464bd06e1ca46d741aa127ac9baf5db5e9a))
* update Dockerfile ([bb47f34](https://github.com/boulardlab/3t-seq/commit/bb47f3467ac9b9ef4008c2ffccdf17cc6a8f48d2))


### Bug Fixes

* add singularity-args and reformat lines in CI ([598de8d](https://github.com/boulardlab/3t-seq/commit/598de8d5011170d031c3a01976147f635f9728d6))
* cleanup Snakemake test profile ([2c2cdf7](https://github.com/boulardlab/3t-seq/commit/2c2cdf7400b66d7450e3479beb2d6591bb5bc406))
* formatting edit_condition_file.py ([e26f101](https://github.com/boulardlab/3t-seq/commit/e26f10193d4003ca94f448687ee180e603e5c45d))
* properly download using UCSC API. Fix all paths and enable caching ([db8da50](https://github.com/boulardlab/3t-seq/commit/db8da5075de206edf14bd603d78d152a89ca40a0))
* switch to wget and write some logging when downloading genome annotations ([c49e4ff](https://github.com/boulardlab/3t-seq/commit/c49e4ff1e13cccd039315bd6700af56ce641042c))


### Miscellaneous Chores

* release 1.1.0 ([0340bbc](https://github.com/boulardlab/3t-seq/commit/0340bbc99ecc270f434f7b545118c5b31e08abd7))
* release 1.1.0 ([74e2339](https://github.com/boulardlab/3t-seq/commit/74e23395390c4df2de2f4c3f04395e613b63383b))

## [1.0.1](https://github.com/boulardlab/3t-seq/compare/v1.0.0...v1.0.1) (2024-03-06)


### Features

* add example Slurm configuration ([96ce038](https://github.com/boulardlab/3t-seq/commit/96ce0383c48844a7691342dca721a9a1bb64440f))
* decouple TE, tRNA and single copy genes pipelines. Allow users to decide which branches to run.[#3](https://github.com/boulardlab/3t-seq/issues/3) ([d3b12e4](https://github.com/boulardlab/3t-seq/commit/d3b12e4017c69df6a91c929088cd6641c8790bbb))
* describe how to run tests ([46e023f](https://github.com/boulardlab/3t-seq/commit/46e023f37e9cfbd96b669695027587fa4888d1f4))
* DESeq2 reference is now separate for each sequencing library ([b4cf102](https://github.com/boulardlab/3t-seq/commit/b4cf10264ffd7bc33ad01ef3e4ae1d12482e3718))


### Bug Fixes

* allow deseq2 ananlysi with different variables and different levels ([76c2562](https://github.com/boulardlab/3t-seq/commit/76c2562501c091b615640a0ec3973a814ec6cca3))
* fix deseq2 scripts to look for correct column names ([ffc4537](https://github.com/boulardlab/3t-seq/commit/ffc4537558e28ba5ee278975d63fb7958d8f3384))
* fix problem with SalmonTE quant input array on single-end libraries ([3fd51ff](https://github.com/boulardlab/3t-seq/commit/3fd51fff87474b1243f7844d4699eb34e445e352))
* fix Snakemake version to be below 8 ([5fa9be3](https://github.com/boulardlab/3t-seq/commit/5fa9be31940270c93024a8cf6a5190215a58ea7d))
* look for correct column name in edit_conditions_file.py ([e60456e](https://github.com/boulardlab/3t-seq/commit/e60456ed4c3e471b21417c2dc4126c46b6d91382))
* make test dataset smaller to let it run in GH Actions ([17ac70d](https://github.com/boulardlab/3t-seq/commit/17ac70d120fd772f369c76d9514e260b2a31746d))
* pin all conda packages versions ([feddd48](https://github.com/boulardlab/3t-seq/commit/feddd48cd06365923fbee22ffc5385a6ec7bebfb))
* remove params.mem_mb from starTE rules ([b6f1848](https://github.com/boulardlab/3t-seq/commit/b6f18484851a60224d1520166b68ed0df9c2c5f0))
* restore test profile ([c50b75e](https://github.com/boulardlab/3t-seq/commit/c50b75ebb4e27459c58e298b1546729b651efc83))
* revert alpine base for SalmonTE Docker ([bfb4a48](https://github.com/boulardlab/3t-seq/commit/bfb4a48f156494a8605dd385c34da6b468c60d1c))
* switch to miniconda3-alpine as base ([7ac38c7](https://github.com/boulardlab/3t-seq/commit/7ac38c79ed9f88fb2a0d93d0cae99302cc752563))
* test config now has library-specific deseq2 reference level ([44d4b43](https://github.com/boulardlab/3t-seq/commit/44d4b43565939785370fa4b1233ec025c1c08614))
* typo ([6e64c84](https://github.com/boulardlab/3t-seq/commit/6e64c848bd9b1765329bac1d8a0a075b98546878))
* update docs ([c811991](https://github.com/boulardlab/3t-seq/commit/c811991b8186ae6aba348964fbbb1c01096292b8))


### Miscellaneous Chores

* release 1.0.1 ([25c11bc](https://github.com/boulardlab/3t-seq/commit/25c11bc0886cf580768fd61ac373ceaca97b997c))

## 1.0.0 (2023-11-17)


### Features

* add file for Snakemake workflow catalog ([cf8723d](https://github.com/boulardlab/3t-seq/commit/cf8723d5883490917ded763e23468eef14529598))
* add pseudocounts to count matrix if cannot estimate size factors ([219d665](https://github.com/boulardlab/3t-seq/commit/219d66534c0b2c7ba8c7b9440d7de7c1d5d7774c))
* explicit all star outputs ([2f2d7b9](https://github.com/boulardlab/3t-seq/commit/2f2d7b9d747fd3e417cc6b70b9e8c8530c4a619b))
* start basic documentation ([1d57ff5](https://github.com/boulardlab/3t-seq/commit/1d57ff57d3e081914db188ac4ea5ca8ca2165fd9))
* update Dockerfile ([877bf88](https://github.com/boulardlab/3t-seq/commit/877bf8846652c89035864a97cc83945bc0217a19))


### Bug Fixes

* add tests ([418381b](https://github.com/boulardlab/3t-seq/commit/418381be91cbba78942f953a11a23d4c819db230))
* add text data ([328781a](https://github.com/boulardlab/3t-seq/commit/328781a3ee4721c6501968f98cf05f703fe70bf5))
* config directory ([bbf4926](https://github.com/boulardlab/3t-seq/commit/bbf492680d9cb1807719818ae2c126506b755715))
* correctly pass genome label to get_rmsk.py ([8157a11](https://github.com/boulardlab/3t-seq/commit/8157a11836f56526c70ef94511ae50ff658196b2))
* edit_condition_file rule ([adb0a5e](https://github.com/boulardlab/3t-seq/commit/adb0a5ec36bd4fdab975d00fba255a7ecf3505de))
* mkdir ([646a3d4](https://github.com/boulardlab/3t-seq/commit/646a3d4fa781622f3abdab9baaf862ab8baddb32))
* mv configuration docs to specific README file ([a8adec4](https://github.com/boulardlab/3t-seq/commit/a8adec4e0c5e1efdfe8cc68cfc326bae7b48c02d))
* pin R version for SalmonTE ([b7c238b](https://github.com/boulardlab/3t-seq/commit/b7c238baafd2510295944cfebe21d6af76a93b17))
* polish ([5d4b3b7](https://github.com/boulardlab/3t-seq/commit/5d4b3b70e62eeec3f9dde11f203349f669773b6d))
* remove useless print ([fa11e58](https://github.com/boulardlab/3t-seq/commit/fa11e58041ab6427b6aea8737e9a3abedd6408ad))
* snakemake variable ([a123c63](https://github.com/boulardlab/3t-seq/commit/a123c63128c31f11cdc41093e1d6adb14dc8829e))
* syntax ([1d20176](https://github.com/boulardlab/3t-seq/commit/1d201761d5f4489ecee120f718d1f1a94460f948))
* update TOKEN in main CI workflow ([3d9faaa](https://github.com/boulardlab/3t-seq/commit/3d9faaa92e4e51a1059d18b73fcba03fdc8a12f1))
