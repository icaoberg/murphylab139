# murphylab139

> **Deprecated:** This repository is no longer maintained. It is archived here for reproducibility purposes only.

## About

This repository contains the MATLAB code and data scripts for the work described in:

**Shann-Ching Chen, Ting Zhao, Geoffrey J. Gordon, and Robert F. Murphy.**
*Automated Image Analysis of Protein Localization in Budding Yeast.*
Bioinformatics 23:i66–i71, 2007.

[https://doi.org/10.1093/bioinformatics/btm207](https://doi.org/10.1093/bioinformatics/btm207)

## Summary

This work presents an automated pipeline for classifying subcellular protein localization patterns in budding yeast (*S. cerevisiae*) from fluorescence microscopy images. The system:

- Segments yeast cells from GFP images
- Extracts cell-level and field-level features
- Classifies protein localization using support vector machines (SVM)
- Evaluates results against UCSF and CYGD localization annotations

## Dependencies

- MATLAB
- [LIBSVM](http://www.csie.ntu.edu.tw/~cjlin/libsvm/) (included under `./SLIC/classify/libsvm`)
- [SLIC (Subcellular Localization Image Classifier)](http://pslid.cbi.cmu.edu/release/) (included under `./SLIC`)

## Data

Yeast GFP images are from the UCSF Yeast GFP Fusion Localization Database:
[http://yeastgfp.ucsf.edu/](http://yeastgfp.ucsf.edu/)

Use `download_data_for_paper.sh` to download the image data.

## Usage

### Single machine

Run `makedata.m` to execute the full pipeline and reproduce all results.

### Cluster

Run scripts in order:

| Step | Script | Description |
|------|--------|-------------|
| 0 | `makedata0.m` | Generate `yeastdata.mat` and `procfiles.mat` |
| 1 | `makedata1_clusterv1.sh` | Image segmentation |
| 2 | `makedata2_clusterv1.sh` | Cell-level feature calculation |
| 3 | `makedata3_clusterv1.sh` | Field-level feature calculation |
| 4 | `makedata4_clusterv1.sh` | Load features for cell-level classification |
| 5 | See below | Cell-level classification (25 seeds) |
| 6 | `makedata6_clusterv1.sh` | Load features for field-level classification |
| 7 | `makedata7_clusterv1.sh` | Field-level classification |
| 8 | `makedata8.m` | Generate tables and figures |

Cell-level classification scripts for step 5 must be run in this order:

1. `makedata5_createnfoldscriptFlag3.sh`
2. `makedata5_classifyscriptFlag3R.sh`
3. `makedata5_combinefoldcriptFlag3.sh`
4. `makedata5_combineseedcriptFlag3.sh`
5. `makedata5_combinefoldcriptFlag41.sh`
6. `makedata5_combinefoldcriptFlag42.sh`
7. `makedata5_combineseedcriptFlag41.sh`
8. `makedata5_combineseedcriptFlag42.sh`

## License

GNU General Public License v2 or later. See [LICENSE](LICENSE).

## Contact

Murphy Lab, Carnegie Mellon University
[http://murphylab.web.cmu.edu](http://murphylab.web.cmu.edu)
