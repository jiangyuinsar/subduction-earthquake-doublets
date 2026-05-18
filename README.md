# subduction-earthquake-doublets

MATLAB code and catalog file for identifying earthquake doublets in global subduction zones.

This repository accompanies the manuscript:

> Yu Jiang, Ross S. Stein, and Daniel T. Trugman.  
> **Earthquake doublets in subduction zones are an Omori process promoted by megathrust-tear fault interaction.**  
> *Communications Earth & Environment*.

## Repository contents

```text
subduction-earthquake-doublets/
│
├── identify_subduction_earthquake_doublets.m
├── gCMT_catalog_19760101_20260101_M6+.mat
└── ThirdPartyScript/
    └── llh2local.m
```

## Description

The MATLAB script identifies candidate earthquake doublets from a preprocessed Global CMT earthquake catalog. A doublet is defined as a pair of earthquakes that satisfy specified criteria for:

1. magnitude similarity,
2. spatial proximity,
3. temporal separation, and
4. exclusion of pairs where both events are likely aftershocks of larger nearby preceding earthquakes.

The script tests multiple parameter combinations, including different magnitude-difference thresholds, distance thresholds, time windows, and 2-D versus 3-D distance calculations.

## Input data

The included MATLAB catalog file is:

```text
gCMT_catalog_19760101_20260101_M6+.mat
```

It contains the variable:

```matlab
earthquake_catalog
```

with columns:

```text
[longitude latitude depth Mw year month day type]
```

where `type` indicates the event mechanism or classification used in the analysis.

## How to run

Open MATLAB, move to the repository folder, and run:

```matlab
identify_subduction_earthquake_doublets
```

The script automatically adds the `ThirdPartyScript` folder to the MATLAB path.

## Output

Detected doublet lists are saved in:

```text
txt_Doublets/
```

Each output file corresponds to one combination of magnitude, distance, time-window, and 2-D/3-D distance criteria.

Each row of the output file contains:

```text
[event_index_1 event_index_2 magnitude_difference time_difference_days distance_km mean_magnitude type_1 type_2]
```

## Required MATLAB function

The script uses:

```matlab
llh2local.m
```

which is included in the `ThirdPartyScript` folder.

## Citation

If you use this code, please cite:

> Jiang, Y., Stein, R. S., and Trugman, D. T.  
> **Earthquake doublets in subduction zones are an Omori process promoted by megathrust-tear fault interaction.**  
> *Communications Earth & Environment*.

After publication, please cite the final DOI of the article and, if available, the archived DOI of this GitHub repository.

## License

This repository is released under the MIT License. See the `LICENSE` file for details.
