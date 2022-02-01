# Benchmark Results

<!-- tip: use https://tabletomarkdown.com/convert-website-table-to-markdown/ -->

## Cotter

### write flags - ord_half_1247842824 none

Cotter, obsid 1247842824 (batch 00) when run on the same machine `scripts/time_cotter.sh` multiple times produces the following results.

| reading         | processing      | writing         |
| --------------- | --------------- | --------------- |
| 00:02:26.659619 | 00:00:56.169862 | 00:00:49.997268 |
| 00:02:26.579678 | 00:00:56.697208 | 00:00:50.511478 |
| 00:02:26.591573 | 00:00:55.746477 | 00:00:51.935500 |
| 00:02:26.314757 | 00:00:53.764417 | 00:00:55.843333 |
| 00:02:26.571668 | 00:00:54.717671 | 00:00:57.614381 |
| 00:02:26.644677 | 00:00:53.283708 | 00:00:59.344707 |
| 00:02:27.158254 | 00:00:57.250386 | 00:00:51.650592 |
| 00:02:28.907182 | 00:00:55.070992 | 00:00:52.034414 |
| 00:02:33.169952 | 00:00:56.878490 | 00:00:52.091193 |
| 00:02:31.007132 | 00:00:55.591012 | 00:00:50.448318 |

### write uvfits - ord_half_1247842824 both

| reading         | processing      | writing         |
| --------------- | --------------- | --------------- |
| 00:02:52.452216 | 00:01:21.827627 | 00:08:35.466152 |

### write uvfits - ord_half_1196175296 both

| reading         | processing      | writing         |
| --------------- | --------------- | --------------- |
| 00:05:35.140382 | 00:02:35.595915 | 00:18:47.143991 |

### write uvfits - ord_small_1254670392 both

| reading         | processing      | writing         |
| --------------- | --------------- | --------------- |
| 00:00:02.467415 | 00:00:00.159839 | 00:00:01.953713 |
| 00:00:02.393905 | 00:00:00.143395 | 00:00:02.388053 |

## [commit a762a27 (v0.1.1)](https://github.com/MWATelescope/Birli/actions/runs/812007480)

total duration: 5h 47m 35s

### a762a27 context_to_baseline_imgsets - ord_half_1196175296

![a762a27 context_to_baseline_imgsets - ord_half_1196175296 pdf](img/a762a27%20context_to_baseline_imgsets%20-%20ord_half_1196175296%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0574858   | 0.0809071 | 0.0629698   |
| Mean      | 1868.3 s    | 1873.5 s  | 1877.9 s    |
| Std. Dev. | 3.6753 s    | 8.2157 s  | 11.304 s    |
| Median    | 1869.1 s    | 1875.4 s  | 1879.4 s    |
| MAD       | 1.7995 s    | 7.4133 s  | 12.325 s    |

## [commit ea16381](https://github.com/MWATelescope/Birli/actions/runs/813672907)

performance changes:

- initial optimization of visibility loading by reducing `pin()`s.
- added mwax benchmark

total duration: 3h 31m 2s

### ea16381 context_to_baseline_imgsets - ord_half_1196175296

![ea16381 context_to_baseline_imgsets - ord_half_1196175296 pdf](img/ea16381%20context_to_baseline_imgsets%20-%20ord_half_1196175296%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.1254305   | 0.1673010 | 0.1302936   |
| Mean      | 725.88 s    | 729.49 s  | 732.81 s    |
| Std. Dev. | 2.7240 s    | 6.0055 s  | 6.8294 s    |
| Median    | 721.67 s    | 732.86 s  | 734.09 s    |
| MAD       | 485.85 ms   | 2.9981 s  | 9.2880 s    |

### ea16381 context_to_baseline_imgsets - mwax_half_1247842824

![ea16381 context_to_baseline_imgsets - mwax_half_1247842824 pdf](img/dbdc0ee%20context_to_baseline_imgsets%20-%20mwax_half_1247842824%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0143499   | 0.0194045 | 0.0141804   |
| Mean      | 371.83 s    | 372.64 s  | 373.46 s    |
| Std. Dev. | 892.73 ms   | 1.4113 s  | 1.7007 s    |
| Median    | 371.33 s    | 372.40 s  | 373.89 s    |
| MAD       | 219.26 ms   | 1.9937 s  | 2.2897 s    |

## [commit dbdc0ee](https://github.com/MWATelescope/Birli/actions/runs/815126545)

![dbdc0ee context_to_baseline_imgsets - ord_half_1196175296 pdf](img/dbdc0ee%20context_to_baseline_imgsets%20-%20ord_half_1196175296%20pdf.svg)

performance changes:

- parallelize visibility loading with crossbeam
- added mwax benchmark

total duration: 2h 8m 33s

### dbdc0ee context_to_baseline_imgsets - ord_half_1196175296

![dbdc0ee context_to_baseline_imgsets - ord_half_1196175296 pdf](img/dbdc0ee%20context_to_baseline_imgsets%20-%20ord_half_1196175296%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0275373   | 0.0383631 | 0.0283015   |
| Mean      | 454.55 s    | 458.65 s  | 462.56 s    |
| Std. Dev. | 4.2645 s    | 6.7653 s  | 8.1761 s    |
| Median    | 453.08 s    | 459.36 s  | 464.87 s    |
| MAD       | 1.3541 s    | 8.6210 s  | 11.241 s    |

### dbdc0ee context_to_baseline_imgsets - mwax_half_1247842824

![dbdc0ee context_to_baseline_imgsets - mwax_half_1247842824 pdf](img/dbdc0ee%20context_to_baseline_imgsets%20-%20mwax_half_1247842824%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0126678   | 0.0172310 | 0.0124256   |
| Mean      | 223.02 s    | 226.30 s  | 229.69 s    |
| Std. Dev. | 3.7627 s    | 5.7011 s  | 6.6885 s    |
| Median    | 220.65 s    | 226.49 s  | 231.33 s    |
| MAD       | 942.88 ms   | 7.1007 s  | 9.4088 s    |

## [commit 09abc69](https://github.com/MWATelescope/Birli/actions/runs/831282207)

![09abc69 context_to_baseline_imgsets - ord_half_1196175296 pdf](img/09abc69%20context_to_baseline_imgsets%20-%20ord_half_1196175296%20pdf.svg)

performance changes:

- parallelize writing to img_bufs
- parallelize over polarizations
- use Vec instead of BTreeMap where possible.
- use chunks_exact instead of chunks

total duration: 52m 52s

### 09abc69 context_to_baseline_imgsets - ord_half_1196175296

![09abc69 context_to_baseline_imgsets - ord_half_1196175296 pdf](img/09abc69%20context_to_baseline_imgsets%20-%20ord_half_1196175296%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0102914   | 0.0144746 | 0.0108148   |
| Mean      | 169.76 s    | 172.69 s  | 175.37 s    |
| Std. Dev. | 1.9382 s    | 4.8114 s  | 6.3284 s    |
| Median    | 168.97 s    | 174.25 s  | 175.44 s    |
| MAD       | 297.86 ms   | 3.1589 s  | 7.8800 s    |

### 09abc69 context_to_baseline_imgsets - mwax_half_1247842824

![09abc69 context_to_baseline_imgsets - mwax_half_1247842824 pdf](img/09abc69%20context_to_baseline_imgsets%20-%20mwax_half_1247842824%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0613268   | 0.0798494 | 0.0575464   |
| Mean      | 85.874 s    | 87.029 s  | 88.337 s    |
| Std. Dev. | 1.0761 s    | 2.1247 s  | 2.7667 s    |
| Median    | 85.245 s    | 86.548 s  | 88.561 s    |
| MAD       | 142.20 ms   | 2.0188 s  | 3.7638 s    |

## [commit 3a5671c](https://github.com/MWATelescope/Birli/actions/runs/1157414447)

performance changes:

- use ndarray (only in context_to...)

total duration: 3h 44m

### 3a5671c context_to_baseline_imgsets - ord_half_1196175296

![3a5671c context_to_baseline_imgsets - ord_half_1196175296 pdf](img/3a5671c%20context_to_baseline_imgsets%20-%20ord_half_1196175296%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0078053   | 0.0100363 | 0.0069116   |
| Mean      | 97.148 s    | 98.537 s  | 100.28 s    |
| Std. Dev. | 943.72 ms   | 2.7233 s  | 3.8568 s    |
| Median    | 96.527 s    | 97.738 s  | 99.662 s    |
| MAD       | 238.80 ms   | 2.0151 s  | 3.6442 s    |

### 3a5671c context_to_baseline_imgsets - mwax_half_1247842824

![3a5671c context_to_baseline_imgsets - mwax_half_1247842824 pdf](img/3a5671c%20context_to_baseline_imgsets%20-%20mwax_half_1247842824%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0124288   | 0.0157458 | 0.0106958   |
| Mean      | 50.970 s    | 51.940 s  | 53.230 s    |
| Std. Dev. | 567.02 ms   | 1.9627 s  | 2.8412 s    |
| Median    | 50.695 s    | 51.221 s  | 52.629 s    |
| MAD       | 40.160 ms   | 1.1959 s  | 2.4080 s    |

### 3a5671c context_to_jones_array - ord_half_1196175296

![3a5671c context_to_jones_array - ord_half_1196175296 pdf](img/3a5671c%20context_to_jones_array%20-%20ord_half_1196175296%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0035676   | 0.0044346 | 0.0029756   |
| Mean      | 63.475 s    | 64.717 s  | 66.482 s    |
| Std. Dev. | 735.02 ms   | 2.6507 s  | 3.9595 s    |
| Median    | 63.293 s    | 64.095 s  | 65.239 s    |
| MAD       | 265.59 ms   | 1.4380 s  | 2.6267 s    |

### 3a5671c context_to_jones_array - mwax_half_1247842824

![3a5671c context_to_jones_array - mwax_half_1247842824 pdf](img/3a5671c%20context_to_jones_array%20-%20mwax_half_1247842824%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0305077   | 0.0449094 | 0.0364725   |
| Mean      | 34.915 s    | 35.321 s  | 35.605 s    |
| Std. Dev. | 103.66 ms   | 608.60 ms | 879.39 ms   |
| Median    | 35.118 s    | 35.531 s  | 35.652 s    |
| MAD       | 57.163 ms   | 194.02 ms | 592.96 ms   |

### 3a5671c correct_cable_lengths - ord_half_1196175296

![3a5671c correct_cable_lengths - ord_half_1196175296 pdf](img/3a5671c%20correct_cable_lengths%20-%20ord_half_1196175296%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0045064   | 0.0066382 | 0.0051820   |
| Mean      | 10.450 s    | 10.543 s  | 10.614 s    |
| Std. Dev. | 35.885 ms   | 141.66 ms | 206.94 ms   |
| Median    | 10.508 s    | 10.581 s  | 10.602 s    |
| MAD       | 9.1543 ms   | 54.572 ms | 165.92 ms   |

### 3a5671c correct_cable_lengths - mwax_half_1247842824

![3a5671c correct_cable_lengths - mwax_half_1247842824 pdf](img/3a5671c%20correct_cable_lengths%20-%20mwax_half_1247842824%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0329479   | 0.0433217 | 0.0308323   |
| Mean      | 6.3688 s    | 6.3826 s  | 6.3982 s    |
| Std. Dev. | 10.892 ms   | 25.335 ms | 30.811 ms   |
| Median    | 6.3630 s    | 6.3694 s  | 6.4105 s    |
| MAD       | 2.9548 ms   | 13.848 ms | 42.508 ms   |

### 3a5671c correct_geometry - ord_half_1196175296

![3a5671c correct_geometry - ord_half_1196175296 pdf](img/3a5671c%20correct_geometry%20-%20ord_half_1196175296%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0365140   | 0.0472116 | 0.0326765   |
| Mean      | 646.13 s    | 650.03 s  | 654.84 s    |
| Std. Dev. | 2.7767 s    | 7.4174 s  | 10.225 s    |
| Median    | 644.78 s    | 648.94 s  | 653.43 s    |
| MAD       | 623.61 ms   | 4.7188 s  | 10.869 s    |

### 3a5671c correct_geometry - mwax_half_1247842824

![3a5671c correct_geometry - mwax_half_1247842824 pdf](img/3a5671c%20correct_geometry%20-%20mwax_half_1247842824%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0221183   | 0.0298131 | 0.0216011   |
| Mean      | 249.90 s    | 252.79 s  | 255.81 s    |
| Std. Dev. | 2.9989 s    | 5.0869 s  | 5.9781 s    |
| Median    | 248.61 s    | 251.05 s  | 258.07 s    |
| MAD       | 1.0386 s    | 5.6726 s  | 8.3296 s    |

### 3a5671c uvfits_output - 1254670392_avg

![3a5671c uvfits_output - 1254670392_avg pdf](img/3a5671c%20uvfits_output%20-%201254670392_avg%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.1386326   | 0.1791628 | 0.1337237   |
| Mean      | 7.7458 s    | 7.7792 s  | 7.8152 s    |
| Std. Dev. | 30.948 ms   | 59.036 ms | 73.729 ms   |
| Median    | 7.7353 s    | 7.7719 s  | 7.8286 s    |
| MAD       | 14.180 ms   | 52.486 ms | 111.55 ms   |

## [commit ee893da](https://github.com/MWATelescope/Birli/actions/runs/1168717144)

performance changes:

- use ndarray everywhere

total duration: 35m 13s

### ee893da context_to_jones_array - ord_half_1196175296

![ee893da context_to_jones_array - ord_half_1196175296 pdf](img/ee893da%20context_to_jones_array%20-%20ord_half_1196175296%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0354699   | 0.0427925 | 0.0288164   |
| Mean      | 69.038 s    | 70.697 s  | 73.240 s    |
| Std. Dev. | 941.93 ms   | 3.7655 s  | 5.7054 s    |
| Median    | 68.634 s    | 69.785 s  | 70.937 s    |
| MAD       | 91.021 ms   | 1.5115 s  | 3.1542 s    |

### ee893da context_to_jones_array - mwax_half_1247842824

![ee893da context_to_jones_array - mwax_half_1247842824 pdf](img/ee893da%20context_to_jones_array%20-%20mwax_half_1247842824%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0122324   | 0.0164283 | 0.0118829   |
| Mean      | 37.796 s    | 37.996 s  | 38.208 s    |
| Std. Dev. | 142.38 ms   | 358.17 ms | 461.30 ms   |
| Median    | 37.773 s    | 37.948 s  | 38.294 s    |
| MAD       | 8.3125 ms   | 204.59 ms | 603.20 ms   |

### ee893da correct_cable_lengths - ord_half_1196175296

![ee893da correct_cable_lengths - ord_half_1196175296 pdf](img/ee893da%20correct_cable_lengths%20-%20ord_half_1196175296%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.3145424   | 0.4039586 | 0.3459765   |
| Mean      | 18.255 s    | 18.495 s  | 18.678 s    |
| Std. Dev. | 124.87 ms   | 365.04 ms | 525.90 ms   |
| Median    | 18.378 s    | 18.583 s  | 18.743 s    |
| MAD       | 31.692 ms   | 202.47 ms | 445.61 ms   |

### ee893da correct_cable_lengths - mwax_half_1247842824

![ee893da correct_cable_lengths - mwax_half_1247842824 pdf](img/ee893da%20correct_cable_lengths%20-%20mwax_half_1247842824%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0003975   | 0.0005287 | 0.0003714   |
| Mean      | 9.2153 s    | 9.2928 s  | 9.3807 s    |
| Std. Dev. | 49.930 ms   | 142.27 ms | 171.98 ms   |
| Median    | 9.1880 s    | 9.2344 s  | 9.4055 s    |
| MAD       | 16.507 ms   | 75.859 ms | 231.62 ms   |

### ee893da correct_geometry - ord_half_1196175296

![ee893da correct_geometry - ord_half_1196175296 pdf](img/ee893da%20correct_geometry%20-%20ord_half_1196175296%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0097952   | 0.0128601 | 0.0089791   |
| Mean      | 9.3852 s    | 9.5513 s  | 9.7464 s    |
| Std. Dev. | 138.49 ms   | 310.87 ms | 419.62 ms   |
| Median    | 9.3199 s    | 9.4600 s  | 9.7426 s    |
| MAD       | 51.771 ms   | 271.98 ms | 492.18 ms   |

### ee893da correct_geometry - mwax_half_1247842824

![ee893da correct_geometry - mwax_half_1247842824 pdf](img/ee893da%20correct_geometry%20-%20mwax_half_1247842824%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0145496   | 0.0195230 | 0.0139053   |
| Mean      | 5.9858 s    | 6.5979 s  | 7.2634 s    |
| Std. Dev. | 643.79 ms   | 1.0928 s  | 1.3857 s    |
| Median    | 5.6613 s    | 6.4340 s  | 7.3634 s    |
| MAD       | 215.60 ms   | 1.2470 s  | 2.0262 s    |

### ee893da uvfits_output - 1254670392_avg

![ee893da uvfits_output - 1254670392_avg pdf](img/ee893da%20uvfits_output%20-%201254670392_avg%20pdf.svg)

|           | Lower bound | Estimate  | Upper bound |
| --------- | ----------- | --------- | ----------- |
| R²        | 0.0129782   | 0.0179183 | 0.0130896   |
| Mean      | 4.8250 s    | 4.8460 s  | 4.8666 s    |
| Std. Dev. | 18.218 ms   | 35.451 ms | 45.378 ms   |
| Median    | 4.8170 s    | 4.8457 s  | 4.8729 s    |
| MAD       | 4.1059 ms   | 16.946 ms | 72.153 ms   |

## Summary

| bench                                              | goal  | a762a27 | ea16381 | dbdc0ee | 09abc69 | 3a5671c | ee893da | v0.1.9  | v0.3.1 |
| -------------------------------------------------- | ----- | ------- | ------- | ------- | ------- | ------- | ------- | ------- | ------ |
| context_to_baseline_imgsets - ord_half_1196175296  |       | 1873.5s | 729.49s | 458.65s | 172.69s | 98.537s |         |         |
| context_to_baseline_imgsets - mwax_half_1247842824 |       |         | 372.64s | 226.30s | 87.029s | 51.940s |         |         |
| context_to_jones_array - ord_half_1196175296       | ~300s |         |         |         |         | 64.717s | 70.697s | 57.96s  | 54.79s |
| context_to_jones_array - mwax_half_1247842824      |       |         |         |         |         | 35.321s | 37.996s | 30.79s  | 28.15s |
| correct_cable_lengths - ord_half_1196175296        |       |         |         |         |         | 10.543s | 18.495s | 17.30s  | 16.51s |
| correct_cable_lengths - mwax_half_1247842824       |       |         |         |         |         | 6.3826s | 9.2928s | 8.78s   | 8.96s  |
| correct_geometry - ord_half_1196175296             |       |         |         |         |         | 650.03s | 9.5513s | 7.80s   | 7.06s  |
| correct_geometry - mwax_half_1247842824            |       |         |         |         |         | 252.79s | 6.5979s | 5.64s   | 5.41s  |
| uvfits_output - 1254670392_avg                     | ~2s   |         |         |         |         | 7.7792s | 4.8460s | 2.48s   | 5.16s  |
| uvfits_output - ord_half_1196175296                |       |         |         |         |         |         |         | 676.46s | -      |
| uvfits_output - mwax_half_1247842824               |       |         |         |         |         |         |         | 328.83s | -      |

## simplified output only

| bench                                | v0.4.0 |
| ------------------------------------ | ------ |
| uvfits_out - 1196175296 10 timesteps | 98.52  |
| ms_out - 1196175296 10 timesteps     | 89.74  |

## MWALib-only testing

|                                                      | mean  | std dev |
| ---------------------------------------------------- | ----- | ------- |
| read_by_baseline - mwax_half_1247842824              | 4.172 | 0.364   |
| read_by_baseline - ord_half_1196175296               | 3.889 | 0.374   |
| read_by_baseline_into_buffer - mwax_half_1247842824  | 3.928 | 0.324   |
| read_by_baseline_into_buffer - ord_half_1196175296   | 3.773 | 0.228   |
| read_by_frequency - mwax_half_1247842824             | 3.124 | 0.198   |
| read_by_frequency - ord_half_1196175296              | 2.950 | 0.180   |
| read_by_frequency_into_buffer - mwax_half_1247842824 | 3.019 | 0.161   |
| read_by_frequency_into_buffer - ord_half_1196175296  | 2.889 | 0.127   |

## on garrawarla

```txt
04:13:13 /astro/mwaops/gsleap/1090008640$ time singularity exec --bind /nvmetmp:/nvmetmp /pawsey/mwa/singularity/birli/birli_latest.sif /app/target/release/birli aoflagger -m 1090008640.metafits -u /nvmetmp/1090008640.uvfits 1090008640*gpubox*_00.fits
coarse_chan 000 : [=================================================================================================================] 53  /53
coarse_chan 001 : [=================================================================================================================] 53  /53
coarse_chan 002 : [=================================================================================================================] 53  /53
coarse_chan 003 : [=================================================================================================================] 53  /53
coarse_chan 004 : [=================================================================================================================] 53  /53
coarse_chan 005 : [=================================================================================================================] 53  /53
coarse_chan 006 : [=================================================================================================================] 53  /53
coarse_chan 007 : [=================================================================================================================] 53  /53
coarse_chan 008 : [=================================================================================================================] 53  /53
coarse_chan 009 : [=================================================================================================================] 53  /53
coarse_chan 010 : [=================================================================================================================] 53  /53
coarse_chan 011 : [=================================================================================================================] 53  /53
coarse_chan 012 : [=================================================================================================================] 53  /53
coarse_chan 013 : [=================================================================================================================] 53  /53
coarse_chan 014 : [=================================================================================================================] 53  /53
coarse_chan 015 : [=================================================================================================================] 53  /53
coarse_chan 016 : [=================================================================================================================] 53  /53
coarse_chan 017 : [=================================================================================================================] 53  /53
coarse_chan 018 : [=================================================================================================================] 53  /53
coarse_chan 019 : [=================================================================================================================] 53  /53
coarse_chan 020 : [=================================================================================================================] 53  /53
coarse_chan 021 : [=================================================================================================================] 53  /53
coarse_chan 022 : [=================================================================================================================] 53  /53
coarse_chan 023 : [=================================================================================================================] 53  /53
loading hdus    : [00:00:01] [===================================================================================================] 100% (0s   )
cable corrections: [00:00:01] [==================================================================================================] 100% (0s   )
flagging b'lines: [00:00:16] [===================================================================================================] 100% (0s   )
geom corrections: [00:00:01] [===================================================================================================] 100% (0s   )
write uv vis    : [00:00:29] [===================================================================================================] 100% (0s   )

real    1m27.251s
user    12m53.060s
04:15:25 /astro/mwaops/gsleap/1090008640$ time singularity exec --bind /nvmetmp:/nvmetmp /pawsey/mwa/singularity/cotter/cotter_latest.sif cotter -m 1090008640.metafits -o /nvmetmp/1090008640.uvfits 1090008640*gpubox*_00.fits
Running Cotter MWA preprocessing pipeline, version 4.5 (2020-08-10).
Flagging is performed by AOFlagger 3.0 (2020-07-21).
Input filenames succesfully parsed: using 24 files covering 1 timeranges from 24 GPU boxes.
Detected 377.6 GB of system memory.
Ignored keyword: ATTEN_DB
Ignored keyword: CHANSEL
Ignored keyword: INSTRUME
Ignored keyword: QUACKTIM
Ignored keyword: GOODTIME
Ignored keyword: COMMENT
Ignored keyword: COMMENT
Ignored keyword: COMMENT
Ignored keyword: COMMENT
Ignored keyword: COMMENT
Ignored keyword: COMMENT
Ignored keyword: COMMENT
Ignored keyword: COMMENT
Ignored keyword: HISTORY
Ignored keyword: HISTORY
Ignored keyword: HISTORY
Observation covers 167.04-197.72 MHz.
Output resolution: 2 s / 40 kHz (time avg: 1x, freq avg: 1x).
The first 2 samples (4 s), last 0 samples (0 s) and 2 edge channels will be flagged.
Using a-priori subband passband with 32 channels.
Using per-input subband gains. Average gains: 1.03369,1.04211,1.04211,1.04553,1.052,1.05872,1.07434,1.08765,1.10815,1.12842,1.14075,1.15637,1.16956,1.18518,1.20654,1.22681,1.2511,1.27039,1.29553,1.31116,1.33044,1.34607,1.36304,1.37866
Subband order: 23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0
Observation's bandwidth is contiguous.
Antenna positions are written in LOCAL MERIDIAN
All 56 scans fit in memory; no partitioning necessary.
=== Processing chunk 1 of 1 ===
There are 54 HDUs in file 1090008640_20140721201027_gpubox24_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox23_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox22_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox21_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox20_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox19_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox18_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox17_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox16_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox15_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox14_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox13_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox12_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox11_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox10_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox09_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox08_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox07_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox06_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox05_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox04_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox03_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox02_00.fits
There are 54 HDUs in file 1090008640_20140721201027_gpubox01_00.fits
Will stop on HDU 54.
Reading GPU files: 0%....10%....20%....30%....40%....50%....60%....70%....80%....90%....100%
WARNING: start time according to raw files is 2014-07-21 20:10:25,
but meta files say 2014-07-21 20:10:24 !
Will use start time from raw file, which should be most accurate.
Warning: header specifies 56 scans, but there are only 53 in the data.
Last 3 scan(s) will be flagged.
RFI detection, conjugations, subband ordering and cable length corrections:
 0%....10%....20%....30%....40%....50%....60%....70%....80%....90%....100%
Writing: 0%....10%....20%....30%....40%....50%....60%....70%....80%....90%....100%
Writing MWA fields to UVFits file...
Wall-clock time in reading: 00:00:24.295526 processing: 00:00:14.817472 writing: 00:01:04.135953

real    1m44.495s
user    12m37.101s
sys     6m17.461s
```

## huge file:

```txt
Running Cotter MWA preprocessing pipeline, version 4.5 (2020-08-10).
Flagging is performed by AOFlagger 3.0 (2020-07-21).
Input filenames succesfully parsed: using 48 files covering 2 timeranges from 24 GPU boxes.
Detected 377.6 GB of system memory.
Ignored keyword: DELAYMOD
Ignored keyword: CABLEDEL
Ignored keyword: GEODEL
Ignored keyword: CALIBDEL
Ignored keyword: DELDESC
Ignored keyword: ATTEN_DB
Ignored keyword: CALIBSRC
Ignored keyword: CHANSEL
Ignored keyword: INSTRUME
Ignored keyword: QUACKTIM
Ignored keyword: GOODTIME
Ignored keyword: COMMENT
Ignored keyword: COMMENT
Ignored keyword: COMMENT
Ignored keyword: COMMENT
Ignored keyword: COMMENT
Ignored keyword: COMMENT
Ignored keyword: COMMENT
Ignored keyword: COMMENT
Ignored keyword: HISTORY
Ignored keyword: HISTORY
Ignored keyword: HISTORY
Observation covers 169.6-200.31 MHz.
Output resolution: 0.5 s / 10 kHz (time avg: 1x, freq avg: 1x).
The first 8 samples (4 s), last 0 samples (0 s) and 8 edge channels will be flagged.
Using a-priori subband passband with 128 channels.
Using per-input subband gains. Average gains: 1.09082,1.0918,1.0957,1.10864,1.12793,1.15283,1.1748,1.19971,1.21509,1.22144,1.23438,1.23804,1.24194,1.25757,1.27686,1.30811,1.32373,1.34521,1.35742,1.3728,1.37646,1.38037,1.38672,1.40601
Subband order: 23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0
Observation's bandwidth is contiguous.
Antenna positions are written in LOCAL MERIDIAN
All 240 scans fit in memory; no partitioning necessary.
=== Processing chunk 1 of 1 ===
There are 121 HDUs in file 1303162952_20210422214216_gpubox24_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox23_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox22_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox21_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox20_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox19_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox18_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox17_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox16_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox15_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox14_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox13_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox12_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox11_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox10_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox09_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox08_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox07_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox06_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox05_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox04_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox03_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox02_00.fits
There are 121 HDUs in file 1303162952_20210422214216_gpubox01_00.fits
Will stop on HDU 121.
Reading GPU files: 0%....10%....20%....30%....40%....50%....60%....70%....80%....90%....100%
WARNING: start time according to raw files is 2021-04-22 21:42:16,
but meta files say 2021-04-22 21:42:14 !
Will use start time from raw file, which should be most accurate.
There are 109 HDUs in file 1303162952_20210422214316_gpubox24_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox23_01.fits
There are 111 HDUs in file 1303162952_20210422214316_gpubox22_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox21_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox20_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox19_01.fits
There are 111 HDUs in file 1303162952_20210422214316_gpubox18_01.fits
There are 111 HDUs in file 1303162952_20210422214316_gpubox17_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox16_01.fits
There are 111 HDUs in file 1303162952_20210422214316_gpubox15_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox14_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox13_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox12_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox11_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox10_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox09_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox08_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox07_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox06_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox05_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox04_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox03_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox02_01.fits
There are 109 HDUs in file 1303162952_20210422214316_gpubox01_01.fits
WARNING: Files had not the same number of HDUs.
Will stop on HDU 109.
Reading GPU files: 0%....10%....20%....30%....40%....50%....60%....70%....80%....90%....100%
Warning: header specifies 240 scans, but there are only 228 in the data.
Last 12 scan(s) will be flagged.
RFI detection, conjugations, subband ordering and cable length corrections:
 0%....10%....20%....30%....40%....50%....60%....70%....80%....90%....100%
Writing: 0%....10%....20%....30%....40%....50%....60%....70%....80%....90%....100%
Writing MWA fields to UVFits file...
Wall-clock time in reading: 00:12:10.804084 processing: 00:05:18.828621 writing: 00:21:31.543545
```

```
[2021-09-24T09:04:48Z INFO  birli] start main
[2021-09-24T09:04:48Z DEBUG birli] args:
    Args { inner: ["/app/target/release/birli", "aoflagger", "-m", "1303162952.metafits", "-u", "/nvmetmp/1303162952.uvfits", "1303162952_20210422214216_gpubox01_00.fits", "1303162952_20210422214216_gpubox02_00.fits", "1303162952_20210422214216_gpubox03_00.fits", "1303162952_20210422214216_gpubox04_00.fits", "1303162952_20210422214216_gpubox05_00.fits", "1303162952_20210422214216_gpubox06_00.fits", "1303162952_20210422214216_gpubox07_00.fits", "1303162952_20210422214216_gpubox08_00.fits", "1303162952_20210422214216_gpubox09_00.fits", "1303162952_20210422214216_gpubox10_00.fits", "1303162952_20210422214216_gpubox11_00.fits", "1303162952_20210422214216_gpubox12_00.fits", "1303162952_20210422214216_gpubox13_00.fits", "1303162952_20210422214216_gpubox14_00.fits", "1303162952_20210422214216_gpubox15_00.fits", "1303162952_20210422214216_gpubox16_00.fits", "1303162952_20210422214216_gpubox17_00.fits", "1303162952_20210422214216_gpubox18_00.fits", "1303162952_20210422214216_gpubox19_00.fits", "1303162952_20210422214216_gpubox20_00.fits", "1303162952_20210422214216_gpubox21_00.fits", "1303162952_20210422214216_gpubox22_00.fits", "1303162952_20210422214216_gpubox23_00.fits", "1303162952_20210422214216_gpubox24_00.fits", "1303162952_20210422214316_gpubox01_01.fits", "1303162952_20210422214316_gpubox02_01.fits", "1303162952_20210422214316_gpubox03_01.fits", "1303162952_20210422214316_gpubox04_01.fits", "1303162952_20210422214316_gpubox05_01.fits", "1303162952_20210422214316_gpubox06_01.fits", "1303162952_20210422214316_gpubox07_01.fits", "1303162952_20210422214316_gpubox08_01.fits", "1303162952_20210422214316_gpubox09_01.fits", "1303162952_20210422214316_gpubox10_01.fits", "1303162952_20210422214316_gpubox11_01.fits", "1303162952_20210422214316_gpubox12_01.fits", "1303162952_20210422214316_gpubox13_01.fits", "1303162952_20210422214316_gpubox14_01.fits", "1303162952_20210422214316_gpubox15_01.fits", "1303162952_20210422214316_gpubox16_01.fits", "1303162952_20210422214316_gpubox17_01.fits", "1303162952_20210422214316_gpubox18_01.fits", "1303162952_20210422214316_gpubox19_01.fits", "1303162952_20210422214316_gpubox20_01.fits", "1303162952_20210422214316_gpubox21_01.fits", "1303162952_20210422214316_gpubox22_01.fits", "1303162952_20210422214316_gpubox23_01.fits", "1303162952_20210422214316_gpubox24_01.fits"] }
[2021-09-24T09:04:48Z DEBUG birli] arg matches:
    ArgMatches { args: {}, subcommand: Some(SubCommand { name: "aoflagger", matches: ArgMatches { args: {"fits-files": MatchedArg { occurs: 48, indices: [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52], vals: ["1303162952_20210422214216_gpubox01_00.fits", "1303162952_20210422214216_gpubox02_00.fits", "1303162952_20210422214216_gpubox03_00.fits", "1303162952_20210422214216_gpubox04_00.fits", "1303162952_20210422214216_gpubox05_00.fits", "1303162952_20210422214216_gpubox06_00.fits", "1303162952_20210422214216_gpubox07_00.fits", "1303162952_20210422214216_gpubox08_00.fits", "1303162952_20210422214216_gpubox09_00.fits", "1303162952_20210422214216_gpubox10_00.fits", "1303162952_20210422214216_gpubox11_00.fits", "1303162952_20210422214216_gpubox12_00.fits", "1303162952_20210422214216_gpubox13_00.fits", "1303162952_20210422214216_gpubox14_00.fits", "1303162952_20210422214216_gpubox15_00.fits", "1303162952_20210422214216_gpubox16_00.fits", "1303162952_20210422214216_gpubox17_00.fits", "1303162952_20210422214216_gpubox18_00.fits", "1303162952_20210422214216_gpubox19_00.fits", "1303162952_20210422214216_gpubox20_00.fits", "1303162952_20210422214216_gpubox21_00.fits", "1303162952_20210422214216_gpubox22_00.fits", "1303162952_20210422214216_gpubox23_00.fits", "1303162952_20210422214216_gpubox24_00.fits", "1303162952_20210422214316_gpubox01_01.fits", "1303162952_20210422214316_gpubox02_01.fits", "1303162952_20210422214316_gpubox03_01.fits", "1303162952_20210422214316_gpubox04_01.fits", "1303162952_20210422214316_gpubox05_01.fits", "1303162952_20210422214316_gpubox06_01.fits", "1303162952_20210422214316_gpubox07_01.fits", "1303162952_20210422214316_gpubox08_01.fits", "1303162952_20210422214316_gpubox09_01.fits", "1303162952_20210422214316_gpubox10_01.fits", "1303162952_20210422214316_gpubox11_01.fits", "1303162952_20210422214316_gpubox12_01.fits", "1303162952_20210422214316_gpubox13_01.fits", "1303162952_20210422214316_gpubox14_01.fits", "1303162952_20210422214316_gpubox15_01.fits", "1303162952_20210422214316_gpubox16_01.fits", "1303162952_20210422214316_gpubox17_01.fits", "1303162952_20210422214316_gpubox18_01.fits", "1303162952_20210422214316_gpubox19_01.fits", "1303162952_20210422214316_gpubox20_01.fits", "1303162952_20210422214316_gpubox21_01.fits", "1303162952_20210422214316_gpubox22_01.fits", "1303162952_20210422214316_gpubox23_01.fits", "1303162952_20210422214316_gpubox24_01.fits"] }, "metafits": MatchedArg { occurs: 1, indices: [2], vals: ["1303162952.metafits"] }, "uvfits-out": MatchedArg { occurs: 1, indices: [4], vals: ["/nvmetmp/1303162952.uvfits"] }}, subcommand: None, usage: Some("USAGE:\n    birli aoflagger [FLAGS] [OPTIONS] <fits-files>... -m <metafits>") } }), usage: Some("USAGE:\n    birli [SUBCOMMAND]") }
[2021-09-24T09:04:50Z DEBUG birli] mwalib correlator context:
    CorrelatorContext (
                Metafits Context:           MetafitsContext (
        obsid:                     1303162952,
        mode:                      HW_LFILES,

        If Correlator Mode:
         fine channel resolution:  10 kHz,
         integration time:         0.50 s
         num fine channels/coarse: 128,

        If Voltage Mode:
         fine channel resolution:  10 kHz,
         num fine channels/coarse: 128,

        Geometric delays applied          : No,
        Cable length corrections applied  : false,
        Calibration delays & gains applied: false,

        Creator:                  andrew,
        Project ID:               G0079,
        Observation Name:         Transient_cal_3C444,
        Receivers:                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
        Delays:                   [3, 8, 13, 18, 2, 7, 12, 17, 1, 6, 11, 16, 0, 5, 10, 15],
        Global attenuation:       1 dB,

        Scheduled start (UNIX)    1619127734,
        Scheduled end (UNIX)      1619127854,
        Scheduled start (GPS)     1303162952,
        Scheduled end (GPS)       1303163072,
        Scheduled start (utc)     2021-04-22 21:42:14 +00:00,
        Scheduled end (utc)       2021-04-22 21:44:14 +00:00,
        Scheduled start (MJD)     59326.90432870371,
        Scheduled end (MJD)       59326.9057175926,
        Scheduled duration        120 s,
        Quack time:               0.5 s,
        Good UNIX start time:     1619127734.5,

        Num timesteps:            240,
        Timesteps:                [unix=1619127734.000, gps=1303162952.000, unix=1619127734.500, gps=1303162952.500, unix=1619127735.000, gps=1303162953.000, unix=1619127735.500, gps=1303162953.500, unix=1619127736.000, gps=1303162954.000, unix=1619127736.500, gps=1303162954.500, unix=1619127737.000, gps=1303162955.000, unix=1619127737.500, gps=1303162955.500, unix=1619127738.000, gps=1303162956.000, unix=1619127738.500, gps=1303162956.500, unix=1619127739.000, gps=1303162957.000, unix=1619127739.500, gps=1303162957.500, unix=1619127740.000, gps=1303162958.000, unix=1619127740.500, gps=1303162958.500, unix=1619127741.000, gps=1303162959.000, unix=1619127741.500, gps=1303162959.500, unix=1619127742.000, gps=1303162960.000, unix=1619127742.500, gps=1303162960.500, unix=1619127743.000, gps=1303162961.000, unix=1619127743.500, gps=1303162961.500, unix=1619127744.000, gps=1303162962.000, unix=1619127744.500, gps=1303162962.500, unix=1619127745.000, gps=1303162963.000, unix=1619127745.500, gps=1303162963.500, unix=1619127746.000, gps=1303162964.000, unix=1619127746.500, gps=1303162964.500, unix=1619127747.000, gps=1303162965.000, unix=1619127747.500, gps=1303162965.500, unix=1619127748.000, gps=1303162966.000, unix=1619127748.500, gps=1303162966.500, unix=1619127749.000, gps=1303162967.000, unix=1619127749.500, gps=1303162967.500, unix=1619127750.000, gps=1303162968.000, unix=1619127750.500, gps=1303162968.500, unix=1619127751.000, gps=1303162969.000, unix=1619127751.500, gps=1303162969.500, unix=1619127752.000, gps=1303162970.000, unix=1619127752.500, gps=1303162970.500, unix=1619127753.000, gps=1303162971.000, unix=1619127753.500, gps=1303162971.500, unix=1619127754.000, gps=1303162972.000, unix=1619127754.500, gps=1303162972.500, unix=1619127755.000, gps=1303162973.000, unix=1619127755.500, gps=1303162973.500, unix=1619127756.000, gps=1303162974.000, unix=1619127756.500, gps=1303162974.500, unix=1619127757.000, gps=1303162975.000, unix=1619127757.500, gps=1303162975.500, unix=1619127758.000, gps=1303162976.000, unix=1619127758.500, gps=1303162976.500, unix=1619127759.000, gps=1303162977.000, unix=1619127759.500, gps=1303162977.500, unix=1619127760.000, gps=1303162978.000, unix=1619127760.500, gps=1303162978.500, unix=1619127761.000, gps=1303162979.000, unix=1619127761.500, gps=1303162979.500, unix=1619127762.000, gps=1303162980.000, unix=1619127762.500, gps=1303162980.500, unix=1619127763.000, gps=1303162981.000, unix=1619127763.500, gps=1303162981.500, unix=1619127764.000, gps=1303162982.000, unix=1619127764.500, gps=1303162982.500, unix=1619127765.000, gps=1303162983.000, unix=1619127765.500, gps=1303162983.500, unix=1619127766.000, gps=1303162984.000, unix=1619127766.500, gps=1303162984.500, unix=1619127767.000, gps=1303162985.000, unix=1619127767.500, gps=1303162985.500, unix=1619127768.000, gps=1303162986.000, unix=1619127768.500, gps=1303162986.500, unix=1619127769.000, gps=1303162987.000, unix=1619127769.500, gps=1303162987.500, unix=1619127770.000, gps=1303162988.000, unix=1619127770.500, gps=1303162988.500, unix=1619127771.000, gps=1303162989.000, unix=1619127771.500, gps=1303162989.500, unix=1619127772.000, gps=1303162990.000, unix=1619127772.500, gps=1303162990.500, unix=1619127773.000, gps=1303162991.000, unix=1619127773.500, gps=1303162991.500, unix=1619127774.000, gps=1303162992.000, unix=1619127774.500, gps=1303162992.500, unix=1619127775.000, gps=1303162993.000, unix=1619127775.500, gps=1303162993.500, unix=1619127776.000, gps=1303162994.000, unix=1619127776.500, gps=1303162994.500, unix=1619127777.000, gps=1303162995.000, unix=1619127777.500, gps=1303162995.500, unix=1619127778.000, gps=1303162996.000, unix=1619127778.500, gps=1303162996.500, unix=1619127779.000, gps=1303162997.000, unix=1619127779.500, gps=1303162997.500, unix=1619127780.000, gps=1303162998.000, unix=1619127780.500, gps=1303162998.500, unix=1619127781.000, gps=1303162999.000, unix=1619127781.500, gps=1303162999.500, unix=1619127782.000, gps=1303163000.000, unix=1619127782.500, gps=1303163000.500, unix=1619127783.000, gps=1303163001.000, unix=1619127783.500, gps=1303163001.500, unix=1619127784.000, gps=1303163002.000, unix=1619127784.500, gps=1303163002.500, unix=1619127785.000, gps=1303163003.000, unix=1619127785.500, gps=1303163003.500, unix=1619127786.000, gps=1303163004.000, unix=1619127786.500, gps=1303163004.500, unix=1619127787.000, gps=1303163005.000, unix=1619127787.500, gps=1303163005.500, unix=1619127788.000, gps=1303163006.000, unix=1619127788.500, gps=1303163006.500, unix=1619127789.000, gps=1303163007.000, unix=1619127789.500, gps=1303163007.500, unix=1619127790.000, gps=1303163008.000, unix=1619127790.500, gps=1303163008.500, unix=1619127791.000, gps=1303163009.000, unix=1619127791.500, gps=1303163009.500, unix=1619127792.000, gps=1303163010.000, unix=1619127792.500, gps=1303163010.500, unix=1619127793.000, gps=1303163011.000, unix=1619127793.500, gps=1303163011.500, unix=1619127794.000, gps=1303163012.000, unix=1619127794.500, gps=1303163012.500, unix=1619127795.000, gps=1303163013.000, unix=1619127795.500, gps=1303163013.500, unix=1619127796.000, gps=1303163014.000, unix=1619127796.500, gps=1303163014.500, unix=1619127797.000, gps=1303163015.000, unix=1619127797.500, gps=1303163015.500, unix=1619127798.000, gps=1303163016.000, unix=1619127798.500, gps=1303163016.500, unix=1619127799.000, gps=1303163017.000, unix=1619127799.500, gps=1303163017.500, unix=1619127800.000, gps=1303163018.000, unix=1619127800.500, gps=1303163018.500, unix=1619127801.000, gps=1303163019.000, unix=1619127801.500, gps=1303163019.500, unix=1619127802.000, gps=1303163020.000, unix=1619127802.500, gps=1303163020.500, unix=1619127803.000, gps=1303163021.000, unix=1619127803.500, gps=1303163021.500, unix=1619127804.000, gps=1303163022.000, unix=1619127804.500, gps=1303163022.500, unix=1619127805.000, gps=1303163023.000, unix=1619127805.500, gps=1303163023.500, unix=1619127806.000, gps=1303163024.000, unix=1619127806.500, gps=1303163024.500, unix=1619127807.000, gps=1303163025.000, unix=1619127807.500, gps=1303163025.500, unix=1619127808.000, gps=1303163026.000, unix=1619127808.500, gps=1303163026.500, unix=1619127809.000, gps=1303163027.000, unix=1619127809.500, gps=1303163027.500, unix=1619127810.000, gps=1303163028.000, unix=1619127810.500, gps=1303163028.500, unix=1619127811.000, gps=1303163029.000, unix=1619127811.500, gps=1303163029.500, unix=1619127812.000, gps=1303163030.000, unix=1619127812.500, gps=1303163030.500, unix=1619127813.000, gps=1303163031.000, unix=1619127813.500, gps=1303163031.500, unix=1619127814.000, gps=1303163032.000, unix=1619127814.500, gps=1303163032.500, unix=1619127815.000, gps=1303163033.000, unix=1619127815.500, gps=1303163033.500, unix=1619127816.000, gps=1303163034.000, unix=1619127816.500, gps=1303163034.500, unix=1619127817.000, gps=1303163035.000, unix=1619127817.500, gps=1303163035.500, unix=1619127818.000, gps=1303163036.000, unix=1619127818.500, gps=1303163036.500, unix=1619127819.000, gps=1303163037.000, unix=1619127819.500, gps=1303163037.500, unix=1619127820.000, gps=1303163038.000, unix=1619127820.500, gps=1303163038.500, unix=1619127821.000, gps=1303163039.000, unix=1619127821.500, gps=1303163039.500, unix=1619127822.000, gps=1303163040.000, unix=1619127822.500, gps=1303163040.500, unix=1619127823.000, gps=1303163041.000, unix=1619127823.500, gps=1303163041.500, unix=1619127824.000, gps=1303163042.000, unix=1619127824.500, gps=1303163042.500, unix=1619127825.000, gps=1303163043.000, unix=1619127825.500, gps=1303163043.500, unix=1619127826.000, gps=1303163044.000, unix=1619127826.500, gps=1303163044.500, unix=1619127827.000, gps=1303163045.000, unix=1619127827.500, gps=1303163045.500, unix=1619127828.000, gps=1303163046.000, unix=1619127828.500, gps=1303163046.500, unix=1619127829.000, gps=1303163047.000, unix=1619127829.500, gps=1303163047.500, unix=1619127830.000, gps=1303163048.000, unix=1619127830.500, gps=1303163048.500, unix=1619127831.000, gps=1303163049.000, unix=1619127831.500, gps=1303163049.500, unix=1619127832.000, gps=1303163050.000, unix=1619127832.500, gps=1303163050.500, unix=1619127833.000, gps=1303163051.000, unix=1619127833.500, gps=1303163051.500, unix=1619127834.000, gps=1303163052.000, unix=1619127834.500, gps=1303163052.500, unix=1619127835.000, gps=1303163053.000, unix=1619127835.500, gps=1303163053.500, unix=1619127836.000, gps=1303163054.000, unix=1619127836.500, gps=1303163054.500, unix=1619127837.000, gps=1303163055.000, unix=1619127837.500, gps=1303163055.500, unix=1619127838.000, gps=1303163056.000, unix=1619127838.500, gps=1303163056.500, unix=1619127839.000, gps=1303163057.000, unix=1619127839.500, gps=1303163057.500, unix=1619127840.000, gps=1303163058.000, unix=1619127840.500, gps=1303163058.500, unix=1619127841.000, gps=1303163059.000, unix=1619127841.500, gps=1303163059.500, unix=1619127842.000, gps=1303163060.000, unix=1619127842.500, gps=1303163060.500, unix=1619127843.000, gps=1303163061.000, unix=1619127843.500, gps=1303163061.500, unix=1619127844.000, gps=1303163062.000, unix=1619127844.500, gps=1303163062.500, unix=1619127845.000, gps=1303163063.000, unix=1619127845.500, gps=1303163063.500, unix=1619127846.000, gps=1303163064.000, unix=1619127846.500, gps=1303163064.500, unix=1619127847.000, gps=1303163065.000, unix=1619127847.500, gps=1303163065.500, unix=1619127848.000, gps=1303163066.000, unix=1619127848.500, gps=1303163066.500, unix=1619127849.000, gps=1303163067.000, unix=1619127849.500, gps=1303163067.500, unix=1619127850.000, gps=1303163068.000, unix=1619127850.500, gps=1303163068.500, unix=1619127851.000, gps=1303163069.000, unix=1619127851.500, gps=1303163069.500, unix=1619127852.000, gps=1303163070.000, unix=1619127852.500, gps=1303163070.500, unix=1619127853.000, gps=1303163071.000, unix=1619127853.500, gps=1303163071.500],

        Num coarse channels:      24,
        Coarse Channels:          [gpu=24 corr=23 rec=133 @ 170.240 MHz, gpu=23 corr=22 rec=134 @ 171.520 MHz, gpu=22 corr=21 rec=135 @ 172.800 MHz, gpu=21 corr=20 rec=136 @ 174.080 MHz, gpu=20 corr=19 rec=137 @ 175.360 MHz, gpu=19 corr=18 rec=138 @ 176.640 MHz, gpu=18 corr=17 rec=139 @ 177.920 MHz, gpu=17 corr=16 rec=140 @ 179.200 MHz, gpu=16 corr=15 rec=141 @ 180.480 MHz, gpu=15 corr=14 rec=142 @ 181.760 MHz, gpu=14 corr=13 rec=143 @ 183.040 MHz, gpu=13 corr=12 rec=144 @ 184.320 MHz, gpu=12 corr=11 rec=145 @ 185.600 MHz, gpu=11 corr=10 rec=146 @ 186.880 MHz, gpu=10 corr=9 rec=147 @ 188.160 MHz, gpu=9 corr=8 rec=148 @ 189.440 MHz, gpu=8 corr=7 rec=149 @ 190.720 MHz, gpu=7 corr=6 rec=150 @ 192.000 MHz, gpu=6 corr=5 rec=151 @ 193.280 MHz, gpu=5 corr=4 rec=152 @ 194.560 MHz, gpu=4 corr=3 rec=153 @ 195.840 MHz, gpu=3 corr=2 rec=154 @ 197.120 MHz, gpu=2 corr=1 rec=155 @ 198.400 MHz, gpu=1 corr=0 rec=156 @ 199.680 MHz],

        Num fine channels:        3072,
        Fine Channels (kHz):      Map { iter: Iter([169600000.0, 169610000.0, 169620000.0, 169630000.0, 169640000.0, 169650000.0, 169660000.0, 169670000.0, 169680000.0, 169690000.0, 169700000.0, 169710000.0, 169720000.0, 169730000.0, 169740000.0, 169750000.0, 169760000.0, 169770000.0, 169780000.0, 169790000.0, 169800000.0, 169810000.0, 169820000.0, 169830000.0, 169840000.0, 169850000.0, 169860000.0, 169870000.0, 169880000.0, 169890000.0, 169900000.0, 169910000.0, 169920000.0, 169930000.0, 169940000.0, 169950000.0, 169960000.0, 169970000.0, 169980000.0, 169990000.0, 170000000.0, 170010000.0, 170020000.0, 170030000.0, 170040000.0, 170050000.0, 170060000.0, 170070000.0, 170080000.0, 170090000.0, 170100000.0, 170110000.0, 170120000.0, 170130000.0, 170140000.0, 170150000.0, 170160000.0, 170170000.0, 170180000.0, 170190000.0, 170200000.0, 170210000.0, 170220000.0, 170230000.0, 170240000.0, 170250000.0, 170260000.0, 170270000.0, 170280000.0, 170290000.0, 170300000.0, 170310000.0, 170320000.0, 170330000.0, 170340000.0, 170350000.0, 170360000.0, 170370000.0, 170380000.0, 170390000.0, 170400000.0, 170410000.0, 170420000.0, 170430000.0, 170440000.0, 170450000.0, 170460000.0, 170470000.0, 170480000.0, 170490000.0, 170500000.0, 170510000.0, 170520000.0, 170530000.0, 170540000.0, 170550000.0, 170560000.0, 170570000.0, 170580000.0, 170590000.0, 170600000.0, 170610000.0, 170620000.0, 170630000.0, 170640000.0, 170650000.0, 170660000.0, 170670000.0, 170680000.0, 170690000.0, 170700000.0, 170710000.0, 170720000.0, 170730000.0, 170740000.0, 170750000.0, 170760000.0, 170770000.0, 170780000.0, 170790000.0, 170800000.0, 170810000.0, 170820000.0, 170830000.0, 170840000.0, 170850000.0, 170860000.0, 170870000.0, 170880000.0, 170890000.0, 170900000.0, 170910000.0, 170920000.0, 170930000.0, 170940000.0, 170950000.0, 170960000.0, 170970000.0, 170980000.0, 170990000.0, 171000000.0, 171010000.0, 171020000.0, 171030000.0, 171040000.0, 171050000.0, 171060000.0, 171070000.0, 171080000.0, 171090000.0, 171100000.0, 171110000.0, 171120000.0, 171130000.0, 171140000.0, 171150000.0, 171160000.0, 171170000.0, 171180000.0, 171190000.0, 171200000.0, 171210000.0, 171220000.0, 171230000.0, 171240000.0, 171250000.0, 171260000.0, 171270000.0, 171280000.0, 171290000.0, 171300000.0, 171310000.0, 171320000.0, 171330000.0, 171340000.0, 171350000.0, 171360000.0, 171370000.0, 171380000.0, 171390000.0, 171400000.0, 171410000.0, 171420000.0, 171430000.0, 171440000.0, 171450000.0, 171460000.0, 171470000.0, 171480000.0, 171490000.0, 171500000.0, 171510000.0, 171520000.0, 171530000.0, 171540000.0, 171550000.0, 171560000.0, 171570000.0, 171580000.0, 171590000.0, 171600000.0, 171610000.0, 171620000.0, 171630000.0, 171640000.0, 171650000.0, 171660000.0, 171670000.0, 171680000.0, 171690000.0, 171700000.0, 171710000.0, 171720000.0, 171730000.0, 171740000.0, 171750000.0, 171760000.0, 171770000.0, 171780000.0, 171790000.0, 171800000.0, 171810000.0, 171820000.0, 171830000.0, 171840000.0, 171850000.0, 171860000.0, 171870000.0, 171880000.0, 171890000.0, 171900000.0, 171910000.0, 171920000.0, 171930000.0, 171940000.0, 171950000.0, 171960000.0, 171970000.0, 171980000.0, 171990000.0, 172000000.0, 172010000.0, 172020000.0, 172030000.0, 172040000.0, 172050000.0, 172060000.0, 172070000.0, 172080000.0, 172090000.0, 172100000.0, 172110000.0, 172120000.0, 172130000.0, 172140000.0, 172150000.0, 172160000.0, 172170000.0, 172180000.0, 172190000.0, 172200000.0, 172210000.0, 172220000.0, 172230000.0, 172240000.0, 172250000.0, 172260000.0, 172270000.0, 172280000.0, 172290000.0, 172300000.0, 172310000.0, 172320000.0, 172330000.0, 172340000.0, 172350000.0, 172360000.0, 172370000.0, 172380000.0, 172390000.0, 172400000.0, 172410000.0, 172420000.0, 172430000.0, 172440000.0, 172450000.0, 172460000.0, 172470000.0, 172480000.0, 172490000.0, 172500000.0, 172510000.0, 172520000.0, 172530000.0, 172540000.0, 172550000.0, 172560000.0, 172570000.0, 172580000.0, 172590000.0, 172600000.0, 172610000.0, 172620000.0, 172630000.0, 172640000.0, 172650000.0, 172660000.0, 172670000.0, 172680000.0, 172690000.0, 172700000.0, 172710000.0, 172720000.0, 172730000.0, 172740000.0, 172750000.0, 172760000.0, 172770000.0, 172780000.0, 172790000.0, 172800000.0, 172810000.0, 172820000.0, 172830000.0, 172840000.0, 172850000.0, 172860000.0, 172870000.0, 172880000.0, 172890000.0, 172900000.0, 172910000.0, 172920000.0, 172930000.0, 172940000.0, 172950000.0, 172960000.0, 172970000.0, 172980000.0, 172990000.0, 173000000.0, 173010000.0, 173020000.0, 173030000.0, 173040000.0, 173050000.0, 173060000.0, 173070000.0, 173080000.0, 173090000.0, 173100000.0, 173110000.0, 173120000.0, 173130000.0, 173140000.0, 173150000.0, 173160000.0, 173170000.0, 173180000.0, 173190000.0, 173200000.0, 173210000.0, 173220000.0, 173230000.0, 173240000.0, 173250000.0, 173260000.0, 173270000.0, 173280000.0, 173290000.0, 173300000.0, 173310000.0, 173320000.0, 173330000.0, 173340000.0, 173350000.0, 173360000.0, 173370000.0, 173380000.0, 173390000.0, 173400000.0, 173410000.0, 173420000.0, 173430000.0, 173440000.0, 173450000.0, 173460000.0, 173470000.0, 173480000.0, 173490000.0, 173500000.0, 173510000.0, 173520000.0, 173530000.0, 173540000.0, 173550000.0, 173560000.0, 173570000.0, 173580000.0, 173590000.0, 173600000.0, 173610000.0, 173620000.0, 173630000.0, 173640000.0, 173650000.0, 173660000.0, 173670000.0, 173680000.0, 173690000.0, 173700000.0, 173710000.0, 173720000.0, 173730000.0, 173740000.0, 173750000.0, 173760000.0, 173770000.0, 173780000.0, 173790000.0, 173800000.0, 173810000.0, 173820000.0, 173830000.0, 173840000.0, 173850000.0, 173860000.0, 173870000.0, 173880000.0, 173890000.0, 173900000.0, 173910000.0, 173920000.0, 173930000.0, 173940000.0, 173950000.0, 173960000.0, 173970000.0, 173980000.0, 173990000.0, 174000000.0, 174010000.0, 174020000.0, 174030000.0, 174040000.0, 174050000.0, 174060000.0, 174070000.0, 174080000.0, 174090000.0, 174100000.0, 174110000.0, 174120000.0, 174130000.0, 174140000.0, 174150000.0, 174160000.0, 174170000.0, 174180000.0, 174190000.0, 174200000.0, 174210000.0, 174220000.0, 174230000.0, 174240000.0, 174250000.0, 174260000.0, 174270000.0, 174280000.0, 174290000.0, 174300000.0, 174310000.0, 174320000.0, 174330000.0, 174340000.0, 174350000.0, 174360000.0, 174370000.0, 174380000.0, 174390000.0, 174400000.0, 174410000.0, 174420000.0, 174430000.0, 174440000.0, 174450000.0, 174460000.0, 174470000.0, 174480000.0, 174490000.0, 174500000.0, 174510000.0, 174520000.0, 174530000.0, 174540000.0, 174550000.0, 174560000.0, 174570000.0, 174580000.0, 174590000.0, 174600000.0, 174610000.0, 174620000.0, 174630000.0, 174640000.0, 174650000.0, 174660000.0, 174670000.0, 174680000.0, 174690000.0, 174700000.0, 174710000.0, 174720000.0, 174730000.0, 174740000.0, 174750000.0, 174760000.0, 174770000.0, 174780000.0, 174790000.0, 174800000.0, 174810000.0, 174820000.0, 174830000.0, 174840000.0, 174850000.0, 174860000.0, 174870000.0, 174880000.0, 174890000.0, 174900000.0, 174910000.0, 174920000.0, 174930000.0, 174940000.0, 174950000.0, 174960000.0, 174970000.0, 174980000.0, 174990000.0, 175000000.0, 175010000.0, 175020000.0, 175030000.0, 175040000.0, 175050000.0, 175060000.0, 175070000.0, 175080000.0, 175090000.0, 175100000.0, 175110000.0, 175120000.0, 175130000.0, 175140000.0, 175150000.0, 175160000.0, 175170000.0, 175180000.0, 175190000.0, 175200000.0, 175210000.0, 175220000.0, 175230000.0, 175240000.0, 175250000.0, 175260000.0, 175270000.0, 175280000.0, 175290000.0, 175300000.0, 175310000.0, 175320000.0, 175330000.0, 175340000.0, 175350000.0, 175360000.0, 175370000.0, 175380000.0, 175390000.0, 175400000.0, 175410000.0, 175420000.0, 175430000.0, 175440000.0, 175450000.0, 175460000.0, 175470000.0, 175480000.0, 175490000.0, 175500000.0, 175510000.0, 175520000.0, 175530000.0, 175540000.0, 175550000.0, 175560000.0, 175570000.0, 175580000.0, 175590000.0, 175600000.0, 175610000.0, 175620000.0, 175630000.0, 175640000.0, 175650000.0, 175660000.0, 175670000.0, 175680000.0, 175690000.0, 175700000.0, 175710000.0, 175720000.0, 175730000.0, 175740000.0, 175750000.0, 175760000.0, 175770000.0, 175780000.0, 175790000.0, 175800000.0, 175810000.0, 175820000.0, 175830000.0, 175840000.0, 175850000.0, 175860000.0, 175870000.0, 175880000.0, 175890000.0, 175900000.0, 175910000.0, 175920000.0, 175930000.0, 175940000.0, 175950000.0, 175960000.0, 175970000.0, 175980000.0, 175990000.0, 176000000.0, 176010000.0, 176020000.0, 176030000.0, 176040000.0, 176050000.0, 176060000.0, 176070000.0, 176080000.0, 176090000.0, 176100000.0, 176110000.0, 176120000.0, 176130000.0, 176140000.0, 176150000.0, 176160000.0, 176170000.0, 176180000.0, 176190000.0, 176200000.0, 176210000.0, 176220000.0, 176230000.0, 176240000.0, 176250000.0, 176260000.0, 176270000.0, 176280000.0, 176290000.0, 176300000.0, 176310000.0, 176320000.0, 176330000.0, 176340000.0, 176350000.0, 176360000.0, 176370000.0, 176380000.0, 176390000.0, 176400000.0, 176410000.0, 176420000.0, 176430000.0, 176440000.0, 176450000.0, 176460000.0, 176470000.0, 176480000.0, 176490000.0, 176500000.0, 176510000.0, 176520000.0, 176530000.0, 176540000.0, 176550000.0, 176560000.0, 176570000.0, 176580000.0, 176590000.0, 176600000.0, 176610000.0, 176620000.0, 176630000.0, 176640000.0, 176650000.0, 176660000.0, 176670000.0, 176680000.0, 176690000.0, 176700000.0, 176710000.0, 176720000.0, 176730000.0, 176740000.0, 176750000.0, 176760000.0, 176770000.0, 176780000.0, 176790000.0, 176800000.0, 176810000.0, 176820000.0, 176830000.0, 176840000.0, 176850000.0, 176860000.0, 176870000.0, 176880000.0, 176890000.0, 176900000.0, 176910000.0, 176920000.0, 176930000.0, 176940000.0, 176950000.0, 176960000.0, 176970000.0, 176980000.0, 176990000.0, 177000000.0, 177010000.0, 177020000.0, 177030000.0, 177040000.0, 177050000.0, 177060000.0, 177070000.0, 177080000.0, 177090000.0, 177100000.0, 177110000.0, 177120000.0, 177130000.0, 177140000.0, 177150000.0, 177160000.0, 177170000.0, 177180000.0, 177190000.0, 177200000.0, 177210000.0, 177220000.0, 177230000.0, 177240000.0, 177250000.0, 177260000.0, 177270000.0, 177280000.0, 177290000.0, 177300000.0, 177310000.0, 177320000.0, 177330000.0, 177340000.0, 177350000.0, 177360000.0, 177370000.0, 177380000.0, 177390000.0, 177400000.0, 177410000.0, 177420000.0, 177430000.0, 177440000.0, 177450000.0, 177460000.0, 177470000.0, 177480000.0, 177490000.0, 177500000.0, 177510000.0, 177520000.0, 177530000.0, 177540000.0, 177550000.0, 177560000.0, 177570000.0, 177580000.0, 177590000.0, 177600000.0, 177610000.0, 177620000.0, 177630000.0, 177640000.0, 177650000.0, 177660000.0, 177670000.0, 177680000.0, 177690000.0, 177700000.0, 177710000.0, 177720000.0, 177730000.0, 177740000.0, 177750000.0, 177760000.0, 177770000.0, 177780000.0, 177790000.0, 177800000.0, 177810000.0, 177820000.0, 177830000.0, 177840000.0, 177850000.0, 177860000.0, 177870000.0, 177880000.0, 177890000.0, 177900000.0, 177910000.0, 177920000.0, 177930000.0, 177940000.0, 177950000.0, 177960000.0, 177970000.0, 177980000.0, 177990000.0, 178000000.0, 178010000.0, 178020000.0, 178030000.0, 178040000.0, 178050000.0, 178060000.0, 178070000.0, 178080000.0, 178090000.0, 178100000.0, 178110000.0, 178120000.0, 178130000.0, 178140000.0, 178150000.0, 178160000.0, 178170000.0, 178180000.0, 178190000.0, 178200000.0, 178210000.0, 178220000.0, 178230000.0, 178240000.0, 178250000.0, 178260000.0, 178270000.0, 178280000.0, 178290000.0, 178300000.0, 178310000.0, 178320000.0, 178330000.0, 178340000.0, 178350000.0, 178360000.0, 178370000.0, 178380000.0, 178390000.0, 178400000.0, 178410000.0, 178420000.0, 178430000.0, 178440000.0, 178450000.0, 178460000.0, 178470000.0, 178480000.0, 178490000.0, 178500000.0, 178510000.0, 178520000.0, 178530000.0, 178540000.0, 178550000.0, 178560000.0, 178570000.0, 178580000.0, 178590000.0, 178600000.0, 178610000.0, 178620000.0, 178630000.0, 178640000.0, 178650000.0, 178660000.0, 178670000.0, 178680000.0, 178690000.0, 178700000.0, 178710000.0, 178720000.0, 178730000.0, 178740000.0, 178750000.0, 178760000.0, 178770000.0, 178780000.0, 178790000.0, 178800000.0, 178810000.0, 178820000.0, 178830000.0, 178840000.0, 178850000.0, 178860000.0, 178870000.0, 178880000.0, 178890000.0, 178900000.0, 178910000.0, 178920000.0, 178930000.0, 178940000.0, 178950000.0, 178960000.0, 178970000.0, 178980000.0, 178990000.0, 179000000.0, 179010000.0, 179020000.0, 179030000.0, 179040000.0, 179050000.0, 179060000.0, 179070000.0, 179080000.0, 179090000.0, 179100000.0, 179110000.0, 179120000.0, 179130000.0, 179140000.0, 179150000.0, 179160000.0, 179170000.0, 179180000.0, 179190000.0, 179200000.0, 179210000.0, 179220000.0, 179230000.0, 179240000.0, 179250000.0, 179260000.0, 179270000.0, 179280000.0, 179290000.0, 179300000.0, 179310000.0, 179320000.0, 179330000.0, 179340000.0, 179350000.0, 179360000.0, 179370000.0, 179380000.0, 179390000.0, 179400000.0, 179410000.0, 179420000.0, 179430000.0, 179440000.0, 179450000.0, 179460000.0, 179470000.0, 179480000.0, 179490000.0, 179500000.0, 179510000.0, 179520000.0, 179530000.0, 179540000.0, 179550000.0, 179560000.0, 179570000.0, 179580000.0, 179590000.0, 179600000.0, 179610000.0, 179620000.0, 179630000.0, 179640000.0, 179650000.0, 179660000.0, 179670000.0, 179680000.0, 179690000.0, 179700000.0, 179710000.0, 179720000.0, 179730000.0, 179740000.0, 179750000.0, 179760000.0, 179770000.0, 179780000.0, 179790000.0, 179800000.0, 179810000.0, 179820000.0, 179830000.0, 179840000.0, 179850000.0, 179860000.0, 179870000.0, 179880000.0, 179890000.0, 179900000.0, 179910000.0, 179920000.0, 179930000.0, 179940000.0, 179950000.0, 179960000.0, 179970000.0, 179980000.0, 179990000.0, 180000000.0, 180010000.0, 180020000.0, 180030000.0, 180040000.0, 180050000.0, 180060000.0, 180070000.0, 180080000.0, 180090000.0, 180100000.0, 180110000.0, 180120000.0, 180130000.0, 180140000.0, 180150000.0, 180160000.0, 180170000.0, 180180000.0, 180190000.0, 180200000.0, 180210000.0, 180220000.0, 180230000.0, 180240000.0, 180250000.0, 180260000.0, 180270000.0, 180280000.0, 180290000.0, 180300000.0, 180310000.0, 180320000.0, 180330000.0, 180340000.0, 180350000.0, 180360000.0, 180370000.0, 180380000.0, 180390000.0, 180400000.0, 180410000.0, 180420000.0, 180430000.0, 180440000.0, 180450000.0, 180460000.0, 180470000.0, 180480000.0, 180490000.0, 180500000.0, 180510000.0, 180520000.0, 180530000.0, 180540000.0, 180550000.0, 180560000.0, 180570000.0, 180580000.0, 180590000.0, 180600000.0, 180610000.0, 180620000.0, 180630000.0, 180640000.0, 180650000.0, 180660000.0, 180670000.0, 180680000.0, 180690000.0, 180700000.0, 180710000.0, 180720000.0, 180730000.0, 180740000.0, 180750000.0, 180760000.0, 180770000.0, 180780000.0, 180790000.0, 180800000.0, 180810000.0, 180820000.0, 180830000.0, 180840000.0, 180850000.0, 180860000.0, 180870000.0, 180880000.0, 180890000.0, 180900000.0, 180910000.0, 180920000.0, 180930000.0, 180940000.0, 180950000.0, 180960000.0, 180970000.0, 180980000.0, 180990000.0, 181000000.0, 181010000.0, 181020000.0, 181030000.0, 181040000.0, 181050000.0, 181060000.0, 181070000.0, 181080000.0, 181090000.0, 181100000.0, 181110000.0, 181120000.0, 181130000.0, 181140000.0, 181150000.0, 181160000.0, 181170000.0, 181180000.0, 181190000.0, 181200000.0, 181210000.0, 181220000.0, 181230000.0, 181240000.0, 181250000.0, 181260000.0, 181270000.0, 181280000.0, 181290000.0, 181300000.0, 181310000.0, 181320000.0, 181330000.0, 181340000.0, 181350000.0, 181360000.0, 181370000.0, 181380000.0, 181390000.0, 181400000.0, 181410000.0, 181420000.0, 181430000.0, 181440000.0, 181450000.0, 181460000.0, 181470000.0, 181480000.0, 181490000.0, 181500000.0, 181510000.0, 181520000.0, 181530000.0, 181540000.0, 181550000.0, 181560000.0, 181570000.0, 181580000.0, 181590000.0, 181600000.0, 181610000.0, 181620000.0, 181630000.0, 181640000.0, 181650000.0, 181660000.0, 181670000.0, 181680000.0, 181690000.0, 181700000.0, 181710000.0, 181720000.0, 181730000.0, 181740000.0, 181750000.0, 181760000.0, 181770000.0, 181780000.0, 181790000.0, 181800000.0, 181810000.0, 181820000.0, 181830000.0, 181840000.0, 181850000.0, 181860000.0, 181870000.0, 181880000.0, 181890000.0, 181900000.0, 181910000.0, 181920000.0, 181930000.0, 181940000.0, 181950000.0, 181960000.0, 181970000.0, 181980000.0, 181990000.0, 182000000.0, 182010000.0, 182020000.0, 182030000.0, 182040000.0, 182050000.0, 182060000.0, 182070000.0, 182080000.0, 182090000.0, 182100000.0, 182110000.0, 182120000.0, 182130000.0, 182140000.0, 182150000.0, 182160000.0, 182170000.0, 182180000.0, 182190000.0, 182200000.0, 182210000.0, 182220000.0, 182230000.0, 182240000.0, 182250000.0, 182260000.0, 182270000.0, 182280000.0, 182290000.0, 182300000.0, 182310000.0, 182320000.0, 182330000.0, 182340000.0, 182350000.0, 182360000.0, 182370000.0, 182380000.0, 182390000.0, 182400000.0, 182410000.0, 182420000.0, 182430000.0, 182440000.0, 182450000.0, 182460000.0, 182470000.0, 182480000.0, 182490000.0, 182500000.0, 182510000.0, 182520000.0, 182530000.0, 182540000.0, 182550000.0, 182560000.0, 182570000.0, 182580000.0, 182590000.0, 182600000.0, 182610000.0, 182620000.0, 182630000.0, 182640000.0, 182650000.0, 182660000.0, 182670000.0, 182680000.0, 182690000.0, 182700000.0, 182710000.0, 182720000.0, 182730000.0, 182740000.0, 182750000.0, 182760000.0, 182770000.0, 182780000.0, 182790000.0, 182800000.0, 182810000.0, 182820000.0, 182830000.0, 182840000.0, 182850000.0, 182860000.0, 182870000.0, 182880000.0, 182890000.0, 182900000.0, 182910000.0, 182920000.0, 182930000.0, 182940000.0, 182950000.0, 182960000.0, 182970000.0, 182980000.0, 182990000.0, 183000000.0, 183010000.0, 183020000.0, 183030000.0, 183040000.0, 183050000.0, 183060000.0, 183070000.0, 183080000.0, 183090000.0, 183100000.0, 183110000.0, 183120000.0, 183130000.0, 183140000.0, 183150000.0, 183160000.0, 183170000.0, 183180000.0, 183190000.0, 183200000.0, 183210000.0, 183220000.0, 183230000.0, 183240000.0, 183250000.0, 183260000.0, 183270000.0, 183280000.0, 183290000.0, 183300000.0, 183310000.0, 183320000.0, 183330000.0, 183340000.0, 183350000.0, 183360000.0, 183370000.0, 183380000.0, 183390000.0, 183400000.0, 183410000.0, 183420000.0, 183430000.0, 183440000.0, 183450000.0, 183460000.0, 183470000.0, 183480000.0, 183490000.0, 183500000.0, 183510000.0, 183520000.0, 183530000.0, 183540000.0, 183550000.0, 183560000.0, 183570000.0, 183580000.0, 183590000.0, 183600000.0, 183610000.0, 183620000.0, 183630000.0, 183640000.0, 183650000.0, 183660000.0, 183670000.0, 183680000.0, 183690000.0, 183700000.0, 183710000.0, 183720000.0, 183730000.0, 183740000.0, 183750000.0, 183760000.0, 183770000.0, 183780000.0, 183790000.0, 183800000.0, 183810000.0, 183820000.0, 183830000.0, 183840000.0, 183850000.0, 183860000.0, 183870000.0, 183880000.0, 183890000.0, 183900000.0, 183910000.0, 183920000.0, 183930000.0, 183940000.0, 183950000.0, 183960000.0, 183970000.0, 183980000.0, 183990000.0, 184000000.0, 184010000.0, 184020000.0, 184030000.0, 184040000.0, 184050000.0, 184060000.0, 184070000.0, 184080000.0, 184090000.0, 184100000.0, 184110000.0, 184120000.0, 184130000.0, 184140000.0, 184150000.0, 184160000.0, 184170000.0, 184180000.0, 184190000.0, 184200000.0, 184210000.0, 184220000.0, 184230000.0, 184240000.0, 184250000.0, 184260000.0, 184270000.0, 184280000.0, 184290000.0, 184300000.0, 184310000.0, 184320000.0, 184330000.0, 184340000.0, 184350000.0, 184360000.0, 184370000.0, 184380000.0, 184390000.0, 184400000.0, 184410000.0, 184420000.0, 184430000.0, 184440000.0, 184450000.0, 184460000.0, 184470000.0, 184480000.0, 184490000.0, 184500000.0, 184510000.0, 184520000.0, 184530000.0, 184540000.0, 184550000.0, 184560000.0, 184570000.0, 184580000.0, 184590000.0, 184600000.0, 184610000.0, 184620000.0, 184630000.0, 184640000.0, 184650000.0, 184660000.0, 184670000.0, 184680000.0, 184690000.0, 184700000.0, 184710000.0, 184720000.0, 184730000.0, 184740000.0, 184750000.0, 184760000.0, 184770000.0, 184780000.0, 184790000.0, 184800000.0, 184810000.0, 184820000.0, 184830000.0, 184840000.0, 184850000.0, 184860000.0, 184870000.0, 184880000.0, 184890000.0, 184900000.0, 184910000.0, 184920000.0, 184930000.0, 184940000.0, 184950000.0, 184960000.0, 184970000.0, 184980000.0, 184990000.0, 185000000.0, 185010000.0, 185020000.0, 185030000.0, 185040000.0, 185050000.0, 185060000.0, 185070000.0, 185080000.0, 185090000.0, 185100000.0, 185110000.0, 185120000.0, 185130000.0, 185140000.0, 185150000.0, 185160000.0, 185170000.0, 185180000.0, 185190000.0, 185200000.0, 185210000.0, 185220000.0, 185230000.0, 185240000.0, 185250000.0, 185260000.0, 185270000.0, 185280000.0, 185290000.0, 185300000.0, 185310000.0, 185320000.0, 185330000.0, 185340000.0, 185350000.0, 185360000.0, 185370000.0, 185380000.0, 185390000.0, 185400000.0, 185410000.0, 185420000.0, 185430000.0, 185440000.0, 185450000.0, 185460000.0, 185470000.0, 185480000.0, 185490000.0, 185500000.0, 185510000.0, 185520000.0, 185530000.0, 185540000.0, 185550000.0, 185560000.0, 185570000.0, 185580000.0, 185590000.0, 185600000.0, 185610000.0, 185620000.0, 185630000.0, 185640000.0, 185650000.0, 185660000.0, 185670000.0, 185680000.0, 185690000.0, 185700000.0, 185710000.0, 185720000.0, 185730000.0, 185740000.0, 185750000.0, 185760000.0, 185770000.0, 185780000.0, 185790000.0, 185800000.0, 185810000.0, 185820000.0, 185830000.0, 185840000.0, 185850000.0, 185860000.0, 185870000.0, 185880000.0, 185890000.0, 185900000.0, 185910000.0, 185920000.0, 185930000.0, 185940000.0, 185950000.0, 185960000.0, 185970000.0, 185980000.0, 185990000.0, 186000000.0, 186010000.0, 186020000.0, 186030000.0, 186040000.0, 186050000.0, 186060000.0, 186070000.0, 186080000.0, 186090000.0, 186100000.0, 186110000.0, 186120000.0, 186130000.0, 186140000.0, 186150000.0, 186160000.0, 186170000.0, 186180000.0, 186190000.0, 186200000.0, 186210000.0, 186220000.0, 186230000.0, 186240000.0, 186250000.0, 186260000.0, 186270000.0, 186280000.0, 186290000.0, 186300000.0, 186310000.0, 186320000.0, 186330000.0, 186340000.0, 186350000.0, 186360000.0, 186370000.0, 186380000.0, 186390000.0, 186400000.0, 186410000.0, 186420000.0, 186430000.0, 186440000.0, 186450000.0, 186460000.0, 186470000.0, 186480000.0, 186490000.0, 186500000.0, 186510000.0, 186520000.0, 186530000.0, 186540000.0, 186550000.0, 186560000.0, 186570000.0, 186580000.0, 186590000.0, 186600000.0, 186610000.0, 186620000.0, 186630000.0, 186640000.0, 186650000.0, 186660000.0, 186670000.0, 186680000.0, 186690000.0, 186700000.0, 186710000.0, 186720000.0, 186730000.0, 186740000.0, 186750000.0, 186760000.0, 186770000.0, 186780000.0, 186790000.0, 186800000.0, 186810000.0, 186820000.0, 186830000.0, 186840000.0, 186850000.0, 186860000.0, 186870000.0, 186880000.0, 186890000.0, 186900000.0, 186910000.0, 186920000.0, 186930000.0, 186940000.0, 186950000.0, 186960000.0, 186970000.0, 186980000.0, 186990000.0, 187000000.0, 187010000.0, 187020000.0, 187030000.0, 187040000.0, 187050000.0, 187060000.0, 187070000.0, 187080000.0, 187090000.0, 187100000.0, 187110000.0, 187120000.0, 187130000.0, 187140000.0, 187150000.0, 187160000.0, 187170000.0, 187180000.0, 187190000.0, 187200000.0, 187210000.0, 187220000.0, 187230000.0, 187240000.0, 187250000.0, 187260000.0, 187270000.0, 187280000.0, 187290000.0, 187300000.0, 187310000.0, 187320000.0, 187330000.0, 187340000.0, 187350000.0, 187360000.0, 187370000.0, 187380000.0, 187390000.0, 187400000.0, 187410000.0, 187420000.0, 187430000.0, 187440000.0, 187450000.0, 187460000.0, 187470000.0, 187480000.0, 187490000.0, 187500000.0, 187510000.0, 187520000.0, 187530000.0, 187540000.0, 187550000.0, 187560000.0, 187570000.0, 187580000.0, 187590000.0, 187600000.0, 187610000.0, 187620000.0, 187630000.0, 187640000.0, 187650000.0, 187660000.0, 187670000.0, 187680000.0, 187690000.0, 187700000.0, 187710000.0, 187720000.0, 187730000.0, 187740000.0, 187750000.0, 187760000.0, 187770000.0, 187780000.0, 187790000.0, 187800000.0, 187810000.0, 187820000.0, 187830000.0, 187840000.0, 187850000.0, 187860000.0, 187870000.0, 187880000.0, 187890000.0, 187900000.0, 187910000.0, 187920000.0, 187930000.0, 187940000.0, 187950000.0, 187960000.0, 187970000.0, 187980000.0, 187990000.0, 188000000.0, 188010000.0, 188020000.0, 188030000.0, 188040000.0, 188050000.0, 188060000.0, 188070000.0, 188080000.0, 188090000.0, 188100000.0, 188110000.0, 188120000.0, 188130000.0, 188140000.0, 188150000.0, 188160000.0, 188170000.0, 188180000.0, 188190000.0, 188200000.0, 188210000.0, 188220000.0, 188230000.0, 188240000.0, 188250000.0, 188260000.0, 188270000.0, 188280000.0, 188290000.0, 188300000.0, 188310000.0, 188320000.0, 188330000.0, 188340000.0, 188350000.0, 188360000.0, 188370000.0, 188380000.0, 188390000.0, 188400000.0, 188410000.0, 188420000.0, 188430000.0, 188440000.0, 188450000.0, 188460000.0, 188470000.0, 188480000.0, 188490000.0, 188500000.0, 188510000.0, 188520000.0, 188530000.0, 188540000.0, 188550000.0, 188560000.0, 188570000.0, 188580000.0, 188590000.0, 188600000.0, 188610000.0, 188620000.0, 188630000.0, 188640000.0, 188650000.0, 188660000.0, 188670000.0, 188680000.0, 188690000.0, 188700000.0, 188710000.0, 188720000.0, 188730000.0, 188740000.0, 188750000.0, 188760000.0, 188770000.0, 188780000.0, 188790000.0, 188800000.0, 188810000.0, 188820000.0, 188830000.0, 188840000.0, 188850000.0, 188860000.0, 188870000.0, 188880000.0, 188890000.0, 188900000.0, 188910000.0, 188920000.0, 188930000.0, 188940000.0, 188950000.0, 188960000.0, 188970000.0, 188980000.0, 188990000.0, 189000000.0, 189010000.0, 189020000.0, 189030000.0, 189040000.0, 189050000.0, 189060000.0, 189070000.0, 189080000.0, 189090000.0, 189100000.0, 189110000.0, 189120000.0, 189130000.0, 189140000.0, 189150000.0, 189160000.0, 189170000.0, 189180000.0, 189190000.0, 189200000.0, 189210000.0, 189220000.0, 189230000.0, 189240000.0, 189250000.0, 189260000.0, 189270000.0, 189280000.0, 189290000.0, 189300000.0, 189310000.0, 189320000.0, 189330000.0, 189340000.0, 189350000.0, 189360000.0, 189370000.0, 189380000.0, 189390000.0, 189400000.0, 189410000.0, 189420000.0, 189430000.0, 189440000.0, 189450000.0, 189460000.0, 189470000.0, 189480000.0, 189490000.0, 189500000.0, 189510000.0, 189520000.0, 189530000.0, 189540000.0, 189550000.0, 189560000.0, 189570000.0, 189580000.0, 189590000.0, 189600000.0, 189610000.0, 189620000.0, 189630000.0, 189640000.0, 189650000.0, 189660000.0, 189670000.0, 189680000.0, 189690000.0, 189700000.0, 189710000.0, 189720000.0, 189730000.0, 189740000.0, 189750000.0, 189760000.0, 189770000.0, 189780000.0, 189790000.0, 189800000.0, 189810000.0, 189820000.0, 189830000.0, 189840000.0, 189850000.0, 189860000.0, 189870000.0, 189880000.0, 189890000.0, 189900000.0, 189910000.0, 189920000.0, 189930000.0, 189940000.0, 189950000.0, 189960000.0, 189970000.0, 189980000.0, 189990000.0, 190000000.0, 190010000.0, 190020000.0, 190030000.0, 190040000.0, 190050000.0, 190060000.0, 190070000.0, 190080000.0, 190090000.0, 190100000.0, 190110000.0, 190120000.0, 190130000.0, 190140000.0, 190150000.0, 190160000.0, 190170000.0, 190180000.0, 190190000.0, 190200000.0, 190210000.0, 190220000.0, 190230000.0, 190240000.0, 190250000.0, 190260000.0, 190270000.0, 190280000.0, 190290000.0, 190300000.0, 190310000.0, 190320000.0, 190330000.0, 190340000.0, 190350000.0, 190360000.0, 190370000.0, 190380000.0, 190390000.0, 190400000.0, 190410000.0, 190420000.0, 190430000.0, 190440000.0, 190450000.0, 190460000.0, 190470000.0, 190480000.0, 190490000.0, 190500000.0, 190510000.0, 190520000.0, 190530000.0, 190540000.0, 190550000.0, 190560000.0, 190570000.0, 190580000.0, 190590000.0, 190600000.0, 190610000.0, 190620000.0, 190630000.0, 190640000.0, 190650000.0, 190660000.0, 190670000.0, 190680000.0, 190690000.0, 190700000.0, 190710000.0, 190720000.0, 190730000.0, 190740000.0, 190750000.0, 190760000.0, 190770000.0, 190780000.0, 190790000.0, 190800000.0, 190810000.0, 190820000.0, 190830000.0, 190840000.0, 190850000.0, 190860000.0, 190870000.0, 190880000.0, 190890000.0, 190900000.0, 190910000.0, 190920000.0, 190930000.0, 190940000.0, 190950000.0, 190960000.0, 190970000.0, 190980000.0, 190990000.0, 191000000.0, 191010000.0, 191020000.0, 191030000.0, 191040000.0, 191050000.0, 191060000.0, 191070000.0, 191080000.0, 191090000.0, 191100000.0, 191110000.0, 191120000.0, 191130000.0, 191140000.0, 191150000.0, 191160000.0, 191170000.0, 191180000.0, 191190000.0, 191200000.0, 191210000.0, 191220000.0, 191230000.0, 191240000.0, 191250000.0, 191260000.0, 191270000.0, 191280000.0, 191290000.0, 191300000.0, 191310000.0, 191320000.0, 191330000.0, 191340000.0, 191350000.0, 191360000.0, 191370000.0, 191380000.0, 191390000.0, 191400000.0, 191410000.0, 191420000.0, 191430000.0, 191440000.0, 191450000.0, 191460000.0, 191470000.0, 191480000.0, 191490000.0, 191500000.0, 191510000.0, 191520000.0, 191530000.0, 191540000.0, 191550000.0, 191560000.0, 191570000.0, 191580000.0, 191590000.0, 191600000.0, 191610000.0, 191620000.0, 191630000.0, 191640000.0, 191650000.0, 191660000.0, 191670000.0, 191680000.0, 191690000.0, 191700000.0, 191710000.0, 191720000.0, 191730000.0, 191740000.0, 191750000.0, 191760000.0, 191770000.0, 191780000.0, 191790000.0, 191800000.0, 191810000.0, 191820000.0, 191830000.0, 191840000.0, 191850000.0, 191860000.0, 191870000.0, 191880000.0, 191890000.0, 191900000.0, 191910000.0, 191920000.0, 191930000.0, 191940000.0, 191950000.0, 191960000.0, 191970000.0, 191980000.0, 191990000.0, 192000000.0, 192010000.0, 192020000.0, 192030000.0, 192040000.0, 192050000.0, 192060000.0, 192070000.0, 192080000.0, 192090000.0, 192100000.0, 192110000.0, 192120000.0, 192130000.0, 192140000.0, 192150000.0, 192160000.0, 192170000.0, 192180000.0, 192190000.0, 192200000.0, 192210000.0, 192220000.0, 192230000.0, 192240000.0, 192250000.0, 192260000.0, 192270000.0, 192280000.0, 192290000.0, 192300000.0, 192310000.0, 192320000.0, 192330000.0, 192340000.0, 192350000.0, 192360000.0, 192370000.0, 192380000.0, 192390000.0, 192400000.0, 192410000.0, 192420000.0, 192430000.0, 192440000.0, 192450000.0, 192460000.0, 192470000.0, 192480000.0, 192490000.0, 192500000.0, 192510000.0, 192520000.0, 192530000.0, 192540000.0, 192550000.0, 192560000.0, 192570000.0, 192580000.0, 192590000.0, 192600000.0, 192610000.0, 192620000.0, 192630000.0, 192640000.0, 192650000.0, 192660000.0, 192670000.0, 192680000.0, 192690000.0, 192700000.0, 192710000.0, 192720000.0, 192730000.0, 192740000.0, 192750000.0, 192760000.0, 192770000.0, 192780000.0, 192790000.0, 192800000.0, 192810000.0, 192820000.0, 192830000.0, 192840000.0, 192850000.0, 192860000.0, 192870000.0, 192880000.0, 192890000.0, 192900000.0, 192910000.0, 192920000.0, 192930000.0, 192940000.0, 192950000.0, 192960000.0, 192970000.0, 192980000.0, 192990000.0, 193000000.0, 193010000.0, 193020000.0, 193030000.0, 193040000.0, 193050000.0, 193060000.0, 193070000.0, 193080000.0, 193090000.0, 193100000.0, 193110000.0, 193120000.0, 193130000.0, 193140000.0, 193150000.0, 193160000.0, 193170000.0, 193180000.0, 193190000.0, 193200000.0, 193210000.0, 193220000.0, 193230000.0, 193240000.0, 193250000.0, 193260000.0, 193270000.0, 193280000.0, 193290000.0, 193300000.0, 193310000.0, 193320000.0, 193330000.0, 193340000.0, 193350000.0, 193360000.0, 193370000.0, 193380000.0, 193390000.0, 193400000.0, 193410000.0, 193420000.0, 193430000.0, 193440000.0, 193450000.0, 193460000.0, 193470000.0, 193480000.0, 193490000.0, 193500000.0, 193510000.0, 193520000.0, 193530000.0, 193540000.0, 193550000.0, 193560000.0, 193570000.0, 193580000.0, 193590000.0, 193600000.0, 193610000.0, 193620000.0, 193630000.0, 193640000.0, 193650000.0, 193660000.0, 193670000.0, 193680000.0, 193690000.0, 193700000.0, 193710000.0, 193720000.0, 193730000.0, 193740000.0, 193750000.0, 193760000.0, 193770000.0, 193780000.0, 193790000.0, 193800000.0, 193810000.0, 193820000.0, 193830000.0, 193840000.0, 193850000.0, 193860000.0, 193870000.0, 193880000.0, 193890000.0, 193900000.0, 193910000.0, 193920000.0, 193930000.0, 193940000.0, 193950000.0, 193960000.0, 193970000.0, 193980000.0, 193990000.0, 194000000.0, 194010000.0, 194020000.0, 194030000.0, 194040000.0, 194050000.0, 194060000.0, 194070000.0, 194080000.0, 194090000.0, 194100000.0, 194110000.0, 194120000.0, 194130000.0, 194140000.0, 194150000.0, 194160000.0, 194170000.0, 194180000.0, 194190000.0, 194200000.0, 194210000.0, 194220000.0, 194230000.0, 194240000.0, 194250000.0, 194260000.0, 194270000.0, 194280000.0, 194290000.0, 194300000.0, 194310000.0, 194320000.0, 194330000.0, 194340000.0, 194350000.0, 194360000.0, 194370000.0, 194380000.0, 194390000.0, 194400000.0, 194410000.0, 194420000.0, 194430000.0, 194440000.0, 194450000.0, 194460000.0, 194470000.0, 194480000.0, 194490000.0, 194500000.0, 194510000.0, 194520000.0, 194530000.0, 194540000.0, 194550000.0, 194560000.0, 194570000.0, 194580000.0, 194590000.0, 194600000.0, 194610000.0, 194620000.0, 194630000.0, 194640000.0, 194650000.0, 194660000.0, 194670000.0, 194680000.0, 194690000.0, 194700000.0, 194710000.0, 194720000.0, 194730000.0, 194740000.0, 194750000.0, 194760000.0, 194770000.0, 194780000.0, 194790000.0, 194800000.0, 194810000.0, 194820000.0, 194830000.0, 194840000.0, 194850000.0, 194860000.0, 194870000.0, 194880000.0, 194890000.0, 194900000.0, 194910000.0, 194920000.0, 194930000.0, 194940000.0, 194950000.0, 194960000.0, 194970000.0, 194980000.0, 194990000.0, 195000000.0, 195010000.0, 195020000.0, 195030000.0, 195040000.0, 195050000.0, 195060000.0, 195070000.0, 195080000.0, 195090000.0, 195100000.0, 195110000.0, 195120000.0, 195130000.0, 195140000.0, 195150000.0, 195160000.0, 195170000.0, 195180000.0, 195190000.0, 195200000.0, 195210000.0, 195220000.0, 195230000.0, 195240000.0, 195250000.0, 195260000.0, 195270000.0, 195280000.0, 195290000.0, 195300000.0, 195310000.0, 195320000.0, 195330000.0, 195340000.0, 195350000.0, 195360000.0, 195370000.0, 195380000.0, 195390000.0, 195400000.0, 195410000.0, 195420000.0, 195430000.0, 195440000.0, 195450000.0, 195460000.0, 195470000.0, 195480000.0, 195490000.0, 195500000.0, 195510000.0, 195520000.0, 195530000.0, 195540000.0, 195550000.0, 195560000.0, 195570000.0, 195580000.0, 195590000.0, 195600000.0, 195610000.0, 195620000.0, 195630000.0, 195640000.0, 195650000.0, 195660000.0, 195670000.0, 195680000.0, 195690000.0, 195700000.0, 195710000.0, 195720000.0, 195730000.0, 195740000.0, 195750000.0, 195760000.0, 195770000.0, 195780000.0, 195790000.0, 195800000.0, 195810000.0, 195820000.0, 195830000.0, 195840000.0, 195850000.0, 195860000.0, 195870000.0, 195880000.0, 195890000.0, 195900000.0, 195910000.0, 195920000.0, 195930000.0, 195940000.0, 195950000.0, 195960000.0, 195970000.0, 195980000.0, 195990000.0, 196000000.0, 196010000.0, 196020000.0, 196030000.0, 196040000.0, 196050000.0, 196060000.0, 196070000.0, 196080000.0, 196090000.0, 196100000.0, 196110000.0, 196120000.0, 196130000.0, 196140000.0, 196150000.0, 196160000.0, 196170000.0, 196180000.0, 196190000.0, 196200000.0, 196210000.0, 196220000.0, 196230000.0, 196240000.0, 196250000.0, 196260000.0, 196270000.0, 196280000.0, 196290000.0, 196300000.0, 196310000.0, 196320000.0, 196330000.0, 196340000.0, 196350000.0, 196360000.0, 196370000.0, 196380000.0, 196390000.0, 196400000.0, 196410000.0, 196420000.0, 196430000.0, 196440000.0, 196450000.0, 196460000.0, 196470000.0, 196480000.0, 196490000.0, 196500000.0, 196510000.0, 196520000.0, 196530000.0, 196540000.0, 196550000.0, 196560000.0, 196570000.0, 196580000.0, 196590000.0, 196600000.0, 196610000.0, 196620000.0, 196630000.0, 196640000.0, 196650000.0, 196660000.0, 196670000.0, 196680000.0, 196690000.0, 196700000.0, 196710000.0, 196720000.0, 196730000.0, 196740000.0, 196750000.0, 196760000.0, 196770000.0, 196780000.0, 196790000.0, 196800000.0, 196810000.0, 196820000.0, 196830000.0, 196840000.0, 196850000.0, 196860000.0, 196870000.0, 196880000.0, 196890000.0, 196900000.0, 196910000.0, 196920000.0, 196930000.0, 196940000.0, 196950000.0, 196960000.0, 196970000.0, 196980000.0, 196990000.0, 197000000.0, 197010000.0, 197020000.0, 197030000.0, 197040000.0, 197050000.0, 197060000.0, 197070000.0, 197080000.0, 197090000.0, 197100000.0, 197110000.0, 197120000.0, 197130000.0, 197140000.0, 197150000.0, 197160000.0, 197170000.0, 197180000.0, 197190000.0, 197200000.0, 197210000.0, 197220000.0, 197230000.0, 197240000.0, 197250000.0, 197260000.0, 197270000.0, 197280000.0, 197290000.0, 197300000.0, 197310000.0, 197320000.0, 197330000.0, 197340000.0, 197350000.0, 197360000.0, 197370000.0, 197380000.0, 197390000.0, 197400000.0, 197410000.0, 197420000.0, 197430000.0, 197440000.0, 197450000.0, 197460000.0, 197470000.0, 197480000.0, 197490000.0, 197500000.0, 197510000.0, 197520000.0, 197530000.0, 197540000.0, 197550000.0, 197560000.0, 197570000.0, 197580000.0, 197590000.0, 197600000.0, 197610000.0, 197620000.0, 197630000.0, 197640000.0, 197650000.0, 197660000.0, 197670000.0, 197680000.0, 197690000.0, 197700000.0, 197710000.0, 197720000.0, 197730000.0, 197740000.0, 197750000.0, 197760000.0, 197770000.0, 197780000.0, 197790000.0, 197800000.0, 197810000.0, 197820000.0, 197830000.0, 197840000.0, 197850000.0, 197860000.0, 197870000.0, 197880000.0, 197890000.0, 197900000.0, 197910000.0, 197920000.0, 197930000.0, 197940000.0, 197950000.0, 197960000.0, 197970000.0, 197980000.0, 197990000.0, 198000000.0, 198010000.0, 198020000.0, 198030000.0, 198040000.0, 198050000.0, 198060000.0, 198070000.0, 198080000.0, 198090000.0, 198100000.0, 198110000.0, 198120000.0, 198130000.0, 198140000.0, 198150000.0, 198160000.0, 198170000.0, 198180000.0, 198190000.0, 198200000.0, 198210000.0, 198220000.0, 198230000.0, 198240000.0, 198250000.0, 198260000.0, 198270000.0, 198280000.0, 198290000.0, 198300000.0, 198310000.0, 198320000.0, 198330000.0, 198340000.0, 198350000.0, 198360000.0, 198370000.0, 198380000.0, 198390000.0, 198400000.0, 198410000.0, 198420000.0, 198430000.0, 198440000.0, 198450000.0, 198460000.0, 198470000.0, 198480000.0, 198490000.0, 198500000.0, 198510000.0, 198520000.0, 198530000.0, 198540000.0, 198550000.0, 198560000.0, 198570000.0, 198580000.0, 198590000.0, 198600000.0, 198610000.0, 198620000.0, 198630000.0, 198640000.0, 198650000.0, 198660000.0, 198670000.0, 198680000.0, 198690000.0, 198700000.0, 198710000.0, 198720000.0, 198730000.0, 198740000.0, 198750000.0, 198760000.0, 198770000.0, 198780000.0, 198790000.0, 198800000.0, 198810000.0, 198820000.0, 198830000.0, 198840000.0, 198850000.0, 198860000.0, 198870000.0, 198880000.0, 198890000.0, 198900000.0, 198910000.0, 198920000.0, 198930000.0, 198940000.0, 198950000.0, 198960000.0, 198970000.0, 198980000.0, 198990000.0, 199000000.0, 199010000.0, 199020000.0, 199030000.0, 199040000.0, 199050000.0, 199060000.0, 199070000.0, 199080000.0, 199090000.0, 199100000.0, 199110000.0, 199120000.0, 199130000.0, 199140000.0, 199150000.0, 199160000.0, 199170000.0, 199180000.0, 199190000.0, 199200000.0, 199210000.0, 199220000.0, 199230000.0, 199240000.0, 199250000.0, 199260000.0, 199270000.0, 199280000.0, 199290000.0, 199300000.0, 199310000.0, 199320000.0, 199330000.0, 199340000.0, 199350000.0, 199360000.0, 199370000.0, 199380000.0, 199390000.0, 199400000.0, 199410000.0, 199420000.0, 199430000.0, 199440000.0, 199450000.0, 199460000.0, 199470000.0, 199480000.0, 199490000.0, 199500000.0, 199510000.0, 199520000.0, 199530000.0, 199540000.0, 199550000.0, 199560000.0, 199570000.0, 199580000.0, 199590000.0, 199600000.0, 199610000.0, 199620000.0, 199630000.0, 199640000.0, 199650000.0, 199660000.0, 199670000.0, 199680000.0, 199690000.0, 199700000.0, 199710000.0, 199720000.0, 199730000.0, 199740000.0, 199750000.0, 199760000.0, 199770000.0, 199780000.0, 199790000.0, 199800000.0, 199810000.0, 199820000.0, 199830000.0, 199840000.0, 199850000.0, 199860000.0, 199870000.0, 199880000.0, 199890000.0, 199900000.0, 199910000.0, 199920000.0, 199930000.0, 199940000.0, 199950000.0, 199960000.0, 199970000.0, 199980000.0, 199990000.0, 200000000.0, 200010000.0, 200020000.0, 200030000.0, 200040000.0, 200050000.0, 200060000.0, 200070000.0, 200080000.0, 200090000.0, 200100000.0, 200110000.0, 200120000.0, 200130000.0, 200140000.0, 200150000.0, 200160000.0, 200170000.0, 200180000.0, 200190000.0, 200200000.0, 200210000.0, 200220000.0, 200230000.0, 200240000.0, 200250000.0, 200260000.0, 200270000.0, 200280000.0, 200290000.0, 200300000.0, 200310000.0]) },

        R.A. (tile_pointing):     330.8820584879053 degrees,
        Dec. (tile_pointing):     -14.70067817554409 degrees,
        R.A. (phase center):      Some(Some(333.607)) degrees,
        Dec. (phase center):      Some(Some(-17.0267)) degrees,
        Azimuth:                  78.6901 degrees,
        Altitude:                 52.8056 degrees,
        Sun altitude:             -11.8554687635075 degrees,
        Sun distance:             64.7231885364707 degrees,
        Moon distance:            173.574994053323 degrees,
        Jupiter distance:         2.2906025553914 degrees,
        LST:                      293.641528461053 degrees,
        Hour angle:                02:31:05.95 degrees,
        Grid name:                sweet,
        Grid number:              82,

        num antennas:             128,
        antennas:                 [Tile011, Tile012, Tile013, Tile014, Tile015, Tile016, Tile017, Tile018, Tile021, Tile022, Tile023, Tile024, Tile025, Tile026, Tile027, Tile028, Tile031, Tile032, Tile033, Tile034, Tile035, Tile036, Tile037, Tile038, Tile041, Tile042, Tile043, Tile044, Tile045, Tile046, Tile047, Tile048, Tile061, Tile062, Tile063, Tile064, Tile065, Tile066, Tile067, Tile068, Tile081, Tile082, Tile083, Tile084, Tile085, Tile086, Tile087, Tile088, Tile091, Tile092, Tile093, Tile094, Tile095, Tile096, Tile097, Tile098, RFIpole, HexE1, HexE2, HexE3, HexE4, HexE5, HexE6, HexE7, HexE8, HexE9, HexE10, HexE11, HexE12, HexE13, HexE14, HexE15, HexE16, HexE17, HexE18, HexE19, HexE20, HexE21, HexE22, HexE23, HexE24, HexE25, HexE26, HexE27, HexE28, HexE29, HexE30, HexE31, HexE32, HexE33, HexE34, HexE35, HexE36, HexS1, HexS2, HexS3, HexS4, HexS5, HexS6, HexS7, HexS8, HexS9, HexS10, HexS11, HexS12, HexS13, HexS14, HexS15, HexS16, HexS17, HexS18, HexS19, HexS20, HexS21, HexS22, HexS23, HexS24, HexS25, HexS26, HexS27, HexS28, HexS29, HexS30, HexS32, HexS33, HexS34, HexS35, HexS36],
        rf_inputs:                [Tile011X, Tile011Y, Tile012X, Tile012Y, Tile013X, Tile013Y, Tile014X, Tile014Y, Tile015X, Tile015Y, Tile016X, Tile016Y, Tile017X, Tile017Y, Tile018X, Tile018Y, Tile021X, Tile021Y, Tile022X, Tile022Y, Tile023X, Tile023Y, Tile024X, Tile024Y, Tile025X, Tile025Y, Tile026X, Tile026Y, Tile027X, Tile027Y, Tile028X, Tile028Y, Tile031X, Tile031Y, Tile032X, Tile032Y, Tile033X, Tile033Y, Tile034X, Tile034Y, Tile035X, Tile035Y, Tile036X, Tile036Y, Tile037X, Tile037Y, Tile038X, Tile038Y, Tile041X, Tile041Y, Tile042X, Tile042Y, Tile043X, Tile043Y, Tile044X, Tile044Y, Tile045X, Tile045Y, Tile046X, Tile046Y, Tile047X, Tile047Y, Tile048X, Tile048Y, Tile061X, Tile061Y, Tile062X, Tile062Y, Tile063X, Tile063Y, Tile064X, Tile064Y, Tile065X, Tile065Y, Tile066X, Tile066Y, Tile067X, Tile067Y, Tile068X, Tile068Y, Tile081X, Tile081Y, Tile082X, Tile082Y, Tile083X, Tile083Y, Tile084X, Tile084Y, Tile085X, Tile085Y, Tile086X, Tile086Y, Tile087X, Tile087Y, Tile088X, Tile088Y, Tile091X, Tile091Y, Tile092X, Tile092Y, Tile093X, Tile093Y, Tile094X, Tile094Y, Tile095X, Tile095Y, Tile096X, Tile096Y, Tile097X, Tile097Y, Tile098X, Tile098Y, RFIpoleX, RFIpoleY, HexE1X, HexE1Y, HexE2X, HexE2Y, HexE3X, HexE3Y, HexE4X, HexE4Y, HexE5X, HexE5Y, HexE6X, HexE6Y, HexE7X, HexE7Y, HexE8X, HexE8Y, HexE9X, HexE9Y, HexE10X, HexE10Y, HexE11X, HexE11Y, HexE12X, HexE12Y, HexE13X, HexE13Y, HexE14X, HexE14Y, HexE15X, HexE15Y, HexE16X, HexE16Y, HexE17X, HexE17Y, HexE18X, HexE18Y, HexE19X, HexE19Y, HexE20X, HexE20Y, HexE21X, HexE21Y, HexE22X, HexE22Y, HexE23X, HexE23Y, HexE24X, HexE24Y, HexE25X, HexE25Y, HexE26X, HexE26Y, HexE27X, HexE27Y, HexE28X, HexE28Y, HexE29X, HexE29Y, HexE30X, HexE30Y, HexE31X, HexE31Y, HexE32X, HexE32Y, HexE33X, HexE33Y, HexE34X, HexE34Y, HexE35X, HexE35Y, HexE36X, HexE36Y, HexS1X, HexS1Y, HexS2X, HexS2Y, HexS3X, HexS3Y, HexS4X, HexS4Y, HexS5X, HexS5Y, HexS6X, HexS6Y, HexS7X, HexS7Y, HexS8X, HexS8Y, HexS9X, HexS9Y, HexS10X, HexS10Y, HexS11X, HexS11Y, HexS12X, HexS12Y, HexS13X, HexS13Y, HexS14X, HexS14Y, HexS15X, HexS15Y, HexS16X, HexS16Y, HexS17X, HexS17Y, HexS18X, HexS18Y, HexS19X, HexS19Y, HexS20X, HexS20Y, HexS21X, HexS21Y, HexS22X, HexS22Y, HexS23X, HexS23Y, HexS24X, HexS24Y, HexS25X, HexS25Y, HexS26X, HexS26Y, HexS27X, HexS27Y, HexS28X, HexS28Y, HexS29X, HexS29Y, HexS30X, HexS30Y, HexS32X, HexS32Y, HexS33X, HexS33Y, HexS34X, HexS34Y, HexS35X, HexS35Y, HexS36X, HexS36Y],

        num antenna pols:         2,
        num baselines:            8256,
        baselines:                0 v 0 to 127 v 127
        num auto-correlations:    128,
        num cross-correlations:   8128,

        num visibility pols:      4,
        visibility pols:          XX, XY, YX, YY,

        metafits FREQCENT key:    184.96 MHz,

        metafits filename:        1303162952.metafits,
    )

                MWA version:                Correlator v1 Legacy,

                num timesteps:              240,
                timesteps:                  [unix=1619127734.000, gps=1303162952.000, unix=1619127734.500, gps=1303162952.500, unix=1619127735.000, gps=1303162953.000, unix=1619127735.500, gps=1303162953.500, unix=1619127736.000, gps=1303162954.000, unix=1619127736.500, gps=1303162954.500, unix=1619127737.000, gps=1303162955.000, unix=1619127737.500, gps=1303162955.500, unix=1619127738.000, gps=1303162956.000, unix=1619127738.500, gps=1303162956.500, unix=1619127739.000, gps=1303162957.000, unix=1619127739.500, gps=1303162957.500, unix=1619127740.000, gps=1303162958.000, unix=1619127740.500, gps=1303162958.500, unix=1619127741.000, gps=1303162959.000, unix=1619127741.500, gps=1303162959.500, unix=1619127742.000, gps=1303162960.000, unix=1619127742.500, gps=1303162960.500, unix=1619127743.000, gps=1303162961.000, unix=1619127743.500, gps=1303162961.500, unix=1619127744.000, gps=1303162962.000, unix=1619127744.500, gps=1303162962.500, unix=1619127745.000, gps=1303162963.000, unix=1619127745.500, gps=1303162963.500, unix=1619127746.000, gps=1303162964.000, unix=1619127746.500, gps=1303162964.500, unix=1619127747.000, gps=1303162965.000, unix=1619127747.500, gps=1303162965.500, unix=1619127748.000, gps=1303162966.000, unix=1619127748.500, gps=1303162966.500, unix=1619127749.000, gps=1303162967.000, unix=1619127749.500, gps=1303162967.500, unix=1619127750.000, gps=1303162968.000, unix=1619127750.500, gps=1303162968.500, unix=1619127751.000, gps=1303162969.000, unix=1619127751.500, gps=1303162969.500, unix=1619127752.000, gps=1303162970.000, unix=1619127752.500, gps=1303162970.500, unix=1619127753.000, gps=1303162971.000, unix=1619127753.500, gps=1303162971.500, unix=1619127754.000, gps=1303162972.000, unix=1619127754.500, gps=1303162972.500, unix=1619127755.000, gps=1303162973.000, unix=1619127755.500, gps=1303162973.500, unix=1619127756.000, gps=1303162974.000, unix=1619127756.500, gps=1303162974.500, unix=1619127757.000, gps=1303162975.000, unix=1619127757.500, gps=1303162975.500, unix=1619127758.000, gps=1303162976.000, unix=1619127758.500, gps=1303162976.500, unix=1619127759.000, gps=1303162977.000, unix=1619127759.500, gps=1303162977.500, unix=1619127760.000, gps=1303162978.000, unix=1619127760.500, gps=1303162978.500, unix=1619127761.000, gps=1303162979.000, unix=1619127761.500, gps=1303162979.500, unix=1619127762.000, gps=1303162980.000, unix=1619127762.500, gps=1303162980.500, unix=1619127763.000, gps=1303162981.000, unix=1619127763.500, gps=1303162981.500, unix=1619127764.000, gps=1303162982.000, unix=1619127764.500, gps=1303162982.500, unix=1619127765.000, gps=1303162983.000, unix=1619127765.500, gps=1303162983.500, unix=1619127766.000, gps=1303162984.000, unix=1619127766.500, gps=1303162984.500, unix=1619127767.000, gps=1303162985.000, unix=1619127767.500, gps=1303162985.500, unix=1619127768.000, gps=1303162986.000, unix=1619127768.500, gps=1303162986.500, unix=1619127769.000, gps=1303162987.000, unix=1619127769.500, gps=1303162987.500, unix=1619127770.000, gps=1303162988.000, unix=1619127770.500, gps=1303162988.500, unix=1619127771.000, gps=1303162989.000, unix=1619127771.500, gps=1303162989.500, unix=1619127772.000, gps=1303162990.000, unix=1619127772.500, gps=1303162990.500, unix=1619127773.000, gps=1303162991.000, unix=1619127773.500, gps=1303162991.500, unix=1619127774.000, gps=1303162992.000, unix=1619127774.500, gps=1303162992.500, unix=1619127775.000, gps=1303162993.000, unix=1619127775.500, gps=1303162993.500, unix=1619127776.000, gps=1303162994.000, unix=1619127776.500, gps=1303162994.500, unix=1619127777.000, gps=1303162995.000, unix=1619127777.500, gps=1303162995.500, unix=1619127778.000, gps=1303162996.000, unix=1619127778.500, gps=1303162996.500, unix=1619127779.000, gps=1303162997.000, unix=1619127779.500, gps=1303162997.500, unix=1619127780.000, gps=1303162998.000, unix=1619127780.500, gps=1303162998.500, unix=1619127781.000, gps=1303162999.000, unix=1619127781.500, gps=1303162999.500, unix=1619127782.000, gps=1303163000.000, unix=1619127782.500, gps=1303163000.500, unix=1619127783.000, gps=1303163001.000, unix=1619127783.500, gps=1303163001.500, unix=1619127784.000, gps=1303163002.000, unix=1619127784.500, gps=1303163002.500, unix=1619127785.000, gps=1303163003.000, unix=1619127785.500, gps=1303163003.500, unix=1619127786.000, gps=1303163004.000, unix=1619127786.500, gps=1303163004.500, unix=1619127787.000, gps=1303163005.000, unix=1619127787.500, gps=1303163005.500, unix=1619127788.000, gps=1303163006.000, unix=1619127788.500, gps=1303163006.500, unix=1619127789.000, gps=1303163007.000, unix=1619127789.500, gps=1303163007.500, unix=1619127790.000, gps=1303163008.000, unix=1619127790.500, gps=1303163008.500, unix=1619127791.000, gps=1303163009.000, unix=1619127791.500, gps=1303163009.500, unix=1619127792.000, gps=1303163010.000, unix=1619127792.500, gps=1303163010.500, unix=1619127793.000, gps=1303163011.000, unix=1619127793.500, gps=1303163011.500, unix=1619127794.000, gps=1303163012.000, unix=1619127794.500, gps=1303163012.500, unix=1619127795.000, gps=1303163013.000, unix=1619127795.500, gps=1303163013.500, unix=1619127796.000, gps=1303163014.000, unix=1619127796.500, gps=1303163014.500, unix=1619127797.000, gps=1303163015.000, unix=1619127797.500, gps=1303163015.500, unix=1619127798.000, gps=1303163016.000, unix=1619127798.500, gps=1303163016.500, unix=1619127799.000, gps=1303163017.000, unix=1619127799.500, gps=1303163017.500, unix=1619127800.000, gps=1303163018.000, unix=1619127800.500, gps=1303163018.500, unix=1619127801.000, gps=1303163019.000, unix=1619127801.500, gps=1303163019.500, unix=1619127802.000, gps=1303163020.000, unix=1619127802.500, gps=1303163020.500, unix=1619127803.000, gps=1303163021.000, unix=1619127803.500, gps=1303163021.500, unix=1619127804.000, gps=1303163022.000, unix=1619127804.500, gps=1303163022.500, unix=1619127805.000, gps=1303163023.000, unix=1619127805.500, gps=1303163023.500, unix=1619127806.000, gps=1303163024.000, unix=1619127806.500, gps=1303163024.500, unix=1619127807.000, gps=1303163025.000, unix=1619127807.500, gps=1303163025.500, unix=1619127808.000, gps=1303163026.000, unix=1619127808.500, gps=1303163026.500, unix=1619127809.000, gps=1303163027.000, unix=1619127809.500, gps=1303163027.500, unix=1619127810.000, gps=1303163028.000, unix=1619127810.500, gps=1303163028.500, unix=1619127811.000, gps=1303163029.000, unix=1619127811.500, gps=1303163029.500, unix=1619127812.000, gps=1303163030.000, unix=1619127812.500, gps=1303163030.500, unix=1619127813.000, gps=1303163031.000, unix=1619127813.500, gps=1303163031.500, unix=1619127814.000, gps=1303163032.000, unix=1619127814.500, gps=1303163032.500, unix=1619127815.000, gps=1303163033.000, unix=1619127815.500, gps=1303163033.500, unix=1619127816.000, gps=1303163034.000, unix=1619127816.500, gps=1303163034.500, unix=1619127817.000, gps=1303163035.000, unix=1619127817.500, gps=1303163035.500, unix=1619127818.000, gps=1303163036.000, unix=1619127818.500, gps=1303163036.500, unix=1619127819.000, gps=1303163037.000, unix=1619127819.500, gps=1303163037.500, unix=1619127820.000, gps=1303163038.000, unix=1619127820.500, gps=1303163038.500, unix=1619127821.000, gps=1303163039.000, unix=1619127821.500, gps=1303163039.500, unix=1619127822.000, gps=1303163040.000, unix=1619127822.500, gps=1303163040.500, unix=1619127823.000, gps=1303163041.000, unix=1619127823.500, gps=1303163041.500, unix=1619127824.000, gps=1303163042.000, unix=1619127824.500, gps=1303163042.500, unix=1619127825.000, gps=1303163043.000, unix=1619127825.500, gps=1303163043.500, unix=1619127826.000, gps=1303163044.000, unix=1619127826.500, gps=1303163044.500, unix=1619127827.000, gps=1303163045.000, unix=1619127827.500, gps=1303163045.500, unix=1619127828.000, gps=1303163046.000, unix=1619127828.500, gps=1303163046.500, unix=1619127829.000, gps=1303163047.000, unix=1619127829.500, gps=1303163047.500, unix=1619127830.000, gps=1303163048.000, unix=1619127830.500, gps=1303163048.500, unix=1619127831.000, gps=1303163049.000, unix=1619127831.500, gps=1303163049.500, unix=1619127832.000, gps=1303163050.000, unix=1619127832.500, gps=1303163050.500, unix=1619127833.000, gps=1303163051.000, unix=1619127833.500, gps=1303163051.500, unix=1619127834.000, gps=1303163052.000, unix=1619127834.500, gps=1303163052.500, unix=1619127835.000, gps=1303163053.000, unix=1619127835.500, gps=1303163053.500, unix=1619127836.000, gps=1303163054.000, unix=1619127836.500, gps=1303163054.500, unix=1619127837.000, gps=1303163055.000, unix=1619127837.500, gps=1303163055.500, unix=1619127838.000, gps=1303163056.000, unix=1619127838.500, gps=1303163056.500, unix=1619127839.000, gps=1303163057.000, unix=1619127839.500, gps=1303163057.500, unix=1619127840.000, gps=1303163058.000, unix=1619127840.500, gps=1303163058.500, unix=1619127841.000, gps=1303163059.000, unix=1619127841.500, gps=1303163059.500, unix=1619127842.000, gps=1303163060.000, unix=1619127842.500, gps=1303163060.500, unix=1619127843.000, gps=1303163061.000, unix=1619127843.500, gps=1303163061.500, unix=1619127844.000, gps=1303163062.000, unix=1619127844.500, gps=1303163062.500, unix=1619127845.000, gps=1303163063.000, unix=1619127845.500, gps=1303163063.500, unix=1619127846.000, gps=1303163064.000, unix=1619127846.500, gps=1303163064.500, unix=1619127847.000, gps=1303163065.000, unix=1619127847.500, gps=1303163065.500, unix=1619127848.000, gps=1303163066.000, unix=1619127848.500, gps=1303163066.500, unix=1619127849.000, gps=1303163067.000, unix=1619127849.500, gps=1303163067.500, unix=1619127850.000, gps=1303163068.000, unix=1619127850.500, gps=1303163068.500, unix=1619127851.000, gps=1303163069.000, unix=1619127851.500, gps=1303163069.500, unix=1619127852.000, gps=1303163070.000, unix=1619127852.500, gps=1303163070.500, unix=1619127853.000, gps=1303163071.000, unix=1619127853.500, gps=1303163071.500],
                num coarse channels,        24,
                coarse channels:            [gpu=24 corr=23 rec=133 @ 170.240 MHz, gpu=23 corr=22 rec=134 @ 171.520 MHz, gpu=22 corr=21 rec=135 @ 172.800 MHz, gpu=21 corr=20 rec=136 @ 174.080 MHz, gpu=20 corr=19 rec=137 @ 175.360 MHz, gpu=19 corr=18 rec=138 @ 176.640 MHz, gpu=18 corr=17 rec=139 @ 177.920 MHz, gpu=17 corr=16 rec=140 @ 179.200 MHz, gpu=16 corr=15 rec=141 @ 180.480 MHz, gpu=15 corr=14 rec=142 @ 181.760 MHz, gpu=14 corr=13 rec=143 @ 183.040 MHz, gpu=13 corr=12 rec=144 @ 184.320 MHz, gpu=12 corr=11 rec=145 @ 185.600 MHz, gpu=11 corr=10 rec=146 @ 186.880 MHz, gpu=10 corr=9 rec=147 @ 188.160 MHz, gpu=9 corr=8 rec=148 @ 189.440 MHz, gpu=8 corr=7 rec=149 @ 190.720 MHz, gpu=7 corr=6 rec=150 @ 192.000 MHz, gpu=6 corr=5 rec=151 @ 193.280 MHz, gpu=5 corr=4 rec=152 @ 194.560 MHz, gpu=4 corr=3 rec=153 @ 195.840 MHz, gpu=3 corr=2 rec=154 @ 197.120 MHz, gpu=2 corr=1 rec=155 @ 198.400 MHz, gpu=1 corr=0 rec=156 @ 199.680 MHz],

                provided timesteps indices:   230: [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233],
                provided coarse chan indices: 24: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23],
                Common timestep indices:    228: [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231],
                Common coarse chan indices: 24: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23],
                Common UNIX start time:     1619127736,
                Common UNIX end time:       1619127850,
                Common GPS start time:      1303162954,
                Common GPS end time:        1303163068,
                Common duration:            114 s,
                Common bandwidth:           30.72 MHz,

                Common/Good timestep indices:    228: [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231],
                Common/Good coarse chan indices: 24: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23],
                Common/Good UNIX start time:     1619127736,
                Common/Good UNIX end time:       1619127850,
                Common/Good GPS start time:      1303162954,
                Common/Good GPS end time:        1303163068,
                Common/Good duration:            114 s,
                Common/Good bandwidth:           30.72 MHz,

                gpubox HDU size:            32.25 MiB,
                Memory usage per scan:      1548 MiB,

                gpubox batches:             [
        batch_number=0 gpubox_files=[filename=1303162952_20210422214216_gpubox01_00.fits channelidentifier=1, filename=1303162952_20210422214216_gpubox02_00.fits channelidentifier=2, filename=1303162952_20210422214216_gpubox03_00.fits channelidentifier=3, filename=1303162952_20210422214216_gpubox04_00.fits channelidentifier=4, filename=1303162952_20210422214216_gpubox05_00.fits channelidentifier=5, filename=1303162952_20210422214216_gpubox06_00.fits channelidentifier=6, filename=1303162952_20210422214216_gpubox07_00.fits channelidentifier=7, filename=1303162952_20210422214216_gpubox08_00.fits channelidentifier=8, filename=1303162952_20210422214216_gpubox09_00.fits channelidentifier=9, filename=1303162952_20210422214216_gpubox10_00.fits channelidentifier=10, filename=1303162952_20210422214216_gpubox11_00.fits channelidentifier=11, filename=1303162952_20210422214216_gpubox12_00.fits channelidentifier=12, filename=1303162952_20210422214216_gpubox13_00.fits channelidentifier=13, filename=1303162952_20210422214216_gpubox14_00.fits channelidentifier=14, filename=1303162952_20210422214216_gpubox15_00.fits channelidentifier=15, filename=1303162952_20210422214216_gpubox16_00.fits channelidentifier=16, filename=1303162952_20210422214216_gpubox17_00.fits channelidentifier=17, filename=1303162952_20210422214216_gpubox18_00.fits channelidentifier=18, filename=1303162952_20210422214216_gpubox19_00.fits channelidentifier=19, filename=1303162952_20210422214216_gpubox20_00.fits channelidentifier=20, filename=1303162952_20210422214216_gpubox21_00.fits channelidentifier=21, filename=1303162952_20210422214216_gpubox22_00.fits channelidentifier=22, filename=1303162952_20210422214216_gpubox23_00.fits channelidentifier=23, filename=1303162952_20210422214216_gpubox24_00.fits channelidentifier=24],
        batch_number=1 gpubox_files=[filename=1303162952_20210422214316_gpubox01_01.fits channelidentifier=1, filename=1303162952_20210422214316_gpubox02_01.fits channelidentifier=2, filename=1303162952_20210422214316_gpubox03_01.fits channelidentifier=3, filename=1303162952_20210422214316_gpubox04_01.fits channelidentifier=4, filename=1303162952_20210422214316_gpubox05_01.fits channelidentifier=5, filename=1303162952_20210422214316_gpubox06_01.fits channelidentifier=6, filename=1303162952_20210422214316_gpubox07_01.fits channelidentifier=7, filename=1303162952_20210422214316_gpubox08_01.fits channelidentifier=8, filename=1303162952_20210422214316_gpubox09_01.fits channelidentifier=9, filename=1303162952_20210422214316_gpubox10_01.fits channelidentifier=10, filename=1303162952_20210422214316_gpubox11_01.fits channelidentifier=11, filename=1303162952_20210422214316_gpubox12_01.fits channelidentifier=12, filename=1303162952_20210422214316_gpubox13_01.fits channelidentifier=13, filename=1303162952_20210422214316_gpubox14_01.fits channelidentifier=14, filename=1303162952_20210422214316_gpubox15_01.fits channelidentifier=15, filename=1303162952_20210422214316_gpubox16_01.fits channelidentifier=16, filename=1303162952_20210422214316_gpubox17_01.fits channelidentifier=17, filename=1303162952_20210422214316_gpubox18_01.fits channelidentifier=18, filename=1303162952_20210422214316_gpubox19_01.fits channelidentifier=19, filename=1303162952_20210422214316_gpubox20_01.fits channelidentifier=20, filename=1303162952_20210422214316_gpubox21_01.fits channelidentifier=21, filename=1303162952_20210422214316_gpubox22_01.fits channelidentifier=22, filename=1303162952_20210422214316_gpubox23_01.fits channelidentifier=23, filename=1303162952_20210422214316_gpubox24_01.fits channelidentifier=24],
    ],
            )

[2021-09-24T09:04:50Z TRACE birli] img_coarse_chan_range: 0..24
[2021-09-24T09:04:50Z TRACE birli] img_timestep_range: 4..234
[2021-09-24T09:04:50Z TRACE birli] antenna_flags: [34, 47, 56, 63]
[2021-09-24T09:05:00Z TRACE birli] start context_to_jones_array
[2021-09-24T09:07:02Z WARN  birli] could not read hdu ts=232, cc=5 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 5 }
[2021-09-24T09:07:02Z WARN  birli] could not read hdu ts=233, cc=5 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 5 }
[2021-09-24T09:07:03Z WARN  birli] could not read hdu ts=232, cc=4 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 4 }
[2021-09-24T09:07:03Z WARN  birli] could not read hdu ts=233, cc=4 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 4 }
[2021-09-24T09:07:03Z WARN  birli] could not read hdu ts=232, cc=17 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 17 }
[2021-09-24T09:07:03Z WARN  birli] could not read hdu ts=233, cc=17 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 17 }
[2021-09-24T09:07:04Z WARN  birli] could not read hdu ts=232, cc=12 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 12 }
[2021-09-24T09:07:04Z WARN  birli] could not read hdu ts=233, cc=12 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 12 }
[2021-09-24T09:07:05Z WARN  birli] could not read hdu ts=232, cc=8 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 8 }
[2021-09-24T09:07:05Z WARN  birli] could not read hdu ts=233, cc=8 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 8 }
[2021-09-24T09:07:05Z WARN  birli] could not read hdu ts=232, cc=21 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 21 }
[2021-09-24T09:07:05Z WARN  birli] could not read hdu ts=233, cc=21 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 21 }
[2021-09-24T09:07:05Z WARN  birli] could not read hdu ts=232, cc=11 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 11 }
[2021-09-24T09:07:05Z WARN  birli] could not read hdu ts=233, cc=11 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 11 }
[2021-09-24T09:07:05Z WARN  birli] could not read hdu ts=232, cc=1 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 1 }
[2021-09-24T09:07:05Z WARN  birli] could not read hdu ts=233, cc=1 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 1 }
[2021-09-24T09:07:05Z WARN  birli] could not read hdu ts=232, cc=13 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 13 }
[2021-09-24T09:07:05Z WARN  birli] could not read hdu ts=233, cc=13 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 13 }
[2021-09-24T09:07:06Z WARN  birli] could not read hdu ts=232, cc=22 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 22 }
[2021-09-24T09:07:06Z WARN  birli] could not read hdu ts=233, cc=22 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 22 }
[2021-09-24T09:07:09Z WARN  birli] could not read hdu ts=232, cc=16 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 16 }
[2021-09-24T09:07:09Z WARN  birli] could not read hdu ts=233, cc=16 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 16 }
[2021-09-24T09:07:09Z WARN  birli] could not read hdu ts=232, cc=14 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 14 }
[2021-09-24T09:07:09Z WARN  birli] could not read hdu ts=233, cc=14 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 14 }
[2021-09-24T09:07:10Z WARN  birli] could not read hdu ts=232, cc=3 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 3 }
[2021-09-24T09:07:10Z WARN  birli] could not read hdu ts=233, cc=3 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 3 }
[2021-09-24T09:07:10Z WARN  birli] could not read hdu ts=232, cc=0 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 0 }
[2021-09-24T09:07:10Z WARN  birli] could not read hdu ts=233, cc=0 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 0 }
[2021-09-24T09:07:11Z WARN  birli] could not read hdu ts=232, cc=20 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 20 }
[2021-09-24T09:07:11Z WARN  birli] could not read hdu ts=233, cc=20 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 20 }
[2021-09-24T09:07:12Z WARN  birli] could not read hdu ts=232, cc=15 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 15 }
[2021-09-24T09:07:12Z WARN  birli] could not read hdu ts=233, cc=15 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 15 }
[2021-09-24T09:07:13Z WARN  birli] could not read hdu ts=232, cc=18 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 18 }
[2021-09-24T09:07:13Z WARN  birli] could not read hdu ts=233, cc=18 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 18 }
[2021-09-24T09:07:14Z WARN  birli] could not read hdu ts=232, cc=19 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 19 }
[2021-09-24T09:07:14Z WARN  birli] could not read hdu ts=233, cc=19 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 19 }
[2021-09-24T09:07:14Z WARN  birli] could not read hdu ts=232, cc=10 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 10 }
[2021-09-24T09:07:14Z WARN  birli] could not read hdu ts=233, cc=10 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 10 }
[2021-09-24T09:07:15Z WARN  birli] could not read hdu ts=232, cc=23 NoDataForTimeStepCoarseChannel { timestep_index: 232, coarse_chan_index: 23 }
[2021-09-24T09:07:15Z WARN  birli] could not read hdu ts=233, cc=23 NoDataForTimeStepCoarseChannel { timestep_index: 233, coarse_chan_index: 23 }
[2021-09-24T09:07:15Z TRACE birli] end context_to_jones_array
[2021-09-24T09:07:15Z DEBUG birli] Applying cable delays. applied: false, desired: true
[2021-09-24T09:07:15Z TRACE birli::corrections] start correct_cable_lengths
[2021-09-24T09:07:33Z TRACE birli::corrections] end correct_cable_lengths
[2021-09-24T09:07:33Z DEBUG birli] flagging with strategy /usr/share/aoflagger/strategies/mwa-default.lua
[2021-09-24T09:07:33Z TRACE birli::flags] start flag_jones_array
[2021-09-24T09:14:05Z TRACE birli::flags] end flag_jones_array
[2021-09-24T09:14:05Z DEBUG birli] Applying geometric delays. applied: No, desired: true
[2021-09-24T09:14:05Z TRACE birli::corrections] start correct_geometry
[2021-09-24T09:14:23Z TRACE birli::corrections] end correct_geometry
[2021-09-24T09:14:23Z TRACE birli::io] start write_uvfits
[2021-09-24T09:30:13Z WARN  birli::io::uvfits] we are using MWA lat / lng / height from mwalib for array position
[2021-09-24T09:54:25Z TRACE birli::io] end write_uvfits
[2021-09-24T09:54:26Z INFO  birli] end main
```
