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

| bench                                              | goal  | a762a27        | ea16381       | dbdc0ee       | 09abc69       | 3a5671c         | ee893da        |
| -------------------------------------------------- | ----- | -------------- | ------------- | ------------- | ------------- | --------------- | -------------- |
| context_to_baseline_imgsets - ord_half_1196175296  |       | 1873.5 ± 8.2 s | 729.49 ± 6.0s | 458.65 ± 6.8s | 172.69 ± 4.8s | 98.537 ± 2.7s   |                |
| context_to_baseline_imgsets - mwax_half_1247842824 |       |                | 372.64 ± 1.4s | 226.30 ± 5.7s | 87.029 ± 2.1s | 51.940 ± 1.9s   |                |
| context_to_jones_array - ord_half_1196175296       | ~300s |                |               |               |               | 64.717 ± 2.6s   | 70.697 ± 3.7s  |
| context_to_jones_array - mwax_half_1247842824      |       |                |               |               |               | 35.321 ± 0.61s  | 37.996 ± 0.36s |
| correct_cable_lengths - ord_half_1196175296        |       |                |               |               |               | 10.543 ± 0.14s  | 18.495 ± 0.36s |
| correct_cable_lengths - mwax_half_1247842824       |       |                |               |               |               | 6.3826 ± 0.025s | 9.2928 ± 0.14s |
| correct_geometry - ord_half_1196175296             |       |                |               |               |               | 650.03 ± 7.4s   | 9.5513 ± 0.31s |
| correct_geometry - mwax_half_1247842824            |       |                |               |               |               | 252.79 ± 5.1s   | 6.5979 ± 1.09s |
| uvfits_output - 1254670392_avg                     | ~2s   |                |               |               |               | 7.7792 ± 0.059s | 4.8460 ± 0.04s |

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
