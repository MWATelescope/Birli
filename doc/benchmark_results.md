# Benchmark Results

## Cotter

Cotter, when run on the same machine with the same data using `scripts/time_cotter.sh` multiple times produces the following results.

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

<!-- tip: use https://tabletomarkdown.com/convert-website-table-to-markdown/ -->

## [commit a762a27 (v0.1.1)](https://github.com/MWATelescope/Birli/actions/runs/812007480)

total duration:  5h 47m 35s

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

total duration:  3h 31m 2s

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

total duration:  2h 8m 33s

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

![dbdc0ee context_to_baseline_imgsets - ord_half_1196175296 pdf](img/dbdc0ee%20context_to_baseline_imgsets%20-%20ord_half_1196175296%20pdf.svg)

performance changes:

- parallelize writing to img_bufs
- parallelize over polarizations
- use Vec instead of BTreeMap where possible.
- use chunks_exact instead of chunks

total duration:  52m 52s

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
