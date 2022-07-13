tb = locals()['tb']
dump_path = tb.name() + ".csv"
num_freqs = tb.getcoldesc('DATA')['shape'].tolist()[1]
row_limit = locals().get('row_limit', 250)

with open(dump_path, 'w') as dump:
    header = ["time", "ant1", "ant2", "u", "v", "w", "type", "pol"]

    dump.write(",".join(header + list(map(str, range(num_freqs)))) + "\n")
    for row_idx in range(0, min(tb.nrows(), row_limit)):
        time = tb.getcell("TIME", row_idx)
        ant1 = tb.getcell("ANTENNA1", row_idx)
        ant2 = tb.getcell("ANTENNA2", row_idx)
        if ant1 == ant2:
            continue
        uvw = tb.getcell("UVW", row_idx).tolist()
        row = [time, ant1, ant2, *uvw]

        vis = tb.getcell("DATA", row_idx).tolist()
        flag = tb.getcell("FLAG", row_idx).tolist()
        weight = tb.getcell("WEIGHT_SPECTRUM", row_idx).tolist()
        for pol_idx, pol in enumerate(["xx", "xy", "yx", "yy"]):
            dump.write(",".join(map(str, row[:] + ["vis", pol] + vis[pol_idx])) + "\n")
            dump.write(",".join(map(str, row[:] + ["weight", pol] + weight[pol_idx])) + "\n")
            dump.write(",".join(map(str, row[:] + ["flag", pol] + flag[pol_idx])) + "\n")
