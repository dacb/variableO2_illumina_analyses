{
	++line;
	if (line == 1) {
		for (i = 1; i <= NF; ++i)
			fields[i] = $i;
	} else {
		for (i = 1; i <= NF; ++i)
			data[line - 1, fields[i]] = $i;
	}
} END {
	for (fi in fields) {
		f = fields[fi];
		select = 0;
		week = 0;
		if (f == "sediment.rep.1")
			select = 1;
		else {
			n = split(f, a, ".");
			if (a[1] == expt) {
				select = 1;
				week = a[2];
			}
		}
		if (select) {
			for (i = 1; i < line; ++i) {
				printf("%s\t%d\t%s\t%s\t%d\t%f\n", expt, week, f, data[i, fields[1]], data[i, f], log10(data[i, f]));
			}
		}
	}
}
function log10(value) {
	return log(value)/log(10)
}
