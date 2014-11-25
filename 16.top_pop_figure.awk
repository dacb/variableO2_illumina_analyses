// {
	++line;
	for (i = 1; i <= NF; ++i) {
		data[line, i] = $i;
	}
	nf = NF;
}
END {
	for (i = 1; i <= nf; ++i) {
		split(data[1, i], a, "_");
		if (a[2] == "5A" || a[2] == "5B" || a[2] == "5C")
			of = "16.top_pop_figure.15uM.tab";
		else if (a[2] == "15A" || a[2] == "15B" || a[2] == "15C")
			of = "16.top_pop_figure.45uM.tab";
		else if (a[2] == "25A" || a[2] == "25B" || a[2] == "25C")
			of = "16.top_pop_figure.75uM.tab";
		else if (a[2] == "50A" || a[2] == "50B" || a[2] == "50C")
			of = "16.top_pop_figure.150uM.tab";
		else if (a[2] == "75A" || a[2] == "75B" || a[2] == "75C")
			of = "16.top_pop_figure.225uM.tab";
		else {
			files[1] = "16.top_pop_figure.15uM.tab";
			files[2] = "16.top_pop_figure.45uM.tab";
			files[3] = "16.top_pop_figure.75uM.tab";
			files[4] = "16.top_pop_figure.150uM.tab";
			files[5] = "16.top_pop_figure.225uM.tab";
			for (f in files) {
				of = files[f]
				for (j = 1; j <= line; ++j) {
					printf("%s\t", data[j, i]) >> of;
				}
				printf("\n") >> of;
			}
			continue;
		}
		for (j = 1; j <= line; ++j) {
			printf("%s\t", data[j, i]) >> of;
		}
		printf("\n") >> of;
	}
}
