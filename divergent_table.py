import re


def cfreq_suffix():
	return "cons_freq"


def gfreq_suffix():
	return "gap_freq"


def get_header(other_groups):
	cf = cfreq_suffix()
	gf = gfreq_suffix()
	h = "\t".join(["ind", "cons", cf, gf])
	for g in other_groups:
		h += "\t" + g + "_" + cf + "\t" + g + "_" + gf
	return h


def get_cfreq_colnames(colnames):
	r = re.compile(".*_cons_freq")
	cols = list(filter(r.match, colnames)) # Read Note below
	return cols



