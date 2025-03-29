# command line interface ------------------------------------------------------#
library(docopt)
library(tidyverse)
'Merge assembly info .tsv and genomad plasmid summary .tsv into a bakta replicon .tsv.

Usage:
  bakta-replicons.R [options] <asm.tsv> <plasmid.tsv> <out.tsv>

Replicon table format: https://github.com/oschwengers/bakta?tab=readme-ov-file#input-and-output

' -> doc

opt <- docopt(doc)
# devel
# opt <- docopt(doc, "analysis/restart/vup37/vup37.tsv analysis/genomad/vup37/vup37_summary/vup37_plasmid_summary.tsv foo.tsv")
#opt <- docopt(doc, "~/Code/projects/detectEVE/issues/17/anopheles-lequime-2017_rvdb-notax/results/APHL01-retro.bed foo.tsv bar.tsv")
#opt <- map_at(opt, ~str_detect(., "score"), as.numeric)

library(tidyverse)
a0 <- read_tsv(opt[["asm.tsv"]])
p0 <- read_tsv(opt[["plasmid.tsv"]])
a1 <- a0 |> 
	transmute(
		seq_id, length,
		new_id="",
		topology=if_else(circular == "Y", "circular", ""),
		type = case_when(
			seq_id %in% p0$seq_name ~ "plasmid",
			length > 1e6 & topology == "circular" ~ "chromosome",
			.default = ""),
		name="") |> 
	select(-length)

print(opt[["out.tsv"]])
write_tsv(a1, opt[["out.tsv"]])