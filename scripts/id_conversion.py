#!/usr/bin/env python3

import pandas as pd
import sys

df1 = pd.read_table(sys.argv[1], header=0)
df2 = pd.read_table(sys.argv[2], header=0)

df1.merge(df2, on='TE_Family')[['Chrom', 'Start', 'End', 'TEfamily', 'Type', 'TypeA', 'TypeB', 'TypeC', 'PASS_FILTER' ]].to_csv(sys.argv[3], sep="\t", header=True, index=False)
