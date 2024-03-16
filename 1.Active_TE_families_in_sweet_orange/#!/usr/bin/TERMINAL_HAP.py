#!/usr/bin/env python3

import pandas as pd
import os
import sys

PRE=sys.argv[1]

RIGHT_ALL=pd.read_table(PRE+'_RIGHT_50bp.align.fasta',names=['ID','R_HAP','R_COUNT']).drop(['R_COUNT'],axis=1)
RIGHT_COUNT = RIGHT_ALL.groupby('R_HAP').count().reset_index().rename(columns={"ID": "R_COUNT"})
RIGHT_COUNT = RIGHT_COUNT.sort_values('R_COUNT',ascending=False)
RIGHT_COUNT = RIGHT_COUNT.reset_index(drop=True)
RIGHT_COUNT['RHAP_ID'] = 'R'+RIGHT_COUNT.index.astype(str)
RIGHT_ALL = RIGHT_ALL.merge(RIGHT_COUNT, on='R_HAP')

LEFT_ALL=pd.read_table(PRE+'_LEFT_50bp.align.fasta',names=['ID','L_HAP','L_COUNT']).drop(['L_COUNT'],axis=1)
LEFT_COUNT = LEFT_ALL.groupby('L_HAP').count().reset_index().rename(columns={"ID": "L_COUNT"})
LEFT_COUNT = LEFT_COUNT.sort_values('L_COUNT',ascending=False)
LEFT_COUNT = LEFT_COUNT.reset_index(drop=True)
LEFT_COUNT['LHAP_ID'] = 'L'+LEFT_COUNT.index.astype(str)
LEFT_ALL = LEFT_ALL.merge(LEFT_COUNT, on='L_HAP')

ALL = LEFT_ALL.merge(RIGHT_ALL,on='ID')
ALL.to_csv(PRE+'_MEMBER_TERMINAL_HAP.tsv', sep="\t", header=True, index=False)
ALL_SUB=ALL.loc[(ALL['L_COUNT']>1) | (ALL['R_COUNT']>1) ]

LEFT=ALL_SUB[['LHAP_ID', 'L_HAP']].drop_duplicates()
LEFT['LHAP_ID']=">"+LEFT['LHAP_ID']
LEFT.to_csv(PRE+'_LEFT_50bp_SUB.fasta', sep="\t", header=False, index=False)

RIGHT=ALL_SUB[['RHAP_ID', 'R_HAP']].drop_duplicates()
RIGHT['RHAP_ID']=">"+RIGHT['RHAP_ID']
RIGHT.to_csv(PRE+'_RIGHT_50bp_SUB.fasta', sep="\t", header=False, index=False)
