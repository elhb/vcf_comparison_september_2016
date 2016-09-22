#!/usr/bin/env python

#
# set input info
#
DEBUG = True
DATA_PATH ='/proj/b2013064/private/ACTIVE/erik/JG_PL.vcf_parsing/data/'
RESULTS_PATH ='/proj/b2013064/private/ACTIVE/erik/JG_PL.vcf_parsing/results/'
VCF_FILE_NAMES = [
    DATA_PATH+'P4107_1003.clean.dedup.recal.bam.raw.annotated.vcf.gz',
    DATA_PATH+'P4728_1004.clean.dedup.recal.bam.raw.annotated.vcf.gz',
    DATA_PATH+'P4728_1005.clean.dedup.recal.bam.raw.annotated.vcf.gz',
    DATA_PATH+'P4728_1006.clean.dedup.recal.bam.raw.annotated.vcf.gz'
]
samples = ['P4107_1003','P4728_1004','P4728_1005','P4728_1006']

#
# import
#
import sys
import time
start_time = time.time()

#
# Connect to database
#
import sqlite3
database = sqlite3.connect(RESULTS_PATH+'variations_and_samples_database.sqlite3.db')

sys.stderr.write( '# INFO: loading variations from sample tables...\n')
variation_ids_by_sample = {}
for sample_name in samples:
    sys.stderr.write( '# INFO: '+sample_name+'...\n')
    variation_ids_by_sample[sample_name] = { tmp[0]:True for tmp in database.cursor().execute('SELECT inc_int_id FROM '+sample_name+'') }

sys.stderr.write( '# INFO: sorting variations within sample...\n')
samples_by_variation_ids = {}
for sample_name in samples:
    sys.stderr.write( '# INFO: '+sample_name+'...\n')
    for variation_id in variation_ids_by_sample[sample_name].keys():
        try: samples_by_variation_ids[variation_id].append( sample_name )
        except KeyError: samples_by_variation_ids[variation_id] = [sample_name]

overlap_counter = {}
for variation_id, sample_overlap in samples_by_variation_ids.iteritems():
    try:
        overlap_counter[','.join(sample_overlap)] += 1
    except KeyError:
        overlap_counter[','.join(sample_overlap)] = 1

for overlap, count in overlap_counter.iteritems():print count,'\t',overlap

database.close()

sys.stderr.write( 'DONE, finished in '+ str(time.time()-start_time)+' seconds\n' )