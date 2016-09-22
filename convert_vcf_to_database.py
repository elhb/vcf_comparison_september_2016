#!/usr/bin/env python

#
# set input info
#
DEBUG = False
DATA_PATH ='/proj/b2013064/private/ACTIVE/erik/JG_PL.vcf_parsing/data/'
RESULTS_PATH ='/proj/b2013064/private/ACTIVE/erik/JG_PL.vcf_parsing/results/'
VCF_FILE_NAMES = [
    DATA_PATH+'P4107_1003.clean.dedup.recal.bam.raw.annotated.vcf.gz',
    DATA_PATH+'P4728_1004.clean.dedup.recal.bam.raw.annotated.vcf.gz',
    DATA_PATH+'P4728_1005.clean.dedup.recal.bam.raw.annotated.vcf.gz',
    DATA_PATH+'P4728_1006.clean.dedup.recal.bam.raw.annotated.vcf.gz'
]

#
# import
#
import gzip
import re
import sys
import time
start_time = time.time()

#
# set initial values
#
variant_positions = {}
samples = []

#
# go through the vcf files
#
for vcf_file_name in VCF_FILE_NAMES:
    vcf_file_start_time = time.time()
    
    # get the sample name and add to list
    sample_name = re.search("P\d{4}\_\d{4}",vcf_file_name).group()
    samples.append(sample_name)
    
    # open file
    vcf_file = gzip.open(vcf_file_name)
    tmp_line_counter = 0
    
    #
    # save each vcf line to variations dictionary
    #
    for line in vcf_file:
        
        if line[0] == '#': continue # skip header
        tmp_line_counter += 1 # increment line counter
        line = line.rstrip().split('\t') # split according to tsv format
        sub_line = [line[2],line[3],line[4],line[5],line[8],line[9]] # trash some info
        
        chromosome_name = line[0]
        position = int(line[1])
        
        # save to in mem dictionary
        try: variant_positions[ chromosome_name ][ position ][sample_name] = sub_line
        except KeyError as err:
            if str(err)[1:-1] == chromosome_name:
                variant_positions[ chromosome_name ] = { position:{sample_name:sub_line} }
            elif str(err) == str(position):
                variant_positions[ chromosome_name ][ position ] = {sample_name:sub_line}
            else:
                print err,chromosome_name,position
                raise err
#
# too keep track of list indexes
#                 0      1        2       3                       4       5 SUBLINE
# 0       1       2      3        4       5       6       7       8       9 LINE
# CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  P4107_1003
        
        # for debugging
        if DEBUG and tmp_line_counter >= 1e5: break
    sys.stderr.write('#INFO: vcf file'+vcf_file_name+' read to mem in '+str(time.time()-vcf_file_start_time)+' seconds.\n')
sys.stderr.write('#INFO: '+str(sum([len(positions) for chromosome_name, positions in variant_positions.iteritems()]))+' variants read in total.\n')    
#if DEBUG:
#    for chromosome_name, positions in variant_positions.iteritems():
#        for position, samples_info in positions.iteritems():
#            print chromosome_name, position, {sample_name:sub_line[2] for sample_name,sub_line in samples_info.iteritems()}

#
# parse through variations in sync
#
variations_database_values = []
per_sample_database_values = {sample_name:[] for sample_name in samples}
id_counter = xrange(int(3e9)).__iter__() # define a counter

for chromosome_name, positions in variant_positions.iteritems(): # go thorug all positions on al chromomes
    for position, samples_info in positions.iteritems():          # each uniq pos is a variation
        
        inc_int_id = id_counter.next()

        variation_id = {}
        reference_allele = {}
        alternative_alleles = {}
    
        for sample_name, sub_line in samples_info.iteritems():
            
            try:             variation_id[ sub_line[0] ].append(sample_name)
            except KeyError: variation_id[ sub_line[0] ] = [sample_name]
            
            try:             reference_allele[ sub_line[1] ].append(sample_name)
            except KeyError: reference_allele[ sub_line[1] ] = [sample_name]
             
            try:             alternative_alleles[ sub_line[2] ].append(sample_name)
            except KeyError: alternative_alleles[ sub_line[2] ] = [sample_name]

        # MIGHT NEED FIXING DUE TO MORE THAN ONE REF AND ALT ALLELES: look at this some time later...
        #if len(variation_id) != 1: print 'problem 1',variation_id,reference_allele,alternative_alleles;sys.exit()
        #if len(reference_allele) != 1:  print 'problem 2',variation_id,reference_allele,alternative_alleles;sys.exit()
        #if len(alternative_alleles.keys()) != 1:  print 'problem 3',variation_id,reference_allele,alternative_alleles;sys.exit()

        for sample_name, sub_line in samples_info.iteritems():
            format_string = sub_line[4]
            info_string   = sub_line[5]
            
            info_dict = {key:value for key,value in zip(format_string.split(':'),info_string.split(':'))}
            
            QUAL = sub_line[3]

            GT = info_dict['GT'] # MIGHT NEED FIXING DUE TO MORE REF AND ALT ALLELES
            AD = info_dict['AD'] # MIGHT NEED FIXING DUE TO MORE REF AND ALT ALLELES
            try: DP = info_dict['DP']
            except KeyError as err:
                assert str(err) == "'DP'", str(err)+'!='+"'DP'"
                DP = sum([int(i) for i in AD.split(',')])
                sys.stderr.write('#WARNING: no DP value found for variant at chromosome='+chromosome_name+' position='+str(position)+' in sample='+sample_name+', using sum of AD instead ('+str(DP)+'='+str(AD)+')\n')
            GQ = info_dict['GQ']
            PL = info_dict['PL'] # MIGHT NEED FIXING DUE TO MORE REF AND ALT ALLELES
            
            per_sample_database_values[sample_name].append( (inc_int_id, QUAL, GT, AD, DP, GQ, PL) )
        
        #print inc_int_id,chromosome_name,position,variation_id,reference_allele,alternative_alleles
        variations_database_values.append( (inc_int_id,chromosome_name,position,str(variation_id),str(reference_allele),str(alternative_alleles)) )
        
sys.stderr.write( '#INFO: '+str(inc_int_id)+ ' variations read into mem and parsed, adding to database...\n')

#
# Create database
#
import sqlite3
database = sqlite3.connect(RESULTS_PATH+'variations_and_samples_database.sqlite3.db')

#
# Create table to hold the variations
#
sys.stderr.write('#INFO: creating variations table in database.\n')
create_variation_table_string = """CREATE TABLE variations (
inc_int_id int,
chromosome_name varchar(10),
position int(9),
variation_id varchar(255),
reference_allele varchar(255),
alternative_alleles TEXT,
PRIMARY KEY (chromosome_name, position))"""
try: database.cursor().execute(create_variation_table_string)
except sqlite3.OperationalError as err:
    if str(err) == 'table variations already exists': pass
    else:
        raise err

#
# Create tables to hold per sample info
#
for sample_name in samples:
    sys.stderr.write('#INFO: creating sample specific table in database for sample '+sample_name+'.\n')
    create_sample_table_string = """CREATE TABLE """+sample_name+""" (
    inc_int_id int,
    QUAL FLOAT,
    GT int(9),
    AD varchar(255),
    DP int,
    GQ int(3),
    PL varchar(255),
    PRIMARY KEY (inc_int_id))"""
    try: database.cursor().execute(create_sample_table_string)
    except sqlite3.OperationalError as err:
        if str(err) == 'table '+sample_name+' already exists': pass
        else:
            raise err

#
# fill the tables
#
sys.stderr.write('#INFO: filling variations table with data.\n')
fill_variations_table_string = 'INSERT INTO variations VALUES (?, ?, ?, ?, ?, ?)'
database.cursor().executemany(fill_variations_table_string, variations_database_values)
for sample_name in samples:
    sys.stderr.write('#INFO: filling sample specific table in database for sample '+sample_name+'.\n')
    fill_samples_table_string = 'INSERT INTO '+sample_name+' VALUES (?, ?, ?, ?, ?, ?, ?)'
    database.cursor().executemany(fill_samples_table_string, per_sample_database_values[sample_name])

database.commit()
database.close()

sys.stderr.write( '#DONE, finished in '+ str(time.time()-start_time)+' seconds\n' )