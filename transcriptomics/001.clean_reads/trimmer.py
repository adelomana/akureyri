###
### usage: time python trimmer.py &> messages.trimmer.txt
###

import os, datetime, sys, re

def printt(label):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(label)))

    return None

def trimmomatic_caller(sample):

    executable='time java -jar {}trimmomatic-0.39.jar PE -threads {} -phred33 '.format(trimmomatic_path,number_threads)
    options=' ILLUMINACLIP:{}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'.format(adapter_file)

    input1 = raw_fastq_dir + sample + '_1.fq.gz'
    input2 = raw_fastq_dir + sample + '_2.fq.gz'

    output_dir = clean_fastq_dir + sample + '/'
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)

    output1 = output_dir + sample + '_R1_clean.fastq.gz'
    output2 = output_dir + sample + '_R2_clean.fastq.gz'

    garbage1 = output_dir + sample + '_R1_garbage.fastq.gz'
    garbage2 = output_dir + sample + '_R2_garbage.fastq.gz'

    input_files = input1 + ' ' + input2
    output_files = output1 + ' ' + garbage1 + ' ' + output2 + ' ' + garbage2

    command = executable + input_files + ' ' + output_files + options

    printt('about to clean {}'.format(sample))
    print('')
    print(command)
    print('')
    os.system(command)
    print('')
    #sys.exit()

    return None

# 0. user defined variablessample
raw_fastq_dir = '/Users/adrian/research/akureyri/data/raw/'
clean_fastq_dir = '/Users/adrian/research/akureyri/data/clean_fastq/'
trimmomatic_path = '/Users/adrian/software/Trimmomatic-0.39/'
adapter_file = trimmomatic_path + 'adapters/TruSeq3-PE-2.fa'
number_threads = 4 # not sure if there is an impact on 4 vs 8. Need to profile.
# @ 8 threads, python trimmer.py &> messages.trimmer.txt  46278.63s user 737.45s system 301% cpu 4:20:09.23 total. it took four hours and a quarter.
# @ 4 threads, python trimmer.py &> messages.trimmer.txt  50630.44s user 797.15s system 302% cpu 4:43:45.32 total. four hours and forty minutes

# 1. recover samples
all_files = os.listdir(raw_fastq_dir)

labels = []
for file in all_files:
    label = '_'.join(file.split('_')[:-1])
    labels.append(label)
unique_labels = list(set(labels))
unique_labels.sort()
print(unique_labels)
print()

# 2. iterate Trimmomatic
for label in unique_labels:
    trimmomatic_caller(label)
