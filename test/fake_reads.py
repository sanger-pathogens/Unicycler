import os
import unicycler.misc
import random


def make_fake_qual_string(length):
    qual_string = ''
    for i in range(length):
        qual_string += chr(random.randint(70, 75))
    return qual_string


def make_fake_reads(seq):
    """
    This function makes super-simple fake reads (no errors and even distribution) for a circular
    genome and saves them to FASTQ files.
    """
    read_length = 100
    insert_length = 300
    looped_seq_forward = seq + seq[0:insert_length]
    looped_seq_reverse = unicycler.misc.reverse_complement(looped_seq_forward)

    out_dir = 'TEST_TEMP_' + str(os.getpid())
    os.makedirs(out_dir)

    reads_1 = os.path.join(out_dir, 'reads_1.fastq')
    reads_2 = os.path.join(out_dir, 'reads_2.fastq')

    read_num = 1
    with open(reads_1, 'wt') as r_1, open(reads_2, 'wt') as r_2:
        for i in range(0, 2):
            for j in range(0, len(looped_seq_forward) - insert_length + 1):
                if i == 0:
                    looped_seq = looped_seq_forward
                else:
                    looped_seq = looped_seq_reverse
                insert = looped_seq[j:j+insert_length]
                read_1 = insert[:read_length]
                read_2 = insert[-read_length:]
                r_1.write('@read_' + str(read_num) + '/1\n')
                r_2.write('@read_' + str(read_num) + '/2\n')
                r_1.write(read_1 + '\n')
                r_2.write(read_2 + '\n')
                r_1.write('+\n')
                r_2.write('+\n')
                r_1.write(make_fake_qual_string(read_length) + '\n')
                r_2.write(make_fake_qual_string(read_length) + '\n')
                read_num += 1
    return out_dir
