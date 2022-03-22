# My charseq pipeline example command sequence:
# date10.02.19
# the plan: trim (fastp) -> check quality (fastqc R1 and R2 seperately) -> PEAR (merge overlapping paired end reads) -> char_bridge tools -> quantify % contaminating free-floating RNA (with Human) -> -> mapp -> ...

## command line execution protocol
# import subprocess
## to simply execute command with arguements
# subprocess.call(["ls", "-l", "/etc/resolv.conf"])
## to retain output
# p = subprocess.Popen(["ls", "-l", "/etc/resolv.conf"], stdout=subprocess.PIPE)
# cmdout, err = p.communicate()
# print "*** Running ls -l command ***\n", cmdout


def seaborn_kde_cumulative(x, outpath):
    # x is list of floats normalized (by DPNII sites and by gene length)
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt

    print('Plotting cumulative distribution')

    fig, axes = plt.subplots(figsize=(8, 8))  # 2, 1, gridspec_kw={'height_ratios': [4, 1]}

    sns.kdeplot(x, cumulative=True, ax=axes)  # vmin=0, vmax=1,
    # axes[1].barh(y=[0,0], width=[0.2,0.4], left=[0,0.5])
    if outpath != '':
        fig.savefig(outpath)
        print('Finished making cumulative distribution: %s' % outpath)
    else:
        print('Finished making cumulative distribution')

    plt.show()
    plt.close()

    return


def retrieve_rna_sequences_from_contacts_inside_features(seqnames, rnafastqfiles, dnabuffer=0, rnabuffer=0):
    import gzip
    import os

    totalrnacount = 0
    byteseqnames = [bytes(k, encoding='utf8') for k in seqnames]
    outrnas, outshortrnas, outiesrnas, outscnrnas = [], [], [], []
    print(rnafastqfiles[:5])
    print(len(rnafastqfiles))
    for filecount, rnafile in enumerate(rnafastqfiles):
        drna = {}  # mainly here to clear memory
        with gzip.open(rnafile, 'r') as FILE:
            for count, line in enumerate(FILE):
                if count % 4 == 0:  # then it is read name line
                    name = line.strip().split()[0][1:]  # remove @ and 1:N:0:CAGATC
                elif count % 4 == 1:  # then it is sequence line
                    drna[name] = line.strip()
                    totalrnacount += 1
        # append sequences to output lists after every file
        for k in byteseqnames:  # unfortunately i have to scuff the solution to save memory
            try:
                drna[k]
            except KeyError:
                pass
            else:
                # the sequence needs to be output
                seqlen = len(drna[k])
                outrnas.append(b'>%s\n%s' % (k, drna[k]))
                if (seqlen <= 30) and (seqlen >= 20):
                    outshortrnas.append(b'>%s\n%s' % (k, drna[k]))
                if (seqlen <= 30) and (seqlen >= 25):
                    outiesrnas.append(b'>%s\n%s' % (k, drna[k]))
                if seqlen == 25:
                    outscnrnas.append(b'>%s\n%s' % (k, drna[k]))
        print(f'Retrieved sequences from {filecount+1} files')

    # print(list(drna.items())[:2])
    print(f'Total Number of RNAs (after pear, before alignment): {totalrnacount}')  # len(list(drna.items()))
    # print(seqnames[:2])
    print(f'Total strength of all contacts in all features: {len(seqnames)}')
    # setnames = set([bytes(i, encoding='utf8') for i in seqnames])
    # setkeys = set(list(drna.keys()))
    # intersect = setnames.intersection(setkeys)
    # un = setnames.union(setkeys)
    # diff = setnames.difference(setkeys) # The difference between two sets results in a third set with the element from the first, that are not present on the second.
    # print(f'intersected: {len(intersect)}, {list(intersect)[:5]}')
    # print(f'union: {len(un)}, {list(un)[:5]}')
    # print(f'difference: {len(diff)}, {list(diff)[:5]}')

    #outrnas = [b'>%s\n%s' % (k, drna[k]) for k in byteseqnames]
    #outshortrnas = [b'>%s\n%s' % (k, drna[k]) for k in byteseqnames if len(drna[k]) <= 30]
    #outiesrnas = [b'>%s\n%s' % (k, drna[k]) for k in byteseqnames if (len(drna[k]) >= 25) and (len(drna[k]) <= 30)]
    #outscnrnas = [b'>%s\n%s' % (k, drna[k]) for k in byteseqnames if len(drna[k]) == 25]

    print(f'Number of unique rna sequences from contacts in features: {len(outrnas)}')
    print(f'Number of unique short rna sequences (<=30bps) from contacts in features: {len(outshortrnas)}')
    print(f'Number of unique short rna sequences (>=25bps, <=30bps) from contacts in features: {len(outiesrnas)}')
    print(f'Number of unique short rna sequences (==25bps) from contacts in features: {len(outscnrnas)}')

    path, file = os.path.split(rnafastqfiles[0])
    outfile = os.path.join(path, file.split('split')[0] + f'RNAs.ContactsNearIES.dnabuffer{dnabuffer}.rnabuffer{rnabuffer}.fa.gz')
    outshortfile = os.path.join(path, file.split('split')[0] + f'RNAs.ContactsNearIES.dnabuffer{dnabuffer}.rnabuffer{rnabuffer}.Short.fa.gz')
    outiesfile = os.path.join(path, file.split('split')[0] + f'RNAs.ContactsNearIES.dnabuffer{dnabuffer}.rnabuffer{rnabuffer}.25-30bp.fa.gz')
    outscnfile = os.path.join(path, file.split('split')[0] + f'RNAs.ContactsNearIES.dnabuffer{dnabuffer}.rnabuffer{rnabuffer}.25bp.fa.gz')
    with gzip.open(outfile, 'w') as OUT:
        OUT.write(b'\n'.join(outrnas))
    with gzip.open(outshortfile, 'w') as OUT:
        OUT.write(b'\n'.join(outshortrnas))
    with gzip.open(outiesfile, 'w') as OUT:
        OUT.write(b'\n'.join(outiesrnas))
    with gzip.open(outscnfile, 'w') as OUT:
        OUT.write(b'\n'.join(outscnrnas))

    return outfile, outscnfile


def list_of_files_from_dir(path, extension='.rna.fastq.gz'):
    import glob
    import os

    #targetPattern = r"C:\Test\*.txt"
    pattern = os.path.join(path, f'*{extension}')
    filepathlist = glob.glob(pattern)

    return filepathlist


def unique_sequence_names(contactfile):
    from natsort import natsorted

    with open(contactfile, 'r') as FILE:
        seqnames = [seqname for line in FILE for seqname in line.strip().split()[-1].split(',')]
        uniquenames = natsorted(list(set(seqnames)))

    return uniquenames


def rpkm(dcounts, dgff3):
    # dcounts is dictionary, key is feature ID, value is feature count
    # dgff3 is dictionar, key is feature ID, value is list of all elements in gff3
    # elem 3 is start coordinate and elem 4 is end coordinate of feature
    total = sum([v for v in list(dcounts.values())])
    # this is 'per million' scaling factor  # corrects seq depth
    scaledtotal = total / 1000000
    # larger values i.e. 1mill 'scale' count values to be larger (makes scaled total smaller)
    drpm = {k: v / scaledtotal for k, v in dcounts.items()}
    # drpkm is reads per kilobase per million reads
    # (int(dgff3[k][4]) - int(dgff3[k][3])) / 1000) == length of feature in kb
    drpkm = {k: v / ((int(dgff3[k][4]) - int(dgff3[k][3])) / 1000) for k, v in drpm.items()}

    return drpkm


def count_contacts_by_feature_norm_length(dgff3, contactfile, minrpkm=5):
    import math
    import os

    print('Collecting Number of connections per feature')
    dcounts = {}
    with open(contactfile, 'r') as FILE:
        for line in FILE:
            l = line.strip().split()
            support = int(math.ceil(float(l[-1])))  # round up to nearest integer # math.ceil returns float()
            rnascaff = 'scaffold51_%d' % int(l[0][2:])
            rnastart = int(l[1])
            dnascaff = 'scaffold51_%d' % int(l[3][2:])

            for k, v in dgff3.items():
                # only keep connections where RNA map to features included in dgff3 (ex: keep only ncRNA or only mRNA)
                if rnascaff == v[0]:
                    if (rnastart >= int(v[3])) and (rnastart <= int(v[4])):
                        dcounts[k] = dcounts.setdefault(k, 0) + support
    print('Normalizing number of connections per feature by feature length')
    drpkm = rpkm(dcounts, dgff3)
    outlines = ['\t'.join([k, str(v)]) for k, v in drpkm.items()]

    path, file = os.path.split(contactfile)
    outfile = os.path.join(path, file + '.GeneLengthNormCounts.tsv')
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(outlines))

    drpkm = {k: v for k, v in drpkm.items() if v > minrpkm}
    print(len(list(drpkm.keys())))
    outlines = ['\t'.join([k, str(v)]) for k, v in drpkm.items()]

    path, file = os.path.split(contactfile)
    outfile = os.path.join(path, file + f'.GeneLengthNormCounts.MinRPKM{minrpkm}.tsv')
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(outlines))

    return drpkm, outfile


def normalize_contact_counts(dgff3, dcounts, contactfile, side='DNA', buffer=0):
    import os
    from natsort import natsorted

    # here dgff3 is not indexed by scaffold
    print('normalize feature counts by \'rpkm\'')  # many IESs are 27bps so only scale by 100
    dcountsrpkm = rpkm(dcounts, dgff3)

    path, file = os.path.split(contactfile)
    outfile = os.path.join(path, '.'.join(file.split('.')[:-1]) + f'.{side}insidefeatures.buffer{buffer}.normcounts.txt')
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(['\t'.join([k, str(round(dcountsrpkm[k], 4))]) for k in natsorted(list(dcounts.keys()))]))

    return

def record_deletion_coordinates_from_cigar(CIGAR, leftmost_coordinate):
    # returns list of coordinates for gaps (splices, or deletions) in sequence
    import re

    gaps = []
    cigar = re.findall(r'\d+[A-Z]', CIGAR)
    if cigar == []:  # if there was a match # sometimes CIGAR string == * # this skips unmapped reads
        print(f'Provided CIGAR string: {CIGAR} does not match CIGAR pattern \\d+[A-Z]')
        rightmost_position = 0  # assumes unmapped read
    else:  # then read should be mapped
        rightmost_position = leftmost_coordinate  # subtract 1 because leftmost base is 1-based
        for i in cigar:
            if i[-1] in ['N', 'D']:
                gapleft = rightmost_position
                gapright = rightmost_position + int(i[:-1])
                gaps.append([gapleft, gapright - 1])  # subtract 1 from right because leftmost base is 1-based
                rightmost_position += int(i[:-1])
            elif i[-1] in ['M', 'X', '=']:
                rightmost_position += int(i[:-1])
            elif i[-1] in ['I', 'S', 'H', 'P']:
                pass
            else:
                pass

    return gaps


def record_match_coordinates_from_cigar(CIGAR, leftmost_coordinate):
    # returns list of coordinates for matching stretches in sequence
    # unfortunately Hisat2 includes mismatching and matching bases as M
    import re

    matches = []
    cigar = re.findall(r'\d+[A-Z]', CIGAR)
    if cigar == []:  # if there was a match # sometimes CIGAR string == * # this skips unmapped reads
        print(f'Provided CIGAR string: {CIGAR} does not match CIGAR pattern \\d+[A-Z]')
        rightmost_position = 0  # assumes unmapped read
    else:  # then read should be mapped
        rightmost_position = leftmost_coordinate  # subtract 1 because leftmost base is 1-based
        for i in cigar:
            if i[-1] in ['M', 'X', '=']:  # M is match or mismatch, X is mismatch, = is match
                left = rightmost_position
                right = rightmost_position + int(i[:-1])
                matches.append([left, right - 1])  # subtract 1 from right because leftmost base is 1-based
                rightmost_position += int(i[:-1])
            elif i[-1] in ['N', 'D']:
                rightmost_position += int(i[:-1])
            elif i[-1] in ['I', 'S', 'H', 'P']:
                pass
            else:
                pass

    return matches


def limit_contacts_by_feature_with_indices(dgff3, contactfile, dnabuffer=0, rnabuffer=0, rnamatch=10, maxgapoverlap=0.75):
    # dgff3 is two nested dictionaries outer dict key is scaffold name, outer dict value is dictionary containing all genes
    # inner dict key is feature ID, inner dict value is list of gff3 entries for gene feature
    # contact file is pt1 10 11 pt5 20 21 3 (delimited by spaces)
    # i.e. rnascaffold, rnaleft, rnaright, rnaorientation, rnacigar, dnascaffold, dnaleft, dnaright, dnaorientation, dnacigar, contactstrength, readname1,readname2,readname3
    # side == 'DNA' or 'RNA': molecule RNA-bridge-DNA, RNA contact position, or DNA contact position inside feature
    # dnabuffer allows extension of feature coordinates by specified amount start - dnabuffer & end + dnabuffer
    # rnabuffer extends the feature size by specified amount start + rnabuffer & end - rnabuffer (should be <20)
    # rnamatch is number of required bases to me match/mismatches to feature (def = 10)
    # maxgapoverlap is float 0-1. (def = 0.75). requires that RNA are NOT gapped over 75% of feature
    import os
    from natsort import natsorted

    dcounts = {}
    outlines = []
    featureids = []
    # check if rna (where RNA originated) is inside features
    rnascaffcolumn = 0  # number use 0-based numbering for columns
    rnaleftcolumn = 1
    rnarightcolumn = 2
    rnaorientationcolumn = 3
    rnacigarcolumn = 4
    # check if dna contact position (where RNA is associated with chromatin) is inside features
    dnascaffcolumn = 5  # number use 0-based numbering for columns
    dnaleftcolumn = 6
    dnarightcolumn = 7
    dnaorientationcolumn = 8
    dnacigarcolumn = 9


    with open(contactfile, 'r') as FILE:
        for count, line in enumerate(FILE):  # for each contact in contactfile
            support = int(line.strip().split()[-2])  # [-1] is comma separated reads names
            scaff = line.strip().split()[dnascaffcolumn][2:]  # isolate dna scaffold number only (keep as string)
            if line.strip().split()[dnaorientationcolumn] == '+':
                dna5prime = int(line.strip().split()[dnaleftcolumn])
            elif line.strip().split()[dnaorientationcolumn] == '-':
                dna5prime = int(line.strip().split()[dnarightcolumn])
            #dnascaff = 'scaffold51_%d' % int(l[3][2:])

            try:
                dgff3[scaff]
            except KeyError:
                pass  # if there are no features for this scaffold then skip it
            else:
                for k, v in dgff3[scaff].items():  # use gff3 index for scaffold to reduce runtime
                    # k is gene id, v is list of all gff3 entries for the feature (limited by scaffold)
                    # only keep connections where RNA (or DNA) contact positions map inside features included in dgff3
                    #if (rnascaff == v[0])
                    # require 5' dna end to be inside feature + buffer
                    featurestart, featureend = int(v[3]), int(v[4])
                    if (dna5prime >= (featurestart - dnabuffer)) and (dna5prime <= (featureend + dnabuffer)):
                        rnascaff = line.strip().split()[rnascaffcolumn][2:]
                        rnaleft = int(line.strip().split()[rnaleftcolumn])
                        rnaright = int(line.strip().split()[rnarightcolumn])
                        # find when rna mapping positions and feature positions overlap (to some extent)
                        if (rnascaff == scaff) and (rnaleft <= (featureend + rnabuffer)) and (rnaright >= (featurestart - rnabuffer)):
                            rnacigar = line.strip().split()[rnacigarcolumn]
                            matches = record_match_coordinates_from_cigar(rnacigar, rnaleft)
                            gaps = record_deletion_coordinates_from_cigar(rnacigar, rnaleft)
                            matchgate = 0  # if > 10 bases match/mismatching feature open the gate (make = 1)
                            gapgate = 1  # if read is gapped across 95% of feature, close gate (make = 0)
                            matchesinsidefeature = 0
                            gapsinsidefeature = 0
                            for match in matches:
                                # check every match in cigar and add up all matches if the overlap feature
                                # below calculates amount of overlap, if they don't overlap, returns 0
                                # len(range(max(x[0], y[0]), min(x[-1], y[-1]) + 1))
                                matchesinsidefeature += len(range(max(featurestart, match[0]), min(featureend, match[1]) + 1))
                            for gap in gaps:  # gap is list of len == 2, i.e. left and right coordinate of gap
                                gapsinsidefeature += len(range(max(featurestart, gap[0]), min(featureend, gap[1]) + 1))
                            if matchesinsidefeature >= rnamatch:
                                matchgate = 1  # open gate, more then 10 bases match/mismatch feature
                            if gapsinsidefeature / (featureend - featurestart + 1) > maxgapoverlap:
                                gapgate = 0  # close gate, gap covers more than 75% of feature. i.e. too many gaps

                            if (matchgate == 1) and (gapgate == 1):
                                dcounts[k] = dcounts.setdefault(k, 0) + support
                                outlines.append(' '.join([k, line.strip()]))
                                featureids.append(k)

            if count % 1000000 == 0:
                print(f'Count: {count} | Current contact: {line.strip()}')
    #dnumcontacts == number of unique contacts to the feature
    uniquefeatureids = natsorted(list(set(featureids)))
    print(f'Number of unique contacts inside features = {len(outlines)}')
    print(f'Number of unique features, with contacts inside = {len(uniquefeatureids)}')
    dnumcontacts = {k: featureids.count(k) for k in uniquefeatureids}

    path, file = os.path.split(contactfile)
    # this file has all contacts inside features
    outfile = os.path.join(path, '.'.join(file.split('.')[:-1]) + f'.insidefeatures.dnabuff{dnabuffer}.rnabuff{rnabuffer}.match{rnamatch}.gap{maxgapoverlap}.txt')
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(outlines))
    # number of lines in .counts.txt file = number of unique features with contacts inside
    outfilecounts = os.path.join(path, '.'.join(file.split('.')[:-1]) + f'.insidefeatures.dnabuff{dnabuffer}.rnabuff{rnabuffer}.match{rnamatch}.gap{maxgapoverlap}.counts.txt')
    with open(outfilecounts, 'w') as OUT:
        # output featureid \t total strength of all contacts in feature \t total unique contacts in feature
        OUT.write('\n'.join(['\t'.join([k, str(round(dcounts[k], 2)), str(dnumcontacts[k])]) for k in uniquefeatureids]))

    # send back dcounts so we can normalize the counts
    return outfile


def trim_features_by_scaffold_and_coordinates(dgff3, names, scaffold, start, end):
    dgff3trim = {}
    namestrim = []
    featurestarts = []
    featurelengths = []
    print('Trimming features to scaffold: %s, start: %d, end: %d' % (scaffold, start, end))
    for id in names:
        start_gff = int(dgff3[id][3])
        end_gff = int(dgff3[id][4])
        scaff = dgff3[id][0]
        # then the feature somehow overlaps scaffold interval
        # if (scaff == scaffold) and (start_gff <= end) and (start <= end_gff):
        # then gfffeature is completely inside scaffold interval
        if (scaff == scaffold) and (start_gff >= start) and (end_gff <= end):
            dgff3trim[id] = dgff3[id]
            namestrim.append(id)
            featurestarts.append(start_gff)
            featurelengths.append(end_gff - start_gff + 1)
            featurelengths.append(end_gff - start_gff + 1)

    return dgff3trim, namestrim, featurestarts, featurelengths


def read_gff3_index_by_scaffold(GFF3file, features=['all']):
    print('Reading Gff3 file and indexing by scaffold: %s' % GFF3file)
    flatfeatures = ', '.join(features)
    print(f'Features kept: {flatfeatures}')
    d = {}
    with open(GFF3file, 'r') as FILE:
        for line in FILE:
            if (line[0] == '#') or (line.strip() == ''):
                pass
            elif features == ['all']:  # keep all lines
                # key = line.strip().split('\t')[0] + ':' +  line.strip().split('\t')[3] + '_' + line.strip().split('\t')[4]
                #keyscaffold = '_'.join(line.strip().split('\t')[0].split('_')[:2])
                keyscaffold = line.strip().split('\t')[0].split('_')[1]  # make keyscaffold number only
                keyid = line.strip().split('\t')[8].split(';')[0][len('ID='):]
                try:
                    d[keyscaffold]
                except KeyError:
                    d[keyscaffold] = {keyid: line.strip().split('\t')}
                else:
                    d[keyscaffold][keyid] = line.strip().split('\t')
            elif line.strip().split('\t')[2] in features:  # keep lines only if in features
                # key = line.strip().split('\t')[0] + ':' + line.strip().split('\t')[3] + '_' + line.strip().split('\t')[4]
                keyscaffold = line.strip().split('\t')[0].split('_')[1]  # make keyscaffold number only
                keyid = line.strip().split('\t')[8].split(';')[0][len('ID='):]
                try:
                    d[keyscaffold]
                except KeyError:
                    d[keyscaffold] = {keyid: line.strip().split('\t')}
                else:
                    d[keyscaffold][keyid] = line.strip().split('\t')
    print('Number of scaffolds: %d' % len(list(d.keys())))
    print(list(d.keys())[:10])
    return d


def read_gff3(GFF3file, features=['all']):
    print('Reading Gff3 file: %s' % GFF3file)
    flatfeatures = ', '.join(features)
    print(f'Features kept: {flatfeatures}')
    d = {}
    names = []
    with open(GFF3file, 'r') as FILE:
        for line in FILE:
            if (line[0] == '#') or (line.strip() == ''):
                pass
            elif features == ['all']:  # keep all lines
                # key = line.strip().split('\t')[0] + ':' +  line.strip().split('\t')[3] + '_' + line.strip().split('\t')[4]
                key = line.strip().split('\t')[8].split(';')[0][len('ID='):]
                d[key] = line.strip().split('\t')
                names.append(key)
            elif line.strip().split('\t')[2] in features:  # keep lines only if in features
                # key = line.strip().split('\t')[0] + ':' + line.strip().split('\t')[3] + '_' + line.strip().split('\t')[4]
                key = line.strip().split('\t')[8].split(';')[0][len('ID='):]
                d[key] = line.strip().split('\t')
                names.append(key)
    print('Number of features: %d' % len(list(d.keys())))
    print(list(d.keys())[:10])
    return d, names


def seaborn_heatmap_with_barh(matrix, outpath, names, ys, starts, lengths):
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt

    print('Plotting heatmap with horizontal bar plot')

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True, gridspec_kw={'height_ratios': [1, 32]})
    sns.heatmap(matrix, ax=axes[1], cbar_kws={'label': 'Number of RNA-DNA Contacts', 'orientation': 'horizontal'})
    # axes[1].set(ylabel=None)
    axes[0].get_yaxis().set_visible(False)
    axes[0].barh(y=ys, width=lengths, left=starts)
    ends = [starts[i] + lengths[i] for i in range(len(starts))]
    for s, e in zip(starts, ends):
        axes[1].axvline(s, 0, len(matrix), linestyle='--', linewidth=0.25)
        axes[1].axvline(e, 0, len(matrix), linestyle='--', linewidth=0.25)
    if outpath != '':
        fig.savefig(outpath)
        print('Finished making heatmap: %s' % outpath)
    else:
        print('Finished making heatmap')

    plt.show()
    plt.close()

    return


def seaborn_heatmap(matrix, outpath=''):
    # matrix is numpy matrix, values in matrix are equal to intensity of all connections in genomic region
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt

    print('Plotting heatmap')

    fig, axes = plt.subplots(figsize=(8, 8))  # 2, 1, gridspec_kw={'height_ratios': [4, 1]}

    sns.heatmap(matrix, xticklabels=False, yticklabels=False, ax=axes)  # vmin=0, vmax=1,
    # axes[1].barh(y=[0,0], width=[0.2,0.4], left=[0,0.5])
    if outpath != '':
        fig.savefig(outpath)
        print('Finished making heatmap: %s' % outpath)
    else:
        print('Finished making heatmap')

    plt.show()
    plt.close()

    return


def np_matrix_sum_number_of_unique_xycoordinate_in_window(xs, ys, maxcoordinate, windowsize=10000, mincoordinate=0):
    import math
    import numpy as np

    # matrix size determined by max coordinate and window size
    # data in each matrix cell represents intensity in heatmap
    # coordinates are 0 based

    print('Maximum coordinate length: %d, Minimum coordinate length: %d, Window size: %d' % (maxcoordinate, mincoordinate, windowsize))
    # genome length / window size
    number_categories = int(math.ceil(float((maxcoordinate - mincoordinate) / windowsize)))  # number of x and y categories in matrix
    print(number_categories)
    # make matrix full of zeros
    matrix = np.zeros((number_categories, number_categories))
    print(matrix.shape)
    # fill in the matrix of zeros with intensity of connection based upon x, y coordinate
    for x, y in zip(xs, ys):
        if (x <= maxcoordinate) and (y <= maxcoordinate) and (x >= mincoordinate) and (y >= mincoordinate):
            xmatrix = math.floor((x - mincoordinate) / windowsize)  # round values down because zero based element numbering for np matrix
            ymatrix = math.floor((y - mincoordinate) / windowsize)  # so if y=500 and window=1000 then 0.5 rounds to zero because in zero cell
            matrix[xmatrix][ymatrix] += 1  # add one because contact has unique cooridnates, orientation, and CIGAR for RNA and DNA. Summed this value equals number of unique contacts in a window.

    matrixlog = np.log(matrix + 1)

    return matrix, matrixlog


def np_matrix_sum_intensity_if_xycoordinate_in_window(xs, ys, intensities, maxcoordinate, windowsize=10000, mincoordinate=0):
    import math
    import numpy as np

    # matrix size determined by max coordinate and window size
    # data in each matrix cell represents intensity in heatmap
    # coordinates are 0 based

    print('Maximum coordinate length: %d, Minimum coordinate length: %d, Window size: %d' % (maxcoordinate, mincoordinate, windowsize))
    # genome length / window size
    number_categories = int(math.ceil(float((maxcoordinate - mincoordinate) / windowsize)))  # number of x and y categories in matrix
    print(number_categories)
    # make matrix full of zeros
    matrix = np.zeros((number_categories, number_categories))
    print(matrix.shape)
    # fill in the matrix of zeros with intensity of connection based upon x, y coordinate
    for x, y, intensity in zip(xs, ys, intensities):
        if (x <= maxcoordinate) and (y <= maxcoordinate) and (x >= mincoordinate) and (y >= mincoordinate):
            xmatrix = math.floor((x - mincoordinate) / windowsize)  # round values down because zero based element numbering for np matrix
            ymatrix = math.floor((y - mincoordinate) / windowsize)  # so if y=500 and window=1000 then 0.5 rounds to zero because in zero cell
            matrix[xmatrix][ymatrix] += intensity

    matrixlog = np.log(matrix+1)

    return matrix, matrixlog


def seaborn_scatter_color_intensity(x, y, intensities, outpath=''):
    # x is list of x int values
    # y is list of y int values
    # intensities is list of float values between 0 and 1
    # outpath is full path to output file

    import seaborn as sns
    import matplotlib.pyplot as plt

    print('Plotting scatter plot colored by intensity')

    fig, axes = plt.subplots(figsize=(8, 8))
    cmap = sns.cubehelix_palette(start=.5, rot=-.5, as_cmap=True)  # lighter values have lower hue values
    sns.scatterplot(x, y, hue=intensities, ax=axes, palette=cmap)

    if outpath != '':
        fig.savefig(outpath)
        print('Finished making scatter plot colored by intensity: %s' % outpath)
    else:
        print('Finished making scatter plot colored by intensity')

    plt.show()
    plt.close()

    return


def plot_kde_heatmap(l1, l2, nbins=100, outpath=''):
    import statistics
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import kde

    x = np.array(l1)
    y = np.array(l2)
    # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
    k = kde.gaussian_kde([x, y])
    xi, yi = np.mgrid[x.min():x.max():nbins * 1j, y.min():y.max():nbins * 1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    # Make the plot
    plt.pcolormesh(xi, yi, zi.reshape(xi.shape))
    # Add color bar
    # plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.Greens_r)
    plt.colorbar()
    if outpath != '':
        plt.savefig(outpath)
    plt.show()
    plt.close()


def convert_scaff_coords_to_continuous_genome_coords(contactfile, dgcoords):
    # dgcoords is sum of all scaffold lengths 'below current scaffold
    # so if on scaff 3 then length of sacff 1 + length scaff 2
    import os
    import math

    print('Converting coordinates for each scaffold to continuous genome coordinates')

    newlines = []
    with open(contactfile, 'r') as FILE:
        for line in FILE:
            # rnascaffold rnaleft rnaright rnaorientation rnacigar dnascaffold dnaleft dnaright dnaorientation dnacigar AANNNAA strength
            l = line.strip().split()
            support = int(math.ceil(float(l[-1])))  # round up to nearest integer # math.ceil returns float()
            rnascaff = 'scaffold51_%d' % int(l[0][2:])
            dnascaff = 'scaffold51_%d' % int(l[5][2:])
            if l[3] == '+':  # if RNA is mapped in forward orientation
                rna3prime = int(l[2])
            elif l[3] == '-':  # if RNA read is mapped in reverse orientation
                rna3prime = int(l[1])
            if l[8] == '+':  # if DNA read is mapped in forward orientation
                dna5prime = int(l[6])
            elif l[8] == '-':  # if DNA read is mapped in reverse orientation
                dna5prime = int(l[7])
            temp = [rna3prime + dgcoords[rnascaff], dna5prime + dgcoords[dnascaff], support]
            newlines.append(temp)

    output = '\n'.join(['\t'.join(map(str, l)) for l in newlines])

    path, file = os.path.split(contactfile)
    outfile = os.path.join(path, '.'.join(contactfile.split('.')[:-1]) + '.gencoords.tsv')
    with open(outfile, 'w') as OUT:
        OUT.write(output)

    # transpose newlines to get list of three tuples: x's in one tuple and y's in another tuple and intensities in another [(xs), (ys), (intensities)]
    tlines = list(zip(*newlines))
    # get lists of x values and y values
    x, y, intensity = list(tlines[0]), list(tlines[1]), list(tlines[2])

    return outfile, x, y, intensity


def sum_fasta_lengths_by_sortlist(dlengths, sortlist):
    print('Summing lengths of all scaffolds')
    total = 0
    dnew = {}
    for k in sortlist:
        dnew[k] = total
        if total == 0:
            firstlength = dlengths[k]
        total += dlengths[k]
        # print('%s %d' % (k, total))
    return dnew, total, firstlength


def length_of_fasta_sequences(genomefile):
    import os

    print('Counting lengths of all scaffolds')
    path, f = os.path.split(genomefile)
    dgenome, names = read_fasta_as_dict(genomefile)
    d = {k: len(v) for k, v in dgenome.items()}

    return d, names


def summarize_raw_rna_dna_contacts_and_names(rawcontactfilenames):
    # input raw contact file with read names as last column
    # line structure rnascaffold, rnastartpos, rnaendpos, dnascaff, dnastart, dnaend, AANNN barcode, read name
    from natsort import natsorted
    import os
    import math
    path, f = os.path.split(rawcontactfilenames)

    print('Starting Summarize RNA DNA contacts, and sequence names')

    # record all sequence names for all unique contacts
    d = {}  # d is dict, key is unique rna_dna contact position 'rnascaff rnastart ...', value is list of all read names
    count = 0
    with open(rawcontactfilenames, 'r') as FILE:
        for l in FILE:  # if memory wasn't a problem i would do this
            # line structure rnascaffold, rnastartpos, rnaendpos, rnaorientation, rnacigar, dnascaff, dnastart, dnaend, dnaorientation, dnacigar, NNNbarcode, read name
            d.setdefault(' '.join(l.strip().split()[:-1]), []).append(l.strip().split()[-1])
            count += 1

    print(f'Number of raw contacts = {count}')
    print(f'Number of unique contacts = {len(list(d.keys()))}')
    print(f'Average contact strength = {count / len(list(d.keys()))}')

    print('writing out unique, but not sorted contacts')
    # write unique contacts, strength, and sequence names to file
    outpathnames = os.path.join(path, 'RNA.DNA.Contacts.Unique.Strength.AllSeqNames.txt')
    with open(outpathnames, 'w') as OUT:
        # line structure: below, AANNNAA column is omitted with PCRDupRemoval = True
        # rnascaffold rnaleft rnaright rnaorientation dnascaffold dnaleft dnaright dnaorientation AANNNAA strength all,read,names,delimited,by,commas
        OUT.write('\n'.join([' '.join([k] + [str(len(v)), ','.join(v)]) for k, v in d.items()]))

    print('sorting contacts')
    # count strength of each unique contact
    #sortcontacts = [i.split() for i in list(d.keys())]  # make list of lists
    #sortcontacts = natsorted(sortcontacts, key=lambda y: (y[0], y[1]))  # sort by first then second column
    #sortcontacts = natsorted([i.split() for i in list(d.keys())], key=lambda y: (y[0], y[1]))  # sort by first then second column
    #sortcontacts = natsorted(d)  # natsorted will sort keys of dict and make list
    #uniquecontacts = None  # clear memory
    #dcount = {}  # dcount is dict, key is unique rna_dna contact position, value is integer, strength of contact.
    #for l in sortcontacts:
    #    k = ' '.join(l)  # count strength of unique contacts
    #    dcount[k] = len(d[k])  # strength = number of read names per contact

    print('contacts sorted')
    # write unique contacts, strength, and sequence names to file
    outpathnames = os.path.join(path, 'RNA.DNA.Contacts.Unique.Strength.AllSeqNames.sort.txt')
    with open(outpathnames, 'w') as OUT:
        # line structure: below, AANNNAA column is omitted with PCRDupRemoval = True
        # rnascaffold rnastartpos rnaendpos dnascaff dnastart dnaend AANNNAA strength all,read,names,delimited,by,commas
        #OUT.write('\n'.join([' '.join(i + [str(len(d[' '.join(i)])), ','.join(d[' '.join(i)])]) for i in sortcontacts]))
        OUT.write('\n'.join([' '.join([i] + [str(len(d[i])), ','.join(d[i])]) for i in natsorted(d)]))  # sortcontacts
    #d = None  # clear memory

    # write unique contacts and strength to file
    outpath = os.path.join(path, 'RNA.DNA.Contacts.Unique.Strength.sort.txt')
    with open(outpath, 'w') as OUT:
        # line structure: below, AANNNAA column is omitted with PCRDupRemoval = True
        # rnascaffold rnastartpos rnaendpos rnaorientation rnacigar dnascaff dnastart dnaend dnaorientation dnacigar AANNNAA strength
        #OUT.write('\n'.join([' '.join(i + [str(len(d[' '.join(i)]))]) for i in sortcontacts]))  # sortcontacts
        OUT.write('\n'.join([' '.join([i] + [str(len(d[i]))]) for i in natsorted(d)]))  # sortcontacts

    # write log2(strength+1) to file for unique contacts & make file for circos (w/thickness)
    # take the log2 the strength of contacts (+1) to reduce high values to <100
    # dcount = {k: math.log(len(v)+1, 2) for k, v in d.items()}
    # the keys of d are now unique and while the thickness value will represent how often the contacts appeared
    # need to split the keys again so we can sort them for output
    outpathlog2 = os.path.join(path, 'RNA.DNA.Contacts.Unique.Strength.log2Plus1.sort.txt')
    with open(outpathlog2, 'w') as OUT:
        #OUT.write('\n'.join([' '.join(i + [str(math.log(len(d[' '.join(i)]) + 1, 2))]) for i in sortcontacts]))  # sortcontacts
        OUT.write('\n'.join([' '.join([i] + [str(math.log(len(d[i]) + 1, 2))]) for i in natsorted(d)]))  # sortcontacts
    outpathcircos = os.path.join(path, 'RNA.DNA.Contacts.Unique.Strength.log2Plus1.sort.circos.txt')
    with open(outpathcircos, 'w') as OUT:
        #OUT.write('\n'.join([' '.join(i + ['thickness=' + str(math.log(len(d[' '.join(i)]) + 1, 2)) + 'p']) for i in sortcontacts]))  # sortcontacts
        OUT.write('\n'.join([' '.join(i.split()[:3] + i.split()[5:8] + ['thickness=' + str(math.log(len(d[i]) + 1, 2)) + 'p']) for i in natsorted(d)]))  # sortcontacts
    #[i + ['thickness=' + str(dcount[' '.join(i)]) + 'p'] for i in sortcontacts]

    print('Finished summarizing RNA DNA contacts')

    return outpathnames, outpath, outpathlog2, outpathcircos


def record_rna_dna_contacts(fDNA, drna, ddna, sn, count=0, PCRDupRemoval=False):
    # fDNA is full path to dna bam file
    from natsort import natsorted
    import os

    sntrim = sn[0].lower() + sn.split(' ')[1][0]
    # rnaorderedkeys = natsorted(drna.keys(), alg=ns.IGNORECASE)
    # print(rnaorderedkeys[:10])
    output = []
    for k in drna.keys():
        try:
            ddna[k]
        except KeyError:
            pass
        else:
            # if the read maps in the RNA and DNA sam files
            for rnapos in drna[k]:  # the read could map to multiple locations 'equally' well (according to MQ cutoff)
                # print(rnapos)
                for dnapos in ddna[k]:  # the read could map to multiple locations 'equally' well (according to MQ)
                    # print(dnapos)
                    if PCRDupRemoval is False:
                        ##### THIS SPLIT STATEMENT IS NOT UNIVERSAL! ##### Used to isolate chromosome number
                        # '_'.join(rnapos[0].split('_')[1:])
                        # '_'.join(dnapos[0].split('_')[1:])
                        # rnascaff rnaleft rnaright rnaorientation rnacigar dnascaff dnaleft dnaright dnaorientation dnacigar readname
                        output.append([sntrim + rnapos[0].split('_')[1], rnapos[1], rnapos[2], rnapos[3], rnapos[4],
                                       sntrim + dnapos[0].split('_')[1], dnapos[1], dnapos[2], dnapos[3], dnapos[4],
                                       k])
                        # , 'thickness=' + str(len(ddna[k])) # str(math.log(len(drna[k]) + len(ddna[k]), 10))
                    elif PCRDupRemoval is True:
                        ##### THIS SPLIT STATEMENT IS NOT UNIVERSAL! ##### Used to isolate chromosome number
                        # k.split(':')[-1]] == the NNN PCR duplicate barcode (I attached to ends of read names with
                        # modified char_bridge_trackall.py script
                        output.append([sntrim + rnapos[0].split('_')[1], rnapos[1], rnapos[2], rnapos[3], rnapos[4],
                                       sntrim + dnapos[0].split('_')[1], dnapos[1], dnapos[2], dnapos[3], rnapos[4],
                                       k.split(':')[-1], k])  # k.split(':')[-1]]
                        # , 'thickness=' + str(len(ddna[k])) # str(math.log(len(drna[k]) + len(ddna[k]), 10))
    sortoutput = natsorted(output, key=lambda y: (y[0], y[1]))
    path, f = os.path.split(fDNA)
    outpath = os.path.join(path, 'RNA.DNA.Contacts.%s.raw.txt' % sntrim)
    outpathnames = os.path.join(path, 'RNA.DNA.Contacts.wNames.%s.raw.txt' % sntrim)
    # use output != [] to prevent writing only '\n' when no contacts present
    if (count == 0) and (output != []):
        with open(outpath, 'w') as OUT:
            OUT.write('\n'.join([' '.join(l[:-1]) for l in sortoutput]) + '\n')  # l[:-1] removes read name from end
        with open(outpathnames, 'w') as OUT:
            OUT.write('\n'.join([' '.join(l) for l in sortoutput]) + '\n')
    elif (count > 0) and (output != []):
        with open(outpath, 'a') as OUT:
            OUT.write('\n'.join([' '.join(l[:-1]) for l in sortoutput]) + '\n')  # l[:-1] removes read name from end
        with open(outpathnames, 'a') as OUT:
            OUT.write('\n'.join([' '.join(l) for l in sortoutput]) + '\n')

    return outpath, outpathnames


def get_rightmost_reference_based_alignment_coordinate(CIGAR, leftmost_coordinate):
    import re
    cigar = re.findall(r'\d+[A-Z]', CIGAR)
    if cigar == []:  # if there was a match # sometimes CIGAR string == * # this skips unmapped reads
        print(f'Provided CIGAR string: {CIGAR} does not match CIGAR pattern \\d+[A-Z]')
        rightmost_position = 0  # assumes unmapped read
    else:  # then read should be mapped
        rightmost_position = leftmost_coordinate - 1  # subtract 1 because leftmost base is 1-based
        for i in cigar:
            if i[-1] in ['M', 'N', 'D', 'X', '=']:
                rightmost_position += int(i[:-1])
            elif i[-1] in ['I', 'S', 'H', 'P']:
                pass
            else:
                pass

    return rightmost_position


def record_read_left_right_orientation(bamfile, MQ, dreads={}):
    #bamfile is full path to bam file
    #MQ = int() = minimum mapping quality threshold
    #dreads is class object utilizing
    # end specifies if the 5' or 3' ends are kept. '5' == 5 prime end and '3' == 3 prime end
    # We want 5' end of DNA sequences and 3' ends of RNA sequences
    # if f'{flag:012b}'[-5] == '1':  # if 1 then it is reverse # '000000010000' # 16 # then the read is mapped in reverse
    # dreads key is readname, value is scaffold leftmost_coordinate rightmost_coordinate orientation
    import pysam

    # print(f'{bamfile}')
    with pysam.AlignmentFile(bamfile, 'rb') as FILE:
        for line in FILE:
            line = line.tostring()  # changes the pysam formatting stuff to the actual format of the sam file
            if line[0] == '@':
                pass
            elif int(line.strip().split('\t')[4]) >= MQ:
                flag = int(line.strip().split('\t')[1])  # sam flag
                # if () and (f'{flag:012b}'[-9] == '0')
                if f'{flag:012b}'[-3] == '0':  # if the read is mapped # '000000000100' == 4 == read unmapped #[-3] != 1
                    if f'{flag:012b}'[-5] == '1':  # if 1 then it is reverse # '000000010000' == 16 read mapped reverse
                        # key is read name, value is scaffold name and then 5' or 3' position of read mapped
                        # read name is reduced to final three fields which should be a unique identifier
                        # scaffold name is reduced to only scaffold number
                        # '_'.join(line.strip().split('\t')[0].split(':')[-3:])
                        orientation = '-'
                    else:  # read is in forward orientation
                        orientation = '+'
                    readname = line.strip().split('\t')[0]
                    scaffold = line.strip().split('\t')[2]
                    leftcoord = int(line.strip().split('\t')[3])
                    cigar = line.strip().split('\t')[5]
                    rightcoord = get_rightmost_reference_based_alignment_coordinate(cigar, leftcoord)
                    dreads.setdefault(readname, []).append([scaffold, str(leftcoord), str(rightcoord), orientation, cigar])

    return dreads


def read_fasta_as_dict(f):
    d = {}  # fasta names are key, values are sequences as string
    namelist = []
    with open(f, 'r') as FILE:
        for line in FILE:
            if line[0] == '>':
                if ' ' in line:
                    name = line.strip().split()[0][1:-len('_with_IES')]
                    namelist.append(name)
                    d[name] = []
                else:
                    name = line.strip()[1:]
                    namelist.append(name)
                    d[name] = []
            elif line.strip() != '':  # else: # trying to prevent some crap happening on the last line
                d[name].append(line.strip())
    for name in namelist:
        d[name] = ''.join(d[name])  # join list of partial sequences. Useful if interleaved fasta

    return d, namelist


def run_samtools(samfile):
    # samfile is full path to sam file
    print('Starting samtools')
    import os
    import subprocess

    path, sam = os.path.split(samfile)
    bamfile = '.'.join(samfile.split('.')[:-1] + ['bam'])
    sortbamfile = '.'.join(samfile.split('.')[:-1] + ['sort', 'bam'])
    # sortsamfile = '.'.join(samfile.split('.')[:-1] + ['sort', 'sam'])
    # flagstatfile = '.'.join(samfile.split('.')[:-1] + ['sort', 'sam', 'flagstat'])

    with open(bamfile, 'w') as OUT:
        # cmd = 'samtools view -h -b %s > %s' % (samfile, bamfile)
        cmd = 'samtools view -h -b %s' % (samfile)
        ps = subprocess.Popen(cmd.split(), stdout=OUT)
        ps.wait()
    # cmd = 'rm %s' % samfile
    os.remove(samfile)  # delete sam file
    with open(sortbamfile, 'w') as OUT:
        # cmd = 'samtools sort %s > %s' % (bamfile, sortbamfile)
        cmd = 'samtools sort %s' % (bamfile)
        ps = subprocess.Popen(cmd.split(), stdout=OUT)
        ps.wait()
    os.remove(bamfile)  # delete bam file
    cmd = 'samtools index %s' % sortbamfile
    ps = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    ps.wait()
    # with open(sortsamfile, 'w') as OUT:
    #     # cmd = 'samtools view -h %s > %s' % (sortbamfile, sortsamfile)
    #     cmd = 'samtools view -h %s' % (sortbamfile)
    #     ps = subprocess.Popen(cmd.split(), stdout=OUT)
    #     ps.wait()
    # with open(flagstatfile, 'w') as OUT:
    #     # cmd = 'samtools flagstat %s > %s' % (sortsamfile, flagstatfile)
    #     cmd = 'samtools flagstat %s' % (sortsamfile)
    #     ps = subprocess.Popen(cmd.split(), stdout=OUT)
    #     ps.wait()

    print('Finished with Samtools\n')
    return sortbamfile  # , sortsamfile


def run_hisat2(aligndatabase, fastqfile, alignargs=''):
    print('Starting Hisat2: aligning\n%s' % fastqfile)
    import os
    import subprocess

    path, f = os.path.split(fastqfile)
    pathminusonedir, dir = os.path.split(path)
    make_directory(os.path.join(pathminusonedir, 'hisat2'))
    outsamfile = os.path.join(pathminusonedir, 'hisat2', '.'.join(f.split('.')[:-2] + ['sam']))  # assuming .fastq.gz
    cmd = 'hisat2 -q %s -x %s -U %s -S %s' % (alignargs, aligndatabase, fastqfile, outsamfile)
    subprocess.call(cmd.split())

    print('Finished with Hisat2\n')
    return outsamfile


def divide_number_of_lines_in_files_mult_100(f1, f2, expectedhuman=1.0):
    # divide number of lines in f1 by f2 multiply by 100 (f1/f2)*100.00
    # expected = expected percent of human contamination
    c1, c2 = 0, 0
    with open(f1, 'r') as FILE:
        for line in FILE:
            c1 += 1.0
    with open(f2, 'r') as FILE:
        for line in FILE:
            c2 += 1.0
    # percent of human RNA from all reads (human + species of interest)
    # only includes RNA that pass char-bridge filters i.e. RNA min length 15 and DNA min length 15 etc.
    percent = (c1 / (c2 / 2)) * 100.00  # divide c2 by 2 because counting line in rna fasta file
    # this is percent of human free-floating RNA that should represent percent free-floating RNA in Species of Interest
    expectedfreefloatingRNA = (percent / expectedhuman) * 100  # if 0.01 / 1.0 * 100 then = 1%

    return percent, expectedfreefloatingRNA


def run_blastn_match_db(fastafile, database, outformat=6, percentidentitiythreshold=98.00, bestbitscoreonly=True):
    # fasta file is full path to .fasta query file
    # database is full path to blastn database
    # calculate % of human RNA spike-in by with BLASTn # this represents expected % free-floating RNA in sample
    print('Start BLASTn on file:\n%s\nTo Database:\n%s\n' % (fastafile, database))
    import subprocess
    import os

    path, f = os.path.split(fastafile)
    pathminusonedir, dir = os.path.split(path)
    outpath = os.path.join(pathminusonedir, 'blastn')
    make_directory(outpath)
    if bestbitscoreonly == True:
        # takes the best BLAST result by bit score
        outfile = 'best_bit_score_per_query.blastn.RNA.tsv'
        fulloutpath = os.path.join(outpath, outfile)
        # full cmd = 'blastn -query %s -db %s -outfmt %d | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > %s' % (fastafile, database, outformat, fulloutpath)
        cmdpipe = ['blastn -query %s -db %s -outfmt %d' % (fastafile, database, outformat), 'sort -k1,1 -k12,12nr -k11,11n %s' % (fulloutpath + '.blastn'), 'sort -u -k1,1 --merge %s' % (fulloutpath + '.sort1')]
        for count, cmd in enumerate(cmdpipe):
            if count == 0:
                with open(fulloutpath + '.blastn', 'w') as OUT:
                    print('Pipe step 1.1')
                    print(cmd)
                    ps = subprocess.Popen(cmd.split(), bufsize=-1, stdout=OUT)
                    print('Pipe step 1.2')
                    ps.wait()
                    print('Pipe step 1.3')
            elif count != len(cmdpipe) - 1:  # if it is not the last command
                with open(fulloutpath + '.sort1', 'w') as OUT:
                    print('Pipe step 2.1')
                    ps = subprocess.Popen(cmd.split(), bufsize=-1, stdout=OUT)
                    print('Pipe step 2.2')
                    ps.wait()
                    print('Pipe step 2.3')
            else:  # it must be the last command
                with open(fulloutpath, 'w') as OUT:
                    print('Pipe step 3.1')
                    ps = subprocess.Popen(cmd.split(), bufsize=-1, stdout=OUT)
                    print('Pipe step 3.2')
                    ps.wait()
                    print('Pipe step 3.3')
        #cmdpipe = ['blastn -query %s -db %s -outfmt %d' % (fastafile, database, outformat), 'sort -k1,1 -k12,12nr -k11,11n', 'sort -u -k1,1 --merge -']
        #for count, cmd in enumerate(cmdpipe):
        #    if count == 0:
        #        print('Pipe step 1.1')
        #        print(cmd)
        #        ps = subprocess.Popen(cmd.split(), bufsize=-1, stdout=subprocess.PIPE)
        #        print('Pipe step 1.2')
        #        ps.wait()
        #        print('Pipe step 1.3')
        #    elif count != len(cmdpipe)-1:  # if it is not the last command
        #        print('Pipe step 2.1')
        #        ps = subprocess.Popen(cmd.split(), stdin=ps.stdout, stdout=subprocess.PIPE)
        #        print('Pipe step 2.2')
        #        ps.wait()
        #        print('Pipe step 2.3')
        #    else:  # it must be the last command
        #        with open(fulloutpath, 'w') as OUT:
        #            print('Pipe step 3.1')
        #            ps = subprocess.Popen(cmd.split(), stdin=ps.stdout, stdout=OUT)
        #            print('Pipe step 3.2')
        #            ps.wait()
        #            print('Pipe step 3.3')

    else:
        ### won't work right now
        outfile = 'All_hits_per_query.blastn.RNA.tsv'
        fulloutpath = os.path.join(outpath, outfile)
        cmd = 'blastn -query %s -db %s -outfmt %d > %s' % (fastafile, database, outformat, fulloutpath)
        cmdlist = cmd.split()
        subprocess.call(cmdlist)
    outdbmatching = os.path.join(outpath, 'human.rna.freefloating.tsv')
    #cmd = 'awk -F \"\t\" \'$3 > %f {print $1}\' %s > %s' % (percentidentitiythreshold, fulloutpath, outdbmatching)
    #subprocess.call(cmd.split())  # doesnt work with ' characters somehow
    with open(fulloutpath, 'r') as FILE:
        output = [line.strip() for line in FILE if float(line.split('\t')[2]) > percentidentitiythreshold]
    with open(outdbmatching, 'w') as OUT:
        OUT.write('\n'.join(output))
    print('Finished BLASTn, output results to:\n%s\n%s\n' % (fulloutpath, outdbmatching))

    return fulloutpath, outdbmatching


def fastq_to_fasta(fullpath):  # fullpath = string, full path to file
    """ Converts a .fastq file to fasta file """
    print('converting fastq to fasta')
    import os

    path, f = os.path.split(fullpath)
    fout = os.path.join(path, '.'.join(f.split('.')[:-1]) + '.fa')
    if os.path.exists(fout):  # if the output file exists delete it before appending
        os.remove(fout)

    FILE = open(fullpath, 'r')
    OUT = open(fout, 'a')

    line = FILE.readline()
    count = 0
    while line:
        if count % 4 == 0:  # if the entire line is only a +
            o = '>' + line[1::].strip() + '\n'
            OUT.write(o)
        elif count % 4 == 1:
            o = line.strip() + '\n'
            OUT.write(o)
        line = FILE.readline()
        count += 1

    FILE.close()
    OUT.close()

    print('Converted fastq to fasta:\n%s\n%s\n' % (fullpath, fout))

    return fout


def run_char_bridge(filegz, rnalen=18, dnalen=18, rnamaxlen=1000):
    # cd /media/sf_LinuxShare/Projects/Lyna/flypipe
    # python char_bridge_trackall.py --FASTQGZ 500_L1_R1R2.trim.AssUnFUnR.fastq.gz --NAME 500_R1R2.trim.AssUnFUnR.DefaultBridge.DNA15RNA15.CodeMod. --minRNA 15 --minDNA 15
    print('Start Charseq Bridge Removal on file:\n%s' % filegz)
    import subprocess
    import os

    path, f = os.path.split(filegz)
    pathminusonedir, dir = os.path.split(path)
    path_to_char_bridge = os.path.join(pathminusonedir, 'char_bridge_trackall_simple.py')
    outprefix = '.'.join(filegz.split('.')[:-2] + ['DNA%dRNA%dRNAmax%d.' % (rnalen, dnalen, rnamaxlen)])
    cmd = "python2 %s --FASTQGZ %s --NAME %s --minRNA %d --minDNA %d --maxRNA %d" % (path_to_char_bridge, filegz, outprefix, rnalen, dnalen, rnamaxlen)
    print(cmd)
    cmdlist = cmd.split()
    subprocess.call(cmdlist)
    rnafile, dnafile = outprefix + 'rna.fastq.gz', outprefix + 'dna.fastq.gz'
    print('Finished Charseq Bridge Removal, output files:\n%s\n%s\n' % (rnafile, dnafile))

    return rnafile, dnafile


def file_splitter(f, chunk=1200000, gz=False):
    # f = str, full path to file
    # chunk = int ## chunk is number of lines per split file, default 400000, multiples of 4 work well for fastq files
    # f is full path to file # limit is line number limit 1st line is 1 NOT 0 # still outputs the linecount <= limit

    if gz is False:
        with open(f, 'r') as FILE:
            splitcount = 0
            splitfilenames = []
            outlines = []
            for count, line in enumerate(FILE, start=1):
                if count % chunk == 0:
                    # then append last line and split file
                    outlines.append(line)
                    outf = '.'.join(f.split('.')[:-1] + [f'split{splitcount}', f.split('.')[-1]])
                    splitfilenames.append(outf)
                    with open(outf, 'w') as OUT:
                        OUT.write(''.join(outlines))
                    outlines = []
                    splitcount += 1
                elif count % chunk < chunk:
                    outlines.append(line)
        if len(outlines) > 0:  # if some lines were not written to file b/c len(outlines) < chunk value
            outf = '.'.join(f.split('.')[:-1] + [f'split{splitcount}', f.split('.')[-1]])
            splitfilenames.append(outf)
            with open(outf, 'w') as OUT:
                OUT.write(''.join(outlines))
            outlines = []

    elif gz is True:
        import gzip
        with gzip.open(f, 'rb') as FILE:
            splitcount = 0
            splitfilenames = []
            outlines = []
            for count, line in enumerate(FILE, start=1):
                if count % chunk == 0:
                    # then append last line and split file
                    outlines.append(line)
                    outf = '.'.join(f.split('.')[:-2] + [f'split{splitcount}', '.'.join(f.split('.')[-2:])])
                    splitfilenames.append(outf)
                    with gzip.open(outf, 'wb') as OUT:
                        OUT.write(b''.join(outlines))
                    outlines = []
                    splitcount += 1
                elif count % chunk < chunk:
                    outlines.append(line)
        if len(outlines) > 0:  # if some lines were not written to file b/c len(outlines) < chunk value
            outf = '.'.join(f.split('.')[:-2] + [f'split{splitcount}', '.'.join(f.split('.')[-2:])])
            splitfilenames.append(outf)
            with gzip.open(outf, 'wb') as OUT:
                OUT.write(b''.join(outlines))
            outlines = []

    else:
        print('gz parameter incorrectly specified, nothing happened')

    print(f)
    print(f'Split into {splitcount+1} files')
    print('Finished')

    return splitfilenames, splitcount+1


def merge_filelist_to_one_file(filelist, mergedoutfile, gzip=False):
    # the forward and reverse reads are still fine to use if they simply have a bigger insert size and can not be merged together, so combine three files to one file
    # cat /media/sf_LinuxShare/Projects/Lyna/DATA/500_L1_R1R2.trim.assembled.fastq > /media/sf_LinuxShare/Projects/Lyna/flypipe/500_L1_R1R2.trim.AssUnFUnR.fastq
    # cat /media/sf_LinuxShare/Projects/Lyna/DATA/500_L1_R1R2.trim.unassembled.forward.fastq >> /media/sf_LinuxShare/Projects/Lyna/flypipe/500_L1_R1R2.trim.AssUnFUnR.fastq
    # cat /media/sf_LinuxShare/Projects/Lyna/DATA/500_L1_R1R2.trim.unassembled.reverse.fastq >> /media/sf_LinuxShare/Projects/Lyna/flypipe/500_L1_R1R2.trim.AssUnFUnR.fastq
    print('Start merging files')
    if gzip is True:
        import gzip

        with gzip.open(mergedoutfile, 'wb') as OUT:
            OUT.write('')  # write nothing to output file to clear any previous information from file
        with gzip.open(mergedoutfile, 'ab') as OUT:
            for f in filelist:
                with gzip.open(f, 'rb') as FILE:
                    outlines = []
                    for count, line in enumerate(FILE):
                        if count % 10000 < 10000:
                            outlines.append(line)  # OUT.write(''.join(FILE.readlines())) # inefficient memory
                        elif count % 10000 == 0:
                            outlines.append(line)
                            OUT.write(b''.join(outlines))
                            outlines = []
                    if len(outlines) > 0:
                        OUT.write(b''.join(outlines))
                        outlines = []
    else:  # assume gz is false
        with open(mergedoutfile, 'w') as OUT:
            OUT.write('')  # write nothing to output file to clear any previous information from file
        with open(mergedoutfile, 'a') as OUT:
            for f in filelist:
                with open(f, 'r') as FILE:
                    outlines = []
                    for count, line in enumerate(FILE):
                        if count % 10000 < 10000:
                            outlines.append(line)  # OUT.write(''.join(FILE.readlines())) # inefficient memory
                        elif count % 10000 == 0:
                            outlines.append(line)
                            OUT.write(''.join(outlines))
                            outlines = []
                    if len(outlines) > 0:
                        OUT.write(''.join(outlines))
                        outlines = []

    print('Finished merging files:\n%s\n' % mergedoutfile)


def run_pear(forwardout, reverseout):
    # pear -f 500_LK_R1.trim.fastq.gz -r 500_LK_R2.trim.fastq.gz -o 500_L1_R1R2.trim
    import os
    import subprocess

    print('Starting PEAR:\n%s\n%s' % (forwardout, reverseout))
    path, f = os.path.split(forwardout)
    pathminusonedir, dir = os.path.split(path)
    outpath = os.path.join(pathminusonedir, 'pear')
    make_directory(outpath)
    pearoutfileprefix = '_'.join(f.split('_')[:-1] + ['R1R2.trim'])
    fulloutpath = os.path.join(outpath, pearoutfileprefix)
    cmd = 'pear -f %s -r %s -o  %s' % (forwardout, reverseout, fulloutpath)
    print(cmd)
    subprocess.call(cmd.split())
    print('Finished PEAR, output at:\n%s\n' % fulloutpath)

    return pearoutfileprefix, fulloutpath


def run_fastqc(fullpath):
    # fullpath is full path to input file
    # '/media/sf_LinuxShare/Programs/FastQC/fastqc -o /media/sf_LinuxShare/Projects/Lyna/DATA/fastqc -f fastq fastq 200107_NB501850_A_L1-4_ADPF-98_R1.fastq'
    print('Starting fastqc')
    import os
    import subprocess

    path, f = os.path.split(fullpath)
    pathminusonedir, dir = os.path.split(path)
    outpath = os.path.join(pathminusonedir, 'fastqc')
    print(outpath)
    make_directory(outpath)
    # /media/sf_LinuxShare/Programs/FastQC/fastqc is path to executable
    cmd = '/media/sf_LinuxShare/Programs/FastQC/fastqc -o %s -f fastq fastq %s' % (outpath, fullpath)
    print(cmd)
    subprocess.call(cmd.split())
    outputfile = os.path.join(outpath, f)
    print('Finished fastqc, output at directory:\n%s\n' % outpath)

    return


def run_gzip_compress(fullpath, cleanup=False):
    import subprocess

    # pack files (using full path)
    print('Start packing files:\n%s' % fullpath)
    if cleanup is True:
        cmd = "gzip -f %s" % fullpath
    else:
        cmd = "gzip -f -k %s" % fullpath
    print(cmd)
    cmdlist = cmd.split()
    subprocess.call(cmdlist)
    filegz = fullpath + '.gz'
    print('Finished packing:\n%s\n' % filegz)

    return filegz


def run_gzip_decompress(fullpath):
    import subprocess
    # unpack fastp trimmed file (using full path)
    print('Start unpacking trimmed.gz file:\n%s' % fullpath)
    cmd = "gzip -f -d -k %s" % fullpath
    print(cmd)
    cmdlist = cmd.split()
    subprocess.call(cmdlist)
    filetrim = fullpath[:-len('.gz')]
    print('Finished unpacking:\n%s\n' % filetrim)

    return filetrim


def run_fastp(forwardreadfile, reversereadfile):
    import subprocess
    import os

    ### start fastp ###
    path1, f1 = os.path.split(forwardreadfile)
    path2, f2 = os.path.split(reversereadfile)

    print('Start trimming of:\n%s\n%s' % (forwardreadfile, reversereadfile))
    make_directory(os.path.join(path1, 'fastp'))
    make_directory(os.path.join(path2, 'fastp'))
    # example: cd /media/sf_LinuxShare/Projects/Lyna/DATA
    #cmd = "cd %s" % path1
    #print(cmd)
    #cmdlist = cmd.split()
    #p = subprocess.call(cmdlist)

    ## fastp -i 500_LK_L1_R1.fastq.gz -I 500_LK_L1_R2.fastq.gz -o 500_LK_R1.trim.fastq.gz -O 500_LK_R2.trim.fastq.gz
    ## fastp -i /media/sf_LinuxShare/Projects/Lyna/TestMyPipe/500_LK_L1_R1.fastq.gz -I /media/sf_LinuxShare/Projects/Lyna/TestMyPipe/500_LK_L1_R2.fastq.gz -o /media/sf_LinuxShare/Projects/Lyna/TestMyPipe/fastp/500_LK_L1_R1.fastq.gz -O /media/sf_LinuxShare/Projects/Lyna/TestMyPipe/fastp/500_LK_L1_R2.fastq.gz
    forwardout = os.path.join(path1, "fastp", '.'.join(f1.split('.')[:-2] + ['trim', 'fastq', 'gz']))  # output file for
    reverseout = os.path.join(path2, "fastp", '.'.join(f2.split('.')[:-2] + ['trim', 'fastq', 'gz']))  # output file rev
    cmd = "fastp -i %s -I %s -o %s -O %s" % (forwardreadfile, reversereadfile, forwardout, reverseout)
    print(cmd)
    cmdlist = cmd.split()
    p = subprocess.Popen(cmdlist, stdout=subprocess.PIPE)
    cmdout, err = p.communicate()
    print(cmdout)
    print('Finished trimming, files output to:\n%s\n%s\n' % (forwardout, reverseout))
    ### end fastp ###

    return forwardout, reverseout


def make_directory(dirName):
    # dirName is directory in cwd or full path to directory
    import os

    if not os.path.exists(dirName):
        os.mkdir(dirName)
        print("Directory ", dirName, " Created ")
    else:
        print("Directory ", dirName, " already exists")


def main(forwardreadfile, reversereadfile, blastdatabase, aligndatabase, altreadfile1='/path/to/altreads1.fastq',
         altreadfile2='/path/to/altreads2.fastq', steps=['fastp', 'fastqc', 'pear', 'mergefiles', 'split2chunks',
         'charbridge', 'fastq2fasta', 'blast', 'hisat2', 'samtools', 'rnadnacontacts', 'heatmap',
         'contacts in features', 'get sequences', 'blast rna', 'cumulative distribution'], rnalen=20, dnalen=20,
         rnamaxlen=1000, REsequence='GATC', windowsize=300, contactcutoff=1, MQ=1, dnabuffer=0, rnabuffer=0,
         rnamatch=10, maxgapoverlap=0.75, sn='Paramecium tetraurelia', alignargs='',
         genomefile='/media/sf_LinuxShare/Ciliates/Genomes/Seqs/ptetraurelia_mac_51.fa', PCRDupRemoval=False,
         annot_file='/media/sf_LinuxShare/Ciliates/Genomes/Annotations/ptetraurelia_mac_51_annotation_v2.0.gff3',
         zoomscaffname='', zoomwinsize=20, scaffcoordmin=[0], scaffcoordmax=[0], features=['mRNA']):

    # forwardreadfile = '/media/sf_LinuxShare/Projects/Lyna/TestMyPipe/500_LK_L1_R1.fastq.gz'
    # reversereadfile = '/media/sf_LinuxShare/Projects/Lyna/TestMyPipe/500_LK_L1_R2.fastq.gz'
    # rnalen = 15
    # dnalen = 15
    # REsequence = 'GATC'
    # windowsize = 300
    # contactcutoff = 2
    # MQ = 10
    # buffer = number of bp to extend feature coordinates
    # sn = 'Paramecium tetraurelia'
    # blastdatabase = '/media/sf_LinuxShare/Humans/Genome/Seqs/GRCh38_top_level.fa'
    # aligndatabase = '/media/sf_LinuxShare/Ciliates/Genomes/Hisat2_Indexes/Pt_51_Mac'
    # 'normalize', normalizing for DPNii sites is removed
    # main(forwardreadfile, reversereadfile, rnalen, dnalen, REsequence, windowsize, contactcutoff, MQ, sn,
    # blastdatabase, aligndatabase)

    # run script from directory with flypipe scripts and data
    # assumes data is paired end
    # forwardreadfile == full path to forward read file .fastq.gz
    # reversereadfile == full path to reverse read file .fastq.gz
    import os
    import sys
    import math
    from namedlist import namedlist
    from natsort import natsorted, ns
    # insert at 1, 0 is the script path (or '' in REPL)
    ## directory, forward = os.path.split(forwardreadfile)
    ## sys.path.insert(1, directory)  # insert into path so we can find the RecordRNADNAContacts script
    ## import RecordRNADNAContacts

    if 'fastp' in steps:
        ### start fastp ###
        forwardout, reverseout = run_fastp(forwardreadfile, reversereadfile)

    if 'fastqc' in steps:
        ### start gzip -dk file.trim.fastq.gz
        forwardtrim = run_gzip_decompress(forwardout)
        reversetrim = run_gzip_decompress(reverseout)

        ### start fastqc ###
        run_fastqc(forwardtrim)
        run_fastqc(reversetrim)

    if 'pear' in steps:
        ### start PEAR ###
        pearfileprefix, pearoutpath = run_pear(forwardout, reverseout)

    if 'mergefiles' in steps:
        ### merge files: merged "assembled" PE reads, unmerged forward reads, and unmerged reverse reads to 1 file / sample
        filelist = [pearoutpath + '.assembled.fastq', pearoutpath + '.unassembled.forward.fastq', pearoutpath + '.unassembled.reverse.fastq']
        mergedoutfile = pearoutpath + '.AssUnFUnR.fastq'
        merge_filelist_to_one_file(filelist, mergedoutfile)

        ### start gzip to compress merged file ###
        filegz = run_gzip_compress(mergedoutfile, cleanup=True)

    if 'split2chunks' in steps:
        if steps[0] == 'split2chunks':
            filegz = altreadfile1  ## assuming user wants to start at this step

        splitgzfiles, numsplitfiles = file_splitter(filegz, chunk=4000000, gz=True)

    if 'charbridge' in steps:
        if steps[0] == 'charbridge':
            splitgzfiles = altreadfile1  # this is list of full paths to split files

        rnagzfiles, dnagzfiles = [], []
        for count, filegz in enumerate(splitgzfiles):
            ### start char_bridge_trackall.py (uses charbridgetools) to remove bridge ###
            rnafile, dnafile = run_char_bridge(filegz, rnalen, dnalen, rnamaxlen)
            rnagzfiles.append(rnafile)
            dnagzfiles.append(dnafile)

    if 'fastq2fasta' in steps:
        if steps[0] == 'fastq2fasta':
            dnafile = altreadfile1  ## assuming user wants to start at this step
            rnafile = altreadfile2
        ### start gzip decompress ###
        rnafiletrim = run_gzip_decompress(rnafile)  # rnafiletrim == decompressed RNA file .fastq
        dnafiletrim = run_gzip_decompress(dnafile)

        ### convert fastq to fasta for BLASTn ###
        rnafastafile = fastq_to_fasta(rnafiletrim)

    if 'blast' in steps:
        ### calculate percent free-floating RNA w/ BLAST & human spike in ###
        print('Start calculations of percent free-floating human RNA contamination')
        print('Expected percent of human-spiked-in of total RNA is 1 percent')
        blastoutpath, blastoutdbmatching = run_blastn_match_db(rnafastafile, blastdatabase, 6, 98.00, True)  # fullpath to fasta file, outformat == 6, percent identity threshold to determine if read is from a species, BestBitScoreOnly?? == True
        percent, expectedfreefloatingRNA = divide_number_of_lines_in_files_mult_100(blastoutdbmatching, rnafastafile)  # blastoutpath
        print('Measured percent of free-floating human RNA results: %f' % percent)
        print('Expected percent of free-floating RNA from %s \"contaminating\" results: %f' % (sn, expectedfreefloatingRNA))

    if 'hisat2' in steps:
        if steps[0] == 'hisat2':
            ## assuming user wants to start at this step ## assumes altreadfile1&2 are lists of files
            dnagzfiles = altreadfile1
            rnagzfiles = altreadfile2

        ### align reads to genome ###
        dnasamfiles, rnasamfiles = [], []
        for dnafile, rnafile in zip(dnagzfiles, rnagzfiles):
            rnasam = run_hisat2(aligndatabase, rnafile, alignargs)
            dnasam = run_hisat2(aligndatabase, dnafile, alignargs)
            rnasamfiles.append(rnasam)
            dnasamfiles.append(dnasam)

    if 'samtools' in steps:
        ### samtools ###
        dnasortbamfiles, rnasortbamfiles = [], []
        for dnasam, rnasam in zip(dnasamfiles, rnasamfiles):
            rnasortbam = run_samtools(rnasam)  # , rnasortsam
            dnasortbam = run_samtools(dnasam)  # , dnasortsam
            rnasortbamfiles.append(rnasortbam)
            dnasortbamfiles.append(dnasortbam)

    if 'rnadnacontacts' in steps:
        if steps[0] == 'rnadnacontacts':
            ## assuming user wants to start at this step
            ## assumes altreadfile1&2 are lists of files
            dnasortbamfiles = altreadfile1
            rnasortbamfiles = altreadfile2

        ### calculate number rna and dna contacts per window ###
        # full path to RNA sam file # full path for DNA sam file # MQ = minimum mapping quality # sn = Species name
        #RecordRNADNAContacts.main(rnasortsam, dnasortsam, MQ, sn)

        count = 0
        print('Starting to count rna and dna read positions and contacts')
        if PCRDupRemoval is True:
            print(f'Remove PCR duplicates because PCRDupRemoval is set to {PCRDupRemoval}')
        # this will be slow but will only keep one file of reads in memory at any one time
        # count = number of files, when 2 or more pairs of files it will concatenate results in record_rna_dna_contacts
        for count, (rnabam, dnabam) in enumerate(zip(rnasortbamfiles, dnasortbamfiles)):
            print(dnabam)
            drna, ddna = {}, {}
            drna = record_read_left_right_orientation(rnabam, MQ, drna)  # key: readname, value: scaff, left, right, orientation
            ddna = record_read_left_right_orientation(dnabam, MQ, ddna)
            rawcontactfile, rawcontactfilenames = record_rna_dna_contacts(dnabam, drna, ddna, sn, count, PCRDupRemoval)

        print('Finished recording 3\' rna and 5\'dna read positions and raw RNA-DNA contacts')

    if 'summarizecontacts' in steps:
        if steps[0] == 'summarizecontacts':
            ## assuming user wants to start at this step
            ## assumes altreadfile1&2 are lists of files
            rawcontactfilenames = altreadfile1  # this is raw contact file with read names (assumed not sorted)
        contactfilenames, contactfile, contactfilelog2, contactfilecircos = summarize_raw_rna_dna_contacts_and_names(rawcontactfilenames)

        # contactfile = summarize_raw_rna_dna_contacts(rawcontactfile, rawcontactfilenames, sn, PCRDupRemoval)


        #for dnasortbam, rnasortbam in zip(dnasortbamfiles, rnasortbamfiles):
        #    drna, ddna = record_read_positions(rnasortbam, dnasortbam, MQ, drna=drna, ddna=ddna)
        # possible that drna and ddna get to be large dictionaries


    #if 'normalize' in steps:
    #    ### Normalize ###
    #    ### calculate number of DPNII restriction sites per window size for the whole genome ###
    #    # GATC is RE site for DPNII
    #    # genomefile = '/media/sf_LinuxShare/Ciliates/Genomes/Seqs/ptetraurelia_mac_51.fa'
    #    if steps[0] == 'normalize':
    #        contactfile = altreadfile1
    #    else:
    #        path, f = os.path.split(dnasortbamfiles[0])
    #        contactfile = os.path.join(path, 'RNA.DNA.Contacts.pt.wthickness.txt')
    #    dpncountsperwindowfile, dRECounts, winnames = calculate_number_RE_per_window(genomefile, REsequence, windowsize)
    #    ### Divide support for each RNA DNA contact by number of DPNII sites surrounding each contact ###
    #    normalizedcontactsfile = normalize_rna_dna_contact_support(contactfile, dRECounts, MQ, contactcutoff, windowsize)


    if 'heatmap' in steps:
        if steps[0] == 'heatmap':
            #normalizedcontactsfile = altreadfile1
            contactfile = altreadfile1
        ## convert scaffold coordinates to whole genome coordinates where we add coordinates from all previous scaffolds
        # calculate lengths of all scaffolds
        dscafflengths, names = length_of_fasta_sequences(genomefile)
        # naturally sort scaffold names, assume naturally sorted names are increasing logically. 1,2,3,4.. not 1,10,11..
        natsortedscaffs = natsorted(names)
        # sum scaffold lengths for all scaffolds 'below' current scaffold
        dscaffgenomecoords, genomelength, firstlength = sum_fasta_lengths_by_sortlist(dscafflengths, natsortedscaffs)

        gencoordcontactsfile, x, y, intensity = convert_scaff_coords_to_continuous_genome_coords(contactfile,
                                                                                                   dscaffgenomecoords)
        print(gencoordcontactsfile)
        # sort x and y coordinates by intensity values (large are last)
        # forces most intense points to be plotted last and therefore always visible
        intensity, x, y = [list(tuple) for tuple in zip(*sorted(zip(intensity, x, y)))]  # sorts second and third list by values of first list
        maxintent = max(intensity)  # avoid division by zero
        intensities = [round(x / maxintent, 3) for x in intensity]
        #logintensities = [math.log2(x) / math.log2(maxintent) for x in intensity]

        # x is list of x int() values, y is list of corresponding y coordinates int() values
        print(intensities[:10])
        print(intensities[-100000:-100010])
        print(intensities[-10:])
        #seaborn_scatter_color_intensity(x[-1000:], y[-1000:], intensities[-1000:])
        #seaborn_scatter_color_intensity(x[-1000:], y[-1000:], logintensities[-1000:])
        ##seaborn_scatter_color_intensity(x[-2500:], y[-2500:], intensities[-2500:])
        ##seaborn_scatter_color_intensity(x[-2500:], y[-2500:], logintensities[-2500:])
        #seaborn_scatter_color_intensity(x[-5000:], y[-5000:], intensities[-5000:])
        # seaborn_scatter_color_intensity(x, y, intensities, outpath=gencoordcontactsfile + '.pdf')
        # seaborn_scatter_color_intensity(x, y, logintensities, outpath=gencoordcontactsfile + '.log2.pdf')

        ### zoom in on IESs ###
        # scaffold51_2_with_IES 355k to 375k
        # zoomscaffname = 'scaffold51_2'
        # zoomwinsize = 20
        # scaffcoordmin, scaffcoordmax = [354999], [374999]
        if zoomscaffname != '':
            dgff3, features = read_gff3(annot_file, features=['all'])
            for scaffmin, scaffmax in zip(scaffcoordmin, scaffcoordmax):
                print(f'Zooming in on {zoomscaffname} from {scaffmin} to {scaffmax}')
                mincoord = dscaffgenomecoords[zoomscaffname] + scaffmin
                maxcoord = dscaffgenomecoords[zoomscaffname] + scaffmax
                dgff3trim, trimfeatures, featurestarts, featurelengths = trim_features_by_scaffold_and_coordinates(dgff3, features, scaffold=zoomscaffname + '_with_IES', start=scaffmin+1, end=scaffmax+1)
                print(f'first five features: {trimfeatures[:5]}')
                print(f'number of features: {len(trimfeatures)}')
                # print(featurestarts[:20])
                # print(featurelengths[:20])
                ys = [0] * len(trimfeatures)
                starts = [(i-scaffmin)/zoomwinsize for i in featurestarts]
                lengths = [i/zoomwinsize for i in featurelengths]
                # print(len(starts))
                # print(len(lengths))
                # print(starts[:5])
                # print(lengths[:5])
                matrix, matrixlog = np_matrix_sum_intensity_if_xycoordinate_in_window(x, y, intensity, maxcoordinate=maxcoord,
                                                                                      windowsize=zoomwinsize,
                                                                                      mincoordinate=mincoord)
                seaborn_heatmap_with_barh(matrixlog, '%s.%s.%d.%d.%d.barh.log.pdf' % (gencoordcontactsfile, zoomscaffname, scaffmin,
                                                                            scaffmax, zoomwinsize), trimfeatures, ys, starts, lengths)
                # seaborn_heatmap_with_barh(matrix, '%s.%s.%d.%d.%d.barh.pdf' % (gencoordcontactsfile, zoomscaffname, scaffmin,
                #                                                         scaffmax, zoomwinsize), trimfeatures, ys, starts, lengths)
        ### ###

        # plot_kde_heatmap(x, y, nbins=100, outpath=gencoordcontactsfile + '.pdf')
        print(natsortedscaffs[0])
        winsize = 20000
        coord = firstlength
        matrix, matrixlog = np_matrix_sum_intensity_if_xycoordinate_in_window(x, y, intensity, maxcoordinate=coord, windowsize=winsize)
        seaborn_heatmap(matrixlog, '%s.%d.%d.log.pdf' % (gencoordcontactsfile, coord, winsize))
        matrix, matrixlog = np_matrix_sum_number_of_unique_xycoordinate_in_window(x, y, maxcoordinate=coord, windowsize=winsize)
        seaborn_heatmap(matrixlog, '%s.%d.%d.log.UniqueContacts.pdf' % (gencoordcontactsfile, coord, winsize))
        # winsize = 10000
        # coord = firstlength
        # matrix, matrixlog = np_matrix_sum_intensity_if_xycoordinate_in_window(x, y, intensity, maxcoordinate=coord, windowsize=winsize)
        # seaborn_heatmap(matrixlog, '%s.%d.%d.log.pdf' % (gencoordcontactsfile, coord, winsize))
        # seaborn_heatmap(matrix, '%s.%d.%d.pdf' % (gencoordcontactsfile, coord, winsize))
        winsize = 5000
        coord = firstlength
        matrix, matrixlog = np_matrix_sum_intensity_if_xycoordinate_in_window(x, y, intensity, maxcoordinate=coord, windowsize=winsize)
        seaborn_heatmap(matrixlog, '%s.%d.%d.log.pdf' % (gencoordcontactsfile, coord, winsize))
        matrix, matrixlog = np_matrix_sum_number_of_unique_xycoordinate_in_window(x, y, maxcoordinate=coord, windowsize=winsize)
        seaborn_heatmap(matrixlog, '%s.%d.%d.log.UniqueContacts.pdf' % (gencoordcontactsfile, coord, winsize))
        # winsize = 2000
        # coord = firstlength
        # matrix, matrixlog = np_matrix_sum_intensity_if_xycoordinate_in_window(x, y, intensity, maxcoordinate=coord, windowsize=winsize)
        # seaborn_heatmap(matrixlog, '%s.%d.%d.log.pdf' % (gencoordcontactsfile, coord, winsize))
        # seaborn_heatmap(matrix, '%s.%d.%d.pdf' % (gencoordcontactsfile, coord, winsize))

        print(natsortedscaffs[:3])
        winsize = 5000
        coord = 3000000
        matrix, matrixlog = np_matrix_sum_intensity_if_xycoordinate_in_window(x, y, intensity, maxcoordinate=coord, windowsize=winsize)
        seaborn_heatmap(matrixlog, '%s.%d.%d.log.pdf' % (gencoordcontactsfile, coord, winsize))
        matrix, matrixlog = np_matrix_sum_number_of_unique_xycoordinate_in_window(x, y, maxcoordinate=coord, windowsize=winsize)
        seaborn_heatmap(matrixlog, '%s.%d.%d.log.UniqueContacts.pdf' % (gencoordcontactsfile, coord, winsize))
        winsize = 50000
        matrix, matrixlog = np_matrix_sum_intensity_if_xycoordinate_in_window(x, y, intensity, maxcoordinate=genomelength, windowsize=winsize)
        seaborn_heatmap(matrixlog, '%s.%d.%d.log.pdf' % (gencoordcontactsfile, genomelength, winsize))
        matrix, matrixlog = np_matrix_sum_number_of_unique_xycoordinate_in_window(x, y, maxcoordinate=genomelength, windowsize=winsize)
        seaborn_heatmap(matrixlog, '%s.%d.%d.log.UniqueContacts.pdf' % (gencoordcontactsfile, coord, winsize))
        matrix, matrixlog, intensities = None, None, None   # clear memory

    if 'contacts in features' in steps:
        if steps[0] == 'contacts in features':
            contactfilenames = altreadfile1
        dgff3 = read_gff3_index_by_scaffold(annot_file, features)
        featurecontactfile = limit_contacts_by_feature_with_indices(dgff3, contactfilenames, dnabuffer, rnabuffer, rnamatch, maxgapoverlap)
        #dgff3 = read_gff3(annot_file, features)
        #normalize_contact_counts(dgff3, dcounts, contactfile, side='DNA', buffer=0)
        dgff3 = None

    if 'get sequences' in steps:
        if steps[0] == 'get sequences':
            featurecontactfile = altreadfile1
        seqnames = unique_sequence_names(featurecontactfile)
        path, file = os.path.split(forwardreadfile)
        peardirectory = os.path.join(path, 'pear')
        rnafastqfiles = list_of_files_from_dir(peardirectory, extension='.rna.fastq.gz')
        rnaseqfile, scnrnaseqfile = retrieve_rna_sequences_from_contacts_inside_features(seqnames, rnafastqfiles, dnabuffer, rnabuffer)

    if 'blast rna' in steps:
        # assumes input is .fa.gz so decompress first before blasting
        if steps[0] == 'blast rna':
            rnaseqfile = altreadfile1
        rnaseqfilefasta = run_gzip_decompress(rnaseqfile)
        run_blastn_match_db(rnaseqfilefasta, blastdatabase, 6, 50.00, True)

    if 'cumulative distribution' in steps:
        if steps[0] == 'cumulative distribution':
            contactfile = altreadfile1
        # cumulative distribution plot
        dgff3, featureids = read_gff3(annot_file, features)
        rnafeaturecontactfile = limit_contacts_by_feature(dgff3, contactfile)
        minrpkm = 0
        dcountnorm, normcountfile = count_contacts_by_feature_norm_length(dgff3, rnafeaturecontactfile, minrpkm)
        normcounts = [v for v in list(dcountnorm.values())]
        seaborn_kde_cumulative(normcounts, normcountfile + f'.cumulative.MinRPKM{minrpkm}.pdf')

    #if 'sequences' in steps:
    #    if steps[0] == 'sequences':
    #        # normalizedcontactsfile = altreadfile1
    #        contactfile = altreadfile1
    #    for feature in features:
    #        print(feature)
    #        dgff3, featureids = read_gff3(annot_file, feature)
    #        rnawindowstart = [dgff3[k][3] for k in featureids]
    #        rnawindowend = [dgff3[k][4] for k in featureids]
    #        # rnafeaturecontactfile = limit_contacts_by_feature(dgff3, normalizedcontactsfile)
    #        rnafeaturecontactfile = limit_contacts_by_feature(dgff3, contactfile)
    #        gencoordcontactsfile2, x2, y2, intensity2 = convert_scaff_coords_to_continuous_genome_coords(rnafeaturecontactfile, dscaffgenomecoords)
    #        winsize = 100000
    #        matrix2, matrixlog2 = np_matrix_sum_intensity_if_xycoordinate_in_window(x2, y2, intensity2, maxcoordinate=genomelength, windowsize=winsize)
    #        seaborn_heatmap(matrixlog2, '%s.%d.%d.%s.log.pdf' % (gencoordcontactsfile2, genomelength, winsize, feature[0]))
    #        seaborn_heatmap(matrix2, '%s.%d.%d.%s.pdf' % (gencoordcontactsfile2, genomelength, winsize, feature[0]))
    #
    #        # cumulative distribution plot
    #        minrpkm = 0
    #        dcountnorm, normcountfile = count_contacts_by_feature_norm_length(dgff3, rnafeaturecontactfile, minrpkm)
    #        normcounts = [v for v in list(dcountnorm.values())]
    #        seaborn_kde_cumulative(normcounts, normcountfile + f'.cumulative.MinRPKM{minrpkm}.pdf')

    print('##########\n    Fin\n##########')

    # calculate # of DPNII sites per window, calculate number of RNA-DNA contacts per window, divide # of contacts by # of DPNII sites
    # calculate # of DPNII sites per FEATURE (based on annotation), calculate number of RNA-DNA contacts per window, divide # of contacts (or RPKM/FPKM/TPM) by number of DPNII sites