"""A set of function definitions that can turn useful when dealing with biological data expressed with
genomic coordinates, such as bed/bedgraph formats"""
import pandas as pd

chromosome_lengths = pd.read_table('/home/mary/internship/hg19.chrom.sizes.txt', header=None)
chr_lengths_dict = dict()
for i in range(len(chromosome_lengths)):
    chr_lengths_dict[chromosome_lengths[0][i]] = chromosome_lengths[1][i]


def is_included(range1, range2):
    """range1 and range2 are tuples of int that contain the coordinates of start and end of the features to compare
    answer to question 'is range 1 included in range 2?' i.e. the smaller interval should be given first
    return values """
    if range1[0] > range2[0] and range1[1] < range2[1]:
        return True
    else:
        return False


def is_overlapping(range1, range2):
    # improvement to implement: need to create a function assessing "left overlap" or "right overlap"
    # so I can establish a policy of including a value in a bin only if the feature overlaps the bin in one side
    if is_included(range1, range2) or is_included(range2, range1):
        return True
    elif range1[0] in range(range2[0], range2[1]) or range1[1] in range(range2[0], range2[1]):
        return True
    else:
        return False


def overlap(range1, range2):
    if not is_overlapping(range1, range2):
        raise ValueError('Disjoint intervals, cannot calculate overlap')
    if is_included(range1, range2):
        return range1[1] - range1[0]
    else:
        if range1[0] in range(range2[0], range2[1]):
            return range2[1] - range1[0]
        if range1[1] in range(range2[0], range2[1]):
            return range1[1] - range2[0]


def extract_tot_gene_coords(intron_coords, exon_coords):
    """intron_coords and exon_coords are lists formatted like [chr, begin, end, ensemble_id, name, strand]
    return dictionary of lists of tuples
    the bad thing with this function is that it only works with the format Miguel gave me"""
    tot_coords = dict()
    for i in intron_coords:
        chrdict = tot_coords.get(i[0], dict())
        genelist = chrdict.get(i[3], list())
        genelist.append((int(i[1]), int(i[2])))
        chrdict[i[3]] = genelist
        tot_coords[i[0]] = chrdict
    for i in exon_coords:
        chrdict = tot_coords.get(i[0], dict())
        genelist = chrdict.get(i[3], list())
        genelist.append((int(i[1]), int(i[2])))
        chrdict[i[3]] = genelist
        tot_coords[i[0]] = chrdict
    return tot_coords


def extreme_coords(tot_coords):
    """tot_coords is already a dictionary
    return dictionary {chrN: {gene1: (start, end), gene2: (start, end), ...}}"""
    gene_coords = dict()
    for i in tot_coords.keys():
        chrdict = gene_coords.get(i, dict())
        for gene in tot_coords[i].keys():
            coord_list = tot_coords[i][gene]
            min_coord = 10000000000000
            max_coord = -1
            for j in coord_list:
                if j[0] < min_coord:
                    min_coord = j[0]
                if j[1] > max_coord:
                    max_coord = j[1]
            chrdict[gene] = (min_coord, max_coord)
        gene_coords[i] = chrdict
    return gene_coords


def extract_bedgraph(bedgraph_list, val_field=3):
    """bedgraph_list is the result of converting the file to list of lists: [chr, start, end, value]
    format of the dict: {chrN: {(start, end): feature, (start, end): feature)}}"""
    bedgraph_features = dict()
    for i in bedgraph_list:
        curr_chr = bedgraph_features.get(i[0], dict())
        curr_chr[(int(i[1]), int(i[2]))] = float(i[val_field])
        bedgraph_features[i[0]] = curr_chr
    return bedgraph_features


def make_bins(bin_width, chr_lengths=chr_lengths_dict):
    """defines start and end of the bins given the total length of each chr and a width
    return dict {chr1: [(1-100000), (100001-200000), ...]}"""
    bins_dict = dict()
    for chromosome in chr_lengths.keys():
        chr_length = chr_lengths[chromosome]
        start = 1
        end = bin_width
        chr_bins = list()
        while start < chr_length and end < chr_length:
            chr_bins.append((start, end))
            start = end + 1
            end = end + bin_width
        # append the "leftover" interval (probably not necessary to check for end >= chr_length)
        if end >= chr_length:
            chr_bins.append((start, chr_length))
        bins_dict[chromosome] = chr_bins
    return bins_dict


def avg_bins(bedgraph_features, bins_dict):
    """average out values in the bedgraph_features dictionary and make a new dictionary with the average value per bin
    structure of bedgraph_features: {chrN: {(start, end): feature, (start, end): feature)}}
    return dict: {chrN: {(start, end): value, (start, end): value, ...}}"""
    to_avg_dict = dict()
    for i in bedgraph_features.keys():  # chromosomes
        to_avg_dict[i] = dict()
        for j in bedgraph_features[i].keys():  # features
            # look for the right bin
            for bin in bins_dict[i]:
                to_avg = to_avg_dict[i].get(bin, list())
                # when the right bin is found, store the value and go to the next feature
                if is_included(j, bin):
                    # improvement to implement: if the feature is included or overlaps the end of the bin
                    to_avg.append(bedgraph_features[i][j])
                    to_avg_dict[i][bin] = to_avg
                    break
    avg_dict = dict()
    for chrom in to_avg_dict.keys():
        avg_dict[chrom] = dict()
        for bin in to_avg_dict[chrom].keys():
            avg_dict[chrom][bin] = sum(to_avg_dict[chrom][bin]) / len(to_avg_dict[chrom][bin])
    return avg_dict


def sum_bins(bedgraph_features, bins_dict):
    """sum all the values in the bedgraph_features dictionary and make a new dictionary with
    the total value per bin"""
    to_sum_dict = dict()
    for i in bedgraph_features.keys():  # chromosomes
        to_sum_dict[i] = dict()
        for j in bedgraph_features[i].keys():  # features
            # look for the right bin
            for bin in bins_dict[i]:
                to_sum = to_sum_dict[i].get(bin, list())
                # when the right bin is found, store the value and go to the next feature
                if is_included(j, bin):
                    # improvement to implement: if the feature is included or overlaps the end of the bin
                    to_sum.append(bedgraph_features[i][j])
                    to_sum_dict[i][bin] = to_sum
                    break
    sum_dict = dict()
    for chrom in to_sum_dict.keys():
        sum_dict[chrom] = dict()
        for bin in to_sum_dict[chrom].keys():
            sum_dict[chrom][bin] = sum(to_sum_dict[chrom][bin])
    return sum_dict


def ensg_to_value(gene_dict, avg_dict):
    """gene_dict is a dictionary resulting from extreme_coords:
    {chrN: {gene1: (start, end), gene2: (start, end), ...}}
    avg_dict is the output of avg_bins:
    {chrN: {(start, end): value, (start, end): value, ...}}
    return {chrN: {gene: value, gene: value, ...}}"""
    gene_binning = dict()
    for chromosome in gene_dict.keys():
        chrdict = gene_binning.get(chromosome, dict())
        for gene in gene_dict[chromosome].keys():
            for chrbin in avg_dict[chromosome].keys():
                if is_overlapping(gene_dict[chromosome][gene], chrbin):
                    chrdict[gene] = avg_dict[chromosome][chrbin]
        gene_binning[chromosome] = chrdict
    return gene_binning


def ensg_as_bins(gene_dict, bed_features):
    """gene_dict is the output of extreme_coords(); bed_features is the output of extract_bedgraph()
    map features to each gene and then make an average of all the values falling in a certain gene
    formats:
    {chrN: {gene1: (start, end), gene2: (start, end), ...}}
    {chrN: {(start, end): feature, (start, end): feature)}}

    output:
    {chrN: {gene1: feature, gene2: feature, ...}}
    I am so sorry for the complexity of this but I don't have the strength to implement it smartly

    (breaking news: I wrote this for nothing because the dnase data were not dense enough. but who knows,
    it might come in handy some other time)"""
    gene_binning = dict()

    # make list of values that pertain to each gene
    for chromosome in gene_dict.keys():
        chrdict = gene_binning.get(chromosome, dict())
        for gene in gene_dict[chromosome].keys():
            genelist = chrdict.get(gene, list())
            for feature in bed_features[chromosome].keys():
                if is_overlapping(feature, gene_dict[chromosome][gene]):
                    genelist.append(bed_features[chromosome][feature])
            chrdict[gene] = genelist
        gene_binning[chromosome] = chrdict

    # make average of all listed values for each gene
    zerogenes = 0
    for chromosome in gene_binning.keys():
        for gene in gene_binning[chromosome].keys():
            if len(gene_binning[chromosome][gene]) == 0:
                gene_binning[chromosome][gene] = 0
                zerogenes += 1
            else:
                gene_binning[chromosome][gene] = sum(gene_binning[chromosome][gene])/len(gene_binning[chromosome][gene])
    print('genes with no associated dnase value: %d' % zerogenes)
    return gene_binning


if __name__ == '__main__':
    # check for interval overlap functions
    # print(is_included((2, 10), (0, 150)))
    # print(is_overlapping((1,10),(5,20)))
    # print(is_overlapping((5,10),(2,8)))
    # print(is_overlapping((1,5),(7,9)))

    # check bin building function
    # print(make_bins(100000))

    # check bedgraph extraction function
    f = open('/home/mary/internship/covariates/replicationdomain/RT_CyT49_ESC_Int73235540_hg19.bedgraph')
    bedlines = f.readlines()[11:]
    for i in range(len(bedlines)):
        bedlines[i] = bedlines[i].strip().split()
    # print(extract_bedgraph(bedlines))

    # check avging function
    avg_rt = avg_bins(extract_bedgraph(bedlines), make_bins(100000))

    # check coord extraction function
    f1 = open('/home/mary/internship/mut_freq/clean_intron_exon_coords/introns_CCDS.bed')
    f2 = open('/home/mary/internship/mut_freq/clean_intron_exon_coords/exons_CCDS.bed')
    lines1 = f1.readlines()
    lines2 = f2.readlines()
    for i in range(len(lines1)):
        lines1[i] = lines1[i].strip().split()
    for i in range(len(lines2)):
        lines2[i] = lines2[i].strip().split()

    gene_coords = extreme_coords(extract_tot_gene_coords(lines1, lines2))

    # check gene binning function
    print(ensg_to_value(gene_coords, avg_rt))
