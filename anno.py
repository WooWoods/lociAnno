import os
import re
import sys
import subprocess as sp

from collections import defaultdict

def read_config(configfile):
    config = {}
    with open(configfile) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            if re.match(r'^\s$', line):
                continue
            line = re.sub(r'[\r\n\s]+','',line)
            arr = line.split('=')
            config[arr[0]] = arr[1]
    #logging.info('configuration info loaded.')
    return config

def create_lib(vcf, outdir='.'):
    library = os.path.join(outdir, 'avinput')
    with open(vcf) as fh, open(library, 'w') as fo:
        for line in fh:
            if line.startswith('#CHROM'):
                fo.write(line)
                continue
            if line.startswith('#'):
                continue
            arr = line.strip().split('\t')
            chrom, start, _, ref, alters = arr[0:5]
            ref2 = ref
            alt_arr = alters.split(',')
            for alt in alt_arr:
                if alt == '*':
                    alt = '-'
                if len(alt) > 1 or len(ref) > 1:
                    if len(alt) > len(ref):
                        alt = re.sub(ref, '', alt)
                        ref2 = '-'
                    else:
                        ref2 = re.sub(alt, '', ref)
                        alt = '-'
                end = int(start) + len(ref2) -1
                samples_gt = '\t'.join(arr[9:])
                fo.write(f'{chrom}\t{start}\t{end}\t{ref2}\t{alt}\t{samples_gt}\n')
    return library

def annovar_filter(cfg, infile, dbtype):
    annovar = os.path.join(cfg['annovar'], 'annotate_variation.pl')
    build_ver = cfg['build_ver']
    humandb = cfg['humandb']
    
    sp.call(' '.join([annovar, '-filter', infile,'-buildver',
                              build_ver,'-dbtype', dbtype, humandb,
                              ]),shell=True)

def region_anno(cfg, infile):
    annovar = os.path.join(cfg['annovar'], 'annotate_variation.pl')
    build_ver = cfg['build_ver']
    humandb = cfg['humandb']
    sp.call(' '.join([annovar, '-geneanno','-buildver',
                              build_ver, infile, humandb
                              ]),shell=True)   

def sample_name(vcf):
    fname = os.path.basename(vcf)
    sname, *tmp = fname.split('.')
    return sname

def annotation(cfg, vcf):
    outdir = '.'
    filter_db = ['avsnp150', '1000g2015aug_all', 'cosmic70']
    library = create_lib(vcf, outdir)

    region_anno(cfg, library)
    for db in filter_db:
        annovar_filter(cfg, library, db)

    anno_info, items = anno_parser(library, filter_db)

    sname = sample_name(vcf)
    output = os.path.join(outdir, f'{sname}.tsv')
    dumper(library, anno_info, items, output, fmt='maf')

def anno_parser(library, dbs, build_ver='hg38'):
    anno_info = defaultdict(lambda : dict())
    anno_exon = library + '.exonic_variant_function'
    anno_nonexon = library + '.variant_function'
    items = ['type', 'function']

    with open(anno_exon) as fh:
        for line in fh:
            arr = line.split('\t')
            key = '_'.join(arr[3:5] + arr[6:7])
            anno_info[key]['function'] = arr[2]
            anno_info[key]['type'] = arr[1]
    with open(anno_nonexon) as fh:
        for line in fh:
            arr = line.split('\t')
            key = '_'.join(arr[2:4] + arr[5:6])
            if key in anno_info:
                continue
            anno_info[key]['function'] = arr[1]
            anno_info[key]['type'] = arr[0]

    for db in dbs:
        if db == '1000g2015aug_all':
            db = 'ALL.sites.2015_08'
        items.append(db)
        fannoed = f'{library}.{build_ver}_{db}_dropped'
        with open(fannoed) as fh:
            for line in fh:
                arr = line.split('\t')
                key = '_'.join(arr[2:4] + arr[5:6])
                anno_info[key][db] = arr[1]

    return anno_info, items

def dumper(library, anno_info, items, output, fmt='maf'):
    with open(library) as fh, open(output, 'w') as fo:
        header = fh.readline().strip().split()
        samples = header[9:]
        maf_title = ['sample', 'genotype', 'ref.cov', 'alt.cov']
        common_title = ['Chr', 'start', 'end', 'ref', 'alt']
        if fmt == 'maf':
            titles = common_title + items + maf_title
            fo.write('\t'.join(titles) + '\n')
        if fmt == 'tsv':
            titles = common_title + items + samples
            fo.write('\t'.join(titles) + '\n')

        for line in fh:
            arr = line.strip().split()
            tmp_arr = arr[0:5]
            key = '_'.join(arr[0:2] + arr[3:4])
            for item in items:
                annot = anno_info[key].get(item, ' ')
                tmp_arr.append(annot)

            point_dict = gt_parser(arr[5:], samples)
            if fmt == 'maf':
                for bc, var in point_dict.items():
                    fo.write('\t'.join(tmp_arr) + '\t' + bc + '\t' + var + '\n')
            if fmt == 'tsv':
                gts = [point_dict.get(bc, './.') for bc in samples]
                #gts = [extract_gt(point_dict.get(bc, './.')) for bc in samples]
                fo.write('\t'.join(tmp_arr + gts) + '\n')

def extract_gt(gt_item):
    arr = gt_item.split('\t')
    return arr[0]

def gt_parser(gt_info, bcs):
    gt = []
    ad = []
    point_dict = dict()
    for idx, item in enumerate(gt_info):
        arr = item.split(':')
        if arr[0] == './.':
            continue
        point_dict[bcs[idx]] = f'{arr[0]}\t{arr[1]}'
    return point_dict


if __name__ == '__main__':
    config, vcf = sys.argv[1:]
    cfg = read_config(config)
    annotation(cfg, vcf)

