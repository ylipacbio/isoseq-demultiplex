#!/usr/bin/env python
import sys
import os
import os.path as op
from argparse import ArgumentParser
from .utils import *
from .cluster_to_consensus_primer import ClusterDict, get_most_common_or_none


def parse_z2c_line(line):
    fs = line.strip().split('\t')
    zmw = fs[0]
    cids = fs[1][1:-1].split(',')
    return (zmw, cids)

def yield_zmw_cids_from_z2c_fn(z2c_fn):
    with open(z2c_fn, 'r') as reader:
        for line in reader:
            zmw, cids = parse_z2c_line(line=line)
            yield (zmw, cids)


def get_z2cp(z2c_iterator, c2cp, min_fraction):
    """
    ...doctest:
        >>> z2c = {'movie/0': ['cid1', 'cid2'], 'movie/1': ['cid2', 'cid3'], 'movie/3': ['cid3', 'cid4'], 'movie/4': ['cid4']}
        >>> c2cp = {'cid1': 0, 'cid2': 0, 'cid3': 1, 'cid4':  None}
        >>> dict(get_z2cp(z2c.iteritems(), c2cp, 0.6))
        {'movie/4': None, 'movie/3': None, 'movie/0': 0, 'movie/1': 1}
    """
    z2cp = defaultdict(lambda: None)
    for zmw, cids in z2c_iterator:
        cprimers = [c2cp[cid] for cid in cids]
        cprimers_wo_none = [cp for cp in cprimers if cp is not None]
        ccprimer = get_most_common_or_none(cprimers_wo_none, min_fraction=min_fraction)
        if len(cids) > 0:
            z2cp[zmw] = c2cp[cid]
        else:
            z2cp[zmw] = None
    return z2cp

def get_c2cp_from_cluster_dict_reader(cluster_dict_reader, min_fraction):
    """
    ...doctest:
        >>> s = 'cid\\t[movie1/1,movie1/2]\\t[movie2/3,movie3/4]\\t[0,0]\\t[1,None]\\t0'
        >>> dict(get_c2cp_from_cluster_dict_reader([s], 0.6))
        {'cid': 0}
    """
    c2cp = defaultdict(lambda: None)
    for s in cluster_dict_reader:
        obj = ClusterDict.fromString(s, min_fraction)
        c2cp[obj.cid] = obj.consensus_primer
    return c2cp


def get_c2cp_from_cluster_dict_fn(cluster_dict_fn, min_fraction):
    return get_c2cp_from_cluster_dict_reader(open(cluster_dict_fn, 'r'), min_fraction)


def get_parser():
    """return arg parser"""
    desc = """Link zmw to consensus primers of its associated clusters."""
    parser = ArgumentParser(description=desc)
    #parser.add_argument("flnc_fa_fn", help="Input FLNC FASTA, e.g., %s" % FLNC_FA_FN)
    parser.add_argument("flnc_z2c_fn", help="Input FLNC zmws to cluster csv, i.e., flnc_z2c.csv generated by cluster-to-consensus-primer")
    parser.add_argument("nfl_z2c_fn", help="Input NFL zmws to cluster csv, i.e., nfl_z2c.csv generated by cluster-to-consensus-primer")
    parser.add_argument("cluster_dict_fn", help="Input cluster to consensus primer csv, i.e., cluster_dict.csv generated by cluster-to-consensus-primer")
    parser.add_argument("out_dir", help="Output directory")
    parser.add_argument("--min_fraction", help="Minimum fraction of consensus primer among all known primers.", default=0.6, type=float)
    return parser

def run(args):
    print 'get_c2cp_from_cluster_dict_fn'
    c2cp = get_c2cp_from_cluster_dict_fn(args.cluster_dict_fn, args.min_fraction)

    print 'get_z2cp(flnc_z2c, c2cp...)'
    flnc_z2c_iterator = yield_zmw_cids_from_z2c_fn(args.flnc_z2c_fn)
    flnc_z2cp = get_z2cp(flnc_z2c_iterator, c2cp, args.min_fraction)
    print 'write_dict(flnc_z2cp)'
    write_dict(dict(flnc_z2cp), o_prefix=op.join(args.out_dir, 'flnc_z2cp'), headers=['flnc_zmw', 'consensus_primer'])

    print 'get_z2cp(nfl_z2c, c2cp...)'
    nfl_z2c_iterator = yield_zmw_cids_from_z2c_fn(args.nfl_z2c_fn)
    nfl_z2cp = get_z2cp(nfl_z2c_iterator, c2cp, args.min_fraction)
    print 'write_dict(nfl_z2cp)'
    write_dict(dict(nfl_z2cp), o_prefix=op.join(args.out_dir, 'nfl_z2cp'), headers=['nfl_zmw', 'consensus_primer'])


def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
