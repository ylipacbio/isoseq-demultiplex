#!/usr/bin/env python

import sys
import os
import os.path as op
from argparse import ArgumentParser
from .utils import *

MOVIES = ['m54006_170729_232022', 'm54026_170727_103805', 'm54200_170721_210832', 'm54200_170722_173443']
IDXS = range(0, 4)
MOVIE2IDX = dict(zip(MOVIES, IDXS))
IDX2MOVIE = dict(zip(IDXS, MOVIES))

def parse_cluster_report_line(s):
    """
    i0_ICE_sample6343b5|c11,m54200_170722_173443/38404462/30_644_CCS,FL
    return (cid, zmw, is_FL)
    """
    cid, read, is_FL = s.strip().split(',')[0:3]
    return (cid, readname2moviezmw(read), True if is_FL == 'FL' else False)


class ReportObj(object):
    def __init__(self, cid, zmw, is_flnc):
        self.cid, self.zmw, self.is_flnc = cid, zmw, is_flnc

    def __hash__(self):
        return (self.cid, self.zmw, self.is_flnc)

    def __str__(self):
        return '(%s, %s, %s)' % (self.cid, self.zmw, self.is_flnc)

    def __repr__(self):
        return str(self)

    def __lt__(self, other):
        return self.__cmp__(other) < 0

    def __cmp__(self, other):
        return cmp((self.cid, self.zmw, self.is_flnc), (other.cid, other.zmw, other.is_flnc))

def yield_cluster_report(fp):
    """
    ...doctest: >>> lines = ['cid1,movie/100,FL', 'cid1,movie/101,FL', 'cid1,movie/102,NonFL', 'cid2,movie/103,FL', 'cid3,movie/103,FL', 'cid4,movie/104,NonFL'] >>> [r for r in yield_cluster_report(lines)]
        [[(cid1, movie/100, True), (cid1, movie/101, True), (cid1, movie/102, False)], [(cid2, movie/103, True)], [(cid3, movie/103, True)], [(cid4, movie/104, False)]]
    """
    prev = None
    ret = []
    for line in fp:
        if not line.startswith('#') and not line.startswith('cluster_id'):
            cid, zmw, is_flnc = parse_cluster_report_line(line)
            if prev is None:
                prev = ReportObj(cid, zmw, is_flnc)
            else:
                if cid == prev.cid:
                    ret.append(prev)
                    prev = ReportObj(cid, zmw, is_flnc)
                else:
                    ret.append(prev)
                    yield ret
                    ret = []
                    prev = ReportObj(cid, zmw, is_flnc)
    ret.append(prev)
    yield ret


def list_to_str(items):
    """
    ...doctest:
        >>> list_to_str([1, 2, 34, 'a'])
        '[1,2,34,a]'
    """
    return '[' + ','.join([str(x) for x in items]) + ']'

def get_most_common_or_none(items, min_fraction):
    """
    ...doctest:
        >>> get_most_common_or_none([1,1,1,1,1,2], 0.6)
        1
        >>> get_most_common_or_none([], 0.6) is None
        True
    """
    if items:
        return get_most_common_item(items, min_fraction)
    else:
        return None

def get_consensus_primer_from_flnc_nfl_primers(flnc_primers, nfl_primers, min_fraction, cid):
    """
    ...doctest:
        >>> get_consensus_primer_from_flnc_nfl_primers([0, 0, 0, 0], [], 0.6, 'cid')
        0
        >>> get_consensus_primer_from_flnc_nfl_primers([0, 0, 1, 1], [0, 0], 0.6, 'cid')
        0
        >>> get_consensus_primer_from_flnc_nfl_primers([0, 0, 1, 1], [0, 1], 0.6, 'cid') is None
        True
    """
    flnc_primers_wo_none = [p for p in flnc_primers if p is not None]
    nfl_primers_wo_none = [p for p in nfl_primers if p is not None]

    from_flnc = get_most_common_or_none(flnc_primers_wo_none, min_fraction=min_fraction)
    from_nfl = get_most_common_or_none(nfl_primers_wo_none, min_fraction=min_fraction)
    from_both = get_most_common_or_none(flnc_primers_wo_none + nfl_primers_wo_none, min_fraction=min_fraction)
    if from_flnc is not None and from_both is not None and from_flnc != from_both:
        print "Warning %s consensus primer from flnc and (flnc+nfl) are different" % (cid)
    return from_flnc if from_flnc is not None else from_both


class ClusterDict(object):
    __sep__ = '\t'
    def __init__(self, cluster_reports, flnc_z2p, nfl_z2p, min_fraction=0.6):
        assert len(cluster_reports) > 0
        self.cid = cluster_reports[0].cid
        self.flnc_zmws = [cluster_report.zmw for cluster_report in cluster_reports if cluster_report.is_flnc is True]
        self.nfl_zmws = [cluster_report.zmw for cluster_report in cluster_reports if cluster_report.is_flnc is False]
        self.flnc_primers = [flnc_z2p[flnc_zmw] for flnc_zmw in self.flnc_zmws]
        self.nfl_primers = [nfl_z2p[nfl_zmw] for nfl_zmw in self.nfl_zmws]
        self.min_fraction = min_fraction

    @property
    def consensus_primer(self):
        return get_consensus_primer_from_flnc_nfl_primers(self.flnc_primers, self.nfl_primers, self.min_fraction, self.cid)

    @property
    def header(self):
        return self.__sep__.join(['cluster_id', 'flnc_zmws', 'nfl_zmws', 'flnc_primers', 'nfl_primers', 'consensus_primer'])

    def to_str(self):
        return self.__sep__.join([str(x) for x in [self.cid, list_to_str(self.flnc_zmws), list_to_str(self.nfl_zmws), list_to_str(self.flnc_primers), list_to_str(self.nfl_primers), self.consensus_primer]])

    @property
    def flnc_z2c(self):
        return {flnc_zmw: self.cid for flnc_zmw in self.flnc_zmws}

    @property
    def nfl_z2c(self):
        return {nfl_zmw: self.cid for nfl_zmw in self.nfl_zmws}

    @classmethod
    def fromString(cls, s, min_fraction=0.6):
        """
        >>> s = 'cid\\t[movie1/1,movie1/2]\\t[movie2/3,movie3/4]\\t[0,0]\\t[1,None]\\t0'
        >>> o = ClusterDict.fromString(s)
        >>> o.cid
        'cid'
        >>> o.flnc_primers
        [0, 0]
        >>> o.nfl_primers
        [1, None]
        >>> o.flnc_z2c
        {'movie1/1': 'cid', 'movie1/2': 'cid'}
        >>> o.nfl_z2c
        {'movie3/4': 'cid', 'movie2/3': 'cid'}
        >>> o.flnc_zmws
        ['movie1/1', 'movie1/2']
        >>> o.nfl_zmws
        ['movie2/3', 'movie3/4']
        """
        def int_or_none(x):
            return None if x == 'None' else int(x)

        def str_or_none(x):
            return None if x == 'None' else str(x)

        fs = s.strip().split(cls.__sep__)
        cid = fs[0]
        flnc_zmws = [str_or_none(x) for x in fs[1][1:-1].split(',')]
        nfl_zmws = [str_or_none(x) for x in fs[2][1:-1].split(',')]
        flnc_primers = [int_or_none(x) for x in fs[3][1:-1].split(',')]
        nfl_primers = [int_or_none(x) for x in fs[4][1:-1].split(',')]
        reports = [ReportObj(cid, zmw, True) for zmw in flnc_zmws] + [ReportObj(cid, zmw, False) for zmw in nfl_zmws]
        flnc_z2p = dict(zip(flnc_zmws, flnc_primers))
        nfl_z2p = dict(zip(nfl_zmws, nfl_primers))
        return ClusterDict(reports, flnc_z2p, nfl_z2p, min_fraction)


def simplified_zmw(zmw, movie2idx):
    """
    ...doctest:
        >>> simplified_zmw('movie1/0', {'movie1':1, 'movie2':2})
        (1, 0)
    """
    movie, zmw_id = zmw.split('/')
    movie_id = movie2idx[movie]
    return (int(movie_id), int(zmw_id))


def simplified_dict(d, movie2idx):
    """
    ...doctest:
        >>> movie2idx = {'movie1': 1, 'movie2': 2}
        >>> d = {'movie1/101': 'cid1', 'movie2/102': 'cid2'}
        >>> simplified_dict(d, movie2idx)
        {(1, 101): 'cid1', (2, 102): 'cid2'}
    """
    ret = {}
    for zmw, cid in d.iteritems():
        ret[simplified_zmw(zmw, movie2idx)] = cid
    return ret

def write_z2c(z2c_dict, z2c_fn, idx2movie):
    z2c_writer = open(z2c_fn, 'w')
    for (movie_id, zmw_id), cids in z2c_dict.iteritems():
        z2c_writer.write('\t'.join(['%s/%s' % (idx2movie[movie_id], zmw_id), list_to_str(cids)]) + '\n')
    z2c_writer.close()


def write_ophan_zmws(zmws, simplified_zmws_in_cluster, movie2idx, o_fn):
    with open(o_fn, 'w') as writer:
        for zmw in zmws:
            if not simplified_zmw(zmw, movie2idx) in simplified_zmws_in_cluster:
                writer.write('%s\n' % zmw)


def update_simplified_z2c(simplified_z2c, z2c, movie2idx):
    """
    ...doctest:
        >>> s_z2c = {(1,1): [], (2,2):['cid2']}
        >>> z2c = {'movie1/1': 'cid1', 'movie2/2': 'cid3'}
        >>> update_simplified_z2c(s_z2c, z2c, {'movie1':1, 'movie2':2})
        {(1, 1): ['cid1'], (2, 2): ['cid2', 'cid3']}
        >>> s_z2c
        {(1, 1): ['cid1'], (2, 2): ['cid2', 'cid3']}
    """
    for zmw, cid in z2c.iteritems():
        s_zmw = simplified_zmw(zmw, movie2idx)
        simplified_z2c[s_zmw].append(cid)
    return simplified_z2c


def cluster_to_consensus_primer(cluster_report_fn, flnc_z2p, nfl_z2p, out_dir):
    flnc_z2c_fn, nfl_z2c_fn = op.join(out_dir, 'flnc.z2c.csv'), op.join(out_dir, 'nfl.z2c.csv')
    o_cluster_dict_csv_fn = op.join(out_dir, 'cluster_dict.csv')
    cluster_report_reader = open(cluster_report_fn, 'r')
    cluster_dict_writer = open(o_cluster_dict_csv_fn, 'w')

    s_flnc_z2c = defaultdict(lambda: [])
    s_nfl_z2c = defaultdict(lambda: [])
    print 'Reading %s' %  (cluster_report_fn)

    for cluster_records in yield_cluster_report(cluster_report_reader):
        cluster_dict = ClusterDict(cluster_records, flnc_z2p, nfl_z2p)
        cluster_dict_writer.write(cluster_dict.to_str() + '\n')
        update_simplified_z2c(s_flnc_z2c, cluster_dict.flnc_z2c, MOVIE2IDX)
        update_simplified_z2c(s_nfl_z2c, cluster_dict.nfl_z2c, MOVIE2IDX)

    cluster_report_reader.close()
    cluster_dict_writer.close()

    print 'Writing z2c %s' %  (flnc_z2c_fn)
    write_z2c(s_flnc_z2c, flnc_z2c_fn, IDX2MOVIE)
    print 'Writing nfl z2c %s' %  (nfl_z2c_fn)
    write_z2c(s_nfl_z2c, nfl_z2c_fn, IDX2MOVIE)


#super slow, ignore
#def write_other(flnc_z2p, flnc_z2c, nfl_z2p, nfl_z2c):
#    print 'Writing orphan flnc zmws %s' %  ('flnc.orphan.csv')
#    write_ophan_zmws(zmws=flnc_z2p.keys(), simplified_zmws_in_cluster=flnc_z2c.keys(), movie2idx=MOVIE2IDX, o_fn='flnc.orphan.csv')
#    print 'Writing orphan nfl zmws %s' %  ('nfl.orphan.csv')
#    write_ophan_zmws(zmws=nfl_z2p.keys(), simplified_zmws_in_cluster=nfl_z2c.keys(), movie2idx=MOVIE2IDX, o_fn='nfl.orphan.csv')


CLUSTER_REPORT_FN = '/pbi/dept/secondary/siv/smrtlink/smrtlink-alpha/jobs-root/020/020643/tasks/pbtranscript.tasks.separate_flnc-0/combined/all.cluster_report.csv'
FLNC_FA_FN, NFL_FA_FN = '/pbi/dept/secondary/siv/yli/isoseq/lima/data/flnc.fasta', '/pbi/dept/secondary/siv/yli/isoseq/lima/data/nfl.fasta'

def get_parser():
    """return arg parser"""
    desc = """Find consensus primers of clusters."""
    parser = ArgumentParser(description=desc)
    parser.add_argument("flnc_fa_fn", help="Input FLNC FASTA, e.g., %s" % FLNC_FA_FN)
    parser.add_argument("nfl_fa_fn", help="Input FLNC FASTA, e.g., %s" % NFL_FA_FN)
    parser.add_argument("cluster_report_fn", help="Input consensus_report_fn, e.g., %s" % CLUSTER_REPORT_FN)
    parser.add_argument("out_dir", help="Output directory")
    return parser

def run(args):
    lazy = False
    if not op.exists(args.out_dir):
        raise ValueError("Must create output directory %s" % args.out_dir)
    flnc_z2p, nfl_z2p = get_all_z2p(args.flnc_fa_fn, args.nfl_fa_fn, o_dir=args.out_dir, lazy=lazy)
    cluster_to_consensus_primer(args.cluster_report_fn, flnc_z2p, nfl_z2p, args.out_dir)


def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
