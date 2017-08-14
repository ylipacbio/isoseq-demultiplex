from pbcore.io import FastaReader
import os
import os.path as op
import json
import pickle
from collections import defaultdict, Counter


class Obj(object):
    def __init__(self, primer):
        self.primer = primer
        self.val = None


def readname2moviezmw(name):
    """
    ...doctest:
        >>> readname2moviezmw('movie/100/0_1000 blala')
        'movie/100'
        >>> readname2moviezmw('movie/100')
        'movie/100'
    """
    return '/'.join(name.split('/')[0:2])

def readname2primer(name):
    """
    ...doctest:
        >>> readname2primer('myname polyA=50;primer=1')
        1
        >>> readname2primer('myname polyA=50;primer=NA') is None
        True
    """
    primer = [x for x in name.split(';') if 'primer=' in x][0][len('primer='):]
    if primer is None or primer == 'NA' or primer == 'None':
        return None
    else:
        return int(primer)

def zmw_to_primer_from_fasta(fasta_fn):
    """Input isoseq_flnc or isoseq_nfl or isoseq_draft fasta file, return
    dict {'movie/zmw': primer index}
    """
    zmw2primer = {}
    for r in FastaReader(fasta_fn):
        zmw2primer[readname2moviezmw(r.name)] = readname2primer(r.name)
    return zmw2primer


class Pcid(object):

    def __init__(self, c_prefix, c_id):
        self.c_prefix = c_prefix
        self.c_id = c_id

    def __hash__(self):
        return (self.c_prefix, self.c_id)

    def __repr__(self):
        return '(%r, %r)' % (self.c_prefix, self.c_id)

    def __cmp__(self, other):
        if self.c_prefix < other.c_prefix:
            return -1
        elif self.c_prefix > other.c_prefix:
            return 1
        else:
            if self.c_id < other.c_id:
                return -1
            elif self.c_id > other.c_id:
                return 1
            else:
                return 0

    def __eq__(self, other):
        return self.__cmp__(other) == 0

    def __lt__(self, other):
        return self.__cmp__(other) < 0

    def __gt__(self, other):
        return self.__cmp__(other) > 0


def zmw_to_cids_from_partial_pickle_fns(c_prefix_to_pickle_fn_dict):
    ret = defaultdict(lambda: [])
    for c_prefix, pickle_fn in c_prefix_to_pickle_fn_dict.iteritems():
         d = zmw_to_cids_from_partial_pickle_fn(c_prefix, pickle_fn)
         for zmw, cids in d.iteritems():
             ret[zmw].extend(cids)
    return ret


def zmw_to_cids_from_partial_pickle_fn(c_prefix, pickle_fn):
    pickle_d = pickle.load(open(pickle_fn, 'r'))['partial_uc']
    return zmw_to_cids_from_partial_pickle_d(c_prefix, pickle_d)


def zmw_to_cids_from_partial_pickle_d(c_prefix, pickle_d):
    """
    pickle_d: {
    ...doctest:
        >>> d = {1855: ['m54006_170729_232022/10486538/0_3405_CCS', 'm54200_170721_210832/45154853/0_5425_CCS'], 1900:['m54006_170729_232022/10486538/0_3405_CCS']}
        >>> dict(zmw_to_cids_from_partial_pickle_d('my_c_prefix', d))
        {'m54200_170721_210832/45154853': [('my_c_prefix', 1855)], 'm54006_170729_232022/10486538': [('my_c_prefix', 1900), ('my_c_prefix', 1855)]}
    """
    zmw2cids = defaultdict(lambda: [])
    for c_id, reads_in_c in pickle_d.iteritems():
        for read in reads_in_c:
            zmw = readname2moviezmw(read)
            zmw2cids[zmw].append((c_prefix, c_id))
    return zmw2cids


def zmw_to_cid_from_flnc_pickle_fns(c_prefix_to_pickle_fn_dict):
    """
    c_prefix_to_pickle_fn_dict: {c_prefix: flnc_pickle_fn}
    """
    ret = {}
    for c_prefix, pickle_fn in c_prefix_to_pickle_fn_dict.iteritems():
        d = dict(zmw_to_cid_from_flnc_pickle_fn(c_prefix, pickle_fn))
        ret.update(d)
    return ret


def zmw_to_cid_from_flnc_pickle_fn(c_prefix, pickle_fn):
    pickle_d = pickle.load(open(pickle_fn, 'r'))['d'] # dict{readname: {cid: weight}}
    return zmw_to_cid_from_flnc_pickle_d(c_prefix, pickle_d)

def zmw_to_cid_from_flnc_pickle_d(c_prefix, pickle_d):
    """
    ...doctest:
        >>> d = {'m54200_170722_173443/71500425/13572_72_CCS': {38: -135.67953402228963}}
        >>> dict(zmw_to_cid_from_flnc_pickle_d('my_c_prefix', d))
        {'m54200_170722_173443/71500425': ('my_c_prefix', 38)}
    """
    zmw2cid = defaultdict(lambda: -1)
    for read, cid2w in pickle_d.iteritems():
        zmw = readname2moviezmw(read)
        for c_id, weight in cid2w.iteritems():
            if zmw2cid[zmw] != -1:
                raise ValueError("FLNC zmw %s maps to multiple cids %s, (%s, %s)" % (zmw, zmw2cid[zmw], c_prefix, c_id))
            zmw2cid[zmw] = (c_prefix, c_id)
    return zmw2cid


def cid_to_zmws_from_flnc_pickle_fns(c_prefix_to_pickle_fn_dict):
    ret = defaultdict(lambda: [])
    for c_prefix, pickle_fn in c_prefix_to_pickle_fn_dict.iteritems():
        d = cid_to_zmws_from_flnc_pickle_fn(c_prefix, pickle_fn)
        for cid, zmws in d.iteritems():
            ret[cid].extend(zmws)
    return ret


def cid_to_zmws_from_flnc_pickle_fn(c_prefix, flnc_pickle_fn):
    """return {(c_prefix, cid): [zmws]}
    """
    pickle_d = pickle.load(open(flnc_pickle_fn, 'r'))['d'] # dict{readname: {cid: weight}}
    return cid_to_zmws_from_flnc_pickle_d(c_prefix, pickle_d)


def cid_to_zmws_from_flnc_pickle_d(c_prefix, pickle_d):
    """
    ...doctest:
        >>> d = {'m54200_170722_173443/71500425/13572_72_CCS': {38: -135.67953402228963}}
        >>> dict(cid_to_zmws_from_flnc_pickle_d('my_c_prefix', d))
        {('my_c_prefix', 38): ['m54200_170722_173443/71500425']}

    """
    cid2zmws = defaultdict(lambda: [])
    for read, cid2w in pickle_d.iteritems():
        for c_id, w in cid2w.iteritems():
            cid2zmws[(c_prefix,c_id)].append(readname2moviezmw(read))
        break
    return cid2zmws


def cid_to_zmws_from_partial_pickle_fn(c_prefix, partial_pickle_fn):
    pickle_d = pickle.load(open(partial_pickle_fn, 'r'))['partial_uc']
    return cid_to_zmws_from_partial_pickle_d(c_prefix, pickle_d)


def cid_to_zmws_from_partial_pickle_d(c_prefix, pickle_d):
    """
    ...doctest:
        >>> d = {1855: ['m54006_170729_232022/10486538/0_3405_CCS', 'm54200_170721_210832/45154853/0_5425_CCS'], 1900:['m54006_170729_232022/10486538/0_3405_CCS']}
        >>> dict(cid_to_zmws_from_partial_pickle_d('my_c_prefix', d))
        {('my_c_prefix', 1855): ['m54006_170729_232022/10486538', 'm54200_170721_210832/45154853'], ('my_c_prefix', 1900): ['m54006_170729_232022/10486538']}
    """
    cid2zmws = defaultdict(lambda: [])
    for c_id, reads_in_c in pickle_d.iteritems():
        for read in reads_in_c:
            cid2zmws[(c_prefix, c_id)].append(readname2moviezmw(read))
    return cid2zmws


#def zmw_to_cids_from_partial_pickle_fn(c_prefix, partial_pickle_fn):
#    pickle_d = pickle.load(open(pickle_fn, 'r'))['partial_uc']
#    return zmw_to_cids_from_partial_pickle_d(c_prefix, pickle_d)
#
#
#def zmw_to_cids_from_partial_pickle_d(c_prefix, pickle_d):
#    """
#    ...doctest:
#        >>> d = {1855: ['m54006_170729_232022/10486538/0_3405_CCS', 'm54200_170721_210832/45154853/0_5425_CCS'], 1900:['m54006_170729_232022/10486538/0_3405_CCS']}
#        >>> dict(zmw_to_cids_from_partial_pickle_d('my_c_prefix', d))
#        {'m54200_170721_210832/45154853': [('my_c_prefix', 1855)], 'm54006_170729_232022/10486538': [('my_c_prefix', 1900), ('my_c_prefix', 1855)]}
#    """
#    zmw2cids = defaultdict(lambda: [])
#    for c_id, reads_in_c in pickle_d.iteritems():
#        for read in reads_in_c:
#            zmw2cids[readname2moviezmw(read)].extend([(c_prefix, c_id)])
#    return  zmw2cids


def cid_to_primer_count(cid2zmws, zmw2primer, default_primers):
    """
    cid2zmws: dict{cid: zmws}
    """
    cid2primer_count = defaultdict(lambda: {primer:0 for primer in default_primers})
    for c_id, zmws in cid2zmws.iteritems():
        for zmw in zmws:
            zmw_primer = zmw2primer[zmw]
            cid2primer_count[c_id][zmw_primer] += 1
    return cid2primer_count

def get_most_common_item(items, min_fraction=0.9):
    """Return the most common item, which >= min_fraction among all items
    ...doctest:
        >>> get_most_common_item([1] * 10 + [2], 0.9)
        1
        >>> get_most_common_item([1] * 10 + [2]*10, 0.9) is None
        True
    """
    counter = Counter(items)
    most_item, most_item_number = counter.most_common()[0]
    if most_item_number >= len(items) * min_fraction:
        return most_item
    return None


def cid_to_consensus_primer(flnc_c2z, z2p, get_consensus_func):
    """
    flnc_c2z: dict{(c_prefix, cid): zmw}
    z2p: dict{zmw: int(primer)}, where primer -1 meaning no primer is detected by isoseq classify
    return {(c_prefix, cid): consensus_primer_of_cid}
    ...doctest:
        >>> def f(items): return items[0]
        >>> c2z = {('p', 100): range(0, 10), ('m', 101): range(5, 15), ('z', 102): range(15, 20)}
        >>> z2p = dict(zip(range(0, 20), [1]*14 + [-1] + [14,15,14,15,14]))
        >>> dict(cid_to_consensus_primer(c2z, z2p, f))
        {('m', 101): 1, ('z', 102): 14, ('p', 100): 1}
        >>> dict(cid_to_consensus_primer(c2z, z2p, get_most_common_item))
        {('m', 101): 1, ('z', 102): None, ('p', 100): 1}
    """
    c2cp = defaultdict(lambda: None) # default is not processed
    for cprefix_cid_tuple, zmws in flnc_c2z.iteritems():
        zmw_primers = [z2p[zmw] for zmw in zmws]
        c2cp[cprefix_cid_tuple] = get_consensus_func(zmw_primers)
    return c2cp

def flnc_zmw_to_consensus_primer(flnc_zmws, flnc_z2c, c2cp):
    """
    flnc_zmws: a list of flnc zmws
    flnc_z2c: defaultdict({flnc_zmw: (c_prefix, cid)})
    c2cp: {(c_prefix, cid): consensus_primer}, if consensus_primer is None, meaning could not get consensus primer.
    return # dict{'moive/zmw': consensus_primer}
    ...doctest:
        >>> zmws = [101, 102, 103]
        >>> z2c = {101: 'c1', 102: 'c2', 103: None}
        >>> c2cp = defaultdict(lambda: None)
        >>> c2cp['c1'] = 'cp0'
        >>> dict(flnc_zmw_to_consensus_primer(zmws, z2c, c2cp))
        {101: 'cp0', 102: 'cid_no_cprimer', 103: 'no_cid_no_primer'}
    """
    flnc_zmw2cp = defaultdict(lambda: None) # {flnc_zmw: consensus_primer}
    for zmw in flnc_zmws:
        if zmw in flnc_z2c.keys() and flnc_z2c[zmw] is not None and len(flnc_z2c[zmw]) != 0: # this flnc zmw has associated clusters
            if c2cp[flnc_z2c[zmw]] is None:
                flnc_zmw2cp[zmw] = 'cid_no_cprimer'
            else:
                flnc_zmw2cp[zmw] = c2cp[flnc_z2c[zmw]]
        else:
            flnc_zmw2cp[zmw] = 'no_cid_no_primer'
    return flnc_zmw2cp


def nfl_zmw_to_consensus_primer(nfl_zmws, nfl_z2c, c2cp, get_consensus_func):
    """
    nfl_zmws: a list of nfl zmws
    nfl_z2c: {nfl_zmw: [(c_prefix, cid),...]}
    c2cp: {(c_prefix, cid): consensus_primer}
    get_consensus_func: function to return a consensus item out of multiple items
    ...doctest:
        >>> zmws = [101, 102, 103, 104, 105]
        >>> z2c = {101: ['c0', 'c1', 'c2'], 102: ['c0', 'c3'], 103: [], 104: ['c4']}
        >>> c2cp = {'c0': 0, 'c1': 0, 'c2': 0, 'c3': 1, 'c4': None}
        >>> d = nfl_zmw_to_consensus_primer(zmws, z2c, c2cp, get_most_common_item)
        >>> [(k, d[k]) for k in sorted(d.keys())]
        [(101, 0), (102, None), (103, 'no_cid_no_cprimer'), (104, 'cid_no_cprimer'), (105, 'no_cid_no_cprimer')]
    """
    nfl_zmw2cp = defaultdict(lambda: None) # {nfl_zmw: consensus_primer}
    for nfl_zmw in nfl_zmws:
        if nfl_zmw in nfl_z2c.keys() and nfl_z2c[nfl_zmw] is not None and len(nfl_z2c[nfl_zmw]) != 0: # this nfl zmw has associated clusters
            cps = [c2cp[cid] for cid in nfl_z2c[nfl_zmw] if c2cp[cid] is not None]
            if cps:
                nfl_zmw2cp[nfl_zmw] = get_consensus_func(cps)
            else:
                nfl_zmw2cp[nfl_zmw] = 'cid_no_cprimer'
        else:
            nfl_zmw2cp[nfl_zmw] = 'no_cid_no_cprimer'
    return nfl_zmw2cp

def yield_nfl_zmw_to_consensus_primer(nfl_zmws, nfl_z2c, c2cp, get_consensus_func):
    """
    nfl_zmws: a list of nfl zmws
    nfl_z2c: {nfl_zmw: [(c_prefix, cid),...]}
    c2cp: {(c_prefix, cid): consensus_primer}
    get_consensus_func: function to return a consensus item out of multiple items
    ...doctest:
        >>> zmws = [101, 102, 103, 104, 105]
        >>> z2c = {101: ['c0', 'c1', 'c2'], 102: ['c0', 'c3'], 103: [], 104: ['c4']}
        >>> c2cp = {'c0': 0, 'c1': 0, 'c2': 0, 'c3': 1, 'c4': None}
        >>> [r for r in yield_nfl_zmw_to_consensus_primer(zmws, z2c, c2cp, get_most_common_item)]
        [(101, 0), (102, None), (103, 'no_cid_no_cprimer'), (104, 'cid_no_cprimer'), (105, 'no_cid_no_cprimer')]
    """
    for nfl_zmw in nfl_zmws:
        if nfl_zmw in nfl_z2c.keys() and nfl_z2c[nfl_zmw] is not None and len(nfl_z2c[nfl_zmw]) != 0: # this nfl zmw has associated clusters
            cps = [c2cp[cid] for cid in nfl_z2c[nfl_zmw] if c2cp[cid] is not None]
            if cps:
                yield (nfl_zmw, get_consensus_func(cps))
            else:
                yield (nfl_zmw, 'cid_no_cprimer')
        else:
            yield (nfl_zmw, 'no_cid_no_cprimer')


def dump_d_to_json(d, o_prefix):
    json_fn = o_prefix + '.json'
    json.dump(d, open(json_fn, 'w'))


def dump_d_to_pickle(d, o_prefix):
    pickle_fn = o_prefix + '.pickle'
    pickle.dump(d, open(pickle_fn, 'w'))

def write_dict(d, o_prefix, headers):
    """write dict to both json and csv files
    ...doctest:
        >>> d = {1: 'a', 2: 'b'}
        >>> p = '/home/UNIXHOME/yli/tmp/test_write_dict'
        >>> write_dict(d, p, ['key', 'val'])
        >>> f1, f2, f3 = p + '.csv', p + '.json', p + '.pickle'
        >>> open(f1, 'r').readlines() == ['key\\tval\\n', '1\\ta\\n', '2\\tb\\n']
        True
        >>> json.load(open(f2, 'r'))
        {u'1': u'a', u'2': u'b'}
        >>> pickle.load(open(f3, 'r'))
        {1: 'a', 2: 'b'}
        >>> os.remove(f1)
        >>> os.remove(f2)
    """
    dump_d_to_json(d, o_prefix)
    dump_d_to_pickle(d, o_prefix)
    csv_fn = o_prefix + '.csv'
    with open(csv_fn, 'w') as writer:
        writer.write("\t".join(headers) + "\n")
        for k, v in d.iteritems():
            writer.write("%s\t%s\n" % (k, v))


def get_all_z2c(c_prefix_to_flnc_pickle_fn_dict, c_prefix_to_partial_pickle_fn_dict, lazy=False):
    if lazy:
        print 'Step 2: lazy zmw_to_cid_from_flnc_pickle_fn'
        flnc_z2c = pickle.load(open('flnc_z2c.json', 'r'))
        print 'Step 3: lazy zmw_to_cids_from_partial_pickle_fn'
        nfl_z2c = pickle.load(open('nfl_z2c.json', 'r'))
        return flnc_z2c, nfl_z2c

    print 'Step 2: zmw_to_cid_from_flnc_pickle_fn'
    flnc_z2c = zmw_to_cid_from_flnc_pickle_fns(c_prefix_to_flnc_pickle_fn_dict)
    print 'Step 3: zmw_to_cids_from_partial_pickle_fn'
    nfl_z2c = zmw_to_cids_from_partial_pickle_fns(c_prefix_to_partial_pickle_fn_dict) # dict{'movie/zmw': [cid]}, from ice_partial pickle file
    write_dict(flnc_z2c, o_prefix='flnc_z2c', headers=['flnc_zmw', 'cid'])
    write_dict(nfl_z2c, o_prefix='nfl_z2c', headers=['nfl_zmw', 'cids'])
    return flnc_z2c, nfl_z2c


def get_all_z2p(lazy):
    """Get all z2p zmw_to_primer dict"""
    if lazy:
        print 'lazy get_all_z2p'
        flnc_z2p = pickle.load(open('flnc_z2p.pickle', 'r'))
        nfl_z2p = pickle.load(open('nfl_z2p.pickle', 'r'))
        return flnc_z2p, nfl_z2p

    draft_fa, flnc_fa_fn, nfl_fa_fn = 'data/draft.fasta', 'data/flnc.fasta', 'data/nfl.fasta'
    flnc_z2p = zmw_to_primer_from_fasta(flnc_fa_fn) # dict{'movie/zmw': int(primer)}, where primer=0, 1 or -1(None)
    nfl_z2p = zmw_to_primer_from_fasta(nfl_fa_fn) # dict{'movie/zmw': int(primer)}, where primer=0, 1 or -1(None)
    write_dict(flnc_z2p, o_prefix='flnc_z2p', headers=['flnc_zmw', 'classify_primer'])
    write_dict(nfl_z2p, o_prefix='nfl_z2p', headers=['nfl_zmw', 'classify_primer'])
    return flnc_z2p, nfl_z2p


def get_c2cp(flnc_c2z, z2p, lazy):
    if lazy:
        c2cp = pickle.load(open('c2cp.pickle', 'r'))
        return c2cp
    c2cp = cid_to_consensus_primer(flnc_c2z, z2p, get_most_common_item) #{(c_prefix, cid): consensus_primer}
    write_dict(c2cp, o_prefix='c2cp', headers=['cid', 'consensus_primer'])
    return c2cp


def get_flnc_c2z(c_prefix_to_flnc_pickle_fn_dict, lazy):
    if lazy:
        return pickle.load(open('flnc_c2z.pickle', 'r'))
    flnc_c2z = cid_to_zmws_from_flnc_pickle_fns(c_prefix_to_flnc_pickle_fn_dict) # dict{(c_prefix, cid): [zmws]} , from flnc pickle
    write_dict(flnc_c2z, o_prefix='flnc_c2z', headers=['cid', 'flnc_zmws'])
    return flnc_c2z


def get_flnc_z2cp(flnc_z2p, flnc_z2c, c2cp, lazy):
    if lazy:
        return pickle.load(open('flnc_z2cp.pickle', 'r'))
    flnc_z2cp = flnc_zmw_to_consensus_primer(flnc_z2p.keys(), flnc_z2c, c2cp) # dict{'moive/zmw': consensus_primer}
    write_dict(flnc_z2cp, o_prefix='flnc_z2cp', headers=['flnc_zmw', 'consensus_primer'])
    return flnc_z2cp


