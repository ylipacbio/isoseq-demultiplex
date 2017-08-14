#!/usr/bin/env python
import sys
import os
import os.path as op
from argparse import ArgumentParser
from .utils import *
from .cluster_to_consensus_primer import *

def z2cp(flnc_z2c, c2cp):
    pass


def get_c2cp_from_cluster_dict_reader(cluster_dict_reader, min_fraction):
    """
    ...doctest:
        >>> s = 'cid\\t[movie1/1,movie1/2]\\t[movie2/3,movie3/4]\\t[0,0]\\t[1,None]\\t0'
        >>> dict(get_c2cp_from_cluster_reader([s], 0.6))
        {'cid': 0}
    """
    c2cp = defaultdict(lambda: None)
    for s in cluster_dict_reader:
        obj = ClusterDict.fromString(s, min_fraction)
        c2cp[obj.cid] = obj.consensus_primer
    return c2cp


def get_c2cp_from_cluster_dict_fn(cluster_dict_fn, min_fraction):
    return get_c2cp_from_cluster_dict_reader(open(cluster_dict_fn, 'r'), min_fraction)

if __name__ == "__main__":
    flnc_z2c_fn = 'flnc_z2c.csv'
    nfl_z2c_fn = 'nfl_z2c.csv'
    cluster_dict_fn = 'cluster_dict.csv'
