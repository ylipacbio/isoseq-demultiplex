"""
    Step 1: map zmw (represented by 'movie/zmw') to primer based on classfied outputs, where primer can be 0, 1 or -1 (None).
    Step 2: map FLNC zmw (aka 'movie/zmw') to cluster id (aka cid) since one FLNC zmw can be assigned to at most one cluster id.
    Step 3: map NFL zmw (aka 'movie/zmw') to cluster ids (aka cids) since one NFL zmw can be assigned to multiple cluster ids.
    Step 4: map cluster id to all supportive FLNC zmws
    Step 5: map cluster id to consensus primer of all supportive FLNC zmws, where consensus primer must be agreed by >= min_fraction supportive FLCN zmws. If there is no consensus, mark this cluster id as ambiguous.
    Step 6: for all FLNC zmws, get predicted primer from their associated cluster:
    Step 6.1 if a FLNC zmw has an associated cluster,
    Step 6.1.1 if the associated cluster has a consensus primer, use the consensus primer as predicted primer
    Step 6.1.2 if the associated cluster has no consensus primer, use 'cid_no_cprimer' as the predicted primer
    Step 6.2 if a FLNC zmw does not have an associated cluster, use 'no_cid_no_cprimer' as the predicted primer

    Step 7: for all NFL czmws, get predicted primer from all their associated clusters:
    Step 7.1.1 if a NFL zmw has a list of associated clusters, use consensus primer of consensus primers as predicted primer, where consensus primer must be agreed by >= min_fraction clusters.
    Step 7.1.2 if a NFL zmw has no associated clusters, use 'no_cid_no_cprimer' as the predicted primer
"""
from .utils import *


def run_all(sep_dir):
    lazy = False
    # Step 1
    print 'Step 1: get_all_z2p'
    flnc_z2p, nfl_z2p = get_all_z2p(lazy=lazy)
    z2p = flnc_z2p.copy()
    z2p.update(nfl_z2p)

    bin_dirs = ['%sto%skb_part0' % (i, i+1) for i in range(0, 15)] + ['5to6kb_part1']
    # ['0to1kb_part0', ..., '14to15kb_part0', '5to6kb_part1']
    nfl_pickle_fns = [op.join(sep_dir, bin_dir, 'cluster_out/output/map_noFL/nfl.all.partial_uc.pickle') for bin_dir in bin_dirs]
    flnc_pickle_fns = [op.join(sep_dir, bin_dir, 'cluster_out/output/final.pickle') for bin_dir in bin_dirs]
    c_prefix_to_partial_pickle_fn_dict = dict(zip(bin_dirs, nfl_pickle_fns))
    c_prefix_to_flnc_pickle_fn_dict = dict(zip(bin_dirs, flnc_pickle_fns))

    print 'Step 2 and 3: get_all_z2c'
    flnc_z2c, nfl_z2c = get_all_z2c(c_prefix_to_flnc_pickle_fn_dict, c_prefix_to_partial_pickle_fn_dict, lazy=lazy)

    print 'Step 4: get_fnc_c2z'
    flnc_c2z = get_flnc_c2z(c_prefix_to_flnc_pickle_fn_dict=c_prefix_to_flnc_pickle_fn_dict, lazy=lazy)

    print 'Step 5: cid_to_consensus_primer'
    c2cp = get_c2cp(flnc_c2z, z2p, lazy=lazy)

    print 'Step 6: flnc_zmw_to_consensus_primer'
    flnc_z2cp = get_flnc_z2cp(flnc_z2p, flnc_z2c, c2cp, lazy=lazy)

    print 'Step 7: nfl_zmw_to_consensus_primer'
    #nfl_z2cp = nfl_zmw_to_consensus_primer(nfl_zmws=nfl_z2p.keys(), nfl_z2c=nfl_z2c, c2cp=c2cp, get_consensus_func=get_most_common_item)
    with open('nfl_z2cp.csv', 'w') as writer:
        writer.write('\t'.join(['flnc_zmw', 'consensus_primer']) + '\n')
        for r in yield_nfl_zmw_to_consensus_primer(nfl_z2p.keys(), nfl_z2c, c2cp, get_most_common_item):
            writer.write('\t'.join(r) + '\n')

            writer.write('\t'.join(r) + '\n')

def run_one_bin(sep_dir): #flnc_pickle_fn, partial_pickle_fn, c_prefix, flnc_z2p, nfl_z2p, z2p):
    flnc_z2p, nfl_z2p = get_all_z2p(lazy=False)
    z2p = flnc_z2p.copy()
    z2p.update(nfl_z2p)

    bin_dir = '13to14kb_part0'
    c_prefix = bin_dir
    job_dir = op.join(sep_dir, bin_dir)
    flnc_pickle_fn = op.join(job_dir, 'cluster_out/output/final.pickle')
    partial_pickle_fn = op.join(job_dir, 'cluster_out/output/map_noFL/nfl.all.partial_uc.pickle')

    #run_one_bin(final_pickle_fn, ice_partial_pickle_fn, c_prefix=bin_dir, flnc_z2p, nfl_z2p, z2p)

    # Step 2
    print 'Step 2: zmw_to_cid_from_flnc_pickle_fn'
    flnc_z2c = zmw_to_cid_from_flnc_pickle_fn(c_prefix, flnc_pickle_fn) # dict{'movie/zmw': cid} from final.pickle
    write_dict(flnc_z2c, o_prefix='flnc_z2c', headers=['flnc_zmw', 'cid'])

    # Step 3
    print 'Step 3: zmw_to_cids_from_partial_pickle_fn'
    nfl_z2c = zmw_to_cids_from_partial_pickle_fn(c_prefix, partial_pickle_fn) # dict{'movie/zmw': [cid]}, from ice_partial pickle file
    write_dict(nfl_z2c, o_prefix='nfl_z2c', headers=['nfl_zmw', 'cids'])

    # Step 4
    print 'Step 4: cid_to_zmws_from_flnc_pickle_fn'
    flnc_c2z = cid_to_zmws_from_flnc_pickle_fn(c_prefix, flnc_pickle_fn) # dict{(c_prefix, cid): [zmws]} , from flnc pickle
    write_dict(flnc_c2z, o_prefix='flnc_c2z', headers=['cid', 'flnc_zmws'])

    # Step 5
    print 'Step 5: cid_to_consensus_primer'
    c2cp = cid_to_consensus_primer(flnc_c2z, z2p, get_most_common_item) #{(c_prefix, cid): consensus_primer}
    write_dict(c2cp, o_prefix='c2cp', headers=['cid', 'consensus_primer'])
    # Step 6
    print 'Step 6: flnc_zmw_to_consensus_primer'
    flnc_z2cp = flnc_zmw_to_consensus_primer(flnc_z2p.keys(), flnc_z2c, c2cp) # dict{'moive/zmw': consensus_primer}
    write_dict(flnc_z2cp, o_prefix='flnc_z2cp', headers=['flnc_zmw', 'consensus_primer'])
    # Step 7
    print 'Step 7: nfl_zmw_to_consensus_primer'
    #nfl_c2z = cid_to_zmws_from_partial_pickle_fn(c_prefix, partial_pickle_fn)
    #nfl_z2cp = nfl_zmw_to_consensus_primer(nfl_zmws=nfl_z2p.keys(), nfl_z2c=nfl_z2c, c2cp=c2cp, get_consensus_func=get_most_common_item)
    with open('nfl_z2cp.csv', 'w') as writer:
        writer.write('\t'.join(['flnc_zmw', 'consensus_primer']) + '\n')
        for r in yield_nfl_zmw_to_consensus_primer(nfl_z2p.keys(), nfl_z2c, c2cp, get_most_common_item):
            writer.write('\t'.join(r) + '\n')


