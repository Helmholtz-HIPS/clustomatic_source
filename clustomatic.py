#!/usr/bin/python3
import scipy
import scipy.cluster
import os
import sys
import numpy as np
import fastcluster
import subprocess
import shutil
from Bio import SeqIO

def parse_input( path, lib ):
    fasta_sequences = SeqIO.parse(open(path),"fasta")
    index = 0
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        (qcid, qgid) = name.split("_")
        if not lib.get(qcid):
            lib[qcid] = { "index": index }
            index = index + 1
        lib[qcid]["length"] = lib[qcid].get("length", 0) + len( sequence )

def prepare_matrix( lib ):
    return np.ones( int( ( ( len(lib) - 1 ) * ( len(lib) ) ) / 2 ) )


def make_blast( path ):
    pid = os.posix_spawnp( "diamond", [
                                        "diamond",
                                        "makedb",
                                        "--quiet",
                                        "--db",
                                        "{0}.db".format( path ),
                                        "--in",
                                        "{0}".format( path )
                                      ], os.environ )
    os.waitpid( pid, 0 )

def run_blast( path ):
    pid = os.posix_spawnp( "diamond", [
                                        "diamond",
                                        "blastp",
                                        "--quiet",
                                        "--db",
                                        "{0}.db".format( path ),
                                        "-q",
                                        "{0}".format( path ),
                                        "--max-target-seqs",
                                        "0",
                                        "--outfmt",
                                        "6",
                                        "qseqid",
                                        "sseqid",
                                        "nident",
                                        "--out",
                                        "{0}.output".format( path )
                                      ], os.environ )
    os.waitpid( pid, 0 )

def parse_blast( path, matrix, lib ):
    n = len(lib)
    with open( "{0}.output".format( path ), "r" ) as output:
        for pool in generate_query_pools(output):
            track = {}
            q = lib[pool[0]['qcid']] # get query index only once per query change
            i = q['index']
            for line in sorted(pool, key=lambda x: -x['score']): # https://stackoverflow.com/questions/34455594/why-is-sorting-a-python-list-of-tuples-faster-when-i-explicitly-provide-the-key
                if not track.setdefault(line['qcid']+line['scid'], {}).get(line['qgid']) and not track.setdefault(line['qcid']+line['scid'], {}).get(line['sgid']):
                    track[line['qcid']+line['scid']][line['qgid']] = True
                    track[line['qcid']+line['scid']][line['sgid']] = True
                    h = lib[line['scid']]
                    j = h['index']
                    if i < j:
                        k = int( n*(n-1)/2 - (n-i)*((n-i)-1)/2 + j - i - 1 )
                        matrix[k] = matrix[k] - 2*line['score'] / ( q['length'] + h['length'] )
                        if matrix[k] < 0:
                            matrix[k] = 0

def generate_query_pools( stream ):
    # generator for query pools
    pool = []
    for line in stream:
        line = line.split()
        (qcid, qgid) = line[0].split('_')
        (scid, sgid) = line[1].split('_')
        if len(pool) > 0 and pool[-1]['qcid'] != qcid: # reset pool with every 'new' query key
            yield pool
            pool = []
        pool.append( { 'qcid': qcid, 'qgid': qgid, 'scid': scid, 'sgid': sgid, 'score': int(line[2]) } )
    yield pool # one final yield

def build_linkage( matrix ):
    # build a linkage
    linkage = fastcluster.linkage( matrix, method='average', preserve_input=False )
    assert scipy.cluster.hierarchy.is_valid_linkage( linkage ), 'reconstructed linkage not valid'
    return linkage

def form_clusters( linkage ):
    return scipy.cluster.hierarchy.fcluster( linkage, t=sys.argv[2], criterion="distance" )

def clustomatic( path ):

    lib = {}

    parse_input( path, lib )

    make_blast( path )

    run_blast( path )

    matrix = prepare_matrix( lib )

    parse_blast( path, matrix, lib )

    linkage = build_linkage( matrix )

    clusters = form_clusters( linkage )

    for (i, key) in enumerate(lib):
        print ( "{0}\t{1}".format( key, clusters[i] ) )

def main():
    if len(sys.argv) < 3:
        raise SystemExit("Usage: clustomatic.py input.fasta threshold\nFasta protein records must have following naming: >clusterID_proteinID, see example_input.fasta")
    if not os.path.isfile(sys.argv[1]):
        raise SystemExit("File not found: {0}".format( sys.argv[1] ))
    if float( sys.argv[2] ) < 0 or float( sys.argv[2] ) > 1:
        raise SystemExit("Threshold has to be a float between 0 and 1")

    if not shutil.which("diamond"):
        raise SystemExit("Please install diamond (https://github.com/bbuchfink/diamond)")

    clustomatic( sys.argv[1] )

main()
