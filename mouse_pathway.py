import sys
from collections import Counter
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact #, chi2_contingency
from os.path import dirname, join

def get_id2sym( homology_fn = None ):
    id2sym = {}

    if not homology_fn  :
        homology_fn = join( dirname( sys.argv[0], 'ensembl_homology.data') )

    for l in homology_fn:
        line = l.split('\t')
        h, m = line[0].strip(), line[7].strip()
        if h and m :
            id2sym[h] = m
    return id2sym

def pathway2fisherp( nlist, npathway, nlistnpathway, ntotal, verbose=False ):
    if verbose :
        print( '# of genes in the list:', nlist )
        print( "# of genes in the pathway:", npathway )
        print( "# of genes in the both:", nlistnpathway )
        print( "# of genes total:", ntotal )

    a = nlistnpathway
    b = npathway - a
    c = nlist - a
    d = ntotal - a - b - c
    if verbose:
        print([a,b],[c,d], sep="\n")
    return fisher_exact([[a,b],[c,d]])[1]


def build_gene2path( pathfn, gid_map={}, minsize_cutoff=1 ) :
    g2p = {}

    for l in open(pathfn):
        line = l.split('\t')
        pid = line[0]
        genes = set(line[9].split())
        if len(genes) < minsize_cutoff :
            continue

        for g in genes :
            if gid_map :
                g = gid_map[g]

            if g in g2p :
                g2p[g].add(pid)
            else :
                g2p[g] = set()
                g2p[g].add(pid)
    return g2p

def build_pathways( pathfn, minsize_cutoff=1 ):
    pid2pathway = {}
    pid2genecounts = {}

    for l in open( pathfn ):
        line = l.split('\t')
        pid = line[0]
        if int(line[8]) < minsize_cutoff :
            continue

        pid2pathway[pid] = l
        pid2genecounts[pid] = len(line[9].split())

    return pid2pathway, pid2genecounts


def print_output( fn, gnames, pvalues, adj_pvalues, pathway_ids, pathways, cutoff=0.001) :
    nlist=len(gnames)
    ofp = open( fn, 'w')
    for i,(p, adj_p, pid) in enumerate(zip( pvalues, adj_pvalues, pathway_ids)):
        #pathways[pid]
        line = pathways[pid].split('\t')
        pnames = set(line[9].split())
        both = gnames.intersection(pnames)


        print(i, p, adj_p,
              len(both)/len(pnames)-nlist/total,
              len(both),
              len(pnames),
              nlist,
              total,
              line[1],
              line[2],
              line[3],
              line[7],
              " ".join(both),
              file=ofp, sep="\t")
        if p > cutoff :
            break
    ofp.close()



if __name__ == '__main__' :
    from os.path import dirname, join
    mouse_pathway_fn = join( dirname(sys.argv[0]), 'warpath.mouse.txt' )

    input_list = sys.argv[1]
    output_fn = sys.argv[2]

    g2p = build_gene2path( mouse_pathway_fn )
    pathways, genecounts = build_pathways(mouse_pathway_fn)
    id2sym = get_id2sym()

    #first try ID to symbol conversion
    genes = set([id2sym[g] for g in [ l.split(".")[0] for l in open(input_list) ] if g in id2sym ])
    #if the conversion fails
    #try them as symbols
    if not genes :
        symbols = set(list(id2sym.values()))
        genes = set( [g in [ l.split(".")[0] for l in open(input_list) ] if g in symbols])

    assert genes #test if the genes are empty or not.

    pathways_in = [g2p[g] for g in genes if g in g2p]

    c = Counter()
    for gs in pathways_in :
        for g in gs :
            c[g]+=1

    nlist = len(pathways_in)
    total = len(g2p)

    results = [ (pathway2fisherp(nlist, genecounts[pid],
                                 both, total), pid ) for pid, both in c.items() ]
    sorted_results = [(p,v) for p,v in sorted(results)]
    pvalues = [p for p,v in sorted_results]
    adj_pvalues = multipletests(pvalues, method='fdr_bh', is_sorted=True)[1]
    pathway_ids = [p for _,p in sorted_results]

    gnames = set(g for g in genes if g in g2p)

    print_output( output_fn, gnames, pvalues, adj_pvalues, pathway_ids, pathways )
