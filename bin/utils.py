import pandas as pd

def annotateByKnn(data,knn_key,ann_key):
    knn = data.obsp[knn_key]
    res = []
    for row in range(data.shape[0]):
        inxs = knn.indices[knn.indptr[row]:knn.indptr[row+1]]
        inxs = [i for i in inxs if  i != row]
        r = data.obs[ann_key][inxs].value_counts()
        res.append([r.index[0],r.values[0],r.sum()])
    res=pd.DataFrame(res,columns=[ann_key,'major_cnt','total'])
    res['confidence'] = res.major_cnt/res.total
    res.index = data.obs_names
    return res

def getMajorAnn(ann,cl):
    d = pd.crosstab(ann,cl)
    res = pd.DataFrame({"cluster" : d.columns,
                    "major_ann" : [d.iloc[:,i].idxmax() for i in range(d.shape[1])],
                    "major_ann_cnt" : [d.iloc[:,i].max() for i in range(d.shape[1])],
                    "total" : d.sum(0)
                   })
    res['major_frac'] = res.major_ann_cnt/res['total']
    return res

def loadGeneFromGTF(gtf_file):
    genes = pd.read_table(gtf_file,comment='#',header=None)
    genes = genes.loc[genes.loc[:,2]=='gene',:]
    genes['gene_id'] = genes.loc[:,8].str.extract('gene_id "([^"]+)";')
    genes['gene_name'] = genes.loc[:,8].str.extract('gene_name "([^"]+)";')
    return genes

def addGeneToPeaks(adata,gtf,upstream,downstream):
    import snapatac2 as snap
    genes = loadGeneFromGTF(gtf)
    genes.index = genes.gene_id

    g2p = snap.tl.init_network_from_annotation(adata.var_names, gtf, upstream=upstream, downstream=downstream, id_type='gene_id', coding_gene_only=False)
    nodes = pd.DataFrame([x.__dict__  for x in g2p.nodes()])
    nodes['inx'] = range(nodes.shape[0])

    regions = nodes.loc[nodes.type=='region',:]
    regions.index = regions.id

    regions['gene_ids'] = [list(nodes.id[g2p.neighbors(inx)]) for inx in regions.inx]
    regions['gene_names'] = [list(genes.loc[gids,'gene_name']) for gids in regions.gene_ids]
    adata.var['gene_ids'] = adata.var['gene_names'] = pd.NA
    adata.var.loc[regions.id,'gene_ids'] =regions.gene_ids.str.join(',')
    adata.var.loc[regions.id,'gene_names'] =regions.gene_names.str.join(',')

