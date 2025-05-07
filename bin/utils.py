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

def getMajorAnn(ann,cl,is_ref):
    d = pd.crosstab(ann[is_ref],cl[is_ref])
    res = pd.DataFrame({"cluster" : d.columns,
                    "major_ann" : [d.iloc[:,i].idxmax() for i in range(d.shape[1])],
                    "major_ann_cnt" : [d.iloc[:,i].max() for i in range(d.shape[1])],
                    "ref_size" : d.sum(0)
                   })
    sizes = pd.crosstab(cl,is_ref)
    res['query_size'] = sizes.loc[res.index,False]
    res['major_frac'] = res.major_ann_cnt/res['ref_size']
    res['query_frac'] = res['query_size']/(res['query_size']+res['ref_size'])
    return res