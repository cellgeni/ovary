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