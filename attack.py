import clt

if __name__=="__main__":
    lam = 10
    k = 5

    params = clt.MMP.set_params(lam, k)
    mmap = clt.MMP(params)