import clt

if __name__=="__main__":
    lam = 10
    k = 5

    params = clt.MMP.set_params(lam, k)
    mmap = clt.MMP(params)

# attack psuedocode

# generate a bunch of level-k encodings of 0 c_i
# multiply all c_i by p_zt yielding w_i for each c_i
# compute lattice L orthogonal to w = (w_1, ..., w_l)
    # run LLL/BKZ on basis B of L
# get all vectors {u} in B that are short enough
# compute lattice L_r orthogonal to L with basis B_r
    # run LLL/BKZ on basis B_r

# for r in reduced_B_r
    # compute s such that <s, r> = 0 (mod x_0)
    # if <s, w> = 0 (mod x_0)
        # continue to next loop iteration
    # compute p = gcd(<s, w>, x_0)
    # p should be a factor of x_0

# repeat algorithm until x_0 is fully factored