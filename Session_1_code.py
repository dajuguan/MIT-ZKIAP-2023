from math import floor
from random import random
import numpy as np
from functools import reduce

PRIME = 273389558745553615023177755634264971227
GEN   = 191981998178538467192271372964660528157
ROUND_OF_VERIFY = 100
# args: x:secret, g: generator, p: large prime, b: random number
# returns: (y: residue, proof: pf)
def dlogProof(x, g, p):
    y = cal_exp(g,x) % p #residue
    r_arr =  np.array([]) #randoms
    for i in range(ROUND_OF_VERIFY):
        r_arr = np.append(r_arr, floor(random()*(p-1)))
    h = np.array([])  #proof 1
    for r in r_arr:
        h = np.append(h,cal_exp(g,r) % p) 
    b = pseudoRandomBits(h,g,p,ROUND_OF_VERIFY)
    s = (r_arr + b*x) % (p -1)    # proof 3
    pf = (h,s)
    return (y, pf)

def cal_exp(x,exp):
    half = x**(floor(exp/2))
    if exp % 2 == 0:
        return half*half
    else:
        return x * half * half

def pseudoRandomBits(h,g,p,ROUND_OF_VERIFY):
    seed = np.sum(h)
    hash = cal_exp(g,seed) % p
    b = np.array([])
    for i in range(0,ROUND_OF_VERIFY):
        temp = hash & i
        bit = 0
        if temp > 0:
            bit = 1
        b = np.append(b, bit)
    return b

# Verifier computes A^s(mod p) which should equal hy^b(mod p).
def verify(y,g,p,pf):
    (h,s) = pf
    b = pseudoRandomBits(h,g,p,ROUND_OF_VERIFY)
    lhs = (g**s) % p
    rhs = (h*(y**b)) % p
    res = reduce(lambda l,r: l and r, lhs == rhs) 
    return res

def test():
    p = 31
    g = 3
    x = 17

    (y,pf) = dlogProof(x,g,p)
    res = verify(y,g,p,pf)
    print("The proof is: ", res)
    res = verify(y-1,g,p,pf)
    print("The proof is: ", res)

    p = 67
    g = 2
    x = 10
    (y,pf) = dlogProof(x,g,p)
    res = verify(y,g,p,pf)
    print("The proof is: ", res)
    res = verify(y-1,g,p,pf)
    print("The proof is: ", res)

test()