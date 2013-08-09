import random

def myran(a,b,Nrand,distrib='uniform'):
    nums = [None] * Nrand
    for i in range(Nrand):
        if distrib == 'gaussian':
            nums[i]=random.gauss(a,b)
        else: 
            nums[i]=random.uniform(a,b)
    return nums
