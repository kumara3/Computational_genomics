import sys
sys.path.append("../lib")

import argparse 
import random
import SBH

def spectrum(s, l):
    s2 = s+s
    S = [s2[i:i+l] for i in range(len(s))]
    random.shuffle(S)
    return S



def random_test(n, l):
    s = "".join([random.choice("ACGT") for i in range(n)])
    S = spectrum(s, l)
    s2 = SBH.SBH(S)
    S2 = spectrum(s2, l)

    return sorted(S) == sorted(S2)


    
def run_test(n, l, t, RNGseed = None):
    if RNGseed:
        random.seed(RNGseed)

    T = [random_test(n,l) for i in range(t)]
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Testing dnaSeq project')
    parser.add_argument('-n', type = int, help = "Length of string", default = 100)
    parser.add_argument('-l', type = int, help = "l value for test", default = 20)
    parser.add_argument('-s', type = int, help = "Seed for random number generator", default = None)
    parser.add_argument('-t', type = int, help = "Number of trials", default = 30)
    args = parser.parse_args()
    run_test(args.n, RNGseed = args.s, l = args.l, t = args.t)




