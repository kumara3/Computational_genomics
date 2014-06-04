import argparse
import random
import SBH

def CircularSpectrum(s, l):
    s2 = s+s
    S = [s2[i:i+l] for i in range(len(s))]
    random.shuffle(S)
    return S

def LinearSpectrum(s, l):
    S = [s[i:i+l] for i in range(len(s)-l+1)]
    random.shuffle(S)
    return S



def random_test(n, l):
    s = "".join([random.choice("ACGT") for i in range(n)])
    S = CircularSpectrum(s, l)
    s2 = SBH.SBH(S)
    S2 = CircularSpectrum(s2, l)

    if sorted(S) == sorted(S2):
        return 1

    # Partial credit if student forgot to account the circular genome
    S = LinearSpectrum(s, l)
    s2 = SBH.SBH(S)
    S2 = LinearSpectrum(s2, l)

    if sorted(S) == sorted(S2):
        return 0.8
    return 0


def avg_test(num_trials, str_len, l):
    T = [random_test(str_len, l) for n in range(num_trials)]
    return sum(T) / float(len(T))

    
def run_test(n, show_results = True, RNGseed = None):
    if RNGseed:
        random.seed(RNGseed)

    T1 = avg_test(n, 20, 5)
    T2 = avg_test(n, 50, 5)
    T3 = avg_test(n, 100, 10)
    T4 = avg_test(n, 500, 15)
    score = (T1 + T2 + T3 + T4)/4.0
    print("SCORE: %5.2f" % (100*score))

    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Testing dnaSeq project')
    parser.add_argument('-n', type = int, help = "Each result will be averaged over n runs", default = 100)
    parser.add_argument('-s', type = int, help = "Seed for random number generator", default = None)
    parser.add_argument('-d', '--display', action = "store_true", help = "display indvidual test results", default = False)
    args = parser.parse_args()
    run_test(args.n, RNGseed = args.s, show_results = args.display)




