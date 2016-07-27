#!/usr/bin/env python

from argparse import ArgumentParser

WAV_MIN = 5100
WAV_MAX = 6300
TARG_MAX = 310

if __name__ == '__main__':
    psr = ArgumentParser(description="Generate script for parallelization")
    psr.add_argument('outpath', type=str, help="Path to output file")
    args = psr.parse_args()

    f = open(args.outpath, "w")

#     for i in range(WAV_MIN, WAV_MAX, 100):
#         for j in range(TARG_MAX):
#             f.write("source ~/.bash_profile; python /home/syee/specmatchemp-working/specmatchemp/tests/library_match_parallel.py \
# '/home/syee/specmatchemp-working/specmatchemp/lib/library_reduced.h5' \
# '/home/syee/specmatchemp-working/specmatchemp/results/{0:d}_results_{1:d}.csv' {1:d} {0:d} 100\n".format(i, j))
#     f.close()
    
    for i in range(WAV_MIN, WAV_MAX, 100):
        for j in range(6, 9):
            f.write("source ~/.bash_profile; python /home/syee/specmatchemp-working/specmatchemp/tests/library_match_lincomb.py " + \
                "'/home/syee/specmatchemp-working/specmatchemp/lib/library_reduced.h5' "+\
                "'/home/syee/specmatchemp-working/specmatchemp/results/{0:d}_results.csv' ".format(i)+\
                "'/home/syee/specmatchemp-working/specmatchemp/results/{0:d}_results_lincomb_{1:d}.csv' ".format(i,j)+\
                "{1:d} {0:d} 100\n".format(i,j))
    f.close()
