from __future__ import absolute_import
from builtins import map
import os
import re
import sys
import glob
import csv
import argparse
import presto.sifting as sifting
from operator import itemgetter, attrgetter

# Note:  You will almost certainly want to adjust
#        the following variables for your particular search

# glob for ACCEL files
globaccel = "*ACCEL_*0"
# globaccel = "*ACCEL*"
# glob for .inf files
globinf = "*DM*.inf"
# In how many DMs must a candidate be detected to be considered "good"
min_num_DMs = 1
# Lowest DM to consider as a "real" pulsar
low_DM_cutoff = 1.0
# Ignore candidates with a sigma (from incoherent power summation) less than this
sifting.sigma_threshold = 4.0
# Ignore candidates with a coherent power less than this
sifting.c_pow_threshold = 100.0

# If the birds file works well, the following shouldn't
# be needed at all...  If they are, add tuples with the bad
# values and their errors.
#                (ms, err)
sifting.known_birds_p = []
#                (Hz, err)
sifting.known_birds_f = []

# The following are all defined in the sifting module.
# But if we want to override them, uncomment and do it here.
# You shouldn't need to adjust them for most searches, though.

# How close a candidate has to be to another candidate to                
# consider it the same candidate (in Fourier bins)
sifting.r_err = 1.1
# Shortest period candidates to consider (s)
sifting.short_period = 0.0005
# Longest period candidates to consider (s)
sifting.long_period = 15.0
# Ignore any candidates where at least one harmonic does exceed this power
sifting.harm_pow_cutoff = 8.0

#--------------------------------------------------------------

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='presto accel sifter for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')

    args = parser.parse_args()

    # Try to read the .inf files first, as _if_ they are present, all of
    # them should be there.  (if no candidates are found by accelsearch
    # we get no ACCEL files...
    path = f"{args.out_dir}/inj_{args.injection_number:06}/processing/PRESTO/ACCEL"
    inffiles = glob.glob(globinf, root_dir=path)
    # accelsearch writes 'root_ACCEL_<zmax>' (the text list we want) alongside
    # 'root_ACCEL_<zmax>.cand'/'.txtcand'; select by excluding those extensions
    # rather than assuming zmax ends in a '0'.
    candfiles = [f for f in glob.glob(globaccel, root_dir=path)
                 if not f.endswith(('.cand', '.txtcand', '.inf'))]

    if not (inffiles and candfiles):
        sys.exit('No PRESTO candidates to sift.')

    # DMs are no longer unique across the search: a segmented, multi-downsample
    # run produces the same DM once per (segment, downsample). Dedupe so a DM is
    # counted once by remove_DM_problems.
    dms = sorted(set(float(x.split("DM")[-1].split(".inf")[0]) for x in inffiles))
    dmstrs = ["%.2f"%x for x in dms]

    # sifting reads the observation length out of the .inf appended to each ACCEL
    # file, so every ACCEL file needs its OWN .inf. Pair them on the full rootname
    # (which carries segment + downsample), never on the DM string alone.
    accel_re = re.compile(r'_ACCEL_\d+(?:_JERK_\d+)?$')

    candfiles_new = []
    for candf in candfiles:
        source_file = f"{path}/{accel_re.sub('', candf)}.inf"
        destination_file = f"{path}/{candf}"
        if os.path.exists(source_file):
            with open(source_file, "r") as src, open(destination_file, "a") as dest:
                for line in src:
                    dest.write(line)
            candfiles_new.append(destination_file)
        else:
            print(f"No .inf found for {candf}, skipping.")

    cands = sifting.read_candidates(candfiles_new) # edit RS

    # Remove candidates that are duplicated in other ACCEL files
    if len(cands):
        cands = sifting.remove_duplicate_candidates(cands)

    # Remove candidates with DM problems
    if len(cands):
        cands = sifting.remove_DM_problems(cands, min_num_DMs, dmstrs, low_DM_cutoff)

    # Remove candidates that are harmonically related to each other
    # Note:  this includes only a small set of harmonics
    if len(cands):
        cands = sifting.remove_harmonics(cands)

    # def to_file(self, candfilenm):
    #     candfile = open(candfilenm, "w")
    #     candfile.write("#" + "file:candnum".center(66) + "DM".center(9) +
    #                     "SNR".center(8) + "sigma".center(8) + "numharm".center(9) +
    #                     "ipow".center(9) + "cpow".center(9) +  "P(ms)".center(14) +
    #                     "r".center(12) + "z".center(8) + "numhits".center(9) + "\n")
    #     for goodcand in self.cands:
    #         candfile.write("%s (%d)\n" % (str(goodcand), len(goodcand.hits)))
    #         if (len(goodcand.hits) > 1):
    #             goodcand.hits.sort(key=lambda cand: float(cand[0]))
    #             for hit in goodcand.hits:
    #                 numstars = int(hit[2]/3.0)
    #                 candfile.write("  DM=%6.2f SNR=%5.2f Sigma=%5.2f   "%hit + \
    #                                 numstars*'*' + '\n')
    #     if candfilenm is not None:
    #         candfile.close()


    def get_jerk_value(filename, candnum):
        # accelsearch's own jerk-search ('w') column isn't exposed by presto.sifting's
        # candidate class (it only tracks r/z), so for jerk-search runs (filenames
        # containing '_JERK_') re-parse the raw candidate line ourselves.
        if '_JERK_' not in filename:
            return 0.0
        try:
            with open(filename, "r") as f:
                for line in f:
                    if sifting.fund_re.match(line):
                        split_line = line.split()
                        if int(split_line[0]) == candnum:
                            return float(split_line[11].split("(")[0])
        except (OSError, IndexError, ValueError):
            pass
        return 0.0

    def to_file(self, candfilenm):
        with open(candfilenm, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)

            # Header (same columns as before, plus 'w' for jerk-search candidates)
            writer.writerow([
                "file",
                "candnum",
                "DM",
                "SNR",
                "sigma",
                "numharm",
                "ipow",
                "cpow",
                "P(ms)",
                "r",
                "z",
                "w",
                "numhits"
            ])

            for goodcand in self.cands:
                fields = str(goodcand).split()
                file, candnum = fields[0].split(':')
                fields[0] = file
                fields.insert(1, candnum)
                fields.append(get_jerk_value(os.path.join(goodcand.path, goodcand.filename), int(candnum)))
                fields.append(len(goodcand.hits))
                writer.writerow(fields)


    cands.sort(key=attrgetter('sigma'), reverse=True)
    to_file(cands, candfilenm=f"{args.out_dir}/inj_{args.injection_number:06}/processing/PRESTO/PRESTO_candidates.csv")


