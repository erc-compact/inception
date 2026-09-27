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

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'PY_general')))
import pipeline_tools as inj_tools

# accelsearch writes 'root_ACCEL_<zmax>' (the text list we want) alongside
# 'root_ACCEL_<zmax>.cand' / '.txtcand'. These are fixed by PRESTO's own output
# naming rather than being a user choice, so they stay in code.
globaccel = "*ACCEL_*"
globinf = "*DM*.inf"

# Everything the sifter thresholds on lives in 'presto_sift_args'. These defaults
# reproduce PRESTO's stock values and are only used when a key is absent.
SIFT_DEFAULTS = {
    # ignore candidates with a sigma (from incoherent power summation) below this
    "sigma_threshold": 4.0,
    # ignore candidates with a coherent power less than this
    "c_pow_threshold": 100.0,
    # ignore candidates where no harmonic exceeds this power
    "harm_pow_cutoff": 8.0,
    # how close (in Fourier bins) two candidates must be to count as the same one
    "r_err": 1.1,
    # shortest / longest period candidates to consider (s)
    "short_period": 0.0005,
    "long_period": 15.0,
    # in how many DMs must a candidate be detected to be considered "good"
    "min_num_DMs": 1,
    # lowest DM to consider as a "real" pulsar
    "low_DM_cutoff": 1.0,
    # known interference to zap, as [value, error] pairs: periods (ms) and freqs (Hz)
    "known_birds_p": [],
    "known_birds_f": [],
}


def apply_sift_args(processing_args):
    """Push the configured thresholds into the sifting module. Returns the few
    that are passed as arguments rather than set as module globals."""
    args = inj_tools.parse_JSON(processing_args).get('presto_sift_args', {})
    sift = dict(SIFT_DEFAULTS, **args)

    for key in ("sigma_threshold", "c_pow_threshold", "harm_pow_cutoff",
                "r_err", "short_period", "long_period"):
        setattr(sifting, key, sift[key])

    # sifting expects (value, error) tuples; JSON can only give nested lists
    sifting.known_birds_p = [tuple(bird) for bird in sift["known_birds_p"]]
    sifting.known_birds_f = [tuple(bird) for bird in sift["known_birds_f"]]

    print(f"sifting with sigma_threshold={sift['sigma_threshold']}, "
          f"c_pow_threshold={sift['c_pow_threshold']}, "
          f"period range {sift['short_period']}-{sift['long_period']} s")

    return sift["min_num_DMs"], sift["low_DM_cutoff"]

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

def to_file(cands, candfilenm):
    with open(candfilenm, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)

        # Header (same columns as before, plus 'w' for jerk-search candidates and
        # 'T', the length of the searched time series straight from its .inf)
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
            "T",
            "numhits"
        ])

        for goodcand in cands:
            fields = str(goodcand).split()
            file, candnum = fields[0].split(':')
            fields[0] = file
            fields.insert(1, candnum)
            fields.append(get_jerk_value(os.path.join(goodcand.path, goodcand.filename), int(candnum)))
            fields.append(goodcand.T)
            fields.append(len(goodcand.hits))
            writer.writerow(fields)

#--------------------------------------------------------------

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='presto accel sifter for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')

    args = parser.parse_args()

    min_num_DMs, low_DM_cutoff = apply_sift_args(args.processing_args)

    # Try to read the .inf files first, as _if_ they are present, all of
    # them should be there.  (if no candidates are found by accelsearch
    # we get no ACCEL files...
    path = f"{args.out_dir}/inj_{args.injection_number:06}/processing/PRESTO/ACCEL"
    out_csv = f"{args.out_dir}/inj_{args.injection_number:06}/processing/PRESTO/PRESTO_candidates.csv"
    inffiles = glob.glob(globinf, root_dir=path)
    # accelsearch writes 'root_ACCEL_<zmax>' (the text list we want) alongside
    # 'root_ACCEL_<zmax>.cand'/'.txtcand'; select by excluding those extensions
    # rather than assuming zmax ends in a '0'.
    candfiles = [f for f in glob.glob(globaccel, root_dir=path)
                 if not f.endswith(('.cand', '.txtcand', '.inf'))]

    if not (inffiles and candfiles):
        # finding nothing is a result, not an error (a non-zero exit would be
        # retried and then end the whole run): write an empty list, which also
        # replaces any stale one from an earlier run, and let the matcher carry on
        print('No PRESTO candidates to sift.')
        to_file([], out_csv)
        sys.exit(0)

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

    # sifting assumes every file it is given searched the same stretch of data:
    # remove_duplicate_candidates merges on Fourier bin r (equal for equal-length
    # segments) and remove_harmonics zaps equal frequencies (its factor 1). Run
    # over all segments at once, only the strongest detection of a pulsar in the
    # whole search survives. So sift each segment on its own - its downsamples and
    # DM trials still merge, as intended - and combine the survivors.
    seg_re = re.compile(r'_SEG_(\d+_\d+)_')
    by_segment = {}
    for candf in candfiles_new:
        seg = seg_re.search(os.path.basename(candf))
        by_segment.setdefault(seg.group(1) if seg else '', []).append(candf)

    goodcands = []
    for segment, seg_files in sorted(by_segment.items()):
        # the same DM appears once per downsample: dedupe so remove_DM_problems
        # counts each DM once
        dms = sorted(set(float(accel_re.sub('', os.path.basename(f)).split("DM")[-1]) for f in seg_files))
        dmstrs = ["%.2f"%x for x in dms]

        cands = sifting.read_candidates(seg_files) # edit RS

        # Remove candidates that are duplicated in other ACCEL files
        if len(cands):
            cands = sifting.remove_duplicate_candidates(cands)

        # Remove candidates with DM problems
        if len(cands):
            cands = sifting.remove_DM_problems(cands, min_num_DMs, dmstrs, low_DM_cutoff)

        # Remove candidates that are harmonically related to each other
        # Note:  this includes only a small set of harmonics. A lone candidate
        # has nothing to be a harmonic of (and remove_harmonics would compare it
        # with itself).
        if len(cands) > 1:
            cands = sifting.remove_harmonics(cands)

        print(f"SEG_{segment}: {len(cands)} candidates after sifting")
        goodcands.extend(cands.cands)

    goodcands.sort(key=attrgetter('sigma'), reverse=True)
    to_file(goodcands, out_csv)


