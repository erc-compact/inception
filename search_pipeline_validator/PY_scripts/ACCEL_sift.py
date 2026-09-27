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

globaccel = "*ACCEL_*"
globinf = "*DM*.inf"

SIFT_DEFAULTS = {
    "sigma_threshold": 4.0,
    "c_pow_threshold": 100.0,
    "harm_pow_cutoff": 8.0,
    "r_err": 1.1,
    "short_period": 0.0005,
    "long_period": 15.0,
    "min_num_DMs": 1,
    "low_DM_cutoff": 1.0,
    "known_birds_p": [],
    "known_birds_f": [],
}


def apply_sift_args(processing_args):
    args = inj_tools.parse_JSON(processing_args).get('presto_sift_args', {})
    sift = dict(SIFT_DEFAULTS, **args)

    for key in ("sigma_threshold", "c_pow_threshold", "harm_pow_cutoff",
                "r_err", "short_period", "long_period"):
        setattr(sifting, key, sift[key])

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

        # Header (same columns as before)
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
    path = f"{args.out_dir}/inj_{args.injection_number:06}/processing/PRESTO/SEARCH"
    out_csv = f"{args.out_dir}/inj_{args.injection_number:06}/processing/PRESTO/PRESTO_candidates.csv"
    inffiles = glob.glob(globinf, root_dir=path)
    candfiles = [f for f in glob.glob(globaccel, root_dir=path)
                 if not f.endswith(('.cand', '.txtcand', '.inf'))]

    if not (inffiles and candfiles):
        print('No PRESTO candidates to sift.')
        to_file([], out_csv)
        sys.exit(0)

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

    seg_re = re.compile(r'_SEG_(\d+_\d+)_')
    by_segment = {}
    for candf in candfiles_new:
        seg = seg_re.search(os.path.basename(candf))
        by_segment.setdefault(seg.group(1) if seg else '', []).append(candf)

    goodcands = []
    for segment, seg_files in sorted(by_segment.items()):
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
        # Note:  this includes only a small set of harmonics
        if len(cands) > 1:
            cands = sifting.remove_harmonics(cands)

        print(f"SEG_{segment}: {len(cands)} candidates after sifting")
        goodcands.extend(cands.cands)

    goodcands.sort(key=attrgetter('sigma'), reverse=True)
    to_file(goodcands, out_csv)


