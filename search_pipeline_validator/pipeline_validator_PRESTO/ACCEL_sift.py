from __future__ import absolute_import
from builtins import map
import re
import glob
import presto.sifting as sifting
from operator import itemgetter, attrgetter

# Note:  You will almost certainly want to adjust
#        the following variables for your particular search

# glob for ACCEL files
# globaccel = "*ACCEL_*0"
globaccel = "*ACCEL_200"
# glob for .inf files
globinf = "*DM*.inf"
# In how many DMs must a candidate be detected to be considered "good"
min_num_DMs = 2
# Lowest DM to consider as a "real" pulsar
low_DM_cutoff = 2.0
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

# Try to read the .inf files first, as _if_ they are present, all of
# them should be there.  (if no candidates are found by accelsearch
# we get no ACCEL files...
path = "/tmp/rsenzel/input_data" # edit RS
inffiles = glob.glob(globinf, root_dir=path) # edit RS
candfiles = glob.glob(globaccel, root_dir=path) # edit RS
# Check to see if this is from a short search
if len(re.findall("_[0-9][0-9][0-9]M_" , inffiles[0])):
    dmstrs = [x.split("DM")[-1].split("_")[0] for x in candfiles]
else:
    dmstrs = [x.split("DM")[-1].split(".inf")[0] for x in inffiles]
dms = list(map(float, dmstrs))
dms.sort()
dmstrs = ["%.2f"%x for x in dms]

# Read in all the candidates
candfiles_new = [path+"/"+candf for candf in candfiles] # edit RS
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

def to_file(self, candfilenm=path+"/cands1.txt"): # edit RS
    candfile = open(candfilenm, "w")
    candfile.write("#" + "file:candnum".center(66) + "DM".center(9) +
                    "SNR".center(8) + "sigma".center(8) + "numharm".center(9) +
                    "ipow".center(9) + "cpow".center(9) +  "P(ms)".center(14) +
                    "r".center(12) + "z".center(8) + "numhits".center(9) + "\n")
    for goodcand in self.cands:
        candfile.write("%s (%d)\n" % (str(goodcand), len(goodcand.hits)))
        if (len(goodcand.hits) > 1):
            goodcand.hits.sort(key=lambda cand: float(cand[0]))
            for hit in goodcand.hits:
                numstars = int(hit[2]/3.0)
                candfile.write("  DM=%6.2f SNR=%5.2f Sigma=%5.2f   "%hit + \
                                numstars*'*' + '\n')
    if candfilenm is not None:
        candfile.close()

# Write candidates to STDOUT
if len(cands):
    cands.sort(key=attrgetter('sigma'), reverse=True)
    to_file(cands) # edit RS
