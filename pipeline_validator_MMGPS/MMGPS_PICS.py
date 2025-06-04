import os
import sys
sys.path.append('/home/psr')

import glob
import cPickle
import optparse

from ubc_AI.data import pfdreader


def get_id_from_cand_file(filename):
    return(filename.split('_')[-1].strip(".ar"))


def extract_and_score(opts):
    path = opts.in_path
    file_type = opts.file_type
    AI_PATH = opts.model_dir
    models = os.listdir(AI_PATH)

    classifiers = []
    for model in models:
        with open(os.path.join(AI_PATH, model), "rb") as f:
            classifiers.append(cPickle.load(f))

    arfiles = sorted(glob.glob("{}/*.{}".format(path,file_type)),key=get_id_from_cand_file)

    scores = []
    readers = [pfdreader(f) for f in arfiles]
    for classifier in classifiers:
        scores.append(classifier.report_score(readers))

    combined  = zip(arfiles, *scores)
    names = ",".join(["{}".format(model.split("/")[-1]) for model in models])
    with open("{}/pics_scores.txt".format(path),"w") as fout:
        fout.write("arfile,{}\n".format(names))
        for row in combined:
            scores = ",".join(map(str,row[1:]))
            fout.write("{},{}\n".format(row[0], scores))


if __name__ == '__main__':

    parser = optparse.OptionParser()
    parser.add_option('--in_path',type=str,help='input path for files to be scored',dest="in_path")
    parser.add_option('--file_type',type=str,help='Type of file (ar/pfd/ar2/clfd)',dest="file_type",default="ar")
    parser.add_option('--model_dir',type=str,help='directory of PICS models')
    parser.add_option('--work_dir',type=str,help='work directory')
    opts,args = parser.parse_args()

    os.environ['THEANO_FLAGS'] = 'base_compiledir={}'.format(opts.work_dir)

    extract_and_score(opts)
