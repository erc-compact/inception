import os
import sys
sys.path.append('/home/psr')

import glob
import cPickle
import argparse

from ubc_AI.data import pfdreader


def extract_and_score(args):

    path = args.cand_dir
    file_type = args.extension
    AI_PATH = args.model_dir
    models = os.listdir(AI_PATH)

    def get_id_from_cand_file(filename):
        return(filename.split('_')[-1].strip(".{}".format(file_type)))

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

    with open(args.output_file, "w") as fout:
        fout.write("file,{}\n".format(names))
        for row in combined:
            scores = ",".join(map(str,row[1:]))
            fout.write("{},{}\n".format(row[0], scores))


if __name__ == '__main__':

    parser = argparse.ArgumentParser('PICS classifier')
    parser.add_argument('--cand_dir', type=str, help='directory of files to be scored')
    parser.add_argument('--model_dir', type=str, help='directory of PICS models')
    parser.add_argument('--extension',type=str, help='Type of file (ar/pfd/ar2/clfd)')
    parser.add_argument('--output_file', type=str, help='output results file')
    args = parser.parse_args()

    os.environ['THEANO_FLAGS'] = 'base_compiledir={}'.format(os.getcwd())
    extract_and_score(args)
