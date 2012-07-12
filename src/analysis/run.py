#!/usr/bin/env python

import biotools.analysis.predict as genepredict
import biotools.analysis.options as options
import threading
from os import sep
try:
    import Queue as queue
except ImportError:
    import queue


def run(infile, strains):
    '''
    Run several instances of `genepredict.run` at once.
    '''

    q = queue.Queue()
    filenames = []

    def run_predict():
        while 1:
            try:
                strainf = q.get(False)
            except queue.Empty:
                break
            strain = strainf.split(sep)[-1]
            pos = strain.rfind('.')
            if pos > 1 or (pos == 1 and strain[0] != '.'):
                strain = strain[:pos]

            options.debug("Predicting for %s." % strain)

            try:
                genepredict.run(infile, strainf, strain, filenames)
            except RuntimeError:
                pass
            q.task_done()

    for strain in strains:
        q.put(strain)

    for i in range(options.NUM_PROCESSES - 1):
        curr = threading.Thread(target=run_predict)
        curr.start()

    run_predict()
    q.join()

    return filenames

if __name__ == "__main__":
    import sys
    try:
        run(sys.argv[1:])
    except IndexError:
        pass
