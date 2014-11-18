#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; -*-
# A work-in-progress attempt at a "computational glue library"
# Jürgen Jänes <jj374@cam.ac.uk>

import hashlib
import cPickle
import pprint
import tempfile
import numpy as np
import atexit
import commands
import datetime
import math
import os
import re
import pickle
import shutil
import subprocess
import sys
import inspect
import __main__

#import numpy as np
#import scipy as sp

def ts(datetime=datetime.datetime.now()):
    return(datetime.strftime("%y%m%d-%H%M%S"))

def ts_res():
    return ts(startTime)

def log(msg, time=datetime.datetime.now(), v=True):
    ts = abs(startTime - datetime.datetime.now())
    if sys.version_info < (2, 6):
        print "%s %s\n" % (ts, msg)        
        return
    if v:
        line_data = inspect.getframeinfo(inspect.stack()[1][0])

        fn = os.path.split(line_data.filename)[1]
        ln = line_data.lineno
        # multi-line strings
        if (isinstance(msg, str) and ("\n" in msg)):
            print "%s %s:%d\n%s" % (ts, fn, ln, msg)
        # pretty-print dictionaries
        elif (isinstance (msg, dict)): 
            print "%s %s:%d\n%s" % (ts, fn, ln, pprint.pformat(msg))
        else:
            #cols = int(os.popen('stty size', 'r').read().split()[1])
            cols=80
            st_l = "%s %s" % (ts, msg)
            st_r = "(%s:%d)" % (fn, ln)
            if (len(st_l) + len(st_r)) < cols:
                print st_l + st_r.rjust(cols - len(st_l))
            else:
                st_l = "%s" % (ts,)
                print st_l + st_r.rjust(cols - len(st_l)) + ("\n%s" % (msg,))

def logList(msg, val):
    log("%s: %s" % (msg, val[:4]))

def start():
    log("Starting.\nFile:\t%s\nTime:\t%s" % (__main__.__file__, ts(startTime)))

# http://jimmyg.org/blog/2009/working-with-python-subprocess.html#without-the-shell
def whereis(program):
    for path in os.environ.get('PATH', '').split(':'):
        if os.path.exists(os.path.join(path, program)) and not os.path.isdir(os.path.join(path, program)):
            return os.path.join(path, program)
    return None

def run(cmd='ls -l', abort=True, work_dir=None, extraPaths=None):
    if not (work_dir is None):
        prev_dir = os.getcwd()
        os.chdir(work_dir)
    if (isinstance(cmd, list)):
        cmd = " ".join(cmd)
    log("Executing: %s" % (cmd,))
    run_start = datetime.datetime.now()
    if extraPaths is None:
        ret_val = subprocess.call(cmd, stderr=subprocess.STDOUT, shell=True)
    else:
        ret_val = subprocess.call(cmd, stderr=subprocess.STDOUT, shell=True, env = {"PATH": ':'.join(extraPaths + [os.getenv('PATH')])})
    log("Finished executing. Return value: %s, Wall clock time: %s" % 
        (ret_val, abs(run_start - datetime.datetime.now()),))
    if not (work_dir is None):
        os.chdir(prev_dir)
    if (abort and (ret_val != 0)):
        log("Non-zero return value. Aborting")
        sys.exit(1)
    return ret_val

def R(rscript="a <- 3; a"):
    """ Run a chunk of R code from python. Useful for e.g. plotting stuff that matplotlib doesn't
      do well.
    """
    log("Executing R code: %s" % (rscript,))
    cmd = 'R --vanilla --slave'
    R_process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    R_process.communicate(input=rscript)

# calculates the mean and uncertainty from a list of measurements
# assumes standard distribution and 95% confidence intervals
#def sprintu(values, format):
#    avg = sp.mean(values)
#    err = 1.96 * sp.std(values) / math.sqrt(len(values)) 
#    fmt = "%s±%s" % (format, format)
#    return fmt % (avg, err)
# calculate row means
#def means(a):
#    m = []
#    for i in range(len(a)):
#        m.append(sp.mean(a[i]))
#    return m

# calculate row errors
#def errors(a):
#    m = []
#    for i in range(len(a)):
#        m.append(1.96 * sp.std(a[i]) / math.sqrt(len(a[i])))
#    return m

# create a directory if it does not exist
def mkresdir(exp_name):
    # TODO add microsecond-resolutions?
    timestamp = datetime.datetime.now().strftime("%y%m%d-%H%M%S")
    dirname = "%s-%s" % (exp_name, timestamp)
    if not os.path.isdir(dirname):
        os.mkdir( dirname )
    return (dirname)
    #print mkresdir("compare_models")
def mkresfile(file_name, exp_name):
    return os.path.join(exp_name, file_name)

#def save_variable(val, exp_name="", file_name):
def save_variable(val, file_name, exp_name="", v=True):
    #cr.save_variable((self.x1, self.y1, self.x2, self.y2), dout, "%d-%d-coordinates" % (iSlice, iROI), v)
    outf = open(mkresfile('%s.pkl' % (file_name), exp_name), 'wb')
    pickle.dump(val, outf)
    outf.close()
    log("save_variable: wrote to %s" % (outf,), v)

def load_variable(file_name, exp_name, v=True):
    #(self.x1, self.y1, self.x2, self.y2) = cr.load_variable(din, "%d-%d-coordinates" % (iSlice, iROI), v)
    inpf = open(mkresfile('%s.pkl' % (file_name), exp_name), 'rb')
    val = pickle.load(inpf)
    if v: log("load_variable: read from %s" % (inpf,), v)
    inpf.close()
    return val

class ResultsWriter:
    def __init__(self, exp_name=None, v=True):
        # use the name of the main script file as the name of the output directory
        if (exp_name is None):
            exp_name = os.path.splitext(__main__.__file__)[0]

        self.exp_name = exp_name
        self.exp_fname = mkresdir(exp_name)

        # ugly hack to duplicate screen output to a file
        #sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
        #tee = subprocess.Popen(["tee", "%s/script-output.txt" % (self.exp_fname)], stdin=subprocess.PIPE)
        #os.dup2(tee.stdin.fileno(), sys.stdout.fileno())
        #os.dup2(tee.stdin.fileno(), sys.stderr.fileno())

        # make a copy of the s    
        log("Making copy of the main script file\nFrom:\t%s\nTo:\t%s" % (__main__.__file__, self.exp_fname))
        shutil.copy(__main__.__file__, self.exp_fname)

    def fp(self, file_name):
        return mkresfile(file_name, self.exp_fname)
    
    def save_variable(self, val, file_name, v=True):
        save_variable(val, file_name, self.exp_fname, v)

    def cdto_expdir(self, subDir = None):
        self.prevdir = os.getcwd()
        newdir = self.exp_fname if (subDir is None) else os.path.join(self.prevdir, self.exp_fname, subDir)
        if not os.path.isdir(newdir):
            os.mkdir(newdir)
        os.chdir(newdir)
        print "ic: now in directory: %s" % (newdir,)

    def cdto_prevdir(self):
        os.chdir(self.prevdir)
        print "ic: now in directory: %s" % (self.prevdir,)

class ResultsReader(ResultsWriter):
    def __init__(self, exp_name=None):
        #self.exp_name = os.path.basename( os.path.splitext(__main__.__file__)[0] )
        #if not (exp_fname is None):
        #    self.exp_fname = exp_fname
        #else:
        self.exp_name = exp_name
        results = filter(os.path.isdir, os.listdir("./"))
        results.sort()
        p = re.compile("^%s" % self.exp_name)
        results = filter(p.match, results)
        self.exp_fname = results[-1]
        #print "setting exp_fname to: %s" % (self.exp_fname)
        log("loading results from: %s" % (self.exp_fname))

    def load_variable(self, file_name, v=True):
        return load_variable(file_name, self.exp_fname, v)
    
    def load_d_fn(self, file_name, v=True):
        #load file name dictionary
        d_fn = self.load_variable(file_name, v)
        #if v: log(d_fn)
        for (par, fn) in d_fn.iteritems():
            d_fn[par] = os.path.join(self.exp_fname, fn)
        if v: log(d_fn)
        return(d_fn)

class NaivelyCachedComputation:
    def __init__(self, recompute=False, resume=False, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.path = self.makePath()
        if ((recompute) and os.path.exists(self.path)): shutil.rmtree(self.path)

        # Remove non-finished computation
        if (os.path.exists(self.path) and not os.path.exists(os.path.join(self.path, "TP-Finished"))): 
            log("Clearing incomplete calculation: %s" % (self.path,))
            shutil.rmtree(self.path)

        # If start fresh, dir e.g. TranscriptPool2_flSubset=True_seedVal=1.res.tp not exit
        if not(os.path.exists(self.path)):
            log("New computation: %s" % (self.path,))
            os.mkdir(self.path)
            self.compute(*self.args, **self.kwargs)
            self.save()
            log("Creating TP-finished file.")
            run("touch %s" %(os.path.join(self.path, "TP-Finished"),))
        else:
            log("Cached computation: %s" % (self.path,))
            self.load()

    def makePath(self): # this can be overloaded(?considering the hashing we're doing) in subclasses for specific cases where arguments become more complex
        # argStr is a string representation of all the arguments
        #argStr = ",".join(map(str, self.args) + map(lambda (key, val): "%s=%s" % (key, val), self.kwargs.iteritems()))
        argList = []
        for arg in self.args:
            if isinstance(arg,NaivelyCachedComputation):
                argList.append("%s" % (arg.path,))
            else:
                argList.append("%s" % (arg,))
        for (key, val) in sorted(self.kwargs.iteritems()):
            if isinstance(val,NaivelyCachedComputation):
                argList.append("%s=%s" % (key,val.path))
            else:
                argList.append("%s=%s" % (key,val))

        m = hashlib.md5()
        for arg in argList: m.update(arg)
        argHash = m.hexdigest()

        hashThreshold = 60
        argStr = ""
        for arg in argList:
            if (len(argStr) + len(arg) > hashThreshold): 
                argStr += "_" + argHash
                break
            if (arg.find('/') == -1): #path name validity check
                argStr += "_" + arg
        #argStr = argStr[1:]

        # path is the path where all the data will be stored
        return("%s%s.res.tp" % (self.__class__.__name__, argStr))

    def compute(self, *args, **kwargs):
        raise NotImplementedError("NaivelyCachedComputation.compute() not implemented.")

    def fp(self, fn): return(os.path.join(self.path, fn))

    def fpPkl(self): return(self.fp("%s.pkl" % (self.path,)))

    # load() and save() from: http://stackoverflow.com/questions/2709800/how-to-pickle-yourself
    def load(self):
        fh = open(self.fpPkl(), 'rb')
        tmp_dict = cPickle.load(fh)
        fh.close()
        self.__dict__.update(tmp_dict) 

    def save(self):
        fh = open(self.fpPkl(), 'wb')
        cPickle.dump(self.__dict__, fh, 2)
        fh.close()

    # TODO class method .rm() removes all results for that specific type
    # TODO class method .iter() iterates over all calls it can find over the arguments provided

@atexit.register
def fin():
    log("Aloha!")#. Total running time: %s" % (abs(startTime - datetime.datetime.now())))
    # Total running time is already visible from the timestamp
    # sys.exit(0)
    # Explicitly trying to exit here occasionally causes this:
    # Exception KeyError: KeyError(140735252056416,) in <module 'threading' from '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/threading.pyc'> ignored

startTime = datetime.datetime.now()
outFname = "%s-%s.log.tp" % (os.path.splitext(__main__.__file__)[0], ts(startTime))

## emulate 'tee' for all output
# http://stackoverflow.com/questions/616645/how-do-i-duplicate-sys-stdout-to-a-log-file-in-python
#sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
#tee = subprocess.Popen(["tee", outFname], stdin=subprocess.PIPE)
#os.dup2(tee.stdin.fileno(), sys.stdout.fileno())
#os.dup2(tee.stdin.fileno(), sys.stderr.fileno())

start()

if __name__ == '__main__':
    ## emulate 'tee' for all output using ipython
    #ip = get_ipython()
    #ip.magic_logstart("-o %s backup" % (outFname,))
    run("This command will fail, but the script will continue.", abort=False)
    run("echo 'Hello, world'")
    run("sleep 0.2")
    #exp = mkresdir("store_variables")
    a = [1, 2, 3]
    b = "Aloha!"
    #save_variable((a, b), "a_and_b", exp)
    #(a_new, b_new) = load_variable("a_and_b", exp)
    #print a_new
    #print b_new
    #rw_transcript_pool = ResultsWriter()
    #rw_transcript_pool.save_variable((a,b), "a_and_b_via_rw")
    log("Aloha!")
    #import IPython as ip
    #ip.embed()