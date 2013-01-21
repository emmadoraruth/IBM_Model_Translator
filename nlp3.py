import math
import re
import sys
import json
import collections

# IBM Model 1
def ibm1(native, foreign, iterations):
    tparams = init1(native, foreign) # Initalize translation parameters hash
    for i in range (iterations): # Repeat for given number of iterations
        nfile = open(native, 'r')
        ffile = open(foreign, 'r')
        # Counts default to 0
        ncounts = collections.defaultdict(default_factory0) # c(n) counts
        fncounts = collections.defaultdict(default_factory0) # c(f,n) counts
        for nline in nfile: # Parse native file by line
            fline = ffile.readline() # Simultaneously parse foreign file by line
            nparsed = str.split(nline) # Parse native line by space
            fparsed = str.split(fline) # Parse foreign line by space
            for m in range(len(fparsed)): # Iterate through foreign sentence positions
                fword = fparsed[m]
                # Calculate sum used in delta value
                sum = tparams['NULL'][fword]
                for l in range(len(nparsed)):
                    sum += tparams[nparsed[l]][fword]
                for l in range(len(nparsed)+1): # Iterate through native sentence positions
                    if l == 0:
                        nword = 'NULL'
                    else:
                        nword = nparsed[l-1]
                    # Calculate delta value and increment counts
                    delta = tparams[nword][fword]/sum
                    ncounts[nword] += delta
                    fncounts[(fword, nword)] += delta
        nfile.close()
        ffile.close()
        # Update translation parameters
        for n in tparams.keys():
            for f in tparams[n].keys():
                tparams[n][f] = float(fncounts[(f,n)])/float(ncounts[n])
    tfile = open('tparams.txt', 'w')
    tfile.write(json.dumps(tparams))
    tfile.close()
    return tparams # Return translation parameters hash

#IBM Model 2
def ibm2(native, foreign, iterations):
    params = init2(native, foreign)
    tparams = params[0] # Initialize translation parameters hash
    qparams = params[1] # Initialize alignment parameters hash
    for i in range (iterations): # Repeat for given number of iterations
        nfile = open(native, 'r')
        ffile = open(foreign, 'r')
        ncounts = collections.defaultdict(default_factory0) # c(n) counts
        fncounts = collections.defaultdict(default_factory0) # c(f,n) counts
        mcounts = collections.defaultdict(default_factory0) # c(i,l,m) counts
        lmcounts = collections.defaultdict(default_factory0) # c(j|i,l,m) counts
        for nline in nfile: # Parse native file by line
            fline = ffile.readline() # Simultaneously parse foreign file by line
            nparsed = str.split(nline) # Parse native line by space
            fparsed = str.split(fline) # Parse foreign line by space
            nlen = len(nparsed)
            flen = len(fparsed)
            for m in range(flen): # Iterate through foreign sentence positions
                fword = fparsed[m]
                # Calculate sum used in delta value
                sum = qparams[(nlen,flen)][(0,m+1)]*tparams['NULL'][fword]
                for l in range(nlen):
                    sum += qparams[(nlen,flen)][(l+1,m+1)]*tparams[nparsed[l]][fword]
                for l in range(nlen+1): # Iterate through native sentence positions
                    if l == 0:
                        nword = 'NULL'
                    else:
                        nword = nparsed[l-1]
                    # Calculate delta value and increment counts
                    delta = tparams[nword][fword]*qparams[(nlen,flen)][(l,m+1)]/sum
                    ncounts[nword] += delta
                    fncounts[(fword, nword)] += delta
                    mcounts[(m+1,nlen,flen)] += delta
                    lmcounts[(l+1,m+1,nlen,flen)] += delta
        nfile.close()
        ffile.close()
        # Update translation parameters
        for n in tparams.keys():
            for f in tparams[n].keys():
                tparams[n][f] = float(fncounts[(f,n)])/float(ncounts[n])
        # Update alignment parameters
        for k in (qparams[(nlen,flen)]).keys():
            qparams[(nlen,flen)][k] = lmcounts[(k[0],k[1],nlen,flen)]/mcounts[(k[1],nlen,flen)]
        qcopy = collections.defaultdict(default_factory)
        for key1 in qparams.keys():
            for key2 in (qparams[key1]).keys():
                qcopy[str(key1)][str(key2)] = qparams[key1][key2]
        qfile = open('qparams.txt', 'w')
        qfile.write(json.dumps(qcopy))
        qfile.close()
        tfile = open('tparams.txt', 'w')
        tfile.write(json.dumps(tparams))
        tfile.close()
        return (tparams, qparams) # Return translation and alignment parameter hashes

# Helper function to initialize dictionary of empty dictionaries
def default_factory():
    return {}

# Helper function to initialize dictionary of with default values of 0
def default_factory0():
    return 0

# Initialize t parameters
def init1(native, foreign):
    tparams = collections.defaultdict(default_factory)
    nfile = open(native, 'r')
    ffile = open(foreign, 'r')
    fwords = collections.defaultdict(default_factory0)
    for nline in nfile:
        fline = ffile.readline()
        nparsed = str.split(nline)
        fparsed = str.split(fline)
        for n in nparsed:
            for f in fparsed:
                fwords[f] = 1
                tparams[n][f] = 0  # Add all native/foregin word pairs that appear in corresponding sentences to possible translations
                tparams['NULL'][f] = 0 # Add alignment of all foreign words to 'NULL' word
    for n in tparams.keys(): # Iterate over all native words
        num = len(tparams[n]) # Number of possible translations
        for f in (tparams[n]).keys(): # Iterate over all possible translations
            tparams[n][f] = 1.0/num # Initialize translation parameters
    nfile.close()
    ffile.close()
    fwordsfile = open('foreign_words.txt', 'w')
    fwordsfile.write(json.dumps(fwords))
    fwordsfile.close()
    return tparams # Return translation parameters

# Initialize t and q parameters
def init2(native, foreign):
    tparams = ibm1(native, foreign, 10) # Initialize translation parameters using ten iterations of IBM Model 1
    qparams = collections.defaultdict(default_factory)
    nfile = open(native, 'r')
    ffile = open(foreign, 'r')
    for nline in nfile:
        fline = ffile.readline()
        nparsed = str.split(nline)
        fparsed = str.split(fline)
        nlen = len(nparsed)
        flen = len(fparsed)
        num = nlen + 1 # Number of possible alignments
        for l in range(nlen+1): # Iterate over all native sentence positions
            for m in range(flen): # Iterate over all foreign sentence positions
                qparams[(nlen,flen)][(l,m+1)] = 1.0/num # Initialize alignment parameters
    nfile.close()
    ffile.close()
    return (tparams, qparams) # Return alignment and translation parameters

# Find most likely foreign translations for each native word in file
def translations(file, num_options):
    tparams = parameters(1)[0]
    nfile = open(file, 'r')
    ffile = open('translations.txt', 'w')
    for line in nfile:
        n = str.split(line)[0] # Iterate over words in the file
        fwords = sorted(tparams[n], key=(tparams[n]).get, reverse=True) # Get all possible translations, sorted by probability descending
        fprobs = sorted((tparams[n]).values(), reverse=True) # Get corresponding probabilities
        ffile.write(n + ':\n')
        for i in range(num_options): # Write specified number of translations and corresponding probabilities to output file
            ffile.write(str(fwords[i]) + ' ' + str(fprobs[i]) + '\n')
        ffile.write('\n')
    nfile.close()
    ffile.close()

# Find the most likely alignments for corresponding pairs of sentences in the native language and foreign translation files
def alignments(native, foreign, model, num_sentences=0): # 0 default value means align all (not used in this assignment)
    temp = parameters(model)
    params = temp[0]
    string = temp[1]
    fwordsfile = open('foreign_words.txt', 'r')
    fwords = json.loads(fwordsfile.readline())
    fwordsfile.close()
    nfile = open(native, 'r')
    ffile = open(foreign, 'r')
    if model == 1 or model == str(1):
        afile = open('alignments1.txt', 'w')
    else:
        afile = open('alignments2.txt', 'w')
    i = 0 # Sentence counter
    for nline in nfile: # Iterate over native sentences
        fline = ffile.readline() # Simultaneously iterate over corresonding foreign sentences
        a = (alignment_helper(nline, fline, params, fwords, model, string))[0] # Find maximum likelihood alignment for sentence pair
        afile.write(str(nline) + str(fline) + str(a) + '\n\n') # Write alignment to output file
        i += 1
        if num_sentences != 0 and i >= num_sentences: # Only align specified number of sentences
            break
    afile.close()
    nfile.close()
    ffile.close()

# Find most likely alignment for given pair of native and foreign sentences
def alignment_helper(nline, fline, params, foreign_words, model, string):
    a = []
    nparsed = str.split(nline)
    fparsed = str.split(fline)
    flen = len(fparsed)
    nlen = len(nparsed)
    prob = 1.0
    for f in range(flen): # Iterate over all foreign sentence positions
        max_prob = 0.0
        max = 0
        if model == 1 or model == str(1):
            tparams = params
            if fparsed[f] in tparams['NULL']:
                max_prob = float(tparams['NULL'][fparsed[f]])
        else:
            tparams = params[0]
            qparams = params[1]
            lens = (nlen, flen)
            inds = (0, f+1)
            if string == 1:
                lens = str(lens)
                inds = str(inds)
            if fparsed[f] in tparams['NULL'] and lens in qparams:
                max_prob = float(tparams['NULL'][fparsed[f]] * qparams[lens][inds])
        for n in range(nlen): # Iterate over all native sentence positions and find maximum probability alignment
            curr = 0.0
            if nparsed[n] in tparams and fparsed[f] in foreign_words:
                if fparsed[f] in tparams[nparsed[n]]:
                    curr = tparams[nparsed[n]][fparsed[f]]
            elif nparsed[n] not in tparams and fparsed[f] not in foreign_words:
                curr = float(10**-18)
            else:
                curr = float(10**-19)
            if model == 2 or model == str(2):
                lens = (nlen, flen)
                inds = (n+1, f+1)
                if string == 1:
                    lens = str(lens)
                    inds = str(inds)
                if lens in qparams:
                    curr *= qparams[lens][inds]
                else:
                    curr *= 0.0
            if curr > max_prob:
                max = n+1
                max_prob = float(curr)
        prob *= max_prob # Update probability of maximum likelihood alignment
        a.append(max) # Update alignment
    return (a, prob) # Return maximum likelihood alignment and its probability

# Get parameters from files or generate if files do not exist
def parameters(model):
    string = 1
    tfile = open('tparams.txt', 'r')
    tparams = json.loads(tfile.readline())
    tfile.close()
    params = tparams
    if model == 2 or model == str(2):
        qfile = open('qparams.txt', 'r')
        qparams = json.loads(qfile.readline())
        qfile.close()
        params = (tparams, qparams)
        if tparams is None or qparams is None:
            params = ibm2('corpus.en', 'corpus.de', 5)
            string = 0
    elif tparams is None:
        string = 0
        params = ibm1('corpus.en', 'corpus.de', 5)
    return (params, string)

# Unscramble native alignments
def unscramble(native, foreign):
    temp = parameters(2)
    params = temp[0]
    string = temp[1]
    ffile = open(foreign, 'r')
    afile = open('unscrambled.en', 'w')
    fwordsfile = open('foreign_words.txt', 'r')
    fwords = json.loads(fwordsfile.readline())
    fwordsfile.close()
    for fline in ffile: # Iterate over all foreign sentences
        nfile = open(native, 'r')
        max = -1.0
        translation = ''
        for nline in nfile: # Iterate over all native sentences and find maximum likelihood translation
            nparsed = str.split(nline)
            fparsed = str.split(fline)
            flen = len(fparsed)
            nlen = len(nparsed)
            curr = (alignment_helper(nline, fline, params, fwords, 2, string))[1] # Find probability of maximum likelihood alignment of sentence pair
            if curr > max:
                max = curr
                translation = nline
        afile.write(str(translation)) # Write maximum likelihood translation to output file
        nfile.close()
    afile.close()
    ffile.close()

def main(args):
    if args[1] == 'ibm1':
        tparams = ibm1(args[2], args[3], 5)
    elif args[1] == 'translations':
        translations(args[2], 10)
    elif args[1] == 'alignments':
        alignments(args[2], args[3], args[4], 20)
    elif args[1] == 'ibm2':
        params = ibm2(args[2], args[3], 5)
    elif args[1] == 'unscramble':
        unscramble(args[2], args[3])

if __name__ == "__main__":
    main(sys.argv)