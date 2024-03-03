#!/Users/russell/miniconda3/envs/azimuth/bin/python

import argparse
import csv
import numpy
import os
import re
import subprocess
from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import Seq
from Bio.Restriction import BsaI
from itertools import islice
from time import time
import azimuth.model_comparison
import pickle


# minimum off-target score of standard off-targets (those that end with NGG)
# This should probably be based on the CFD score these days
# But for now, I'll let the user do the filtering
MINSCORE = 0.0

# minimum off-target score for alternative PAM off-targets
# There is not a lot of data to support this cutoff, but it seems
# reasonable to have at least some cutoff, as otherwise we would show
# NAG and NGA like NGG and the data shows clearly that the alternative
# PAMs are not recognized as well as the main NGG PAM.
# so for now, I just filter out very degenerative ones. the best solution
# would be to have a special penalty on the CFD score, but CFS does not
# support non-NGG PAMs (is this actually true?)
ALTPAMMINSCORE = 1.0

designs = (((BsaI,), 'TGTGGAAAGGACGAAACACCG{}GTTTAAGAGCTATGCTGGAAGGAGA'),)


def design_test(n):
    def _test(seq):
        as_design = Seq(design.format(seq), IUPACAmbiguousDNA)
        is_cut = any(enzyme.search(as_design) for enzyme in enzymes)
        return not (pattern.search(seq) or is_cut)
    pattern = re.compile('(^[Gg]{5})|([Tt]{4})')
    enzymes, design = designs[n]
    return _test

def calcMitGuideScore(hitSum):
    """ Sguide defined on http://crispr.mit.edu/about 
    Input is the sum of all off-target hit scores. Returns the specificity of the guide.
    """
    score = 100 / (100+hitSum)
    score = int(round(score*100))
    return score

hitScoreM = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]


def fz_score(guideSeq, otSeqs, maxMm=None, minHitScore=None, minAltHitScore=None):
    """
    adapted from https://snippets.siftie.com/harijay/fz-score-d1324dab/
    linked in CRISPOR paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4934014/
    """
    # n = -len(b)
    # # The Patrick Hsu weighting scheme
    # M = [
    #     0.000, 0.000, 0.014, 0.000, 0.000, 0.395, 0.317, 0.000, 0.389, 0.079,
    #     0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]
    # if n > -19:
    #     return 0
    # I = [i for i, (j, k) in islice(enumerate(zip(a[n:], b)), 20)  if j != k]
    # dists = [i - j for i, j in zip(I[1:], I[:-1])]
    # denom = 1. if not dists else (1 - sum(dists) / len(dists) / 19) * 4 + 1
    # score = 1. if not I else 1. / len(I) ** 2
    # for i in I:
    #     score *= 1 - M[i + n]
    # return score / denom
    if maxMm:
        newOtSeqs = []
        for otSeq in otSeqs:
            mmCount, diffLogo = countMms(otSeq, guideSeq)
            if mmCount <= maxMm:
                newOtSeqs.append(otSeq)
        otSeqs = newOtSeqs
    altOts = []
    mainOts = []
    for ot in otSeqs:
        if ot.endswith("AG") or ot.endswith("GA"):
            altOts.append(ot)
        if ot.endswith("GG"):
            mainOts.append(ot)

    # calc and filter hit scores
    mainHitScores = [calcHitScore(guideSeq, ot) for ot in mainOts]
    mainHitScores = [h for h in mainHitScores if h > minHitScore]
    altHitScores  = [calcHitScore(guideSeq, ot) for ot in altOts]
    altHitScores  = [h for h in altHitScores if h > minAltHitScore]

    mainHitScores.extend(altHitScores)
    scoreSum = sum(mainHitScores[1:])
    score = calcMitGuideScore(scoreSum)
    return score

def calcHitScore(string1,string2, startPos=0):
    """ 
    The MIT off-target score
    see 'Scores of single hits' on http://crispr.mit.edu/about
    startPos can be used to feed sequences longer than 20bp into this function
    the most likely off-targets have a score of 100
    >>> int(calcHitScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
    100
    >>> int(calcHitScore("AGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
    100
    >>> int(calcHitScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGA"))
    41
    """
    # The Patrick Hsu weighting scheme
    #print(string1)
    #print(string2)
    
    if len(string1)==len(string2)==23:
        string1 = string1[:20]
        string2 = string2[:20]
    
    

    dists = [] # distances between mismatches, for part 2
    mmCount = 0 # number of mismatches, for part 3
    lastMmPos = None # position of last mismatch, used to calculate distance

    score1 = 1.0
    
    for pos in range(min(20, len(string2))):
        if (string1[pos]!=string2[pos]):
            mmCount+=1
            if lastMmPos!=None:
                dists.append(pos-lastMmPos)
            score1 *= 1-hitScoreM[pos]
            lastMmPos = pos
    # 2nd part of the score
    if mmCount<2: # special case, not shown in the paper
        score2 = 1.0
    else:
        avgDist = sum(dists)/len(dists)
        score2 = 1.0 / (((19-avgDist)/19.0) * 4 + 1)
    # 3rd part of the score
    if mmCount==0: # special case, not shown in the paper
        score3 = 1.0
    else:
        score3 = 1.0 / (mmCount**2)

    score = score1 * score2 * score3 * 100
    return score

############################################################################
#CFD SCORE CALCULATION
#
############################################################################

#Reverse complements a given string
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

#Unpickle mismatch scores and PAM scores
pipline = os.path.dirname(os.path.abspath(__file__))
def get_mm_pam_scores():
    try:
        mm_scores = pickle.load(open(pipline+'/mismatch_score.pkl','rb'))
        pam_scores = pickle.load(open(pipline+'/pam_scores.pkl','rb'))
        return (mm_scores,pam_scores)
    except: 
        raise Exception("Could not find file with mismatch scores or PAM scores")

#Calculates CFD score

def calc_cfd(wt,sg,pam):
    #mm_scores,pam_scores = get_mm_pam_scores()
    score = 1
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    s_list = list(sg)
    wt_list = list(wt)
    for i,sl in enumerate(s_list):
        if wt_list[i] == sl:

            score*=1
        else:
            key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
            score*= mm_scores[key]
    if len(pam) == 2:
        score*=pam_scores[pam]
    else:
        score*= 0
    return (score)

mm_scores, pam_scores = None, None


def calcCfdScore(guideSeq, otSeq):
    """ based on source code provided by John Doench
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGAAAGGG")
    0.4635989007074176
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGGGGGGG")
    1.0
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "aaaaGaGaGGGGGGGGGGGGGGG")
    0.5140384614450001
    # mismatches:      *               !!
    >>> calcCfdScore("ATGGTCGGACTCCCTGCCAGAGG", "ATGGTGGGACTCCCTGCCAGAGG")
    0.5
    # mismatches:    *  ** *
    >>> calcCfdScore("ATGGTCGGACTCCCTGCCAGAGG", "ATGATCCAAATCCCTGCCAGAGG")
    0.53625000020625
    >>> calcCfdScore("ATGTGGAGATTGCCACCTACCGG", "ATCTGGAGATTGCCACCTACAGG")
    """
    global mm_scores, pam_scores
    if mm_scores is None:
        mm_scores,pam_scores = get_mm_pam_scores()
    wt = guideSeq.upper()
    off = otSeq.upper()
    m_wt = re.search('[^ATCG]',wt)
    m_off = re.search('[^ATCG]',off)
    if (m_wt is None) and (m_off is None):
        pam = off[-2:]
        sg = off[:20]
        cfd_score = calc_cfd(wt,sg,pam)
        return cfd_score
# ==== END CFD score source provided by John Doench



def get_pams(seq):
    for match in re.finditer(r'(?=([ACGT]{25}GG[ACGT]{3}))', seq, re.IGNORECASE):
        #yield match.group(1)
        yield match.start(), match.group(1), 'sense', match.start()+22
    for match in re.finditer(r'(?=([ACGT]{3}CC[ACGT]{25}))', seq, re.IGNORECASE):
        yield match.end(), str(Seq(match.group(1)).reverse_complement()), 'antisense', match.end() + 10
    

def parse_bwa(refs, seqs, exclude, nontargeting):
    #print('# Running bwa...')
    run_bwa(refs, seqs)
    #print('# Reading genome...')
    with open(refs, 'rt') as handle:
        genome = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
    #print('# Reading bwa .sam file...')
    ot = re.compile(r'([^:;]+),(\S)(\d+),\w+,\d+')
    xa = re.compile(r'\s+XA:\S+')
    with open(temp_path+'/temp.sam', 'rt') as handle:
        lines = (line.strip().split('\t', 10) for line in handle if line[0] != '@')
        scores = []
        for line in lines:
            seq = line[0]
            alts = xa.search(line[-1])
            # print(seq)
            # print(alts)
            if not alts:
                yield seq, 0, 0,{}
                continue
            hits = [
                (x, y, z)
                for x, y, z in [(line[2], '-+'[seq == line[9].upper()], line[3])] + ot.findall(alts.group())
                if not re.match(exclude, x)]
            # print(hits)
            hit_list = extract_seq(hits, genome, n=23)
            # print(hit_list)
            MIT_score= fz_score(seq, hit_list)
            idx_list = alignment_list(seq, hit_list)
            location_list = [ value for idx, value in enumerate(hits) if idx in idx_list]
            seq_list = [value for idx, value in enumerate(hit_list) if idx in idx_list]
            cfdScores = [ calcCfdScore(seq, x) for x in hit_list[1:]]
            cfdScores = [0 if x== None else x for x in cfdScores]
            # print(sum(cfdScores))
            guideCfdScore = calcMitGuideScore(sum(cfdScores))

            target_dict = {}
            for x in range(len(location_list)):
                if seq_list[x] not in target_dict:
                    target_dict.setdefault(seq_list[x], []).append(location_list[x])
                else:
                    target_dict[seq_list[x]].append(location_list[x])
    
            yield seq, MIT_score, guideCfdScore, target_dict
            #print('test################')
            #print(scores)
            # if nontargeting != (max(scores) == 1.):
            #     yield seq, int(round(100 / (1 + sum(scores + [0]))))


def run_bwa(refs, seqs):
    if not all([os.path.exists('{}.{}'.format(refs, x)) for x in ('amb', 'ann', 'bwt', 'pac', 'sa')]):
        subprocess.call('bwa index "{}"'.format(refs), shell=True)
    with open(temp_path+'/temp.fna', 'w') as handle:
        handle.write('\n'.join('>{0}\n{0}'.format(seq.upper()) for n, seq in enumerate(seqs)))
    command = 'bwa aln -k 4 -l 20 -n 4 -o 0 -N -t 4 "{}" "{}"/temp.fna > "{}"/temp.sai'.format(refs, temp_path, temp_path)
    subprocess.call(command, shell=True)
    command = 'bwa samse -n 60000 "{}" "{}"/temp.sai "{}"/temp.fna > "{}"/temp.sam'.format(refs, temp_path, temp_path, temp_path)
    subprocess.call(command, shell=True)




# def score_guide(seq, alts, genome, n):
#     pams = (
#         dict(
#             GG=1., gg=1., Gg=1., gG=1., 
#             AG=.2, ag=.2, Ag=.2, aG=.2, 
#             GA=.2, ga=.2, Ga=.2, gA=.2), 
#         dict(
#             CC=1., cc=1., Cc=1., cC=1., 
#             CT=.2, ct=.2, Ct=.2, cT=.2, 
#             TC=.2, tc=.2, Tc=.2, tC=.2))

#     # print("############seq")
#     # print(seq)
#     # print("############alts")
#     # print(alts)
#     # print("############genome")
#     # print(genome)
#     # print("############n")
#     # print(n)
#     for seqid, strand, pos in alts:
# 	# print("############seqid")
# 	# print(seqid)
# 	# print("############strand")
# 	# print(strand)
# 	# print("############pos")
# 	# print(pos)
#         x = int(pos) - 1
#         if strand == '+':
#             score = pams[0].get(str(genome[seqid][x + n + 1: x + n + 3].seq), 0.)
#             if score:
#                 alt = genome[seqid][x: x + n].upper()
#                 if len(alt) == len(seq):
#                     yield score * fz_score(seq, alt)
#         if strand == '-':
#             score = pams[1].get(str(genome[seqid][x - 3: x - 1].seq), 0.)
#             if score:
#                 alt = genome[seqid][x: x + n].upper().reverse_complement()
#                 if len(alt) == len(seq):
#                     yield score * fz_score(seq, alt)

def extract_seq(alts, genome, n):
    # print("############seq")
    # print(seq)
    # print("############alts")
    # print(alts)
    # print("############genome")
    # print(genome)
    # print("############n")
    # print(n)
    seq_list = []
    for seqid, strand, pos in alts:
	# print("############seqid")
	# print(seqid)
	# print("############strand")
	# print(strand)
	# print("############pos")
	# print(pos)
        x = int(pos) - 1
        if strand == '+':
            alt = genome[seqid][x: x + n].upper()
            alt = str(alt.seq)
            seq_list.append(alt)

        if strand == '-':
            alt = genome[seqid][x: x + n].upper().reverse_complement()
            alt = str(alt.seq)
            seq_list.append(alt)
    #print(seq_list)
    return seq_list

def alignment_list(guideSeq, otSeqs, maxMm=None, minHitScore=None, minAltHitScore=None):
    """
    adapted from https://snippets.siftie.com/harijay/fz-score-d1324dab/
    linked in CRISPOR paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4934014/
    """
    # n = -len(b)
    # # The Patrick Hsu weighting scheme
    # M = [
    #     0.000, 0.000, 0.014, 0.000, 0.000, 0.395, 0.317, 0.000, 0.389, 0.079,
    #     0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]
    # if n > -19:
    #     return 0
    # I = [i for i, (j, k) in islice(enumerate(zip(a[n:], b)), 20)  if j != k]
    # dists = [i - j for i, j in zip(I[1:], I[:-1])]
    # denom = 1. if not dists else (1 - sum(dists) / len(dists) / 19) * 4 + 1
    # score = 1. if not I else 1. / len(I) ** 2
    # for i in I:
    #     score *= 1 - M[i + n]
    # return score / denom
    altOts = []
    mainOts = []
    for ot in otSeqs:
        if ot.endswith("AG") or ot.endswith("GA"):
            altOts.append(ot)
        if ot.endswith("GG"):
            mainOts.append(ot)

    # calc and filter hit scores
    mainHitScores = [calcHitScore(guideSeq, ot) for ot in mainOts]
    mainHitScores = [h for h in mainHitScores if h > minHitScore]
    altHitScores  = [calcHitScore(guideSeq, ot) for ot in altOts]
    altHitScores  = [h for h in altHitScores if h > minAltHitScore]

    mainHitScores.extend(altHitScores)
    index_list = [idx for idx, value in enumerate(mainHitScores) if value == 100]
    return index_list

def dict_reformat(target_dict, seq):
    """used to extract perfect match location"""
    new_target_dict = {}
    for key, values in target_dict.items():
        if key[:20] == seq[:20] and key[21:]=='GG':
            new_target_dict[key] = values
    return new_target_dict


def off_target_format(specificities, seq):
    """used for return sgRNA position, off-target position"""
    new_dict = dict_reformat(specificities.get(seq, -1)[-1], seq)
    if seq in new_dict:
        perfect_match_seq_list = new_dict[seq]
        perfect_match_seq_list = ["_".join(list(x))for x in perfect_match_seq_list]
        perfect_match_num = len(perfect_match_seq_list)
        off_target_list = [item for sublist in new_dict.values() for item in sublist ]
        off_target_list = [ "_".join(list(x))for x in off_target_list]
        if set(off_target_list) == set(perfect_match_seq_list):
            off_target_list = 'None'
            off_target_num = 0
        else:
            off_target_list = list(set(off_target_list) - set(perfect_match_seq_list))
            off_target_num = len(off_target_list)
    else:
        perfect_match_seq_list = []
        perfect_match_num = 0
        off_target_list = []
        off_target_num = 0
    return [perfect_match_seq_list, perfect_match_num, off_target_list, off_target_num]
        
    

def main(targets, reference, exclude=None, nontargeting=False, number=10):

    #print(exclude)
    #print(nontargeting)
    mm_scores,pam_scores = get_mm_pam_scores()
    test = design_test(0)
    yield '#locus', 'locusID', 'position', 'strand', 'guide', 'PAM', 'cloneable', 'MIT_score', 'CFD_score','Azimuth', 'cds_len', 'perfect_position_list', 'perfect_match_count', 'off_targets_list', 'off_targets_number'
    with open(targets, 'rt') as handle:
        for n, target in enumerate(SeqIO.parse(handle, 'fasta')):
            print('# Target=', target.id, 'length=', len(target.seq), 'bp')
            hits = list(get_pams(str(target.seq)))
            #print(hits)
            if len(hits) > 0:
                activities = azimuth.model_comparison.predict(
                    numpy.array([y.upper() for x, y, a, b in hits]), aa_cut=None, percent_peptide=None)
                #print(activities)
                specificities = dict((x, [y, z, l])for x, y, z, l in parse_bwa(reference, [y[4:27] for x, y, a, b in hits], exclude, nontargeting))
                # specificities = dict(list(parse_bwa(reference, [y[4:27] for x, y in hits], exclude, nontargeting)))
                # print(specificities)
                result = [(v,w, x[:20], y, test(x[:20]), specificities.get(x, -1)[0], specificities.get(x, -1)[1], z, len(target.seq), off_target_format(specificities, x)[0], off_target_format(specificities, x)[1], off_target_format(specificities, x)[2], off_target_format(specificities, x)[3]) for v, w, x, y, z in (
                    (i[3], i[2], i[1][4:27].upper(), i[1][24:27].upper(), j) for i, j in zip(hits, activities))]
                print('# Obtained', len(result), 'candidate guides.')
                result.sort(key=lambda x: (x[3], x[4], x[5]), reverse=True)
                for row in result[:number]:
                    yield (n + 1, target.id) + row
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('targets', help='path to FASTA file with target sequences')
    parser.add_argument('reference', help='path to FASTA file with reference genome sequences')
    parser.add_argument('save', help='path to save results')
    parser.add_argument(
        '--exclude',
        default='^$',
        help='regex pattern for seqids to exclude from references (e.g. alternate chromosomes)')
    parser.add_argument('--nontargeting', action='store_true', help='guides should not occur in genome')
    parser.add_argument('--number', type=int, default=10, help='number of guides to pick per target record')
    args = parser.parse_args()
    start = time()
    save_path = args.save
    temp_path = ('/').join(save_path.split('/')[:-1])
    with open(args.save, 'wt') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerows(main(
            args.targets, args.reference, exclude=args.exclude, nontargeting=args.nontargeting, number=args.number))
    print('# Done in', time() - start, 'seconds')
