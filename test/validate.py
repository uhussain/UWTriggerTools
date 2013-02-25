#!/usr/bin/env python

# Stupid check to ensure the validation files haven't changed
#
# Usage:
# ./validate.py mt_numEvent200.root validation/mt_numEvent200.root \
        # rlxTauEfficiency/Ntuple

import sys
import math

file1, file2, ntuple = sys.argv[1], sys.argv[2], sys.argv[3]

print file1, file2, ntuple

# So ROOT doesn't monkey with them
sys.argv[:] = []

import ROOT

f1 = ROOT.TFile(file1)
f2 = ROOT.TFile(file2)

ntuple1 = f1.Get(ntuple)
ntuple2 = f2.Get(ntuple)

branch_names = [x.GetName() for x in ntuple1.GetListOfBranches()]

entries = ntuple1.GetEntries()

num_errors = 0

num_above_15 = 0


def deltaPhi(x1, x2):
    diff = x1 - x2
    while diff > math.pi:
        diff = diff - 2 * math.pi
    if diff < -math.pi:
        diff += math.pi
    return diff

for i in range(entries):
    ntuple1.GetEntry(i)
    ntuple2.GetEntry(i)
    if ntuple1.l1gPt > 15:
        num_above_15 += 1

    for branch in branch_names:
        res1 = getattr(ntuple1, branch)
        res2 = getattr(ntuple2, branch)
        #print i, branch, res1, res2
        delta = abs(res1 - res2)
        if 'Phi' in branch:
            delta = deltaPhi(res1, res2)
        if 'UIC' in branch:
            delta = 0
        # we don't care about these branches in the tau files
        if 'MIP' in branch and 'Tau' in ntuple:
            delta = 0
        if 'TauVeto' in branch and 'Tau' in ntuple:
            delta = 0

        if delta > 1e-5:
            print ntuple1.l1gPt, ntuple1.l1gEta, ntuple1.l1gPhi
            print "mismatch event %i branch %s (%s %s != %s %s)" % (
                i, branch, file1, str(res1), file2, str(res2))
            num_errors += 1

print "found %i errors in %i entries, with %i high-pt events" % (
    num_errors, entries, num_above_15)
