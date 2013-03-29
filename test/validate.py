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
eff_mode = any('l1' in x for x in branch_names)

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


def iterable(x):
    try:
        iter(x)
    except TypeError:
        return False
    return True

for i in range(entries):
    ntuple1.GetEntry(i)
    ntuple2.GetEntry(i)

    if eff_mode:
        if ntuple1.l1gPt > 15:
            num_above_15 += 1
    else:
        if ntuple1.pt > 15:
            num_above_15 += 1

    for branch in branch_names:
        if branch == 'type':
            continue
        results1 = getattr(ntuple1, branch)
        results2 = getattr(ntuple2, branch)

        if not iterable(results1):
            results1 = [results1]
            results2 = [results2]

        for j, (res1, res2) in enumerate(zip(results1, results2)):
            delta = abs(res1 - res2)
            if 'phi' in branch.lower():
                delta = deltaPhi(res1, res2)
            if 'uic' in branch.lower():
                delta = 0
            # we don't care about these branches in the tau files
            if 'mip' in branch.lower() and 'Tau' in ntuple:
                delta = 0
            if 'tauveto' in branch.lower() and 'Tau' in ntuple:
                delta = 0

            if delta > 1e-5:
                #print ntuple1.l1gPt, ntuple1.l1gEta, ntuple1.l1gPhi
                print "mismatch event %i-%i branch %s (%s %s != %s %s)" % (
                    i, j, branch, file1, str(res1), file2, str(res2))
                num_errors += 1

print "found %i errors in %i entries, with %i high-pt events" % (
    num_errors, entries, num_above_15)
