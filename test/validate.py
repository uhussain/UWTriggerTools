#!/usr/bin/env python

# Stupid check to ensure the validation files haven't changed
#
# Usage:
# ./validate.py mt_numEvent200.root validation/mt_numEvent200.root rlxTauEfficiency/Ntuple

import sys

file1, file2, ntuple = sys.argv[1], sys.argv[2], sys.argv[3]

print file1, file2, ntuple

# So ROOT doesn't monkey with them
sys.argv[:] = []

import ROOT

f1 = ROOT.TFile(file1)
f2 = ROOT.TFile(file2)

ntuple1 = f1.Get(ntuple)
ntuple2 = f2.Get(ntuple)

branch_names = (x.GetName() for x in ntuple1.GetListOfBranches())

entries = ntuple1.GetEntries()

num_errors = 0

num_above_15 = 0

for i in range(entries):
    ntuple1.GetEntry(i)
    ntuple2.GetEntry(i)
    if ntuple1.l1gPt > 15:
        num_above_15 += 1

    for branch in branch_names:
        res1 = getattr(ntuple1, branch)
        res2 = getattr(ntuple2, branch)
        if res1 != res2:
            print ntuple1.l1gPt
            print ntuple2.l1gPt
            print "mismatch event %i branch %s (%s %s != %s %s)" % (
                i, branch, file1, str(res1), file2, str(res2))
            num_errors += 1

print "found %i errors in %i entries, with %i high-pt events" % (
    num_errors, entries, num_above_15)
