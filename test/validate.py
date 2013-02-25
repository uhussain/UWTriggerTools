#!/usr/bin/env python

# Stupid check to ensure the validation files haven't changed

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

for i in range(entries):
    ntuple1.GetEntry(i)
    ntuple2.GetEntry(i)
    for branch in branch_names:
        res1 = getattr(ntuple1, branch)
        res2 = getattr(ntuple2, branch)
        if res1 != res2:
            print "mismatch event %i branch %s (%s != %s)" % (
                i, branch, str(res1), str(res2))

print "found no errors in %i entries" % entries
