#!/usr/bin/env python

import ROOT
import sys

ROOT.gStyle.SetOptFit(111111)

file = ROOT.TFile(sys.argv[1])
ntuple = file.Get("%s/Ntuple" % sys.argv[2])

canvas = ROOT.TCanvas("asdf", "asdf", 800, 800)

def fit(x, y):
    htemp = ROOT.gDirectory.Get('htemp')
    prof = htemp.ProfileX('blah')
    prof.Fit('pol1')
    htemp.SetFillColor(16)
    htemp.Draw('box')
    prof.SetLineColor(ROOT.EColor.kBlue)
    prof.Draw('same')
    htemp.GetXaxis().SetTitle(x)
    htemp.GetYaxis().SetTitle(y)
    return htemp, prof


ntuple.Draw("l1gPt:l1gRegionEt>>htemp(20, 0, 100, 20, 0, 100)", "l1gPt > 0", "box")

x = fit('4x4 p_{T}', '3x3 p_{T}')

canvas.SaveAs(sys.argv[2] + "_" +"3x3_vs_center4x4.png")

ntuple.Draw("l1gPt:l1g2ndRegionEt>>htemp(20, 0, 100, 20, 0, 100)", "l1gPt > 0", "box")

x = fit('4x4 #2 p_{T}', '3x3 p_{T}')

canvas.SaveAs(sys.argv[2] + "_" +"3x3_vs_next4x4.png")

ntuple.Draw("l1gPt:recoPt>>htemp(20, 0, 100, 20, 0, 100)", "l1gPt > 0", "box")
x = fit('RECO p_{T}', '3x3 p_{T}')

canvas.SaveAs(sys.argv[2] + "_" +"3x3_vs_reco.png")

ntuple.Draw("l1gPt:l1gJetPt>>htemp(20, 0, 100, 20, 0, 100)", "l1gPt > 0", "box")
x = fit('12x12 p_{T}', '3x3 p_{T}')

canvas.SaveAs(sys.argv[2] + "_" +"3x3_vs_12x12.png")
