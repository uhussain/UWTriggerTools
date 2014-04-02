#!/usr/bin/env python
import ROOT
import random

eff_file = ROOT.TFile("/scratch/efriis/uct_tau_efficiency.root")
rate_file = ROOT.TFile("/scratch/efriis/uct_rates_eic3.root")

tau_eff_ntuple = eff_file.Get("rlxTauEfficiency/Ntuple")
tau_rate_ntuple = rate_file.Get("rlxTauUCTRate/Ntuple")

tau_eff_ntuple.SetAlias("pt1", "l1gPt")
tau_eff_ntuple.SetAlias("pt2", "max(l1gPt, l1gRegionEt + 0)")
tau_eff_ntuple.SetAlias("pt3", "max(l1gPt, l1g2ndRegionEt + l1gRegionEt)")
tau_eff_ntuple.SetAlias("iso2", "l1gJetPt - pt2")
tau_eff_ntuple.SetAlias("iso3", "l1gJetPt - pt3")

tau_rate_ntuple.SetAlias("pt1", "pt")
tau_rate_ntuple.SetAlias("pt2", "max(pt, regionPt + 0)")
tau_rate_ntuple.SetAlias("pt3", "max(pt, regionPt + secondRegionPt)")
tau_rate_ntuple.SetAlias("iso2", "jetPt - pt2")
tau_rate_ntuple.SetAlias("iso3", "jetPt - pt3")

coeffs = (0.235, 0.116, 0.236, -0.098, 14.649)  # pt 20

#coeffs = (0.244, 0.068, 0.350, -0.033, 19.311)  # pt 30

tau_eff_ntuple.SetAlias(
    "pt5",
    "%0.3f * (l1gRegionEt + l1g2ndRegionEt + l1gPt)"
    "+ %0.3f * (l1gRegionEt - l1g2ndRegionEt + l1gPt)"
    "+ %0.3f * (l1gRegionEt + l1g2ndRegionEt - l1gPt)"
    "+ %0.3f * l1gPUUIC + %0.3f" % tuple(coeffs))

tau_rate_ntuple.SetAlias(
    "pt5",
    "%0.3f * (regionPt + secondRegionPt + pt)"
    "+ %0.3f * (regionPt - secondRegionPt + pt)"
    "+ %0.3f * (regionPt + secondRegionPt - pt)"
    "+ %0.3f * pu + %0.3f" % tuple(coeffs))


def get_pt95(formula):
    for i in range(100):
        if formula.Eval(i) / formula.GetParameter(0) > 0.90:
            return i
    return 100

canvas = ROOT.TCanvas("asdf", "adsf", 600, 600)


def make_turnon(l1_cut, color, name, fix_plateau=None):
    tau_eff_ntuple.Draw("recoPt>>denom(20, 0, 100)", "", "goff")
    tau_eff_ntuple.Draw("recoPt>>num(20, 0, 100)", l1_cut, "goff")
    denom = ROOT.gDirectory.Get("denom")
    num = ROOT.gDirectory.Get("num")
    graph = ROOT.TGraphAsymmErrors(num, denom)
    graph.SetMarkerColor(color)
    graph.SetMarkerStyle(20)
    formula = ROOT.TF1('fitter' + str(random.randint(0, 1000)),
                       "[0]*(0.5*(tanh([1]*(x-[2]))+1))", 0, 100)
    formula.SetParName(0, 'plateau')
    formula.SetParName(1, 'width')
    formula.SetParName(2, 'threshold')
    formula.SetParameter(0, 0.95)
    formula.SetParLimits(0, 0.1, 1.01)
    if fix_plateau is not None:
        print 'fix plat'
        formula.FixParameter(0, fix_plateau)
    formula.SetParameter(1, 1. / 10)
    formula.SetParameter(2, 30)
    formula.SetLineColor(color)
    graph.Fit(formula)
    # refit, not using the inefficient junk at start of curve
    formula.SetRange(formula.GetParameter(2), 100)
    graph.Fit(formula, "R")
    graph.Draw("ape")
    graph.GetHistogram().SetMaximum(1.1)
    graph.GetHistogram().SetMinimum(0)
    pt95 = get_pt95(formula)
    graph.GetHistogram().GetXaxis().SetTitle('pt95 = %i' % pt95)
    canvas.SaveAs(name + ".png")
    return graph, formula, pt95

pt5_95_map = {}
pt5_iso3_95_map = {}
pt5_iso3rel_95_map = {}

pt2_95_map = {}
pt2_iso2_95_map = {}

mapping = {}

colors = [ROOT.EColor.kRed, ROOT.EColor.kBlue, ROOT.EColor.kOrange,
          ROOT.EColor.kViolet, ROOT.EColor.kGreen, ROOT.EColor.kBlack,
          ROOT.EColor.kTeal, ROOT.EColor.kPink]

print 'enhancing'
temp_file = ROOT.TFile("/scratch/efriis/tmpfile.root", "RECREATE")
enhanced_rate_ntuple = tau_rate_ntuple.CopyTree(
    "pt > 19 || pt2 > 19 || pt5 > 19")


def get_rate(l1_var):
    return enhanced_rate_ntuple.GetEntries(
        "MaxIf$(1, %s) > 0.5" % (l1_var))


def map_cut(the_cut, l1cut, the_color, the_name):
    print "mapping", the_cut, l1cut
    cut_dict = mapping.setdefault(the_name, {})
    fix_plateau = None
    if 25 in cut_dict:
        fix_plateau = cut_dict[25][0][1].GetParameter(0)
    result = make_turnon(the_cut % l1cut, the_color,
                         ('turnon_%i_' % l1cut) + the_name,
                         fix_plateau=fix_plateau)
    rate = get_rate(the_cut % l1cut)
    print "rate", rate, 'pt95', result[2]
    cut_dict[l1cut] = (result, rate)
    print "\n"

l1cuts = [25, 30, 35, 40, 45, 50, 60]

for l1cut, color in zip(l1cuts, colors):
    map_cut('pt5 > %i', l1cut, color, 'pt5')
    map_cut('pt5 > %i && pt1/pt3 > 0.4', l1cut, color, 'pt5_pt1pt3')
    map_cut('pt5 > %i && pt1/pt3 > 0.6', l1cut, color, 'pt5_pt1pt3_2')
    map_cut('pt5 > %i && pt1/pt3 > 0.45 && iso3/pt5 < 0.07', l1cut, color,
            'pt5_pt1pt3_iso')
    map_cut('pt2 > %i', l1cut, color, 'pt2')
    map_cut('pt2 > %i && iso2/pt2 < 0.2', l1cut, color, 'pt2_iso2')
    map_cut('pt5 > %i && iso3/pt5 < 0.065', l1cut, color, 'pt5_iso3rel')
    #map_cut('pt5 > %i && iso3 < 3', l1cut, color, 'pt5_iso3')

rate_multigraph = ROOT.TMultiGraph()
legend = ROOT.TLegend(0.4, 0.7, 0.89, 0.89, '', 'brNDC')
legend.SetFillStyle(0)
rates = []

name_map = {
    'pt5_iso3': "NewPT AbsIso",
    'pt5_iso3rel': "NewPT NewIso",
    'pt5': "NewPT",
    'pt5_pt1pt3': "NewPT LdTrkLoose",
    'pt5_pt1pt3_2': "NewPT LdTrkTight",
    'pt5_pt1pt3_iso': "NewPT LdTrkLoose Iso",
    'pt2': "OldPT",
    'pt2_iso2': "OldPT OldIso",
}


for type, color in zip(mapping.keys(), colors):
    rate_curve = ROOT.TGraph(len(l1cuts))
    avg_plateau = 0
    for i, cut in enumerate(l1cuts):
        result, rate = mapping[type][cut]
        avg_plateau += result[1].GetParameter(0)
        rate_curve.SetPoint(i, result[2], rate)
    avg_plateau /= len(l1cuts)
    rate_curve.SetLineColor(color)
    rate_curve.SetLineWidth(2)
    rate_curve.SetMarkerStyle(20)
    rate_curve.SetMarkerSize(1)
    rate_curve.SetMarkerColor(color)
    rates.append(rate_curve)
    legend.AddEntry(
        rate_curve,
        name_map[type] + ": plateau eff: %0.0f%%" % (100 * avg_plateau), "lp")
    rate_multigraph.Add(rate_curve, "lp")


canvas.Clear()
canvas.SetLogy(True)
rate_multigraph.Draw('alp')
rate_multigraph.GetXaxis().SetTitle("Offline p_{T} @ 90% of plateau")
legend.Draw()
canvas.SaveAs("rates_multi.png")
