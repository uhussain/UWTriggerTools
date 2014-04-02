#!/usr/bin/env python
import ROOT
import numpy

eff_file = ROOT.TFile("../data/LSB50/uct_tau_efficiency.root")
rate_file = ROOT.TFile("../data/LSB50/uct_rates_eic3.root")

tau_eff_ntuple = eff_file.Get("rlxTauEfficiency/Ntuple")
tau_rate_ntuple = rate_file.Get("rlxTauUCTRate/Ntuple")

# Setup nice aliases
#
# 1) just pt
# 2) max(pt, regionPt)
# 2) regionPt + 2ndRegionPt

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

offline_cut = 30


def get_resolution(variable, cut=offline_cut):
    tau_eff_ntuple.Draw(
        "1-(%s/recoPt)>>htemp(100, -5, 5)" % variable,
        "l1gMatch && recoPt > %i" % cut)
    histo = ROOT.gDirectory.Get("htemp")
    return "%0.2f %0.2f" % (histo.GetMean(), histo.GetRMS())

print "pt1 20", get_resolution("pt1", 20)
print "pt1 50", get_resolution("pt1", 50)
print "pt2 20", get_resolution("pt2", 20)
print "pt2 50", get_resolution("pt2", 50)
print "pt3 20", get_resolution("pt3", 20)
print "pt3 50", get_resolution("pt3", 50)

# we want to find the best combo of
#
# A*pt1 + B*pt2 + C*pt3 + D*pu + E

coefficients = []
targets = []

for row in tau_eff_ntuple:
    if row.recoPt > offline_cut and row.l1gMatch:
        coefficients.append([
            #row.l1gPt,
            row.l1gRegionEt + row.l1g2ndRegionEt - row.l1gPt,
            row.l1gRegionEt,
            row.l1g2ndRegionEt,
            row.l1gPUUIC, 1
        ]
        )
        targets.append(row.recoPt)


coeffs, residuals, rank, s = numpy.linalg.lstsq(coefficients, targets)

print "pt4 coeff (" + ", ".join("%0.3f" % x for x in coeffs) + ")"

tau_eff_ntuple.SetAlias(
    "pt4",
    "%0.3f * l1gPt"
    #"%0.3f * (l1gRegionEt + l1g2ndRegionEt - l1gPt)"
    "+ %0.3f * l1gRegionEt"
    "+ %0.3f * l1g2ndRegionEt"
    "+ %0.3f * l1gPUUIC + %0.3f" % tuple(coeffs))

tau_rate_ntuple.SetAlias(
    "pt4",
    "%0.3f * pt"
    #"%0.3f * (regionPt + secondRegionPt - pt)"
    "+ %0.3f * regionPt"
    "+ %0.3f * secondRegionPt"
    "+ %0.3f * puUIC + %0.3f" % tuple(coeffs))

print "pt4 20", get_resolution("pt4", 20)
print "pt4 50", get_resolution("pt4", 50)

# let's try another method (more stable)

coefficients = []
targets = []

for row in tau_eff_ntuple:
    if row.recoPt > offline_cut and row.l1gMatch:
        coefficients.append([
            #row.l1gPt,
            row.l1gRegionEt + row.l1g2ndRegionEt + row.l1gPt,
            row.l1gRegionEt - row.l1g2ndRegionEt + row.l1gPt,
            row.l1gRegionEt + row.l1g2ndRegionEt - row.l1gPt,
            row.l1gPUUIC, 1
        ]
        )
        targets.append(row.recoPt)


coeffs, residuals, rank, s = numpy.linalg.lstsq(coefficients, targets)

print "pt4 coeff (" + ", ".join("%0.3f" % x for x in coeffs) + ")"

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
    "+ %0.3f * puUIC + %0.3f" % tuple(coeffs))

print "pt5 20", get_resolution("pt5", 20)
print "pt5 50", get_resolution("pt5", 50)

# isolation calibration
iso_coefficients = []
iso_targets = []

for row in tau_eff_ntuple:
    if row.recoPt > offline_cut and row.l1gMatch and row.l1gJetPt > 0:
        pt4 = sum(x * y for (x, y) in zip(
            coeffs,
            [row.l1gPt, row.l1gRegionEt, row.l1g2ndRegionEt, row.l1gPUUIC, 1]
        ))
        iso_coefficients.append(
            [row.l1gPUUIC,  1]
        )
        #iso_targets.append(row.l1gJetPt - pt4)
        iso_targets.append(row.l1gJetPt - row.l1gRegionEt - row.l1g2ndRegionEt)

iso_coeffs, iso_residuals, iso_rank, iso_s = numpy.linalg.lstsq(
    iso_coefficients, iso_targets)

print ", ".join("%0.3f" % x for x in iso_coeffs)

tau_eff_ntuple.SetAlias(
    "isoCorr",
    "l1gJetPt - l1gRegionEt - l1g2ndRegionEt"
    "- %0.3f * l1gPUUIC - %0.3f" % tuple(iso_coeffs))

tau_rate_ntuple.SetAlias(
    "isoCorr",
    "jetPt - regionPt - secondRegionPt"
    "- %0.3f * puUIC - %0.3f" % tuple(iso_coeffs))


def make_turnon(l1_cut, color):
    tau_eff_ntuple.Draw("recoPt>>denom(25, 0, 100)", "", "goff")
    tau_eff_ntuple.Draw("recoPt>>num(25, 0, 100)", l1_cut, "goff")
    denom = ROOT.gDirectory.Get("denom")
    num = ROOT.gDirectory.Get("num")
    graph = ROOT.TGraphAsymmErrors(num, denom)
    graph.SetLineColor(color)
    graph.SetLineWidth(2)
    graph.Fit("tanh")
    return graph


def get_rate(l1_var, cut):
    return tau_rate_ntuple.GetEntries(
        "MaxIf$(%s, pt > 0) > %i" % (l1_var, cut))


canvas = ROOT.TCanvas("asdf", "adsf", 600, 600)

frame = ROOT.TH1F('frame', 'frame', 50, 0, 100)
frame.SetMaximum(1.1)

frame.Draw()

pt1_turnon = make_turnon("pt1 > 30", ROOT.EColor.kRed)
#pt1_turnon.Draw('lp')

pt2_turnon = make_turnon("pt2 > 22", ROOT.EColor.kBlue)
pt2_turnon.Draw('lx')

pt3_turnon = make_turnon("pt3 > 30", ROOT.EColor.kGreen)
pt3_turnon.Draw('lx')

pt4_turnon = make_turnon("pt4 > 30", ROOT.EColor.kBlack)
pt4_turnon.Draw('lx')

pt5_turnon = make_turnon("pt5 > 30", ROOT.EColor.kOrange)
pt5_turnon.Draw('lx')

legend = ROOT.TLegend(0.6, 0.2, 0.89, 0.5, '', 'brNDC')
legend.SetFillStyle(0)
legend.AddEntry(pt2_turnon, "Current - %0.f" % get_rate("pt2", 22), "l")
#legend.AddEntry(pt3_turnon, "Add 2nd - %0.f" % get_rate("pt3", 30), "l")
#legend.AddEntry(pt4_turnon, "Combo1 - %0.f" % get_rate("pt4", 30), "l")
legend.AddEntry(pt5_turnon, "Combo2 - %0.f" % get_rate("pt5", 30), "l")

legend.Draw()

canvas.SaveAs("turnons.png")

print "enhancing tree"

temp_file = ROOT.TFile("/scratch/efriis/tmpfile.root", "RECREATE")
enhanced_rate_ntuple = tau_rate_ntuple.CopyTree("pt5 > 30")

print "checking discriminants"

COLOR = 2
roc_graphs = []
roc_legend = ROOT.TLegend(0.2, 0.6, 0.5, 0.89, '', 'brNDC')
roc_legend.SetFillStyle(0)
roc_multigraph = ROOT.TMultiGraph()


def compare_discriminant(discS, discB, cutS, cutB, binning, name,
                         invert=False):
    global COLOR
    COLOR += 1
    tau_eff_ntuple.Draw("%s>>signal(%s)" % (discS, binning), cutS)
    enhanced_rate_ntuple.Draw("%s>>bkg(%s)" % (discB, binning), cutB)
    signal = ROOT.gDirectory.Get("signal")
    bkg = ROOT.gDirectory.Get("bkg")
    print name, signal.Integral(), bkg.Integral()
    signal.Scale(1.0 / signal.Integral())
    bkg.Scale(1.0 / bkg.Integral())
    bkg.SetLineColor(ROOT.EColor.kRed)
    signal.Draw()
    bkg.Draw('same')
    signal_int = signal.GetIntegral()
    bkg_int = bkg.GetIntegral()
    roc_curve = []
    wp_90_found = False
    wp_95_found = False
    wp_85_found = False
    wp_80_found = False
    wp_60_found = False
    for i in range(signal.GetNbinsX()):
        sig_eff = signal_int[i]
        bkg_eff = bkg_int[i]
        if invert:
            sig_eff = 1 - sig_eff
            bkg_eff = 1 - bkg_eff

        if sig_eff > 0.60 and not wp_60_found:
            print "0.60 WP", signal.GetBinCenter(i), \
                "eff:", "%0.2f" % sig_eff, "rate:", "%0.2f" % bkg_eff
            wp_60_found = True

        if sig_eff > 0.80 and not wp_80_found:
            print "0.80 WP", signal.GetBinCenter(i), \
                "eff:", "%0.2f" % sig_eff, "rate:", "%0.2f" % bkg_eff
            wp_80_found = True

        if sig_eff > 0.85 and not wp_85_found:
            print "0.85 WP", signal.GetBinCenter(i), \
                "eff:", "%0.2f" % sig_eff, "rate:", "%0.2f" % bkg_eff
            wp_85_found = True

        if sig_eff > 0.90 and not wp_90_found:
            print "0.90 WP", signal.GetBinCenter(i), \
                "eff:", "%0.2f" % sig_eff, "rate:", "%0.2f" % bkg_eff
            wp_90_found = True

        if sig_eff > 0.95 and not wp_95_found:
            print "0.95 WP", signal.GetBinCenter(i), \
                "eff:", "%0.2f" % sig_eff, "rate:", "%0.2f" % bkg_eff
            wp_95_found = True

        if not roc_curve or roc_curve[-1] != (sig_eff, bkg_eff):
            roc_curve.append((sig_eff, bkg_eff))
    roc_graph = ROOT.TGraph(len(roc_curve))
    for i, (x, y) in enumerate(roc_curve):
        roc_graph.SetPoint(i, x, y)
    roc_graph.SetLineColor(COLOR)
    roc_graphs.append(roc_graph)
    roc_multigraph.Add(roc_graph, "l")
    roc_legend.AddEntry(roc_graph, name, "l")
    canvas.SaveAs(name + ".png")
    return roc_graph

iso2_rel = compare_discriminant(
    "iso2/pt2", "iso2/pt2",
    "pt5 > 30 && l1gJetPt > 0", "pt5 > 30 && jetPt > 0",
    "100, -2, 3", "iso2_rel")

iso2_pt5_rel = compare_discriminant(
    "iso2/pt5", "iso2/pt5",
    "pt5 > 30 && l1gJetPt > 0", "pt5 > 30 && jetPt > 0",
    "100, -2, 3", "iso2_pt5_rel")

iso2 = compare_discriminant(
    "iso2", "iso2",
    "pt5 > 30 && l1gJetPt > 0", "pt5 > 30 && jetPt > 0",
    "100, -20, 100", "iso2", invert=False)

iso3 = compare_discriminant(
    "iso3", "iso3",
    "pt5 > 30 && l1gJetPt > 0", "pt5 > 30 && jetPt > 0",
    "100, -20, 100", "iso3", invert=False)

iso3 = compare_discriminant(
    "iso3", "iso3",
    "pt5 > 30 && l1gJetPt > 0 && pt1/pt3 > 0.45",
    "pt5 > 30 && jetPt > 0 && pt1/pt3 > 0.45",
    "100, -20, 100", "iso3_pt1pt3", invert=False)

iso3_rel = compare_discriminant(
    "iso3/pt5", "iso3/pt5",
    "pt5 > 30 && l1gJetPt > 0", "pt5 > 30 && jetPt > 0",
    "100, -2, 3", "iso3_rel")

isoCorr = compare_discriminant(
    "isoCorr", "isoCorr",
    "pt5 > 30 && l1gJetPt > 0", "pt5 > 30 && jetPt > 0",
    "100, -20, 50", "isoCorr", invert=False)

isoCorr_rel = compare_discriminant(
    "isoCorr/pt5", "isoCorr/pt5",
    "pt5 > 30 && l1gJetPt > 0", "pt5 > 30 && jetPt > 0",
    "100, -2, 2", "isoCorr_rel")

#pt1pt5 = compare_discriminant(
    #"pt1/pt5", "pt1/pt5",
    #"pt5 > 30 && l1gJetPt > 0", "pt5 > 30 && jetPt > 0",
    #"100, -2, 2", "pt1pt5", invert=True)

pt1pt3 = compare_discriminant(
    "pt1/pt3", "pt1/pt3",
    "pt5 > 30 && l1gJetPt > 0", "pt5 > 30 && jetPt > 0",
    "100, -2, 2", "pt1pt3", invert=True)

pt1pt3 = compare_discriminant(
    "pt2/pt3", "pt2/pt3",
    "pt5 > 30 && l1gJetPt > 0", "pt5 > 30 && jetPt > 0",
    "100, -2, 2", "pt2pt3", invert=True)

pt1pt2 = compare_discriminant(
    "pt1/pt2", "pt1/pt2",
    "pt5 > 30 && l1gJetPt > 0", "pt5 > 30 && jetPt > 0",
    "100, -2, 2", "pt1pt2", invert=True)

pt1 = compare_discriminant(
    "pt1", "pt1",
    "pt5 > 30 && l1gJetPt > 0", "pt5 > 30 && jetPt > 0",
    "100, 0, 200", "pt1", invert=True)

canvas.Clear()

roc_multigraph.Draw("al")
roc_legend.Draw()
canvas.SaveAs("roc_curves.png")
