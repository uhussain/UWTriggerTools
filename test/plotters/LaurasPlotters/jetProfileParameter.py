'''
Usage:python jetEfficiencyPlot.py RootFile.root label[optional]

Script to make some quick efficiency plots to test ntuple generation.


Author: L. Dodd, UW Madison

'''

from subprocess import Popen
from sys import argv, exit, stdout, stderr

import ROOT

# So things don't look like crap.
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

######## File #########
if len(argv) < 2:
   print 'Usage:python jetEffCorrAppPlot.py RootFile.root label[optional]'
   exit()

infile = argv[1]
ntuple_file = ROOT.TFile(infile)

######## LABEL & SAVE WHERE #########

if len(argv)>2:
   saveWhere='~/www/'+argv[2]+'_'
else:
   saveWhere='~/www/'
####### Calibration factor ####
L1_CALIB_FACTOR = 1.0
L1G_CALIB_FACTOR = 1.0


#####################################
#Get Jet Eff NTUPLE                 #
#####################################


jet_ntuple = ntuple_file.Get("corrjetEfficiency/Ntuple")
jet_ntuple_uncorr = ntuple_file.Get("jetEfficiency/Ntuple")

canvas = ROOT.TCanvas("asdf", "adsf", 800, 800)

def make_plot(tree, variable, selection, binning, xaxis='', title=''):
    ''' Plot a variable using draw and return the histogram '''
    draw_string = "%s>>htemp(%s)" % (variable, ", ".join(str(x) for x in binning))
    tree.Draw(draw_string, selection, "BOX")
    output_histo = ROOT.gDirectory.Get("htemp").Clone()
    output_histo.GetXaxis().SetTitle(xaxis)
    output_histo.SetTitle(title)
    return output_histo





def profile_ptoffset(ntuple, variable, l1PtCut_low, l1PtCut_high, binning, filename,
                         title='', xaxis=''):
    l1g = make_plot(
        ntuple, variable,
        "recoEta<3.0 && l1gMatch && %0.2f * recoPt > %0.2f && %0.2f * recoPt< %0.2f " % (L1G_CALIB_FACTOR, l1PtCut_low,L1G_CALIB_FACTOR,l1PtCut_high),
        binning
    )

    frame = ROOT.TH2F("frame", "frame", *binning)
    frame.SetMaximum(30)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.Draw()
    l1g.Draw('BOX')
    profilel1g=l1g.ProfileX(l1g.GetName() + "_profilel1g")
    profilel1g.SetMarkerStyle(23)
    l1g.GetXaxis().SetTitle("nPVs")
    l1g.GetYaxis().SetTitle("Reco Pt-Trigger Pt")
    profilel1g.SetLineWidth(2)
    return profilel1g
    #saveas = saveWhere+filename+'.png'
    #print saveas
    #canvas.SaveAs(saveas)

def profile_fit_parameter(ntuple, variable, l1PtCut_low, l1PtCut_high, binning, filename,
                         title='', xaxis=''):
    profile=profile_ptoffset(ntuple, variable, l1PtCut_low, l1PtCut_high,binning, filename,title='', xaxis='')
    profile.Fit("p1")
    p0=profile.GetParameter(0)
    return p0


def compare_profiles(ntuple, variable, ptbinsize, ptbinstart, ptbinnumber, binning, filename,title='', xaxis=''):

    bin1 = profile_ptoffset(ntuple, variable, ptbinstart, ptbinstart+ptbinsize, binning, filename,title='', xaxis='')
    bin2 = profile_ptoffset(ntuple, variable, ptbinstart+ptbinsize, ptbinstart+2*ptbinsize,binning, filename,title='', xaxis='')
    bin3 = profile_ptoffset(ntuple, variable, ptbinstart+2*ptbinsize, ptbinstart+3*ptbinsize,binning, filename,title='', xaxis='')
    bin4 = profile_ptoffset(ntuple, variable, ptbinstart+3*ptbinsize, ptbinstart+4*ptbinsize,binning, filename,title='', xaxis='')

    print bin1, bin2, bin3, bin4
    bin1.Draw()
    bin1.GetYaxis().SetRangeUser(-30,30.0)
    bin2.Draw("same")
    bin3.Draw("same")
    bin4.Draw("same")
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)


compare_profiles(jet_ntuple,'(recoPt-l1gPt):nPVs',10,30,4,[40,0,40,100,-50,50],'ptoffset_reco', "Jet Pt 3040","nPVs")
