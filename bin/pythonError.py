from sys import argv, exit, stdout, stderr
import sys
import ROOT


# So things don't look like crap.
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(1011)
######## File #########
if len(argv) < 2:
   print 'Usage:python pythontest.py RootFile.root SaveLabel[optional]'
   exit()

infile = argv[1]
ntuple_file = ROOT.TFile(infile)


ntuple = ntuple_file.Get("treeInclHF/Ntuple")

canvas = ROOT.TCanvas("asdf", "adsf", 800, 800)

######## LABEL & SAVE WHERE #########

if len(argv)>2:
   saveWhere='~/www/Research/'+argv[2]+'_'
else:
   saveWhere='~/www/Research/'


hist_eta = [] 
#
#ntuple.Draw("npvs:puMult0>>plotforbin(18,0,396)", "", "goff")
#plotforbin = ROOT.gDirectory.Get("plotforbin")
#numbins = plotforbin.GetXaxis().GetNbinsX()
#numbins = plotforbin.GetNbinsX()
#print numbins
# Create a 22x22 array of histograms
histos = []
for i in range(0,22):
    hname_eta = "hist_eta%d" %(i)
    hist_eta.append(ROOT.TH1F(hname_eta,"",18,0,17))
    histos.append([])
    for j in range(0,18):
        hname = "histos%d_%d" % (i, j) # Each histogram must have a unique name
        histos[i].append( ROOT.TH1F(hname,"",400,0,200) )

for event in ntuple:
    nentries = len(event.regionPt)
    nentriesEta = len(event.regionEta)
#    print 'nentries'
#    print nentries
    for pu in range(0,18):
        if event.puMult0 in range(pu*22,22+pu*22):
           #print 'pumult'
           #print event.puMult0
           for eta in range(0,22):
               for i in range(0,396):
                   if event.regionEta[i] == eta:
                      histos[eta][pu].Fill(event.regionPt[i])


print 'mean'
for i in range(0,22):
    print '\t'
    for j in range(0,18):
        Mean =histos[i][j].GetMean()
        MeanError =histos[i][j].GetMeanError()
        hist_eta[i].SetBinContent(j,Mean)
        hist_eta[i].SetBinError(j,MeanError)
        print '%d \t %d \t %f' %(i,j,Mean)
    save = 'hist_eta%d.png' % i  
    hist_eta[i].Draw("p")
    hist_eta[i].Fit("pol1")
    canvas.SaveAs(save) 


file=ROOT.TFile("outfile.root","RECREATE")
file.cd()

for i in range(0,22):
  hist_eta[i].Write()
  for j in range(0,18):
      histos[i][j].Write()                 


#for j in range(0,18):
#    print ' pu: %d' % j
#    for i in range (0,22):
#        Mean = histos[i][j].GetMean()
#        print 'eta: %d Mean: %f' %(i,Mean) 


