import ROOT
from array import array
ROOT.gROOT.SetBatch(True)

# '''
# 256 -             826
# 128 -             281.3
# 64 -            - 288.4
# 32 -            - 490.2
# 16              - 890.9
# '''
# '''data
# time=array('d', [890.9, 490.2, 288.4, 281.3, 826])
# cores=array('d', [16, 32, 64, 128, 256])
# ne=462515184
# '''

# '''
# WJets

# 52538760 events processed in 1169.1 s rate 558147.264571 histograms written:  6351
# runwj_16.log:652538760 events processed in 3732.9 s rate 174807.063997 histograms written:  6351
# runwj_32.log:652538760 events processed in 2057.2 s rate 317192.80834 histograms written:  6351
# runwj_64.log:652538760 events processed in 1304.4 s rate 500265.672499 histograms written:  6351
# runwj_128.log:652538760 events processed in 1156.3 s rate 564343.067005 histograms written:  6351
# runwj_256.log:652538760 events processed in 2615.6 s rate 249479.614437 histograms written:  6351
# '''
# time=array('d', [3732.9, 2057.2, 1304.4, 1156.3, 2615.6])
# cores=array('d', [16, 32, 64, 128, 256])
# ne=462515184

# rates=array('d')

# for i in range(0, 5):
#     print i, cores[i], '\t', time[i], '\t', ne/time[i]
#     rates.append((0.001*ne/time[i]))

# print "validation"
# for i in range(0, 5):
#     print i, cores[i], '\t', time[i], '\t', rates[i]

# sc=TGraph( 5, cores, rates)
# sc.SetName("scaling")
# sc.SetTitle( ";Ncores;Event processing rate [kHz] ")
# sc.SetMarkerStyle(21)
# sc.SetMarkerSize(1.2)
# c=TCanvas("scCanvas", "")
# c.cd()
# sc.Draw("ap")
# c.Update()
# c.Modified()
# c.SaveAs("scaling.root")
# c.SaveAs("scaling.png")


timeDict ={
    'Data' : {
        'time' : array('d', [890.9, 490.2, 288.4, 281.3, 826]),
        'cores' : array('d', [16, 32, 64, 128, 256]),
        'nev' : 462515184,
        'color' : ROOT.kBlack,
        'marker' : 20,
        'rates' : array('d'),
        'histo' : '13',
        'size' : '55.9 GB'
    },
    'Signal' : {
        'time' : array('d', [3732.9, 2057.2, 1304.4, 1156.3, 2615.6]),
        'cores' : array('d', [16, 32, 64, 128, 256]),
        'nev' : 652538760,
        'color' : ROOT.kRed+1,
        'marker' : 21,
        'rates' : array('d'),
        'histo' : '6351',
        'size' : '397 GB',

        
    }
}

cmslab = "#bf{CMS} #scale[0.7]{#it{Work in progress}}"
lumilab = " #scale[0.7]{13 TeV}"
cmsLatex = ROOT.TLatex()


graphDict = {}

# time=array('d', [890.9, 490.2, 288.4, 281.3, 826])
# ne=462515184


# time=array('d', [3732.9, 2057.2, 1304.4, 1156.3, 2615.6])
# cores=array('d', [16, 32, 64, 128, 256])
# ne=652538760

leg = ROOT.TLegend(0.5,0.75,0.9,0.91)
leg.SetFillStyle(0)
leg.SetLineWidth(0)

for k, val in timeDict.items() : 

    for i in range(0, len(val['cores'])):
        print(i, val['cores'][i], '\t', val['time'][i], '\t', val['nev']/val['time'][i])
        val['rates'].append((0.001*val['nev']/val['time'][i])) #in kHz

    graphDict[k]=ROOT.TGraph( len(val['cores']), val['cores'], val['rates'])
    graphDict[k].SetName("scaling")
    # graphDict[k].SetTitle( "Compared processing rate scaling ;N_{cores};Event processing rate [kHz] ")
    graphDict[k].SetTitle( " ;N_{cores};Event processing rate [kHz] ")
    graphDict[k].GetXaxis().SetTitleSize(0.06)
    graphDict[k].GetXaxis().SetTitleOffset(0.7)
    graphDict[k].GetYaxis().SetTitleSize(0.05)
    graphDict[k].GetYaxis().SetTitleOffset(0.93)
    graphDict[k].SetMarkerStyle(val['marker'])
    graphDict[k].SetMarkerColor(val['color'])
    graphDict[k].SetLineColor(val['color'])
    graphDict[k].SetLineWidth(2)
    graphDict[k].SetMarkerSize(2)
    graphDict[k].SetFillStyle(0)
    graphDict[k].SetLineStyle(2)
    # graphDict[k].Draw("ap")

    leg.AddEntry(graphDict[k], k+', N_{hist}='+val['histo']+', size='+val['size'])
 
c=ROOT.TCanvas("fw_time_scaling", "fw_time_scaling", 1600, 1200)
c.cd()
c.SetGridx()
c.SetGridy()
graphDict['Data'].Draw('apL')
graphDict['Signal'].Draw('pL same')
graphDict['Data'].GetYaxis().SetRangeUser(0,1800)
leg.Draw("same")

cmsLatex.SetNDC()
cmsLatex.SetTextFont(42)
cmsLatex.SetTextColor(ROOT.kBlack)
cmsLatex.SetTextAlign(31) 
cmsLatex.DrawLatex(1-c.GetRightMargin(),1-0.8*c.GetTopMargin(),lumilab)
cmsLatex.SetTextAlign(11) 
cmsLatex.DrawLatex(c.GetLeftMargin(),1-0.8*c.GetTopMargin(),cmslab)

   
c.Update()
c.Modified()
c.SaveAs("FW_time_scaling.pdf")
c.SaveAs("FW_time_scaling.png")
c.SaveAs("FW_time_scaling.root")




