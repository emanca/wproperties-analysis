lumiperEra = { "B" : 5.824235614,
               "C" : 2.621295973,
               "D" : 4.285851496,
               "E" : 4.065974639,
               "F" : 2.717344923,
               "F_postVFP" : 0.418120616,
               "G" : 7.652808366,
               "H" : 8.739883636
}

lumi_preVFP=float(lumiperEra['B']+lumiperEra['C']+lumiperEra['D']+lumiperEra['E']+lumiperEra['F'])
lumi_postVFP=float(lumiperEra['F_postVFP']+lumiperEra['G']+lumiperEra['H'])
lumi_total2016=lumi_preVFP+lumi_postVFP

# print(lumi_preVFP,lumi_postVFP)
