from module import *
import wremnants

class defineWeight(module):
   
    def __init__(self, dataset,isWlike,eventCount,pileup_helper=None,muon_efficiency_helper=None,muon_prefiring_helper=None):

        self.dataset = dataset
        self.eventCount = eventCount
        self.pileup_helper = pileup_helper
        self.muon_efficiency_helper = muon_efficiency_helper
        self.muon_prefiring_helper = muon_prefiring_helper
        self.isWlike = isWlike

    def run(self,d):

        self.d = d

        if self.isWlike:
            flag = "trigMuons"
        else:
            flag = "goodMuons"
        if not self.dataset.is_data:
            if not self.eventCount:
                
                self.d = self.d.Define("weight_pu", self.pileup_helper, ["Pileup_nTrueInt"])
                self.d = self.d.Define("weight_fullMuonSF_withTrackingReco", self.muon_efficiency_helper, [f"{flag}_pt0", f"{flag}_eta0", f"{flag}_SApt0", f"{flag}_SAeta0", f"{flag}_charge0","passIso"])
                self.d = self.d.Define("weight_newMuonPrefiringSF", self.muon_prefiring_helper, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId"])

                self.d = self.d.Define("nominal_weight", "weight*weight_pu*weight_fullMuonSF_withTrackingReco*weight_newMuonPrefiringSF*L1PreFiringWeight_ECAL_Nom")
            else:
                self.d = self.d.Define("weight", "std::copysign(1.0, genWeight)")
                self.d = self.d.Define("unity", "1.0")
        else:
            if self.eventCount:
                self.d = self.d.DefinePerSample("weight", "1.0")
        
        return self.d
