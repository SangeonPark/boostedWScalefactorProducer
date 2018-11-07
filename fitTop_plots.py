import ROOT as r,sys,math,array,os
from optparse import OptionParser
from ROOT import *
import tdrstyle
import numpy as np
tdrstyle.setTDRStyle()

pt        = np.array([375.,437.5,500.,575.])
pt_unc    = np.array([0.,0.,0.,0.])
scale     = np.array([0.983,0.991,0.982,0.984])
scale_unc = np.array([0.004,0.003,0.004,0.004])
res       = np.array([1.086,1.192,1.091,1.167])
res_unc   = np.array([0.03,0.06,0.17,0.05])

n = 4


c1 = r.TCanvas("c1", "c1",800,600)
p1 = r.TGraphErrors(n,pt,scale,pt_unc,scale_unc)
f1 = p1.Fit("pol1")
#f1.SetLineWidth().SetLineWidth(2)
f1.Draw()
p1.SetTitle("Top mass scale")
p1.GetYaxis().SetRangeUser(0.96,1.01)
p1.GetXaxis().SetTitle("pT (GeV)")
p1.GetYaxis().SetTitle("Scale (#mu_{Data}/#mu_{MC})")
p1.Draw("ap")
c1.Update()
c1.SaveAs("scale.png")
c1.SaveAs("scale.pdf")


c2 = r.TCanvas("c2", "c2",800,600)
p2 = r.TGraphErrors(n,pt,res,pt_unc,res_unc)
f2 = p2.Fit("pol1")
f2.Draw()
p2.SetTitle("Top mass resolution")
p2.GetYaxis().SetRangeUser(0.8,1.4)
p2.GetXaxis().SetTitle("pT (GeV)")
p2.GetYaxis().SetTitle("Resolution (#sigma_{Data}/#sigma_{MC})")
p2.Draw("ap")
c2.Update()
c2.SaveAs("res.png")
c2.SaveAs("res.pdf")

