#!/usr/bin/env python
import ROOT as r,sys,math,array,os
from optparse import OptionParser
from ROOT import *
import tdrstyle
import numpy as np 
import re
tdrstyle.setTDRStyle()

fVar="Puppijet0_msd"
fPt="Puppijet0_pt"
#fCut=["(Puppijet0_pt>350 && Puppijet0_pt<400)"]
fPt_bins = [350,400,475,550,100000]
fWeight="weight"
fLumi=41.1
fBinWidth =5
fDir="/afs/cern.ch/user/c/cmantill/work/public/forJeff/sklimWtag"
fPlotDir="plots/"
fOutput="top.root"

fFiles = {'data': 'Data',
          #'wlnu_mc' : 'WJets',
          #'wlnu_data': 'WJets',
          #'vv':   'VV', # ignore diboson for now
          'st_mc':   'ST',
          'st_data':   'ST',
          #'tt_mc':   'TT',
          #'tt_data':   'TT',
          'tt_fakeW_mc': 'TT',
          'tt_realW_mc': 'TT',
          'tt_fakeW_data': 'TT',
          'tt_realW_data': 'TT',
          'mc': 'PseudoData'
          }

fPDFs = {
         #'wlnu_mc':       'Exp_mc',
         #'wlnu_data':     'Exp_data',
         #'st_mc':        'GausErfExp_tt',
         'st_mc':         'PolyGaus_st_mc',
         'st_data':       'PolyGaus_st_data',
         #'tt':            'GausExp_tt',
         #'tt_mc':         'GausErfExp_tt_mc',
         #'tt_data':       'GausErfExp_tt_data',
         #'tt_data':       'GausGaus_tt_data',
         'tt_fakeW_mc':   'GausExp_tt_mc',
         'tt_realW_mc':   'GausErfExp_tt_mc',
         'tt_fakeW_data': 'GausExp_tt_data',
         'tt_realW_data': 'GausErfExp_tt_data',
         }

r.gSystem.Load("./PDFs/HWWLVJRooPdfs_cxx.so")
r.gSystem.Load("./PDFs/PdfDiagonalizer_cc.so")
r.gSystem.Load("./PlotStyle/Util_cxx.so")
r.gSystem.Load("./PlotStyle/PlotUtils_cxx.so")

def parser():
    parser = OptionParser()
    parser.add_option('--tag',     action='store', type='string', dest='tag',   default='CA15TopSelTau32', help='samples tag')
    parser.add_option('--xMin',    action='store', type='float',  dest='xmin',  default=40,                help='x-min')
    parser.add_option('--xMax',    action='store', type='float',  dest='xmax',  default=250,               help='x-max')
    parser.add_option('-b',        action='store_true',           dest='noX',   default=True,              help='no X11 windows')
    parser.add_option('--match',   action='store_true',           dest='match', default=False,             help='divide tt sample in matched and unmatched W')
    (options,args) = parser.parse_args()
    return options

# can do a couple of things here: colors, errors, add legends!
def drawFrame(iFrame,iData,iFuncs,iLegend,iColor,isPoisson=False):
    if isPoisson:
        iData.plotOn(iFrame,r.RooFit.DataError(r.RooAbsData.Poisson),r.RooFit.XErrorSize(0))
    else:
        iData.plotOn(iFrame)
    #iColor=50
    for pFunc in iFuncs:
        pFunc.plotOn(iFrame,r.RooFit.LineColor(iColor),r.RooFit.LineStyle(iColor != 50+1)) 
        #iLegend.AddEntry(pFunc, "Data", "l")
        pFunc.plotOn(iFrame,r.RooFit.Name("ErfExp"),r.RooFit.Components("erfExp*"),r.RooFit.LineColor(4))
        pFunc.plotOn(iFrame,r.RooFit.Name("Gaus"),r.RooFit.Components("gaus*"),r.RooFit.LineColor(3))
        pFunc.plotOn(iFrame,r.RooFit.Name("Poly"),r.RooFit.Components("poly*"),r.RooFit.LineColor(40))
        pFunc.plotOn(iFrame,r.RooFit.Name("Poly"),r.RooFit.Components("exp*"),r.RooFit.LineColor(6))

    print iFrame.Print("")

def fitSeparate(par):
    pt_bin_centers = [375.,437.5,500.,575.] 
    mc_mean   = [pt_bin_centers[i]*par["slope_mc"] + par["offset_mc"] for i in range(len(pt_bin_centers))]
    data_mean = [pt_bin_centers[i]*par["slope_data"] + par["offset_data"] for i in range(len(pt_bin_centers))]
    mc_mean_e = [np.sqrt(par["offset_mc_e"]**2 * pt_bin_centers[i]**2 + par["slope_mc_e"]**2) for i in range(len(pt_bin_centers))]
    data_mean_e = [np.sqrt(par["offset_data_e"]**2 * pt_bin_centers[i]**2 + par["slope_data_e"]**2) for i in range(len(pt_bin_centers))]
    scale     = np.asarray([data_mean[i]/mc_mean[i] for i in range(len(pt_bin_centers))])
    scale_e   = np.asarray([scale[i]*np.sqrt(mc_mean_e[i]**2/mc_mean[i]**2 + data_mean_e[i]**2/data_mean[i]**2) for i in range(len(pt_bin_centers))])
    print par
    print scale
    print scale_e
    pt_bin_centers = np.asarray(pt_bin_centers)
    pt_e      =np.asarray([0.,0.,0.,0.])
    lCan   = r.TCanvas("lCan","lCan",800,600)
    lP     = r.TGraphErrors(4,pt_bin_centers,scale)
    #lP     = r.TGraphErrors(4,pt_bin_centers,scale,pt_e,scale_e)
    pF     = lP.Fit("pol1")
    pF.Draw()
    lP.SetTitle("Top mass scale")
    lP.GetYaxis().SetRangeUser(0.96,1.01)
    lP.GetXaxis().SetTitle("pT (GeV)")
    lP.GetYaxis().SetTitle("Scale (#mu_{Data}/#mu_{MC})")
    lP.Draw("ap")
    lCan.Update()
    lCan.SaveAs("scale_pseudosimultaneous.png")
    lCan.SaveAs("scale_pseudosimultaneous.pdf")
    

def draw(iVar,iData,iFuncs,iRooFitResult,iLabel="A"):
    lCan   = r.TCanvas(str(iLabel),str(iLabel),800,600)
    lFrame = iVar.frame()
    lLegend = getLegend()
    drawFrame(lFrame,iData,iFuncs,lLegend,50)
    lFrame.GetYaxis().SetRangeUser(0,lFrame.GetMaximum()*1.2)
    lFrame.GetYaxis().SetTitle(" Events / %.1f GeV"%fBinWidth);
    lFrame.GetXaxis().SetTitle("PUPPI Softdrop Jet Mass (GeV) ");
    print iData,iFuncs,iRooFitResult,iVar
    lchisq,lndof = getChi2NDOF(iData,iFuncs,iRooFitResult,iVar)
    addInfo = getPavetext()
    addInfo.AddText("#chi^{2}/nDOF = %.3f/%i"%(lchisq,lndof))
    print "chi^{2}/nDOF = %.3f/%i"%(lchisq,lndof)
    lFrame.Draw()
    addInfo.Draw()
    lCan.Modified()
    lCan.Update()
    lCan.SaveAs(fPlotDir+iLabel+".pdf")
    lCan.SaveAs(fPlotDir+iLabel+".png")

def getLegend():
    legend = TLegend(0.8,0.75,0.95,0.9)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.040)
    legend.SetTextAlign(12)
    return legend

def getPavetext():
    addInfo = TPaveText(0.2,0.75,0.65,0.9,"NDC")
    addInfo.SetFillColor(1)
    addInfo.SetLineColor(1)
    addInfo.SetFillStyle(1)
    addInfo.SetBorderSize(0)
    addInfo.SetTextFont(42)
    addInfo.SetTextSize(0.040)
    addInfo.SetTextAlign(12)
    return addInfo

def drawDataMc(iVar,iData,iFuncsData,iMc,iFuncsMc,iRooFitResult_mc,iRooFitResult_data,params,iLabel='A'):
    lCan   = r.TCanvas(str(iLabel),str(iLabel),800,600)
    lFrame = iVar.frame()
    lLegend = getLegend()
    drawFrame(lFrame,iData,iFuncsData,lLegend,4,True)
    data_chisq,data_ndof = getChi2NDOF(iData,iFuncsData,iRooFitResult_data,iVar)
    mc_chisq,mc_ndof     = getChi2NDOF(iMc,iFuncsMc,iRooFitResult_mc,iVar)
    data_chi2Prob        = r.TMath.Prob(data_chisq,data_ndof)
    mc_chi2Prob          = r.TMath.Prob(mc_chisq,mc_ndof)
    #print 'scale factor + resolution for t sample: '
    #print 'scale =', params["data_mean_val"]/params["mc_mean_val"], "+/-", params["data_mean_val"]/params["mc_mean_val"]*np.sqrt((params["mc_mean_err"]/params["data_mean_val"])**2 + (params["mc_mean_err"]/params["mc_mean_val"])**2)
    #print 'resolution =', params["data_sigma_val"]/params["mc_sigma_val"], "+/-", params["data_sigma_val"]/params["mc_sigma_val"]*np.sqrt((params["data_sigma_err"]/params["data_sigma_val"])**2 + (params["mc_sigma_err"]/params["mc_sigma_val"])**2)
    #print  
    #print 'chisq/ndof for data =', data_chisq,'/',data_ndof
    #print 'chisq/ndof for mc =', mc_chisq,'/',mc_ndof
    drawFrame(lFrame,iMc,iFuncsMc,lLegend,2)
    lFrame.GetYaxis().SetRangeUser(0,lFrame.GetMaximum()*1.3)
    lFrame.GetYaxis().SetTitle(" Events / %.1f GeV"%fBinWidth);
    lFrame.GetXaxis().SetTitle( "PUPPI Softdrop Jet Mass (GeV)")
    lLegend.AddEntry(lFrame.getCurve("model_mc_Norm[Puppijet0_msd]"),"MC","l")
    lLegend.AddEntry(lFrame.getCurve("model_data_Norm[Puppijet0_msd]"),"Data","l")
    addInfo = getPavetext()
    addInfo.AddText("Data #chi^{2}/nDOF = %.1f/%i, p = %.2f"%(data_chisq,data_ndof,data_chi2Prob))
    addInfo.AddText("MC #chi^{2}/nDOF = %.1f/%i, p = %.2f"%(mc_chisq,mc_ndof,mc_chi2Prob))
    lFrame.Draw()
    addInfo.Draw() 
    lLegend.Draw()
    lCan.Modified() 
    lCan.Update()
    lCan.SaveAs(fPlotDir+iLabel+".pdf")
    lCan.SaveAs(fPlotDir+iLabel+".png")

def getChi2NDOF(iData,iModel,iRooFitResult,iVar):
    lDataHist = iData.binnedClone(iData.GetName()+"_binnedClone",iData.GetName()+"_binnedClone");
    pChi2     = iModel[0].createChi2(lDataHist,r.RooFit.Range(100.,250.),r.RooFit.Extended(r.kTRUE),r.RooFit.DataError(r.RooAbsData.Poisson))
    #print '#chi sq',pChi2.getVal()
    #print 'number of bins', iVar.getBins()
    #print 'number of free parameters', iRooFitResult.floatParsFinal().getSize()
    #return pChi2.getVal()/(int(iVar.getBins() - iRooFitResult.floatParsFinal().getSize())); 
    return pChi2.getVal(), int(iVar.getBins() - iRooFitResult.floatParsFinal().getSize())

def getPull(iVar,iPlot):
    lHPull = iPlot.pullHist();
    pPull  = iVar.frame(r.RooFit.Title("Pull"), r.RooFit.Bins(int(iVar.getBins())));
    lLine  = r.TLine(iVar.getMin(),0,iVar.getMax(),0);
    lLine.SetLineWidth(2); lLine.SetLineColor(r.kRed);
    pPull.addObject(lLine);
    pPull.addPlotable(lHPull,"P");
    pPull.SetTitle("");
    pPull.GetXaxis().SetTitle("");
    pPull.GetYaxis().SetTitle("pull");
    pPull.GetYaxis().SetRangeUser(-4,4);
    pPull.GetXaxis().SetTitleSize(0.10);
    pPull.GetXaxis().SetLabelSize(0.10);
    pPull.GetYaxis().SetTitleSize(0.10);
    pPull.GetYaxis().SetLabelSize(0.10);
    return pPull

# add gaussian constraint                                                                                                   
def addConstraint(iWorkspace,iVar,iMean,iSigma,iList):
    print '---- Adding gaussian constraint to %s'%iVar.GetName()
    lMean = r.RooRealVar("%s_mean"%iVar.GetName(),"%s_mean"%iVar.GetName(),iMean);
    lSigma = r.RooRealVar("%s_sigma"%iVar.GetName(),"%s_sigma"%iVar.GetName(),iSigma);
    lConstraint = r.RooGaussian("constraint_pdf_%s"%iVar.GetName(),"constraint_pdf_%s"%iVar.GetName(),iVar,lMean,lSigma)
    lConstraint.Print()
    iList.append(lConstraint.GetName())
    getattr(iWorkspace,"import")(lConstraint,r.RooFit.RecycleConflictNodes())


# make Pdf from model, probably should put parameters in dictionary
def makePdf(iWorkspace,iLabel,iModel,iMc=False):
    print '---- Making pdf' 
    lVar = iWorkspace.var(fVar);
    lModelPdf = None
    lTag = "%s_%s"%(iModel,iLabel)
    pVarHigh       = r.RooRealVar("%sHi_%s"%(lVar.GetName(),lTag),"%sHi_%s"%(lVar.GetName(),lTag),0.5,0.,1.);
    pVarHigh1      = r.RooRealVar("%sHi1_%s"%(lVar.GetName(),lTag),"%sHi1_%s"%(lVar.GetName(),lTag),0.5,0.,1.);
    lConstraints = []

    if iModel == "Exp_mc":
        lC_Exp_mc        = r.RooRealVar("c_mc"+lTag,"c_mc"+lTag,-0.03,-2.,0.05)
        lModelPdf = r.RooExponential("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_Exp_mc);
        
    if iModel == "Exp_data":
        lC_Exp_data      = r.RooRealVar("c_data"+lTag,"c_data"+lTag,-0.03,-2.,0.05)
        lModelPdf = r.RooExponential("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_Exp_data);
        
    if "PolyGaus_st" in iModel:
        lMean_Gaus       = r.RooRealVar("mean"+lTag,"mean"+lTag,175.,120.,190.)
        lSigma_Gaus      = r.RooRealVar("sigma"+lTag,"sigma"+lTag,15.,0.,100.)
        pGaus            = r.RooGaussian("gaus"+lTag,"gaus"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        lMean2_Gaus      = r.RooRealVar("mean2"+lTag,"mean2"+lTag,150.,140.,175.)
        lSigma2_Gaus     = r.RooRealVar("sigma2"+lTag,"sigma2"+lTag,80.,40.,100.)
        pGaus2           = r.RooGaussian("gaus2"+lTag,"gaus2"+lTag,lVar,lMean2_Gaus,lSigma2_Gaus)
        la1              = r.RooRealVar("a1"+lTag,"a1"+lTag,-0.003,-0.005,0.)
        #la1             = r.RooRealVar("a1"+lTag,"a1"+lTag,-0.04,-.1,0.)
        #la2             = r.RooRealVar("a2"+lTag,"a2"+lTag,0.1,0.,1.)
        #pPoly           = r.RooExponential("poly"+lTag,"poly"+lTag,lVar,la1)
        pPoly            = r.RooPolynomial("poly"+lTag,"poly"+lTag,lVar,r.RooArgList(la1),1)
        lModelPdf = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus,pPoly),r.RooArgList(pVarHigh))

    if "ErfExpGaus_st" in iModel:
        lC_ErfExp        = r.RooRealVar("c"+lTag,"c"+lTag,-0.05,-0.2,0.0)
        lOffSet_ErfExp   = r.RooRealVar("offset"+lTag,"offset"+lTag,175.,140.,300.)
        lWidth_ErfExp    = r.RooRealVar("width"+lTag,"width"+lTag,30.,10.,300.)
        lMean_Gaus       = r.RooRealVar("mean"+lTag,"mean"+lTag,175.,140.,200.)
        lSigma_Gaus      = r.RooRealVar("sigma"+lTag,"sigma"+lTag,7.,0.,40.)
        pGaus            = r.RooGaussian("gaus"+lTag,"gaus"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        #add W 
        #lMean_Gaus2     = r.RooRealVar("mean2"+lTag,"mean2"+lTag,85.,75.,95.)
        #lSigma_Gaus2    = r.RooRealVar("sigma2"+lTag,"sigma2"+lTag,10.,2.,20.)
        #pGaus2          = r.RooGaussian("gaus2"+lTag,"gaus2"+lTag,lVar,lMean_Gaus2,lSigma_Gaus2)
        
        pErfExp          = r.RooErfExpPdf("erfExp"+lTag,"erfExp"+lTag,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        #lModelPdf       = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pErfExp,pGaus,pGaus2),r.RooArgList(pVarHigh,pVarHigh1));
        lModelPdf = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pErfExp,pGaus),r.RooArgList(pVarHigh));

    if iModel == "ErfExp_st_mc":
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-5,1.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,175.,50.,1000.)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,50.,20.,1000.)
        lModelPdf = r.RooErfExpPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);

    if iModel == "ErfExp_st_data":
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-5,1.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,175.,50.,1000.)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,50.,20.,1000.)
        lModelPdf = r.RooErfExpPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);

    if "GausExp_tt" in iModel:
        if iLabel == "tt_mc" or iLabel == "tt_data":
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,175.,130.,200.)
        elif "fakeW" in iLabel:
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,130.,100.,150.)
        elif "realW" in iLabel:
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,175.,160.,200.)
        lSigma_Gaus        = r.RooRealVar("sigma_"+lTag,"sigma_"+lTag,7.,0.,40.)
        lC_Exp             = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-1.,1.)
        #la_Poly           = r.RooRealVar("a_"+lTag,"a_"+lTag,-0.005,0.,0.005)
        #la2_Poly           = r.RooRealVar("a2_"+lTag,"a2_"+lTag,-0.005,0.,0.005)
        pGaus              = r.RooGaussian("gaus_"+lTag,"gaus_%s"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pExp               = r.RooExponential("exp_"+lTag,"exp_%s"+lTag,lVar,lC_Exp)
        #pPoly              = r.RooPolynomial("poly_"+lTag,"poly_%s"+lTag,lVar,r.RooArgList(la_Poly,la2_Poly),1)
        lModelPdf          = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus,pExp),r.RooArgList(pVarHigh),1)
        #if iLabel == 'tt_mc' or iLabel == 'tt_data' or 'realW' in iLabel:
        #    addConstraint(iWorkspace,lMean_Gaus,175,8,lConstraints)
        #    addConstraint(iWorkspace,lSigma_Gaus,17,5,lConstraints)
        #elif ('tt_data' in iLabel or 'tt_mc' in iLabel) and 'fakeW' in iLabel:
        #    addConstraint(iWorkspace,lMean_Gaus,130,8,lConstraints)
        #    addConstraint(iWorkspace,lSigma_Gaus,17,6,lConstraints)

    if  "GausErfExp_tt" in iModel:
        ptbin = [int(s) for s in re.findall('\d+',iLabel)]
        print ptbin
        print iWorkspace.var("offset_mc")
        print iWorkspace.var("slope_mc")
        print iLabel
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-10,10.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,20.,0.,300.)
        if iLabel == "tt_mc" or iLabel == "tt_data":
            #lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,20.,0.,200.)
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,175.,130.,200.)
        elif "fakeW" in iLabel:  
            #lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,100.,0.,150.)
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,130.,100.,150.)
        elif "realW" in iLabel:
            #lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,200.,150.,250.)
            #lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,175.,160.,200.)
            #if iLabel.startswith("mc"):
            pt = (ptbin[0]+ptbin[1])/2
            if ptbin[1] == 100000:
                pt = 575
            elif ptbin[1] == 550:
                pt = 500
            if "mc" in iLabel: 
                lMean_Gaus     = r.RooFormulaVar("mean_"+lTag,"mean_"+lTag,"@0 + @1*@2", r.RooArgList(iWorkspace.var("offset_mc"),iWorkspace.var("slope_mc"),r.RooFit.RooConst(pt)))
            elif "data" in iLabel:
                lMean_Gaus     = r.RooFormulaVar("mean_"+lTag,"mean_"+lTag,"@0 + @1*@2", r.RooArgList(iWorkspace.var("offset_data"),iWorkspace.var("slope_data"),r.RooFit.RooConst(pt)))
            #elif iLabel.startswith("data"):
            #    lMean_Gaus     = r.RooFormulaVar("mean_"+lTag,"mean_"+lTag,"@0 + @1*@2", r.RooArgList(iWorkspace.var("offset_data"),iWorkspace.var("slope_data"),r.RooFit.RooConst((ptbin[0]+ptbin[1])/2)))

        a1             = r.RooRealVar("a1_"+lTag,"a1_"+lTag,0.0001,-0.01,0.01)
        a2             = r.RooRealVar("a2_"+lTag,"a2_"+lTag,0.00001,-0.0001,0.01)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,20.,0.,60.) #WAS[ 0,50]
        lSigma_Gaus    = r.RooRealVar("sigma_"+lTag,"sigma_"+lTag,7.,0.,60.)
        pPoly          = r.RooPolynomial("poly_"+lTag,"poly_%s"+lTag,lVar,r.RooArgList(a1,a2),1)
        pErfExp        = r.RooErfExpPdf("erfExp_"+lTag,"erfExp_%s"+lTag,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        pGaus          = r.RooGaussian("gaus_"+lTag,"gaus_%s"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pVarFrac       = r.RooRealVar("%sFrac_%s"%(lVar.GetName(),iLabel),"%sFrac_%s"%(lVar.GetName(),iLabel),1.0); # is this needed?
        #lC_ErfExp.setConstant(r.kTRUE);
        #lOffSet_ErfExp.setConstant(r.kTRUE);
        #lWidth_ErfExp.setConstant(r.kTRUE);
        pPoly          = r.RooPolynomial("poly_"+lTag,"poly_%s"+lTag,lVar,r.RooArgList(a1),1)
        pGaus          = r.RooGaussian("gaus_"+lTag,"gaus_%s"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pErfExp        = r.RooErfExpPdf("erfExp_%s"%iLabel,"erfExp_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        lModelPdf = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus,pErfExp,pPoly),r.RooArgList(pVarHigh,pVarHigh),1);
        #if iLabel == 'tt_mc' or iLabel == 'tt_data' or 'realW' in iLabel:
            #addConstraint(iWorkspace,lMean_Gaus,178,4,lConstraints)
            #addConstraint(iWorkspace,lSigma_Gaus,19,3,lConstraints)
            #addConstraint(iWorkspace,lWidth_ErfExp,50,4,lConstraints)
            #addConstraint(iWorkspace,lOffSet_ErfExp,164,4,lConstraints)
            #addConstraint(iWorkspace,lC_ErfExp,-1.6e-2,1e-1,lConstraints)

            #addConstraint(iWorkspace,lC_ErfExp,1e-4.,5e-4,lConstraints)
        #elif ('tt_data' in iLabel or 'tt_mc' in iLabel) and 'fakeW' in iLabel:
        #    addConstraint(iWorkspace,lMean_Gaus,110,8,lConstraints)
        #    addConstraint(iWorkspace,lSigma_Gaus,24,6,lConstraints)
        #    addConstraint(iWorkspace,lWidth_ErfExp,13,5,lConstraints)
        #    addConstraint(iWorkspace,lOffSet_ErfExp,128,8,lConstraints)
        #    addConstraint(iWorkspace,lC_ErfExp,-7e-3,2e-2,lConstraints)

    if iModel == "ErfExp_tt_mc":
        lModelPdf      = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus,pErfExp),r.RooArgList(pVarHigh),1);
        addConstraint(iWorkspace,lMean_Gaus,175,3,lConstraints)
        addConstraint(iWorkspace,lSigma_Gaus,7,2,lConstraints)
        if 'realW' in iLabel:
            addConstraint(iWorkspace,lOffSet_ErfExp,175,10,lConstraints)
            addConstraint(iWorkspace,lWidth_ErfExp,30,5,lConstraints) 

    if "GausGaus_tt" in iModel:
        #narrow component
        lMean_Gaus1     = r.RooRealVar("mean1"+lTag,"mean1"+lTag,175.,140.,200.)
        lSigma_Gaus1    = r.RooRealVar("sigma1"+lTag,"sigma1"+lTag,7.,0.,40.)
        pGaus1          = r.RooGaussian("gaus1"+lTag,"gaus1"+lTag,lVar,lMean_Gaus1,lSigma_Gaus1)
        pVarFrac1       = r.RooRealVar("%sFrac_%s"%(lVar.GetName(),iLabel),"%sFrac_%s"%(lVar.GetName(),iLabel),1.0); # is this needed?
        pGaus1          = r.RooGaussian("gaus1"+lTag,"gaus1"+lTag,lVar,lMean_Gaus1,lSigma_Gaus1)
        #broad component
        lMean_Gaus2     = r.RooRealVar("mean2"+lTag,"mean2"+lTag,120.,100.,175.) #offset left
        lSigma_Gaus2    = r.RooRealVar("sigma2"+lTag,"sigma2"+lTag,60.,0.,80.) #wider
        pGaus2          = r.RooGaussian("gaus2"+lTag,"gaus2"+lTag,lVar,lMean_Gaus2,lSigma_Gaus2)
        #W
        lMean_Gaus3     = r.RooRealVar("mean3"+lTag,"mean3"+lTag,80.,75.,85.)
        lSigma_Gaus3    = r.RooRealVar("sigma3"+lTag,"sigma3"+lTag,10.,0.,20.)
        pGaus3          = r.RooGaussian("gaus3"+lTag,"gaus3"+lTag,lVar,lMean_Gaus3,lSigma_Gaus3)
        
        #try poly
        a1              = r.RooRealVar("a1"+lTag,"a1"+lTag,0.004,0.,0.007)
        a2              = r.RooRealVar("a2"+lTag,"a2"+lTag,-1e-6,-1e-9,0.)
        a3              = r.RooRealVar("a3"+lTag,"a3"+lTag,1e-7,0.,0.005)
        lExp2             = r.RooPolynomial("lExp2"+lTag,"lExp2"+lTag,lVar,r.RooArgList(a1,a2,a3),1)
        
        #a3              = r.RooRealVar("a3"+lTag,"a3"+lTag,1e-8,-0.005,0.005)
        #lPoly           = r.RooPolynomial("lPoly"+lTag,"lPoly"+lTag,lVar, r.RooArgList(a1,a2),1)
        
        #lModelPdf      = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus1,pGaus2),r.RooArgList(pVarFrac1),1);
        #lModelPdf1      = r.RooAddPdf("model_pdf1_%s"%iLabel,"model_pdf1_%s"%iLabel,r.RooArgList(pGaus1,pGaus2),r.RooArgList(pVarHigh));
        addConstraint(iWorkspace,lMean_Gaus1,175,3,lConstraints)
        addConstraint(iWorkspace,lSigma_Gaus1,7,2,lConstraints)
        lModelPdf       = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus1,lExp2),r.RooArgList(pVarHigh))

    getattr(iWorkspace,"import")(lModelPdf,r.RooFit.RecycleConflictNodes())
    return iWorkspace.pdf("model_pdf_%s"%iLabel),lConstraints

# return RooExtendPdf and Constraints
def makeModel(iWorkspace,iLabel,iModel):
    print '---- Making model'
    lNumber = r.RooRealVar("number_%s"%iLabel,"number_%s"%iLabel,500,0.,1e7); # check if needed
    lModelPdf,lConstraints = makePdf(iWorkspace,iLabel,iModel)
    lModel = r.RooExtendPdf("model_%s"%iLabel,"model_%s"%iLabel,lModelPdf,lNumber)
    getattr(iWorkspace,"import")(lNumber,r.RooFit.RecycleConflictNodes())
    getattr(iWorkspace,"import")(lModel,r.RooFit.RecycleConflictNodes())
    iWorkspace.pdf("model_%s"%iLabel).Print()
    return iWorkspace.pdf("model_%s"%iLabel),lConstraints

# fit single rooDataset to model
# return list of constraints
def fitSingleMC(iWorkspace,iLabel,iModel):
    print '---- Fitting for %s'%iLabel
    lVar   = iWorkspace.var(fVar);
    lData  = iWorkspace.data(iLabel+"_D");
    lModel,lConstraints = makeModel(iWorkspace,iLabel,iModel)

    #pRooFitResult = lModel.fitTo(lData,r.RooFit.Save(1),                             r.RooFit.Extended(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.Verbose(r.kFALSE))
    pRooFitResult = lModel.fitTo(lData,r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Extended(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.Verbose(r.kFALSE));
    #pRooFitResult = lModel.fitTo(lData,r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Extended(r.kTRUE),r.RooFit.Minimizer("Minuit2"));
    getattr(iWorkspace,"import")(lModel,r.RooFit.RecycleConflictNodes())
    pRooFitResult.Print()
    getattr(iWorkspace,"import")(pRooFitResult)
    
    draw(lVar,lData,[lModel],pRooFitResult,iLabel+"_only")

    return lConstraints

class TopPeak():
    def __init__(self,options):

        self._lMSD_lo = fXMin
        self._lMSD_hi = fXMax
        self._lMSD_blind_lo = 140;
        self._lMSD_blind_hi = 200;

        self._lMatch  = options.match

        self._lW      = r.RooWorkspace("w","w")

        self._lMSD    = r.RooRealVar(fVar,fVar,self._lMSD_lo,self._lMSD_hi)
        self._lMSD.setRange('Low',self._lMSD_lo,self._lMSD_blind_lo)
        self._lMSD.setRange('Blind',self._lMSD_blind_lo,self._lMSD_blind_hi)
        self._lMSD.setRange('High',self._lMSD_blind_hi,self._lMSD_hi)
        self._lMSD.setBins(fNBins)
        getattr(self._lW,"import")(self._lMSD)

        Pt_cats = r.RooCategory("Pt_cats","Pt_cats")
        for Pt_it in range(0,len(fPt_bins)-1):
            Pt_cats.defineType("pT_" + str(fPt_bins[Pt_it]) + "_" +str(fPt_bins[Pt_it+1]))

        self._lWeight = r.RooRealVar(fWeight,fWeight,-1e+25,1e+25)
        self._lPt     = r.RooRealVar(fPt,fPt,0,2000)

        self._lPDatas = {}
        self._lHPdfs  = {}
        self._lNPdfs  = {}
        self._lModels = {}
        self._lConstraints = {}
        self._Scale_val = {}
        self._Scale_unc = {}
        self._Resolution = {}
        self._combData_data = {}
        self._combData_MC = {}

        Slope_mc  = r.RooRealVar("slope_mc","slope_mc",0.0,-0.1,0.1)
        Offset_mc = r.RooRealVar("offset_mc","offset_mc",175.0,165.,200.)
        getattr(self._lW,"import")(Slope_mc,r.RooFit.RecycleConflictNodes())
        getattr(self._lW,"import")(Offset_mc,r.RooFit.RecycleConflictNodes())

        Slope_data  = r.RooRealVar("slope_data","slope_data",0.0,-0.1,0.1)
        Offset_data = r.RooRealVar("offset_data","offset_data",175.,165.,200.)
        getattr(self._lW,"import")(Slope_data,r.RooFit.RecycleConflictNodes())
        getattr(self._lW,"import")(Offset_data,r.RooFit.RecycleConflictNodes())

        for iLabel,iFile in fFiles.iteritems():
          for Pt_it in range(len(fPt_bins)-1): 

            print 'preparing sample %s in pT bin[%s,%s]'%(iLabel,fPt_bins[Pt_it],fPt_bins[Pt_it+1])
            pFile = fDir+fTag+"/"+iFile+"_"+fTag+".root"
            fCut = "(Puppijet0_pt>" + str(fPt_bins[Pt_it]) + " && Puppijet0_pt<" + str(fPt_bins[Pt_it+1]) + ")"
            lLabel = iLabel + '_pT_' + str(fPt_bins[Pt_it]) + '_' + str(fPt_bins[Pt_it+1])
            Pt_label = "pT_" + str(fPt_bins[Pt_it]) + '_' + str(fPt_bins[Pt_it+1])
            self._lPDatas[lLabel] = self.getDataset(iLabel,lLabel,Pt_label,pFile,fCut,True)
            if iLabel != "mc" and iLabel != "data":
                self._lConstraints[lLabel] = fitSingleMC(self._lW,lLabel,fPDFs[iLabel])

            # check why are these needed
            self._lHPdfs[lLabel] = self.histPdf(lLabel)
            self._lNPdfs[lLabel] = r.RooRealVar("number_%s"%lLabel,"number_%s"%lLabel,self._lW.data(lLabel+"_D").sumEntries()) # *sf? ttsf, what is that?
            self._lNPdfs[lLabel].Print()
            print self._lW.data(lLabel+"_D").sumEntries()
        
        self._lW.Print()

    def getDataset(self,iLabel,lLabel,Pt_label,iFile,iCut="(1==1)",iWeight=False,iTree="otree2"):

        print 'preparing dataset for ',iFile,' with cut ',iCut
        lFile   = r.TFile(iFile)
        lTree   = lFile.Get(iTree)
        #lData   = r.RooDataSet(iLabel+"_D",iLabel+"_D",lTree,r.RooAet(self._lMSD,self._lPt,self._lWeight),iCut)
        lData = r.RooDataSet(lLabel+"_D",lLabel+"_D",r.RooArgSet(self._lMSD,self._lWeight),r.RooFit.WeightVar(self._lWeight))
        lCut = r.TTreeFormula("lCut",iCut,lTree)  
        lTree.SetNotify(lCut)       
        for i0 in range(lTree.GetEntriesFast()):
        #for i0 in range(100000): 
            lTree.LoadTree(i0) 
            lSel = False 
            for i1 in range(lCut.GetNdata()):  
                if (lCut.EvalInstance(i1)):  
                    lSel=True; break;
            if not lSel: continue     
            lTree.GetEntry(i0) 
            if iWeight:
                lWeight = getattr(lTree,fWeight)*fLumi
            else:
                lWeight = 1
            lMass = getattr(lTree,fVar)
            lMatched = 0
            jmatched = getattr(lTree,"Puppijet0_vMatching");
            jhadronic = getattr(lTree,"Puppijet0_isHadronicV");
            if jhadronic == 1.0 and jmatched < 1. and jmatched > 0.:
                lMatched = 1;
            if 'realW' in iLabel and lMatched == 0: continue
            if 'fakeW' in iLabel and lMatched == 1: continue
            if lMass < self._lMSD_hi and lMass > self._lMSD_lo:
                self._lMSD.setVal(lMass)
                lData.add(r.RooArgSet(self._lMSD), lWeight)
        getattr(self._lW,'import')(lData,r.RooFit.RecycleConflictNodes())
        lData.Print()
        lFile.Close()
        return self._lW.data(lLabel+"_D")
                     
    # shape and normalization parameters fixed
    def getModel(self, iLabel):
        print 'getting model for ', iLabel
        lModel    = self._lW.pdf("model_"+iLabel)
        lModel.Print()
        self._lW.data(iLabel+"_D").Print()
        # lPars     = lModel.getParameters(self._lW.data(iLabel+"_D"))
        # pParsIter = lPars.createIterator()
        # pParsIter.Reset()
        # pPar = pParsIter.Next()
        # while (pPar):
        #     pPar.setConstant(r.kTRUE)
        #     pPar = pParsIter.Next()
        return lModel

    # return roohistpdf
    def histPdf(self, iLabel):
        lDataSet = self._lW.data(iLabel+"_D")
        lVar = self._lW.var(fVar)
        lReduced  = lDataSet.reduce(r.RooArgSet(lVar))
        lDataHist = lReduced.binnedClone(lReduced.GetName()+"_binnedClone",lReduced.GetName()+"_binnedClone");
        lHPdf     = r.RooHistPdf(lReduced.GetName()+"P",lReduced.GetName()+"P",r.RooArgSet(lVar),lDataHist,0)
        print lReduced.sumEntries()
        getattr(self._lW,"import")(lHPdf,r.RooFit.RecycleConflictNodes())
        return lHPdf

    # fit mc and data
    def fit(self):

      for Pt_it in range(len(fPt_bins)-1):
        print 'processing combined fit for pt bin [%s,%s]'%(fPt_bins[Pt_it],fPt_bins[Pt_it+1])
        lVar   = self._lW.var(fVar);  
        print self._lPDatas
        lSData = self._lPDatas["data_pT_" + str(fPt_bins[Pt_it]) + "_" + str(fPt_bins[Pt_it+1])]
        lSMc   = self._lPDatas["mc_pT_"   + str(fPt_bins[Pt_it]) + "_" + str(fPt_bins[Pt_it+1])]

        # impose normalization constraint
        norm_fakeW = self._lW.arg("number_tt_fakeW_mc_pT_" + str(fPt_bins[Pt_it]) + '_' + str(fPt_bins[Pt_it+1])).getVal()
        norm_realW = self._lW.arg("number_tt_realW_mc_pT_" + str(fPt_bins[Pt_it]) + '_' + str(fPt_bins[Pt_it+1])).getVal()
        norm_st    = self._lW.arg("number_st_mc_pT_" + str(fPt_bins[Pt_it]) + '_' + str(fPt_bins[Pt_it+1])).getVal()
        norm_ratio = (norm_fakeW + norm_realW) / norm_st
        print 'norm ratio: %s'%norm_ratio

        # TODO define ratio && add constraint
        '''
        norm_constraint_mc = r.RooFormulaVar("norm_constraint_mc","(@0 + @1)/@2",r.RooArgList(self._lW.arg("number_tt_fakeW_mc"),self._lW.arg("number_tt_realW_mc"), self._lW.arg("number_st_mc")))
        norm_constraint_data = r.RooFormulaVar("norm_constraint_data","(@0 + @1)/@2",r.RooArgList(self._lW.arg("number_tt_fakeW_data"),self._lW.arg("number_tt_realW_data"), self._lW.arg("number_st_data")))
        print '---- Adding gaussian constraint to normalization'
        lMean_mc = r.RooRealVar("mean_mc","mean_mc",norm_ratio);
        lSigma_mc = r.RooRealVar("sigma_mc","sigma_mc",norm_ratio*.4);
        lMean_data = r.RooRealVar("mean_data","mean_data",norm_ratio);
        lSigma_data = r.RooRealVar("sigma_data","sigma_data",norm_ratio*.4);
        
        lConstraint_mc = r.RooGaussian("constraint_norm_mc","constraint_norm_mc",norm_constraint_mc,lMean_mc,lSigma_mc)
        lConstraint_data = r.RooGaussian("constraint_norm_data","constraint_norm_data",norm_constraint_data,lMean_data,lSigma_data)
        lConstraint_mc.Print()
        lConstraint_data.Print()
        '''
        # constraints
        pPdfConstraints_Mc   = r.RooArgSet("pPdfConstraints_Mc")
        pPdfConstraints_Data = r.RooArgSet("pPdfConstraints_Data")

        #pPdfConstraints_Mc.add(lConstraint_mc)
        #pPdfConstraints_Data.add(lConstraint_data)

        pSamples = ['st_mc','st_data','tt_mc','tt_data']
        if self._lMatch:
            pSamples = ['st_mc','st_data','tt_fakeW_mc','tt_fakeW_data','tt_realW_mc','tt_realW_data']
        for iLabel in pSamples:
            lLabel = iLabel + "_pT_" + str(fPt_bins[Pt_it]) + '_' + str(fPt_bins[Pt_it+1])
            self._lModels[lLabel] = self.getModel(lLabel)
            self._lModels[lLabel].Print()
            for i0 in range(len(self._lConstraints[lLabel])):
                if 'mc' in lLabel:
                    pPdfConstraints_Mc.add(self._lW.pdf(self._lConstraints[lLabel][i0]))
                if 'data' in lLabel:
                    pPdfConstraints_Data.add(self._lW.pdf(self._lConstraints[lLabel][i0]))

        args = self._lW.allVars()
        args.Print()
        iter = args.createIterator()
        var = iter.Next()
        while var:
            print var.GetName()
            #put floating variables here
            if "offset_mc" in var.GetName() or "slope_mc" in var.GetName() or "offset_data" in var.GetName() or "slope_data" in var.GetName() or "mean_GausErfExp_tt" in var.GetName() or "c_GausErfExp_tt_data_tt_realW_data" in var.GetName() or "offset_GausErfExp_tt_data_tt_realW_data" in var.GetName() or "width_GausErfExp_tt_data_tt_realW_data" in var.GetName()  or "a1_GausExp_tt_data" in var.GetName()  or "number_st" in var.GetName() or "sigma_GausErfExp_tt" in var.GetName():
                pass
            else:
                var.setConstant(r.kTRUE)
                print 'set constant'
            var = iter.Next()

        Pt_label = "_pT_" + str(fPt_bins[Pt_it]) + '_' + str(fPt_bins[Pt_it+1])

        st_norm_mc = self._lW.arg("number_st_mc"+Pt_label)
        norm_ratio_mc = r.RooRealVar("norm_ratio_mc"+Pt_label,"norm_ratio_mc"+Pt_label,norm_ratio) 
        tt_norm_mc = r.RooFormulaVar("tt_norm_mc"+Pt_label,"(@0*@1)",r.RooArgList(norm_ratio_mc,st_norm_mc))

        st_norm_data = self._lW.arg("number_st_data"+Pt_label)
        norm_ratio_data = r.RooRealVar("norm_ratio_data"+Pt_label,"norm_ratio_data"+Pt_label,norm_ratio)
        tt_norm_data = r.RooFormulaVar("tt_norm_data"+Pt_label,"(@0*@1)",r.RooArgList(norm_ratio_data,st_norm_data))

        print 'made normalizations'
        # models
        if self._lMatch:
            #self._lModels['Mc'] = r.RooAddPdf("model_mc","model_mc",r.RooArgList(self._lModels['st_mc'],self._lModels['tt_fakeW_mc'],self._lModels['tt_realW_mc']))
            #self._lModels['Data'] = r.RooAddPdf("model_data","model_data",r.RooArgList(self._lModels['st_data'],self._lModels['tt_fakeW_data'],self._lModels['tt_realW_data']))
            #st_model_mc     = r.RooExtendPdf("model_mc_st","model_mc_st",r.RooArgList(self._lModels['st_mc']))
            tt_model_mc     = r.RooAddPdf("model_mc_tt"+Pt_label,"model_mc_tt"+Pt_label,r.RooArgList(self._lW.pdf("model_tt_fakeW_mc"+Pt_label),self._lW.pdf("model_tt_realW_mc"+Pt_label)))
            tt_model_mc_ext = r.RooExtendPdf("ext_model_mc_tt"+Pt_label,"ext_model_mc_tt"+Pt_label,tt_model_mc,tt_norm_mc)
            st_model_mc_ext = r.RooExtendPdf("ext_model_mc_st"+Pt_label,"ext_model_mc_st"+Pt_label,self._lW.pdf("model_st_mc"+Pt_label),st_norm_mc)
            self._lModels["Mc"+Pt_label] = r.RooAddPdf("model_mc"+Pt_label,"model_mc"+Pt_label,r.RooArgList(tt_model_mc_ext,st_model_mc_ext))


            tt_model_data     = r.RooAddPdf("model_data_tt"+Pt_label,"model_data_tt"+Pt_label,r.RooArgList(self._lW.pdf("model_tt_fakeW_data"+Pt_label),self._lW.pdf("model_tt_realW_data"+Pt_label)))
            tt_model_data_ext = r.RooExtendPdf("ext_model_data_tt"+Pt_label,"ext_model_data_tt"+Pt_label,tt_model_data,tt_norm_data)
            st_model_data_ext = r.RooExtendPdf("ext_model_data_st"+Pt_label,"ext_model_data_st"+Pt_label,self._lW.pdf("model_st_data"+Pt_label),st_norm_data)
            self._lModels["Data"+Pt_label] = r.RooAddPdf("model_data"+Pt_label,"model_data"+Pt_label,r.RooArgList(tt_model_data_ext,st_model_data_ext))

        else:
            self._lModels['Mc'] = r.RooAddPdf("model_mc","model_mc",r.RooArgList(self._lModels['st_mc'],self._lModels['tt_mc']))
            self._lModels['Data'] = r.RooAddPdf("model_data","model_data",r.RooArgList(self._lModels['st_data'],self._lModels['tt_data']))
        

        # fit to mc
        print self._lModels
        pRooFitResult_Mc   = self._lModels["Mc"+Pt_label].fitTo(lSMc,r.RooFit.Strategy(1),r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.Verbose(r.kFALSE))
        #pRooFitResult_Mc   = self._lModels['Mc'].fitTo(lSMc,r.RooFit.Strategy(1),r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.ExternalConstraints(pPdfConstraints_Mc))

        getattr(self._lW,"import")(self._lModels["Mc"+Pt_label],r.RooFit.RecycleConflictNodes())
         
        pRooFitResult_Mc.Print()
        draw(lVar,lSMc,[self._lModels["Mc"+Pt_label]],pRooFitResult_Mc,"mc_only"+Pt_label)
 
        # print mc parameters
        x1 = self._lModels["Mc"+Pt_label].getVariables()
        print x1.Print("v")

        # fit to data
        pRooFitResult_Data = self._lModels["Data"+Pt_label].fitTo(lSData,r.RooFit.Strategy(1),r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.Verbose(r.kFALSE))
        #pRooFitResult_Data = self._lModels['Data'].fitTo(lSData,r.RooFit.Strategy(2),r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.ExternalConstraints(pPdfConstraints_Data),r.RooFit.Verbose(r.kFALSE))
        #pRooFitResult_Data = self._lModels['Data'].fitTo(lSData,r.RooFit.Strategy(1),r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.ExternalConstraints(pPdfConstraints_Data))

        getattr(self._lW,"import")(self._lModels["Data"+Pt_label],r.RooFit.RecycleConflictNodes())

        pRooFitResult_Data.Print()
        draw(lVar,lSData,[self._lModels["Data"+Pt_label]],pRooFitResult_Data,"data_only"+Pt_label) 

        # print data parameters
        x2 = self._lModels["Data"+Pt_label].getVariables()
        print x2.Print("v")
        params = {}
        if self._lMatch:
            lTagMc = fPDFs["tt_realW_mc"]+"_tt_realW_mc"+Pt_label
            lTagData = fPDFs["tt_realW_data"]+"_tt_realW_data"+Pt_label
        else:
            lTagMc = fPDFs['tt_mc']+'_tt_mc'
            lTagData = fPDFs['tt_data']+'_tt_data'
        """
        params["mc_mean_val"+Pt_label]    = x1.find("mean_%s"%lTagMc).getVal()
        params["mc_mean_err"+Pt_label]    = x1.find("mean_%s"%lTagMc).getError()
        params["mc_sigma_val"+Pt_label]   = x1.find("sigma_%s"%lTagMc).getVal() 
        params["mc_sigma_err"+Pt_label]   = x1.find("sigma_%s"%lTagMc).getError()

        params["data_mean_val"+Pt_label]    = x2.find("mean_%s"%lTagData).getVal()
        params["data_mean_err"+Pt_label]    = x2.find("mean_%s"%lTagData).getError()
        params["data_sigma_val"+Pt_label]   = x2.find("sigma_%s"%lTagData).getVal()
        params["data_sigma_err"+Pt_label]   = x2.find("sigma_%s"%lTagData).getError()
        print params
        #self._Scale["scale"+Pt_label] = r.RooFormulaVar("fScale%s","fScale%s","@0/@1",r.RooArgList(x2.find("mean_%s"%lTagData),x1.find("mean_%s"%lTagMc)))
        self._Scale_val["scale_val"+Pt_label] = x2.find("mean_%s"%lTagData).getVal() / x1.find("mean_%s"%lTagMc).getVal()
        self._Scale_unc["scale_unc"+Pt_label] = self._Scale_val["scale_val"+Pt_label] * np.sqrt(x2.find("mean_%s"%lTagData).getError()**2/x2.find("mean_%s"%lTagData).getVal()**2 + x1.find("mean_%s"%lTagMc).getError()**2/x1.find("mean_%s"%lTagMc).getVal()**2)

        print "scale += unc. =%.3f +/- %.3f"%(self._Scale_val["scale_val"+Pt_label], self._Scale_unc["scale_unc"+Pt_label])
        """
        # draw data and mc
        drawDataMc(lVar,lSData,[self._lModels["Data"+Pt_label]],lSMc,[self._lModels["Mc"+Pt_label]],pRooFitResult_Mc,pRooFitResult_Data,params,"data_mc"+Pt_label)
        # write workspace
        self._lW.writeToFile(fOutput)
   
      self._lW.Print()
      print self._Scale_val
      print self._Scale_unc
      print "DATASETS:", self._lPDatas
      print "MODELS:", self._lModels

      Pt_cats_mc   = r.RooCategory("Pt_cats_mc"  ,"Pt_cats_mc")
      Pt_cats_data = r.RooCategory("Pt_cats_data","Pt_cats_data")

      for Pt_it in range(0,len(fPt_bins)-1):
          Pt_label = "pT_" + str(fPt_bins[Pt_it]) + "_" +str(fPt_bins[Pt_it+1])
          Pt_cats_mc  .defineType("mc_"  +Pt_label,Pt_it)
          Pt_cats_data.defineType("data_"+Pt_label,Pt_it)

      #Pt_cats_mc.Print("v")
      #Pt_cats_data.Print("v")

      combData_data = r.RooDataSet("combData_data","combData_data",r.RooArgSet(self._lMSD,self._lWeight),r.RooFit.WeightVar(self._lWeight),r.RooFit.Index(Pt_cats_data),r.RooFit.Import("data_pT_350_400",self._lPDatas["data_pT_350_400"]),r.RooFit.Import("data_pT_400_475",self._lPDatas["data_pT_400_475"]), r.RooFit.Import("data_pT_475_550",self._lPDatas["data_pT_475_550"]),r.RooFit.Import("data_pT_550_100000",self._lPDatas["data_pT_550_100000"]))
      combData_mc = r.RooDataSet("combData_mc","combData_mc",r.RooArgSet(self._lMSD,self._lWeight),r.RooFit.Index(Pt_cats_mc),r.RooFit.Import("mc_pT_350_400",self._lPDatas["mc_pT_350_400"]),r.RooFit.Import("mc_pT_400_475",self._lPDatas["mc_pT_400_475"]), r.RooFit.Import("mc_pT_475_550",self._lPDatas["mc_pT_475_550"]),r.RooFit.Import("mc_pT_550_100000",self._lPDatas["mc_pT_550_100000"]))

      combData_data.Print("v")
      combData_mc.Print("v")

      simPdf_data = r.RooSimultaneous("simPdf_data","simPdf_data",Pt_cats_data)
      simPdf_mc   = r.RooSimultaneous("simPdf_mc"  ,"simPdf_mc"  ,Pt_cats_mc)

      for Pt_it in range(0,len(fPt_bins)-1):
         Pt_label = "pT_" + str(fPt_bins[Pt_it]) + "_" +str(fPt_bins[Pt_it+1])
         simPdf_data.addPdf(self._lW.pdf("model_data_"+Pt_label),"data_"+Pt_label)
         simPdf_mc  .addPdf(self._lW.pdf("model_mc_"  +Pt_label),"mc_"  +Pt_label)
      """
      offset_mc = r.RooRealVar("offset_mc","offset_mc",0.5,0.50,1.10)
      slope_mc  = r.RooRealVar("slope_mc" ,"slope_mc" ,-0.5,-0.10,0.10)
      scale_pt_mc_cat0 = r.RooFormulaVar("scale_pt_mc_cat0","scale_pt_mc_cat0","@0 = @1 + @2*@3",r.RooArgList(self._lW.arg("mean_GausErfExp_tt_mc_tt_realW_mc_pT_350_400"),offset_mc,slope_mc,r.RooFit.RooConst(375)))
      scale_pt_mc_cat1 = r.RooFormulaVar("scale_pt_mc_cat1","scale_pt_mc_cat1","@0 = @1 + @2*@3",r.RooArgList(self._lW.arg("mean_GausErfExp_tt_mc_tt_realW_mc_pT_400_475"),offset_mc,slope_mc,r.RooFit.RooConst(437)))
      scale_pt_mc_cat2 = r.RooFormulaVar("scale_pt_mc_cat2","scale_pt_mc_cat2","@0 = @1 + @2*@3",r.RooArgList(self._lW.arg("mean_GausErfExp_tt_mc_tt_realW_mc_pT_475_550"),offset_mc,slope_mc,r.RooFit.RooConst(512)))
      scale_pt_mc_cat3 = r.RooFormulaVar("scale_pt_mc_cat3","scale_pt_mc_cat3","@0 = @1 + @2*@3",r.RooArgList(self._lW.arg("mean_GausErfExp_tt_mc_tt_realW_mc_pT_550_100000"),offset_mc,slope_mc,r.RooFit.RooConst(600)))

      offset_data = r.RooRealVar("offset_data","offset_data",0.5,0.90,1.10)
      slope_data  = r.RooRealVar("slope_data" ,"slope_data" ,-0.05,-0.10,0.10)
      scale_pt_data_cat0 = r.RooFormulaVar("scale_pt_data_cat0","scale_pt_data_cat0","@0 = @1 + @2*@3",r.RooArgList(self._lW.arg("mean_GausErfExp_tt_data_tt_realW_data_pT_350_400"),offset_data,slope_data,r.RooFit.RooConst(375)))
      scale_pt_data_cat1 = r.RooFormulaVar("scale_pt_data_cat1","scale_pt_data_cat1","@0 = @1 + @2*@3",r.RooArgList(self._lW.arg("mean_GausErfExp_tt_data_tt_realW_data_pT_400_475"),offset_data,slope_data,r.RooFit.RooConst(437)))
      scale_pt_data_cat2 = r.RooFormulaVar("scale_pt_data_cat2","scale_pt_data_cat2","@0 = @1 + @2*@3",r.RooArgList(self._lW.arg("mean_GausErfExp_tt_data_tt_realW_data_pT_475_550"),offset_data,slope_data,r.RooFit.RooConst(512)))
      scale_pt_data_cat3 = r.RooFormulaVar("scale_pt_data_cat3","scale_pt_data_cat3","@0 = @1 + @2*@3",r.RooArgList(self._lW.arg("mean_GausErfExp_tt_data_tt_realW_data_pT_550_100000"),offset_data,slope_data,r.RooFit.RooConst(600)))
      """

      print 'run simultaneous fit...'
      simFit_mc   = simPdf_mc  .fitTo(combData_mc,r.RooFit.Save(r.kTRUE),r.RooFit.Verbose(r.kFALSE),r.RooFit.Minimizer("Minuit2"),r.RooFit.SumW2Error(r.kTRUE))
      simFit_mc.Print()
      simFit_data = simPdf_data.fitTo(combData_data,r.RooFit.Save(r.kTRUE),r.RooFit.Verbose(r.kFALSE),r.RooFit.Minimizer("Minuit2"),r.RooFit.SumW2Error(r.kTRUE))
      simFit_data.Print()

      par = {}
      par["slope_mc"]    = self._lW.var("slope_mc") .getValV()
      par["offset_mc"]   = self._lW.var("offset_mc").getValV()
      par["slope_data"]  = self._lW.var("slope_data") .getValV()
      par["offset_data"] = self._lW.var("offset_data").getValV()

      par["slope_mc_e"]  = self._lW.var("slope_mc") .getError()
      par["offset_mc_e"]  = self._lW.var("offset_mc") .getError()
      par["slope_data_e"]  = self._lW.var("slope_data") .getError()
      par["offset_data_e"] = self._lW.var("offset_data").getError()


      fitSeparate(par)

if __name__ == "__main__":
    options = parser()
    print options
    global fTag,fXMin,fXMax,fNBins
    fTag   = options.tag
    fXMin  = options.xmin
    fXMax  = options.xmax
    fNBins = int( (fXMax - fXMin) / fBinWidth)
    # get roodataset and make individual fits
    lTop = TopPeak(options);
    # make combined fit
    lTop.fit();
