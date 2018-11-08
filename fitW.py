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
#fPt_bins = [200,400,475,550,100000]
fPt_bins = [200,100000]
fCats = ['pass','fail']
fWCut = {'pass': "Puppijet0_N2DDT<0",
         'fail': "Puppijet0_N2DDT>=0"}
fWeight="weight"
fLumi=41.4
fBinWidth =5
fDir="/afs/cern.ch/user/c/cmantill/work/public/forJeff/sklimWtag"
fPlotDir="plotsWtag/"
fOutput="w.root"

fFiles = {'data': 'Data',
          'wlnu_mc' : 'WJets',
          'wlnu_data': 'WJets',
          'st_mc':   'ST',
          'st_data':   'ST',
          'tt_fakeW_mc': 'TT',
          'tt_realW_mc': 'TT',
          'tt_fakeW_data': 'TT',
          'tt_realW_data': 'TT',
          "tt_signal_mc": 'TT',
          "tt_signal_data": 'TT',
          "tt_bkg_mc": 'TT',
          "tt_bkg_data": 'TT',
          'mc': 'PseudoData'
          }

fPDFs = {
         'wlnu_mc':       'Exp_mc',
         'wlnu_data':     'Exp_data',
         'st_mc':         'ErfExpGaus_st_mc',
         'st_data':       'ErfExpGaus_st_data',
         'tt_fakeW_mc':   'ErfExp_tt_mc',
         'tt_realW_mc':   'GausErfExp_tt_mc',
         'tt_fakeW_data': 'ErfExp_tt_data',
         'tt_realW_data': 'GausErfExp_tt_data'}

fPDFs['tt_signal_mc'] = fPDFs['tt_realW_mc']
fPDFs['tt_signal_data'] = fPDFs['tt_realW_data']
fPDFs['tt_bkg_mc'] = fPDFs['tt_fakeW_mc']
fPDFs['tt_bkg_data'] = fPDFs['tt_fakeW_data']

r.gSystem.Load("./PDFs/HWWLVJRooPdfs_cxx.so")
r.gSystem.Load("./PDFs/PdfDiagonalizer_cc.so")
r.gSystem.Load("./PlotStyle/Util_cxx.so")
r.gSystem.Load("./PlotStyle/PlotUtils_cxx.so")

def parser():
    parser = OptionParser()
    parser.add_option('--tag',     action='store', type='string', dest='tag',   default='AK8v42017',      help='samples tag')
    parser.add_option('--xMin',    action='store', type='float',  dest='xmin',  default=50,                help='x-min')
    parser.add_option('--xMax',    action='store', type='float',  dest='xmax',  default=130,               help='x-max')
    parser.add_option('-b',        action='store_true',           dest='noX',   default=True,              help='no X11 windows')
    (options,args) = parser.parse_args()
    return options

def drawFrame(iFrame,iData,iFuncs,iLegend,iColor,isPoisson=False):
    if isPoisson:
        iData.plotOn(iFrame,r.RooFit.DataError(r.RooAbsData.Poisson),r.RooFit.XErrorSize(0))
    else:
        iData.plotOn(iFrame)
    for pFunc in iFuncs:
        pFunc.plotOn(iFrame,r.RooFit.LineColor(iColor),r.RooFit.LineStyle(iColor != 50+1))
        pFunc.plotOn(iFrame,r.RooFit.Name("ErfExp"),r.RooFit.Components("erfExp*"),r.RooFit.LineColor(4))
        pFunc.plotOn(iFrame,r.RooFit.Name("Gaus"),r.RooFit.Components("gaus*"),r.RooFit.LineColor(3))
        pFunc.plotOn(iFrame,r.RooFit.Name("Poly"),r.RooFit.Components("poly*"),r.RooFit.LineColor(40))
        pFunc.plotOn(iFrame,r.RooFit.Name("Poly"),r.RooFit.Components("exp*"),r.RooFit.LineColor(6))

    print iFrame.Print("")

def draw(iVar,iData,iFuncs,iRooFitResult,iLabel="A"):
    lCan   = r.TCanvas(str(iLabel),str(iLabel),800,600)
    lFrame = iVar.frame()
    lLegend = getLegend()
    drawFrame(lFrame,iData,iFuncs,lLegend,50)
    lFrame.GetYaxis().SetRangeUser(0,lFrame.GetMaximum()*1.2)
    lFrame.GetYaxis().SetTitle(" Events / %.1f GeV"%fBinWidth);
    lFrame.GetXaxis().SetTitle("PUPPI Softdrop Jet Mass (GeV) ");
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
    addInfo = TPaveText(0.2,0.75,0.5,0.9,"NDC")
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
    pChi2     = iModel[0].createChi2(lDataHist,r.RooFit.Range(fXMin,fXMax),r.RooFit.Extended(r.kTRUE),r.RooFit.DataError(r.RooAbsData.Poisson))
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
    print '---- Making pdf for %s with model %s'%(iLabel,iModel) 
    lVar = iWorkspace.var(fVar);
    lModelPdf = None
    lTag = "%s_%s"%(iModel,iLabel)
    pVarHigh       = r.RooRealVar("%sHi_%s"%(lVar.GetName(),lTag),"%sHi_%s"%(lVar.GetName(),lTag),0.5,0.,1.);
    pVarHigh1      = r.RooRealVar("%sHi1_%s"%(lVar.GetName(),lTag),"%sHi1_%s"%(lVar.GetName(),lTag),0.5,0.,1.);
    lConstraints = []

    if iModel == "Exp_mc":
        lC_Exp_mc        = r.RooRealVar("c_mc"+lTag,"c_mc"+lTag,-0.01,-2.,0.05)
        lModelPdf = r.RooExponential("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_Exp_mc);
        
    if iModel == "Exp_data":
        lC_Exp_data      = r.RooRealVar("c_data"+lTag,"c_data"+lTag,-0.01,-2.,0.05)
        lModelPdf = r.RooExponential("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_Exp_data);
        
    if "ErfExpGaus_st" in iModel:
        lC_ErfExp        = r.RooRealVar("c"+lTag,"c"+lTag,-0.04,-0.2,0.0)
        lOffSet_ErfExp   = r.RooRealVar("offset"+lTag,"offset"+lTag,82.,75.,90.)
        lWidth_ErfExp    = r.RooRealVar("width"+lTag,"width"+lTag,30.,10.,300.)
        lMean_Gaus       = r.RooRealVar("mean"+lTag,"mean"+lTag,82.,75.,90.)
        lSigma_Gaus      = r.RooRealVar("sigma"+lTag,"sigma"+lTag,7.,0.,40.)
        pGaus            = r.RooGaussian("gaus"+lTag,"gaus"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pErfExp          = r.RooErfExpPdf("erfExp"+lTag,"erfExp"+lTag,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        lModelPdf = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pErfExp,pGaus),r.RooArgList(pVarHigh));

    if iModel == "ErfExp_st_mc":
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-5,1.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,82.,50.,1000.)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,50.,20.,1000.)
        lModelPdf = r.RooErfExpPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);

    if iModel == "ErfExp_st_data":
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-5,1.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,82.,50.,1000.)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,50.,20.,1000.)
        lModelPdf = r.RooErfExpPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);

    if "GausExp_tt" in iModel:
        if iLabel == "tt_mc" or iLabel == "tt_data":
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,80.,75.,90.)
        elif "fakeW" in iLabel:
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,80.,75.,90.)
        elif "realW" in iLabel:
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,80.,75.,90.)
        lSigma_Gaus        = r.RooRealVar("sigma_"+lTag,"sigma_"+lTag,7.,0.,40.)
        lC_Exp             = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-1.,1.)
        pGaus              = r.RooGaussian("gaus_"+lTag,"gaus_%s"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pExp               = r.RooExponential("exp_"+lTag,"exp_%s"+lTag,lVar,lC_Exp)
        lModelPdf          = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus,pPoly,pExp),r.RooArgList(pVarHigh,pVarHigh1),1)

    if  "GausErfExp_tt" in iModel:
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-10,10.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,20.,0.,300.)
        if "realW" in iLabel:
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,82.,75.,90.)
        elif "fakeW" in iLabel:
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,82.,75.,90.)
        else:
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,82.,75.,90.)
        a1             = r.RooRealVar("a1_"+lTag,"a1_"+lTag,0.0001,-0.01,0.01)
        a2             = r.RooRealVar("a2_"+lTag,"a2_"+lTag,0.00001,-0.0001,0.01)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,20.,0.,60.)
        lSigma_Gaus    = r.RooRealVar("sigma_"+lTag,"sigma_"+lTag,7.,0.,60.)
        pPoly          = r.RooPolynomial("poly_"+lTag,"poly_%s"+lTag,lVar,r.RooArgList(a1,a2),1)
        pErfExp        = r.RooErfExpPdf("erfExp_"+lTag,"erfExp_%s"+lTag,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        pGaus          = r.RooGaussian("gaus_"+lTag,"gaus_%s"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pPoly          = r.RooPolynomial("poly_"+lTag,"poly_%s"+lTag,lVar,r.RooArgList(a1),1)
        pGaus          = r.RooGaussian("gaus_"+lTag,"gaus_%s"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pErfExp        = r.RooErfExpPdf("erfExp_%s"%iLabel,"erfExp_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        lModelPdf = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus,pErfExp,pPoly),r.RooArgList(pVarHigh,pVarHigh),1);

    if iModel == "ErfExp_tt_mc" or iModel == "ErfExp_tt_data":
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-10,10.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,20.,0.,300.)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,20.,0.,60.)
        pErfExp        = r.RooErfExpPdf("erfExp_%s"%iLabel,"erfExp_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        lModelPdf = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pErfExp),r.RooArgList(pVarHigh),1);
        #addConstraint(iWorkspace,lMean_Gaus,175,3,lConstraints)
        #addConstraint(iWorkspace,lSigma_Gaus,7,2,lConstraints)
        #if 'realW' in iLabel:
        #    addConstraint(iWorkspace,lOffSet_ErfExp,175,10,lConstraints)
        #    addConstraint(iWorkspace,lWidth_ErfExp,30,5,lConstraints) 

    getattr(iWorkspace,"import")(lModelPdf,r.RooFit.RecycleConflictNodes())
    return iWorkspace.pdf("model_pdf_%s"%iLabel),lConstraints

# return RooExtendPdf and Constraints
def makeModel(iWorkspace,iLabel,iModel):
    print '---- Making model'
    lNumber = r.RooRealVar("number_%s"%iLabel,"number_%s"%iLabel,500,0.,1e7); 
    lModelPdf,lConstraints = makePdf(iWorkspace,iLabel,iModel)
    lModel = r.RooExtendPdf("model_%s"%iLabel,"model_%s"%iLabel,lModelPdf,lNumber)
    getattr(iWorkspace,"import")(lNumber,r.RooFit.RecycleConflictNodes())
    getattr(iWorkspace,"import")(lModel,r.RooFit.RecycleConflictNodes())
    iWorkspace.pdf("model_%s"%iLabel).Print()
    return iWorkspace.pdf("model_%s"%iLabel),lConstraints

# switch to MC
# make TT model with eff => is this only "sig" i.e. realWs?
# can have: tt_mc_pass,tt_mc_fail,tt_data_pass,tt_data_fail
# for pass = eff*number(real+fake)
# for fail = (1-eff)*number(real+fake)
def makeTTModel(iWorkspace,iLabel,iPtLabel,iCat,iModel):
    print '---- Making TT model for label %s, cat %s and model %s'%(iLabel,iPtLabel,iModel)
    nrealW_pass = iWorkspace.var('number_tt_realW_mc_'+iPtLabel+"_pass").getVal()
    nrealW_fail = iWorkspace.var('number_tt_realW_mc_'+iPtLabel+"_fail").getVal()
    nfakeW_pass = iWorkspace.var('number_tt_fakeW_mc_'+iPtLabel+"_pass").getVal()
    nfakeW_fail = iWorkspace.var('number_tt_fakeW_mc_'+iPtLabel+"_fail").getVal()
    effS      = nrealW_pass/(nrealW_pass+nfakeW_pass+nrealW_fail+nfakeW_fail)
    print ""; "Real MC pass / all = ",effS;
    effS      = nrealW_pass/(nrealW_pass+nrealW_fail)
    print "Real MC signal efficiency = " ,effS; print ""; 
    
    lLabel = iLabel+'_'+iPtLabel+'_'+iCat
    # tt_data
    if iLabel.find("tt_signal_data")!=-1:
        if iPtLabel.find("fail")==-1:
            lNumberTotal = r.RooRealVar("numbertotal_%s"%iLabel+'_'+iPtLabel,"numbertotal_%s"%iLabel+'_'+iPtLabel,500,0.,1e7);
            lEff = r.RooRealVar("eff_%s"%iLabel+'_'+iPtLabel,"eff_%s"%lLabel,effS,effS*0.8,effS*1.2);  
            lNumber = r.RooFormulaVar("number_%s_%s_pass"%(iLabel,iPtLabel), "@0*@1", r.RooArgList(lEff,lNumberTotal));
        else:
            lNumberTotal = iWorkspace.var("numbertotal_%s"%iLabel+'_'+iPtLabel)
            lEff = iWorkspace.var("eff_%s"%iLabel+'_'+iPtLabel)
            lNumber = r.RooFormulaVar("number_%s_%s_fail"%(iLabel,iPtLabel), "(1-@0)*@1", r.RooArgList(lEff,lNumberTotal));

    # tt_mc
    if iLabel.find("tt_signal_mc")!=-1:
        if iPtLabel.find("fail")==-1:
            lNumberTotal = r.RooRealVar("numbertotal_%s"%iLabel+'_'+iPtLabel,"numbertotal_%s"%iLabel+'_'+iPtLabel,500,0.,1e7);
            lEff = r.RooRealVar("eff_%s"%iLabel+'_'+iPtLabel,"eff_%s"%iLabel+'_'+iPtLabel,effS,effS*0.8,effS*1.2);
            lNumber = r.RooFormulaVar("number_%s_%s_pass"%(iLabel,iPtLabel), "@0*@1", r.RooArgList(lEff,lNumberTotal));
        else:
            lNumberTotal = iWorkspace.var("numbertotal_%s"%lLabel)
            lEff = iWorkspace.var("eff_%s"%iLabel+'_'+iPtLabel)
            lNumber = r.RooFormulaVar("number_%s_%s_fail"%(iLabel,iPtLabel), "(1-@0)*@1", r.RooArgList(lEff,lNumberTotal));

    lModelPdf,lConstraints = makePdf(iWorkspace,lLabel,iModel)
    lModel = r.RooExtendPdf("model_%s"%lLabel,"model_%s"%lLabel,lModelPdf,lNumber)  
    getattr(iWorkspace,"import")(lNumber,r.RooFit.RecycleConflictNodes())
    getattr(iWorkspace,"import")(lModel,r.RooFit.RecycleConflictNodes())
    iWorkspace.pdf("model_%s"%lLabel).Print()
    return iWorkspace.pdf("model_%s"%lLabel),lConstraints

# fit single rooDataset to model
# return list of constraints
def fitSingleMC(iWorkspace,iLabel,iModel):
    print '---- Fitting for %s'%iLabel
    lVar   = iWorkspace.var(fVar);
    lData  = iWorkspace.data(iLabel+"_D");
    lModel,lConstraints = makeModel(iWorkspace,iLabel,iModel)

    pRooFitResult = lModel.fitTo(lData,r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Extended(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.Verbose(r.kFALSE));
    getattr(iWorkspace,"import")(lModel,r.RooFit.RecycleConflictNodes())
    pRooFitResult.Print()
    getattr(iWorkspace,"import")(pRooFitResult)
    
    draw(lVar,lData,[lModel],pRooFitResult,iLabel+"_only")
    return lConstraints

class WPeak():
    def __init__(self,options):

        self._lMSD_lo = fXMin
        self._lMSD_hi = fXMax

        self._lW      = r.RooWorkspace("w","w")

        self._lMSD    = r.RooRealVar(fVar,fVar,self._lMSD_lo,self._lMSD_hi)
        self._lMSD.setBins(fNBins)
        getattr(self._lW,"import")(self._lMSD)

        self._lWeight = r.RooRealVar(fWeight,fWeight,-1e+25,1e+25)
        self._lPt     = r.RooRealVar(fPt,fPt,0,2000)

        self._lPDatas = {}
        self._lHPdfs  = {}
        self._lNPdfs  = {}
        self._lModels = {}
        self._lConstraints = {}
        self._lScaleNumber = {}

        self._TTsf = {};

        # get dataset for individual processes except for mc and tt_bkg(fakew), tt_signal(realW)
        for iLabel,iFile in fFiles.iteritems():
            if iLabel == "mc" or  "bkg" in iLabel or "signal" in iLabel: continue
            for Pt_it in range(len(fPt_bins)-1):
                pFile = fDir+fTag+"/"+iFile+"_"+fTag+".root"
                fCut = "(Puppijet0_pt>" + str(fPt_bins[Pt_it]) + " && Puppijet0_pt<" + str(fPt_bins[Pt_it+1]) + ")"
                Pt_label = "pT_" + str(fPt_bins[Pt_it]) + '_' + str(fPt_bins[Pt_it+1])
                for Cat_it in fCats:
                    self._TTsf[Pt_label+ '_' + Cat_it] = 1
                    print 'preparing sample %s in pT bin[%s,%s] and %s category'%(iLabel,fPt_bins[Pt_it],fPt_bins[Pt_it+1],Cat_it)
                    pName = iLabel + '_' + Pt_label + '_' + Cat_it
                    pCut = "(%s&&%s)"%(fCut,fWCut[Cat_it])
                    if 'wlnu' in iLabel: pCut = "(%s&&weight<0.2)"%pCut
                    if iLabel != "data":
                        self._lPDatas[pName] = self.getDataset(iLabel,pName,Pt_label+ '_' + Cat_it,pFile,pCut,True)
                    else:
                        self._lPDatas[pName] = self.getDataset(iLabel,pName,Pt_label+ '_' + Cat_it,pFile,pCut,False)
                    if iLabel != "data":
                        self._lConstraints[pName] = fitSingleMC(self._lW,pName,fPDFs[iLabel])
                    print self._lW.data(pName+"_D").sumEntries()

        # get normalization SF for tt for pass and fail
        for Pt_it in range(len(fPt_bins)-1):
            Pt_label = "pT_" + str(fPt_bins[Pt_it]) + '_' + str(fPt_bins[Pt_it+1])
            for Cat_it in fCats:
                self._TTsf[Pt_label+'_'+Cat_it] = self.getTTSF(Pt_label,Cat_it)
        
        # get dataset for mc pass and fail and scale it by tt_sf norm
        iLabel = "mc"
        print 'SF!! ',iLabel
        for Pt_it in range(len(fPt_bins)-1):
            pFile = fDir+fTag+"/"+fFiles[iLabel]+"_"+fTag+".root"
            Pt_label = "pT_" + str(fPt_bins[Pt_it]) + '_' + str(fPt_bins[Pt_it+1])
            for Cat_it in fCats:
                print 'preparing sample %s in pT bin[%s,%s] and %s category'%(iLabel,fPt_bins[Pt_it],fPt_bins[Pt_it+1],Cat_it)
                pName = iLabel + '_' + Pt_label + '_' + Cat_it
                pCut = "(%s&&%s)"%(fCut,fWCut[Cat_it])
                self._lPDatas[pName] = self.getDataset(iLabel,pName,Pt_label+ '_' + Cat_it,pFile,pCut,True)
                print self._lW.data(pName+"_D").sumEntries()

        print self._lPDatas
        self._lW.Print()
        
    # get dataset for each process from tree
    def getDataset(self,iLabel,iName,Pt_label,iFile,iCut="(1==1)",iWeight=False,iTree="otree2"):
        print 'preparing dataset for ',iFile,' with cut ',iCut, ' and TTsf ',self._TTsf[Pt_label]
        lFile   = r.TFile(iFile)
        lTree   = lFile.Get(iTree)
        lData = r.RooDataSet(iName+"_D",iName+"_D",r.RooArgSet(self._lMSD,self._lWeight),r.RooFit.WeightVar(self._lWeight))
        lCut = r.TTreeFormula("lCut",iCut,lTree)  
        lTree.SetNotify(lCut)       
        for i0 in range(lTree.GetEntriesFast()):
            lTree.LoadTree(i0) 
            lSel = False 
            for i1 in range(lCut.GetNdata()):  
                if (lCut.EvalInstance(i1)):  
                    lSel=True; break;
            if not lSel: continue     
            lTree.GetEntry(i0) 
            if iWeight:
                lWeight = getattr(lTree,fWeight)*fLumi*self._TTsf[Pt_label]
                #if 'wlnu' in iLabel: 
                #    print getattr(lTree,fWeight)
            else:
                lWeight = 1
            lMass = getattr(lTree,fVar)
            lMatched = 0
            jmatched = getattr(lTree,"Puppijet0_vMatching");
            jhadronic = getattr(lTree,"Puppijet0_isHadronicV");
            if jhadronic == 1.0 and jmatched < 0.8 and jmatched > 0.:
                lMatched = 1;
            if 'realW' in iLabel and lMatched == 0: continue
            if 'fakeW' in iLabel and lMatched == 1: continue
            if lMass < self._lMSD_hi and lMass > self._lMSD_lo:
                self._lMSD.setVal(lMass)
                lData.add(r.RooArgSet(self._lMSD), lWeight)
        getattr(self._lW,'import')(lData,r.RooFit.RecycleConflictNodes())
        lData.Print()
        lFile.Close()
        return self._lW.data(iName+"_D")
                     
    # getModel
    def getModel(self, iLabel):
        print 'getting model for ', iLabel
        lModel    = self._lW.pdf("model_"+iLabel)
        lModel.Print()
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

    # get TT normalization SF for pt label
    def getTTSF(self,iPtLabel,Cat_it):
        nData = nTT = nMinorBKG = 0;
        for Data_it in self._lPDatas:
            if Cat_it not in Data_it: continue
            if iPtLabel not in Data_it: continue
            nEntries = self._lW.data(Data_it+"_D").sumEntries()
            print Data_it,nEntries
            if Data_it == "data_"+iPtLabel+"_"+Cat_it:
                nData += nEntries
                print 'adding to data'
            if "_mc" in Data_it:
                if "realW" in Data_it or "fakeW" in Data_it:
                    nTT += nEntries
                    print 'adding to tt'
                else:
                    nMinorBKG += nEntries
                    print 'adding to bkg'
        print 'cat ', Cat_it, ' data ',nData,' tt ',nTT,' bkg ',nMinorBKG
        ttScalefactor = (nData-nMinorBKG)/nTT
        #ttScalefactor = 1
        scalefactor = r.RooRealVar("tt_scalefactor_"+Cat_it,"tt_scalefactor_"+Cat_it,ttScalefactor)
        getattr(self._lW,'import')(scalefactor,r.RooFit.RecycleConflictNodes())
        return scalefactor.getVal()
        #return 1

    # get Wtag SFs
    def getWtagSFs(self,iPtLabel,iModelLabel,iTTLabel):
        print 'ttlabel ',iTTLabel
        print 'pt ',iPtLabel
        # efficiency: realW tt evts:  pass / pass+fail
        # N2 sf = efficiency in data / efficiency in mc 
        # mass sf = mean in data / mean in mc
        # mass shift = mean in data - mean in mc
        # mass smear = sigma in data / sigma in mc
        # sigma sf = sigma in data / sigma in mc
        pEff = {}; pMean = {}; pSigma = {}
        pEffErr = {}; pMeanErr = {}; pSigmaErr = {}

        for i0 in ['mc','data']:
            lTTLabel = iModelLabel.replace('mc',i0)+"_"+iTTLabel.replace('mc',i0)
            if "signal" in lTTLabel:
                pEff[i0] = self._lW.var("eff_"+iTTLabel.replace('mc',i0)+"_"+iPtLabel).getVal(); 
                pEffErr[i0] = self._lW.var("eff_"+iTTLabel.replace('mc',i0)+"_"+iPtLabel).getError();
            if "realW" in lTTLabel:
                nrealW_pass = self._lW.var("number_"+iTTLabel.replace('mc',i0)+"_"+iPtLabel+"_pass").getVal(); 
                nrealW_fail = self._lW.var("number_"+iTTLabel.replace('mc',i0)+"_"+iPtLabel+"_fail").getVal()
                nrealW_pass_Err = self._lW.var("number_"+iTTLabel.replace('mc',i0)+"_"+iPtLabel+"_pass").getError();
                nrealW_fail_Err = self._lW.var("number_"+iTTLabel.replace('mc',i0)+"_"+iPtLabel+"_fail").getError();
                pEff[i0] = nrealW_pass/(nrealW_pass+nrealW_fail)
                pEffErr[i0] = pEff[i0] * math.sqrt( (nrealW_pass_Err/nrealW_pass)**2 + 
                                                    (math.sqrt( nrealW_pass_Err**2+nrealW_fail_Err**2 )/(nrealW_pass+nrealW_fail))**2 )

            print "mean_"+lTTLabel+"_"+iPtLabel+"_pass"
            pMean[i0] = self._lW.var("mean_"+lTTLabel+"_"+iPtLabel+"_pass").getVal();
            pMeanErr[i0] = self._lW.var("mean_"+lTTLabel+"_"+iPtLabel+"_pass").getError();
            
            pSigma[i0] = self._lW.var("sigma_"+lTTLabel+"_"+iPtLabel+"_pass").getVal(); 
            pSigmaErr[i0] = self._lW.var("sigma_"+lTTLabel+"_"+iPtLabel+"_pass").getError(); 

        pEffSf = pEff['data']/pEff['mc'];
        pEffSfErr = pEffSf * math.sqrt( (pEffErr['data']/pEff['data'])**2 + (pEffErr['mc']/pEff['mc'])**2 )

        pMeanShift = pMean['data'] - pMean['mc']
        pMeanShiftErr = math.sqrt(pMeanErr['data']**2 +pMeanErr['mc']**2)

        pSigmaSmear = pSigma['data']/pSigma['mc']
        pSigmaSmearErr = pSigmaSmear * math.sqrt( (pSigmaErr['data']/pSigma['data'])**2 + (pSigmaErr['mc']/pSigma['mc'])**2 )

        pMeanSf = pMean['data']/pMean['mc']
        pMeanSfErr = pMeanSf * math.sqrt( (pMeanErr['data']/pMean['data'])**2 + (pMeanErr['mc']/pMean['mc'])**2 )

        pSigmaSf = pSigmaSmear
        pSigmaSfErr =  pSigmaSmearErr

        print 'W-tagging SF: %.3f +/- %.3f '%(pEffSf,pEffSfErr)
        print 'Mass shift [GeV]: %.3f +/- %.3f '%(pMeanShift,pMeanShiftErr)
        print 'Mass SF: %.3f +/- %.3f '%(pMeanSf,pMeanSfErr)
        print 'Resolution SF: %.3f +/- %.3f '%(pSigmaSf,pSigmaSfErr)
        print '\n'
        print 'eff data: %.3f +/- %.3f, mc: %.3f +/- %.3f'%(pEff['data'],pEffErr['data'],pEff['mc'],pEffErr['mc'])
        print '<m> data: %.3f +/- %.3f, mc: %.3f +/- %.3f'%(pMean['data'],pMeanErr['data'],pMean['mc'],pMeanErr['mc'])
        print 'res data: %.3f +/- %.3f, mc: %.3f +/- %.3f'%(pSigma['data'],pSigmaErr['data'],pSigma['mc'],pSigmaErr['mc'])

    # fit mc and data
    def fit(self):

        for Pt_it in range(len(fPt_bins)-1):
            
            # define categories
            pCats  = r.RooCategory("pCats"  ,"pCats")
            for Cat_it in fCats:
                pCats.defineType(Cat_it)

            pPtLabel = "pT_" + str(fPt_bins[Pt_it]) + "_" + str(fPt_bins[Pt_it+1])

            # nData
            print self._lPDatas
            nData =self._lPDatas["data_"+pPtLabel+"_pass"].sumEntries()+self._lPDatas["data_"+pPtLabel+"_fail"].sumEntries()
            print 'nData ',nData

            # get models
            for Cat_it in fCats:
                print 'processing combined fit for pT bin[%s,%s] and %s category'%(fPt_bins[Pt_it],fPt_bins[Pt_it+1],Cat_it)
                print 'and with pass sf: %.3f and fail sf: %.3f'%(self._lW.var("tt_scalefactor_pass").getVal(),self._lW.var("tt_scalefactor_fail").getVal())

                pPtLabelCat = pPtLabel + "_"+Cat_it

                lVar   = self._lW.var(fVar);
                lSData = self._lPDatas["data_"+pPtLabelCat]
                lSMc   = self._lPDatas["mc_"+pPtLabelCat]
                
                # constraints
                pPdfConstraints_Mc   = r.RooArgSet("pPdfConstraints_Mc")
                pPdfConstraints_Data = r.RooArgSet("pPdfConstraints_Data")

                nBkg = 0
                pSamples = ['wlnu_mc','wlnu_data','st_mc','st_data','tt_fakeW_mc','tt_fakeW_data','tt_realW_mc','tt_realW_data']
                for iLabel in pSamples:
                    lLabel = iLabel + '_' + pPtLabelCat

                    if "_data" not in iLabel:
                        if 'tt' in iLabel:
                            nBkg += self._lW.data(lLabel+"_D").sumEntries()*self._TTsf[pPtLabelCat]
                        else:
                            nBkg += self._lW.data(lLabel+"_D").sumEntries()

                    self._lHPdfs[lLabel] = self.histPdf(lLabel)
                    self._lModels[lLabel] = self.getModel(lLabel)

                    for i0 in range(len(self._lConstraints[lLabel])):
                        if 'mc' in lLabel:
                            pPdfConstraints_Mc.add(self._lW.pdf(self._lConstraints[lLabel][i0]))
                        if 'data' in lLabel:
                            pPdfConstraints_Data.add(self._lW.pdf(self._lConstraints[lLabel][i0]))
                            
                # build tt model for realw(signal) so that eff is included in parameter of fit
                xConstraints = {}
                for iLabel in ['tt_signal_mc','tt_signal_data']:
                    lLabel = iLabel + '_' + pPtLabelCat
                    self._lModels[lLabel],xConstraints[lLabel]  =  makeTTModel(self._lW,iLabel,pPtLabel,Cat_it,fPDFs[iLabel])

                # build model for tt bkg (fakeW) for mc and data                                                                                                                 
                for iLabel in ["tt_bkg_mc","tt_bkg_data"]:
                    lLabel = iLabel + '_' + pPtLabelCat
                    lModel,lConstraints = makeModel(self._lW,lLabel,fPDFs[iLabel])
                    self._lModels[lLabel] = self.getModel(lLabel)


                # constraints
                pPdfConstraints_Mc.Print()
                pPdfConstraints_Data.Print()

                # fix variables from MC and only leave floating normalization, mean and sigma from TT
                # variables floating: eff(realW), numbertotal(realW)
                args = self._lW.allVars()
                args.Print()
                iter = args.createIterator()
                var = iter.Next()
                lFloating = ["mean_GausErfExp_tt_",
                             "c_GausErfExp_tt_",
                             "c_ErfExp_tt_",
                             "offset_GausErfExp_tt_",
                             "offset_ErfExp_tt_",
                             "sigma_GausErfExp_tt_",
                             "width_ErfExp_tt_",
                             "width_GausErfExp_tt_",
                             #"a1_GausErfExp_tt_",
                             "number_tt",
                             "eff_tt_signal_",
                             "numbertotal_",
                             ]
                while var:
                    #print var.GetName()
                    if any(f in var.GetName() for f in lFloating):
                        print 'float: ',var.GetName()
                        pass
                    else:
                        var.setConstant(r.kTRUE)
                        print 'set constant'
                    var = iter.Next()

                # build model for tt bkg (fakeW) for mc and data
                for iLabel in ["tt_bkg_mc","tt_bkg_data"]:
                    lLabel = iLabel + '_' + pPtLabelCat
                    lModel,lConstraints = makeModel(self._lW,lLabel,fPDFs[iLabel])
                    self._lModels[lLabel] = self.getModel(lLabel)

                self._lModels['TotalMc_'+pPtLabelCat] = r.RooAddPdf(("model_total_mc_"+pPtLabelCat),("model_totalmc_"+pPtLabelCat),
                                                                    r.RooArgList(self._lModels["tt_signal_mc_"+pPtLabelCat],
                                                                                 self._lModels["tt_bkg_mc_"+pPtLabelCat],
                                                                                 self._lModels['wlnu_mc_'+pPtLabelCat],
                                                                                 self._lModels['st_mc_'+pPtLabelCat]))
                self._lModels['TotalData_'+pPtLabelCat] = r.RooAddPdf("model_total_data_"+pPtLabelCat,"model_totaldata_"+pPtLabelCat,
                                                                      r.RooArgList(self._lModels["tt_signal_data_"+pPtLabelCat],
                                                                                   self._lModels["tt_bkg_data_"+pPtLabelCat],
                                                                                   self._lModels['wlnu_data_'+pPtLabelCat],
                                                                                   self._lModels['st_data_'+pPtLabelCat]))

                getattr(self._lW,"import")(self._lModels["TotalMc_"+pPtLabelCat],r.RooFit.RecycleConflictNodes())
                getattr(self._lW,"import")(self._lModels["TotalData_"+pPtLabelCat],r.RooFit.RecycleConflictNodes())

                # fit to data (cat)
                pRooFitResult_Data = self._lModels['TotalData_'+pPtLabelCat].fitTo(lSData,r.RooFit.Strategy(1),r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Minimizer("Minuit2"),
                                                                                   r.RooFit.ExternalConstraints(pPdfConstraints_Data))
                pRooFitResult_Data.Print()
                lTagData = 'data_only_'+pPtLabelCat
                lTagData+= '_ttmatched'
                draw(lVar,lSData,[self._lModels['TotalData_'+pPtLabelCat]],pRooFitResult_Data,lTagData)
                x2 = self._lModels['TotalData_'+pPtLabelCat].getVariables()
                print x2.Print("v")
                
                # fit to mc (cat)
                pRooFitResult_Mc   = self._lModels['TotalMc_'+pPtLabelCat].fitTo(lSMc,r.RooFit.Strategy(1),r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Minimizer("Minuit2"),
                                                                                 r.RooFit.ExternalConstraints(pPdfConstraints_Mc))
                lTagMc = 'mc_only_'+pPtLabelCat
                lTagMc += '_ttmatched'
                draw(lVar,lSMc,[self._lModels['TotalMc_'+pPtLabelCat]],pRooFitResult_Mc,lTagMc)
                x1 = self._lModels['TotalMc_'+pPtLabelCat].getVariables()
                print x1.Print("v")

                # draw data and mc (for that cat)
                lTagDataMc = "data_mc_"+pPtLabelCat
                lTagDataMc += "_ttmatched"
                params = {}
                drawDataMc(lVar,lSData,[self._lModels['TotalData_'+pPtLabelCat]],lSMc,[self._lModels['TotalMc_'+pPtLabelCat]],pRooFitResult_Mc,pRooFitResult_Data,params,lTagDataMc)


            # combined data (pass and fail)
            combData_data = r.RooDataSet("combData_data","combData_data",r.RooArgSet(self._lMSD,self._lWeight),r.RooFit.WeightVar(self._lWeight),
                                         RooFit.Index(pCats),
                                         RooFit.Import("data_"+pPtLabel+"_pass",self._lPDatas["data_"+pPtLabel+"_pass"]),
                                         RooFit.Import("data_"+pPtLabel+"_fail",self._lPDatas["data_"+pPtLabel+"_fail"]))
            # combined mc (pass and fail)
            combData_mc = r.RooDataSet("combData_mc","combData_mc",r.RooArgSet(self._lMSD,self._lWeight),r.RooFit.WeightVar(self._lWeight),
                                       RooFit.Index(pCats),
                                       RooFit.Import("mc_"+pPtLabel+"_pass",self._lPDatas["mc_"+pPtLabel+"_pass"]),
                                       RooFit.Import("mc_"+pPtLabel+"_fail",self._lPDatas["mc_"+pPtLabel+"_fail"]))

            # simultaneous fit with tt model
            simPdf_total_data = r.RooSimultaneous("simPdf_total_data","simPdf_total_data",pCats)
            simPdf_total_mc   = r.RooSimultaneous("simPdf_total_mc"  ,"simPdf_total_mc"  ,pCats)

            for Cat_it in fCats:
                simPdf_total_data.addPdf(self._lW.pdf("model_total_data_"+pPtLabel+"_"+Cat_it),"data_"+pPtLabel+"_"+Cat_it)
                simPdf_total_mc  .addPdf(self._lW.pdf("model_total_mc_"  +pPtLabel+"_"+Cat_it),"mc_"  +pPtLabel+"_"+Cat_it)

            # do simult fit with tt model
            print "simultaneous pass and fail fit with tt model"
            simFit_total_mc   = simPdf_total_mc  .fitTo(combData_mc,r.RooFit.Save(r.kTRUE),r.RooFit.Verbose(r.kFALSE),r.RooFit.Minimizer("Minuit2"),r.RooFit.SumW2Error(r.kTRUE))
            simFit_total_mc.Print()
            simFit_total_data = simPdf_total_data.fitTo(combData_data,r.RooFit.Save(r.kTRUE),r.RooFit.Verbose(r.kFALSE),r.RooFit.Minimizer("Minuit2"),r.RooFit.SumW2Error(r.kTRUE))
            simFit_total_data.Print()

            # get Wtag with tt model
            self.getWtagSFs(pPtLabel,fPDFs["tt_signal_mc"],"tt_signal_mc");

            # draw simult fit with tt model
            drawDataMc(lVar,lSData,
                       [self._lW.pdf("model_total_data_"+pPtLabel+"_pass")],lSMc,[self._lW.pdf("model_total_mc_"+pPtLabel+"_pass")],
                       simFit_total_mc,simFit_total_data,params,"data_pass_simult_tt")
            drawDataMc(lVar,lSData,
                       [self._lW.pdf("model_total_data_"+pPtLabel+"_fail")],lSMc,[self._lW.pdf("model_total_mc_"+pPtLabel+"_fail")],
                       simFit_total_mc,simFit_total_data,params,"data_fail_simult_tt")

            # write workspace
            self._lW.writeToFile(fOutput)

if __name__ == "__main__":
    options = parser()
    print options
    global fTag,fXMin,fXMax,fNBins
    fTag   = options.tag
    fXMin  = options.xmin
    fXMax  = options.xmax
    fNBins = int( (fXMax - fXMin) / fBinWidth)
    # get roodataset and make individual fits
    lW = WPeak(options);
    # make combined fit
    lW.fit();
