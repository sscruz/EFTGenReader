import ROOT as r 
from math import sqrt
import numpy as np 
r.gROOT.SetBatch(True)

from inputs import inputs


def getIndexCoef( c1, c2, fit):
    pairs=fit.getPairs()
    names=fit.getNames()

    thepair=fit.getIndexPair(c1,c2)
    for i,pair in enumerate(pairs):
        if pair.first == thepair.first and pair.second==thepair.second:
            idx=i
            break
    return idx

def getRawCoeff(c1,c2, fit):
    idx=getIndexCoef(c1,c2, fit)
    return fit.getCoefficient(idx)

def getRawCovTerm( term1, term2, fit):
    idx1=getIndexCoef(term1[0],term1[1], fit)
    idx2=getIndexCoef(term2[0],term2[1], fit)
    errpairs=fit.getErrorPairs()
    errcoefs=fit.getErrorCoefficients()

    for i,errcoef in enumerate(errcoefs):
        errpair=errpairs.at(i)
        if errpair.first==idx1 and errpair.second==idx2:
            # idx_pair=pairs.at(errpair.first)
            # print names.at(idx_pair.first)
            # print names.at(idx_pair.second)
            # idx_pair=pairs.at(errpair.second)
            # print names.at(idx_pair.first)
            # print names.at(idx_pair.second)
            # print errcoef
            return errcoef
            break

def normTermAndError(c1, c2, fit):
    rawterm=getRawCoeff(c1,c2,fit)
    rawsm=getRawCoeff('sm','sm', fit)
    norm_coef=rawterm/rawsm
    norm_err=getRawCovTerm((c1,c2), (c1,c2), fit)/rawsm**2+(rawterm**2/rawsm**4)*getRawCovTerm(('sm','sm'),('sm','sm'), fit)-2*(rawterm/rawsm**3)*getRawCovTerm((c1,c2),('sm','sm'),fit)
    if norm_err<0:
        norm_err=0

    return norm_coef, sqrt(norm_err)

def printOperatorList(oplist, fit):
    for coef in oplist:
        print "{0: <15}    {1: <15}    {2: <15}".format('sm*sm', '%s*sm'%coef, '%s*%s'%(coef,coef))
        sm_sm=normTermAndError('sm','sm',fit)
        eft_sm=normTermAndError(coef,'sm',fit)
        eft_eft=normTermAndError(coef,coef,fit)
        print "{0:.4f}+/-{1:.4f}    {2:.4f}+/-{3:.4f}    {4:.4f}+/-{5:.4f}".format(sm_sm[0],sm_sm[1], eft_sm[0],eft_sm[1], eft_eft[0],eft_eft[1])
        print ''

def plotOperatorDependence(operator, fit, process):
    rang=np.linspace(-10,10,20)
    gr=r.TGraphErrors(len(rang)-1)
    sm=r.WCPoint()
    c=r.TCanvas()
    for i, xval in enumerate(rang):
        gr.SetPoint(i,xval,fit.evalPoint(operator, xval)/fit.evalPoint(sm))
        gr.SetPointError(i,0,fit.evalPointError(operator,xval)/fit.evalPoint(sm))
    gr.SetFillColor(6);
    gr.SetFillStyle(3005);
    gr.Draw("a3")
    c.SaveAs("%s_%s.png"%(process, operator))

# for op in ["cpQM","ctp","ctW","ctZ","ctG","cbW","cpQ3","cptb","cpt","cQl3i","cQlMi","cQei","ctli","ctei","ctlSi","ctlTi"]+["cQQ1","cQQ8","cQt1","cQt8","ctt1","ctb1","cQtQb1","cQtQb8"]+["cQq13","cQq83","cQq11","cQq81","cQu1","cQu8","cQd1","cQd8","ctq1","ctq8","ctu1","ctu8","ctd1","ctd8"]:
#     plotOperatorDependence(op)

def summaryPlot(operators, fit, process):
    h=r.TH1F("summary_%s"%process,"", 2*len(operators), -0.5, 2*len(operators)-0.5)

    for i,coef in enumerate(operators):
        h.GetXaxis().SetBinLabel(2*i+1, coef + ', lin')
        eft_sm=normTermAndError(coef,'sm',fit)
        h.SetBinContent(2*i+1,eft_sm[0])
        h.SetBinError  (2*i+1,eft_sm[1])

        h.GetXaxis().SetBinLabel(2*i+2, coef + ', quad')
        eft_eft=normTermAndError(coef,coef,fit)
        h.SetBinContent(2*i+2,eft_eft[0])
        h.SetBinError  (2*i+2,eft_eft[1])
    h.SetDirectory(0)
    return h
#    tf=r.TFile.Open("summary_%s.root"%process,'recreate')
#    tf.WriteTObject(h,'summary')
#    tf.Close()


r.gROOT.SetBatch(True)
r.gROOT.ProcessLine(".x tdrstyle.cc")
r.gStyle.SetOptStat(0)
r.gStyle.SetOptTitle(0)

hists=[]
legend=r.TLegend(0.2,0.5,0.5,0.8)
for process in inputs:

    tf=r.TFile.Open(process['file'])
    hist=tf.Get('EFTGenHistsWithCuts/h_eventsumEFT')
    fit=hist.GetBinFit(1)



    idx=-1

    

#    hist = summaryPlot( ["ctp","cpQM","ctW","ctZ","ctG","cbW","cpQ3","cptb","cpt","cQl3i","cQlMi","cQei","ctli","ctei","ctlSi","ctlTi"] + 
#                        ["cQQ1","cQQ8","cQt1","cQt8","ctt1","ctb1","cQtQb1","cQtQb8"]+ 
#                        ["cQq13","cQq83","cQq11","cQq81","cQu1","cQu8","cQd1","cQd8","ctq1","ctq8","ctu1","ctu8","ctd1","ctd8"], fit, process['name'])

    hist = summaryPlot( ["cQQ1","cQQ8","cQt1","cQt8","cQb1","cQb8","ctt1","ctb1","cQtQb1","cQtQb8","cQq13","cQq83","cQq11","ctq1","cQq81","ctq8"] , fit, process['name'])
    hist.SetLineColor( process['color'])
    hist.SetLineWidth(2)
    legend.AddEntry( hist, process['name'], 'l' )
    hists.append(hist)

hists.sort( key = lambda x : x.GetMaximum(), reverse=True)

c = r.TCanvas('canvas', '', 2560, 1440)
for i, h in enumerate(hists): 
    if i: 
        h.Draw('hist,same')
    else: 
        h.Draw('hist')
legend.Draw('same')
c.SaveAs("summary.pdf")    
c.SaveAs("summary.C")    

    # print "Two quarks and bosons, two quarks and two leptons"     
    # printOperatorList(["ctp","cpQM","ctW","ctZ","ctG","cbW","cpQ3","cptb","cpt","cQl3i","cQlMi","cQei","ctli","ctei","ctlSi","ctlTi"],fit)
    
    # print "Four heavy quarks" 
    # printOperatorList(["cQQ1","cQQ8","cQt1","cQt8","ctt1","ctb1","cQtQb1","cQtQb8"]) # "cQb1","cQb8"

    # print "Two heavy-two light quarks" 
    # printOperatorList(["cQq13","cQq83","cQq11","cQq81","cQu1","cQu8","cQd1","cQd8","ctq1","ctq8","ctu1","ctu8","ctd1","ctd8"])

    


#    getCoeffAndError(coef,coef), getCoeffAndError(coef,'sm')
#print getRawCovTerm(("cpQM","cpQM"),("cpQM","cpQM"))
#print getRawCoeff("cpQM","cpQM")
#print normTermAndError('cpQM','cpQM')



