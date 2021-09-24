#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Tue Mar 21 17:27:35 2017
# @author: stefano.magni@uni.lu

print("Starting...")
from Script_data import *  # STE
from Script_information import *  # STE
from genes_1 import *  # STE
from scipy.stats.stats import pearsonr
import Script_data as sd  # STE
import numpy # STE
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import matplotlib.ticker as plticker
import statsmodels.stats.weightstats as wstats
import seaborn as sns
import copy
import matplotlib as mpl
import math
import os


# Tags to decide what to run
GetDataFromFiles = True # True     # NECESSARY = True FOR JOCHEN ANALYSYS!!!
ComputeKneePlots = True # True
PlotKneePlots = True # True
CutCellNumber = True # True         # NECESSARY = True FOR JOCHEN ANALYSYS!!!

SaveTransformedData = False # False
GetFromFileTransformedData = False # False

DimensionalityReductionColorStem = False # False

DRwithList_stemGenesOnly = False # False
DRwithList_CellCycleGenesOnly = False # False
DRwithList_MitoGenesOnly = False # False
DRwithList_ROSGenesOnly = False # False

DoComparisonPreliminarySteps = True # True   # NECESSARY = 1 FOR JOCHEN ANALYSYS!!!
DoComparison = False # False 
ProduceListOfDifferentiallyExpressedGenes = True # True

PlotHistogramsIndividualGenes = False # False
PlotHistogramsLRRK2 = False # False
PlotHistogramsPAX6 = False # False

PlotHistogramsCumulative = True # True
PlotHistogramsCumulativeCellDeath = True # True

PlotHistogramsCumulativeLinnarsonDAneurons = False # False False
LoadLuisData = False # False False
LuisDataPreliminarySteps = False # False False

ReDoCumulativeWithLinnarsonLists = False # False False

PlotHistogramsIndividualGenesCoreRegulatCircuitLasse = False # False
PlotHistogramsIndividualGenesCoreRegulatCircuitLasseExtendedGenes = False # False
PlotHistogramsIndividualGenesCoreRegulatCircuitLasse_5_TOP_GENES = True # True
PlotHistogramsIndividualGenes_SRR = False # False
DoComparison_Lasse = False # False
DoCorrelationCRC_Genes_Lasse_ALL_cells = True # True
DoCorrelationCRC_Genes_Lasse_LowHigh_NR2F1_cells = False # False

Plots_ZIC_12345_Lasse_Jocken = False # False False
DoCorrelation_ZIC_12345_Lasse_Jocken =  False # False False

#RemoveDoubleNamesFromGeneNamesLists = 1
ComputeIntersectionGenesLists = False # False
HeatMapOfRowExpressionData =  False # False False

JochenAnalysis = False # False False               # NECESSARY =1 FOR JOCHEN ANALYSYS!!!
CumulativeGeneExpressionCRC = False  # False False # NECESSARY =1 FOR JOCHEN ANALYSYS!!! But do it alternatively to subsequents!!!
UsePCApreviousTotSNE = False  # False False # True = use PCA before tSNE, False = do not use it.
NumberOfPCAbeforetSNE = 6
FortSNEUseAllGenesInsteadOfListDAneurons = False  # False False
DRwithList_CRCgenesOnly_Mutant = False  # False False # NECESSARY =1 FOR JOCHEN ANALYSYS!!! But do it alternatively to subsequents!!!
DRwithList_CRCgenesOnly_Control = False  # False False # NECESSARY =1 FOR JOCHEN ANALYSYS!!! But do it alternatively to subsequents!!!
#GaussianMixtureModelOnALLtSNE = 1  # False False # NECESSARY =1 FOR JOCHEN ANALYSYS!!! But do it alternatively to subsequents!!!
GaussianMixtureModelOnIndividualtSNEwithCumGenExpr = False # False False # NECESSARY =1 FOR JOCHEN ANALYSYS!!! But do it alternatively to subsequents!!!
NumberOfClustersChosen = 8 

ComputeIntersectionGenesListsDEG_MitoOnly_Days = True # True # THIS FOR MITO PAPER: Walter et al., 2019, Stem Cell Reports. 2019 May 14;12(5):878-889. doi: 10.1016/j.stemcr.2019.03.004. Epub 2019 Apr 11. 

PlotHistogramsERN1 = False # False False


MYDIR = ("Results")
CHECK_FOLDER = os.path.isdir(MYDIR)
# If folder doesn't exist, then create it.
if not CHECK_FOLDER:
    os.makedirs(MYDIR)
    print("created folder : ", MYDIR)
else:
    print(MYDIR, "folder already exists.")


if GetDataFromFiles == 1:
    print("Getting data from file, Day 0...")
    A_0_genes, A_0_expr, A_0_expr_np = sd.get_from_txt('./Data_SingleCellRNASeq_DigitalExpressionMatrices/A_0.txt')
    B_0_genes, B_0_expr, B_0_expr_np = sd.get_from_txt('./Data_SingleCellRNASeq_DigitalExpressionMatrices/B_0.txt')
    print("Getting data from file, Day 10...")
  
    A_10_genes, A_10_expr, A_10_expr_np = sd.get_from_txt('./Data_SingleCellRNASeq_DigitalExpressionMatrices/A_10.txt')
    B_10_genes, B_10_expr, B_10_expr_np = sd.get_from_txt('./Data_SingleCellRNASeq_DigitalExpressionMatrices/B_10.txt')
    print("Getting data from file, Day 14...")
    
    A_14_genes, A_14_expr, A_14_expr_np = sd.get_from_txt('./Data_SingleCellRNASeq_DigitalExpressionMatrices/A_14.txt')
    B_14_genes, B_14_expr, B_14_expr_np = sd.get_from_txt('./Data_SingleCellRNASeq_DigitalExpressionMatrices/B_14.txt')
    print("Getting data from file, Day 42...")
    
    A_42_genes, A_42_expr, A_42_expr_np = sd.get_from_txt('./Data_SingleCellRNASeq_DigitalExpressionMatrices/A_42.txt')
    B_42_genes, B_42_expr, B_42_expr_np = sd.get_from_txt('./Data_SingleCellRNASeq_DigitalExpressionMatrices/B_42.txt')
    print("Getting data from file, all days done.")

if ComputeKneePlots == 1:
    print("Compuing knee plots...")
    a0,b0 = sd.get_knee_plot(A_0_expr_np)  # STE
    a10,b10 = sd.get_knee_plot(A_10_expr_np)  # STE
    a14,b14 = sd.get_knee_plot(A_14_expr_np)  # STE
    a42,b42 = sd.get_knee_plot(A_42_expr_np)  # STE
    
    a0B,b0B = sd.get_knee_plot(B_0_expr_np)  # STE
    a10B,b10B = sd.get_knee_plot(B_10_expr_np)  # STE
    a14B,b14B = sd.get_knee_plot(B_14_expr_np)  # STE
    a42B,b42B = sd.get_knee_plot(B_42_expr_np)  # STE
    
    print("Computing knee plots done.")

if PlotKneePlots == 1:
    plt.clf()    
    plt.xlim([0,2000])  # STE
    plt.xlabel('Cell Number')  # STE
    plt.ylabel('Cumulative fraction of transcripts')  # STE
    plt.plot(a0, label='Day 0')  # STE
    plt.plot(a10, label='Day 10')  # STE
    plt.plot(a14, label='Day 14')  # STE
    plt.plot(a42, label='Day 42')  # STE
    plt.legend(loc='center right')
    plt.title('G2019S')
    plt.savefig('./Results/KneePlotCumulative_G2019S.pdf')
    
    plt.clf() 
    plt.xlim([-10,1000])  # STE
    plt.xlabel('Cell Number')  # STE
    plt.ylabel('Number of transcripts')  # STE   
    plt.plot(b0, label='Day 0')  # STE
    plt.plot(b10, label='Day 10')  # STE
    plt.plot(b14, label='Day 14')  # STE
    plt.plot(b42, label='Day 42')  # STE
    plt.legend(loc='center right')
    plt.title('G2019S')
    plt.savefig('./Results/KneePlotAbsolute_G2019S.pdf')
    
    plt.clf()    
    plt.xlim([0,2000])  # STE
    plt.xlabel('Cell Number')  # STE
    plt.ylabel('Cumulative fraction of transcripts')  # STE
    plt.plot(a0B, label='Day 0')  # STE
    plt.plot(a10B, label='Day 10')  # STE
    plt.plot(a14B, label='Day 14')  # STE
    plt.plot(a42B, label='Day 42')  # STE
    plt.legend(loc='center right')
    plt.title('GC')
    plt.savefig('./Results/KneePlotCumulative_GC.pdf')
    
    plt.clf()  
    plt.xlim([-10,1000])  # STE
    plt.xlabel('Cell Number')  # STE
    plt.ylabel('Number of transcripts')  # STE  
    plt.plot(b0B, label='Day 0')  # STE
    plt.plot(b10B, label='Day 10')  # STE
    plt.plot(b14B, label='Day 14')  # STE
    plt.plot(b42B, label='Day 42')  # STE
    plt.legend(loc='center right')
    plt.title('GC')
    plt.savefig('./Results/KneePlotAbsolute_GC.pdf')
    
    print("Knee plots plotted.")

if CutCellNumber == 1:
    print("Cutting cell numbers, Day 0...")
    A_0 = sd.get_top_cells_matrix(A_0_expr_np,500)[0] # Tomasz = 500
    B_0 = sd.get_top_cells_matrix(B_0_expr_np,500)[0] # Tomasz = 500
    A_0_norm, A_0_log, A_0_sca = sd.get_normalization(A_0)
    B_0_norm, B_0_log, B_0_sca = sd.get_normalization(B_0)
    
    print("Cutting cell numbers, Day 10...")
    A_10 = sd.get_top_cells_matrix(A_10_expr_np,250)[0] # Tomasz = 250
    B_10 = sd.get_top_cells_matrix(B_10_expr_np,250)[0] # Tomasz = 250
    A_10_norm, A_10_log, A_10_sca = sd.get_normalization(A_10)
    B_10_norm, B_10_log, B_10_sca = sd.get_normalization(B_10)


    print("Cutting cell numbers, Day 14...")    
    A_14 = sd.get_top_cells_matrix(A_14_expr_np,505)[0]#250 # Tomasz = 505
    B_14 = sd.get_top_cells_matrix(B_14_expr_np,350)[0]#250 # Tomasz = 350
    #Here, I would take 250 cells from each subpopulation. The numbers 505 and 350 are the
    #ones used by Suaresh. I may be too conservative, I'm not sure. 
    #There's no big difference in reults after all.
    A_14_norm, A_14_log, A_14_sca = sd.get_normalization(A_14)
    B_14_norm, B_14_log, B_14_sca = sd.get_normalization(B_14)
    
    print("Cutting cell numbers, Day 42...")
    A_42 = sd.get_top_cells_matrix(A_42_expr_np,400)[0] # Tomasz = 400
    B_42 = sd.get_top_cells_matrix(B_42_expr_np,300)[0] # Tomasz = 300
    A_42_norm, A_42_log, A_42_sca = sd.get_normalization(A_42)
    B_42_norm, B_42_log, B_42_sca = sd.get_normalization(B_42)
    
    print("Cutting cell numbers, done.")
    
    #C_0 = Compare((A_0_norm,A_0_genes),(B_0_norm,B_0_genes)) 
    #C_10 = Compare((A_10_norm,A_10_genes),(B_10_norm,B_10_genes)) 
    #C_14 = Compare((A_14_norm,A_14_genes),(B_14_norm,B_14_genes)) 
    #C_42 = Compare((A_42_norm,A_42_genes),(B_42_norm,B_42_genes)) 
    

if SaveTransformedData == 1:
    
    print("Here I am 18")
    scaAlist = [A_0_sca,A_10_sca,A_14_sca,A_42_sca]
    genesAlist = [A_0_genes,A_10_genes,A_14_genes,A_42_genes]
    normAlist = [A_0_norm,A_10_norm,A_14_norm,A_42_norm]
    logAlist = [A_0_log,A_10_log,A_14_log,A_42_log]
    Alist = [A_0,A_10,A_14,A_42]
    
    scaBlist = [B_0_sca,B_10_sca,B_14_sca,B_42_sca]
    genesBlist = [B_0_genes,B_10_genes,B_14_genes,B_42_genes]
    normBlist = [B_0_norm,B_10_norm,B_14_norm,B_42_norm]
    logBlist = [B_0_log,B_10_log,B_14_log,B_42_log]
    Blist = [B_0,B_10,B_14,B_42]
    
    i = 0
    for DAY in [0,10,14,42]: 
        
        sca = scaAlist[i]
        genes = genesAlist[i]  
        norm = normAlist[i]
        log = logAlist[i]
        A = Alist[i]
        #You can use it to save transformed data
        sd.save_gene_list(genes,'./Results/IntermediateFiles/G2019S_day' + str(DAY) + '_genes.txt')
        sd.save_array(norm,'./Results/IntermediateFiles/G2019S_day' + str(DAY) + '_norm.csv')
        sd.save_array(log,'./Results/IntermediateFiles/G2019S_day' + str(DAY) + '_log.csv')
        sd.save_array(sca,'./Results/IntermediateFiles/G2019S_day' + str(DAY) + '_sca.csv')
        sd.save_array(A,'./Results/IntermediateFiles/G2019S_day' + str(DAY) + '.csv')
        print("Here I am 19a")
        
        sca = scaBlist[i]
        genes = genesBlist[i]  
        norm = normBlist[i]
        log = logBlist[i]
        b = Blist[i]
        #You can use it to save transformed data
        sd.save_gene_list(genes,'./Results/IntermediateFiles/GC_day' + str(DAY) + '_genes.txt')
        sd.save_array(norm,'./Results/IntermediateFiles/GC_day' + str(DAY) + '_norm.csv')
        sd.save_array(log,'./Results/IntermediateFiles/GC_day' + str(DAY) + '_log.csv')
        sd.save_array(sca,'./Results/IntermediateFiles/GC_day' + str(DAY) + '_sca.csv')
        sd.save_array(A,'./Results/IntermediateFiles/GC_day' + str(DAY) + '.csv')
        print("Here I am 19b")
        
        i = i + 1
        
if GetFromFileTransformedData == 1:

    A_0_genes = sd.get_from_txt('./Results/IntermediateFiles/G2019S_day0_genes.txt')[0]
    A_0_norm = sd.get_from_csv('./Results/IntermediateFiles/G2019S_day0_norm.csv',';')
    A_0_log = sd.get_from_csv('./Results/IntermediateFiles/G2019S_day0_log.csv',';')
    A_0_sca = sd.get_from_csv('./Results/IntermediateFiles/G2019S_day0_sca.csv',';')
    A_0 = sd.get_from_csv('./Results/IntermediateFiles/G2019S_day0.csv',';')
    print("Here I am 20a")
    A_10_genes = sd.get_from_txt('./Results/IntermediateFiles/G2019S_day10_genes.txt')[0]
    A_10_norm = sd.get_from_csv('./Results/IntermediateFiles/G2019S_day10_norm.csv',';')
    A_10_log = sd.get_from_csv('./Results/IntermediateFiles/G2019S_day10_log.csv',';')
    A_10_sca = sd.get_from_csv('./Results/IntermediateFiles/G2019S_day10_sca.csv',';')
    A_10 = sd.get_from_csv('./Results/IntermediateFiles/G2019S_day10.csv',';')
    print("Here I am 20b")    
    A_14_genes = sd.get_from_txt('./Results/IntermediateFiles/G2019S_day14_genes.txt')[0]
    A_14_norm = sd.get_from_csv('./Results/IntermediateFiles/G2019S_day14_norm.csv',';')
    A_14_log = sd.get_from_csv('./Results/IntermediateFiles/G2019S_day14_log.csv',';')
    A_14_sca = sd.get_from_csv('./Results/IntermediateFiles/G2019S_day14_sca.csv',';')
    A_14 = sd.get_from_csv('./Results/IntermediateFiles/G2019S_day14.csv',';')
    print("Here I am 20c")    
    A_42_genes = sd.get_from_txt('./Results/IntermediateFiles/G2019S_day42_genes.txt')[0]
    A_42_norm = sd.get_from_csv('./Results/IntermediateFiles/G2019S_day42_norm.csv',';')
    A_42_log = sd.get_from_csv('./Results/IntermediateFiles/G2019S_day42_log.csv',';')
    A_42_sca = sd.get_from_csv('./Results/IntermediateFiles/G2019S_day42_sca.csv',';')
    A_42 = sd.get_from_csv('./Results/IntermediateFiles/G2019S_day42.csv',';')
    print("Here I am 20d")
    
    B_0_genes = sd.get_from_txt('./Results/IntermediateFiles/GC_day0_genes.txt')[0]
    B_0_norm = sd.get_from_csv('./Results/IntermediateFiles/GC_day0_norm.csv',';')
    B_0_log = sd.get_from_csv('./Results/IntermediateFiles/GC_day0_log.csv',';')
    B_0_sca = sd.get_from_csv('./Results/IntermediateFiles/GC_day0_sca.csv',';')
    B_0 = sd.get_from_csv('./Results/IntermediateFiles/GC_day0.csv',';')
    print("Here I am 20a")
    B_10_genes = sd.get_from_txt('./Results/IntermediateFiles/GC_day10_genes.txt')[0]
    B_10_norm = sd.get_from_csv('./Results/IntermediateFiles/GC_day10_norm.csv',';')
    B_10_log = sd.get_from_csv('./Results/IntermediateFiles/GC_day10_log.csv',';')
    B_10_sca = sd.get_from_csv('./Results/IntermediateFiles/GC_day10_sca.csv',';')
    B_10 = sd.get_from_csv('./Results/IntermediateFiles/GC_day10.csv',';')
    print("Here I am 20b")    
    B_14_genes = sd.get_from_txt('./Results/IntermediateFiles/GC_day14_genes.txt')[0]
    B_14_norm = sd.get_from_csv('./Results/IntermediateFiles/GC_day14_norm.csv',';')
    B_14_log = sd.get_from_csv('./Results/IntermediateFiles/GC_day14_log.csv',';')
    B_14_sca = sd.get_from_csv('./Results/IntermediateFiles/GC_day14_sca.csv',';')
    B_14 = sd.get_from_csv('./Results/IntermediateFiles/GC_day14.csv',';')
    print("Here I am 20c")    
    B_42_genes = sd.get_from_txt('./Results/IntermediateFiles/GC_day42_genes.txt')[0]
    B_42_norm = sd.get_from_csv('./Results/IntermediateFiles/GC_day42_norm.csv',';')
    B_42_log = sd.get_from_csv('./Results/IntermediateFiles/GC_day42_log.csv',';')
    B_42_sca = sd.get_from_csv('./Results/IntermediateFiles/GC_day42_sca.csv',';')
    B_42 = sd.get_from_csv('./Results/IntermediateFiles/GC_day42.csv',';')
    print("Here I am 20d")

if DimensionalityReductionColorStem == 1:
    auxAlist = [A_0_norm,A_10_norm,A_14_norm,A_42_norm]
    genesAlist = [A_0_genes,A_10_genes,A_14_genes,A_42_genes]
    auxBlist = [B_0_norm,B_10_norm,B_14_norm,B_42_norm]
    genesBlist = [B_0_genes,B_10_genes,B_14_genes,B_42_genes]
    
    i = 0
    for DAY in [0,10,14,42]:
        
        print("Here I am 12")
        #Dimensionality reduction, like in Tuorial_DR
        #SOX2 looks good, it's a stemnes marker
        aux = auxAlist[i] # A_10_norm
        genes = genesAlist[i]
        t1 = 'G2019S, day ' + str(DAY) + ' - PCA\npoints colored by expression of SOX2'
        t2 = 'G2019S, day ' + str(DAY) + ' - tSNE\npoints colored by expression of SOX2'
        SaveAs1 = './Results/ProducedPlots/1-Stemness/G2019S_day' + str(DAY) + '_PCA_SOX2.pdf'
        SaveAs2 = './Results/ProducedPlots/1-Stemness/G2019S_day' + str(DAY) + '_tSNE_SOX2.pdf'
        a,b = get_DR(aux,genes.index('SOX2'),2,[],t1,t2,SaveAs1,SaveAs2) # STE
        
        aux = auxBlist[i] # B_42_norm ERROR BY TOMASZ??!?!?!!?
        genes = genesBlist[i]
        t1 = 'GC, day ' + str(DAY) + ' -PCA\npoints colored by expression of SOX2'
        t2 = 'GC, day ' + str(DAY) + ' - tSNE\npoints colored by expression of SOX2'
        SaveAs1 = './Results/ProducedPlots/1-Stemness/GC_day' + str(DAY) + '_PCA_SOX2.pdf'
        SaveAs2 = './Results/ProducedPlots/1-Stemness/GC_day' + str(DAY) + '_tSNE_SOX2.pdf'
        a,b = get_DR(aux,genes.index('SOX2'),2,[],t1,t2,SaveAs1,SaveAs2)
        del(a,b)
        print("Here I am 13")
        
        i = i+1
    
    ### print("Here I am 14")
    #### Here, Iwas just playing around,
    #### Let's keep it here, it may be useful
    ###C = Compare((A_10,A_10_genes),(B_10,B_10_genes))
    ###aux = C.merge()[0]
    ###C_norm, C_log, C_sca = sd.get_normalization(aux)
    ###a,b = get_DR(C_sca,C.c_names.index('MAP2'))
    ### a,b = get_DR_bis(C_sca[:,0:250],C_sca[:,250:500])
    ### a,b = get_DR_bis(C_norm[:,0:250],C_norm[:,250:500])
    ### print("Here I am 15")

if DRwithList_stemGenesOnly == 1:
    print("Here I am 16")
    #DR with genes from list_stem ONLY
    L_stem=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_stem.txt')[0] #L_kamil stands for listr from Kamil
    
    scaAlist = [A_0_sca,A_10_sca,A_14_sca,A_42_sca]
    genesAlist = [A_0_genes,A_10_genes,A_14_genes,A_42_genes]
    scaBlist = [B_0_sca,B_10_sca,B_14_sca,B_42_sca]
    genesBlist = [B_0_genes,B_10_genes,B_14_genes,B_42_genes]
                     
    i = 0
    for DAY in [0,10,14,42]:
        
        sca = scaAlist[i]
        #print(len(sca))
        genes = genesAlist[i]
        #print(len(genes))
        M = sd.get_submatrix(L_stem, sca, genes)
        t1 = 'G2019S,' + str(DAY) + ',sca - PCA\npoints colored by expression of SOX2'
        t2 = 'G2019S,' + str(DAY) + ',sca - tSNE\npoints colored by expression of SOX2'
        SaveAs1 = './Results/ProducedPlots/1-Stemness/G2019S_day' + str(DAY) + '_PCA_StemGenesOnly_SOX2.pdf'
        SaveAs2 = './Results/ProducedPlots/1-Stemness/G2019S_day' + str(DAY) + '_tSNE_StemGenesOnly_SOX2.pdf'
        y_pca, y_tsne = get_DR(M,genes.index('SOX2'),2,[],t1,t2,SaveAs1,SaveAs2)
        print("Here I am 17a")
        
        sca = scaBlist[i]
        genes = genesBlist[i]
        M = sd.get_submatrix(L_stem, sca, genes)
        t1 = 'GC' + str(DAY) + ',sca - PCA\npoints colored by expression of SOX2'
        t2 = 'GC' + str(DAY) + ',sca - tSNE\npoints colored by expression of SOX2'
        SaveAs1 = './Results/ProducedPlots/1-Stemness/GC_day' + str(DAY) + '_PCA_StemGenesOnly_SOX2.pdf'
        SaveAs2 = './Results/ProducedPlots/1-Stemness/GC_day' + str(DAY) + '_tSNE_StemGenesOnly_SOX2.pdf'
        y_pca, y_tsne = get_DR(M,genes.index('SOX2'),2,[],t1,t2,SaveAs1,SaveAs2)
        print("Here I am 17b")
        
        i = i+1

if DRwithList_CellCycleGenesOnly == 1:

    print("Here I am 21")
    #DR with genes from list_CellCycleTomaszSte ONLY
    L_CellCycle=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_CellCycle.txt')[0] #L_kamil stands for listr from Kamil
    print(L_CellCycle)
    scaAlist = [A_0_sca,A_10_sca,A_14_sca,A_42_sca]
    genesAlist = [A_0_genes,A_10_genes,A_14_genes,A_42_genes]
    scaBlist = [B_0_sca,B_10_sca,B_14_sca,B_42_sca]
    genesBlist = [B_0_genes,B_10_genes,B_14_genes,B_42_genes]
    
    i = 0
    for DAY in [0,10,14,42]: 
        
        sca = scaAlist[i]
        genes = genesAlist[i]
        M = sd.get_submatrix(L_CellCycle, sca, genes)
        t1 = 'G2019S day' + str(DAY) + ', sca - PCA\npoints colored by expression of TOP2A'
        t2 = 'G2019S day' + str(DAY) + ', sca - tSNE\npoints colored by expression of TOP2A'
        SaveAs1 = './Results/ProducedPlots/2-CellCycle/G2019S_day' + str(DAY) + '_PCA_CellCycleGenesOnly_TOP2A.pdf'
        SaveAs2 = './Results/ProducedPlots/2-CellCycle/G2019S_day' + str(DAY) + '_tSNE_CellCycleGenesOnly_TOP2A.pdf'
        y_pca, y_tsne = get_DR(M,genes.index('TOP2A'),2,[],t1,t2,SaveAs1,SaveAs2)
        print("Here I am 22")
        
        sca = scaBlist[i]
        genes = genesBlist[i]                   
        M = sd.get_submatrix(L_CellCycle, sca, genes)
        t1 = 'GC day' + str(DAY) + ', sca - PCA\npoints colored by expression of TOP2A'
        t2 = 'GC day' + str(DAY) + ', sca - tSNE\npoints colored by expression of TOP2A'
        SaveAs1 = './Results/ProducedPlots/2-CellCycle/GC_day' + str(DAY) + '_PCA_CellCycleGenesOnly_TOP2A.pdf'
        SaveAs2 = './Results/ProducedPlots/2-CellCycle/GC_day' + str(DAY) + '_tSNE_CellCycleGenesOnly_TOP2A.pdf'
        y_pca, y_tsne = get_DR(M,genes.index('TOP2A'),2,[],t1,t2,SaveAs1,SaveAs2)
        print("Here I am 23")

        i = i+1

if DRwithList_MitoGenesOnly == 1:
    print("Here I am 16")
    #DR with genes from list_stem ONLY
    L_Mito=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Mito.txt')[0] #L_kamil stands for listr from Kamil
    
    scaAlist = [A_0_sca,A_10_sca,A_14_sca,A_42_sca]
    genesAlist = [A_0_genes,A_10_genes,A_14_genes,A_42_genes]
    scaBlist = [B_0_sca,B_10_sca,B_14_sca,B_42_sca]
    genesBlist = [B_0_genes,B_10_genes,B_14_genes,B_42_genes]
                     
    i = 0
    for DAY in [0,10,14,42]:
        
        sca = scaAlist[i]
        #print(len(sca))
        genes = genesAlist[i]
        #print(len(genes))
        M = sd.get_submatrix(L_Mito, sca, genes)
        t1 = 'G2019S,' + str(DAY) + ',sca - PCA\npoints colored by expression of CYC1'
        t2 = 'G2019S,' + str(DAY) + ',sca - tSNE\npoints colored by expression of CYC1'
        SaveAs1 = './Results/ProducedPlots/3-Mito/G2019S_day' + str(DAY) + '_PCA_MitoGenesOnly_CYC1.pdf'
        SaveAs2 = './Results/ProducedPlots/3-Mito/G2019S_day' + str(DAY) + '_tSNE_MitoGenesOnly_CYC1.pdf'
        y_pca, y_tsne = get_DR(M,genes.index('CYC1'),2,[],t1,t2,SaveAs1,SaveAs2)
        print("Here I am 17a")
        
        sca = scaBlist[i]
        genes = genesBlist[i]
        M = sd.get_submatrix(L_Mito, sca, genes)
        t1 = 'GC' + str(DAY) + ',sca - PCA\npoints colored by expression of CYC1'
        t2 = 'GC' + str(DAY) + ',sca - tSNE\npoints colored by expression of CYC1'
        SaveAs1 = './Results/ProducedPlots/3-Mito/GC_day' + str(DAY) + '_PCA_MitoGenesOnly_CYC1.pdf'
        SaveAs2 = './Results/ProducedPlots/3-Mito/GC_day' + str(DAY) + '_tSNE_MitoGenesOnly_CYC1.pdf'
        y_pca, y_tsne = get_DR(M,genes.index('CYC1'),2,[],t1,t2,SaveAs1,SaveAs2)
        print("Here I am 17b")
        
        i = i+1
        
if DRwithList_ROSGenesOnly == 1:
    print("Here I am 16")
    #DR with genes from list_stem ONLY
    L_ROS=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_ROS.txt')[0] #L_kamil stands for listr from Kamil
    
    scaAlist = [A_0_sca,A_10_sca,A_14_sca,A_42_sca]
    genesAlist = [A_0_genes,A_10_genes,A_14_genes,A_42_genes]
    scaBlist = [B_0_sca,B_10_sca,B_14_sca,B_42_sca]
    genesBlist = [B_0_genes,B_10_genes,B_14_genes,B_42_genes]
                     
    i = 0
    for DAY in [0,10,14,42]:
        
        sca = scaAlist[i]
        #print(len(sca))
        genes = genesAlist[i]
        #print(len(genes))
        M = sd.get_submatrix(L_ROS, sca, genes)
        t1 = 'G2019S,' + str(DAY) + ',sca - PCA\npoints colored by expression of ANXA6'
        t2 = 'G2019S,' + str(DAY) + ',sca - tSNE\npoints colored by expression of ANXA6'
        SaveAs1 = './Results/ProducedPlots/4-ROS/G2019S_day' + str(DAY) + '_PCA_ROSGenesOnly_ANXA6.pdf'
        SaveAs2 = './Results/ProducedPlots/4-ROS/G2019S_day' + str(DAY) + '_tSNE_ROSGenesOnly_ANXA6.pdf'
        y_pca, y_tsne = get_DR(M,genes.index('ANXA6'),2,[],t1,t2,SaveAs1,SaveAs2)
        print("Here I am 17a")
        
        sca = scaBlist[i]
        genes = genesBlist[i]
        M = sd.get_submatrix(L_ROS, sca, genes)
        t1 = 'GC' + str(DAY) + ',sca - PCA\npoints colored by expression of ANXA6'
        t2 = 'GC' + str(DAY) + ',sca - tSNE\npoints colored by expression of ANXA6'
        SaveAs1 = './Results/ProducedPlots/4-ROS/GC_day' + str(DAY) + '_PCA_ROSGenesOnly_ANXA6.pdf'
        SaveAs2 = './Results/ProducedPlots/4-ROS/GC_day' + str(DAY) + '_tSNE_ROSGenesOnly_ANXA6.pdf'
        y_pca, y_tsne = get_DR(M,genes.index('ANXA6'),2,[],t1,t2,SaveAs1,SaveAs2)
        print("Here I am 17b")
        
        i = i+1
        
if DoComparisonPreliminarySteps == 1:
    print("Doing preliminary steps for comparison...")
    A_0_250_genes, A_0_250_expr, A_0_250_expr_np = sd.get_from_txt('./Data_SingleCellRNASeq_DigitalExpressionMatrices/A_0.txt')
    B_0_250_genes, B_0_250_expr, B_0_250_expr_np = sd.get_from_txt('./Data_SingleCellRNASeq_DigitalExpressionMatrices/B_0.txt')  
    A_10_250_genes, A_10_250_expr, A_10_250_expr_np = sd.get_from_txt('./Data_SingleCellRNASeq_DigitalExpressionMatrices/A_10.txt')
    B_10_250_genes, B_10_250_expr, B_10_250_expr_np = sd.get_from_txt('./Data_SingleCellRNASeq_DigitalExpressionMatrices/B_10.txt')
    A_14_250_genes, A_14_250_expr, A_14_250_expr_np = sd.get_from_txt('./Data_SingleCellRNASeq_DigitalExpressionMatrices/A_14.txt')
    B_14_250_genes, B_14_250_expr, B_14_250_expr_np = sd.get_from_txt('./Data_SingleCellRNASeq_DigitalExpressionMatrices/B_14.txt')    
    A_42_250_genes, A_42_250_expr, A_42_250_expr_np = sd.get_from_txt('./Data_SingleCellRNASeq_DigitalExpressionMatrices/A_42.txt')
    B_42_250_genes, B_42_250_expr, B_42_250_expr_np = sd.get_from_txt('./Data_SingleCellRNASeq_DigitalExpressionMatrices/B_42.txt')
    
    print("Normalizing Day 0...")
    A_0_250 = sd.get_top_cells_matrix(A_0_250_expr_np,250)[0]
    B_0_250 = sd.get_top_cells_matrix(B_0_250_expr_np,250)[0]
    A_0_250_norm, A_0_250_log, A_0_250_sca = sd.get_normalization(A_0_250)
    B_0_250_norm, B_0_250_log, B_0_250_sca = sd.get_normalization(B_0_250)
    
    print("Normalizing Day 10...")
    A_10_250 = sd.get_top_cells_matrix(A_10_250_expr_np,250)[0]
    B_10_250 = sd.get_top_cells_matrix(B_10_250_expr_np,250)[0]
    A_10_250_norm, A_10_250_log, A_10_250_sca = sd.get_normalization(A_10_250)
    B_10_250_norm, B_10_250_log, B_10_250_sca = sd.get_normalization(B_10_250)

    print("Normalizing Day 14...") 
    A_14_250 = sd.get_top_cells_matrix(A_14_250_expr_np,250)[0]
    B_14_250 = sd.get_top_cells_matrix(B_14_250_expr_np,250)[0]
    A_14_250_norm, A_14_250_log, A_14_250_sca = sd.get_normalization(A_14_250)
    B_14_250_norm, B_14_250_log, B_14_250_sca = sd.get_normalization(B_14_250)

    print("Normalizing Day 42...")
    A_42_250 = sd.get_top_cells_matrix(A_42_250_expr_np,250)[0]
    B_42_250 = sd.get_top_cells_matrix(B_42_250_expr_np,250)[0]
    A_42_250_norm, A_42_250_log, A_42_250_sca = sd.get_normalization(A_42_250)
    B_42_250_norm, B_42_250_log, B_42_250_sca = sd.get_normalization(B_42_250)

    print("Concatenating matrices...")
    C0 = Compare((A_0_250,A_0_250_genes),(B_0_250,B_0_250_genes)) 
    C10 = Compare((A_10_250,A_10_250_genes),(B_10_250,B_10_250_genes)) 
    C14 = Compare((A_14_250,A_14_250_genes),(B_14_250,B_14_250_genes)) 
    C42 = Compare((A_42_250,A_42_250_genes),(B_42_250,B_42_250_genes)) 

    print("Merging matrices...")
    aux0 = C0.merge()[0]
    aux10 = C10.merge()[0]
    aux14 = C14.merge()[0]
    aux42 = C42.merge()[0]

    print("Normalizing Comparison Matrices...")         
    C0_norm, C0_log, C0_sca = sd.get_normalization(aux0)
    C10_norm, C10_log, C10_sca = sd.get_normalization(aux10)
    C14_norm, C14_log, C14_sca = sd.get_normalization(aux14)
    C42_norm, C42_log, C42_sca = sd.get_normalization(aux42)
    print("Preliminaries for Comparison done.")

if DoComparison == 1:
    print("Doing actual comparison.")

    normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
    scaClist = [C0_sca,C10_sca,C14_sca,C42_sca]
    Clist = [C0,C10,C14,C42]
    genesClist = []
                     
    i = 0
    for DAY in [0,10,14,42]:
        print("Starting Day " + str(DAY) + "...")
        C_norm = normClist[i]
        C_sca = scaClist[i]
        C = Clist[i]
        
        #print("PCA/tSNE colored wih MAP2")
        #Myt1 = 'COMPARISON,' + str(DAY) + ',sca - PCA\npoints colored by expression of MAP2'
        #Myt2 = 'COMPARISON,' + str(DAY) + ',sca - tSNE\npoints colored by expression of MAP2'
        #MySaveAs1 = './Results/COMPARISON_day' + str(DAY) + '_PCA_MAP2.pdf'
        #MySaveAs2 = './Results/COMPARISON_day' + str(DAY) + '_tSNE_MAP2.pdf'
        #a,b = get_DR(C_sca,C.c_names.index('MAP2'),2,[],Myt1,Myt2,MySaveAs1,MySaveAs2)
    
        print("PCA/tSNE with sca normalizaion")
        Myt1 = 'COMPARISON,' + str(DAY) + ',sca - PCA'
        Myt2 = 'COMPARISON,' + str(DAY) + ',sca - tSNE'
        MySaveAs1 = './Results/COMPARISON_day' + str(DAY) + '_PCA_sca.pdf'
        MySaveAs2 = './Results/COMPARISON_day' + str(DAY) + '_tSNE_sca.pdf'
        a,b = get_DR_bis(C_sca[:,0:250],C_sca[:,250:500],t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)
    
        print("PCA/tSNE with norm normalizaion")
        Myt1 = 'COMPARISON,' + str(DAY) + ',norm - PCA'
        Myt2 = 'COMPARISON,' + str(DAY) + ',norm - tSNE'
        MySaveAs1 = './Results/COMPARISON_day' + str(DAY) + '_PCA_norm.pdf'
        MySaveAs2 = './Results/COMPARISON_day' + str(DAY) + '_tSNE_norm.pdf'
        a,b = get_DR_bis(C_norm[:,0:250],C_norm[:,250:500],t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)
        
        print("Other steps comparison...")

        L_stem=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_stem.txt')[0]
        
        C_genes = C.comm_genes()
        M = sd.get_submatrix(L_stem, C_sca, C_genes)
        Myt1 = 'COMPARISON,' + str(DAY) + ',sca, ONLY STEM GENES - PCA'
        Myt2 = 'COMPARISON,' + str(DAY) + ',sca, ONLY STEM GENES- tSNE'
        MySaveAs1 = './Results/COMPARISON_day' + str(DAY) + '_PCA_sca_STEMgenes.pdf'
        MySaveAs2 = './Results/COMPARISON_day' + str(DAY) + '_tSNE_sca_STEMgenes.pdf'
        a,b = get_DR_bis(M[:,0:250],M[:,250:500],t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)
        
        #        print("Comp HERE")
        #        Myt1 = 'COMPARISON,' + str(DAY) + ',sca, ONLY STEM GENES - PCA\npoints colored by expression of MAP2'
        #        Myt2 = 'COMPARISON,' + str(DAY) + ',sca, ONLY STEM GENES- tSNE\npoints colored by expression of MAP2'
        #        MySaveAs1 = './Results/ProducedPlots/5-Compare/COMPARISON_day' + str(DAY) + '_PCA_sca_STEMgenes_MAP2.pdf'
        #        MySaveAs2 = './Results/ProducedPlots/5-Compare/COMPARISON_day' + str(DAY) + '_tSNE_sca_STEMgenes_MAP2.pdf'
        #        a,b = get_DR(M,C.c_names.index('MAP2'),2,[],Myt1,Myt2,MySaveAs1,MySaveAs2)
        #        print("Comp THERE")
        
        L_CellCycle=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_CellCycle.txt')[0]

        C_genes = C.comm_genes()
        M = sd.get_submatrix(L_CellCycle, C_sca, C_genes)
        Myt1 = 'COMPARISON,' + str(DAY) + ',sca, ONLY CELLC CYCLE GENES - PCA'
        Myt2 = 'COMPARISON,' + str(DAY) + ',sca, ONLY CELLC CYCLE - tSNE'
        MySaveAs1 = './Results/COMPARISON_day' + str(DAY) + '_PCA_sca_CELL_CYCLE_genes.pdf'
        MySaveAs2 = './Results/COMPARISON_day' + str(DAY) + '_tSNE_sca_CELL_CYCLE_genes.pdf'
        a,b = get_DR_bis(M[:,0:250],M[:,250:500],t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)
        
        L_Mito=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Mito.txt')[0]
        
        C_genes = C.comm_genes()
        M = sd.get_submatrix(L_Mito, C_sca, C_genes)
        Myt1 = 'COMPARISON,' + str(DAY) + ',sca, ONLY MITO GENES - PCA'
        Myt2 = 'COMPARISON,' + str(DAY) + ',sca, ONLY MITO GENES- tSNE'
        MySaveAs1 = './Results/COMPARISON_day' + str(DAY) + '_PCA_sca_MITOgenes.pdf'
        MySaveAs2 = './Results/COMPARISON_day' + str(DAY) + '_tSNE_sca_MITOgenes.pdf'
        a,b = get_DR_bis(M[:,0:250],M[:,250:500],t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)
        
        L_ROS=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_ROS.txt')[0]
        
        C_genes = C.comm_genes()
        M = sd.get_submatrix(L_ROS, C_sca, C_genes)
        Myt1 = 'COMPARISON,' + str(DAY) + ',sca, ONLY ROS GENES - PCA'
        Myt2 = 'COMPARISON,' + str(DAY) + ',sca, ONLY ROS GENES- tSNE'
        MySaveAs1 = './Results/COMPARISON_day' + str(DAY) + '_PCA_sca_ROSgenes.pdf'
        MySaveAs2 = './Results/COMPARISON_day' + str(DAY) + '_tSNE_sca_ROSgenes.pdf'
        a,b = get_DR_bis(M[:,0:250],M[:,250:500],t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)
    
        i = i+1


if ProduceListOfDifferentiallyExpressedGenes == 1:
    #for pval in []:
    Clist = [C0,C10,C14,C42]
    i = 0
    for DAY in [0,10,14,42]:
        C = Clist[i]
        
        print("Compute Differentially Expressed Genes...")
        #Saves list of genes that are different between A and B on a given day.
        aux=C.avg() #Here, choose the day
        a = np.array(aux) #some transformations
        a = np.transpose(a) #some transformations
        L_min, L_min_stats = C.avg_list(0.01/20000, aux[3], a) #you need to change this
        sd.save_array(L_min_stats, './Results/TabS2_DiffExprGenes_day' + str(DAY) + '.csv') #and here
        
        #L_min_pval005, L_min_stats_pval005 = C.avg_list(0.05, aux[3], a) #you need to change this
        #sd.save_array(L_min_stats_pval005, './Results/jonas_DiffExprGenes_day' + str(DAY) + '_pval005.csv') #and here
        
        print("Compare DEGs with Lists of Relevant Genes...this gives the numbers needed to compute the percentages of FigS4B!!!")
        #produces list of genes common for the two lists
        ListOfDifferentiallyExpressedSTEMNESSGenes = []
        ListOfDifferentiallyExpressedCELLCYCLEGenes = []
        ListOfDifferentiallyExpressedSMITOGenes = []
        ListOfDifferentiallyExpressedROSGenes = []
        ListOfDifferentiallyExpressedD0BulkAnalysisGenes = []
        ListOfDifferentiallyExpressedDOPAMINNEURONSGenes = []
        ListOfDifferentiallyExpressedProApoptotic = []
        ListOfDifferentiallyExpressedAntiApoptotic = []
        ListOfDifferentiallyExpressedCaspases = []

        L_stem=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_stem.txt')[0]
        L_CellCycle=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_CellCycle.txt')[0]
        L_Mito=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Mito.txt')[0]
        L_ROS=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_ROS.txt')[0]
        L_D0BulkAnalysis=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/listGenes20degD0BulkAnalysisEnrico.txt')[0]
        L_DopaminNeurons=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_DopaminergicNeurons.txt')[0]
              
        L_ProApoptotic = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Pro-Apoptotic.txt')[0]
        L_AntiApoptotic = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Anti-Apoptotic.txt')[0]
        L_Caspases = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Caspases.txt')[0]

        for GeneName in L_min:
        
            if GeneName in L_stem:
                ListOfDifferentiallyExpressedSTEMNESSGenes.append(GeneName)
            if GeneName in L_CellCycle:
                ListOfDifferentiallyExpressedCELLCYCLEGenes.append(GeneName)
            if GeneName in L_Mito:
                ListOfDifferentiallyExpressedSMITOGenes.append(GeneName)
            if GeneName in L_ROS:
                ListOfDifferentiallyExpressedROSGenes.append(GeneName)
            if GeneName in L_D0BulkAnalysis:
                ListOfDifferentiallyExpressedD0BulkAnalysisGenes.append(GeneName)
            if GeneName in L_DopaminNeurons:
                ListOfDifferentiallyExpressedDOPAMINNEURONSGenes.append(GeneName)
            if GeneName in L_ProApoptotic:
                ListOfDifferentiallyExpressedProApoptotic.append(GeneName)            
            if GeneName in L_AntiApoptotic:
                ListOfDifferentiallyExpressedAntiApoptotic.append(GeneName)            
            if GeneName in L_Caspases:
                ListOfDifferentiallyExpressedCaspases.append(GeneName)
                
        # sd.save_array(ListOfDifferentiallyExpressedSTEMNESSGenes, './Results/jonas_DiffExprGenes_STEMNESS_day' + str(DAY) + '.csv') #and here
        # sd.save_array(ListOfDifferentiallyExpressedCELLCYCLEGenes, './Results/jonas_DiffExprGenes_CELLCYCLE_day' + str(DAY) + '.csv') #and here
        # sd.save_array(ListOfDifferentiallyExpressedSMITOGenes, './Results/jonas_DiffExprGenes_MITO_day' + str(DAY) + '.csv') #and here
        # sd.save_array(ListOfDifferentiallyExpressedROSGenes, './Results/jonas_DiffExprGenes_ROS_day' + str(DAY) + '.csv') #and here

        sd.save_from_list(ListOfDifferentiallyExpressedSTEMNESSGenes, './Results/TabS2_DiffExprGenes_STEMNESS_day' + str(DAY) + '.txt') #and here
        sd.save_from_list(ListOfDifferentiallyExpressedCELLCYCLEGenes, './Results/TabS2_DiffExprGenes_CELLCYCLE_day' + str(DAY) + '.txt') #and here
        sd.save_from_list(ListOfDifferentiallyExpressedSMITOGenes, './Results/TabS2_DiffExprGenes_MITO_day' + str(DAY) + '.txt') #and here
        sd.save_from_list(ListOfDifferentiallyExpressedROSGenes, './Results/TabS2_DiffExprGenes_ROS_day' + str(DAY) + '.txt') #and here
        sd.save_from_list(ListOfDifferentiallyExpressedD0BulkAnalysisGenes, './Results/TabS2_DiffExprGenes_D0BulkAnalysis_day' + str(DAY) + '.txt') #and here
        sd.save_from_list(ListOfDifferentiallyExpressedDOPAMINNEURONSGenes, './Results/TabS2_DiffExprGenes_DopaminergicNeurons_day' + str(DAY) + '.txt') #and here
        sd.save_from_list(ListOfDifferentiallyExpressedProApoptotic, './Results/TabS2_DiffExprGenes_ProApoptotic_day' + str(DAY) + '.txt') #and here
        sd.save_from_list(ListOfDifferentiallyExpressedAntiApoptotic, './Results/TabS2_DiffExprGenes_AntiApoptotic_day' + str(DAY) + '.txt') #and here
        sd.save_from_list(ListOfDifferentiallyExpressedCaspases, './Results/TabS2_DiffExprGenes_Caspases_day' + str(DAY) + '.txt') #and here

        print(len(ListOfDifferentiallyExpressedSTEMNESSGenes))
        print(len(ListOfDifferentiallyExpressedCELLCYCLEGenes))
        print(len(ListOfDifferentiallyExpressedDOPAMINNEURONSGenes))
        print(len(ListOfDifferentiallyExpressedSMITOGenes))
        print(len(ListOfDifferentiallyExpressedROSGenes))
        print(len(ListOfDifferentiallyExpressedD0BulkAnalysisGenes))
        print(len(ListOfDifferentiallyExpressedProApoptotic))
        print(len(ListOfDifferentiallyExpressedAntiApoptotic))
        print(len(ListOfDifferentiallyExpressedCaspases))

        i = i+1
        
if PlotHistogramsIndividualGenes == 1:
    
    Clist = [C0,C10,C14,C42]    
    normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
        
    i = 0
    for DAY in [0,10,14,42]:
        C = Clist[i]
        C_norm = normClist[i]
        
#        plt.plot(C_sca[C.c_names.index('NR2F1'),1:250])
#        plt.plot(C_sca[C.c_names.index('NR2F1'),251:500])
#
#        plt.plot(C_sca[C.c_names.index('TOP2A'),1:250])
#        plt.plot(C_sca[C.c_names.index('TOP2A'),251:500])
#        
#        plt.plot(C_sca[C.c_names.index('CHCHD2'),1:250])
#        plt.plot(C_sca[C.c_names.index('CHCHD2'),251:500])
        
        # PREVIOUSLY NR2F1
        fig1 = plt.figure()
        plt.hist(C_norm[C.c_names.index('TH'),0:250],bins=20, histtype='stepfilled', color='b', label='G2019S', alpha=0.5)        
        plt.hist(C_norm[C.c_names.index('TH'),250:500],bins=20, histtype='stepfilled', color='r', label='GC', alpha=0.5)
        plt.title("Differential expression of TH")
        plt.xlabel("Gene expression (norm)")
        plt.ylabel("Number of Cells")
        plt.legend()
        plt.show()
        plt.savefig('./Results/jonas_DiffExprTH_day' + str(DAY) + '.pdf')
        
        fig2 = plt.figure()
        plt.hist(C_norm[C.c_names.index('TOP2A'),10:250],bins=20, histtype='stepfilled',  color='b', label='G2019S', alpha=0.5)        
        plt.hist(C_norm[C.c_names.index('TOP2A'),250:500],bins=20, histtype='stepfilled', color='r', label='GC', alpha=0.5)
        plt.title("Differential expression of TOP2A")
        plt.xlabel("Gene expression (norm)")
        plt.ylabel("Number of Cells")
        plt.legend()
        plt.show()
        plt.savefig('./Results/jonas_DiffExprTOP2A_day' + str(DAY) + '.pdf')

        fig3 = plt.figure()
        plt.hist(C_norm[C.c_names.index('CHCHD2'),0:250],bins=20, histtype='stepfilled', color='b', label='G2019S', alpha=0.5)        
        plt.hist(C_norm[C.c_names.index('CHCHD2'),250:500],bins=20, histtype='stepfilled', color='r', label='GC', alpha=0.5)
        plt.title("Differential expression of CHCHD2")
        plt.xlabel("Gene expression (norm)")
        plt.ylabel("Number of Cells")
        plt.legend()
        plt.show()
        plt.savefig('./Results/jonas_DiffExprCHCHD2_day' + str(DAY) + '.pdf')
        
        fig4 = plt.figure()
        plt.hist(C_norm[C.c_names.index('EN1'),0:250],bins=20, histtype='stepfilled', color='b', label='G2019S', alpha=0.5)        
        plt.hist(C_norm[C.c_names.index('EN1'),250:500],bins=20, histtype='stepfilled', color='r', label='GC', alpha=0.5)
        plt.title("Differential expression of EN1")
        plt.xlabel("Gene expression (norm)")
        plt.ylabel("Number of Cells")
        plt.legend()
        plt.show()
        plt.savefig('./Results/jonas_DiffExprEN1_day' + str(DAY) + '.pdf')
        
if PlotHistogramsLRRK2 == 1:
    
    Clist = [C0,C10,C14,C42]    
    normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
    
#        fig1 = plt.figure()
#        plt.hist(C_norm[C.c_names.index('LRRK2'),0:250],bins=20, histtype='stepfilled', color='r', label='G2019S', alpha=0.5)        
#        plt.hist(C_norm[C.c_names.index('LRRK2'),250:500],bins=20, histtype='stepfilled', color='gray', label='GC', alpha=0.5)
#        plt.title("Differential expression of LRRK2")
#        plt.xlabel("Gene expression (norm)")
#        plt.ylabel("Number of Cells")
#        plt.legend()
#        plt.show()
#        plt.savefig('./Results/jonas_DiffExprLRRK2_day' + str(DAY) + '.pdf')
#        
    Nlines = 2
    Ncols = 4
    fig, axes = plt.subplots(nrows=Nlines, ncols=Ncols)
    fig.subplots_adjust(hspace=0.6, wspace=0.15)

    for ax in axes.flat:
        # Hide all ticks and labels
        #ax.xaxis.set_visible(False)
        #ax.yaxis.set_visible(False)

        # Set up ticks only on one side for the "edge" subplots...
        if ax.is_first_col():
            ax.yaxis.set_ticks_position('left')
            ax.yaxis.set_visible(True)

    j = 0
    for DAY in [0,10,14,42]:
        print('Day ' + str(DAY))
        print()
        C = Clist[j]
        C_norm = normClist[j]
        
        #plt.title("Differential expression of LRRK2")

        axes[0,0].set_ylabel("Number of cells")
        axes[1,0].set_ylabel("Number of cells")
        axes[0,j].set_xlabel("LRRK2 expression (norm)")
        axes[1,j].set_xlabel("LRRK2 expression (norm)")

        axes[0,j].hist(C_norm[C.c_names.index('LRRK2'),0:250],bins=12, histtype='stepfilled', color=(255.0/256,0.,0.), label='PAR2-G2019S', alpha=0.5)   #  normed=0,     
        axes[0,j].hist(C_norm[C.c_names.index('LRRK2'),250:500],bins=12, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), label='PAR2-GC', alpha=0.5) #  normed=0,

        axes[1,j].hist(C_norm[C.c_names.index('LRRK2'),0:250],bins=12, histtype='stepfilled', color=(255.0/256,0.,0.), label='PAR2-G2019S', alpha=0.5)      #  normed=0,  
        axes[1,j].hist(C_norm[C.c_names.index('LRRK2'),250:500],bins=12, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), label='PAR2-GC', alpha=0.5) #  normed=0,
       
        axes[0,0].set_xlim(0,0.5)
        axes[0,1].set_xlim(0,1.5)        
        axes[0,2].set_xlim(0,5.0)        
        axes[0,3].set_xlim(0,4.0)
        
        axes[0,0].set_ylim(0,260)
        axes[0,1].set_ylim(0,260)        
        axes[0,2].set_ylim(0,260)        
        axes[0,3].set_ylim(0,260)
        
        axes[1,0].set_xlim(0,0.5)
        axes[1,1].set_xlim(0,1.5)        
        axes[1,2].set_xlim(0,5.0)        
        axes[1,3].set_xlim(0,4.0)
        
        axes[1,0].set_ylim(0,11)
        axes[1,1].set_ylim(0,11)        
        axes[1,2].set_ylim(0,11)        
        axes[1,3].set_ylim(0,11)
        ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
        # of 2 populations are statistically significantly different ######
        
        #DataPopulation1 = numpy.array([1,2,3,4,5,2,3,4,3])
        DataPopulation1 = C_norm[C.c_names.index('LRRK2'),0:250]
        #DataPopulation2 = numpy.array([1,2,3.001,4,5,2,3,4,3])
        DataPopulation2 = C_norm[C.c_names.index('LRRK2'),250:500]
        DataPopulation1_Object = wstats.DescrStatsW(DataPopulation1)
        DataPopulation2_Object = wstats.DescrStatsW(DataPopulation2)
        MyComparison = wstats.CompareMeans(DataPopulation1_Object,DataPopulation2_Object)
        TestStatistics, pvalue = MyComparison.ztest_ind(alternative='two-sided', usevar='unequal', value=0)

        # Remember Bonferroni correction for multiple testing, 
        # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
        # by the number of repetitions of the test, 4 days x 7 lists = 28
        # The 7 lists are Stem, DAneur, Mito, Cell Cycle, Pro-apoptosis, Anti-apoptosis, Caspases
        print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
        print(' ')

        plt.legend()

        j=j+1
        plt.show()
        plt.savefig('./Results/jonas_DiffExprLRRK2.pdf')
        
if PlotHistogramsPAX6 == 1:
    
    Clist = [C0,C10,C14,C42]    
    normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
    
#        fig1 = plt.figure()
#        plt.hist(C_norm[C.c_names.index('LRRK2'),0:250],bins=20, histtype='stepfilled', color='r', label='G2019S', alpha=0.5)        
#        plt.hist(C_norm[C.c_names.index('LRRK2'),250:500],bins=20, histtype='stepfilled', color='gray', label='GC', alpha=0.5)
#        plt.title("Differential expression of LRRK2")
#        plt.xlabel("Gene expression (norm)")
#        plt.ylabel("Number of Cells")
#        plt.legend()
#        plt.show()
#        plt.savefig('./Results/jonas_DiffExprLRRK2_day' + str(DAY) + '.pdf')
#        
    Nlines = 2
    Ncols = 4
    fig, axes = plt.subplots(nrows=Nlines, ncols=Ncols)
    fig.subplots_adjust(hspace=0.6, wspace=0.15)

    for ax in axes.flat:
        # Hide all ticks and labels
        #ax.xaxis.set_visible(False)
        #ax.yaxis.set_visible(False)

        # Set up ticks only on one side for the "edge" subplots...
        if ax.is_first_col():
            ax.yaxis.set_ticks_position('left')
            ax.yaxis.set_visible(True)

    j = 0
    for DAY in [0,10,14,42]:
        print('Day ' + str(DAY))
        print()
        C = Clist[j]
        C_norm = normClist[j]
        
        #plt.title("Differential expression of LRRK2")

        axes[0,0].set_ylabel("Number of cells")
        axes[1,0].set_ylabel("Number of cells")
        axes[0,j].set_xlabel("PAX6 expression (norm. to max)")
        axes[1,j].set_xlabel("PAX6 expression (norm. to max)")

        MaxX = max(max(C_norm[C.c_names.index('PAX6'),0:250]),max(C_norm[C.c_names.index('PAX6'),250:500]))
        GenesExpressionA_NORM = C_norm[C.c_names.index('PAX6'),0:250] / MaxX
        GenesExpressionB_NORM = C_norm[C.c_names.index('PAX6'),250:500] / MaxX
        binwidth = 0.035
        axes[0,j].hist(GenesExpressionA_NORM, bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(255.0/256,0.,0.), label='PAR2-G2019S', alpha=0.5)    #  normed=0,    
        axes[0,j].hist(GenesExpressionB_NORM, bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), label='PAR2-GC', alpha=0.5) #  normed=0,

        axes[0,j].set_xlim(0,1.1)
        axes[0,j].set_ylim(0,250)

#        axes[0,j].hist(C_norm[C.c_names.index('PAX6'),0:250],bins=12, histtype='stepfilled', color=(255.0/256,0.,0.), normed=0, label='PAR2-G2019S', alpha=0.5)        
#        axes[0,j].hist(C_norm[C.c_names.index('PAX6'),250:500],bins=12, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), normed=0, label='PAR2-GC', alpha=0.5)

#        axes[1,j].hist(C_norm[C.c_names.index('PAX6'),0:250],bins=12, histtype='stepfilled', color=(255.0/256,0.,0.), normed=0, label='PAR2-G2019S', alpha=0.5)        
#        axes[1,j].hist(C_norm[C.c_names.index('PAX6'),250:500],bins=12, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), normed=0, label='PAR2-GC', alpha=0.5)
#               
        axes[1,j].hist(GenesExpressionA_NORM, bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(255.0/256,0.,0.), label='PAR2-G2019S', alpha=0.5)     #  normed=0,   
        axes[1,j].hist(GenesExpressionB_NORM, bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), label='PAR2-GC', alpha=0.5) #  normed=0, 

        axes[1,j].set_xlim(0,1.1)
        axes[1,j].set_ylim(0,30)
        
#        axes[0,0].set_xlim(0,0.5)
#        axes[0,1].set_xlim(0,1.5)        
#        axes[0,2].set_xlim(0,5.0)        
#        axes[0,3].set_xlim(0,4.0)
#        
#        axes[0,0].set_ylim(0,260)
#        axes[0,1].set_ylim(0,260)        
#        axes[0,2].set_ylim(0,260)        
#        axes[0,3].set_ylim(0,260)
#        
#        axes[1,0].set_xlim(0,0.5)
#        axes[1,1].set_xlim(0,1.5)        
#        axes[1,2].set_xlim(0,5.0)        
#        axes[1,3].set_xlim(0,4.0)
#        
#        axes[1,0].set_ylim(0,11)
#        axes[1,1].set_ylim(0,11)        
#        axes[1,2].set_ylim(0,11)        
#        axes[1,3].set_ylim(0,11)
        ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
        # of 2 populations are statistically significantly different ######
        
        import statsmodels.stats.weightstats as wstats
        #DataPopulation1 = numpy.array([1,2,3,4,5,2,3,4,3])
        DataPopulation1 = C_norm[C.c_names.index('PAX6'),0:250]
        #DataPopulation2 = numpy.array([1,2,3.001,4,5,2,3,4,3])
        DataPopulation2 = C_norm[C.c_names.index('PAX6'),250:500]
        DataPopulation1_Object = wstats.DescrStatsW(DataPopulation1)
        DataPopulation2_Object = wstats.DescrStatsW(DataPopulation2)
        MyComparison = wstats.CompareMeans(DataPopulation1_Object,DataPopulation2_Object)
        TestStatistics, pvalue = MyComparison.ztest_ind(alternative='two-sided', usevar='unequal', value=0)

        # Remember Bonferroni correction for multiple testing, 
        # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
        # by the number of repetitions of the test, 4 days x 7 lists = 28
        # The 7 lists are Stem, DAneur, Mito, Cell Cycle, Pro-apoptosis, Anti-apoptosis, Caspases
        print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
        print(' ')

        plt.legend()

        j=j+1
        plt.show()
        plt.savefig('./Results/jonas_DiffExprPAX6.pdf')
        
if PlotHistogramsCumulative == 1:

    print("Plotting histograms of cumulative gene expression scores for relevant lists of genes, namely stemness, DA neurons, cell cycle, mitochondria and ROS...")
    
    normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
    scaClist = [C0_sca,C10_sca,C14_sca,C42_sca]
    Clist = [C0,C10,C14,C42]
    
    """
        L_stem=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_stem.txt')[0]
        L_CellCycle=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_CellCycle.txt')[0]
        L_Mito=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Mito.txt')[0]
        L_DopaminNeurons=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_DopaminergicNeurons.txt')[0]
        L_ROS=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_ROS.txt')[0]
        
        CumulativeGenesExpressionA_Stemness = np.zeros(250)
        CumulativeGenesExpressionB_Stemness = np.zeros(250)
        for GeneName in L_stem:
            CumulativeGenesExpressionA_Stemness = CumulativeGenesExpressionA_Stemness + C_norm[C.c_names.index(GeneName),0:250]
            CumulativeGenesExpressionB_Stemness = CumulativeGenesExpressionB_Stemness + C_norm[C.c_names.index(GeneName),250:500]
  
        CumulativeGenesExpressionA_CellCycle = np.zeros(250)
        CumulativeGenesExpressionB_CellCycle = np.zeros(250)
        for GeneName in L_CellCycle:
            if GeneName in C.c_names:
                CumulativeGenesExpressionA_CellCycle = CumulativeGenesExpressionA_CellCycle + C_norm[C.c_names.index(GeneName),0:250]
                CumulativeGenesExpressionB_CellCycle = CumulativeGenesExpressionB_CellCycle + C_norm[C.c_names.index(GeneName),250:500]
            
        CumulativeGenesExpressionA_Mito = np.zeros(250)
        CumulativeGenesExpressionB_Mito = np.zeros(250)
        for GeneName in L_Mito:
            if GeneName in C.c_names:
                CumulativeGenesExpressionA_Mito = CumulativeGenesExpressionA_Mito + C_norm[C.c_names.index(GeneName),0:250]
                CumulativeGenesExpressionB_Mito = CumulativeGenesExpressionB_Mito + C_norm[C.c_names.index(GeneName),250:500
                                                                                                ]
        CumulativeGenesExpressionA_DopaminNeurons = np.zeros(250)
        CumulativeGenesExpressionB_DopaminNeurons = np.zeros(250)
        for GeneName in L_DopaminNeurons:
            if GeneName in C.c_names:
                CumulativeGenesExpressionA_DopaminNeurons = CumulativeGenesExpressionA_DopaminNeurons + C_norm[C.c_names.index(GeneName),0:250]
                CumulativeGenesExpressionB_DopaminNeurons = CumulativeGenesExpressionB_DopaminNeurons + C_norm[C.c_names.index(GeneName),250:500]
            
        CumulativeGenesExpressionA_ROS = np.zeros(250)
        CumulativeGenesExpressionB_ROS = np.zeros(250)
        for GeneName in L_stem:
            CumulativeGenesExpressionA_ROS = CumulativeGenesExpressionA_ROS + C_norm[C.c_names.index(GeneName),0:250]
            CumulativeGenesExpressionB_ROS = CumulativeGenesExpressionB_ROS + C_norm[C.c_names.index(GeneName),250:500]
  
        fig1 = plt.figure()
        plt.hist(CumulativeGenesExpressionA_Stemness,bins=15, histtype='stepfilled', color='b', normed=1, label='G2019S', alpha=0.5)        
        plt.hist(CumulativeGenesExpressionB_Stemness,bins=15, histtype='stepfilled', color='r', normed=1, label='GC', alpha=0.5)
        plt.title("Differential expression of Stemness Genes (Cumulative)")
        plt.xlabel("Gene expression (norm)")
        plt.ylabel("Frequency")
        plt.xlim(0, 100)
        plt.ylim(0, 0.12)
        plt.legend()
        plt.show()
        plt.savefig('./Results/jonas_DiffExprStemness_day' + str(DAY) + '.pdf')
        
        fig2 = plt.figure()
        plt.hist(CumulativeGenesExpressionA_CellCycle,bins=15, histtype='stepfilled',  color='b', normed=1, label='G2019S', alpha=0.5)        
        plt.hist(CumulativeGenesExpressionB_CellCycle,bins=15, histtype='stepfilled', color='r', normed=1, label='GC', alpha=0.5)
        plt.title("Differential expression of Cell Cycle Genes (Cumulative)")
        plt.xlabel("Gene expression (norm)")
        plt.ylabel("Frequency")
        plt.xlim(0, 400)
        plt.ylim(0, 0.06)
        plt.legend()
        plt.show()
        plt.savefig('./Results/jonas_DiffExprCellCycle_day' + str(DAY) + '.pdf')

        fig3 = plt.figure()
        plt.hist(CumulativeGenesExpressionA_Mito,bins=15, histtype='stepfilled', color='b', normed=1, label='G2019S', alpha=0.5)        
        plt.hist(CumulativeGenesExpressionB_Mito,bins=15, histtype='stepfilled', color='r', normed=1, label='GC', alpha=0.5)
        plt.title("Differential expression of Mitocondrial Genes (Cumulative)")
        plt.xlabel("Gene expression (norm)")
        plt.ylabel("Frequency")
        plt.xlim(0, 1500)
        plt.ylim(0, 0.006)
        plt.legend()
        plt.show()
        plt.savefig('./Results/jonas_DiffExprMito_day' + str(DAY) + '.pdf')
        
        fig4 = plt.figure()
        plt.hist(CumulativeGenesExpressionA_DopaminNeurons,bins=15, histtype='stepfilled', color='b', normed=1, label='G2019S', alpha=0.5)        
        plt.hist(CumulativeGenesExpressionB_DopaminNeurons,bins=15, histtype='stepfilled', color='r', normed=1, label='GC', alpha=0.5)
        plt.title("Differential expression of Dopaminergic Neurons Genes (Cumulative)")
        plt.xlabel("Gene expression (norm)")
        plt.ylabel("Frequency")
        plt.xlim(0, 35)
        plt.ylim(0, 0.5)
        plt.legend()
        plt.show()
        plt.savefig('./Results/jonas_DiffExprDopaminNeurons_day' + str(DAY) + '.pdf')
        
        fig5 = plt.figure()
        plt.hist(CumulativeGenesExpressionA_ROS,bins=15, histtype='stepfilled', color='b', normed=1, label='G2019S', alpha=0.5)        
        plt.hist(CumulativeGenesExpressionB_ROS,bins=15, histtype='stepfilled', color='r', normed=1, label='GC', alpha=0.5)
        plt.title("Differential expression of ROS Genes (Cumulative)")
        plt.xlabel("Gene expression (norm)")
        plt.ylabel("Frequency")
        plt.xlim(0, 100)
        plt.ylim(0, 0.11)
        plt.legend()
        plt.show()
        plt.savefig('./Results/jonas_DiffExprROS_day' + str(DAY) + '.pdf')
        
        i=i+1
    """
        
    ### Start cumulative plot
    
    Nlines = 5
    Ncols = 4
    fig, axes = plt.subplots(nrows=Nlines, ncols=Ncols, figsize=(20, 15))
    fig.subplots_adjust(hspace=0.3, wspace=0.05)

    for ax in axes.flat:
        # Hide all ticks and labels
        ax.xaxis.set_visible(True)
        ax.yaxis.set_visible(False)
        #ax.tick_params(axis=u'x', which=u'x',length=0)
        #fig.setp(ax.get_xticklabels(),visible=False)

        # Set up ticks only on one side for the "edge" subplots...
        if ax.is_first_col():
            ax.yaxis.set_ticks_position('left')
            ax.yaxis.set_visible(True)
#        if ax.is_last_col():
#            ax.yaxis.set_ticks_position('right')
#            ax.yaxis.set_label_position('right')
#        if ax.is_first_row():
#            ax.xaxis.set_ticks_position('top')
#            ax.xaxis.set_label_position('top')
#        if ax.is_last_row():
#           #ax.xaxis.set_ticks_position('bottom')
#           ax.xaxis.set_visible(True)

    L_stem=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_stem.txt')[0]
    L_CellCycle=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_CellCycle.txt')[0]
    L_Mito=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Mito.txt')[0]
    L_DopaminNeurons=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_DopaminergicNeurons.txt')[0]
    L_ROS=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_ROS.txt')[0]
    
    j = 0
    for DAY in [0,10,14,42]:
        print('Day ' + str(DAY))
        print()
        C = Clist[j]
        C_norm = normClist[j]
        
#        ListOfMeansB = []
#        ListOfStdDevB = []
#        ListOfMeansA = []
#        ListOfStdDevA = []
#        ListOfCoeffOfVarB = []
#        ListOfCoeffOfVarA = []

        MatrixOfCoeffOfVarB = np.ndarray(shape=(Nlines,Ncols))
        MatrixOfCoeffOfVarA = np.ndarray(shape=(Nlines,Ncols))
        
        ListSelectedGenes = []
        for i in range(Nlines):
            CumulativeGenesExpressionA = np.zeros(250)
            CumulativeGenesExpressionB = np.zeros(250)            
            if i == 0:
                ListSelectedGenes = L_stem
                #xmax = 99
                ymax = 7.2#0.11
                if j == 3:
                    plt.legend()
            elif i == 1:
                ListSelectedGenes = L_DopaminNeurons
                #xmax = 34
                ymax = 8.5#0.50
            elif i == 2:
                ListSelectedGenes = L_CellCycle
                #xmax = 395
                ymax = 24#0.058
            elif i == 3:
                ListSelectedGenes = L_Mito
                #xmax = 1495
                ymax = 8#0.0067
            elif i == 4:                
                ListSelectedGenes = L_ROS
                #xmax = 88
                ymax = 4.7#0.12
 
            axes[i,0].set_ylabel("Frequency")
            axes[4,j].set_xlabel("Cumulative Gene Expression (norm. to max)")
            
            for GeneName in ListSelectedGenes:
                if GeneName in C.c_names:
                    CumulativeGenesExpressionA = CumulativeGenesExpressionA + C_norm[C.c_names.index(GeneName),0:250]
                    CumulativeGenesExpressionB = CumulativeGenesExpressionB + C_norm[C.c_names.index(GeneName),250:500]
            
            ######### COMPUTE COEFF OF VARIATION ##########
            
            MeanB = numpy.mean(CumulativeGenesExpressionB)
            StdDevB = numpy.std(CumulativeGenesExpressionB)
            MeanA = numpy.mean(CumulativeGenesExpressionA)
            StdDevA = numpy.std(CumulativeGenesExpressionA)
            
#            ListOfMeansB.append(MeanB)
#            ListOfStdDevB.append(StdDevB)
#            ListOfMeansA.append(MeanA)
#            ListOfStdDevA.append(StdDevA)
#            
#            ListOfCoeffOfVarB.append(StdDevB/MeanB)
#            ListOfCoeffOfVarA.append(StdDevA/MeanA)

            MatrixOfCoeffOfVarB[i,j] = StdDevB/MeanB
            MatrixOfCoeffOfVarA[i,j] = StdDevA/MeanA
                
#            plt.title("Differential expression of Stemness Genes (Cumulative)")
#            plt.xlabel("Gene expression (norm)")
#            plt.ylabel("Frequency")
#            plt.savefig('./Results/jonas_DiffExprStemness_day' + str(DAY) + '.pdf')

            ########## MAKE HISTOGRAMS ##########
            MaxX = max(max(CumulativeGenesExpressionA),max(CumulativeGenesExpressionB))
            CumulativeGenesExpressionA_NORM = CumulativeGenesExpressionA / MaxX
            CumulativeGenesExpressionB_NORM = CumulativeGenesExpressionB / MaxX
            axes[i,j].hist(CumulativeGenesExpressionA_NORM,bins=15, histtype='stepfilled', color=(255.0/256,0.,0.), density=True, label='PAR2-G2019S', alpha=0.5)  #  normed=1,      
            axes[i,j].hist(CumulativeGenesExpressionB_NORM,bins=15, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), density=True, label='PAR2-GC', alpha=0.5) #  normed=1,

            axes[i,j].set_xlim(0,1.1)
            axes[i,j].set_ylim(0,ymax)
            
            plt.legend()

#            axes[0,j].xaxis.set_visible(True)
#            axes[0,j].locator_params(axis='x',nbins=5)
#            plt.setp(axes[0,j].get_xticklabels(), rotation=90, horizontalalignment='left')
#            if i == 0:
#                axes[0,j].set_xlabel(ListOfParametersNamesForLegend[k])
#            axes[0,j].get_xaxis().set_tick_params(direction='out')
#
#            axes[Nlines-1,j].xaxis.set_visible(True)
#            axes[Nlines-1,j].locator_params(axis='x',nbins=5)
#            plt.setp(axes[Nlines-1,j].get_xticklabels(), rotation=90, horizontalalignment='right')
#            if i == 1:
#                axes[Nlines-1,j].set_xlabel(ListOfParametersNamesForLegend[k])
#            axes[Nlines-1,j].get_xaxis().set_tick_params(direction='out')
#
#            axes[i,0].yaxis.set_visible(True)
#            axes[i,0].locator_params(axis='y',nbins=8)
#            axes[i,0].set_ylabel("Root Mean Square (adim.)")
#            axes[i,0].get_yaxis().set_tick_params(direction='out')
#
#            axes[i,Ncols-1].yaxis.set_visible(True)
#            axes[i,Ncols-1].locator_params(axis='y',nbins=8)
#            axes[i,Ncols-1].set_ylabel("Root Mean Square (adim.)")
#            axes[i,Ncols-1].get_yaxis().set_tick_params(direction='out')
#            print("k is ")
#            print(k)

            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            
            import statsmodels.stats.weightstats as wstats
            #DataPopulation1 = numpy.array([1,2,3,4,5,2,3,4,3])
            DataPopulation1 = CumulativeGenesExpressionA
            #DataPopulation2 = numpy.array([1,2,3.001,4,5,2,3,4,3])
            DataPopulation2 = CumulativeGenesExpressionB
            DataPopulation1_Object = wstats.DescrStatsW(DataPopulation1)
            DataPopulation2_Object = wstats.DescrStatsW(DataPopulation2)
            MyComparison = wstats.CompareMeans(DataPopulation1_Object,DataPopulation2_Object)
            TestStatistics, pvalue = MyComparison.ztest_ind(alternative='two-sided', usevar='unequal', value=0)

            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test, 4 days x 7 lists = 28
            # The 7 lists are Stem, DAneur, Mito, Cell Cycle, Pro-apoptosis, Anti-apoptosis, Caspases
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')

        j=j+1

    plt.show()
    plt.savefig('./Results/Fig2B_Fig3C_Fig3E_MitoPaperFig1C_CumulativeExpressionsScores_Stem_DaNeur_CellCyc_Mito_ROS.pdf')

    ### Start Coefficient of Variation plot
    
    Nlines = 5
    Ncols = 1
    fig2, axes2 = plt.subplots(nrows=Nlines, ncols=Ncols)
    fig2.subplots_adjust(hspace=0.3, wspace=0.05)

    for ax in axes2.flat:
        ax.yaxis.set_visible(False)
        if ax.is_first_col():
            ax.yaxis.set_ticks_position('left')
            ax.yaxis.set_visible(True)
 
    for i in range(Nlines):    
        
        if i == 0:
            ymax = 1.1
        elif i == 1:
            ymax = 1.1
        elif i == 2:
            ymax = 2.1
        elif i == 3:
            ymax = 0.21
        elif i == 4:                
            ymax = 0.5
        
        axes2[i].set_ylabel(r"$c_v$")
        axes2[4].set_xlabel("time (days)")                                                 

        axes2[i].plot([0,10,14,42],MatrixOfCoeffOfVarB[i,:], marker='o', color='blue', label='GC')
        axes2[i].plot([0,10,14,42],MatrixOfCoeffOfVarA[i,:], marker='o', color='r', label='G2019S')        
            
        #axes[i].set_xlim(0,1)
        axes2[i].set_ylim(0,ymax)
            
        plt.legend(loc="center right")
    plt.show()

#if RemoveDoubleNamesFromGeneNamesLists == 1:
#
#        L_Aall=sd.get_from_txt('./ListsOfGenesJonas/list_stem.txt')[0]
#        L_Ball=sd.get_from_txt('./ListsOfGenesJonas/list_DopaminergicNeurons.txt')[0]
#        L_Call=sd.get_from_txt('./ListsOfGenesJonas/list_CellCycle.txt')[0]
#        L_Dall=sd.get_from_txt('./ListsOfGenesJonas/list_Mito.txt')[0]
#        L_Eall=sd.get_from_txt('./ListsOfGenesJonas/list_ROS.txt')[0]
#        
#        L_Acleaned = []
#        L_Bcleaned = []
#        L_Ccleaned = []
#        L_Dcleaned = []
#        L_Ecleaned = []
#        
#        L_all = []
#        for GeneName in L_Aall:
#            if not(GeneName in L_Acleaned):
#                L_Acleaned.append(GeneName)
#        for GeneName in L_Ball:
#            if not(GeneName in L_Bcleaned):
#                L_Bcleaned.append(GeneName)
#        for GeneName in L_Call:
#            if not(GeneName in L_Ccleaned):
#                L_Ccleaned.append(GeneName)
#        for GeneName in L_Dall:
#            if not(GeneName in L_Dcleaned):
#                L_Dcleaned.append(GeneName)
#        for GeneName in L_Eall:
#            if not(GeneName in L_Ecleaned):
#                L_Ecleaned.append(GeneName)
#                
#        print(len(L_Aall))
#        print(len(L_Ball))
#        print(len(L_Call))
#        print(len(L_Dall))
#        print(len(L_Eall))
#        print(" ")        
#        print(len(L_Acleaned))
#        print(len(L_Bcleaned))
#        print(len(L_Ccleaned))
#        print(len(L_Dcleaned))
#        print(len(L_Ecleaned))
#        print(" ")

if PlotHistogramsCumulativeCellDeath == 1:
    
    print("Plotting histograms of cumulative gene expression scores for cell death (pro-apoptosis, anti-apoptosis, caspases)...")    
    
    normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
    scaClist = [C0_sca,C10_sca,C14_sca,C42_sca]
    Clist = [C0,C10,C14,C42]
        
    ### Start cumulative plot
    
    Nlines = 3
    Ncols = 4
    fig, axes = plt.subplots(nrows=Nlines, ncols=Ncols, figsize=(20, 15))
    fig.subplots_adjust(hspace=0.3, wspace=0.05)

    for ax in axes.flat:
        # Hide all ticks and labels
        ax.xaxis.set_visible(True)
        ax.yaxis.set_visible(False)
        #ax.tick_params(axis=u'x', which=u'x',length=0)
        #fig.setp(ax.get_xticklabels(),visible=False)

        # Set up ticks only on one side for the "edge" subplots...
        if ax.is_first_col():
            ax.yaxis.set_ticks_position('left')
            ax.yaxis.set_visible(True)
#        if ax.is_last_col():
#            ax.yaxis.set_ticks_position('right')
#            ax.yaxis.set_label_position('right')
#        if ax.is_first_row():
#            ax.xaxis.set_ticks_position('top')
#            ax.xaxis.set_label_position('top')
#        if ax.is_last_row():
#           #ax.xaxis.set_ticks_position('bottom')
#           ax.xaxis.set_visible(True)

    L_ProApoptotic = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Pro-Apoptotic.txt')[0]
    L_AntiApoptotic = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Anti-Apoptotic.txt')[0]
    L_Caspases = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Caspases.txt')[0]

    
    j = 0
    for DAY in [0,10,14,42]:
        print('Day ' + str(DAY))
        print()
        C = Clist[j]
        C_norm = normClist[j]
        
        ListSelectedGenes = []
        for i in range(Nlines):
            CumulativeGenesExpressionA = np.zeros(250)
            CumulativeGenesExpressionB = np.zeros(250)            
            if i == 0:
                ListSelectedGenes = L_ProApoptotic
                ymax = 3.5#0.18
            elif i == 1:
                ListSelectedGenes = L_AntiApoptotic
                ymax = 4#0.2
            elif i == 2:
                ListSelectedGenes = L_Caspases
                ymax = 10#0.7
                if j == 3:
                    plt.legend()

            axes[i,0].set_ylabel("Frequency")
            axes[2,j].set_xlabel("Cumulative Gene Expression (norm. to max)")
            
            for GeneName in ListSelectedGenes:
                #print(GeneName)
                if GeneName in C.c_names:
                    CumulativeGenesExpressionA = CumulativeGenesExpressionA + C_norm[C.c_names.index(GeneName),0:250]
                    CumulativeGenesExpressionB = CumulativeGenesExpressionB + C_norm[C.c_names.index(GeneName),250:500]
            
            ########## MAKE HISTOGRAMS ##########
            MaxX = max(max(CumulativeGenesExpressionA),max(CumulativeGenesExpressionB))
            CumulativeGenesExpressionA_NORM = CumulativeGenesExpressionA / MaxX
            CumulativeGenesExpressionB_NORM = CumulativeGenesExpressionB / MaxX
            axes[i,j].hist(CumulativeGenesExpressionA_NORM,bins=15, histtype='stepfilled', color=(255.0/256,0.,0.),  density=True, label='PAR2-G2019S', alpha=0.5)   # normed = 1,     
            axes[i,j].hist(CumulativeGenesExpressionB_NORM,bins=15, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256),  density=True, label='PAR2-GC', alpha=0.5)    # normed = 1,  

            axes[i,j].set_xlim(0,1.1)
            axes[i,j].set_ylim(0,ymax)

            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            
            #DataPopulation1 = numpy.array([1,2,3,4,5,2,3,4,3])
            DataPopulation1 = CumulativeGenesExpressionA
            #DataPopulation2 = numpy.array([1,2,3.001,4,5,2,3,4,3])
            DataPopulation2 = CumulativeGenesExpressionB
            DataPopulation1_Object = wstats.DescrStatsW(DataPopulation1)
            DataPopulation2_Object = wstats.DescrStatsW(DataPopulation2)
            MyComparison = wstats.CompareMeans(DataPopulation1_Object,DataPopulation2_Object)
            TestStatistics, pvalue = MyComparison.ztest_ind(alternative='two-sided', usevar='unequal', value=0)

            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test, 4 days x 7 lists = 28
            # The 7 lists are Stem, DAneur, Mito, Cell Cycle, Pro-apoptosis, Anti-apoptosis, Caspases
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')

        j=j+1
        
    plt.legend(loc="center right")
    plt.show()
    plt.savefig('./Results/FigS4C_CumulativeExpressionsScores_ProApop_AntiApop_Casp.pdf')
    
    
if PlotHistogramsCumulativeLinnarsonDAneurons == 1:
    
    if LoadLuisData == 1:
    
        a_0_genes, a_0_expr, a_0_expr_np = sd.get_from_txt('./AllData_Pheno_and_GeneExpr_13022018_Luis/A_0_Red_GenExprMat.txt')
        b_0_genes, b_0_expr, b_0_expr_np = sd.get_from_txt('./AllData_Pheno_and_GeneExpr_13022018_Luis/B_0_Red_GenExprMat.txt')
    
        a_10_genes, a_10_expr, a_10_expr_np = sd.get_from_txt('./AllData_Pheno_and_GeneExpr_13022018_Luis/A_10_Red_GenExprMat.txt')
        b_10_genes, b_10_expr, b_10_expr_np = sd.get_from_txt('./AllData_Pheno_and_GeneExpr_13022018_Luis/B_10_Red_GenExprMat.txt')
    
        a_14_genes, a_14_expr, a_14_expr_np = sd.get_from_txt('./AllData_Pheno_and_GeneExpr_13022018_Luis/A_14_Red_GenExprMat.txt')
        b_14_genes, b_14_expr, b_14_expr_np = sd.get_from_txt('./AllData_Pheno_and_GeneExpr_13022018_Luis/B_14_Red_GenExprMat.txt')
    
        a_42_genes, a_42_expr, a_42_expr_np = sd.get_from_txt('./AllData_Pheno_and_GeneExpr_13022018_Luis/A_42_Red_GenExprMat.txt')
        b_42_genes, b_42_expr, b_42_expr_np = sd.get_from_txt('./AllData_Pheno_and_GeneExpr_13022018_Luis/B_42_Red_GenExprMat.txt')
    
    if LuisDataPreliminarySteps == 1:

        a_0_cut = sd.get_top_cells_matrix(a_0_expr_np,500)[0]
        b_0_cut = sd.get_top_cells_matrix(b_0_expr_np,500)[0]
        a_0_norm, a_0_log, a_0_sca = sd.get_normalization(a_0_cut)
        b_0_norm, b_0_log, b_0_sca = sd.get_normalization(b_0_cut)
        print("c1")
     
        a_10_cut = sd.get_top_cells_matrix(a_10_expr_np,250)[0]
        b_10_cut = sd.get_top_cells_matrix(b_10_expr_np,250)[0]
        a_10_norm, a_10_log, a_10_sca = sd.get_normalization(a_10_cut)
        b_10_norm, b_10_log, b_10_sca = sd.get_normalization(b_10_cut)
        print("c2")
     
        a_14_cut = sd.get_top_cells_matrix(a_14_expr_np,505)[0]
        b_14_cut = sd.get_top_cells_matrix(b_14_expr_np,350)[0]
        a_14_norm, a_14_log, a_14_sca = sd.get_normalization(a_14_cut)
        b_14_norm, b_14_log, b_14_sca = sd.get_normalization(b_14_cut)
        print("c3")
    
        a_42_cut = sd.get_top_cells_matrix(a_42_expr_np,400)[0]
        b_42_cut = sd.get_top_cells_matrix(b_42_expr_np,300)[0]
        a_42_norm, a_42_log, a_42_sca = sd.get_normalization(a_42_cut)
        b_42_norm, b_42_log, b_42_sca = sd.get_normalization(b_42_cut)
        print("c4")
    
        c0 = Compare((a_0_cut,a_0_genes),(b_0_cut,b_0_genes)) 
        c10 = Compare((a_10_cut,a_10_genes),(b_10_cut,b_10_genes)) 
        c14 = Compare((a_14_cut,a_14_genes),(b_14_cut,b_14_genes)) 
        c42 = Compare((a_42_cut,a_42_genes),(b_42_cut,b_42_genes)) 
        print("c5")
     
        LUISaux0 = c0.merge()[0]
        LUISaux10 = c10.merge()[0]
        LUISaux14 = c14.merge()[0]
        LUISaux42 = c42.merge()[0]
        print("c6")
        
        c0_norm, c0_log, c0_sca = sd.get_normalization(LUISaux0)
        c10_norm, c10_log, c10_sca = sd.get_normalization(LUISaux10)
        c14_norm, c14_log, c14_sca = sd.get_normalization(LUISaux14)
        c42_norm, c42_log, c42_sca = sd.get_normalization(LUISaux42)
    
        norm_c_list = [c0_norm,c10_norm,c14_norm,c42_norm]
        sca_c_list = [c0_sca,c10_sca,c14_sca,c42_sca]
        c_list = [c0,c10,c14,c42]
        
        ########
        
        ArrayOfDANeurons_A_Day0 = sd.get_from_txt('./AllData_Pheno_and_GeneExpr_13022018_Luis/A_0_Highest_Pheno.txt')[1]
        ArrayOfDANeurons_B_Day0 = sd.get_from_txt('./AllData_Pheno_and_GeneExpr_13022018_Luis/B_0_Highest_Pheno.txt')[1]
        
        indexesOfDANeurons_A_Day0 = [i for i, j in enumerate(ArrayOfDANeurons_A_Day0[0]) if j == 7]
        indexesOfOtherCells_A_Day0 = [i for i, j in enumerate(ArrayOfDANeurons_A_Day0[0]) if j != 7]
        
        indexesOfDANeurons_B_Day0 = [i for i, j in enumerate(ArrayOfDANeurons_B_Day0[0]) if j == 7]
        indexesOfOtherCells_B_Day0 = [i for i, j in enumerate(ArrayOfDANeurons_B_Day0[0]) if j != 7]
        TempArray1 = numpy.array(indexesOfDANeurons_B_Day0)
        TempArray2 = numpy.array(indexesOfOtherCells_B_Day0)
        TranslatedindexesOfDANeurons_B_Day0 = TempArray1 + 500
        TranslatedindexesOfOtherCells_B_Day0 = TempArray2 + 500
        
#        DANeuronsExpression_A_Day0 = c0_norm[:,indexesOfDANeurons_A_Day0]
#        OtherCellsExpression_A_Day0 = c0_norm[:,indexesOfOtherCells_A_Day0]
#        DANeuronsExpression_B_Day0 = c0_norm[:,TranslatedindexesOfDANeurons_B_Day0]
#        OtherCellsExpression_B_Day0 = c0_norm[:,TranslatedindexesOfOtherCells_B_Day0]
        
        ########
        
        ArrayOfDANeurons_A_Day10 = sd.get_from_txt('./AllData_Pheno_and_GeneExpr_13022018_Luis/A_10_Highest_Pheno.txt')[1]
        ArrayOfDANeurons_B_Day10 = sd.get_from_txt('./AllData_Pheno_and_GeneExpr_13022018_Luis/B_10_Highest_Pheno.txt')[1]
        
        indexesOfDANeurons_A_Day10 = [i for i, j in enumerate(ArrayOfDANeurons_A_Day10[0]) if j == 7]
        indexesOfOtherCells_A_Day10 = [i for i, j in enumerate(ArrayOfDANeurons_A_Day10[0]) if j != 7]
        
        indexesOfDANeurons_B_Day10 = [i for i, j in enumerate(ArrayOfDANeurons_B_Day10[0]) if j == 7]
        indexesOfOtherCells_B_Day10 = [i for i, j in enumerate(ArrayOfDANeurons_B_Day10[0]) if j != 7]
        TempArray1 = numpy.array(indexesOfDANeurons_B_Day10)
        TempArray2 = numpy.array(indexesOfOtherCells_B_Day10)
        TranslatedindexesOfDANeurons_B_Day10 = TempArray1 + 250
        TranslatedindexesOfOtherCells_B_Day10 = TempArray2 + 250
        
        ########
        
        ArrayOfDANeurons_A_Day14 = sd.get_from_txt('./AllData_Pheno_and_GeneExpr_13022018_Luis/A_14_Highest_Pheno.txt')[1]
        ArrayOfDANeurons_B_Day14 = sd.get_from_txt('./AllData_Pheno_and_GeneExpr_13022018_Luis/B_14_Highest_Pheno.txt')[1]
        
        indexesOfDANeurons_A_Day14 = [i for i, j in enumerate(ArrayOfDANeurons_A_Day14[0]) if j == 7]
        indexesOfOtherCells_A_Day14 = [i for i, j in enumerate(ArrayOfDANeurons_A_Day14[0]) if j != 7]
        
        indexesOfDANeurons_B_Day14 = [i for i, j in enumerate(ArrayOfDANeurons_B_Day14[0]) if j == 7]
        indexesOfOtherCells_B_Day14 = [i for i, j in enumerate(ArrayOfDANeurons_B_Day14[0]) if j != 7]
        TempArray1 = numpy.array(indexesOfDANeurons_B_Day14)
        TempArray2 = numpy.array(indexesOfOtherCells_B_Day14)
        TranslatedindexesOfDANeurons_B_Day14 = TempArray1 + 505
        TranslatedindexesOfOtherCells_B_Day14 = TempArray2 + 505
        
        ########
        
        ArrayOfDANeurons_A_Day42 = sd.get_from_txt('./AllData_Pheno_and_GeneExpr_13022018_Luis/A_42_Highest_Pheno.txt')[1]
        ArrayOfDANeurons_B_Day42 = sd.get_from_txt('./AllData_Pheno_and_GeneExpr_13022018_Luis/B_42_Highest_Pheno.txt')[1]
        
        indexesOfDANeurons_A_Day42 = [i for i, j in enumerate(ArrayOfDANeurons_A_Day42[0]) if j == 7]
        indexesOfOtherCells_A_Day42 = [i for i, j in enumerate(ArrayOfDANeurons_A_Day42[0]) if j != 7]
        
        indexesOfDANeurons_B_Day42 = [i for i, j in enumerate(ArrayOfDANeurons_B_Day42[0]) if j == 7]
        indexesOfOtherCells_B_Day42 = [i for i, j in enumerate(ArrayOfDANeurons_B_Day42[0]) if j != 7]
        TempArray1 = numpy.array(indexesOfDANeurons_B_Day42)
        TempArray2 = numpy.array(indexesOfOtherCells_B_Day42)
        TranslatedindexesOfDANeurons_B_Day42 = TempArray1 + 400
        TranslatedindexesOfOtherCells_B_Day42 = TempArray2 + 400
        
        ########
        
        DAneurons_Indexes_List_Mutant = [indexesOfDANeurons_A_Day0, indexesOfDANeurons_A_Day10, indexesOfDANeurons_A_Day14, indexesOfDANeurons_A_Day42]
        DAneurons_Indexes_List_Control = [TranslatedindexesOfDANeurons_B_Day0, TranslatedindexesOfDANeurons_B_Day10, TranslatedindexesOfDANeurons_B_Day14, TranslatedindexesOfDANeurons_B_Day42]
        Others_Indexes_List_Mutant = [indexesOfOtherCells_A_Day0, indexesOfOtherCells_A_Day10, indexesOfOtherCells_A_Day14, indexesOfOtherCells_A_Day42]
        Others_Indexes_List_Control = [TranslatedindexesOfOtherCells_B_Day0, TranslatedindexesOfOtherCells_B_Day10, TranslatedindexesOfOtherCells_B_Day14, TranslatedindexesOfOtherCells_B_Day42]
        
    ### Start cumulative plot
    
    Nlines = 4
    Ncols = 4
    fig, axes = plt.subplots(nrows=Nlines, ncols=Ncols)
    fig.subplots_adjust(hspace=0.3, wspace=0.05)

    for ax in axes.flat:
        # Hide all ticks and labels
        ax.xaxis.set_visible(True)
        ax.yaxis.set_visible(False)
        #ax.tick_params(axis=u'x', which=u'x',length=0)
        #fig.setp(ax.get_xticklabels(),visible=False)

        # Set up ticks only on one side for the "edge" subplots...
        if ax.is_first_col():
            ax.yaxis.set_ticks_position('left')
            ax.yaxis.set_visible(True)
#        if ax.is_last_col():
#            ax.yaxis.set_ticks_position('right')
#            ax.yaxis.set_label_position('right')
#        if ax.is_first_row():
#            ax.xaxis.set_ticks_position('top')
#            ax.xaxis.set_label_position('top')
#        if ax.is_last_row():
#           #ax.xaxis.set_ticks_position('bottom')
#           ax.xaxis.set_visible(True)

    L_Linnarson_Long = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLinnarson/list_ph_mDANeurons_Linnarson.txt')[0]
    L_Linnarson_Short = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLinnarson/list_ph_mDANeurons_Linnarson2.txt')[0]
    
    j = 0
    for DAY in [0,10,14,42]:
        print('Day ' + str(DAY))
        print()
        C = c_list[j]
        c_norm = norm_c_list[j]
        
#        ListOfMeansB = []
#        ListOfStdDevB = []
#        ListOfMeansA = []
#        ListOfStdDevA = []
#        ListOfCoeffOfVarB = []
#        ListOfCoeffOfVarA = []

        MatrixOfCoeffOfVarB = np.ndarray(shape=(Nlines,Ncols))
        MatrixOfCoeffOfVarA = np.ndarray(shape=(Nlines,Ncols))

        Red = (255.0/256,0.,0.)
        Gray = (128.0/256,128.0/256,128.0/256)
        RedNEURONS = (128.0/256, 255.0/256, 0.0/256)
        GrayNEURONS = (0/256, 191/256, 255/256)

        Color1 = Red
        Color2 = Gray

        ListSelectedGenes = []
        for i in range(Nlines):          
            if i == 0:
                ListSelectedGenes = L_Linnarson_Long
                #xmax = 99
                ymax = 7.2#0.11
                Color1 = Red
                Color2 = RedNEURONS
                if j == 3:
                    plt.legend()
                IndexesOthers = Others_Indexes_List_Mutant[j]
                IndexesDANeurons = DAneurons_Indexes_List_Mutant[j]
                CumulativeGenesExpressionOthers = np.zeros(len(IndexesOthers))
                CumulativeGenesExpressionDANeurons = np.zeros(len(IndexesDANeurons))  
                labelOthers = 'Other Cells (G2019S - Linnarson Long)'
                labelDANeurons = 'DANeurons (G2019S - Linnarson Long)'
            elif i == 1:
                ListSelectedGenes = L_Linnarson_Long
                #xmax = 34
                ymax = 8.5#0.50
                Color1 = Gray
                Color2 = GrayNEURONS
                IndexesOthers = Others_Indexes_List_Control[j]
                IndexesDANeurons = DAneurons_Indexes_List_Control[j]
                CumulativeGenesExpressionOthers = np.zeros(len(IndexesOthers))
                CumulativeGenesExpressionDANeurons = np.zeros(len(IndexesDANeurons)) 
                labelOthers = 'Other Cells (GC - Linnarson Long)'
                labelDANeurons = 'DANeurons (GC - Linnarson Long)'
            elif i == 2:
                ListSelectedGenes = L_Linnarson_Short
                #xmax = 395
                ymax = 24#0.058
                Color1 = Red
                Color2 = RedNEURONS
                IndexesOthers = Others_Indexes_List_Mutant[j]
                IndexesDANeurons = DAneurons_Indexes_List_Mutant[j]
                CumulativeGenesExpressionOthers = np.zeros(len(IndexesOthers))
                CumulativeGenesExpressionDANeurons = np.zeros(len(IndexesDANeurons))  
                labelOthers = 'Other Cells (G2019S - Linnarson Short)'
                labelDANeurons = 'DANeurons (G2019S - Linnarson Short)'
            elif i == 3:
                ListSelectedGenes = L_Linnarson_Short
                #xmax = 1495
                ymax = 8#0.0067
                Color1 = Gray
                Color2 = GrayNEURONS
                IndexesOthers = Others_Indexes_List_Control[j]
                IndexesDANeurons = DAneurons_Indexes_List_Control[j]
                CumulativeGenesExpressionOthers = np.zeros(len(IndexesOthers))
                CumulativeGenesExpressionDANeurons = np.zeros(len(IndexesDANeurons))  
                labelOthers = 'Other Cells (GC - Linnarson Short)'
                labelDANeurons = 'DANeurons (GC - Linnarson Short)'
#            elif i == 4:                
#                ListSelectedGenes = L_ROS
#                #xmax = 88
#                ymax = 4.7#0.12
 
            axes[i,0].set_ylabel("Number of cells")
            axes[3,j].set_xlabel("Cumulative Gene Expression (norm. to max)")

#
#if GetDataFromFiles == 1:
#    print("Here I am 1")
#    A_0_genes, A_0_expr, A_0_expr_np = sd.get_from_txt('./DataJonas/A_0.txt')
#    B_0_genes, B_0_expr, B_0_expr_np = sd.get_from_txt('./DataJonas/B_0.txt')
#    print("Here I am 2")
#  
#    A_10_genes, A_10_expr, A_10_expr_np = sd.get_from_txt('./DataJonas/A_10.txt')
#    B_10_genes, B_10_expr, B_10_expr_np = sd.get_from_txt('./DataJonas/B_10.txt')
#    print("Here I am 3")
#    
#    A_14_genes, A_14_expr, A_14_expr_np = sd.get_from_txt('./DataJonas/A_14.txt')
#    B_14_genes, B_14_expr, B_14_expr_np = sd.get_from_txt('./DataJonas/B_14.txt')
#    print("Here I am 4")
#    
#    A_42_genes, A_42_expr, A_42_expr_np = sd.get_from_txt('./DataJonas/A_42.txt')
#    B_42_genes, B_42_expr, B_42_expr_np = sd.get_from_txt('./DataJonas/B_42.txt')
#    print("Here I am 4b")

            for GeneName in ListSelectedGenes:
                if GeneName in C.c_names:
                    CumulativeGenesExpressionOthers = CumulativeGenesExpressionOthers + c_norm[C.c_names.index(GeneName),IndexesOthers]
                    CumulativeGenesExpressionDANeurons = CumulativeGenesExpressionDANeurons + c_norm[C.c_names.index(GeneName),IndexesDANeurons]

#            ######### COMPUTE COEFF OF VARIATION ##########
#            
#            MeanB = numpy.mean(CumulativeGenesExpressionB)
#            StdDevB = numpy.std(CumulativeGenesExpressionB)
#            MeanA = numpy.mean(CumulativeGenesExpressionA)
#            StdDevA = numpy.std(CumulativeGenesExpressionA)
#            
##            ListOfMeansB.append(MeanB)
##            ListOfStdDevB.append(StdDevB)
##            ListOfMeansA.append(MeanA)
##            ListOfStdDevA.append(StdDevA)
##            
##            ListOfCoeffOfVarB.append(StdDevB/MeanB)
##            ListOfCoeffOfVarA.append(StdDevA/MeanA)
#
#            MatrixOfCoeffOfVarB[i,j] = StdDevB/MeanB
#            MatrixOfCoeffOfVarA[i,j] = StdDevA/MeanA
#                
##            plt.title("Differential expression of Stemness Genes (Cumulative)")
##            plt.xlabel("Gene expression (norm)")
##            plt.ylabel("Frequency")
##            plt.savefig('./Results/jonas_DiffExprStemness_day' + str(DAY) + '.pdf')

            ########## MAKE HISTOGRAMS ##########
            MaxX = max(max(CumulativeGenesExpressionDANeurons),max(CumulativeGenesExpressionOthers))
            CumulativeGenesExpressionDANeurons_NORM = CumulativeGenesExpressionDANeurons / MaxX
            CumulativeGenesExpressionOthers_NORM = CumulativeGenesExpressionOthers / MaxX
            binwidth = 0.05

            axes[i,j].hist(CumulativeGenesExpressionOthers_NORM, bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=Color1, normed=0, label=labelOthers, alpha=0.5)        
            axes[i,j].hist(CumulativeGenesExpressionDANeurons_NORM, bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=Color2, normed=0, label=labelDANeurons, alpha=0.5)

            axes[i,j].set_xlim(0,1.1)
            axes[0,j].set_ylim(0,100)
            axes[1,j].set_ylim(0,100)
            axes[2,j].set_ylim(0,100)
            axes[3,j].set_ylim(0,110)
            
            #plt.legend()

#            axes[0,j].xaxis.set_visible(True)
#            axes[0,j].locator_params(axis='x',nbins=5)
#            plt.setp(axes[0,j].get_xticklabels(), rotation=90, horizontalalignment='left')
#            if i == 0:
#                axes[0,j].set_xlabel(ListOfParametersNamesForLegend[k])
#            axes[0,j].get_xaxis().set_tick_params(direction='out')
#
#            axes[Nlines-1,j].xaxis.set_visible(True)
#            axes[Nlines-1,j].locator_params(axis='x',nbins=5)
#            plt.setp(axes[Nlines-1,j].get_xticklabels(), rotation=90, horizontalalignment='right')
#            if i == 1:
#                axes[Nlines-1,j].set_xlabel(ListOfParametersNamesForLegend[k])
#            axes[Nlines-1,j].get_xaxis().set_tick_params(direction='out')
#
#            axes[i,0].yaxis.set_visible(True)
#            axes[i,0].locator_params(axis='y',nbins=8)
#            axes[i,0].set_ylabel("Root Mean Square (adim.)")
#            axes[i,0].get_yaxis().set_tick_params(direction='out')
#
#            axes[i,Ncols-1].yaxis.set_visible(True)
#            axes[i,Ncols-1].locator_params(axis='y',nbins=8)
#            axes[i,Ncols-1].set_ylabel("Root Mean Square (adim.)")
#            axes[i,Ncols-1].get_yaxis().set_tick_params(direction='out')
#            print("k is ")
#            print(k)

            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            
            #DataPopulation1 = numpy.array([1,2,3,4,5,2,3,4,3])
            DataPopulation1 = CumulativeGenesExpressionOthers
            #DataPopulation2 = numpy.array([1,2,3.001,4,5,2,3,4,3])
            DataPopulation2 = CumulativeGenesExpressionDANeurons
            DataPopulation1_Object = wstats.DescrStatsW(DataPopulation1)
            DataPopulation2_Object = wstats.DescrStatsW(DataPopulation2)
            MyComparison = wstats.CompareMeans(DataPopulation1_Object,DataPopulation2_Object)
            TestStatistics, pvalue = MyComparison.ztest_ind(alternative='two-sided', usevar='unequal', value=0)

            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test, 4 days x 7 lists = 28
            # The 7 lists are Stem, DAneur, Mito, Cell Cycle, Pro-apoptosis, Anti-apoptosis, Caspases
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')
            
            #plt.legend()

        j=j+1
    
    axes[0,3].legend()
    axes[1,3].legend()    
    axes[2,3].legend()
    axes[3,3].legend()
    
    plt.show()
    
###############################################################################
if ReDoCumulativeWithLinnarsonLists == 1:
    import matplotlib
    
    normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
    scaClist = [C0_sca,C10_sca,C14_sca,C42_sca]
    Clist = [C0,C10,C14,C42]
   
    ### Start cumulative plot
    
    Nlines = 3
    Ncols = 4
    fig, axes = plt.subplots(nrows=Nlines, ncols=Ncols)
    fig.subplots_adjust(hspace=0.3, wspace=0.05)

    for ax in axes.flat:
        # Hide all ticks and labels
        ax.xaxis.set_visible(True)
        ax.yaxis.set_visible(False)
        #ax.tick_params(axis=u'x', which=u'x',length=0)
        #fig.setp(ax.get_xticklabels(),visible=False)

        # Set up ticks only on one side for the "edge" subplots...
        if ax.is_first_col():
            ax.yaxis.set_ticks_position('left')
            ax.yaxis.set_visible(True)
#        if ax.is_last_col():
#            ax.yaxis.set_ticks_position('right')
#            ax.yaxis.set_label_position('right')
#        if ax.is_first_row():
#            ax.xaxis.set_ticks_position('top')
#            ax.xaxis.set_label_position('top')
#        if ax.is_last_row():
#           #ax.xaxis.set_ticks_position('bottom')
#           ax.xaxis.set_visible(True)

    #L_stem=sd.get_from_txt('./ListsOfGenesJonas/list_stem.txt')[0]
    #L_CellCycle=sd.get_from_txt('./ListsOfGenesJonas/list_CellCycle.txt')[0]
    #L_Mito=sd.get_from_txt('./ListsOfGenesJonas/list_Mito.txt')[0]
    L_DopaminNeurons=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_DopaminergicNeurons.txt')[0]
    #L_ROS=sd.get_from_txt('./ListsOfGenesJonas/list_ROS.txt')[0]
    
    L_Linnarson_Long = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLinnarson/list_ph_mDANeurons_Linnarson.txt')[0]
    L_Linnarson_Short = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLinnarson/list_ph_mDANeurons_Linnarson2.txt')[0]
    
    j = 0
    for DAY in [0,10,14,42]:
        print('Day ' + str(DAY))
        print()
        C = Clist[j]
        C_norm = normClist[j]
        

        MatrixOfCoeffOfVarB = np.ndarray(shape=(Nlines,Ncols))
        MatrixOfCoeffOfVarA = np.ndarray(shape=(Nlines,Ncols))
        
        ListSelectedGenes = []
        for i in range(Nlines):
            CumulativeGenesExpressionA = np.zeros(250)
            CumulativeGenesExpressionB = np.zeros(250)            
            if i == 0:
                ListSelectedGenes = L_DopaminNeurons
                #xmax = 99
#                ymax = 7.2#0.11
                if j == 3:
                    plt.legend()
            elif i == 1:
                ListSelectedGenes = L_Linnarson_Long
                #xmax = 34
#                ymax = 8.5#0.50
            elif i == 2:
                ListSelectedGenes = L_Linnarson_Short
                #xmax = 395
#                ymax = 24#0.058
#            elif i == 3:
#                ListSelectedGenes = L_Mito
#                #xmax = 1495
#                ymax = 8#0.0067
#            elif i == 4:                
#                ListSelectedGenes = L_ROS
#                #xmax = 88
#                ymax = 4.7#0.12
 
            axes[i,0].set_ylabel("Frequency")
            axes[2,j].set_xlabel("Cumulative Gene Expression (norm. to max)")
            
            for GeneName in ListSelectedGenes:
                if GeneName in C.c_names:
                    CumulativeGenesExpressionA = CumulativeGenesExpressionA + C_norm[C.c_names.index(GeneName),0:250]
                    CumulativeGenesExpressionB = CumulativeGenesExpressionB + C_norm[C.c_names.index(GeneName),250:500]
            
#            ######### COMPUTE COEFF OF VARIATION ##########
#            
#            MeanB = numpy.mean(CumulativeGenesExpressionB)
#            StdDevB = numpy.std(CumulativeGenesExpressionB)
#            MeanA = numpy.mean(CumulativeGenesExpressionA)
#            StdDevA = numpy.std(CumulativeGenesExpressionA)
#            
##            ListOfMeansB.append(MeanB)
##            ListOfStdDevB.append(StdDevB)
##            ListOfMeansA.append(MeanA)
##            ListOfStdDevA.append(StdDevA)
##            
##            ListOfCoeffOfVarB.append(StdDevB/MeanB)
##            ListOfCoeffOfVarA.append(StdDevA/MeanA)
#
#            MatrixOfCoeffOfVarB[i,j] = StdDevB/MeanB
#            MatrixOfCoeffOfVarA[i,j] = StdDevA/MeanA
#                
##            plt.title("Differential expression of Stemness Genes (Cumulative)")
##            plt.xlabel("Gene expression (norm)")
##            plt.ylabel("Frequency")
##            plt.savefig('./Results/jonas_DiffExprStemness_day' + str(DAY) + '.pdf')

            ########## MAKE HISTOGRAMS ##########
            MaxX = max(max(CumulativeGenesExpressionA),max(CumulativeGenesExpressionB))
            CumulativeGenesExpressionA_NORM = CumulativeGenesExpressionA / MaxX
            CumulativeGenesExpressionB_NORM = CumulativeGenesExpressionB / MaxX
            axes[i,j].hist(CumulativeGenesExpressionA_NORM,bins=15, histtype='stepfilled', color=(255.0/256,0.,0.), normed=1, label='PAR2-G2019S', alpha=0.5)        
            axes[i,j].hist(CumulativeGenesExpressionB_NORM,bins=15, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), normed=1, label='PAR2-GC', alpha=0.5)

            axes[i,j].set_xlim(0,1.1)
            #axes[i,j].set_ylim(0,ymax)
            
            #plt.legend()

            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            
            import statsmodels.stats.weightstats as wstats
            #DataPopulation1 = numpy.array([1,2,3,4,5,2,3,4,3])
            DataPopulation1 = CumulativeGenesExpressionA
            #DataPopulation2 = numpy.array([1,2,3.001,4,5,2,3,4,3])
            DataPopulation2 = CumulativeGenesExpressionB
            DataPopulation1_Object = wstats.DescrStatsW(DataPopulation1)
            DataPopulation2_Object = wstats.DescrStatsW(DataPopulation2)
            MyComparison = wstats.CompareMeans(DataPopulation1_Object,DataPopulation2_Object)
            TestStatistics, pvalue = MyComparison.ztest_ind(alternative='two-sided', usevar='unequal', value=0)

            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test, 4 days x 7 lists = 28
            # The 7 lists are Stem, DAneur, Mito, Cell Cycle, Pro-apoptosis, Anti-apoptosis, Caspases
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')

        j=j+1
    plt.show()
    
    
###############################################################################
######## Histograms Individual genes for Core Regulatory Circuit Lasse ########
###############################################################################

NBonferroni = (16.0 + 51.0) * 4.0

font05 = {'family': 'serif',
    'color':  'darkred',
    'weight': 'normal',
    'size': 10,
    }
                
font01 = {'family': 'serif',
    'color':  'orange',
    'weight': 'normal',
    'size': 10,
    }
                
font001 = {'family': 'serif',
    'color':  'green',
    'weight': 'normal',
    'size': 10,
    }

font0001 = {'family': 'serif',
    'color':  'darkgreen',
    'weight': 'normal',
    'size': 10,
    }

font2 = {'family': 'serif',
    'color':  'blue',
    'weight': 'normal',
    'size': 15,
    }

font3 = {'family': 'serif',
    'color':  'black',
    'weight': 'normal',
    'size': 20,
    }

font3s = {'family': 'serif',
    'color':  'black',
    'weight': 'normal',
    'size': 10,
    }


if PlotHistogramsIndividualGenesCoreRegulatCircuitLasse == 1:
    
    normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
    scaClist = [C0_sca,C10_sca,C14_sca,C42_sca]
    Clist = [C0,C10,C14,C42]

    L_CoreRegulatoryGenesLasse = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLasse/20180427-CRC-CoreTFs-unique-AllCellTypes.txt')[0]
    List_Part_I = L_CoreRegulatoryGenesLasse[0:4]
    List_Part_II = L_CoreRegulatoryGenesLasse[4:8]
    List_Part_III = L_CoreRegulatoryGenesLasse[8:12]
    List_Part_IV = L_CoreRegulatoryGenesLasse[12:16]
    
#    l = [0,0,0,0]
#    TableOfMeansComparingMutCtr = []
#    for i in range(16):
#        TableOfMeansComparingMutCtr.append(l)
    
    ijkijk = 0
    for CurrentList in (List_Part_I, List_Part_II, List_Part_III, List_Part_IV):
            ### Start cumulative plot
    
        Nlines = 4 # 16
        Ncols = 4
        fig, axes = plt.subplots(nrows=Nlines, ncols=Ncols)
        fig.subplots_adjust(hspace=0.3, wspace=0.05)
    
        for ax in axes.flat:
            # Hide all ticks and labels
            ax.xaxis.set_visible(True)
            ax.yaxis.set_visible(False)
            #ax.tick_params(axis=u'x', which=u'x',length=0)
            #fig.setp(ax.get_xticklabels(),visible=False)
    
            # Set up ticks only on one side for the "edge" subplots...
            if ax.is_first_col():
                ax.yaxis.set_ticks_position('left')
                ax.yaxis.set_visible(True)
    #        if ax.is_last_col():
    #            ax.yaxis.set_ticks_position('right')
    #            ax.yaxis.set_label_position('right')
    #        if ax.is_first_row():
    #            ax.xaxis.set_ticks_position('top')
    #            ax.xaxis.set_label_position('top')
    #        if ax.is_last_row():
    #           #ax.xaxis.set_ticks_position('bottom')
    #           ax.xaxis.set_visible(True)
        ijkijk = ijkijk + 1

        j = 0
        ListOfAveragesDays = []
        for DAY in [0,10,14,42]:
            print('Day ' + str(DAY))
            print()

            C = Clist[j]
            C_norm = normClist[j]
            
            ListOfAveragesUpDown = []
            ListSelectedGenes = []
            for i in range(Nlines):
                CumulativeGenesExpressionA = np.zeros(250)
                CumulativeGenesExpressionB = np.zeros(250) 
                ListSelectedGenes = CurrentList[i]
    #
    #            if i == 0:
    #                ymax = 3.5#0.18
    #            elif i == 1:
    #                ListSelectedGenes = L_AntiApoptotic
    #                ymax = 4#0.2
    #            elif i == 2:
    #                ListSelectedGenes = L_Caspases
    #                ymax = 10#0.7

                axes[i,0].set_ylabel("Number of Cells")
                axes[3,j].set_xlabel("Gene Expression (norm. to max)")
                
    #            for GeneName in ListSelectedGenes:
    #                print(GeneName)
    #                if GeneName in C.c_names:
    #                    print('IN')
                GeneName = ListSelectedGenes
    
                CumulativeGenesExpressionA = CumulativeGenesExpressionA + C_norm[C.c_names.index(GeneName),0:250]
                CumulativeGenesExpressionB = CumulativeGenesExpressionB + C_norm[C.c_names.index(GeneName),250:500]
                
                ##### COMPUTE TABLE OF AVERAGES #####
                MeanA = sum(CumulativeGenesExpressionA)/len(CumulativeGenesExpressionA)
                MeanB = sum(CumulativeGenesExpressionB)/len(CumulativeGenesExpressionB)
                if MeanA > MeanB:
                    # MyColorForMeans = MeanA - MeanB
                    MyColorForMeans = 1
                elif MeanA < MeanB:
                    # MyColorForMeans = MeanA - MeanB
                    MyColorForMeans = -1
#                print(MyColorForMeans)
#                kkk = (ijkijk-1) * 4 + i
#                print(kkk)
#                print(j)
#                TableOfMeansComparingMutCtr[kkk][j] = MyColorForMeans
                #####################################
                
                ########## MAKE HISTOGRAMS ##########
                MaxX = max(max(CumulativeGenesExpressionA),max(CumulativeGenesExpressionB))
    #            if MaxX == 0:
    #                MaxX = 1
    #                print('Max Gene Expression for gene ' + str(GeneName) + ' is 0')
                CumulativeGenesExpressionA_NORM = CumulativeGenesExpressionA / MaxX
                CumulativeGenesExpressionB_NORM = CumulativeGenesExpressionB / MaxX
                
                ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
                # of 2 populations are statistically significantly different ######
                
                import statsmodels.stats.weightstats as wstats
                #DataPopulation1 = numpy.array([1,2,3,4,5,2,3,4,3])
                DataPopulation1 = CumulativeGenesExpressionA
                #DataPopulation2 = numpy.array([1,2,3.001,4,5,2,3,4,3])
                DataPopulation2 = CumulativeGenesExpressionB
                DataPopulation1_Object = wstats.DescrStatsW(DataPopulation1)
                DataPopulation2_Object = wstats.DescrStatsW(DataPopulation2)
                MyComparison = wstats.CompareMeans(DataPopulation1_Object,DataPopulation2_Object)
                TestStatistics, pvalue = MyComparison.ztest_ind(alternative='two-sided', usevar='unequal', value=0)
    
                # Remember Bonferroni correction for multiple testing, 
                # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
                # by the number of repetitions of the test, 4 days x 7 lists = 28
                # The 7 lists are Stem, DAneur, Mito, Cell Cycle, Pro-apoptosis, Anti-apoptosis, Caspases
                print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
                print(' ')
                
#                ##########   MAKE HISTOGRAMS   ##########
#                axes[i,j].hist(CumulativeGenesExpressionA_NORM,bins=15, histtype='stepfilled', color=(255.0/256,0.,0.), normed=1, label='PAR2-G2019S', alpha=0.5)        
#                axes[i,j].hist(CumulativeGenesExpressionB_NORM,bins=15, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), normed=1, label='PAR2-GC', alpha=0.5)
#                
#                axes[i,j].hist(CumulativeGenesExpressionA_NORM,bins=15, histtype='stepfilled', color=(255.0/256,0.,0.), normed=0, label='PAR2-G2019S', alpha=0.5)        
#                axes[i,j].hist(CumulativeGenesExpressionB_NORM,bins=15, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), normed=0, label='PAR2-GC', alpha=0.5)
    
                binwidth = 0.05
                axes[i,j].hist(CumulativeGenesExpressionA_NORM,bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(255.0/256,0.,0.), density=True, label='PAR2-G2019S', alpha=0.5)    #  normed=0,    
                axes[i,j].hist(CumulativeGenesExpressionB_NORM,bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), density=True, label='PAR2-GC', alpha=0.5) #  normed=0,
    
    
                if i == 3 and j == 3 :
                    plt.legend()
                                
                BONF = 'For each distribution 250 cells were considered \n Bonferroni correction for multiple hypothesis testing on ' + str(round(NBonferroni)) + ' z-tests (assessing significance of difference of the means) was applied'
                plt.suptitle(BONF)
                    
                if pvalue < 0.001 / NBonferroni:
                    PValueLabel = 'p-val < 0.001'
                    font = font001
                elif pvalue < 0.01 / NBonferroni:
                    PValueLabel = 'p-val < 0.01'
                    font = font01
                elif pvalue < 0.05 / NBonferroni:
                    PValueLabel = 'p-val < 0.05'
                    font = font05
                else:
                    PValueLabel = ' '
                    font = font05
                    MyColorForMeans = 0

                ListOfAveragesUpDown.append(MyColorForMeans)
                
                axes[i,j].text(0.65, 120, PValueLabel, fontdict=font)
                axes[i,3].text(1.1, 0, GeneName, fontdict=font2)

#                if ijkijk == 1:
#                    axes[0,j].set_ylim(0,30)
#                    axes[1,j].set_ylim(0,20)
#                    axes[2,j].set_ylim(0,70)
#                    axes[3,j].set_ylim(0,45)
#                elif ijkijk == 2:
#                    axes[0,j].set_ylim(0,25)
#                    axes[1,j].set_ylim(0,50)
#                    axes[2,j].set_ylim(0,50)
#                    axes[3,j].set_ylim(0,80)                
#                elif ijkijk == 3:
#                    axes[0,j].set_ylim(0,70)
#                    axes[1,j].set_ylim(0,25)
#                    axes[2,j].set_ylim(0,20)
#                    axes[3,j].set_ylim(0,15)                
#                elif ijkijk == 4:
#                    axes[0,j].set_ylim(0,15)
#                    axes[1,j].set_ylim(0,30)
#                    axes[2,j].set_ylim(0,60)
#                    axes[3,j].set_ylim(0,100)  
                    
                if ijkijk == 1:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)
                elif ijkijk == 2:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)                
                elif ijkijk == 3:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)                
                elif ijkijk == 4:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)  
                    
            ListOfAveragesDays.append(ListOfAveragesUpDown)
            j=j+1
        plt.show()
        
        fig, ax = plt.subplots()
        TransposedListOfAverages = np.array(ListOfAveragesDays).T.tolist()
        heatmap = ax.pcolor(TransposedListOfAverages, cmap=plt.cm.seismic, vmin=-1, vmax=1) # RdGy_r
        ax.invert_yaxis()
        ax.text(-0.5, 0.5, CurrentList[0], fontdict=font3)
        ax.text(-0.5, 1.5, CurrentList[1], fontdict=font3)
        ax.text(-0.5, 2.5, CurrentList[2], fontdict=font3)
        ax.text(-0.5, 3.5, CurrentList[3], fontdict=font3)
        ax.text(0.3, -0.3, 'Day 0', fontdict=font3)
        ax.text(1.3, -0.3, 'Day 10', fontdict=font3)
        ax.text(2.3, -0.3, 'Day 14', fontdict=font3)
        ax.text(3.3, -0.3, 'Day 42', fontdict=font3)
        #Spacing between each line
        intervals = 1#float(sys.argv[1])
        loc = plticker.MultipleLocator(base=intervals)
        ax.xaxis.set_major_locator(loc)
        ax.yaxis.set_major_locator(loc)
        ax.grid(which="major", color="black", linestyle='-', linewidth=1)
        ax.tick_params(labelbottom='off')    
        ax.tick_params(labelleft='off')    
        plt.show()
    
#        CountsOfElements = {x:ListOfAveragesDays.count(x) for x in ListOfAveragesDays}
#        
#        CounterControl = 
#        CounterMutant = 
#        CounterNotSignific = 
        

if PlotHistogramsIndividualGenesCoreRegulatCircuitLasseExtendedGenes == 1:
    
    normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
    scaClist = [C0_sca,C10_sca,C14_sca,C42_sca]
    Clist = [C0,C10,C14,C42]

    List_Part_1 = L_CoreRegulatoryGenesLasse[0:4]
    List_Part_2 = L_CoreRegulatoryGenesLasse[4:8]
    List_Part_3 = L_CoreRegulatoryGenesLasse[8:12]
    List_Part_4 = L_CoreRegulatoryGenesLasse[12:16]
    List_Part_5 = L_CoreRegulatoryGenesLasse[16:20]
    List_Part_6 = L_CoreRegulatoryGenesLasse[20:24]
    List_Part_7 = L_CoreRegulatoryGenesLasse[24:28]
    List_Part_8 = L_CoreRegulatoryGenesLasse[28:32]
    List_Part_9 = L_CoreRegulatoryGenesLasse[32:36]
    List_Part_10 = L_CoreRegulatoryGenesLasse[36:40]
    List_Part_11 = L_CoreRegulatoryGenesLasse[40:44]
    List_Part_12 = L_CoreRegulatoryGenesLasse[44:48]
    List_Part_13 = L_CoreRegulatoryGenesLasse[48:52]
#    List_Part_14 = L_CoreRegulatoryGenesLasse[52:56]
#    List_Part_15 = L_CoreRegulatoryGenesLasse[56:60]
#    List_Part_16 = L_CoreRegulatoryGenesLasse[60:64]
#    List_Part_17 = L_CoreRegulatoryGenesLasse[64:68]
    
    ijkijk = 0
    for CurrentList in (List_Part_1, List_Part_2, List_Part_3, List_Part_4, List_Part_5, List_Part_6, List_Part_7, List_Part_8, List_Part_9, List_Part_10, List_Part_11, List_Part_12, List_Part_13):
        ### Start cumulative plot
    
        Nlines = 4 # 16
        if len(CurrentList) == 3:
            Nlines = 3 # 16
        Ncols = 4
        fig, axes = plt.subplots(nrows=Nlines, ncols=Ncols)
        fig.subplots_adjust(hspace=0.3, wspace=0.05)
    
        for ax in axes.flat:
            # Hide all ticks and labels
            ax.xaxis.set_visible(True)
            ax.yaxis.set_visible(False)
            #ax.tick_params(axis=u'x', which=u'x',length=0)
            #fig.setp(ax.get_xticklabels(),visible=False)
    
            # Set up ticks only on one side for the "edge" subplots...
            if ax.is_first_col():
                ax.yaxis.set_ticks_position('left')
                ax.yaxis.set_visible(True)
    #        if ax.is_last_col():
    #            ax.yaxis.set_ticks_position('right')
    #            ax.yaxis.set_label_position('right')
    #        if ax.is_first_row():
    #            ax.xaxis.set_ticks_position('top')
    #            ax.xaxis.set_label_position('top')
    #        if ax.is_last_row():
    #           #ax.xaxis.set_ticks_position('bottom')
    #           ax.xaxis.set_visible(True)
        ijkijk = ijkijk + 1

        j = 0
        ListOfAveragesDays = []
        for DAY in [0,10,14,42]:
            print('Day ' + str(DAY))
            print()

            C = Clist[j]
            C_norm = normClist[j]
            
            ListOfAveragesUpDown = []
            ListSelectedGenes = []
            for i in range(Nlines):
                CumulativeGenesExpressionA = np.zeros(250)
                CumulativeGenesExpressionB = np.zeros(250) 
                ListSelectedGenes = CurrentList[i]
    #
    #            if i == 0:
    #                ymax = 3.5#0.18
    #            elif i == 1:
    #                ListSelectedGenes = L_AntiApoptotic
    #                ymax = 4#0.2
    #            elif i == 2:
    #                ListSelectedGenes = L_Caspases
    #                ymax = 10#0.7

                axes[i,0].set_ylabel("Number of Cells")
                axes[Nlines-1,j].set_xlabel("Gene Expression (norm. to max)")
                
    #            for GeneName in ListSelectedGenes:
    #                print(GeneName)
    #                if GeneName in C.c_names:
    #                    print('IN')
                GeneName = ListSelectedGenes
    
                print(GeneName)
    
                CumulativeGenesExpressionA = CumulativeGenesExpressionA + C_norm[C.c_names.index(GeneName),0:250]
                CumulativeGenesExpressionB = CumulativeGenesExpressionB + C_norm[C.c_names.index(GeneName),250:500]
                
                ##### COMPUTE TABLE OF AVERAGES #####
                MeanA = sum(CumulativeGenesExpressionA)/len(CumulativeGenesExpressionA)
                MeanB = sum(CumulativeGenesExpressionB)/len(CumulativeGenesExpressionB)
                if MeanA > MeanB:
                    # MyColorForMeans = MeanA - MeanB
                    MyColorForMeans = 1
                elif MeanA < MeanB:
                    # MyColorForMeans = MeanA - MeanB
                    MyColorForMeans = -1
#                print(MyColorForMeans)
#                kkk = (ijkijk-1) * 4 + i
#                print(kkk)
#                print(j)
#                TableOfMeansComparingMutCtr[kkk][j] = MyColorForMeans
                #####################################
                
                ########## MAKE HISTOGRAMS ##########
                MaxX = max(max(CumulativeGenesExpressionA),max(CumulativeGenesExpressionB))
    #            if MaxX == 0:
    #                MaxX = 1
    #                print('Max Gene Expression for gene ' + str(GeneName) + ' is 0')
                CumulativeGenesExpressionA_NORM = CumulativeGenesExpressionA / MaxX
                CumulativeGenesExpressionB_NORM = CumulativeGenesExpressionB / MaxX
                
                ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
                # of 2 populations are statistically significantly different ######
                
                import statsmodels.stats.weightstats as wstats
                #DataPopulation1 = numpy.array([1,2,3,4,5,2,3,4,3])
                DataPopulation1 = CumulativeGenesExpressionA
                #DataPopulation2 = numpy.array([1,2,3.001,4,5,2,3,4,3])
                DataPopulation2 = CumulativeGenesExpressionB
                DataPopulation1_Object = wstats.DescrStatsW(DataPopulation1)
                DataPopulation2_Object = wstats.DescrStatsW(DataPopulation2)
                MyComparison = wstats.CompareMeans(DataPopulation1_Object,DataPopulation2_Object)
                TestStatistics, pvalue = MyComparison.ztest_ind(alternative='two-sided', usevar='unequal', value=0)
    
                # Remember Bonferroni correction for multiple testing, 
                # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
                # by the number of repetitions of the test, 4 days x 7 lists = 28
                # The 7 lists are Stem, DAneur, Mito, Cell Cycle, Pro-apoptosis, Anti-apoptosis, Caspases
                print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
                print(' ')
                
                ##########   MAKE HISTOGRAMS   ##########
#                axes[i,j].hist(CumulativeGenesExpressionA_NORM,bins=15, histtype='stepfilled', color=(255.0/256,0.,0.), normed=1, label='PAR2-G2019S', alpha=0.5)        
#                axes[i,j].hist(CumulativeGenesExpressionB_NORM,bins=15, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), normed=1, label='PAR2-GC', alpha=0.5)

                binwidth = 0.05
                axes[i,j].hist(CumulativeGenesExpressionA_NORM,bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(255.0/256,0.,0.),  density=True, label='PAR2-G2019S', alpha=0.5)  #  normed=0,      
                axes[i,j].hist(CumulativeGenesExpressionB_NORM,bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256),  density=True,  label='PAR2-GC', alpha=0.5) # normed=0,
    
                if i == 3 and j == 3 :
                    plt.legend()
    
                BONF = 'For each distribution 250 cells were considered \n Bonferroni correction for multiple hypothesis testing on ' + str(round(NBonferroni)) + ' z-tests (assessing significance of difference of the means) was applied'
                plt.suptitle(BONF)
                    
                if pvalue < 0.001 / NBonferroni:
                    PValueLabel = 'p-val < 0.001'
                    font = font001
                elif pvalue < 0.01 / NBonferroni:
                    PValueLabel = 'p-val < 0.01'
                    font = font01
                elif pvalue < 0.05 / NBonferroni:
                    PValueLabel = 'p-val < 0.05'
                    font = font05
                else:
                    PValueLabel = ' '
                    font = font05
                    MyColorForMeans = 0

                ListOfAveragesUpDown.append(MyColorForMeans)
                axes[i,j].text(0.65, 120, PValueLabel, fontdict=font)
                axes[i,3].text(1.1, 0, GeneName, fontdict=font2)
                
#                if ijkijk == 1:
#                    axes[0,j].set_ylim(0,50)
#                    axes[1,j].set_ylim(0,20)
#                    axes[2,j].set_ylim(0,30)
#                    axes[3,j].set_ylim(0,45)
#                elif ijkijk == 2:
#                    axes[0,j].set_ylim(0,40)
#                    axes[1,j].set_ylim(0,15)
#                    axes[2,j].set_ylim(0,250)
#                    axes[3,j].set_ylim(0,25)                
#                elif ijkijk == 3:
#                    axes[0,j].set_ylim(0,25)
#                    axes[1,j].set_ylim(0,30)
#                    axes[2,j].set_ylim(0,30)
#                    axes[3,j].set_ylim(0,30)                
#                elif ijkijk == 4:
#                    axes[0,j].set_ylim(0,30)
#                    axes[1,j].set_ylim(0,25)
#                    axes[2,j].set_ylim(0,55)
#                    axes[3,j].set_ylim(0,20)    
#                elif ijkijk == 5:
#                    axes[0,j].set_ylim(0,15)
#                    axes[1,j].set_ylim(0,25)
#                    axes[2,j].set_ylim(0,30)
#                    axes[3,j].set_ylim(0,20)                
#                elif ijkijk == 6:
#                    axes[0,j].set_ylim(0,40)
#                    axes[1,j].set_ylim(0,60)
#                    axes[2,j].set_ylim(0,20)
#                    axes[3,j].set_ylim(0,40)                
#                elif ijkijk == 7:
#                    axes[0,j].set_ylim(0,30)
#                    axes[1,j].set_ylim(0,90)
#                    axes[2,j].set_ylim(0,70)
#                    axes[3,j].set_ylim(0,50)   
#                elif ijkijk == 8:
#                    axes[0,j].set_ylim(0,30)
#                    axes[1,j].set_ylim(0,30)
#                    axes[2,j].set_ylim(0,30)
#                    axes[3,j].set_ylim(0,15)                
#                elif ijkijk == 9:
#                    axes[0,j].set_ylim(0,25)
#                    axes[1,j].set_ylim(0,25)
#                    axes[2,j].set_ylim(0,40)
#                    axes[3,j].set_ylim(0,30)                
#                elif ijkijk == 10:
#                    axes[0,j].set_ylim(0,15)
#                    axes[1,j].set_ylim(0,30)
#                    axes[2,j].set_ylim(0,20)
#                    axes[3,j].set_ylim(0,40)   
#                elif ijkijk == 11:
#                    axes[0,j].set_ylim(0,25)
#                    axes[1,j].set_ylim(0,25)
#                    axes[2,j].set_ylim(0,40)
#                    axes[3,j].set_ylim(0,25)                
#                elif ijkijk == 12:
#                    axes[0,j].set_ylim(0,40)
#                    axes[1,j].set_ylim(0,20)
#                    axes[2,j].set_ylim(0,30)
#                    axes[3,j].set_ylim(0,75)                
#                elif ijkijk == 13:
#                    axes[0,j].set_ylim(0,15)
#                    axes[1,j].set_ylim(0,110)
#                    axes[2,j].set_ylim(0,40)
                    
                if ijkijk == 1:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)
                elif ijkijk == 2:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)                
                elif ijkijk == 3:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)                
                elif ijkijk == 4:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)    
                elif ijkijk == 5:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)                
                elif ijkijk == 6:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)                
                elif ijkijk == 7:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)   
                elif ijkijk == 8:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)                
                elif ijkijk == 9:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)                
                elif ijkijk == 10:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)   
                elif ijkijk == 11:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)                
                elif ijkijk == 12:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)                
                elif ijkijk == 13:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    
            ListOfAveragesDays.append(ListOfAveragesUpDown)

            j=j+1
        plt.show()

        fig, ax = plt.subplots()
        TransposedListOfAverages = np.array(ListOfAveragesDays).T.tolist()
        heatmap = ax.pcolor(TransposedListOfAverages, cmap=plt.cm.seismic, vmin=-1, vmax=1) # RdGy_r
        ax.invert_yaxis()
        ax.text(-0.5, 0.5, CurrentList[0], fontdict=font3)
        ax.text(-0.5, 1.5, CurrentList[1], fontdict=font3)
        ax.text(-0.5, 2.5, CurrentList[2], fontdict=font3)
        if ijkijk != 13:
            ax.text(-0.5, 3.5, CurrentList[3], fontdict=font3)
        ax.text(0.3, -0.3, 'Day 0', fontdict=font3)
        ax.text(1.3, -0.3, 'Day 10', fontdict=font3)
        ax.text(2.3, -0.3, 'Day 14', fontdict=font3)
        ax.text(3.3, -0.3, 'Day 42', fontdict=font3)
        #Spacing between each line
        intervals = 1#float(sys.argv[1])
        loc = plticker.MultipleLocator(base=intervals)
        ax.xaxis.set_major_locator(loc)
        ax.yaxis.set_major_locator(loc)
        ax.grid(which="major", color="black", linestyle='-', linewidth=1)
        ax.tick_params(labelbottom='off')    
        ax.tick_params(labelleft='off')    
        plt.show()


if PlotHistogramsIndividualGenesCoreRegulatCircuitLasse_5_TOP_GENES == 1:
    
    print("Plot distributions of expression of 5 individual genes of the CRC, including NR2F1, and their fold changes...")
    
    normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
    scaClist = [C0_sca,C10_sca,C14_sca,C42_sca]
    Clist = [C0,C10,C14,C42]

    L_CoreRegulatoryGenesLasse = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLasse/20180427-CRC-CoreTFs-unique-AllCellTypes_5_TOP_GENES.txt')[0]
    
#    l = [0,0,0,0]
#    TableOfMeansComparingMutCtr = []
#    for i in range(16):
#        TableOfMeansComparingMutCtr.append(l)
    
    ijkijk = 0
    for CurrentList in (L_CoreRegulatoryGenesLasse, L_CoreRegulatoryGenesLasse):
            ### Start cumulative plot
    
        Nlines = 5 # 16
        Ncols = 4
        fig, axes = plt.subplots(nrows=Nlines, ncols=Ncols, figsize=(16, 12))
        fig.subplots_adjust(hspace=0.3, wspace=0.05)
    
        for ax in axes.flat:
            # Hide all ticks and labels
            ax.xaxis.set_visible(True)
            ax.yaxis.set_visible(False)
            #ax.tick_params(axis=u'x', which=u'x',length=0)
            #fig.setp(ax.get_xticklabels(),visible=False)
    
            # Set up ticks only on one side for the "edge" subplots...
            if ax.is_first_col():
                ax.yaxis.set_ticks_position('left')
                ax.yaxis.set_visible(True)
    #        if ax.is_last_col():
    #            ax.yaxis.set_ticks_position('right')
    #            ax.yaxis.set_label_position('right')
    #        if ax.is_first_row():
    #            ax.xaxis.set_ticks_position('top')
    #            ax.xaxis.set_label_position('top')
    #        if ax.is_last_row():
    #           #ax.xaxis.set_ticks_position('bottom')
    #           ax.xaxis.set_visible(True)
        ijkijk = ijkijk + 1

        j = 0
        ListOfAveragesDays = []
        MatrixOfFoldChanges = []
        MatrixOfSEMforFoldChanges = []
        for DAY in [0,10,14,42]:
            print('Day ' + str(DAY))
            print()

            C = Clist[j]
            C_norm = normClist[j]
            
            ListOfAveragesUpDown = []
            ListSelectedGenes = []
            ListOfFoldChangesThisDay = []
            ListOfSTDEVforFoldChangesThisDay = []
            for i in range(Nlines):
                CumulativeGenesExpressionA = np.zeros(250) 
                CumulativeGenesExpressionB = np.zeros(250) 
                ListSelectedGenes = CurrentList[i]
    #
    #            if i == 0:
    #                ymax = 3.5#0.18
    #            elif i == 1:
    #                ListSelectedGenes = L_AntiApoptotic
    #                ymax = 4#0.2
    #            elif i == 2:
    #                ListSelectedGenes = L_Caspases
    #                ymax = 10#0.7

                axes[i,0].set_ylabel("Number of Cells")
                axes[4,j].set_xlabel("Gene Expression (norm. to max)")
                
    #            for GeneName in ListSelectedGenes:
    #                print(GeneName)
    #                if GeneName in C.c_names:
    #                    print('IN')
                GeneName = ListSelectedGenes
    
                print(GeneName)
                CumulativeGenesExpressionA = CumulativeGenesExpressionA + C_norm[C.c_names.index(GeneName),0:250]
                CumulativeGenesExpressionB = CumulativeGenesExpressionB + C_norm[C.c_names.index(GeneName),250:500]
                
                ##### COMPUTE TABLE OF AVERAGES #####
                MeanA = sum(CumulativeGenesExpressionA)/len(CumulativeGenesExpressionA)
                MeanB = sum(CumulativeGenesExpressionB)/len(CumulativeGenesExpressionB)
                
                StdDevA = numpy.std(CumulativeGenesExpressionA)
                StdDevB = numpy.std(CumulativeGenesExpressionB)
                
                StandardErrorOfTheMeanA = StdDevA / math.sqrt(250)
                StandardErrorOfTheMeanB = StdDevB / math.sqrt(250)
                
                if MeanA > MeanB:
                    # MyColorForMeans = MeanA - MeanB
                    MyColorForMeans = 1
                elif MeanA < MeanB:
                    # MyColorForMeans = MeanA - MeanB
                    MyColorForMeans = -1
#                print(MyColorForMeans)
#                kkk = (ijkijk-1) * 4 + i
#                print(kkk)
#                print(j)
#                TableOfMeansComparingMutCtr[kkk][j] = MyColorForMeans
                #####################################
                
                RatioMutOnCtrl = MeanA/MeanB
                MyBase = 2
                FoldChangeCurrentGene = math.log(RatioMutOnCtrl, MyBase) 
                ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
                
                Temp1 = StandardErrorOfTheMeanA ** 2 / (MeanA ** 2 * math.log(2) ** 2)
                Temp2 = StandardErrorOfTheMeanB ** 2 / (MeanB ** 2 * math.log(2) ** 2)
                
                STDDEVofFoldChangeCurrentGene = math.sqrt(Temp1 + Temp2)
                ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)

                ########## MAKE HISTOGRAMS ##########
                MaxX = max(max(CumulativeGenesExpressionA),max(CumulativeGenesExpressionB))
    #            if MaxX == 0:
    #                MaxX = 1
    #                print('Max Gene Expression for gene ' + str(GeneName) + ' is 0')
                CumulativeGenesExpressionA_NORM = CumulativeGenesExpressionA / MaxX
                CumulativeGenesExpressionB_NORM = CumulativeGenesExpressionB / MaxX
                
                ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
                # of 2 populations are statistically significantly different ######
                
                #DataPopulation1 = numpy.array([1,2,3,4,5,2,3,4,3])
                DataPopulation1 = CumulativeGenesExpressionA
                #DataPopulation2 = numpy.array([1,2,3.001,4,5,2,3,4,3])
                DataPopulation2 = CumulativeGenesExpressionB
                DataPopulation1_Object = wstats.DescrStatsW(DataPopulation1)
                DataPopulation2_Object = wstats.DescrStatsW(DataPopulation2)
                MyComparison = wstats.CompareMeans(DataPopulation1_Object,DataPopulation2_Object)
                TestStatistics, pvalue = MyComparison.ztest_ind(alternative='two-sided', usevar='unequal', value=0)
    
                # Remember Bonferroni correction for multiple testing, 
                # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
                # by the number of repetitions of the test, 4 days x 7 lists = 28
                # The 7 lists are Stem, DAneur, Mito, Cell Cycle, Pro-apoptosis, Anti-apoptosis, Caspases
                print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
                print(' ')
                
#                ##########   MAKE HISTOGRAMS   ##########
#                axes[i,j].hist(CumulativeGenesExpressionA_NORM,bins=15, histtype='stepfilled', color=(255.0/256,0.,0.), normed=1, label='PAR2-G2019S', alpha=0.5)        
#                axes[i,j].hist(CumulativeGenesExpressionB_NORM,bins=15, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), normed=1, label='PAR2-GC', alpha=0.5)
#                
#                axes[i,j].hist(CumulativeGenesExpressionA_NORM,bins=15, histtype='stepfilled', color=(255.0/256,0.,0.), normed=0, label='PAR2-G2019S', alpha=0.5)        
#                axes[i,j].hist(CumulativeGenesExpressionB_NORM,bins=15, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), normed=0, label='PAR2-GC', alpha=0.5)
    
                binwidth = 0.05
                axes[i,j].hist(CumulativeGenesExpressionA_NORM,bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(255.0/256,0.,0.), label='PAR2-G2019S', alpha=0.5)    #  normed=0,    
                axes[i,j].hist(CumulativeGenesExpressionB_NORM,bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), label='PAR2-GC', alpha=0.5) # normed=0,
    
    
                if i == 4 and j == 3 :
                    plt.legend()
                                
                BONF = 'For each distribution 250 cells were considered \n Bonferroni correction for multiple hypothesis testing on ' + str(round(NBonferroni)) + ' z-tests (assessing significance of difference of the means) was applied'
                plt.suptitle(BONF)
                    
                if pvalue < 0.001 / NBonferroni:
                    PValueLabel = 'p-val < 0.001'
                    font = font001
                elif pvalue < 0.01 / NBonferroni:
                    PValueLabel = 'p-val < 0.01'
                    font = font01
                elif pvalue < 0.05 / NBonferroni:
                    PValueLabel = 'p-val < 0.05'
                    font = font05
                else:
                    PValueLabel = ' '
                    font = font05
                    MyColorForMeans = 0

                ListOfAveragesUpDown.append(MyColorForMeans)
                
                axes[i,3].text(1.1, 0, GeneName, fontdict=font2)

#                if ijkijk == 1:
#                    axes[0,j].set_ylim(0,30)
#                    axes[1,j].set_ylim(0,20)
#                    axes[2,j].set_ylim(0,70)
#                    axes[3,j].set_ylim(0,45)
#                elif ijkijk == 2:
#                    axes[0,j].set_ylim(0,25)
#                    axes[1,j].set_ylim(0,50)
#                    axes[2,j].set_ylim(0,50)
#                    axes[3,j].set_ylim(0,80)                
#                elif ijkijk == 3:
#                    axes[0,j].set_ylim(0,70)
#                    axes[1,j].set_ylim(0,25)
#                    axes[2,j].set_ylim(0,20)
#                    axes[3,j].set_ylim(0,15)                
#                elif ijkijk == 4:
#                    axes[0,j].set_ylim(0,15)
#                    axes[1,j].set_ylim(0,30)
#                    axes[2,j].set_ylim(0,60)
#                    axes[3,j].set_ylim(0,100)  
                    
                if ijkijk == 1:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)
                    axes[4,j].set_ylim(0,250)
                    axes[i,j].text(0.65, 120, PValueLabel, fontdict=font)

                elif ijkijk == 2:
                    axes[0,j].set_ylim(0,60)
                    axes[1,j].set_ylim(0,35)
                    axes[2,j].set_ylim(0,60)
                    axes[3,j].set_ylim(0,60)
                    axes[4,j].set_ylim(0,50)
                    axes[i,j].text(0.65, 20, PValueLabel, fontdict=font)
                    
            MatrixOfFoldChanges.append(ListOfFoldChangesThisDay)
            MatrixOfSEMforFoldChanges.append(ListOfSTDEVforFoldChangesThisDay)

            ListOfAveragesDays.append(ListOfAveragesUpDown)
            j=j+1
        plt.show()
        if ijkijk == 1:
            MyLabel = ''
        elif ijkijk == 2:
            MyLabel = '_ZOOM'
        plt.savefig('./Results/Fig5C_Expression_CRCGenes_NR2F1' + MyLabel + '.pdf')

        fig2, ax2 = plt.subplots(figsize=(8, 6))
#        ax.xaxis.set_visible(True)
#        ax.yaxis.set_visible(True)
        TransposedListOfAverages = np.array(ListOfAveragesDays).T.tolist()
        heatmap = ax2.pcolor(TransposedListOfAverages, cmap = mpl.colors.ListedColormap([(128.0/256,128.0/256,128.0/256), 'white', (255.0/256,0.,0.)]), vmin=-1, vmax=1, edgecolors='k') # cmap=plt.cm.RdGy_r # cmap=plt.cm.seismic
        ax2.invert_yaxis()
        ax2.text(-0.5, 0.5, CurrentList[0], fontdict=font3s)
        ax2.text(-0.5, 1.5, CurrentList[1], fontdict=font3s)
        ax2.text(-0.5, 2.5, CurrentList[2], fontdict=font3s)
        ax2.text(-0.5, 3.5, CurrentList[3], fontdict=font3s)
        ax2.text(-0.5, 4.5, CurrentList[4], fontdict=font3s)
        ax2.text(0.3, -0.3, 'Day 0', fontdict=font3s)
        ax2.text(1.3, -0.3, 'Day 10', fontdict=font3s)
        ax2.text(2.3, -0.3, 'Day 14', fontdict=font3s)
        ax2.text(3.3, -0.3, 'Day 42', fontdict=font3s)
        #Spacing between each line
        intervals = 1#float(sys.argv[1])
        loc = plticker.MultipleLocator(base=intervals)
        ax2.xaxis.set_major_locator(loc)
        ax2.yaxis.set_major_locator(loc)
        # plt.grid(True)
        # plt.rc('grid', linestyle="-", color='black')
        #ax2.grid(which="major", color="black", linestyle='-', linewidth=1)
        ax2.tick_params(labelbottom='off')    
        ax2.tick_params(labelleft='off')   
        plt.show()
        # plt.grid(False)
    
#        CountsOfElements = {x:ListOfAveragesDays.count(x) for x in ListOfAveragesDays}
#        
#        CounterControl = 
#        CounterMutant = 
#        CounterNotSignific = 

        ########## NOW LET'S PLOT THESE FOLD-CHANGES ##########
        
        N = 4 # Number of days
        
        ind = np.arange(N)  # the x locations for the groups
        width = 0.15       # the width of the bars
        
        fig = plt.figure( figsize=(8, 6))
        ax = fig.add_subplot(111)
        
        Gene0FoldChanges = (MatrixOfFoldChanges[0][0], MatrixOfFoldChanges[1][0], MatrixOfFoldChanges[2][0], MatrixOfFoldChanges[3][0])
        # Gene0Std = (0.1, 0.1, 0.1, 0.1)
        Gene0Std = (MatrixOfSEMforFoldChanges[0][0], MatrixOfSEMforFoldChanges[1][0], MatrixOfSEMforFoldChanges[2][0], MatrixOfSEMforFoldChanges[3][0])
        rects1 = ax.bar(ind, Gene0FoldChanges, width, color='r', yerr=Gene0Std)

        Gene1FoldChanges = (MatrixOfFoldChanges[0][1], MatrixOfFoldChanges[1][1], MatrixOfFoldChanges[2][1], MatrixOfFoldChanges[3][1])
        # Gene1Std = (0.1, 0.1, 0.1, 0.1)
        Gene1Std =  (MatrixOfSEMforFoldChanges[0][1], MatrixOfSEMforFoldChanges[1][1], MatrixOfSEMforFoldChanges[2][1], MatrixOfSEMforFoldChanges[3][1]) #  (0.1, 0.1, 0.1, 0.1)
        rects2 = ax.bar(ind+width, Gene1FoldChanges, width, color='y', yerr=Gene1Std)

        Gene2FoldChanges = (MatrixOfFoldChanges[0][2], MatrixOfFoldChanges[1][2], MatrixOfFoldChanges[2][2], MatrixOfFoldChanges[3][2])
        # Gene2Std = (0.1, 0.1, 0.1, 0.1)
        Gene2Std = (MatrixOfSEMforFoldChanges[0][2], MatrixOfSEMforFoldChanges[1][2], MatrixOfSEMforFoldChanges[2][2], MatrixOfSEMforFoldChanges[3][2])
        rects3 = ax.bar(ind+width+width, Gene2FoldChanges, width, color='b', yerr=Gene2Std)
        
        Gene3FoldChanges = (MatrixOfFoldChanges[0][3], MatrixOfFoldChanges[1][3], MatrixOfFoldChanges[2][3], MatrixOfFoldChanges[3][3])
        # Gene3Std = (0.1, 0.1, 0.1, 0.1)
        Gene3Std = (MatrixOfSEMforFoldChanges[0][3], MatrixOfSEMforFoldChanges[1][3], MatrixOfSEMforFoldChanges[2][3], MatrixOfSEMforFoldChanges[3][3])
        rects4 = ax.bar(ind+width+width+width, Gene3FoldChanges, width, color='c', yerr=Gene3Std)
        
        Gene4FoldChanges = (MatrixOfFoldChanges[0][4], MatrixOfFoldChanges[1][4], MatrixOfFoldChanges[2][4], MatrixOfFoldChanges[3][4])
        # Gene4Std = (0.1, 0.1, 0.1, 0.1)
        Gene4Std = (MatrixOfSEMforFoldChanges[0][4], MatrixOfSEMforFoldChanges[1][4], MatrixOfSEMforFoldChanges[2][4], MatrixOfSEMforFoldChanges[3][4])
        rects5 = ax.bar(ind+width+width+width+width, Gene4FoldChanges, width, color='k', yerr=Gene4Std)
        
        # add some
        ax.set_ylabel(r'Fold change ($log_2(\langle mutant \rangle / \langle wt \rangle)$)')
        # ax.set_title('Fold changes for CRC core genes')
        ax.set_xticks(ind+width)
        ax.set_xticklabels( ('Day 0', 'Day 10', 'Day 14', 'Day 42') )

        ax.legend( (rects1[0], rects2[0], rects3[0], rects4[0], rects5[0]), ('NR2F1', 'NR2F2', 'POU3F2', 'POU3F3', 'SOX2') )

        plt.show()

        def autolabel(rects, values, GeneFoldChanges, GeneStd):
            # attach some text labels
            i = 0
            for rect in rects:
                if GeneFoldChanges[i] > 0:
                    height = GeneFoldChanges[i]
                    Displacement = GeneStd[i] - 0.1
                elif GeneFoldChanges[i] < 0:
                    height = GeneFoldChanges[i]
                    Displacement = - GeneStd[i] - 0.35
                ax.text(rect.get_x()+rect.get_width()/2., height + Displacement, values[i], #'%d'%int(height),
                        ha='center', va='bottom')
                i=i+1
        
        autolabel(rects1,('****','****','****','****'), Gene0FoldChanges, Gene0Std)
        autolabel(rects2,(' ','****','****','****'), Gene1FoldChanges, Gene1Std)
        autolabel(rects3,(' ',' ','****','**'), Gene2FoldChanges, Gene2Std)
        autolabel(rects4,('****','****',' ',' '), Gene3FoldChanges, Gene3Std)
        autolabel(rects5,('**','*',' ','****'), Gene4FoldChanges, Gene4Std)
        print(' ')
        print('ATTENTION!!! SIGNIFICANCES INDICATED WITH ASTERISKS IN THE FOLD CHANGE PLOTS ARE HARD CODED, SO THEY WILL NOT CHANGE IF YOU CHANGE SOME PARAMETERS OF THE ANALYSIS!!!')
        print(' ')

        plt.savefig('./Results/Fig5A_FoldChanges_CRCGenes_NR2F1.pdf')


if PlotHistogramsIndividualGenes_SRR == 1:
    
    normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
    scaClist = [C0_sca,C10_sca,C14_sca,C42_sca]
    Clist = [C0,C10,C14,C42]

    #L_CoreRegulatoryGenesLasse = 
    
#    l = [0,0,0,0]
#    TableOfMeansComparingMutCtr = []
#    for i in range(16):
#        TableOfMeansComparingMutCtr.append(l)
    
    ijkijk = 0
    #for CurrentList in (L_CoreRegulatoryGenesLasse, L_CoreRegulatoryGenesLasse):
            ### Start cumulative plot

    Nlines = 1 # 16
    Ncols = 4
    fig, axes = plt.subplots(nrows=Nlines, ncols=Ncols)
    fig.subplots_adjust(hspace=0.3, wspace=0.05)

    for ax in axes.flat:
        # Hide all ticks and labels
        ax.xaxis.set_visible(True)
        ax.yaxis.set_visible(False)
        #ax.tick_params(axis=u'x', which=u'x',length=0)
        #fig.setp(ax.get_xticklabels(),visible=False)

        # Set up ticks only on one side for the "edge" subplots...
        if ax.is_first_col():
            ax.yaxis.set_ticks_position('left')
            ax.yaxis.set_visible(True)
#        if ax.is_last_col():
#            ax.yaxis.set_ticks_position('right')
#            ax.yaxis.set_label_position('right')
#        if ax.is_first_row():
#            ax.xaxis.set_ticks_position('top')
#            ax.xaxis.set_label_position('top')
#        if ax.is_last_row():
#           #ax.xaxis.set_ticks_position('bottom')
#           ax.xaxis.set_visible(True)
    ijkijk = ijkijk + 1

    j = 0
    ListOfAveragesDays = []
    for DAY in [0,10,14,42]:
        print('Day ' + str(DAY))
        print()

        C = Clist[j]
        C_norm = normClist[j]
        
        ListOfAveragesUpDown = []
        ListSelectedGenes = []
        for i in range(Nlines):
            CumulativeGenesExpressionA = np.zeros(250)
            CumulativeGenesExpressionB = np.zeros(250) 
            #ListSelectedGenes = CurrentList[i]
#
#            if i == 0:
#                ymax = 3.5#0.18
#            elif i == 1:
#                ListSelectedGenes = L_AntiApoptotic
#                ymax = 4#0.2
#            elif i == 2:
#                ListSelectedGenes = L_Caspases
#                ymax = 10#0.7

            axes[0].set_ylabel("Number of Cells")
            axes[j].set_xlabel("Gene Expression (norm. to max)")
            
#            for GeneName in ListSelectedGenes:
#                print(GeneName)
#                if GeneName in C.c_names:
#                    print('IN')
            GeneName = 'SRR'#ListSelectedGenes

            CumulativeGenesExpressionA = CumulativeGenesExpressionA + C_norm[C.c_names.index(GeneName),0:250]
            CumulativeGenesExpressionB = CumulativeGenesExpressionB + C_norm[C.c_names.index(GeneName),250:500]
            
            ##### COMPUTE TABLE OF AVERAGES #####
            MeanA = sum(CumulativeGenesExpressionA)/len(CumulativeGenesExpressionA)
            MeanB = sum(CumulativeGenesExpressionB)/len(CumulativeGenesExpressionB)
            if MeanA > MeanB:
                # MyColorForMeans = MeanA - MeanB
                MyColorForMeans = 1
            elif MeanA < MeanB:
                # MyColorForMeans = MeanA - MeanB
                MyColorForMeans = -1
#                print(MyColorForMeans)
#                kkk = (ijkijk-1) * 4 + i
#                print(kkk)
#                print(j)
#                TableOfMeansComparingMutCtr[kkk][j] = MyColorForMeans
            #####################################
            
            ########## MAKE HISTOGRAMS ##########
            MaxX = max(max(CumulativeGenesExpressionA),max(CumulativeGenesExpressionB))
#            if MaxX == 0:
#                MaxX = 1
#                print('Max Gene Expression for gene ' + str(GeneName) + ' is 0')
            CumulativeGenesExpressionA_NORM = CumulativeGenesExpressionA # / MaxX
            CumulativeGenesExpressionB_NORM = CumulativeGenesExpressionB # / MaxX
            
            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            
            import statsmodels.stats.weightstats as wstats
            #DataPopulation1 = numpy.array([1,2,3,4,5,2,3,4,3])
            DataPopulation1 = CumulativeGenesExpressionA
            #DataPopulation2 = numpy.array([1,2,3.001,4,5,2,3,4,3])
            DataPopulation2 = CumulativeGenesExpressionB
            DataPopulation1_Object = wstats.DescrStatsW(DataPopulation1)
            DataPopulation2_Object = wstats.DescrStatsW(DataPopulation2)
            MyComparison = wstats.CompareMeans(DataPopulation1_Object,DataPopulation2_Object)
            TestStatistics, pvalue = MyComparison.ztest_ind(alternative='two-sided', usevar='unequal', value=0)

            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test, 4 days x 7 lists = 28
            # The 7 lists are Stem, DAneur, Mito, Cell Cycle, Pro-apoptosis, Anti-apoptosis, Caspases
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')
            
#                ##########   MAKE HISTOGRAMS   ##########
#                axes[i,j].hist(CumulativeGenesExpressionA_NORM,bins=15, histtype='stepfilled', color=(255.0/256,0.,0.), normed=1, label='PAR2-G2019S', alpha=0.5)        
#                axes[i,j].hist(CumulativeGenesExpressionB_NORM,bins=15, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), normed=1, label='PAR2-GC', alpha=0.5)
#                
#                axes[i,j].hist(CumulativeGenesExpressionA_NORM,bins=15, histtype='stepfilled', color=(255.0/256,0.,0.), normed=0, label='PAR2-G2019S', alpha=0.5)        
#                axes[i,j].hist(CumulativeGenesExpressionB_NORM,bins=15, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), normed=0, label='PAR2-GC', alpha=0.5)

            binwidth = 0.05
            axes[j].hist(CumulativeGenesExpressionA_NORM,bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(255.0/256,0.,0.), density = True, label='PAR2-G2019S', alpha=0.5)  # normed=0,       
            axes[j].hist(CumulativeGenesExpressionB_NORM,bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), density = True, label='PAR2-GC', alpha=0.5) #  normed=0,


            if i == 4 and j == 3 :
                plt.legend()
                            
            BONF = 'For each distribution 250 cells were considered \n Bonferroni correction for multiple hypothesis testing on ' + str(round(NBonferroni)) + ' z-tests (assessing significance of difference of the means) was applied'
            plt.suptitle(BONF)
                
            if pvalue < 0.001 / NBonferroni:
                PValueLabel = 'p-val < 0.001'
                font = font001
            elif pvalue < 0.01 / NBonferroni:
                PValueLabel = 'p-val < 0.01'
                font = font01
            elif pvalue < 0.05 / NBonferroni:
                PValueLabel = 'p-val < 0.05'
                font = font05
            else:
                PValueLabel = ' '
                font = font05
                MyColorForMeans = 0

            ListOfAveragesUpDown.append(MyColorForMeans)
            
            axes[3].text(1.1, 0, GeneName, fontdict=font2)

#                if ijkijk == 1:
#                    axes[0,j].set_ylim(0,30)
#                    axes[1,j].set_ylim(0,20)
#                    axes[2,j].set_ylim(0,70)
#                    axes[3,j].set_ylim(0,45)
#                elif ijkijk == 2:
#                    axes[0,j].set_ylim(0,25)
#                    axes[1,j].set_ylim(0,50)
#                    axes[2,j].set_ylim(0,50)
#                    axes[3,j].set_ylim(0,80)                
#                elif ijkijk == 3:
#                    axes[0,j].set_ylim(0,70)
#                    axes[1,j].set_ylim(0,25)
#                    axes[2,j].set_ylim(0,20)
#                    axes[3,j].set_ylim(0,15)                
#                elif ijkijk == 4:
#                    axes[0,j].set_ylim(0,15)
#                    axes[1,j].set_ylim(0,30)
#                    axes[2,j].set_ylim(0,60)
#                    axes[3,j].set_ylim(0,100)  
                
            axes[j].set_ylim(0,250)
            axes[j].text(0.65, 120, PValueLabel, fontdict=font)

#            if ijkijk == 1:
#                axes[0,j].set_ylim(0,250)
#                axes[1,j].set_ylim(0,250)
#                axes[2,j].set_ylim(0,250)
#                axes[3,j].set_ylim(0,250)
#                axes[4,j].set_ylim(0,250)
#                axes[i,j].text(0.65, 120, PValueLabel, fontdict=font)
#
#            elif ijkijk == 2:
#                axes[0,j].set_ylim(0,60)
#                axes[1,j].set_ylim(0,35)
#                axes[2,j].set_ylim(0,60)
#                axes[3,j].set_ylim(0,60)
#                axes[4,j].set_ylim(0,50)
#                axes[i,j].text(0.65, 20, PValueLabel, fontdict=font)

        ListOfAveragesDays.append(ListOfAveragesUpDown)
        j=j+1
        plt.show()
        
#        fig, ax = plt.subplots()
#        TransposedListOfAverages = np.array(ListOfAveragesDays).T.tolist()
#        heatmap = ax.pcolor(TransposedListOfAverages, cmap=plt.cm.seismic, vmin=-1, vmax=1) # RdGy_r
#        ax.invert_yaxis()
#        ax.text(-0.5, 0.5, CurrentList[0], fontdict=font3)
##        ax.text(-0.5, 1.5, CurrentList[1], fontdict=font3)
##        ax.text(-0.5, 2.5, CurrentList[2], fontdict=font3)
##        ax.text(-0.5, 3.5, CurrentList[3], fontdict=font3)
##        ax.text(-0.5, 4.5, CurrentList[4], fontdict=font3)
#        ax.text(0.3, -0.3, 'Day 0', fontdict=font3)
#        ax.text(1.3, -0.3, 'Day 10', fontdict=font3)
#        ax.text(2.3, -0.3, 'Day 14', fontdict=font3)
#        ax.text(3.3, -0.3, 'Day 42', fontdict=font3)
#        #Spacing between each line
#        intervals = 1#float(sys.argv[1])
#        loc = plticker.MultipleLocator(base=intervals)
#        ax.xaxis.set_major_locator(loc)
#        ax.yaxis.set_major_locator(loc)
#        ax.grid(which="major", color="black", linestyle='-', linewidth=1)
#        ax.tick_params(labelbottom='off')    
#        ax.tick_params(labelleft='off')    
#        plt.show()
    
#        CountsOfElements = {x:ListOfAveragesDays.count(x) for x in ListOfAveragesDays}
#        
#        CounterControl = 
#        CounterMutant = 
#        CounterNotSignific = 

        
if DoComparison_Lasse == 1:

    #C0_sca(1)
    #C0_norm(1)

    normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
    scaClist = [C0_sca,C10_sca,C14_sca,C42_sca]
    Clist = [C0,C10,C14,C42]
    genesClist = []
                     
    i = 0
    for DAY in [0,10,14,42]:
        print("Preliminary steps for comparison, day " + str(DAY) + "...")
        C_norm = normClist[i]
        C_sca = scaClist[i]
        C = Clist[i]
        i = i+1
        
        #L_CoreRegulatoryGenesLasse = sd.get_from_txt('./ListsOfGenesLasse/20180427-CRC-ExtGenes-unique-AllCellTypes_MODIFIED.txt')[0]
        #L_CoreRegulatoryGenesLasse = sd.get_from_txt('./ListsOfGenesLasse/20180427-CRC-CoreTFs-unique-AllCellTypes.txt')[0]
        L_CoreRegulatoryGenesLasse = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLasse/20180427-CRC-CoreTFs-unique-AllCellTypes_5_TOP_GENES.txt')[0]

        C_genes = C.comm_genes()
        M = sd.get_submatrix(L_CoreRegulatoryGenesLasse, C_sca, C_genes)
        Myt1 = 'COMPARISON,' + str(DAY) + ',sca, ONLY CRC CoreTFs GENES - PCA'
        Myt2 = 'COMPARISON,' + str(DAY) + ',sca, ONLY CRC CoreTFs GENES - tSNE'
        MySaveAs1 = './COMPARISON_day' + str(DAY) + '_PCA_sca_CRC_Lasse_genes.pdf'
        MySaveAs2 = './COMPARISON_day' + str(DAY) + '_tSNE_sca_CRC_Lasse_genes.pdf'
        #a,b,pca_Lasse = get_DR_tris(M[:,0:250],M[:,250:500],com=3,t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)
        #a,b,pca_Lasse = get_DR_tris(M[:,0:250],M[:,250:500],com=2,t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)
        a,b,pca_Lasse = get_DR_tris(M[:,0:250],M[:,250:500],com=2,
                                    t1=Myt1,
                                    t2=Myt2,
                                    SaveAs1=MySaveAs1,
                                    SaveAs2=MySaveAs2, 
                                    Color1=(255.0/256,0.,0.), 
                                    Color2=(128.0/256,128.0/256,128.0/256), 
                                    Label1='G2019S', 
                                    Label2='GC')

        ExplainedVariance = pca_Lasse.explained_variance_ratio_ # Fraction, total over components should be 1.
        print(' ')
        print(ExplainedVariance)
        print(' ')
        Loadings = pca_Lasse.components_ # This extract the loadings, 
        
        #i.e. the Principal axes in feature space, representing the directions 
        # of maximum variance in the data. The components are sorted by explained_variance_.

        #binwidth = 0.05
        #plt.hist(Loadings,bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(255.0/256,0.,0.), normed=0, label='XXX', alpha=0.5)        


if DoCorrelationCRC_Genes_Lasse_ALL_cells == 1:
    
    print("Computing gene-gene correlation over time for the 5 CRC genes including NR2F1...")
    
    WhiteDiagonal = 1 # 1 --> mask values (1) from the diagonal, 0 --> do not mask
    RemoveNotSignificantCorr = 1
    N_Bonf_Corr = 80 # Divide by Bonferroni correction N=80?
    ThrasoldForPval = 0.05/N_Bonf_Corr
    
    ######### ALL CELLS TOGETHER #########
    
    normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
    scaClist = [C0_sca,C10_sca,C14_sca,C42_sca]
    Clist = [C0,C10,C14,C42]
    genesClist = []
    
    #L_CoreRegulatoryGenesLasse = sd.get_from_txt('./ListsOfGenesLasse/20180427-CRC-ExtGenes-unique-AllCellTypes_MODIFIED.txt')[0]
    #L_CoreRegulatoryGenesLasse = sd.get_from_txt('./ListsOfGenesLasse/20180427-CRC-CoreTFs-unique-AllCellTypes.txt')[0]
    L_CoreRegulatoryGenesLasse = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLasse/20180427-CRC-CoreTFs-unique-AllCellTypes_5_TOP_GENES.txt')[0]
    
    MyAverageCorrVector = []
    for ControlMutant in [0, 250]:
        k = 0
        for DAY in [0,10,14,42]:
            C_norm = normClist[k]
            C_sca = scaClist[k]
            C = Clist[k]
            k = k+1
            
            CurrentList = []
            N = len(L_CoreRegulatoryGenesLasse)
            pvalPearsonCorr_Matrix = 100 * np.ones([N, N])
            MyCorrelationMatrix_Full = 100 * np.ones([N, N])
            MyCorrelationMatrix_Masked = 100 * np.ones([N, N])

            i = 0
            for GeneName1 in L_CoreRegulatoryGenesLasse:
                j = 0
                CurrentList.append(GeneName1)
                for GeneName2 in L_CoreRegulatoryGenesLasse:
                    
                    Gene_1_Expression = C_norm[C.c_names.index(GeneName1),0+ControlMutant:250+ControlMutant]
        #            Gene_1_ExpressionB = C_norm[C.c_names.index(GeneName1),250:500]
                    Gene_2_Expression = C_norm[C.c_names.index(GeneName2),0+ControlMutant:250+ControlMutant]
        #            Gene_2_ExpressionB = C_norm[C.c_names.index(GeneName2),250:500]
                    
                    PearsonCorrCoeff, pvalPearsonCorr = pearsonr(Gene_1_Expression, Gene_2_Expression)
                    pvalPearsonCorr_Matrix[i,j] = pvalPearsonCorr
                    
                    if RemoveNotSignificantCorr == 0:
                        MyCorrelationMatrix_Full[i,j] = PearsonCorrCoeff
                                      
                        if WhiteDiagonal == 1:
                            MyCorrelationMatrix_Masked[i,j] = PearsonCorrCoeff
                            if i == j:
                                MyCorrelationMatrix_Masked[i,j] = 0
                    elif RemoveNotSignificantCorr == 1:
                        if pvalPearsonCorr <= ThrasoldForPval:
                            MyCorrelationMatrix_Full[i,j] = PearsonCorrCoeff
                        elif pvalPearsonCorr > ThrasoldForPval:
                            MyCorrelationMatrix_Full[i,j] = 0
                                      
                        if WhiteDiagonal == 1:
                            MyCorrelationMatrix_Masked[i,j] = MyCorrelationMatrix_Full[i,j]
                            if i == j:
                                MyCorrelationMatrix_Masked[i,j] = 0
                                                      
                                                      
                                                      
#                    print(GeneName1)
#                    print(GeneName2)
#                    print(PearsonCorrCoeff)
#                    print(' ')
                    
                    j = j + 1
                i = i + 1
        
            fig, ax = plt.subplots(figsize=(16,12))
            #TransposedListOfAverages = np.array(ListOfAveragesDays).T.tolist()
            #sns.heatmap(MyCorrelationMatrix, annot=False,  linewidths=.5, cmap=plt.cm.seismic, vmin=-1, vmax=1, ax = ax)
            #ax.set_xticks() sns.heatmap(data, annot=True, linewidths=.5)
            
            if WhiteDiagonal == 0:
                MyCorrelationMatrix = MyCorrelationMatrix_Full
            elif WhiteDiagonal == 1:
                MyCorrelationMatrix = MyCorrelationMatrix_Masked
            else:
                print("Error in Masking Correlations!")
            
            MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix)
            for i in [0,1,2,3,4]:
                for j in [0,1,2,3,4]:
                    if j>=i:
                        MyCorrelationMatrix_Triangle[i,j] = 0
            
            heatmap = ax.pcolor(MyCorrelationMatrix_Triangle, cmap=plt.cm.seismic_r, vmin=-1, vmax=1) # RdGy_r  # PiYG  # seismic   # PRGn_r
            ax.invert_yaxis()
            ax.text(-0.7, 0.5, CurrentList[0], fontdict=font3)
            ax.text(-0.7, 1.5, CurrentList[1], fontdict=font3)
            ax.text(-0.7, 2.5, CurrentList[2], fontdict=font3)
            ax.text(-0.7, 3.5, CurrentList[3], fontdict=font3)
            ax.text(-0.7, 4.5, CurrentList[4], fontdict=font3)
            ax.text(0.3, -0.3, CurrentList[0], fontdict=font3)
            ax.text(1.3, -0.3, CurrentList[1], fontdict=font3)
            ax.text(2.3, -0.3, CurrentList[2], fontdict=font3)
            ax.text(3.3, -0.3, CurrentList[3], fontdict=font3)
            ax.text(4.3, -0.3, CurrentList[4], fontdict=font3)
            #Spacing between each line
            intervals = 1#float(sys.argv[1])
            loc = plticker.MultipleLocator(base=intervals)
            ax.xaxis.set_major_locator(loc)
            ax.yaxis.set_major_locator(loc)
            ax.grid(which="major", color="black", linestyle='-', linewidth=1)
            ax.tick_params(labelbottom='off')    
            ax.tick_params(labelleft='off')    
            fig.colorbar(heatmap)
            i = 0
            j = 0
            for i in [0,1,2,3,4]:
                for j in [0,1,2,3,4]:
                    if j>=i:
                        if i != j:
                            if MyCorrelationMatrix[i,j] != 0:
                                Current_pvalPearsonCorr = pvalPearsonCorr_Matrix[i,j]
                                if Current_pvalPearsonCorr <= 0.001/N_Bonf_Corr:
                                    Significance = ' ***'
                                elif Current_pvalPearsonCorr <= 0.01/N_Bonf_Corr:
                                    Significance = ' **'
                                elif Current_pvalPearsonCorr <= 0.05/N_Bonf_Corr:
                                    Significance = ' *'
                                MyText = str(round(MyCorrelationMatrix[i,j],2)) + Significance
                                ax.text(i+0.5, j+0.5, MyText, 
                                    horizontalalignment='center', verticalalignment='center', fontdict=font3)
                            elif MyCorrelationMatrix[i,j] == 0:
                                ax.text(i+0.5, j+0.5, r'$\approx 0$', 
                                    horizontalalignment='center', verticalalignment='center', fontdict=font3)
                        elif i == j:
                            if WhiteDiagonal == 0:
                                ax.text(i+0.5, j+0.5, round(MyCorrelationMatrix[i,j],2), 
                                        horizontalalignment='center', verticalalignment='center', fontdict=font3)
                            elif WhiteDiagonal == 1:
                                ax.text(i+0.5, j+0.5,' ', 
                                        horizontalalignment='center', verticalalignment='center', fontdict=font3)
            plt.show()
            
            Dim = len(MyCorrelationMatrix)
            if WhiteDiagonal == 0:
                AverageCorr = (sum(sum(abs(MyCorrelationMatrix)))-Dim)/(Dim*Dim-Dim)
            elif WhiteDiagonal == 1: # Note that MyCorrelationMatrix contains the diagonal here above, doesn't here below.
                AverageCorr = (sum(sum(abs(MyCorrelationMatrix))))/(Dim*Dim-Dim)
            print(AverageCorr)
            print(' ')
            MyAverageCorrVector.append(AverageCorr)            
            
            
    plt.figure(figsize=(10,4))      
    plt.plot([0,10,14,42], MyAverageCorrVector[0:4], marker='o', linestyle='--', color=(255.0/256,0.,0.), label='G2019S')
    plt.plot([0,10,14,42], MyAverageCorrVector[4:8], marker='s', linestyle='-', color=(128.0/256,128.0/256,128.0/256), label='GC')
    plt.xlabel('Time (Days after differentiation)')
    plt.ylabel(r'Average Gene-Gene Correlation')
    plt.ylim(0,0.12)#0.2)
    plt.title('Average of the absolute value of the Gene-Gene Pearson Correlation Coefficient')
    plt.legend(loc='upper left')
    plt.show()
    plt.savefig('./Results/Fig5B_AverageGeneGeneCorrelation_CRC.pdf')

            
if DoCorrelationCRC_Genes_Lasse_LowHigh_NR2F1_cells == 1:
    
    ######### SEPARATE CELLS WITH LOW / HIGH NR2F1 #########
    
    normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
    scaClist = [C0_sca,C10_sca,C14_sca,C42_sca]
    Clist = [C0,C10,C14,C42]
    genesClist = []
    
    #L_CoreRegulatoryGenesLasse = sd.get_from_txt('./ListsOfGenesLasse/20180427-CRC-ExtGenes-unique-AllCellTypes_MODIFIED.txt')[0]
    #L_CoreRegulatoryGenesLasse = sd.get_from_txt('./ListsOfGenesLasse/20180427-CRC-CoreTFs-unique-AllCellTypes.txt')[0]
    L_CoreRegulatoryGenesLasse = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLasse/20180427-CRC-CoreTFs-unique-AllCellTypes_5_TOP_GENES.txt')[0]
    
    MyAverageCorrVector = []
    for ControlMutant in [0, 250]:
        k = 0
        for DAY in [0,10,14,42]:
            C_norm = normClist[k]
            C_sca = scaClist[k]
            C = Clist[k]
            k = k+1
            
            CurrentList = []
            N = len(L_CoreRegulatoryGenesLasse)
            MyCorrelationMatrix = 100 * np.ones([N, N])
            i = 0
            for GeneName1 in L_CoreRegulatoryGenesLasse:
                j = 0
                CurrentList.append(GeneName1)
                for GeneName2 in L_CoreRegulatoryGenesLasse:
                    
                    Gene_1_Expression = C_norm[C.c_names.index(GeneName1),0+ControlMutant:250+ControlMutant]
        #            Gene_1_ExpressionB = C_norm[C.c_names.index(GeneName1),250:500]
                    Gene_2_Expression = C_norm[C.c_names.index(GeneName2),0+ControlMutant:250+ControlMutant]
        #            Gene_2_ExpressionB = C_norm[C.c_names.index(GeneName2),250:500]
                    
                    PearsonCorrCoeff, pvalPearsonCorr = pearsonr(Gene_1_Expression, Gene_2_Expression)
                    MyCorrelationMatrix[i,j] = PearsonCorrCoeff
        
#                    print(GeneName1)
#                    print(GeneName2)
#                    print(PearsonCorrCoeff)
#                    print(' ')
                    
                    j = j + 1
                i = i + 1
        
            fig, ax = plt.subplots()
            #TransposedListOfAverages = np.array(ListOfAveragesDays).T.tolist()
            #sns.heatmap(MyCorrelationMatrix, annot=False,  linewidths=.5, cmap=plt.cm.seismic, vmin=-1, vmax=1, ax = ax)
            #ax.set_xticks() sns.heatmap(data, annot=True, linewidths=.5)
            
            heatmap = ax.pcolor(MyCorrelationMatrix, cmap=plt.cm.seismic, vmin=-1, vmax=1) # RdGy_r
            ax.invert_yaxis()
            ax.text(-0.7, 0.5, CurrentList[0], fontdict=font3)
            ax.text(-0.7, 1.5, CurrentList[1], fontdict=font3)
            ax.text(-0.7, 2.5, CurrentList[2], fontdict=font3)
            ax.text(-0.7, 3.5, CurrentList[3], fontdict=font3)
            ax.text(-0.7, 4.5, CurrentList[4], fontdict=font3)
            ax.text(0.3, -0.3, CurrentList[0], fontdict=font3)
            ax.text(1.3, -0.3, CurrentList[1], fontdict=font3)
            ax.text(2.3, -0.3, CurrentList[2], fontdict=font3)
            ax.text(3.3, -0.3, CurrentList[3], fontdict=font3)
            ax.text(4.3, -0.3, CurrentList[4], fontdict=font3)
            #Spacing between each line
            intervals = 1#float(sys.argv[1])
            loc = plticker.MultipleLocator(base=intervals)
            ax.xaxis.set_major_locator(loc)
            ax.yaxis.set_major_locator(loc)
            ax.grid(which="major", color="black", linestyle='-', linewidth=1)
            ax.tick_params(labelbottom='off')    
            ax.tick_params(labelleft='off')    
            fig.colorbar(heatmap)
            i = 0
            j = 0
            for i in [0,1,2,3,4]:
                for j in [0,1,2,3,4]:
                    ax.text(i+0.5, j+0.5, round(MyCorrelationMatrix[i,j],2), 
                            horizontalalignment='center', verticalalignment='center', fontdict=font3)
            plt.show()
            
            Dim = len(MyCorrelationMatrix)
            AverageCorr = (sum(sum(abs(MyCorrelationMatrix)))-Dim)/(Dim*Dim-Dim)
            print(AverageCorr)
            print(' ')
            MyAverageCorrVector.append(AverageCorr)
            
    plt.figure()      
    plt.plot([0,10,14,42], MyAverageCorrVector[0:4], marker='o', linestyle='--', color='r', label='G2019S')
    plt.plot([0,10,14,42], MyAverageCorrVector[4:8], marker='s', linestyle='-', color='gray', label='GC')
    plt.xlabel('Time (Days after differentiation)')
    plt.ylabel(r'Average Gene-Gene Correlation')
    plt.ylim(0,0.2)
    plt.title('Average of the absolute value of the Gene-Gene Pearson Correlation Coefficient')
    plt.legend()
    plt.show()            



if Plots_ZIC_12345_Lasse_Jocken == 1: 
    
    normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
    scaClist = [C0_sca,C10_sca,C14_sca,C42_sca]
    Clist = [C0,C10,C14,C42]

    L_CoreRegulatoryGenesLasse = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJochen/list_ZIC_GENES.txt')[0]
    
#    l = [0,0,0,0]
#    TableOfMeansComparingMutCtr = []
#    for i in range(16):
#        TableOfMeansComparingMutCtr.append(l)
    
    ijkijk = 0
    for CurrentList in (L_CoreRegulatoryGenesLasse, L_CoreRegulatoryGenesLasse):
            ### Start cumulative plot
    
        Nlines = 5 # 16
        Ncols = 4
        fig, axes = plt.subplots(nrows=Nlines, ncols=Ncols)
        fig.subplots_adjust(hspace=0.3, wspace=0.05)
    
        for ax in axes.flat:
            # Hide all ticks and labels
            ax.xaxis.set_visible(True)
            ax.yaxis.set_visible(False)
            #ax.tick_params(axis=u'x', which=u'x',length=0)
            #fig.setp(ax.get_xticklabels(),visible=False)
    
            # Set up ticks only on one side for the "edge" subplots...
            if ax.is_first_col():
                ax.yaxis.set_ticks_position('left')
                ax.yaxis.set_visible(True)
    #        if ax.is_last_col():
    #            ax.yaxis.set_ticks_position('right')
    #            ax.yaxis.set_label_position('right')
    #        if ax.is_first_row():
    #            ax.xaxis.set_ticks_position('top')
    #            ax.xaxis.set_label_position('top')
    #        if ax.is_last_row():
    #           #ax.xaxis.set_ticks_position('bottom')
    #           ax.xaxis.set_visible(True)
        ijkijk = ijkijk + 1

        j = 0
        ListOfAveragesDays = []
        MatrixOfFoldChanges = []
        MatrixOfSEMforFoldChanges = []
        for DAY in [0,10,14,42]:
            print('Day ' + str(DAY))
            print()

            C = Clist[j]
            C_norm = normClist[j]
            
            ListOfAveragesUpDown = []
            ListSelectedGenes = []
            ListOfFoldChangesThisDay = []
            ListOfSTDEVforFoldChangesThisDay = []
            for i in range(Nlines):
                CumulativeGenesExpressionA = np.zeros(250) 
                CumulativeGenesExpressionB = np.zeros(250) 
                ListSelectedGenes = CurrentList[i]
    #
    #            if i == 0:
    #                ymax = 3.5#0.18
    #            elif i == 1:
    #                ListSelectedGenes = L_AntiApoptotic
    #                ymax = 4#0.2
    #            elif i == 2:
    #                ListSelectedGenes = L_Caspases
    #                ymax = 10#0.7

                axes[i,0].set_ylabel("Number of Cells")
                axes[4,j].set_xlabel("Gene Expression (norm. to max)")
                
    #            for GeneName in ListSelectedGenes:
    #                print(GeneName)
    #                if GeneName in C.c_names:
    #                    print('IN')
                GeneName = ListSelectedGenes
                if not (DAY == 14 and GeneName == 'ZIC3'):
                    CumulativeGenesExpressionA = CumulativeGenesExpressionA + C_norm[C.c_names.index(GeneName),0:250]
                    CumulativeGenesExpressionB = CumulativeGenesExpressionB + C_norm[C.c_names.index(GeneName),250:500]
                
                ##### COMPUTE TABLE OF AVERAGES #####
                if not (DAY == 14 and GeneName == 'ZIC3'):
                    MeanA = sum(CumulativeGenesExpressionA)/len(CumulativeGenesExpressionA)
                    MeanB = sum(CumulativeGenesExpressionB)/len(CumulativeGenesExpressionB)
                    
                if not (DAY == 14 and GeneName == 'ZIC3'):
                    StdDevA = numpy.std(CumulativeGenesExpressionA)
                    StdDevB = numpy.std(CumulativeGenesExpressionB)
                
                
                if not (DAY == 14 and GeneName == 'ZIC3'):
                    StandardErrorOfTheMeanA = StdDevA / math.sqrt(250)
                    StandardErrorOfTheMeanB = StdDevB / math.sqrt(250)
                
#                if MeanA > MeanB:
#                    # MyColorForMeans = MeanA - MeanB
#                    MyColorForMeans = 1
#                elif MeanA < MeanB:
#                    # MyColorForMeans = MeanA - MeanB
#                    MyColorForMeans = -1
##                print(MyColorForMeans)
##                kkk = (ijkijk-1) * 4 + i
##                print(kkk)
##                print(j)
##                TableOfMeansComparingMutCtr[kkk][j] = MyColorForMeans
                #####################################
                
                RatioMutOnCtrl = MeanA/MeanB
                MyBase = 2
                print()
                print(MeanA)
                print(MeanB)
                print(RatioMutOnCtrl)
                print()
                
                if not (DAY == 14 and GeneName == 'ZIC3'):
                    Temp1 = StandardErrorOfTheMeanA ** 2 / (MeanA ** 2 * math.log(2) ** 2)
                    Temp2 = StandardErrorOfTheMeanB ** 2 / (MeanB ** 2 * math.log(2) ** 2)
                
                    STDDEVofFoldChangeCurrentGene = math.sqrt(Temp1 + Temp2)
                else:
                    STDDEVofFoldChangeCurrentGene = 0
                ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)

                ########## MAKE HISTOGRAMS ##########
                if not (DAY == 14 and GeneName == 'ZIC3'):
                    MaxX = max(max(CumulativeGenesExpressionA),max(CumulativeGenesExpressionB))
                else:
                    MaxX = max(CumulativeGenesExpressionB)

    #            if MaxX == 0:
    #                MaxX = 1
    #                print('Max Gene Expression for gene ' + str(GeneName) + ' is 0')
                if not (DAY == 14 and GeneName == 'ZIC3'):
                     CumulativeGenesExpressionA_NORM = CumulativeGenesExpressionA # / MaxX
                     CumulativeGenesExpressionB_NORM = CumulativeGenesExpressionB # / MaxX
                
                ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
                # of 2 populations are statistically significantly different ######
                
                if not (DAY == 14 and GeneName == 'ZIC3'):
                   
                    import statsmodels.stats.weightstats as wstats
                    #DataPopulation1 = numpy.array([1,2,3,4,5,2,3,4,3])
                    DataPopulation1 = CumulativeGenesExpressionA
                    #DataPopulation2 = numpy.array([1,2,3.001,4,5,2,3,4,3])
                    DataPopulation2 = CumulativeGenesExpressionB
                    DataPopulation1_Object = wstats.DescrStatsW(DataPopulation1)
                    DataPopulation2_Object = wstats.DescrStatsW(DataPopulation2)
                    MyComparison = wstats.CompareMeans(DataPopulation1_Object,DataPopulation2_Object)
                    TestStatistics, pvalue = MyComparison.ztest_ind(alternative='two-sided', usevar='unequal', value=0)
        
                    # Remember Bonferroni correction for multiple testing, 
                    # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
                    # by the number of repetitions of the test, 4 days x 7 lists = 28
                    # The 7 lists are Stem, DAneur, Mito, Cell Cycle, Pro-apoptosis, Anti-apoptosis, Caspases
                else:
                    TestStatistics = 999 
                    pvalue =999
                                        
                print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
                print(' ')
                    
                if MeanA == 0.0:
                    FoldChangeCurrentGene = 0 
                elif DAY == 14 and GeneName == 'ZIC3':
                    FoldChangeCurrentGene = 0 
                else:
                    FoldChangeCurrentGene = math.log(RatioMutOnCtrl, MyBase) 
                ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
                
#                ##########   MAKE HISTOGRAMS   ##########
#                axes[i,j].hist(CumulativeGenesExpressionA_NORM,bins=15, histtype='stepfilled', color=(255.0/256,0.,0.), normed=1, label='PAR2-G2019S', alpha=0.5)        
#                axes[i,j].hist(CumulativeGenesExpressionB_NORM,bins=15, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), normed=1, label='PAR2-GC', alpha=0.5)
#                
#                axes[i,j].hist(CumulativeGenesExpressionA_NORM,bins=15, histtype='stepfilled', color=(255.0/256,0.,0.), normed=0, label='PAR2-G2019S', alpha=0.5)        
#                axes[i,j].hist(CumulativeGenesExpressionB_NORM,bins=15, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), normed=0, label='PAR2-GC', alpha=0.5)
    
                if not (DAY == 14 and GeneName == 'ZIC3'):
                    TestStatistics, pvalue
                    axes[i,j].hist(CumulativeGenesExpressionA_NORM,bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(255.0/256,0.,0.), density = True, label='PAR2-G2019S', alpha=0.5)  #  normed=0,      
                    axes[i,j].hist(CumulativeGenesExpressionB_NORM,bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), density = True, label='PAR2-GC', alpha=0.5) #  normed=0,
                    
    
                if i == 4 and j == 3 :
                    plt.legend()
                                
                BONF = 'For each distribution 250 cells were considered \n Bonferroni correction for multiple hypothesis testing on ' + str(round(NBonferroni)) + ' z-tests (assessing significance of difference of the means) was applied'
                plt.suptitle(BONF)
                
                NBonferroni = 20
                if pvalue < 0.0001 / NBonferroni:
                    PValueLabel = 'p-val < 0.0001'
                    font = font0001
                elif pvalue < 0.001 / NBonferroni:
                    PValueLabel = 'p-val < 0.001'
                    font = font001
                elif pvalue < 0.01 / NBonferroni:
                    PValueLabel = 'p-val < 0.01'
                    font = font01
                elif pvalue < 0.05 / NBonferroni:
                    PValueLabel = 'p-val < 0.05'
                    font = font05
                else:
                    PValueLabel = ' '
                    font = font05
                    MyColorForMeans = 0

                ListOfAveragesUpDown.append(MyColorForMeans)
                
                axes[i,3].text(1.1, 0, GeneName, fontdict=font2)

#                if ijkijk == 1:
#                    axes[0,j].set_ylim(0,30)
#                    axes[1,j].set_ylim(0,20)
#                    axes[2,j].set_ylim(0,70)
#                    axes[3,j].set_ylim(0,45)
#                elif ijkijk == 2:
#                    axes[0,j].set_ylim(0,25)
#                    axes[1,j].set_ylim(0,50)
#                    axes[2,j].set_ylim(0,50)
#                    axes[3,j].set_ylim(0,80)                
#                elif ijkijk == 3:
#                    axes[0,j].set_ylim(0,70)
#                    axes[1,j].set_ylim(0,25)
#                    axes[2,j].set_ylim(0,20)
#                    axes[3,j].set_ylim(0,15)                
#                elif ijkijk == 4:
#                    axes[0,j].set_ylim(0,15)
#                    axes[1,j].set_ylim(0,30)
#                    axes[2,j].set_ylim(0,60)
#                    axes[3,j].set_ylim(0,100)  
                    
                if ijkijk == 1:
                    axes[0,j].set_ylim(0,250)
                    axes[1,j].set_ylim(0,250)
                    axes[2,j].set_ylim(0,250)
                    axes[3,j].set_ylim(0,250)
                    axes[4,j].set_ylim(0,250)
                    axes[i,j].text(0.65, 120, PValueLabel, fontdict=font)

                elif ijkijk == 2:
                    axes[0,j].set_ylim(0,10)
                    axes[1,j].set_ylim(0,10)
                    axes[2,j].set_ylim(0,10)
                    axes[3,j].set_ylim(0,10)
                    axes[4,j].set_ylim(0,10)
                    axes[i,j].text(0.65, 8, PValueLabel, fontdict=font)
                    
            MatrixOfFoldChanges.append(ListOfFoldChangesThisDay)
            MatrixOfSEMforFoldChanges.append(ListOfSTDEVforFoldChangesThisDay)

            ListOfAveragesDays.append(ListOfAveragesUpDown)
            j=j+1
        plt.show()

#        fig2, ax2 = plt.subplots()
##        ax.xaxis.set_visible(True)
##        ax.yaxis.set_visible(True)
#        TransposedListOfAverages = np.array(ListOfAveragesDays).T.tolist()
#        heatmap = ax2.pcolor(TransposedListOfAverages, cmap = mpl.colors.ListedColormap([(128.0/256,128.0/256,128.0/256), 'white', (255.0/256,0.,0.)]), vmin=-1, vmax=1, edgecolors='k') # cmap=plt.cm.RdGy_r # cmap=plt.cm.seismic
#        ax2.invert_yaxis()
#        ax2.text(-0.5, 0.5, CurrentList[0], fontdict=font3s)
#        ax2.text(-0.5, 1.5, CurrentList[1], fontdict=font3s)
#        ax2.text(-0.5, 2.5, CurrentList[2], fontdict=font3s)
#        ax2.text(-0.5, 3.5, CurrentList[3], fontdict=font3s)
#        ax2.text(-0.5, 4.5, CurrentList[4], fontdict=font3s)
#        ax2.text(0.3, -0.3, 'Day 0', fontdict=font3s)
#        ax2.text(1.3, -0.3, 'Day 10', fontdict=font3s)
#        ax2.text(2.3, -0.3, 'Day 14', fontdict=font3s)
#        ax2.text(3.3, -0.3, 'Day 42', fontdict=font3s)
#        #Spacing between each line
#        intervals = 1#float(sys.argv[1])
#        loc = plticker.MultipleLocator(base=intervals)
#        ax2.xaxis.set_major_locator(loc)
#        ax2.yaxis.set_major_locator(loc)
#        # plt.grid(True)
#        # plt.rc('grid', linestyle="-", color='black')
#        #ax2.grid(which="major", color="black", linestyle='-', linewidth=1)
#        ax2.tick_params(labelbottom='off')    
#        ax2.tick_params(labelleft='off')   
#        plt.show()
        # plt.grid(False)
    
#        CountsOfElements = {x:ListOfAveragesDays.count(x) for x in ListOfAveragesDays}
#        
#        CounterControl = 
#        CounterMutant = 
#        CounterNotSignific = 

        ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
        
        N = 4 # Number of days
        
        ind = np.arange(N)  # the x locations for the groups
        width = 0.15       # the width of the bars
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        Gene0FoldChanges = (MatrixOfFoldChanges[0][0], MatrixOfFoldChanges[1][0], MatrixOfFoldChanges[2][0], MatrixOfFoldChanges[3][0])
        # Gene0Std = (0.1, 0.1, 0.1, 0.1)
        Gene0Std = (MatrixOfSEMforFoldChanges[0][0], MatrixOfSEMforFoldChanges[1][0], MatrixOfSEMforFoldChanges[2][0], MatrixOfSEMforFoldChanges[3][0])
        rects1 = ax.bar(ind, Gene0FoldChanges, width, color='r', yerr=Gene0Std)

        Gene1FoldChanges = (MatrixOfFoldChanges[0][1], MatrixOfFoldChanges[1][1], MatrixOfFoldChanges[2][1], MatrixOfFoldChanges[3][1])
        # Gene1Std = (0.1, 0.1, 0.1, 0.1)
        Gene1Std =  (MatrixOfSEMforFoldChanges[0][1], MatrixOfSEMforFoldChanges[1][1], MatrixOfSEMforFoldChanges[2][1], MatrixOfSEMforFoldChanges[3][1]) #  (0.1, 0.1, 0.1, 0.1)
        rects2 = ax.bar(ind+width, Gene1FoldChanges, width, color='y', yerr=Gene1Std)

        Gene2FoldChanges = (MatrixOfFoldChanges[0][2], MatrixOfFoldChanges[1][2], MatrixOfFoldChanges[2][2], MatrixOfFoldChanges[3][2])
        # Gene2Std = (0.1, 0.1, 0.1, 0.1)
        Gene2Std = (MatrixOfSEMforFoldChanges[0][2], MatrixOfSEMforFoldChanges[1][2], MatrixOfSEMforFoldChanges[2][2], MatrixOfSEMforFoldChanges[3][2])
        rects3 = ax.bar(ind+width+width, Gene2FoldChanges, width, color='b', yerr=Gene2Std)
        
        Gene3FoldChanges = (MatrixOfFoldChanges[0][3], MatrixOfFoldChanges[1][3], MatrixOfFoldChanges[2][3], MatrixOfFoldChanges[3][3])
        # Gene3Std = (0.1, 0.1, 0.1, 0.1)
        Gene3Std = (MatrixOfSEMforFoldChanges[0][3], MatrixOfSEMforFoldChanges[1][3], MatrixOfSEMforFoldChanges[2][3], MatrixOfSEMforFoldChanges[3][3])
        rects4 = ax.bar(ind+width+width+width, Gene3FoldChanges, width, color='c', yerr=Gene3Std)
        
        Gene4FoldChanges = (MatrixOfFoldChanges[0][4], MatrixOfFoldChanges[1][4], MatrixOfFoldChanges[2][4], MatrixOfFoldChanges[3][4])
        # Gene4Std = (0.1, 0.1, 0.1, 0.1)
        Gene4Std = (MatrixOfSEMforFoldChanges[0][4], MatrixOfSEMforFoldChanges[1][4], MatrixOfSEMforFoldChanges[2][4], MatrixOfSEMforFoldChanges[3][4])
        rects5 = ax.bar(ind+width+width+width+width, Gene4FoldChanges, width, color='k', yerr=Gene4Std)
        
        # add some
        ax.set_ylabel(r'Fold change ($log_2(\langle mutant \rangle / \langle wt \rangle)$)')
        # ax.set_title('Fold changes for CRC core genes')
        ax.set_xticks(ind+width)
        ax.set_xticklabels( ('Day 0', 'Day 10', 'Day 14', 'Day 42') )

        ax.legend( (rects1[0], rects2[0], rects3[0], rects4[0], rects5[0]), ('ZIC1', 'ZIC2', 'ZIC3', 'ZIC4', 'ZIC5') )

        plt.show()

        def autolabel(rects, values, GeneFoldChanges, GeneStd):
            # attach some text labels
            i = 0
            for rect in rects:
                if GeneFoldChanges[i] > 0:
                    height = GeneFoldChanges[i]
                    Displacement = GeneStd[i] - 0.1
                elif GeneFoldChanges[i] < 0:
                    height = GeneFoldChanges[i]
                    Displacement = - GeneStd[i] - 0.35
                ax.text(rect.get_x()+rect.get_width()/2., height + Displacement, values[i], #'%d'%int(height),
                        ha='center', va='bottom')
                i=i+1
        
        autolabel(rects1,('*',' ',' ','****'), Gene0FoldChanges, Gene0Std)
        autolabel(rects2,(' ',' ',' ','****'), Gene1FoldChanges, Gene1Std)
        autolabel(rects3,(' ',' ',' ','**'), Gene2FoldChanges, Gene2Std) # # = No Data, $ = MeanA = 0
        autolabel(rects4,(' ',' ',' ',' '), Gene3FoldChanges, Gene3Std)
        autolabel(rects5,(' ',' ',' ','***'), Gene4FoldChanges, Gene4Std)
        print(' ')
        print('ATTENTION!!! SIGNIFICANCES INDICATED WITH ASTERISKS IN THE FOLD CHANGE PLOTS ARE HARD CODED, SO THEY WILL NOT CHANGE IF YOU CHANGE SOME PARAMETERS OF THE ANALYSIS!!!')
        print(' ')



if DoCorrelation_ZIC_12345_Lasse_Jocken == 1:
    WhiteDiagonal = 1 # 1 --> mask values (1) from the diagonal, 0 --> do not mask
    RemoveNotSignificantCorr = 1
    N_Bonf_Corr = 80 # Divide by Bonferroni correction N=80?
    ThrasoldForPval = 0.05/N_Bonf_Corr
    
    ######### ALL CELLS TOGETHER #########
    
    normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
    scaClist = [C0_sca,C10_sca,C14_sca,C42_sca]
    Clist = [C0,C10,C14,C42]
    genesClist = []
    
    #L_CoreRegulatoryGenesLasse = sd.get_from_txt('./ListsOfGenesLasse/20180427-CRC-ExtGenes-unique-AllCellTypes_MODIFIED.txt')[0]
    #L_CoreRegulatoryGenesLasse = sd.get_from_txt('./ListsOfGenesLasse/20180427-CRC-CoreTFs-unique-AllCellTypes.txt')[0]
    L_CoreRegulatoryGenesLasse = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJochen/list_ZIC_GENES.txt')[0]
    
    MyAverageCorrVector = []
    for ControlMutant in [0, 250]:
        k = 0
        for DAY in [0,10,14,42]:
            C_norm = normClist[k]
            C_sca = scaClist[k]
            C = Clist[k]
            k = k+1
            
            CurrentList = []
            N = len(L_CoreRegulatoryGenesLasse)
            pvalPearsonCorr_Matrix = 100 * np.ones([N, N])
            MyCorrelationMatrix_Full = 100 * np.ones([N, N])
            MyCorrelationMatrix_Masked = 100 * np.ones([N, N])

            i = 0
            for GeneName1 in L_CoreRegulatoryGenesLasse:
                j = 0
                CurrentList.append(GeneName1)
                for GeneName2 in L_CoreRegulatoryGenesLasse:
                    if not (DAY == 14 and (GeneName1 == 'ZIC3' or GeneName2 == 'ZIC3')):
                        Gene_1_Expression = C_norm[C.c_names.index(GeneName1),0+ControlMutant:250+ControlMutant]
                        Gene_2_Expression = C_norm[C.c_names.index(GeneName2),0+ControlMutant:250+ControlMutant]
                    else:
                        Gene_1_Expression = numpy.zeros(250)
                        Gene_2_Expression = numpy.zeros(250)

        #            Gene_1_ExpressionB = C_norm[C.c_names.index(GeneName1),250:500]
        #            Gene_2_ExpressionB = C_norm[C.c_names.index(GeneName2),250:500]
                    
                    PearsonCorrCoeff, pvalPearsonCorr = pearsonr(Gene_1_Expression, Gene_2_Expression)
                    pvalPearsonCorr_Matrix[i,j] = pvalPearsonCorr
                    
                    if RemoveNotSignificantCorr == 0:
                        MyCorrelationMatrix_Full[i,j] = PearsonCorrCoeff
                                      
                        if WhiteDiagonal == 1:
                            MyCorrelationMatrix_Masked[i,j] = PearsonCorrCoeff
                            if i == j:
                                MyCorrelationMatrix_Masked[i,j] = 0
                    elif RemoveNotSignificantCorr == 1:
                        if pvalPearsonCorr <= ThrasoldForPval:
                            MyCorrelationMatrix_Full[i,j] = PearsonCorrCoeff
                        elif pvalPearsonCorr > ThrasoldForPval:
                            MyCorrelationMatrix_Full[i,j] = 0
                                      
                        if WhiteDiagonal == 1:
                            MyCorrelationMatrix_Masked[i,j] = MyCorrelationMatrix_Full[i,j]
                            if i == j:
                                MyCorrelationMatrix_Masked[i,j] = 0
                                                      
                                                      
                                                      
#                    print(GeneName1)
#                    print(GeneName2)
#                    print(PearsonCorrCoeff)
#                    print(' ')
                    
                    j = j + 1
                i = i + 1
        
            fig, ax = plt.subplots()
            #TransposedListOfAverages = np.array(ListOfAveragesDays).T.tolist()
            #sns.heatmap(MyCorrelationMatrix, annot=False,  linewidths=.5, cmap=plt.cm.seismic, vmin=-1, vmax=1, ax = ax)
            #ax.set_xticks() sns.heatmap(data, annot=True, linewidths=.5)
            
            if WhiteDiagonal == 0:
                MyCorrelationMatrix = MyCorrelationMatrix_Full
            elif WhiteDiagonal == 1:
                MyCorrelationMatrix = MyCorrelationMatrix_Masked
            else:
                print("Error in Masking Correlations!")
            
            MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix)
            for i in [0,1,2,3,4]:
                for j in [0,1,2,3,4]:
                    if j>=i:
                        MyCorrelationMatrix_Triangle[i,j] = 0
            
            heatmap = ax.pcolor(MyCorrelationMatrix_Triangle, cmap=plt.cm.seismic_r, vmin=-1, vmax=1) # RdGy_r  # PiYG  # seismic   # PRGn_r
            ax.invert_yaxis()
            ax.text(-0.7, 0.5, CurrentList[0], fontdict=font3)
            ax.text(-0.7, 1.5, CurrentList[1], fontdict=font3)
            ax.text(-0.7, 2.5, CurrentList[2], fontdict=font3)
            ax.text(-0.7, 3.5, CurrentList[3], fontdict=font3)
            ax.text(-0.7, 4.5, CurrentList[4], fontdict=font3)
            ax.text(0.3, -0.3, CurrentList[0], fontdict=font3)
            ax.text(1.3, -0.3, CurrentList[1], fontdict=font3)
            ax.text(2.3, -0.3, CurrentList[2], fontdict=font3)
            ax.text(3.3, -0.3, CurrentList[3], fontdict=font3)
            ax.text(4.3, -0.3, CurrentList[4], fontdict=font3)
            #Spacing between each line
            intervals = 1#float(sys.argv[1])
            loc = plticker.MultipleLocator(base=intervals)
            ax.xaxis.set_major_locator(loc)
            ax.yaxis.set_major_locator(loc)
            ax.grid(which="major", color="black", linestyle='-', linewidth=1)
            ax.tick_params(labelbottom='off')    
            ax.tick_params(labelleft='off')    
            fig.colorbar(heatmap)
            i = 0
            j = 0
            for i in [0,1,2,3,4]:
                for j in [0,1,2,3,4]:
                    if j>=i:
                        if i != j:
                            if MyCorrelationMatrix[i,j] != 0:
                                Current_pvalPearsonCorr = pvalPearsonCorr_Matrix[i,j]
                                if Current_pvalPearsonCorr <= 0.001/N_Bonf_Corr:
                                    Significance = ' ***'
                                elif Current_pvalPearsonCorr <= 0.01/N_Bonf_Corr:
                                    Significance = ' **'
                                elif Current_pvalPearsonCorr <= 0.05/N_Bonf_Corr:
                                    Significance = ' *'
                                MyText = str(round(MyCorrelationMatrix[i,j],2)) + Significance
                                ax.text(i+0.5, j+0.5, MyText, 
                                    horizontalalignment='center', verticalalignment='center', fontdict=font3)
                            elif MyCorrelationMatrix[i,j] == 0:
                                ax.text(i+0.5, j+0.5, r'$\approx 0$', 
                                    horizontalalignment='center', verticalalignment='center', fontdict=font3)
                        elif i == j:
                            if WhiteDiagonal == 0:
                                ax.text(i+0.5, j+0.5, round(MyCorrelationMatrix[i,j],2), 
                                        horizontalalignment='center', verticalalignment='center', fontdict=font3)
                            elif WhiteDiagonal == 1:
                                ax.text(i+0.5, j+0.5,' ', 
                                        horizontalalignment='center', verticalalignment='center', fontdict=font3)
            plt.show()
            
            Dim = len(MyCorrelationMatrix)
            if WhiteDiagonal == 0:
                AverageCorr = (sum(sum(abs(MyCorrelationMatrix)))-Dim)/(Dim*Dim-Dim)
            elif WhiteDiagonal == 1: # Note that MyCorrelationMatrix contains the diagonal here above, doesn't here below.
                AverageCorr = (sum(sum(abs(MyCorrelationMatrix))))/(Dim*Dim-Dim)
            print(AverageCorr)
            print(' ')
            MyAverageCorrVector.append(AverageCorr)            
            
            
    plt.figure()      
    plt.plot([0,10,14,42], MyAverageCorrVector[0:4], marker='o', linestyle='--', color=(255.0/256,0.,0.), label='G2019S')
    plt.plot([0,10,14,42], MyAverageCorrVector[4:8], marker='s', linestyle='-', color=(128.0/256,128.0/256,128.0/256), label='GC')
    plt.xlabel('Time (Days after differentiation)')
    plt.ylabel(r'Average Gene-Gene Correlation')
    plt.ylim(0,0.3)
    plt.title('Average of the absolute value of the Gene-Gene Pearson Correlation Coefficient')
    plt.legend()
    plt.show()



###############################################################################
###############################################################################
    

if ComputeIntersectionGenesLists == 1:
        
        print("Compute intersection of lists of relevant genes...")
    
        L_Aall=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_stem.txt')[0]
        L_Ball=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_DopaminergicNeurons.txt')[0]
        L_Call=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_CellCycle.txt')[0]
        L_Dall=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Mito.txt')[0]
        L_Eall=sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_ROS.txt')[0]
        
        L_all = []
        for GeneName in L_Aall:
            L_all.append(GeneName)
        for GeneName in L_Ball:
            L_all.append(GeneName)
        for GeneName in L_Call:
            L_all.append(GeneName)
        for GeneName in L_Dall:
            L_all.append(GeneName)
        for GeneName in L_Eall:
            L_all.append(GeneName)

        L_A = []
        L_B = []
        L_C = []
        L_D = []
        L_E = []
        
        L_AB = []
        L_AC = [] 
        L_AD = []
        L_AE = []  
        L_BC = []
        L_BD = []
        L_BE = []
        L_CD = []
        L_CE = []
        L_DE = []

        L_ABC = []
        L_ABD = [] 
        L_ABE = []
        L_ACD = []###
        L_ACE = []###
        L_ADE = []###
        L_BCD = []
        L_BCE = []
        L_BDE = []
        L_CDE = []
                
        L_ABCD = []#
        L_ABCE = []#
        L_ABDE = []# ###
        L_ACDE = []#
        L_BCDE = []# ###
        
        L_ABCDE = []#
        
        for GeneName in L_all:
            if (GeneName in L_Aall) and (GeneName in L_Ball) and (GeneName in L_Call) and (GeneName in L_Dall) and (GeneName in L_Eall):
                L_ABCDE.append(GeneName)
            elif (GeneName in L_Aall) and (GeneName in L_Call) and (GeneName in L_Dall) and (GeneName in L_Eall):
                L_ACDE.append(GeneName)
            elif (GeneName in L_Aall) and (GeneName in L_Ball) and (GeneName in L_Call) and (GeneName in L_Eall):
                L_ABCE.append(GeneName)
            elif (GeneName in L_Aall) and (GeneName in L_Ball) and (GeneName in L_Call) and (GeneName in L_Dall):
                L_ABCD.append(GeneName)
            elif (GeneName in L_Aall) and (GeneName in L_Ball) and (GeneName in L_Dall) and (GeneName in L_Eall):
                L_ABDE.append(GeneName)
            elif (GeneName in L_Ball) and (GeneName in L_Call) and (GeneName in L_Dall) and (GeneName in L_Eall):
                L_BCDE.append(GeneName)
            elif (GeneName in L_Call) and (GeneName in L_Dall) and (GeneName in L_Eall):
                L_CDE.append(GeneName)
            elif (GeneName in L_Ball) and (GeneName in L_Dall) and (GeneName in L_Eall):
                L_BDE.append(GeneName)
            elif (GeneName in L_Ball) and (GeneName in L_Call) and (GeneName in L_Eall):
                L_BCE.append(GeneName)
            elif (GeneName in L_Ball) and (GeneName in L_Call) and (GeneName in L_Dall):
                L_BCD.append(GeneName)
            elif (GeneName in L_Aall) and (GeneName in L_Ball) and (GeneName in L_Eall):
                L_ABE.append(GeneName)
            elif (GeneName in L_Aall) and (GeneName in L_Ball) and (GeneName in L_Dall):
                L_ABD.append(GeneName)
            elif (GeneName in L_Aall) and (GeneName in L_Ball) and (GeneName in L_Call):
                L_ABC.append(GeneName)
            elif (GeneName in L_Aall) and (GeneName in L_Call) and (GeneName in L_Dall):
                L_ACD.append(GeneName)
            elif (GeneName in L_Aall) and (GeneName in L_Call) and (GeneName in L_Eall):
                L_ACE.append(GeneName)
            elif (GeneName in L_Aall) and (GeneName in L_Dall) and (GeneName in L_Eall):
                L_ADE.append(GeneName)
            elif (GeneName in L_Dall) and (GeneName in L_Eall):
                L_DE.append(GeneName)
            elif (GeneName in L_Call) and (GeneName in L_Eall):
                L_CE.append(GeneName)
            elif (GeneName in L_Call) and (GeneName in L_Dall):
                L_CD.append(GeneName)
            elif (GeneName in L_Ball) and (GeneName in L_Eall):
                L_BE.append(GeneName)
            elif (GeneName in L_Ball) and (GeneName in L_Dall):
                L_BD.append(GeneName)
            elif (GeneName in L_Ball) and (GeneName in L_Call):
                L_BC.append(GeneName)
            elif (GeneName in L_Aall) and (GeneName in L_Eall):
                L_AE.append(GeneName)
            elif (GeneName in L_Aall) and (GeneName in L_Dall):
                L_AD.append(GeneName)
            elif (GeneName in L_Aall) and (GeneName in L_Call):
                L_AC.append(GeneName)
            elif (GeneName in L_Aall) and (GeneName in L_Ball):
                L_AB.append(GeneName)
            elif (GeneName in L_Aall):
                L_A.append(GeneName)
            elif (GeneName in L_Ball):
                L_B.append(GeneName)
            elif (GeneName in L_Call):
                L_C.append(GeneName)
            elif (GeneName in L_Dall):
                L_D.append(GeneName)
            elif (GeneName in L_Eall):
                L_E.append(GeneName)
                
#        for MyList in [L_AB,L_AC,L_AD,L_AE,L_BC,L_BD,L_BE,
#                     L_CD,L_CE,L_DE,L_ABC,L_ABD,L_ABE,L_ACD,L_ACE,L_ADE,L_BCD,
#                     L_BCE,L_BDE,L_CDE,L_ABCD,L_ABCE,L_ABDE,L_ACDE,L_BCDE,L_ABCDE]:
     
        if not(len(L_CD)==0):
            ListAux = L_CD
            L_CD = []
            for GeneName in ListAux:
                if not(GeneName in L_CD):
                    #print(GeneName)
                    #rint(L_CD)
                    L_CD.append(GeneName)
       
        if not(len(L_DE)==0):
            ListAux = L_DE
            L_DE = []
            for GeneName in ListAux:
                if not(GeneName in L_DE):
                    #print(GeneName)
                    #print(L_DE)
                    L_DE.append(GeneName)
                
        print("This computes the numbers of genes in the intersections between lists of relevant genes (Stem, DA neuro, cell cycle, mito ,ROS).")
        print(" ")                
        print(len(L_all))
        print(" ")
        print(len(L_Aall))
        print(len(L_Ball))
        print(len(L_Call))
        print(len(L_Dall))
        print(len(L_Eall))
        print(" ")
        print(len(L_A))
        print(len(L_B))
        print(len(L_C))
        print(len(L_D))
        print(len(L_E))
        print(" ")
        print(len(L_AB))
        print(len(L_AC)) 
        print(len(L_AD))
        print(len(L_AE)) 
        print(len(L_BC))
        print(len(L_BD))
        print(len(L_BE))
        print(len(L_CD))
        print(len(L_CE))
        print(len(L_DE))
        print(" ")
        print(len(L_ABC))
        print(len(L_ABD)) 
        print(len(L_ABE))
        print(len(L_BCD))
        print(len(L_BCE))
        print(len(L_BDE))
        print(len(L_CDE))
        print(len(L_ACD))
        print(len(L_ACE))
        print(len(L_ADE))
        print(" ")
        print(len(L_ABCD))
        print(len(L_ABCE))
        print(len(L_ACDE))
        print(len(L_ABDE))
        print(len(L_BCDE))
        print(" ")
        print(len(L_ABCDE))


if ComputeIntersectionGenesListsDEG_MitoOnly_Days == 1:

    print("Compute intersection of genes lists for differentially expressed genes and mio genes, for mito paper...")
    
    L_DEG_Mito_day0 = sd.get_from_txt('./Results/TabS2_DiffExprGenes_MITO_day0.txt')[0]
    L_DEG_Mito_day10 = sd.get_from_txt('./Results/TabS2_DiffExprGenes_MITO_day10.txt')[0]
    L_DEG_Mito_day14 = sd.get_from_txt('./Results/TabS2_DiffExprGenes_MITO_day14.txt')[0]
    L_DEG_Mito_day42 = sd.get_from_txt('./Results/TabS2_DiffExprGenes_MITO_day42.txt')[0]
    
    L_Mito_all = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Mito.txt')[0]

    print(len(L_DEG_Mito_day0))
    print(len(L_DEG_Mito_day10))
    print(len(L_DEG_Mito_day14))
    print(len(L_DEG_Mito_day42))
    print(" ")
    print(len(L_Mito_all))
    print(" ")
    print(len(L_DEG_Mito_day0)/len(L_Mito_all))
    print(len(L_DEG_Mito_day10)/len(L_Mito_all))
    print(len(L_DEG_Mito_day14)/len(L_Mito_all))
    print(len(L_DEG_Mito_day42)/len(L_Mito_all))
    
    L_Aall = L_DEG_Mito_day0
    L_Ball = L_DEG_Mito_day10
    L_Call = L_DEG_Mito_day14
    L_Dall = L_DEG_Mito_day42
    
    L_all = []
    for GeneName in L_Aall:
        L_all.append(GeneName)
    for GeneName in L_Ball:
        L_all.append(GeneName)
    for GeneName in L_Call:
        L_all.append(GeneName)
    for GeneName in L_Dall:
        L_all.append(GeneName)

    L_all = list(set(L_all))

    L_A = []
    L_B = []
    L_C = []
    L_D = []
    
    L_AB = []
    L_AC = [] 
    L_AD = []
    L_BC = []
    L_BD = []
    L_CD = []

    L_ABC = []
    L_ABD = [] 
    L_ACD = []
    L_BCD = []
            
    L_ABCD = []

    for GeneName in L_all:
        if (GeneName in L_Aall) and (GeneName in L_Ball) and (GeneName in L_Call) and (GeneName in L_Dall):
            L_ABCD.append(GeneName)
        elif (GeneName in L_Ball) and (GeneName in L_Call) and (GeneName in L_Dall):
            L_BCD.append(GeneName)
        elif (GeneName in L_Aall) and (GeneName in L_Ball) and (GeneName in L_Dall):
            L_ABD.append(GeneName)
        elif (GeneName in L_Aall) and (GeneName in L_Ball) and (GeneName in L_Call):
            L_ABC.append(GeneName)
        elif (GeneName in L_Aall) and (GeneName in L_Call) and (GeneName in L_Dall):
            L_ACD.append(GeneName)
        elif (GeneName in L_Call) and (GeneName in L_Dall):
            L_CD.append(GeneName)
        elif (GeneName in L_Ball) and (GeneName in L_Dall):
            L_BD.append(GeneName)
        elif (GeneName in L_Ball) and (GeneName in L_Call):
            L_BC.append(GeneName)
        elif (GeneName in L_Aall) and (GeneName in L_Dall):
            L_AD.append(GeneName)
        elif (GeneName in L_Aall) and (GeneName in L_Call):
            L_AC.append(GeneName)
        elif (GeneName in L_Aall) and (GeneName in L_Ball):
            L_AB.append(GeneName)
        elif (GeneName in L_Aall):
            L_A.append(GeneName)
        elif (GeneName in L_Ball):
            L_B.append(GeneName)
        elif (GeneName in L_Call):
            L_C.append(GeneName)
        elif (GeneName in L_Dall):
            L_D.append(GeneName)
            
    print("This computes the numbers of genes in the intersections between Mito genes list and DEGs at various days, it is used to generate the numbers in Mito paper (https://doi.org/10.1016/j.stemcr.2019.03.004) Fig1D!!!")

    print(" ")                
    print(len(L_all))
    print(" ")
    print(len(L_Aall))
    print(len(L_Ball))
    print(len(L_Call))
    print(len(L_Dall))
    print(" ")
    print(len(L_A))
    print(len(L_B))
    print(len(L_C))
    print(len(L_D))
    print(" ")
    print(len(L_AB))
    print(len(L_AC)) 
    print(len(L_AD))
    print(len(L_BC))
    print(len(L_BD))
    print(len(L_CD))
    print(" ")
    print(len(L_ABC))
    print(len(L_ABD)) 
    print(len(L_BCD))
    print(len(L_ACD))
    print(" ")
    print(len(L_ABCD))
    

if HeatMapOfRowExpressionData == 1:
    
    # REQUIRES: "DoComparisonPreliminarySteps"
    
    """
    print("aaa 1")
    A_0_250_genes, A_0_250_expr, A_0_250_expr_np = sd.get_from_txt('./DataJonas/A_0.txt')
    B_0_250_genes, B_0_250_expr, B_0_250_expr_np = sd.get_from_txt('./DataJonas/B_0.txt')  
    A_10_250_genes, A_10_250_expr, A_10_250_expr_np = sd.get_from_txt('./DataJonas/A_10.txt')
    B_10_250_genes, B_10_250_expr, B_10_250_expr_np = sd.get_from_txt('./DataJonas/B_10.txt')
    A_14_250_genes, A_14_250_expr, A_14_250_expr_np = sd.get_from_txt('./DataJonas/A_14.txt')
    B_14_250_genes, B_14_250_expr, B_14_250_expr_np = sd.get_from_txt('./DataJonas/B_14.txt')    
    A_42_250_genes, A_42_250_expr, A_42_250_expr_np = sd.get_from_txt('./DataJonas/A_42.txt')
    B_42_250_genes, B_42_250_expr, B_42_250_expr_np = sd.get_from_txt('./DataJonas/B_42.txt')
    print("aaa 2")
    A_0_250 = sd.get_top_cells_matrix(A_0_250_expr_np,250)[0]
    B_0_250 = sd.get_top_cells_matrix(B_0_250_expr_np,250)[0]
    A_0_250_norm, A_0_250_log, A_0_250_sca = sd.get_normalization(A_0_250)
    B_0_250_norm, B_0_250_log, B_0_250_sca = sd.get_normalization(B_0_250)
    print("aaa 3")
 
    A_10_250 = sd.get_top_cells_matrix(A_10_250_expr_np,250)[0]
    B_10_250 = sd.get_top_cells_matrix(B_10_250_expr_np,250)[0]
    A_10_250_norm, A_10_250_log, A_10_250_sca = sd.get_normalization(A_10_250)
    B_10_250_norm, B_10_250_log, B_10_250_sca = sd.get_normalization(B_10_250)
    print("aaa 4")
 
    A_14_250 = sd.get_top_cells_matrix(A_14_250_expr_np,250)[0]
    B_14_250 = sd.get_top_cells_matrix(B_14_250_expr_np,250)[0]
    A_14_250_norm, A_14_250_log, A_14_250_sca = sd.get_normalization(A_14_250)
    B_14_250_norm, B_14_250_log, B_14_250_sca = sd.get_normalization(B_14_250)
    print("aaa 5")

    A_42_250 = sd.get_top_cells_matrix(A_42_250_expr_np,250)[0]
    B_42_250 = sd.get_top_cells_matrix(B_42_250_expr_np,250)[0]
    A_42_250_norm, A_42_250_log, A_42_250_sca = sd.get_normalization(A_42_250)
    B_42_250_norm, B_42_250_log, B_42_250_sca = sd.get_normalization(B_42_250)
    print("aaa 6")
    """
#    import numpy
#    column_labels = list('ABCD')
#    row_labels = list('WXYZ')
#    data = numpy.random.rand(4,4)
#    from matplotlib import pyplot as plt
#    heatmap = plt.pcolor(data)

    
    #Clist = [C0,C10,C14,C42]    
    scaClist = [C0_sca,C10_sca,C14_sca,C42_sca]
        
    i = 0
    for DAY in [0,10,14,42]:
        for ControlMutant in [0,250]:
            #C = Clist[i]
            C_sca = scaClist[i]
        
            #column_labels = list('ABCD')
            #row_labels = list('WXYZ')
            data = C_sca[:,0+ControlMutant:250+ControlMutant] # A_0_250_sca # np.random.rand(20000,250)
            fig, ax = plt.subplots()
            heatmap = ax.pcolor(data, cmap=plt.cm.seismic)
            
        #    C_norm[C.c_names.index(GeneName),0:250]
        #    C_norm[C.c_names.index(GeneName),250:500]
            
        #    # put the major ticks at the middle of each cell
        #    ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
        #    ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)
        #    
        #    # want a more natural, table-like display
        #    ax.invert_yaxis()
        #    ax.xaxis.tick_top()
        #    
        #    ax.set_xticklabels(row_labels, minor=False)
        #    ax.set_yticklabels(column_labels, minor=False)
            plt.show()
            if ControlMutant == 0:
                CM = 'G2019S'
            elif ControlMutant == 250:
                CM = 'GC'
            plt.savefig('./Results/jonas_ExpressionMatrix_' + CM + '_day' + str(DAY) + '.jpeg')
        i = i+1
        
        
        
        
        
        
######################################################################
########### Now after Jonas analysis perform Jochen analysis#########
######################################################################
        
if JochenAnalysis == 1:
    if CumulativeGeneExpressionCRC == 1:
        print("Here I am JOCHEN 1")

        import matplotlib
    
        normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
        scaClist = [C0_sca,C10_sca,C14_sca,C42_sca]
        Clist = [C0,C10,C14,C42]
            
        ### Start cumulative plot
        Nlines = 3
        Ncols = 4
        fig, axes = plt.subplots(nrows=Nlines, ncols=Ncols)
        fig.subplots_adjust(hspace=0.3, wspace=0.05)
    
        for ax in axes.flat:
            # Hide all ticks and labels
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            #ax.tick_params(axis=u'x', which=u'x',length=0)
            #fig.setp(ax.get_xticklabels(),visible=False)
    
            # Set up ticks only on one side for the "edge" subplots...
#            if ax.is_last_col():
#                if ax.is_first_row():
#                    #ax.set_title("Title 1")
#                    ax.annotate("DA neurons\nJonas",xy=(5, 0.20))
#                elif ax.is_last_row():
#                    ax.annotate("core genes\nJochen",xy=(15, 0.40))
#                else:
#                    ax.annotate("CRC\nJochen", xy=(0.01, 0.01))
            if ax.is_first_col():
                ax.yaxis.set_ticks_position('left')
                ax.yaxis.set_visible(True)
    #        if ax.is_last_col():
    #            ax.yaxis.set_ticks_position('right')
    #            ax.yaxis.set_label_position('right')
    #        if ax.is_first_row():
    #            ax.xaxis.set_ticks_position('top')
    #            ax.xaxis.set_label_position('top')
            if ax.is_last_row():
               #ax.xaxis.set_ticks_position('bottom')
               ax.xaxis.set_visible(True)
               ax.set_xticklabels([])
               ax.xaxis.set_ticks_position('none')
               #fig.setp(ax.get_xticklabels(),visible=True)
    
        L_JonasDA = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_DopaminergicNeurons.txt')[0]
        L_JochenCRC = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJochen/list_core_genes_JO_DAneurons_CRC.txt')[0]
        L_JochenInnerCore = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJochen/list_core_genes_JO_DAneurons_innercore.txt')[0]

        j = 0
        for DAY in [0,10,14,42]:
            print('Day ' + str(DAY))
            print()
            C = Clist[j]
            C_norm = normClist[j]
            
            ListSelectedGenes = []
            for i in range(Nlines):
                CumulativeGenesExpressionA = np.zeros(250)
                CumulativeGenesExpressionB = np.zeros(250)            
                if i == 0:
                    ListSelectedGenes = L_JonasDA
                    ymax = 0.5
                elif i == 1:
                    ListSelectedGenes = L_JochenCRC
                    ymax = 0.035
                elif i == 2:

                    ListSelectedGenes = L_JochenInnerCore
                    ymax = 1.0
                    if j == 3:
                        plt.legend()
    
                axes[i,0].set_ylabel("Frequency")
                axes[2,j].set_xlabel("Cumulative Gene Expression")
                    
                for GeneName in ListSelectedGenes:
                    if GeneName in C.c_names:
                        CumulativeGenesExpressionA = CumulativeGenesExpressionA + C_norm[C.c_names.index(GeneName),0:250]
                        CumulativeGenesExpressionB = CumulativeGenesExpressionB + C_norm[C.c_names.index(GeneName),250:500]

                axes[0,0].set_title("Day 0")
                axes[0,1].set_title("Day 10")
                axes[0,2].set_title("Day 14")
                axes[0,3].set_title("Day 42")

                axes[0,3].annotate("DA Neuron Genes\nJonas",xy=(1.05, 0.4),xycoords='axes fraction', size='large')
                axes[1,3].annotate("CRC Genes\nJochen", xy=(1.05, 0.4),xycoords='axes fraction', size='large')
                axes[2,3].annotate("Core Genes\nJochen",xy=(1.05, 0.4),xycoords='axes fraction', size='large')

#                axes[0,3].annotate("DA neurons\nJonas",xy=(5, 0.20))
#                axes[1,3].annotate("CRC\nJochen", xy=(150, 0.01))
#                axes[2,3].annotate("core genes\nJochen",xy=(15, 0.40))
                            
                ########## MAKE HISTOGRAMS ##########
                axes[i,j].hist(CumulativeGenesExpressionA,bins=15, histtype='stepfilled', color=(255.0/256,0.,0.), normed=1, label='PAR2-G2019S', alpha=0.5)        
                axes[i,j].hist(CumulativeGenesExpressionB,bins=15, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), normed=1, label='PAR2-GC', alpha=0.5)
        
                #axes[i,j].set_xlim(0,xmax)
                xmax = max(max(CumulativeGenesExpressionA),max(CumulativeGenesExpressionB))
                axes[i,j].set_ylim(0,ymax)
                axes[i,j].set_xlim(0,xmax)
    
                ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
                # of 2 populations are statistically significantly different ######
                
                import statsmodels.stats.weightstats as wstats
                #DataPopulation1 = numpy.array([1,2,3,4,5,2,3,4,3])
                DataPopulation1 = CumulativeGenesExpressionA
                #DataPopulation2 = numpy.array([1,2,3.001,4,5,2,3,4,3])
                DataPopulation2 = CumulativeGenesExpressionB
                DataPopulation1_Object = wstats.DescrStatsW(DataPopulation1)
                DataPopulation2_Object = wstats.DescrStatsW(DataPopulation2)
                MyComparison = wstats.CompareMeans(DataPopulation1_Object,DataPopulation2_Object)
                TestStatistics, pvalue = MyComparison.ztest_ind(alternative='two-sided', usevar='unequal', value=0)
    
                # Remember Bonferroni correction for multiple testing, 
                # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
                # by the number of repetitions of the test, 4 days x 7 lists = 28
                # The 7 lists are Stem, DAneur, Mito, Cell Cycle, Pro-apoptosis, Anti-apoptosis, Caspases
                print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
                print(' ')
        
                plt.legend()
        
            j=j+1
        plt.show()
            
    if DRwithList_CRCgenesOnly_Mutant == 1:
        ##### MUTANT #####        
        
        ### for normalization in [norm, sca, log] 
        
        print("Here I am JOCHEN 2")
        
        L_JonasDA = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_DopaminergicNeurons.txt')[0]
        L_JochenCRC = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJochen/list_core_genes_JO_DAneurons_CRC.txt')[0]
            
        CurrentNormalizationAlist = [A_0_norm,A_10_norm,A_14_norm,A_42_norm] # [A_0_log,A_10_log,A_14_log,A_42_log] #[A_0_sca,A_10_sca,A_14_sca,A_42_sca]
        genesAlist = [A_0_genes,A_10_genes,A_14_genes,A_42_genes]
    #    scaBlist = [B_0_norm,B_10_norm,B_14_norm,B_42_norm]#[B_0_sca,B_10_sca,B_14_sca,B_42_sca]
    #    genesBlist = [B_0_genes,B_10_genes,B_14_genes,B_42_genes]
        
        import matplotlib
            
        Nlines = 2
        Ncols = 4
            
        fig, axes = plt.subplots(nrows=Nlines, ncols=Ncols)
        fig.subplots_adjust(hspace=0.3, wspace=0.05)
        plt.suptitle('tSNE of G2019S (the mutant)')
    
        for ax in axes.flat:
            # Hide all ticks and labels
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            
    #        if ax.is_first_col():
    #            ax.yaxis.set_ticks_position('left')
    #            ax.yaxis.set_visible(True)
    #        if ax.is_last_row():
    #           ax.xaxis.set_visible(True)
    #           ax.set_xticklabels([])
    #           ax.xaxis.set_ticks_position('none')         
        
        ListOftSNE_Mutant_DAneuronsJonasList = []
        ListOftSNE_Mutant_CRCJochenList = []
        j = 0
        for DAY in [0,10,14,42]:
            print()
            print('Day ' + str(DAY))
            print()
    
            for i in range(Nlines):         
                print('List ' + str(i))
                if i == 0:
                    ListSelectedGenes = L_JonasDA
                    ListName = 'Jonas'
                    GeneToColor = 'CALB2'
                elif i == 1:
                    ListSelectedGenes = L_JochenCRC
                    ListName = 'Jochen'
                    GeneToColor = 'MSI2'
    
                #axes[i,0].set_ylabel("X label")
                #axes[1,j].set_xlabel("Y label")
    
                axes[0,0].set_title("Day 0")
                axes[0,1].set_title("Day 10")
                axes[0,2].set_title("Day 14")
                axes[0,3].set_title("Day 42")
                print('1')
                if FortSNEUseAllGenesInsteadOfListDAneurons == True:
                    print('2')
                    axes[0,3].annotate("ALL Genes",xy=(1.05, 0.4),xycoords='axes fraction', size='large')
                    axes[1,3].annotate("ALL Genes", xy=(1.05, 0.4),xycoords='axes fraction', size='large')
                    #axes2[2,3].annotate("Core Genes\nJochen",xy=(1.05, 0.4),xycoords='axes fraction', size='large')
                elif FortSNEUseAllGenesInsteadOfListDAneurons == False:
                    axes[0,3].annotate("DA Neuron Genes\nJonas",xy=(1.05, 0.4),xycoords='axes fraction', size='large')
                    axes[1,3].annotate("CRC Genes\nJochen", xy=(1.05, 0.4),xycoords='axes fraction', size='large')
                    #axes2[2,3].annotate("Core Genes\nJochen",xy=(1.05, 0.4),xycoords='axes fraction', size='large')
                print('3')
    
                ########## MAKE tSNE ##########
                MyAdata = CurrentNormalizationAlist[j]
                print('4')
    
                genes = genesAlist[j]
                if FortSNEUseAllGenesInsteadOfListDAneurons == True:
                    print('5')
                    M = MyAdata
                    ListName = 'ALL'
                    print('6')
                elif FortSNEUseAllGenesInsteadOfListDAneurons == False:
                    M = sd.get_submatrix(ListSelectedGenes, MyAdata, genes)   
                    
                    
                if UsePCApreviousTotSNE == True:
                    print('7')
    
                    # DO PCA BEFORE tSNE ADN USE RESULTS AS IMPUT FOR SUBSEQUENT tSNE!!!
                    y_pca_temp, y_tsne_temp = get_DR(M,g_ind=-1, d=NumberOfPCAbeforetSNE, off=[])
                    print('8')
    
                    Myt1 = 'G2019S,'+str(DAY)+',norm - PCA\npoints colored by expression of '+GeneToColor+'_List '+ListName+'(tSNE PRECEDEED BY PCA)'
                    Myt2 = 'G2019S,'+str(DAY)+',norm - tSNE\npoints colored by expression of '+GeneToColor+'_List '+ListName+'(tSNE PRECEDEED BY PCA)'
                    MySaveAs1 = './Results/Jochen/PCA/Mutant/G2019S_day' + str(DAY) + '_PCA_DAneurons_'+GeneToColor+'_List ' + ListName + '_PCAbeforetSNE.pdf'
                    MySaveAs2 = './Results/Jochen/tSNE/Mutant/G2019S_day' + str(DAY) + '_tSNE_DAneurons_'+GeneToColor+'_List ' + ListName + '_PCAbeforetSNE.pdf'
                    print('9')
                                     
                    y_pca, y_tsne = get_DR(np.transpose(y_pca_temp),g_ind=-1, d=2, off=[],t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)
                    axes[i,j].scatter(y_tsne[:,0],y_tsne[:,1])  
                    print('10')
    
                elif UsePCApreviousTotSNE == False:
                    t1 = 'G2019S,'+str(DAY)+',norm - PCA\npoints colored by expression of '+GeneToColor+'_List '+ListName
                    t2 = 'G2019S,'+str(DAY)+',norm - tSNE\npoints colored by expression of '+GeneToColor+'_List '+ListName
                    SaveAs1 = './Results/Jochen/PCA/Mutant/G2019S_day' + str(DAY) + '_PCA_DAneurons_'+GeneToColor+'_List ' + ListName + '.pdf'
                    SaveAs2 = './Results/Jochen/tSNE/Mutant/G2019S_day' + str(DAY) + '_tSNE_DAneurons_'+GeneToColor+'_List ' + ListName + '.pdf'
                                                                                  
                    y_pca, y_tsne = get_DR(M,genes.index(GeneToColor),2,[],t1,t2,SaveAs1,SaveAs2)
                    axes[i,j].scatter(y_tsne[:,0],y_tsne[:,1])  
                
                if i == 0:
                    ListOftSNE_Mutant_DAneuronsJonasList.append(y_tsne)
                elif i == 1:
                    ListOftSNE_Mutant_CRCJochenList.append(y_tsne)
                    
            j=j+1
            
        plt.show()
        print("Io pure pure!")
    
        
    if DRwithList_CRCgenesOnly_Control == 1:
        ##### CONTROL #####
    
        print("Here I am JOCHEN 3")
        
        L_JonasDA = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_DopaminergicNeurons.txt')[0]
        L_JochenCRC = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJochen/list_core_genes_JO_DAneurons_CRC.txt')[0]
            
    #    scaAlist = [A_0_norm,A_10_norm,A_14_norm,A_42_norm]#[A_0_sca,A_10_sca,A_14_sca,A_42_sca]
    #    genesAlist = [A_0_genes,A_10_genes,A_14_genes,A_42_genes]
        CurrentNormalizationBlist = [B_0_norm,B_10_norm,B_14_norm,B_42_norm]#[B_0_sca,B_10_sca,B_14_sca,B_42_sca]
        genesBlist = [B_0_genes,B_10_genes,B_14_genes,B_42_genes]
                
        Nlines = 2
        Ncols = 4
        
        fig2, axes2 = plt.subplots(nrows=Nlines, ncols=Ncols)
        fig2.subplots_adjust(hspace=0.3, wspace=0.05)
        plt.suptitle('tSNE of GC (the control)')
    
        for ax in axes2.flat:
            # Hide all ticks and labels
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            
    #        if ax.is_first_col():
    #            ax.yaxis.set_ticks_position('left')
    #            ax.yaxis.set_visible(True)
    #        if ax.is_last_row():
    #           ax.xaxis.set_visible(True)
    #           ax.set_xticklabels([])
    #           ax.xaxis.set_ticks_position('none') 
        
        ListOftSNE_Control_DAneuronsJonasList = []
        ListOftSNE_Control_CRCJochenList = []
        j = 0
        for DAY in [0,10,14,42]:
            print()
            print('Day ' + str(DAY))
            print()
    
            for i in range(Nlines):         
                print('List ' + str(i))
                if i == 0:
                    ListSelectedGenes = L_JonasDA
                    ListName = 'Jonas'
                    GeneToColor = 'CALB2'
                elif i == 1:
                    ListSelectedGenes = L_JochenCRC
                    GeneToColor = 'MSI2'
    
                #axes2[i,0].set_ylabel("X label")
                #axes2[1,j].set_xlabel("Y label")
    
                axes2[0,0].set_title("Day 0")
                axes2[0,1].set_title("Day 10")
                axes2[0,2].set_title("Day 14")
                axes2[0,3].set_title("Day 42")
    
                if FortSNEUseAllGenesInsteadOfListDAneurons == True:
                    axes2[0,3].annotate("ALL Genes",xy=(1.05, 0.4),xycoords='axes fraction', size='large')
                    axes2[1,3].annotate("ALL Genes", xy=(1.05, 0.4),xycoords='axes fraction', size='large')
                    #axes2[2,3].annotate("Core Genes\nJochen",xy=(1.05, 0.4),xycoords='axes fraction', size='large')
                elif FortSNEUseAllGenesInsteadOfListDAneurons == False:
                    axes2[0,3].annotate("DA Neuron Genes\nJonas",xy=(1.05, 0.4),xycoords='axes fraction', size='large')
                    axes2[1,3].annotate("CRC Genes\nJochen", xy=(1.05, 0.4),xycoords='axes fraction', size='large')
                    #axes2[2,3].annotate("Core Genes\nJochen",xy=(1.05, 0.4),xycoords='axes fraction', size='large')        
                
                
                ########## MAKE tSNE ##########
                MyBdata = CurrentNormalizationBlist[j]
                genes = genesBlist[j]
                if FortSNEUseAllGenesInsteadOfListDAneurons == True:
                    M = MyBdata
                    ListName = 'ALL'
                elif FortSNEUseAllGenesInsteadOfListDAneurons == False:
                    M = sd.get_submatrix(ListSelectedGenes, MyBdata, genes)
                
                if UsePCApreviousTotSNE == True:
    
                    # DO PCA BEFORE tSNE ADN USE RESULTS AS IMPUT FOR SUBSEQUENT tSNE!!!
                    y_pca_temp, y_tsne_temp = get_DR(M,g_ind=-1, d=NumberOfPCAbeforetSNE, off=[])
    
                    Myt1 = 'G2019S,'+str(DAY)+',norm - PCA\npoints colored by expression of '+GeneToColor+'_List '+ListName+'(tSNE PRECEDEED BY PCA)'
                    Myt2 = 'G2019S,'+str(DAY)+',norm - tSNE\npoints colored by expression of '+GeneToColor+'_List '+ListName+'(tSNE PRECEDEED BY PCA)'
                    MySaveAs1 = './Results/Jochen/PCA/Mutant/G2019S_day' + str(DAY) + '_PCA_DAneurons_'+GeneToColor+'_List ' + ListName + '_PCAbeforetSNE.pdf'
                    MySaveAs2 = './Results/Jochen/tSNE/Mutant/G2019S_day' + str(DAY) + '_tSNE_DAneurons_'+GeneToColor+'_List ' + ListName + '_PCAbeforetSNE.pdf'
                                                     
                    y_pca, y_tsne = get_DR(np.transpose(y_pca_temp),g_ind=-1, d=2, off=[],t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)
                    axes2[i,j].scatter(y_tsne[:,0],y_tsne[:,1])  
    
                elif UsePCApreviousTotSNE == False:
                    
                    t1 = 'GC,'+str(DAY)+', norm - PCA\npoints colored by expression of '+GeneToColor+'_List '+ListName
                    t2 = 'GC,'+str(DAY)+', norm - tSNE\npoints colored by expression of '+GeneToColor+'_List '+ListName
                    SaveAs1 = './Results/Jochen/PCA/Control/GC_day' + str(DAY) + '_PCA_DAneurons_'+GeneToColor+'_List ' + ListName + '.pdf'
                    SaveAs2 = './Results/Jochen/tSNE/Control/GC_day' + str(DAY) + '_tSNE_DAneurons_'+GeneToColor+'_List ' + ListName + '.pdf'
                                                                                  
                    y_pca, y_tsne = get_DR(M,genes.index(GeneToColor),2,[],t1,t2,SaveAs1,SaveAs2)
                    axes2[i,j].scatter(y_tsne[:,0],y_tsne[:,1])  
                
                if i == 0:
                    ListOftSNE_Control_DAneuronsJonasList.append(y_tsne)
                elif i == 1:
                    ListOftSNE_Control_CRCJochenList.append(y_tsne)
                    
            j=j+1
            
        plt.show()
      
            
    if GaussianMixtureModelOnIndividualtSNEwithCumGenExpr == 1:
    
        L_JonasDA = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_DopaminergicNeurons.txt')[0]
        L_JochenCRC = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJochen/list_core_genes_JO_DAneurons_CRC.txt')[0]
        
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.mlab as mlab
        import math
        import copy
        from mpl_toolkits.mplot3d import Axes3D
        from sklearn import mixture
        from scipy import linalg
        from matplotlib.patches import Ellipse
    
        def ClusterDataWithGaussianMixtureModel(samples, Ncomponents, MyTitle, SaveAs='DoNotSave'):
            #fit the gaussian mixture model to the data
            gmix = mixture.GaussianMixture(n_components=Ncomponents, covariance_type='full')
            gmix.fit(samples)
            #print(gmix.means_)
            #print(' ')
            #print(gmix.covariances_)
            
            # plot the data coloured differently for each cluster
            colors = ['r' if i==0 else ('g' if i==1 else ('b' if i==2 else 
                                                          ('y' if i==3 else 
                                                           ('c' if i==4 else 
                                                            ('m' if i==5 else 
                                                             ('k' if i==6 else
                                                              ('tab:orange' if i==7 else
                                                               ('tab:purple' if i==8 else
                                                                ('tab:brown' if i==9 else
                                                                 ('tab:pink' if i==10 else
                                                                  ('tab:gray'))))))))))) 
                                                                for i in gmix.predict(samples)]
            figGMM = plt.figure()
            plt.title(MyTitle)
            ax = plt.gca()
            ax.scatter(samples[:,0], samples[:,1], c=colors, alpha=0.8)
            
            # Make ellipses around each claster
            ListOfAngles = []
            ListOfVs = []
            for Covar in gmix.covariances_:
                v, w = linalg.eigh(Covar)
                v = 2. * np.sqrt(2.) * np.sqrt(v)
                u = w[0] / linalg.norm(w[0])
                angle = np.arctan(u[1] / u[0])
                angle = 180. * angle / np.pi  # convert to degrees
                ListOfAngles.append(angle)
                ListOfVs.append(v)
            ColorsList = ['r','g','b','y','c','m','k','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray']
            i=0
            for Center in gmix.means_:
                v=ListOfVs[i]
                ellipse = Ellipse(xy=Center, width=v[0], height=v[1], angle=180+ListOfAngles[i], edgecolor='black')
                ellipse.set_alpha(0.3)
                ellipse.set_facecolor(ColorsList[i])
                ax.add_artist(ellipse)
                i=i+1
                
            plt.show()
            if SaveAs != 'DoNotSave':
                plt.savefig(SaveAs)
            return gmix, figGMM
        
        
        for IndexDay in [0 ,1 ,2 ,3]:
            for IndexControlMutant in [0, 1]:
                for IndexGenesList in [0, 1]:
        
                    ThisDay = IndexDay # 0 = Day0, 1 = Day 10, 2 = Day 14, 3 = Day 42
                    CurrentListIndex = IndexControlMutant # 0 = Control, 1 = Mutant
                    CurrentListOfDAneuronsGenes = IndexGenesList # 0 = Jonas DA neurons, 1 = Jochen CRC genes
                    
    #    ThisDay = 2 # 0 = Day0, 1 = Day 10, 2 = Day 14, 3 = Day 42
    #    CurrentListIndex = 1 # 0 = Control, 1 = Mutant
    #    CurrentListOfDAneuronsGenes = 0 # 0 = Jonas DA neurons, 1 = Jochen CRC genes
                                    
                    if ThisDay == 0:
                        ThisDayString = "Day 0"
                    elif ThisDay == 1:
                        ThisDayString = "Day 10"
                    elif ThisDay == 2:
                        ThisDayString = "Day 14"
                    elif ThisDay == 3:
                        ThisDayString = "Day 42"
                        
                    ThisDayStringNoSpaces = ThisDayString.replace(" ", "")
                
                
                    if CurrentListIndex == 0:
                        ThisListString = "GC"
                        normABlist = [B_0_norm,B_10_norm,B_14_norm,B_42_norm]
                        genesABlist = [B_0_genes,B_10_genes,B_14_genes,B_42_genes]
                        if CurrentListOfDAneuronsGenes == 0:
                            L_UsedForCumulative = L_JonasDA
                            GenesListString = "Jonas DA neurons genes list"
                            CurrentListOftSNEs = ListOftSNE_Control_DAneuronsJonasList
                        elif CurrentListOfDAneuronsGenes == 1:
                            L_UsedForCumulative = L_JochenCRC
                            GenesListString = "Jochen CRC genes list"
                            CurrentListOftSNEs = ListOftSNE_Control_CRCJochenList
                    elif CurrentListIndex == 1:
                        ThisListString = "G2019S"
                        normABlist = [A_0_norm,A_10_norm,A_14_norm,A_42_norm]
                        genesABlist = [A_0_genes,A_10_genes,A_14_genes,A_42_genes]
                        if CurrentListOfDAneuronsGenes == 0:
                            L_UsedForCumulative = L_JonasDA
                            GenesListString = "Jonas DA neurons genes list"
                            CurrentListOftSNEs = ListOftSNE_Mutant_DAneuronsJonasList
                        elif CurrentListOfDAneuronsGenes == 1:
                            L_UsedForCumulative = L_JochenCRC
                            GenesListString = "Jochen CRC genes list"
                            CurrentListOftSNEs = ListOftSNE_Mutant_CRCJochenList
                        
                    GenesListStringNoSpaces = GenesListString.replace(" ", "")
                    
                    MyData = CurrentListOftSNEs[ThisDay]
                    
                    # Plot Row data
                #    fig1 = plt.figure()
                #    plt.title("Original data")
                #    plt.scatter(MyData[:,0],MyData[:,1])
                #    plt.show()
                #    
                #    ########## NOW DO THE GAUSSIAN MIXED MODEL ##########
                #    
                #    # Train a Gaussian Mixture Model with 2 components
                    MyTitle = "Unsupervised clustering via Gaussian Mixture Model\nof tSNE (" + GenesListString + ") of "  + ThisListString + " at " + ThisDayString
                    SaveAsGMM = './Results/Jochen/GaussianMixtureModels/IndividualPlotsClusters/GMMclustering_' + ThisListString + '_' + GenesListStringNoSpaces + '_' + ThisDayStringNoSpaces + '.pdf'
                    MyGMMmodel, FigureGMM = ClusterDataWithGaussianMixtureModel(MyData,NumberOfClustersChosen,MyTitle,SaveAsGMM)
                    
                    # Take out from the model the means of the gaussians and the flags telling you wich datapoints belong to wich cluster
                #    MeansOfGMMs = MyGMMmodel.means_
                #    print(MeansOfGMMs)
                    WhichClusterDataLinesBelongTo = MyGMMmodel.predict(MyData)
                    
                    ########## PLOT CUMULATIVE GENE EXPRESSION CLUSTER BY CLUSTER ##########    
                    
                    fig, axes = plt.subplots(nrows=1, ncols=NumberOfClustersChosen)
                    #plt.subplots(1,2, gridspec_kw = {'width_ratios':[3, 1]})
                    fig.subplots_adjust(hspace=0.3, wspace=0.1)
                    plt.suptitle('Cumulative gene expression (' + GenesListString + ') from clusters\nof tSNE of '  + ThisListString + " at " + ThisDayString)
                
                    for ax in axes.flat:
                        # Hide all ticks and labels
                        ax.xaxis.set_visible(True)
                        ax.yaxis.set_visible(True)
                    
                    j=0
                    ListOfCumulativeGeneExpressions = []
                    ListOfYmax = []
                    for Cluster in range(NumberOfClustersChosen):
                
                #        temp = (WhichClusterDataLinesBelongTo != Cluster)
                #        temp2 = np.array(range(len(temp)))
                #        IndexesOfNOTcurrentCluster = temp2[temp]
                #
                #        import copy
                #        CurrentClasterScaA = copy.deepcopy(scaAlist[ThisDay])
                #        for i in IndexesOfNOTcurrentCluster:
                #            CurrentClasterScaA[:,i]=np.zeros(len(CurrentClasterScaA))
                #
                        temp = (WhichClusterDataLinesBelongTo == Cluster)
                        temp2 = np.array(range(len(temp)))
                        IndexesOfcurrentCluster = temp2[temp]
                
                        CurrentClasterScaA = copy.deepcopy(normABlist[ThisDay])
                #        for i in IndexesOfNOTcurrentCluster:
                #            CurrentClasterScaA[:,i]=np.zeros(len(CurrentClasterScaA))
                
                        CumulativeGenesExpressionAB = np.zeros(len(IndexesOfcurrentCluster))
                
                        for GeneName in L_UsedForCumulative:
                            if GeneName in genesABlist[ThisDay]:
                                CumulativeGenesExpressionAB = CumulativeGenesExpressionAB + CurrentClasterScaA[genesABlist[ThisDay].index(GeneName),IndexesOfcurrentCluster]
                        
                #        NumberOfCellsInThisCluster = len(IndexesOfcurrentCluster)
                #        AverageGeneExpression = CumulativeGenesExpressionA / NumberOfCellsInThisCluster
                #       axes[0,0].set_title("Day 0")
                #       axes[0,1].set_title("Day 10")
                #       axes[0,2].set_title("Day 14")
                #       axes[0,3].set_title("Day 42")
                #
                #       axes[0,3].annotate("DA Neuron Genes\nJonas",xy=(1.05, 0.4),xycoords='axes fraction', size='large')
                #       axes[1,3].annotate("CRC Genes\nJochen", xy=(1.05, 0.4),xycoords='axes fraction', size='large')
                #       axes[2,3].annotate("Core Genes\nJochen",xy=(1.05, 0.4),xycoords='axes fraction', size='large')
                
                #       axes[0,3].annotate("DA neurons\nJonas",xy=(5, 0.20))
                #       axes[1,3].annotate("CRC\nJochen", xy=(150, 0.01))
                #       axes[2,3].annotate("core genes\nJochen",xy=(15, 0.40))
                
                        CurrentClusterColor = ['r' if Cluster==0 else ('g' if Cluster==1 else ('b' if Cluster==2 else 
                                                  ('y' if Cluster==3 else 
                                                   ('c' if Cluster==4 else 
                                                    ('m' if Cluster==5 else 
                                                     ('k' if Cluster==6 else
                                                      ('tab:orange' if Cluster==7 else
                                                       ('tab:purple' if Cluster==8 else
                                                        ('tab:brown' if Cluster==9 else
                                                         ('tab:pink' if Cluster==10 else
                                                          ('tab:gray')))))))))))]
                        ########## MAKE HISTOGRAMS ##########
                        BinsList = np.arange(math.ceil(max(CumulativeGenesExpressionAB))+1)
                        y, x, _ = axes[j].hist(CumulativeGenesExpressionAB, BinsList, histtype='stepfilled', color=CurrentClusterColor, normed=0, label='PAR2-G2019S', alpha=0.8)        
                        #axes[i,j].hist(CumulativeGenesExpressionB,bins=15, histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), normed=1, label='PAR2-GC', alpha=0.5)
                
                        ListOfYmax.append(y.max())
                
                #        plt.legend()
                        ListOfCumulativeGeneExpressions.append(CumulativeGenesExpressionAB)
                        j=j+1
                    
                    MaxY = max(ListOfYmax)
                    ListOfMax = []
                    for k in range(Cluster+1):
                        ListOfMax.append(max(ListOfCumulativeGeneExpressions[k]))
                    xmax = max(ListOfMax)
                    
                    for k in range(Cluster+1):
                        axes[k].set_xlim(0,math.ceil(xmax))  
                        axes[k].set_ylim(0,MaxY+1) # min(MaxY+1,70)
                        axes[k].set_xlabel("Cumulative\nGene Expression")
                        if k != 0:
                            axes[k].set_yticklabels([])
                    axes[0].set_ylabel("Number of Cells")
                
                    plt.show()
                    SaveAs = './Results/Jochen/GaussianMixtureModels/CumulativeGenesExpressions/CumulGeneExprGMM_' + ThisListString + '_' + GenesListStringNoSpaces + '_' + ThisDayStringNoSpaces + '.pdf'
                    plt.savefig(SaveAs)
                    
                    
#if GaussianMixtureModelOnALLtSNE == 1:
#        ### THIS DOES NOT WORK!!! ###
#
#    Nlines = 2
#    Ncols = 4
#    
#    fig3, axes3 = plt.subplots(nrows=Nlines, ncols=Ncols)
#    fig3.subplots_adjust(hspace=0.3, wspace=0.05)
#    plt.suptitle('GMM Clustering of tSNEs performed of GC (the control) XXX')
#
#    for ax in axes3.flat:
#        # Hide all ticks and labels
#        ax.xaxis.set_visible(False)
#        ax.yaxis.set_visible(False)
#        
##        if ax.is_first_col():
##            ax.yaxis.set_ticks_position('left')
##            ax.yaxis.set_visible(True)
##        if ax.is_last_row():
##           ax.xaxis.set_visible(True)
##           ax.set_xticklabels([])
##           ax.xaxis.set_ticks_position('none') 
#    
#    j = 0
#    for DAY in [0,10,14,42]:
#        print()
#        print('Day ' + str(DAY))
#        print()
#
#        for i in range(Nlines):         
#            print('List ' + str(i))
#            if i == 0:
#                ListSelectedGenes = L_JonasDA
#            elif i == 1:
#                ListSelectedGenes = L_JochenCRC
#
#            #axes2[i,0].set_ylabel("X label")
#            #axes2[1,j].set_xlabel("Y label")
#
#            axes3[0,0].set_title("Day 0")
#            axes3[0,1].set_title("Day 10")
#            axes3[0,2].set_title("Day 14")
#            axes3[0,3].set_title("Day 42")
#
#            axes3[0,3].annotate("DA Neuron Genes\nJonas",xy=(1.05, 0.4),xycoords='axes fraction', size='large')
#            axes3[1,3].annotate("CRC Genes\nJochen", xy=(1.05, 0.4),xycoords='axes fraction', size='large')
#            #axes2[2,3].annotate("Core Genes\nJochen",xy=(1.05, 0.4),xycoords='axes fraction', size='large')
#                        
#            ########## MAKE tSNE ##########
##            sca = scaBlist[j]
##            genes = genesBlist[j]
##            M = sd.get_submatrix(ListSelectedGenes, sca, genes)
##            t1 = 'GC,'+str(DAY)+',sca - PCA\npoints colored by expression of '+GeneToColor+'_List '+str(i)
##            t2 = 'GC,'+str(DAY)+',sca - tSNE\npoints colored by expression of '+GeneToColor+'_List '+str(i)
##            SaveAs1 = './Results/Jochen/PCA/Control/GC_day' + str(DAY) + '_PCA_DAneurons_'+GeneToColor+'_List ' + str(i) + '.pdf'
##            SaveAs2 = './Results/Jochen/tSNE/Control/GC_day' + str(DAY) + '_tSNE_DAneurons_'+GeneToColor+'_List ' + str(i) + '.pdf'
##                                                                          
##            y_pca, y_tsne = get_DR(M,genes.index(GeneToColor),2,[],t1,t2,SaveAs1,SaveAs2)
#
#            #MyTitle = "Unsupervised clustering via Gaussian Mixture Model\nof tSNE (" + GenesListString + ") of "  + ThisListString + " at " + ThisDayString
#            MyGMMmodel, FigureGMM = ClusterDataWithGaussianMixtureModel(MyData,NumberOfClustersChosen,"ciao")
#            FigureGMM.show()
#            axes3[i,j].scatter([1,2,3],[1,2,1])  
#            
#        j=j+1
#        
#    plt.show()
#            



###### KEEP THESE COMMENTS SINCE THEY ARE OBSERVATIONS!!! ######
# NORMALIZATIONS: norm gives rise to better tSNE than sca, log is as good as norm and does not change much
# NUMBER OF CELLS: The weird numbers used by Suresh are so far unjustified.
# Tomasz suggested everything to 250.
# In any case, the higher the cell number, the more cells clusters, 
# so comparing tSNE from days with different cell numbers may be misleading!!!

# PROBLEM: tSNE gives for mutant much better clusters separation at Day 0 than at Day 10, 
# against bio and previous Luis results!!!
# VERIFICATIONS:
# 1) NAIVE ERRORS IN CODE? 
# NO.
# 2) DOES CELL NUMBER METTERS? 
# YES. DAY 10 --> 250, DAY 0 --> 500, IF BOTH INCREASED TO 500 CELL CLUSTERING IS SIMILAR.
# 3) DOES NUMBER OF PCA COMPONENT USED AS IMPUT FOR tSNE MATTERS? 
# 4) DOES tSNE done with ALL genes (and then 40 PCAs as input) looks more meaningfull??? 
# YES. The configuration as Luis did, with ALL genes, 40 PCA inputs and cell numbers from Suresh, 
# gives a tSNE which separate clusters much better at days 10 and 14, and no separation at day 0, bad separation at day 42.


if PlotHistogramsERN1 == 1:
    
    Clist = [C0,C10,C14,C42]    
    normClist = [C0_norm,C10_norm,C14_norm,C42_norm]
    
#        fig1 = plt.figure()
#        plt.hist(C_norm[C.c_names.index('LRRK2'),0:250],bins=20, histtype='stepfilled', color='r', label='G2019S', alpha=0.5)        
#        plt.hist(C_norm[C.c_names.index('LRRK2'),250:500],bins=20, histtype='stepfilled', color='gray', label='GC', alpha=0.5)
#        plt.title("Differential expression of LRRK2")
#        plt.xlabel("Gene expression (norm)")
#        plt.ylabel("Number of Cells")
#        plt.legend()
#        plt.show()
#        plt.savefig('./Results/jonas_DiffExprLRRK2_day' + str(DAY) + '.pdf')
#        
    Nlines = 2
    Ncols = 4
    fig, axes = plt.subplots(nrows=Nlines, ncols=Ncols)
    fig.subplots_adjust(hspace=0.6, wspace=0.15)

    for ax in axes.flat:
        # Hide all ticks and labels
        #ax.xaxis.set_visible(False)
        #ax.yaxis.set_visible(False)

        # Set up ticks only on one side for the "edge" subplots...
        if ax.is_first_col():
            ax.yaxis.set_ticks_position('left')
            ax.yaxis.set_visible(True)

    j = 0
    for DAY in [0,10,14,42]:
        print('Day ' + str(DAY))
        print()
        C = Clist[j]
        C_norm = normClist[j]
        
        #plt.title("Differential expression of LRRK2")

        axes[0,0].set_ylabel("Number of cells")
        axes[1,0].set_ylabel("Number of cells")
        axes[0,j].set_xlabel("ERN1 expression (norm. to max)")
        axes[1,j].set_xlabel("ERN1 expression (norm. to max)")

        MaxX = max(max(C_norm[C.c_names.index('ERN1'),0:250]),max(C_norm[C.c_names.index('ERN1'),250:500]))
        GenesExpressionA_NORM = C_norm[C.c_names.index('ERN1'),0:250] / MaxX
        GenesExpressionB_NORM = C_norm[C.c_names.index('ERN1'),250:500] / MaxX
        binwidth = 0.035
        axes[0,j].hist(GenesExpressionA_NORM, bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(255.0/256,0.,0.), normed=0, label='PAR2-G2019S', alpha=0.5)        
        axes[0,j].hist(GenesExpressionB_NORM, bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), normed=0, label='PAR2-GC', alpha=0.5)

        axes[0,j].set_xlim(0,1.1)
        axes[0,j].set_ylim(0,250)

        axes[1,j].hist(GenesExpressionA_NORM, bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(255.0/256,0.,0.), normed=0, label='PAR2-G2019S', alpha=0.5)        
        axes[1,j].hist(GenesExpressionB_NORM, bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=(128.0/256,128.0/256,128.0/256), normed=0, label='PAR2-GC', alpha=0.5)

        axes[1,j].set_xlim(0,1.1)
        axes[1,j].set_ylim(0,30)

        ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
        # of 2 populations are statistically significantly different ######
        
        #DataPopulation1 = numpy.array([1,2,3,4,5,2,3,4,3])
        DataPopulation1 = C_norm[C.c_names.index('ERN1'),0:250]
        #DataPopulation2 = numpy.array([1,2,3.001,4,5,2,3,4,3])
        DataPopulation2 = C_norm[C.c_names.index('ERN1'),250:500]
        DataPopulation1_Object = wstats.DescrStatsW(DataPopulation1)
        DataPopulation2_Object = wstats.DescrStatsW(DataPopulation2)
        MyComparison = wstats.CompareMeans(DataPopulation1_Object,DataPopulation2_Object)
        TestStatistics, pvalue = MyComparison.ztest_ind(alternative='two-sided', usevar='unequal', value=0)

        # Remember Bonferroni correction for multiple testing, 
        # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
        # by the number of repetitions of the test, 4 days x 7 lists = 28
        # The 7 lists are Stem, DAneur, Mito, Cell Cycle, Pro-apoptosis, Anti-apoptosis, Caspases
        print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
        print(' ')

        plt.legend()

        j=j+1
        plt.show()
        plt.savefig('./Results/jonas_DiffExprERN1.pdf')