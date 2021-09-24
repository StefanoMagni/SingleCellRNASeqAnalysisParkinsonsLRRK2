#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Mon Dec 02 2019
# @authors: stefano.magni@uni.lu

from Script_data import *
from Script_information import *
from genes_1 import *
from SingleCellLibrariesStefano import *

from scipy.stats.stats import pearsonr
import Script_data as sd
import numpy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import matplotlib.ticker as plticker
import seaborn as sns
import copy
import matplotlib as mpl
import math
import matplotlib.mlab as mlab
from mpl_toolkits.mplot3d import Axes3D
from scipy import linalg
from matplotlib.patches import Ellipse
import matplotlib.patches as mpatches
import os

from sklearn import mixture
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


# Pre Processing
GetDataFromFiles = True # True                                                            # PRE-PROCESSING
ComputeKneePlots = True # True                                                            # PRE-PROCESSING
PlotKneePlots = True # True                                                               # PRE-PROCESSING
CutCellNumber = True # True # BOTH CutCellNumber AND FilterAndCutCellNumber CAN BE DONE   # PRE-PROCESSING
# RUN THE FOLLOWING ONLY IF YOU WISH TO ADDITIONALLY FILTER OUT CELLS ACCORDING TO QUALITY CRITERIA
FilterAndCutCellNumber = False # False                                                     # PRE-PROCESSING
Normalize = True # True                                                                   # PRE-PROCESSING
ComparePreliminaries = True # True                                                        # PRE-PROCESSING
# RUN THE FOLLOWING ONLY IF YOU WISH TO ADDITIONALLY FILTER OUT GENES ACCORDING TO QUALITY CRITERIA
FilterOutGenesExpressedInFewCells = False # False  # ALWAYS RUN TOGETHER WITH FilterAndCutCellNumber  # PRE-PROCESSING
SaveListsOfCommonGenes = True # True                                                      # PRE-PROCESSING

# SaveCutFilteredNormalizedData = 0 # NOT FULLY WORKING YET                     # SAVE DATA
# GetFromFileTransformedData = 0 # NOT FULLY WORKING YET                        # LOAD DATA 

ProduceListsOfDifferentiallyExpressedGenes = False # False                                  # DEGs

HistoAndFoldChangesGeneListsGeneralCellProcesses = True # True 
DoGeneGeneCorrelation_CellDeath_VS_CellTypesGenes = False # False 


ListOfDays = ['0','10','14','42']
N_Days = len(ListOfDays)
ColorMU = (255.0/256,0.,0.)
ColorCT = (128.0/256,128.0/256,128.0/256)

# 250 Cells is a good choice, I also considered alternatively 800 which seems to be about the max we could consider according to knee plot (plotted below)
N_cells = 250 # 800


MYDIR = ("Results")
CHECK_FOLDER = os.path.isdir(MYDIR)
# If folder doesn't exist, then create it.
if not CHECK_FOLDER:
    os.makedirs(MYDIR)
    print("created folder : ", MYDIR)
else:
    print(MYDIR, "folder already exists.")


if GetDataFromFiles == 1:
    A_genes, A_expr, A_expr_np = [None] * N_Days, [None] * N_Days, [None] * N_Days
    B_genes, B_expr, B_expr_np = [None] * N_Days, [None] * N_Days, [None] * N_Days
    for IndexDay in range(N_Days):
        Day = ListOfDays[IndexDay]
        print('Loading data for Day ' + Day + '...')
        A_genes[IndexDay], A_expr[IndexDay], A_expr_np[IndexDay] = sd.get_from_txt(
                './Data_SingleCellRNASeq_DigitalExpressionMatrices/A_' + str(Day) + '.txt')
        B_genes[IndexDay], B_expr[IndexDay], B_expr_np[IndexDay] = sd.get_from_txt(
                './Data_SingleCellRNASeq_DigitalExpressionMatrices/B_' + str(Day) + '.txt')
    print('Done.')


if ComputeKneePlots == 1:
    print('Computing Knee Plots...')
    aA, bA = [None] * N_Days, [None] * N_Days
    aB, bB = [None] * N_Days, [None] * N_Days
    for IndexDay in range(N_Days):
        aA[IndexDay],bA[IndexDay] = sd.get_knee_plot(A_expr_np[IndexDay])
        aB[IndexDay],bB[IndexDay] = sd.get_knee_plot(B_expr_np[IndexDay])
    print('Done.')


if PlotKneePlots == 1:
    
    print('Plotting Knee Plots...')
    fig = plt.figure(facecolor="white")
    plt.clf()    
    plt.xlabel('Cell Number')
    plt.ylabel('Cumulative fraction of transcripts')
    linestyles = ['-','--','-.',':']
    for IndexDay in range(N_Days):
        Day = ListOfDays[IndexDay]
        plt.plot(aA[IndexDay], linestyle=linestyles[IndexDay], c=ColorMU, label='Day ' + Day + ', G2019S')
        plt.plot(aB[IndexDay], linestyle=linestyles[IndexDay], c=ColorCT, label='Day ' + Day + ', GC')
    plt.legend()
    plt.title('Knee Plot, Cumulative')
    plt.show()
    plt.savefig('./Results/KneePlotCumulative_G2019S_GC.pdf')
    
    fig = plt.figure(facecolor="white")
    plt.clf() 
    plt.xlabel('Cell Number')  
    plt.ylabel('Number of transcripts')   
    for IndexDay in range(N_Days):
        Day = ListOfDays[IndexDay]
        plt.plot(bA[IndexDay], linestyle=linestyles[IndexDay], c=ColorMU, label='Day ' + Day + ', G2019S')
        plt.plot(bB[IndexDay], linestyle=linestyles[IndexDay], c=ColorCT, label='Day ' + Day + ', GC')
    plt.legend()
    plt.title('Knee Plot, Absolute')
    plt.show()
    plt.savefig('./Results/KneePlotAbsolute_G2019S_GC.pdf')
    print('Done.')


if CutCellNumber == 1:
    ListOfAs_cut = [None] * N_Days
    ListOfBs_cut = [None] * N_Days
    for IndexDay in range(N_Days):
        print('Cutting Cell Number for Day ' + ListOfDays[IndexDay] + '...')
        ListOfAs_cut[IndexDay] = sd.get_top_cells_matrix(A_expr_np[IndexDay],N_cells)[0]
        ListOfBs_cut[IndexDay] = sd.get_top_cells_matrix(B_expr_np[IndexDay],N_cells)[0]
    print('Done.')
    ListOfAs_Filtered = ListOfAs_cut    
    ListOfBs_Filtered = ListOfBs_cut
    

if FilterAndCutCellNumber == 1:
    
    ListOfAs_Filtered = [None] * N_Days
    ListOfBs_Filtered = [None] * N_Days
    List_N_cells_A_AfterRNAfilters = [None] * N_Days
    for IndexDay in range(N_Days):
        print('Filtering out cells for Day ' + ListOfDays[IndexDay] + '...')
        print('')
        Day = ListOfDays[IndexDay]
        
        ##### CUT 2: Cut Cells Expressing too much the MITO Genes, in percantage w.r.t. Total Expression #####
        
        TotMITOReadCountsAsPercentage_A, TotMITOReadCounts_A, TotReadCounts_A = ComputeTotalMITOCHONDRIALGenesExpression(ListOfAs_CutAndOrFiltered[IndexDay], A_genes[IndexDay])
        TotMITOReadCountsAsPercentage_B, TotMITOReadCounts_B, TotReadCounts_B = ComputeTotalMITOCHONDRIALGenesExpression(ListOfBs_CutAndOrFiltered[IndexDay], B_genes[IndexDay])
    
        fig = plt.figure(facecolor="white")
        # Remove nans for plotting only
        Aux_A = TotMITOReadCountsAsPercentage_A[numpy.logical_not(numpy.isnan(TotMITOReadCountsAsPercentage_A))]
        Aux_B = TotMITOReadCountsAsPercentage_B[numpy.logical_not(numpy.isnan(TotMITOReadCountsAsPercentage_B))]
        plt.hist(Aux_A,bins=np.arange(0,1,0.05), histtype='stepfilled', color=ColorMU, normed=0, alpha=0.5)
        plt.hist(Aux_B,bins=np.arange(0,1,0.05), histtype='stepfilled', color=ColorCT, normed=0, alpha=0.5)
        plt.xlabel('Total MITOCHONDRIAL Read Counts per Cell')
        plt.ylabel('Cells')
        plt.title('Distributions of cells across Total MITOCHONDRIAL Read Counts, Day ' + Day)
        plt.savefig('./Results/DistributionsTotalMITOCHONDRIALReadCountsDay' + Day +'.pdf')
        
        fig = plt.figure(facecolor="white")
        MaxXLim = max(max(TotReadCounts_A),max(TotReadCounts_B))
        plt.hist(TotReadCounts_A,bins=np.arange(0,MaxXLim,50), histtype='stepfilled', color=ColorMU, normed=0, alpha=0.5)
        plt.hist(TotReadCounts_B,bins=np.arange(0,MaxXLim,50), histtype='stepfilled', color=ColorCT, normed=0, alpha=0.5)
        plt.xlabel('Total Read Counts per Cell')
        plt.ylabel('Cells')
        plt.title('Distributions of cells across Total Read Counts, Day ' + Day)
        plt.savefig('./Results/DistributionsTotalReadCounts_1_KneePlotCutDay' + Day +'.pdf')
        
        fig = plt.figure(facecolor="white")
        plt.scatter(TotReadCounts_A, TotMITOReadCountsAsPercentage_A, color=ColorMU, alpha=0.5)
        plt.scatter(TotReadCounts_B, TotMITOReadCountsAsPercentage_B, color=ColorCT, alpha=0.5)
        plt.xlabel('Total Read Counts')
        plt.ylabel('Total MITOCHONDRIAL Read Counts')
        plt.title('Total MITOCHONDRIAL Read Counts VS Total Read Counts, Day ' + Day)
        plt.savefig('./Results/DistributionsTotalMITOCHONDRIALReadCountsVSTotalReadCountsDay' + Day +'.pdf')
        
        MeanMITOtotexpr = np.mean(np.concatenate((Aux_A,Aux_B)))
        StdMITOtotexpr = np.std(np.concatenate((Aux_A,Aux_B)))
        ThrasholdCutCellsMITOReadCounts = MeanMITOtotexpr + 2 * StdMITOtotexpr
    
        print('WE START WITH ' + str(np.shape(ListOfAs_cut[IndexDay])[1]) + ' CELLS')
        print("The mean Mitocondrial Expression is:")
        print(MeanMITOtotexpr)
        print("The std Mitocondrial Expression is:")
        print(StdMITOtotexpr)
        print("Thus the thrashold for Mitocondrial Expression is:")
        print(ThrasholdCutCellsMITOReadCounts)
                
        # Keep only cells for which the total expression of Mito genes is below Thrashold
        IndexesOfAbelowMITOThrashold = np.where(np.array(TotMITOReadCountsAsPercentage_A) <= ThrasholdCutCellsMITOReadCounts)[0]
        IndexesOfBbelowMITOThrashold = np.where(np.array(TotMITOReadCountsAsPercentage_B) <= ThrasholdCutCellsMITOReadCounts)[0]
        A_cut = ListOfAs_cut[IndexDay]
        B_cut = ListOfBs_cut[IndexDay]
        A_tempI = A_cut[:,IndexesOfAbelowMITOThrashold]
        B_tempI = B_cut[:,IndexesOfBbelowMITOThrashold]
    
        print('AFTER MITO CUT WE HAVE FOR A ' + str(np.shape(A_tempI)[1]) 
        + ' CELLS, AND FOR B ' +  str(np.shape(B_tempI)[1])  + ' CELLS')
        print('')
        
        # To ensure we performed the filtering we intended to, let's plot the cells left after filtering
        fig = plt.figure(facecolor="white")
        TotMITOReadCountsAsPercentage_A_AFTERFILTERING, TotMITOReadCounts_A_AFTERFILTERING, TotReadCounts_A_AFTERFILTERING = ComputeTotalMITOCHONDRIALGenesExpression(
                A_tempI, A_genes[IndexDay])
        TotMITOReadCountsAsPercentage_B_AFTERFILTERING, TotMITOReadCounts_B_AFTERFILTERING, TotReadCounts_B_AFTERFILTERING = ComputeTotalMITOCHONDRIALGenesExpression(
                B_tempI, B_genes[IndexDay])
        plt.scatter(TotReadCounts_A_AFTERFILTERING, TotMITOReadCountsAsPercentage_A_AFTERFILTERING, color=ColorMU, alpha=0.5)
        plt.scatter(TotReadCounts_B_AFTERFILTERING, TotMITOReadCountsAsPercentage_B_AFTERFILTERING, color=ColorCT, alpha=0.5)
        plt.xlabel('Total Read Counts AFTER MITO FILTERING')
        plt.ylabel('Total MITOCHONDRIAL Read Counts AFTER MITO FILTERING')
        plt.title('Total MITOCHONDRIAL Read Counts VS Total Read Counts AFTER MITO FILTERING, Day ' + Day)
        plt.savefig('./Results/DistributionsTotalMITOCHONDRIALReadCountsVSTotalReadCounts_AFTER_MITO_FILTERINGDay' + Day +'.pdf')
    
        fig = plt.figure(facecolor="white")
        MaxXLim = max(max(TotReadCounts_A_AFTERFILTERING),max(TotReadCounts_B_AFTERFILTERING))
        plt.hist(TotReadCounts_A_AFTERFILTERING,bins=np.arange(0,MaxXLim,50), histtype='stepfilled', color=ColorMU, normed=0, alpha=0.5)
        plt.hist(TotReadCounts_B_AFTERFILTERING,bins=np.arange(0,MaxXLim,50), histtype='stepfilled', color=ColorCT, normed=0, alpha=0.5)
        plt.xlabel('Total Read Counts per Cell')
        plt.ylabel('Cells')
        plt.title('Distributions of cells across Total Read Count, Day ' + Day)
        plt.savefig('./Results/DistributionsTotalReadCounts_2_AfterMITOFilteringDay' + Day + '.pdf')
        
        ##### CUT 3: Cut Cells Expressing too much the RIBO Genes, in percantage w.r.t. Total Expression #####
    
        TotRIBOReadCountsAsPercentage_A, TotRIBOReadCounts_A, TotReadCounts_A = ComputeTotalRIBOSOMALGenesExpression(A_tempI, A_genes[IndexDay])
        TotRIBOReadCountsAsPercentage_B, TotRIBOReadCounts_B, TotReadCounts_B = ComputeTotalRIBOSOMALGenesExpression(B_tempI, B_genes[IndexDay])
    
        fig = plt.figure(facecolor="white")
        # Remove nans for plotting only
        Aux_A = TotRIBOReadCountsAsPercentage_A[numpy.logical_not(numpy.isnan(TotRIBOReadCountsAsPercentage_A))]
        Aux_B = TotRIBOReadCountsAsPercentage_B[numpy.logical_not(numpy.isnan(TotRIBOReadCountsAsPercentage_B))]
        plt.hist(Aux_A,bins=np.arange(0,1,0.05), histtype='stepfilled', color=ColorMU, normed=0, alpha=0.5)
        plt.hist(Aux_B,bins=np.arange(0,1,0.05), histtype='stepfilled', color=ColorCT, normed=0, alpha=0.5)
        plt.xlabel('Total RIBOSOMAL Read Counts per Cell')
        plt.ylabel('Cells')
        plt.title('Distributions of cells across Total RIBOSOMAL Read Counts, Day ' + Day)
        plt.savefig('./Results/DistributionsTotalRIBOSOMALReadCountsDay' + Day +'.pdf')
        
        fig = plt.figure(facecolor="white")
        plt.scatter(TotReadCounts_A, TotRIBOReadCountsAsPercentage_A, color=ColorMU, alpha=0.5)
        plt.scatter(TotReadCounts_B, TotRIBOReadCountsAsPercentage_B, color=ColorCT, alpha=0.5)
        plt.xlabel('Total Read Counts')
        plt.ylabel('Total RIBOSOMAL Read Counts')
        plt.title('Total RIBOSOMAL Read Counts VS Total Read Counts, Day ' + Day)
        plt.savefig('./Results/DistributionsTotalRIBOSOMALReadCountsVSTotalReadCountsDay' + Day +'.pdf')
        
        fig = plt.figure(facecolor="white")
        plt.scatter(TotMITOReadCountsAsPercentage_A_AFTERFILTERING, TotRIBOReadCountsAsPercentage_A, color=ColorMU, alpha=0.5)
        plt.scatter(TotMITOReadCountsAsPercentage_B_AFTERFILTERING, TotRIBOReadCountsAsPercentage_B, color=ColorCT, alpha=0.5)
        plt.xlabel('Total MITOCONDRIAL Read Counts')
        plt.ylabel('Total RIBOSOMAL Read Counts')
        plt.title('Total RIBOSOMAL Reads VS Total MITOCONDRIAL Reads, Day ' + Day)
        plt.savefig('./Results/DistributionsTotalRIBOSOMALReadCountsVSTotalMITOCONDRIALReadCountsDay' + Day +'.pdf')
        
        MeanRIBOtotexpr = np.mean(np.concatenate((Aux_A,Aux_B)))
        StdRIBOtotexpr = np.std(np.concatenate((Aux_A,Aux_B)))
        ThrasholdCutCellsRIBOReadCounts = MeanRIBOtotexpr + 2 * StdRIBOtotexpr # 3 * StdRIBOtotexpr
    
        print("The mean Ribosomal Expression is:")
        print(MeanRIBOtotexpr)
        print("The std Ribosomal Expression is:")
        print(StdRIBOtotexpr)
        print("Thus the thrashold for Ribosomal Expression is:")
        print(ThrasholdCutCellsRIBOReadCounts)
                
        # Keep only cells for which the total expression of Mito genes is below Thrashold
        IndexesOfAbelowRIBOThrashold = np.where(np.array(TotRIBOReadCountsAsPercentage_A) <= ThrasholdCutCellsRIBOReadCounts)[0]
        IndexesOfBbelowRIBOThrashold = np.where(np.array(TotRIBOReadCountsAsPercentage_B) <= ThrasholdCutCellsRIBOReadCounts)[0]
        A_tempII = A_tempI[:,IndexesOfAbelowRIBOThrashold]
        B_tempII = B_tempI[:,IndexesOfBbelowRIBOThrashold]
        
        print('AFTER RIBO CUT WE HAVE FOR A ' + str(np.shape(A_tempII)[1]) 
        + ' CELLS, AND FOR B ' +  str(np.shape(B_tempII)[1])  + ' CELLS')
        print('')
        
        # To ensure we performed the filtering we intended to, let's plot the cells left after filtering
        fig = plt.figure(facecolor="white")
        TotRIBOReadCountsAsPercentage_A_AFTERFILTERING, TotRIBOReadCounts_A_AFTERFILTERING, TotReadCounts_A_AFTERFILTERING = ComputeTotalRIBOSOMALGenesExpression(A_tempII, A_genes[IndexDay])
        TotRIBOReadCountsAsPercentage_B_AFTERFILTERING, TotRIBOReadCounts_B_AFTERFILTERING, TotReadCounts_B_AFTERFILTERING = ComputeTotalRIBOSOMALGenesExpression(B_tempII, B_genes[IndexDay])
        plt.scatter(TotReadCounts_A_AFTERFILTERING, TotRIBOReadCountsAsPercentage_A_AFTERFILTERING, color=ColorMU, alpha=0.5)
        plt.scatter(TotReadCounts_B_AFTERFILTERING, TotRIBOReadCountsAsPercentage_B_AFTERFILTERING, color=ColorCT, alpha=0.5)
        plt.xlabel('Total Read Counts AFTER RIBO FILTERING')
        plt.ylabel('Total RIBOSOMAL Read Counts AFTER RIBO FILTERING')
        plt.title('Total RIBOSOMAL Read Counts VS Total Read Counts AFTER RIBO FILTERING, Day ' + Day)
        plt.savefig('./Results/DistributionsTotalRIBOSOMALReadCountsVSTotalReadCounts_AFTER_RIBO_FILTERINGDay' + Day +'.pdf')
        
        fig = plt.figure(facecolor="white")
        MaxXLim = max(max(TotReadCounts_A_AFTERFILTERING),max(TotReadCounts_B_AFTERFILTERING))
        plt.hist(TotReadCounts_A_AFTERFILTERING,bins=np.arange(0,MaxXLim,50), histtype='stepfilled', color=ColorMU, normed=0, alpha=0.5)
        plt.hist(TotReadCounts_B_AFTERFILTERING,bins=np.arange(0,MaxXLim,50), histtype='stepfilled', color=ColorCT, normed=0, alpha=0.5)
        plt.xlabel('Total Read Counts per Cell')
        plt.ylabel('Cells')
        plt.title('Distributions of cells across Total Read Counts, Day ' + Day)
        plt.savefig('./Results/DistributionsTotalReadCounts_3_AfterRIBOFilteringDay' + Day +'.pdf')
    
        ##### CUT 4: Cut further Cells Expressing too much or too few of Total Transcripts #####
    
        MeanTotExpr = np.mean(np.concatenate((TotReadCounts_A_AFTERFILTERING,TotReadCounts_B_AFTERFILTERING)))
        StdTotExpr = np.std(np.concatenate((TotReadCounts_A_AFTERFILTERING,TotReadCounts_B_AFTERFILTERING)))
        ThrasholdCutCellsHigherTotalReadCounts = MeanTotExpr + 2 * StdTotExpr
        ThrasholdCutCellsLowerTotalReadCounts = MeanTotExpr - 2 * StdTotExpr
    
        print("The mean Total Expression is:")
        print(MeanTotExpr)
        print("The std Total Expression is:")
        print(StdTotExpr)
        print("Thus the upper thrashold for Total Expression is:")
        print(ThrasholdCutCellsHigherTotalReadCounts)
        print("Thus the lower thrashold for Total Expression is:")
        print(ThrasholdCutCellsLowerTotalReadCounts)
    
        IndexesOfALowerThanThrasholdsMax = np.where(np.array(TotReadCounts_A_AFTERFILTERING) <= ThrasholdCutCellsHigherTotalReadCounts)[0]
        IndexesOfAHigherThanThrasholdsMix = np.where(np.array(TotReadCounts_A_AFTERFILTERING) >= ThrasholdCutCellsLowerTotalReadCounts)[0]
        IndexesOfAwithinThrasholds = np.intersect1d(IndexesOfALowerThanThrasholdsMax, IndexesOfAHigherThanThrasholdsMix)
    
        IndexesOfBLowerThanThrasholdsMax = np.where(np.array(TotReadCounts_B_AFTERFILTERING) <= ThrasholdCutCellsHigherTotalReadCounts)[0]
        IndexesOfBHigherThanThrasholdsMix = np.where(np.array(TotReadCounts_B_AFTERFILTERING) >= ThrasholdCutCellsLowerTotalReadCounts)[0]
        IndexesOfBwithinThrasholds = np.intersect1d(IndexesOfBLowerThanThrasholdsMax, IndexesOfBHigherThanThrasholdsMix)
    
        A_tempIII = A_tempII[:,IndexesOfAwithinThrasholds]
        B_tempIII = B_tempII[:,IndexesOfBwithinThrasholds]
        print('AFTER TOTAL mRNA CUT WE HAVE FOR A ' + str(np.shape(A_tempIII)[1]) 
        + ' CELLS, AND FOR B ' +  str(np.shape(B_tempIII)[1])  + ' CELLS')
        print('')
        
        # To ensure we performed the filtering we intended to, let's plot the cells left after filtering
        fig = plt.figure(facecolor="white")
        TotRIBOReadCountsAsPercentage_A_AFTERFILTERING, TotRIBOReadCounts_A_AFTERFILTERING, TotReadCounts_A_AFTERFILTERING = ComputeTotalRIBOSOMALGenesExpression(A_tempIII, A_genes[IndexDay])
        TotRIBOReadCountsAsPercentage_B_AFTERFILTERING, TotRIBOReadCounts_B_AFTERFILTERING, TotReadCounts_B_AFTERFILTERING = ComputeTotalRIBOSOMALGenesExpression(B_tempIII, B_genes[IndexDay])
        plt.scatter(TotReadCounts_A_AFTERFILTERING, TotRIBOReadCountsAsPercentage_A_AFTERFILTERING, color=ColorMU, alpha=0.5)
        plt.scatter(TotReadCounts_B_AFTERFILTERING, TotRIBOReadCountsAsPercentage_B_AFTERFILTERING, color=ColorCT, alpha=0.5)
        plt.xlabel('Total Read Counts AFTER TOTAL READS FILTERING')
        plt.ylabel('Total RIBOSOMAL Read Counts AFTER TOTAL READS FILTERING')
        plt.title('Total RIBOSOMAL Read Counts VS Total Read Counts AFTER TOTAL READS FILTERING, Day ' + Day)
        plt.savefig('./Results/DistributionsTotalRIBOSOMALReadCountsVSTotalReadCounts_AFTER_TOTAL_READS_FILTERINGDay' + Day +'.pdf')
    
        fig = plt.figure(facecolor="white")
        MaxXLim = max(max(TotReadCounts_A_AFTERFILTERING),max(TotReadCounts_B_AFTERFILTERING))
        plt.hist(TotReadCounts_A_AFTERFILTERING,bins=np.arange(0,MaxXLim,50), histtype='stepfilled', color=ColorMU, normed=0, alpha=0.5)
        plt.hist(TotReadCounts_B_AFTERFILTERING,bins=np.arange(0,MaxXLim,50), histtype='stepfilled', color=ColorCT, normed=0, alpha=0.5)
        plt.xlabel('Total Read Counts per Cell')
        plt.ylabel('Cells')
        plt.title('Distributions of cells across Total Read Counts, Day ' + Day)
        plt.savefig('./Results/DistributionsTotalReadCounts_4_AfterTotReadsFilteringDay' + Day +'.pdf')
        
        A_Filtered = A_tempIII
        B_Filtered = B_tempIII
        
        List_N_cells_A_AfterRNAfilters[IndexDay] = len(IndexesOfAwithinThrasholds)

    ListOfAs_Filtered = A_Filtered
    ListOfBs_Filtered = B_Filtered
    
    print('Done.')
         
    
if Normalize == 1:
    List_A_norm, List_A_log, List_A_sca = [None] * N_Days, [None] * N_Days, [None] * N_Days
    List_B_norm, List_B_log, List_B_sca = [None] * N_Days, [None] * N_Days, [None] * N_Days
    for IndexDay in range(N_Days):
        Day = ListOfDays[IndexDay]
        print('Normalizing Cell Matrices for Day ' + Day + '...')
        List_A_norm[IndexDay], List_A_log[IndexDay], List_A_sca[IndexDay] = sd.get_normalization(ListOfAs_Filtered[IndexDay])
        List_B_norm[IndexDay], List_B_log[IndexDay], List_B_sca[IndexDay] = sd.get_normalization(ListOfBs_Filtered[IndexDay])
    print('Done.')
    
    
if ComparePreliminaries == 1:
    listOfCobjects = [None] * N_Days
    C_list = [None] * N_Days
    C_norm_list, C_log_list, C_sca_list = [None] * N_Days, [None] * N_Days, [None] * N_Days
    for IndexDay in range(N_Days):
        print('Do Compare preliminaries for Day ' + ListOfDays[IndexDay])
        C_list[IndexDay] = Compare((List_A_norm[IndexDay],A_genes[IndexDay]),(List_B_norm[IndexDay],B_genes[IndexDay]))
        aux_temp = C_list[IndexDay].merge()[0]
        C_norm_list[IndexDay], C_log_list[IndexDay], C_sca_list[IndexDay] = sd.get_normalization(aux_temp)
    print('Done.')
    
    
if FilterOutGenesExpressedInFewCells == 1:
    
    ThrasholdMinFractionOfCellsAGeneMustBeExprIn = 0.01 # Alternative choice: 0.05 

    for IndexDay in range(N_Days):
        Day = ListOfDays[IndexDay]
        print('')
        print('Filter out genes not expressed in at least ' + str(
                ThrasholdMinFractionOfCellsAGeneMustBeExprIn) + ' fraction of cells, for Day ' + Day)
        ListOfNnonZeroCells = []
        ListOfNnonZeroCells_A = []
        ListOfNnonZeroCells_B = []
        ListOfGenesPassingTheTest = []
        for line in C_norm_list[IndexDay]:
            # First Compute Number of Cells Expressing each gene, for MU and CT separately
            N_NonZeroCells_A = np.count_nonzero(line[0:len(ListOfAs_Filtered[IndexDay][0])])
            N_NonZeroCells_B = np.count_nonzero(line[len(ListOfAs_Filtered[IndexDay][0]):len(line)])
            ListOfNnonZeroCells_A.append(N_NonZeroCells_A)
            ListOfNnonZeroCells_B.append(N_NonZeroCells_B)
            # Second, Compute Number of Cells Expressing each gene, for MU and CT separately
            N_NonZeroCells = np.count_nonzero(line)
            ListOfNnonZeroCells.append(N_NonZeroCells)
            # Third, verify which genes are expressed in at least the desired fraction of cells
            if N_NonZeroCells >= ThrasholdMinFractionOfCellsAGeneMustBeExprIn*len(line):
                Tag = True
            else:
                Tag = False
            ListOfGenesPassingTheTest.append(Tag)
    
        fig = plt.figure()
        plt.hist2d(ListOfNnonZeroCells_A, ListOfNnonZeroCells_B, cmap=plt.cm.jet, bins=(len(ListOfAs_Filtered[IndexDay][0]), len(ListOfBs_Filtered[IndexDay][0]))) # ZOOM AFTERWARD
        plt.xlabel('Number of Cells Expressing a Gene in MU')
        plt.ylabel('Number of Cells Expressing a Gene in CT')
        plt.title('Colorscale: Number of genes being expressed in that number of cells')
        plt.show()
    
        fig = plt.figure(facecolor="white")
        plt.hist(ListOfNnonZeroCells, bins=np.arange(0,len(C_norm_list[IndexDay][0]),1), histtype='stepfilled', color='r', normed=0, alpha=0.5)
        plt.xlabel('Number of Cells Expressing a Gene')
        plt.ylabel('Genes')
        plt.title('Distributions of Genes Across Number of Cells Expressing them,\n for MU and its CT together')
    
        # Actually perform the filtering of genes wich do not fulfill the quality criterion above
        NewCdataList = []
        NewCGeneNames = []
        C_common_genes = C_list[IndexDay].comm_genes()
        index = 0
        for TestPassedTag in ListOfGenesPassingTheTest:
            if TestPassedTag == True:
                NewCdataList.append(C_norm_list[IndexDay][index,:])
                NewCGeneNames.append(C_common_genes[index])
            else:
                print('Gene ' + str(C_common_genes[index]) + ' did not pass the filtering, as present in only ' + str(ListOfNnonZeroCells[index]) + ' cells') 
            index = index + 1
        NewCdataArray = np.array(NewCdataList)
        C_list[IndexDay] = Compare((NewCdataArray[:,0:len(ListOfAs_Filtered[IndexDay][0])],NewCGeneNames),
                    (NewCdataArray[:,len(ListOfAs_Filtered[IndexDay][0]):len(ListOfAs_Filtered[IndexDay][0])+len(ListOfBs_Filtered[IndexDay][0])],NewCGeneNames))
#        C_list[IndexDay] = Compare((NewCdataArray[:,0:len(ListOfAs_Filtered[IndexDay][0])],NewCGeneNames),
#                    (NewCdataArray[:,len(A_Filtered[0]):len(ListOfAs_Filtered[IndexDay][0])+len(ListOfBs_Filtered[IndexDay][0])],NewCGeneNames))
        aux = C_list[IndexDay].merge()[0]
        C_norm_list[IndexDay], C_log_list[IndexDay], C_sca_list[IndexDay] = sd.get_normalization(aux)
    print('Done.')
    
    
if SaveListsOfCommonGenes == 1:
    C_common_genes = [None] * N_Days
    for IndexDay in range(N_Days):
        Day = ListOfDays[IndexDay]
        C_common_genes[IndexDay] = C_list[IndexDay].comm_genes()    

        ListOfTrueFalseRPGenes = []
        for GeneName in C_common_genes[IndexDay]:
            if GeneName.startswith('RP'):
                ListOfTrueFalseRPGenes.append(True)
            else:
                ListOfTrueFalseRPGenes.append(False) 
        print('The total number of genes in the common dataset, for day ' + Day + ' , after cells and genes filtering (depending on the options) is: ' + str(len(C_common_genes[IndexDay])) )
        print('The number of genes in the common dataset which are ribosomal (i.e. starts with RP) is: ' + str(sum(ListOfTrueFalseRPGenes)) )
        print('So the fraction of the two above is: ' + str(sum(ListOfTrueFalseRPGenes)/len(C_common_genes[IndexDay])) )
        sd.save_array(np.transpose(np.array([np.transpose(C_common_genes[IndexDay]),np.transpose(C_common_genes[IndexDay])])), './Results/NamesOfCommonGenesBetweenControlAndMutand_AfterGenesFiltering_Day'+Day+'.csv')
        
    
#if SaveCutFilteredNormalizedData == 1:
#    print('Saving Data to Files for Future Use...')
#    AcellsNumber = [None] * N_Days
#    BcellsNumber = [None] * N_Days
#    for IndexDay in range(N_Days):
#        Day = ListOfDays[IndexDay]
#        print('Saving Data to Files for Future Use, Day ' + Day + '...')
#        genes = C_list[IndexDay].comm_genes() 
#        norm = C_norm_list[IndexDay]
#        log = C_log_list[IndexDay]
#        sca = C_sca_list[IndexDay]
#        #You can use it to save transformed data
#        sd.save_gene_list(genes,'./Results/G2019S_GC_day' + str(Day) + '_genes.txt')
#        sd.save_array(norm,'./Results/G2019S_GC_day' + str(Day) + '_norm.csv')
#        sd.save_array(log,'./Results/G2019S_GC_day' + str(Day) + '_log.csv')
#        sd.save_array(sca,'./Results/G2019S_GC_day' + str(Day) + '_sca.csv')
#           
#        AcellsNumber[IndexDay] = np.shape(List_A_norm[IndexDay][1])[0]
#        BcellsNumber[IndexDay] = np.shape(List_B_norm[IndexDay][1])[0]
#    sd.save_array([AcellsNumber,BcellsNumber],'./Results/G2019S_GC_CellsNumbers.csv')
#    print('Done.')
#
#    
#if GetFromFileTransformedData == 1:
#
#    D_genes_list = [None] * N_Days
#    D_norm_list, D_log_list, D_sca_list = [None] * N_Days, [None] * N_Days, [None] * N_Days
#    for IndexDay in range(N_Days):
#        Day = ListOfDays[IndexDay]
#        print('Loading Data from Files for Future Use, Day ' + Day + '...')
#        D_genes_list[IndexDay] = sd.get_from_txt('./Results/G2019S_GC_day' + str(Day) + '_genes.txt')[0]
#        D_norm_list[IndexDay] = sd.get_from_csv('./Results/G2019S_GC_day' + str(Day) + '_norm.csv',';')
#        D_log_list[IndexDay] = sd.get_from_csv('./Results/G2019S_GC_day' + str(Day) + '_log.csv',';')
#        D_sca_list[IndexDay] = sd.get_from_csv('./Results/G2019S_GC_day' + str(Day) + '_sca.csv',';')
#    AAA = sd.get_from_csv('./Results/G2019S_GC_CellsNumbers.csv',';')
#    print('Done.')
    
    
if ProduceListsOfDifferentiallyExpressedGenes == 1:

    for IndexDay in range(N_Days):
        Day = ListOfDays[IndexDay]
        
        C = C_list[IndexDay]
                
        # Saves list of genes that are different between A and B on a given day.
        aux = C.avg() 
        a = np.array(aux)
        a = np.transpose(a)
    
        BONFERRONI_DEG = len(aux[0]) 
        
        L_min_pval005, L_min_stats_pval005 = C.avg_list(0.05/BONFERRONI_DEG, aux[3], a)
        sd.save_array(L_min_stats_pval005, './DiffExprGenes_pval_0p05_PvalFromBestOf3Tests_' + Day + '.csv')
        # aux[2] means to use the p-value coming from the MI test only
        L_min_pval005, L_min_stats_pval005 = C.avg_list(0.05/BONFERRONI_DEG, aux[2], a)
        sd.save_array(L_min_stats_pval005, './DiffExprGenes_pval_0p05_from_MI_' + Day + '.csv')
    
        ListOfTrueFalseRPGenesInDEGs = []
        for GeneName in L_min_pval005:
            if GeneName.startswith('RP'):
                ListOfTrueFalseRPGenesInDEGs.append(True)
            else:
                ListOfTrueFalseRPGenesInDEGs.append(False) 
                
        print('The total number of DEGs at Day ' + Day + ' is: ' + str(len(L_min_pval005)) )
        print('The number of genes in the DEGs list which are ribosomal (i.e. starts with RP) is: ' + str(sum(ListOfTrueFalseRPGenesInDEGs)) )
        print('So the fraction of the two above is: ' + str(sum(ListOfTrueFalseRPGenesInDEGs)/len(L_min_pval005)) )
    
        # Now order DEGs by decreasing p-values
    
        OrderedIndexes = np.argsort(aux[3])
        PvaluesArray = np.array(aux[3])
        OrderedPvalues = PvaluesArray[OrderedIndexes]
        GeneNamesArray = np.array(C.c_names)
        OrderedGeneNames = GeneNamesArray[OrderedIndexes]
        
        GenesExpMeansA = np.array(aux[5])
        GenesExpMeansB = np.array(aux[6])
        OrderedMeansA = GenesExpMeansA[OrderedIndexes]
        OrderedMeansB = GenesExpMeansB[OrderedIndexes]
        
        N_Genes = len(L_min_stats_pval005)
        OrderedGeneNamesCUT = OrderedGeneNames[0:N_Genes]
        OrderedPvaluesCUT = OrderedPvalues[0:N_Genes]
        OrderedMeansACUT = OrderedMeansA[0:N_Genes]
        OrderedMeansBCUT = OrderedMeansB[0:N_Genes]
        OrderedRATIOofMeansAoverB = OrderedMeansACUT/OrderedMeansBCUT
        
        OrderedNamesPvalsRatios = np.transpose(np.array([OrderedGeneNamesCUT, OrderedPvaluesCUT, OrderedRATIOofMeansAoverB]))
    
        # Now order DEGs by decreasing fold change MU / CT
        
        ArrayData = np.array(OrderedNamesPvalsRatios)
        IndexesOrderedByFC = np.argsort(ArrayData[:,2].astype(np.float))
            
        OrderedNames = ArrayData[IndexesOrderedByFC,0]
        OrderedPvals = ArrayData[IndexesOrderedByFC,1]
        OrderedFCs = ArrayData[IndexesOrderedByFC,2]
    
        OrderedNames_reversed = OrderedNames[::-1]
        OrderedPvals_reversed = OrderedPvals[::-1]
        OrderedFCs_reversed = OrderedFCs[::-1]
            
        OrderedNamesPvalsRatios_FC = np.transpose(np.array([OrderedNames_reversed, OrderedPvals_reversed, OrderedFCs_reversed]))
    
        sd.save_array(OrderedNamesPvalsRatios_FC, './DiffExprGenes_pval_0p05_from_MI_' + Day + '_Ordered_by_FC.csv') # and here
    
        ListDEGsAllCellsForLater = L_min_pval005
        
        
        
        
        
        
        
        
if HistoAndFoldChangesGeneListsGeneralCellProcesses == True: 
    
    print("HistoAndFoldChangesGeneListsGeneralCellTypesProcesses")

    L_ProApoptotic = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Pro-Apoptotic.txt')[0]
    L_AntiApoptotic = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Anti-Apoptotic.txt')[0]
    L_Caspases = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Caspases.txt')[0]

    ListsOfListsOfGenesCellProcesses = [L_ProApoptotic, L_AntiApoptotic, L_Caspases]

    ListOfListsNames_CellProcesses = [
                        'List of\nGenes\nPro Apoptotic',
                        'List of\nGenes\nAnti Apoptotic',
                        'List of\nGenes\nCaspases']

    ListOfListsNames_CellProcesses_Short = [
                        'ProApoptotic',
                        'AntiApoptotic',
                        'Caspases']

    ListOfListsOfGenesForCumulatives = ListsOfListsOfGenesCellProcesses

    NBonferroni = len(ListOfListsOfGenesForCumulatives) * len(ListOfDays)
    
    j = 0
    MatrixOfFoldChanges = []
    MatrixOfSEMforFoldChanges = []
    MatrixOfPvalues = []
    MatrixOfCumulatives_A = []
    MatrixOfCumulatives_B = []
    for DAY in ListOfDays:
        print('Day ' + str(DAY))
        C = C_list[j]
        C_norm = C_norm_list[j]
            
        ListOfFoldChangesThisDay = []
        ListOfSTDEVforFoldChangesThisDay = []
        ListOfPvalues = []
        
        ListOfCumulativesA = []
        ListOfCumulativesB = []
    
        for CurrentList in ListOfListsOfGenesForCumulatives:
    
            CumulativeGeneExpressionA, CumulativeGeneExpressionB = ComputeCumulativeGeneExpressionFromListAndCompareObject(
                    C, C_norm, N_cells, CurrentList)
            
            ListOfCumulativesA.append(CumulativeGeneExpressionA)
            ListOfCumulativesB.append(CumulativeGeneExpressionB)
            
            FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                    CumulativeGeneExpressionA,CumulativeGeneExpressionB)
            
            ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
            ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
            
            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            pvalue = ApplyZtestOnMeansOfDistributions(CumulativeGeneExpressionA, CumulativeGeneExpressionB)
    
            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')
            ListOfPvalues.append(pvalue)
        MatrixOfCumulatives_A.append(ListOfCumulativesA)
        MatrixOfCumulatives_B.append(ListOfCumulativesB)
        MatrixOfFoldChanges.append(ListOfFoldChangesThisDay)
        MatrixOfSEMforFoldChanges.append(ListOfSTDEVforFoldChangesThisDay)
        MatrixOfPvalues.append(ListOfPvalues)
        j = j + 1

    ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
    ListOfDays = ['Day 0','Day 10','Day 14','Day 42']
    LabelUpBy = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    LabelDownBy = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
    Ncols = 4#len(ListOfDays)
    My_ylabel_Fold_Changes = r'Fold Change G2019S / GC'
    PlotFoldChangesManyGenesOneComparison(ListOfListsNames_CellProcesses_Short, 
                                              MatrixOfFoldChanges, 
                                              MatrixOfSEMforFoldChanges, 
                                              MatrixOfPvalues, 
                                              NBonferroni, 
                                              My_ylabel_Fold_Changes,
                                              LabelUpBy,
                                              LabelDownBy,
                                              Ncols,
                                              ListOfDays,
                                              MyTitle = "Fold Change in Expression of Genes from Cell Death lists", 
                                              BarsColors = "Rainbow",
                                              NamesOnEachRectangle = 'No')
    
    plt.savefig('./Results/FigS4I_FoldChanges_CellDeath_ListsOfGenes.pdf')

                
    ########## MAKE HISTOGRAMS ##########     
    Nlines = len(ListOfListsNames_CellProcesses_Short) 
    Ncols = len(ListOfDays)
    ylimVector = [120,120,120,120]
    PlotHistogramsCumulatives(Nlines, Ncols, MatrixOfCumulatives_A, MatrixOfCumulatives_B, 
                              ['','',''], MatrixOfPvalues, NBonferroni, 
                              ['LRRK2G2019S','LRRK2WT'], 'Cumulative Gene Expression (norm. to max)',
                              ylimVector, font3s, Colors=[ColorMU,ColorCT])
    

    print("Done.")
        
    
    
    
    
    
    
    
    
    
if DoGeneGeneCorrelation_CellDeath_VS_CellTypesGenes == 1:
    print("DoGeneGeneCorrelation")
    
    L_stem = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_stem.txt')[0] 
    L_Neurons = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLisa/EnrichedNeuronList.txt')[0]
    L_DA_Neuro = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_DopaminergicNeurons.txt')[0]

    L_Gluta_Neuro = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLisa/L_Glutamatergic_Neurons.txt')[0]
    L_GABA_Neuro = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLisa/L_GABAergic_Neurons.txt')[0]
    L_Seroton_Neuro = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLisa/L_Serotonergic_Neurons.txt')[0]

    WhiteDiagonal = True # True --> mask values (1s) from the diagonal, False --> do not mask
    RemoveNotSignificantCorr = False
    RemoveNANs = True
    DoubleFace = True
    
    # Make the list of genes to be used for correlation
    List_Gene_Names_for_Correlation_BeforeCleaning = L_ProApoptotic  + L_Caspases + L_stem + L_Neurons + L_DA_Neuro + L_Gluta_Neuro + L_GABA_Neuro + L_Seroton_Neuro

#  + L_AntiApoptotic
    # Remove Duplicate Genes
    List_Gene_Names_for_Correlation_BeforeCleaning = list(dict.fromkeys(List_Gene_Names_for_Correlation_BeforeCleaning))
    
    Nsamples = 2  # Here 2 samples, in fact the samples in C object are 0 = MU, 1 = CT

    ##### DETERMINE WHICH GENES ARE IN COMMON BETWEEN MU70, MU35, WT70, WT35.
    ijk = 0
    ListOfListsOfGenesForCorrelationPresentInData = []
    for DAY in ListOfDays:        
        C = C_list[ijk]
        ListGenesC = C.comm_genes() 
        List_Gene_Names_for_Correlation_PresentInData = []
        for GeneInMyList in List_Gene_Names_for_Correlation_BeforeCleaning:
            if GeneInMyList in ListGenesC:
                List_Gene_Names_for_Correlation_PresentInData.append(GeneInMyList)
        ListOfListsOfGenesForCorrelationPresentInData.append(List_Gene_Names_for_Correlation_PresentInData)
        ijk = ijk + 1
        
    List_MuWt_Day0 = ListOfListsOfGenesForCorrelationPresentInData[0]
    List_MuWt_Day10 = ListOfListsOfGenesForCorrelationPresentInData[1]
    List_MuWt_Day14 = ListOfListsOfGenesForCorrelationPresentInData[2]
    List_MuWt_Day42 = ListOfListsOfGenesForCorrelationPresentInData[3]

    Commongenes_Mu_Wt_Days_All = []
    for GeneName in List_MuWt_Day0:
        if GeneName in List_MuWt_Day10 and GeneName in List_MuWt_Day14 and GeneName in List_MuWt_Day42:
            Commongenes_Mu_Wt_Days_All.append(GeneName)

    List_Gene_Names_for_Correlation = Commongenes_Mu_Wt_Days_All

    Ngenes = len(Commongenes_Mu_Wt_Days_All)
    N_Bonf_Corr = (Ngenes*Ngenes-Ngenes)/2 * Nsamples * len(ListOfDays)
    ThrasoldForPval = 0.05/N_Bonf_Corr

    kkk = 0
    for DAY in ListOfDays:        
        print('Day ' + str(DAY))
        
        C = C_list[kkk]
        C_norm = C_norm_list[kkk]

        MyAverageCorrVector = []
        for SampleIndex in range(Nsamples): # Here samples in C object are 0 = MU, 1 = CT
            
            ######### COMPUTE GENE GENE CORRELATION MATRIX #########
            MyCorrelationMatrix_Full, pvalPearsonCorr_Matrix, MyCorrelationMatrix_Masked = ComputeGeneGeneCorrelationOn2SamplesOfCompareObject(C,
                                                                C_norm,
                                                                List_Gene_Names_for_Correlation,
                                                                N_cells,
                                                                ThrasoldForPval,
                                                                RemoveNotSignificantCorr,
                                                                WhiteDiagonal,
                                                                SampleIndex)
            
            ######### PLOT GENE GENE CORRELATION MATRIX #########
            figX, ax = plt.subplots(facecolor="white")
            
            if DoubleFace == False:
                
                if WhiteDiagonal == False:
                    MyCorrelationMatrix = MyCorrelationMatrix_Full
                elif WhiteDiagonal == True:
                    MyCorrelationMatrix = MyCorrelationMatrix_Masked
                else:
                    print("Error in Masking Correlations!")
                
                MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix)
                for i in range(len(List_Gene_Names_for_Correlation)):
                    for j in range(len(List_Gene_Names_for_Correlation)):
                        if j>i:
                            MyCorrelationMatrix_Triangle[i,j] = 0
                        if RemoveNANs and numpy.isnan(MyCorrelationMatrix_Triangle[i,j]):
                            print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                            MyCorrelationMatrix_Triangle[i,j] = 0
                            
            elif DoubleFace == True:
                MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix_Full)
                for i in range(len(List_Gene_Names_for_Correlation)):
                    for j in range(len(List_Gene_Names_for_Correlation)):
                        if j>i:
                            MyCorrelationMatrix_Triangle[i,j] = MyCorrelationMatrix_Masked[i,j]
                        if RemoveNANs and np.isnan(MyCorrelationMatrix_Triangle[i,j]):
                            print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                            MyCorrelationMatrix_Triangle[i,j] = 0
            
            heatmap = ax.pcolor(MyCorrelationMatrix_Triangle, cmap=plt.cm.seismic_r, vmin=-1, vmax=1)# RdGy_r  # PiYG  # seismic   # PRGn_r
            
            ax.invert_yaxis()
            for i in range(len(List_Gene_Names_for_Correlation)):
                for j in range(len(List_Gene_Names_for_Correlation)):
                    ax.text(- 3, 0.8 + i, List_Gene_Names_for_Correlation[i], fontdict=font3ss) # ss
                    ax.text(0.1 + j, - 3, List_Gene_Names_for_Correlation[j], fontdict=font3ss, rotation=90)
                    
            #Spacing between each line
            intervals = 1
            loc = plticker.MultipleLocator(base=intervals)
            ax.xaxis.set_major_locator(loc)
            ax.yaxis.set_major_locator(loc)
            ax.set_xlim(0,len(List_Gene_Names_for_Correlation))
            ax.set_ylim(len(List_Gene_Names_for_Correlation),0)
            #ax.grid(which="major", color="black", linestyle='-', linewidth=1)
            ax.tick_params(labelbottom='off')    
            ax.tick_params(labelleft='off')    
            cbar = figX.colorbar(heatmap)
            cbar.set_label('Pearson Correlation', rotation=90)
            if SampleIndex == 0:
                SampleName = "G2019S"
            elif SampleIndex == 1:
                SampleName = "WT"
            # plt.suptitle("Gene-gene Pearson Correlation Coefficient, " + SampleName + ' at Day ' + str(DAY))
            
            plt.show()

        kkk = kkk + 1
    print("Done.")