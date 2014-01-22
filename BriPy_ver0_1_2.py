#!/usr/bin/env python
#-------------------------------------------------------------------------------
# Name:        BriPy: NGS data analysis tools for BRIC-seq
# Purpose:     Calculating RNA Half-life from BRIC-seq data
# Version:     0.1.2-beta

# Author:      Naoto Imamachi

# Created:     19/01/2014
# Copyright:   (c) Imamachi 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------

#import_built-in_modules
import sys, os
import argparse
from time import strftime
from time import clock
import shutil

#import_extention_modules
import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

#changing_history_to_this_module
__author__ = "Naoto Imamachi"
__copyright__ = "Copyright 2014. All rights reserved"
__version__ = "0.1.2-beta"
__status__ = "Testing"

def printlog (mesg):
    '''print progress into stderr'''
    mesg = "[" + strftime("%Y-%m-%d %H:%M:%S") + "] " + mesg
    print >>sys.stderr, mesg

def linear(t,a):
    return a*t

def func1(t,a):
    return np.exp(-a*t)

def func2(t,a,b):
    return np.exp(-a*t)+b

def func3(t,a,b,c):
    return a*np.exp(-b*t)+(1-a)*np.exp(-c*t)

def linear_residual(param,t,y):
    a = param[0]
    return y - linear(param,a)

def func1_residual(param,t,y):
    a = param[0]
    return y - func1(t,a)

def func2_residual(param,t,y):
    a = param[0]
    b = param[1]
    return y - func2(t,a,b)

def func3_residual(param,t,y):
    a = param[0]
    b = param[1]
    c = param[2]
    return y - func3(t,a,b,c)

def half_linear(t,a):
    return 0.5 - linear(t,a)

def half1(t,a):
    return 0.5 - func1(t,a)

def half2(t,a,b):
    return 0.5 - func2(t,a,b)

def half3(t,a,b,c):
    return 0.5 - func3(t,a,b,c)

def logplot(xfit,yfit,args_output_prefix,symbol):
    plt.plot(xfit,yfit,'b-')
    plt.title(symbol)
    plt.xlabel('Time')
    plt.ylabel('Relative RNA remaining (0hr = 1)')
    #plt.yscale('log')
    plt.savefig(args_output_prefix + '/' + 'images' + '/' + symbol + '.png')
    plt.clf()

def main():
    #start_time
    StartTime = clock()

    #command_options
    parser = argparse.ArgumentParser(description='BriPy 0.1.2-beta: NGS data analysis tools for BRIC-seq')
    parser.add_argument('-i','--input-file',action='store',dest='input_file',help='RPKM file [required] (default: Cuffdiff output format) ')
    parser.add_argument('-o','--out-prefix',action='store',dest='output_prefix',help='Prefix of output directory')
    parser.add_argument('-t','--time-course',action='store',dest='time_course',help='Time course [required] For example, 0,4,8,12')
    parser.add_argument('--norm-genes',action='store',dest='norm_genes',help='reference genes for normalization (default: GAPDH)')
    parser.add_argument('--min-rpkm-value',action='store',dest='minimum_RPKM_value',help='set minimum_RPKM_value')
    parser.add_argument('--min-time-zero-rpkm-value',action='store',dest='minimum_RPKM_value_time_zero',help='set minimum_RPKM_value in time zero')
    parser.add_argument('--version',action='version',version='%(prog)s ver.0.1.2-beta',help='show this program version')
    args = parser.parse_args()

    #beginning_BriPy_script
    printlog('Beginning BriPy run (v0.1.2-beta)')
    print >>sys.stderr, '-----------------------------------------------'

    #checking_command_errors
    if not (args.input_file and args.time_course):
        print >>sys.stderr, "ERROR: Please set input file and time course!"
        print >>sys.stderr, "e.g.) python BriPy.py -i /path/to/genes.fpkm_tracking -t 0,4,8,12"
        sys.exit(0)
    if not os.path.exists(args.input_file):
        print >>sys.stderr, "ERROR: Input file does not exist!"
        sys.exit(0)

    #prepare_directory
    if not args.output_prefix:
        try:
            args.output_prefix = 'BriPy_out'
            os.mkdir(args.output_prefix)
            os.mkdir(args.output_prefix + '/' + 'images')

        except OSError:
            shutil.rmtree(args.output_prefix)
            os.mkdir(args.output_prefix)
            os.mkdir(args.output_prefix + '/' + 'images')
    else:
        try:
            os.mkdir(args.output_prefix)
            os.mkdir(args.output_prefix + '/' + 'images')

        except OSError:
            shutil.rmtree(args.output_prefix)
            os.mkdir(args.output_prefix)
            os.mkdir(args.output_prefix + '/' + 'images')

    '''get time course data'''
    #time_course
    timeList = (args.time_course.split(','))
    time = map(float,timeList)
    firstTime = float(time[0])
    lastTime = float(time[-1])
    time = np.array(time)
    N = time.size

    '''get reference genes for normalization'''
    #reference genes
    if args.norm_genes:
        geneNormList = (args.norm_genes.split(','))
    else:
        geneNorm = 'GAPDH'

    '''get minimum RPKM value'''
    #minimum RPKM value
    if args.minimum_RPKM_value: #--min-rpkm-value => true
        minRPKM = float(args.minimum_RPKM_value)
    else:                       #--min-rpkm-value => false
        minRPKM = 0

    #minimum RPKM value (time_zero)
    if args.minimum_RPKM_value_time_zero: #--min-time-zero-rpkm-value => true
        minRPKM_zero = float(args.minimum_RPKM_value_time_zero)
    else:                       #--min-time-zero-rpkm-value => false
        minRPKM_zero = 0

    '''get weight factor for normalization'''
    #file_input
    INfile = open(args.input_file,'r')

    #file_output
    OUTwf = open(args.output_prefix + '/' + 'BriPy_weight_factor.log','w')

    #get_weight_factor (house-keeping genes [default: GAPDH])
    wf = []
    printlog("Estimating weight factor for normalization...")
    for line in INfile:
        line = line.rstrip()
        data = line.split("\t")
        if(data[0] == geneNorm):
            #time-exp_length
            dataLength = len(data)
            timeLength = 9 + N*4
            if not(dataLength == timeLength):
                print >>sys.stderr, "ERROR: Time course is not the same as the sample"
                sys.exit(0)

            #extract_each_RPKM_for_estimating_weight_factor
            expList = []
            for point in range(N):
                expList.append(data[9+point*4])
            exp = map(float,expList)
            exp_0h = exp[0]
            wf = [line/exp_0h for line in exp]
            time_print = '\t'.join(timeList)
            wf_print = '\t'.join(map(str,wf))
            print >>OUTwf, 'Time:', '\t' , time_print
            print >>OUTwf, 'Weight_factor:', '\t' , wf_print

    '''calculate RNA half-lives'''
    #file_input
    INfile = open(args.input_file,'r')

    #file_output
    OUTLog = open(args.output_prefix + '/' + 'BriPy_param.log','w')
    OUTPoints = open(args.output_prefix + '/' + 'BriPy_points.log','w')
    OUTResult = open(args.output_prefix + '/' + 'RNA_half-life_calc.result','w')

    printlog("Calculating RNA half-life...")

    #header
    header1_log = "symbol\ta_1\tHalfLife_1\tR2_1\tAIC_1\t"
    header2_log = "a_2\tb_2\tHalfLife_2\tR2_2\tAIC_2\t"
    header3_log = "a_3\tb_3\tc_3\tHalfLife_3\tR2_3\tAIC_3"
    header_result = "symbol\tmodel\tR2\thalf-life"
    print >>OUTLog, header1_log, header2_log, header3_log
    print >>OUTPoints, 'Time: ', '\t' , time_print
    print >>OUTResult, header_result

    for line in INfile:
        try:
            expList = []
            modelSelect = {}
            line = line.rstrip("\n")
            data = line.split("\t")

            #remove_header
            if data[0] == 'tracking_id': continue

            symbol = data[0]
            print >>OUTLog, symbol, '\t',
            print >>OUTPoints, symbol, '\t',
            print >>OUTResult, symbol, '\t',

            #time-exp_length
            dataLength = len(data)
            timeLength = 9 + N*4
            if not(dataLength == timeLength):
                print >>sys.stderr, "ERROR: Time course is not the same as the sample"
                sys.exit(0)

            #extract_exp
            for point in range(N):
                expList.append(data[9+point*4])
            exp = map(float,expList)
            '''skip too few points in exp'''
            if not len(exp) >= 3:
                print >>OUTLog, "NA","\t","NA","\t","NA","\t","NA","\t",
                print >>OUTLog, "NA","\t","NA","\t","NA","\t","NA","\t", "NA","\t",
                print >>OUTLog, "NA","\t","NA","\t","NA","\t","NA","\t", "NA","\t","NA"
                print >>OUTPoints, "too_few_points"
                print >>OUTResult, "too_few_points","\t","NA","\t","NA"
                continue
            elif not max(exp) > minRPKM:
                print >>OUTLog, "NA","\t","NA","\t","NA","\t","NA","\t",
                print >>OUTLog, "NA","\t","NA","\t","NA","\t","NA","\t", "NA","\t",
                print >>OUTLog, "NA","\t","NA","\t","NA","\t","NA","\t", "NA","\t","NA"
                print >>OUTPoints, "low_expression"
                print >>OUTResult, "low_expression","\t","NA","\t","NA"
                continue
            elif not float(exp[0]) > minRPKM_zero:
                print >>OUTLog, "NA","\t","NA","\t","NA","\t","NA","\t",
                print >>OUTLog, "NA","\t","NA","\t","NA","\t","NA","\t", "NA","\t",
                print >>OUTLog, "NA","\t","NA","\t","NA","\t","NA","\t", "NA","\t","NA"
                print >>OUTPoints, "low_expression"
                print >>OUTResult, "low_expression","\t","NA","\t","NA"
                continue
            exp_0h = exp[0]
            exp = [line/exp_0h for line in exp]
            expNorm = [x/y for (x,y) in zip(exp,wf)]
            exp_print = '\t'.join(map(str,expNorm))
            print >>OUTPoints, exp_print
            expNorm = np.array(expNorm)

        except ZeroDivisionError:
            print >>OUTLog, "NA","\t","NA","\t","NA","\t","NA","\t",
            print >>OUTLog, "NA","\t","NA","\t","NA","\t","NA","\t", "NA","\t",
            print >>OUTLog, "NA","\t","NA","\t","NA","\t","NA","\t", "NA","\t","NA"
            print >>OUTPoints, "low_expression"
            print >>OUTResult, "low_expression","\t","NA","\t","NA"
            continue

        #select_parameter
        p0_1 = [1]
        p0_2 = [1,1]
        p0_3 = [1,1,1]
        p0_half = [1]

        #MODEL1 - exp(-a*t)
        try:
            #Fitting_curve
            fitParams1, fitCovariances1 = leastsq(func1_residual,p0_1,args=(time,expNorm))
            a_1 = fitParams1[0]
            print >>OUTLog, a_1, "\t",

            #Calculate_RNA_half-life
            HalfLife1, HalfCovariances1 = leastsq(half1,p0_half,args=(a_1))
            print >>OUTLog, float(HalfLife1), "\t",

            #Calculate_R2
            SSReg1 = np.dot((expNorm-func1(time,*fitParams1)),(expNorm-func1(time,*fitParams1)))
            expMean1 = float(np.mean(expNorm))
            SSTet1 = np.dot((expNorm-expMean1),(expNorm-expMean1))
            if not SSTet1 == 0:
                R2_1 = 1-SSReg1/SSTet1
            else:
                R2_1 = 0
            print >>OUTLog, R2_1, "\t",

            #Akaike Information Criterion(AIC) - Parameter: 1
            p1 = len(fitParams1)
            logLik1 = 0.5*(-N*(np.log2(2*np.pi)+1-np.log2(N)+np.log2(SSReg1)))
            AIC1 = 2*p1 - 2*logLik1
            modelSelect['model1'] = [AIC1,R2_1,HalfLife1]
            print >>OUTLog, AIC1, "\t",

        except RuntimeError:
            print >>sys.stderr, "ERROR: Curve fitting failed for MODEL1."
            print >>OUTLog, symbol,"\t","NA","\t","NA","\t","NA","\t","NA","\t",

        #MODEL2 - exp(-a*t)+b
        try:
            #Fitting_curve
            fitParams2, fitCovariances2 = leastsq(func2_residual,p0_2,args=(time,expNorm))
            a_2 = fitParams2[0]
            b_2 = fitParams2[1]
            print >>OUTLog, a_2, "\t", b_2, "\t",

            #Calculate_RNA_half-life
            HalfLife2, HalfCovariances2 = leastsq(half2,p0_half,args=(a_2,b_2))
            print >>OUTLog, float(HalfLife2), "\t",

            #Calculate_R2
            SSReg2 = np.dot((expNorm-func2(time,*fitParams2)),(expNorm-func2(time,*fitParams2)))
            expMean2 = float(np.mean(expNorm))
            SSTet2 = np.dot((expNorm-expMean2),(expNorm-expMean2))
            if not SSTet2 == 0:
                R2_2 = 1-SSReg2/SSTet2
            else:
                R2_2 = 0
            print >>OUTLog, R2_2, "\t",

            #Akaike Information Criterion(AIC) - Parameter: 2
            p2 = len(fitParams2)
            logLik2 = 0.5*(-N*(np.log2(2*np.pi)+1-np.log2(N)+np.log2(SSReg2)))
            AIC2 = 2*p2 - 2*logLik2
            modelSelect['model2'] = [AIC2,R2_2,HalfLife2]
            print >>OUTLog, AIC2, "\t",

        except RuntimeError:
            print >>sys.stderr, "ERROR: Curve fitting failed for MODEL2."
            print >>OUTLog, "NA","\t","NA","\t","NA","\t","NA","\t" "NA","\t",

        #MODEL3 - a*exp(-b*t)+(1-a)*exp(-c*t)
        try:
            #Fitting_curve
            fitParams3, fitCovariances3 = leastsq(func3_residual,p0_3,args=(time,expNorm))
            a_3 = fitParams3[0]
            b_3 = fitParams3[1]
            c_3 = fitParams3[2]
            print >>OUTLog, a_3, "\t", b_3, "\t", c_3, "\t",

            #Calculate_RNA_half-life
            HalfLife3, HalfCovariances3 = leastsq(half3,p0_half,args=(a_3,b_3,c_3))
            print >>OUTLog, float(HalfLife3), "\t",

            #Calculate_R2
            SSReg3 = np.dot((expNorm-func3(time,*fitParams3)),(expNorm-func3(time,*fitParams3)))
            expMean3 = float(np.mean(expNorm))
            SSTet3 = np.dot((expNorm-expMean3),(expNorm-expMean3))
            if not SSTet3 == 0:
                R2_3 = 1-SSReg3/SSTet3
            else:
                R2_3 = 0
            print >>OUTLog, R2_3, "\t",

            #Akaike Information Criterion(AIC) - Parameter: 3
            p3 = len(fitParams3)
            logLik3 = 0.5*(-N*(np.log2(2*np.pi)+1-np.log2(N)+np.log2(SSReg3)))
            AIC3 = 2*p3 - 2*logLik3
            modelSelect['model3'] = [AIC3,R2_3,HalfLife3]
            print >>OUTLog, AIC3

        except RuntimeError:
            print >>sys.stderr, "ERROR: Curve fitting failed for MODEL3."
            print >>OUTLog, "NA","\t","NA","\t","NA","\t","NA","\t" "NA","\t","NA"

        #Model_Selection (AIC)
        if modelSelect: #checking if dictionary is empty or not.
            bestModel = min(modelSelect.items(),key=lambda x:x[1])[0]
            bestR2 = min(modelSelect.items(),key=lambda x:x[1])[1][1]
            bestHalf = min(modelSelect.items(),key=lambda x:x[1])[1][2]
            print >>OUTResult, bestModel, "\t",
            print >>OUTResult, bestR2, "\t",
            print >>OUTResult, float(bestHalf)
        else:
            print >>OUTResult, symbol,"\t", "Curve_fitting_failed","\t","NA","\t","NA"

        #plot_fitting_curve
        plt.plot(time,expNorm,'bo  ')
        xfit = np.linspace(firstTime,lastTime)
        if bestModel == 'model1':
            yfit = func1(xfit,a_1)
            logplot(xfit,yfit,args.output_prefix,symbol)
        elif bestModel == 'model2':
            yfit = func2(xfit,a_2,b_2)
            logplot(xfit,yfit,args.output_prefix,symbol)
        elif bestModel == 'model3':
            yfit = func3(xfit,a_3,b_3,c_3)
            logplot(xfit,yfit,args.output_prefix,symbol)
        else:
            pass

    INfile.close()
    OUTLog.close()
    OUTResult.close()

    #end_time
    EndTime = clock()
    TimeLog = EndTime - StartTime
    TimeLog = float(TimeLog)/60
    TimeLog = str(TimeLog)
    printlog('Run complete: ' + TimeLog + 'min elapsed')


if __name__ == '__main__':
    main()
