#-------------------------------------------------------------------------------
# Name:            GAMETES_Archive_Maker.py
# Author:          Tuan Nguyen, Ryan Urbanowicz 
# Created:         07/01/2016
# Description: Designed to make single models at a time, and make models and datasets separately.
#                Assumes use of EDM and not COR.
#GAMETES_use options: model, data, hetdata
#Example Command: python ./gametes_archive.py --use model 
# Status: First Test Format

#-------------------------------------------------------------------------------
#!/usr/bin/env python
import sys
import os
import argparse
import time


def main(argv):
    #Ryan's Hard Coding (Ignore if you are not Ryan) ----------------------------------------------
    wrtPath = '/home/ryanurb/idata/'

    #ARGUMENTS:------------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='')
    #On/Off Analysis Options:----------------------------------------------------------------------
    parser.add_argument('--use', dest='GAMETES_use', help='', type=str, default ='model') #defaults to model generation
    parser.add_argument('--write-path', dest='writePath', help='', type=str, default = 'RyansHardCode') #full path/filename
    parser.add_argument('--namemod', dest='addFolderNameText', help='', type=str, default ='') #defaults to model generation
     
    options=parser.parse_args(argv[1:])
    #On/Off Analysis Options:----------------------------------------------------------------------
    GAMETES_use = options.GAMETES_use
    #Required Command Line: (No Default)-----------------------------------------------------------
    if options.writePath == 'RyansHardCode':
        wrtPath = wrtPath
    else:
        wrtPath = options.writePath + '/'
    addFolderNameText = options.addFolderNameText
    
    modelPath = wrtPath+'models/' 
    dataPath = wrtPath+'datasets/'
    scratchPath = wrtPath+'scratch/'  
    logPath = wrtPath+'logs/'  
    
    #Folder Management------------------------------
    #Model Path-----------------
    if not os.path.exists(modelPath):
        os.mkdir(modelPath)  
    #Dataset Path--------------------
    if not os.path.exists(dataPath):
        os.mkdir(dataPath)    
    #Scratch Path-------------------- 
    if not os.path.exists(scratchPath):
        os.mkdir(scratchPath) 
    #LogFile Path--------------------
    if not os.path.exists(logPath):
        os.mkdir(logPath) 
    #----------------------------------------------

    # General Parameters ##################################################################################
    locus = 2
    heritability = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4]
    minorAF = [0.2] 
    setK = True #if True, then K value will be specified as a constraint - if False, then K will be allowed to vary.
    K = 0.3  #population prevelance 
    pop_count = 100000 #
    try_count = 10000000 
    quantiles = 2

    #Dataset attributes
    samplesize = [200, 400, 800, 1600] #[200, 400, 800, 1600, 3200, 6400] #assumes balanced datasets (#cases = #controls)
    ratio = [50,75]  # have to do the math for both ratio, X and 100-X
    numberofattributes = [20] # [20, 100, 1000, 10000, 100000] #[200, 100]
    AF_Min = 0.01
    AF_Max = 0.5
    replicates = 100 #100
    
    
    ###################################################################################
    gametesFileName = "gametes_"+str(version)+".jar"
    model_folder_name = "GAMETES_"+str(version)+"_"+str(addFolderNameText)+"_Models_Loc_"+str(locus)+"_Qnt_"+str(quantiles)+"_Pop_"+str(pop_count)
    if GAMETES_use == 'data':
        data_folder_name = "GAMETES_"+str(version)+"_"+str(addFolderNameText)+"_Datasets_Loc_"+str(locus)+"_Qnt_"+str(quantiles)+"_Pop_"+str(pop_count)
    elif GAMETES_use == 'hetdata':
        data_folder_name = "GAMETES_"+str(version)+"_"+str(addFolderNameText)+"_Datasets_2Het_Loc_"+str(locus)+"_Qnt_"+str(quantiles)+"_Pop_"+str(pop_count)
        
    if GAMETES_use == "model":
        print("GAMETES: Generating Models")
        jobCounter = 0
        if not os.path.exists(modelPath+model_folder_name):
            os.mkdir(modelPath+model_folder_name)  
        for h in heritability:
            print(h)
            for m in minorAF:           
                print(m) 
                genModelName = modelPath+model_folder_name+"/her_"+str(h)+"__maf_"+str(m)
                #MAKE CLUSTER JOBS###################################################################
                jobName = scratchPath+'gametes_'+"her_"+str(h)+"__maf_"+str(m)+'_run.sh'
                pbsFile = open(jobName, 'w')

                pbsFile.write('#!/bin/bash -l\n') #NEW
                pbsFile.write('#BSUB -q moore_long\n')
                pbsFile.write('#BSUB -J '+GAMETES_use+'_'+str(h)+'_'+str(m)+'\n')
                pbsFile.write('#BSUB -u ryanurb@upenn.edu\n\n')
                pbsFile.write('#BSUB -o ' + logPath+GAMETES_use+'_'+str(h)+'_'+str(m)+'_'+str(time.time())+'.o\n')
                pbsFile.write('#BSUB -e ' + logPath+GAMETES_use+'_'+str(h)+'_'+str(m)+'_'+str(time.time())+'.e\n\n')
                if setK:
                    filewrite = 'time java -jar '+gametesPath+gametesFileName+' -M " -h '+str(h)+' -p '+str(K)
                else:
                    filewrite = 'time java -jar '+gametesPath+gametesFileName+' -M " -h '+str(h) 
                for i in range(locus):
                    filewrite = filewrite +' -a '+str(m)
                filewrite = filewrite +' -o '+str(genModelName)+'.txt'+'" -q '+str(quantiles)+' -p '+str(pop_count)+' -t '+str(try_count)
                pbsFile.write(filewrite)
                pbsFile.close()
                os.system('bsub < '+jobName)
                #####################################################################################  
                jobCounter +=1
        print((str(jobCounter)+ " jobs submitted."))   
         
        # Write model generation report #############################################################
        reportName = modelPath+model_folder_name+"/GAMETES_"+str(version)+"_"+str(addFolderNameText)+"_Report_Loc_"+str(locus)+"_Qnt_"+str(quantiles)+"_Pop_"+str(pop_count)+".txt"
        report = open(reportName, 'w')
        report.write("GAMETES Version = "+str(version)+'\n'+'\n')  
        report.write("Number of loci = "+str(locus)+'\n')  
        report.write("Heritabilities Generated = "+str(heritability)+'\n')  
        report.write("Minor Allele Frequencies Generated = "+str(minorAF)+'\n')  
        report.write("Used a Fixed K = "+str(setK)+'\n')  
        report.write("K (makes no difference if above is False) = "+str(K)+'\n')  
        report.write("Model population size = "+str(pop_count)+'\n') 
        report.write("Tries to reach population size = "+str(try_count)+'\n') 
        report.write("Quantiles (number of models selected from population) = "+str(quantiles)+'\n')                     
                                    
    elif GAMETES_use == "data":
        print("GAMETES: Generating Datasets")
        jobCounter = 0
        if not os.path.exists(dataPath+data_folder_name):
            os.mkdir(dataPath+data_folder_name) 
        for n in numberofattributes:
            print(n)
            if not os.path.exists(dataPath+data_folder_name+'/'+"a_"+str(n)):
                os.mkdir(dataPath+data_folder_name+'/'+"a_"+str(n))  
            for s in samplesize: 
                print(s)
                if not os.path.exists(dataPath+data_folder_name+'/'+"a_"+str(n)+'/'+"s_"+str(s)):
                    os.mkdir(dataPath+data_folder_name+'/'+"a_"+str(n)+'/'+"s_"+str(s))  
                for h in heritability:
                    print(h)
                    for m in minorAF:
                        modelName = "her_"+str(h)+"__maf_"+str(m)
                        dataName = "a_"+str(n)+"s_"+str(s)+str(modelName)
                        if not os.path.exists(dataPath+data_folder_name+'/'+"a_"+str(n)+'/'+"s_"+str(s)+'/'+str(modelName)):
                            os.mkdir(dataPath+data_folder_name+'/'+"a_"+str(n)+'/'+"s_"+str(s)+'/'+str(modelName)) 
                            
                        modelFile = modelPath+model_folder_name+"/her_"+str(h)+"__maf_"+str(m)+"_Models.txt"
                        genModelName = dataPath+data_folder_name+'/'+"a_"+str(n)+'/'+"s_"+str(s)+'/'+str(modelName)+'/'+dataName
                        jobName = scratchPath+'gametes_'+str(n)+'_'+str(s)+'_'+str(modelName)+'_run.sh'
                        pbsFile = open(jobName, 'w')
                        pbsFile.write('#!/bin/bash -l\n') #NEW
                        pbsFile.write('#BSUB -q moore_long\n')
                        pbsFile.write('#BSUB -J '+GAMETES_use+'_'+str(h)+'_'+str(m)+'_'+str(s)+'_'+str(n)+'\n')
                        pbsFile.write('#BSUB -u ryanu@upenn.edu\n\n')
                        pbsFile.write('#BSUB -o ' + logPath+GAMETES_use+'_'+str(h)+'_'+str(m)+'_'+str(s)+'_'+str(n)+'_'+str(time.time())+'.o\n')
                        pbsFile.write('#BSUB -e ' + logPath+GAMETES_use+'_'+str(h)+'_'+str(m)+'_'+str(s)+'_'+str(n)+'_'+str(time.time())+'.e\n\n')
                        pbsFile.write('time java -jar '+gametesPath+gametesFileName+' -i '+str(modelFile)+' -D " -n '+str(AF_Min)+' -x '+str(AF_Max)+' -a '+str(n)+' -s '+str(int(s/2))+' -w '+str(int(s/2))+' -r '+str(replicates)+' -o '+str(genModelName)+'"')
                        pbsFile.close()
                        os.system('bsub < '+jobName)
                        jobCounter +=1
        print((str(jobCounter)+ " jobs submitted."))  

        # Write data generation report #############################################################
        reportName = dataPath+data_folder_name+"/GAMETES_"+str(version)+"_"+str(addFolderNameText)+"_Report_Loc_"+str(locus)+"_Qnt_"+str(quantiles)+"_Pop_"+str(pop_count)+".txt"
        report = open(reportName, 'w')
        report.write("GAMETES Version = "+str(version)+'\n'+'\n')  
        report.write("Model Archive Used = "+str(model_folder_name)+'\n'+'\n')  

        report.write("Sample Sizes Generated = "+str(samplesize)+'\n')   
        report.write("Number of Attributes Generated = "+str(numberofattributes)+'\n')  
        report.write("Minimum Minor Allele Frequency for Noise Atributes = "+str(AF_Min)+'\n')  
        report.write("Maximum Minor Allele Frequency for Noise Atributes = "+str(AF_Max)+'\n')  
        report.write("Dataset Replicates Generated = "+str(replicates)+'\n') 
    
    elif GAMETES_use == "hetdata": #Script only set up to use the same model twice in the same dataset.
        print("GAMETES: Generating Heterogeneous Datasets")
        jobCounter = 0
        if not os.path.exists(dataPath+data_folder_name):
            os.mkdir(dataPath+data_folder_name) 
        for n in numberofattributes:
            print(n)
            if not os.path.exists(dataPath+data_folder_name+'/'+"a_"+str(n)):
                os.mkdir(dataPath+data_folder_name+'/'+"a_"+str(n))  
            for s in samplesize: 
                print(s)
                if not os.path.exists(dataPath+data_folder_name+'/'+"a_"+str(n)+'/'+"s_"+str(s)):
                    os.mkdir(dataPath+data_folder_name+'/'+"a_"+str(n)+'/'+"s_"+str(s))  
                for h in heritability:
                    print(h)
                    for m in minorAF:
                        print(m)
                        modelName = "h_"+str(h)+"MAF_"+str(m)
                        if not os.path.exists(dataPath+data_folder_name+'/'+"a_"+str(n)+'/'+"s_"+str(s)+'/'+str(modelName)):
                            os.mkdir(dataPath+data_folder_name+'/'+"a_"+str(n)+'/'+"s_"+str(s)+'/'+str(modelName)) 
                        for r in ratio:
                            print(r)
                            if not os.path.exists(dataPath+data_folder_name+'/'+"a_"+str(n)+'/'+"s_"+str(s)+'/'+str(modelName)+'/'+"r_"+str(r)):
                                os.mkdir(dataPath+data_folder_name+'/'+"a_"+str(n)+'/'+"s_"+str(s)+'/'+str(modelName)+'/'+"r_"+str(r)) 
                            
                            dataName = "a_"+str(n)+"s_"+str(s)+'_Het_'+str(modelName)+'_'+"r_"+str(r)
                            modelFile = modelPath+model_folder_name+"/her_"+str(h)+"__maf_"+str(m)+"_Models.txt"
                            genDataName = dataPath+data_folder_name+'/'+"a_"+str(n)+'/'+"s_"+str(s)+'/'+str(modelName)+'/'+"r_"+str(r)+'/'+dataName
                            jobName = scratchPath+'gametes_'+str(n)+'_'+str(s)+'_'+str(modelName)+'_run.sh'
                            pbsFile = open(jobName, 'w')
                            pbsFile.write('#!/bin/bash -l\n') #NEW
                            pbsFile.write('#BSUB -q moore_long\n')
                            pbsFile.write('#BSUB -N '+GAMETES_use+'_'+str(h)+'_'+str(m)+'_'+str(s)+'_'+str(n)+'_'+str(r)+'\n') 
                            pbsFile.write('#BSUB -u ryanurb@upenn.edu\n\n')
                            pbsFile.write('#BSUB -o ' + logPath+GAMETES_use+'_'+str(h)+'_'+str(m)+'_'+str(s)+'_'+str(n)+'_'+str(r)+'_'+str(time.time())+'.o\n')
                            pbsFile.write('#BSUB -e ' + logPath+GAMETES_use+'_'+str(h)+'_'+str(m)+'_'+str(s)+'_'+str(n)+'_'+str(r)+'_'+str(time.time())+'.e\n\n')
                            pbsFile.write('time java -jar '+gametesPath+gametesFileName+' -i '+str(modelFile)+' -f '+str(r)+' -i '+str(modelFile)+' -f '+str(100-r)+ ' -D " -n '+str(AF_Min)+' -x '+str(AF_Max)+' -a '+str(n)+' -s '+str(int(s/2))+' -w '+str(int(s/2))+' -r '+str(replicates)+' -o '+str(genDataName)+'"')
    
                            pbsFile.close()
                            os.system('bsub < '+jobName)
                            jobCounter +=1
        print((str(jobCounter)+ " jobs submitted."))  

        # Write data generation report #############################################################
        reportName = dataPath+data_folder_name+"/GAMETES_"+str(version)+"_"+str(addFolderNameText)+"_Report_Loc_"+str(locus)+"_Qnt_"+str(quantiles)+"_Pop_"+str(pop_count)+".txt"
        report = open(reportName, 'w')
        report.write("GAMETES Version = "+str(version)+'\n'+'\n')  
        report.write("Model Archive Used = "+str(model_folder_name)+'\n')  
        report.write("Heterogeneous datasets generated using the same model to generate independent samples for two different pairs of attributes in the dataset."+'\n')  
        report.write("Sample Sizes Generated = "+str(samplesize)+'\n')  
        report.write("Ratios Generated = "+str(ratio)+'\n')  
        report.write("Number of Attributes Generated = "+str(numberofattributes)+'\n')  
        report.write("Minimum Minor Allele Frequency for Noise Atributes = "+str(AF_Min)+'\n')  
        report.write("Maximum Minor Allele Frequency for Noise Atributes = "+str(AF_Max)+'\n')  
        report.write("Dataset Replicates Generated = "+str(replicates)+'\n') 
        
    else:
        print("Specified function not found:  Please enter either 'model' or 'data'. ")
    
    
if __name__=="__main__":
    
    wdrPath = os.path.dirname(os.path.realpath(__file__)) #Grabs the path to the location of this file from whereever it is run from. 
    wdrPath = wdrPath.split('/')
    gametesPath = ''
    for i in range(0,len(wdrPath)-1):
        gametesPath = gametesPath+wdrPath[i]+'/'

    version = wdrPath[-2] #location of the version folder name
    
    sys.exit(main(sys.argv))
