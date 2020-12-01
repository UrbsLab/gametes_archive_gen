#-------------------------------------------------------------------------------
# Name:            GAMETES_Archive_Maker.py
# Author:          Tuan Nguyen, Ryan Urbanowicz 
# Created:         07/01/2016
# Description: Designed to make single models at a time, and make models and datasets separately.
#                Assumes use of EDM and not COR.
#GAMETES_use options: model, data, hetdata
#
# Status: First Test Format

#-------------------------------------------------------------------------------
#!/usr/bin/env python

def main(GAMETES_use):
    # General Parameters ##################################################################################
    gametesVersion = 2.2
    addFolderNameText = ""  #optional added text to folder\file names

    loci = [1, 2, 3, 4, 5, 6]
    heritability = [0.4]
    minorAF = [0.2] 
    setK = False #if True, then K value will be specified as a constraint - if False, then K will be allowed to vary.
    K = 0.3  #population prevelance 
    pop_count = 100000 #
    try_count = 10000000 
    quantiles = 1 #only easiest model 

    #Dataset attributes
    samplesize = [1600] #[200, 400, 800, 1600, 3200, 6400] #assumes balanced datasets (#cases = #controls)
    ratio = [50,75]  # have to do the math for both ratio, X and 100-X
    numberofattributes = [20] #[200, 100]
    AF_Min = 0.01
    AF_Max = 0.5
    replicates = 20 #100
    
    
    ###################################################################################
    gametesFileName = "GAMETES_"+str(gametesVersion)+".jar"
    for locus in loci:
        model_folder_name = "GAMETES_"+str(gametesVersion)+"_"+str(addFolderNameText)+"_Models_Loc_"+str(locus)+"_Qnt_"+str(quantiles)+"_Pop_"+str(pop_count)
        if GAMETES_use == 'data':
            data_folder_name = "GAMETES_"+str(gametesVersion)+"_"+str(addFolderNameText)+"_Datasets_Loc_"+str(locus)+"_Qnt_"+str(quantiles)+"_Pop_"+str(pop_count)
        elif GAMETES_use == 'hetdata':
            data_folder_name = "GAMETES_"+str(gametesVersion)+"_"+str(addFolderNameText)+"_Datasets_2Het_Loc_"+str(locus)+"_Qnt_"+str(quantiles)+"_Pop_"+str(pop_count)
            
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
                    jobName = scratchPath+"her_"+str(h)+"__maf_"+str(m)+'_run.sh'
                    pbsFile = open(jobName, 'w')

                    pbsFile.write('#!/bin/bash -l\n') #NEW
                    pbsFile.write('#BSUB -q moore_long\n')
                    pbsFile.write('#BSUB -J '+GAMETES_use+'_'+str(h)+'_'+str(m)+'\n')
                    pbsFile.write('#BSUB -u tnguyen4@swarthmore.edu\n\n')
                    pbsFile.write('#BSUB -o /project/gametes_simulated_data_set/models\n\n')
                    pbsFile.write('#BSUB -e /project/gametes_simulated_data_set/models\n\n')
                    if setK:
                        filewrite = 'time java -jar '+myHome+'/gametes/'+gametesFileName+' -M " -h '+str(h)+' -p '+str(K)
                    else:
                        filewrite = 'time java -jar '+myHome+'/gametes/'+gametesFileName+' -M " -h '+str(h) 
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
            reportName = modelPath+model_folder_name+"/GAMETES_"+str(gametesVersion)+"_"+str(addFolderNameText)+"_Report_Loc_"+str(locus)+"_Qnt_"+str(quantiles)+"_Pop_"+str(pop_count)+".txt"
            report = open(reportName, 'w')
            report.write("GAMETES Version = "+str(gametesVersion)+'\n'+'\n')  
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
                            jobName = scratchPath+str(n)+'_'+str(s)+'_'+str(modelName)+'_run.sh'
                            pbsFile = open(jobName, 'w')
                            pbsFile.write('#!/bin/bash -l\n') #NEW
                            pbsFile.write('#BSUB -q moore_long\n')
                            pbsFile.write('#BSUB -J '+GAMETES_use+'_'+str(h)+'_'+str(m)+'_'+str(s)+'_'+str(n)+'\n')
                            pbsFile.write('#BSUB -u tnguyen4@swarthmore.edu\n\n')
                            pbsFile.write('#BSUB -o /project/gametes_simulated_data_set/models\n\n')
                            pbsFile.write('#BSUB -e /project/gametes_simulated_data_set/models\n\n')
                            pbsFile.write('time java -jar '+myHome+'/gametes/'+gametesFileName+' -i '+str(modelFile)+' -D " -n '+str(AF_Min)+' -x '+str(AF_Max)+' -a '+str(n)+' -s '+str(int(s/2))+' -w '+str(int(s/2))+' -r '+str(replicates)+' -o '+str(genModelName)+'"')

                            pbsFile.close()
                            os.system('bsub < '+jobName)
                            jobCounter +=1
            print((str(jobCounter)+ " jobs submitted."))  

            # Write data generation report #############################################################
            reportName = dataPath+data_folder_name+"/GAMETES_"+str(gametesVersion)+"_"+str(addFolderNameText)+"_Report_Loc_"+str(locus)+"_Qnt_"+str(quantiles)+"_Pop_"+str(pop_count)+".txt"
            report = open(reportName, 'w')
            report.write("GAMETES Version = "+str(gametesVersion)+'\n'+'\n')  
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
                                jobName = scratchPath+str(n)+'_'+str(s)+'_'+str(modelName)+'_run.sh'
                                pbsFile = open(jobName, 'w')
                                pbsFile.write('#!/bin/bash -l\n') #NEW
                                pbsFile.write('#BSUB -q moore_long\n')
                                pbsFile.write('#BSUB -N '+GAMETES_use+'_'+str(h)+'_'+str(m)+'_'+str(s)+'_'+str(n)+'\n')
                                pbsFile.write('#BSUB -u tnguyen4@swarthmore.edu\n\n')
                                pbsFile.write('#BSUB -o /project/gametes_simulated_data_set/models\n\n')
                                pbsFile.write('#BSUB -e /project/gametes_simulated_data_set/models\n\n')
                                pbsFile.write('time java -jar '+myHome+'/gametes/'+gametesFileName+' -i '+str(modelFile)+' -f '+str(r)+' -i '+str(modelFile)+' -f '+str(100-r)+ ' -D " -n '+str(AF_Min)+' -x '+str(AF_Max)+' -a '+str(n)+' -s '+str(int(s/2))+' -w '+str(int(s/2))+' -r '+str(replicates)+' -o '+str(genDataName)+'"')
        
                                pbsFile.close()
                                os.system('bsub < '+jobName)
                                jobCounter +=1
            print((str(jobCounter)+ " jobs submitted."))  

            # Write data generation report #############################################################
            reportName = dataPath+data_folder_name+"/GAMETES_"+str(gametesVersion)+"_"+str(addFolderNameText)+"_Report_Loc_"+str(locus)+"_Qnt_"+str(quantiles)+"_Pop_"+str(pop_count)+".txt"
            report = open(reportName, 'w')
            report.write("GAMETES Version = "+str(gametesVersion)+'\n'+'\n')  
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
    import sys
    import os
    GAMETES_use = sys.argv[1] 
    myHome = os.environ.get('HOME')
    
    modelPath = '/project/gametes_simulated_data_set/models/' 
    dataPath = '/project/gametes_simulated_data_set/datasets/'
    scratchPath = '/project/gametes_simulated_data_set/scratch/tnguyen/'    
    main(GAMETES_use)
