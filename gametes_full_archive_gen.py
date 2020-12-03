"""
Author: Ryan Urbanowicz
Created: 11/30/20
Description: Script to apply GAMETES to generate and organize a variety of SNP simulated models and a corresponding datasets
"""

import sys
import os
import argparse
import time

def main(argv):
    #Parse arguments
    parser = argparse.ArgumentParser(description="")
    #No defaults
    parser.add_argument('--output-path',dest='output_path',type=str,help='path to output directory')
    parser.add_argument('--archive-name', dest='archive_name',type=str, help='name of archive output folder (no spaces)')
    parser.add_argument('--run-parallel',dest='run_parallel',type=str,help='path to directory containing datasets',default="True")
    parser.add_argument('--use', dest='use', help='', type=str, default ='model') #defaults to model generation

    options = parser.parse_args(argv[1:])
    output_path = options.output_path
    archive_name = options.archive_name
    run_parallel = options.run_parallel
    use = options.use

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    this_file_path = os.path.dirname(os.path.realpath(__file__))

    model_dest = output_path+'/'+archive_name+'/models'
    job_dest = output_path+'/temporary'+'/jobs'
    log_dest = output_path+'/temporary'+'/logs'

    #Create folders
    if not os.path.exists(output_path+'/'+archive_name):
        os.mkdir(output_path+'/'+archive_name)
    if not os.path.exists(model_dest):
        os.mkdir(model_dest)
    if not os.path.exists(output_path+'/temporary'):
        os.mkdir(output_path+'/temporary')
    if not os.path.exists(job_dest):
        os.mkdir(job_dest)
    if not os.path.exists(log_dest):
        os.mkdir(log_dest)

    if use == 'model':
        #Generate core main effect models
        univariate_core_model(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path)
        #Generate core 2-way epistasis models
        epistasis_2_locus_core_model(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path)
        #Generate 3-way epistasis models
        epistasis_3_locus_model(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path)

    elif use == 'data':
        #Generate core main effect data
        univariate_core_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path)
        #Generate core epistasis data
        epistasis_2_locus_core_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path)
        #Generate 3-way epistasis data
        epistasis_3_locus_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path)
        #Generate heterogeneous data (2 subgroups of 2-way epistasis)
        epistasis_2_locus_hetero_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path)
        #Generate additive data (2 additively combined 2-way epistasis models, yielding 'impure' epistasis)
        epistasis_2_locus_additive_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path)
        #Generate heterogeneous data (2 subgroups of univariate efects)
        univariate_2_locus_hetero_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path)
        #Generate additive data (2 subgroups of univariate efects)
        univariate_2_locus_additive_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path)
        #Generate heterogeneous data (4 subgroups of univariate efects)
        univariate_4_locus_hetero_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path)
        #Generate additive data (4 subgroups of univariate efects)
        univariate_4_locus_additive_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path)
        #Generate imbalanced dataset (with 2-way epistasis)
        epistasis_2_locus_imbalanced_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path)
        #Generate continuous endpoint data (with 2-way epistasis)
        epistasis_2_locus_quantitative_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path)
        #Generate increasing feature count datasets (with 2-way epistasis)
        epistasis_2_locus_numfeatures_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path)
    else:
        print("GAMETES use not recognized.")

def univariate_core_model(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path):
    #Define model parameters
    locus = 1
    heritability = [0.05, 0.1, 0.2, 0.4]
    minorAF = [0.2]
    setK = True #if True, then K value will be specified as a constraint - if False, then K will be allowed to vary.
    K = 0.3  #population prevelance
    pop_count = 100000 #
    try_count = 10000000
    quantiles = 2

    #Generate models
    for h in heritability:
        for m in minorAF:
            model_path_name = model_dest+"/L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)
            #Create gametes run command
            if setK:
                model_path_name = model_path_name+'_K_'+str(K)
                filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar '+'-M " -h '+str(h)+' -p '+str(K)
            else:
                filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar '+'-M " -h '+str(h)
            for i in range(locus):
                filewrite = filewrite +' -a '+str(m)
            filewrite = filewrite +' -o '+model_path_name+'.txt'+'" -q '+str(quantiles)+' -p '+str(pop_count)+' -t '+str(try_count)

            if run_parallel:
                job_ref = str(time.time())
                job_path_name = job_dest+'/gametes_'+"L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_'+job_ref+'_run.sh'
                sh_file = open(job_path_name,'w')
                sh_file.write('#!/bin/bash\n')
                sh_file.write('#BSUB -q i2c2_normal'+'\n')
                sh_file.write('#BSUB -J '+job_ref+'\n')
                sh_file.write('#BSUB -R "rusage[mem=4G]"'+'\n')
                sh_file.write('#BSUB -M 15GB'+'\n')
                sh_file.write('#BSUB -o ' + log_dest+'/gametes_'+"L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_'+job_ref+'.o\n')
                sh_file.write('#BSUB -e ' + log_dest+'/gametes_'+"L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_'+job_ref+'.e\n')
                sh_file.write(filewrite)
                sh_file.close()
                os.system('bsub < '+job_path_name)
                pass
            else:
                os.system(filewrite)

def epistasis_2_locus_core_model(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path):
    #Define model parameters
    locus = 2
    heritability = [0.05, 0.1, 0.2, 0.4]
    minorAF = [0.2]
    setK = True #if True, then K value will be specified as a constraint - if False, then K will be allowed to vary.
    K = 0.3  #population prevelance
    pop_count = 100000 #
    try_count = 10000000
    quantiles = 2

    #Generate models
    for h in heritability:
        for m in minorAF:
            model_path_name = model_dest+"/L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)
            #Create gametes run command
            if setK:
                model_path_name = model_path_name+'_K_'+str(K)
                filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar '+'-M " -h '+str(h)+' -p '+str(K)
            else:
                filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar '+'-M " -h '+str(h)
            for i in range(locus):
                filewrite = filewrite +' -a '+str(m)
            filewrite = filewrite +' -o '+model_path_name+'.txt'+'" -q '+str(quantiles)+' -p '+str(pop_count)+' -t '+str(try_count)

            if run_parallel:
                job_ref = str(time.time())
                job_path_name = job_dest+'/gametes_'+"L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_'+job_ref+'_run.sh'
                sh_file = open(job_path_name,'w')
                sh_file.write('#!/bin/bash\n')
                sh_file.write('#BSUB -q i2c2_normal'+'\n')
                sh_file.write('#BSUB -J '+job_ref+'\n')
                sh_file.write('#BSUB -R "rusage[mem=4G]"'+'\n')
                sh_file.write('#BSUB -M 15GB'+'\n')
                sh_file.write('#BSUB -o ' + log_dest+'/gametes_'+"L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_'+job_ref+'.o\n')
                sh_file.write('#BSUB -e ' + log_dest+'/gametes_'+"L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_'+job_ref+'.e\n')
                sh_file.write(filewrite)
                sh_file.close()
                os.system('bsub < '+job_path_name)
                pass
            else:
                os.system(filewrite)

def epistasis_3_locus_model(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path):
    #Define model parameters
    locus = 3
    heritability = [0.2]
    minorAF = [0.2]
    setK = True #if True, then K value will be specified as a constraint - if False, then K will be allowed to vary.
    K = 0.3  #population prevelance
    pop_count = 100000 #
    try_count = 100000000
    quantiles = 2

    #Generate models
    for h in heritability:
        for m in minorAF:
            model_path_name = model_dest+"/L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)
            #Create gametes run command
            if setK:
                model_path_name = model_path_name+'_K_'+str(K)
                filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar '+'-M " -h '+str(h)+' -p '+str(K)
            else:
                filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar '+'-M " -h '+str(h)
            for i in range(locus):
                filewrite = filewrite +' -a '+str(m)
            filewrite = filewrite +' -o '+model_path_name+'.txt'+'" -q '+str(quantiles)+' -p '+str(pop_count)+' -t '+str(try_count)

            if run_parallel:
                job_ref = str(time.time())
                job_path_name = job_dest+'/gametes_'+"L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_'+job_ref+'_run.sh'
                sh_file = open(job_path_name,'w')
                sh_file.write('#!/bin/bash\n')
                sh_file.write('#BSUB -q i2c2_normal'+'\n')
                sh_file.write('#BSUB -J '+job_ref+'\n')
                sh_file.write('#BSUB -R "rusage[mem=4G]"'+'\n')
                sh_file.write('#BSUB -M 15GB'+'\n')
                sh_file.write('#BSUB -o ' + log_dest+'/gametes_'+"L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_'+job_ref+'.o\n')
                sh_file.write('#BSUB -e ' + log_dest+'/gametes_'+"L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_'+job_ref+'.e\n')
                sh_file.write(filewrite)
                sh_file.close()
                os.system('bsub < '+job_path_name)
                pass
            else:
                os.system(filewrite)

def univariate_core_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path):
    #Model parameters needed
    locus = 1
    heritability = [0.05, 0.1, 0.2, 0.4]
    minorAF = [0.2]
    K = 0.3  #population prevelance
    #Define dataset parameters
    data_name = 'gametes_univariate'
    samplesize = [200, 400, 800, 1600] #[200, 400, 800, 1600, 3200, 6400] #assumes balanced datasets (#cases = #controls)
    numberofattributes = [100] # [20, 100, 1000, 10000, 100000] #[200, 100]
    AF_Min = 0.01
    AF_Max = 0.5
    replicates = 30 #100

    #Make dataset folder
    if not os.path.exists(output_path+'/'+archive_name+'/'+data_name):
        os.mkdir(output_path+'/'+archive_name+'/'+data_name)

    #Generate datasets and folders
    for n in numberofattributes:
        for s in samplesize:
            for h in heritability:
                for m in minorAF:
                    modelName = "L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_K_'+str(K)
                    modelFile = model_dest+'/'+modelName+"_Models.txt"
                    genDataName = output_path+'/'+archive_name+'/'+data_name+'/'+data_name+'_A_'+str(n)+'_S_'+str(s)+'_'+str(modelName)
                    #Create gametes run command
                    filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar -i '+modelFile+' -D "-n '+str(AF_Min)+' -x '+str(AF_Max)+' -a '+str(n)+' -s '+str(int(s/2))+' -w '+str(int(s/2))+' -r '+str(replicates)+' -o '+str(genDataName)+'"'

                    if run_parallel:
                        job_ref = str(time.time())
                        job_path_name = job_dest+'/gametes_'+data_name+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'_run.sh'
                        sh_file = open(job_path_name,'w')
                        sh_file.write('#!/bin/bash\n')
                        sh_file.write('#BSUB -q i2c2_normal'+'\n')
                        sh_file.write('#BSUB -J '+job_ref+'\n')
                        sh_file.write('#BSUB -R "rusage[mem=4G]"'+'\n')
                        sh_file.write('#BSUB -M 15GB'+'\n')
                        sh_file.write('#BSUB -o ' + log_dest+'/gametes_'+data_name+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.o\n')
                        sh_file.write('#BSUB -e ' + log_dest+'/gametes_'+data_name+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.e\n')
                        sh_file.write(filewrite)
                        sh_file.close()
                        os.system('bsub < '+job_path_name)
                        pass
                    else:
                        os.system(filewrite)

def epistasis_2_locus_core_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path):
    #Model parameters needed
    locus = 2
    heritability = [0.05, 0.1, 0.2, 0.4]
    minorAF = [0.2]
    K = 0.3  #population prevelance
    #Define dataset parameters
    data_name = 'gametes_2way_epistasis'
    samplesize = [200, 400, 800, 1600] #[200, 400, 800, 1600, 3200, 6400] #assumes balanced datasets (#cases = #controls)
    numberofattributes = [100] # [20, 100, 1000, 10000, 100000] #[200, 100]
    AF_Min = 0.01
    AF_Max = 0.5
    replicates = 30 #100

    #Make dataset folder
    if not os.path.exists(output_path+'/'+archive_name+'/'+data_name):
        os.mkdir(output_path+'/'+archive_name+'/'+data_name)

    #Generate datasets and folders
    for n in numberofattributes:
        for s in samplesize:
            for h in heritability:
                for m in minorAF:
                    modelName = "L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_K_'+str(K)
                    modelFile = model_dest+'/'+modelName+"_Models.txt"
                    genDataName = output_path+'/'+archive_name+'/'+data_name+'/'+data_name+'_A_'+str(n)+'_S_'+str(s)+'_'+str(modelName)
                    #Create gametes run command
                    filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar -i '+modelFile+' -D "-n '+str(AF_Min)+' -x '+str(AF_Max)+' -a '+str(n)+' -s '+str(int(s/2))+' -w '+str(int(s/2))+' -r '+str(replicates)+' -o '+str(genDataName)+'"'

                    if run_parallel:
                        job_ref = str(time.time())
                        job_path_name = job_dest+'/gametes_'+data_name+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'_run.sh'
                        sh_file = open(job_path_name,'w')
                        sh_file.write('#!/bin/bash\n')
                        sh_file.write('#BSUB -q i2c2_normal'+'\n')
                        sh_file.write('#BSUB -J '+job_ref+'\n')
                        sh_file.write('#BSUB -R "rusage[mem=4G]"'+'\n')
                        sh_file.write('#BSUB -M 15GB'+'\n')
                        sh_file.write('#BSUB -o ' + log_dest+'/gametes_'+data_name+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.o\n')
                        sh_file.write('#BSUB -e ' + log_dest+'/gametes_'+data_name+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.e\n')
                        sh_file.write(filewrite)
                        sh_file.close()
                        os.system('bsub < '+job_path_name)
                        pass
                    else:
                        os.system(filewrite)

def epistasis_3_locus_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path):
    #Model parameters needed
    locus = 3
    heritability = [0.2]
    minorAF = [0.2]
    K = 0.3  #population prevelance
    #Define dataset parameters
    data_name = 'gametes_3way_epistasis'
    samplesize = [1600] #[200, 400, 800, 1600, 3200, 6400] #assumes balanced datasets (#cases = #controls)
    numberofattributes = [100] # [20, 100, 1000, 10000, 100000] #[200, 100]
    AF_Min = 0.01
    AF_Max = 0.5
    replicates = 30 #100

    #Make dataset folder
    if not os.path.exists(output_path+'/'+archive_name+'/'+data_name):
        os.mkdir(output_path+'/'+archive_name+'/'+data_name)

    #Generate datasets and folders
    for n in numberofattributes:
        for s in samplesize:
            for h in heritability:
                for m in minorAF:
                    modelName = "L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_K_'+str(K)
                    modelFile = model_dest+'/'+modelName+"_Models.txt"
                    genDataName = output_path+'/'+archive_name+'/'+data_name+'/'+data_name+'_A_'+str(n)+'_S_'+str(s)+'_'+str(modelName)
                    #Create gametes run command
                    filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar -i '+modelFile+' -D "-n '+str(AF_Min)+' -x '+str(AF_Max)+' -a '+str(n)+' -s '+str(int(s/2))+' -w '+str(int(s/2))+' -r '+str(replicates)+' -o '+str(genDataName)+'"'

                    if run_parallel:
                        job_ref = str(time.time())
                        job_path_name = job_dest+'/gametes_'+data_name+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'_run.sh'
                        sh_file = open(job_path_name,'w')
                        sh_file.write('#!/bin/bash\n')
                        sh_file.write('#BSUB -q i2c2_normal'+'\n')
                        sh_file.write('#BSUB -J '+job_ref+'\n')
                        sh_file.write('#BSUB -R "rusage[mem=4G]"'+'\n')
                        sh_file.write('#BSUB -M 15GB'+'\n')
                        sh_file.write('#BSUB -o ' + log_dest+'/gametes_'+data_name+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.o\n')
                        sh_file.write('#BSUB -e ' + log_dest+'/gametes_'+data_name+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.e\n')
                        sh_file.write(filewrite)
                        sh_file.close()
                        os.system('bsub < '+job_path_name)
                        pass
                    else:
                        os.system(filewrite)

def epistasis_2_locus_hetero_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path):
    #Model parameters needed
    locus = 2
    heritability = [0.4]
    minorAF = [0.2]
    K = 0.3  #population prevelance
    #Define dataset parameters
    data_name = 'gametes_2way_epi_2het'
    samplesize = [1600] #[200, 400, 800, 1600, 3200, 6400] #assumes balanced datasets (#cases = #controls)
    numberofattributes = [100] # [20, 100, 1000, 10000, 100000] #[200, 100]
    weight = [50,75]  # have to do the math for both ratio, X and 100-X
    AF_Min = 0.01
    AF_Max = 0.5
    replicates = 30 #100

    #Make dataset folder
    if not os.path.exists(output_path+'/'+archive_name+'/'+data_name):
        os.mkdir(output_path+'/'+archive_name+'/'+data_name)

    #Generate datasets and folders
    for n in numberofattributes:
        for s in samplesize:
            for h in heritability:
                for m in minorAF:
                    modelName = "L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_K_'+str(K)
                    modelFile = model_dest+'/'+modelName+"_Models.txt"
                    for w in weight:
                        genDataName = output_path+'/'+archive_name+'/'+data_name+'/'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_'+str(modelName)
                        #Create gametes run command
                        filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar -i '+modelFile+' -w '+str(w)+' -i '+modelFile+' -w '+str(100-w)+' -D "-h heterogeneous -b -n '+str(AF_Min)+' -x '+str(AF_Max)+' -a '+str(n)+' -s '+str(int(s/2))+' -w '+str(int(s/2))+' -r '+str(replicates)+' -o '+str(genDataName)+'"'

                        if run_parallel:
                            job_ref = str(time.time())
                            job_path_name = job_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'_run.sh'
                            sh_file = open(job_path_name,'w')
                            sh_file.write('#!/bin/bash\n')
                            sh_file.write('#BSUB -q i2c2_normal'+'\n')
                            sh_file.write('#BSUB -J '+job_ref+'\n')
                            sh_file.write('#BSUB -R "rusage[mem=4G]"'+'\n')
                            sh_file.write('#BSUB -M 15GB'+'\n')
                            sh_file.write('#BSUB -o ' + log_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.o\n')
                            sh_file.write('#BSUB -e ' + log_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.e\n')
                            sh_file.write(filewrite)
                            sh_file.close()
                            os.system('bsub < '+job_path_name)
                            pass
                        else:
                            os.system(filewrite)

def epistasis_2_locus_additive_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path):
    #Model parameters needed
    locus = 2
    heritability = [0.4]
    minorAF = [0.2]
    K = 0.3  #population prevelance
    #Define dataset parameters
    data_name = 'gametes_2way_epi_2add'
    samplesize = [1600] #[200, 400, 800, 1600, 3200, 6400] #assumes balanced datasets (#cases = #controls)
    numberofattributes = [100] # [20, 100, 1000, 10000, 100000] #[200, 100]
    weight = [50,75]  # have to do the math for both ratio, X and 100-X
    AF_Min = 0.01
    AF_Max = 0.5
    replicates = 30 #100

    #Make dataset folder
    if not os.path.exists(output_path+'/'+archive_name+'/'+data_name):
        os.mkdir(output_path+'/'+archive_name+'/'+data_name)

    #Generate datasets and folders
    for n in numberofattributes:
        for s in samplesize:
            for h in heritability:
                for m in minorAF:
                    modelName = "L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_K_'+str(K)
                    modelFile = model_dest+'/'+modelName+"_Models.txt"
                    for w in weight:
                        genDataName = output_path+'/'+archive_name+'/'+data_name+'/'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_'+str(modelName)
                        #Create gametes run command
                        filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar -i '+modelFile+' -w '+str(w)+' -i '+modelFile+' -w '+str(100-w)+' -D "-h hierarchical -n '+str(AF_Min)+' -x '+str(AF_Max)+' -a '+str(n)+' -s '+str(int(s/2))+' -w '+str(int(s/2))+' -r '+str(replicates)+' -o '+str(genDataName)+'"'

                        if run_parallel:
                            job_ref = str(time.time())
                            job_path_name = job_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'_run.sh'
                            sh_file = open(job_path_name,'w')
                            sh_file.write('#!/bin/bash\n')
                            sh_file.write('#BSUB -q i2c2_normal'+'\n')
                            sh_file.write('#BSUB -J '+job_ref+'\n')
                            sh_file.write('#BSUB -R "rusage[mem=4G]"'+'\n')
                            sh_file.write('#BSUB -M 15GB'+'\n')
                            sh_file.write('#BSUB -o ' + log_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.o\n')
                            sh_file.write('#BSUB -e ' + log_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.e\n')
                            sh_file.write(filewrite)
                            sh_file.close()
                            os.system('bsub < '+job_path_name)
                            pass
                        else:
                            os.system(filewrite)

def univariate_2_locus_hetero_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path):
    #Model parameters needed
    locus = 1
    heritability = [0.4]
    minorAF = [0.2]
    K = 0.3  #population prevelance
    #Define dataset parameters
    data_name = 'gametes_uni_2het'
    samplesize = [1600] #[200, 400, 800, 1600, 3200, 6400] #assumes balanced datasets (#cases = #controls)
    numberofattributes = [100] # [20, 100, 1000, 10000, 100000] #[200, 100]
    weight = [50]  # have to do the math for both ratio, X and 100-X
    AF_Min = 0.01
    AF_Max = 0.5
    replicates = 30 #100

    #Make dataset folder
    if not os.path.exists(output_path+'/'+archive_name+'/'+data_name):
        os.mkdir(output_path+'/'+archive_name+'/'+data_name)

    #Generate datasets and folders
    for n in numberofattributes:
        for s in samplesize:
            for h in heritability:
                for m in minorAF:
                    modelName = "L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_K_'+str(K)
                    modelFile = model_dest+'/'+modelName+"_Models.txt"
                    for w in weight:
                        genDataName = output_path+'/'+archive_name+'/'+data_name+'/'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_'+str(modelName)
                        #Create gametes run command
                        filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar -i '+modelFile+' -w '+str(w)+' -i '+modelFile+' -w '+str(w)+' -D "-h heterogeneous -b -n '+str(AF_Min)+' -x '+str(AF_Max)+' -a '+str(n)+' -s '+str(int(s/2))+' -w '+str(int(s/2))+' -r '+str(replicates)+' -o '+str(genDataName)+'"'

                        if run_parallel:
                            job_ref = str(time.time())
                            job_path_name = job_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'_run.sh'
                            sh_file = open(job_path_name,'w')
                            sh_file.write('#!/bin/bash\n')
                            sh_file.write('#BSUB -q i2c2_normal'+'\n')
                            sh_file.write('#BSUB -J '+job_ref+'\n')
                            sh_file.write('#BSUB -R "rusage[mem=4G]"'+'\n')
                            sh_file.write('#BSUB -M 15GB'+'\n')
                            sh_file.write('#BSUB -o ' + log_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.o\n')
                            sh_file.write('#BSUB -e ' + log_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.e\n')
                            sh_file.write(filewrite)
                            sh_file.close()
                            os.system('bsub < '+job_path_name)
                            pass
                        else:
                            os.system(filewrite)

def univariate_2_locus_additive_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path):
    #Model parameters needed
    locus = 1
    heritability = [0.4]
    minorAF = [0.2]
    K = 0.3  #population prevelance
    #Define dataset parameters
    data_name = 'gametes_uni_2add'
    samplesize = [1600] #[200, 400, 800, 1600, 3200, 6400] #assumes balanced datasets (#cases = #controls)
    numberofattributes = [100] # [20, 100, 1000, 10000, 100000] #[200, 100]
    weight = [50]  # have to do the math for both ratio, X and 100-X
    AF_Min = 0.01
    AF_Max = 0.5
    replicates = 30 #100

    #Make dataset folder
    if not os.path.exists(output_path+'/'+archive_name+'/'+data_name):
        os.mkdir(output_path+'/'+archive_name+'/'+data_name)

    #Generate datasets and folders
    for n in numberofattributes:
        for s in samplesize:
            for h in heritability:
                for m in minorAF:
                    modelName = "L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_K_'+str(K)
                    modelFile = model_dest+'/'+modelName+"_Models.txt"
                    for w in weight:
                        genDataName = output_path+'/'+archive_name+'/'+data_name+'/'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_'+str(modelName)
                        #Create gametes run command
                        filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar -i '+modelFile+' -w '+str(w)+' -i '+modelFile+' -w '+str(w)+' -D "-h hierarchical -b -n '+str(AF_Min)+' -x '+str(AF_Max)+' -a '+str(n)+' -s '+str(int(s/2))+' -w '+str(int(s/2))+' -r '+str(replicates)+' -o '+str(genDataName)+'"'

                        if run_parallel:
                            job_ref = str(time.time())
                            job_path_name = job_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'_run.sh'
                            sh_file = open(job_path_name,'w')
                            sh_file.write('#!/bin/bash\n')
                            sh_file.write('#BSUB -q i2c2_normal'+'\n')
                            sh_file.write('#BSUB -J '+job_ref+'\n')
                            sh_file.write('#BSUB -R "rusage[mem=4G]"'+'\n')
                            sh_file.write('#BSUB -M 15GB'+'\n')
                            sh_file.write('#BSUB -o ' + log_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.o\n')
                            sh_file.write('#BSUB -e ' + log_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.e\n')
                            sh_file.write(filewrite)
                            sh_file.close()
                            os.system('bsub < '+job_path_name)
                            pass
                        else:
                            os.system(filewrite)

def univariate_4_locus_hetero_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path):
    #Model parameters needed
    locus = 1
    heritability = [0.4]
    minorAF = [0.2]
    K = 0.3  #population prevelance
    #Define dataset parameters
    data_name = 'gametes_uni_4het'
    samplesize = [1600] #[200, 400, 800, 1600, 3200, 6400] #assumes balanced datasets (#cases = #controls)
    numberofattributes = [100] # [20, 100, 1000, 10000, 100000] #[200, 100]
    weight = [50]  # have to do the math for both ratio, X and 100-X
    AF_Min = 0.01
    AF_Max = 0.5
    replicates = 30 #100

    #Make dataset folder
    if not os.path.exists(output_path+'/'+archive_name+'/'+data_name):
        os.mkdir(output_path+'/'+archive_name+'/'+data_name)

    #Generate datasets and folders
    for n in numberofattributes:
        for s in samplesize:
            for h in heritability:
                for m in minorAF:
                    modelName = "L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_K_'+str(K)
                    modelFile = model_dest+'/'+modelName+"_Models.txt"
                    for w in weight:
                        genDataName = output_path+'/'+archive_name+'/'+data_name+'/'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_'+str(modelName)
                        #Create gametes run command
                        filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar -i '+modelFile+' -w '+str(w)+' -i '+modelFile+' -w '+str(w)+' -i '+modelFile+' -w '+str(w)+' -i '+modelFile+' -w '+str(w)+' -D "-h heterogeneous -b -n '+str(AF_Min)+' -x '+str(AF_Max)+' -a '+str(n)+' -s '+str(int(s/2))+' -w '+str(int(s/2))+' -r '+str(replicates)+' -o '+str(genDataName)+'"'

                        if run_parallel:
                            job_ref = str(time.time())
                            job_path_name = job_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'_run.sh'
                            sh_file = open(job_path_name,'w')
                            sh_file.write('#!/bin/bash\n')
                            sh_file.write('#BSUB -q i2c2_normal'+'\n')
                            sh_file.write('#BSUB -J '+job_ref+'\n')
                            sh_file.write('#BSUB -R "rusage[mem=4G]"'+'\n')
                            sh_file.write('#BSUB -M 15GB'+'\n')
                            sh_file.write('#BSUB -o ' + log_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.o\n')
                            sh_file.write('#BSUB -e ' + log_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.e\n')
                            sh_file.write(filewrite)
                            sh_file.close()
                            os.system('bsub < '+job_path_name)
                            pass
                        else:
                            os.system(filewrite)

def univariate_4_locus_additive_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path):
    #Model parameters needed
    locus = 1
    heritability = [0.4]
    minorAF = [0.2]
    K = 0.3  #population prevelance
    #Define dataset parameters
    data_name = 'gametes_uni_4add'
    samplesize = [1600] #[200, 400, 800, 1600, 3200, 6400] #assumes balanced datasets (#cases = #controls)
    numberofattributes = [100] # [20, 100, 1000, 10000, 100000] #[200, 100]
    weight = [50]  # have to do the math for both ratio, X and 100-X
    AF_Min = 0.01
    AF_Max = 0.5
    replicates = 30 #100

    #Make dataset folder
    if not os.path.exists(output_path+'/'+archive_name+'/'+data_name):
        os.mkdir(output_path+'/'+archive_name+'/'+data_name)

    #Generate datasets and folders
    for n in numberofattributes:
        for s in samplesize:
            for h in heritability:
                for m in minorAF:
                    modelName = "L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_K_'+str(K)
                    modelFile = model_dest+'/'+modelName+"_Models.txt"
                    for w in weight:
                        genDataName = output_path+'/'+archive_name+'/'+data_name+'/'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_'+str(modelName)
                        #Create gametes run command
                        filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar -i '+modelFile+' -w '+str(w)+' -i '+modelFile+' -w '+str(w)+' -i '+modelFile+' -w '+str(w)+' -i '+modelFile+' -w '+str(w)+' -D "-h hierarchical -b -n '+str(AF_Min)+' -x '+str(AF_Max)+' -a '+str(n)+' -s '+str(int(s/2))+' -w '+str(int(s/2))+' -r '+str(replicates)+' -o '+str(genDataName)+'"'

                        if run_parallel:
                            job_ref = str(time.time())
                            job_path_name = job_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'_run.sh'
                            sh_file = open(job_path_name,'w')
                            sh_file.write('#!/bin/bash\n')
                            sh_file.write('#BSUB -q i2c2_normal'+'\n')
                            sh_file.write('#BSUB -J '+job_ref+'\n')
                            sh_file.write('#BSUB -R "rusage[mem=4G]"'+'\n')
                            sh_file.write('#BSUB -M 15GB'+'\n')
                            sh_file.write('#BSUB -o ' + log_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.o\n')
                            sh_file.write('#BSUB -e ' + log_dest+'/gametes_'+data_name+'_W_'+str(w)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.e\n')
                            sh_file.write(filewrite)
                            sh_file.close()
                            os.system('bsub < '+job_path_name)
                            pass
                        else:
                            os.system(filewrite)

def epistasis_2_locus_imbalanced_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path):
    #Model parameters needed
    locus = 2
    heritability = [0.4]
    minorAF = [0.2]
    K = 0.3  #population prevelance
    #Define dataset parameters
    data_name = 'gametes_2way_epistasis_inbal'
    samplesize = [1600] #[200, 400, 800, 1600, 3200, 6400] #assumes balanced datasets (#cases = #controls)
    numberofattributes = [100] # [20, 100, 1000, 10000, 100000] #[200, 100]
    balance = [.6,.9]
    AF_Min = 0.01
    AF_Max = 0.5
    replicates = 30 #100

    #Make dataset folder
    if not os.path.exists(output_path+'/'+archive_name+'/'+data_name):
        os.mkdir(output_path+'/'+archive_name+'/'+data_name)

    #Generate datasets and folders
    for n in numberofattributes:
        for s in samplesize:
            for h in heritability:
                for m in minorAF:
                    for b in balance:
                        #Calculate case and control counts
                        controlCount = int(float(s)*b)
                        caseCount = int(s-controlCount)

                        modelName = "L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_K_'+str(K)
                        modelFile = model_dest+'/'+modelName+"_Models.txt"
                        genDataName = output_path+'/'+archive_name+'/'+data_name+'/'+data_name+'_'+str(b)+'_A_'+str(n)+'_S_'+str(s)+'_'+str(modelName)
                        #Create gametes run command
                        filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar -i '+modelFile+' -D "-n '+str(AF_Min)+' -x '+str(AF_Max)+' -a '+str(n)+' -s '+str(caseCount)+' -w '+str(controlCount)+' -r '+str(replicates)+' -o '+str(genDataName)+'"'

                        if run_parallel:
                            job_ref = str(time.time())
                            job_path_name = job_dest+'/gametes_'+data_name+'_'+str(b)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'_run.sh'
                            sh_file = open(job_path_name,'w')
                            sh_file.write('#!/bin/bash\n')
                            sh_file.write('#BSUB -q i2c2_normal'+'\n')
                            sh_file.write('#BSUB -J '+job_ref+'\n')
                            sh_file.write('#BSUB -R "rusage[mem=4G]"'+'\n')
                            sh_file.write('#BSUB -M 15GB'+'\n')
                            sh_file.write('#BSUB -o ' + log_dest+'/gametes_'+data_name+'_'+str(b)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.o\n')
                            sh_file.write('#BSUB -e ' + log_dest+'/gametes_'+data_name+'_'+str(b)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.e\n')
                            sh_file.write(filewrite)
                            sh_file.close()
                            os.system('bsub < '+job_path_name)
                            pass
                        else:
                            os.system(filewrite)

def epistasis_2_locus_quantitative_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path):
    #Model parameters needed
    locus = 2
    heritability = [0.4]
    minorAF = [0.2]
    K = 0.3  #population prevelance
    #Define dataset parameters
    data_name = 'gametes_2way_epistasis_quant'
    samplesize = [1600] #[200, 400, 800, 1600, 3200, 6400] #assumes balanced datasets (#cases = #controls)
    numberofattributes = [100] # [20, 100, 1000, 10000, 100000] #[200, 100]
    standardDev = [0.2,0.5,0.8]
    AF_Min = 0.01
    AF_Max = 0.5
    replicates = 30 #100

    #Make dataset folder
    if not os.path.exists(output_path+'/'+archive_name+'/'+data_name):
        os.mkdir(output_path+'/'+archive_name+'/'+data_name)

    #Generate datasets and folders
    for n in numberofattributes:
        for s in samplesize:
            for h in heritability:
                for m in minorAF:
                    for d in standardDev:
                        modelName = "L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_K_'+str(K)
                        modelFile = model_dest+'/'+modelName+"_Models.txt"
                        genDataName = output_path+'/'+archive_name+'/'+data_name+'/'+data_name+'_'+str(d)+'_A_'+str(n)+'_S_'+str(s)+'_'+str(modelName)
                        #Create gametes run command
                        filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar -i '+modelFile+' -D "-c -d '+ str(d) + ' -t '+ str(s) + ' -n '+str(AF_Min)+' -x '+str(AF_Max)+' -a '+str(n)+' -r '+str(replicates)+' -o '+str(genDataName)+'"'

                        if run_parallel:
                            job_ref = str(time.time())
                            job_path_name = job_dest+'/gametes_'+data_name+'_'+str(d)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'_run.sh'
                            sh_file = open(job_path_name,'w')
                            sh_file.write('#!/bin/bash\n')
                            sh_file.write('#BSUB -q i2c2_normal'+'\n')
                            sh_file.write('#BSUB -J '+job_ref+'\n')
                            sh_file.write('#BSUB -R "rusage[mem=4G]"'+'\n')
                            sh_file.write('#BSUB -M 15GB'+'\n')
                            sh_file.write('#BSUB -o ' + log_dest+'/gametes_'+data_name+'_'+str(d)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.o\n')
                            sh_file.write('#BSUB -e ' + log_dest+'/gametes_'+data_name+'_'+str(d)+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.e\n')
                            sh_file.write(filewrite)
                            sh_file.close()
                            os.system('bsub < '+job_path_name)
                            pass
                        else:
                            os.system(filewrite)

def epistasis_2_locus_numfeatures_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path):
    #Model parameters needed
    locus = 2
    heritability = [0.4]
    minorAF = [0.2]
    K = 0.3  #population prevelance
    #Define dataset parameters
    data_name = 'gametes_2way_epistasis_numfeat'
    samplesize = [1600] #[200, 400, 800, 1600, 3200, 6400] #assumes balanced datasets (#cases = #controls)
    numberofattributes = [1000,10000,100000] # [20, 100, 1000, 10000, 100000] #[200, 100]
    AF_Min = 0.01
    AF_Max = 0.5
    replicates = 30 #100

    #Make dataset folder
    if not os.path.exists(output_path+'/'+archive_name+'/'+data_name):
        os.mkdir(output_path+'/'+archive_name+'/'+data_name)

    #Generate datasets and folders
    for n in numberofattributes:
        for s in samplesize:
            for h in heritability:
                for m in minorAF:
                    modelName = "L_"+str(locus)+"_H_"+str(h)+"_F_"+str(m)+'_K_'+str(K)
                    modelFile = model_dest+'/'+modelName+"_Models.txt"
                    genDataName = output_path+'/'+archive_name+'/'+data_name+'/'+data_name+'_A_'+str(n)+'_S_'+str(s)+'_'+str(modelName)
                    #Create gametes run command
                    filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar -i '+modelFile+' -D "-n '+str(AF_Min)+' -x '+str(AF_Max)+' -a '+str(n)+' -s '+str(int(s/2))+' -w '+str(int(s/2))+' -r '+str(replicates)+' -o '+str(genDataName)+'"'

                    if run_parallel:
                        job_ref = str(time.time())
                        job_path_name = job_dest+'/gametes_'+data_name+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'_run.sh'
                        sh_file = open(job_path_name,'w')
                        sh_file.write('#!/bin/bash\n')
                        sh_file.write('#BSUB -q i2c2_normal'+'\n')
                        sh_file.write('#BSUB -J '+job_ref+'\n')
                        sh_file.write('#BSUB -R "rusage[mem=4G]"'+'\n')
                        sh_file.write('#BSUB -M 15GB'+'\n')
                        sh_file.write('#BSUB -o ' + log_dest+'/gametes_'+data_name+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.o\n')
                        sh_file.write('#BSUB -e ' + log_dest+'/gametes_'+data_name+'_A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.e\n')
                        sh_file.write(filewrite)
                        sh_file.close()
                        os.system('bsub < '+job_path_name)
                        pass
                    else:
                        os.system(filewrite)

######################################
if __name__ == '__main__':
    sys.exit(main(sys.argv))
