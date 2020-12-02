"""
Author: Ryan Urbanowicz
Created: 11/30/20
Description: Script to apply GAMETES to generate and organize a variety of SNP 2-way interaction models and a corresponding Datasets
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
        #Generate core epistasis models
        epistasis_2_locus_core_model(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path)

    elif use == 'data':
        #Generate core epistasis data
        epistasis_2_locus_core_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path)

    else:
        print("GAMETES use not recognized.")

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
            model_path_name = model_dest+"/H_"+str(h)+"_F_"+str(m)
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
                job_path_name = job_dest+'/gametes_'+"H_"+str(h)+"_F_"+str(m)+'_'+job_ref+'_run.sh'
                sh_file = open(job_path_name,'w')
                sh_file.write('#!/bin/bash\n')
                sh_file.write('#BSUB -q i2c2_normal'+'\n')
                sh_file.write('#BSUB -J '+job_ref+'\n')
                sh_file.write('#BSUB -R "rusage[mem=4G]"'+'\n')
                sh_file.write('#BSUB -M 15GB'+'\n')
                sh_file.write('#BSUB -o ' + log_dest+'/gametes_'+"H_"+str(h)+"_F_"+str(m)+'_'+job_ref+'.o\n')
                sh_file.write('#BSUB -e ' + log_dest+'/gametes_'+"H_"+str(h)+"_F_"+str(m)+'_'+job_ref+'.e\n')
                sh_file.write(filewrite)
                sh_file.close()
                os.system('bsub < '+job_path_name)
                pass
            else:
                os.system(filewrite)

def epistasis_2_locus_core_data(output_path,archive_name,model_dest,job_dest,log_dest,run_parallel,this_file_path):
    #Model parameters needed
    heritability = [0.05, 0.1, 0.2, 0.4]
    minorAF = [0.2]
    K = 0.3  #population prevelance
    #Define dataset parameters
    data_name = 'gametes_2way_epistasis'
    samplesize = [200, 400, 800, 1600] #[200, 400, 800, 1600, 3200, 6400] #assumes balanced datasets (#cases = #controls)
    numberofattributes = [20,100] # [20, 100, 1000, 10000, 100000] #[200, 100]
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
                    modelName = "H_"+str(h)+"_F_"+str(m)+'_K_'+str(K)
                    modelFile = model_dest+'/'+modelName+"_Models.txt"
                    genDataName = output_path+'/'+archive_name+'/'+data_name+'/'+data_name+'_A_'+str(n)+'_S_'+str(s)+'_'+str(modelName)
                    #Create gametes run command
                    filewrite = 'java -jar '+this_file_path+'/gametes_2.2_dev.jar -i '+modelFile+' -D "-n '+str(AF_Min)+' -x '+str(AF_Max)+' -a '+str(n)+' -s '+str(int(s/2))+' -w '+str(int(s/2))+' -r '+str(replicates)+' -o '+str(genDataName)+'"'

                    if run_parallel:
                        job_ref = str(time.time())
                        job_path_name = job_dest+'/gametes_'+'A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'_run.sh'
                        sh_file = open(job_path_name,'w')
                        sh_file.write('#!/bin/bash\n')
                        sh_file.write('#BSUB -q i2c2_normal'+'\n')
                        sh_file.write('#BSUB -J '+job_ref+'\n')
                        sh_file.write('#BSUB -R "rusage[mem=4G]"'+'\n')
                        sh_file.write('#BSUB -M 15GB'+'\n')
                        sh_file.write('#BSUB -o ' + log_dest+'/gametes_'+'A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.o\n')
                        sh_file.write('#BSUB -e ' + log_dest+'/gametes_'+'A_'+str(n)+'_S_'+str(s)+'_H_'+str(h)+'_F_'+str(m)+'_'+job_ref+'.e\n')
                        sh_file.write(filewrite)
                        sh_file.close()
                        os.system('bsub < '+job_path_name)
                        pass
                    else:
                        os.system(filewrite)

######################################
if __name__ == '__main__':
    sys.exit(main(sys.argv))
