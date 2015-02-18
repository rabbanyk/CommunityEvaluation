
"""

/*
    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    *                                                                               *
    *    This program is free software; you can redistribute it and/or modify         *
    *  it under the terms of the GNU General Public License as published by         *
    *  the Free Software Foundation; either version 2 of the License, or            *
    *  (at your option) any later version.                                          *
    *                                                                               *
    *  This program is distributed in the hope that it will be useful,              *
    *  but WITHOUT ANY WARRANTY; without even the implied warranty of               *
    *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
    *  GNU General Public License for more details.                                 *
    *                                                                               *
    *  You should have received a copy of the GNU General Public License            *
    *  along with this program; if not, write to the Free Software                  *
    *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA    *
    *                                                                               *
    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    *                                                                               *
    *  Created by Andrea Lancichinetti  (email: arg.lanci@gmail.com)                *
    *                                                                               *
    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    */
    
    
"""


"""
    this module contains all the wrappers
    for the clustering programs and the layouts

"""


import os
import glob

from methods import *
import time


def print_statement():
    
    print "python select.py -n [network_file]  -p [program_number] -f [folder] -c [consensus_runs] -l [label_file] -t [threshold] -slices [slice_file]" 
    print "Example:\npython select.py -n datasets/wordnet.dat  -l datasets/wordlabels.dat -p 4 -f test_folder -c 1 \n"    
    print "Example:\npython select.py -n datasets/example.dat -p 2 -f test_folder -c 10 \n"
    print "Example:\npython select.py -slices datasets/slice_file -p 4 -f slice_test -c 5\n"
    print "*** [network_file] is the edge list network " 
    print "*** [folder]: a directory where all the results will be written. If [folder] already exists, the code will exit."
    print "*** [consensus_runs] is the number of consensus runs. Set it to 1 if you want to skip this step."  
    print "*** [program_number] is to select which algorithm you would like to run:"  
    print "   0: oslom undirected"  
    print "   1: oslom directed"  
    print "   2: infomap undirected"  
    print "   3: infomap directed"  
    print "   4: louvain"  
    print "   5: label propagation method"  
    print "   6: hierarchical infomap undirected"  
    print "   7: hierarchical infomap directed"
    print "   8: modularity optimization (simulated annealing zero temperature)"
    print "+++ [label_file] is in the format \"label id\" where 'label' is the name of the node 'id' (optional)."
    print "+++ [threshold] is the paramter for pruning the consensus matrix. Defaulf value is 0.5 (optional)"
    print "/// [slice_file] is a file containing different slices of the network (snapshots, for instance). If you use this, do not use option -n."
    exit()


def parse_argv(args):
    
    
    netfile=''
    folder_name=''
    program_number=-1
    
    runs=-1
    threshold=0.5
    label_file='NULL'
    slice_file=''
    
    
    for k,s in enumerate(args):
        
        if s=='-n':
            netfile=args[k+1]
        if s=='-l':
            label_file=args[k+1]
        if s=='-p':
            program_number=int(args[k+1])
        if s=='-f':
            folder_name=args[k+1]+'/'
        if s=='-c':
            runs=int(args[k+1])
        if s=='-t':
            threshold=float(args[k+1])
        if s=='-slices':
            slice_file=args[k+1]
                

    if program_number==-1:
        print 'option -p was not selected. Which algorithm would you like to run?'
        print_statement()

    if netfile=='' and slice_file=='':
        print 'option -n was not selected. You need to specify a network or a slice file'
        print_statement()
        
    if folder_name=='':
        print 'option -f was not selected. You need to specify a folder where to write the results.'
        print_statement()

    if runs==-1:
        print 'option -c was not selected. How many runs for the running the consensus? Use \'-c 1\' if you want to skip this.'
        print_statement()


                
    return (netfile, label_file, folder_name, program_number, runs, threshold, slice_file)




def run_method_given_folder_and_network(folder_name, program_number, seed, run_number):
        
    current_dir = os.path.dirname(os.path.realpath(__file__))
    os.chdir(folder_name)
        
    if program_number==0 or program_number==1:
        exec_oslom(program_number, seed, run_number)
    if program_number==2 or program_number==3 or program_number==5 or program_number==8:
        flat_partition(program_number, seed, run_number)
    if program_number==4 or program_number==6 or program_number==7:
        hiermap_and_louvain(program_number, seed, run_number)
    
    os.system("mv network_file_cvis/cvis_*.gdf results_"+str(run_number))
    os.system("rm -r network_file_cvis")

    os.chdir(current_dir)



def make_program_number_undirected(pn):

    if pn==1:
        pn=0
    if pn==3:
        pn=2
    if pn==7:
        pn=6
    
    return pn


def run_consensus_method(glob_string, folder_name, threshold, program_number):
    
    pts=glob.glob(glob_string)
    f=open(folder_name+'/pts_file', 'w')
    for gf in pts:
        f.write(gf+'\n')
    f.close()
        
    command_line=os.path.dirname(os.path.realpath(__file__)) + '/bin/pr_consensus '+folder_name+'/pts_file '+str(threshold)
    print command_line
    os.system(command_line)
        
    os.system('mv '+folder_name+'/pts_file_net.dat '+folder_name+'/network_file_consensus')
    run_method_given_folder_and_network(folder_name, make_program_number_undirected(program_number), str(10), 'consensus')


def run_all(runs, folder_name, program_number, netfile, label_file, threshold, run_consensus):
    
    print "input:", '\nnetwork:', netfile, '\nlabel file: ', label_file, '\nfolder: ', folder_name, '\nprogram: ', program_number, '\nconsensus runs:', runs, '\nthreshold:', threshold
    
    if(os.path.exists(folder_name)):
        print folder_name +" already exists. exiting..."
        exit()
    
    #os.system("mkdir " +  folder_name)
    os.mkdir(folder_name)
    print 'checking file format...'
    check_netfile(netfile, folder_name+"network_file")
    
    if(os.path.exists(label_file)):
        os.system("cp "+label_file+" "+folder_name+"label_file_copy")
    
    print 'NB. do not worry if the program tries to delete files which do not exist: it\'s just a double check.'
    print 'To check what the program is doing, please open the log files '+folder_name+'cluot and '+folder_name+'cvis_out'
    time.sleep(2)
    
    for k in range(runs):
        run_method_given_folder_and_network(folder_name, program_number, str(k+1), str(k+1))
        print 'run ', k+1, 'done!'
        
    if runs>1 and run_consensus:
        print 'running consensus clustering'
        run_consensus_method(folder_name+'results_*/tp', folder_name, threshold, program_number)






if __name__ == '__main__':
    
    
    
    
    
    netfile, label_file, folder_name, program_number, runs, threshold, slice_file = parse_argv(sys.argv)
        
    if slice_file=='':
        run_all(runs, folder_name, program_number, netfile, label_file, threshold, True)
    else:
        f=open(slice_file)
        files=[]
        for l in f.readlines():
            files+=l.split()
        f.close()

        print 'slices found:', files

        if(os.path.exists(folder_name)):
            print folder_name +" already exists. exiting..."
            exit()
        #os.system("mkdir " +  folder_name)
        os.mkdir(folder_name)
        
        for k,f in enumerate(files):
            run_all(runs, folder_name+'/slice_'+str(k)+'/', program_number, f, label_file, threshold, False)
            
        try:
            os.system('cp '+folder_name+'/slice_'+str(len(files)-1)+'/network_file '+folder_name+'/network_file')
        except:
            pass

                
        print 'running consensus clustering for the slices'
        run_consensus_method(folder_name+'/slice_*/results_*/tp', folder_name+'/', threshold, program_number)

        try:
            print 'the network used for the visualization in ', folder_name+'results_consensus is ', files[-1]
            print 'which is the last in file ', slice_file
        except:
            pass
    
            
    print 'Done!\nIf you are planning to open the cvis files with gephi, you might want to use the expansion layout a little bit.'













        
        
