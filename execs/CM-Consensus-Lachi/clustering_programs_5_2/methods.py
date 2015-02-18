

import time
import os
import sys
import glob

clus_bin = os.path.dirname(os.path.realpath(__file__)) + '/bin/'
clus_bin_just = os.path.dirname(os.path.realpath(__file__)) + '/bin/'


def exec_oslom(program_number, seed, run_number, run_cvis=True):
	
	command_line=""
	if(program_number==0):
		command_line= clus_bin+"oslom_undir ";
	else: 
		command_line= clus_bin+"oslom_dir ";
	
	if run_number=='consensus':
		command_line+= " -f  network_file_consensus -uw -r 3 -hr 3 -seed "+seed+" > cluout "
	else:
		command_line+= " -f  network_file -uw -r 3 -hr 3 -seed "+seed+" > cluout "
			
	print 'running',  command_line
	os.system(command_line)

	if run_number=='consensus':

		os.system("rm -r network_file_oslo_files")		
		os.system("mv network_file_consensus_oslo_files network_file_oslo_files")

	
	if run_cvis:
		if program_number==0:
			command_line= clus_bin+"cvis_undir -f  network_file -folder network_file_oslo_files -labels label_file_copy -seed 1 > cvis_out"
		else:
			command_line= clus_bin+"cvis_dir -f  network_file -folder network_file_oslo_files -labels label_file_copy -seed 1 > cvis_out"
		
		os.system(command_line)
	
	
	os.system("mkdir results")
	os.system("cp network_file_oslo_files/tp* results/")
	os.system("mv results results_"+run_number)



def flat_partition(program_number, seed, run_number, run_cvis=True):
	
	
	netwfile="network_file"
	if run_number=="consensus":
		netwfile="network_file_consensus"
	
	if(program_number==2):
		command_line= clus_bin+"infomap_undir_script "+netwfile+" "+seed+" 10 "+ clus_bin_just +" > cluout"
		partfile="infomap.part"
	if(program_number==3):
		command_line= clus_bin+"infomap_dir_script "+netwfile+" "+seed+" 10 "+ clus_bin_just +" > cluout"
		partfile="infomap.part"
	if(program_number==5):
		command_line= clus_bin+"lpm -f "+netwfile+" -r 10 -seed "+seed+" > cluout"
		partfile=netwfile+"part"
	if(program_number==8):
		command_line= clus_bin+"modopt "+netwfile+" "+seed+" 1 5 0. 0.  > cluout"
		partfile=netwfile+"part"

	
	print 'running',  command_line

	os.system(command_line)
	os.system("mkdir results")
	os.system("mv "+partfile+" results/tp")
	
	
	if run_cvis:
		
		if program_number==3:
			command_line= clus_bin+"cvis_dir -f  network_file -folder results -labels label_file_copy -seed 1 > cvis_out"
		else:
			command_line= clus_bin+"cvis_undir -f  network_file -folder results -labels label_file_copy -seed 1 > cvis_out"
		
		os.system(command_line)
	
	os.system("mv results results_"+run_number)



def hiermap_and_louvain(program_number, seed, run_number):
	
	
	netwfile="network_file"
	if run_number=="consensus":
		netwfile="network_file_consensus"

	
	if(program_number==6):
		command_line= clus_bin+"infohiermap_undir_script "+netwfile+" "+seed+" 10 "+ clus_bin_just +" > cluout"
	if(program_number==7):
		command_line= clus_bin+"infohiermap_dir_script "+netwfile+" "+seed+" 10 "+ clus_bin_just +" > cluout"
	if(program_number==4):
		
		tp_files=glob.glob('tp*')
		for tf in tp_files:
			os.system('rm '+tf)

		stp_files=glob.glob('short_tp*')
		for tf in stp_files:
			os.system('rm '+tf)

		
		command_line= clus_bin+"louvain_method -f "+netwfile+" -r 10 -seed "+seed+" > cluout"
	
	print 'running',  command_line
	os.system(command_line)
	
	
	if(program_number==7):
		command_line= clus_bin+"cvis_dir -f network_file -tree infomap_net.tree -labels label_file_copy  -seed 1 > cvis_out"
	if(program_number==6): 
		command_line= clus_bin+"cvis_undir -f network_file -tree infomap_net.tree -labels label_file_copy  -seed 1 > cvis_out"
	if(program_number==4): 
		command_line= clus_bin+"cvis_undir -f  network_file -folder . -labels label_file_copy -seed 1 > cvis_out"
	
	
	os.system(command_line)
	
	os.system("mkdir results")
	
	pi=0
	while (True):
		
		file=""
		if program_number==6 or program_number==7:
			file= "network_file_cvis/original_pt_"+str(pi)+".dat"
		else:
			if pi==0:
				file="tp"
			else:
				file="tp"+str(pi)
		
		if(os.path.exists(file)==False):
			break;
		if(pi==0):
			os.system("cp "+file+" results/tp")
		else:
			os.system("cp "+file+" results/tp"+str(pi))
		pi+=1

	os.system("mv results results_"+run_number)



def check_int(num):
	try:
		s=int(num)
	except:
		print num, "is not an integer"
		return False
	
	if s <0:
		print num, "is negative"
		return False
	
	return True

def check_float(num):
	try:
		s=float(num)
	except:
		print num, "is not a float number"
		return False
	
	if s <0:
		print num, "is negative"
		return False
	
	return True


def check_netfile(netfile,new_netfile):
	
	fcleaned=open(new_netfile, "w")
	
	f=open(netfile)
	for h in f.readlines():
		l=h.split()
		
		skip=False
		
		if len(l)> 3:
			print "some rows have more than three numbers: [id] [id] [weight] is expected instead"
			skip=True
		if len(l)==1:
			print "some rows have more only one number: [id] [id] [weight] is expected instead"
			skip=True
		
		if len(l)>1:
			if check_int(l[0])==False:
				skip=True
			if check_int(l[1])==False:
				skip=True
		
		if len(l)==3:
			if check_float(l[2])==False:
				skip=True
		
		if skip==False:
			
			if len(l)==3:
				fcleaned.write(l[0]+' '+l[1]+' '+l[2]+'\n')
			if len(l)==2:
				fcleaned.write(l[0]+' '+l[1]+'\n')
				


	fcleaned.close()
	f.close()



