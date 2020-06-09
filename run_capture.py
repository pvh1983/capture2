import os 


'''
Run capture analyes

'''

#
os_sys = 'windows' # 'windows' or 'linux'
ncores = 2

# get curret directory
#proj_dir = 'f:/projects/capture/codes/' # Windows
#proj_dir = '/home/hpham/Documents/projects/capture/codes/' # Ubuntu
proj_dir = os.getenv('proj_dir') + '/'

model_files_prefix = 'Capture_May2020'
print(f'Project directory is: {proj_dir}')

# Create a temp direction
temp_dir = proj_dir + 'tmp'
if not os.path.exists(temp_dir):
    os.makedirs(temp_dir)

# Load list of pumping locatins
ifile_pmploc = '../input/rnum_test.dat'
count = 0
with open(ifile_pmploc) as f:
    for line in f:
        #print(f'Reading line: {line}')
        val = line.split()
        rid,lay,row,col = int(val[0]),int(val[1]),int(val[2]),int(val[3])
        #print(f'{rid,lay,row,col}')
        
        # Create a new to folder for a single model run
        run_dir = proj_dir + 'tmp/run_' + str(rid).rjust(6, '0')
        if not os.path.exists(run_dir):            
            os.mkdir(run_dir)
        
        # Go to the new created folder
        os.chdir(run_dir)
	
        # write the current run loc to file
        with open('pmploc.dat', 'w') as out:
            out.write(f'{rid} {lay} {row} {col}')

        # Link model files        
        cmd = 'ln -s ' + proj_dir + 'model_final/' + model_files_prefix + '.* .' 
        #print(f'{cmd}')
        os.system(cmd)

        cmd = 'ln -s ' + proj_dir + 'model_final/' + 'f*.inp .' 
        os.system(cmd)	

        cmd = 'ln -s ' + proj_dir + 'model_final/' + 'mfnwt .' 
        os.system(cmd)	

        cmd = 'ln -s ' + proj_dir + 'scripts/' + 'main2.py .' 
        os.system(cmd)	

        cmd = 'rm -f *.out *.hed *.ccf *._os *.wel *.drw'               
        os.system(cmd)

        # run capture
        os.system('python main2.py &')

        # Change to the script folder
        os.chdir(proj_dir + 'scripts')
        count+=1
        
        # If all cores are filled, wait ...
        if count % ncores == 0:
            count = 0
            os.system('wait')



