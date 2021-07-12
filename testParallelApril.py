
from multiprocessing import Process
import os
import time
# git remote set-url origin https://mgrecu35@github.com/mgrecu35/cmbv7.git

def info(title):
    print(title)
    print('module name:', __name__)
    print('parent process:', os.getppid())
    print('process id:', os.getpid())
    
def fsh(fname_p):
    fname=fname_p[0]
    p=fname_p[1]
    lines=open(fname,"r").readlines()
    print(lines[-3])
    logName="log."+lines[-5][9:]
    cmd='bpsh %i ./combAlg2.exe junk %s >&out/out.%s'%(189-p,fname,logName)
    print(cmd)
    os.system(cmd)
    #time.sleep(1)
import glob
t11=time.time()
if __name__ == '__main__':
    fs=glob.glob("ite755ParamApril2018/*")
    fs=sorted(fs)
    nf=len(fs)
    n1=int(len(fs)/5)
    n1=0
    for i in range(0,n1+1):#range(13,n1+1):#,n1+2):
        jobs=[]
        t1=time.time()
        for k in range(5):
            if 10*i+k<nf:
                p = Process(target=fsh, args=((fs[10*i+k],k),))
                jobs.append(p)
                p.start()
        for j in jobs:
            j.join()
        print('all done')
        print(time.time()-t1)
    print(time.time()-t11)
