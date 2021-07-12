vars=open('oeVars.txt','r').readlines()
#call copylognws1_fs( log10NwMean,dPRRet%n9(:,i,j),i-1)
#call copyenvsfqvs1t(dPRData%envQv(:,i,j),i-1)

fout=open('src_c/oe_output_kugmi.c','w')
foutF=open('oe_output.f90','w')
lines=open('header.txt','r').readlines()
for l in lines:
    fout.write(l)

def write_3(vname,sname,vtype,swath_type,n1):
    if 'integer' in vtype:
        argType='int *'
    if 'real' in vtype:
        argType="float *"
    sL=[]
    s='void copy_%s_%s_(%s%s, int *i)\n'%(vname.lower(),swath_type.lower(),argType,vname)
    sL.append(s)
    sL.append('{');
    sL.append('int k;\n');
    sL.append('extern L2BCMB_SWATHS swathx;\n');
    sL.append('for(k=0;k<%s;k++)\n'%n1)
    sL.append('{\n')
    sL.append('swathx.%s.OptEst.OE%s[*i][k]=%s[k];\n'%(swath_type,sname,vname))
    sL.append('}\n')
    sL.append('}\n\n')
    for s in sL:
        print(s)
        fout.write(s)

def write_3_call(vname,sname,vtype,swath_type):
    if 'integer' in vtype:
        argType='int *'
    if 'real' in vtype:
        argType="float *"
    sL=[]
    s='call copy_%s_%s(dPRData%%%s(:,i,j), i-1)\n'%(vname.lower(),swath_type.lower(),vname)
    foutF.write(s)

def write_2(vname,sname,vtype,swath_type):
    if 'integer' in vtype:
        argType='int *'
    if 'real' in vtype:
        argType="float *"
    sL=[]
    s='void copy_%s_%s_(%s%s, int *i)\n'%(vname.lower(),swath_type.lower(),argType,vname)
    sL.append(s)
    sL.append('{');
    sL.append('int k;\n');
    sL.append('extern L2BCMB_SWATHS swathx;\n');
    #sL.append('for(k=0;k<%s;k++)\n'%n1)
    #sL.append('{\n')
    sL.append('swathx.%s.OptEst.OE%s[*i]=*%s;\n'%(swath_type,sname,vname))
    #sL.append('}\n')
    sL.append('}\n\n')
    for s in sL:
        print(s)
        fout.write(s)

def write_2_call(vname,sname,vtype,swath_type):
    sL=[]
    s='call copy_%s_%s(dPRData%%%s(i,j), i-1)\n'%(vname.lower(),swath_type.lower(),vname)
    print(s)
    foutF.write(s)

def write_n9(vname,sname,vtype,swath_type):
    if 'integer' in vtype:
        argType='int *'
    if 'real' in vtype:
        argType="float *"
    sL=[]
    s='void copy_%s_%s_(%s%s,int *n9,int *i)\n'%(vname.lower(),swath_type.lower(),argType,vname)
    sL.append(s)
    sL.append('{');
    sL.append('int k;\n');
    sL.append('extern L2BCMB_SWATHS swathx;\n');
    sL.append('for(k=0;k<9;k++)\n')
    sL.append('{\n')
    sL.append('swathx.%s.OptEst.OE%s[*i][n9[k]]=%s[k];\n'%(swath_type,sname,vname))
    sL.append('}\n')
    sL.append('}\n\n')
    for s in sL:
        print(s)
        fout.write(s)

def write_n10(vname,sname,vtype,swath_type):
    if 'integer' in vtype:
        argType='int *'
    if 'real' in vtype:
        argType="float *"
    sL=[]
    s='void copy_%s_%s_(%s%s,short *n10,int *i)\n'%(vname.lower(),swath_type.lower(),argType,vname)
    sL.append(s)
    sL.append('{');
    sL.append('int k;\n');
    sL.append('extern L2BCMB_SWATHS swathx;\n');
    sL.append('for(k=0;k<10;k++)\n')
    sL.append('{\n')
    sL.append('swathx.%s.OptEst.OE%s[*i][k]=%s[n10[k]-1];\n'%(swath_type,sname,vname))
    sL.append('}\n')
    sL.append('}\n\n')
    for s in sL:
        print(s)
        fout.write(s)
        
def write_9_call(vname,sname,vtype,swath_type):
    if 'integer' in vtype:
        argType='int *'
    if 'real' in vtype:
        argType="float *"
    sL=[]
    s='call copy_%s_%s(dPRData%%%s(:,i,j),dPRRet%%n9(:,i,j), i-1)\n'%(vname.lower(),swath_type.lower(),vname)
    foutF.write(s)

def write_10_call(vname,sname,vtype,swath_type):
    if 'integer' in vtype:
        argType='int *'
    if 'real' in vtype:
        argType="float *"
    sL=[]
    s='call copy_%s_%s(dPRData%%%s(:,i,j),env_nodes(:,i), i-1)\n'%(vname.lower(),swath_type.lower(),vname)
    foutF.write(s)

d={"OEQv":'n10',"OETemp":'n10',"OECloud":'nbins',"OEsimTbNonRain":"n13",\
       "OEemis":"n13", "OEemisSigma":"n13", "OEemisA":"n13",  "OEpiaNonRain":"n2"}


#OEpiaNonRain piaNoPrecip

for v in vars[:20]:
    vs=v.split()
    #print(vs[5],vs[6])
    vtype=vs[1]
    ndims=vs[2].count(':')
    sname=vs[6][8:]
    vname=vs[5]
    if ndims==3:
        print(vname,sname)
        swath_type='KuGMI'
        if vname in d.keys():
            if d[vname]=="n13":
                n1=str(13)
                write_3(vname,sname,vtype,swath_type,n1)
                write_3_call(vname,sname,vtype,swath_type)
            if d[vname]=="nbins":
                write_3(vname,sname,vtype,swath_type,"nbins")
                write_3_call(vname,sname,vtype,swath_type)
            if d[vname]=="n10":
                write_n10(vname,sname,vtype,swath_type)
                write_10_call(vname,sname,vtype,swath_type)
            if d[vname]=="n2":
                n1=str(2)
                write_3(vname,sname,vtype,swath_type,n1)
                write_3_call(vname,sname,vtype,swath_type)
    if ndims==2:
         write_2(vname,sname,vtype,swath_type)
         write_2_call(vname,sname,vtype,swath_type)

for v in vars[:20]:
    vs=v.split()
    #print(vs[5],vs[6])
    vtype=vs[1]
    ndims=vs[2].count(':')
    sname=vs[6][8:]
    vname=vs[5]
    if ndims==3:
        print(vname,sname)
        swath_type='KuKaGMI'
        if vname in d.keys():
            if d[vname]=="n13":
                n1=str(13)
                write_3(vname,sname,vtype,swath_type,n1)
                write_3_call(vname,sname,vtype,swath_type)
            if d[vname]=="nbins":
                write_3(vname,sname,vtype,swath_type,"nbins")
                write_3_call(vname,sname,vtype,swath_type)
            if d[vname]=="n10":
                write_n10(vname,sname,vtype,swath_type)
                write_10_call(vname,sname,vtype,swath_type)
            if d[vname]=="n2":
                n1=str(2)
                write_3(vname,sname,vtype,swath_type,n1)
                write_3_call(vname,sname,vtype,swath_type)
    if ndims==2:
         write_2(vname,sname,vtype,swath_type)
         write_2_call(vname,sname,vtype,swath_type)
fout.close()
foutF.close()
