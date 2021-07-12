lines=open("out_s2_fs")
for l in lines:
    ind=l.index("_fs")
    l=l[ind+4:]
    ind2=l.index(")")
    if "(" in  l[:ind2+1]:
        s="%s=missing_r4"%l[:ind2+1]
    else:
        ind2=l.index(",")
        s="%s=missing_r4"%l[:ind2]
    print(s)
    
