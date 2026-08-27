
import sys

def IndexConverter(MSAseq):
    poscolIdxs=[]
    colposIdxs=[]
    Pc=0
    Cc=-1
    lastalpha=None
    for col in range(len(MSAseq)):
        if MSAseq[col].isalpha():
            Cc+=1
            if(lastalpha!=None):
                for i in range(((col + lastalpha)//2+1),col):
                    colposIdxs[i]=Cc
            lastalpha=col
            poscolIdxs.append(Pc+Cc)
        else:
            Pc+=1
        
        colposIdxs.append(Cc)
       
    colposIdxs=[None if x == -1 else x for x in colposIdxs]
    return poscolIdxs, colposIdxs

def MuliIndexConverter(MSAseqs):
    posColIdx=[]
    for MSAseq in MSAseqs:
        posColIdx.append(IndexConverter(MSAseq))
    return posColIdx


def main():
    sequences=["FPIKWTAPEAALYGRFTIKSDVWSFGILLTKGRVPYPGMVNREVLDQVERG","FPIKWTAPEAALYGRFTIKSDVWSFGILLTELVTKGRVVMVNREVLEQVERG"]
    MSAseqs=["-FPIKWTAPEAALY----GRFTIKSDVWSFGILL----TKGRVPYPGMVNR-EVLDQVERG","FPIKWTAPEAALY---GRFTIKSDVWSFGILLTELVTKGRV--TV-MVNR-EVLEQVERG"]
    
    posIdxs,colIdxs= IndexConverter(MSAseqs[0])
    print(posIdxs, colIdxs)

if __name__ == '__main__':
    sys.exit(main()) 