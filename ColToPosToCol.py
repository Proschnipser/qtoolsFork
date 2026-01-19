
import sys

def IndexConverter(MSAseq):
    poscolIdxs=[]
    colposIdxs=[]
    Pc=0
    Cc=-1

    for col in range(len(MSAseq)):
        if MSAseq[col].isalpha():
            Cc+=1
            poscolIdxs.append(Pc+Cc)
        else:
            Pc+=1
        colposIdxs.append(Cc)

    return poscolIdxs, colposIdxs

def main():
    sequences=["FPIKWTAPEAALYGRFTIKSDVWSFGILLTKGRVPYPGMVNREVLDQVERG","FPIKWTAPEAALYGRFTIKSDVWSFGILLTELVTKGRVVMVNREVLEQVERG"]
    MSAseqs=["-FPIKWTAPEAALY---GRFTIKSDVWSFGILL----TKGRVPYPGMVNR-EVLDQVERG","FPIKWTAPEAALY---GRFTIKSDVWSFGILLTELVTKGRV--TV-MVNR-EVLEQVERG"]
    
    posIdxs,colIdxs= IndexConverter(MSAseqs[0])
    print(posIdxs, colIdxs)

if __name__ == '__main__':
    sys.exit(main()) 