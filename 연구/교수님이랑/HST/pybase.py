import numpy as np
import pandas as pd
import rpy2
import rpy2.robjects as ro
# ro.r('library(devtools)') ## to use source_url 
# ro.r('library(tidyverse)')

############################## rpy2 functions ##############################
# def pull(r):
#     return r2p(ro.globalenv[r])

# def push(py,rname=None):
#     import inspect
#     def retrieve_name(var):
#         for fi in reversed(inspect.stack()):
#             names = [var_name for var_name, var_val in fi.frame.f_locals.items() if var_val is var]
#             if len(names) > 0:
#                 return names[0]
#     if rname==None: rname = retrieve_name(py)
#     ro.globalenv[rname]=p2r(py)
    
# def p2r(A):
#     from rpy2.robjects.vectors import FloatVector 
#     from rpy2.robjects.vectors import StrVector as s2r_temp

#     def a2r_temp(a):
#         if type(a) in {float,int,bool}: a=[a]
#         a=list(a)
#         rtn=FloatVector(a)
#         return rtn

#     def m2r_temp(A):
#         A=np.matrix(A)
#         Acopy=A.T.copy()
#         nrow=Acopy.shape[0]
#         Acopy.shape=(np.prod(Acopy.shape),1)
#         rtn=ro.r.matrix(a2r_temp(list(Acopy)),ncol=nrow)
#         return rtn

#     from rpy2.robjects import pandas2ri
#     from rpy2.robjects.conversion import localconverter

#     def pd2r_temp(A):
#         with localconverter(ro.default_converter + pandas2ri.converter):
#             rtn = ro.conversion.py2rpy(A)
#         return rtn
    
#     if type(A)==type(pd.DataFrame(np.zeros([2,2]))):
#         rtn=pd2r_temp(A) 
#     elif type(A)==type(np.matrix(np.zeros([2,2]))):
#         rtn=m2r_temp(A)
#     elif type(A)==type(np.zeros([2,2])):
#         if len(A.shape)==1: 
#             rtn=a2r_temp(A)
#         else:
#             rtn=m2r_temp(A)
#     elif type(A)==str: 
#         rtn=s2r_temp(A)        
#     elif type(pd.DataFrame(np.matrix(A)).iloc[0,0])==str:
#         rtn=s2r_temp(pd.DataFrame(np.matrix(A)).T.iloc[:,0])
#     else:
#         rtn=a2r_temp(A)
#     return rtn

# def r2p(A):
#     from rpy2.robjects import pandas2ri
#     from rpy2.robjects.conversion import localconverter

#     def r2a_temp(a):
#         return list(a)
    
#     def r2m_temp(A):
#         return np.matrix(A)
    
#     def r2pd_temp(A):
#         with localconverter(ro.default_converter + pandas2ri.converter):
#             rtn = ro.conversion.rpy2py(A)
#         return rtn   
    
#     ro.globalenv['temp']=A
#     if ro.r('is.null(dim(temp))')[0]==False: ## in the cases of matrix or dataframe
#         if ro.r('is.data.frame(temp)')[0]: 
#             rtn=r2pd_temp(A)
#         elif ro.r('is.matrix(temp)')[0]:
#             rtn=r2m_temp(A)
#         else:
#             print('I don\`t know which type of this data in R.')
#     else:
#         rtn=r2a_temp(A)
#     ro.r('rm("temp")')
#     return rtn

def cbind(*Mat):
    lenofMat=len(Mat)
    if lenofMat==1: 
        print("You must enter two or more input objects.")
        rtn=Mat[0]
    elif lenofMat==2: 
        rtn=cbindtemp(Mat[0],Mat[1])
    else: 
        rtn=cbindtemp(Mat[0],Mat[1])
        for i in np.arange(2,lenofMat):
            rtn=cbindtemp(rtn,Mat[i])
    return rtn

def rbind(*Mat):
    lenofMat=len(Mat)
    if lenofMat==1: 
        print("You must enter two or more input objects.")
        rtn=Mat[0]
    elif lenofMat==2: 
        rtn=rbindtemp(Mat[0],Mat[1])
    else: 
        rtn=rbindtemp(Mat[0],Mat[1])
        for i in np.arange(2,lenofMat):
            rtn=rbindtemp(rtn,Mat[i])
    return rtn

def cbindtemp(A,B):
    typ=['matrix','matrix']
    
    if isinstance(A, pd.core.series.Series): 
        A=a2c(A)
    if isinstance(B, pd.core.series.Series): 
        B=a2c(B)
    A=np.asmatrix(A)
    B=np.asmatrix(B)

    # row-vector에 대한 처리 
    if A.shape[0]==1: typ[0]='rowvec'
    if B.shape[0]==1: typ[1]='rowvec'

    # col-vector에 대한 처리 
    if A.shape[1]==1: typ[0]='colvec'
    if B.shape[1]==1: typ[1]='colvec'    
    
    # 스칼라에 대한 처리 
    if A.shape==(1,1): typ[0]='scala'
    if B.shape==(1,1): typ[1]='scala'
        
    if typ==['scala','scala']:  A=np.array(A); B=np.array(B)
    if typ==['scala','rowvec']: A=np.array(A); 
    if typ==['scala','colvec']: A=np.full(B.shape,A[0,0]); 
    if typ==['scala','matrix']: A=np.full((B.shape[0],1),A[0,0]); 

    if typ==['rowvec','scala']: B=np.array(B)
    #if typ==['rowvec','rowvec']:
    if typ==['rowvec','colvec']: A=A.T
    if typ==['rowvec','matrix']: A=A.T
        
    if typ==['colvec','scala']:  B=np.full(A.shape,B[0,0])
    if typ==['colvec','rowvec']: B=B.T
    #if typ==['colvec','colvec']: 
    #if typ==['colvec','matrix']: 
    
    if typ==['matrix','scala']:  B=np.full((A.shape[0],1),B[0,0])
    if typ==['matrix','rowvec']: B=B.T
    #if typ==['matrix','colvec']: 
    #if typ==['matrix','matrix']:
    
    return np.hstack([A,B])
    
def rbindtemp(A,B):
    typ=['matrix','matrix']
    
    A=np.asmatrix(A)
    B=np.asmatrix(B)

    # row-vector에 대한 처리 
    if A.shape[0]==1: typ[0]='rowvec'
    if B.shape[0]==1: typ[1]='rowvec'

    # col-vector에 대한 처리 
    if A.shape[1]==1: typ[0]='colvec'
    if B.shape[1]==1: typ[1]='colvec'    
    
    # 스칼라에 대한 처리 
    if A.shape==(1,1): typ[0]='scala'
    if B.shape==(1,1): typ[1]='scala'
        
    if typ==['scala','scala']:  A=np.array(A); B=np.array(B)
    if typ==['scala','rowvec']: A=np.full(B.shape,A[0,0]); 
    if typ==['scala','colvec']: A=np.array(A);
    if typ==['scala','matrix']: A=np.full((1,B.shape[1]),A[0,0]); 

    if typ==['rowvec','scala']: B=np.full((1,A.shape[1]),B[0,0]); 
    #if typ==['rowvec','rowvec']:
    if typ==['rowvec','colvec']: B=B.T
    #if typ==['rowvec','matrix']: 
        
    #if typ==['colvec','scala']:  
    if typ==['colvec','rowvec']: A=A.T
    #if typ==['colvec','colvec']: 
    if typ==['colvec','matrix']: A=A.T
    
    if typ==['matrix','scala']:  B=np.full((1,A.shape[1]),B[0,0])
    #if typ==['matrix','rowvec']: 
    if typ==['matrix','colvec']: B=B.T
    #if typ==['matrix','matrix']:
    
    return np.vstack([A,B])

def ids(pddata):
    push(pddata.columns,"vname")
    print(r2p(ro.r("str_c(str_c('(',str_c(1:length(vname)-1),') ',vname),collapse='\n')"))[0])

def l2distance(X): #X:=n*p ndarray
    X=np.array(X)
    n=len(X)
    rtn=np.array(np.zeros([n,n]))
    try: 
        rtn=np.sum((X[:,np.newaxis,:]-X[np.newaxis,:,:])**2,axis=-1)
    except MemoryError:
        for i in np.arange(0,n):
            rtn[i,:]=np.sum((X[i,:]-X[:,:])**2,axis=1)
    return rtn
