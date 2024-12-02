### 1. hst 
def hst(f,W,V,τ,b,γ=1,T=999999999): #samefunction with hst1realization except print
    #from random import sample

    n=len(f)
    E=W>0
    trajectory = [0]*(τ+1)#cc(0,τ)*0
    flowcount= [0]*(τ+1)
    
    # 1. calculate π0, i.e.,  
    d=np.sum(np.array(W),0)
    π0=d/np.sum(d) ## stationary dist
    
    # 2. Define hst_onewalk function 
    def hst_onewalk(h,W,b,currentnode,flowcount): #supporting hst
        hnext=h.copy()
        # 1. h(u) <- h(u)+b : update current node when flowcound > 0 
        hnext[currentnode]=hnext[currentnode]+b
        # 2. check that: are there any nodes to which snow can flow from u. 
        neighbor=np.where(W[currentnode,:]>0)[1] ## Nu is np.array
        if len(neighbor)==0: downstream=np.array([])
        else: downstream=neighbor[list((np.where(hnext[neighbor]<=hnext[currentnode]))[0])]
        # 3. check flowable 
        flowable = len(downstream)>0 and flowcount < T 
        # 4. determine flow or block
        if flowable==0: # block!
            #hnext[currentnode]=hnext[currentnode]+b ### important!! update current node again
            _node0=np.asscalar(np.random.choice(n, 1, p=π0))
            _neighbor=np.where(W[_node0,:]>0)[1] 
            _downstream=_neighbor[list((np.where(hnext[_neighbor]<=hnext[_node0]))[0])]
            if len(_downstream)==0: 
                nextnode=_node0
            else: 
                nextnode=np.asscalar(np.random.choice(list(_downstream),1))
            flowcount=0
        else: #flow
            nextnode=np.asscalar(np.random.choice(list(downstream),1))
            flowcount=flowcount+1
        return [hnext,flowcount,nextnode]
   
    # 2. h^0
    print('hst start (' +'τ='+str(τ)+', b='+str(b)+')')
    trajectory[0]=np.asscalar(np.random.choice(n,1,p=π0)) # first node 
    flowcount[0]=0
    
    hst_results=pd.DataFrame(data={'Nodename(=v)':V,'h0':f}) 
    
    # 3. h^1 ~ h^τ
    for ℓ in np.arange(1,τ+1): 
        print('\r'+str(ℓ)+'/'+str(τ),sep='',end='')
        Wthreshed=W>np.asmatrix(np.random.random((n,n)))
        hstwalkrslt=hst_onewalk(h=hst_results['h'+str(ℓ-1)],W=Wthreshed,b=b*γ**flowcount[ℓ-1],
                                currentnode=trajectory[ℓ-1],flowcount=flowcount[ℓ-1])
        hst_results['h'+str(ℓ)]=hstwalkrslt[0]
        flowcount[ℓ]=hstwalkrslt[1]
        trajectory[ℓ]=hstwalkrslt[2]
    print('\n'+'hst end')
    
    return hst_results

# def hmat(hstresult):
#     τ=int((hstresult.shape[1]-2))
#     rtn=np.asmatrix(hstresult[sprod('h',np.arange(0,τ+1))])
#     return rtn

def hmat(hst__):
    τ=int((hst__.shape[1]-2))
    rtn=np.asmatrix(hst__[sprod('h',np.arange(0,τ+1))])
    rtn=rtn-np.matrix(np.apply_along_axis(np.mean,1,rtn)).T
    return rtn

def sigmamat(hmat__): #supporting snowdist, #hh:=n*p 
    hmat__=np.array(hmat__)
    n=len(hmat__)
    rtn=np.array(np.zeros([n,n]))
    try: 
        rtn=np.sum((hmat__[:,np.newaxis,:]-hmat__[np.newaxis,:,:])**2,axis=-1)
    except MemoryError:
        for i in np.arange(0,n):
            rtn[i,:]=np.sum((hmat__[i,:]-hmat__[:,:])**2,axis=1)
    return np.asmatrix(rtn)

# def sigmamat_v2(hmat_,τmax=None): 
#     if τmax==None: rtn=l2dist(hmat_)
#     else: rtn=l2dist(hmat_[:,0:(τmax+1)])
#     return rtn

def weightmat(sigmamat__,θ=1): 
    n=len(sigmamat__)
    rtn=np.zeros([n,n])
    for i in np.arange(0,n):
        for j in np.arange(0,n):
            if i==j: rtn[i,j]=0
            else: rtn[i,j]=np.exp(-sigmamat__[i,j]**2/(2*θ**2))
    return rtn

# def total_degree(W):
#     return np.sum(W)

# ### 2. hst: visualization 
# def datavis4ts(f,nodename=None,groupindex=None,
#            figname='temp',figsize=(1,1),dpi=1,cex=1,
#            text=None,fade=0.5):
#     n=len(f)
    
#     if groupindex==None: colors=[0]*n
#     elif groupindex=='continuous': colors=cm.rainbow(np.linspace(1, 0, n))
#     else: colors=np.array(groupindex)
    
#     f=np.array(f)
#     figsize=(10*figsize[0],6*figsize[1])
#     dpi=150*dpi
#     cex=50*cex
#     if text!=None: text=10*text
#     fade=fade
    
#     Fig=plt.figure(figsize=figsize, dpi=dpi) # Make figure object 
#     ax=plt.axes()
#     plt.scatter(np.arange(1,n+1),f,c=colors,s=cex,alpha=fade)
#     style=dict(size=text,color='k')
    
#     if nodename==None:
#         for i in np.arange(1,n+1): 
#             if text!=None: ax.text(i,f[i-1],'%s'% str(i), **style) # numbering index of nodes 
#         rtn=Fig 
#     else: 
#         for i in np.arange(1,n+1): 
#             if text!=None: ax.text(i,f[i-1],'%s'% nodename[i-1], **style) # numbering index of nodes 
#         rtn=Fig 
#     rtn.savefig(figname+'.pdf')

# def datavis4sct(v1,v2,nodename=None,groupindex=None,
#            figname='temp',figsize=(1,1),dpi=1,cex=1,
#            text=None,fade=0.5):
#     n=len(v1)
    
#     if groupindex==None: colors=[0]*n
#     elif groupindex=='continuous': colors=cm.rainbow(np.linspace(1, 0, n))
#     else: colors=np.array(groupindex)
    
#     v1=np.array(v1)
#     v2=np.array(v2)
#     figsize=(10*figsize[0],6*figsize[1])
#     dpi=150*dpi
#     cex=50*cex
#     if text!=None: text=10*text
#     fade=fade
    
#     Fig=plt.figure(figsize=figsize, dpi=dpi) # Make figure object 
#     ax=plt.axes()
#     plt.scatter(v1,v2,c=colors,s=cex,alpha=fade)
#     style=dict(size=text,color='k')
    
#     if nodename==None:
#         for i in np.arange(1,n+1): 
#             if text!=None: ax.text(v1[i-1],v2[i-1],'%s'% str(i), **style) # numbering index of nodes 
#         rtn=Fig 
#     else: 
#         for i in np.arange(1,n+1): 
#             if text!=None: ax.text(v1[i-1],v2[i-1],'%s'% nodename[i-1], **style) # numbering index of nodes 
#         rtn=Fig 
#     rtn.savefig(figname+'.pdf')

# def pca4vis3d(Σ,nodename=None,groupindex=None,
#            figname='temp',title=None,figsize=(1,1),dpi=1,cex=1,
#            xlim=None,ylim=None,zlim=None,text=None,fade=0.5,
#            prnt=False): # size=(size of obs representation, size of text which represent obs index)

#     from sklearn.decomposition import PCA 
#     from mpl_toolkits import mplot3d
#     if prnt==True: print('PCA start')
#     pca=PCA(n_components=3) # PCA start 
#     n=Σ.shape[0]
#     pca.fit(Σ) 
#     pcarslt=pca.transform(Σ) # PCA end 
#     if prnt==True: print('end')

#     if groupindex==None: colors=['gray']*n
#     elif groupindex=='continuous': colors=cm.rainbow(np.linspace(1, 0, n))
#     else: colors=np.array(groupindex)
    
#     figsize=(20*figsize[0],6*figsize[1])
#     dpi=200*dpi
#     cex=50*cex
#     if text!=None: text=10*text
#     fade=fade
    
#     Fig=plt.figure(figsize=figsize, dpi=dpi)  # Make figure object 
#     ax=plt.axes(projection='3d') # define type of axes: 3d plot 
#     if title!=None: ax.text(0,0,zlim[1],title,fontsize=30)
#     if xlim!=None: ax.set_xlim(xlim[0],xlim[1])
#     if ylim!=None: ax.set_ylim(ylim[0],ylim[1])
#     if zlim!=None: ax.set_zlim(zlim[0],zlim[1])
#     ax.scatter3D(pcarslt[:,0],pcarslt[:,1],pcarslt[:,2],s=cex,c=colors,alpha=fade) # drawing each obs by scatter in 3d axes   
#     ax.ticklabel_format(style='sci', axis='x',scilimits=(0,0))
#     ax.ticklabel_format(style='sci', axis='y',scilimits=(0,0))
#     ax.ticklabel_format(style='sci', axis='z',scilimits=(0,0))
#     ax.xaxis.set_major_locator(plt.MaxNLocator(9))
#     ax.yaxis.set_major_locator(plt.MaxNLocator(5))
#     ax.zaxis.set_major_locator(plt.MaxNLocator(3))
#     ax.w_xaxis.pane.fill = False
#     ax.w_yaxis.pane.fill = False
#     ax.w_zaxis.pane.fill = False
#     if prnt==True: print('labeling (observation-wise)')
#     if nodename==None:
#         for i in np.arange(1,n+1): 
#             if prnt==True: print('\r'+str(i),'/'+str(n),sep='',end='')
#             if text!=None: ax.text(pcarslt[i-1,0],pcarslt[i-1,1],pcarslt[i-1,2],'%s'% (str(i)), size=text, zorder=1,color='k') # numbering index of nodes 
#         if prnt==True: print('\n'+'end')
#     else: 
#         for i in np.arange(1,n+1): 
#             if prnt==True: print('\r'+str(i),'/'+str(n),sep='',end='')
#             if text!=None: ax.text(pcarslt[i-1,0],pcarslt[i-1,1],pcarslt[i-1,2],'%s'% (nodename[i-1]), size=text, zorder=1,color='k') # numbering index of nodes 
#         if prnt==True: print('\n'+'end')
#     Fig.savefig(figname+'.pdf')
#     rtn=Fig

# def pca4msvis3d(hstresult,τlist,
#               nodename=None,groupindex=None,
#               figname='temp',figsize=(1,1),dpi=1,cex=1,text=None,fade=1,
#               prnt=False,logscale=(False,False,False)): # size=(size of obs representation, size of text which represent obs index)
#     Σresult=τlist.copy()
#     hh=hmat(hstresult)
#     M=len(τlist)
#     from sklearn.decomposition import PCA 
#     pca=PCA(n_components=3) # PCA start 
#     pca.fit(Sigma(hh)) 
#     pcarslt=pca.transform(Sigma(hh))
#     x0=min(pcarslt[:,0]);x1=max(pcarslt[:,0])
#     y0=min(pcarslt[:,1]);y1=max(pcarslt[:,1])
#     z0=min(pcarslt[:,2]);z1=max(pcarslt[:,2])
#     if prnt==True: print('obtain snowdist')
#     for m in np.arange(0,M):
#         if prnt==True: print('\r'+str(m),'/'+str(M),sep='',end='')
#         Σresult[m]=Sigma(hh,τmax=τlist[m])
#         pca4vis3d(Σresult[m],nodename=nodename,groupindex=groupindex,
#             title='τ='+str(τlist[m]),
#             figname=figname+str(m+1),figsize=figsize,dpi=dpi,
#             cex=cex,text=text,fade=fade,
#             xlim=(x0,x1),ylim=(y0,y1),zlim=(z0,z1))      
#     if prnt==True: print('\n'+'end')    

