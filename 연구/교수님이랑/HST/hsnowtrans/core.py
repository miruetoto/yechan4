import numpy as np
#from plotnine import * 
import tqdm 

class GraphSignal: 
    def __init__(self,V,W,f):
        self.f=np.array(f) 
        self.V=V 
        self.W=np.array(W)
        self.n=len(self.f)
        self.degree=np.sum(np.array(self.W),0)
        self.initdist=self.degree/np.sum(self.degree) 
        
class Heavysnow:
    def __init__(self,graphsignal):
        self.n=graphsignal.n
        self.f=graphsignal.f
        self.V=graphsignal.V 
        self.nodeindex=np.arange(0,self.n,dtype='int')
        self.graphweight=graphsignal.W
        self.initdist=graphsignal.initdist
        self._initialize()
    
    def _initialize(self):
        self.tau=None
        self.b=None
        self.trajectory=[np.random.choice(self.n,1,p=self.initdist).item()]
        self.flowcount=[0]
        self.snowygrounds=None 
        self.eucliddistance=(np.array(self.f)[:,np.newaxis]-np.array(self.f)[np.newaxis,:])**2
        self.snowdistance=np.zeros((self.n,self.n))
        self.euclidweight=np.zeros((self.n,self.n))
        self.snowweight=np.zeros((self.n,self.n))
        self.status=['init']
        self.theta=None
        
    def _snowonce(self,ell,maxflow):
        b=self.b
        flowcount=self.flowcount[-1] 
        snowyground=self.snowygrounds[...,ell-1].copy()
        currentnode=self.trajectory[-1] 
        neighbor=self.nodeindex[(self.graphweight>0)[currentnode]]
        transitionprob=None
        nextnode=None 
        downstream=None
        
        if flowcount > maxflow: # reset 
            nextnode=np.random.choice(self.n,1,p=self.initdist).item()
            snowyground[nextnode]=snowyground[nextnode]+b
            flowcount=0
            self.status.append('reset')
        elif sum(neighbor)==0: # empty neighborhood 
            nextnode=np.random.choice(self.n,1,p=self.initdist).item()
            snowyground[nextnode]=snowyground[nextnode]+b
            flowcount=0
            self.status.append('empty_neighborhood')
        else:
            _transitionprob=self.graphweight[currentnode]/sum(self.graphweight[currentnode])
            nextnode=np.random.choice(self.n,1,p=_transitionprob).item()
            downstream=neighbor[np.where(snowyground[currentnode]>=snowyground[neighbor])]
            if nextnode in set(downstream): # flow 
                snowyground[nextnode]=snowyground[nextnode]+b
                flowcount=flowcount+1
                self.status.append('flow')
            else: # block 
                nextnode=np.random.choice(self.n,1,p=self.initdist).item()
                snowyground[currentnode]=snowyground[currentnode]+b
                flowcount=0
                self.status.append('block')
        self.snowygrounds[...,ell]=snowyground
        self.flowcount=self.flowcount+[flowcount]
        self.trajectory=self.trajectory+[nextnode]
        
    def _getdegreematrix(self,matrix):
        return np.diag(np.sum(matrix,1))
    
    def _normalize(self,matrix): 
        return np.sqrt(self._getdegreematrix(matrix))@matrix@np.sqrt(self._getdegreematrix(matrix))
    
    def _updateeuclidweight(self):
        self.euclidweight=np.exp(-self.eucliddistance/(2*self.theta**2))-np.eye(self.n)
        
    def _updatesnowdistance(self):
        try: 
            self.snowdistance=np.sqrt(np.sum((self.snowygrounds[:,np.newaxis,:]-self.snowygrounds[np.newaxis,:,:])**2,axis=-1))
        except MemoryError:
            print("Due to insufficient memory, the distance is calculated using a for loop.") 
            for u in tqdm.tqdm(range(self.n)):                
                for v in range(self.n):                    
                    self.snowdistance[u,v]=np.sqrt(np.sum((self.snowygrounds[u,:] - self.snowygrounds[v,:])**2))
            self.snowdistance = self.snowdistance + self.snowdistance.T 
             
    def _updatesnowweight(self):
        self.snowweight=np.exp(-self.snowdistance/(2*self.tau*self.theta**2))-np.eye(self.n)
    
    def snow(self,tau,b=1,maxflow=100000,prnt=True):
        self._initialize()
        self.b=b
        self.theta=np.sqrt(1/2)*self.b
        self.tau=tau
        self.snowygrounds=np.repeat(self.f,self.tau+1).reshape(self.n,self.tau+1)

        if prnt: 
            print('HST (tau= %s, b=%s)' % (self.tau,self.b))
            for ell in tqdm.tqdm(range(1,self.tau+1)):
                self._snowonce(ell,maxflow)
            print('Calculate distance and weights')
            self._updatesnowdistance()
            self._updatesnowweight()
            self._updateeuclidweight()  
            print('HST completed and all history is recorded.')
        else:
            for ell in range(1,self.tau+1):
                self._snowonce(ell,maxflow)
            self._updatesnowdistance()
            self._updatesnowweight()
            self._updateeuclidweight()  

    def snow2(self): ## 시뮬레이션용 (tau, b등이 미리 저장되어있어야함) 
        maxflow=100000 
        self.theta=np.sqrt(1/2)*self.b
        self.snowygrounds=np.repeat(self.f,self.tau+1).reshape(self.n,self.tau+1)
        for ell in range(1,self.tau+1):
            self._snowonce(ell,maxflow)
        self._updatesnowdistance()
        self._updatesnowweight()
        self._updateeuclidweight()                 

    def _snow_for_simulation(self,tau,b=1,maxflow=100000): # 시뮬레이션용 (시간측정만)
        self._initialize()
        self.b=b
        self.theta=np.sqrt(1/2)*self.b
        self.tau=tau
        self.snowygrounds=np.repeat(self.f,self.tau+1).reshape(self.n,self.tau+1)
        for ell in tqdm.tqdm(range(1,self.tau+1)):
            self._snowonce(ell,maxflow)
        
    def adjustingtheta(self,theta): 
        self.theta=theta
        self._updatesnowweight()
        self._updateeuclidweight()                

class SpectralAnalysis:        
    def __init__(self,graphsignal):
        self.f=graphsignal.f
        self.V=graphsignal.V
        self.W=graphsignal.W
        self.n=graphsignal.n
        self.degree=graphsignal.degree        
        self._initialize()
    
    def _initialize(self):
        self.D=np.diag(self.degree)
        self.L=self.D-self.W
        _drootinv=np.zeros_like(self.degree)
        _drootinv[self.degree>1e-2]= 1/np.sqrt(self.degree[self.degree>1e-2])
        self.Lz=np.diag(_drootinv)@(self.L)@np.diag(_drootinv)
        self.Psi= None 
        self.lamb= None
        self.Lamb= None
        self.fbar= None
        self.components=None
    
    def graphFouriertransform(self):
        self.Psi,self.lamb,_ = np.linalg.svd(self.Lz)
        self.Lamb = np.diag(self.lamb)
        self.fbar = self.Psi.T @ self.f 
        
    def decompose(self):
        _complist = [(np.outer(self.Psi[i],self.Psi[i])@self.f).tolist() for i in range(self.n)] 
        self.components = np.array(_complist).T

    def _specplot(self):
        thm = theme_classic()
        aes1= aes(x=self.lamb,y=np.abs(self.fbar),xend=self.lamb,yend=0)
        aes2= aes(x=self.lamb,y=np.abs(self.fbar))
        g1= geom_segment(aes1)  
        g2= geom_point(aes2)    
        fig = ggplot()+thm+g1+g2+xlim(0,2.1)+ xlab("")+ ylab("") + ggtitle("Periodogram")
        return fig 
    
#     def _decomplot(self):       
# #         ggplot(data=tibble(V=V,f=f,Vtext=Vtext),aes(x=V,y=f,label=Vtext))+
# #         geom_col(aes(fill=(f>0)),width=0.1)+geom_hline(aes(yintercept=0),col="gray60",lty=2)+
# #         geom_text(fontface = 4,size=8)+
# #         xlab("")+ylab("")+guides(fill=FALSE)+theme(plot.title=element_text(face="bold.italic"))+
# #         theme_bw()+
# #         theme(strip.text.x = element_text(size = 20, color = "black", face = "bold.italic"))+
# #         theme(strip.text.y = element_text(size = 15, color = "black", face = "bold.italic"))+
# #         ylim(-1.2,1.2)+
# #         theme(plot.title=element_text(face="bold.italic"))
        
        
#         df=
#         thm = theme_classic()
#         aes1= aes(x=self.V,y=np.abs(self.fbar),xend=self.lamb,yend=0)
#         aes2= aes(x=self.lamb,y=np.abs(self.fbar))
#         g1= geom_col(aes1)  
#         g2= geom_text(aes2)    
#         fig = ggplot()+thm+g1+g2+xlim(0,2.1)+ xlab("")+ ylab("") + ggtitle("Periodogram")
#         return fig 
        