import scipy as sp

class base(object):
    _instances = []
    def __init__(self):
        self._instances.append(self) #Track all instances derived from this class for kicks
    
    def find_objects_by_name(self,name):
        for item in self._instances:
            if item.name == name:
                obj = item
        return obj
    
    def print_dicts(self):
        print('Data dictionaries:')
        for item in self._data.keys():
            print('  '+item)
        print('Info dictionaries:')
        for item in self._info.keys():
            print('  '+item)

    def show_associations(self):
        if self.__class__.__name__ == 'geometry':
            try: print(self._net.name+' <'+self._net.__class__.__name__+'>')
            except: pass
        try: 
            for item in self._fluid:
                print(item.name+' <'+item.__class__.__name__+'>')
        except: 
            try: print(self._fluid.name+' <'+self._fluid.__class__.__name__+'>')
            except: pass
        try: 
            for item in self._geometry: 
                print(item.name+' <'+item.__class__.__name__+'>')
        except:
            try: print(self._geometry.name+' <'+self._geometry.__class__.__name__+'>')
            except: pass
        try: 
            for item in self._physics: 
                print(item.name+' <'+item.__class__.__name__+'>')
        except:
            try: print(self._physics.name+' <'+self._physics.__class__.__name__+'>')
            except: pass
    
class tools(object):
    def __init__(self):
        pass
    
    def get_num_pores(self,subdomain=['all']):
        if type(subdomain) == str: subdomain = [subdomain] #convert string to list, if necessary
        Np = sp.shape(self._data['nums'])[0]
        if subdomain == ['all']: #return all pores
            return Np
        else:
            if type(subdomain) == str: subdomain = [subdomain] #convert string to list, if necessary
            temp = sp.zeros((Np,),dtype=bool)
            for item in subdomain: #iterate over subdomain list and accumulate Trues
                temp = temp + self._info[item]
            return sp.sum(temp) #return sum of Trues
            
    def get_pore_indices(self,subdomain=['all'],indices=True,mode='union'):
        r'''
        '''
        if type(subdomain) == str: subdomain = [subdomain] #convert string to list, if necessary
        if subdomain == ['all']: #Return full index; easier than get_data(prop='nums')
            if indices:
                ind = sp.r_[0:self.get_num_pores()]
            else:
                ind = sp.ones((self.get_num_pores(),),dtype=bool)
        else:
            if mode == 'union':
                union = sp.zeros((self.get_num_pores(),),dtype=bool)
                for item in subdomain: #iterate over subdomain list and collect all indices
                    union = union + self.get_info(prop=item)
                ind = union
            elif mode == 'intersection':
                intersect = sp.ones((self.get_num_pores(),),dtype=bool)
                for item in subdomain: #iterate over subdomain list and collect all indices
                    intersect = intersect*self.get_info(prop=item)
                ind = intersect
            if indices: ind = sp.where(ind==True)[0]
        return ind
    
    def get_data(self,subdomain='',phase='',prop=''):
        r''' 
        '''
        try: subdomain = subdomain.name #allow passing of geometry objects
        except: pass #Otherwise, accept string
        if phase and not subdomain: return phase._data[prop] #Get fluid prop
        elif subdomain and not phase: #Get geometry property
            ind = self.get_pore_indices(subdomain)
            return self._data[prop][ind]
        elif fluid and subdomain: #Get physics property
            ind = self.get_pore_indices(subdomain)
            return phase._data[prop][ind]
        elif not (phase or subdomain): return self._data[prop] #Get topology property
            
    def set_data(self,subdomain='',phase='',prop='',data=''):
        r'''
        '''
        try: subdomain = subdomain.name #allow passing of geometry objects
        except: pass #Otherwise, accept string
        if phase and not subdomain: phase._data[prop] = sp.array(data,ndmin=1) #Set fluid property
        elif subdomain and not phase: #Set geometry property
            ind = self.get_pore_indices(subdomain)
            try: self._data[prop] #Test existance of prop
            except: self._data[prop] = sp.zeros((self.get_num_pores(),))*sp.nan
            if sp.shape(ind) == sp.shape(data): self._data[prop][ind] = data
            elif sp.shape(sp.array(data,ndmin=1))[0] == 1: self._data[prop][ind] = data
            else: print('data is the wrong size!')
        elif phase and subdomain: #Set pore scale physics property
            ind = self.get_pore_indices(subdomain)
            try: phase._data[prop]
            except: phase._data[prop] = sp.zeros((self.get_num_pores(),))
            phase._data[prop][ind] = sp.array(data,ndmin=1)
        elif not (phase or subdomain): self._data[prop] = sp.array(data,ndmin=1) #Set topology property
        
    def set_info(self,prop='',data='',indices=False):
        r'''
        '''
        if indices:
            try: self._info[prop]
            except: self._info[prop] = sp.zeros((self.get_num_pores(),),dtype=bool)
            self._info[prop][data] = True
        else:
            self._info[prop] = sp.array(data,dtype=bool,ndmin=1)
            
    def get_info(self,prop='',indices=False):
        r'''
        '''
        if indices:
            return sp.where(self._info[prop]==True)[0]
        else:
            return self._info[prop]

class generic_network(base,tools): #Generic network class, inherits from base and tools
    _data = {}
    _info = {}
    def __init__(self):
        self._geometry = []
        
    def generate(self,**kwargs):
        pass
        
    def visualize_model(self):
        try: nx
        except: 
            import networkx as nx
        G = nx.Graph()
        #Scan connections
        for item1 in self._geometry:
            G.add_edge('network',item1.name)
            for item2 in item1._physics:
                G.add_edge(item1.name,item2.name)
                G.add_edge(item2.name,item2._fluid.name)
        ps=nx.spring_layout(G,iterations=1000)
        nx.draw_networkx_edges(G,ps)
        nx.draw_networkx_nodes(G,ps,node_size=1000)
        nx.draw_networkx_labels(G,ps,font_size=10,font_family='sans-serif')
        plt.axis('off')        
        plt.show() # display
        return G
        
class specific_network(generic_network): #Subclass of generic_network, ovrerloading generate
    def __init__(self,name):
        print('initializing network instance')
        super(specific_network,self).__init__()
        self.name = name

    def generate(self,Np): #The network can create itself
        print('creating topology in new instance')
        self._generate_pores(Np)
        return self
        
    def _generate_pores(self,Np):
        self.set_data(prop='nums',data=sp.r_[0:Np])
        
class geometry(base): #Size calculation methods
    def __init__(self,network,name='main'):
        print('initializing geometry instance for '+name)
        super(geometry,self).__init__()
        self.name = name
        self.subdomain = name
        self._net = network
        self._net._geometry.append(self)
        self._physics = []
        
    def calc_diameter(self,f=1):
        print('appying diameters to '+self.subdomain)
        Np = self._net.get_num_pores(self.name)
        self._net.set_data(subdomain=self.name, prop='diameter',data=sp.rand(Np)*f)

class fluid(base,tools): #Inherits from base, but also tools so it can access info about itself
    _data = {}
    _info = {}
    def __init__(self,network,name):
        print('initializing fluid instance for '+name)
        super(fluid,self).__init__()
        self._physics = []
        self.name = name
        self._net = network
        self._net.set_data(phase=self,prop='temperature',data=298)
        self.set_data(prop='nums',data=self._net.get_pore_indices()) #This is necessary for the methods from 'tools' to work.  They must know network size.
        
    def diffusivity(self,value):
        self._net.set_data(phase=self,prop='diffusivity',data=value)

class physics(base): #Pore scale physics calculations
    def __init__(self,network,fluid,geometry,name):
        super(physics,self).__init__()
        print('initializing physics instance for '+fluid.name+' in '+geometry.name)
        self.name = name
        self._fluid = fluid #Attach fluid object
        self._net = network #Attach network object
        self._geometry = geometry #Attach geometry object
        self._fluid._physics.append(self)
        self._geometry._physics.append(self)
        self.subdomain = self._geometry.name
        
    def diffusive_conductance(self):
        T= self._net.get_data(phase=self._fluid,prop='temperature')
        DAB = self._net.get_data(phase=self._fluid,prop='diffusivity')
        A = self._net.get_data(subdomain=self._geometry,prop='diameter')
        self._net.set_data(subdomain=self._geometry,phase=self._fluid,prop='diffusive_conductance', data=(A)*(DAB)/(T))
        
if __name__ == '__main__':
    #Create base network object
    Np = 10
    pn = specific_network(name='net1').generate(Np)
    
    #Create base fluid objects
    water = fluid(network=pn,name='water')
    air = fluid(network=pn,name='air')
    water.diffusivity(value=1)
    air.diffusivity(value=100)
    
    #Define geometry for each subdomain
    GDL = geometry(pn,name='GDL')
    MPL = geometry(pn,name='MPL')

    #Define two domains GDL and MPL
    subdomain_GDL = sp.array([1,1,1,1,1,1,1,0,0,0])
    pn.set_info(prop='GDL',data=subdomain_GDL)
    
    subdomain_MPL = [7,8,9]
    pn.set_info(prop='MPL',data=subdomain_MPL,indices=True)
    
    #Set some pores as top, front and right
    pn.set_info(prop='top',data=[1,2,3],indices=True)
    pn.set_info(prop='front',data=[1,4,5],indices=True)
    pn.set_info(prop='side',data=[1,2,6],indices=True)
    #Now find corners
    corners = pn.get_pore_indices(['top','front','side'],mode='intersection',indices=False)
    
    GDL.calc_diameter()
    MPL.calc_diameter(f=0.001)
    
    #Define physics for each domain and fluid
    phys1 = physics(network=pn,fluid=water,geometry=GDL,name='phys1')
    phys2 = physics(network=pn,fluid=air,geometry=GDL,name='phys2')
    phys3 = physics(network=pn,fluid=water,geometry=MPL,name='phys3')
    phys4 = physics(network=pn,fluid=air,geometry=MPL,name='phys4')
    
    #Calculate some pore scale physics properties
    phys1.diffusive_conductance()
    phys2.diffusive_conductance()
    phys3.diffusive_conductance()
    phys4.diffusive_conductance()
    