import struct 
import numpy as np
from anytree import NodeMixin, RenderTree
from anytree import LevelOrderIter
lookup =       {'00464542': 'FEB',
                '01000000': 'ROOT',
                '01010000': 'HEADER',
                '01010001': 'VERSION',
                '01010002': 'NODES',
                '01010003': 'MAX_FACET_NODES',
                '01010004': 'COMPRESSION',
                '01010005': 'AUTHOR',
                '01010006': 'SOFTWARE',
                '01020000': 'DICTIONARY',
                '01021000': 'GLOBAL_VAR',
                '01022000': 'MATERIAL_VAR',
                '01023000': 'NODESET_VAR',
                '01024000': 'DOMAIN_VAR',
                '01025000': 'SURFACE_VAR',
                '01020001': 'DICTIONARY_ITEM',
                '01020002': 'ITEM_TYPE',
                '01020003': 'ITEM_FORMAT',
                '01020004': 'ITEM_NAME',
                '01020005': 'ITEM_UNKNOWN',
                '01030000': 'MATERIALS',
                '01030001': 'MATERIAL',
                '01030002': 'MATERIAL_ID',
                '01030003': 'MATERIAL_NAME',
                '01040000': 'MESH',
                '01041000': 'NODE_SECTION',
                '01041100': 'NODE_HEADER',
                '01041101': 'NODES',
                '01041102': 'DIM',
                '01041103': 'NAME',
                '01041200': 'NODE_COORDS',
                '01042000': 'DOMAIN_SECTION',
                '01042100': 'DOMAIN',
                '01042101': 'DOMAIN_HEADER',
                '01042102': 'ELEM_TYPE',
                '01042103': 'PART_ID',
                '01032104': 'ELEMENTS',
                '01032105': 'NAME',
                '01042200': 'ELEMENT_LIST',
                '01042201': 'ELEMENT',
                '01043000': 'SURFACE_SECTION',
                '01043100': 'SURFACE',
                '01043101': 'SURFACE_HEADER',
                '01043102': 'SURFACE_ID',
                '01043103': 'FACES',
                '01043200': 'FACE_LIST',
                '01043201': 'FACE',
                '01044000': 'NODESET_SECTION',
                '01044100': 'NODESET',
                '01044101': 'NODESET_HEADER',
                '01044102': 'NODESET_ID',
                '01044103': 'NODESET_NAME',
                '01044104': 'NODESET_SIZE',
                '01044200': 'NODELIST',
                '01045000': 'PARTS_SECTION',
                '01045100': 'PART',
                '01045101': 'PART_ID',
                '01045102': 'PART_NAME',
                '02000000': 'STATE_SECTION',
                '02010000': 'STATE_HEADER',
                '02010002': 'STATE_TIME',
                '02010003': 'STATE_UNKNOWN',
                '02020000': 'STATE_DATA',
                '02020001': 'STATE_VAR',
                '02020002': 'VARIABLE_ID',
                '02020003': 'VARIABLE_DATA',
                '02020100': 'GLOBAL_DATA',
                '02020200': 'MATERIAL_DATA',
                '02020300': 'NODE_DATA',
                '02020400': 'DOMAIN_DATA',
                '02020500': 'SURFACE_DATA',
                '02030000': 'STATE_UNKNOWN2',
                }

def asHexString(x):
    return "{:08x}".format(x)

class xplt_node(NodeMixin):
    def __init__(self,name,size=None,read=False,parent=None,children=None):
        self.name = name
        self.size = size
        self.read = read
        self.parent = parent
        self.debug = False
        if children:
            self.children = children

    def set(self):
        self.read = True
        if self.parent is not None:
            self.parent.set()

    def read_file(self,fid):
        chunk = np.fromfile(fid,dtype=np.uint32,count=1)
        if not lookup[asHexString(chunk[0])]==self.name:
            if self.debug:
                print("ID not met",self.name,lookup[asHexString(chunk[0])])
                print("Assuming",self.name,"not present and continuing")
            fid.seek(-4,1)
            #return
            raise ValueError("ID not met",self.name,lookup[asHexString(chunk[0])])
        if self.debug:
            print('Reading',self.name)
        if self.name != 'FEB':
            self.size = np.fromfile(fid,dtype=np.uint32,count=1)[0]
        if self.debug:
            print('Reading',self.name,self.size)
        if self.read == True:
            if self.is_leaf:
                if self.size>0:
                    globals()["CALL_"+self.name](fid,self)
                else:
                    if self.debug:
                        print(self.name, 'has a size',self.size,'skipping reading it')
            else:
                for c in self.children:
                    try:
                        c.read_file(fid)
                    except:
                        c.parent=None #essentially delete the node
        else:
            fid.seek(self.size,1)

    def rank(self):
        i=0
        for node in LevelOrderIter(self.parent.parent,filter_=lambda n: n.depth==self.depth and n.name==self.name):
            if node is self:
                return i
            i += 1
        return i

def find_nstates(fid):
    nstate = 1
    feb.read_file(fid)
    while(True):
        try:
            state_section[nstate-1].read_file(fid)
            nstate += 1
            state_section.append(xplt_node('STATE_SECTION',parent=feb))
        except:
            print("Number of states is",nstate)
            return nstate

def read1int32(fid):
    return np.fromfile(fid,dtype=np.uint32,count=1)[0]

def readchar(fid,n):
    strbits = np.fromfile(fid, dtype=np.int8, count=n)
    a_string = ''
    for i in strbits:
        if i>0:
            a_string += chr(i)
        else:
            return a_string
    return a_string

def CALL_DIM(fid,node):
    nDim = read1int32(fid)
    print("Number of dimensions is",nDim)

def CALL_SOFTWARE(fid,node):
    chunk = np.fromfile(fid,dtype=np.int32,count=1)
    print('Software is',readchar(fid,chunk[0]))

def CALL_ITEM_TYPE(fid,node):
    print('Item type is',read1int32(fid))

def CALL_ELEMENTS(fid,node):
    nElem = read1int32(fid)
    print('Number of elements is',nElem)

def CALL_ELEM_TYPE(fid,node):
    print('Element type is',read1int32(fid))

def CALL_ITEM_FORMAT(fid,node):
    print('Item format is',read1int32(fid))

def CALL_ITEM_UNKNOWN(fid,node):
    chunk = np.fromfile(fid,dtype=np.int32,count=1)
    print('Unknown is',readchar(fid,chunk[0]))

def CALL_ITEM_NAME(fid,node):
    print('Item name is',readchar(fid,64))

def CALL_NODES(fid,node):
    nNodes = read1int32(fid)
    print("Number of nodes is",nNodes)

def CALL_NODE_COORDS(fid,node):
    nNodes = 1111
    nDim = 3
    for i in range(nNodes):
        print(np.fromfile(fid,dtype=np.uint32,count=1)[0])
        coords = np.fromfile(fid,dtype=np.float32,count=nDim)
        print(coords)

def CALL_STATE_TIME(fid,node):
    #state_time = struct.unpack('f',fid.read(4))[0]
    state_time = np.fromfile(fid,dtype=np.float32,count=1)[0]
    print("State time is",state_time)

def CALL_VARIABLE_ID(fid,node):
    print('Variable_ID is',read1int32(fid))

def CALL_VARIABLE_DATA(fid,node):
    print(node.parent.parent.name, node.rank())
    chunk = np.fromfile(fid,dtype=np.uint32,count=1) #unsure what this refers to
    chunk2= np.fromfile(fid,dtype=np.uint32,count=1) #another data size 
    #print(chunk,chunk2)
    if node.parent.parent.name == 'NODE_DATA':
        chunk = np.fromfile(fid,dtype=np.float32,count=1111*3)
    elif node.parent.parent.name == 'DOMAIN_DATA':
        #print(node.rank())
        chunk = np.fromfile(fid,dtype=np.float32,count=1061*6)
    print(chunk)

feb = xplt_node('FEB')
feb.set()
root = xplt_node('ROOT',parent=feb)
mesh = xplt_node('MESH',parent=feb)
state_section = []
state_section.append(xplt_node('STATE_SECTION',parent=feb))

xpltFile = 'resample.xplt'
fid = open(xpltFile, 'rb')
fid.seek(0,2)
file_size=fid.tell()
fid.seek(0)
nstates = find_nstates(fid)
fid.seek(0)

header = xplt_node('HEADER',parent=root)
dictionary = xplt_node('DICTIONARY',parent=root)

version = xplt_node('VERSION',parent=header)
max_facet_nodes = xplt_node('MAX_FACET_NODES',parent=header)
compression = xplt_node('COMPRESSION',parent=header)
author = xplt_node('AUTHOR',parent=header)
software = xplt_node('SOFTWARE',parent=header)

global_var  = xplt_node('GLOBAL_VAR',parent=dictionary)
nodeset_var  = xplt_node('NODESET_VAR',parent=dictionary)
domain_var  = xplt_node('DOMAIN_VAR',parent=dictionary)
surface_var  = xplt_node('SURFACE_VAR',parent=dictionary)

dictionary_item_global = []
dictionary_item_nodeset = []
dictionary_item_domain = []
dictionary_item_surface = []
item_type = []
item_format = []
item_unknown = []
item_name = []
MAXVAR = 2
for i in range(MAXVAR): #assume there are maximum MAXVAR dictionary items for each (global, nodeset, domain, and surface)
    dictionary_item_global.append(xplt_node('DICTIONARY_ITEM',parent=global_var))
    item_type.append(xplt_node('ITEM_TYPE',parent=dictionary_item_global[i]))
    item_format.append(xplt_node('ITEM_FORMAT',parent=dictionary_item_global[i]))
    item_unknown.append(xplt_node('ITEM_UNKNOWN',parent=dictionary_item_global[i]))
    item_name.append(xplt_node('ITEM_NAME',parent=dictionary_item_global[i]))
    item_type[-1].set()
    item_name[-1].set()
    item_format[-1].set()

    dictionary_item_nodeset.append(xplt_node('DICTIONARY_ITEM',parent=nodeset_var))
    item_type.append(xplt_node('ITEM_TYPE',parent=dictionary_item_nodeset[i]))
    item_format.append(xplt_node('ITEM_FORMAT',parent=dictionary_item_nodeset[i]))
    item_unknown.append(xplt_node('ITEM_UNKNOWN',parent=dictionary_item_nodeset[i]))
    item_name.append(xplt_node('ITEM_NAME',parent=dictionary_item_nodeset[i]))
    item_type[-1].set()
    item_name[-1].set()
    item_format[-1].set()

    dictionary_item_domain.append(xplt_node('DICTIONARY_ITEM',parent=domain_var))
    item_type.append(xplt_node('ITEM_TYPE',parent=dictionary_item_domain[i]))
    item_format.append(xplt_node('ITEM_FORMAT',parent=dictionary_item_domain[i]))
    item_unknown.append(xplt_node('ITEM_UNKNOWN',parent=dictionary_item_domain[i]))
    item_name.append(xplt_node('ITEM_NAME',parent=dictionary_item_domain[i]))
    item_type[-1].set()
    item_name[-1].set()
    item_format[-1].set()

    dictionary_item_surface.append(xplt_node('DICTIONARY_ITEM',parent=surface_var))
    item_type.append(xplt_node('ITEM_TYPE',parent=dictionary_item_surface[i]))
    item_format.append(xplt_node('ITEM_FORMAT',parent=dictionary_item_surface[i]))
    item_unknown.append(xplt_node('ITEM_UNKNOWN',parent=dictionary_item_surface[i]))
    item_name.append(xplt_node('ITEM_NAME',parent=dictionary_item_surface[i]))
    item_type[-1].set()
    item_name[-1].set()
    item_format[-1].set()

node_section    = xplt_node('NODE_SECTION',parent=mesh)
domain_section  = xplt_node('DOMAIN_SECTION',parent=mesh)
surface_section = xplt_node('SURFACE_SECTION',parent=mesh)
nodeset_section = xplt_node('NODESET_SECTION',parent=mesh)
parts_section   = xplt_node('PARTS_SECTION',parent=mesh)

node_header = xplt_node('NODE_HEADER',parent=node_section)
node_coords = xplt_node('NODE_COORDS',parent=node_section)

domain = xplt_node('DOMAIN',parent=domain_section)

nodes = xplt_node('NODES',parent=node_header)
dim   = xplt_node('DIM',parent=node_header)
name  = xplt_node('NAME',parent=node_header) 

domain_header = xplt_node('DOMAIN_HEADER',parent=domain)
element_list = xplt_node('ELEMENT_LIST',parent=domain)

elem_type = xplt_node('ELEM_TYPE',parent=domain_header)
part_id = xplt_node('PART_ID',parent=domain_header)
elements = xplt_node('ELEMENTS',parent=domain_header)
name = xplt_node('NAME',parent=domain_header) 

state_header = []
state_unknown2 = []
state_data = []
state_time = []
state_unknown = []
global_data = []
node_data = []
domain_data = []
surface_data = []
state_var_node = []
variable_id_node = []
variable_data_node = []
state_var_domain = []
variable_id_domain = []
variable_data_domain = []

for i in range(nstates):
    state_header.append(xplt_node('STATE_HEADER',parent=state_section[i]))
    state_unknown2.append(xplt_node('STATE_UNKNOWN2',parent=state_section[i]))
    state_data.append(xplt_node('STATE_DATA',parent=state_section[i]))

    state_time.append(xplt_node('STATE_TIME',parent=state_header[i]))
    state_unknown.append(xplt_node('STATE_UNKNOWN',parent=state_header[i]))

    global_data.append(xplt_node('GLOBAL_DATA',parent=state_data[i]))
    node_data.append(xplt_node('NODE_DATA',parent=state_data[i]))
    domain_data.append(xplt_node('DOMAIN_DATA',parent=state_data[i]))
    surface_data.append(xplt_node('SURFACE_DATA',parent=state_data[i]))

    for j in range(MAXVAR):
        state_var_node.append(xplt_node('STATE_VAR',parent=node_data[i]))
        variable_id_node.append(xplt_node('VARIABLE_ID',parent=state_var_node[-1]))
        variable_data_node.append(xplt_node('VARIABLE_DATA',parent=state_var_node[-1]))
        variable_data_node[-1].set()
        variable_id_node[-1].set()

        state_var_domain.append(xplt_node('STATE_VAR',parent=domain_data[i]))
        variable_id_domain.append(xplt_node('VARIABLE_ID',parent=state_var_domain[-1]))
        variable_data_domain.append(xplt_node('VARIABLE_DATA',parent=state_var_domain[-1]))
        variable_id_domain[-1].set()
    variable_data_domain[-2].set()
    variable_data_domain[-1].set()

nodeset_var.set()
domain_var.set()
software.set()
nodes.set()
dim.set()
domain.set()
#node_coords.set()
elem_type.set()
elements.set()

for i in range(nstates):
    state_time[i].set()

#for pre,fill,node in RenderTree(feb):
#    treestr = u"%s%s" %(pre,node.name)
#    print(treestr.ljust(8),node.name,node.read)

feb.read_file(fid)

for pre,fill,node in RenderTree(feb):
    treestr = u"%s%s" %(pre,node.name)
    print(treestr.ljust(8),node.name,node.read)

#print([node.name for node in LevelOrderIter(feb,filter_=lambda n: n.depth==3)])
#print([node.name for node in LevelOrderIter(state_section[0],filter_=lambda n: n.depth==variable_data_node[0].depth and n.name==variable_data_node[0].name)])
#for v in variable_data_node:
#    v.set()
#fid.seek(0)
#feb.read_file(fid)
