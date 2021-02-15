import numpy as np
from SMILES import SMILES
from BigSMILES_BigSmilesObj import BigSMILES
from BigSmilesPattern import BigSmilesPattern
from BigSMILES_Bond import BigSMILES_Bond
from BigSMILES_StoObj import BigSMILES_StoObj
from utility import errorMsg, flatten_list
from error import BigSMILESError,BigSMILES_BondInconsistencyError,\
    BigSMILES_StoObjMissingTerminalError, BigSMILES_StoObjError
import networkx as nx
from collections import deque
import csv


### Creating a class to facilitate condition writing later in code
class Atom():
    def __init__(self,atom_type,hybrid,num_H,num_noneH,position):
        self.atom_type = atom_type
        self.hybrid = hybrid
        self.num_H = num_H
        self.position = position
        self.num_noneH = num_noneH


### Generates vector/list of all possible singlets from listed atoms.
def make_singlets(atoms_considered):
    singles = []
    singlist = []

    atoms=['C','c','N','n','O','o','Si','S','F','Cl','Br']
    hybridization = ['sp3', 'sp2', 'sp1', 'sp3d3', '--']
    number_H = [3,2,1,0]
    position = ['bb', 'sg']


    ########## CARBONS
    ##### non aromatic
    if 'C' in atoms_considered:
        for i in range(len(hybridization[:3])):
            for j in range(i,len(number_H)):
                if i == j:
                    singles.append(Atom(atoms[0], hybridization[i],\
                                        number_H[j], 4-i-number_H[j],\
                                        'sg'))
                    singlist.append([atoms[0], hybridization[i],\
                                     number_H[j], 4-i-number_H[j],\
                                     'sg'])
                if i!=j:
                    for p in position:
                        singles.append(Atom(atoms[0], hybridization[i],\
                                            number_H[j], 4-i-number_H[j],\
                                            p))
                        singlist.append([atoms[0], hybridization[i],\
                                         number_H[j], 4-i-number_H[j],\
                                         p])

    ##### aromatic
    if 'c' in atoms_considered:
        for p in position:
            for numH in number_H[2:]:
                singles.append(Atom(atoms[1], hybridization[1],\
                                    numH, 2+(1-numH), p))
                singlist.append([atoms[1], hybridization[1],\
                                 numH, 2+(1-numH), p])



    ########## NITROGEN
    ##### non-aromatic
    if 'N' in atoms_considered:
        for i in range(len(hybridization[:-2])):
            for j in range(i+1,len(number_H)):
                if i+1==j:
                    singles.append(Atom(atoms[2], hybridization[i],\
                                        number_H[j], 3-i-number_H[j],\
                                        'sg'))
                    singlist.append([atoms[2], hybridization[i],\
                                     number_H[j], 3-i-number_H[j],\
                                     'sg'])
                if i+1!=j:
                    for p in position:
                        singles.append(Atom(atoms[2], hybridization[i],\
                                            number_H[j], 3-i-number_H[j],\
                                            p))
                        singlist.append([atoms[2], hybridization[i],\
                                         number_H[j], 3-i-number_H[j],\
                                         p])

    ##### aromatic
    if 'n' in atoms_considered:
        for p in position:
            singles.append(Atom(atoms[3], hybridization[1],\
                                number_H[3], 2, p))
            singlist.append([atoms[3], hybridization[1],\
                             number_H[3], 2, p])


    ########## OXYGENS
    ##### non-aromatic
    if 'O' in atoms_considered:
        for i in range(len(hybridization[:2])):
            for j in range(i+2,len(number_H)):
                if i+2==j:
                    singles.append(Atom(atoms[4], hybridization[i],\
                                        number_H[j], 2-i-number_H[j],\
                                        'sg'))
                    singlist.append([atoms[4], hybridization[i],\
                                     number_H[j], 2-i-number_H[j],\
                                     'sg'])
                if i+2!=j:
                    for p in position:
                        singles.append(Atom(atoms[4], hybridization[i],\
                                            number_H[j], 2-i-number_H[j],\
                                            p))
                        singlist.append([atoms[4], hybridization[i],\
                                         number_H[j], 2-i-number_H[j],\
                                         p])

    ##### aromatic

    ########## SILICON
    ##### non aromatic

    if 'Si' in atoms_considered:
        for i in range(len(hybridization[:-2])):
            for j in range(i,len(number_H)):
                if i == j:
                    singles.append(Atom('Si', hybridization[i],\
                                        number_H[j], 4-i-number_H[j],\
                                        'sg'))
                    singlist.append(['Si', hybridization[i],\
                                     number_H[j], 4-i-number_H[j],\
                                     'sg'])
                if i!=j:
                    for p in position:
                        singles.append(Atom('Si', hybridization[i],\
                                            number_H[j], 4-i-number_H[j],\
                                            p))
                        singlist.append(['Si', hybridization[i],\
                                         number_H[j], 4-i-number_H[j],\
                                         p])



    ########## SULFUR
    if 'S' in atoms_considered:
        singles.append(Atom('S', hybridization[0], number_H[2],\
                            2-number_H[j], 'sg'))
        singlist.append(['S', hybridization[0], number_H[2],\
                         2-number_H[j], 'sg'])
        for p in position:
            singles.append(Atom('S', hybridization[0], number_H[3],\
                                2-number_H[j], p))
            singlist.append(['S', hybridization[0], number_H[3],\
                             2-number_H[j], p])
        for p in position:
            singles.append(Atom('S', 'sp3d3', 0, 4, p))
            singlist.append(['S', 'sp3d3', 0, 4, p])



    ########## HALOGENS
    for a in atoms[-3:]:
        if a in atoms_considered:
            singles.append(Atom(a, hybridization[-1], 0, 1, 'sg'))
            singlist.append([a, hybridization[-1], 0, 1, 'sg'])

    return singles,singlist





##################### GETTING DOUBLES ######################

### From list of sinlgets, generates all chemically possible covalent pair
### Here are listed all the different combinations possible
def make_doublets(singles,singleslist):
    doublets = []
    doubletl = []
    aromatic = ['c', 'n']
    for i in range(len(singles)-3): # no need to go through the halogens at the end
        for j in range(i, len(singles)):
            a = singles[i]
            b = singles[j]
            al = singleslist[i]
            bl = singleslist[j]

            ### condition for single bonding
            condsp2 = (b.hybrid == 'sp2') and\
                ((b.atom_type in aromatic) and (b.num_noneH>2))\
                or ((b.atom_type not in aromatic) and (b.num_noneH>1))
            condsp1 = (b.hybrid == 'sp1') and (b.num_noneH >1)
            condsp3 = (b.hybrid == 'sp3') or (b.hybrid == '--')
            condd3 = (b.hybrid == 'sp3d3')
            condsinglebond = condsp1 or condsp2 or condsp3 or condd3

            if (a.num_noneH != 1) or (b.num_noneH != 1):
                if (a.hybrid == 'sp1'):
                    if (b.hybrid == 'sp1'):
                        doublets.append([a,'#',b])
                        doubletl.append([al,'#',bl])
                    if (a.num_noneH >1) and condsinglebond:
                        doublets.append([a,'-',b])
                        doubletl.append([al,'-',bl])

                if (a.hybrid == 'sp2') and (a.atom_type in aromatic):
                    if (b.atom_type in aromatic):
                        doublets.append([a,':',b])
                        doubletl.append([al,':',bl])
                    if (a.num_noneH == 3):
                        if condsinglebond:
                            doublets.append([a,'-',b])
                            doubletl.append([al,'-',bl])

                if (a.hybrid == 'sp2') and (a.atom_type not in aromatic):
                    if ((b.hybrid == 'sp2') and (b.atom_type not in aromatic))\
                            or (b.hybrid == 'sp3d3'):
                        doublets.append([a,'=',b])
                        doubletl.append([al,'=',bl])
                    if (a.num_noneH >1):
                        if condsinglebond:
                            doublets.append([a,'-',b])
                            doubletl.append([al,'-',bl])

                if (a.hybrid == 'sp3'):
                    if condsinglebond:
                        doublets.append([a,'-',b])
                        doubletl.append([al,'-',bl])

                if (a.hybrid == 'sp3d3'):
                    if (b.hybrid == 'sp2') and (b.atom_type not in aromatic):
                        doublets.append([a,'=',b])
                        doubletl.append([al,'=',bl])
                    if condsinglebond:
                        doublets.append([a,'-',b])
                        doubletl.append([al,'-',bl])

    return doublets, doubletl



##################### GETTING TRIPLETS ######################

### Same idea as for doublets
def make_triplets(singles, singleslist):
    triplets = []
    tripletl = []
    aromatic = ['c', 'n']
    for i in range(len(singles)-3):
        for j in range(len(singles)):
            for k in range(j,len(singles)):
                a = singles[j]
                b = singles[i]
                c = singles[k]
                al = singleslist[j]
                bl = singleslist[i]
                cl = singleslist[k]

                ### condition single bonding for a
                acondsp2 = (a.hybrid == 'sp2') and\
                    ((a.atom_type in aromatic) and (a.num_noneH>2)) or\
                    ((a.atom_type not in aromatic) and (a.num_noneH>1))
                acondsp1 = (a.hybrid == 'sp1') and (a.num_noneH >1)
                acondsp3 = (a.hybrid == 'sp3') or (a.hybrid == '--')
                acondd3 = (b.hybrid == 'sp3d3')
                asinglebond = acondsp1 or acondsp2 or acondsp3 or acondd3
                ### condition single bonding for c
                ccondsp2 = (c.hybrid == 'sp2') and\
                    ((c.atom_type in aromatic) and (c.num_noneH>2)) or\
                    ((c.atom_type not in aromatic) and (c.num_noneH>1))
                ccondsp1 = (c.hybrid == 'sp1') and (c.num_noneH >1)
                ccondsp3 = (c.hybrid == 'sp3') or (c.hybrid == '--')
                ccondd3 = (b.hybrid == 'sp3d3')
                csinglebond = ccondsp1 or ccondsp2 or ccondsp3 or ccondd3

                ### basic conditions for triplet to exist
                cond1 = (a.position =='sg') and (c.position == 'sg') and\
                    (b.position == 'bb') and (b.num_noneH >2)
                cond2 = (b.num_noneH >1) and\
                    ((b.position == a.position) or (b.position == c.position))

                if cond1 or cond2:

                    if (b.hybrid == 'sp1'):
                        if (a.hybrid == 'sp1'):
                            if csinglebond:
                                triplets.append([a,'#',b,'-',c])
                                tripletl.append([al,'#',bl,'-',cl])
                        if (c.hybrid == 'sp1') and (c != a):
                            if asinglebond:
                                triplets.append([a,'-',b,'#',c])
                                tripletl.append([al,'-',bl,'#',cl])

                    if (b.atom_type in aromatic):
                        if (a.atom_type in aromatic):
                            if (c.atom_type in aromatic):
                                triplets.append([a,':',b,':',c])
                                tripletl.append([al,':',bl,':',cl])
                            if (b.num_noneH ==3):
                                if csinglebond:
                                    triplets.append([a,':',b,'-',c])
                                    tripletl.append([al,':',bl,'-',cl])
                        if (c.atom_type in aromatic) and (c != a):
                            if (b.num_noneH ==3):
                                if asinglebond:
                                    triplets.append([a,'-',b,':',c])
                                    tripletl.append([al,'-',bl,':',cl])

                    if (b.hybrid == 'sp2') and (b.atom_type not in aromatic):
                        if ((a.hybrid == 'sp2') and\
                                (a.atom_type not in aromatic)) or\
                                acondd3:
                            if csinglebond:
                                triplets.append([a,'=',b,'-',c])
                                tripletl.append([al,'=',bl,'-',cl])
                        if (c != a):
                            if ((c.hybrid == 'sp2') and\
                                    (c.atom_type not in aromatic)) or\
                                    acondd3:
                                if asinglebond:
                                    triplets.append([a,'-',b,'=',c])
                                    tripletl.append([al,'-',bl,'=',cl])
                        if (b.num_noneH ==3):
                            if asinglebond and csinglebond:
                                triplets.append([a,'-',b,'-',c])
                                tripletl.append([al,'-',bl,'-',cl])

                    if (b.hybrid == 'sp3'):
                        if asinglebond and csinglebond:
                            triplets.append([a,'-',b,'-',c])
                            tripletl.append([al,'-',bl,'-',cl])

                    if (b.hybrid == 'sp3d3'):
                        if ((a.hybrid == 'sp2') and\
                                (a.atom_type not in aromatic)) or\
                                acondd3:
                            if ((c.hybrid == 'sp2') and\
                                    (c.atom_type not in aromatic)) or\
                                    acondd3:
                                triplets.append([a,'=',b,'=',c])
                                tripletl.append([al,'=',bl,'=',cl])
                            if csinglebond:
                                triplets.append([a,'=',b,'-',c])
                                tripletl.append([al,'=',bl,'-',cl])
                        if ((c.hybrid == 'sp2') and\
                                (c.atom_type not in aromatic)) or\
                                acondd3:
                            if asinglebond:
                                triplets.append([a,'-',b,'=',c])
                                tripletl.append([al,'-',bl,'=',cl])
                        if asinglebond and csinglebond:
                            triplets.append([a,'-',b,'-',c])
                            tripletl.append([al,'-',bl,'-',cl])

    return triplets, tripletl




##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################


###### GENERATING THE LIST OF TUPLETS OF A REPEAT UNIT




aromatic = ['c', 'n']
types_atom = ['', 'bracket_atom', 'organic']
types_bond = ['', '/', '\\']


### lists all the hibridity conditions possible per atom type
def find_hybridization(atom, bonds):
    if atom in ['F', 'Cl', 'Br']:
        return '--'
    if (atom == 'S') and (len(bonds)>2):
        return 'sp3d3'
    else:
        if ('=' in bonds) or (atom in aromatic):
            return 'sp2'
        if '#' in bonds:
            return 'sp1'
        if ('=' not in bonds) and ('#' not in bonds):
            return 'sp3'


### finds number of H neighbors from hybridization and # of neighbors from BigSMILES string
def find_numH(atom,hybridization,noneH):
    if atom in ['C','Si','c']:
        return int(hybridization[-1])+1 - noneH
    if hybridization == 'sp3d3':
        return 4 - noneH
    if atom in ['O', 'S']:
        return int(hybridization[-1])-1 - noneH
    if atom in ['N','n']:
        return int(hybridization[-1]) - noneH
    if atom in ['F', 'Cl', 'Br']:
        return 0


### Get all the details of a single atom
def find_singlet(graphe, backbone,  node_index):
    node = graphe.nodes[node_index]
    atom = node['atom']
    noneH = 0
    for i in node['neighList']:
        if graphe.nodes[i]['atom'] != 'H':
            noneH +=1
    bond_list = []
    for i in node['neighList']:
        bond_list.append(graphe.edges[(node_index,i)]['type'])
    hybrid = find_hybridization(atom,bond_list)
    numH = find_numH(atom,hybrid,noneH)
    if node_index in backbone:
        position = 'bb'
    else:
        position = 'sg'
    return [atom,hybrid,numH,noneH,position]



### From the large vector of all possible singlets, doublets and triplets
### generated at the begining, find the index that corresponds to the
### tuplet we are looking at
def find_index(description, listlet):
    i=0
    while description != listlet[i]:
        i+=1
    return i

### Generate list of singlets in a repeat unit
def get_all_singlets(graphe,singlist):
    node = graphe.nodes()
    Gsingles = []
    backbone = find_backbone(graphe)
    for i in range(1,1+len(graphe)):
        if node[i]['_type'] in types_atom:
            description = find_singlet(graphe,backbone,i)
            Gsingles.append([description,find_index(description,singlist)])
    return Gsingles

### Determins if atom is along backbone or sidegroup
def find_backbone(graphe):
    types_atom = ['', 'bracket_atom', 'organic']
    index = []
    for i in graphe:
        if graphe.nodes[i]['_type'] not in types_atom:
            index.append(i)
    return nx.shortest_path(graphe,source=index[0],target=index[1])

### Finds all the doublets in a repeat unit
def get_all_doublets(graphe,singlist,doublist):
    node = graphe.nodes()
    edge = graphe.edges()
    Gdoubles = []
    Gconnect = []
    backbone = find_backbone(graphe)
    for (i,j) in edge():
        if  node[i]['_type'] in types_atom:
            a = find_singlet(graphe, backbone,i)
            aind = find_index(a,singlist)
        if  node[i]['_type'] not in types_atom:
            a = node[i]
            aind = -1
        if  node[j]['_type'] in types_atom:
            b = find_singlet(graphe,backbone, j)
            bind = find_index(b,singlist)
        if  node[j]['_type'] not in types_atom:
            b = node[j]
            bind = -1
        if aind<0:
            Gconnect.append([a,b])
        if bind<0:
            Gconnect.append([b,a])
        if (aind>=0) and (bind>=0):
            bond = edge[i,j]['type']
            if (a[0] in aromatic) and (b[0] in aromatic):
                bond = ':'
            if bond in types_bond:
                bond = '-'
            if aind>bind:
                double = [b,bond,a]
            else:
                double = [a,bond,b]
            Gdoubles.append([double,find_index(double,doublist)])
    return Gdoubles, Gconnect

### Previous code generates doublets such as ['<', 'C sp3 ...]
### This function takes in all the doublets that contain connectors
### such as $, < or > and creats all the connections possible
### and return the list of these generated doublets
def make_double_connect(connectlist,singlist,doublist):
    Gdoubles = []
    for i in range(len(connectlist)-1):
        a = connectlist[i]
        abond_num = a[0]['BigSMILES_Bond']._id
        asymb = a[0]['rawStr']
        if '>' in asymb:
            symb = '<'
        if '<' in asymb:
            symb = '>'
        if '$' in asymb:
            symb = '$'
        for j in range(i+1,len(connectlist)):
            b = connectlist[j]
            bbond_num = b[0]['BigSMILES_Bond']._id
            bsymb = b[0]['rawStr']
            if (symb in bsymb) and (bbond_num == abond_num):
                aind = find_index(a[1],singlist)
                bind = find_index(b[1],singlist)
                if (a[1][0] in aromatic) and (b[1][0] in aromatic):
                    if aind<bind:
                        double = [a[1],':',b[1]]
                    else:
                        double = [b[1],':',a[1]]
                else:
                    if aind<bind:
                        double = [a[1],'-',b[1]]
                    else:
                        double = [b[1],'-',a[1]]
                Gdoubles.append([double,find_index(double,doublist)])
    return Gdoubles

### Generates all triplets
def get_all_triplets(graphe,singlist,triplist):
    Gtriplets = []
    Gconnect = []
    Gtwiceconnect = []
    node = graphe.nodes()
    edge = graphe.edges()
    backbone = find_backbone(graphe)
    for i in node:
        if node[i]['_type'] in types_atom:
            neigh = node[i]['neighList']
            if len(neigh)>1:
                b = find_singlet(graphe,backbone,i)
                for j in range(len(neigh)-1):
                    a = node[neigh[j]]
                    aind = -1
                    if a['_type'] in types_atom:
                        a = find_singlet(graphe,backbone,neigh[j])
                        aind = find_index(a,singlist)
                        abond = edge[i,neigh[j]]['type']
                        if (a[0] in aromatic) and (b[0] in aromatic):
                            abond = ':'
                        if abond in types_bond:
                            abond = '-'
                    for k in range(j+1, len(neigh)):
                        c = node[neigh[k]]
                        cind = -1
                        if c['_type'] in types_atom:
                            c = find_singlet(graphe,backbone,neigh[k])
                            cind = find_index(c,singlist)
                            cbond = edge[i,neigh[k]]['type']
                            if (c[0] in aromatic) and (b[0] in aromatic):
                                cbond = ':'
                            if cbond in types_bond:
                                cbond = '-'
                        if (aind<0) and (cind>=0):
                            Gconnect.append([a,[b,cbond,c]])
                        if (cind<0) and (aind>=0):
                            Gconnect.append([c,[b,abond,a]])
                        if (aind<0) and (cind<0):
                            Gtwiceconnect.append([a,b,c])
                        if (aind>=0) and (cind>=0):
                            if aind<cind:
                                triplet = [a,abond,b,cbond,c]
                            else:
                                triplet = [c,cbond,b,abond,a]
                            Gtriplets.append([triplet,find_index(triplet,triplist)])
    return Gtriplets, Gconnect, Gtwiceconnect

### As before, the previous funciton can generate triplets such as ['$', 'C', 'C']
### This function return a list of triplets with all possible connecitons, uses the
### previously generated list of doublets to do so.
def make_triple_connect(connectTriples,connectDoubles,singlist,triplist):
    Gtriples = []
    for i in range(len(connectTriples)):
        a = connectTriples[i]
        abond_num = a[0]['BigSMILES_Bond']._id
        asymb = a[0]['rawStr']
        if '>' in asymb:
            symb = '<'
        if '<' in asymb:
            symb = '>'
        if '$' in asymb:
            symb = '$'
        for j in range(len(connectDoubles)):
            b = connectDoubles[j]
            bbond_num = b[0]['BigSMILES_Bond']._id
            bsymb = b[0]['rawStr']
            if (symb in bsymb) and (abond_num == bbond_num) :
                aind = find_index(a[1][2],singlist)
                bind = find_index(b[1],singlist)
                if (a[1][0][0] in aromatic) and (b[1][0] in aromatic):
                    if aind<bind:
                        triple = [a[1][2],a[1][1],a[1][0],':',b[1]]
                    else:
                        triple = [b[1],':',a[1][0],a[1][1],a[1][2]]
                else:
                    if aind<bind:
                        triple = [a[1][2],a[1][1],a[1][0],'-',b[1]]
                    else:
                        triple = [b[1],'-',a[1][0],a[1][1],a[1][2]]

                Gtriples.append([triple,find_index(triple,triplist)])
    return Gtriples

### This is in case a triplet such as ['$', 'C', '$'] is genreated, does as above
def make_triple_twiceconnect(twiceconnectTriples,connectDoubles,singlist,triplist):
    Gtripleconnect = []
    for i in range(len(twiceconnectTriples)):
        a = twiceconnectTriples[i]
        abond_num = a[0]['BigSMILES_Bond']._id
        asymb = a[0]['rawStr']
        if '>' in asymb:
            symb = '<'
        if '<' in asymb:
            symb = '>'
        if '$' in asymb:
            symb = '$'
        for j in range(len(connectDoubles)):
            b = connectDoubles[j]
            bbond_num = b[0]['BigSMILES_Bond']._id
            bsymb = b[0]['rawStr']
            if (symb in bsymb) and (abond_num == bbond_num) :
                if (a[1][0] in aromatic) and (b[1][0] in aromatic):
                    triple = [a[2],[a[1],':',b[1]]]
                else:
                    triple = [a[2],[a[1],'-',b[1]]]
            Gtripleconnect.append(triple)
    return make_triple_connect(Gtripleconnect,connectDoubles,singlist,triplist)



### Returns the singlets, doublets and triplets vectors from a BigSMILES string
def get_vectorization(string,singlist,doublist,triplist):
    polymer = BigSMILES(string)
    repUnits = polymer[0]
    allSingles = []
    allDoubles = []
    connectDoubles = []
    allTriples = []
    connectTriples = []
    twiceconnectTriples = []
    for unit in repUnits:
        allSingles = allSingles+get_all_singlets(unit.G,singlist)
        Gdoubles, Gconnect = get_all_doublets(unit.G,singlist,doublist)
        allDoubles = allDoubles + Gdoubles
        connectDoubles = connectDoubles + Gconnect
        Gtriples, Gconnect, Gtwiceconnect = get_all_triplets(unit.G,singlist,triplist)
        allTriples = allTriples + Gtriples
        connectTriples = connectTriples + Gconnect
        twiceconnectTriples = twiceconnectTriples + Gtwiceconnect
    allDoubles = allDoubles + make_double_connect(connectDoubles,singlist,doublist)
    allTriples = allTriples + make_triple_connect(connectTriples,connectDoubles,singlist,triplist)
    allTriples = allTriples + make_triple_twiceconnect(twiceconnectTriples,connectDoubles,singlist,triplist)
    return allSingles, allDoubles, allTriples



##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################


### from the list of tuplets, generate the vector that simply adds +1 at the
### index corresponding to that tuplet in the giant all possible tuplet vector
def make_vector(vector,vectorlist):
    vect = np.zeros(len(vectorlist))
    for i in vector:
        vect[i[-1]] +=1
    return vect


### Takes away all tuplets that doesn't appear in our dataset
def reduce_dim(feature_vector):
    shape = feature_vector.shape
    feature_index = np.arange(0,shape[0])
    zero = np.zeros(shape[1])
    for i in range(shape[0]):
        if np.all(feature_vector[i,:] == zero):
            feature_index[i] = -1
    indices = np.where(feature_index>-1)[0]
    return indices

### enter a csv file with the bigSMILES strings listed, and
### the directory of that file (ex: './' = this folder, or
### /Users/yourname/Desktop/ ...)
### Must end wiht a / !!!!
### This will generate files in your directory with the vectors and their corresponding tuplet
def vectorize(filename,directory):
    atoms=['C','c','N','n','O','o','Si','S','F','Cl','Br']
    dataset = filename[:-4]
    singlets, singlist = make_singlets(atoms)
    doublets, doublist = make_doublets(singlets, singlist)
    triplets, triplist = make_triplets(singlets, singlist)

    with open(directory+filename, newline='') as f:
        reader = csv.reader(f)
        data = list(reader)

    string = data[0][0]
    print(string)
    singles,doubles,triples = get_vectorization(string,singlist,doublist,triplist)
    Singles = np.array([make_vector(singles,singlist)]).T/len(singles)
    Doubles = np.array([make_vector(doubles,doublist)]).T/len(singles)
    Triples = np.array([make_vector(triples,triplist)]).T/len(singles)


    for i in range(1,len(data)):
        string = data[i][0]
        print(string)
        singles,doubles,triples = get_vectorization(string,singlist,doublist,triplist)
        Single = np.array([make_vector(singles,singlist)]).T/len(singles)
        Double = np.array([make_vector(doubles,doublist)]).T/len(singles)
        Triple = np.array([make_vector(triples,triplist)]).T/len(singles)

        Singles = np.concatenate((Singles,Single),axis=1)
        Doubles = np.concatenate((Doubles,Double),axis=1)
        Triples = np.concatenate((Triples,Triple),axis=1)

    singind = reduce_dim(Singles)
    doubind = reduce_dim(Doubles)
    tripind = reduce_dim(Triples)
    Singredu = Singles[singind]
    Doubredu = Doubles[doubind]
    Tripredu = Triples[tripind]

    singlistredu=[]
    doublistredu=[]
    triplistredu=[]
    listredu = [singlistredu,doublistredu,triplistredu]
    lists = [singlist,doublist,triplist]
    indices = [singind,doubind,tripind]
    vector = ['singles','doubles','triples']
    for n in range(3):
        for i in indices[n]:
            listredu[n].append([lists[n][i], i])
        with open(directory+dataset+'_reduced_list_'+vector[n]+'.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerows(listredu[n])

    full_vectors = [Singles,Doubles,Triples]
    reduced_vectors = [Singredu,Doubredu,Tripredu]
    for n in range(3):
 #       np.savetxt(directory+dataset+'_full_vector_'+vector[n]+'.csv',full_vectors[n], delimiter=',')
        np.savetxt(directory+dataset+'_reduced_vector_'+vector[n]+'.csv',reduced_vectors[n], delimiter=',')

### Example
#vectorize('BIGsmiles_Dataset.csv', './')

