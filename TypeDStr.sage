#Need to implement:
#simplify in class TypeDStr
#def infty_type_D_bd(k):
#def zero_type_D_bd(k):
#Improve mor_to_d
#DElt.__hash__

def paths(graph):
    "Returns a list of all oriented paths in the directed graph. Raises exception if graph has loop."
    if graph.has_loops():
        raise Exception("Graph has directed loops; set of paths not finite.")
    answer = list()
    for v in graph.vertices():
        answer.extend(paths_from_v(graph,v))
    return answer

def paths_from_v(graph,v):
    "Returns a list of all oriented paths in the directed graph starting at v. Does not check this is finie."
    answer = list()
    for w in graph.neighbors_out(v):
        for p in paths_from_v(graph,w):
            answer.append([v]+p)
    answer.append([v])
    return answer

def longest_path_len(graph):
    "Returns the length of the longest path in directed graph. Raises exception if graph has directed loop.)"
    if graph.has_loops():
        raise Exception("Graph has directed loops; set of paths not finite.")
    path_list = paths(graph)
    l=0
    for x in path_list:
        if len(x)>l:
            l=len(x)
    return l-1

class DGen(object):
    """A class for generators of a type D structure.

    Inputs:
    name: name for this generator.
    pmc: underlying pointed matched circle.
    idem: idempotent for this generator.

    Both pmc and idem should be given EXCEPT if just want to obtain an element with a given name.

    __eq__ JUST LOOKS AT THE NAME.  Similarly __hash__.

    Example:
    x = DGen('x',PMC.split_matching(1),[(0,2)])
    """
    def __init__(self, name, pmc=[], idem=[]):
        #Should do some sanity checking.
        self.name = name
        self.pmc = pmc
        self.idem =idem
        self.idem.sort()

    def __repr__(self):
        return self.name.__repr__()
    
    def __eq__(self, other):
        "Checks if the NAMES of self and other are the same."
        if type(self) != type(other):
            return False
        if self.name == other.name:
            return True
        return False

    def __ne__(self, other):
        return not (self.__eq__(other))

    def __hash__(self):
        return self.name.__hash__()

    def __rmul__(self, other):
        "Multiply self by alg elt other to give a DElt."
        coeff = AlgElt([Strand_Diagram(self.pmc, [], left_idem=self.idem)])
        coeff = other*coeff
        return DElt({self:coeff})

from UserDict import UserDict
class DElt(UserDict):
    """A class for elements of a left Type D structure.

    Input:
    A dict 'data' with keys x_i DGen's and values algebra elements a_i.
    Alternatively, a_i may be a list of Strand_Diagram's. Also, it's okay for x_i just to be a name, not a DGen: in that case, the idempotent will be extracted from a_i[0].
    pmc: optional argument, giving the pointed matched circle for this element. Must be given if dict is empty.

    Examples:
    small_pmc = PMC.split_matching(1)
    x = DGen('x',small_pmc,[(0,2)])
    y = DGen('y',small_pmc,[(1,3)])
    rho1=Strand_Diagram(small_pmc,[(0,1)],left_idem=[(0,2)],name='rho1')
    rho2=Strand_Diagram(small_pmc, [(1,2)],left_idem=[(1,3)],name='rho2')
    rho3=Strand_Diagram(small_pmc, [(2,3)],left_idem=[(0,2)],name='rho3')
    ax = DElt({x:AlgElt([rho2])})
    by = DElt({y:AlgElt([rho1,rho3])})
    bigelt = DElt({'x':[rho2],'y':[rho1,rho3]})
    """
    def __init__(self, data, pmc=[]):
        gooddata = dict()
        for x in data.keys():
            if isinstance(x,DGen):
                goodx = x
            else:
                goodx = DGen(x,data[x][0].pmc,data[x][0].right_idem)
            if isinstance(data[x],AlgElt):
                gooddata[goodx]=data[x]
            else:
                gooddata[goodx]=AlgElt(data[x])
        UserDict.__init__(self, gooddata)
        if pmc:
            self.pmc = pmc
        else:
            self.pmc = self.data.values()[0].pmc
        self.reduce()

#    def keys(self):
#        return self.data.keys()

    def coeff(self, x):
        return self.data[x]

    def diff_coeffs(self):
        "Returns the element obtained by differentiating the coefficients of self."
        answer = dict()
        for x in self.data.keys():
            answer[x]=list()
            for a in self.data[x]:
                answer[x].extend(a.differential())
        return DElt(answer, pmc=self.pmc)

    def basis_expansion(self):
        "Returns a list of pairs [(rho_i,x_i)], with rho_i a Strand_Diagram and x_i a DGen so that self is sum_i rho_i \otimes x_i."
        self.reduce()
        answer = list()
        for x in self.data.keys():
            for a in self.data[x]:
                answer.append((a,x))
        return answer

    def __repr__(self):
        return repr(self.data)

    def opposite(self):
        "Return result of taking opposite() of algebra elements."
        revdata = dict()
        for x in self.data.keys():
            revdata[x]=self.data[x].opposite()
        return DElt(revdata, self.pmc.opposite())

    def reduce(self):
        "Delete duplicates from self, and also keys with values the empty list."
        for x in self.data.keys():
            self.data[x].reduce_mod_2()
        for x in self.data.keys():
            if self.data[x]==[]:
                del self.data[x]

    def __add__(self, other):
        if self.pmc != other.pmc:
            raise Exception("Can't add elements over different pointed matched circles.")
        ansdict = dict()
        for x in self.data.keys():
            if not (x in other.data.keys()):
                ansdict[x]=self.data[x]
            else:
                ansdict[x] = self.data[x]+other.data[x]
        for x in other.data.keys():
            if not (x in self.data.keys()):
                ansdict[x]=other.data[x]
        return DElt(ansdict, self.pmc)
    
    def __rmul__(self, other):
        "Left multiply an element by an algebra element."
        ansdict = dict()
        for x in self.data.keys():
            ansdict[x] = other*self.data[x]
        return DElt(ansdict, self.pmc)

    def __eq__(self, other):
        if type(self) != type(other):
            return False
        if Set(self.data.keys()) != Set(other.data.keys()):
            return False
        for x in self.data.keys():
            if self.data[x] != other.data[x]:
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

#    def __hash__(self):
#        pass

class TypeDStr(object):
    """A type D structure.

    Inputs: 
    --pmc: the underlying pointed matched circle.
    --basis: a dictionary with keys names for basis elements (arbitrary) and values the corresponding idempotents in pmc 
    OR 
    a list of DGen's
    --structure_consts: a dictionary with keys the basis elements and values lists of pairs (a,x) where:
    -----a is an algebra element (or, alternatively, list of ChordDiagrams for pmc).
    -----x is a basis element.
    OR
    a dictionary with keys basis elements and values DElt's

    Currently, to give structure_consts in first form, must also give basis in first form.

    #OLD:
    Example:
    my_pmc = PMC.split_matching(1)
    rho2=Strand_Diagram(my_pmc, [(1,2)],left_idem=[(1,3)])
    rho3=Strand_Diagram(my_pmc, [(2,3)],left_idem=[(0,2)])
    i02=Strand_Diagram(my_pmc, [], left_idem=[(0,2)])
    gens = {'x':[(1,3)],'a':[(0,2)], 'b':[(0,2)]}
    str_consts={'x':[([rho2],'a')], 'a':[], 'b':[([rho3],'x'),([i02],'a')]}
    CFD = TypeDStr(my_pmc, basis=gens, structure_consts=str_consts)
    CFD
    CFD.display()
    CFD.is_bounded()
    
    rho23=Strand_Diagram(my_pmc, [(1,3)], left_idem=[(1,3)])
    new_gens = {'x':[(1,3)]}
    new_str_consts={'x':[([rho23],'x')]}
    new_CFD = TypeDStr(my_pmc,basis=new_gens, structure_consts=new_str_consts)
    """

    def __init__(self, pmc, basis, structure_consts={}):
        self.pmc=pmc
        #Coerce first form of basis to second form.
        if type(basis) == type({0:1,2:3}):
            self.basis = list()
            for x in basis.keys():
                self.basis.append(DGen(x,self.pmc,basis[x]))
        else:
            if type(basis) == type([1,2,3,4]):
                self.basis = basis
            else:
                raise Exception("Basis must be either a dict or a list.")
        #Coerce first form of structure_consts to second.
        self.structure_consts = dict()
        for x in structure_consts.keys():
            if type(structure_consts[x])==type([1,2,3,4]):
                dx = DElt({},self.pmc)
                for (a,y) in structure_consts[x]:
                    dx = dx+DElt({y:a})
                if type(basis)==type({0:1,2:3}):
                    self.structure_consts[DGen(x,self.pmc,basis[x])]=dx
                else:
                    self.structure_consts[x]=dx
            else:
                if type(basis)==type({0:1,2:3}):
                    self.structure_consts[DGen(x,self.pmc,basis[x])]=structure_consts[x]
                else:
                    self.structure_consts[x]=structure_consts[x]
        #Currently does almost no sanity checking; improve.

    def __repr__(self):
        return 'Type D structure with generators '+repr(self.basis)+'\nDelta1 given by '+repr(self.structure_consts)

    def graph(self):
        "Returns a graph representing self."
        verts = self.basis
        edges = dict()
        for v in verts:
            targets = dict()
            for e in self.structure_consts[v].keys():
                targets[repr(e)]=repr(self.structure_consts[v].coeff(e))
#                targets[e]=repr(self.structure_consts[v].coeff(e))
#            edges[v]=targets #Fails with prev two lines for mysterious reasons -- something about __trunc__
            edges[repr(v)]=targets
        G=DiGraph(edges)
        return G


    def show(self):
        self.graph().show(edge_labels=True)

    def is_bounded(self):
        "Returns True if self is bounded and False otherwise."
        return (not self.graph().has_loops())

    def bounding_num(self):
        "Returns the highest n so that deltan is nonzero."
        return longest_path_len(self.graph())

    def delta1(self, elt):
        if elt in self.structure_consts.keys():
            return self.structure_consts[elt]
        if DGen(elt) in self.structure_consts.keys():
            return self.structure_consts[DGen(elt)]

    def olddelta1(self, elt):
        "Returns delta_1(x), as a list (sum) of tuples (a,y)."
        newdel = self.delta1(elt)
        answer = list()
        for x in newdel.keys():
            for a in newdel.coeffs(x):
                answer.append((a,x))
        return answer

    def deltan(self, elt, n):
        "Returns delta_n(x), as a list (sum) of tuples (a_1,a_2,...,a_n,y)."
        if n==0:
            return [([],elt)]
        if n==1:
            return self.olddelta1(elt)
        else:
            almost = self.olddeltan(elt, n-1)
            answer = list()
            for s in almost:
                for t in self.olddelta1(s[n-1]):
                    answer.append(s[:n-1]+t)
        return answer

    def delta(self, elt):
        if not self.is_bounded():
            raise Exception("Type D structure not bounded, so delta not a finite sum.")
        answer = list()
        for i in range(self.bounding_num()+1):
            answer.extend(self.deltan(elt,i))
        return answer

    def simplify(self, keep_bounded=False):
        """Cancel out naked differentias in delta_1. If keep_bounded is True, do this only as much as preserves boundedness.

        Note: keep_bounded option not yet implemented.
        """
        while self.cancel_an_edge(keep_bounded):
            pass

    def cancel_an_edge(self, keep_bounded=False):
        """Find an edge and cancel it. Return True if an edge was found, False if not.

        Note: keep_bounded option not yet implemented.
        """
        for x in self.basis:
            for y in self.structure_consts[x].keys():
                if self.structure_consts[x][y]*self.structure_consts[x][y]==self.structure_consts[x][y]:
                    self.cancel_edge(x,y)
                    return True
        return False

    def cancel_edge(self, source, target):
        """Cancel the edge from source to tagert.

        Should be able to do this if label is a unit, but currently can only cancel if label is an idempotent.
        """
        if self.delta1(source)[target] != self.delta1(source)[target]*self.delta1(source)[target]:
            raise Exception("Can only cancel an edge labeled by an idempotent.")
        #Find the edges into target
        updated_consts = dict(self.structure_consts)
        del updated_consts[source]
        del updated_consts[target]
        otherverts = list(self.basis)
        otherverts.remove(source)
        otherverts.remove(target)
        for x in otherverts:
            if target in self.structure_consts[x].keys():
                updated_consts[x]=updated_consts[x]+self.structure_consts[x][target]*self.structure_consts[source]
        #Finally, delete all references to edges we've deleted
        for x in otherverts:
            if source in updated_consts[x].keys():
                del updated_consts[x][source]
            if target in updated_consts[x].keys():
                del updated_consts[x][target]
        self.basis = otherverts
        self.structure_consts=updated_consts

    def d_sq_zero(self):
        """Returns True is self satisfies d^2=0 and False otherwise.

        Example:
        sage: small_pmc = PMC.split_matching(1)
        sage: rho2=Strand_Diagram(small_pmc, [(1,2)],left_idem=[(1,3)])
        sage: rho3=Strand_Diagram(small_pmc, [(2,3)],left_idem=[(0,2)])
        sage: iota02=Strand_Diagram(small_pmc, [], left_idem=[(0,2)])
        sage: gens = {'x':[(1,3)],'a':[(0,2)], 'b':[(0,2)]}
        sage: bad_str_consts={'x':[([rho2],'a')], 'a':[], 'b':[([rho1],'x'),([iota02],'a')]}
        sage: bad_CFD = TypeDStr(small_pmc, basis=gens, structure_consts=bad_str_consts)
        sage: bad_CFD.d_sq_zero()
        False
        sage: my_pmc = PMC.split_matching(2)
        sage: newbad_consts={'x':[([Strand_Diagram(my_pmc,[(0,6),(1,3)],[(0,2),(1,3)])],'x')]}
        sage: newbad_CFD = TypeDStr(my_pmc, basis={'x':[(0,2),(1,3)]},structure_consts=newbad_consts)
        """
        #Hasn't been tested much; probably doesn't work quite right.
        for x in self.basis:
            first_part = self.delta1(x).diff_coeffs()
            second_part = DElt({},self.pmc)
            for y in self.delta1(x).keys():
                second_part = second_part + self.delta1(x).coeff(y)*self.delta1(y)
            if first_part + second_part != DElt({},self.pmc):
                return False
        return True

    def shorten_names(self):
        "Replace the names in basis by short ones."
        newbasis = list()
        old_to_new = dict()
        new_to_old = dict()
        for i in range(len(self.basis)):
            ithnewelt = DGen(i,pmc=self.pmc,idem=self.basis[i].idem)
            newbasis.append(ithnewelt)
            old_to_new[self.basis[i]]=ithnewelt
            new_to_old[ithnewelt]=self.basis[i]
        new_str_consts = dict()
        for x in self.structure_consts.keys():
            new_dx_data = dict()
            for y in self.structure_consts[x].keys():
                new_dx_data[old_to_new[y]]=self.structure_consts[x][y]
            new_str_consts[old_to_new[x]]=DElt(new_dx_data,self.pmc)
        self.basis = newbasis
        self.structure_consts = new_str_consts

    def mor_to_d(self, other):
        """Returns the chain complex of homomorphisms from self to other, where other is a type D structure.
        
        Examples:
        sage: infty_type_D(1).mor_to_d(infty_type_D(1))
        Digraph on 2 vertices
        sage: infty_type_D(2).mor_to_d(m_one_type_D(2))
        Digraph on 9 vertices
        sage: infty_type_D(2).mor_to_d(m_one_type_D(2)).homology()
        1
        """
        if self.pmc != other.pmc:
            raise Exception("Can only compute Mor's between type D structures over the same algebra.")
        gens = list()
        alg_gens = self.pmc.alg_basis()
        for x in self.basis:
            #Next loop quite inefficient.
            for a in alg_gens:
                for y in other.basis:
                    if (x.idem==a.left_idem) and (y.idem==a.right_idem):
                        gens.append((x,a,y))
                        #The tuple (x,a,y) stands for the map x -> ay, z->0 if z != x.

        diffs = dict()
        for f in gens:
            df = list()
            #d(f(x))
            #First differentiate y in (x,a,y)
            deltay = other.delta1(f[2])
            ady = AlgElt([f[1]])*deltay
            for (b,z) in ady.basis_expansion():
                df.append((f[0],b,z))
            #Now differentiate a.
            for b in f[1].differential():
                df.append((f[0],b,f[2]))
            #f(dx)
            for xp in self.basis:
                for (b,w) in self.delta1(xp).basis_expansion():
                    if w==f[0]:
                        if b*f[1]:
                            df.append((xp,b*f[1],f[2]))
            diffs[f]=reduce_mod_2(df)
        return ChainCx(diffs)

    def hom_to_d(self, other):
        "Returns a basis for the module of type D structure homomorphisms from self to other."
        pass

