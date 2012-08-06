#Note that TypeDDStr is a left-left module.


class DDGen(object):
    """A class for generators of a type DD structure.

    Inputs:
    name: name for this generator. Probably needs to be hashable, and unique.
    pmc_1: underlying first pointed matched circle.
    pmc_2: underlying second pointed matched circle. Defaults to pmc_1 if not given.
    idem_1: first idempotent for this generator.
    idem_2: second idempotent for this generator.

    All data should be given EXCEPT if just want to obtain an element with a given name.

    __eq__ JUST LOOKS AT THE NAME.  Similarly __hash__.

    Example:
    x = DDGen('x',split_matching(2), split_matching(2), [(0,2),(1,3)],[(0,2),(4,6)])
    """
    def __init__(self, name, pmc_1=[], pmc_2=[], idem_1=[], idem_2=[]):
        #Should do some sanity checking.
        self.name = name
        self.pmc_1 = pmc_1
        if not pmc_2:
            self.pmc_2=pmc_1
        else:
            self.pmc_2=pmc_2
        self.idem_1=idem_1
        self.idem_2=idem_2
        self.idem_1.sort()
        self.idem_2.sort()

    def __repr__(self):
        return self.name.__repr__()
    
    def __eq__(self, other):
        "Checks if the NAMES of self and other are the same."
        if not isinstance(other, DDGen):
            return False
        if self.name == other.name:
            return True
        return False

    def __ne__(self, other):
        return not (self.__eq__(other))

    def __hash__(self):
        return self.name.__hash__()

    def exchangeLR(self):
        "Returns result of exchanging pmc_1 and pmc_2, i.e., a\otimes b becomes b\otimes a."
        return DDGen(self.name, self.pmc_2, self.pmc_1, self.idem_2, self.idem_1)

class DDElt(UserDict):
    """A class for elements of a left Type DD structure.

    Input:
    A dict 'data' with keys x_i DDGen's and values AlgBlgElt's (the coefficients of the DDGen's).
    pmc_1: optional argument, giving the first pointed matched circle for this element. Must be given if dict is empty.
    pmc_2: optional argument, giving the second pointed matched circle for this element. Must be given if dict is empty.

    Examples:
    """
    def __init__(self, data, pmc_1=[], pmc_2=[]):
        UserDict.__init__(self, data)
        if pmc_1:
            self.pmc_1 = pmc_1
        else:
            self.pmc_1 = self.data.keys()[0].pmc_1
        if pmc_2:
            self.pmc_2 = pmc_2
        else:
            self.pmc_2 = self.data.keys()[0].pmc_2
        self.reduce()

    def coeff(self, x):
        return self.data[x]

    def diff_1_coeffs(self):
        "Returns the element obtained by differentiating the coefficients on the first side of self."
        answer = dict()
        for x in self.data.keys():
            answer[x] = self.data[x].diff_1()
        return DDElt(answer, pmc_1=self.pmc_1,pmc_2=self.pmc_2)

    def diff_2_coeffs(self):
        "Returns the element obtained by differentiating the coefficients on the second side of self."
        answer = dict()
        for x in self.data.keys():
            answer[x] = self.data[x].diff_2()
        return DDElt(answer, pmc_1=self.pmc_1,pmc_2=self.pmc_2)

    def basis_expansion(self):
        "Returns a list of triples [(rho_i,sigma_i,x_i)], with rho_i and sigma_i are Strand_Diagram and x_i a DDGen so that self is sum_i rho_i\otimes sigma_i \otimes x_i."
        self.reduce()
        answer = list()
        for x in self.data.keys():
            for (a,b) in self.data[x].basis_expansion():
                answer.append((a,b,x))
        return answer

    def __repr__(self):
        return repr(self.data)

    def reduce(self):
        "Delete duplicates from self, and also keys with values the empty list."
        for x in self.data.keys():            
        #Eliminate terms where idempotents don't match up.
            self.data[x] = self.data[x]*AlgBlgElt({Strand_Diagram(self.pmc_1,[],x.idem_1):AlgElt([Strand_Diagram(self.pmc_2,[],x.idem_2)])})
            #Now delete pairs of same term.
            self.data[x].reduce_mod_2()
            if self.data[x]==[]:
                del self.data[x]

    def __add__(self, other):
        if self.pmc_1 != other.pmc_1:
            raise Exception("Can't add DDElts with different pmc_1's.")
        if self.pmc_2 != other.pmc_2:
            raise Exception("Can't add DDElts with different pmc_2's.")
        ansdict = dict()
        for x in self.data.keys():
            if not (x in other.data.keys()):
                ansdict[x]=self.data[x]
            else:
                ansdict[x] = self.data[x]+other.data[x]
        for x in other.data.keys():
            if not (x in self.data.keys()):
                ansdict[x]=other.data[x]
        return DDElt(ansdict, self.pmc_1, self.pmc_2)
    
    def __rmul__(self, other):
        "Left multiply an element by an AlgBlgElt."
        ansdict = dict()
        for x in self.data.keys():
            ansdict[x] = other*self.data[x]
        return DDElt(ansdict, self.pmc_1,self.pmc_2)

    def __eq__(self, other):
        if not isinstance(other, DDElt):
            return False
        if Set(self.data.keys()) != Set(other.data.keys()):
            return False
        for x in self.data.keys():
            if self.data[x] != other.data[x]:
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def exchangeLR(self):
        "Returns result of exchanging pmc_1 and pmc_2, i.e., a\otimes b becomes b\otimes a."
        answer = dict()
        for x in self.data.keys():
            answer[x.exchangeLR()]=self.data[x].exchangeLR()
        return DDElt(answer, self.pmc_2, self.pmc_1)

#    def __hash__(self):
#        pass

class TypeDDStr(object):
    """A type DD structure.

    Inputs: 
    --pmc_1: first underlying pointed matched circle.
    --pmc_2: second underlying pointed matched circle.
    --basis: a list of DDGen's
    --structure_consts: a dictionary with keys basis elements and values DDElt's

    Examples:
    #CFDD(Id) for torus boundary:
    p = DDGen('p',small_pmc, small_pmc, [(0,2)], [(0,2)])
    q = DDGen('q',small_pmc, small_pmc, [(1,3)], [(1,3)])
    chords = AlgBlgElt({rho1:AlgElt([rho1]),rho2:AlgElt([rho2]),rho3:AlgElt([rho3]),rho123:AlgElt([rho123])})
    i02chords = AlgBlgElt({iota02:AlgElt([iota02])})*chords
    i13chords = AlgBlgElt({iota13:AlgElt([iota13])})*chords
    TorusId = TypeDDStr(small_pmc, small_pmc, [p,q], {p:DDElt({p:i02chords,q:i02chords}),q:DDElt({p:i13chords,q:i13chords})})
    """

    def __init__(self, pmc_1, pmc_2, basis, structure_consts={}):
        self.pmc_1=pmc_1
        self.pmc_2=pmc_2
        self.basis = basis
        self.structure_consts = structure_consts
        #Should do some sanity checking.

    def __repr__(self):
        parsedstrconsts = str()
        for x in self.structure_consts.keys():
            parsedstrconsts=parsedstrconsts + repr(x)+" -> "+repr(self.structure_consts[x])+"\n"
        return 'Type DD structure with generators '+repr(self.basis)+'\nDelta1 given by\n'+parsedstrconsts
#        return 'Type DD structure with generators '+repr(self.basis)+'\nDelta1 given by '+repr(self.structure_consts)

    def graph(self):
        "Returns a graph representing self."
        verts = self.basis
        edges = dict()
        for v in verts:
            targets = dict()
            for e in self.structure_consts[v].keys():
                targets[repr(e)]=repr(self.structure_consts[v].data[e])
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

    def is_left_bounded(self):
        "Returns True if self is left bounded and False otherwise."
        pass

    def is_right_bounded(self):
        "Returns True if self is right bounded and False otherwise."
        pass

    def bounding_num(self):
        "Returns the highest n so that deltan is nonzero."
        return longest_path_len(self.graph())

    def delta1(self, elt):
        if elt in self.structure_consts.keys():
            return self.structure_consts[elt]
        if DGen(elt) in self.structure_consts.keys():
            return self.structure_consts[DGen(elt)]

    def d_sq_zero(self):
        "Returns True is self satisfies d^2=0 and False otherwise."
        for x in self.basis:
            dx = self.delta1(x)
            dsqx = DDElt({},self.pmc_1,self.pmc_2)
            #differentiate algebra elements on two sides of delta^1(x)
            dsqx = dsqx+dx.diff_1_coeffs()
            dsqx = dsqx+dx.diff_2_coeffs()
            for y in dx.keys():
                dy = self.delta1(y)
                dsqx = dsqx+dx[y]*dy
            dsqx.reduce()
            if dsqx:
                return False
        return True

    def why_d_sq_not_zero(self):
        "Prints an explanation of why d^2 not zero."
        for x in self.basis:
            dx = self.delta1(x)
            dsqx = DDElt({},self.pmc_1,self.pmc_2)
            #differentiate algebra elements on two sides of delta^1(x)
            dsqx = dsqx+dx.diff_1_coeffs()
            dsqx = dsqx+dx.diff_2_coeffs()
            for y in dx.keys():
                dy = self.delta1(y)
                dsqx = dsqx+dx[y]*dy
            dsqx.reduce()
            if dsqx:
                print "Starting with generator: "+repr(x)
                print "d^2 of generator: "+repr(dsqx)
                print "Differential of generator: "+repr(dx)
                print "Differential of first coeffs of df: "+repr(dx.diff_1_coeffs())
                print "Differential of second coeffs of df: "+repr(dx.diff_2_coeffs())
                for y in dx.keys():
                    print "Differential through "+repr(y)+" contributes: "+repr(dx[y]*self.delta1(y))

    def exchangeLR(self):
        "Returns result of exchanging pmc_1 and pmc_2, a DD module over (pmc_2, pmc_1)."
        answer_basis = list()
        for x in self.basis:
            answer_basis.append(x.exchangeLR())
        answer_consts = dict()
        for x in self.basis:
            answer_consts[x.exchangeLR()]=self.structure_consts[x].exchangeLR()
        return TypeDDStr(self.pmc_2, self.pmc_1, answer_basis, answer_consts)

    #Currently assumes we're in middle SpinC structure.
    def mor_to_d(self, other):
        """Returns the chain complex of homomorphisms over A(pmc_1) from self to other, where other is a type D structure.
        
        Examples:
        """
        if self.pmc_1 != other.pmc:
            raise Exception("Can only compute Mor's between structures over the same algebra.")
        #Compile a basis for Mor
        mor_basis = list()
        for m in self.basis:
            for n in other.basis:
                for a in self.pmc_1.alg_basis():
                    if (Set(m.idem_1)==Set(a.left_idem)) and (Set(n.idem)==Set(a.right_idem)):
                        mor_basis.append(DGen((m,a,n),pmc=self.pmc_2,idem=m.idem_2))
        diffs = dict()
#delta^1 on Mor in three parts.
#First: differentiating the algebra element
        for f in mor_basis:
            diffs[f]=DElt({},pmc=self.pmc_2)
            (m,a,n)=f.name
            for b in a.differential():
                diffs[f]=diffs[f]+DElt({DGen((m,b,n),pmc=self.pmc_2,idem=m.idem_2):AlgElt([Strand_Diagram(self.pmc_2,[],m.idem_2)])})
#Second: differentiating the image (element of other)
#        print diffs
        for f in mor_basis:
            (m,a,n)=f.name
            adn = AlgElt([a])*other.delta1(n)
            for (b,p) in adn.basis_expansion():
                diffs[f] = diffs[f]+DElt({DGen((m,b,p),pmc=self.pmc_2,idem=m.idem_2):AlgElt([Strand_Diagram(self.pmc_2,[],m.idem_2)])})
#Third: differentiating the source (element of self). This is the only part which outputs non-idempotent algebra elements.
#        print diffs
        for f in mor_basis:
            (m,a,n)=f.name
            for l in self.basis:
                dl=self.delta1(l)
                if m in dl.keys():
                    for b in dl[m].keys():
                        if b*a:
                            diffs[f] = diffs[f]+DElt({DGen((l,b*a,n),pmc=self.pmc_2,idem=l.idem_2):dl[m][b]})
#        print diffs
        #Now view PMC as over opposite algebra:
        #First, reverse the idempotents.
        rev_basis = list()
        for f in mor_basis:
            (m,a,n)=f.name
            rev_idem = list()
            for (i,j) in m.idem_2:
                rev_idem.append((4*self.pmc_2.genus-1-j,4*self.pmc_2.genus-1-i))
            rev_basis.append(DGen((m,a,n),pmc=self.pmc_2,idem=rev_idem))
        rev_diffs = dict()
        #Now, reverse the algebra element outputs.
        for f in diffs.keys():
            rev_diffs[f]=diffs[f].opposite()
        return TypeDStr(self.pmc_2.opposite(), rev_basis, rev_diffs)


def change_framing(module, bimodules):
    "Twist the framing on module by taking mor_to_d from the bimodules in list bimodules."
    answer = TypeDStr(module.pmc, module.basis, module.structure_consts)
    for x in bimodules:
        answer = x.mor_to_d(answer)
        answer.simplify()
        answer.shorten_names()
    return answer

def dd_identity(pmc):
    "Returns the type DD identity bimodule associated to pmc."
    generators = list()
    oppositepmc = pmc.opposite()
    #Find a list of generators
    idemps = pmc.idempotents()
    for x in idemps:
        generators.append(DDGen(idemps.index(x),pmc, oppositepmc,x,pmc.complementary_idem(x)))
    #Compile a basis for the diagonal subalgebra
    diag_basis = pmc.zero()**oppositepmc.zero()
    for i in range(4*pmc.genus):
        for j in range(i+1,4*pmc.genus):
            for idem in pmc.idempotents():
                if l_idem_compat(pmc, [(i,j)],idem) and l_idem_compat(oppositepmc, [(4*pmc.genus-j-1,4*pmc.genus-i-1)],pmc.complementary_idem(idem)):
                    diag_basis = diag_basis + AlgElt(Strand_Diagram(pmc,[(i,j)],idem))**AlgElt(Strand_Diagram(oppositepmc,[(4*pmc.genus-j-1,4*pmc.genus-i-1)],pmc.complementary_idem(idem)))
    #Now, differentials...
    differentials = dict()
#    print "Diag basis is"+repr(diag_basis)
    for x in generators:
        dx = DDElt({}, pmc, oppositepmc)
        for y in generators:
            dx = dx + DDElt({y:(AlgElt(Strand_Diagram(pmc,[],x.idem_1))**AlgElt(Strand_Diagram(oppositepmc, [], x.idem_2))*diag_basis*AlgElt(Strand_Diagram(pmc,[],y.idem_1))**AlgElt(Strand_Diagram(oppositepmc, [], y.idem_2)))}, pmc, oppositepmc)
        differentials[x]=dx
    return TypeDDStr(pmc, oppositepmc, generators, differentials)

class Underslide(object):
    """Class for computing DD modules of underslides.

    Input:
    pmc: The pointed matched circle you're arc-sliding.
    i: the moving endpoint of the arc-slide
    j: the endpoint i is sliding over

    Methods:
    dd_mod: generates and returns the type DD module for this underslide. The DD module is also stored in self.ddstr. Calling dd_mod() only generates it the first time; afterwards, the existing dd structure is returned. (So, subsequent calls are quick.)
    show(): display this arc-slide graphically.
    Lots of others used in the construction (and debugging) of dd_mod.

    Other useful data:
    self.pmc_2: the pointed matched circle resulting from the arc-slide.
    self.pmc_2rev: the orientation-reverse of self.pmc_2.

    Caveats:
    Currently restricts to the central SpinC structure, though this would be relatively easy to change.
    """
    def __init__(self, pmc, i, j=i+1):
        if not pmc.is_underslide(i,j):
            raise Exception("Sliding i over j is not an underslide.")
        self.i = i
        self.j = j
        self.u1chords = None
        self.u2chords = None
        self.u3chords = None
        self.u4chords = None
        self.u5chords = None
        self.u6chords = None
        self.allchords = None
        self.gen_idems = None
        self.generators = None
        self.ddstr = None
    #The PMC resulting from slide:
        self.pmc_1 = pmc
        self.pmc_2 = pmc.arcslide(i,j)
        self.pmc_2rev = self.pmc_2.opposite()
        self.genus = self.pmc_1.genus
    #Notation from HFa:
        self.c1 = j
        self.c2 = pmc.matched_point(self.c1)
        self.b1 = i
        self.b2 = pmc.matched_point(self.b1)
        self.c1p = self.c1
        self.c2p = self.pmc_2.matched_point(self.c1p)
        if self.b1 > self.c1:
            self.b1p = self.c2p-1
        if self.b1 < self.c1:
            self.b1p = self.c2p+1
        self.b2p = self.pmc_2.matched_point(self.b1p)
        
    def is_up_underslide(self):
        if self.j==self.i+1:
            return True
        return False

    def is_down_underslide(self):
        if self.j==self.i-1:
            return True
        return False

    def r(self, x):
        "Take a point a in self.pmc_1 to the corresponding point in self.pmc_2, assuming a is not b_1."
        if x==self.b1:
            raise Exception("r not defined on b1.")
        if (x < self.b1) and (x<self.c2):
            return x
        if (x > self.b1) and (x>self.c2):
            return x
        if self.b1<x<self.c2:
            return x-1
        if self.c2<x<self.b1:
            return x+1
        if x==self.c2:
            return self.c2p

    def rinv(self, x):
        "The inverse of the map r."
        if x==self.b1p:
            raise Exception("r inverse not defined on b1'.")
        if (x < self.b1p) and (x<self.c1p):
            return x
        if (x > self.b1p) and (x>self.c1p):
            return x
        if self.b1p<x<self.c1p:
            return x-1
        if self.c1p<x<self.b1p:
            return x+1
        if x==self.c1p:
            return self.c1

    def generate_U1(self):
        "Generates a list of chords of type U1. Stored in self.u1chords, in the form (xi, xi'), where xi is a list of chords (containing exactly one chord), and similarly xi'."
        answer = list()
        for x in range(4*self.pmc_1.genus):
            for y in range(x+1, 4*self.pmc_1.genus):
                if (x != self.b1) and (y!= self.b1) and (not (x,y) in self.pmc_1.matching) and (not (y,x) in self.pmc_1.matching):
                    answer.append(([(x,y)],[(self.r(x),self.r(y))]))
        self.u1chords = answer

    def generate_U2(self):
        "Generates a list of chords of type U2. Stored in self.u2chords"
        answer = list()
        if self.b1 < self.c1:
            answer.append(([(self.b1,self.c1)],[]))
        if self.b1 > self.c1:
            answer.append(([(self.c1,self.b1)],[]))
        if self.b1p < self.c2p:
            answer.append(([],[(self.b1p,self.c2p)]))
        if self.b1p > self.c2p:
            answer.append(([],[(self.c2p,self.b1p)]))
        self.u2chords = answer

    def generate_U3(self):
        "Generates a list of chords of type U3. Stored in self.u3chords"
        answer = list()
        #Extra sigma
        if self.b1 < self.c1:
            for x in range(self.c1+1,4*self.genus):
                answer.append(([(self.b1,x)],[(self.c1p,self.r(x))]))
        if self.b1 > self.c1:
            for x in range(self.c1):
                answer.append(([(x,self.b1)],[(self.r(x),self.c1p)]))
        #Extra sigma'
        if self.b1p < self.c2p:
            for x in range(self.c2p+1,4*self.genus):
                answer.append(([(self.c2,self.rinv(x))],[(self.b1p,x)]))
        if self.b1p > self.c2p:
            for x in range(self.c2p):
                answer.append(([(self.rinv(x),self.c2)],[(x,self.b1p)]))
        self.u3chords = answer

    def generate_U4(self):
        "Generates a list of chords of type U4. Stored in self.u4chords, in the form (xi, xi'), where xi is a list of chords, and similarly xi'."
        answer = list()
        #Two connected chords
        #Missing sigma:
        if self.b1 < self.c1:
            for x in range(self.b1):
                answer.append(([(x,self.b1)],[(self.r(x),self.c1p)]))
        if self.b1 > self.c1:
            for x in range(self.b1+1, 4*self.genus):
                answer.append(([(self.b1,x)],[(self.c1p,self.r(x))]))
        #Missing sigma':
        if self.b1p < self.c2p:
            for x in range(self.b1p):
                answer.append(([(self.rinv(x),self.c2)],[(x,self.b1p)]))
        if self.b1p > self.c2p:
            for x in range(self.b1p+1, 4*self.genus):
                answer.append(([(self.c2, self.rinv(x))],[(self.b1p,x)]))
        #Three connected chords total, two on left.
        if self.b1 < self.c1:
            for x in range(self.b1):
                for y in range(self.c1+1,4*self.genus):
                    answer.append(([(x,self.b1),(self.c1,y)],[(self.r(x),self.r(y))]))
        if self.b1 > self.c1:
            for x in range(self.c1):
                for y in range(self.b1+1, 4*self.genus):
                    answer.append(([(x,self.c1),(self.b1,y)],[(self.r(x),self.r(y))]))
        #Three connected chords total, two on right.
        if self.b1p < self.c2p:
            for x in range(self.b1p):
                for y in range(self.c2p+1, 4*self.genus):
                    answer.append(([(self.rinv(x),self.rinv(y))],[(x,self.b1p),(self.c2p,y)]))
        if self.b1p > self.c2p:
            for x in range(self.c2p):
                for y in range(self.b1p+1, 4*self.genus):
                    answer.append(([(self.rinv(x),self.rinv(y))],[(x,self.c2p),(self.b1p,y)]))
        self.u4chords = answer

    def generate_U5(self):
        "Generates a list of chords of type U5. Stored in self.u5chords"
        answer = list()
        smaller = min(self.c1, self.c2)
        bigger = max(self.c1,self.c2)
        rbigger = self.r(bigger)
        rsmaller = self.r(smaller)
        for x in range(smaller):
            for y in range(bigger+1, 4*self.genus): #Used to be range(bigger,4*self.genus). Bohua caught the bug.
                answer.append(([(x,smaller),(bigger,y)],[(self.r(x),rsmaller),(rbigger,self.r(y))]))
        self.u5chords = answer

    def generate_U6(self):
        "Generates a list of chords of type U6. Stored in self.u6chords"
        answer = list()
        #Extra sigma:
        if self.b1 < self.c1:
            for x in range(self.c2+1, 4*self.genus):
                if (x != self.b1) and (x != self.c1):
                    answer.append(([(self.c2,x), (self.b1,self.c1)],[(self.b1p, self.r(x))]))
        if self.b1 > self.c1:
            for x in range(self.c2):
                if (x != self.b1) and (x != self.c1):
                    answer.append(([(x, self.c2), (self.c1,self.b1)],[(self.r(x),self.b1p)]))
        #Extra sigma':
        if self.b1p < self.c2p:
            for x in range(self.c1p+1, 4*self.genus):
                if (x != self.b1p) and (x != self.c2p):
                    answer.append(([(self.b1, self.rinv(x))],[(self.c1p,x),(self.b1p,self.c2p)]))
        if self.b1p > self.c2p:
            for x in range(self.c1p):
                if (x != self.b1p) and (x != self.c2p):
                    answer.append(([(self.rinv(x),self.b1)],[(x,self.c1p),(self.c2p,self.b1p)]))
        self.u6chords = answer

    def generate_chords(self):
        "Generates a list of near-chords for self."
        regenallchords = False
        if self.allchords == None:
            self.regenallchords = True
        if self.u1chords == None:
            self.generate_U1()
            regenallchords = True
        if self.u2chords == None:
            self.generate_U2()
            regenallchords = True
        if self.u3chords == None:
            self.generate_U3()
            regenallchords = True
        if self.u4chords == None:
            self.generate_U4()
            regenallchords = True
        if self.u5chords == None:
            self.generate_U5()
            regenallchords = True
        if self.u6chords == None:
            self.generate_U6()
            regenallchords = True
        if regenallchords:
            self.allchords = self.u1chords + self.u2chords + self.u3chords + self.u4chords + self.u5chords + self.u6chords

    def complementary_idem(self, idem):
        "Given an idempotent of self.pmc_1, return the complementary idempotent of self.pmc_2."
        answer = list(self.pmc_2.matching)
        for (x,y) in self.pmc_2.matching:
            if (x != self.b1p) and (y != self.b1p ):
                if ( (self.rinv(x), self.rinv(y)) in idem ) or ( (self.rinv(y), self.rinv(x))  in idem):
                    answer.remove((x,y))
        if ( (self.b1, self.b2) in idem ) or ( (self.b2,self.b1) in idem ):
            if self.b1p < self.b2p:
                answer.remove((self.b1p, self.b2p))
            else:
                answer.remove((self.b2p, self.b1p))
        answer.sort()
        return answer

    def generate_generators(self):
        """Generates a list of generators of the DD module for this arcslide. Stored in self.generators
        Currently restricts to the middle SpinC structure.
        """
        answer = list()
        #Generate gen_idems, if necessary
        if self.gen_idems == None:
            self.generate_gen_idems()
        for (x,y) in self.gen_idems:
            yrev = list()
            for (m,n) in y:
                yrev.append((4*self.genus-n-1, 4*self.genus-m-1))
            answer.append(DDGen(repr(x)+'|'+repr(yrev),self.pmc_1, self.pmc_2rev,x,yrev))
        self.generators = answer

    def generate_gen_idems(self):
        """Generates a list of pairs of idempotents for generators for the DD module. Stored in self.gen_idems
        Currently restricts to the middle SpinC structure.
        """
        #Generators of type X (complementary)
        xidems = list()
        pmc1idems = self.pmc_1.idempotents()
        for idem in pmc1idems:
            xidems.append((idem, self.complementary_idem(idem)))
        #Generators of type Y (sub-complementary)
        yidems = list()
        for idem in pmc1idems:
            if ( (self.c1, self.c2) in idem ) or ( (self.c2, self.c1) in idem):
                if ( not ( (self.b1, self.b2) in idem)) and ( not ( (self.b2, self.b1) in idem) ):
                    scompidem = self.complementary_idem(idem)
                    if self.b1p<self.b2p:
                        scompidem.remove((self.b1p,self.b2p))
                    if self.b2p<self.b1p:
                        scompidem.remove((self.b2p,self.b1p))
                    if self.c1p < self.c2p:
                        scompidem.append((self.c1p,self.c2p))
                    if self.c2p < self.c1p:
                        scompidem.append((self.c2p,self.c1p))
                    scompidem.sort()
                    yidems.append((idem, scompidem))                        
        answer = xidems + yidems
        self.gen_idems = answer

    def dd_mod(self):
        "Return the type DD module for this arcslide."
        if self.ddstr != None:
            return self.ddstr
        if self.generators == None:
            self.generate_generators()
        self.generate_chords()
        diffs = dict()
        for x in self.generators:
            dx = DDElt({}, self.pmc_1, self.pmc_2rev)
            for y in self.generators:
                for (a,b) in self.allchords:
                    #If a is compatible with x's left idempotents and y's left idempotents, and b is compatible with x's right idempotents and y's right idempotents, add corresponding elt to dx.
                    if l_idem_compat(self.pmc_1, a, x.idem_1) and l_idem_compat(self.pmc_2rev, self.opposite_strands(b), x.idem_2):
                        aalgelt = AlgElt(Strand_Diagram(self.pmc_1, a, x.idem_1))*AlgElt(Strand_Diagram(self.pmc_1, [], y.idem_1))
                        balgelt = AlgElt(Strand_Diagram(self.pmc_2rev, self.opposite_strands(b), x.idem_2))*AlgElt(Strand_Diagram(self.pmc_2rev, [], y.idem_2))
                        abalgelt = aalgelt**balgelt
                        if abalgelt:
                            dx = dx + DDElt({y:abalgelt}, self.pmc_1, self.pmc_2rev)
            diffs[x] = dx
        self.ddstr = TypeDDStr(self.pmc_1, self.pmc_2rev, self.generators, diffs)
        return self.ddstr

    def opposite_strands(self, strands):
        answer = list()
        for (x,y) in strands:
            answer.append((4*self.genus-y-1, 4*self.genus-x-1))
        return answer

    def show(self):
        "Display graphically the pair (lchord, rchord). Intended for testing generate_U1() ... generate_U6()."
        #First show the arcslide
        pict = line([(0,-1), (0,4*self.genus)], rgbcolor = (0,0,0))
        pict = pict + line([(2,-1), (2,4*self.genus)], rgbcolor = (0,0,0))
        for i in range(4*self.genus):
            if i != self.b1:
                pict = pict + line([(0,i),(2,self.r(i))], rgbcolor = (1,0,0))
        if self.b1 < self.c1:
            pict = pict + line([(0,self.b1), (1,self.b1+1)], rgbcolor = (1,.25,.25)) + line([(1,self.b1p-1), (2,self.b1p)], rgbcolor = (1,.25,.25))
        if self.b1 > self.c1:
            pict = pict + line([(0,self.b1), (1,self.b1-1)], rgbcolor = (1,.25,.25)) + line([(1,self.b1p+1), (2,self.b1p)], rgbcolor = (1,.25,.25))
        #Show the matchings on the two sides:
        #left:
        for l in range(len(self.pmc_1.matching)):
            pict = pict + line([(-0.2, self.pmc_1.matching[l][0]), (-0.3-l/(8*self.genus), self.pmc_1.matching[l][0])], rgbcolor = (.2,.2,.2))
            pict = pict + line([(-0.2, self.pmc_1.matching[l][1]), (-0.3-l/(8*self.genus), self.pmc_1.matching[l][1])], rgbcolor = (.2,.2,.2))
            pict = pict + line([(-0.3-l/(8*self.genus), self.pmc_1.matching[l][0]), (-0.3-l/(8*self.genus), self.pmc_1.matching[l][1])], rgbcolor = (.2,.2,.2))
        #right:
        for l in range(len(self.pmc_2.matching)):
            pict = pict + line([(2.2, self.pmc_2.matching[l][0]), (2.3+l/(8*self.genus), self.pmc_2.matching[l][0])], rgbcolor = (.2,.2,.2))
            pict = pict + line([(2.2, self.pmc_2.matching[l][1]), (2.3+l/(8*self.genus), self.pmc_2.matching[l][1])], rgbcolor = (.2,.2,.2))
            pict = pict + line([(2.3+l/(8*self.genus), self.pmc_2.matching[l][0]), (2.3+l/(8*self.genus), self.pmc_2.matching[l][1])], rgbcolor = (.2,.2,.2))
        pict.show(axes=False) #, xmin=-2,xmax=2, ymin=-.5, ymax=4*self.genus-.5)


    def display_chord(self, chords):
        "Display graphically the pair (lchord, rchord). Intended for testing generate_U1() ... generate_U6()."
        lchord = chords[0]
        rchord = chords[1]
        #First show the arcslide
        pict = line([(0,-1), (0,4*self.genus)], rgbcolor = (0,0,0))
        pict = pict + line([(2,-1), (2,4*self.genus)], rgbcolor = (0,0,0))
        for i in range(4*self.genus):
            if i != self.b1:
                pict = pict + line([(0,i),(2,self.r(i))], rgbcolor = (1,0,0))
        if self.b1 < self.c1:
            pict = pict + line([(0,self.b1), (1,self.b1+1)], rgbcolor = (1,.25,.25)) + line([(1,self.b1p-1), (2,self.b1p)], rgbcolor = (1,.25,.25))
        if self.b1 > self.c1:
            pict = pict + line([(0,self.b1), (1,self.b1-1)], rgbcolor = (1,.25,.25)) + line([(1,self.b1p+1), (2,self.b1p)], rgbcolor = (1,.25,.25))
        #light lines going farther out, for clarity.
        for x in range(4*self.genus):
            pict = pict + line([(-1,x),(0,x)], rgbcolor = (.8,.8,.8))
            pict = pict + line([(2,x),(3,x)], rgbcolor = (.8,.8,.8))
        #Show the matchings on the two sides:
        #left:
        for l in range(len(self.pmc_1.matching)):
            pict = pict + line([(-0.2, self.pmc_1.matching[l][0]), (-0.3-l/(8*self.genus), self.pmc_1.matching[l][0])], rgbcolor = (.2,.2,.2))
            pict = pict + line([(-0.2, self.pmc_1.matching[l][1]), (-0.3-l/(8*self.genus), self.pmc_1.matching[l][1])], rgbcolor = (.2,.2,.2))
            pict = pict + line([(-0.3-l/(8*self.genus), self.pmc_1.matching[l][0]), (-0.3-l/(8*self.genus), self.pmc_1.matching[l][1])], rgbcolor = (.2,.2,.2))
        #right:
        for l in range(len(self.pmc_2.matching)):
            pict = pict + line([(2.2, self.pmc_2.matching[l][0]), (2.3+l/(8*self.genus), self.pmc_2.matching[l][0])], rgbcolor = (.2,.2,.2))
            pict = pict + line([(2.2, self.pmc_2.matching[l][1]), (2.3+l/(8*self.genus), self.pmc_2.matching[l][1])], rgbcolor = (.2,.2,.2))
            pict = pict + line([(2.3+l/(8*self.genus), self.pmc_2.matching[l][0]), (2.3+l/(8*self.genus), self.pmc_2.matching[l][1])], rgbcolor = (.2,.2,.2))
        #Now draw the chord on the two sides.
        for (x,y) in lchord:
            pict = pict + line([(-1,x),(-0.1, y)], rgbcolor = (0,0,1))
        for (x,y) in rchord:
            pict = pict + line([(2.1,x),(3, y)], rgbcolor = (0,0,1))
        pict.show(axes=False) #, xmin=-2,xmax=2, ymin=-.5, ymax=4*self.genus-.5)

class Overslide(object):
    """Class for computing DD modules of overslides.

    Input:
    pmc: The pointed matched circle you're arc-sliding.
    i: the moving endpoint of the arc-slide
    j: the endpoint i is sliding over

    Methods:
    dd_mod: generates and returns the type DD module for this underslide. The DD module is also stored in self.ddstr. Calling dd_mod() only generates it the first time; afterwards, the existing dd structure is returned. (So, subsequent calls are quick.)
    show(): display this arc-slide graphically.
    Lots of others used in the construction (and debugging) of dd_mod.

    Other useful data:
    self.pmc_2: the pointed matched circle resulting from the arc-slide.
    self.pmc_2rev: the orientation-reverse of self.pmc_2.

    Caveats:
    Currently restricts to the central SpinC structure, though this would be relatively easy to change.
    """
    def __init__(self, pmc, i, j=i+1):
        if not pmc.is_overslide(i,j):
            raise Exception("Sliding i over j is not an overslide.")
        self.i = i
        self.j = j
        self.o1DetChords = None
        self.o2DetChords = None
        self.o3DetChords = None
        self.o4DetChords = None
        self.o5DetChords = None
        self.o6DetChords = None
        self.o3IndetChords = None
        self.o4IndetChords = None
        self.o7IndetChords = None
        self.o8IndetChords = None
        self.DetChords = None
        self.IndetChords = None
        self.gen_idems = None
        self.generators = None
        self.ddstr = None
    #The PMC resulting from slide:
        self.pmc_1 = pmc
        self.pmc_2 = pmc.arcslide(i,j)
        self.pmc_2rev = self.pmc_2.opposite()
        self.genus = self.pmc_1.genus
    #Notation from HFa:
        self.c1 = j
        self.c2 = pmc.matched_point(self.c1)
        self.b1 = i
        self.b2 = pmc.matched_point(self.b1)
        if self.b1>self.c1:
            self.c1p = self.c1+1
        if self.b1<self.c1:
            self.c1p=self.c1-1
        self.c2p = self.pmc_2.matched_point(self.c1p)
        if self.b1 > self.c1:
            self.b1p = self.c2p-1
        if self.b1 < self.c1:
            self.b1p = self.c2p+1
        self.b2p = self.pmc_2.matched_point(self.b1p)
        
    def is_up_overslide(self):
        if self.j==self.i+1:
            return True
        return False

    def is_down_overslide(self):
        if self.j==self.i-1:
            return True
        return False

    def r(self, x):
        "Take a point a in self.pmc_1 to the corresponding point in self.pmc_2, assuming a is not b_1."
        if x==self.b1:
            raise Exception("r not defined on b1.")
        if (x < self.b1) and (x<self.c2):
            return x
        if (x > self.b1) and (x>self.c2):
            return x
        if self.b1<x<self.c2:
            return x-1
        if self.c2<x<self.b1:
            return x+1
        if x==self.c2:
            return self.c2p

    def rinv(self, x):
        "The inverse of the map r."
        if x==self.b1p:
            raise Exception("r inverse not defined on b1'.")
        if (x < self.b1p) and (x<self.c1p):
            return x
        if (x > self.b1p) and (x>self.c1p):
            return x
        if self.b1p<x<self.c1p:
            return x-1
        if self.c1p<x<self.b1p:
            return x+1
        if x==self.c1p:
            return self.c1

    def generate_O1(self):
        "Generates a list of chords of type O1. Stored in self.o1DetChords (since they're all determinate), in the form (xi, xi'), where xi is a list of chords (containing exactly one chord), and similarly xi'."
        answer = list()
        for x in range(4*self.pmc_1.genus):
            for y in range(x+1, 4*self.pmc_1.genus):
                if (x != self.b1) and (y!= self.b1) and (not (x,y) in self.pmc_1.matching) and (not (y,x) in self.pmc_1.matching):
                    answer.append(([(x,y)],[(self.r(x),self.r(y))]))
        self.o1DetChords = answer

    def generate_O2(self):
        "Generates a list of chords of type O2. Stored in self.O2DetChords"
        answer = list()
        if self.b1 < self.c1:
            answer.append(([(self.b1,self.c1)],[]))
        if self.b1 > self.c1:
            answer.append(([(self.c1,self.b1)],[]))
        if self.b1p < self.c2p:
            answer.append(([],[(self.b1p,self.c2p)]))
        if self.b1p > self.c2p:
            answer.append(([],[(self.c2p,self.b1p)]))
        self.o2DetChords = answer

    def generate_O3(self):
        "Generates a list of chords of type O3. Stored in self.o3DetChords and self.o3IndetChords, depending on whether the chord is determinate or not."
        #The determinate ones:
        answer = list()
        #Extra sigma
        if self.b1 < self.c1:
            for x in range(self.c1+1,4*self.genus):
                if x != self.c2:
                    answer.append(([(self.b1,x)],[(self.c1p,self.r(x))]))
        if self.b1 > self.c1:
            for x in range(self.c1):
                if x!= self.c2:
                    answer.append(([(x,self.b1)],[(self.r(x),self.c1p)]))
        #Extra sigma'
        if self.b1p < self.c2p:
            for x in range(self.c2p+1,4*self.genus):
                if x!= self.c1p:
                    answer.append(([(self.c2,self.rinv(x))],[(self.b1p,x)]))
        if self.b1p > self.c2p:
            for x in range(self.c2p):
                if x!= self.c1p:
                    answer.append(([(self.rinv(x),self.c2)],[(x,self.b1p)]))
        self.o3DetChords = answer
        #Now, the indeterminate ones
        indet = list()
        if self.b1<self.c1:
            indet.append(([(self.b1, self.c2)],[(self.c1p, self.c2p)]))
            indet.append(([(self.c1, self.c2)],[(self.c1p, self.b1p)]))
        if self.b1>self.c1:
            indet.append(([(self.c2, self.b1)],[(self.c2p, self.c1p)]))
            indet.append(([(self.c2, self.c1)],[(self.b1p, self.c1p)]))
        self.o3IndetChords = indet

    def generate_O4(self):
        "Generates a list of chords of type U4. Stored in self.o4DetChords and self.o4IndetChords, in the form (xi, xi'), where xi is a list of chords, and similarly xi'."
        detanswer = list()
        indetanswer = list()
        #Two connected chords
        #Missing sigma:
        if self.b1 < self.c1:
            for x in range(self.b1):
                detanswer.append(([(x,self.b1)],[(self.r(x),self.c1p)]))
        if self.b1 > self.c1:
            for x in range(self.b1+1, 4*self.genus):
                detanswer.append(([(self.b1,x)],[(self.c1p,self.r(x))]))
        #Missing sigma':
        if self.b1p < self.c2p:
            for x in range(self.b1p):
                detanswer.append(([(self.rinv(x),self.c2)],[(x,self.b1p)]))
        if self.b1p > self.c2p:
            for x in range(self.b1p+1, 4*self.genus):
                detanswer.append(([(self.c2, self.rinv(x))],[(self.b1p,x)]))
        #Three connected chords total, two on left.
        if self.b1 < self.c1:
            for x in range(self.b1):
                for y in range(self.c1+1,4*self.genus):
                    #Is this indeterminate?
                    if y==self.c2:
                        indetanswer.append(([(x,self.b1),(self.c1,y)],[(self.r(x),self.r(y))]))
                    else:
                        detanswer.append(([(x,self.b1),(self.c1,y)],[(self.r(x),self.r(y))]))                        
        if self.b1 > self.c1:
            for x in range(self.c1):
                for y in range(self.b1+1, 4*self.genus):
                    #Is this indeterminate?
                    if x==self.c2:
                        indetanswer.append(([(x,self.c1),(self.b1,y)],[(self.r(x),self.r(y))]))                        
                    else:
                        detanswer.append(([(x,self.c1),(self.b1,y)],[(self.r(x),self.r(y))]))
        #Three connected chords total, two on right.
        if self.b1p < self.c2p:
            for x in range(self.b1p):
                for y in range(self.c2p+1, 4*self.genus):
                    #Is this indeterminate?
                    if y==self.c1p:
                        indetanswer.append(([(self.rinv(x),self.rinv(y))],[(x,self.b1p),(self.c2p,y)]))
                    else:
                        detanswer.append(([(self.rinv(x),self.rinv(y))],[(x,self.b1p),(self.c2p,y)]))
        if self.b1p > self.c2p:
            for x in range(self.c2p):
                for y in range(self.b1p+1, 4*self.genus):
                    #Is this indeterminate?
                    if x==self.c1p:
                        indetanswer.append(([(self.rinv(x),self.rinv(y))],[(x,self.c2p),(self.b1p,y)]))
                    else:
                        detanswer.append(([(self.rinv(x),self.rinv(y))],[(x,self.c2p),(self.b1p,y)]))
        self.o4DetChords = detanswer
        self.o4IndetChords = indetanswer

    def generate_O5(self):
        "Generates a list of chords of type O5. Stored in self.o5DetChords"
        answer = list()
        smaller = min(self.c1, self.c2)
        bigger = max(self.c1,self.c2)
        rbigger = self.r(bigger)
        rsmaller = self.r(smaller)
        #Chords can be disjoint:
        for x in range(smaller+1,bigger):
            for y in range(x+1,bigger):
                answer.append(([(smaller,x),(y,bigger)],[(rsmaller,self.r(x)),(self.r(y),rbigger)]))
        #Or one can be contained in the other:
        for x in range(bigger,4*self.genus):
            for y in range(smaller+1,bigger):
                if (x != self.b1) and (x != self.b2):
                    answer.append(([(smaller,x),(y,bigger)],[(rsmaller,self.r(x)),(self.r(y),rbigger)]))
        for x in range(smaller+1,bigger):
            for y in range(smaller):
                if (y != self.b1) and (y != self.b2):
                    answer.append(([(smaller,x),(y,bigger)],[(rsmaller,self.r(x)),(self.r(y),rbigger)]))
        self.o5DetChords = answer

    def generate_O6(self):
        "Generates a list of chords of type O6. Stored in self.o6DetChords"
        answer = list()
        #Extra sigma:
        if self.b1 > self.c1:
            for x in range(self.c2):
                answer.append(([(x, self.c2),(self.c1,self.b1)],[(self.r(x),self.b1p)]))
        if self.b1 < self.c1:
            for x in range(self.c2+1,4*self.genus):
                answer.append(([(self.b1,self.c1),(self.c2,x)],[(self.b1p,self.r(x))]))

        #Extra sigma'
        if self.b1 > self.c1:
            for x in range(self.b1+1,4*self.genus):
                answer.append(([(self.b1, x)],[(self.c1p,self.r(x)),(self.b1p,self.c2p)]))
        if self.b1 < self.c1:
            for x in range(self.b1):
                answer.append(([(x,self.b1)],[(self.r(x),self.c1p),(self.c2p,self.b1p)]))
        self.o6DetChords = answer

    def generate_O7(self):
        "Generates a list of chords of type O7. Stored in self.o7IndetChords"   
        answer = list()
        #Break on left.
        if self.c1<self.c2:
            for x in range(self.c1+1,self.c2):
                answer.append(([(self.c1,x),(x,self.c2)],[(self.c1p,self.c2p)]))
        if self.c2<self.c1:
            for x in range(self.c2+1,self.c1):
                answer.append(([(self.c2,x),(x,self.c1)],[(self.c2p,self.c1p)]))
        #Break on right.
        if self.c1p<self.c2p:
            for x in range(self.c1p+1,self.c2p):
                answer.append(([(self.c1,self.c2)],[(self.c1p,x),(x,self.c2p)]))
        if self.c2p<self.c1p:
            for x in range(self.c2p+1,self.c1p):
                answer.append(([(self.c2,self.c1)],[(self.c2p,x),(x,self.c1p)]))
        self.o7IndetChords = answer

    def generate_O8(self):
        "Generates a list of chords of type O7. Stored in self.o7IndetChords"   
        answer = list()
        #xi contained in [c1,c2]
        smaller = min(self.c1, self.c2)
        bigger = max(self.c1, self.c2)
        rsmaller = min(self.c1p, self.c2p)
        rbigger = max(self.c1p, self.c2p)
        for x in range(smaller+1,bigger):
            for y in range(x+1,bigger):
                if (not (x,y) in self.pmc_1.matching) and (not (y,x) in self.pmc_1.matching):
                    answer.append(([(x,y),(smaller,bigger)],[(self.r(x),self.r(y)),(rsmaller,rbigger)]))
        #xi disjoint from [c1,c2]
        #xi below than c1:
        for x in range(smaller):
            for y in range(x+1,smaller):
                if (not (x,y) in self.pmc_1.matching) and (not (y,x) in self.pmc_1.matching) and (x!=self.b1) and (y!=self.b1):
                    answer.append(([(x,y),(smaller,bigger)],[(self.r(x),self.r(y)),(rsmaller,rbigger)]))
        #xi above c2:
        for x in range(bigger+1, 4*self.genus):
            for y in range(x+1, 4*self.genus):
                if (not (x,y) in self.pmc_1.matching) and (not (y,x) in self.pmc_1.matching) and (x != self.b1) and (y !=self.b1):
                    answer.append(([(x,y),(smaller,bigger)],[(self.r(x),self.r(y)),(rsmaller,rbigger)]))
        #NOTE: for idempotent reasons, chords of this type with an endpoint on b1 or b2 wont count. But they have not been explicitly disallowed above.
        self.o8IndetChords = answer

    def generate_chords(self):
        "Generates a list of near-chords for self."
        regenallchords = False
        if (self.DetChords == None) or (self.Indetchords == None) :
            regenallchords =  True
        if self.o1DetChords == None:
            self.generate_O1()
            regenallchords = True
        if self.o2DetChords == None:
            self.generate_O2()
            regenallchords = True
        if self.o3DetChords == None:
            self.generate_O3()
            regenallchords = True
        if self.o4DetChords == None:
            self.generate_O4()
            regenallchords = True
        if self.o5DetChords == None:
            self.generate_O5()
            regenallchords = True
        if self.o6DetChords == None:
            self.generate_O6()
            regenallchords = True
        if self.o7IndetChords == None:
            self.generate_O7()
            regenallchords = True
        if self.o8IndetChords == None:
            self.generate_O8()
            regenallchords = True
        if regenallchords:
            self.DetChords = self.o1DetChords+self.o2DetChords+self.o3DetChords+self.o4DetChords+self.o5DetChords+self.o6DetChords
            self.IndetChords = self.o3IndetChords+self.o4IndetChords+self.o7IndetChords+self.o8IndetChords

    def basic_choice(self):
        "Generate the indeterminate chords corresponding to a particular basic chioce."
        answer = list()
        #The indeterminate O3 chords: those using sigma', not those using sigma.
        if self.b1<self.c1:
            answer.append(([(self.c1, self.c2)],[(self.c1p, self.b1p)]))
        if self.b1>self.c1:
            answer.append(([(self.c2, self.c1)],[(self.b1p, self.c1p)]))
        #The indeterminate O4 chords
        #Three connected chords total, two on left.
        if self.b1 < self.c1:
            for x in range(self.b1):
                answer.append(([(x,self.b1),(self.c1,self.c2)],[(self.r(x),self.c2p)]))
        if self.b1 > self.c1:
            for y in range(self.b1+1, 4*self.genus):
                answer.append(([(self.c2,self.c1),(self.b1,y)],[(self.c2p,self.r(y))]))                        
        #The O7 chords: those with a break on the left.
        if self.c1<self.c2:
            for x in range(self.c1+1,self.c2):
                answer.append(([(self.c1,x),(x,self.c2)],[(self.c1p,self.c2p)]))
        if self.c2<self.c1:
            for x in range(self.c2+1,self.c1):
                answer.append(([(self.c2,x),(x,self.c1)],[(self.c2p,self.c1p)]))
        #No O8 chords.
        self.choice_chords = answer

    def complementary_idem(self, idem):
        "Given an idempotent of self.pmc_1, return the complementary idempotent of self.pmc_2."
        answer = list(self.pmc_2.matching)
        for (x,y) in self.pmc_2.matching:
            if (x != self.b1p) and (y != self.b1p ):
                if ( (self.rinv(x), self.rinv(y)) in idem ) or ( (self.rinv(y), self.rinv(x))  in idem):
                    answer.remove((x,y))
        if ( (self.b1, self.b2) in idem ) or ( (self.b2,self.b1) in idem ):
            if self.b1p < self.b2p:
                answer.remove((self.b1p, self.b2p))
            else:
                answer.remove((self.b2p, self.b1p))
        answer.sort()
        return answer

    def generate_generators(self):
        """Generates a list of generators of the DD module for this arcslide. Stored in self.generators
        Currently restricts to the middle SpinC structure.
        """
        answer = list()
        #Generate gen_idems, if necessary
        if self.gen_idems == None:
            self.generate_gen_idems()
        for (x,y) in self.gen_idems:
            yrev = list()
            for (m,n) in y:
                yrev.append((4*self.genus-n-1, 4*self.genus-m-1))
            answer.append(DDGen(repr(x)+'|'+repr(yrev),self.pmc_1, self.pmc_2rev,x,yrev))
        self.generators = answer

    def generate_gen_idems(self):
        """Generates a list of pairs of idempotents for generators for the DD module. Stored in self.gen_idems
        Currently restricts to the middle SpinC structure.
        """
        #Generators of type X (complementary)
        xidems = list()
        pmc1idems = self.pmc_1.idempotents()
        for idem in pmc1idems:
            xidems.append((idem, self.complementary_idem(idem)))
        #Generators of type Y (sub-complementary)
        yidems = list()
        for idem in pmc1idems:
            if ( (self.c1, self.c2) in idem ) or ( (self.c2, self.c1) in idem):
                if ( not ( (self.b1, self.b2) in idem)) and ( not ( (self.b2, self.b1) in idem) ):
                    scompidem = self.complementary_idem(idem)
                    if self.b1p<self.b2p:
                        scompidem.remove((self.b1p,self.b2p))
                    if self.b2p<self.b1p:
                        scompidem.remove((self.b2p,self.b1p))
                    if self.c1p < self.c2p:
                        scompidem.append((self.c1p,self.c2p))
                    if self.c2p < self.c1p:
                        scompidem.append((self.c2p,self.c1p))
                    scompidem.sort()
                    yidems.append((idem, scompidem))                        
        answer = xidems + yidems
        self.gen_idems = answer

    def dd_mod(self):
        "Return the type DD module for this arcslide."
        if self.ddstr != None:
            return self.ddstr
        if self.generators == None:
            self.generate_generators()
        self.generate_chords()
        self.basic_choice()
        diffs = dict()
        for x in self.generators:
            dx = DDElt({}, self.pmc_1, self.pmc_2rev)
            for y in self.generators:
                #Process the determinate chords.
                for (a,b) in self.DetChords+self.choice_chords:
                    #If a is compatible with x's left idempotents and y's left idempotents, and b is compatible with x's right idempotents and y's right idempotents, add corresponding elt to dx.
                    if l_idem_compat(self.pmc_1, a, x.idem_1) and l_idem_compat(self.pmc_2rev, self.opposite_strands(b), x.idem_2):
                        aalgelt = AlgElt(Strand_Diagram(self.pmc_1, a, x.idem_1))*AlgElt(Strand_Diagram(self.pmc_1, [], y.idem_1))
                        balgelt = AlgElt(Strand_Diagram(self.pmc_2rev, self.opposite_strands(b), x.idem_2))*AlgElt(Strand_Diagram(self.pmc_2rev, [], y.idem_2))
                        abalgelt = aalgelt**balgelt
                        if abalgelt:
                            dx = dx + DDElt({y:abalgelt}, self.pmc_1, self.pmc_2rev)
                #Now process the indeterminate chords -- write.
            diffs[x] = dx
        self.ddstr = TypeDDStr(self.pmc_1, self.pmc_2rev, self.generators, diffs)
        return self.ddstr

    def opposite_strands(self, strands):
        answer = list()
        for (x,y) in strands:
            answer.append((4*self.genus-y-1, 4*self.genus-x-1))
        return answer

#Method show() has been adjusted for overslides.
    def show(self):
        "Display graphically the pair (lchord, rchord). Intended for testing generate_U1() ... generate_U6()."
        pict = line([(0,-1), (0,4*self.genus)], rgbcolor = (0,0,0))
        pict = pict + line([(2,-1), (2,4*self.genus)], rgbcolor = (0,0,0))
        for i in range(4*self.genus):
            if i != self.b1:
                pict = pict + line([(0,i),(2,self.r(i))], rgbcolor = (1,0,0))
        if self.b1 < self.c1:
            pict = pict + line([(0,self.b1), (1,self.b1+1/2)], rgbcolor = (1,.25,.25)) + line([(1,self.b1p-1/2), (2,self.b1p)], rgbcolor = (1,.25,.25))
        if self.b1 > self.c1:
            pict = pict + line([(0,self.b1), (1,self.b1-1/2)], rgbcolor = (1,.25,.25)) + line([(1,self.b1p+1/2), (2,self.b1p)], rgbcolor = (1,.25,.25))
        #Show the matchings on the two sides:
        #left:
        for l in range(len(self.pmc_1.matching)):
            pict = pict + line([(-0.2, self.pmc_1.matching[l][0]), (-0.3-l/(8*self.genus), self.pmc_1.matching[l][0])], rgbcolor = (.2,.2,.2))
            pict = pict + line([(-0.2, self.pmc_1.matching[l][1]), (-0.3-l/(8*self.genus), self.pmc_1.matching[l][1])], rgbcolor = (.2,.2,.2))
            pict = pict + line([(-0.3-l/(8*self.genus), self.pmc_1.matching[l][0]), (-0.3-l/(8*self.genus), self.pmc_1.matching[l][1])], rgbcolor = (.2,.2,.2))
        #right:
        for l in range(len(self.pmc_2.matching)):
            pict = pict + line([(2.2, self.pmc_2.matching[l][0]), (2.3+l/(8*self.genus), self.pmc_2.matching[l][0])], rgbcolor = (.2,.2,.2))
            pict = pict + line([(2.2, self.pmc_2.matching[l][1]), (2.3+l/(8*self.genus), self.pmc_2.matching[l][1])], rgbcolor = (.2,.2,.2))
            pict = pict + line([(2.3+l/(8*self.genus), self.pmc_2.matching[l][0]), (2.3+l/(8*self.genus), self.pmc_2.matching[l][1])], rgbcolor = (.2,.2,.2))
        pict.show(axes=False) #, xmin=-2,xmax=2, ymin=-.5, ymax=4*self.genus-.5)


    def display_chord(self, chords):
        "Display graphically the pair (lchord, rchord). Intended for testing generate_U1() ... generate_U6()."
        lchord = chords[0]
        rchord = chords[1]
        #First show the arcslide
        pict = line([(0,-1), (0,4*self.genus)], rgbcolor = (0,0,0))
        pict = pict + line([(2,-1), (2,4*self.genus)], rgbcolor = (0,0,0))
        for i in range(4*self.genus):
            if i != self.b1:
                pict = pict + line([(0,i),(2,self.r(i))], rgbcolor = (1,0,0))
        if self.b1 < self.c1:
            pict = pict + line([(0,self.b1), (1,self.b1+1/2)], rgbcolor = (1,.25,.25)) + line([(1,self.b1p-1/2), (2,self.b1p)], rgbcolor = (1,.25,.25))
        if self.b1 > self.c1:
            pict = pict + line([(0,self.b1), (1,self.b1-1/2)], rgbcolor = (1,.25,.25)) + line([(1,self.b1p+1/2), (2,self.b1p)], rgbcolor = (1,.25,.25))
        #light lines going farther out, for clarity.
        for x in range(4*self.genus):
            pict = pict + line([(-1,x),(0,x)], rgbcolor = (.8,.8,.8))
            pict = pict + line([(2,x),(3,x)], rgbcolor = (.8,.8,.8))
        #Show the matchings on the two sides:
        #left:
        for l in range(len(self.pmc_1.matching)):
            pict = pict + line([(-0.2, self.pmc_1.matching[l][0]), (-0.3-l/(8*self.genus), self.pmc_1.matching[l][0])], rgbcolor = (.2,.2,.2))
            pict = pict + line([(-0.2, self.pmc_1.matching[l][1]), (-0.3-l/(8*self.genus), self.pmc_1.matching[l][1])], rgbcolor = (.2,.2,.2))
            pict = pict + line([(-0.3-l/(8*self.genus), self.pmc_1.matching[l][0]), (-0.3-l/(8*self.genus), self.pmc_1.matching[l][1])], rgbcolor = (.2,.2,.2))
        #right:
        for l in range(len(self.pmc_2.matching)):
            pict = pict + line([(2.2, self.pmc_2.matching[l][0]), (2.3+l/(8*self.genus), self.pmc_2.matching[l][0])], rgbcolor = (.2,.2,.2))
            pict = pict + line([(2.2, self.pmc_2.matching[l][1]), (2.3+l/(8*self.genus), self.pmc_2.matching[l][1])], rgbcolor = (.2,.2,.2))
            pict = pict + line([(2.3+l/(8*self.genus), self.pmc_2.matching[l][0]), (2.3+l/(8*self.genus), self.pmc_2.matching[l][1])], rgbcolor = (.2,.2,.2))
        #Now draw the chord on the two sides.
        for (x,y) in lchord:
            pict = pict + line([(-1,x),(-0.1, y)], rgbcolor = (0,0,1))
        for (x,y) in rchord:
            pict = pict + line([(2.1,x),(3, y)], rgbcolor = (0,0,1))
        pict.show(axes=False) #, xmin=-2,xmax=2, ymin=-.5, ymax=4*self.genus-.5)

class Arcslide(object):
    """Class for computing DD modules of arcslides. Contains an instance of either Overslide or Underslide, as appropriate.

    Input:
    pmc: The pointed matched circle you're arc-sliding.
    i: the moving endpoint of the arc-slide
    j: the endpoint i is sliding over

    Methods:
    dd_mod: generates and returns the type DD module for this underslide. The DD module is also stored; Calling dd_mod() only generates it the first time; afterwards, the existing dd structure is returned. (So, subsequent calls are quick.)
    show(): display this arc-slide graphically.

    Other useful data:
    self.pmc_2: the pointed matched circle resulting from the arc-slide.
    self.pmc_2rev: the orientation-reverse of self.pmc_2.

    Caveats:
    Currently restricts to the central SpinC structure, though this would be relatively easy to change.
    """
    def __init__(self, pmc, i, j=i+1):
        self.pmc_1 = pmc
        self.pmc_2 = pmc.arcslide(i,j)
        self.i=i
        self.j=j
        if self.pmc_1.is_overslide(i,j):
            self.arcslide = Overslide(self.pmc_1,i,j)
        if self.pmc_1.is_underslide(i,j):
            self.arcslide = Underslide(self.pmc_1,i,j)
        self.pmc_2rev = self.arcslide.pmc_2rev
        self.b1 = self.arcslide.b1
        self.b1p = self.arcslide.b1p
        self.c1 = self.arcslide.c1
        self.c1p = self.arcslide.c1p
        self.c2 = self.arcslide.c2
        self.c2p = self.arcslide.c2p

    def __repr__(self):
        if self.is_overslide():
            answer = "Overslide of "
        else:
            answer = "Underslide of "
        answer = answer + repr(self.i)+" over "+repr(self.j)+" of PMC "+repr(self.pmc_1)+"\n"
        return answer

    def dd_mod(self):
        return self.arcslide.dd_mod()

    def dd_mod_bar(self):
        return self.inverse().dd_mod().exchangeLR()
    
    def show(self):
        self.arcslide.show()

    def is_overslide(self):
        "Returns True if this is an overslide, False if an underslide."
        if self.pmc_1.is_overslide(self.i,self.j):
            return True
        return False

    def is_underslide(self):
        "Returns True if this is an underslide, False if an overslide."
        if self.pmc_1.is_overslide(self.i,self.j):
            return False
        return True

    def slide_type(self):
        "Returns 'Overslide' if this is an overslide, 'Underslide' if this is an underslide."
        if self.is_overslide():
            return "Overslide"
        return "Underslide"

    def inverse(self):
        "Returns the inverse arcslide to self."
        return Arcslide(self.pmc_2, self.b1p, self.c2p)
