#Different in this version: algebra elements now have a parameter "spinc", the number of strands - g.

#Still to implement:
#Strand_Diagram: small_grading
#SmallGradingGroup: everything.
#Document AlgElt

class PMC(object):
    """A pointed matched circle.

    Inputs:
    genus: the genus of the underlying surface. Number of points is 4*genus.
    matching: A list of pairs (two tuples) of elements of range(4*genus).

    Data:
    genus: same as input 'genus.'
    matching: same as input 'matching', but stored as a tuple.
    matching_inv: a function which takes points of the pmc to their matched pairs.
    basis_computed: True if we have already computed a basis for this PMC, false otherwise.
    basis[i]: the basis for the algebras associated to this PMC, in idempotent i. Is *not* computed when pmc is initialized, since it's slow and takes lots of memory.
    trunc_basis_computed: True if we have already computed a truncated basis for this PMC, false otherwise.
    trunc_basis: those elements of basis with multiplicity at most 1 everywhere.
    """

    def __init__(self,genus,matching):
        self.genus=genus
        #Sanitize the matching, so that each pair is smaller number first.
        better_matching = list()
        for (i,j) in matching:
            if i<j:
                better_matching.append((i,j))
            if j<i:
                better_matching.append((j,i))
            if i==j:
                raise Exception("Can't have a point matched to itself!")
        better_matching.sort()
        self.matching=tuple(better_matching)
        self.basis_computed = dict()
        for i in range(-genus,genus+1):
            self.basis_computed[i]=False
        self.basis = dict()
        self.trunc_basis_computed = dict()
        for i in range(-genus,genus+1):
            self.trunc_basis_computed[i]=False
        self.trunc_basis = dict()
        self.matching_inv=dict()
        for s in self.matching:
            self.matching_inv[s[0]]=s[1]
            self.matching_inv[s[1]]=s[0]

    def __repr__(self):
        return 'Genus '+repr(self.genus)+' pointed matched circle with matching '+repr(self.matching)

    def __eq__(self,other):
        #Slightly wrong: order in each pair of matching shouldn't matter, either.
        if (other.genus==self.genus) and (Set(other.matching) == Set(self.matching)):
            return True
        return False
    
    def __ne__(self,other):
        if self.__eq__(other):
            return False
        return True

    def __hash__(self):
        return hash(self.matching)

    def selftest(self):
        "Test whether given data is, in fact, a pointed matched circle."
        #Check if length is right.
        if len(self.matching) != 2*self.genus:
            raise Exception("Matching not the right length.")
        #Check each strand appears exactly once in matching.
        points_list = range(4*self.genus)
        for m in self.matching:
            for i in m:
                if i not in points_list:
                    raise Exception("Point "+i+"appears twice in matching, or is outside matching range.")
                points_list.remove(i)
        #Check matching leads to a connected surface. - WRITE
        return True

    def which_pair(self, arc):
        "Return which pair (element of self.matching) the arc belongs to."
        for m in self.matching:
            if arc in m:
                return m
        raise Exception("No idempotent seems to contain arc.")

    def matched_point(self, point):
        "Return the point matched with the given one."
        pair = self.which_pair(point)
        if point == pair[0]:
            return pair[1]
        return pair[0]

    def idempotents(self, spinc=0):
        """Returns a list of the idempotents for this pointed matched circle with g+i strands.
        Idempotents are given as lists of elements of matching."""
        return Combinations(self.matching,self.genus+spinc).list()

    def compute_basis(self, spinc=0):
        "Compues a basis for Alg(PMC, i=spinc). Stores in self.basis[spinc]"
        idems=self.idempotents(spinc)
        self.basis[spinc] = list()
        for [i,j] in CartesianProduct(idems,idems).list():
            ijcompat = list()
            ijcompat.append([i,j,[]])
            for t in range(len(i)):
                ijcompatnew=list()
                for x in ijcompat:
                    for termPair in x[1]:
                        for termPoint in termPair:
                            if x[0][0][0]<termPoint:
                                newi = list(x[0])
                                newi.remove(x[0][0])
                                newj = list(x[1])
                                newj.remove(termPair)
                                newstrands = list(x[2])
                                newstrands.append((x[0][0][0],termPoint))
                                ijcompatnew.append([newi, newj, newstrands])
                            if x[0][0][0]==termPoint:
                                newi = list(x[0])
                                newi.remove(x[0][0])
                                newj = list(x[1])
                                newj.remove(termPair)
                                newstrands = list(x[2])
                                ijcompatnew.append([newi, newj, newstrands])
                            if x[0][0][1]<termPoint:
                                newi = list(x[0])
                                newi.remove(x[0][0])
                                newj = list(x[1])
                                newj.remove(termPair)
                                newstrands = list(x[2])
                                newstrands.append((x[0][0][1],termPoint))
                                ijcompatnew.append([newi, newj, newstrands])
                    ijcompat = ijcompatnew
                ijcompat=list(ijcompatnew)
            for x in ijcompat:
                self.basis[spinc].append(Strand_Diagram(self, x[2], left_idem=i, right_idem=j))                  
        self.basis_computed[spinc] = True

    def alg_basis(self, spinc=0):
        "Returns a basis for Alg(PMC)."
        if self.basis_computed[spinc] == False:
            self.compute_basis(spinc)
        return list(self.basis[spinc])

    def compute_trunc_basis(self, spinc=0):
        "Computes a basis for the truncated algebra (no multiplicities > 1). Stores in self.trunc_basis. Currently very inefficient."
        if self.basis_computed[spinc] == False:
            self.compute_basis(spinc)
        self.trunc_basis=list()
        for x in self.basis[spinc]:
            if x.is_in_trunc():
                self.trunc_basis[spinc].append(x)
        self.trunc_basis_computed[spinc] = True

    def alg_trunc_basis(self,spinc=0):
        "Returns a basis for the truncated algebra (no multiplicities > 1)."
        if self.trunc_basis_computed[spinc] == False:
            self.compute_trunc_basis(spinc)
        return list(self.trunc_basis[spinc])

    def opposite(self):
        "Returns the PMC obtained by orientation-reversing self."
        revmatch = list()
        for (i,j) in self.matching:
            revmatch.append((4*self.genus-j-1,4*self.genus-i-1))
        return PMC(self.genus, revmatch)

    def complementary_idem(self, idem):
        "Returns the idempotent in self.opposite() complementary to i."
        answer = list()
        for (i,j) in self.matching:
            if not (i,j) in idem:
                answer.append((4*self.genus-j-1,4*self.genus-i-1))
        return answer

    def propagated_left_equals_right(self, strands, left_idem, right_idem):
        ri = list(left_idem)
        for rho in strands:
            if self.which_pair(rho[0]) not in ri:
                raise Exception("Unable to propagate left idempotent: missing element.")
            ri.remove(self.which_pair(rho[0]))
        for rho in strands:
            if self.which_pair(rho[1]) in ri:
                raise Exception("Unable to propagate left idempotent: repeated element.")
            ri.append(self.which_pair(rho[1]))
        if (right_idem !=[]) and (right_idem != ri):
            return False
        return True

    def show(self):
        circ = line([(0,-.5),(0,4*self.genus-.5)], rgbcolor=(0,0,0))
        pict = circ
        for i in range(4*self.genus):
            pict = pict + line([(-.25,i),(.25,i)], rgbcolor=(.1,.1,.1))
        for pair in self.matching:
            j=self.matching.index(pair)
            offset = j/(2*self.genus)
            pict = pict+line([(-.3,pair[0]),(-.5-offset,pair[0]),(-.5-offset,pair[1]),(-.3,pair[1])], rgbcolor=(.3-j/(5*self.genus),.6+j/(5*self.genus),.3-j/(5*self.genus)))
        pict.show(axes=False, xmin=-2,xmax=2, ymin=-.5, ymax=4*self.genus-.5)

    def zero(self):
        "Return the zero element of self."
        return AlgElt([],self)

    def chords(self, spinc=0):
        "Return sum_(xi a chord)a(xi), in specified strands grading. Defaults to middle strands grading."
        arc_list = list()
        for i in range(4*self.genus):
            for j in range(i+1,4*self.genus):
                arc_list.append((i,j))
        answer = self.zero()
        for xi in arc_list:
            answer = answer + alg_element(self,[xi],spinc)
        return answer

    def arcslide(self, i, j):
        "Return the PMC obtained from self by doing an arcslide of arc end i over arc ending in j.  i and j must be consecutive."
        if abs(i-j) != 1:
            raise Exception("Can only slide a handle over a consecutive one.")
        answer_matching = list()
        jp = self.matched_point(j)
        if i<jp:
            for (a,b) in self.matching:
                if (a < i) or (a>jp):
                    newa = a
                if i<a<jp:
                    newa = a-1
                if (b < i) or (b>jp):
                    newb = b
                if i<b<jp:
                    newb = b-1
                if self.is_overslide(i,j):
                    if a == i:
                        newa = jp
                    if a == jp:
                        newa = jp-1
                    if b == i:
                        newb = jp
                    if b == jp:
                        newb = jp-1
                if self.is_underslide(i,j):
                    if a == i:
                        newa = jp-1
                    if a == jp:
                        newa = jp
                    if b == i:
                        newb = jp-1
                    if b == jp:
                        newb = jp                
                answer_matching.append((newa,newb))
        if i>jp:
            for (a,b) in self.matching:
                if (a > i) or (a<jp):
                    newa = a
                if jp<a<i:
                    newa = a+1
                if (b > i) or (b<jp):
                    newb = b
                if jp<b<i:
                    newb = b+1
                if self.is_overslide(i,j):
                    if a == i:
                        newa = jp
                    if a == jp:
                        newa = jp+1
                    if b == i:
                        newb = jp
                    if b == jp:
                        newb = jp+1
                if self.is_underslide(i,j):
                    if a == i:
                        newa = jp+1
                    if a == jp:
                        newa = jp
                    if b == i:
                        newb = jp+1
                    if b == jp:
                        newb = jp
                answer_matching.append((newa,newb))
        return PMC(self.genus, answer_matching)

    def is_underslide(self, i, j):
        "Returns True if sliding i over j is an overslide, and False if it's an underslide."
        jp = self.matched_point(j)
        if ((jp < i) and (i<j)) or ( (j<i) and (i<jp)):
            return True
        return False

    def is_overslide(self, i, j):
        "Returns True if sliding i over j is an overslide, and False if it's an underslide."
        jp = self.matched_point(j)
        if ((jp < i) and (j<i)) or ( (i<jp) and (i<j)):
            return True
        return False

    #Some particular pmc's.
    @classmethod
    def split_matching(cls, g):
        match = list()
        for i in range(g):
            match.append((4*i,4*i+2))
            match.append((4*i+1,4*i+3))
        return PMC(g,match)

    @classmethod
    def antipodal_matching(cls, g):
        match = list()
        for i in range(2*g):
            match.append((i,i+2*g))
        return PMC(g,match)

    
    def infty_type_D(self):
        "Returns the standard type D structure for an infinity-framed handlebody of genus k."
        xidems = list()
        for i in range(self.genus):
            xidems.append((4*i+1,4*i+3))
        diffs = list()
        for i in range(self.genus):
            diffs.append(Strand_Diagram(self, [(4*i+1,4*i+3)], left_idem=xidems))
        return TypeDStr(self, {'x':xidems}, {'x':[(diffs,'x')]})

    
    def infty_type_D_bd(cls, k):
        "Returns a bounded type D structure for an infinity-framed handlebody of genus k."
        pass


    def zero_type_D(self):
        "Returns the standard type D structure for an zero-framed handlebody of genus k."
        xidems = list()
        for i in range(self.genus):
            xidems.append((4*i,4*i+2))
        diffs = list()
        for i in range(self.genus):
            diffs.append(Strand_Diagram(self, [(4*i,4*i+2)], left_idem=xidems))
        return TypeDStr(self, {'y':xidems}, {'y':[(diffs,'y')]})

    
    def zero_type_D_bd(cls, k):
        "Returns a bounded type D structure for an zero-framed handlebody of genus k."
        pass

    
    def m_one_type_D(self):
        "Returns the standard type D structure for a -1-framed handlebody of genus k."
        gens = [[]]
        for i in range(self.genus):
            new_gens=list()
            for x in gens:
                new_gens.append(x+[(4*i,4*i+2)])
                new_gens.append(x+[(4*i+1,4*i+3)])
            gens = new_gens
        diffs = dict()
        for x in gens:
            ans = list()
            for i in range(self.genus):
                if x[i]==(4*i,4*i+2):
                    appendme = list(x)
                    appendme[i] = (4*i+1,4*i+3)##EDIT HERE.
                    ans.append(([Strand_Diagram(self, [(4*i,4*i+1)], left_idem=x), Strand_Diagram(self,[(4*i+2,4*i+3)],left_idem=x)],gens.index(appendme)))
            diffs[gens.index(x)]=ans
        basis = dict()
        for x in gens:
            basis[gens.index(x)] = x
        return TypeDStr(self,basis,diffs)



class Strand_Diagram(object):
    """A strand diagram for a pointed matched circle.

    Elements of the algebra are to be stored as lists of Strand_Diagram's.

    Inputs:
    pmc: the pointed matched circle in which this is a strand diagram.
    strands: and list of pairs in points in the pmc (range(4*genus)).
    left_idem: the left idempotent of the diagram. Leave blank to have it computed from the right idempotent.
    right_idem: the right idempotent of the diagram. Leave blank to have it computed from the left idempotent.
    name: optionally, a name for this strand diagram; if given, the name is printed; if not, the key data is printed.
    spinc: the spinc structure. Number of elements of idempotent is pmc.genus + spinc. Defaults to 0.

    Example:
    my_pmc = split_matching(1)
    my_idemp = my_pmc.which_pair(1)
    my_chord = Strand_Diagram(my_pmc, [(1,2)], left_idem=[my_idemp])
    """
    def __init__(self, pmc, strands=[], left_idem=[], right_idem=[], name=''):
        self.pmc=pmc; self.strands=strands; 
        #Propogate idempotents.
        if left_idem != []:
            self.name = name
            self.left_idem=left_idem
            ri = list(left_idem)
            for rho in strands:
                if self.pmc.which_pair(rho[0]) not in ri:
                    raise Exception("Unable to propagate left idempotent: missing element.")
                ri.remove(self.pmc.which_pair(rho[0]))
            for rho in strands:
                if self.pmc.which_pair(rho[1]) in ri:
                    raise Exception("Unable to propagate left idempotent: repeated element.")
                ri.append(self.pmc.which_pair(rho[1]))
            if (right_idem !=[]) and (Set(right_idem) != Set(ri)):
                print "Propogated: ",ri,"  Given: ",right_idem
                raise Exception("Propogated right idempotent, does not agree with given right idempotent.")
            self.right_idem=ri
        else: 
            self.name = name
            self.right_idem=right_idem
            li = list(right_idem)
            reverse_strands=list(strands)
            reverse_strands.reverse()
            for rho in reverse_strands:
                if self.pmc.which_pair(rho[1]) not in li:
                    raise Exception("Unable to propagate right idempotent: missing element.")
                li.remove(self.pmc.which_pair(rho[1]))
                if self.pmc.which_pair(rho[0]) in li:
                    raise Exception("Unable to propagate right idempotent: repeated element.")
                li.append(self.pmc.which_pair(rho[0]))
            self.left_idem=li
        self.strands.sort()
        self.left_idem.sort()
        self.right_idem.sort()
        self.spinc=len(self.left_idem)
            
    def __repr__(self):
        if self.name:
            return self.name
        return ' | LI:'+repr(self.left_idem)+' S:'+repr(self.strands)+' RI:'+repr(self.right_idem)+' | '

    def __eq__(self, other):
        if not isinstance(other, Strand_Diagram):
            return False
        if (self.pmc==other.pmc) and (Set(self.strands)==Set(other.strands)) and (Set(self.left_idem)==Set(other.left_idem)):
            return True
        return False
    
    def __ne__(self,other):
        if self.__eq__(other):
            return False
        return True

    #This seems to work, but RL a little concerned this goes beyond his coding skills.
    def __hash__(self):
        return hash(tuple(self.strands)).__xor__(hash(tuple(self.left_idem)))

    def __mul__(self, other):
        return self.r_multiply(other)
    
    def name_me(self, name):
        self.name=name

    def multiplicities(self):
        "Returns the list of local multiplicities of self."
        local_mults = (4*self.pmc.genus-1)*[0]
        for rho in self.strands:
            for i in range(rho[0],rho[1]):
                local_mults[i]+=1
        return local_mults

    def is_in_trunc(self):
        "Returns True if none of the local multiplicities of self are >1 and False otherwise."
        local_mults = self.multiplicities()
        for i in range(4*self.pmc.genus-1):
            if local_mults[i]>1:
                return False
        return True

    def inv(self):
        "Returns the number of inversions of self. Used by big_grading()."
        answer=0
        for rho in self.strands:
            for sigma in self.strands[self.strands.index(rho):]:
                if rho[0]<sigma[0] and rho[1]>sigma[1]:
                    answer+=1
        return answer

    def iota(self):
        "Returns iota(self), as in LOT Definition 3.18. Used by big_grading."
        local_mults = self.multiplicities()
        m = 0/1
        reduced_starting_idem = list()
        for s in self.strands:
            reduced_starting_idem.append(s[0])
        for i in reduced_starting_idem:
            m = m + local_mults[i]/2
            if i != 0:
                m = m + local_mults[i-1]/2
        return self.inv()-m

    def big_grading(self):
        "Returns G' (big grading group) grading of self."
        answer = list()
        answer.append(self.iota())
        answer.extend(self.multiplicities())
        return answer

    def small_grading(self):
        "Returns G (small grading group) grading of self."
        pass

    #RL August 5 2011: fixed bug: differential() used to not check for double crossings. (Thanks, Bohua Zhan.)
    def differential(self):
        "Return the differential of self."
        answer = list()
        grgrp = BigGradingGroup(self.pmc)
        #First run over all of the moving strands.
        for x in self.strands:
            for y in self.strands[self.strands.index(x):]:
                if (x[0]<y[0]) and (x[1]>y[1]):
                    new_strands = list(self.strands)
                    new_strands[self.strands.index(x)]=(x[0],y[1])
                    new_strands[self.strands.index(y)]=(y[0],x[1])
                    append_me = Strand_Diagram(self.pmc, new_strands,left_idem=self.left_idem)
                    #Test whether this element had a double crossing.
                    if grgrp.multiply(append_me.big_grading(),grgrp.central()) == self.big_grading():
                        answer.append(append_me)
        #Now find the unmoving strands:
        moving_pts = list()
        for s in self.strands:
            moving_pts.append(s[0])
        unmoving_idem=list(self.left_idem)
        for i in list(unmoving_idem):
            for p in i:
                if (p in moving_pts) and (i in unmoving_idem):
                    unmoving_idem.remove(i)
        #And consider their contribution to differential
        for i in unmoving_idem:
            for x in self.strands:
                for p in i:
                    if (x[0]<p) and (p<x[1]):
                        new_strands=list(self.strands)
                        new_strands.remove(x)
                        new_strands.append((p,x[1]))
                        new_strands.append((x[0],p))
                        append_me = Strand_Diagram(self.pmc,new_strands,left_idem=self.left_idem)
                    #Test whether this element had a double crossing.
                        if grgrp.multiply(append_me.big_grading(),grgrp.central()) == self.big_grading():
                            answer.append(append_me)
        return answer

    def r_multiply(self, other):
        "Return the product self times other, where other is a Strand_Diagram."
        #Do the idempotents agree?
        if self.right_idem != other.left_idem:
            return []
        right_endpoints = list()
        for s in self.strands:
            right_endpoints.append(s[1])
        unmoving_arcs = list()
        for x in self.right_idem:
            if (not (x[0] in right_endpoints)) and (not (x[1] in right_endpoints)):
                unmoving_arcs.append(x)
        unmoving_points = list()
        for x in unmoving_arcs:
            unmoving_points.append(x[0])
            unmoving_points.append(x[1])
        #Do the idempotents agree but the strands not match up?
        for s in other.strands:
            if (not (s[0] in right_endpoints)) and (not (s[0] in unmoving_points)):
                return []
        #Join the strands that can be joined; take union of rest.
        answer_strands = list(self.strands)
        for s in other.strands:
            if s[0] in right_endpoints:
                for x in self.strands:
                    if x[1]==s[0]:
                        answer_strands.remove(x)
                        answer_strands.append((x[0],s[1]))
            else:
                answer_strands.append(s)
        #Now check: were there double crossings?
        answer=Strand_Diagram(self.pmc, answer_strands, self.left_idem)
        grgrp = BigGradingGroup(self.pmc)
        if grgrp.multiply(self.big_grading(),other.big_grading()) != answer.big_grading():
            return []
        return answer
    
    def l_multiply(self, other):
        "Return the product other times self."
        return other.r_multiply(self)

    def augmentation(self):
        "Returns augmentation(self), i.e., self if self is an idempotent, empty list otherwise."
        if self.strands == []:
            return self
        return []

    def id_minus_aug(self):
        "Returns (Id - augmentation)(self), i.e., self if self is not an idempotent, empty list otherwise."
        if self.strands != []:
            return self
        return []

    def show(self):
        #First draw the pointed matched circle.
        circ = line([(0,-.5),(0,4*self.pmc.genus-.5)], rgbcolor=(0,0,0))
        pict = circ
        for i in range(4*self.pmc.genus):
            pict = pict + line([(-.25,i),(.25,i)], rgbcolor=(.1,.1,.1))
        for pair in self.pmc.matching:
            j=self.pmc.matching.index(pair)
            offset = j/(2*self.pmc.genus)
            pict = pict+line([(-.3,pair[0]),(-.5-offset,pair[0]),(-.5-offset,pair[1]),(-.3,pair[1])], rgbcolor=(.3-j/(5*self.pmc.genus),.6+j/(5*self.pmc.genus),.3-j/(5*self.pmc.genus)))
        #Now draw the moving strands
        for s in self.strands:
            i=self.strands.index(s)
            pict = pict+line([(.3,s[0]),(1.7,s[1])],rgbcolor=(1-i/(6*self.pmc.genus),.1+i/(6*self.pmc.genus),.1+i/(6*self.pmc.genus)))
        #Find unmoving strands
        moving_pts = list()
        for s in self.strands:
            moving_pts.append(s[0])
        unmoving_idem=list(self.left_idem)
        for i in list(unmoving_idem):
            for p in i:
                if (p in moving_pts) and (i in unmoving_idem):
                    unmoving_idem.remove(i)
        #Now draw the non-moving strands
        for i in unmoving_idem:
            j=unmoving_idem.index(i)
            pict = pict+line([(.3,i[0]), (1.7,i[0])], linestyle='--', rgbcolor=(.1+j/(6*self.pmc.genus),.1+j/(6*self.pmc.genus),1-j/(6*self.pmc.genus)))
            pict = pict+line([(.3,i[1]), (1.7,i[1])], linestyle='--', rgbcolor=(.1+j/(6*self.pmc.genus),.1+j/(6*self.pmc.genus),1-j/(6*self.pmc.genus)))
        #Also draw circle on RHS
        pict = pict+line([(2,-.5),(2,4*self.pmc.genus-.5)], rgbcolor=(0,0,0))
        for i in range(4*self.pmc.genus):
            pict = pict + line([(1.75,i),(2.25,i)], rgbcolor=(.1,.1,.1))
        for pair in self.pmc.matching:
            j=self.pmc.matching.index(pair)
            offset = j/(2*self.pmc.genus)
            pict = pict+line([(2.3,pair[0]),(2.5+offset,pair[0]),(2.5+offset,pair[1]),(2.3,pair[1])], rgbcolor=(.3-j/(5*self.pmc.genus),.6+j/(5*self.pmc.genus),.3-j/(5*self.pmc.genus)))
        pict.show(axes=False)

    def opposite(self):
        "Returns the Strand_Diagram with the same multiplicities and idempotents as self but viewed as in opposite PMC."
        rev_right_idem = list()
        for (i,j) in self.left_idem:
            rev_right_idem.append((4*self.pmc.genus-1-j,4*self.pmc.genus-1-i))
        rev_left_idem = list()
        for (i,j) in self.right_idem:
            rev_left_idem.append((4*self.pmc.genus-1-j,4*self.pmc.genus-1-i))
        rev_strands = list()
        for (i,j) in self.strands:
            rev_strands.append((4*self.pmc.genus-1-j,4*self.pmc.genus-1-i))
        return Strand_Diagram(self.pmc.opposite(),rev_strands,rev_left_idem,rev_right_idem)

from UserList import UserList
class AlgElt(UserList):
    """An element of Alg(PMC). Roughly, a list of Strand_Diagram's.
    Inputs:
    data: a list or tuple of Strand_Diagram's. Should all be over the same PMC.
    pmc (optional): pointed matched cirlce underlying this AlgElt. If data != [] then this is irrelevant. If data==[] then it's good (but not forced) to specify pmc.

    Alternately, data can be a single Strand_Diagram. In this case, pmc argument is ignored.

    Methods: 

    Examples:
    """
    def __init__(self, data, pmc=[]):
        if isinstance(data, Strand_Diagram):
            self.pmc = data.pmc
            UserList.__init__(self, [data])
        else:
            UserList.__init__(self, data)
            if pmc:
                self.pmc=pmc
            else:
                if data:
                    self.pmc=data[0].pmc
            if data:
                for x in data:
                    if self.pmc != x.pmc:
                        raise Exception("All elements in AlgElt must come from same PMC.")
            self.reduce_mod_2()

    def __repr__(self):
        return "<"+repr(self.data)+">"

    def __mul__(self, other):
        if not isinstance(other, AlgElt):
            return other.__rmul__(self)
        answer = list()
        for x in self:
            for y in other:
                if x.r_multiply(y):
                    answer.append(x.r_multiply(y))
        return AlgElt(answer)

    def __add__(self, other):
        answer = self.data + other.data
        return AlgElt(reduce_mod_2(answer))

    def __eq__(self, other):
        #For convenience, the element 0 is also equal to the empty list.
        if not self.data:
            if not other:
                return True
            return False
        if not isinstance(other, AlgElt):
            return False
        if Set(self.data)==Set(other.data):
            return True
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    #Probably works... though mutable objects shouldn't be hashable...
    def __hash__(self):
        return Set(self.data).__hash__()

    def __pow__(self, other):
        "a**b now gives a\otimes b, an AlgBlgElt."
        if not isinstance(other, AlgElt):
            raise Exception("What used to be the power operation takes two algebra elements and gives their external tensor product, and AlgBlgElt.")
        answer = AlgBlgElt({},self.pmc, other.pmc)
        for x in self.data:
            answer = answer + AlgBlgElt({x:other})
        return answer

    def diff(self):
        "Returns differential(self)."
        answer = list()
        for x in self:
            answer.extend(x.differential())
        return AlgElt(answer, self.pmc)

    def differential(self):
        "Returns differential(self)."
        return self.diff()

    def reduce_mod_2(self):
        "Deletes pairs of duplicates. Does not try to preserve order."
        answer = list()
        for x in self.data:
            if x in answer:
                answer.remove(x)
            else:
                answer.append(x)
        self.data = answer

    def pmc(self):
        if not self.data:
            return False
    def opposite(self):
        "Returns the AlgElt obtained by taking opposite of all Strand_Diagram's in self."
        rev_data = list()
        for x in self.data:
            rev_data.append(x.opposite())
        return AlgElt(rev_data, self.pmc.opposite())

    def is_invertible(self):
        "Returns True is self is an idempotent and False otherwise. (SHOULD REVISE.)"
        if len(self.data)==0:
            return False
        for x in self.data:
            if x != x*x:
                return False
        return True
    
    def inverse(self):
        "Currently only works for idempotents."
        if not self.is_invertible():
            raise Exception("Can not invert non-invertible element.")
        return AlgElt(self.data, self.pmc)


#Some useful operations for algebras associated to pmc's.
def alg_element(pmc, strands, spinc=0):
    "Returns list (sum) of all algebra elements in A(pmc) with the given strands, and any consistent idempotents. Naively written at present."
    g = pmc.genus+spinc
    answer = list()
    if not isinstance(strands, list):
        raise Exception("Strands input to alg_element should be a list.")
    for i in pmc.idempotents(spinc):
        if l_idem_compat(pmc,strands,i):
            answer.append(Strand_Diagram(pmc,strands,left_idem=i))

    return AlgElt(answer, pmc) #Used to just return answer.

def l_idem_compat(pmc, strands, idem):
    "Check if strands are compatible with the left idempotent idem of pmc."
#    if idem not in pmc.idempotents():
#        raise Exception("Given idmpotent idem is not an idempotent for given matched circle pmc.")
    delete_list = list(idem)
    add_list = list()
    for s in strands:
        if (pmc.which_pair(s[0]) not in delete_list) or (pmc.which_pair(s[1]) in add_list):
            return False
        delete_list.remove(pmc.which_pair(s[0]))
        add_list.append(pmc.which_pair(s[1]))
    for i in delete_list:
        if i in add_list:
            return False
    return True

def r_idem_compat(pmc, strands, idem):
    "Check if strands are compatible with the right idempotent idem of pmc."
#    if idem not in pmc.idempotents():
#        raise Exception("Given idmpotent idem is not an idempotent for given matched circle pmc.")
    delete_list = list(idem)
    add_list = list()
    for s in strands:
        if (pmc.which_pair(s[1]) not in delete_list) or (pmc.which_pair(s[0]) in add_list):
            return False
        delete_list.remove(pmc.which_pair(s[1]))
        add_list.append(pmc.which_pair(s[0]))
    for i in delete_list:
        if i in add_list:
            return False
    return True


def mul_list(list1, list2):
    "Multiply two lists (sums) of Strand_Diagrams."
    answer = list()
    for x in list1:
        for y in list2:
            if x.r_multiply(y):
                answer.append(x.r_multiply(y))
    return answer

def alg_differentiate_list(x):
    "Differentiate all the elements in the given list x."
    answer = list()
    for s in x:
        answer.extend(s.differential())
    return answer

def aug_list(x):
    "Apply the augmentation to all elements of list; return result."
    answer = list()
    for s in x:
        if s.augmentation():
            answer.append(s.augmentation())
    return answer

def id_minus_aug_list(x):
    "Apply (Id-Augmentation) to all elements of list; return result."
    answer = list()
    for s in x:
        if s.id_minus_aug():
            answer.append(s.id_minus_aug())
    return answer

#Three self-tests for the algebra, in specified spinc structure.
def alg_d_squared_zero(pmc, spinc=0):
    "Return True if d-squared is zero for Alg(pmc), and False otherwise."
    for a in pmc.alg_basis(spinc):
        if reduce_mod_2(alg_differentiate_list(a.differential())):
            return False
    return True

def alg_is_leibniz(pmc, spinc=0):
    "Return True if Alg(pmc) satisfies the Leibniz rule, and False otherwise."
    for a in pmc.alg_basis(spinc):
        for b in pmc.alg_basis(spinc):
            test_list = mul_list(a.differential(),[b])
            other_order = mul_list([a],b.differential())
            if other_order:
                test_list.extend(other_order)
            if a.r_multiply(b):
                test_list.extend((a.r_multiply(b)).differential())
            if reduce_mod_2(test_list):
                return False
    return True

def alg_is_assoc(pmc,spinc=0):
    "Return True if Alg(pmc) is associative, and False otherwise."
    #Written very badly.
    for a in pmc.alg_basis(spinc):
        for b in pmc.alg_basis(spinc):
            for c in pmc.alg_basis(spinc):
                if (b.r_multiply(c)):
                    if a.r_multiply(b):
                        if not ((a.r_multiply(b.r_multiply(c))) == ((a.r_multiply(b)).r_multiply(c))):
                            return False
                    else:
                        if a.r_multiply(b.r_multiply(c)):
                            return False
                else:
                    if a.r_multiply(b):
                        if (a.r_multiply(b)).r_multiply(c):
                            return False
    return True

def alg_homology(pmc,spinc=0):
    "Returns the dimension of the homology of A(pmc,0)"
    vert = pmc.alg_basis(spinc)
#    print vert
#    for i in list(vert):
#        if (i.left_idem != [(0,2),(1,3)]) or (i.right_idem !=[(0,2),(1,3)]):
#            vert.remove(i)
    edges = dict()
    for v in vert:
        edges[v]=v.differential()
    chcx=ChainCx(edges)
    return chcx.homology()

def reduce_mod_2(elts):
    "Takes a list and deletes pairs of duplicates. Does not try to preserve order."
    answer = list()
    for x in elts:
        if x in answer:
            answer.remove(x)
        else:
            answer.append(x)
    return answer

class BigGradingGroup(object):
    """The big grading group associated to a pointed matched circle.
    Elements are lists [m,a_1,...] of length 4*genus where:
    m is the homological (Maslov) component of the grading and the
    a_i are the homological components of the gradings.
    The m and a_i should be Rational's or Integer's.
    Inputs: a pointed matched circle (class PMC)"""
    
    def __init__(self, pmc):
        self.pmc = pmc

    def is_element(self, x):
        "Check if x is an element of self."
        if not isinstance(x,list):
            return False
        if len(x) != 4*self.pmc.genus:
            return False
        for t in x:
            if (not isinstance(t,sage.rings.rational.Rational)) and (not isinstance(t,sage.rings.integer.Integer)):
                return False
        return True
    
    def multiply(self, x, y):
        if (not self.is_element(x)):
            raise Exception("First input is not an element of this BigGradingGroup.")
        if (not self.is_element(y)):
            raise Exception("Second input is not an element of this BigGradingGroup.")
        z=list(x)
        for i in range(4*self.pmc.genus):
            z[i]=x[i]+y[i]
        for i in range(1,4*self.pmc.genus-1):
            z[0]=z[0]+(x[i]*y[i+1]-x[i+1]*y[i])/2
        return z

    def inverse(self, x):
        if (not self.is_element(x)):
            raise Exception("Input is not an element of this BigGradingGroup.")
        z=list()
        for i in range(4*self.pmc.genus):
            z.append(-x[i])
        return z


    def multiply_list(self, list_of_elts):
        "Given a list [x_1,x_2,...] of elements of G' returns the element x_1*x_2*..."
        z=self.identity()
        for x in list_of_elts:
            z=self.multiply(z,x)
        return z

    def power(self, x, n):
        "Return x^n"
        if n==0:
            return self.identity()
        if n<0:
            to_raise=self.inverse(x)
            pow = -n
        if n>0:
            to_raise=x
            pow = n
        return self.multiply_list(pow*[to_raise])

    def identity(self):
        return (4*self.pmc.genus)*[0]

    def central(self):
        "Return the distinguished central element of this BigGradingGroup"
        return [1]+(4*self.pmc.genus-1)*[0]

class SmallGradingGroup(object):
    """The small grading group associated to a pointed matched circle.
    Inputs: a pointed matched circle (class PMC)"""
    
    def __init__(self, pmc):
        pass
    
    def multiply(self, x, y):
        pass

    def inverse(self, x):
        pass
    
    def multiply_list(self, list_of_elts):
        pass

    def central(self):
        pass



from UserDict import UserDict
class AlgBlgElt(UserDict):
    """An element of Alg(PMC1)\otimes Alg(PMC2). 
    Stored as a dictionary with keys Strand_Diagram's for PMC1 and values AlgElt's for PMC2.

    Inputs:
    data: A dictionary with keys Strand_Diagram's for PMC1 and values AlgElts for PMC2.
    pmc_1 (optional): first pointed matched cirlce underlying this AlgElt. If data != [] then this is irrelevant. If data==[] then it's good (but not forced) to specify pmc.
    pmc_2 (optional): second pointed matched cirlce underlying this AlgElt. If data != [] then this is irrelevant. If data==[] then it's good (but not forced) to specify pmc.

    Methods: 

    Examples:
    sage:  ab = AlgBlgElt({rho1:AlgElt([rho1,rho3]),rho2:AlgElt([rho2])},small_pmc,small_pmc)

    sage:  myelt = AlgBlgElt({rho1:AlgElt([rho1,rho2,rho3]),rho2:AlgElt([rho1,rho23])})
    sage:  myelt.basis_expansion()
    [(rho1, rho2), (rho1, rho3), (rho2, rho1), (rho2, rho23)]
    sage:  myelt*myelt
    {LI:[(0, 2)] S:[(0, 2)] RI:[(0, 2)]: [LI:[(0, 2)] S:[(0, 3)] RI:[(1, 3)]]}
    """
    def __init__(self, data, pmc_1=[], pmc_2=[]):
        UserDict.__init__(self, data)
        if pmc_1:
            self.pmc_1=pmc_1
        else:
            if data:
                self.pmc_1=data.keys()[0].pmc
        if pmc_2:
            self.pmc_2=pmc_2
        else:
            if data:
                self.pmc_2=data.values()[0].pmc
        self.reduce_mod_2()
        #Should do more sanity checking.

    def __add__(self, other):
        if (self.pmc_1 != other.pmc_1):
            raise Exception("Can't add elements with different pmc_1's.")
        if (self.pmc_2 != other.pmc_2):
            raise Exception("Can't add elements with different pmc_2's.")
        answer = dict()
        for x in self.data.keys():
            if x in other.data.keys():
                answer[x]=self.data[x]+other.data[x]
            else:
                answer[x]=self.data[x]
        for x in other.data.keys():
            if not (x in self.data.keys()):
                answer[x]=other.data[x]
        return AlgBlgElt(answer, self.pmc_1, self.pmc_2)

    def __mul__(self, other):
        if not (isinstance(other,AlgBlgElt)):
            return other.__rmul__(self)
        if (self.pmc_1 != other.pmc_1):
            raise Exception("Can't multiply elements with different pmc_1's.")
        if (self.pmc_2 != other.pmc_2):
            raise Exception("Can't multiply elements with different pmc_2's.")
        answer = dict()
        for x in self.data.keys():
            for y in other.data.keys():
                if x*y:
                    answer[x*y]=AlgElt([],self.pmc_2)
        for x in self.data.keys():
            for y in other.data.keys():
                if x*y:
                    answer[x*y]=answer[x*y]+self.data[x]*other.data[y]
        return AlgBlgElt(answer, self.pmc_1,self.pmc_2)

    def __eq__(self, other):
        #For convenience, the element 0 is also equal to the empty list.
        if not self.data:
            if not other:
                return True
            return False
        if not isinstance(other, AlgBlgElt):
            return False
        for x in self.data.keys():
            if not (x in other.data.keys()):
                if self.data[x]:
                    return False
            else:
                if self.data[x]!=other.data[x]:
                    return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

#    #Probably works... though mutable objects shouldn't be hashable...
#    def __hash__(self):
#        return Set(self.data).__hash__()

    def diff_1(self):
        "Differentiate part of self in PMC_1; returns result."
        #A bit slow: no need to initialize AlgBlgElt's each time.
        answer = AlgBlgElt({},self.pmc_1, self.pmc_2)
        for x in self.data.keys():
            for y in x.differential():
                answer = answer + AlgBlgElt({y:self.data[x]},self.pmc_1,self.pmc_2)
        return answer

    def diff_2(self):
        "Differentiate part of self in PMC_2; returns result."
        #A bit slow: no need to initialize AlgBlgElt's each time.
        answer = AlgBlgElt({},self.pmc_1, self.pmc_2)
        for x in self.data.keys():
            answer = answer + AlgBlgElt({x:self.data[x].diff()},self.pmc_1,self.pmc_2)
        return answer

    def reduce_mod_2(self):
        "Deletes pairs of duplicates."
        answer = dict()
        for x in self.data.keys():
            self.data[x].reduce_mod_2()
        for x in self.data.keys():
            if not self.data[x]:
                del self.data[x]

    def basis_expansion(self):
        "Returns a list of pairs (a_i,b_i) of strand diagrams so that sum_i a_i\otimes b_i = self"
        self.reduce_mod_2()
        answer = list()
        for x in self.data.keys():
            for y in self.data[x]:
                answer.append((x,y))
        return answer
    
    def show(self):
        "Display self graphically."
        pass

    def exchangeLR(self):
        "Returns result of exchanging pmc_1 and pmc_2, i.e., a\otimes b becomes b\otimes a."
        answer = AlgBlgElt([],self.pmc_2, self.pmc_1)
        for (a,b) in self.basis_expansion():
            answer = answer + (AlgElt(b)**AlgElt(a))
        return answer
