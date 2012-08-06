def linear_hb(g, print_status=False):
    "Returns the type D structure for the 'linear handlebody' of genus g, by doing a sequence of arc-slides on the 0-franed split handlebody."
    hb = dict()
    hb[0]=infty_type_D(g)
#    hb[0]=zero_type_D(g)
    slideseq = list()
    for i in range(g-1):
        slideseq.extend([(4*i+3,4*i+4),(4*i+6,4*i+7),(4*i+5,4*i+6)])
    arcslideseq = list()
    for i in range(len(slideseq)):
        if print_status:
            print "Performing slide number "+repr(i+1)+" of "+repr(3*(g-1))+"."
        hb[i+1] = Arcslide(hb[i].pmc,slideseq[i][0],slideseq[i][1]).dd_mod().mor_to_d(hb[i])
        hb[i+1].simplify()
        hb[i+1].shorten_names()
    return hb[3*(g-1)]

def linear_pmc(g):
    "Returns the linear pointed matched circle of genus g."
    matching = list()
    matching.append((0,2))
    for i in range(2*g-2):
        matching.append((2*i+1,2*i+4))
    matching.append((4*g-3,4*g-1))
    return PMC(g,matching)


class BraidGrp(object):
    """A class for doing computations of branched double covers of braids on num_strands strands. num_strands must be even.

    Methods: 
    ----generator(n): returns a list of Arcslides corresponding go sigma_n if n is positive, and sigma_{|n|}^{-1} if n is negative.
    ----hf_of_knot(braidword, print_status=False): returns the total rank of HF-hat of the branched double cover of the plat closure of braidword.  Here, braidword is a list of sigma_i's.  Passing print_status=True will result in verbose output, saying what it's currently computing.

    Example:
    sage: mygp = BraidGrp(6)
    sage: mygp.hf_of_knot(5*[1]+3*[3]+2*[-5])
    1
    #Computes HF-hat of the branched double cover of the (-2,3,5) pretzel knot (which gives the Poincare homology sphere.)    
    """
    def __init__(self, num_strands):
        self.num_strands = num_strands
        self.genus = (num_strands - 2)/2
        self.computed_gens = dict()
        self.pmc = linear_pmc(self.genus)
        self.handlebody = None

    def generator(self, n):
        """Return a list of Arcslides for the generator n of self. If n is positive, this is sigma_n.  If n is negative, this is sigma_n^{-1}."""
        if n in self.computed_gens.keys():
            return self.computed_gens[n]
        absn = abs(n)
        if absn == 1:
            abslidelist = [(1,0)]
        if 1<absn<self.num_strands-2:
            abslidelist = [(2*(absn-1),2*(absn-1)-1), (2*(absn-1),2*(absn-1)-1)]
        if absn == self.num_strands-2: #Just added this case; thanks Bohua.
            abslidelist = [(2*(absn-1),2*(absn-1)-1)]
        if absn == self.num_strands-1:
            temp = []
            for i in range(self.genus):
                temp = [(4*i+1,4*i)]+temp
            for i in range(2*self.genus-2):  
                temp.append((2*i+2,2*i+1))
            abslidelist = temp
        absarcslidelist = list()
        absarcslidelist.append(Arcslide(self.pmc,abslidelist[0][0],abslidelist[0][1]))
        for i in range(1,len(abslidelist)):
            absarcslidelist.append(Arcslide(absarcslidelist[i-1].pmc_2,abslidelist[i][0], abslidelist[i][1]))
        if n>0:
            self.computed_gens[n]=absarcslidelist
            return absarcslidelist
        if n<0:
            slidelist = list()
            for slide in absarcslidelist:
                slidelist = [slide.inverse()]+slidelist
            self.computed_gens[n]=slidelist
            return slidelist

    def hf_of_knot(self, braidword, print_status=False):
        """Returns the rank of HF-hat of the branched double cover of the plat closure of braidword.

        Braidword: a list of integers [a,b,c,d,...]. a stands for sigma_a if a is positive, and sigma_a^{-1} if a is negative.
        """
    #Keep track of which slides we've already used, so we don't re-compute their DD modules.
        if self.handlebody == None:
            if print_status:
                print "Generating handlebody for capping."
            self.handlebody = linear_hb(self.genus, print_status)
        current_hb = self.handlebody
        k=0
        for b in braidword:
            k=k+1
            i=0
            for s in self.generator(b):
                if print_status:
                    i=i+1
                    print "Performing handleslide "+repr(i)+" of "+repr(len(self.generator(b)))+" for generator "+repr(k)+" of "+repr(len(braidword))+ " (type s_"+repr(b)+")."
                current_hb = s.dd_mod().mor_to_d(current_hb)
                current_hb.simplify()
                current_hb.shorten_names()
        if print_status:
            print "Finishing by capping off."
        return self.handlebody.mor_to_d(current_hb).homology()



