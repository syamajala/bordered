##Much of this code written by Jin Woo Jang, Rachel Vishnepolsky, Xuran Wang and Tom Peters
# Documentation
class ChainCx(DiGraph):
    """
    A class for chain complexes over Z/2.

    Extends DiGraph.

    EXAMPLES:
    sage: egCx=ChainCx({0:[1,2], 1:[3], 2:[3], 3:[]})
    sage: egCx.is_complex()
    True
    sage: egCx.homology()
    0
    
    Fails if d^2 is not zero:

    sage: ChainCx({0:[1,2], 1:[], 2:[3], 3:[]})
    Traceback (click to the left for traceback)
    ...
    Exception: Does not define a chain complex -- d squared not zero.
    """

    def __init__(self, data=None, pos=None, loops=None, format=None, boundary=[], weighted=None, implementation='networkx', sparse=True, vertex_labels=True, **kwds):
        DiGraph.__init__(self, data, pos, loops, format, boundary, weighted, implementation, sparse, vertex_labels, **kwds)
        if self.is_complex() != True:
#            print "d squared not zero."
            raise Exception("Does not define a chain complex -- d squared not zero.")

    def grandchildren(self,v):
        "Compiles a list of grandchildren of a vertex v"
        list=[]
        for w in self.neighbors_out(v):
            for y in self.neighbors_out(w):
                list.append(y)
        return list    	    


    def is_complex(self):
        """
        Checks if the graph is a complex
        """
        for v in self.vertices():
            if v in self.grandchildren(v):
                return False
            else:
                for w in self.grandchildren(v):
                    L = self.all_paths(v,w)
                    if len(L)%2!=0:
                        return False
        return True

    def homology(self):
        """
        Simplifies the graph by isolating an edge, directed from the vertex 'tail' to the vertex 'head:'
        Returns the dimension of the homology.
        """
    # Maintains a dynamic list of edges:
        edges=self.edges()
        while len(edges)>0:
            e=edges[0]
    # Indicates the starting and ending points of the edge we wish to isolate:
            head=e[1]
            tail=e[0]    
            
# Now we make a list of all the other children of the tail, and delete them. 
            sib=self.neighbors_out(tail)
            
            sib.remove(head)

            for y in sib:
                self.delete_edge(tail,y)
                
                edges.remove((e[0], y, None))

# Then we make a list of all the other parents of the head, and delete the edges from them to the head.
            par=self.neighbors_in(head)
            par.remove(tail)

            for y in par:
                self.delete_edge(y,head)
                edges.remove((y, e[1], None))
# Next we delete the edges from the grandparents to the tail.
            gp=self.neighbors_in(tail)
            for y in gp:
                self.delete_edge(y,tail)
                edges.remove((y, e[0], None))
# Now we delete the edges from the head to the grandchildren. 
            gc=self.neighbors_out(head)
            for y in gc:
                self.delete_edge(head,y)
                edges.remove((e[1], y, None))
# Now we either delete or add an edge depending on whether there is an edge between the other children and the other parents.
            for x in par:
                for y in sib:
                    if y in self.neighbors_out(x):
                        self.delete_edge(x,y)
                        edges.remove((x, y, None))
                    else: 
                        self.add_edge(x,y)
                        edges.append((x, y, None))
#FINALLY... we delete the edge between the tail and the head!
            self.delete_vertices([head,tail])
            edges.remove((e[0], e[1], None))
        return self.connected_components_number()
#        print 'Homology is ',self.connected_components_number(),' dimensional'

    def __repr__(self):
        return "Chain complex with "+repr(len(self.vertices()))+" generators."

    def tensor(self, other):
        """
        Returns a graph representing the tensor product of the two complexes self and other
        """
        result = ChainCx()
        for v in self.vertices():
            for w in other.vertices():
                result.add_vertex((v,w))
        for v in result.vertices():
            for w in self.neighbors_out(v[0]):
                result.add_edge( (v,(w,v[1])) )
            for w in other.neighbors_out(v[1]):
                result.add_edge( (v,(v[0],w)) )
        return result
