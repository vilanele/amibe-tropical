import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox, Button
from sage.misc.parser import Parser, Tokenizer
from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
from matplotlib.patches import Circle, Polygon, FancyArrowPatch, ArrowStyle
from matplotlib.lines import Line2D
import random as rd


AMIBE = 0
SUBDIV = 1
TEXTBOX = 2

def my_set_val(textbox, val):
    newval = str(val)
    if textbox.text == newval:
        return
    textbox.text = newval
    textbox.text_disp.remove()
    textbox.text_disp = textbox._make_text_disp(textbox.text)
    textbox._rendercursor()
    textbox._notify_change_observers()
    textbox._notify_submit_observers()


def parse_tropw( trop, w ):

    p = Parser()
    vertices = [ list(x) for x in p.parse_sequence(trop) ]
    w = p.parse_sequence(w)

    return vertices, w

def proj_gon3d( lifted_polytope, i):

    E = lifted_polytope.inequalities_list()[i]


    G = []
    for v in lifted_polytope.vertices():
        if vector(E[1:]).dot_product(vector(v)) == -E[0]:
            G.append( list(v[:2]) )


    return G


def order_gon2d( vertices ):
    sorted_vertices = LatticePolytope_PPL(vertices).ordered_vertices()
    return sorted_vertices


def intersect_ray_z( ray ):
    return [(1/ray[2])*x for x in ray[:2] ]

def intersect_plan_z( x, y ):

    normal = vector(x).cross_product(vector(y))
    A = Matrix(QQ,[ normal, [0,0,1]])
    direction = A.transpose().kernel().basis()[0]

    C = Cone([x,y])
    if not C.contains(direction):
        direction = -1*direction

    return direction[:2]

def find_edge( lifted_polytope, i, j):

    L = []
    A = lifted_polytope.inequalities_list()[i]
    B = lifted_polytope.inequalities_list()[j]

    for v in lifted_polytope.vertices_list():
        if (vector(A[1:]).dot_product(vector(v)) == -A[0]) and (vector(B[1:]).dot_product(vector(v)) == -B[0]):
            L.append(v)

    return L[0],L[1]


def integer_points( A, B):

    p = Polyhedron([A,B])
    return p.integral_points()

def integer_first_vector(vector):

    m = lcm(QQ(vector[0]).denominator(),QQ(vector[1]).denominator())
    g = gcd( vector[0]*m, vector[1] * m)

    return [ (vector[0]* m)/g, (vector[1] *m)/g]

def distance(A,B):
    return RR(norm(vector([A[0]-B[0],A[1]-B[1]]))).nearby_rational(max_denominator=100)

def distance_exact(A,B):
    return RR(norm(vector([A[0]-B[0],A[1]-B[1]])))

def get_min_distance(V):

    distances = []

    for (i,v) in enumerate(V[:-1]):

        W = V[i+1:]
        D = [ distance_exact(v,x) for x in W]
        distances.append(min(D))

    return min(distances)

def make_limits(V, f):

    if len(V) == 1:
        t = 10
        center=V[0]
        return [ center[0] -t, center[0] + t], [center[1] - t, center[1] + t]

    center = Polyhedron( vertices=V).center()
    size = max([distance(center,v) for v in V])
    total = (size + f*size)*2

    return [center[0] - total/2, center[0] + total/2], [center[1] - total/2, center[1] + total/2]

class DualityID:

    def __init__(self,id,color=None):
        self.id = id
        self.color = color

    def set_id(self,id):
        self.id = id

    def get_id(self):
        return id


class Point2D(DualityID):

    def __init__(self,coord,ax=None,id=None):
        DualityID.__init__(self,id)
        self.x = coord[0]
        self.y = coord[1]
        self.coord = coord
        self.ax = ax

    def plot(self):
        self.circle = Circle( self.coord, self.radius,alpha=0.75,color=self.color,fill=True,linewidth=None,picker=True)
        #self.ax.plot(self.x,self.y,marker='o',color=self.color, markersize=12)
        self.ar = self.ax.add_patch(self.circle)

    def set_radius(self,radius):
        self.radius= radius



class Segment2D(DualityID):

    ray = False

    def __init__(self,P1,P2,ax=None,id=None):
        DualityID.__init__(self,id)
        self.P1 = P1
        self.P2 = P2
        self.X = [ P1.x,P2.x ]
        self.Y = [ P1.y,P2.y ]
        self.ax = ax

    def plot(self):
        line = Line2D(self.X,self.Y,linewidth=2,color='black',antialiased=True)
        self.ar = self.ax.add_line(line)


class Segment2D_int(Segment2D):

    def __init__(self,P1,P2,ax=None,id=None):
        DualityID.__init__(self,id)
        Segment2D.__init__( self,P1,P2,ax,id)
        self.integers = integer_points(P1.coord,P2.coord)
        self.integer_length = len( self.integers) - 1

    def plot(self):


        for P in self.integers:
            self.ax.plot(P[0],P[1],marker='o', markeredgewidth=2, markeredgecolor='black',fillstyle='full',markersize='6',color='white',zorder=2)

        self.ar = self.ax.plot(self.X,self.Y,linewidth=2,color='black',zorder=1)



class Ray2D(DualityID):

    ray = True

    def __init__(self,base,dir,ax=None,id=None):
        self.ax=ax
        self.id=id
        self.base=base
        self.dir=dir

    def plot(self):
        t = 30
        line = Line2D([self.base.x,self.base.x + t*self.dir.x], [self.base.y,self.base.y + t*self.dir.y], color='black', linewidth=2,antialiased=True)
        self.ar = self.ax.add_line(line)

class Gon2D(DualityID):

    def __init__(self,vertices,ax=None,id=None):
        DualityID.__init__(self,id)
        self.vertices = list(vertices)
        self.ax = ax


    def plot(self):

         pgon = Polygon( self.vertices, fill=True,color=self.color,alpha=0.75,zorder=0)
         self.ar = self.ax.add_patch(pgon)

class Amoeba:

    def __init__(self):
        self.dic_indices_id = dict()
        self.points = set()
        self.segments = set()
        self.rays = set()
        self.dic = dict()
        self.dic_id_balance = dict()

def find_lower_faces( lifted_polytope ):

    lower_faces = set()

    for (i,ray) in enumerate(lifted_polytope.inequalities_list()):

        if ray[3] > 0:
            lower_faces.add(i)

    return lower_faces

def find_segments_rays(lifted_polytope, all_facets, lower_facets):

    segments = set()
    rays = set()

    m = lifted_polytope.facet_adjacency_matrix()

    tested = set()
    for i in lower_facets:
        for j in all_facets.difference(tested):

            if m[i][j]:

                if i in lower_facets and j in lower_facets:
                    segments.add( (i,j) )
                else:
                    if i in lower_facets:
                        rays.add( (i,j) )
                    else:
                        rays.add( (j,i) )

            tested.add(i)

    return segments, rays



class duality_manager:

    def __init__(self,N, vertices, w):
        self.vertices = vertices
        self.w = w
        self.ids = set(range(1,N+1))
        self.dic = dict()
        self.colormap = dict()
        self.balance_view = False
        self.over=False
        self.over_id=None
        self.artists = None


    def add(self, duality_block ):

        id = self.ids.pop()
        self.dic[ id ] = duality_block

        for e in duality_block.values():
            e.set_id( id )

        return id

    def get(self, id, section):
        return self.dic[id][section]

    def get_block( self, id):
        return self.dic[id]

    def set_amibe(self, amibe):
        self.amibe = amibe

    def set_figure(self, figure):
        self.figure = figure

    def make_color_map(self):

        self.color_map = dict()
        for (id,block) in self.dic.items():

            c = ( rd.random(), rd.random(), rd.random() )
            self.color_map[id] = c
            for e in block.values():
                e.color = c

    def build_axes(self):
        self.figure = plt.figure()
        self.figure.set_size_inches( 10,10 )
        cid2 = self.figure.canvas.mpl_connect('motion_notify_event', self.onmousemove)

        ax_amibe = self.figure.add_axes([0.025,0.35,0.3,0.6] )
        ax_balance = self.figure.add_axes([0.35, 0.35, 0.3,0.6])
        ax_subdiv = self.figure.add_axes([ 0.675, 0.35, 0.3, 0.6])
        ax_textbox = self.figure.add_axes([ 0.25,0.1,0.6,0.05] )

        ax_amibe.set_aspect(aspect=1.0)
        ax_subdiv.set_aspect(aspect=1.0)
        ax_balance.set_aspect(aspect=1.0)

        TextBox( ax_textbox, "polynome")

    def set_axis(self, ax_amibe, ax_subdiv, ax_balance ):
        self.axis = [ ax_amibe, ax_subdiv, ax_balance ]


    def start(self):

        self.make_color_map()
        V = [ self.get(id,AMIBE).coord for id in self.amibe.points]
        if len(V) > 1:

            amibe_xlim, amibe_ylim  = make_limits(V,0.6)



            r = (amibe_xlim[1] - amibe_xlim[0])/20
            min_radius = get_min_distance( V )/2.2
            radius = min( r, min_radius)
        else:
            amibe_xlim, amibe_ylim = [-5,5],[-5,5]
            #subdiv_xlim, subdiv_ylim = [-5,5],[-5,5]
            radius = 0.5

        subdiv_xlim, subdiv_ylim = make_limits(self.vertices, 0.1)

        for id in self.amibe.points:
            self.get(id,AMIBE).set_radius( radius)

        self.radius = radius

        ax_amibe = self.axis[0]
        ax_balance = self.axis[2]
        ax_subdiv = self.axis[1]

        ax_amibe.clear()
        ax_balance.clear()
        ax_subdiv.clear()

        ax_amibe.set_xlim(amibe_xlim[0],amibe_xlim[1])
        ax_amibe.set_ylim(amibe_ylim[0],amibe_ylim[1])
        ax_subdiv.set_xlim(subdiv_xlim[0],subdiv_xlim[1])
        ax_subdiv.set_ylim(subdiv_ylim[0],subdiv_ylim[1])
        self.amibe_limits = [amibe_xlim, amibe_ylim]
        ax_balance.grid(True,color='black',alpha=0.25,linewidth=1)

        self.plot()
        #plt.show()


    def plot(self):



        for block in self.dic.values():

            e = block[AMIBE]
            e.ax = self.axis[AMIBE]
            e.plot()


            f = block[SUBDIV]
            f.ax = self.axis[SUBDIV]
            f.plot()

        self.figure.canvas.draw()

    def plot_balance(self, id_found):

        things = self.amibe.dic_id_balance[id_found]
        xlim = things[0]
        ylim = things[1]
        circle = things[2]
        lines = things[3]
        arrows = things[4]

        ax = self.axis[2]
        ax.clear()
        ax.set_xlim(xlim[0],xlim[1])
        ax.set_ylim(ylim[0],ylim[1])

        if xlim[1] - xlim[0] < 35 and ylim[1] - ylim[0] < 35:
            ax.grid(True,color='black',alpha=0.15,linewidth=1)
            ax.set_xticks( range(RR(xlim[0]).round(),RR(xlim[1]).round()+1))
            ax.set_yticks( range(round(ylim[0]),round(ylim[1])+1))
        else:
            ax.set_yticklabels([])
            ax.set_xticklabels([])

        circle.set_color(self.color_map[id_found])
        ax.add_patch(circle)
        for l in lines:
            ax.add_line(l)
        for a in arrows:
            ax.add_patch(a)


    def onclick(self,event):


        ax = event.inaxes

        if not self.balance_view:

            # On test si on entre dans la balance balance_view
            if ax == self.axis[AMIBE]:

                Found=False
                for id in self.amibe.points:

                    P = self.get(id,AMIBE)
                    if norm(vector([P.x - event.xdata,P.y -event.ydata])) < self.radius:
                        self.balance_view = True

                        break


        else:

            if ax == self.axis[0]:
                self.balance_view = False


    def onmousemove(self,event):

        if not self.balance_view:
            ax = event.inaxes
            if ax == self.axis[AMIBE]:

                Found=False
                for id in self.amibe.points:

                    P = self.get(id,AMIBE)
                    if norm(vector([P.x - event.xdata,P.y -event.ydata])) < self.radius:
                        Found=True
                        id_found = id
                        self.over_id=id
                        break

                if  self.over_id:
                    P = self.get(self.over_id,AMIBE)
                if Found:

                    if not self.over:
                        self.over = True
                        P.ar.set_alpha(1.0)
                        self.get(P.id, SUBDIV).ar.set_alpha(1.0)

                        for id in self.amibe.points:
                            if id != id_found:
                                self.dic[id][AMIBE].ar.set_alpha(0.075)
                                self.dic[id][SUBDIV].ar.set_alpha(0.075)

                        self.plot_balance(id_found)
                        self.figure.canvas.draw()

                else:

                    if self.over:
                        self.over=False
                        P.ar.set_alpha(0.75)
                        self.get(P.id, SUBDIV).ar.set_alpha(0.75)

                        for id in self.amibe.points:
                            if id != P.id:
                                self.dic[id][AMIBE].ar.set_alpha(0.75)
                                self.dic[id][SUBDIV].ar.set_alpha(0.75)

                        self.axis[2].clear()

                        self.figure.canvas.draw()


def build( vertices,w ):

    lifted_polytope = Polyhedron( vertices=[ v + [i] for (v,i) in zip(vertices,w)] )

    all_facets = set( range(0,lifted_polytope.n_facets()) )
    lower_facets = find_lower_faces( lifted_polytope )

    points = lower_facets
    segments, rays = find_segments_rays( lifted_polytope, all_facets, lower_facets )

    manager = duality_manager(40,vertices, w )
    amibe = Amoeba()

    id = []
    for i in points:

        duality_block = dict()

        x = lifted_polytope.inequalities_list()[i][1:]
        duality_block[AMIBE] = Point2D( intersect_ray_z(x) )

        gon2d_vertices = order_gon2d( proj_gon3d( lifted_polytope, i) )
        duality_block[SUBDIV] = Gon2D( gon2d_vertices  )

        id.append( manager.add( duality_block ) )

        amibe.points.add( id[-1] )
        amibe.dic[ id[-1] ] = set()
        amibe.dic_indices_id[i] = id[-1]

    for c in segments:

        duality_block = dict()

        i = c[0]
        j = c[1]

        # Amibe
        P1 = lifted_polytope.inequalities_list()[i][1:]
        P2 = lifted_polytope.inequalities_list()[j][1:]
        duality_block[AMIBE] = Segment2D( Point2D( intersect_ray_z(P1)), Point2D(intersect_ray_z(P2)) )

        # Subdivision
        A, B = find_edge(lifted_polytope, i,j)
        duality_block[SUBDIV] = Segment2D_int( Point2D(A[:2]),Point2D(B[:2]) )

        # add to manager
        id.append(manager.add( duality_block))

        #build amibe structure
        amibe.segments.add(id[-1])
        amibe.dic[ amibe.dic_indices_id[i] ].add( id[-1] )
        amibe.dic[ amibe.dic_indices_id[j] ].add( id[-1])


    for r in rays:

        duality_block = dict()

        # amibe
        L = lifted_polytope.inequalities_list()
        base = intersect_ray_z( L[r[0]][1:] )
        direction = intersect_plan_z( L[r[0]][1:4], L[r[1]][1:4] )
        duality_block[AMIBE] = Ray2D(Point2D(base),Point2D(direction))

        # subdivision
        P1, P2 = find_edge(lifted_polytope, r[0], r[1])
        duality_block[SUBDIV] = Segment2D_int(Point2D(P1[:2]), Point2D(P2[:2]))

        # add to manager
        id.append( manager.add( duality_block ))

        #build amibe structure
        amibe.rays.add(id[-1])
        amibe.dic[ amibe.dic_indices_id[r[0]]].add( id[-1] )


    for (id,indices) in amibe.dic.items():

        circle = Circle([0,0], 0.2,fill=True,linewidth=None,zorder=2)
        point = manager.get(id,AMIBE)

        lines = []
        arrows = []
        L = []
        for i in indices:
            P = manager.get(i,AMIBE)
            if P.ray:
                dir = P.dir.coord
                #line = Line2D([0,200*P.dir.x],[0,200*P.dir.y],color='black',linewidth=0.5,zorder=1,alpha=0.5)
            else:
                if P.P1.x == point.x and P.P1.y == point.y:
                    other = P.P2
                    base = P.P1
                else:
                    other = P.P1
                    base = P.P2

                dir = [ other.x - base.x, other.y - base.y]

            line = Line2D([0,200*dir[0]],[0,200*dir[1]],color='black',linewidth=0.5,zorder=1,alpha=0.5)
            lines.append(line)

            v = integer_first_vector( dir )
            weight = manager.get(i,SUBDIV).integer_length


            style1= ArrowStyle("Fancy", head_length=6, head_width=6, tail_width=1)
            style2=ArrowStyle("Fancy", head_length=6, head_width=3, tail_width=1)
            arrow1 = FancyArrowPatch((0,0),(v[0],v[1]),color='black',alpha=1.0,arrowstyle=style1,linewidth=3,mutation_scale=1.5)
            arrow2 = FancyArrowPatch((0,0),(weight*v[0],weight*v[1]),color='black',alpha=0.75,linewidth=1.5,arrowstyle=style2,mutation_scale=1.5)
            L.append( [weight*v[0],weight*v[1]])
            arrows.append(arrow1)
            arrows.append(arrow2)

        xlim,ylim = make_limits(L,0.2)

        amibe.dic_id_balance[id] = [xlim,ylim,circle,lines,arrows]

    manager.set_amibe(amibe)

    return manager

# It is simpler to treat the case of equal weights (i.e the amibe is just the normal fan) apart
def fan( vertices,w ):

    manager = duality_manager(40,vertices,w )
    amibe = Amoeba()


    duality_block = dict()

    duality_block[AMIBE] = Point2D([0,0])
    duality_block[SUBDIV] = Gon2D( order_gon2d(vertices) )

    id_point = manager.add( duality_block )

    amibe.points.add( id_point)
    amibe.dic_id_balance[id_point] = set()
    poly = Polyhedron(vertices)
    m = poly.incidence_matrix()
    lines = []
    arrows = []
    L = []
    circle = Circle([0,0], 0.2,fill=True,linewidth=None,zorder=2)
    for (i,r) in enumerate( poly.normal_fan().rays() ):

        V = []
        for j in range(0,poly.n_vertices()):
            if m[j][i]:
                V.append( poly.vertices_list()[j] )

        A = V[0]
        B = V[1]

        duality_block=dict()
        duality_block[AMIBE] = Ray2D(Point2D([0,0]), Point2D(r) )
        duality_block[SUBDIV] = Segment2D_int( Point2D(A), Point2D(B) )

        id_seg = manager.add(duality_block )

        line = Line2D([0,200*r[0]],[0,200*r[1]],color='black',linewidth=0.5,zorder=1,alpha=0.5)


        v = integer_first_vector( r)
        weight = manager.get(id_seg,SUBDIV).integer_length
        L.append( [weight*v[0],weight*v[1]])
        style1= ArrowStyle("Fancy", head_length=6, head_width=6, tail_width=1)
        style2=ArrowStyle("Fancy", head_length=6, head_width=3, tail_width=1)
        arrow1 = FancyArrowPatch((0,0),(v[0],v[1]),color='black',alpha=1.0,arrowstyle=style1,linewidth=3,mutation_scale=1.5)
        arrow2 = FancyArrowPatch((0,0),(weight*v[0],weight*v[1]),color='black',alpha=0.75,linewidth=1.5,arrowstyle=style2,mutation_scale=1.5)

        lines.append( line )
        arrows.append( arrow1)
        arrows.append( arrow2 )

    xlim,ylim = make_limits(L,0.2)
    amibe.dic_id_balance[id_point] = [ xlim,ylim,circle,lines,arrows ]

    manager.set_amibe(amibe)
    return manager

class tropicalview_app:

    def __init__(self):

        figure = plt.figure()
        figure.set_size_inches( 11,11 )

        ax_amibe = figure.add_axes([0.025,0.35,0.3,0.6] )
        ax_balance = figure.add_axes([0.35, 0.35, 0.3,0.6])
        ax_subdiv = figure.add_axes([ 0.675, 0.35, 0.3, 0.6])

        ax_trop = figure.add_axes([ 0.1, 0.18  ,0.6,0.05])
        ax_w = figure.add_axes([ 0.7,   0.18 ,0.2,0.05])
        ax_button_tropw = figure.add_axes( [0.95, 0.18, 0.04,0.04])
        self.secondary_axis = [  ax_trop, ax_w]

        ax_amibe.set_aspect(aspect=1.0)
        ax_subdiv.set_aspect(aspect=1.0)
        ax_balance.set_aspect(aspect=1.0)

        trop_initial_text = "[4,0],[6,5],[3,4],[3,0],[0,2],[1,1],[1,0],[0,1],[0,0]"
        w_initial_text = " 0,1,0,-1,0,0,0,1,0"


        textbox_trop = TextBox( ax_trop, "Trop =", initial=trop_initial_text)
        textbox_w = TextBox( ax_w, "Poids =", initial=w_initial_text)
        button_tropw = Button( ax_button_tropw, "Ok" )

        self.refs = [ textbox_trop, textbox_w, button_tropw]

        button_tropw.on_clicked( self.onclick_tropw )

        self.figure = figure
        self.axis = [ ax_amibe, ax_subdiv, ax_balance ]


        if len(sys.argv) == 1:

            vertices = [[8,0],[0,2],[1,1],[0,1],[1,0],[0,0],[3,0],[3,4],[6,5]]
            weights = [0,1,0,-1,0,0,0,1,0]

        else:

            ax_button_next = figure.add_axes( [0.52, 0.05, 0.08, 0.04] )
            ax_button_previous = figure.add_axes( [0.4, 0.05, 0.08, 0.04] )
            self.secondary_axis += [ ax_button_next, ax_button_previous]
            button_next = Button( ax_button_next, "Next" )
            button_previous = Button( ax_button_previous, "Previous" )
            self.refs += [ button_next, button_previous]
            button_next.on_clicked( self.onclick_next )
            button_previous.on_clicked( self.onclick_previous )

            filename = sys.argv[1]
            file = open(filename, 'r')

            s = []

            for line in file:
                if len(line):
                    if line[0] != "#" and line.find(':') != -1:
                        i = line.find(':')
                        s.append( parse_tropw( line[:-1][:i], line[:-1][i+1:]) )

            self.polynomials = s
            self.current_polynomial_indice = 0
            self.polynomials_len = len(s)
            vertices = self.polynomials[0][0]
            weights = self.polynomials[0][1]
            self.set_values( vertices, weights)

        self.launch_manager(vertices, weights)

    def run(self):

        plt.show()

    def launch_manager(self, vertices, w ):

        w_set = set(w)
        if len(w_set) == 1:
            self.manager = fan(vertices,w)
        else:
            self.manager = build(vertices, w)

        cid = self.figure.canvas.mpl_connect( 'motion_notify_event', self.manager.onmousemove)
        cid2 = self.figure.canvas.mpl_connect('button_press_event', self.manager.onclick)
        self.manager.set_axis( self.axis[ 0], self.axis[1], self.axis[2] )
        self.manager.set_figure( self.figure )
        self.manager.start()


    def onclick_tropw(self,event):

        p = Parser()
        vertices = [ list(x) for x in p.parse_sequence(self.refs[0].text) ]
        w = p.parse_sequence(self.refs[1].text)

        self.launch_manager( vertices,w )


    def onclick_next(self, event):



        if self.current_polynomial_indice < self.polynomials_len - 1:
            self.current_polynomial_indice += 1
            vertices, w = self.polynomials[ self.current_polynomial_indice ]

            self.set_values( vertices, w)

            self.launch_manager( vertices, w)

    def onclick_previous(self, event):



        if self.current_polynomial_indice > 0:
            self.current_polynomial_indice -= 1
            vertices, w = self.polynomials[ self.current_polynomial_indice ]

            self.set_values( vertices, w)

            self.launch_manager( vertices, w)

    def set_values( self, vertices, w):

        my_set_val( self.refs[0], str(vertices)[1:-1] )
        my_set_val( self.refs[1], str(w)[1:-1] )



#####
#
#
#####
#
#
#####

app = tropicalview_app()
app.run()
