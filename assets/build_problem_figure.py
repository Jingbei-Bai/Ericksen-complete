"""Draw representatives of the six families from the article's exact maps.

Run: python assets/build_problem_figure.py
Requirements: NumPy, SymPy, Matplotlib. No external images are used.
Each plotted map is also differentiated symbolically to check its signed
Cartesian Jacobian. The examples follow Section 2.2, equations (2.3)--(2.8).
"""
from pathlib import Path
from dataclasses import dataclass
import argparse
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LightSource
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import proj3d
import numpy as np
import sympy as sp

ASSETS = Path(__file__).resolve().parent
TITLE = "A local solution of Ericksen's problem"
SHORT_NAMES = ['Homogeneous','Block bending','Sector straightening',
               'Inflation and torsion','Spherical shell','Logarithmic cylindrical']
for filename in ['lmroman17-regular.otf','lmroman10-regular.otf']:
    font_manager.fontManager.addfont(ASSETS/'fonts'/filename)
TITLE_FONT = font_manager.FontProperties(fname=ASSETS/'fonts/lmroman17-regular.otf')
TEXT_FONT = font_manager.FontProperties(fname=ASSETS/'fonts/lmroman10-regular.otf')
u, v, w = sp.symbols('u v w', real=True)
Q = (u, v, w)
RAT = sp.Rational
INK = '#512D1F'
COPPER = '#BB7044'
GRID = '#7D4228'
IVORY = '#FFF9F0'


@dataclass
class Family:
    number: int
    name: str
    reference: sp.Matrix
    current: sp.Matrix
    bounds: tuple
    samples: tuple
    coordinates: str
    parameters: dict
    density: sp.Expr


def cylindrical(radius, angle, height):
    return sp.Matrix([radius*sp.cos(angle), radius*sp.sin(angle), height])


def families():
    cartesian = sp.Matrix(Q)
    cylinder = cylindrical(u, v, w)
    wedge = ((RAT(7,10), RAT(11,10)), (RAT(1,4), sp.Integer(5)),
             (-RAT(11,20), RAT(11,20)))
    sphere = sp.Matrix([u*sp.sin(w)*sp.cos(v), u*sp.sin(w)*sp.sin(v), u*sp.cos(w)])
    radial_sphere = (u**3+RAT(3,5))**RAT(1,3)
    return [
        Family(0, 'Homogeneous deformation', cartesian,
               sp.Matrix([u+RAT(9,20)*v, v, w]),
               ((-RAT(3,5),RAT(3,5)),(-RAT(3,5),RAT(3,5)),(-RAT(9,20),RAT(9,20))),
               (9,11,8), '(X,Y,Z)', {'L':[['1','9/20','0'],['0','1','0'],['0','0','1']]}, sp.Integer(1)),
        Family(1, 'Bending of a rectangular block', cartesian,
               cylindrical(sp.sqrt(2*u+1),v,w),
               ((sp.Integer(0),RAT(3,5)),(-RAT(11,10),RAT(11,10)),(-RAT(9,20),RAT(9,20))),
               (11,36,9), '(X,Y,Z)', {'alpha':'1','beta':'1','d':'1','gamma':'0'}, sp.Integer(1)),
        Family(2, 'Straightening of a cylindrical sector', cylinder,
               sp.Matrix([(u**2-1)/2,v,w]),
               ((sp.Integer(1),sp.sqrt(RAT(11,5))),(-RAT(11,10),RAT(11,10)),(-RAT(9,20),RAT(9,20))),
               (11,36,9), '(R,Theta,Z)', {'alpha':'1','beta':'1','x0':'-1/2','gamma':'0'},u),
        Family(3, 'Inflation and torsion', cylinder,
               cylindrical(sp.sqrt(RAT(4,5)*u**2+RAT(2,5)),v+RAT(3,5)*w,RAT(5,4)*w),
               wedge, (13,65,17), '(R,Theta,Z)',
               {'A':'4/5','B0':'2/5','k':'1','delta':'3/5','kappa':'0','gamma':'5/4'},u),
        Family(4, 'Spherical shell deformation', sphere,
               sphere*radial_sphere/u,
               ((RAT(16,25),sp.Integer(1)),(sp.Integer(0),sp.pi),(sp.Integer(0),sp.pi)),
               (11,49,41), '(R,Theta,Phi)', {'epsilon':'1','Q':'I','B0':'3/5'},u**2*sp.sin(w)),
        Family(5, 'Logarithmic cylindrical deformation', cylinder,
               cylindrical(RAT(6,5)*u,v+RAT(7,5)*sp.log(u),RAT(25,36)*w),
               ((RAT(13,20),RAT(11,10)),(RAT(1,4),sp.Integer(5)),(-RAT(1,2),RAT(1,2))),
               (13,65,15), '(R,Theta,Z)',
               {'alpha':'6/5','beta':'7/5','k':'1','gamma':'25/36'},u),
    ]


def evaluate(expressions, coordinates):
    """Evaluate the same symbolic expressions used in the Jacobian checks."""
    shape = np.broadcast_shapes(*(np.shape(x) for x in coordinates))
    functions = [sp.lambdify(Q, expression, 'numpy') for expression in expressions]
    return np.stack([np.broadcast_to(f(*coordinates),shape) for f in functions],axis=-1)


def surfaces(family, expressions):
    ranges = [np.linspace(float(a),float(b),n) for (a,b),n in zip(family.bounds,family.samples)]
    patches = []
    for fixed in range(3):
        # The spherical coordinate poles are edges, not two-dimensional faces.
        if family.number==4 and fixed==2:
            continue
        remaining = [j for j in range(3) if j!=fixed]
        s,t = np.meshgrid(ranges[remaining[0]],ranges[remaining[1]],indexing='ij')
        for end in [0,-1]:
            coordinates = [None,None,None]
            coordinates[fixed] = np.full(s.shape,ranges[fixed][end])
            coordinates[remaining[0]] = s
            coordinates[remaining[1]] = t
            patch = evaluate(expressions,coordinates)
            assert np.isfinite(patch).all(), (family.number,fixed,end)
            # Orient every parameter face outwards for surface lighting.
            # The (R,Theta,Phi) spherical chart has negative orientation.
            orientation = (1 if end==-1 else -1)*((-1)**fixed)*(-1 if family.number==4 else 1)
            if orientation==1:
                patch = patch[::-1]
            patches.append(patch)
    return patches


def verify(family):
    reference_jacobian = family.reference.jacobian(Q).det()
    current_jacobian = family.current.jacobian(Q).det()
    determinant = sp.simplify(sp.trigsimp(current_jacobian/reference_jacobian))
    assert determinant==1, (family.number,determinant)
    volume = family.density
    for coordinate,bounds in zip(Q,family.bounds):
        volume = sp.integrate(volume,(coordinate,*bounds))
    return {
        'family':f'F{family.number}', 'paper_label':f'intro:family{family.number}',
        'equation':f'2.{family.number+3}', 'name':family.name,
        'reference_coordinates':family.coordinates,
        'reference_cartesian_map':[sp.sstr(x) for x in family.reference],
        'current_cartesian_map':[sp.sstr(x) for x in family.current],
        'parameters':family.parameters,
        'reference_domain':[[sp.sstr(a),sp.sstr(b)] for a,b in family.bounds],
        'signed_cartesian_jacobian':sp.sstr(determinant),
        'reference_and_current_volume':sp.sstr(sp.simplify(volume)),
    }


def show(ax,patches,centre,span,number):
    ax.set_facecolor(IVORY)
    ax.patch.set_alpha(0)
    light = LightSource(azdeg=135,altdeg=45)
    polygons = []
    for patch in patches:
        points = patch-centre
        polygons.append(np.stack([points[:-1,:-1],points[:-1,1:],
                                  points[1:,1:],points[1:,:-1]],axis=2).reshape(-1,4,3))
    # Sort all boundary facets together, including the inside of each shell.
    collection = Poly3DCollection(np.concatenate(polygons),facecolors=[COPPER],
                                 edgecolors=[GRID],linewidth=.23,antialiased=True,
                                 shade=True,lightsource=light,zsort='average')
    ax.add_collection3d(collection)
    for direction,limit in zip('xyz',span/2):
        getattr(ax,f'set_{direction}lim')(-limit,limit)
    # Matplotlib versions may normalize an array argument in place.
    # Keep the same physical axis scale for the two configurations.
    ax.set_box_aspect(span.copy(),zoom=1.0)
    ax.view_init(elev=19 if number==4 else 24,azim=-75 if number==4 else -57)
    ax.set_proj_type('ortho')
    ax.set_axis_off()


def draw(family,fig,x,y,width,height):
    reference = surfaces(family,family.reference)
    current = surfaces(family,family.current)
    all_points = [np.concatenate([p.reshape(-1,3) for p in group]) for group in [reference,current]]
    centres = [(p.max(axis=0)+p.min(axis=0))/2 for p in all_points]
    span = np.maximum(*[np.ptp(p,axis=0) for p in all_points])*1.16
    fig.add_artist(FancyBboxPatch((x,y),width,height,boxstyle='round,pad=0.002,rounding_size=.006',
                                 transform=fig.transFigure,facecolor=IVORY,
                                 edgecolor='#E6B58B',linewidth=.75,zorder=-5))
    fig.text(x+width/2,y+height*.912,f'$\\mathcal{{F}}_{family.number}$   {SHORT_NAMES[family.number]}',
             ha='center',va='center',fontsize=17,color=INK,fontproperties=TEXT_FONT)
    axis_width = width*.47
    axis_y = y+height*.085
    axis_height = height*.745
    views = []
    for i,(patches,centre) in enumerate(zip([reference,current],centres)):
        axis_x = x+width*(.015 if i==0 else .515)
        ax = fig.add_axes([axis_x,axis_y,axis_width,axis_height],projection='3d')
        show(ax,patches,centre,span,family.number)
        views.append((ax,all_points[i]-centre,span.copy(),family.number))
    fig.add_artist(FancyArrowPatch((x+width*.478,y+height*.485),
                                  (x+width*.522,y+height*.485),
                                  transform=fig.transFigure,arrowstyle='-|>',
                                  mutation_scale=17,linewidth=1.15,color='#945331'))
    return views


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=ASSETS)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    cases = families()
    checks = [verify(case) for case in cases]
    assert [c['signed_cartesian_jacobian'] for c in checks]==['1']*6
    # Check the two geometric distinctions displayed in the corresponding panels.
    assert sp.simplify(cases[1].current[0]**2+cases[1].current[1]**2-(2*u+1))==0
    assert sp.simplify(cases[3].current[0]**2+cases[3].current[1]**2-(RAT(4,5)*u**2+RAT(2,5)))==0
    assert all(sp.simplify(sp.diff(component,w))==0 for component in cases[5].current[:2])
    matplotlib.rcParams.update({'font.family':TEXT_FONT.get_name(),'mathtext.fontset':'cm',
                               'pdf.fonttype':3,'savefig.facecolor':IVORY})
    height_scale = 1.0
    fig = plt.figure(figsize=(16,10.7),facecolor=IVORY)
    fig.text(.5,.953,TITLE,ha='center',va='center',fontsize=39,color=INK,fontproperties=TITLE_FONT)
    statement = 'Determine all deformations admitting equilibrium for every isotropic hyperelastic law.'
    fig.text(.5,.894,statement,ha='center',va='center',fontsize=20,color=INK,fontproperties=TEXT_FONT)
    scope = 'Homogeneous incompressible solids · No body forces'
    fig.text(.5,.858,scope,ha='center',va='center',fontsize=14,color='#845E48',fontproperties=TEXT_FONT)
    views = []
    for family in cases:
        column = family.number%3
        row = family.number//3
        views.extend(draw(family,fig,.021+column*.326,(.52-row*.342)*height_scale,
                          .305,.319*height_scale))
    fig.add_artist(FancyBboxPatch((.022,.018),.956,.126,
                                 boxstyle='round,pad=0.004,rounding_size=.012',
                                 transform=fig.transFigure,facecolor='#FFF0DC',
                                 edgecolor='none',zorder=-5))
    classification = 'Complete local classification: the six classical families.'
    fig.text(.5,.117,classification,ha='center',va='center',fontsize=24,color=INK,fontproperties=TITLE_FONT)
    explanation = [
        'Equilibrium for every isotropic hyperelastic law constrains principal direction fields in three dimensions.',
        'For distinct constant principal stretches, one principal direction is normal to a family of surfaces.',
    ]
    for line,position in zip(explanation,[.076,.046]):
        fig.text(.5,position,line,ha='center',va='center',fontsize=16,color=INK,fontproperties=TEXT_FONT)
    fig.canvas.draw()
    fitted_zoom = {}
    for start in range(0,len(views),2):
        pair = views[start:start+2]
        permitted = []
        for ax,points,span,number in pair:
            xp,yp,_ = proj3d.proj_transform(*points.T,ax.get_proj())
            display = ax.transData.transform(np.column_stack([xp,yp]))
            origin = ax.transData.transform([[0,0]])[0]
            offsets = display-origin
            frame = ax.get_window_extent()
            permitted.extend([
                (origin[0]-frame.x0-6)/(-offsets[:,0].min()),
                (frame.x1-origin[0]-6)/offsets[:,0].max(),
                (origin[1]-frame.y0-6)/(-offsets[:,1].min()),
                (frame.y1-origin[1]-6)/offsets[:,1].max(),
            ])
        zoom = min(1.30,.97*min(permitted))
        fitted_zoom[f'F{number}'] = zoom
        for ax,points,span,number in pair:
            ax.set_box_aspect(span.copy(),zoom=zoom)
    fig.canvas.draw()
    for ax,points,span,number in views:
        xp,yp,_ = proj3d.proj_transform(*points.T,ax.get_proj())
        display = ax.transData.transform(np.column_stack([xp,yp]))
        frame = ax.get_window_extent()
        assert display[:,0].min()>frame.x0+2 and display[:,0].max()<frame.x1-2
        assert display[:,1].min()>frame.y0+2 and display[:,1].max()<frame.y1-2
    fig.savefig(args.output_dir/'ericksen_problem.png',dpi=190)
    fig.savefig(args.output_dir/'ericksen_problem.pdf',metadata={
        'Title':TITLE,'Author':'Jingbei Bai'})
    plt.close(fig)
    report = {'status':'PASS','source':ASSETS.name+'/'+Path(__file__).name,
              'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'paper_section':'2.2','geometry':'images of the same parameter-domain faces under the exact Cartesian maps',
              'paper_title':TITLE,'problem_statement':statement,
              'physical_scope':scope,
              'classification_statement':classification,'principal_direction_explanation':explanation,
              'font':'Latin Modern Roman 2.005',
              'versions':{'numpy':np.__version__,'sympy':sp.__version__,'matplotlib':matplotlib.__version__},
              'families':checks,
              'fitted_pair_zooms':fitted_zoom,
              'cross_checks':{'F1_and_F2_inverse_coordinate_examples':True,
                              'F3_radial_coordinate_independent_of_height':True,
                              'F5_planar_coordinates_independent_of_height':True,
                              'all_projected_mesh_vertices_inside_frames':True,
                              'same_axis_scale_camera_and_zoom_within_each_pair':True}}
    # F2 in cylindrical coordinates is the inverse of the chosen F1 map.
    inverse = cases[2].current.subs(u,sp.sqrt(2*u+1),simultaneous=True)
    assert all(sp.simplify(a-b)==0 for a,b in zip(inverse,sp.Matrix(Q)))
    assert cases[1].parameters['alpha']==cases[2].parameters['alpha']=='1'
    (args.output_dir/'geometry_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','families':6,'signed_cartesian_jacobians':['1']*6,
                      'outputs':[str(args.output_dir/name) for name in
                                 ['ericksen_problem.png','ericksen_problem.pdf','geometry_checks.json']]},indent=2))


if __name__=='__main__':
    main()
