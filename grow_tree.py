# Copyright (C) 2024 Warren Usui, MIT License
"""
Create tree of pieces
"""
from copy import deepcopy
from itertools import chain

def coord_to_id(coord):
    """
    Convert [y, x] coordinates to a unique id number
    """
    def idc_inner(dist):
        return dist * (dist - 1) + dist + coord[0]

    return idc_inner(abs(coord[0]) + abs(coord[1]))

def id_to_dict(numb):
    """
    Convert id number to an intermediate dictionary

    dict values:
        'dfo' -- distance from origin
        'sv' -- y coordinate starting value
    """
    def lval(id_dict):
        return id_dict['sv'] > 0
    def gdict(rangev):
        return {'dfo': rangev + 1, 'sv': numb - rangev * (rangev + 1)}
    return list(filter(lval, map(gdict, range(0, 4))))[-1]

def id_dconv(idict):
    """
    Convert intermediate dictionary (defined in id_to_dict) to
    [y, x] coordinate values.
    """
    def idd_inner(yval):
        return [yval, idict['dfo'] - abs(yval)]
    return idd_inner(idict['sv'] - idict['dfo'])

def id_to_coord(numb):
    """
    Convert id number to [x, y] coordinates
    """
    return id_dconv(id_to_dict(numb))

def max_indx(indx):
    """
    Highest id for the next indx
    """
    return (indx + 1) * (indx + 2)

def get_val(node_v):
    """
    Convert a set of coordinates that represent a figure to one
    unique numerical value
    """
    def exponent2(val):
        return 2 ** val
    return sum(map(exponent2, map(coord_to_id, node_v)))

def comp_val(node_v):
    """
    Wrap a node's info into a set of points to get a unique numerical value
    """
    return get_val(node_v['lineage'] + [node_v['coord']])

def get_off_coords(p_node):
    """
    p_node is a list consisting of a level number and a node.
    Generate a list of coordinates within the range of possible
    open squares.  This list is in reverse order so that points on
    the y-axis get checked first
    """
    def goff_get_sq(sq_id):
        return id_to_coord(sq_id)
    return map(goff_get_sq, range(
                max_indx(len(p_node[1]['lineage'])), 0, -1))

def find_n(p_node):
    """
    Check if a point is a neighbor to p_node.  Return that point
    as a new node if it is.
    """
    def find_n_inner(off_coord):
        def ntest(tcoord):
            return abs(tcoord[0] - off_coord[0]) + \
                    abs(tcoord[1] - off_coord[1]) == 1
        if off_coord == p_node[1]['coord'] or off_coord in \
                p_node[1]['lineage']:
            return []
        return list(filter(ntest, [p_node[1]['coord']] +
                           p_node[1]['lineage']))
    return find_n_inner

def get_off(p_node):
    """
    Get list of neighboring points for a given node
    """
    return filter(find_n(p_node), get_off_coords(p_node))

def glin(n_lin_data):
    """
    Combine lineage and coord points into one list
    """
    return n_lin_data['lineage'] + [n_lin_data['coord']]

def next_level(prev_nodes):
    """
    Add a new level onto the tree.  Prev_nodes is ia list of
    previous levels containing the nodes at that level
    """
    def nl_wnp(nxt_row):
        def hndl_grp(this_grp):
            def set_grp(new_node):
                return {'coord': new_node,
                        'lineage': glin(prev_nodes[-1][this_grp[0]]),
                        'pindx': this_grp[0]}
            return map(set_grp, this_grp[1])
        return map(hndl_grp, nxt_row)
    return list(chain.from_iterable(nl_wnp(list(
                                enumerate(list(map(get_off,
                                enumerate(prev_nodes[-1]))))))))

def gen_tree(prev_nodes):
    """
    Recursive routine that adds new layer to tree and terminates when
    the leaves represent pentominos
    """
    if len(prev_nodes) >= 5:
        return prev_nodes
    return gen_tree(prev_nodes + [rm_dups(next_level(prev_nodes))])

def rm_dups(tree_lev):
    """
    Scan a tree level and remove duplicate figures
    """
    def zindx(indx):
        def zinner(n_data):
            return n_data[indx]
        return zinner
    def rd_proc_zip_data(zip_data):
        return list(map(zindx(0), filter(zindx(1), zip_data)))
    def rd_tf_dup_list(tf_list):
        return zip(tree_lev, tf_list)
    def rd_wvals(nvals):
        def lchk(nv_ind):
            return nvals[nv_ind] not in nvals[0:nv_ind]
        return rd_tf_dup_list(list(map(lchk, range(0, len(nvals)))))
    return rd_proc_zip_data(rd_wvals(list(map(comp_val, tree_lev))))

def add_olinks(rtree):
    """
    Add the offspring links to the tree
    """
    def ol_inn(row_num):
        def add_link(link_data):
            def add_link_inner(inode):
                return dict(inode[1]) | {'offspring': link_data[inode[0]]}
            return add_link_inner
        def get_oindx(node):
            return node['pindx']
        def set_oindx(node):
            return list(map(get_oindx, rtree[row_num + 1])).index(node)
        def gp_index():
            return rtree[row_num + 1][-1]['pindx']
        def olinks_set():
            return list(map(set_oindx, range(0, gp_index() + 1)))
        return list(map(add_link(olinks_set()), enumerate(rtree[row_num])))
    return list(map(ol_inn, range(0, 4))) + [rtree[-1][:]]

def gen_rots(points):
    """
    Compute rotations of a figure to find mappings for the same pentomino.
    """
    def get_dim(dimv):
        def gd_inner(idata):
            return idata[dimv]
        return gd_inner
    def flipf():
        def do_flip(indvp):
            return [indvp[1], indvp[0]]
        return list(map(do_flip, points))
    def gr_inner(rot_fctrs):
        def shift_up(spoints):
            def su_inner(vlim):
                def su_pt(pts):
                    return [pts[0], pts[1] - vlim]
                return list(map(su_pt, spoints))
            return su_inner(min(list(map(get_dim(1), spoints))))
        def shift_right(spoints):
            def sr_inner(hlim):
                def sr_pt(pts):
                    return [pts[0] - hlim, pts[1]]
                return list(map(sr_pt, spoints))
            return sr_inner(min(list(map(get_dim(0), spoints))))
        def cmp_fig(ipoints):
            def rotato(ipoint):
                return [ipoint[0] * rot_fctrs[0], ipoint[1] * rot_fctrs[1]]
            return get_val(shift_up(shift_right(list(map(rotato, ipoints)))))
        return [cmp_fig(points),
                cmp_fig(flipf())]
    return chain.from_iterable(list(map(gr_inner,
                    [[1, 1], [1, -1], [-1, 1], [-1, -1]])))

def get_figure(figure):
    """
    Calculate the figure value for a leaf node.  The figure value
    is unique for all rotations of a figure.
    """
    def gff(fignum):
        return {55: 'P', 87: 'V', 118: 'W', 535: 'L',
                563: 'Y', 566: 'N', 1047: 'U', 1125: 'T',
                1126: 'F', 1140: 'Z', 3110: 'X', 66067: 'I'}[fignum]
    return gff(min(gen_rots(figure['lineage'] + [figure['coord']])))

def get_fnumbs(figures):
    """
    Generate new nodes for the leaves of the tree that have figure
    numbers added
    """
    def gf_inner(figure):
        return {'coord': figure['coord'],
                'lineage': figure['lineage'],
                'pindx': figure['pindx'],
                'figure': get_figure(figure)}
    return list(map(gf_inner, figures))

def add_pent_numbs(tree):
    """
    Add figure numbers onto the leaves of the tree
    """
    return deepcopy(tree[0:4]) + [get_fnumbs(tree[4])]

def grow_raw_tree():
    """
    Wrapper that sets the [0, 0] square and generates tree
    """
    return gen_tree([[{'coord': [0, 0], 'lineage': [], 'pindx':-1}]])

def grow_tree():
    """
    Return node tree with offspring links and shape values set
    """
    return add_pent_numbs(add_olinks(grow_raw_tree()))
