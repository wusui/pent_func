# Copyright (C) 2024 Warren Usui, MIT License
"""
Solve pentomino rectangles
"""
from itertools import chain
from grow_tree import grow_tree

def get_x_x_range(ydim):
    """
    Given a height, find the maximum height (y axis value) of the center
    of the X pentomino
    """
    return (7 - ydim) * 2 + ydim // 3

def mk_board(board_info):
    """
    Make a board (with x coordinate set)
    """
    def gen_xrow(xpattern):
        return ['-'] * xpattern[0] + ['X'] * xpattern[1] + \
                    ['-'] * xpattern[2]
    def mk_indiv_sz(board):
        def mk_b_rows(row_num):
            def gen_board(col_info):
                def mk_row(b_row_num):
                    if b_row_num in [row_num, row_num + 2]:
                        return gen_xrow([col_info + 1, 1,
                                       board['xdim'] - col_info - 2])
                    if b_row_num == row_num + 1:
                        return gen_xrow([col_info, 3,
                                       board['xdim'] - col_info - 3])
                    return ['-'] * board['xdim']
                return list(map(mk_row, range(0, board['ydim'])))
            return list(map(gen_board, range(board['xinfo'][row_num],
                    get_x_x_range(board['ydim']))))
        return [f"{board['ydim']}x{board['xdim']}",
                    list(chain.from_iterable(map(mk_b_rows,
                    range(0, len(board['xinfo'])))))]
    return dict(list(map(mk_indiv_sz, board_info)))

def get_boards():
    """
    Get boards -- Each different dimension of board is indexed in a dict.
    """
    def get_box_parms(ydim):
        def gbp_x(xdim):
            def gbp_rcnt(row_cnt):
                def gbp_rows(xstart):
                    return {'ydim': ydim, 'xdim': xdim,
                            'xinfo': xstart,
                            'xmax': get_x_x_range(ydim)}
                return gbp_rows(range(1, 1 - row_cnt, -1))
            return gbp_rcnt((ydim - 1) // 2)
        return gbp_x(60 // ydim)
    return mk_board(list(map(get_box_parms, range(3, 7))))

def find_origin(board):
    """
    Given a board, return the first open square scanning from left to
    right, top to bottom
    """
    def brd_coord(number):
        return [number % len(board), number // len(board)]
    return brd_coord(list(chain.from_iterable(zip(*board))).index('-'))

def set_points_in_board(b_info):
    """
    b_info is a dict containing three values:
        'points' -- list of points to get set
        'board' -- list of board rows
        'value' -- value to set points to
    """
    def spib_rows(row):
        def spibr_in(col):
            if [row[0], col[0]] in b_info['points']:
                if b_info['board'][row[0]][col[0]] == '-':
                    return b_info['value']
            return b_info['board'][row[0]][col[0]]
        return list(map(spibr_in, enumerate(row[1])))
    return list(map(spib_rows, enumerate(b_info['board'])))

def add_fig_brd(f_data):
    """
    f_data is a dict containing:
        'figure' -- a bottom-level tree node on the structure created
                    by grow_tree()
        'origin' -- coordinates of the first open square (see find_origin)
    """
    def afb_inner(origin):
        def pt_add(a_point):
            return [origin[0] + a_point[0], origin[1] + a_point[1]]
        return list(map(pt_add, f_data['figure']['lineage'])) + [
                    pt_add(f_data['figure']['coord'])]
    return afb_inner(f_data['origin'])

def chk_open_area(board):
    """
    Group the squares marked by 'A' while finding open spaces.
    """
    def coa_y(y_dim):
        def coa_x(x_dim):
            return [y_dim, x_dim]
        def increase_a_range(coord):
            if board[coord[0]][coord[1]] != '-':
                return board[coord[0]][coord[1]]
            if coord[0] > 0 and board[coord[0] - 1][coord[1]] == 'A':
                return 'A'
            if coord[1] > 0 and board[coord[0]][coord[1] - 1] == 'A':
                return 'A'
            if coord[0] < len(board) - 1 and \
                            board[coord[0] + 1][coord[1]] == 'A':
                return 'A'
            if coord[1] < len(board[0]) - 1 and \
                            board[coord[0]][coord[1] + 1] == 'A':
                return 'A'
            return '-'
        return list(map(increase_a_range,
                        map(coa_x, range(0, len(board[0])))))
    return list(map(coa_y, range(0, len(board))))

def count_open_area(board):
    """
    Count the size of the next open area on the board.
    """
    def cnt_oa(pbc):
        def cnt_next(newb):
            if list(chain.from_iterable(newb)).count('A') == pbc:
                return pbc
            return count_open_area(newb)
        return cnt_next(chk_open_area(board))
    return cnt_oa(list(chain.from_iterable(board)).count('A'))

def fragment_ok(board):
    """
    Return True if the next open area size is divisible by 5 (meaning
    pentominos can fit there).
    """
    return count_open_area(set_points_in_board(
                                    {'points': [find_origin(board)],
                                     'board': board,
                                     'value': 'A'})) % 5 == 0

def make_layout_cmd(sdata):
    """
    Generate a board layout with X placed
    """
    def mk_lyt_board_set(board_indx):
        def setup_brd(board):
            if fragment_ok(board):
                return [board]
            return []
        return [board_indx, list(chain.from_iterable(
                map(setup_brd, sdata[board_indx])))]
    return dict(list(map(mk_lyt_board_set, sdata.keys())))

def no_symmetry_issue(board):
    """
    Check for symmetry duplicates.  If there is a odd dimension and the X
    pentomino is centered on that dimension, flipping patterns can cause
    the same pattern to be counted more than once.  This makes sure that
    these duplicates are not counted by not counting layouts where most of
    the Y pentomino is on the left/top side.
    """
    def nsi_comp(nboard):
        def nsi_comp_v(mid_point):
            def nsi_last_comp(ynumbs):
                def ycount(ycnt):
                    return list(chain.from_iterable(ycnt)).count('Y')
                if ycount(ynumbs['top']) > ycount(ynumbs['bot']):
                    return False
                return True
            if nboard[mid_point].count('X') != 3:
                return True
            return nsi_last_comp({'top': nboard[0:mid_point],
                           'bot': nboard[mid_point + 1:]})
        return nsi_comp_v(len(nboard) // 2)
    if len(board[0]) == 10:
        return True
    if len(board[0]) % 2 == 1:
        return nsi_comp(list(zip(*board)))
    return nsi_comp(board)

def stringify(board):
    """
    Convert a board into one string
    """
    def sf_inner(linev):
        return ''.join(linev)
    return '\n'.join(list(map(sf_inner, board)))

def solve():
    """
    Recursively try to add tree representations of pentominos into a
    board pattern with the X pentomino already placed.
    """
    def tsolve(tree):
        def handle_one_layout_set(board):
            def is_bad(pstate):
                def is_bad_inner(f_node):
                    if f_node[0] >= len(pstate['board']):
                        return True
                    if f_node[1] >= len(pstate['board'][0]):
                        return True
                    if f_node[0] < 0:
                        return True
                    if f_node[1] < 0:
                        return True
                    if pstate['board'][f_node[0]][f_node[1]] != '-':
                        return True
                    return False
                return is_bad_inner([pstate['origin'][0] +
                                         pstate['coord'][0],
                                         pstate['origin'][1] +
                                         pstate['coord'][1]])
            def find_tiles(bstate):
                def f_gnodes():
                    def get_nnumb(en_data):
                        return en_data[0]
                    def chk_a_node(node_info):
                        if node_info[1]['pindx'] not in bstate['good_p']:
                            return False
                        return not is_bad({'origin': bstate['origin'],
                                               'coord': node_info[1]['coord'],
                                               'board': bstate['board']})
                    return list(map(get_nnumb, filter(chk_a_node,
                                    enumerate(tree[bstate['level']]))))
                def ft_inner(node_info):
                    if bstate['level'] == 4:
                        return node_info
                    return find_tiles({'level': bstate['level'] + 1,
                                    'good_p': node_info,
                                    'origin': find_origin(bstate['board']),
                                    'board': bstate['board']})
                return ft_inner(f_gnodes())
            def chk_ps(sboard):
                def dup_pc(pnumb):
                    if tree[4][pnumb]['figure'] in list(
                                    chain.from_iterable(sboard)):
                        return False
                    return True
                def do_recurs(fnumb):
                    def w_newboard(newboard):
                        if list(chain.from_iterable(newboard)).count('-') == 0:
                            if no_symmetry_issue(newboard):
                                return [stringify(newboard)]
                            return []
                        if fragment_ok(newboard) and no_symmetry_issue(
                                        newboard):
                            return chk_ps(newboard)
                        return []
                    return w_newboard(set_points_in_board(
                                {'points': add_fig_brd(
                                    {'figure': tree[4][fnumb],
                                    'origin': find_origin(sboard)}),
                                'board': sboard,
                                'value': tree[4][fnumb]['figure']}))
                def ldft():
                    return find_tiles({'level': 1, 'good_p': [0],
                                'origin': find_origin(sboard),
                                'board': sboard})
                return list(chain.from_iterable(
                    list(map(do_recurs, filter(dup_pc, ldft())))))
            return chk_ps(board)
        def ssolve(layout):
            def wrap_layouts(keyv):
                def layout_res(one_lset):
                    return list(map(handle_one_layout_set, one_lset))
                return [keyv, list(chain.from_iterable(
                                layout_res(layout[keyv])))]
            return list(map(wrap_layouts, layout.keys()))
        return dict(ssolve(make_layout_cmd(get_boards())))
    return tsolve(grow_tree())
