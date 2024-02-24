# Copyright (C) 2024 Warren Usui, MIT License
"""
See how long solve takes, and stash data in a json file.  Then check
all solutions found against a website and see that our solutions
match.
"""
from json import dumps
from itertools import chain
from datetime import datetime, timedelta
import requests
from solve import solve

def get_sq_list():
    """
    Return list of strings expressing the dimensions of the boards
    being checked.
    """
    def gen_sq_name(numb):
        return f"{numb}x{60 // numb}"

    return list(map(gen_sq_name, range(3, 7)))

def proc_sol(sstring):
    """
    Generate all four possible orientations of a board
    """
    def rflip(rows):
        return list(map(lambda a: a[::-1], rows))
    def ps_inner(rows):
        return [rows, rows[::-1], rflip(rows), rflip(rows[::-1])]

    return ps_inner(sstring.split())

def sol_org(webdata):
    """
    Take possible orientations and make all data extracted linear
    """
    return list(chain.from_iterable(list(map(proc_sol, webdata))))

def rejoin(sdata):
    """
    Convert list to string with parts separated by newline
    """
    return list(map('\n'.join, sdata))

def fix_info(data):
    """
    Wrap rejoin and sol_org calls
    """
    return sorted(rejoin(sol_org(data)))

def add_brk(sdata):
    """
    set breaks between boards while parsing scraped data
    """
    if "," in sdata:
        return 'xxxxx'
    return sdata

def verify(vpacket):
    """
    Compare our solutions with one found on the Internet.
        we -- our solutions
        outside -- Internet solutions
    """
    def parse_web(wdata):
        def pw_inner(webslv):
            def pwi_2(genslv):
                def comp_l_vs_i(pinfo):
                    return webslv[pinfo[0]] != pinfo[1]
                if len(list(set(genslv))) != len(genslv):
                    return 'Duplicate boards found'
                if list(filter(comp_l_vs_i, enumerate(genslv))):
                    return 'Error in pentomino generation'
                return []
            return pwi_2(fix_info(vpacket['we']))
        def ladd_brk():
            return '\n'.join(list(map(add_brk, wdata))).split('xxxxx')
        return pw_inner(fix_info(list(map(lambda a: a.strip(),
                         list(filter(None, ladd_brk()))))))
    return parse_web(vpacket['outside'].split())

if __name__ == "__main__":
    start = datetime.now()
    opents = solve()
    elapsed = datetime.now() - start
    secnds = elapsed / timedelta(seconds=1)
    print(f'solve() elapsed time = {secnds} seconds')
    with open('pentominos.json', 'w', encoding='utf-8') as outfile:
        outfile.write(dumps(opents, indent=4))
    for rectangle in get_sq_list():
        txt = ['https://isomerdesign.com/Pentomino/',
               f'{rectangle}', '/solutions.txt']
        pdata = requests.get(''.join(txt), timeout=30).text
        response = verify({'level': rectangle, 'we': opents[rectangle],
                      'outside': pdata})
        if response:
            print(response)
