#! /usr/bin/env python3
"""Take symbolic measure matrix and reconstruct possible isotopomer/cumomer/emu groups (as elementary as possible, idealy each group composed of only one i/c/e) for each combination of measurement methods
"""
from time import strftime, localtime, process_time as tproc
globals()["_T0"]=tproc()
def timeme(s="", dig=2):
    "if a global variable TIMEME is True and another global variable _T0 is set, print current CPU time relative to _T0. This printing is preceded by a message from 's'"
    if TIMEME:
        if "_T0" in globals():
            print(s, ":\tCPU=", round(tproc()-_T0, dig), "s", sep="")
        else:
            globals()["_T0"]=tproc()

TIMEME=False


import csv
import numpy as np

tol=np.finfo(float).eps*2**7 # 2.842170943040401e-14

from scipy import stats
import re
import sys, os
from io import StringIO
import argparse
import pandas as pa
from itertools import combinations, chain
from operator import itemgetter

#import isosolve
fve=os.path.join(os.path.dirname(__file__), "isosolve", "version.txt")
with open(fve, "r") as f:
    ver=f.read().rstrip()
    __version__=ver

from sympy import symbols, Symbol, rcollect, Matrix, factor, pprint, Eq, solve, Poly, Add, Mul, Pow, lambdify, IndexedBase, simplify, nan

from mdutils.mdutils import MdUtils
from mdutils import Html
import markdown as md

timeme("imports")
#print(sys.argv)
TEMPLATE = """<!DOCTYPE html>
<html>
<head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
    <meta name="referrer" content="no-referrer" />
    <meta name="referrer" content="unsafe-url" />
    <meta name="referrer" content="origin" />
    <meta name="referrer" content="no-referrer-when-downgrade" />
    <meta name="referrer" content="origin-when-cross-origin" />
    <title>SymNmrMs</title>
    <link href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" rel="stylesheet">
    <style>
        body {
            font-family: Helvetica,Arial,sans-serif;
        }
        code, pre {
            font-family: monospace;
        }
        table {
          border-collapse: collapse;
        }
        table, th, td {
          border: 1px solid black;
          padding: 5px;
        }
        th {
          background-color: #e0e0e0;
          text-align: center;
        }
        tr:nth-child(even) {background-color: #f2f2f2;}
    </style>
</head>
<body>
<div class="container">
{{content}}
</div>
<hr>
<h5>Legal information</h5>
<p style="font-size: small; padding: 5px">This content is produced by <tt>IsoSolve</tt> script ({{version}}) on {{date}} in {{duration}}.<br>
The script is written by Serguei Sokol &lt;sokol [at] insa-toulouse [dot] fr&gt; and it is a result of a collaboration between Mathematics Cell and MetaSys team, namely Pierre Millard at TBI (Toulouse Biotechnology Institute, Bio & Chemical Engineering).<br>
The copyright notice hereafter is relative to the script IsoSolve, not to the content of this page:<br>
Copyright © 2021, INRAE/INSA/CNRS.<br>
<tt>IsoSolve</tt> is released under <a href="https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html">GPL-2 license</a>.<br>
Follow this link for <a href="https://github.com/sgsokol/isosolve/issues">issue tracking</a>.
</p>
</body>
</html>
"""
class Object(object):
    "Object to which user can add fields"
    pass
def isstr(s):
    "return True if 's' is of type 'str'"
    return type(s) == str
def useq(seq):
    "return iterator of unique elements in a sequence preserving initial order"
    seen=set()
    return (x for x in seq if not (x in seen or seen.add(x)))
def is_intstr(s):
    "check if a string 's' represents exactly an integer number"
    try:
        i=int(s)
    except:
        return False
    return str(i)==s
def assign(name, value, d):
    "assign a 'value' to 'name' in dictionary 'd', e.g. 'd' can be locals(). Return 'value'."
    d[name] = value
    return value
def revenum(x):
    "reversed enumerated iterable x"
    return reversed(list(enumerate(x)))
def pdiff(expr, vs):
    "partial derivatives of expr by vars in vs"
    return [expr.diff(v) for v in vs]
def partcomb(x, p):
    "return iterator of iterator tuples ([p distinct elements of x], [the rest len(x)-p distinct elements of x])"
    from itertools import combinations, compress
    import numpy as np
    xn=len(x)
    if p > xn or p < 0:
        return iter(())
    elif p == xn:
        return iter([(iter(x), iter(()))])
    elif p == 0:
        return iter([(iter(()), iter(x))])
    xlog=np.zeros(xn, bool)
    ix=np.array(range(xn))
    for ixp in combinations(ix, p):
        xlog[:]=False
        xlog[list(ixp)]=True
        yield (compress(x, xlog), compress(x, ~xlog))
def limi(li, negative=False):
    "apply or not sign '-' to all terms of the list li"
    if not negative:
        return li
    return [-v for v in li]
def e2li(e):
    "transform expression 'e' in a nested list. Each nested level corresponds to Add, Mul and Pow arguments, Add being the highest level."
    from sympy import Add, Mul, Matrix, Pow, Atom
    if isinstance(e, Atom):
        return [[(e,1)]]
    at=[]
    for v in Add.make_args(e):
        mli=[]
        for u in Mul.make_args(v):
            if type(u) is Mul:
                u=Mul.make_args(term)
            else:
                u=[u]
            # split each factor in couples (a,b) for a^b
            u=[q.args if type(q) == Pow else (q,1) for q in u]
            mli.extend(u)
        at.append(mli) # each add-term is a list of factors
    return at
def lfact(e, esub={}, cache={}, scache={}, pr=False, nsub={}, fast=False):
    "factorize terms in expression 'e'"
    from sympy import Add, Mul, Matrix, Pow, Atom
    #if pr:
    #    import pdb; pdb.set_trace()
    if isinstance(e, Atom):
        return e
    #if str(e) in ("(a*(b + c) + d)/(a - 1)",):
    #    import pdb; pdb.set_trace()
    cl=type(e)
    if cl == Matrix:
        return Matrix([[lfact(v, esub, cache, scache, pr, nsub, fast=fast) for v in li] for li in e.tolist()])
    if cl != list:
        res=cache.get(e, None)
        if res:
            try:
                if nsub and np.abs(e.subs(nsub).evalf() - res.subs(nsub).evalf()) > 1.e-14:
                    print("num e=", e, "\nres_cache=", res)
                    import pdb; pdb.set_trace()
            except:
                print("no num e=", e, "\nres_cache=", res)
                import pdb; pdb.set_trace()
            if pr:
                print("from cache e=", e, "\nres=", res)
            return res
    if cl != Pow and cl != Mul and cl != Add and cl != list:
        cache[e]=e
        return e
    if cl == Mul:
        # recursive call on cofactors
        if fast:
            res=Mul(*[lfact(v, esub, cache, scache, pr, nsub, fast=fast) for v in Mul.make_args(e)])
        else:
            res=Mul(*[lfact(v, esub, cache, scache, pr, nsub, fast=fast) for v in Mul.make_args(e.cancel())])
        cache[e]=res
        return res
    if cl == Pow:
        # call on base
        res=Pow(lfact(e.args[0], esub, cache, scache, pr, nsub, fast=fast), e.args[1])
        cache[e]=res
        return res
    if cl == list:
        # already decomposed => extract the most frequent factor
        
        # make the base unique in each factor group
        ee=[]
        for fli in e:
            dbase={} # base dictionary: base-> power
            for ft in fli:
                if ft[0] in dbase:
                    dbase[ft[0]] += ft[1]
                elif -ft[0] in dbase:
                    po=dbase[-ft[0]]+ft[1]
                    del dbase[-ft[0]]
                    dbase[ft[0]]=po
                    dbase[-1]=dbase.get(-1, 0)+1
                else:
                    dbase[ft[0]]=ft[1]
            ee.append(list(dbase.items()))
        e=ee
        n=len(e)
        ee=Add(*[Mul(*[Pow(*v) for v in li]) for li in e])
        if ee in cache:
            return cache[ee]

        # register factors in each add term -> cnt={fct: [indexes in at]}
        cnt={}
        #import pdb; pdb.set_trace()
        for ia,li in enumerate(e):
            for im,t in enumerate(li):
                if t[0] == 1 or t[0] == -1:
                    #import pdb; pdb.set_trace()
                    continue
                cnt[t[0]]=cnt.get(t[0], [])
                cnt[-(t[0])]=cnt.get(-(t[0]), [])
                cnt[t[0]].append((ia,im,t[1],False))
                cnt[-(t[0])].append((ia,im,t[1],True))
        le=[(t,len(li)) for t,li in cnt.items()]
        if not le:
            # only 1s
            cache[ee]=ee
            return ee
        pw=dict((t, min(li, key=lambda i: i[2])[2]) for t,li in cnt.items()) # base-> minimal power
        wh,n=max(le, key=lambda i: i[1])
        cnt_wh=cnt[wh]
        pw=pw[wh] # keep only minimal power
        wh=Pow(wh, pw)
        # indexes of added terms in fterm
        iawh=set(tu[0] for tu in cnt_wh)
        rest=[li for i,li in enumerate(e) if i not in iawh]
        #print("wh=", wh, "; n=", n, "; cnt_wh=", cnt_wh)
        if n == 1:
            # nothing to factorize => check if any Add-type cofactor is in the rest of add list
            ali=[Mul(*[Pow(*v) for v in ml]) for ml in e] # add list
            alim=[-v for v in ali] # -ali
            found=False
            fminus=False
            for ia,fli in enumerate(e):
                if type(ali[ia]) != Mul:
                    continue
                for ifa,(b,p) in enumerate(fli):
                    if type(b) != Add or p < 0:
                        continue
                    # check if b terms are in the ali (nb b cannot coincide with ali[ia] as b is at best an element of add)
                    ib=iainb(b.args, ali)
                    if ib:
                        found=True
                        break
                    # try -ali
                    ib=iainb(b.args, alim)
                    if ib:
                        found=True
                        fminus=True
                        break
                if found:
                    break
            if found:
                cof=b*(Mul(*[Pow(v[0], v[1] if iv != ifa else v[1]-1) for iv,v in enumerate(fli)])+(-1 if fminus else 1))+Add(*[v for iv,v in enumerate(ali) if iv not in ib and iv != ia])
                cof=rsubs(cof, esub, scache)
            else:
                cof=rsubs(ee, esub, scache)
            #print("e=", ee, "; res=", cof)
            if nsub and np.abs(ee.subs(nsub).evalf() - cof.subs(nsub).evalf()) > 1.e-14:
                print("num e=", ee, "\nres_cof=", cof)
                import pdb; pdb.set_trace()
            if pr:
                print("e=", ee, "\nres_cof=", cof)
            cache[ee]=cof
            return cof
        # factorized term
        cof=[]
        for ia,im,pm,fmin in cnt_wh:
            tpw=e[ia][im][1]-pw # power which rests after simplification
            tpw=[(e[ia][im][0], tpw)] if tpw else [] # wh factor which rests
            tmp=e[ia][:im]+tpw+e[ia][(im+1):] or [(1,1)]
            tmp=Mul(*(Pow(*v) for v in tmp))
            if fmin:
                tmp=-tmp
            cof.append(tmp)
        #print("ee=", ee, "\nwh=", wh, ", cof list=", cof)
        cof=rsubs(Add(*cof), esub, scache)
        if type(cof) == Add:
            cof=lfact(cof, esub, cache, scache, pr, nsub, fast=fast)
        #print("cof factorized=", cof)
        cof=cof.args if type(cof) == Pow else (cof, 1)
        wh=rsubs(wh, esub, scache)
        wh=wh.args if type(wh) == Pow else (wh, 1)
        if wh[0] == cof[0]:
            fterm=[(wh[0],wh[1]+cof[1])]
        elif wh[0] == -cof[0]:
            fterm=[(-1,1), (wh[0],wh[1]+cof[1])]
        elif ('real' in dir(cof[0]) and cof[0] < 0) or ('could_extract_minus_sign' in dir(cof[0]) and cof[0].could_extract_minus_sign()):
            if wh[0].could_extract_minus_sign():
                fterm=[(-wh[0],wh[1]), (-cof[0],cof[1])]
            else:
                fterm=[(-1,1), wh, (-cof[0],cof[1])]
        else:
            if wh[0].could_extract_minus_sign():
                fterm=[(-1,1), (-wh[0],wh[1]), cof]
            else:
                fterm=[wh, cof]
        fterm=Mul(*[Pow(rsubs(v[0], esub, scache), v[1]) for v in fterm])
        #print("after wh=", wh, ", cof=", cof, ", fterm=", fterm)
        #print("rest=", rest, "; fterm=", fterm, "e=", e)
        if rest:
            #import inspect
            #if len(inspect.stack()) > 100:
            #    print("hm, stack is too big on fterm=", fterm, "; and rest=", rest)
            #    if len(inspect.stack()) > 103:
            #        raise Exception("duh")
            if fterm == 0:
                # can happen after simplifications by esub
                res=lfact(rest, esub, cache, scache, pr, nsub, fast=fast)
                try:
                    if nsub and np.abs(ee.subs(nsub).evalf() - res.subs(nsub).evalf()) > 1.e-14:
                        print("num e=", ee, "\nres_rest0=", res)
                        import pdb; pdb.set_trace()
                except:
                    print("no num e=", ee, "\nres_rest0=", res)
                    import pdb; pdb.set_trace()
                if pr:
                    print("e=", ee, "\nres_rest0=", res)
                cache[ee]=res
                return res
            else:
                #import inspect
                #if len(inspect.stack()) > 40:
                #    print("hm, stack is too big on fterm=", fterm, "\nrest=", rest, "\nee=", ee, "\nres=", Add(fterm, lfact(rest, esub, cache, scache, pr, nsub, fast=fast)))
                #    #raise Exception("duh")
                #    #import pdb; pdb.set_trace()
                #if len(inspect.stack()) > 46:
                #    #raise Exception("duh")
                #    import pdb; pdb.set_trace()
                res=Add(fterm, lfact(rest, esub, cache, scache, pr, nsub, fast=fast))
                if res != ee and res not in lfact.estack:
                    res=lfact(res, esub, cache, scache, pr, nsub, fast=fast)
                try:
                    if nsub and np.abs(ee.subs(nsub).evalf() - res.subs(nsub).evalf()) > 1.e-14:
                        print("num e=", ee, "\nres_rest1=", res)
                        import pdb; pdb.set_trace()
                except:
                    print("no num e=", ee, "\nres_rest1=", res)
                    import pdb; pdb.set_trace()
                if pr:
                    print("e=", ee, "\nres_rest1=", res)
                cache[ee]=res
                return res
        else:
            #import inspect
            #if len(inspect.stack()) > 10:
            #    print("hm, stack is too big on fterm=", fterm, "; and e=", e)
            #    raise Exception("duh")
            res=fterm
            if nsub and np.abs(ee.subs(nsub).evalf() - res.subs(nsub).evalf()) > 1.e-14:
                print("num e=", ee, "\nres_fterm=", res)
                import pdb; pdb.set_trace()
            if pr:
                print("e=", ee, "\nres_fterm=", res)
            cache[ee]=res
            return res
    # here e is expr => decompose in list of lists and recall lfact()
    lfact.estack=getattr(lfact, "estack", {})
    if e in lfact.estack:
        import pdb; pdb.set_trace()
        raise Exception(f"cycling call is detected for e='{e}'")
    lfact.estack[e]=None
    
    if type(esub) != type({}):
        esub=dict(esub)
    at=Add.make_args(e)
    li=[]
    if len(at) <= 1:
        res=e
    else:
        for term in at:
            li.extend(e2li(lfact(term, esub, cache, scache, pr, nsub, fast=fast)))
        res=lfact(li, esub, cache, scache, pr, nsub, fast=fast)
    cache[e]=res
    try:
        if nsub and np.abs(e.subs(nsub).evalf() - res.subs(nsub).evalf()) > 1.e-12:
            print("e=", e, "\nres=", res, "\n")
            print("e.num=", e.subs(nsub), "\nres.num=", res.subs(nsub), "\n")
            import pdb; pdb.set_trace()
    except:
        print("no num e=", e, "\nres=", res, "\n")
        import pdb; pdb.set_trace()
    if pr:
        print("e=", e, "\nres=", res, "\n")
    #print("cache=", cache)
    del(lfact.estack[e])
    return res
def iainb(a, b):
    "return a list of indexes of a in b or None if a is not strictly in b. An empty list is returned for empty a. Repeated values in a must be repeated at least the same number of times in b."
    if len(a) > len(b):
        return None
    # dictionary of indexes in b -- {item: list_of_indexes}
    dib={}
    for i,v in enumerate(b):
        dib[v]=dib.get(v, [])
        dib[v].append(i)
    ilast={}
    res=[None]*len(a) # place-holder for results
    for ia,x in enumerate(a):
        #ib=b.index(x, ilast.get(x, 0))
        if x not in dib:
            return None # x is not in b at all
        ilast[x]=ilast.get(x, 0)
        if ilast[x] == len(dib[x]):
            return None # x is more present in a than in b
        res[ia]=dib[x][ilast[x]]
        ilast[x] += 1
    return res
def rsubs(e, d, c={}):
    "restricted substitution of keys from 'd' found in 'e' by d[key]. c is an optional cache dict"
    if not d:
        return e
    from sympy import Add, Mul, Pow, Matrix, Atom
    cl=type(e)
    if type(d) != type({}):
        d=dict(d)
    if type(e) == Matrix:
        #import pdb; pdb.set_trace()
        return Matrix([[rsubs(v, d, c) for v in li] for li in e.tolist()])
    res=c.get(e, None)
    if res:
        return res
    if e in d:
        res=d[e]
        c[e]=res
        return res
    elif -e in d:
        res=-d[-e]
        c[e]=res
        return res
    elif isinstance(e, Atom):
        res=e
        c[e]=res
        return res
    elif isinstance(-e, Atom):
        res=e
        c[e]=res
        return res
    elif cl == Pow:
        ve=e.args
        res=Pow(rsubs(ve[0], d, c), ve[1])
        c[e]=res
        return res
    elif cl == Add or cl == Mul:
        ve=list(rsubs(v, d, c) for v in cl.make_args(e))
        vme=limi(ve, True)
        vli=[(v,cl.make_args(k)) for k,v in d.items()]
        #print("ve=", ve, ", vme=", vme, ", vli=", vli)
        for v,li in vli:
            ie=iainb(li, ve)
            if ie:
                #import pdb; pdb.set_trace()
                ie=sorted(ie, reverse=True)
                for i in ie:
                    ve.pop(i)
                    vme.pop(i)
                ve.append(v)
                vme.append(-v)
                continue
            ie=iainb(li, vme)
            if ie:
                ie=sorted(ie, reverse=True)
                for i in ie:
                    ve.pop(i)
                    vme.pop(i)
                #print("ie=", ie, ", ve cut=", ve)
                ve.append(-v)
                vme.append(v)
                #print("ve new=", ve)
        res=cl(*ve)
        #print("e=", e, "\nres=", res, "\n")
    else:
        res=e
    c[e]=res
    return res
def setsub(l, it, v):
    "in a list l set values with indexes from it to v, i.e. l[it]=v where it is iterator. v is recycled if needed"
    for ii,i in enumerate(it):
        l[i]=v[ii%len(v)]
def hascumo(s1, s2):
    "return a cumulated string (or None if not conformant) for binary cumomers in s1 and s2. E.g. 01x+11x => x1x"
    if len(s1) != len(s2) or len(s1) == 0:
        return None
    res=list(s1)
    ndif=0
    for i,c1 in enumerate(s1):
        c2=s2[i]
        if c1 == c2:
            res[i]=c1
        elif c1 in "01" and c2 in "01" and ndif == 0:
            ndif += 1
            res[i]="x"
        else:
            return None
    return "".join(res)
def lico2tb(lico, mdFile=None, text_align="center"):
    "transform list of columns in lico to a flat list of strings suitable for md table writing. Write md table if mdFile is set."
    # check that all column have the same non zero length
    if not lico:
        return ""
    n0=len(lico[0])
    if n0 == 0:
        return ""
    if not all(len(li)==n0 for li in lico):
        raise Exception("Not all columns have the same length in lico")
    tb=list(chain(*map(list, zip(*[[str(it) for it in li] for li in lico]))))
    if mdFile:
        mdFile.new_table(columns=len(lico), rows=n0, text=tb, text_align=text_align)
    return tb
def ech(m, d={}, c={}, sc={}, nd={}, ld={}, fast=False):
    "reduce augmented matrix m to echelon form applying substitutions from d. Return echelon matrix and permutation list of col-index couples. No check is made for non zero terms in the last column in all zero rows."
    if type(m) != Matrix:
        raise Exception("expected Matrix type on input")
    nr,nc=m.shape
    mres=m.copy()
    #print("mres init")
    #pprint(mres)
    pe=[] # column permutations
    per=[] # row permutations
    if nr*nc == 0 or nr*(nc-1) == 0:
        return (mres, pe)
    nro=min(nr, nc-1)
    if not nd:
        nd=dict((v, np.random.rand(1)[0]) for v in mres.free_symbols)
    mn=np.array(mres[:,:-1].subs(nd)).astype(np.double)
    #print("mn init=", mn)
    # find pivots on numeric matrix
    for i in range(nro):
        #if i == 15:
        #    import pdb; pdb.set_trace()
        # search for shortest non zero rows
        nz1d=list(zip(np.count_nonzero(np.round(mn[i:,i:], 10), axis=1), np.min(np.abs(mn[i:,i:]-1), axis=1)))
        io=sorted(range(len(nz1d)),key=nz1d.__getitem__)
        ip=[ir for ir in io if nz1d[ir][0] != 0]
        if not ip:
            # non zero was not found at all => the end
            break
        ip=ip[0]
        jp=np.nonzero(np.round(mn[i+ip,i:], 10))[0]
        if jp.size > 1:
            # choose the column closest to 1
            jp=jp[np.argmin(np.abs(mn[i+ip,i+jp]-1))]
        else:
            jp=jp[0]
        #print("i=", i, "ip=", ip, "jp=", jp, "mn=\n", np.round(mn, 2))
        #if i == 5:
        #    import pdb; pdb.set_trace()
        #i1=np.where(np.round(np.abs(mn[i:,i:]), 12) == 1)
        #if not len(i1[0]):
            # 1 was not found, search for max.abs non zero
            #mi=mn[i:,i:].min()
            #ma=mn[i:,i:].max()
            #mima=mi if -mi > ma else ma
            #if round(mima, 12) == 0:
            #    # non zero was not found at all => the end
            #    break
            #i1=np.where(mn[i:,i:] == mima)
        if ip != 0:
            mn[[i, ip+i]]=mn[[ip+i, i]]
            per.append((i, ip+i))
            mres.row_swap(i, ip+i)
        if jp != 0:
            mn[:,[i, jp+i]]=mn[:,[jp+i, i]]
            pe.append((i, jp+i))
            mres.col_swap(i, jp+i)
        if mn[i, i] != 1:
            mn[i, i:] *= 1./mn[i, i]
        for ii in range(i+1, nr):
            if round(mn[ii,i], 12):
                mn[ii,(i+1):] -= mn[ii,i]*mn[i,(i+1):]
        #if any(np.abs(mn[i, :]) > 1.e10):
        #    import pdb; pdb.set_trace()
        mn[i+1:,i]=0.
        #print("i=", i, "mn perm=\n", np.round(mn, 2))
        #print("pe=", pe, "\nper=", per)
    if round(mn[i,i], 12) == 0:
        i -= 1
    #import pdb; pdb.set_trace()
    #timeme("mn")
    #import pdb; pdb.set_trace()
    nr=min(i+1, nro)
    mres=mres[:nr,:]
    #print("mres perm=")
    #pprint(mres)
    for i in range(nr):
        if mres[i, i] != 1:
            mres[i, (i+1):] *= 1/mres[i, i]
            #if any(v==nan for v in mres[i, (i+1):]):
            #    import pdb; pdb.set_trace()
            mres[i, (i+1):]=lfact(mres[i, (i+1):], d, c, sc, fast=fast)
            mres[i, i]=1
        vr=mres[i, (i+1):]
        #print("vr=", vr)
        for ii in range(i+1, nr):
            if mres[ii,i]:
                mres[ii,(i+1):] -= mres[ii,i]*vr
                mres[ii,(i+1):]=lfact(mres[ii,(i+1):], d, c, sc, fast=fast)
        mres[i+1:, i]=Matrix.zeros(nr-i-1, 1)
        #mres=lfact(mres, d, c, sc, fast=fast)
        #print("mres i=", i)
        #pprint(mres)
    #xgr=Matrix(ld["xgr"])
    #print("m*x-b=", m[:,:-1].subs(nd)*xgr-m[:,-1].subs(nd))
    #print("mres*x-b=", mres[:,:-1].subs(nd)*xgr.permute(pe)-mres[:,-1].subs(nd))
    return (mres, pe)
def echsol(m, *unk, pe=[], d={}, c={}, sc={}, nsub={}, fast=False):
    """solves augmented linear system in echelon form resulted from ech(...).
    Unknown symbols are taken from \*unk. pe is a possible permutation list of ij-tuples from ech(...).
    Return a symbolic column vector"""
    if type(m) != Matrix:
        raise Exception("expected Matrix type on input")
    nr,nc=m.shape
    sol=Matrix(range(nc-1)) # place-holder for solution
    if nc == 1:
        # empty matrix
        return sol
    m=m.copy()
    #import pdb; pdb.set_trace()
    if nr < nc-1:
        # under-determined system
        sol[nr:,0]=unk[nr:]
        m[:,-1] -= m[:,nr:-1]*sol[nr:,0]
        m[:,-1]=lfact(m[:,-1], d, c, sc, fast=fast)
    sol[nr-1,0]=m[nr-1,-1]
    for i in range(nr-2, -1, -1):
        m[:(i+1),-1] -= m[:(i+1),i+1]*sol[i+1,0]
        m[:(i+1),-1]=lfact(m[:(i+1),-1], d, c, sc, fast=fast)
        sol[i,0]=m[i,-1]
    for i,j in reversed(pe):
        sol.row_swap(i,j)
    return sol
def sumi2sd(covmat, i, tol=tol):
    "calculate sd for a new variable defined as sum of components in 'i'"
    i=np.array(i)
    val=np.sum(covmat[i[:,None],i])
    val=val if val >= 0. else 0. if val >= -tol else val
    return np.sqrt(val)
def new_header(level=None, title="", mdFile=None):
    if mdFile is None or level is None:
        return
    mdFile.new_line('<a name="'+title.lower().replace(" ", "-")+'"></a>')
    mdFile.new_header(level=level, title=title)
    return
def dhms(f, dig=2):
    "Convert float number 'f' representing seconds to a 'D days Hh Mm Ss"
    f=round(f, dig)
    if f < 60:
        return f"{f}s"
    elif f < 3600:
        return f"{int(f//60)}m {round(f%60, dig)}s"
    elif f < 24*3600:
        return f"{int(f//3600)}h {int(f%3600//60)}m {round(f%60, dig)}s"
    else:
        d=f//(24*3600)
        return f"{int(d)} day{'s' if d > 1 else ''} {int(f%(24*3600)//3600)}h {int(f%3600//60)}m {round(f%60, dig)}s"
def main_arg():
    "Wrapper for `main(li=sys.argv[1:])` call when used from executable script"
    main(li=sys.argv[1:])
def main(li=None, mm=None, colsel=None, tim=None, data=None, w=None, s=None, rd=None, path=None, fast=False, vers=None, inchi=False):
    """Resolve symbolically and numerically iso-, cumo- and EMU-mers from isotope labeling measurements (NMR, MS, ...)
    
    :param li: list of argument to be processed, typically a ``sys.argv[1:]``, cf. ``isosolve -h`` or :ref:`cli` section for parameter significance. Defaults to None.
    :type li: list, optional
    :param mm: mapping matrix describes symbolically the mapping of each measurements method to isotopic space. Can be a file name in tsv format or a pandas data.frame. Either this parameter must be set or the only positional argument in li. Defaults to None.
    :type mm: str or pandas.FataFrame, optional
    :param colsel: if str, coma separated column selection from mm. The syntax is described in ``isosolve -h``. If iterable, a collection of items describing column selection. Each item can be an integer (positive or negative for exclusion), slice or string with column name. Column numbers are 1-based (not 0-based as usual for Python). Defaults to None which is equivalent to proceeding all the columns from mm.
    :type colsel: str or iterable, optional
    :param tim: if True, print CPU time after some key steps. Defaults to None.
    :type tim: logical, optional
    :param data: if DataFrame, numeric values and sd of measurements described in mm. If str, a file name in tsv format with such data in three columns: name, value, sd. Defaults to None.
    :type data: str or pandas.DataFrame, optional
    :param w: if True, write formatted results in .md (plain text, MarkDown format) and .html files. The basis of file name is the same as in mm when mm is string, or ``mm`` if mm is a DataFrame. The name ``mm`` can be overwritten with path parameter (cf. hereafter). Defaults to None which is equivalent to True when used in command line and to False when used programmatically.
    :type w: logical, optional
    :param s: seed for pseudo random generation. When used, randomly drawn values become reproducible. Only useful for debugging or issue reporting. Defaults to None.
    :type s: int, optional
    :param rd: if True, use randomly generated values for checking obtained formulas. In the output HTML file, numerical values that are not coinciding are typeset in **bold**. Only useful for debugging and issue reporting. Defaults to None.
    :type rd: logical, optional
    :param path: directory (when ends with "/") of file name (otherwise) for formatted output. Can be useful when mm is a DataFrame. Defaults to None which is equivalent to the current working directory.
    :type path: str, optional
    :param fast: if True, skip calculation of elementary measurable combinations.
    :type fast: logical, optional
    :param vers: if True, print the version number and return without any result.
    :type vers: logical, optional
    :param inchi: if True, writes TSV files with International Chemical Identifiers (InChI) for isotopomers, cumomers, EMUs and elementary measurable combinations.
    :type vers: logical, optional
    :return: dictionary with the following entries:
    
        :xbc: 1-column DataFrame, randomly drawn values for isotopomers;
        :xv: dict, key=measurement name and value=numeric value deduced from xbc;
        :am: sympy Matrix, augmented matrix of system equations. The last column is the right hand side;
        :rdn_mes: list, names of redundant measurements;
        :smm: DataFrame, mapping matrix after column selection;
        :mmd: dict, key=measurement name, value=numpy array with isotopomer indexes mapped onto this measurements;
        :uvals: dict, key=method name (i.e. column name in mm), value=list of unique measurement names involved in this method;
        :sy_meas: list, symbolic equations defining measurements;
        :sbc: dict, key=isotopomer name (i.e. row name in mm), value=isotopomer symbolic variable;
        :sols: list, symbolic solution for isotopomers. The order is the same as in mm;
        :vsubs: dict, key=full symbolic expression, value=simplified symbolic expression. Used for simplification procedures;
        :metrics: dict, the fields are:
        
            :idef: list, defined isotopomer indexes;
            :idef_c: set, indexes of defined cumomers;
            :idef_e: dict, key=EMU name, value=set of defined M+i mass isotopologues in this EMU;
            
        :emusol: dict, key=EMU name, value=list of symbolic solutions for M+i;
        :cusol: dict, key=cumomer name, value=symbolic solution;
        :ls: dict, least square solution. The fields are:
        
            :iso: list, isotopomer values and sd;
            :cov: isotopomer covariance matrix;
            :cumo: list, cumomer values and sd;
            :emu: list of dicts. Each dict has a key=EMU name and value=list of values (the first dict) or sd (the second one);
            :mec: list, measurable elementary combinations values and sd;

        :inchi: dict, InChI information. The fields are:
        
            :iso: DataFrame for isotopomers;
            :cumo: DataFrame for cumomers;
            :emu: DataFrame for EMUs;
            :mcomb: DataFrame for measurable elementary combinations;
            
    :rtype: dict
    """
    #print("li=", li)
    #timeme("init")
    globals()["_T0"]=tproc()
    if vers:
        print(ver)
        return
    if type(li) == list:
        class DefaultHelpParser(argparse.ArgumentParser):
            def error(self, message):
                sys.stderr.write('error: %s\n' % message)
                self.print_help()
                sys.exit(2)
        parser = DefaultHelpParser(description='calculate isotopmer/cumomer/EMU expressions (symbolic and optionally numeric) by a combination of isotope labeling measurements (NMR, MS, ...)', formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('mm',
            help="""file name (tsv) providing measure matrix which is organized as follows:
    - each row corresponds to a given binary cumomer, e.g. '010x'
    - each column corresponds to a given experiment type, e.g. 'HSQC-C_beta'
    - each cell contains a variable name, e.g. 'a', 'b' etc. or empty, i.e. NA
    - each variable name can appear multiple times in a given
     column but not in different columns. If a variable appears
     in several rows, these rows are considered as mixed in
     equal parts in a given measurement method. E.g. in HSQC-C_beta,
     contributions of patterns '001x' and '101x' are mixed
     in one variable, say 'k', and contributions of '011x'
     and '111x' are mixed in another one, say 'l'. Other rows
     in this column must be left empty as they don't
     contribute to any measurement.""")
        parser.add_argument("-c", "--colsel", help="""column selection (full set by default).
    Can be a slice, comma separated list of integers or names,
    regular expression, or a mix of all this. In a slice,
    a negative number means counting from the end. In a list,
    if given negative numbers or names starting with '-',
    the corresponding columns are excluded from treatment.
    Names can be given as Python regular expressions.
    The order of column selection is irrelevant.

    Examples:
    '1:3' - first 3 columns;
    ':3' - the same;
    ':-1' - all columns but the last;
    '2::2' - even columns;
    '1,3,6' - first, third and sixth columns;
    'HSQC.*' - columns with names started by 'HSQC';
    '-HN.*' - all columns except starting with 'HN'""")
        parser.add_argument("-t", "--TIMEME", help="activate or not (default) CPU time printing. Useful only for debugging or issue reporting.", action='store_true')
        parser.add_argument("-d", "--data", help="file name with 3 columns: name, value and sd. Numeric values in columns 'value', and 'sd' must be non-negative and positive respectively. Fields are tab-separated, comments start with sharp-sign '#'")
        parser.add_argument("-w", "--write", help="force .md and .html file writing even in non CLI mode", action='store_true')
        parser.add_argument("-s", "--seed", help="integer value used as a seed for pseudo-random drawing. Useful only for debugging or issue reporting.", type=int)
        parser.add_argument("-r", "--rand", help="make random draws for numerical tests of formulas. Useful only for debugging or issue reporting.", action='store_true')
        parser.add_argument("-p", "--path", help="path for .md and html files. If it ends by '/' it is interpreted as a directory path which is created if nonexistent. Otherwise, it is interpreted as a base before .md and .html extensions. If not given, .md and .html are written in the same directory with the same basename (before the extension) as 'MM' file ")
        parser.add_argument("-f", "--fast", help="skip calculations of measurable combinations", action='store_true')
        parser.add_argument("-i", "--inchi", help="write InChI files", action='store_true')
        parser.add_argument("-v", "--version", help="print version number on stdout and exit. Useful only for debugging or issue reporting.", action='version', version=f'%(prog)s {ver}')
        
        args = parser.parse_args(li)
        fmm=args.mm
        fdata=args.data
        jc=args.colsel
        wri=args.write or __name__ == "__main__" or __name__[:8] == "isosolve"
        np.random.seed(args.seed)
        globals()["TIMEME"]=args.TIMEME
        rdraw=args.rand
        path=args.path
        fast=args.fast
        inchi=args.inchi
        #timeme("parser")
    # overwrite options with explicit parameters
    alist=[]
    smm=None
    dmm=None
    if type(mm) in (str, pa.core.frame.DataFrame) or not li:
        if isstr(mm):
            fmm=mm
            alist.append(f"mm='{mm}'")
            smm=pa.read_csv(fmm, delimiter="\t", comment="#", dtype=str)
            smm.set_index(smm.columns[0], inplace=True)
        elif type(mm) == pa.core.frame.DataFrame:
            # mm is supposed to be data.frame
            fmm="mm.par"
            smm=mm
            alist.append(f"mm='data.frame {'x'.join(str(v) for v in smm.shape)}'")
        elif type(li) == list and len(li) == 0:
            parser.print_help()
            return 1
        else:
            raise Exception(f"unknown type for 'mm'. Expected 'str' or pandas data.frame got {type(mm)}")
    if li and smm is None:
        smm=pa.read_csv(fmm, delimiter="\t", comment="#", dtype=str)
        smm.set_index(smm.columns[0], inplace=True)
    if colsel or not li:
        jc=colsel
        alist.append(f"colsel='{colsel}'")
    if not tim is None or not li:
        globals()["TIMEME"]=tim
        alist.append(f"tim={tim}")
        # get numeric measurement matrix
    if type(data) in (str, pa.core.frame.DataFrame) or not li:
        if isstr(data):
            fdata=data
            alist.append(f"data='{data}'")
            dmm=pa.read_csv(fdata, delimiter="\t", comment="#")
        elif type(data) == pa.core.frame.DataFrame:
            # data is supposed to be data.frame
            fdata="mm.data"
            dmm=data
            alist.append(f"data='data.frame {'x'.join(str(v) for v in dmm.shape)}'")
        elif data is None:
            fdata=None
        else:
            raise Exception(f"unknown type for 'data'. Expected None, 'str' or pandas data.frame got {type(data)}")
    if li and dmm is None and fdata:
        dmm=pa.read_csv(fdata, delimiter="\t", comment="#")

    if not w is None or not li:
        wri=w
        alist.append(f"w={w}")
    if s or not li:
        np.random.seed(s)
        alist.append(f"s={s}")
    if rd or not li:
        rdraw=rd
        alist.append(f"rd={rd}")

    # md: prepare markdown output
    if wri or inchi:
        nm_md=".".join(fmm.split(".")[:-1])
        if path:
            os.makedirs(path, exist_ok=True)
            if path[-1:] == "/":
                nm_md=os.path.join(path, os.path.basename(nm_md))
            else:
                nm_md=path
    if wri:
        mdFile = MdUtils(file_name=nm_md, title="Symbolic Resolution of Labeled Measurements")
        new_header(level=1, title='SymNmrMs', mdFile=mdFile)
        # md: h2 problem setup
        new_header(level=2, title='Problem setup', mdFile=mdFile)
        # md: h3 arguments
        if li:
            new_header(level=3, title='Used arguments', mdFile=mdFile)
            mdFile.insert_code(" ".join("'"+v+"'" for v in li))
        if alist:
            new_header(level=3, title='Explicit Parameters', mdFile=mdFile)
            mdFile.insert_code(" ".join(alist))
    
    # get symbolic measurement matrix
    nr,nc=smm.shape
    nm_bc=smm.index
    nm_exp=smm.columns
    if fdata:
        if dmm.shape[1] != 3:
            raise Exception(f"Expected 3 columns in data file '{fdata}' instead got {dmm.shape[1]}")
        dmm.set_index(dmm.columns[0], inplace=True)
        dmm.columns=["value", "sd"]
    
    ljc=[]
    allj=list(range(nc))
    if isstr(jc) and jc:
        # jc can be like "1:3,4,-1,:5,:-2,-MS.*"
        tmp={}
        ljc=[tmp["vs"] for v in jc.split(",") if assign("vs", v.strip(), tmp)]
    elif jc:
        ljc=[str(v) for v in jc]
    if ljc:
        # get litteral column names and stripe them off
        jc=[nm_exp.get_loc(ljc.pop(i)) for i,v in revenum(ljc) if v in nm_exp]
        # get regexp column names and strip them off
        for i,v in revenum(ljc.copy()):
            # list of matched columns by v
            lima=[nm_exp.get_loc(vma) for vma in nm_exp if re.match(v, vma)]
            jc += lima
            if lima:
                ljc.pop(i)
        # get integers convertible and in good range
        jc += [int(ljc.pop(i))-1 for i,v in revenum(ljc) if is_intstr(v) and int(v) > 0 and int(v) <= nc]
        # get slices
        for i,it in revenum(ljc):
            if it.count(":") == 0:
                continue
            # get begin, end, step of the slice
            bes=it.split(":")
            if len(bes) > 3:
                raise Exception(f"slice '{it}' is badly formed: too many ':'")
            s="" if len(bes) < 3 else bes[2]
            b,e=bes[:2]
            try:
                b=int(b) if b else 0
                e=int(e) if e else nc
                s=int(s) if s else 1
            except:
                raise Exception(f"slice '{it}' is badly formed: not integer values")
            if b > 0:
                b -= 1
            li=allj[slice(b, e, s)]
            jc += li
            ljc.pop(i)
        jc=sorted(set(jc))
        if not jc:
            jc=allj # default: all
        # remove "negative" litteral, regexp and integer indexes
        for i,it in revenum(ljc):
            if not it.startswith("-"):
                continue
            #import pdb; pdb.set_trace()
            it=it[1:]
            # litteral
            if it in nm_exp:
                j=nm_exp.get_loc(it)
                if j in jc:
                    jc.remove(j)
                ljc.pop(i)
                continue
            # regexp
            lima=[nm_exp.get_loc(vma) for vma in nm_exp if re.match(it, vma)]
            for j in lima:
                if j in jc:
                    jc.remove(j)
            if lima:
                ljc.pop(i)
                continue
            # integer
            if is_intstr(it):
                j=int(it)-1
                if j in jc:
                    jc.remove(j)
                ljc.pop(i)
                continue
    else:
        jc=allj
    if ljc:
        sys.stderr.write(f"Warning: following identifiers were not used as not matching any column:\n\t"+"\n\t".join(ljc)+"\n")
    smm=smm.iloc[:,jc]
    nr,nc=smm.shape
    nm_bc=smm.index
    nm_exp=smm.columns
    flen=len(nm_bc[0]) # total fragment length
    clen=nm_bc[0].count("1")+nm_bc[0].count("0") # account for 01 atoms
    if nr != 2**clen:
        raise Exception(f"not sufficient row number in '{fmm}': expected {2**clen}, got {nr}") 
        
    # unique values per column -> dict uvals
    uvals=dict((nm,list(v for v in useq(col) if isstr(v))) for nm,col in smm.items())
    # measurement name to column translator
    v2col=dict((v, nm) for nm,li in uvals.items() for v in li)
    # make dict measure matrix mmd: 'a'->integer indexes of '010x',...
    mmd=dict()
    for nmc, col in smm.items():
        for v in uvals[nmc]:
            if v in mmd:
                raise Exception(f"variable '{v}' appears in different columns")
            mmd[v]=np.where(col==v)[0]
    # MS list of bc indexes
    # cumomer list of bc indexes
    i01=np.where([c=="1" or c=="0" for c in nm_bc[0]])[0]
    cumo=dict()
    cumo_i=dict()
    maskx=list("x"*len(nm_bc[0]))
    for w in range(1,clen+1): # cumomer weight 1:clen
        for i in combinations(i01, w):
            cu=maskx.copy()
            setsub(cu, i, "1")
            nm_cu="".join(cu)
            cumo_i[nm_cu],cumo[nm_cu]=zip(*((ii,nm_bc[ii]) for ii in range(nr) if all(nm_bc[ii][ic] == "1" for ic in i)))
    
    # build a list of linear equations => sy_all
    # dict of symbolic vars 'a': symbol('a'), ...
    sv=dict((v,symbols(v, real=True)) for v in mmd.keys())
    # dict of symbolic bcumos '000x': symbol('bc000x')
    sbc=dict((bc,symbols("i"+bc, real=True)) for bc in nm_bc)
    sbc2i=dict((sbc[bc], i) for i,bc in enumerate(nm_bc))
    bc2i=dict((bc, i) for i,bc in enumerate(nm_bc))
    # expression: sum of bcN (will be 1)
    e_sumbc=sum(sbc.values())
    # expression: sum of bc by v
    e_bcv=dict((v, sum(sbc[nm_bc[i]] for i in inx)) for v,inx in mmd.items())
    # ex: value's sum of each colum in smm
    e_colv=dict((n, sum(sv[v] for v in uvals[n])) for n,col in smm.items())
    sy_col1=[v-1 for v in e_colv.values()]
    # ex: each measure without sf
    #import pdb; pdb.set_trace()
    e_colbc=dict((n, sum(sbc[nm_bc[i]] for i,m in enumerate(col) if isstr(m))) if sum(isstr(v) for v in col) < nr else (n, 1) for n,col in smm.items())
    # ex: cumo 1xx: 10x+11x -> dict e_cubc
    e_cubc=dict((k, sum(sbc[bc] for bc in v)) for k,v in cumo.items())
    # ex: EMU {EEx: list(M+0, M+1, M+2)}, for M+1 we have 10x+01x -> dict e_emu
    e_emu=dict()
    emu_i=dict() # index list for each MID element of each emu
    for le in range(1, clen+1):
        for isub in combinations(i01, le):
            cemu="".join("E" if i in isub else "x" for i in range(flen))
            e_emu[cemu]=[None]*(le+1)
            emu_i[cemu]=[None]*(le+1)
            for i in range(le+1):
                r=itemgetter(*isub)
                emu_i[cemu][i]=[ibc for ibc,bc in enumerate(nm_bc) if r(bc).count("1") == i]
                e_emu[cemu][i]=sum(sbc[nm_bc[ii]] for ii in emu_i[cemu][i])
    sy_meas=[sv[v]*e_colbc[n]-e_bcv[v] for n,col in smm.items() for v in uvals[n]]        
    # generate random bc values => xbc
    xbc=np.random.rand(nr)
    xbc=pa.DataFrame(xbc/sum(xbc), index=nm_bc)
    nbcsub=dict((bc,xbc[0][nm]) for nm,bc in sbc.items())
    # dict of measures a, b, etc. => xv
    xv=dict((v, (e_bcv[v]/e_colbc[v2col[v]]).subs(nbcsub)) for v in mmd)
    nvsub=dict((sv[v],xv[v]) for v in mmd)
    nsub=nvsub.copy()
    nsub.update(nbcsub)
    # dict of provided data
    if fdata:
        #import pdb; pdb.set_trace()
        import nlsic
        setval=set(mmd.keys())
        setdata=set(dmm.index)
        if setval > setdata:
            bad=sorted(setval-setdata)
            raise Exception(f"following measurements are not provided in '{fdata}':\n\t"+"\n\t".join(bad))
        # normalize values and sd by column
        ndata={}
        #print(dmm)
        #print(dmm.index)
        for nm in nm_exp:
            vli=uvals[nm]
            ix=np.in1d(dmm.index, vli)
            s=np.sum(dmm["value"][ix])
            #print("vli=", vli, "; s=", s, "; ix=", ix)
            #import pdb; pdb.set_trace()
            if s == 0:
                raise Exception(f"the following values from '{fdata}' sum up to 0:\n\t"+"\n\t".join(sorted(vli)))
            dmm.iloc[ix, 0] /= s
            dmm.iloc[ix, 1] /= s
            ndata.update((sv[v], dmm.loc[v, "value"]) for v in vli)
            #print(dmm)
        # prepare fresid() for nlsic optimization
        e_sim=dict((v,e_bcv[v]/e_colbc[v2col[v]]) for v in mmd)
        # replace iXXX by par[]
        ipar=IndexedBase("par")
        dtmp=dict((sbc[bc], ipar[i]) for i,bc in enumerate(nm_bc))
        e_sim=Matrix([v.subs(dtmp) for v in e_sim.values()])
        fsim=lambdify(ipar, e_sim, 'numpy')
        fjac=lambdify(ipar, e_sim.jacobian([ipar[i] for i in range(nr)]), 'numpy')
        def fresid(par, cjac=False, *kargs, **kwargs):
            res=Object()
            res.sim=np.array(fsim(par.flat), float)
            res.sim=res.sim.reshape(res.sim.size)
            res.res=(res.sim-kwargs["meas"])/kwargs["vsd"]
            res.res=res.res.reshape(res.res.size)
            if cjac:
                res.jacobian=np.array(fjac(par.flat), float)/kwargs["vsd"][:,None]
            return res
        # prepare inequalities iXXX >= 0
        u=np.eye(nr)
        cu=np.zeros((nr,1))
        # prepare sum(iXXX)=1
        em=np.ones((1, nr))
        eco=np.ones((1,1))
        # least sq --ln solution -> bcls
        bcls=np.ones(nr)/nr # init
        #import pdb; pdb.set_trace()
        meas=np.asarray(dmm.loc[mmd.keys(), "value"])
        vsd=np.asarray(dmm.loc[mmd.keys(), "sd"])
        lsres=nlsic.nlsic(bcls, fresid, u=u, co=cu, control={"trace": False, "tolx": min(dmm.loc[:,"sd"])/10./nr}, e=em, eco=eco, flsi=nlsic.lsi_ln, meas=meas, vsd=vsd)
        bcls=lsres.par
        fres=fresid(bcls, meas=meas, vsd=vsd) # by order of mmd.keys()
        nbcls=dict((sbc[bc], bcls[i,0]) for i,bc in enumerate(nm_bc))
        # estimate covbc
        u,sva,vt=np.linalg.svd(lsres.a, full_matrices=False)
        ra=sum(sva > sva[0]*tol)
        racov=sum(sva > sva[0]*np.sqrt(tol))
        svi,vt=1./sva[:racov],vt[:racov,:]
        #import warnings
        #warnings.simplefilter("error")
        sdfact=np.sqrt(racov/(racov-1.) if racov > 1 else 2.)
        #covbc=lsres.nte.dot(vt.T.dot(np.einsum("ij,i->ij", u.T,sdfact*svi))) # sd(resid) is supposed to be 1 because of diff[i]/sd[i]
        covbc=lsres.nte.dot((np.asarray(vt.T)*(sdfact*svi)))
        covbc=covbc.dot(covbc.T)
        bcsd=np.sqrt(np.diag(covbc))
        #import pdb; pdb.set_trace()
        timeme("least sq")

    # make the rows unique in smm. Indexes of repeated rows will be concatenated with '_': '01_10' etc.
    o=np.array(sorted(range(nr), key=lambda i: list(map(str, list(smm.iloc[i,:])))))
    io=np.arange(nr)
    io[o]=io.copy()# inverse ordering
    dup=smm.iloc[o, :].duplicated()
    # indexes of duplicated rows -> idr
    idr=[]
    for i,v in enumerate(dup):
        if not v:
            # start new group
            idr.append([i])
        else:
            # add to the current group (the last one)
            idr[-1].append(i)
    # index of unique rows in original smm -> iur
    idr=sorted(idr, key=lambda v: o[v[0]])
    # old index to grouped index -> ibc2gr
    ibc2gr=dict((o[v], i) for i,li in enumerate(idr) for v in li)
    iur=[o[v[0]] for v in idr]
    # names of grouped rows
    nmgr=["_".join(smm.index[o[v]]) for v in idr]
    ngr=len(nmgr)
    smm=smm.iloc[iur, :]
    smm.index=nmgr
    xgr=np.array([np.sum(xbc.iloc[o[v],0]) for v in idr])
    # md: h3 table with input measurements
    if wri:
        new_header(3, "Available measurements", mdFile=mdFile)
        mdFile.new_line()
        tb=[[""]+[v.replace("_", "+") for v in nmgr]]+[[nm]+[it if isstr(it) else ""  for it in c] for nm,c in smm.items()]
        lico2tb(tb, mdFile)
    vsubs=dict(\
        [(e_colv[nmc], 1) for nmc in smm]+\
        [(sum(sv[b] for b in uvals[nmc] if b != a)-1, -sv[a]) for nmc in smm if len(uvals[nmc]) == 2 for a in uvals[nmc]]+\
        [(sum(sv[b] for b in uvals[nmc] if b != a), 1-sv[a]) for nmc in smm if len(uvals[nmc]) >= 3 for a in uvals[nmc]]+\
        [(sum(sv[b] for b in it2)-1, -sum(sv[b] for b in it1)) for nmc in smm if len(uvals[nmc]) >= 4 for it1, it2 in partcomb(uvals[nmc], 2)]+\
        [(sum(sv[b] for b in it2), 1-sum(sv[b] for b in it1)) for nmc in smm if len(uvals[nmc]) >= 5 for it1, it2 in partcomb(uvals[nmc], 2)]+\
        [(e_bcv[v], sv[v]) for nmc in smm if e_colbc[nmc] == 1 for v in uvals[nmc]]
        )
    vcache={} # for factorizations
    scache={} # for substitutions
    am=Matrix.ones(1,ngr+1)
    nm_r=["all1"]
    nmlong=sorted(nm_exp, key=lambda nm: str(e_colbc[nm]))
    for n in nmlong:
        col=smm[n]
        for iv,v in enumerate(uvals[n]):
            ve=Matrix.zeros(1,ngr+1)
            # first put plain 1s
            if e_colbc[n] == 1:
                for j in mmd[v]:
                    ve[0,ibc2gr[j]]=1
                ve[0,ngr]=sv[v]
            else:
                for j in range(ngr):
                    if isstr(col[j]):
                        ve[0,j] = 1
                for j in mmd[v]:
                    ve[0,ibc2gr[j]]=1-1/sv[v]
            am=am.vstack(am, ve)
            nm_r.append(v)
    
    timeme("pb setup")
    if wri:
        # md: h3 equations
        new_header(level=3, title="Equation system", mdFile=mdFile)
        mdFile.new_line("i"+" + i".join(nm_bc)+" = 1")
        for v in sy_meas:
            mdFile.new_line(repr(v).replace("*", "·")+" = 0")
    
    ame,pe=ech(am, vsubs, vcache, scache, nsub, ld=locals(), fast=False)
    timeme("echelon")

    unk=[sbc[v] if v in sbc else symbols("i"+v) for i,v in enumerate(nmgr)]
    if pe:
        # apply permutations
        for i,j in pe:
            unk[j], unk[i]=unk[i], unk[j]
    sol=echsol(ame, *unk, pe=pe, d=vsubs, c=vcache, sc=scache, nsub=nsub, fast=False)
    if ngr < nr:
        grsub={}
        for v in unk:
            grsub[v]=sum(sbc[u] for u in str(v)[1:].split("_"))
        sols=[lfact(v, vsubs, vcache, scache, fast=False).subs(grsub) for v in sol]
    else:
        sols=[lfact(v, vsubs, vcache, scache, fast=False) for v in sol]

    # detect redundant measurements
    #rdn_meas=set(sv)-set(str(s) for v in sols for s in v.free_symbols)
    # NB. rdn_meas may be updated after the final solution is found
    # expand grouped terms
    sole=[None]*nr
    for i,nm in enumerate(nmgr):
        if "_" in nm:
            vnm=nm.split("_")
            sole[bc2i[vnm[0]]]=sols[i]-sum(sbc[v] for v in vnm[1:])
            for v in vnm[1:]:
                sole[bc2i[v]]=sbc[v]
        else:
            sole[bc2i[nm]]=sols[i]
    sols=sole
    bcset=set(sbc.values())
    idef=[i for i in range(nr) if not (sols[i].free_symbols & bcset)]
    timeme("solution")

    bcsub=[(sbc[bc], sols[i]) for i,bc in enumerate(nm_bc)]
    if rdraw:
        nsol=[v.subs(nsub) for v in sols]
    # calculate formulas for redundant measurements
    rdn_meas=list(mmd)
    redfrm=dict((v,lfact((e_bcv[v]/e_colbc[v2col[v]]).subs(bcsub), vsubs, vcache, scache, fast=False)) for v in rdn_meas)
    # if a formula has iXXX in it or includes any of measurements of the same column, eliminate the corresponding measurements
    todel=[]
    for k,v in redfrm.items():
        if "free_symbols" in dir(v):
            fsy=v.free_symbols
        else:
            continue
        if fsy & bcset:
            rdn_meas.remove(k)
            todel.append(k)
            continue
        fsy=set(str(v) for v in fsy)
        if fsy & set(uvals[v2col[k]]):
            rdn_meas.remove(k)
            todel.append(k)
    rdn_meas=sorted(rdn_meas)
    for k in todel:
        del(redfrm[k])
    rdnsub=dict((v,sv[k]) for k,v in redfrm.items() if type(v) in (Add, Mul, Pow))
    sols=[v.subs(rdnsub) for v in sols]

    # md: h2 Solution
    if wri:
        new_header(level=2, title="Solution", mdFile=mdFile)
        # md: h3 basic solution
        new_header(level=3, title="Isotopomers", mdFile=mdFile)
        tb=[["variable"]+list(nm_bc)]+ \
            [["formula"]+[re.sub("(i[01x]+)", "**\\1**", str(s).replace("**", "^").replace("*", "·")) for s in sols]]
        if rdraw:
            rd=[str(round(v, 3)).rstrip(".0") for v in xbc[0]]
            fa=[str(round(v, 3)).rstrip(".0") for v in nsol]
            fa=[(v if v == rd[i] else "**"+v+"**") for i,v in enumerate(fa)]

            tb += [["random draw"]+rd]
            tb += [["formula applied"]+fa]
            timeme("num test")
        lico2tb(tb, mdFile)
    
    res={"xbc": xbc, "xv": xv, "nsub": nsub, "am": am, "ame": ame, "rdn_meas": rdn_meas}
    if rdraw:
        res["nsol"]=nsol
    timeme("rdn_meas")
    # substitute bcs in e_emu by sol
    emusol=dict()
    idef_e=dict()
    for k,li in e_emu.items():
        emusol[k]=[None]*len(li)
        idef_e[k]=[]
        for i,v in enumerate(li):
            term=lfact(v.subs(bcsub), vsubs, vcache, scache, fast=False)
            emusol[k][i]=term.subs(rdnsub)
            if not term.free_symbols & bcset:
                idef_e[k].append(i)
        idef_e[k]=set(idef_e[k])
    timeme("emu")
    # substitute bcs in e_cubc by sol
    cusol=dict()
    idef_c=[]
    for i,(k,v) in enumerate(e_cubc.items()):
        term=lfact(v.subs(bcsub), vsubs, vcache, scache, fast=False)
        cusol[k]=term.subs(rdnsub)
        if not (term.free_symbols & bcset):
            idef_c.append(i)
    idef_c=set(idef_c)

    # md: h3 EMU solution
    if wri:
        new_header(level=3, title="EMU", mdFile=mdFile)
        # tb is a list of columns (clen+2 length) which is transposed and chained and the end
        tb=[["variable"]+["M+"+str(i) for i in range(clen+1)]] # first column
        for k,li in emusol.items():
            tb.append([k]+[(re.sub("(i[01x]+)", "**\\1**", str(li[i]).replace("**", "^").replace("*", "·")) if i < len(li) else "") for i in range(clen+1)])
        lico2tb(tb, mdFile)
            
        # md: h3 cumomer solution
        new_header(level=3, title="Cumomers", mdFile=mdFile)
        tb=[["variable"]+list(cusol.keys())]+[["formula"]+[re.sub("(i[01x]+)", "**\\1**", str(it).replace("**", "^").replace("*", "·")) for it in cusol.values()]]
        lico2tb(tb, mdFile)
    timeme("cumo")

    # decide about just-, over- or under-determined
    status="undefined"
    alls=set(it for e in sols for it in e.free_symbols)
    bc_free=sorted(str(it) for it in (bcset & alls))
    if bc_free:
        status="under-determined"
        reason=f"following {len(bc_free)} isotopomers could not be determined: *"+"*, *".join([v[1:] for v in bc_free])+"*"
    elif rdn_meas:
        status="over-determined"
        reason="following measurements were not used: *"+"*, *".join(rdn_meas)+"*"
    else:
        status="just-determined"
        reason=""
    if wri:
        # md: h2 status
        new_header(level=2, title="System Status", mdFile=mdFile)
        mdFile.new_paragraph("The system is *"+status+"*.")
        if reason:
            mdFile.new_paragraph("The reason is that "+reason+".")
        mdFile.new_paragraph("Number of redundant measurements: "+str(len(rdn_meas))+".")
        if rdn_meas:
            # md: h2 Redundant measurements
            new_header(level=2, title="Redundant measurements", mdFile=mdFile)
            tb=[["Measurement"]+rdn_meas]+ \
            [["Formula"]+[str(redfrm[k]).replace("**", "^").replace("*", "·") for k in rdn_meas]]
            text_align=["center", "center"]
            if rdraw:
                vr=[str(round(sv[v].subs(nsub), 3)).rstrip(".0") for v in rdn_meas]
                fa=[str(round(redfrm[v].subs(nsub) if "subs" in dir(redfrm[v]) else redfrm[v], 3)).rstrip(".0") for v in rdn_meas]
                fa=["**"+v+"**" if v != vr[i] else v for i,v in enumerate(fa)]
                tb += [["Random draw"]+vr]
                tb += [["Formula applied"]+fa]
                text_align += ["center", "center"]
            tb += [["Methods"]+[v2col[v]+": *"+", ".join(sorted(set(v2col[str(vv)] for vv in (redfrm[v].free_symbols if "free_symbols" in dir(redfrm[v]) else []))))+"*" for v in rdn_meas]]
            text_align += ["left"]
            lico2tb(tb, mdFile, text_align=text_align)
    
    # for under-determined, system find bc sum that can be measured.
    if not fast:
        compact=list()
        if status == "under-determined":
            # partial derivatives -> pder
            pder=Matrix([pdiff(s, sbc.values()) for s in sols])
            try:
                pn=np.array(pder).astype(np.byte) # see if it in {-1,0,1}
            except TypeError:
                pn=np.array(pder.subs(nsub)).astype(np.float32)
            inz=np.where(np.any(pn, axis=1))[0]
            ifr=np.array(sorted(bc2i[it[1:]] for it in bc_free))
            # remove all zero rows and cols, as well as ifr rows
            inzf=np.delete(inz, np.where(np.isin(inz, ifr)))
            pn=pn[inzf,:][:,ifr]
            # run through combinations summing columns to 0 or -1. If -1, add corresponding ifr index to combi
            zcomb=[]
            for le in range(1, pn.shape[0]):
                for isub in combinations(range(pn.shape[0]), le):
                    cs=pn[isub,:].sum(axis=0) # column sums
                    if not np.all(np.logical_or(cs == 0, cs == -1)):
                        continue
                    # check that isub cannot be decomposed in sets already in zcom
                    im1=np.where(cs == -1)[0]
                    isub=np.sort(np.hstack((inzf[list(isub)], ifr[im1]))) if len(im1) else inzf[list(isub)]
                    found=False
                    for iold in zcomb:
                        if len(iold) >= len(isub):
                            continue
                        if np.all(np.isin(iold, isub, True)):
                            # if previous iold in it => invalidate the whole sequence as non elementary sum
                            found=True
                            break
                    if not found:
                        #print("append=", i)
                        zcomb.append(isub)
                #timeme(f"combi2 le={le}")
            zcomb=sorted(zcomb, key=lambda v: len(v))
            # remove non elementary combs
            n=len(zcomb)
            for i in range(n):
                for j in range(n-1, i, -1):
                    if len(zcomb[i]) == len(zcomb[j]):
                        continue # as i cannot be in j, they are necessary different
                    if np.all(np.isin(zcomb[i], zcomb[j], True)):
                        # remove j-th entry
                        del(zcomb[j])
                        n -= 1
                if i >= n:
                    break
            # print measurable combinations
            if zcomb:
                # symbolic combinations
                scomb=[lfact(sum(sols[i] for i in co).expand(), vsubs, vcache, scache, pr=False, fast=False) for co in zcomb]
                # compacted combinations: 001x+101x=>x01x
                for li in zcomb:
                    lic=[nm_bc[i] for i in li]
                    while len(lic) > 1:
                        #print("lic=", lic)
                        cres=list()
                        ifound=[]
                        for i in range(len(lic)-1):
                            if i in ifound:
                                continue
                            for j in range(i+1,len(lic)):
                                found=hascumo(lic[i], lic[j])
                                if found:
                                    cres.append(found)
                                    ifound.append(i)
                                    ifound.append(j)
                                    break
                            if found:
                                continue
                        if cres:
                            cres.extend(lic[i] for i in range(len(lic)) if i not in ifound)
                            #print("cres=", cres)
                            lic=cres
                        else:
                            break
                    compact.append(cres or lic)
                #print(compact)
                
                for ie,e in enumerate(scomb):
                    #print(("starting e=", e))
                    for bc in sbc.values():
                        e=e.apart(bc)
                        co=e.coeff(bc)
                        if co == 0:
                            continue
                        if np.abs(co.subs(nvsub)) < 1.e-14:
                            cp=Poly(e, bc).coeffs()
                            cp[0]=0
                            e=Poly(cp, bc).as_expr()
                        scomb[ie]=lfact(e, vsubs, vcache, scache, fast=False)
                ncomb=[e.subs(nvsub) for e in scomb]
                timeme("combi")
            else:
                scomb=[]
                ncomb=[]
        else:
            idef=list(range(nr))
            zcomb=[]
            scomb=[]
            ncomb=[]
        if wri and (zcomb or len(idef)):
            # md: h2 Measurable elementary combinations
            #import pdb; pdb.set_trace()
            new_header(level=2, title="Measurable elementary combinations", mdFile=mdFile)
            tb=[["N°"]+list(range(1, len(idef)+1))+list(range(len(idef)+1, len(idef)+len(zcomb)+1))]
            tb += [["Combination"]+[nm_bc[i] for i in idef]+[" + ".join(nm_bc[i] for i in co) for co in zcomb]]+ \
                [["Accumulated"]+[nm_bc[i] for i in idef]+[" + ".join(li) for li in compact]]+ \
                [["Formula"]+[str(sols[i]).replace("**", "^").replace("*", "·") for i in idef]+[str(it).replace("**", "^").replace("*", "·") for it in scomb]]
            if rdraw:
                rd=[str(round(xbc[0][i], 3)).rstrip(".0") for i in idef]+[str(round(sum(xbc[0][i] for i in co), 3)).rstrip(".0") for co in zcomb]
                fa=[str(round(sols[i].subs(nsub), 3)).rstrip(".0") for i in idef]+[str(round(it, 3) if it.is_constant() else it.evalf(3)).rstrip(".0") for it in ncomb]
                fa=[(v if v == rd[i] else "**"+v+"**") for i,v in enumerate(fa)]
                tb += [["Randomly Drawn Values"]+rd]
                tb += [["Formula Applied"]+fa]
            tb += [["Methods"]+[", ".join(sorted(set(v2col[str(v)] for v in sols[i].free_symbols))) for i in idef]+[", ".join(sorted(set(v2col[str(v)] for v in it.free_symbols))) for it in scomb]]
            lico2tb(tb, mdFile)
    if fdata:
        cumols=[np.sum(bcls[li,0]) for li in cumo_i.values()]
        cumosd=[sumi2sd(covbc, li) for li in cumo_i.values()]
        emuls=dict()
        emusd=dict()
        for k,li in emu_i.items():
            emuls[k]=[np.sum(bcls[v,0]) for v in li]
            emusd[k]=[sumi2sd(covbc, v) for v in li]
        if len(idef) or zcomb:
            mecls=[bcls[i,0] for i in idef]+[np.sum(bcls[co]) for co in zcomb]
            mecsd=[bcsd[i] for i in idef]+[sumi2sd(covbc, co) for co in zcomb]
    if fdata:
        nfree=len(mmd)-racov
        chi2=sum(fres.res**2)
        tb=["Item", "Value", "Comment",
            "number of measurements", str(len(mmd)), "m",
            "system rank (number of statistically defined parameters)", str(racov), "p",
            "degree of freedom", str(nfree), "dof=m - p",
            f"χ²({nfree})", str(round(chi2, 3)), "Σ((estimated-measured)/sd)²"+"<br/>95% confidence interval is [0; "+(str(round(stats.chi2.ppf(0.95, nfree))) if nfree else "NaN")+"]",
            f"χ²({nfree})/{nfree}", str(round(chi2/nfree, 3)) if nfree else "NaN", "Reduced χ², χ²(dof)/dof: if ≫ 1 then poor fitting, ~ 1 means 'good' fitting and ≪ 1 means over-fitting or overestimating sd",
            f"p-value of one-tail χ² test, i.e.<br/>P(χ²({nfree}) > {round(chi2, 3)})", str(round(1-stats.chi2.cdf(chi2, nfree), 3)) if nfree else "NaN","Value close to 0 (e.g. under 0.05) means poor fitting. Value close to 1 can be an evidence for over-fitting or that sd are overestimated. It can be NaN (not a number) for dof=0"
        ]
        tb3=[tb[i*3:(i+1)*3] for i in range(len(tb)//3)]
        res["chi2"]=pa.DataFrame(tb3[1:], columns=tb3[0])
        res["chi2"]["Value"]=pa.to_numeric(res["chi2"]["Value"])
        res["chi2"].at[3, "Value"]=chi2
        res["chi2"].at[4, "Value"]=chi2/nfree
        res["chi2"].at[5, "Value"]=1-stats.chi2.cdf(chi2, nfree)
    if wri and fdata:
        # md: h2 least squares
        new_header(level=2, title=("Minimum Norm " if len(idef) < nr else "")+"Least Squares", mdFile=mdFile)
        # md: h3 Summary
        new_header(3, "Problem summary", mdFile=mdFile)
        mdFile.new_table(columns=3, rows=len(tb)//3, text=tb)
        
        # md: h3 table with data
        new_header(3, "Measured values", mdFile=mdFile)
        mdFile.new_line()
        nmv=list(mmd.keys())
        tb=[["name"]+nmv]
        tb += [["sd"]+[round(dmm.loc[it, "sd"], 5) for it in nmv]]
        tb += [["value"]+[round(dmm.loc[it, "value"], 3) for it in nmv]]
        tb += [["estimated"]+[round(fres.sim[i], 3) for i in range(len(nmv))]]
        tb += [["(estimated - value)/sd"]+[np.round(fres.res[i], 5) for i in range(len(nmv))]]
        #tb += [["value"]+[round(xv[it], 5) for it in nmv]]
        #tb += [["sd"]+[round(0.001, 5) for it in nmv]]
        lico2tb(tb, mdFile)
        # md: h3 Redundant measurements
        if rdn_meas:
            new_header(level=3, title="Redundant measurements", mdFile=mdFile)
            tb=[["Measurement"]+rdn_meas]
            vm=[dmm.loc[v,"value"] for v in rdn_meas]
            fa=[redfrm[v].subs(ndata) if "subs" in dir(redfrm[v]) else redfrm[v] for v in rdn_meas]
            tb += [["Measured Value"]+[round(v, 3) for v in vm]]
            tb += [["Formula Applied"]+[round(v, 3) for v in fa]]
            tb += [["Δ=FA-MV"]+[round(fa[i]-v, 3) for i,v in enumerate(vm)]]
            lico2tb(tb, mdFile)
        # md: h3 iso
        new_header(level=3, title="Isotopomers", mdFile=mdFile)
        tb=[["variable"]+list(nm_bc)]
        tb += [["value"]+[np.round(bcls[i,0], 3) if i in idef else "-" for i in range(nr)]]
        tb += [["sd"]+[np.round(bcsd[i], 5) if i in idef else "-" for i in range(nr)]]
        lico2tb(tb, mdFile)
        # md: h3 cumo
        new_header(level=3, title="Cumomers", mdFile=mdFile)
        tb=[["variable"]+list(cumo_i.keys())]
        tb += [["value"]+[np.round(cumols[i], 3) if i in idef_c else "-" for i,li in enumerate(cumo_i.values())]]
        tb += [["sd"]+[np.round(cumosd[i], 5) if i in idef_c else "-" for i,li in enumerate(cumo_i.values())]]
        lico2tb(tb, mdFile)
        # md: h3 emu val
        new_header(level=3, title="EMU values", mdFile=mdFile)
        tb=[["variable"]+["M+"+str(i) for i in range(clen+1)]] # first column
        for k,li in emu_i.items():
            tb += [[k]+[round(emuls[k][i], 3) if i in idef_e[k] else "-" for i,v in enumerate(li)]+[""]*(clen+1-len(li))]
        lico2tb(tb, mdFile)
        # md: h3 emu sd
        new_header(level=3, title="EMU standard deviations", mdFile=mdFile)
        tb=[["variable"]+["M+"+str(i) for i in range(clen+1)]] # first column
        for k,li in emu_i.items():
            tb += [[k]+[round(emusd[k][i], 5) if i in idef_e[k] else "-" for i,v in enumerate(li)]+[""]*(clen+1-len(li))]
        lico2tb(tb, mdFile)
        # md: h3 Measurable elementary combinations
        if len(idef) or zcomb:
            new_header(level=3, title="Measurable combinations", mdFile=mdFile)
            tb=[["Combination"]+[nm_bc[i] for i in idef]+[" + ".join(nm_bc[i] for i in co) for co in zcomb]]
            tb += [["Accumulated"]+[nm_bc[i] for i in idef]+[" + ".join(li) for li in compact]]
            tb += [["value"]+[str(round(v, 3)).rstrip(".0") for v in mecls]]
            tb += [["sd"]+[str(round(v, 5)).rstrip(".0") for v in mecsd]]
            lico2tb(tb, mdFile)
    if inchi:
        res["inchi"]={}
        tb=[["isotopomer"]+list(nm_bc)]
        tb += [["InChI"]+["/i"+",".join(str(k+1)+"+"+sy for k,sy in enumerate(i) if sy in "01") for i in nm_bc]]
        res["inchi"]["iso"]=pa.DataFrame(tb[1][1:], index=tb[0][1:], columns=["inchi"])
        if fdata:
            tb += [["value"]+[str(bcls[i,0]) if i in idef else "" for i in range(nr)]]
            tb += [["sd"]+[str(bcsd[i]) if i in idef else "" for i in range(nr)]]
            res["inchi"]["iso"]["value"]=pa.to_numeric(tb[2][1:])
            res["inchi"]["iso"]["sd"]=pa.to_numeric(tb[3][1:])
        with open(nm_md+"_inchi_iso.tsv", "w") as f:
            f.write("\n".join("\t".join(row) for row in zip(*tb)))
            f.write("\n")
        
        tb=[["cumomer"]+list(cumo_i.keys())]
        tb += [["InChI"]+["/i"+",".join(str(k+1)+"+"+sy for k,sy in enumerate(cu) if sy == "1") for cu in cumo_i.keys()]]
        res["inchi"]["cumo"]=pa.DataFrame(tb[1][1:], index=tb[0][1:], columns=["inchi"])
        if fdata:
            tb += [["value"]+[str(cumols[i]) if i in idef_c else "" for i in range(len(cumo_i))]]
            tb += [["sd"]+[str(cumosd[i]) if i in idef_c else "" for i in range(len(cumo_i))]]
            res["inchi"]["cumo"]["value"]=pa.to_numeric(tb[2][1:])
            res["inchi"]["cumo"]["sd"]=pa.to_numeric(tb[3][1:])
        with open(nm_md+"_inchi_cumo.tsv", "w") as f:
            f.write("\n".join("\t".join(row) for row in zip(*tb)))
            f.write("\n")

        tb=[["EMU"]]+[["InChI"]]
        if fdata:
            tb += [["value"]]+[["sd"]]
        for k in emu_i:
            nw=k.count("E")
            ili=range(nw+1)
            tb[0] += [k+"+"+str(i) for i in ili] # EEE+0 etc.
            a="/a(C"
            atoms=","+",".join(str(i+1) for i,sy in enumerate(k) if sy == "E")
            tb[1] += [a+str(nw if i == 0 else i)+"+"+str((not not i)+0)+atoms+")" for i in ili] # /a(C3+0,1,2,3) etc
            if fdata:
                tb[2] += [str(emuls[k][i]) if i in idef_e[k] else "" for i in ili]
                tb[3] += [str(emusd[k][i]) if i in idef_e[k] else "" for i in ili]
        res["inchi"]["emu"]=pa.DataFrame(tb[1][1:], index=tb[0][1:], columns=["inchi"])
        if fdata:
            res["inchi"]["emu"]["value"]=pa.to_numeric(tb[2][1:])
            res["inchi"]["emu"]["sd"]=pa.to_numeric(tb[3][1:])
        with open(nm_md+"_inchi_emu.tsv", "w") as f:
            f.write("\n".join("\t".join(row) for row in zip(*tb)))
            f.write("\n")
        if len(idef) or zcomb:
            tb=[["Accumulated combination"]+[nm_bc[i] for i in idef]+[i for li in compact for i in li]]
            tb += [["Group id"]+["i"+str(i+1) for i in range(len(idef))]+["c"+str(i+1) for i in range(len(compact)) for j in range(len(compact[i]))]]
            tb += [["InChI"]+["/i"+",".join(str(k+1)+"+"+sy for k,sy in enumerate(i) if sy in "01") for i in tb[0][1:]]]
            res["inchi"]["mcomb"]=pa.DataFrame(tb[0][1:], columns=["accu"])
            res["inchi"]["mcomb"]["group"]=tb[1][1:]
            res["inchi"]["mcomb"]["inchi"]=tb[2][1:]
            if fdata:
                tb += [["value"]+[str(v) for v in mecls[:len(idef)]]+[str(v) for v,li in zip(mecls[len(idef):], compact) for i in range(len(li))]]
                tb += [["sd"]+[str(v) for v in mecsd[:len(idef)]]+[str(v) for v,li in zip(mecsd[len(idef):], compact) for i in range(len(li))]]
                res["inchi"]["mcomb"]["value"]=pa.to_numeric(tb[3][1:])
                res["inchi"]["mcomb"]["sd"]=pa.to_numeric(tb[4][1:])
            with open(nm_md+"_inchi_mcomb.tsv", "w") as f:
                f.write("\n".join("\t".join(row) for row in zip(*tb)))
                f.write("\n")
        #pprint(res["inchi"])
    # md: finalize markdown/html
    if wri:
        mdFile.new_table_of_contents(table_title='Contents', depth=2)
        mdFile.create_md_file()
        md.markdownFromFile(input=nm_md+".md", output=nm_md+".html", extensions=['markdown.extensions.tables', 'extra', 'smarty'], output_format='html5')
        with open(nm_md+".html", "r") as f:
            htm=f.read()
        doc=TEMPLATE.replace('{{content}}', htm);
        doc=doc.replace('{{version}}', ver);
        now=strftime("%a, %d %b %Y %H:%M:%S", localtime())
        doc=doc.replace('{{date}}', now);
        doc=doc.replace('{{duration}}', dhms(tproc()-_T0));
        with open(nm_md+".html", "w") as f:
            print(doc, file=f)
        if __name__ == "__main__" or __name__[:8] == "isosolve" :
            try:
                import webbrowser
                webbrowser.open("file://"+os.path.realpath(f.name))
            except:
                pass
    
    res.update({"smm": smm, "mmd": mmd, "uvals": uvals, "sy_meas": sy_meas, "sv": sv,\
        "sbc": sbc, "sols": sols, "vsubs": vsubs, "metrics": {"idef": idef}})
    res.update({"emusol": emusol, "cusol": cusol})
    res["metrics"].update({"idef_c": idef_c, "idef_e": idef_e})

    if fdata:
        res["ls"]={"iso": [bcls, bcsd], "cov": covbc}
        res["ls"].update({"cumo": [cumols, cumosd], "emu": [emuls, emusd]})
        if not fast and (len(idef) or zcomb):
            res["ls"]["mec"]=[mecls, mecsd]

    timeme("end")
    return res
if __name__ == "__main__":
    main(li=sys.argv[1:])
