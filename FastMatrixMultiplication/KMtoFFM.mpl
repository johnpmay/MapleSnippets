# This Maple script verifies the GF2 Matrix multiplication formulas of
# Kauers and Moosbauer https://arxiv.org/pdf/2210.04045.pdf

# use these for local files
#file := "s47.exp"; n:=4;
#file := "s95.exp"; n := 5;
#filestring := FileTools:-Text:-ReadFile(file);

#filestring := URL:-Get("http://www.algebra.uni-linz.ac.at/people/mkauers/matrix-mult/s47.exp"); n := 4;
filestring := URL:-Get("http://www.algebra.uni-linz.ac.at/people/mkauers/matrix-mult/s95.exp"); n := 5;
if type(filestring, Array) then
    filestring := StringTools:-FromByteArray(filestring);
end if;
dims := [n,n,n];

# code to parse the exp file into the format used by
# the FMM catalog https://fmm.univ-lille.fr/

prods := map(parse,StringTools:-Split(StringTools:-Trim(filestring), "\n"));

prodtotriad := proc(expr::`*`, dims)
local d,tmp, a,b,c;
    (d, tmp) := selectremove(type, expr, integer);
    if nops(tmp)<3 then
        error "expected 3 factors in %1", tmp;
     end if;
    (a,b,c) := op(tmp);
    return [d*sumtomatrix(a,dims[1],dims[2]), sumtomatrix(b,dims[2],dims[3]), sumtomatrix(c,dims[3],dims[1])];
end proc:

sumtomatrix := proc(expr::{`+`,name},n,m)
local am, s, a, t, i, j;
    am := Matrix(n,m):
    for s in [op(expr)] do
        a := coeffs(s);
        t := s/a;
        i := sscanf(substring(t,2),"%d")[];
        j := sscanf(substring(t,3),"%d")[];
        am[i,j] := a;
    end do:
    return am;
end proc:

Tensor := TriadSet([seq(Triad(prodtotriad(p, dims)),p in prods)]);

idx := cat(dims[1],",",dims[2],",",dims[3]);
rank := nops(prods);
A := Matrix(dims[1], dims[2], (i,j)->cat(`a__`,i,`,`,j));
B := Matrix(dims[2], dims[3], (i,j)->cat(`b__`,i,`,`,j));
C := Matrix(dims[1], dims[3], (i,j)->cat(`c__`,i,`,`,j));

chk := add(
        LinearAlgebra:-Trace(LinearAlgebra:-Transpose(op([1,i,1,1],Tensor)).A)
    *LinearAlgebra:-Trace(LinearAlgebra:-Transpose(op([1,i,1,2],Tensor)).B)
    *LinearAlgebra:-Transpose(op([1,i,1,3],Tensor)),
    i=1..rank):
print("tensor check");
print(idx=(map(expand,A.B - chk) mod 2));

tmp := seq(
    LinearAlgebra:-Trace(LinearAlgebra:-Transpose(op([1,i,1,1],Tensor)).A)
    *LinearAlgebra:-Trace(LinearAlgebra:-Transpose(op([1,i,1,2],Tensor)).B)
    *LinearAlgebra:-Transpose(op([1,i,1,3],Tensor)),
    i=1..rank):

MUL := [seq(
        cat(`m__`,i)=remove(type,convert(tmp[i],set)[-1],numeric),
    i=1..rank)]:

ADD := [seq(seq(
    cat(`c__`,i,`,`,j)=
        add(
            ifelse(tmp[k][i,j]=0, 0,
                select(type,tmp[k][i,j],numeric)*cat(`m__`,k)),
        k=1..rank),
        j=1..dims[3]),i=1..dims[1])]:

print(MUL);
print(ADD);

print(`*`=nops(MUL));
print(`+`=`+`(seq(nops(rhs(e))-1, e in ADD)) + `+`(seq(nops(op(1,rhs(e)))+nops(op(2,rhs(e)))-2, e in MUL)) );

print("formula check");
print(idx=(map(expand,subs(subs(MUL,ADD),C-A.B)) mod 2));
