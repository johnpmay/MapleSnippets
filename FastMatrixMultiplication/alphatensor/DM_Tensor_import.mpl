#filename := "factorizations_f2.npz";
#filename := "factorizations_r.npz";
filename := "factorizations_f2.mpl";
#filename := "factorizations_r.mpl";

# use Maple's python to unpack the npz
if filename[-3..-1] = "npz" then
    Python:-Import("numpy");
    import := Python:-EvalString(cat("dict( numpy.load('", filename, "',allow_pickle=True))")):
    factorizations := convert(import, table):
else
    read filename;
end if;

for idx in sort([indices(factorizations)]) do
    print(idx);
    if filename[-3..-1] = "npz" then
        tmp := convert(factorizations[idx[]], Array);
        if type(tmp[1], python) then
            tmp := convert(tmp, list);
            tmp := map(convert, tmp, Array);
        else
            tmp := [tmp[1,...,...], tmp[2,...,...], tmp[3,...,...]];
        end if;
        factorizations[idx[]] := tmp;
    end if;

    dims := sscanf(idx[], "%d,%d,%d");

    A:=Matrix(dims[1], dims[2], (i,j)->cat(`a__`,i,`,`,j));
    B:=Matrix(dims[2], dims[3], (i,j)->cat(`b__`,i,`,`,j));
    C := Matrix(dims[1], dims[3], (i,j)->cat(`c__`,i,`,`,j));

    (u,v,w) := factorizations[idx[]][];
    rank := numelems(u[1]);
    print(idx[]=rank);

# this is the Tensor format used by FMM https://fmm.univ-lille.fr/
    Tensor := TriadSet([seq(Triad([
            Matrix(dims[1], dims[2], convert(convert(u[..,i],list),rational)),
            Matrix(dims[2], dims[3], convert(convert(v[..,i],list),rational)),
            Matrix(dims[3], dims[1], convert(convert(w[..,i],list),rational)) ])
        ,i=1..rank)]);

#check
    chk := add(
         LinearAlgebra:-Trace(LinearAlgebra:-Transpose(op([1,i,1,1],Tensor)).A)
        *LinearAlgebra:-Trace(LinearAlgebra:-Transpose(op([1,i,1,2],Tensor)).B)
        *LinearAlgebra:-Transpose(op([1,i,1,3],Tensor)),
        i=1..rank);

    if filename[-6..-5] = "f2" then
        print(idx=(map(expand,A.B - chk) mod 2));
    else
        print(idx=(map(expand,A.B - chk) ));
    end if;

# evaluation formulas

    tmp := seq(
        LinearAlgebra:-Trace(LinearAlgebra:-Transpose(op([1,i,1,1],Tensor)).A)
        *LinearAlgebra:-Trace(LinearAlgebra:-Transpose(op([1,i,1,2],Tensor)).B)
        *LinearAlgebra:-Transpose(op([1,i,1,3],Tensor)),
        i=1..rank):

    MUL := [seq(
            cat(`m__`,i)=remove(type,convert(tmp[i],set)[-1],numeric),
        i=1..rank)];

    ADD := [seq(seq(
        cat(`c__`,i,`,`,j)=
            add(
                ifelse(tmp[k][i,j]=0, 0,
                    select(type,tmp[k][i,j],numeric)*cat(`m__`,k)),
            k=1..rank),
            j=1..dims[3]),i=1..dims[1])];

    #print(MUL);
    #print(ADD);

    print(`*`=nops(MUL));
    print(`+` = `+`(seq(nops(rhs(e))-1, e in ADD))
              + `+`(seq(nops(op(1,rhs(e)))+nops(op(2,rhs(e)))-2, e in MUL)) );

#check
    if filename[-6..-5] = "f2" then
        print(idx=(map(expand,subs(subs(MUL,ADD),C-A.B)) mod 2));
    else
        print(idx=(map(expand,subs(subs(MUL,ADD),C-A.B))));
    end if;

end do:


if filename[-3..-1] = "npz" then
    save factorizations, cat(filename[1..-4],"mpl");
end if;
