with(Groebner):
with(StringTools):

eqs_dirs := FileTools[ListDirectory]("Models");

for dir in eqs_dirs do

    eqs_file := cat("Models/", dir, "/odes.txt"):

    eqs_raw := FileTools[Text][ReadFile](eqs_file);
    eval(parse(cat("eqs := ", eqs_raw, ";")));
    
    # Plugging random numbers
    state_vars := map(e -> op(lhs(e))[1], eqs);
    all_vars := foldl(`union`, op(map(e -> indets(rhs(e)), eqs)));
    params := all_vars minus {op(state_vars)};
    roll := rand(1..10):
    subs_params := map(p -> p = roll(), params);
    eqs := map(e -> lhs(e) = subs(subs_params, rhs(e)), eqs);
    
    # converting to exponent + coeff vectors
    exp_vectors := [];
    max_deg := 0;
    for e in eqs do
        poly := rhs(e);
        print(poly, state_vars);
        monomials := Groebner[Support](expand(poly), state_vars);
        print(monomials);
        res := [];
        for m in monomials do
            exps := map(x -> degree(m, x), state_vars);
            res := [op(res), exps]:
            max_deg := max(max_deg, add(exps));
        end do;
        exp_vectors := [op(exp_vectors), res]:
    end do:
    
    if max_deg > 2 then
        print("System ", dir, " of degree ", max_deg, " and dimension ", nops(eqs));
        fd := FileTools[Text][Open](cat("parsed/", dir, "_res.txt"), create):
        FileTools[Text][WriteString](fd, cat(convert(nops(state_vars), string), "\n")):
        # Loosing coefficients !
        for r in exp_vectors do
            str_monoms := map(p -> cat("(", Join(map(d -> convert(d, string), p), ", "), ")"), r);
            FileTools[Text][WriteString](fd, cat(Join(str_monoms, " + "), "\n"));
        end do;
        FileTools[Text][Close](fd);
    end if;
end do:
