npts = 1000
c1 = 0.01
c2 = 0.1

funcs = (:bp, :br, :hp, :lp)
types = (:bu, :bu, :bu, :bu)
cs = ((c1, c2), (c1, c2), c1, c1)
passes = (2, 1, 1, 1)
npoles = (4, 6, 5, 3)

open(pipeline(`sac`, stdin), "w") do f
    println(f, """
        echo on
        fg impulse npts $npts
        w impulse.sac
        """)
    for (func, typ, c, p, np) in zip(funcs, types, cs, passes, npoles)
        corners = join(c, "-")
        outfile = join(("impulse", func, typ, "c", corners, "npoles", np,
                        "passes", p), "_") * ".sac"
        println(f, """
            fg impulse npts $npts
            $func $typ c $(join(c, " ")) npoles $np passes $p
            w $outfile
            """)
    end
end
