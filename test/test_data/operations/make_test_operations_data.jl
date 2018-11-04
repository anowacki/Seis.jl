# Differentiation
open(pipeline(`sac`, stdin), "w") do f
    println(f, "echo on")
    for (k, v) in (2=>"two", 3=>"three", 5=>"five")
        outfile = "seis_diff_points_$k.sac"
        println(f, """
            fg seis
            dif $v
            w $outfile
            * Workaround for a bug in dif where npts is decreased in the operation, but the
            * full 1000 points are written to disk, giving a SAC file which is self-inconsistent.
            r $outfile
            w over
            """)
    end
end

# Integration
open(pipeline(`sac`, stdin), "w") do f
    println(f, "echo on")
    # N.B. Need to do trapezium first due to memory bug in MacSAC
    for (k, v) in (:trapezium=>"trapezoidal", :rectangle=>"rectangular")
        outfile = "seis_int_$k.sac"
        println(f, """
            fg seis
            rmean
            int $v
            w $outfile
            """)
    end
    println(f, "quit")
end
