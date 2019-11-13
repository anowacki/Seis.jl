using Seis, Seis.Plot, DSP, Plots
using PyCall

obspy = pyimport("obspy")

ndec = 3

fs = 5
delta = 1/fs
nyquist = fs/2/ndec
npts = round(Int, 1000/fs)
s0 = Trace{Float32}(0, delta, npts)
tim = times(s0)

s0.t .= zeros(Float32, npts)
s0.t[end÷2] = 1
s0.t .= rand(Float32, npts) .- 0.5
s0.t .= [sin(π*t/2) + sin(3π*t)/3 for t in times(s0)]

# Write data
sdec = deepcopy(s0)
write_sac(s0, "/tmp/imp_orig.sac")

# SAC decimation
open(pipeline(`sac`, stdout=devnull), "w", stdin) do f
    println(f, """
        r /tmp/imp_orig.sac
        decimate $ndec
        w /tmp/imp_dec.sac
        q
        """)
end
ssac = read_sac("/tmp/imp_dec.sac")

# DSP.resample decimation
sdec.t = resample(s0.t, 1//ndec)
sdec.delta *= ndec

# Obspy decimation
odec = get(obspy.read("/tmp/imp_orig.sac"), 0)
odec.decimate(ndec)
odec = begin
    t = deepcopy(sdec)
    t.t = odec.data
    t
end

# Compare results
labels = hcat("Orig", "SAC dec", "Seis dec", "Obspy dec")

p = [periodogram(ss.t, fs=1/ss.delta) for ss in [s0, ssac, sdec, odec]]

pl = plot(freq.(p), power.(p), label=labels, legend=:topright,
    lw=[4 3 2 2], la=[0.5 0.5 1 1])
vline!(pl, [nyquist], label="", lw=2, lc=:black)

plot(
    pl,
    plot([s0, ssac, sdec, odec], label=vec(labels), ylim=(-1,1)),
    layout=(2,1),
    size=(600,600),
    )