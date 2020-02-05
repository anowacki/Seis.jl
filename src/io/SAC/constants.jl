# Contains constants and constructors for SAC file precision, and defaults
# for routines and endianness

# SAC types
const SACFloat = Float32
const SACChar = String
const SACInt = Int32
const SACBool = Bool
# Constructors
sacfloat(x) = map(Float32, x)
sacint(x) = map(Int32, x)
const sacchar = ascii
sacbool(x) = x != 0
# Length of SAC floats and ints
const SAC_BYTE_LEN = 4
# Length of SAC character headers (except kevnm, which is twice the length)
const SACCHARLEN = 8
# SAC file version number
const SAC_VER_NUM = SACInt(6)
# Whether this machine is big- or little endian.  SAC files are meant to be big-endian,
# so this determines whether a file is 'native-endian' or not.
const MACHINE_IS_LITTLE_ENDIAN = Base.ENDIAN_BOM == 0x04030201

# Convert a number into a SACChar
sacstring(x, maxlen=SACCHARLEN) = sacchar(string(x)[1:minimum((length(string(x)),maxlen))]*" "^(maximum((0,maxlen-length(string(x))))))

# SAC unset values
const SAC_RNULL = -12345.
const SAC_INULL = -12345
const SAC_CNULL = "-12345"

# Default values for filtering
const SAC_NPOLES = 2
const SAC_PASSES = 1

# For SAC/BRIS (MacSAC), files are always big-endian, so set this appropriately
const SAC_FORCE_SWAP = MACHINE_IS_LITTLE_ENDIAN

# Flag for verbosity
const SAC_VERBOSE = Ref(true)

# Lists of SAC headers as symbols
const SAC_FLOAT_HDR = [:delta, :depmin, :depmax, :scale, :odelta, :b,
    :e, :o, :a, :internal0, :t0, :t1,
    :t2, :t3, :t4, :t5, :t6, :t7,
    :t8, :t9, :f, :resp0, :resp1, :resp2,
    :resp3, :resp4, :resp5, :resp6, :resp7, :resp8,
    :resp9, :stla, :stlo, :stel, :stdp, :evla,
    :evlo, :evel, :evdp, :mag, :user0, :user1,
    :user2, :user3, :user4, :user5, :user6, :user7,
    :user8, :user9, :dist, :az, :baz, :gcarc,
    :internal1, :internal2, :depmen, :cmpaz, :cmpinc, :xminimum,
    :xmaximum, :yminimum, :ymaximum, :unused1, :unused2, :unused3,
    :unused4, :unused5, :unused6, :unused7]
const SAC_INT_HDR = [:nzyear, :nzjday, :nzhour, :nzmin, :nzsec, :nzmsec,
    :nvhdr, :norid, :nevid, :npts, :internal3, :nwfid,
    :nxsize, :nysize, :unused8, :iftype, :idep, :iztype,
    :unused9, :iinst, :istreg, :ievreg, :ievtyp, :iqual,
    :isynth, :imagtyp, :imagsrc, :unused10, :unused11, :unused12,
    :unused13, :unused14, :unused15, :unused16, :unused17]
const SAC_BOOL_HDR = [:leven, :lpspol, :lovrok, :lcalda, :unused18]
const SAC_CHAR_HDR = [:kstnm, :kevnm, :khole, :ko, :ka, :kt0,
    :kt1, :kt2, :kt3, :kt4, :kt5, :kt6,
    :kt7, :kt8, :kt9, :kf, :kuser0, :kuser1,
    :kuser2, :kcmpnm, :knetwk, :kdatrd, :kinst]
const SAC_ALL_HDR = [SAC_FLOAT_HDR; SAC_INT_HDR; SAC_BOOL_HDR; SAC_CHAR_HDR]

# Where in the file the NVHDR value is
const SAC_NVHDR_POS = length(SAC_FLOAT_HDR) + findfirst(SAC_INT_HDR .== :nvhdr) - 1

# Length in bytes of total SAC header, accounting for double-length kevnm
const SAC_HEADER_LEN = SAC_BYTE_LEN*(length(SAC_FLOAT_HDR) + length(SAC_INT_HDR) +
    length(SAC_BOOL_HDR)) + SACCHARLEN*(length(SAC_CHAR_HDR) + 1)

# Some of the enumerated header values (see https://ds.iris.edu/files/sac-manual/manual/file_format.html)
"Time series file"
const SAC_ITIME = 1
"Spectral file---real and imaginary"
const SAC_IRLIM = 2
"Spectral file---amplitude and phase"
const SAC_IAMPH = 3
"General x versus y data"
const SAC_IXY = 4
"Unknown"
const SAC_IUNKN = 5
