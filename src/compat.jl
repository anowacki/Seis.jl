# Compatibility with older Julia versions

@static if VERSION < v"1.3"
    sincosd(x) = (sind(x), cosd(x))
end
