using Printf
using Dates


macro flt()
    t = Dates.format(now(), "HH:MM:SS")
    return :( @sprintf("%s %s:%d |", $t, $(basename(string(__source__.file))), $(__source__.line)) )
end

@doc"""
    @infoe message
"""
macro infoe(message)
    fl = string(basename(string(__source__.file)), ":", (__source__.line))
    return :( printstyled(string("[ Info: ", $fl, " | ", $(esc(message)), "\n"), color = :light_cyan, bold = true) )
end

@doc """
     @warne message
 """
macro warne(message)
    fl = string(basename(string(__source__.file)), ":", (__source__.line))
    return :( printstyled(string("[ Warn: ", $fl, " | ", $(esc(message)), "\n"), color = :green, bold = true) )
end

@doc """
     @errore message
 """
macro errore(message)
    fl = string(basename(string(__source__.file)), ":", (__source__.line))
    return :( printstyled(string("[ Error: ", $fl, " | ", $(esc(message), "\n")), color = :red, bold = true) )
end

export @infoe, @warne, @errore


