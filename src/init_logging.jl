using Printf, Dates
using Logging, LoggingExtras

const colour = Dict(Debug => :blue, Info => :light_cyan, Warn => :green, Error => :red)

logger = FormatLogger() do io, args
    printstyled(io, "[ ", args.level, ": ", basename(args.file), ":", args.line, " | ", args.message, "\n", color = colour[args.level], bold = true)
end;

global_logger(logger);

disable_logging(Debug)
#disable_logging(Logging.BelowMinLevel)

#=
macro flt()
	t = Dates.format(now(), "HH:MM:SS")
	return :( @sprintf("%s %s:%d |", $t, $(basename(string(__source__.file))), $(__source__.line)) )
end
macro info(message)
	fl = string(basename(string(__source__.file)), ":", (__source__.line))
	return :( printstyled(string("[ Info: ", $fl, " | ", $(esc(message)), "\n"), color = :light_cyan, bold = true) )
end
macro warn(message)
	fl = string(basename(string(__source__.file)), ":", (__source__.line))
	return :( printstyled(string("[ Warn: ", $fl, " | ", $(esc(message)), "\n"), color = :green, bold = true) )
end
macro error(message)
	fl = string(basename(string(__source__.file)), ":", (__source__.line))
	return :( printstyled(string("[ Error: ", $fl, " | ", $(esc(message), "\n")), color = :red, bold = true) )
end
=#

nothing
