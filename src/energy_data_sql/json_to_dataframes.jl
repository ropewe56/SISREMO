@inline function normalize_colname(colname)
    colname = replace(colname, " " => "_")
    colname = replace(colname, "/" => "")
    colname = replace(colname, "__" => "_")
    colname = replace(colname, "-" => "_")
    colname = replace(colname, "_(incl._self_consumption)" => "")
    colname = replace(colname, "(" => "")
    colname = replace(colname, ")" => "")
    colname = replace(colname, "," => "")
    colname
end

function normalize_keys(din)
    dout = OrderedDict()
    pk = collect(keys(din))
    for k in pk
        kn = normalize_colname(k)
        dout[kn] = din[k]
    end
    dout
end

function json_to_dataframe(jsonpath)

    cont = open(jsonpath, "r") do io
        JSON3.read(io)
    end

    tk = if haskey(cont, :unix_seconds)
        :unix_seconds
    elseif haskey(cont, :time)
        :time
    end
    
    mk = if haskey(cont, :production_types)
        :production_types
    elseif haskey(cont, :countries)
        :countries
    end

    uts = cont[tk]
    din = OrderedDict()
    for pt in cont[mk]
        data = @. ifelse(pt[:data] === nothing, 0.0, pt[:data])
        din[pt[:name]] = Float64.(data)
    end
    dout = normalize_keys(din)

    colnames = Vector{String}(undef,0)
    values   = []
    push!(colnames, string(tk))
    push!(values, uts)
    for (k,v) in dout
        push!(colnames, k)
        push!(values, v)
    end
    DataFrame(colnames .=> values)
end

 