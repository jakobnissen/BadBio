module BadSeq

using BioSequences
using FASTX
using BioAlignments
using DataStructures

######## NamedSeq ################
struct NamedSeq{T <: BioSequence}
    name::String
    seq::T
end

Base.length(s::NamedSeq) = length(s.seq)
Base.setindex!(x::NamedSeq, i::Integer, v) = (x.seq[i] = v)
Base.getindex(x::NamedSeq, i::Integer) = x.seq[i]

function Base.print(io::IO, s::NamedSeq)
    print(io, '>')
    println(io, s.name)
    print(io, s.seq)
end

function Base.read(reader::FASTA.Reader, ::Type{NamedSeq})
    record = read!(reader, FASTA.Record())
    header = FASTA.identifier(record)
    FASTA.hasdescription(record) && (header *= (" " * FASTA.description(record)))
    seq = FASTA.sequence(record)
    return NamedSeq(header, seq)
end

function read_fasta(io::IO)
    io = io isa FASTA.Reader ? io : FASTA.Reader(io)
    result = NamedSeq[]
    while !eof(io)
        push!(result, read(io, NamedSeq))
    end
    result
end

macro seq_str(str, name)
    :(NamedSeq($(name), LongDNASeq($(str))))
end

"MAFFT align sequences. By default prints a log to /tmp"
function mafft(seqs, stderr=tempname("/tmp"))
    inpath = tempname("/tmp")
    
    open(inpath, "w") do file
        for seq in seqs
            println(file, seq)
        end
    end
    out = IOBuffer()
    if stderr isa String
        open(stderr, "w") do file
            run(pipeline(pipeline(`mafft $inpath`; stderr=file), out))
        end
    else
        run(pipeline(pipeline(`mafft $inpath`; stderr=stderr), out))
    end
    rm(inpath)
    return read_fasta(seekstart(out))
end

"Trim '-' from both ends of an iterable of Seqs"
function trim(seqs)
    strings = [String(s.seq) for s in seqs]
    leftpad  = maximum([length(s) - length(lstrip(s, ['-'])) for s in strings])
    rightpad = maximum([length(s) - length(rstrip(s, ['-'])) for s in strings])
    return [NamedSeq(s.name, typeof(s.seq)(st[leftpad+1:end-rightpad])) for (s, st) in zip(seqs, strings)]
end
#########################

######### DNA - AA alignment ############

function findstarts(seq::LongSequence{DNAAlphabet{2}})
	starts = Int[]
	query = dna"ATG"
	p = 1
	loc = findfirst(query, seq, 1)
	while loc !== nothing
		push!(starts, first(loc))
		p = last(loc) + 1
		loc = findfirst(query, seq, p)
	end
	return starts
end

function orflength(seq::LongSequence{DNAAlphabet{2}}, start)
	it = each(DNACodon, seq[start:end], 3)
	stop = 0
	for (pos, fw, rv) in it
		stop = pos-1
		if (fw == mer"TAG") | (fw == mer"TGA") | (fw == mer"TAA")
			break
		end
	end
	return stop
end

function likely_orfs_forward(seq::LongSequence{DNAAlphabet{2}}, len=1)
	starts = findstarts(seq)
	likely_orfs = UnitRange{Int}[]
	cutoff = round(Int, 0.8 * len)
	for start in starts
		length(seq) - start < cutoff && break
		orflen = orflength(seq, start)
		orflen >= cutoff && push!(likely_orfs, start:start+orflen-1)
	end
	return likely_orfs
end

function best_alignment_forward(seq::LongSequence{DNAAlphabet{2}}, gene::AminoAcidSeq)
	likely_orfs = likely_orfs_forward(seq, 3*length(gene))
	isempty(likely_orfs) && return (1:0, typemin(Int))
	
	scoremodel = AffineGapScoreModel(BLOSUM62, gap_open=-5, gap_extend=-1)
	alignments = []
	
	for orf in likely_orfs
		orfgene = translate(seq[orf])
		push!(alignments, pairalign(GlobalAlignment(), orfgene, gene, scoremodel))
	end
	
	bestindex = argmax(map(score, alignments))
	return (likely_orfs[bestindex], score(alignments[bestindex]))
end

function findbest(gene::AminoAcidSeq, seq::LongSequence{<:NucleicAcidAlphabet})
	seq_ = convert(LongSequence{DNAAlphabet{2}}, seq)
	orf_fw, score_fw = best_alignment_forward(seq_, gene)
	orf_rv, score_rv = best_alignment_forward(reverse_complement(seq_), gene)

    reversed = score_rv > score_fw
    return reversed, ifelse(reversed, orf_rv, orf_fw)
end

###########
"Get all maximum-length orfs separated by a stop codon, if range satisfies f"
function get_orfs(seq::BioSequence{<:NucleicAcidAlphabet})
    orfs = UnitRange{Int}[]
    for frame in 1:3
        lastindex = length(seq) - (length(seq) - frame + 1) % 3
        aas = translate(seq[frame:lastindex])
        start = frame
        local stop = 0
        for i in eachindex(aas)
            index = 3*(i-1) + frame
            if aas[i] == AA_Term
                stop = index - 1
                stop > start && push!(orfs, start:stop)
                start = index + 3
            end
        end
        # If no stop was detected at end.
        start < lastindex && push!(orfs, start:lastindex)
    end
    return orfs
end

function get_met_orfs(seq::BioSequence{<:NucleicAcidAlphabet})
    orfs = UnitRange{Int}[]
    start = typemax(Int)
    stop = 0
    for frame in 1:3
        lastindex = length(seq) - (length(seq) - frame + 1) % 3
        aas = translate(seq[frame:lastindex])
        for i in eachindex(aas)
            index = 3*(i-1) + frame
            if aas[i] == AA_Term
                stop = index - 1
                stop > start && push!(orfs, start:stop)
                start = typemax(Int)
            elseif (start == typemax(Int)) & (aas[i] == AA_M)
                start = index
            end
        end
        stop > start && push!(orfs, start:stop)
    end
    return orfs
end

longest(xs) = maximum(((length(i), i) for i in xs))[2]

"Get (is_rc, largest ORF) in the sequence. `met` is whether it must start with Met"
function largest_orf(sequence::BioSequence{<:NucleicAcidAlphabet}, met=false)
    seq = copy(sequence)
    f = ifelse(met, get_met_orfs, get_orfs)
    fw = longest(f(seq))
    rv = longest(f(reverse_complement!(seq)))
    reversed = length(rv) > length(fw)
    return reversed, ifelse(reversed, rv, fw)
end

function largest_aa(sequence::BioSequence{<:NucleicAcidAlphabet}, met=false)
    seq = copy(sequence)
    isrv, orfrange = largest_orf(sequence, met)
    return translate((isrv ? reverse_complement(sequence) : sequence)[orfrange])
end

end # module
