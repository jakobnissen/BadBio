# BadBio
_Convenience code for bioinformatics, written badly in Julia_

BioJulia is a collection of well-designed, generic and fast packages well suited for laying the foundation of a solid, trustworthy computational science ecosystem.

This is not one of those packages.

Sometimes, you don't care if your code doesn't scale. You don't care about the future implications of your questionable API decisions. You just want stuff done *right now*, interactively. And often, you find yourself writing the same convenience functions over and over again. If only there was a way to save some time by not repeating yourself...

So do you write a package with this functionality? Of course not! Writing quality code takes time as you settle on the design and API over several iterations. Instead, you dump your crap code in a repo like this, so that others may copy it, and it doesn't get lost over time.

## Features:
* No tests, whatsoever.
* It's probably pretty slow, IDK.
* It's not even a package. Like, you need to use `include`.
* Blatantly ignores SemVer 2.0; will push breaking changes straight to master when you least expect it.
* Pirates other packages' types with abandon

## Contributing
Is this repo missing something? Make a PR, and i'll include it. Actually, I'll probably grant you write access to this repo. After all, this is all about making the least amount of effort.

## Content:
__Non-biological functionality__

* `@indir directory expression`: `cd` to `directory`, execute `expression`, then move back
* `nested_dict(K=Any, V=Any)`: Create an infinitely nested `DefaultDict{K,V}`, such that you can do e.g. `d[1]["myname"][2.1][false] = "foo"`
* `Path`: An object representing a file path
    * `Base.:*(::Path, ::Union{Path, String})`: Join the paths using joinpath
    * `split(::Path)`: Split along file delimiter

__Biological functionality__

* `NamedSeq`: A BioSequence with a name attached.
    * `read_fasta(::IO)`: Return a vector of `NamedSeq` from a FASTA-formatted IO
    * `print(::IO, ::NamedSeq)`: Print in FASTA format
* `mafft(seqs)` align a collection of `NamedSeq` in memory using `mafft` from command line, and returns a vector of `NamedSeq`.
* `findbest(gene::AminoAcidSeq, seq::NucleotideSeq)`: Find best-matching position of `gene` in the nucleotide seq `seq`
