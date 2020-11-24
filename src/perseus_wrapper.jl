
##########################################################################################

####	PERSEUS WRAPPER

##########################################################################################

#=
This file is a simple Julia wrapper for the Perseus library for persistent
homology, which is produced and maintained by Vidit Nanda:

    http://people.maths.ox.ac.uk/nanda/perseus/

Note that keyword argument <filepath> should be a filepath to an executable
form of Perseus.

Important keyword argument/values are

model
- 	"vr" (input is a distance matrix)
- 	"pc2vr" (input is a point cloud; Eirene will compute and pass a distance matrix to perseus)
- 	"pc" (eirene will pass a pointcloud, calling the brips option)

=#

function perseusjl(
					datum;					# 	filepaths should end with .txt
					model					= 	"vr",
					rowsare 				= 	"dimensions",
					datapath				= 	testfp("prsip"),
					outpath					= 	testfp("prsop"),
					maxdim 					= 	1,
					minrad					= 	0,
					stepsz					= 	0.1,
					nsteps					= 	10,
					fr						= 	false, # stands for full resolution; used only for "vr" model
					scalefactor				= 	1,
					pointbirths				= 	[],
					perseusfilepath 		= 	"/Users/gh10/Google Drive/gregtaylorgoogledrive_nongooglefiles/GregDirectory/julia_gd/gdc_agora/gdc_a_peresuswrappers/perseusMac"
					)

	if in(model, ["perseusdistmat","perseusbrips"])
		datapath = datum
	else
		writelog 			=
		writeperseusfile(	datum;					# 	filepaths should end with .txt
							model					= 	model,
							rowsare 				= 	rowsare,
							datapath				= 	datapath,
							outpath					= 	outpath,
							maxdim 					= 	maxdim,
							minrad					= 	minrad,
							stepsz					= 	stepsz,
							nsteps					= 	nsteps,
							fr						= 	fr, # stands for full resolution
							scalefactor				= 	scalefactor,
							pointbirths				= 	pointbirths)
	end

	nsteps 	= 	writelog["nsteps"]
	ocg2rad =	writelog["ocg2rad"] # this will only matter when fr == true

	if 		in(model,["vr","perseusdistmat"])
		command  	= `$perseusfilepath distmat $datapath $outpath`
	elseif 	model 	== 	"pc"
		command  	= `$perseusfilepath rips $datapath $outpath`
	end
	run(`$(command)`)

	D = Dict(
			:barcodes 			=> 	Array{Array{Float64,2},1}(undef,maxdim+1),
			:betti   			=> 	Array{Int64}(undef,0,maxdim+1),
			:model 				=>	model,
			:maxdim 			=> 	maxdim,
			:perseusjlversion 	=>  "0.0.0"
			)

	# si stands for shifted index; since peresus ouptput starts indexing at 0 and arrays are 1-indexed, the true barcode of dimension r is D[:filtvalssi][D[:barcodes][r]+1]
	if fr
		D[:filtvalssi] 	=	reverse(ocg2rad,dims=1)
	else
		D[:filtvalssi] 	=	minrad:stepsz:(minrad+(1+nsteps)*stepsz)
	end

	g = open(outpath*"_betti.txt")  	# i check that the file is not empty in order to avoid throwing an error with readdlm
	seekend(g)

	#### the following commented code may be used to encode betti statistics;
	#### however, it may require modification
	#
	# if position(g) != 0
	# 	D[:betti]	= 	readdlm(outpath*"_betti.txt")
	# else
	# 	D[:betti] 	=	Array{Int64,2}(0,maxdim+1)
	# end

	for p = 1:maxdim+1
		g = open(outpath*"_$(p-1).txt")  	# i check that the file is not empty in order to avoid throwing an error with readdlm
		seekend(g)
		if position(g) != 0
			D[:barcodes][p] = readdlm(outpath*"_$(p-1).txt")
		else
			D[:barcodes][p] = Array{Float64,2}(0,2)
		end
	end

	return D
end

function writeperseusfile(
					datum;					# 	filepaths should end with .txt
					model					= 	"vr",
					rowsare 				= 	"dimensions",
					datapath				= 	joinpath(@__DIR__,"perseusreservoir/perseusreservoir_input.txt"),
					outpath					= 	joinpath(@__DIR__,"perseusreservoir/perseusreservoir_output"),
					maxdim 					= 	1,
					minrad					= 	0,
					stepsz					= 	0.1,
					nsteps					= 	10,
					fr						= 	false, # stands for full resolution
					scalefactor				= 	1,
					pointbirths				= 	[],
					)

	if typeof(datum)==String
		filename = datum
		if typeof(readdlm(filename,','))<:Array{Float64}
			s = readdlm(filename,',')
		elseif typeof(readdlm(filename,' '))<:Array{Float64}
			s = readdlm(filename,' ')
		else
			print("Error reading text file.  datum files of this format must be either comma or space delimited.")
			return
		end
	elseif typeof(datum)<:Array{Int64} || typeof(datum)<:Array{Float64}
		s = datum
	else
		print("Error: datum must be a string or an array of type Array{Float64} or Array{Int64}")
	end

	# the value assigned to ocg2rad only matters in specific cases; here
	# we assign it an arbitrary default value
	ocg2rad 	= 	[]

	if in(model,["vr","pc2vr"])

		if 	model 	== 	"pc2vr"
			if rowsare == "dimensions"
				d = Distances.pairwise(Euclidean(),s,dims=2)
			elseif rowsare == "points"
				d = Distances.pairwise(Euclidean(),s',dims=2)
			end
		end

		if model == "vr"
			d = convert(Array{Float64},s)
			if !issymmetric(d)
				println()
				println("Error: when the \"vr\" keyword value is passed the input array should be symmetric.")
				return
			end
		end

		if fr
			d,ocg2rad 	= 	ordercanonicalform_4(
							d;
							maxrad=Inf,
							minrad=-Inf,
							numrad=Inf,
							vscale="diagonal",
							)
			d 			= 	maximum(d)-d
			nsteps		= 	Inf
			stepsz		= 	1
			minrad		= 	0
			maxrad		= 	Inf
		end

		if nsteps == Inf
			nsteps = 1 + ceil(Int64,maximum(d)/stepsz)
		end

		samplesize = size(d,1)

		datastream = open(datapath,"w+")
		close(datastream)  				# this clears the current data file
		datastream = open(datapath,"a+")
		str = "$(samplesize)\n$(minrad) $(stepsz) $(nsteps) $(maxdim)"
		write(datastream,str)
		for p = 1:samplesize
			str = string(d[p,:])
			str = replace(str,"["=>"")
			str = replace(str,"]"=>"")
			str = replace(str,","=>"")
			str = "\n"*str
			write(datastream,str)
		end
		close(datastream)
	elseif model 	== 	"brips"
		println()
		println("NB: functionality for the \"brips\" keyword value is currently in beta.")
		if rowsare == "dimensions"
			s 	= 	s';
		end
		ambdim,numpts 		= 	size(s)
		if nsteps == Inf
			nsteps = 1 + ceil(Int64,maximum(Distances.pairwise(Euclidean(),s',dims=2))/stepsz)
		end
		if 	isempty(pointbirths)
			pointbirths 	= 	zeros(numpts)
		else
			pointbirths 	= 	pointbirths(:)
		end
		d 					= 	[s pointbirths]

		datastream = open(datapath,"w+")
		close(datastream)  				# this clears the current data file
		datastream = open(datapath,"a+")
		str = "$(ambdim)\n$(scalefactor) $(stepsz) $(nsteps)"
		write(datastream,str)
		for p = 1:samplesize
			str = string(d[p,:])
			str = replace(str,"["=>"")
			str = replace(str,"]"=>"")
			str = replace(str,","=>"")
			str = "\n"*str
			write(datastream,str)
		end
		close(datastream)
	end
	writelog 	= 	Dict("nsteps" => nsteps, "ocg2rad" => ocg2rad)
	return writelog
end

function barcode_perseus(D;dim=1)
	p = dim+1
	if !haskey(D,:barcodes)
		print("Error: input dictionary does not appear to contain barcode data.")
		return
	elseif p > D[:maxdim]
		print("Error: input dictionary does not appear to contain barcode data in dimensions greater than $(D[:maxdim]).  The user requested barcodes in dimensions up through $(p).")
		return
	else
		B 					=	D[:barcodes][p]
		B					=	2 .+round.(Int64,B)
		translator 			=	Array{Float64,1}(undef,1+length(D[:filtvalssi]))
		translator[2:end]	=	D[:filtvalssi]
		translator[1]		=	Inf
		return					translator[B]
	end
end

function saveperseustestdata()
	E 					= 	generateperseusvrdata()
	filepath			= 	testfp("prsjd")
	JLD.save(filepath,"E",E)
end


function generateperseusvrdata()
	numits 				= 	1
	ambdim 				= 	40
	maxdim 				= 	2
	numsteps 			= 	10000
	minrad 				= 	0
	stepsize 			= 	1
	calibrationdata 	= 	Array{Any}(undef,numits)

	for p 	= 	1:numits
		numpts 			= 	rand(50:70,1)
		numpts 			= 	numpts[1]
		vrmat 			= 	vertexlifemat(	numpts,
											model 		= 	"rand",
											scale 		= 	0) # we want zeros on the diagonal
		vrmat 			= 	ceil2grid(		vrmat*numsteps;
											origin=minrad,
											stepsize=stepsize,
											numsteps=numsteps)

		# birthtimes 		= 	birthtimesraw./2

		D 				=
		perseusjl(
		vrmat;			# 	filepaths should end with .txt
		model			= 	"vr",
		# rowsare 		= 	"dimensions",
		datapath		= 	testfp("prsip"),
		outpath			= 	testfp("prsop"),
		maxdim 			= 	maxdim,
		minrad			= 	minrad,
		stepsz			= 	stepsize,
		nsteps			= 	numsteps,
		# pointbirths		= 	birthtimes,
		perseusfilepath = 	"/Users/gh10/a/c/j/gdc_agora/gdc_a_peresuswrappers/perseusMac"
		)

		E 				= 	Dict(
							"perseusdict" 	=>	D,
							# "pcloud" 		=>	pcloud,
							"minrad"		=> 	minrad,
							"stepsize"		=> 	stepsize,
							"maxdim" 		=>  maxdim,
							"vrmat"			=>  vrmat)

		calibrationdata[p] 					= 	E
	end
	return calibrationdata
end


function ordercanonicalform_4(
	d;
	minrad=-Inf,
	maxrad=Inf,
	numrad=Inf, # note we have already performed the necessary rounding by this point
	fastop=true,
	vscale="diagonal",
	verbose = verbose)

	# round as necessary
	if 	minrad 		== 	"minedge"
		minrad 		= 	minimum(offdiagmin(d))
	end
	d 				= 	minmaxceil(d,minrad=minrad,maxrad=maxrad,numrad=numrad)
	d[d.>maxrad] 	= 	maxrad+1; # we make this reassignment b/c ordercanonicalform_3 takes only arguments with finite entries

	(t,ocg2rad) = ordercanonicalform_3(
		d;
		minrad=minrad,
		maxrad=maxrad,
		numrad=Inf, # note we have already performed the necessary rounding by this point
		fastop=fastop,
		vscale=vscale,
		verbose = verbose)
	return t, ocg2rad
end
