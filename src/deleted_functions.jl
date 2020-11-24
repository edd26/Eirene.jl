
function handcalcsolution()

	K 			= 	Dict()

    # SPHERE / MACHINE PRECISION

	maxdim 		= 	2
	model 		= 	"complex"
	maxrad 		= 	100
	numrad 		= 	Inf
	space		=   "sphere"
	solkey 		= 	solutionkey(
					model 		= 	"complex",
					maxrad 		= 	maxrad,
					numrad		= 	numrad,
					space		= 	space,
					problemset 	= 	"hand")

	K[solkey]					=	Dict()
	K[solkey][:barcodes]		= 	Array{Any,1}(undef,maxdim+1)
	K[solkey][:cyclerep]		= 	Array{Any,1}(undef,maxdim+1)  # this will be filled in manually

	K[solkey][:barcodes][1] 	=	[20.0 30.0;10.0 Inf]
	K[solkey][:barcodes][2] 	=	[40.0 50.0]
	K[solkey][:barcodes][3] 	=	[60.0 70.0]

	K[solkey][:cyclerep][1]	= 	[1; 2]
	K[solkey][:cyclerep][2]	= 	[1; 2]
	K[solkey][:cyclerep][3]	= 	[1; 2]

    # SPHERE / 10:100 PRECISION
    # (identical to sphere with machine precision, just a different nli & key value)

	model 		= 	"complex"
	maxrad 		= 	100
	numrad 		= 	10
	space		=   "sphere"
	solkey 		= 	solutionkey(
					model 		= 	"complex",
					maxrad 		= 	maxrad,
					numrad		= 	numrad,
					space		= 	space,
					problemset 	= 	"hand")

	K[solkey]					=	Dict()
	K[solkey][:barcodes]		= 	Array{Any,1}(undef,maxdim+1)
	K[solkey][:cyclerep]		= 	Array{Any,1}(undef,maxdim+1)  # this will be filled in manually

	K[solkey][:barcodes][1] 	=	[20.0 30.0;10.0 Inf]
	K[solkey][:barcodes][2] 	=	[40.0 50.0]
	K[solkey][:barcodes][3] 	=	[60.0 70.0]

	K[solkey][:cyclerep][1]	= 	[1; 2]
	K[solkey][:cyclerep][2]	= 	[1; 2]
	K[solkey][:cyclerep][3]	= 	[1; 2]

    # EMPTY SPACE / MACHINE PRECISION

	model 		= 	"complex"
	maxrad 		= 	100
	numrad 		= 	Inf
	space		=   "empty"
	solkey 		= 	solutionkey(
					model 		= 	"complex",
					maxrad 		= 	maxrad,
					numrad		= 	numrad,
					space		= 	space,
					problemset 	= 	"hand")

	K[solkey]					=	Dict()
	K[solkey][:barcodes]		= 	Array{Any,1}(undef,maxdim+1)
	K[solkey][:cyclerep]		= 	Array{Any,1}(undef,maxdim+1)  # this will be filled in manually

	for p = 1:maxdim+1
		K[solkey][:barcodes][p] 	=	Array{Int64,2}(undef,0,2)
	end

    # EMPTY SPACE / 10:100 PRECISION

	model 		= 	"complex"
	maxrad 		= 	100
	numrad 		= 	10
	space		=   "empty"
	solkey 		= 	solutionkey(
					model 		= 	"complex",
					maxrad 		= 	maxrad,
					numrad		= 	numrad,
					space		= 	space,
					problemset 	= 	"hand")

	K[solkey]					=	Dict()
	K[solkey][:barcodes]		= 	Array{Any,1}(undef,maxdim+1)
	K[solkey][:cyclerep]		= 	Array{Any,1}(undef,maxdim+1)  # this will be filled in manually

	for p = 1:maxdim+1
		K[solkey][:barcodes][p] 	=	Array{Int64,2}(undef,0,2)
	end

	# SKRABA TRIANGLE / MACHINE PRECISION

	maxdim 		= 	2
	model 		= 	"complex"
	maxrad 		= 	100
	numrad 		= 	Inf
	space		=   "skrabatriangle"
	solkey 		= 	solutionkey(
					model 		= 	"complex",
					maxrad 		= 	maxrad,
					numrad		= 	numrad,
					space		= 	space,
					problemset 	= 	"hand")

	K[solkey]					=	Dict()
	K[solkey][:barcodes]		= 	Array{Any,1}(undef,maxdim+1)
	K[solkey][:cyclerep]		= 	Array{Any,1}(undef,maxdim+1)

	K[solkey][:barcodes][1] 	=	[1.0 Inf]
	K[solkey][:barcodes][2] 	=	[3.0 7.0]
	K[solkey][:barcodes][3] 	=	zeros(Float64,0,2)

	K[solkey][:cyclerep][1]		= 	[1]
	K[solkey][:cyclerep][2]		= 	[1,2,3]
	K[solkey][:cyclerep][3]		= 	zeros(Float64,0,2)

	# SKRABA TRIANGLE / DUMMY PRECISION

	maxdim 		= 	2
	model 		= 	"complex"
	maxrad 		= 	100
	numrad 		= 	10
	space		=   "skrabatriangle"
	solkey 		= 	solutionkey(
					model 		= 	"complex",
					maxrad 		= 	maxrad,
					numrad		= 	numrad,
					space		= 	space,
					problemset 	= 	"hand")

	K[solkey]					=	Dict()
	K[solkey][:barcodes]		= 	Array{Any,1}(undef,maxdim+1)
	K[solkey][:cyclerep]		= 	Array{Any,1}(undef,maxdim+1)

	K[solkey][:barcodes][1] 	=	[1.0 Inf]
	K[solkey][:barcodes][2] 	=	[3.0 7.0]
	K[solkey][:barcodes][3] 	=	zeros(Float64,0,2)

	K[solkey][:cyclerep][1]		= 	[1]
	K[solkey][:cyclerep][2]		= 	[1,2,3]
	K[solkey][:cyclerep][3]		= 	zeros(Float64,0,2)

	return K
end




function savesolutions()
	print(
	"""

	====================================================

	Please note: the file

	/Users/gh10/Google Drive/gregtaylorgoogledrive_nongooglefiles/GregDirectory/julia_gd/gdc_agora/gdc_a_peresuswrappers/gdc_a_a_perseuswrapper.jl

	must be loaded for this function to work properly.

	====================================================
	""")
	K 	= 	generatecrosscheckdata_perseus()
	manualcycleadditions(K)
	filepath = joinpath(@__DIR__,"testsolutions/testsolutions.jld")
	JLD.save(filepath,"K",K)
	# JLD.save("/Users/gh10/Google Drive/gregtaylorgoogledrive_nongooglefiles/GregDirectory/julia_gd/gdb_eirene/gdb_data/gdb_da_a_calib/gdb_da_a_d_solns/gdb_da_a_dsoln_a.jld","K",K)
	# JLD.save("/Users/greghenselman/Google Drive/GregDirectory/julia_gd/gdb_eirene/gdb_data/gdb_da_a_calib/gdb_da_a_d_solns/gdb_da_a_dsoln_a.jld","K",K)
	return K
end



# 12/28/2017
# This function is meant to generate crosscheck data for the version of <unittest>
# defined 12/30/2017.
function generatecrosscheckdata_perseus()

	K 			= 	Dict()			# K for "calibrate" (letter C was taken)

	for entryformat 	= ["textfile"]
		for model 		= ["vr" "pc"] # cellular examples must be handled separately
			for maxdim 		= [0 1 2]
				for minrad 		= [-Inf 0]
					for maxrad 		= [Inf, 100]
						for NUMRAD 		= [1 10 0]
							for fastop 		= [true,false]
								for vscale 		= [[]]
									for record 		= ["all" "cyclerep" "none"]
										for pointlabels	= [[]]
											for iteration	= [1 2]

												# the last three of these are arbitrary choices
												proceed 		= 	maxdim 	== 	2 			&&
																	model 	== 	"vr"		&&
																	record 	== 	"all"		&&
																	minrad 	== 	0			&&
																	fastop
												if proceed

													filepath 	= 	modit2filepath(model,iteration-1)
													numrad 		= 	NUMRAD2numrad(NUMRAD)
													solutionkey(	model 		= 	"complex",
																	maxrad 		= 	Inf,
																	numrad		= 	Inf,
																	space		= 	0,
																	problemset 	= 	"checkparameters")
													stepsz 		= 	nr2ss(numrad)

													E	= 	perseusjl(
															filepath;				# 	filepaths should end with .txt
															model					= 	model,
															rowsare 				= 	"distances",
															maxdim 					= 	maxdim,
															minrad					= 	0,
															stepsz					= 	stepsz,
															nsteps					= 	Inf,
															fr 						= 	numrad == Inf)

													K[solkey]						=	Dict()
													K[solkey][:barcodes]			= 	Array{Any,1}(undef,maxdim+1)
													K[solkey][:cyclerep]			= 	Array{Any,1}(undef,maxdim+1)  # this will be filled in manually
													for r 	= 	1:(maxdim+1)
														K[solkey][:barcodes][r] = barcode_perseus(E,dim=r-1)
													end
												end
											end
										end
									end
								end
							end
						end
					end
				end
			end
		end
	end
	return K
end

function nr2ss(numrad)
	if 		numrad 	== 	1
		ss	= 	100
	elseif numrad 	== 	10
		ss 	= 	10
	elseif 	numrad 	== 	Inf
		ss 	= 	1
	else
		println()
		println("error: <nr2ss> only takes arguments 1, 10, and 0")
		println("argument passed was: $(numrad)")
	end
	return ss
end

function NUMRAD2numrad(NUMRAD)
	if NUMRAD == 0
		numrad 	= 	Inf
	else
		numrad 	= 	NUMRAD
	end
	return numrad
end



function modit2filepath(model,iteration)
	suffix 	= 	string("testdata/",model,"/",model,"$(iteration)_input.csv")
	joinpath(@__DIR__,suffix)
end


function diagonalentries(x)
	if 	size(x,1) != size(x,2)
		println()
		println("error: d should be square")
		return
	end
	m 			= 	size(x,2)
	v 			= 	zeros(m)
	for 	p 	= 	1:m
		v[p]	= 	x[p,p]
	end
	return 		v
end

function ceil2grid_overflow( 	N;
								minrad 		=	-Inf,
								maxrad		=	Inf,
								mingrid 	= 	Inf,
								maxgrid 	= 	Inf,
								numgrid 	= 	Inf
								)
	S 							= 	truncatearray(N,minrad,maxrad)
	fi 							= 	findall(isfinite,S) # stands for finite indices
	fv 							= 	S[fi]
	ran 						= 	makegrid(mingrid,maxgrid,numgrid)
	S[fi] 						= 	ceil2grid(fv,ran)
	return S
end

function makegrid(mingrid,maxgrid,numrad)
	if numrad < 1
		println("error in <makegrid>: argument numrad must be a positive integer")
		return
	end
	mingrid 		= 	Float64(mingrid)
	maxgrid 		= 	Float64(maxgrid)

	if numrad 		== 	1
		grid 		= 	[maxgrid]
	elseif numrad 	== 	Inf
		grid 		= 	"all"
	else
		grid 		= 	range(mingrid,stop=maxgrid,length=numrad)
		grid 		= 	Array{Float64}(grid)
		grid[end] 	= 	maxgrid # this identity might not hold if we did not enforce it, due to numerical imprecision
	end
	return grid
end


##########################################################################################

####	BENCHMARKS

##########################################################################################

function roadmapbenchmarks(benchmarkdirectory)
	names 	= 	[
				"celegans_weighted_undirected_reindexed_for_matlab_maxdist_2.6429_SP_distmat.txt",
				"klein_bottle_pointcloud_new_400.txt_distmat.txt",
				"HIV1_2011.all.nt.concat.fa_hdm.txt",
				# "dragon_vrip.ply.txt_2000_.txt_distmat.txt",
				# "random_point_cloud_50_16_.txt_distmat.txt",
				"fractal_9_5_2_random_edge_list.txt_0.19795_distmat.txt"
				]
	dims  	=	Dict(
				"celegans_weighted_undirected_reindexed_for_matlab_maxdist_2.6429_SP_distmat.txt" => 1,
				"klein_bottle_pointcloud_new_400.txt_distmat.txt" =>1,
				"HIV1_2011.all.nt.concat.fa_hdm.txt" =>1,
				"dragon_vrip.ply.txt_2000_.txt_distmat.txt" =>1,
				"random_point_cloud_50_16_.txt_distmat.txt" =>7,
				"fractal_9_5_2_random_edge_list.txt_0.19795_distmat.txt"=>2
				)
	# precompile
	C 				= 	eirene(rand(10,10),model="pc")
	C 				= 	eirene(rand(10,10),model="pc",maxdim=3)
	for name in names
		println(name)
		fp 			= 	benchmarkdirectory*"/"*name
		x 			= 	readdlm(fp)
		x 			= 	(x+x')/2
		@time 	eirene(x,model="vr",maxdim=dims[name])
	end
end

function construction_sanitycheck(;numtrials = 10,samplesize = 50,sd=4)
	for i = 1:numtrials
		pcloud = rand(20,samplesize)
		d = Distances.pairwise(Euclidean(),pcloud,dims=2)
		(t,ocg2rad) = ordercanonicalform(d;fastop=false)
		construction_sanitycheck_subroutine(t,sd,samplesize)
		#gc()
	end
end

function construction_sanitycheck_subroutine(t,n,samplesize)
	N = n+2
	D2 = buildcomplex3(t,N,dictionaryoutput = true)
	D3 = buildcomplex3(t,n,dictionaryoutput = true)
	D2supp = trues(length(D2["farfaces"][n]))
	D2supp[D2["farfaces"][n+1][D2["prepairs"][n+1]]]=false
	checkthesearch = false
	for i = 1:samplesize
		D2ran = cran(D2["firstv"][n],i)
		D2ran = D2ran[D2supp[D2ran]]
		D2ff = D2["farfaces"][n][D2ran]
		D3ff = crows(D3["firstv"][n],D3["farfaces"][n],i)
		if sort(D2ff) != sort(D3ff)
			println("	sort(D2ff) != sort(D3ff) ")
			checkthesearch = true
			break
		end
		if checkthesearch == true
			break
		end
	end
	if sort(D2["farfaces"][n][D2["prepairs"][n]])!=sort(D3["farfaces"][n][D3["prepairs"][n]])
		println("please check prepairs")
		sleep(10)
	elseif checkthesearch
		println("checkthesearch evaluated to true")
		sleep(10)
	else
		println("ok so far as these checks are concerned :)")
	end
end

function cellcheck()
	c = 0
	maxdim 	= 	3
	for i = 1:20
		D 	= eirene(rand(20,50),model="pc",maxdim=maxdim)
		N 	= ff2complex(D["farfaces"],D["firstv"])
		Nf 	= ocff2of(D["grain"],D["ocg2rad"])
		Nrv = copy(N[1]);
		Ncp = copy(N[2]);
		Nf_copy = copy(Nf)
		F = persistf2complex(rv=N[1],cp=N[2],fv=Nf,maxdim=maxdim,record="all")
		if Nrv != N[1] || Ncp != N[2]
			print("changed N")
			break
		elseif Nf_copy != Nf
			print("changed Nf")
			break
		end
		for k = 1:maxdim
			if sortslices(barcode(D,dim=k),dims=1)!= sortslices(barcode(F,dim=k),dims=1)
				c+=1
			end
		end
	end
	return c
end


function persistencemultistats(x;numtrials = 10)
	trials = Array{Any,1}(numtrials)
	for i = 1:numtrials
		trials[i] = persistencestats(x)
	end
	return trials
end
